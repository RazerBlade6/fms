//! The ocean: a two-layer slab after Codron (2012) and Bhatnagar et al.
//! (2026), with sea ice and wind-driven transport to come (M2 steps 11–12).
//!
//! Each cell with any ocean has a mixed layer, which exchanges heat with the
//! atmosphere, over a deep layer, which exchanges heat only with the mixed
//! layer and its neighbours. Neither moves yet. Heat travels by:
//!
//! - **Horizontal diffusion** in both layers, as a flux across every edge
//!   between two cells that have ocean, so the ocean's heat is conserved to
//!   round-off and none crosses a coast. An edge's opening is the smaller of
//!   its two cells' ocean fractions, which also keeps a nearly dry cell from
//!   taking more heat than it can hold.
//! - **Convective adjustment**: a mixed layer colder than the deep layer
//!   below it overturns, and both take their mass-weighted mean.
//!
//! The diffusivity is the one the paper's runs actually have. The paper gives
//! 8000 m² s⁻¹, but its code scales that by 50 m over the ocean's total depth
//! (`ocean_slab_mod.F90`), which is 2000 m² s⁻¹ for 50 m over 150 m.
//!
//! The crate knows nothing about the atmosphere: fluxes come in as numbers and
//! temperatures go out.

use fms_grid::Operators;
use fms_grid::operators::MIN_CHUNK;
use rayon::prelude::*;

/// The ocean's layers and how heat moves through them. [`Default`] is the
/// two-layer ocean of Bhatnagar et al. (2026, Table A1).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OceanParams {
    /// Density of seawater (kg m⁻³).
    pub density: f64,
    /// Specific heat of seawater (J kg⁻¹ K⁻¹).
    pub specific_heat: f64,
    /// Depth of the mixed layer (m).
    pub mixed_depth: f64,
    /// Depth of the deep layer (m); 0 for none.
    pub deep_depth: f64,
    /// Horizontal diffusivity in both layers (m² s⁻¹); 0 for none.
    pub diffusivity: f64,
    /// Whether a mixed layer colder than the deep layer overturns.
    pub convective_adjustment: bool,
    /// The deep layer starts at the mixed layer's temperature, but no warmer
    /// than this (K). Without upwelling (M2 step 12) the deep layer can only
    /// lose heat, so a warmer start would come back as a false warm spell.
    pub deep_start_limit: f64,
}

impl Default for OceanParams {
    fn default() -> Self {
        Self {
            density: 1026.0,
            specific_heat: 3994.0,
            mixed_depth: 50.0,
            deep_depth: 150.0,
            diffusivity: 2000.0,
            convective_adjustment: true,
            deep_start_limit: 277.15,
        }
    }
}

impl OceanParams {
    /// A motionless slab of fresh water `depth` metres deep, with no deep
    /// layer: the ocean before M2 step 11, kept for short spin-ups, tuning
    /// runs and the tests pinned to it.
    pub fn slab(depth: f64) -> Self {
        Self {
            density: 1000.0,
            specific_heat: 4186.0,
            mixed_depth: depth,
            deep_depth: 0.0,
            diffusivity: 0.0,
            convective_adjustment: false,
            deep_start_limit: f64::INFINITY,
        }
    }

    /// Heat capacity of the mixed layer (J m⁻² K⁻¹).
    pub fn mixed_capacity(&self) -> f64 {
        self.density * self.specific_heat * self.mixed_depth
    }

    /// Heat capacity of the deep layer (J m⁻² K⁻¹).
    pub fn deep_capacity(&self) -> f64 {
        self.density * self.specific_heat * self.deep_depth
    }
}

/// The ocean's temperatures at one moment, for [`Ocean::jump`].
#[derive(Clone, Debug)]
pub struct OceanSnapshot {
    mixed: Vec<f64>,
    deep: Vec<f64>,
}

/// An edge between two cells that both have ocean.
#[derive(Clone, Copy, Debug)]
struct Link {
    cells: [u32; 2],
    /// Edge length over the distance between the cell centres, times the
    /// edge's opening.
    conductance: f64,
}

/// Largest share of the difference from its neighbours that a cell may close
/// in one diffusion step. Explicit diffusion is monotonic up to 1.
const MAX_DIFFUSION_NUMBER: f64 = 0.5;

/// The ocean's state on one mesh.
#[derive(Clone, Debug)]
pub struct Ocean {
    params: OceanParams,
    /// Temperature of each cell's mixed layer (K); stays at its start where
    /// the cell has no ocean.
    mixed: Vec<f64>,
    /// Temperature of each cell's deep layer (K); likewise.
    deep: Vec<f64>,
    /// Ocean area of each cell (m²).
    area: Vec<f64>,
    links: Vec<Link>,
    /// Largest `Σ conductance / area` of any cell (m⁻²).
    largest_rate: f64,
    /// Heat convergence per cell, reused between steps.
    convergence: Vec<f64>,
}

impl Ocean {
    /// An ocean on the mesh of `ops`, covering what `land_fraction` leaves
    /// of each cell, with its mixed layer at `surface` (K).
    ///
    /// # Panics
    /// If the lengths differ from the mesh's cells, a land fraction is
    /// outside 0–1, or the mixed layer has no depth.
    pub fn new(
        params: OceanParams,
        ops: &Operators,
        land_fraction: &[f64],
        surface: Vec<f64>,
    ) -> Self {
        let (mesh, geometry) = (ops.mesh(), ops.geometry());
        let n = mesh.n_cells();
        assert_eq!(land_fraction.len(), n);
        assert_eq!(surface.len(), n);
        assert!(
            land_fraction.iter().all(|f| (0.0..=1.0).contains(f)),
            "land fractions must be within 0–1"
        );
        assert!(params.mixed_depth > 0.0, "the mixed layer needs a depth");
        let r2 = ops.radius() * ops.radius();
        let open: Vec<f64> = land_fraction.iter().map(|f| 1.0 - f).collect();
        let area: Vec<f64> = (0..n)
            .map(|i| open[i] * geometry.cell_areas[i] * r2)
            .collect();
        let links: Vec<Link> = mesh
            .edges()
            .iter()
            .enumerate()
            .filter_map(|(e, edge)| {
                let opening = open[edge.cells[0] as usize].min(open[edge.cells[1] as usize]);
                (opening > 0.0).then(|| Link {
                    cells: edge.cells,
                    conductance: opening * geometry.edge_lengths[e] / geometry.dual_edge_lengths[e],
                })
            })
            .collect();
        let mut sum = vec![0.0; n];
        for link in &links {
            for c in link.cells {
                sum[c as usize] += link.conductance;
            }
        }
        let largest_rate = sum
            .iter()
            .zip(&area)
            .filter(|(_, a)| **a > 0.0)
            .map(|(s, a)| s / a)
            .fold(0.0, f64::max);
        let deep = surface
            .iter()
            .map(|t| t.min(params.deep_start_limit))
            .collect();
        Self {
            params,
            mixed: surface,
            deep,
            area,
            links,
            largest_rate,
            convergence: vec![0.0; n],
        }
    }

    /// The parameters.
    pub fn params(&self) -> &OceanParams {
        &self.params
    }

    /// Temperature of each cell's mixed layer (K); meaningless where the
    /// cell has no ocean.
    pub fn mixed_temperature(&self) -> &[f64] {
        &self.mixed
    }

    /// Temperature of each cell's deep layer (K); meaningless where the cell
    /// has no ocean or the ocean has no deep layer.
    pub fn deep_temperature(&self) -> &[f64] {
        &self.deep
    }

    /// Heat held in cell `i`'s two layers (J per m² of its ocean).
    pub fn heat(&self, i: usize) -> f64 {
        self.params.mixed_capacity() * self.mixed[i] + self.params.deep_capacity() * self.deep[i]
    }

    /// Heat the mixed layers gain over `dt` seconds from `flux(i)`, the net
    /// flux into cell `i`'s ocean (W per m² of it).
    pub fn add_surface_heat(&mut self, dt: f64, flux: impl Fn(usize) -> f64 + Sync) {
        let capacity = self.params.mixed_capacity();
        let area = &self.area;
        self.mixed
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .filter(|(i, _)| area[*i] > 0.0)
            .for_each(|(i, t)| *t += dt * flux(i) / capacity);
    }

    /// Move heat for `dt` seconds: diffusion in both layers, then convective
    /// adjustment.
    pub fn transport(&mut self, dt: f64) {
        if self.params.diffusivity > 0.0 {
            let number = self.params.diffusivity * dt * self.largest_rate;
            let steps = (number / MAX_DIFFUSION_NUMBER).ceil().max(1.0);
            let amount = self.params.diffusivity * dt / steps;
            for _ in 0..steps as usize {
                self.diffuse_mixed(amount);
                if self.params.deep_depth > 0.0 {
                    self.diffuse_deep(amount);
                }
            }
        }
        if self.params.convective_adjustment && self.params.deep_depth > 0.0 {
            self.adjust();
        }
    }

    fn diffuse_mixed(&mut self, amount: f64) {
        diffuse(
            &mut self.mixed,
            &self.links,
            &self.area,
            &mut self.convergence,
            amount,
        );
    }

    fn diffuse_deep(&mut self, amount: f64) {
        diffuse(
            &mut self.deep,
            &self.links,
            &self.area,
            &mut self.convergence,
            amount,
        );
    }

    /// The temperatures now, for a later [`jump`](Self::jump).
    pub fn snapshot(&self) -> OceanSnapshot {
        OceanSnapshot {
            mixed: self.mixed.clone(),
            deep: self.deep.clone(),
        }
    }

    /// Spin-up only: repeat each layer's change since `from` for `periods`
    /// more of the same span, moving no temperature by more than `limit` (K),
    /// then overturn what that left unstable. Taken over a whole year, the
    /// change has no seasonal part. It doesn't conserve heat, by design.
    pub fn jump(&mut self, from: &OceanSnapshot, periods: f64, limit: f64) {
        let extrapolate = |now: &mut [f64], then: &[f64]| {
            for ((t, t0), a) in now.iter_mut().zip(then).zip(&self.area) {
                if *a > 0.0 {
                    *t += (periods * (*t - t0)).clamp(-limit, limit);
                }
            }
        };
        extrapolate(&mut self.mixed, &from.mixed);
        if self.params.deep_depth > 0.0 {
            extrapolate(&mut self.deep, &from.deep);
            if self.params.convective_adjustment {
                self.adjust();
            }
        }
    }

    /// Mix each column whose mixed layer is colder than its deep layer.
    fn adjust(&mut self) {
        let (hm, hd) = (self.params.mixed_depth, self.params.deep_depth);
        self.mixed
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(self.deep.par_iter_mut())
            .zip(self.area.par_iter())
            .filter(|((mixed, deep), area)| **area > 0.0 && **deep > **mixed)
            .for_each(|((mixed, deep), _)| {
                let mean = (hm * *mixed + hd * *deep) / (hm + hd);
                *mixed = mean;
                *deep = mean;
            });
    }
}

/// One explicit diffusion step of `field`, where `amount` is the diffusivity
/// times the step (m²). Flux form, so `Σ area · field` is conserved.
fn diffuse(field: &mut [f64], links: &[Link], area: &[f64], convergence: &mut [f64], amount: f64) {
    convergence.fill(0.0);
    for link in links {
        let [a, b] = link.cells.map(|c| c as usize);
        let flux = link.conductance * (field[b] - field[a]);
        convergence[a] += flux;
        convergence[b] -= flux;
    }
    for ((t, c), a) in field.iter_mut().zip(convergence.iter()).zip(area) {
        if *a > 0.0 {
            *t += amount * c / a;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fms_grid::{LloydOptions, MeshGeometry, SphereMesh};

    const RADIUS: f64 = 6.371e6;

    struct Grid {
        mesh: SphereMesh,
        geometry: MeshGeometry,
    }

    impl Grid {
        fn new(k: u32) -> Self {
            let (mesh, _) = SphereMesh::centroidal(k, LloydOptions::default());
            let geometry = MeshGeometry::new(&mesh);
            Self { mesh, geometry }
        }

        fn ops(&self) -> Operators<'_> {
            Operators::new(&self.mesh, &self.geometry, RADIUS)
        }

        fn n(&self) -> usize {
            self.mesh.n_cells()
        }

        fn each(&self, f: impl Fn(f64) -> f64) -> Vec<f64> {
            self.mesh.cell_centers().iter().map(|p| f(p.z)).collect()
        }
    }

    /// A number in 0–1 that looks random but is fixed by `i`.
    fn scatter(i: usize, seed: u64) -> f64 {
        let x = (i as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15 ^ seed);
        (x >> 40) as f64 / (1u64 << 24) as f64
    }

    /// Land fractions with whole-land, whole-ocean and coastal cells.
    fn coasts(n: usize) -> Vec<f64> {
        (0..n)
            .map(|i| match i % 5 {
                0 => 0.0,
                1 => 1.0,
                _ => scatter(i, 1),
            })
            .collect()
    }

    /// Total heat of the ocean (J).
    fn total(ocean: &Ocean) -> f64 {
        (0..ocean.area.len())
            .map(|i| ocean.area[i] * ocean.heat(i))
            .sum()
    }

    #[test]
    fn the_defaults_are_the_papers_ocean() {
        let p = OceanParams::default();
        assert_eq!(p.mixed_capacity(), 1026.0 * 3994.0 * 50.0);
        assert_eq!(p.deep_capacity(), 1026.0 * 3994.0 * 150.0);
        // The paper's 8000 m²/s as its code applies it.
        assert_eq!(p.diffusivity, 8000.0 * 50.0 / (50.0 + 150.0));
    }

    #[test]
    fn a_slab_only_takes_up_surface_heat() {
        let grid = Grid::new(8);
        let start = grid.each(|z| 280.0 + 10.0 * z);
        let params = OceanParams::slab(10.0);
        let land = vec![0.0; grid.n()];
        let mut ocean = Ocean::new(params, &grid.ops(), &land, start.clone());
        ocean.transport(3600.0);
        assert_eq!(ocean.mixed_temperature(), &start[..]);
        ocean.add_surface_heat(1800.0, |i| i as f64);
        for (i, (t, s)) in ocean.mixed_temperature().iter().zip(&start).enumerate() {
            assert_eq!(*t, s + 1800.0 * i as f64 / params.mixed_capacity());
        }
    }

    #[test]
    fn the_deep_layer_starts_no_warmer_than_its_limit() {
        let grid = Grid::new(8);
        let start = grid.each(|z| 288.0 - 20.0 * z * z);
        let land = vec![0.0; grid.n()];
        let ocean = Ocean::new(OceanParams::default(), &grid.ops(), &land, start.clone());
        for (d, s) in ocean.deep_temperature().iter().zip(&start) {
            assert_eq!(*d, s.min(277.15));
        }
    }

    #[test]
    fn surface_heat_goes_only_into_ocean_and_is_counted() {
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        let start = grid.each(|z| 285.0 + 12.0 * z);
        let mut ocean = Ocean::new(OceanParams::default(), &grid.ops(), &land, start.clone());
        let before = total(&ocean);
        let (dt, flux) = (1800.0, |i: usize| 300.0 * (scatter(i, 2) - 0.5));
        ocean.add_surface_heat(dt, flux);
        let expected: f64 = (0..grid.n()).map(|i| ocean.area[i] * dt * flux(i)).sum();
        assert!(
            (total(&ocean) - before - expected).abs() < 1e-14 * before,
            "gained {}, expected {expected}",
            total(&ocean) - before
        );
        for i in 0..grid.n() {
            if land[i] == 1.0 {
                assert_eq!(ocean.mixed_temperature()[i], start[i]);
            }
        }
    }

    #[test]
    fn diffusion_conserves_heat_to_round_off_and_stays_within_its_range() {
        // Rough temperatures and rough coasts, at a normal step and at one
        // long enough to need sub-steps.
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        let start: Vec<f64> = (0..grid.n())
            .map(|i| 271.0 + 30.0 * scatter(i, 3))
            .collect();
        let params = OceanParams {
            convective_adjustment: false,
            ..OceanParams::default()
        };
        for dt in [1800.0, 3.0e10] {
            let mut ocean = Ocean::new(params, &grid.ops(), &land, start.clone());
            let open = |t: &[f64]| -> (f64, f64) {
                t.iter()
                    .zip(&land)
                    .filter(|(_, f)| **f < 1.0)
                    .fold((f64::MAX, f64::MIN), |(a, b), (t, _)| {
                        (a.min(*t), b.max(*t))
                    })
            };
            let (low, high) = open(&start);
            let before = total(&ocean);
            for _ in 0..20 {
                ocean.transport(dt);
            }
            let after = total(&ocean);
            assert!(
                ((after - before) / before).abs() < 1e-14,
                "relative drift {} at dt {dt}",
                (after - before) / before
            );
            let (new_low, new_high) = open(ocean.mixed_temperature());
            assert!(new_low >= low && new_high <= high, "{new_low}–{new_high}");
        }
    }

    #[test]
    fn a_uniform_ocean_stays_exactly_uniform() {
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        let start = vec![275.0; grid.n()];
        let mut ocean = Ocean::new(OceanParams::default(), &grid.ops(), &land, start.clone());
        ocean.transport(1800.0);
        assert_eq!(ocean.mixed_temperature(), &start[..]);
        assert_eq!(ocean.deep_temperature(), &start[..]);
    }

    #[test]
    fn diffusion_damps_a_zonal_wave_at_the_analytic_rate() {
        // P₂(z) = (3z² − 1)/2 has ∇²P₂ = −6 P₂ / R².
        let grid = Grid::new(32);
        let land = vec![0.0; grid.n()];
        let p2 = grid.each(|z| 0.5 * (3.0 * z * z - 1.0));
        let params = OceanParams {
            convective_adjustment: false,
            deep_start_limit: f64::INFINITY,
            ..OceanParams::default()
        };
        let start: Vec<f64> = p2.iter().map(|x| 20.0 * x).collect();
        let mut ocean = Ocean::new(params, &grid.ops(), &land, start);
        let time = 1.0e9;
        let steps = 200;
        for _ in 0..steps {
            ocean.transport(time / f64::from(steps));
        }
        let project = |t: &[f64]| -> f64 {
            let a = &grid.geometry.cell_areas;
            let top: f64 = (0..t.len()).map(|i| a[i] * t[i] * p2[i]).sum();
            let bottom: f64 = (0..t.len()).map(|i| a[i] * p2[i] * p2[i]).sum();
            top / bottom
        };
        let exact = 20.0 * fms_math::exp(-6.0 * params.diffusivity * time / (RADIUS * RADIUS));
        for (name, t) in [
            ("mixed", ocean.mixed_temperature()),
            ("deep", ocean.deep_temperature()),
        ] {
            let amplitude = project(t);
            // 20 K decays to about 18.3 K; the decay itself is right to 2%.
            assert!(
                ((20.0 - amplitude) / (20.0 - exact) - 1.0).abs() < 0.02,
                "{name}: amplitude {amplitude}, exact {exact}"
            );
        }
    }

    #[test]
    fn no_heat_crosses_a_continent() {
        // Land around the equator splits the ocean in two; each side keeps
        // its own heat.
        let grid = Grid::new(16);
        let land = grid.each(|z| if z.abs() < 0.2 { 1.0 } else { 0.0 });
        let start = grid.each(|z| 280.0 + 15.0 * z);
        let mut ocean = Ocean::new(OceanParams::default(), &grid.ops(), &land, start);
        let side = |o: &Ocean, north: bool| -> f64 {
            let centers = grid.mesh.cell_centers();
            (0..grid.n())
                .filter(|&i| (centers[i].z > 0.0) == north)
                .map(|i| o.area[i] * o.heat(i))
                .sum()
        };
        let (north, south) = (side(&ocean, true), side(&ocean, false));
        for _ in 0..100 {
            ocean.transport(86_400.0);
        }
        assert!((side(&ocean, true) / north - 1.0).abs() < 1e-14);
        assert!((side(&ocean, false) / south - 1.0).abs() < 1e-14);
    }

    #[test]
    fn a_cold_surface_overturns_and_a_warm_one_does_not() {
        let grid = Grid::new(8);
        let land = vec![0.0; grid.n()];
        let params = OceanParams {
            diffusivity: 0.0,
            deep_start_limit: 280.0,
            ..OceanParams::default()
        };
        // The deep layer starts at 280 K everywhere. The north's surface is
        // warmer than that; the south's is chilled 8 K below it.
        let start = grid.each(|z| if z < 0.0 { 280.0 } else { 290.0 });
        let mut ocean = Ocean::new(params, &grid.ops(), &land, start);
        let centers = grid.mesh.cell_centers();
        let cooling = 8.0 * params.mixed_capacity();
        ocean.add_surface_heat(1.0, |i| if centers[i].z < 0.0 { -cooling } else { 0.0 });
        let before = total(&ocean);
        ocean.transport(1800.0);
        for (i, center) in centers.iter().enumerate() {
            let (m, d) = (ocean.mixed_temperature()[i], ocean.deep_temperature()[i]);
            if center.z < 0.0 {
                // A quarter of the water at 272 K, the rest at 280 K.
                assert_eq!(m, d);
                assert!((m - 278.0).abs() < 1e-9, "{m}");
            } else {
                assert_eq!(m, 290.0);
                assert_eq!(d, 280.0);
            }
        }
        assert!(((total(&ocean) - before) / before).abs() < 1e-14);
    }

    #[test]
    fn a_jump_repeats_the_change_within_its_limit() {
        let grid = Grid::new(4);
        let mut land = vec![0.0; grid.n()];
        land[0] = 1.0;
        let params = OceanParams {
            diffusivity: 0.0,
            ..OceanParams::default()
        };
        let mut ocean = Ocean::new(params, &grid.ops(), &land, vec![290.0; grid.n()]);
        let start = ocean.snapshot();
        let c = params.mixed_capacity();
        ocean.add_surface_heat(1.0, |i| if i == 1 { 2.0 * c } else { 0.5 * c });
        ocean.jump(&start, 3.0, 3.0);
        let mixed = ocean.mixed_temperature();
        assert_eq!(mixed[0], 290.0, "no ocean here");
        assert_eq!(mixed[1], 295.0, "2 K, plus 6 K limited to 3 K");
        assert_eq!(mixed[2], 292.0, "0.5 K, plus 3 × 0.5 K");
    }

    #[test]
    fn a_jump_overturns_what_it_leaves_unstable() {
        let grid = Grid::new(4);
        let land = vec![0.0; grid.n()];
        let params = OceanParams {
            diffusivity: 0.0,
            ..OceanParams::default()
        };
        let mut ocean = Ocean::new(params, &grid.ops(), &land, vec![277.0; grid.n()]);
        let start = ocean.snapshot();
        ocean.add_surface_heat(1.0, |_| -params.mixed_capacity());
        ocean.jump(&start, 1.0, 10.0);
        // A quarter of the water at 275 K, the rest at 277 K.
        assert_eq!(ocean.mixed_temperature()[0], 276.5);
        assert_eq!(ocean.deep_temperature()[0], 276.5);
    }

    #[test]
    fn overturning_releases_the_deep_layers_heat() {
        let grid = Grid::new(8);
        let land = vec![0.0; grid.n()];
        let params = OceanParams {
            diffusivity: 0.0,
            ..OceanParams::default()
        };
        let mut ocean = Ocean::new(params, &grid.ops(), &land, vec![277.0; grid.n()]);
        let before = total(&ocean);
        // Strong cooling chills the mixed layer below the deep layer's 277 K.
        ocean.add_surface_heat(86_400.0, |_| -400.0);
        let cooled = params.mixed_capacity();
        ocean.transport(0.0);
        let drop = 86_400.0 * 400.0 / cooled;
        let expected = 277.0 - 0.25 * drop;
        for (m, d) in ocean
            .mixed_temperature()
            .iter()
            .zip(ocean.deep_temperature())
        {
            assert!((m - expected).abs() < 1e-9 && m == d, "{m}, {d}");
        }
        let removed = 86_400.0 * 400.0 * ocean.area.iter().sum::<f64>();
        assert!(((total(&ocean) - (before - removed)) / before).abs() < 1e-14);
    }
}
