//! Grid imprinting: does a long time mean carry a pattern locked to the
//! icosahedral grid?
//!
//! Under zonally symmetric forcing without terrain (Held–Suarez), the time
//! mean of any field should depend on latitude only. What is left after
//! removing the zonal mean is either sampling noise from the weather, which
//! differs between independent stretches of the run, or a pattern fixed to
//! the grid, which does not. [`GridImprint`] samples a run into two halves of
//! alternating blocks; [`ImprintAnalysis`] splits each field's residual into
//! what the halves share (the persistent pattern) and what they don't (the
//! noise).
//!
//! A pattern fixed to the grid must also share the grid's symmetry: the
//! icosahedral mesh repeats under a half-turn about the axis (`R`,
//! longitude + 180°) and a mirror in the equator (`M`). The residual splits
//! into four parts, even or odd under each; weather the halves happen to
//! share spreads over all four, a grid-locked pattern sits in the part even
//! under both. The analysis also reports the zonal wavenumbers and the 12
//! pentagons.
//!
//! The analysis is a pure function of the mesh and the halves' means, so it
//! works on synthetic fields as well as on a run.

use fms_grid::coords::to_lat_lon;
use fms_grid::zonal::east_north;
use fms_grid::{CellLocator, LatitudeBands, MeshGeometry, SphereMesh};
use glam::DVec3;

use crate::layered::{LayeredModel, LayeredState};

/// Zonal wavenumbers analysed: `1..=WAVENUMBERS`.
pub const WAVENUMBERS: usize = 12;
/// Latitude belts (15° each) the wavenumbers are measured in.
pub const BELTS: usize = 12;
/// Half-width of the local zonal-mean fit, in mean grid spacings.
pub const FIT_SPACINGS: f64 = 2.5;

/// The fields sampled, as `(name, units, mirror)`, in the order of
/// [`GridImprint::results`]: surface pressure, the lowest layer's
/// temperature, the eastward wind of the upper layer (the one holding
/// 250 hPa) and the northward wind of the lowest layer. `mirror` is the
/// field's sign under a mirror in the equator: −1 for the northward wind.
pub const FIELDS: [(&str, &str, f64); 4] = [
    ("surface_pressure", "Pa", 1.0),
    ("temperature_lower", "K", 1.0),
    ("east_wind_upper", "m/s", 1.0),
    ("north_wind_lower", "m/s", -1.0),
];

/// The imprint measures of one field. Amplitudes are in the field's units;
/// every mean is area-weighted over the sphere unless it says otherwise.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Imprint {
    /// The field's name.
    pub name: String,
    /// Its units.
    pub units: String,
    /// RMS of the full time mean's residual (the mean minus its local
    /// zonal mean, [`ImprintAnalysis::residual`]).
    pub rms: f64,
    /// RMS of the band means of that residual: the part that is a function
    /// of latitude, which the zonal-mean fit should leave near zero.
    pub zonal: f64,
    /// RMS of half the difference between the halves' residuals: the
    /// sampling noise left in the full mean.
    pub noise: f64,
    /// Mean of the product of the halves' residuals.
    pub covariance: f64,
    /// `√max(covariance, 0)`: the amplitude of the pattern both halves share.
    pub persistent: f64,
    /// Correlation between the halves' residuals.
    pub correlation: f64,
    /// The cross-half covariance of each symmetry part of the residual
    /// ([`ImprintAnalysis::symmetry_parts`]): even under the half-turn and
    /// the mirror, even and odd, odd and even, odd under both. They sum to
    /// `covariance`.
    pub symmetric: [f64; 4],
    /// `√max(symmetric[0] − mean(symmetric[1..]), 0)`: the amplitude of the
    /// shared pattern with the grid's symmetry, beyond what the other parts
    /// share by chance.
    pub locked: f64,
    /// Mean over cells of the standard deviation in time of every sample.
    pub transient: f64,
    /// `√max(·, 0)` of the plain mean over the 12 pentagon cells of the
    /// product of the halves' residuals.
    pub pentagon: f64,
    /// The same over the pentagons and their neighbours.
    pub ring: f64,
    /// Signed cross power between the halves of zonal wavenumbers
    /// `1..=12` (`waves[m − 1]`): a pattern `a cos mλ` in both contributes
    /// its variance `a²/2`; noise averages out.
    pub waves: [f64; WAVENUMBERS],
    /// The full time mean's residual, per cell.
    pub residual: Vec<f64>,
    /// Its part even under both symmetries, per cell.
    pub locked_residual: Vec<f64>,
}

/// The mesh-dependent parts of the imprint analysis: the local zonal-mean
/// fit, the grid's symmetries, latitude bands and belts, and the pentagons.
#[derive(Clone, Debug)]
pub struct ImprintAnalysis {
    bands: LatitudeBands,
    /// Cell areas (unit sphere).
    areas: Vec<f64>,
    /// Each cell's share of the sphere's area.
    weights: Vec<f64>,
    /// Cell-centre latitude and longitude (radians).
    latitude: Vec<f64>,
    longitude: Vec<f64>,
    /// Cells in order of latitude, and their latitudes.
    by_latitude: Vec<usize>,
    sorted_latitude: Vec<f64>,
    /// Half-width of the zonal-mean fit (radians).
    window: f64,
    /// The cell at each cell's image under the half-turn and the mirror.
    turn: Vec<usize>,
    mirror: Vec<usize>,
    /// The largest distance (radians) from a cell centre's image to the
    /// centre of the cell it maps to.
    symmetry_error: f64,
    belt_of_cell: Vec<usize>,
    /// Each belt's share of the sphere's area.
    belt_weight: [f64; BELTS],
    pentagons: Vec<usize>,
    /// The pentagons and their neighbours.
    ring: Vec<usize>,
}

impl ImprintAnalysis {
    /// The analysis on `mesh`; `bands` equal-width latitude bands measure
    /// what latitude-only residual is left ([`Imprint::zonal`]).
    pub fn new(mesh: &SphereMesh, geometry: &MeshGeometry, bands: usize) -> Self {
        let n = mesh.n_cells();
        let bands = LatitudeBands::new(mesh, geometry, bands);
        let areas = geometry.cell_areas.clone();
        let total: f64 = areas.iter().sum();
        let weights: Vec<f64> = areas.iter().map(|a| a / total).collect();
        let (latitude, longitude): (Vec<f64>, Vec<f64>) =
            mesh.cell_centers().iter().map(|&p| to_lat_lon(p)).unzip();
        let mut by_latitude: Vec<usize> = (0..n).collect();
        by_latitude.sort_by(|&a, &b| latitude[a].total_cmp(&latitude[b]).then(a.cmp(&b)));
        let sorted_latitude = by_latitude.iter().map(|&c| latitude[c]).collect();
        let window = FIT_SPACINGS * (4.0 * std::f64::consts::PI / n as f64).sqrt();
        // The cell holding each centre's image.
        let locator = CellLocator::new(mesh);
        let centers = mesh.cell_centers();
        let mut symmetry_error = 0.0f64;
        let mut image = |f: fn(DVec3) -> DVec3| -> Vec<usize> {
            centers
                .iter()
                .map(|&p| {
                    let q = f(p);
                    let c = locator.locate(mesh, q) as usize;
                    let d = fms_math::acos(q.dot(centers[c]).clamp(-1.0, 1.0));
                    symmetry_error = symmetry_error.max(d);
                    c
                })
                .collect()
        };
        let turn = image(|p| DVec3::new(-p.x, -p.y, p.z));
        let mirror = image(|p| DVec3::new(p.x, p.y, -p.z));
        let width = std::f64::consts::PI / BELTS as f64;
        let belt_of_cell: Vec<usize> = latitude
            .iter()
            .map(|&lat| (((lat + std::f64::consts::FRAC_PI_2) / width) as usize).min(BELTS - 1))
            .collect();
        let mut belt_weight = [0.0; BELTS];
        for (&b, &w) in belt_of_cell.iter().zip(&weights) {
            belt_weight[b] += w;
        }
        let pentagons: Vec<usize> = (0..n)
            .filter(|&c| mesh.cell_neighbors(c as u32).len() == 5)
            .collect();
        let mut ring: Vec<usize> = pentagons
            .iter()
            .flat_map(|&c| {
                std::iter::once(c).chain(mesh.cell_neighbors(c as u32).iter().map(|&n| n as usize))
            })
            .collect();
        ring.sort_unstable();
        ring.dedup();
        Self {
            bands,
            areas,
            weights,
            latitude,
            longitude,
            by_latitude,
            sorted_latitude,
            window,
            turn,
            mirror,
            symmetry_error,
            belt_of_cell,
            belt_weight,
            pentagons,
            ring,
        }
    }

    /// The cells with five neighbours (12 on an icosahedral mesh).
    pub fn pentagons(&self) -> &[usize] {
        &self.pentagons
    }

    /// The pentagons and their neighbours, in cell order.
    pub fn ring(&self) -> &[usize] {
        &self.ring
    }

    /// How far the mesh is from repeating under the half-turn and the
    /// mirror: the largest distance (radians) from a cell centre's image to
    /// the centre of the cell holding it. Near zero on an icosahedral mesh;
    /// the symmetry parts mean nothing otherwise.
    pub fn symmetry_error(&self) -> f64 {
        self.symmetry_error
    }

    /// Whether the half-turn and the mirror each map the cells one to one.
    pub fn symmetries_are_permutations(&self) -> bool {
        [&self.turn, &self.mirror].iter().all(|map| {
            let mut hit = vec![false; map.len()];
            map.iter().all(|&c| !std::mem::replace(&mut hit[c], true))
        })
    }

    /// Half-width of the zonal-mean fit (radians).
    pub fn window(&self) -> f64 {
        self.window
    }

    /// `values` minus its zonal mean at each cell's latitude.
    ///
    /// The zonal mean at cell `i` is the constant term of a weighted
    /// least-squares quadratic in `φ − φᵢ`, fitted to every cell within the
    /// window `w` ([`FIT_SPACINGS`] mean grid spacings) of `φᵢ`, with
    /// weights `Aⱼ (1 − |φⱼ − φᵢ| / w)`. The quadratic follows a profile's
    /// curvature, which a piecewise-linear zonal mean leaves behind as a
    /// residual of latitude only. If the fit is singular it falls back to a
    /// straight line, then to the weighted mean.
    pub fn residual(&self, values: &[f64]) -> Vec<f64> {
        assert_eq!(values.len(), self.areas.len(), "one value per cell");
        let w = self.window;
        (0..values.len())
            .map(|i| {
                let phi = self.latitude[i];
                let lo = self.sorted_latitude.partition_point(|&l| l <= phi - w);
                let hi = self.sorted_latitude.partition_point(|&l| l < phi + w);
                // Moments Σ k xᵏ (x = Δφ / w) and Σ k xᵏ v.
                let (mut s, mut t) = ([0.0f64; 5], [0.0f64; 3]);
                for &j in &self.by_latitude[lo..hi] {
                    let x = (self.latitude[j] - phi) / w;
                    let k = self.areas[j] * (1.0 - x.abs());
                    let (x2, v) = (x * x, values[j]);
                    s[0] += k;
                    s[1] += k * x;
                    s[2] += k * x2;
                    s[3] += k * x2 * x;
                    s[4] += k * x2 * x2;
                    t[0] += k * v;
                    t[1] += k * x * v;
                    t[2] += k * x2 * v;
                }
                values[i] - fit_constant(&s, &t)
            })
            .collect()
    }

    /// Area-weighted mean of a cell field over each band.
    fn band_means(&self, values: &[f64]) -> Vec<f64> {
        assert_eq!(values.len(), self.areas.len(), "one value per cell");
        let mut sums = vec![0.0; self.bands.len()];
        for ((&b, &a), &v) in self
            .bands
            .band_of_cell()
            .iter()
            .zip(&self.areas)
            .zip(values)
        {
            sums[b as usize] += a * v;
        }
        sums.iter()
            .zip(self.bands.areas())
            .map(|(s, a)| s / a)
            .collect()
    }

    /// The four symmetry parts of a cell field `e` whose sign under the
    /// mirror is `mirror` (±1): with `eR`, `eM`, `eRM` its values at the
    /// image cells and `σ = mirror`,
    /// `[(e + eR + σeM + σeRM), (e + eR − σeM − σeRM), (e − eR + σeM − σeRM),
    /// (e − eR − σeM + σeRM)] / 4`: even under the half-turn and the
    /// mirror, even and odd, odd and even, odd under both. They sum to `e`.
    pub fn symmetry_parts(&self, e: &[f64], mirror: f64) -> [Vec<f64>; 4] {
        assert_eq!(e.len(), self.areas.len(), "one value per cell");
        let mut parts: [Vec<f64>; 4] = std::array::from_fn(|_| Vec::with_capacity(e.len()));
        for i in 0..e.len() {
            let (m, rm) = (self.mirror[i], self.turn[self.mirror[i]]);
            let (a, r) = (e[i], e[self.turn[i]]);
            let (sm, srm) = (mirror * e[m], mirror * e[rm]);
            parts[0].push(0.25 * (a + r + sm + srm));
            parts[1].push(0.25 * (a + r - sm - srm));
            parts[2].push(0.25 * (a - r + sm - srm));
            parts[3].push(0.25 * (a - r - sm + srm));
        }
        parts
    }

    /// The imprint measures of a field `(name, units, mirror)` (as in
    /// [`FIELDS`]) from its time means over half A and half B (`counts`
    /// samples each, which weight them in the full mean) and its standard
    /// deviation in time per cell.
    pub fn analyze(
        &self,
        (name, units, mirror): (&str, &str, f64),
        means: [&[f64]; 2],
        counts: [u64; 2],
        transient: &[f64],
    ) -> Imprint {
        let n = self.areas.len();
        assert!(
            means.iter().all(|m| m.len() == n) && transient.len() == n,
            "one value per cell"
        );
        let [ra, rb] = means.map(|m| self.residual(m));
        let total = (counts[0] + counts[1]).max(1) as f64;
        let (wa, wb) = (counts[0] as f64 / total, counts[1] as f64 / total);
        let residual: Vec<f64> = ra.iter().zip(&rb).map(|(a, b)| wa * a + wb * b).collect();
        let mean = |f: &dyn Fn(usize) -> f64| (0..n).map(|c| self.weights[c] * f(c)).sum::<f64>();
        let rms = mean(&|c| residual[c] * residual[c]).sqrt();
        let band_means = self.band_means(&residual);
        let band = self.bands.band_of_cell();
        let zonal = mean(&|c| band_means[band[c] as usize].powi(2)).sqrt();
        let noise = mean(&|c| 0.25 * (ra[c] - rb[c]).powi(2)).sqrt();
        let covariance = mean(&|c| ra[c] * rb[c]);
        let (va, vb) = (mean(&|c| ra[c] * ra[c]), mean(&|c| rb[c] * rb[c]));
        let correlation = if va > 0.0 && vb > 0.0 {
            covariance / (va * vb).sqrt()
        } else {
            0.0
        };
        let (pa, pb) = (
            self.symmetry_parts(&ra, mirror),
            self.symmetry_parts(&rb, mirror),
        );
        let symmetric: [f64; 4] = std::array::from_fn(|k| mean(&|c| pa[k][c] * pb[k][c]));
        let others = (symmetric[1] + symmetric[2] + symmetric[3]) / 3.0;
        let shared = |cells: &[usize]| {
            let m = cells.iter().map(|&c| ra[c] * rb[c]).sum::<f64>() / cells.len().max(1) as f64;
            m.max(0.0).sqrt()
        };
        let [locked_residual, ..] = self.symmetry_parts(&residual, mirror);
        Imprint {
            name: name.into(),
            units: units.into(),
            rms,
            zonal,
            noise,
            covariance,
            persistent: covariance.max(0.0).sqrt(),
            correlation,
            symmetric,
            locked: (symmetric[0] - others).max(0.0).sqrt(),
            transient: mean(&|c| transient[c]),
            pentagon: shared(&self.pentagons),
            ring: shared(&self.ring),
            waves: self.cross_power(&ra, &rb),
            residual,
            locked_residual,
        }
    }

    /// `P_m = Σ_belts w_b · 2 Re(c_m^A conj(c_m^B))` for m = 1..=12, with
    /// `c_m = Σ Aᵢ rᵢ e^(−imλᵢ) / Σ Aᵢ` over a belt's cells and `w_b` the
    /// belt's share of the sphere.
    fn cross_power(&self, ra: &[f64], rb: &[f64]) -> [f64; WAVENUMBERS] {
        let mut out = [0.0; WAVENUMBERS];
        for (m, p) in out.iter_mut().enumerate() {
            let m = (m + 1) as f64;
            // Per belt: Σ w r cos mλ and Σ w r sin mλ for each half.
            let mut sums = [[0.0f64; 4]; BELTS];
            for c in 0..self.areas.len() {
                let (s, co) = fms_math::sin_cos(m * self.longitude[c]);
                let w = self.weights[c];
                let sum = &mut sums[self.belt_of_cell[c]];
                sum[0] += w * ra[c] * co;
                sum[1] += w * ra[c] * s;
                sum[2] += w * rb[c] * co;
                sum[3] += w * rb[c] * s;
            }
            *p = sums
                .iter()
                .zip(&self.belt_weight)
                .filter(|(_, bw)| **bw > 0.0)
                .map(|(s, &bw)| {
                    // c = (Σ w r cos − i Σ w r sin) / bw; Re(cA conj cB).
                    bw * 2.0 * (s[0] * s[2] + s[1] * s[3]) / (bw * bw)
                })
                .sum();
        }
        out
    }
}

/// The constant term `c₀` of the weighted least-squares fit
/// `v ≈ c₀ + c₁x + c₂x²`, from the moments `s[k] = Σ kⱼ xⱼᵏ` and
/// `t[k] = Σ kⱼ xⱼᵏ vⱼ`; a straight line if that is singular, then the
/// weighted mean.
fn fit_constant(s: &[f64; 5], t: &[f64; 3]) -> f64 {
    const TINY: f64 = 1e-10;
    // Cramer's rule on the normal equations.
    let det3 = s[0] * (s[2] * s[4] - s[3] * s[3]) - s[1] * (s[1] * s[4] - s[3] * s[2])
        + s[2] * (s[1] * s[3] - s[2] * s[2]);
    if det3.abs() > TINY * s[0] * s[2] * s[4] {
        let num = t[0] * (s[2] * s[4] - s[3] * s[3]) - s[1] * (t[1] * s[4] - s[3] * t[2])
            + s[2] * (t[1] * s[3] - s[2] * t[2]);
        return num / det3;
    }
    let det2 = s[0] * s[2] - s[1] * s[1];
    if det2.abs() > TINY * s[0] * s[2] {
        return (t[0] * s[2] - s[1] * t[1]) / det2;
    }
    if s[0] > 0.0 { t[0] / s[0] } else { 0.0 }
}

/// Per-cell sums of one field: per half, and of squares over all samples.
#[derive(Clone, Debug)]
struct Sums {
    halves: [Vec<f64>; 2],
    squares: Vec<f64>,
}

/// Samples a run for [`ImprintAnalysis`]: the [`FIELDS`] per cell, summed
/// separately over half A (even blocks of samples) and half B (odd blocks),
/// with the sum of squares over all samples for the transient standard
/// deviation.
#[derive(Clone, Debug)]
pub struct GridImprint {
    analysis: ImprintAnalysis,
    upper: usize,
    lower: usize,
    block: u64,
    samples: u64,
    counts: [u64; 2],
    sums: Vec<Sums>,
}

impl GridImprint {
    /// Samples `model`'s runs in blocks of `block` samples, reading the
    /// eastward wind of layer `upper` and the temperature and northward wind
    /// of layer `lower`; the zonal mean is taken over `bands` latitude bands.
    ///
    /// # Panics
    /// If `block` is zero or a layer index is out of range.
    pub fn new(model: &LayeredModel, upper: usize, lower: usize, block: u64, bands: usize) -> Self {
        assert!(block > 0, "a block holds at least one sample");
        assert!(
            upper < model.layers() && lower < model.layers(),
            "no such layer"
        );
        let ops = model.operators();
        let n = ops.mesh().n_cells();
        let sums = Sums {
            halves: [vec![0.0; n], vec![0.0; n]],
            squares: vec![0.0; n],
        };
        Self {
            analysis: ImprintAnalysis::new(ops.mesh(), ops.geometry(), bands),
            upper,
            lower,
            block,
            samples: 0,
            counts: [0, 0],
            sums: vec![sums; FIELDS.len()],
        }
    }

    /// Add one sample.
    pub fn sample(&mut self, model: &LayeredModel, state: &LayeredState) {
        let ops = model.operators();
        let (mesh, geometry) = (ops.mesh(), ops.geometry());
        let centers = mesh.cell_centers();
        let component = |k: usize, north: bool| -> Vec<f64> {
            let vectors = geometry.reconstruct_vectors(mesh, &state.u[k]);
            centers
                .iter()
                .zip(&vectors)
                .map(|(&p, &v)| {
                    let (e, n) = east_north(p, v);
                    if north { n } else { e }
                })
                .collect()
        };
        let temperature: Vec<f64> = state
            .ps
            .iter()
            .zip(&state.theta_ps[self.lower])
            .map(|(&ps, &t)| t / ps * model.layer_exner(ps)[self.lower])
            .collect();
        let fields = [
            state.ps.clone(),
            temperature,
            component(self.upper, false),
            component(self.lower, true),
        ];
        let half = ((self.samples / self.block) % 2) as usize;
        for (sums, values) in self.sums.iter_mut().zip(&fields) {
            for ((s, q), &v) in sums.halves[half]
                .iter_mut()
                .zip(&mut sums.squares)
                .zip(values)
            {
                *s += v;
                *q += v * v;
            }
        }
        self.counts[half] += 1;
        self.samples += 1;
    }

    /// Samples taken.
    pub fn samples(&self) -> u64 {
        self.samples
    }

    /// Samples per block.
    pub fn block(&self) -> u64 {
        self.block
    }

    /// Complete blocks sampled.
    pub fn complete_blocks(&self) -> u64 {
        self.samples / self.block
    }

    /// Samples in half A and half B.
    pub fn counts(&self) -> [u64; 2] {
        self.counts
    }

    /// The analysis the results come from.
    pub fn analysis(&self) -> &ImprintAnalysis {
        &self.analysis
    }

    /// The time means of field `field` (an index into [`FIELDS`]) over half
    /// A and half B (zero for an empty half).
    pub fn means(&self, field: usize) -> [Vec<f64>; 2] {
        let sums = &self.sums[field];
        [0, 1].map(|h| {
            let n = self.counts[h].max(1) as f64;
            sums.halves[h].iter().map(|s| s / n).collect()
        })
    }

    /// The standard deviation in time of field `field` over all samples, per
    /// cell.
    pub fn transient(&self, field: usize) -> Vec<f64> {
        let sums = &self.sums[field];
        let n = self.samples.max(1) as f64;
        sums.squares
            .iter()
            .zip(&sums.halves[0])
            .zip(&sums.halves[1])
            .map(|((q, a), b)| {
                let mean = (a + b) / n;
                (q / n - mean * mean).max(0.0).sqrt()
            })
            .collect()
    }

    /// The imprint measures of every field, in the order of [`FIELDS`].
    pub fn results(&self) -> Vec<Imprint> {
        FIELDS
            .iter()
            .enumerate()
            .map(|(f, &field)| {
                let [a, b] = self.means(f);
                self.analysis
                    .analyze(field, [&a, &b], self.counts, &self.transient(f))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_fit_recovers_a_quadratic_and_falls_back_when_singular() {
        // v = 2 − x + 3x² at x = −1, −0.5, 0, 0.5, 1 with weights 1.
        let xs = [-1.0f64, -0.5, 0.0, 0.5, 1.0];
        let moments = |xs: &[f64], v: &dyn Fn(f64) -> f64| {
            let (mut s, mut t) = ([0.0; 5], [0.0; 3]);
            for &x in xs {
                for (k, sk) in s.iter_mut().enumerate() {
                    *sk += x.powi(k as i32);
                }
                for (k, tk) in t.iter_mut().enumerate() {
                    *tk += x.powi(k as i32) * v(x);
                }
            }
            (s, t)
        };
        let (s, t) = moments(&xs, &|x| 2.0 - x + 3.0 * x * x);
        assert!((fit_constant(&s, &t) - 2.0).abs() < 1e-12);
        // Two points: a straight line through them.
        let (s, t) = moments(&[-1.0, 1.0], &|x| 4.0 + x);
        assert!((fit_constant(&s, &t) - 4.0).abs() < 1e-12);
        // One point: its value.
        let (s, t) = moments(&[0.3], &|_| 7.0);
        assert!((fit_constant(&s, &t) - 7.0).abs() < 1e-12);
    }
}
