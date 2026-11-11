//! The ocean: a two-layer slab after Codron (2012) and Bhatnagar et al.
//! (2026) with sea ice ([`IceParams`]), wind-driven Ekman transport
//! ([`EkmanParams`]) and Gent–McWilliams eddy transport ([`EddyParams`]).
//!
//! Each cell with any ocean has a mixed layer, which exchanges heat with the
//! atmosphere, over a deep layer, which exchanges heat only with the mixed
//! layer and its neighbours. Heat travels by, in this order:
//!
//! - **Ekman transport** (if asked for): the wind's stress drives a mass
//!   flux in the mixed layer (the paper's Eq. 4) across every edge between
//!   two cells that have ocean, and the deep layer carries the same mass
//!   back. Where the mixed layer's flow diverges, deep water wells up into
//!   it; where it converges, mixed water sinks. Each layer's temperature is
//!   carried upstream, which on hexagons replaces the paper's centred
//!   meridional flux.
//! - **Eddy transport** (if asked for): a second mass flux in the mixed
//!   layer, down the slope of the isotherms between the two layers, returned
//!   and made up from below in the same way and in the same step. It
//!   flattens the isotherms: heat moves towards the cold side at the surface
//!   and the deep layer cools on the warm side.
//! - **Horizontal diffusion** in both layers, as a flux across every edge
//!   between two cells that have ocean, so the ocean's heat is conserved to
//!   round-off and none crosses a coast. An edge's opening is the smaller of
//!   its two cells' ocean fractions, which also keeps a nearly dry cell from
//!   taking more heat than it can hold. Ekman and eddy transport cross the
//!   same openings.
//! - **Convective adjustment**: a mixed layer colder than the deep layer
//!   below it overturns, and both take their mass-weighted mean.
//!
//! The diffusivity is the one the paper's runs actually have. The paper gives
//! 8000 m² s⁻¹, but its code scales that by 50 m over the ocean's total depth
//! (`ocean_slab_mod.F90`), which is 2000 m² s⁻¹ for 50 m over 150 m.
//!
//! Both eddy coefficients, the diffusivity and the eddy transport's, are
//! Earth's. On another planet they scale with [`eddy_scale`].
//!
//! The crate knows nothing about the atmosphere: fluxes and stresses come in
//! as numbers and temperatures go out.

mod ice;

use fms_grid::Operators;
use fms_grid::operators::MIN_CHUNK;
use glam::DVec3;
pub use ice::{IceForcing, IceParams};
use ice::{SeaIce, Water};
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
    /// Horizontal diffusivity in both layers at Earth's rotation (m² s⁻¹);
    /// 0 for none.
    pub diffusivity: f64,
    /// Whether a mixed layer colder than the deep layer overturns.
    pub convective_adjustment: bool,
    /// The deep layer starts at the mixed layer's temperature, but no warmer
    /// than this (K). Without Ekman transport the deep layer can only lose
    /// heat, so a warmer start would come back as a false warm spell.
    pub deep_start_limit: f64,
    /// Sea ice; `None` lets the water cool past freezing.
    pub ice: Option<IceParams>,
    /// Wind-driven Ekman transport, which needs a deep layer; `None` keeps
    /// the water still.
    pub ekman: Option<EkmanParams>,
    /// Gent–McWilliams eddy transport, which needs a deep layer; `None` for
    /// none.
    pub eddies: Option<EddyParams>,
    /// Whether the diffusivity and the eddy transport's coefficient grow on
    /// a planet that turns more slowly than Earth and shrink on a faster
    /// one, by [`eddy_scale`].
    pub rotation_scaling: bool,
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
            ice: Some(IceParams::default()),
            ekman: None,
            eddies: None,
            rotation_scaling: true,
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
            ice: None,
            ekman: None,
            eddies: None,
            rotation_scaling: false,
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

/// Wind-driven Ekman transport after Codron (2012) and Bhatnagar et al.
/// (2026, Sect. 2.3.1). [`Default`] is the paper's.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EkmanParams {
    /// Frictional damping ε of the mixed layer's flow (s⁻¹). Where the
    /// Coriolis parameter `|f| ≫ ε` the flow is the Ekman transport across
    /// the wind, `−k × τ / f`; near the equator, where `|f| ≲ ε`, it turns
    /// to run with the wind, `τ / ε`.
    pub friction: f64,
    /// Near the equator, replace the meridional part of the frictional flow
    /// with the Sverdrup transport `curl τ / β`, as the paper's code does by
    /// default (`slab_sverdrup` and `slab_tyeq`). With the weight `w =
    /// exp(−f²/ε²)` the transport per unit width becomes
    ///
    /// `M = (1 − w)/f (−k × τ) + (w/ε) τ + w B (k · ∇ × τ) north`,
    ///
    /// which is the Ekman transport where `|f| ≫ ε` and is carried by the
    /// friction and the wind's curl within a few degrees of the equator
    /// (about 4° on Earth). The code keeps only the `−∂τₓ/∂y` half of the
    /// curl; this is the whole curl. `B = β / (β² + (ε/R)²)`, with `β = 2Ω
    /// cos φ / R`, is our own: it is `1/β` wherever `β ≫ ε/R` (0.5% less on
    /// Earth's equator) but never more than `R/(2ε)`, and vanishes on a
    /// planet that doesn't turn, where the paper notes the bare `1/β` makes
    /// the transport unrealistically large. `false` is the paper's Eq. 4
    /// everywhere.
    pub sverdrup: bool,
}

impl Default for EkmanParams {
    fn default() -> Self {
        Self {
            friction: 1e-5,
            sverdrup: true,
        }
    }
}

/// `(1 − exp(−f²/ε²)) / f`, the Ekman weight of the Sverdrup blend (s): `1/f`
/// where `|f| ≫ ε`, `f/ε²` as `f → 0`, and exactly zero at `f = 0`.
fn ekman_factor(f: f64, epsilon: f64) -> f64 {
    let x = (f / epsilon) * (f / epsilon);
    if x == 0.0 {
        f / (epsilon * epsilon)
    } else {
        -fms_math::exp_m1(-x) / f
    }
}

/// `B = β / (β² + (ε/R)²)`, the limited `1/β` of the Sverdrup transport
/// (m s), at most `R/(2ε)`.
fn sverdrup_factor(beta: f64, epsilon: f64, radius: f64) -> f64 {
    let damping = epsilon / radius;
    beta / (beta * beta + damping * damping)
}

/// Eddy transport after Gent & McWilliams (1990), as Bhatnagar et al. (2026,
/// Sect. 2.3.2) apply it between two layers: the mixed layer carries
/// `−ρ κ s` per unit width, where `s` is the slope of the isotherms, and the
/// deep layer carries it back.
///
/// The slope is the horizontal gradient of the two layers' mean temperature
/// over the vertical gradient `(T_m − T_d) / separation`. The heat this
/// carries, `ρ c κ · separation · ∇T`, doesn't depend on the layers'
/// difference until the slope reaches its limit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EddyParams {
    /// Thickness diffusivity κ at Earth's rotation (m² s⁻¹). The paper's
    /// 2000; its code's own default (`slab_gmkappa`) is 1000.
    pub diffusivity: f64,
    /// Largest slope of the isotherms the flux follows; the paper's 0.002.
    pub max_slope: f64,
    /// Depth the two layers' temperature difference is taken to span (m).
    /// The code's 500 (`slab_heat_transp_mod.F90`), which stands for a
    /// thermocline: the layers' centres are only 100 m apart.
    pub separation: f64,
    /// Smallest difference between the layers the slope is taken over (K),
    /// which keeps it finite where they have mixed; the code's 0.25.
    pub min_contrast: f64,
}

impl Default for EddyParams {
    fn default() -> Self {
        Self {
            diffusivity: 2000.0,
            max_slope: 2e-3,
            separation: 500.0,
            min_contrast: 0.25,
        }
    }
}

/// Earth's rotation rate relative to the stars (rad s⁻¹), exactly as
/// `fms-planet` derives it from the solar day and the year, so that
/// [`eddy_scale`] is exactly 1 on Earth.
pub const EARTH_ROTATION: f64 = 7.292_115_854_615_987e-5;

/// Most that slow rotation may multiply the eddy coefficients by.
pub const MAX_EDDY_SCALE: f64 = 4.0;

/// What the eddy coefficients, which are given for Earth, are multiplied by
/// on a planet turning at `rotation` (rad s⁻¹): Earth's rate over the
/// planet's, at most [`MAX_EDDY_SCALE`].
///
/// Eddies are about as wide as the deformation radius, which grows as the
/// rotation slows, and a mixing length that wide gives a diffusivity in
/// proportion. The paper leaves this open (Sect. 5.6); the limit is ours,
/// because on a slow rotator the eddies reach the size of the planet and the
/// scaling stops meaning anything.
pub fn eddy_scale(rotation: f64) -> f64 {
    (EARTH_ROTATION / rotation.abs()).min(MAX_EDDY_SCALE)
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
    /// The mesh edge.
    edge: u32,
    /// Its cells, as the mesh orders them.
    cells: [u32; 2],
    /// Edge length over the distance between the cell centres, times the
    /// edge's opening.
    conductance: f64,
    /// Edge length times its opening (m).
    width: f64,
    /// Distance between the cell centres (m).
    distance: f64,
}

/// Largest share of the difference from its neighbours that a cell may close
/// in one diffusion step. Explicit diffusion is monotonic up to 1.
const MAX_DIFFUSION_NUMBER: f64 = 0.5;

/// Largest share of a layer's mass that may leave a cell in one advection
/// step. Donor-cell advection is monotonic up to 1.
const MAX_ADVECTION_NUMBER: f64 = 0.5;

/// What Ekman transport needs to turn the wind into mass fluxes.
#[derive(Clone, Debug)]
struct Ekman {
    friction: f64,
    /// Unit east and north vectors at each cell's centre; zero at a pole,
    /// where east is undefined (as in `fms_grid::zonal::east_north`).
    frames: Vec<[DVec3; 2]>,
    /// One per link.
    faces: Vec<Face>,
    /// Eastward and northward stress on each cell's water, per m² of its
    /// ocean: the stress on open water times the share that is open (N m⁻²).
    stress: Vec<[f64; 2]>,
    /// The stress's curl, for the Sverdrup blend.
    curl: Option<StressCurl>,
}

/// The vertical component of the stress's curl at each vertex: the
/// circulation of the stress around the triangle joining the vertex's three
/// cell centres, as [`Operators::curl`] takes it of an edge field, the edge
/// field being each link's mean stress along its normal.
#[derive(Clone, Debug)]
struct StressCurl {
    /// Per vertex, the links of its three edges and the weight of each,
    /// `t_{e,v} dₑ / (A_v R)` (m⁻¹); `None` where any of its cells has no
    /// ocean, where the stress falls to nothing at the coast and that isn't
    /// the wind's curl.
    stencil: Vec<Option<[(u32, f64); 3]>>,
    /// Each link's mean stress along its normal (N m⁻²).
    along: Vec<f64>,
    /// The curl at each vertex (N m⁻³).
    vertex: Vec<f64>,
}

/// Where a link's Ekman flux is evaluated.
#[derive(Clone, Copy, Debug)]
struct Face {
    /// The edge point, which is the local vertical `k`.
    point: DVec3,
    /// Unit normal at the edge point, from the link's `cells[0]` to
    /// `cells[1]`.
    normal: DVec3,
    /// Coriolis parameter `2Ω sin φ` at the edge point (s⁻¹).
    coriolis: f64,
    /// For the Sverdrup blend: the edge's two vertices, the northward unit
    /// vector at the edge point (zero on a pole), and the weights of the
    /// Ekman, frictional and Sverdrup terms, `(1 − w)/f` (s), `w/ε` (s) and
    /// `w B` (m s).
    vertices: [u32; 2],
    north: DVec3,
    blend: [f64; 3],
}

/// The mean stress of a link's two cells, each a tangent vector at its
/// centre from its east and north `frames` (N m⁻²).
fn link_stress(frames: &[[DVec3; 2]], stress: &[[f64; 2]], link: &Link) -> DVec3 {
    let vector = |i: usize| {
        let ([east, north], [e, n]) = (frames[i], stress[i]);
        e * east + n * north
    };
    let [a, b] = link.cells.map(|c| c as usize);
    0.5 * (vector(a) + vector(b))
}

impl Ekman {
    /// The mixed layer's mass flux across each link (kg s⁻¹, positive from
    /// `cells[0]` to `cells[1]`), with the stress the mean of the two cells'
    /// and `f` and `k` those of the edge point, across the opened width:
    /// the paper's Eq. 4, `M = (ε τ − f k × τ) / (ε² + f²)` per unit width,
    /// or its blend with the Sverdrup transport ([`EkmanParams::sverdrup`]),
    /// where an edge's curl is the mean of its two vertices'. Work per link
    /// and per vertex only, so the result doesn't depend on the thread
    /// count.
    fn mass_flux(&mut self, links: &[Link], flux: &mut [f64]) {
        let Self {
            friction: epsilon,
            frames,
            faces,
            stress,
            curl,
        } = self;
        let (epsilon, frames, stress) = (*epsilon, &*frames, &*stress);
        let link_stress = |link: &Link| link_stress(frames, stress, link);
        let Some(curl) = curl else {
            flux.par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .zip(links)
                .zip(&*faces)
                .for_each(|((flux, link), face)| {
                    let tau = link_stress(link);
                    let f = face.coriolis;
                    let m =
                        (epsilon * tau - f * face.point.cross(tau)) / (epsilon * epsilon + f * f);
                    *flux = m.dot(face.normal) * link.width;
                });
            return;
        };
        curl.along
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(links)
            .zip(&*faces)
            .for_each(|((along, link), face)| *along = link_stress(link).dot(face.normal));
        let along = &curl.along;
        curl.vertex
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(&curl.stencil)
            .for_each(|(out, stencil)| {
                *out = stencil.map_or(0.0, |terms| {
                    terms.iter().map(|&(l, w)| w * along[l as usize]).sum()
                });
            });
        let vertex = &curl.vertex;
        flux.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(links)
            .zip(&*faces)
            .for_each(|((flux, link), face)| {
                let tau = link_stress(link);
                let [ekman, friction, sverdrup] = face.blend;
                let m = -ekman * face.point.cross(tau) + friction * tau;
                let [v0, v1] = face.vertices.map(|v| vertex[v as usize]);
                let curl = 0.5 * (v0 + v1);
                *flux = (m.dot(face.normal) + sverdrup * curl * face.north.dot(face.normal))
                    * link.width;
            });
    }
}

/// The eddy transport's mass flux in the mixed layer across each link (kg
/// s⁻¹, positive from `cells[0]` to `cells[1]`): `−ρ κ s` across the opened
/// width, with `κ` the scaled `diffusivity` and `s` the limited slope of the
/// isotherms ([`EddyParams`]). Work per link only, so the result doesn't
/// depend on the thread count.
fn eddy_flux(
    eddies: &EddyParams,
    diffusivity: f64,
    density: f64,
    links: &[Link],
    mixed: &[f64],
    deep: &[f64],
    flux: &mut [f64],
) {
    flux.par_iter_mut()
        .with_min_len(MIN_CHUNK)
        .zip(links)
        .for_each(|(flux, link)| {
            let [a, b] = link.cells.map(|c| c as usize);
            let gradient = 0.5 * ((mixed[b] + deep[b]) - (mixed[a] + deep[a])) / link.distance;
            let contrast =
                (0.5 * ((mixed[a] - deep[a]) + (mixed[b] - deep[b]))).max(eddies.min_contrast);
            let slope = (gradient * eddies.separation / contrast)
                .clamp(-eddies.max_slope, eddies.max_slope);
            *flux = -density * diffusivity * slope * link.width;
        });
}

/// The mass fluxes that overturn the ocean, by what drives them.
const EKMAN: usize = 0;
const EDDY: usize = 1;

/// Two-layer advection by mass fluxes across links: the mixed layer carries
/// each link's flux and the deep layer carries it back, and what the mixed
/// layer's flow leaves a cell short of is made up from below. Holds what the
/// latest call left, and its scratch. Arrays of two are indexed by
/// [`EKMAN`] and [`EDDY`], then by layer, mixed first.
#[derive(Clone, Debug)]
struct Overturning {
    /// Net horizontal outflow of each cell's mixed-layer mass by each flux,
    /// which the deep layer makes up (kg s⁻¹; positive is upwelling).
    rise: [Vec<f64>; 2],
    /// Mass leaving each cell's mixed and deep layer (kg s⁻¹).
    leaving: [Vec<f64>; 2],
    /// Heat convergence by each flux in each layer in the current sub-step,
    /// per unit of specific heat (kg K s⁻¹).
    convergence: [[Vec<f64>; 2]; 2],
    /// Mean heating by each flux of each layer over the latest call (W per
    /// m² of ocean).
    heating: [[Vec<f64>; 2]; 2],
    /// Sub-steps the latest call took.
    steps: usize,
}

impl Overturning {
    fn new(n: usize) -> Self {
        let pair = || [vec![0.0; n], vec![0.0; n]];
        Self {
            rise: pair(),
            leaving: pair(),
            convergence: [pair(), pair()],
            heating: [pair(), pair()],
            steps: 0,
        }
    }

    /// Advect both layers for `dt` seconds by the sum of `fluxes`, the mixed
    /// layer's mass fluxes across each link (kg s⁻¹, positive from
    /// `cells[0]` to `cells[1]`), held fixed. Temperatures are taken
    /// upstream of the net flow in each layer and in the vertical, as the
    /// paper's code takes the vertical, and the step is split so that no
    /// cell loses more than [`MAX_ADVECTION_NUMBER`] of either layer's mass
    /// in one sub-step, which keeps every new temperature a mean of old
    /// ones. Every flux moves heat between two stores, so `Σ area (H_m T_m +
    /// H_d T_d)` is conserved to round-off, by each flux on its own. Serial
    /// over links, so the result doesn't depend on the thread count.
    #[allow(clippy::too_many_arguments)] // Each is one of the ocean's fields.
    fn advect(
        &mut self,
        mixed: &mut [f64],
        deep: &mut [f64],
        links: &[Link],
        fluxes: &[Vec<f64>; 2],
        area: &[f64],
        params: &OceanParams,
        dt: f64,
    ) {
        for field in self.rise.iter_mut().chain(&mut self.leaving) {
            field.fill(0.0);
        }
        for (l, link) in links.iter().enumerate() {
            let [a, b] = link.cells.map(|c| c as usize);
            for (rise, flux) in self.rise.iter_mut().zip(fluxes) {
                rise[a] += flux[l];
                rise[b] -= flux[l];
            }
            // The mixed layer leaves the cell `f` leaves; the deep layer,
            // carrying `−f`, leaves the other.
            let f = fluxes[EKMAN][l] + fluxes[EDDY][l];
            let (upstream, downstream) = if f > 0.0 { (a, b) } else { (b, a) };
            self.leaving[0][upstream] += f.abs();
            self.leaving[1][downstream] += f.abs();
        }
        // Mass per m² of ocean (kg m⁻²).
        let mass = [
            params.density * params.mixed_depth,
            params.density * params.deep_depth,
        ];
        let mut largest = 0.0f64;
        for (i, &a) in area.iter().enumerate() {
            if a > 0.0 {
                // Upwelling leaves the deep layer; downwelling the mixed one.
                let w = self.rise[EKMAN][i] + self.rise[EDDY][i];
                let mixed_out = self.leaving[0][i] + (-w).max(0.0);
                let deep_out = self.leaving[1][i] + w.max(0.0);
                largest = largest
                    .max(dt * mixed_out / (mass[0] * a))
                    .max(dt * deep_out / (mass[1] * a));
            }
        }
        let steps = (largest / MAX_ADVECTION_NUMBER).ceil().max(1.0);
        let h = dt / steps;
        for heating in self.heating.iter_mut().flatten() {
            heating.fill(0.0);
        }
        for _ in 0..steps as usize {
            for convergence in self.convergence.iter_mut().flatten() {
                convergence.fill(0.0);
            }
            for (l, link) in links.iter().enumerate() {
                let [a, b] = link.cells.map(|c| c as usize);
                let f = fluxes[EKMAN][l] + fluxes[EDDY][l];
                let (tm, td) = if f > 0.0 {
                    (mixed[a], deep[b])
                } else {
                    (mixed[b], deep[a])
                };
                for ([mixed_in, deep_in], flux) in self.convergence.iter_mut().zip(fluxes) {
                    let f = flux[l];
                    mixed_in[a] -= f * tm;
                    mixed_in[b] += f * tm;
                    deep_in[a] += f * td;
                    deep_in[b] -= f * td;
                }
            }
            for (i, &a) in area.iter().enumerate() {
                if a > 0.0 {
                    let w = self.rise[EKMAN][i] + self.rise[EDDY][i];
                    let t = if w > 0.0 { deep[i] } else { mixed[i] };
                    for (c, [mixed_in, deep_in]) in self.convergence.iter_mut().enumerate() {
                        let w = self.rise[c][i];
                        mixed_in[i] += w * t;
                        deep_in[i] -= w * t;
                        self.heating[c][0][i] += h * mixed_in[i] / a;
                        self.heating[c][1][i] += h * deep_in[i] / a;
                    }
                    let [ekman, eddy] = &self.convergence;
                    mixed[i] += h * (ekman[0][i] + eddy[0][i]) / (mass[0] * a);
                    deep[i] += h * (ekman[1][i] + eddy[1][i]) / (mass[1] * a);
                }
            }
        }
        let scale = if dt > 0.0 {
            params.specific_heat / dt
        } else {
            0.0
        };
        for heating in self.heating.iter_mut().flatten() {
            heating.iter_mut().for_each(|q| *q *= scale);
        }
        self.steps = steps as usize;
    }
}

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
    ice: Option<SeaIce>,
    ekman: Option<Ekman>,
    /// What the planet's rotation multiplies the eddy coefficients by.
    eddy_scale: f64,
    /// The mixed layer's mass fluxes across each link in the latest
    /// [`transport`](Ocean::transport), by [`EKMAN`] and [`EDDY`] (kg s⁻¹,
    /// positive from `cells[0]` to `cells[1]`).
    mass_flux: [Vec<f64>; 2],
    overturning: Overturning,
}

impl Ocean {
    /// An ocean on the mesh of `ops`, covering what `land_fraction` leaves
    /// of each cell, with its mixed layer at `surface` (K), on a planet
    /// turning at `rotation` (Ω relative to the stars, rad s⁻¹), which sets
    /// the Coriolis parameter `2Ω sin φ` of Ekman transport and scales the
    /// eddy coefficients.
    ///
    /// # Panics
    /// If the lengths differ from the mesh's cells, a land fraction is
    /// outside 0–1, the mixed layer has no depth, Ekman or eddy transport is
    /// asked for without a deep layer, or Ekman transport without friction.
    pub fn new(
        params: OceanParams,
        ops: &Operators,
        land_fraction: &[f64],
        surface: Vec<f64>,
        rotation: f64,
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
        if let Some(ekman) = params.ekman {
            assert!(
                params.deep_depth > 0.0,
                "Ekman transport needs a deep layer to return its flow"
            );
            assert!(
                ekman.friction > 0.0,
                "Ekman transport needs friction to stay finite at the equator"
            );
        }
        assert!(
            params.eddies.is_none() || params.deep_depth > 0.0,
            "eddy transport needs a deep layer to return its flow"
        );
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
                    edge: e as u32,
                    cells: edge.cells,
                    conductance: opening * geometry.edge_lengths[e] / geometry.dual_edge_lengths[e],
                    width: opening * geometry.edge_lengths[e] * ops.radius(),
                    distance: geometry.dual_edge_lengths[e] * ops.radius(),
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
        let floor = params.ice.map_or(f64::NEG_INFINITY, |p| p.freezing);
        let deep = surface
            .iter()
            .map(|t| t.min(params.deep_start_limit).max(floor))
            .collect();
        let ekman = params.ekman.map(|p| {
            let frames = mesh
                .cell_centers()
                .iter()
                .map(|&x| {
                    let east = DVec3::Z.cross(x);
                    let length = east.length();
                    if length < 1e-12 {
                        return [DVec3::ZERO; 2];
                    }
                    let east = east / length;
                    [east, x.cross(east)]
                })
                .collect();
            let radius = ops.radius();
            let epsilon = p.friction;
            let faces = links
                .iter()
                .map(|link| {
                    let e = link.edge as usize;
                    let point = geometry.edge_points[e];
                    let coriolis = 2.0 * rotation * point.z;
                    let east = DVec3::Z.cross(point);
                    let length = east.length();
                    let north = if length < 1e-12 {
                        DVec3::ZERO
                    } else {
                        point.cross(east / length)
                    };
                    let weight = fms_math::exp(-(coriolis / epsilon) * (coriolis / epsilon));
                    let beta = 2.0 * rotation * length / radius;
                    Face {
                        point,
                        normal: geometry.edge_normals[e],
                        coriolis,
                        vertices: mesh.edges()[e].vertices,
                        north,
                        blend: [
                            ekman_factor(coriolis, epsilon),
                            weight / epsilon,
                            weight * sverdrup_factor(beta, epsilon, radius),
                        ],
                    }
                })
                .collect();
            let curl = p.sverdrup.then(|| {
                let mut link_of = vec![u32::MAX; mesh.edges().len()];
                for (l, link) in links.iter().enumerate() {
                    link_of[link.edge as usize] = l as u32;
                }
                let stencil = (0..mesh.n_vertices())
                    .map(|v| {
                        let wet = mesh.triangles()[v].iter().all(|&c| area[c as usize] > 0.0);
                        let v = v as u32;
                        let (edges, signs) = (mesh.vertex_edges(v), geometry.vertex_edge_signs(v));
                        let scale = geometry.vertex_areas[v as usize] * radius;
                        wet.then(|| {
                            std::array::from_fn(|k| {
                                let e = edges[k] as usize;
                                debug_assert_ne!(
                                    link_of[e],
                                    u32::MAX,
                                    "a wet vertex's edges are links"
                                );
                                (link_of[e], signs[k] * geometry.dual_edge_lengths[e] / scale)
                            })
                        })
                    })
                    .collect();
                StressCurl {
                    stencil,
                    along: vec![0.0; links.len()],
                    vertex: vec![0.0; mesh.n_vertices()],
                }
            });
            Ekman {
                friction: epsilon,
                frames,
                faces,
                stress: vec![[0.0; 2]; n],
                curl,
            }
        });
        Self {
            params,
            mixed: surface,
            deep,
            area,
            eddy_scale: if params.rotation_scaling {
                eddy_scale(rotation)
            } else {
                1.0
            },
            mass_flux: [vec![0.0; links.len()], vec![0.0; links.len()]],
            links,
            largest_rate,
            convergence: vec![0.0; n],
            ice: params.ice.map(|p| SeaIce::new(p, n)),
            ekman,
            overturning: Overturning::new(n),
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
        let water = self.params.mixed_capacity() * self.mixed[i]
            + self.params.deep_capacity() * self.deep[i];
        self.ice.as_ref().map_or(water, |ice| water + ice.heat(i))
    }

    /// Heat the mixed layers gain over `dt` seconds from `flux(i)`, the net
    /// flux into cell `i`'s open water (W per m² of it).
    pub fn add_surface_heat(&mut self, dt: f64, flux: impl Fn(usize) -> f64 + Sync) {
        let capacity = self.params.mixed_capacity();
        let area = &self.area;
        let ice = self.ice.as_ref();
        self.mixed
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .filter(|(i, _)| area[*i] > 0.0)
            .for_each(|(i, t)| {
                let open = 1.0 - ice.map_or(0.0, |ice| ice.fraction[i]);
                *t += dt * (open * flux(i)) / capacity;
            });
    }

    /// The wind's stress on each cell's open water, `stress(i)` eastward and
    /// northward (N per m² of open water), for the Ekman transport of the
    /// following [`transport`](Self::transport) calls. The ocean keeps it
    /// scaled by the share of the cell's ocean that is open now. Does nothing
    /// without Ekman transport.
    pub fn set_wind_stress(&mut self, stress: impl Fn(usize) -> [f64; 2] + Sync) {
        let Some(ekman) = &mut self.ekman else {
            return;
        };
        let (area, ice) = (&self.area, self.ice.as_ref());
        ekman
            .stress
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, s)| {
                *s = if area[i] > 0.0 {
                    let open = 1.0 - ice.map_or(0.0, |ice| ice.fraction[i]);
                    stress(i).map(|t| open * t)
                } else {
                    [0.0; 2]
                };
            });
    }

    /// Move heat for `dt` seconds: Ekman transport under the latest
    /// [wind stress](Self::set_wind_stress) together with eddy transport,
    /// diffusion in both layers, then convective adjustment.
    pub fn transport(&mut self, dt: f64) {
        if let Some(ekman) = &mut self.ekman {
            ekman.mass_flux(&self.links, &mut self.mass_flux[EKMAN]);
        }
        if let Some(eddies) = &self.params.eddies {
            eddy_flux(
                eddies,
                eddies.diffusivity * self.eddy_scale,
                self.params.density,
                &self.links,
                &self.mixed,
                &self.deep,
                &mut self.mass_flux[EDDY],
            );
        }
        if self.ekman.is_some() || self.params.eddies.is_some() {
            self.overturning.advect(
                &mut self.mixed,
                &mut self.deep,
                &self.links,
                &self.mass_flux,
                &self.area,
                &self.params,
                dt,
            );
        }
        let diffusivity = self.params.diffusivity * self.eddy_scale;
        if diffusivity > 0.0 {
            let number = diffusivity * dt * self.largest_rate;
            let steps = (number / MAX_DIFFUSION_NUMBER).ceil().max(1.0);
            let amount = diffusivity * dt / steps;
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

    /// What the planet's rotation multiplies the diffusivity and the eddy
    /// transport's coefficient by ([`eddy_scale`], or 1 without
    /// [`rotation_scaling`](OceanParams::rotation_scaling)).
    pub fn eddy_scale(&self) -> f64 {
        self.eddy_scale
    }

    /// Heating of cell `i`'s column, both layers, by the latest
    /// [`transport`](Self::transport)'s Ekman transport (W per m² of its
    /// ocean); zero without it. It only moves heat, so its integral over the
    /// ocean is zero.
    pub fn ekman_heating(&self, i: usize) -> f64 {
        let [mixed, deep] = &self.overturning.heating[EKMAN];
        mixed[i] + deep[i]
    }

    /// Heating of cell `i`'s mixed layer alone by the latest
    /// [`transport`](Self::transport)'s Ekman transport (W per m² of its
    /// ocean); zero without it.
    pub fn ekman_mixed_heating(&self, i: usize) -> f64 {
        self.overturning.heating[EKMAN][0][i]
    }

    /// Heating of cell `i`'s column, both layers, by the latest
    /// [`transport`](Self::transport)'s eddy transport (W per m² of its
    /// ocean); zero without it. Its integral over the ocean is zero.
    pub fn eddy_heating(&self, i: usize) -> f64 {
        let [mixed, deep] = &self.overturning.heating[EDDY];
        mixed[i] + deep[i]
    }

    /// Heating of cell `i`'s mixed layer alone by the latest
    /// [`transport`](Self::transport)'s eddy transport (W per m² of its
    /// ocean); zero without it.
    pub fn eddy_mixed_heating(&self, i: usize) -> f64 {
        self.overturning.heating[EDDY][0][i]
    }

    /// Deep water welling up into cell `i`'s mixed layer in the latest
    /// [`transport`](Self::transport), by Ekman and eddy transport together
    /// (kg per m² of its ocean per second; negative where the mixed layer
    /// sinks); zero without either.
    pub fn upwelling(&self, i: usize) -> f64 {
        let a = self.area[i];
        if a > 0.0 {
            (self.overturning.rise[EKMAN][i] + self.overturning.rise[EDDY][i]) / a
        } else {
            0.0
        }
    }

    /// The mesh edge of each link between two cells that both have ocean,
    /// in the order of [`ekman_mass_flux`](Self::ekman_mass_flux).
    pub fn link_edges(&self) -> impl ExactSizeIterator<Item = u32> + '_ {
        self.links.iter().map(|link| link.edge)
    }

    /// The mixed layer's Ekman mass flux across each link in the latest
    /// [`transport`](Self::transport) (kg s⁻¹, positive along the edge's
    /// normal, from its `cells[0]` to its `cells[1]`); the deep layer
    /// carries the opposite. Zero without Ekman transport.
    pub fn ekman_mass_flux(&self) -> &[f64] {
        &self.mass_flux[EKMAN]
    }

    /// The mixed layer's eddy mass flux across each link in the latest
    /// [`transport`](Self::transport), as
    /// [`ekman_mass_flux`](Self::ekman_mass_flux). Zero without eddy
    /// transport.
    pub fn eddy_mass_flux(&self) -> &[f64] {
        &self.mass_flux[EDDY]
    }

    /// The sea-ice parameters, if the ocean has ice.
    pub fn ice_params(&self) -> Option<&IceParams> {
        self.params.ice.as_ref()
    }

    /// Share of cell `i`'s ocean under ice.
    pub fn ice_fraction(&self, i: usize) -> f64 {
        self.ice.as_ref().map_or(0.0, |ice| ice.fraction[i])
    }

    /// Mean thickness of cell `i`'s ice (m).
    pub fn ice_thickness(&self, i: usize) -> f64 {
        self.ice.as_ref().map_or(0.0, |ice| ice.thickness[i])
    }

    /// Surface temperature of cell `i`'s ice (K); the freezing point where
    /// it has none.
    pub fn ice_temperature(&self, i: usize) -> f64 {
        self.ice.as_ref().map_or(0.0, |ice| ice.temperature[i])
    }

    /// Temperature of the top of cell `i`'s ice, under any snow (K).
    pub fn ice_interface_temperature(&self, i: usize) -> f64 {
        self.ice.as_ref().map_or(0.0, |ice| ice.interface[i])
    }

    /// Snow on cell `i`'s ice (kg per m² of ice).
    pub fn snow_on_ice(&self, i: usize) -> f64 {
        self.ice.as_ref().map_or(0.0, |ice| ice.snow[i])
    }

    /// Heat the limit on thickness has removed from cell `i` so far (J per
    /// m² of its ocean); [`heat`](Self::heat) still counts it.
    pub fn ice_limit_heat(&self, i: usize) -> f64 {
        self.ice.as_ref().map_or(0.0, |ice| ice.limited[i])
    }

    /// Albedo of cell `i`'s ice in each band, rising with thickness from
    /// `open_water`'s towards its limit.
    pub fn ice_albedo(&self, i: usize, open_water: [f64; 2]) -> [f64; 2] {
        let Some(ice) = &self.ice else {
            return open_water;
        };
        let p = ice.params;
        let thin = fms_math::exp(-ice.thickness[i] / p.albedo_depth);
        let cover = (ice.snow[i] / p.snow_critical).clamp(0.0, 1.0);
        std::array::from_fn(|b| {
            let bare = p.albedo_limit[b] - (p.albedo_limit[b] - open_water[b]) * thin;
            cover * p.snow_albedo + (1.0 - cover) * bare
        })
    }

    /// What the surface of cell `i`'s ice exchanges with its interior over
    /// `dt` seconds: a conductance (W m⁻² K⁻¹) and the temperature it pulls
    /// towards (K).
    pub fn ice_balance(&self, i: usize, dt: f64) -> (f64, f64) {
        self.ice
            .as_ref()
            .expect("this ocean has no ice")
            .balance(i, dt)
    }

    /// Advance the ice's temperature `dt` seconds under `forcing(i)`, and
    /// take the heat the ocean gives it. Does nothing without ice.
    pub fn step_ice(&mut self, dt: f64, forcing: impl Fn(usize) -> IceForcing) {
        if let Some(ice) = &mut self.ice {
            let mut water = Water {
                mixed: &mut self.mixed,
                area: &self.area,
                capacity: self.params.mixed_capacity(),
                density: self.params.density,
                specific_heat: self.params.specific_heat,
            };
            ice.step(dt, &mut water, forcing);
        }
    }

    /// Grow, melt and spread the ice with the heat recorded since the last
    /// call, and freeze open water below the freezing point. Ice areas
    /// change only here. Does nothing without ice.
    pub fn rearrange_ice(&mut self) {
        let capacity = self.params.mixed_capacity();
        if let Some(ice) = &mut self.ice {
            ice.rearrange(&mut self.mixed, &self.area, capacity);
        }
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
        // No ice comes of a jump.
        let floor = self.params.ice.map_or(f64::NEG_INFINITY, |p| p.freezing);
        let extrapolate = |now: &mut [f64], then: &[f64]| {
            for ((t, t0), a) in now.iter_mut().zip(then).zip(&self.area) {
                if *a > 0.0 {
                    *t = (*t + (periods * (*t - t0)).clamp(-limit, limit)).max(floor);
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
    const ROTATION: f64 = EARTH_ROTATION;

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
        let mut ocean = Ocean::new(params, &grid.ops(), &land, start.clone(), ROTATION);
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
        let params = OceanParams {
            ice: None,
            ..OceanParams::default()
        };
        let ocean = Ocean::new(params, &grid.ops(), &land, start.clone(), ROTATION);
        for (d, s) in ocean.deep_temperature().iter().zip(&start) {
            assert_eq!(*d, s.min(277.15));
        }
    }

    #[test]
    fn surface_heat_goes_only_into_ocean_and_is_counted() {
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        let start = grid.each(|z| 285.0 + 12.0 * z);
        let mut ocean = Ocean::new(
            OceanParams::default(),
            &grid.ops(),
            &land,
            start.clone(),
            ROTATION,
        );
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
            let mut ocean = Ocean::new(params, &grid.ops(), &land, start.clone(), ROTATION);
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
        let mut ocean = Ocean::new(
            OceanParams::default(),
            &grid.ops(),
            &land,
            start.clone(),
            ROTATION,
        );
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
            ice: None,
            ..OceanParams::default()
        };
        let start: Vec<f64> = p2.iter().map(|x| 20.0 * x).collect();
        let mut ocean = Ocean::new(params, &grid.ops(), &land, start, ROTATION);
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
        let mut ocean = Ocean::new(OceanParams::default(), &grid.ops(), &land, start, ROTATION);
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
        let mut ocean = Ocean::new(params, &grid.ops(), &land, start, ROTATION);
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

    /// An ocean starting `cold` K below freezing, already turned into ice.
    fn iced(grid: &Grid, cold: f64) -> (Ocean, OceanParams) {
        let params = OceanParams {
            diffusivity: 0.0,
            ..OceanParams::default()
        };
        let freezing = params.ice.unwrap().freezing;
        let land = vec![0.0; grid.n()];
        let start = vec![freezing - cold; grid.n()];
        let mut ocean = Ocean::new(params, &grid.ops(), &land, start, ROTATION);
        ocean.rearrange_ice();
        (ocean, params)
    }

    #[test]
    fn open_water_below_freezing_turns_into_ice_and_keeps_its_heat() {
        let grid = Grid::new(4);
        let params = OceanParams {
            diffusivity: 0.0,
            ..OceanParams::default()
        };
        let p = params.ice.unwrap();
        let land = vec![0.0; grid.n()];
        let start = vec![p.freezing - 0.1; grid.n()];
        let mut ocean = Ocean::new(params, &grid.ops(), &land, start, ROTATION);
        assert_eq!(
            ocean.deep_temperature()[0],
            p.freezing,
            "the deep stays liquid"
        );
        let before = total(&ocean);
        ocean.rearrange_ice();
        assert!(((total(&ocean) - before) / before).abs() < 1e-14);
        assert_eq!(ocean.mixed_temperature()[0], p.freezing);
        let volume = 0.1 * params.mixed_capacity() / p.melting_energy(p.freezing);
        assert!((ocean.ice_thickness(0) - p.thin).abs() < 1e-12);
        assert!((ocean.ice_fraction(0) - volume / p.thin).abs() < 1e-12);
    }

    #[test]
    fn ice_grows_from_below_by_stefans_law() {
        let grid = Grid::new(4);
        let (mut ocean, params) = iced(&grid, 0.1);
        let p = params.ice.unwrap();
        // A surface held 18.2 K below freezing, losing what it conducts.
        let skin = p.freezing - 18.2;
        let ice = ocean.ice.as_mut().unwrap();
        ice.temperature.fill(skin);
        ice.interface.fill(skin);
        let (h0, f0) = (ocean.ice_thickness(0), ocean.ice_fraction(0));
        let dt = 900.0;
        let steps = 8 * 30 * 12;
        for n in 0..steps {
            let h = ocean.ice_thickness(0);
            let net = p.conductivity * (skin - p.freezing) / h;
            ocean.step_ice(dt, |_| IceForcing {
                net,
                skin,
                ..IceForcing::default()
            });
            if n % 12 == 11 {
                ocean.rearrange_ice();
            }
        }
        let time = dt * steps as f64;
        let per_volume = p.melting_energy(skin);
        let analytic = (h0 * h0 + 2.0 * p.conductivity * 18.2 * time / per_volume).sqrt();
        let h = ocean.ice_thickness(0);
        assert!(
            ((h - analytic) / analytic).abs() < 0.01,
            "{h} vs {analytic}"
        );
        assert!(h > 0.5 && h < p.thick, "{h}");
        assert_eq!(
            ocean.ice_fraction(0),
            f0,
            "below h_thick, only thickness grows"
        );
    }

    #[test]
    fn heat_is_conserved_through_freezing_and_melting() {
        let grid = Grid::new(4);
        let (mut ocean, params) = iced(&grid, 0.4);
        let dt = 300.0;
        let n = grid.n();
        let (mut melted, mut formed) = (0, 0);
        // Each cell swings between strong heating and cooling.
        for step in 0..12_000 {
            let phase = step as f64 / 2000.0 * std::f64::consts::TAU;
            let open: Vec<f64> = (0..n)
                .map(|i| 400.0 * fms_math::sin(phase + 6.0 * scatter(i, 1)))
                .collect();
            let net: Vec<f64> = (0..n)
                .map(|i| {
                    let swing = 400.0 * fms_math::sin(phase + 6.0 * scatter(i, 2));
                    swing - 30.0 * (ocean.ice_temperature(i) - 272.0)
                })
                .collect();
            let frac: Vec<f64> = (0..n).map(|i| ocean.ice_fraction(i)).collect();
            let before = total(&ocean);
            ocean.add_surface_heat(dt, |i| open[i]);
            ocean.step_ice(dt, |i| IceForcing {
                net: net[i],
                skin: 271.35 + (net[i] / 40.0).clamp(-15.0, 1.8),
                stress: 0.1,
                sublimation: 0.0,
                snowfall: 2e-5 * (1.0 + fms_math::sin(phase + 6.0 * scatter(i, 3))),
            });
            ocean.transport(dt);
            let gained: f64 = (0..n)
                .map(|i| ocean.area[i] * dt * ((1.0 - frac[i]) * open[i] + frac[i] * net[i]))
                .sum::<f64>()
                - (0..n)
                    .map(|i| {
                        let snowfall = 2e-5 * (1.0 + fms_math::sin(phase + 6.0 * scatter(i, 3)));
                        ocean.area[i] * dt * frac[i] * snowfall * params.ice.unwrap().fusion
                    })
                    .sum::<f64>();
            let error = (total(&ocean) - before - gained).abs();
            assert!(error < 1e-11 * before, "step {step}: {error}");
            if step % 12 == 11 {
                let before = total(&ocean);
                let had: Vec<bool> = (0..n).map(|i| ocean.ice_fraction(i) > 0.0).collect();
                ocean.rearrange_ice();
                assert!(((total(&ocean) - before) / before).abs() < 1e-12);
                for (i, had) in had.iter().enumerate() {
                    match (*had, ocean.ice_fraction(i) > 0.0) {
                        (true, false) => melted += 1,
                        (false, true) => formed += 1,
                        _ => {}
                    }
                }
            }
        }
        assert!(melted > 0 && formed > 0, "{melted} melted, {formed} formed");
        let p = params.ice.unwrap();
        for i in 0..n {
            assert!(ocean.ice_fraction(i) <= p.max_fraction);
            assert!(ocean.ice_thickness(i) <= p.max_thickness);
        }
    }

    #[test]
    fn melting_from_below_takes_area_then_warms_the_ocean() {
        let grid = Grid::new(4);
        let (mut ocean, params) = iced(&grid, 0.2);
        let p = params.ice.unwrap();
        let (f, h) = (ocean.ice_fraction(0), ocean.ice_thickness(0));
        let e = p.melting_energy(ocean.ice_temperature(0));
        // Three times what it takes to melt all of it, through the base.
        ocean.ice.as_mut().unwrap().bottom.fill(3.0 * e * h);
        let before = total(&ocean);
        ocean.rearrange_ice();
        assert_eq!(ocean.ice_fraction(0), 0.0);
        let warmed = 2.0 * e * h * f / params.mixed_capacity();
        assert!((ocean.mixed_temperature()[0] - (p.freezing + warmed)).abs() < 1e-9);
        assert!(((total(&ocean) - before) / before).abs() < 1e-14);
    }

    #[test]
    fn the_ocean_gives_the_ice_heat_only_from_above_freezing() {
        let grid = Grid::new(4);
        let (mut ocean, params) = iced(&grid, 0.2);
        let p = params.ice.unwrap();
        let hold = IceForcing {
            skin: p.freezing,
            ..IceForcing::default()
        };
        ocean.step_ice(1800.0, |_| hold);
        assert_eq!(ocean.mixed_temperature()[0], p.freezing);
        // Warm the water; a long step takes at most what is above freezing.
        ocean.mixed.fill(p.freezing + 0.5);
        ocean.step_ice(1e9, |_| hold);
        let t = ocean.mixed_temperature()[0];
        assert!(t >= p.freezing && t < p.freezing + 0.5, "{t}");
    }

    #[test]
    fn thickness_stops_at_the_limit_and_the_heat_stays_on_record() {
        let grid = Grid::new(4);
        let (mut ocean, params) = iced(&grid, 0.2);
        let p = params.ice.unwrap();
        let e = p.melting_energy(ocean.ice_temperature(0));
        let ice = ocean.ice.as_mut().unwrap();
        ice.fraction.fill(0.5);
        ice.thickness.fill(p.max_thickness - 0.01);
        ice.bottom.fill(-20.0 * e);
        let before = total(&ocean);
        ocean.rearrange_ice();
        assert!(ocean.ice_thickness(0) <= p.max_thickness);
        assert!(ocean.ice_limit_heat(0) > 0.0);
        assert!(((total(&ocean) - before) / before).abs() < 1e-14);
    }

    #[test]
    fn ice_albedo_rises_with_thickness_from_the_waters_in_each_band() {
        let grid = Grid::new(4);
        let (mut ocean, params) = iced(&grid, 0.2);
        let p = params.ice.unwrap();
        let open = [0.08, 0.08];
        ocean.ice.as_mut().unwrap().thickness.fill(0.0);
        let bare = ocean.ice_albedo(0, open);
        assert!((bare[0] - open[0]).abs() < 1e-12 && (bare[1] - open[1]).abs() < 1e-12);
        ocean.ice.as_mut().unwrap().thickness.fill(100.0);
        let thick = ocean.ice_albedo(0, open);
        assert!((thick[0] - p.albedo_limit[0]).abs() < 1e-12);
        assert!((thick[1] - p.albedo_limit[1]).abs() < 1e-12);
    }

    /// An ice-covered ocean with 0.5 m of ice and `depth` metres of snow, at
    /// surface and interface temperatures of `t` K.
    fn snowy(grid: &Grid, depth: f64, t: f64) -> (Ocean, IceParams, OceanParams) {
        let (mut ocean, params) = iced(grid, 0.2);
        let p = params.ice.unwrap();
        let ice = ocean.ice.as_mut().unwrap();
        ice.fraction.fill(0.9);
        ice.thickness.fill(0.5);
        ice.snow.fill(depth * p.snow_density);
        ice.temperature.fill(t);
        ice.interface.fill(t);
        (ocean, p, params)
    }

    #[test]
    fn heat_conducts_through_snow_and_ice_in_series() {
        let grid = Grid::new(4);
        let (mut ocean, p, _) = snowy(&grid, 0.2, 260.0);
        let (h, d) = (0.5, 0.2);
        let skin = 250.0;
        let (ri, rs) = (h / p.conductivity, d / p.snow_conductivity);
        // A surface held at `skin`: it loses what conducts up to it.
        for _ in 0..20_000 {
            let net = (skin - ocean.ice_interface_temperature(0)) / rs;
            ocean.step_ice(1800.0, |_| IceForcing {
                net,
                skin,
                ..IceForcing::default()
            });
        }
        let expected = p.freezing - (p.freezing - skin) * ri / (ri + rs);
        let interface = ocean.ice_interface_temperature(0);
        assert!(
            (interface - expected).abs() < 1e-6,
            "{interface} vs {expected}"
        );
    }

    #[test]
    fn melting_takes_the_snow_then_the_ice_then_warms_the_ocean() {
        let grid = Grid::new(4);
        let (mut ocean, p, params) = snowy(&grid, 0.1, 272.0);
        let (s, f, h) = (0.1 * p.snow_density, 0.9, 0.5);
        let snow = p.snow_melting_energy(272.0);
        let ice = p.melting_energy(272.0);
        // Half the snow.
        ocean.ice.as_mut().unwrap().top.fill(0.5 * s * snow);
        let before = total(&ocean);
        ocean.rearrange_ice();
        assert!(((total(&ocean) - before) / before).abs() < 1e-14);
        assert!((ocean.snow_on_ice(0) - 0.5 * s).abs() < 1e-9);
        assert_eq!(ocean.ice_thickness(0), h);
        // The rest of it and 0.1 m of ice.
        let extra = 0.5 * s * snow + 0.1 * ice;
        ocean.ice.as_mut().unwrap().top.fill(extra);
        ocean.rearrange_ice();
        assert_eq!(ocean.snow_on_ice(0), 0.0);
        assert!((ocean.ice_thickness(0) - (h - 0.1)).abs() < 1e-9);
        assert_eq!(ocean.mixed_temperature()[0], p.freezing);
        // All the ice and more: a tenth of a metre's energy warms the water.
        let top = (h - 0.1) * ice + 0.1 * ice;
        ocean.ice.as_mut().unwrap().top.fill(top);
        let before = total(&ocean);
        ocean.rearrange_ice();
        assert_eq!(ocean.ice_fraction(0), 0.0);
        assert!(((total(&ocean) - before) / before).abs() < 1e-14);
        let warmed = f * 0.1 * ice / params.mixed_capacity();
        let t = ocean.mixed_temperature()[0];
        assert!((t - p.freezing - warmed).abs() < 1e-3 * warmed, "{t}");
    }

    #[test]
    fn snow_lands_with_its_heat_and_part_blows_into_the_water() {
        let grid = Grid::new(4);
        let (mut ocean, p, params) = snowy(&grid, 0.0, 260.0);
        let rate = 1e-4;
        let dt = 3600.0;
        let start = ocean.mixed_temperature()[0];
        let before = total(&ocean);
        ocean.step_ice(dt, |_| IceForcing {
            skin: 260.0,
            net: p.conductivity * (260.0 - p.freezing) / 0.5,
            snowfall: rate,
            ..IceForcing::default()
        });
        let blown = p.snow_blown * 0.1;
        let cooled = dt * 0.9 * blown * rate * p.fusion / params.mixed_capacity();
        assert!((start - ocean.mixed_temperature()[0] - cooled).abs() < 1e-12);
        let expected_heat = -0.9 * dt * rate * p.fusion;
        let change: f64 = total(&ocean) - before;
        let net = 0.9 * dt * p.conductivity * (260.0 - p.freezing) / 0.5;
        let area: f64 = ocean.area.iter().sum();
        assert!(
            (change / area - expected_heat - net).abs() < 1e-3,
            "{change}"
        );
        let before = total(&ocean);
        ocean.rearrange_ice();
        assert!(((total(&ocean) - before) / before).abs() < 1e-14);
        let laid = ocean.ice_fraction(0) * ocean.snow_on_ice(0);
        assert!(
            (laid - 0.9 * dt * rate * (1.0 - blown)).abs() < 1e-9,
            "{laid}"
        );
    }

    #[test]
    fn snow_below_the_waterline_floods_into_ice_and_keeps_the_heat() {
        let grid = Grid::new(4);
        // 0.5 m of ice bears 54.5 kg/m² of snow; this is 1 m of snow.
        let (mut ocean, p, _) = snowy(&grid, 1.0, 265.0);
        let before = total(&ocean);
        ocean.rearrange_ice();
        assert!(((total(&ocean) - before) / before).abs() < 1e-14);
        let (h, s) = (ocean.ice_thickness(0), ocean.snow_on_ice(0));
        assert!(h > 0.5 && s < p.snow_density, "{h} {s}");
        let freeboard = s + p.density * h - p.sea_density * h;
        assert!(freeboard.abs() < 1e-9, "{freeboard}");
        assert!(
            ocean.ice_interface_temperature(0) > 265.0,
            "freezing warms it"
        );
    }

    #[test]
    fn a_thin_snow_layer_merges_with_the_ice_and_keeps_the_heat() {
        let grid = Grid::new(4);
        let (mut ocean, p, _) = snowy(&grid, 0.1, 265.0);
        ocean.ice.as_mut().unwrap().temperature.fill(255.0);
        // Melt most of the snow, leaving it below the critical amount.
        let melted = (0.1 * p.snow_density - 5.0) * p.snow_melting_energy(260.0);
        ocean.ice.as_mut().unwrap().top.fill(melted);
        let before = total(&ocean);
        ocean.rearrange_ice();
        assert!(((total(&ocean) - before) / before).abs() < 1e-14);
        assert_eq!(ocean.ice_temperature(0), ocean.ice_interface_temperature(0));
        let t = ocean.ice_temperature(0);
        assert!(t > 255.0 && t < 265.0, "{t}");
    }

    #[test]
    fn snow_hides_the_ices_albedo() {
        let grid = Grid::new(4);
        let (ocean, p, _) = snowy(&grid, 0.2, 260.0);
        let covered = ocean.ice_albedo(0, [0.08; 2]);
        assert_eq!(covered, [p.snow_albedo; 2]);
    }

    #[test]
    fn a_jump_makes_no_ice() {
        let grid = Grid::new(4);
        let (mut ocean, params) = iced(&grid, 0.2);
        let p = params.ice.unwrap();
        let start = ocean.snapshot();
        ocean.mixed.fill(p.freezing);
        ocean.add_surface_heat(1.0, |_| -params.mixed_capacity());
        ocean.jump(&start, 5.0, 10.0);
        assert_eq!(ocean.mixed_temperature()[0], p.freezing);
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
        let mut ocean = Ocean::new(params, &grid.ops(), &land, vec![290.0; grid.n()], ROTATION);
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
        let mut ocean = Ocean::new(params, &grid.ops(), &land, vec![277.0; grid.n()], ROTATION);
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
        let mut ocean = Ocean::new(params, &grid.ops(), &land, vec![277.0; grid.n()], ROTATION);
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

    /// The two-layer ocean with only Ekman transport moving heat.
    fn ekman_only() -> OceanParams {
        OceanParams {
            diffusivity: 0.0,
            convective_adjustment: false,
            deep_start_limit: f64::INFINITY,
            ice: None,
            ekman: Some(EkmanParams::default()),
            ..OceanParams::default()
        }
    }

    /// A rough wind stress (N m⁻²) of up to 0.15 in each direction.
    fn rough_wind(i: usize) -> [f64; 2] {
        [0.3 * (scatter(i, 5) - 0.5), 0.3 * (scatter(i, 6) - 0.5)]
    }

    /// Lowest and highest temperature of either layer where there is ocean.
    fn layer_range(ocean: &Ocean) -> (f64, f64) {
        (0..ocean.area.len())
            .filter(|&i| ocean.area[i] > 0.0)
            .flat_map(|i| [ocean.mixed[i], ocean.deep[i]])
            .fold((f64::MAX, f64::MIN), |(a, b), t| (a.min(t), b.max(t)))
    }

    #[test]
    fn ekman_transport_conserves_heat_to_round_off_and_stays_within_its_range() {
        // Rough winds, coasts and temperatures, at a normal step and at one
        // long enough to need sub-steps.
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        let mixed: Vec<f64> = (0..grid.n())
            .map(|i| 271.0 + 30.0 * scatter(i, 3))
            .collect();
        let deep: Vec<f64> = (0..grid.n())
            .map(|i| 271.0 + 30.0 * scatter(i, 4))
            .collect();
        for (sverdrup, dt, calls) in [
            (true, 1800.0, 20),
            (true, 1.0e7, 3),
            (false, 1800.0, 20),
            (false, 1.0e7, 3),
        ] {
            let params = OceanParams {
                ekman: Some(EkmanParams {
                    sverdrup,
                    ..EkmanParams::default()
                }),
                ..ekman_only()
            };
            let mut ocean = Ocean::new(params, &grid.ops(), &land, mixed.clone(), ROTATION);
            ocean.deep.copy_from_slice(&deep);
            let (low, high) = layer_range(&ocean);
            let before = total(&ocean);
            for _ in 0..calls {
                ocean.set_wind_stress(rough_wind);
                ocean.transport(dt);
            }
            // Measured 2.5e-15, which is the sum in `total` itself: the
            // change in the sum of each cell's heat above 271 K is 6e-17 of
            // the total.
            let drift = (total(&ocean) - before) / before;
            let steps = ocean.overturning.steps;
            eprintln!("Sverdrup {sverdrup}, dt {dt}: relative drift {drift:e}, {steps} sub-steps");
            assert!(drift.abs() < 1e-14, "relative drift {drift} at dt {dt}");
            if dt > 1e6 {
                assert!(steps > 1, "{steps} sub-steps");
            }
            let (new_low, new_high) = layer_range(&ocean);
            assert!(
                new_low >= low && new_high <= high,
                "{new_low}–{new_high} from {low}–{high}"
            );
            assert!(ocean.mass_flux[EKMAN].iter().any(|f| *f != 0.0));
        }
    }

    #[test]
    fn a_uniform_ocean_stays_uniform_under_wind() {
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        let start = vec![275.0; grid.n()];
        let mut ocean = Ocean::new(
            OceanParams {
                ekman: Some(EkmanParams::default()),
                ..OceanParams::default()
            },
            &grid.ops(),
            &land,
            start.clone(),
            ROTATION,
        );
        for dt in [1800.0, 1.0e7] {
            ocean.set_wind_stress(rough_wind);
            ocean.transport(dt);
        }
        // Each cell's mass budget closes only to round-off, so the heat
        // carried in and out of a uniform cell can differ by round-off, a few
        // units in the last place of 275 K (6e-14 K) per sub-step. Measured
        // exactly 0.
        let error = (0..grid.n())
            .filter(|&i| land[i] < 1.0)
            .flat_map(|i| [ocean.mixed[i], ocean.deep[i]])
            .map(|t| (t - 275.0).abs())
            .fold(0.0, f64::max);
        eprintln!("largest departure {error:e} K");
        assert!(error < 1e-10, "{error} K");
    }

    #[test]
    fn a_uniform_eastward_stress_drives_the_analytic_ekman_flux() {
        // The paper's Eq. 4 everywhere, without the Sverdrup blend.
        let grid = Grid::new(16);
        let land = vec![0.0; grid.n()];
        let params = OceanParams {
            ekman: Some(EkmanParams {
                sverdrup: false,
                ..EkmanParams::default()
            }),
            ..ekman_only()
        };
        let mut ocean = Ocean::new(params, &grid.ops(), &land, vec![290.0; grid.n()], ROTATION);
        let tau = 0.1;
        ocean.set_wind_stress(|_| [tau, 0.0]);
        ocean.transport(1800.0);
        let epsilon = params.ekman.unwrap().friction;
        // Each link against Eq. 4 at its edge point. The scheme's stress
        // there is the mean of the two cells' east vectors, which is short of
        // the edge point's by the square of their longitude difference: the
        // error grows towards the poles, to 9e-3 at 75°.
        let mut worst = 0.0f64;
        for (k, &e) in ocean.link_edges().collect::<Vec<_>>().iter().enumerate() {
            let e = e as usize;
            let x = grid.geometry.edge_points[e];
            let latitude = fms_math::asin(x.z).to_degrees().abs();
            if !(15.0..=60.0).contains(&latitude) {
                continue;
            }
            let east = DVec3::Z.cross(x).normalize();
            let stress = tau * east;
            let f = 2.0 * ROTATION * x.z;
            let m = (epsilon * stress - f * x.cross(stress)) / (epsilon * epsilon + f * f);
            let width = grid.geometry.edge_lengths[e] * RADIUS;
            let exact = m.dot(grid.geometry.edge_normals[e]) * width;
            worst = worst.max((ocean.mass_flux[EKMAN][k] - exact).abs() / (m.length() * width));
        }
        // Measured 2.8e-3.
        eprintln!("largest error {worst:e} of |M| × width");
        assert!(worst < 5e-3, "{worst}");
        // The net northward flux across 45° N, per metre of the parallel,
        // is the Ekman transport −τ/f.
        let z0 = fms_math::sin(45f64.to_radians());
        let centers = grid.mesh.cell_centers();
        let crossing: f64 = ocean
            .link_edges()
            .zip(&ocean.mass_flux[EKMAN])
            .map(|(e, f)| {
                let [a, b] = grid.mesh.edges()[e as usize].cells;
                match (centers[a as usize].z > z0, centers[b as usize].z > z0) {
                    (false, true) => *f,
                    (true, false) => -f,
                    _ => 0.0,
                }
            })
            .sum();
        let per_metre = crossing / (std::f64::consts::TAU * RADIUS * (1.0 - z0 * z0).sqrt());
        let ekman = -tau / (2.0 * ROTATION * z0);
        // Measured 2.1% short, 0.9% of it the friction's ε²/f².
        eprintln!("{per_metre} kg/m/s against −τ/f = {ekman}");
        assert!(
            (per_metre / ekman - 1.0).abs() < 0.04,
            "{per_metre} vs {ekman}"
        );
    }

    #[test]
    fn tropical_easterlies_upwell_cold_water_on_the_equator() {
        let grid = Grid::new(16);
        let land = vec![0.0; grid.n()];
        let params = OceanParams {
            deep_start_limit: 285.0,
            ..ekman_only()
        };
        let mut ocean = Ocean::new(params, &grid.ops(), &land, vec![300.0; grid.n()], ROTATION);
        // Easterlies within 30° of the equator, strongest on it.
        let centers = grid.mesh.cell_centers();
        let easterly = |i: usize| {
            let z: f64 = centers[i].z;
            [-0.1 * (1.0 - (z / 0.5).powi(2)).max(0.0), 0.0]
        };
        ocean.set_wind_stress(easterly);
        ocean.transport(86_400.0);
        let mut equatorial = 0;
        for (i, x) in centers.iter().enumerate() {
            if x.z.abs() < 0.04 {
                equatorial += 1;
                assert!(ocean.upwelling(i) > 0.0, "cell {i}: {}", ocean.upwelling(i));
                assert!(ocean.mixed[i] < 300.0, "cell {i}: {}", ocean.mixed[i]);
                assert!(ocean.ekman_mixed_heating(i) < 0.0);
            }
        }
        assert!(equatorial > 0);
        let a = &ocean.area;
        let net: f64 = (0..grid.n()).map(|i| a[i] * ocean.ekman_heating(i)).sum();
        let gross: f64 = (0..grid.n())
            .map(|i| a[i] * ocean.ekman_heating(i).abs())
            .sum();
        eprintln!("net {net:e} W of gross {gross:e} W");
        assert!(
            gross > 0.0 && (net / gross).abs() < 1e-12,
            "{net} of {gross}"
        );
    }

    #[test]
    fn no_stress_moves_no_water() {
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        let start: Vec<f64> = (0..grid.n())
            .map(|i| 272.0 + 30.0 * scatter(i, 3))
            .collect();
        let params = OceanParams {
            ice: Some(IceParams::default()),
            ..ekman_only()
        };
        // Full ice cover takes all of the wind's stress; so does no wind.
        let mut iced = Ocean::new(params, &grid.ops(), &land, start.clone(), ROTATION);
        iced.ice.as_mut().unwrap().fraction.fill(1.0);
        let mut calm = Ocean::new(params, &grid.ops(), &land, start.clone(), ROTATION);
        let (mixed, deep) = (iced.mixed.clone(), iced.deep.clone());
        iced.set_wind_stress(rough_wind);
        calm.set_wind_stress(|_| [0.0; 2]);
        for ocean in [&mut iced, &mut calm] {
            ocean.transport(1800.0);
            assert!(ocean.mass_flux[EKMAN].iter().all(|f| *f == 0.0));
            assert_eq!(ocean.mixed, mixed);
            assert_eq!(ocean.deep, deep);
            assert!(
                (0..grid.n()).all(|i| ocean.upwelling(i) == 0.0 && ocean.ekman_heating(i) == 0.0)
            );
        }
    }

    #[test]
    fn without_ekman_the_wind_changes_nothing() {
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        let start: Vec<f64> = (0..grid.n())
            .map(|i| 272.0 + 30.0 * scatter(i, 3))
            .collect();
        let new = || {
            Ocean::new(
                OceanParams::default(),
                &grid.ops(),
                &land,
                start.clone(),
                ROTATION,
            )
        };
        let (mut windy, mut still) = (new(), new());
        for _ in 0..5 {
            windy.set_wind_stress(rough_wind);
            windy.transport(1800.0);
            still.transport(1800.0);
        }
        assert_eq!(windy.mixed, still.mixed);
        assert_eq!(windy.deep, still.deep);
    }

    #[test]
    #[should_panic(expected = "Ekman transport needs a deep layer")]
    fn ekman_transport_needs_a_deep_layer() {
        let grid = Grid::new(4);
        let params = OceanParams {
            ekman: Some(EkmanParams::default()),
            ..OceanParams::slab(50.0)
        };
        let land = vec![0.0; grid.n()];
        Ocean::new(params, &grid.ops(), &land, vec![290.0; grid.n()], ROTATION);
    }

    /// The net northward mass flux across the parallel at `z0` per metre of
    /// it (kg m⁻¹ s⁻¹), summing the links whose cells lie on either side.
    fn northward_per_metre(grid: &Grid, ocean: &Ocean, z0: f64) -> f64 {
        northward_flux_per_metre(grid, ocean, &ocean.mass_flux[EKMAN], z0)
    }

    /// As [`northward_per_metre`], of any flux on the ocean's links.
    fn northward_flux_per_metre(grid: &Grid, ocean: &Ocean, flux: &[f64], z0: f64) -> f64 {
        let centers = grid.mesh.cell_centers();
        let crossing: f64 = ocean
            .link_edges()
            .zip(flux)
            .map(|(e, f)| {
                let [a, b] = grid.mesh.edges()[e as usize].cells;
                match (centers[a as usize].z > z0, centers[b as usize].z > z0) {
                    (false, true) => *f,
                    (true, false) => -f,
                    _ => 0.0,
                }
            })
            .sum();
        crossing / (std::f64::consts::TAU * RADIUS * (1.0 - z0 * z0).sqrt())
    }

    /// An eastward stress `T sin φ cos φ` (N m⁻²): none on the equator, where
    /// its curl is `−T/R`.
    fn sverdrup_wind(grid: &Grid, t: f64) -> impl Fn(usize) -> [f64; 2] + Sync + '_ {
        move |i| {
            let z = grid.mesh.cell_centers()[i].z;
            [t * z * (1.0 - z * z).sqrt(), 0.0]
        }
    }

    #[test]
    fn the_ekman_factor_is_finite_and_exact_at_the_equator() {
        let epsilon = 1e-5;
        assert_eq!(ekman_factor(0.0, epsilon), 0.0);
        for f in [1e-12, 1e-9, -3e-8] {
            let small = ekman_factor(f, epsilon) / (f / (epsilon * epsilon));
            assert!((small - 1.0).abs() < 1e-5, "{f}: {small}");
        }
        for f in [1e-3, -2e-4] {
            assert_eq!(ekman_factor(f, epsilon), 1.0 / f);
        }
    }

    #[test]
    fn the_sverdrup_factor_is_one_over_beta_but_bounded() {
        let (epsilon, beta_earth) = (1e-5, 2.0 * ROTATION / RADIUS);
        let earth = sverdrup_factor(beta_earth, epsilon, RADIUS) * beta_earth;
        // 0.5% below 1/β on Earth's equator.
        assert!((0.99..0.996).contains(&earth), "{earth}");
        let bound = RADIUS / (2.0 * epsilon);
        for k in 0..60 {
            let beta = beta_earth * fms_math::pow(10.0, -0.1 * f64::from(k));
            let b = sverdrup_factor(beta, epsilon, RADIUS);
            assert!(b.is_finite() && b <= bound * (1.0 + 1e-15), "{beta}: {b}");
        }
        assert_eq!(sverdrup_factor(0.0, epsilon, RADIUS), 0.0);
    }

    #[test]
    fn the_stress_curl_is_the_grids_curl() {
        let grid = Grid::new(16);
        let land = vec![0.0; grid.n()];
        let mut ocean = Ocean::new(
            ekman_only(),
            &grid.ops(),
            &land,
            vec![290.0; grid.n()],
            ROTATION,
        );
        ocean.set_wind_stress(rough_wind);
        ocean.transport(1800.0);
        let curl = ocean.ekman.as_ref().unwrap().curl.as_ref().unwrap();
        // All ocean: every edge is a link.
        let mut along = vec![0.0; grid.mesh.edges().len()];
        for (e, a) in ocean.link_edges().zip(&curl.along) {
            along[e as usize] = *a;
        }
        let mut expected = vec![0.0; grid.mesh.n_vertices()];
        grid.ops().curl(&along, &mut expected);
        let scale = expected.iter().fold(0.0f64, |a, c| a.max(c.abs()));
        let error = curl
            .vertex
            .iter()
            .zip(&expected)
            .fold(0.0f64, |a, (c, e)| a.max((c - e).abs()));
        eprintln!("largest difference {error:e} of {scale:e} N/m³");
        assert!(scale > 0.0 && error <= 1e-14 * scale, "{error} of {scale}");
    }

    #[test]
    fn a_coastal_vertex_has_no_curl() {
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        let mut ocean = Ocean::new(
            ekman_only(),
            &grid.ops(),
            &land,
            vec![290.0; grid.n()],
            ROTATION,
        );
        ocean.set_wind_stress(rough_wind);
        ocean.transport(1800.0);
        let curl = ocean.ekman.as_ref().unwrap().curl.as_ref().unwrap();
        let (mut coastal, mut wet) = (0, 0);
        for (v, cells) in grid.mesh.triangles().iter().enumerate() {
            if cells.iter().any(|&c| land[c as usize] == 1.0) {
                coastal += 1;
                assert!(curl.stencil[v].is_none(), "vertex {v}");
                assert_eq!(curl.vertex[v], 0.0, "vertex {v}");
            } else {
                wet += 1;
                assert!(curl.stencil[v].is_some(), "vertex {v}");
            }
        }
        assert!(coastal > 0 && wet > 0, "{coastal} coastal, {wet} wet");
        assert!(curl.vertex.iter().any(|c| *c != 0.0));
    }

    #[test]
    fn the_wind_curl_drives_the_sverdrup_transport_across_the_equator() {
        let grid = Grid::new(16);
        let land = vec![0.0; grid.n()];
        let params = ekman_only();
        let mut ocean = Ocean::new(params, &grid.ops(), &land, vec![290.0; grid.n()], ROTATION);
        let t = 0.1;
        ocean.set_wind_stress(sverdrup_wind(&grid, t));
        ocean.transport(1800.0);
        let epsilon = params.ekman.unwrap().friction;
        let beta = 2.0 * ROTATION / RADIUS;
        let damping = epsilon / RADIUS;
        let expected = -t / (2.0 * ROTATION) * beta * beta / (beta * beta + damping * damping);
        let measured = northward_per_metre(&grid, &ocean, 0.0);
        let error = measured / expected - 1.0;
        eprintln!("{measured} kg/m/s across the equator against {expected}: {error:+.4}");
        assert!(error.abs() < 0.03, "{measured} vs {expected}");
        // Measured 0.4% short. Without the blend only the Ekman and
        // frictional flow of the edges a degree or two either side of the
        // equator cross it, and much less does (measured 13%).
        let mut plain = Ocean::new(
            OceanParams {
                ekman: Some(EkmanParams {
                    sverdrup: false,
                    ..EkmanParams::default()
                }),
                ..params
            },
            &grid.ops(),
            &land,
            vec![290.0; grid.n()],
            ROTATION,
        );
        plain.set_wind_stress(sverdrup_wind(&grid, t));
        plain.transport(1800.0);
        let without = northward_per_metre(&grid, &plain, 0.0);
        eprintln!("without Sverdrup: {without} kg/m/s");
        assert!(without.abs() < 0.25 * expected.abs(), "{without}");
    }

    #[test]
    fn away_from_the_equator_the_blend_is_pure_ekman_transport() {
        let grid = Grid::new(16);
        let land = vec![0.0; grid.n()];
        let mut ocean = Ocean::new(
            ekman_only(),
            &grid.ops(),
            &land,
            vec![290.0; grid.n()],
            ROTATION,
        );
        let tau = 0.1;
        ocean.set_wind_stress(|_| [tau, 0.0]);
        ocean.transport(1800.0);
        let mut worst = 0.0f64;
        for (k, e) in ocean.link_edges().enumerate() {
            let e = e as usize;
            let x = grid.geometry.edge_points[e];
            let latitude = fms_math::asin(x.z).to_degrees().abs();
            if !(15.0..=60.0).contains(&latitude) {
                continue;
            }
            let stress = tau * DVec3::Z.cross(x).normalize();
            let m = -x.cross(stress) / (2.0 * ROTATION * x.z);
            let width = grid.geometry.edge_lengths[e] * RADIUS;
            let exact = m.dot(grid.geometry.edge_normals[e]) * width;
            worst = worst.max((ocean.mass_flux[EKMAN][k] - exact).abs() / (m.length() * width));
        }
        // Measured 2.7e-3, from the same averaging of east vectors as Eq. 4's.
        eprintln!("largest error {worst:e} of |M| × width");
        assert!(worst < 5e-3, "{worst}");
    }

    #[test]
    fn without_rotation_only_the_friction_moves_the_water() {
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        let mixed: Vec<f64> = (0..grid.n())
            .map(|i| 271.0 + 30.0 * scatter(i, 3))
            .collect();
        let mut ocean = Ocean::new(ekman_only(), &grid.ops(), &land, mixed, 0.0);
        let before = total(&ocean);
        for _ in 0..10 {
            ocean.set_wind_stress(rough_wind);
            ocean.transport(86_400.0);
        }
        let ekman = ocean.ekman.as_ref().unwrap();
        let epsilon = ekman.friction;
        for (k, link) in ocean.links.iter().enumerate() {
            let tau = link_stress(&ekman.frames, &ekman.stress, link);
            let face = &ekman.faces[k];
            let friction = (1.0 / epsilon * tau).dot(face.normal) * link.width;
            assert_eq!(ocean.mass_flux[EKMAN][k], friction, "link {k}");
        }
        assert!(ocean.mixed.iter().chain(&ocean.deep).all(|t| t.is_finite()));
        let drift = (total(&ocean) - before) / before;
        assert!(drift.abs() < 1e-14, "{drift}");
    }

    #[test]
    fn slow_rotation_keeps_the_sverdrup_transport_bounded() {
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        let params = ekman_only();
        let epsilon = params.ekman.unwrap().friction;
        let bound = RADIUS / (2.0 * epsilon);
        for rotation in [ROTATION, 1e-5, 1e-6, 1e-7, 1e-8] {
            let mut ocean = Ocean::new(params, &grid.ops(), &land, vec![290.0; grid.n()], rotation);
            ocean.deep.iter_mut().for_each(|t| *t = 280.0);
            let before = total(&ocean);
            ocean.set_wind_stress(sverdrup_wind(&grid, 0.1));
            ocean.transport(86_400.0);
            let ekman = ocean.ekman.as_ref().unwrap();
            let curl = ekman.curl.as_ref().unwrap();
            let mut largest = 0.0f64;
            for face in &ekman.faces {
                let edge_curl = 0.5
                    * face
                        .vertices
                        .iter()
                        .map(|&v| curl.vertex[v as usize])
                        .sum::<f64>();
                let part = (face.blend[2] * edge_curl * face.north.dot(face.normal)).abs();
                assert!(part <= bound * edge_curl.abs(), "Ω {rotation}: {part}");
                largest = largest.max(face.blend[2] / bound);
            }
            assert!(ocean.mass_flux[EKMAN].iter().all(|f| f.is_finite()));
            let drift = (total(&ocean) - before) / before;
            eprintln!("Ω {rotation:e}: largest w B is {largest:.3} of R/(2ε), drift {drift:e}");
            assert!(drift.abs() < 1e-14, "Ω {rotation}: {drift}");
        }
    }

    /// Eddy transport alone, on an ocean that otherwise keeps still.
    fn eddies_only() -> OceanParams {
        OceanParams {
            ekman: None,
            eddies: Some(EddyParams::default()),
            ..ekman_only()
        }
    }

    #[test]
    fn eddy_transport_conserves_heat_to_round_off_and_stays_within_its_range() {
        // Rough coasts and temperatures, alone and with a rough wind's Ekman
        // transport, at a normal step and at one long enough for sub-steps.
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        let mixed: Vec<f64> = (0..grid.n())
            .map(|i| 271.0 + 30.0 * scatter(i, 3))
            .collect();
        let deep: Vec<f64> = (0..grid.n())
            .map(|i| 271.0 + 30.0 * scatter(i, 4))
            .collect();
        for (ekman, dt, calls) in [
            (false, 1800.0, 20),
            (false, 1.0e7, 3),
            (true, 1800.0, 20),
            (true, 1.0e7, 3),
        ] {
            let params = OceanParams {
                ekman: ekman.then(EkmanParams::default),
                ..eddies_only()
            };
            let mut ocean = Ocean::new(params, &grid.ops(), &land, mixed.clone(), ROTATION);
            ocean.deep.copy_from_slice(&deep);
            let (low, high) = layer_range(&ocean);
            let before = total(&ocean);
            for _ in 0..calls {
                ocean.set_wind_stress(rough_wind);
                ocean.transport(dt);
            }
            let drift = (total(&ocean) - before) / before;
            let steps = ocean.overturning.steps;
            eprintln!("Ekman {ekman}, dt {dt}: relative drift {drift:e}, {steps} sub-steps");
            assert!(drift.abs() < 1e-14, "relative drift {drift} at dt {dt}");
            if dt > 1e6 {
                assert!(steps > 1, "{steps} sub-steps");
            }
            let (new_low, new_high) = layer_range(&ocean);
            assert!(
                new_low >= low && new_high <= high,
                "{new_low}–{new_high} from {low}–{high}"
            );
            assert!(ocean.mass_flux[EDDY].iter().any(|f| *f != 0.0));
            assert_eq!(ocean.mass_flux[EKMAN].iter().any(|f| *f != 0.0), ekman);
        }
    }

    #[test]
    fn level_isotherms_drive_no_eddy_transport() {
        // Each layer uniform, the mixed one warmer: nothing to flatten.
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        let mut ocean = Ocean::new(
            eddies_only(),
            &grid.ops(),
            &land,
            vec![290.0; grid.n()],
            ROTATION,
        );
        ocean.deep.fill(278.0);
        ocean.transport(1800.0);
        assert!(ocean.mass_flux[EDDY].iter().all(|f| *f == 0.0));
        assert!(ocean.mixed.iter().all(|t| *t == 290.0));
        assert!(ocean.deep.iter().all(|t| *t == 278.0));
        assert!((0..grid.n()).all(|i| ocean.eddy_heating(i) == 0.0 && ocean.upwelling(i) == 0.0));
    }

    #[test]
    fn the_eddy_flux_runs_down_the_isotherm_slope() {
        // Both layers warm towards the south by `A sin φ`, 10 K apart: the
        // mean temperature's northward gradient is `−A cos φ / R` and the
        // slope that over `10 K / separation`, well under its limit.
        let grid = Grid::new(32);
        let land = vec![0.0; grid.n()];
        let (amplitude, contrast) = (20.0, 10.0);
        let mixed = grid.each(|z| 290.0 - amplitude * z);
        let mut ocean = Ocean::new(eddies_only(), &grid.ops(), &land, mixed, ROTATION);
        for i in 0..grid.n() {
            ocean.deep[i] = ocean.mixed[i] - contrast;
        }
        ocean.transport(1800.0);
        let (p, e) = (ocean.params, EddyParams::default());
        // kg m⁻¹ s⁻¹ per unit of the gradient's cosine.
        let strength = p.density * e.diffusivity * (amplitude / RADIUS) * e.separation / contrast;
        assert!(amplitude / RADIUS * e.separation / contrast < 0.1 * e.max_slope);
        let mut worst = 0.0f64;
        for (k, edge) in ocean.link_edges().enumerate() {
            let edge = edge as usize;
            let (x, n) = (
                grid.geometry.edge_points[edge],
                grid.geometry.edge_normals[edge],
            );
            // −ρ κ s along the normal, with ∇T = −(A/R) ẑ's tangent part.
            let exact = strength * n.z * grid.geometry.edge_lengths[edge] * RADIUS;
            let scale = strength * grid.geometry.edge_lengths[edge] * RADIUS;
            let _ = x;
            worst = worst.max((ocean.mass_flux[EDDY][k] - exact).abs() / scale);
        }
        eprintln!("largest error {worst:e} of ρ κ s × width");
        assert!(worst < 2e-3, "{worst}");
        // Towards the cold north in the mixed layer, at 45° N.
        let z0 = std::f64::consts::FRAC_1_SQRT_2;
        let across = northward_flux_per_metre(&grid, &ocean, &ocean.mass_flux[EDDY], z0);
        let expected = strength * z0;
        eprintln!("{across} kg/m/s across 45° N, expected {expected}");
        assert!((across / expected - 1.0).abs() < 0.03, "{across}");
        // Heat leaves the warm hemisphere's columns for the cold one's, at
        // ρ c κ · separation · ∇T whatever the layers' difference. Measured
        // 6.4% more: the mixed layer's upstream cell is warmer and the deep
        // layer's colder by the gradient across one cell, 0.64 K of the 10.
        let south: f64 = (0..grid.n())
            .filter(|&i| grid.mesh.cell_centers()[i].z < 0.0)
            .map(|i| ocean.area[i] * ocean.eddy_heating(i))
            .sum();
        let equator = std::f64::consts::TAU * RADIUS;
        let transport =
            p.density * p.specific_heat * e.diffusivity * e.separation * amplitude / RADIUS;
        eprintln!(
            "{south:e} W leaves the south, expected {:e}",
            -transport * equator
        );
        let excess = south / (-transport * equator) - 1.0;
        assert!((0.0..0.1).contains(&excess), "{south}");
    }

    #[test]
    fn the_eddy_flux_is_limited_by_the_slope_and_the_contrast() {
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        // Steep everywhere, and layers that have mixed.
        let mixed: Vec<f64> = (0..grid.n())
            .map(|i| 271.0 + 30.0 * scatter(i, 3))
            .collect();
        let mut ocean = Ocean::new(eddies_only(), &grid.ops(), &land, mixed, ROTATION);
        ocean.deep.copy_from_slice(&ocean.mixed.clone());
        ocean.transport(1800.0);
        let (p, e) = (ocean.params, EddyParams::default());
        let most = p.density * e.diffusivity * e.max_slope;
        let mut limited = 0;
        for (link, f) in ocean.links.iter().zip(&ocean.mass_flux[EDDY]) {
            assert!(f.is_finite() && f.abs() <= most * link.width, "{f}");
            limited += usize::from(f.abs() == most * link.width);
        }
        assert!(
            limited > ocean.links.len() / 2,
            "{limited} links at the limit"
        );
        // A gentle slope over layers closer than the floor is taken over the
        // floor: 0.25 K, not the 0.1 K between them.
        let grid = Grid::new(32);
        let land = vec![0.0; grid.n()];
        let mixed = grid.each(|z| 290.0 - 0.02 * z);
        let mut ocean = Ocean::new(eddies_only(), &grid.ops(), &land, mixed, ROTATION);
        for i in 0..grid.n() {
            ocean.deep[i] = ocean.mixed[i] - 0.1;
        }
        ocean.transport(1800.0);
        let z0 = std::f64::consts::FRAC_1_SQRT_2;
        let across = northward_flux_per_metre(&grid, &ocean, &ocean.mass_flux[EDDY], z0);
        let expected =
            p.density * e.diffusivity * (0.02 / RADIUS) * e.separation / e.min_contrast * z0;
        assert!(
            (across / expected - 1.0).abs() < 0.03,
            "{across} of {expected}"
        );
    }

    #[test]
    fn eddies_flatten_the_isotherms() {
        // A warm equator and cold poles in both layers, the deep one colder.
        let grid = Grid::new(16);
        let land = vec![0.0; grid.n()];
        let mixed = grid.each(|z| 300.0 - 25.0 * z * z);
        let mut ocean = Ocean::new(eddies_only(), &grid.ops(), &land, mixed.clone(), ROTATION);
        let deep = grid.each(|z| 282.0 - 8.0 * z * z);
        ocean.deep.copy_from_slice(&deep);
        let column = |ocean: &Ocean, i: usize| 0.25 * ocean.mixed[i] + 0.75 * ocean.deep[i];
        let start: Vec<f64> = (0..grid.n()).map(|i| column(&ocean, i)).collect();
        let before = total(&ocean);
        for _ in 0..360 {
            ocean.transport(86_400.0);
        }
        assert!(((total(&ocean) - before) / before).abs() < 1e-14);
        let centers = grid.mesh.cell_centers();
        let mean = |f: &dyn Fn(usize) -> f64, keep: &dyn Fn(f64) -> bool| {
            let (mut sum, mut area) = (0.0, 0.0);
            for i in (0..grid.n()).filter(|&i| keep(centers[i].z.abs())) {
                sum += ocean.area[i] * f(i);
                area += ocean.area[i];
            }
            sum / area
        };
        let (polar, tropical) = (|z: f64| z > 0.8, |z: f64| z < 0.3);
        // Heat has moved from the tropics' columns to the high latitudes'.
        let high = mean(&|i| column(&ocean, i) - start[i], &polar);
        let low = mean(&|i| column(&ocean, i) - start[i], &tropical);
        // The surface warms at high latitudes and the deep layer cools in
        // the tropics, as in the paper's Fig. 2.
        let surface = mean(&|i| ocean.mixed[i] - mixed[i], &polar);
        let depth = mean(&|i| ocean.deep[i] - deep[i], &tropical);
        eprintln!(
            "columns {high:+.3} K polar, {low:+.3} K tropical; polar surface {surface:+.3} K, tropical depth {depth:+.3} K"
        );
        assert!(high > 0.05 && low < -0.01, "{high}, {low}");
        // Measured +0.79 K and −0.006 K: the deep layer's own gradient is
        // small in the tropics.
        assert!(surface > 0.1 && depth < -1e-3, "{surface}, {depth}");
    }

    #[test]
    fn each_flux_accounts_for_its_own_heating() {
        // With both on, each one's heating sums to nothing over the ocean,
        // and together they are the change in every column's heat.
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        let mixed: Vec<f64> = (0..grid.n())
            .map(|i| 285.0 + 15.0 * scatter(i, 3))
            .collect();
        let params = OceanParams {
            ekman: Some(EkmanParams::default()),
            ..eddies_only()
        };
        let mut ocean = Ocean::new(params, &grid.ops(), &land, mixed, ROTATION);
        for i in 0..grid.n() {
            ocean.deep[i] = 275.0 + 8.0 * scatter(i, 4);
        }
        let dt = 3600.0;
        let before: Vec<f64> = (0..grid.n()).map(|i| ocean.heat(i)).collect();
        ocean.set_wind_stress(rough_wind);
        ocean.transport(dt);
        let (mut ekman, mut eddy, mut gross) = (0.0, 0.0, 0.0);
        for i in (0..grid.n()).filter(|&i| ocean.area[i] > 0.0) {
            let (a, b) = (ocean.ekman_heating(i), ocean.eddy_heating(i));
            ekman += ocean.area[i] * a;
            eddy += ocean.area[i] * b;
            gross += ocean.area[i] * (a.abs() + b.abs());
            let change = (ocean.heat(i) - before[i]) / dt;
            assert!(a != 0.0 && b != 0.0, "cell {i}: {a}, {b}");
            // The heat is 2e11 J m⁻², so its change is known to about
            // 1e-5 W m⁻² × 1e11 / dt.
            assert!(
                (change - (a + b)).abs() < 1e-2 + 1e-9 * change.abs(),
                "cell {i}"
            );
        }
        eprintln!("net {ekman:e} and {eddy:e} W of {gross:e} gross");
        assert!(ekman.abs() < 1e-12 * gross && eddy.abs() < 1e-12 * gross);
    }

    #[test]
    fn eddy_coefficients_scale_with_the_rotation() {
        assert_eq!(eddy_scale(EARTH_ROTATION), 1.0);
        assert_eq!(eddy_scale(-EARTH_ROTATION), 1.0);
        assert_eq!(eddy_scale(EARTH_ROTATION / 2.0), 2.0);
        assert_eq!(eddy_scale(4.0 * EARTH_ROTATION), 0.25);
        assert_eq!(eddy_scale(EARTH_ROTATION / 100.0), MAX_EDDY_SCALE);
        assert_eq!(eddy_scale(0.0), MAX_EDDY_SCALE);
        // A planet turning half as fast diffuses and stirs as Earth would
        // with both coefficients doubled, and not at all differently if the
        // scaling is off.
        let grid = Grid::new(16);
        let land = coasts(grid.n());
        let mixed: Vec<f64> = (0..grid.n())
            .map(|i| 285.0 + 15.0 * scatter(i, 3))
            .collect();
        let run = |params: OceanParams, rotation: f64| {
            let mut ocean = Ocean::new(params, &grid.ops(), &land, mixed.clone(), rotation);
            for i in 0..grid.n() {
                ocean.deep[i] = 275.0 + 8.0 * scatter(i, 4);
            }
            for _ in 0..5 {
                ocean.transport(86_400.0);
            }
            (ocean.mixed, ocean.deep)
        };
        let base = OceanParams {
            diffusivity: 2000.0,
            ..eddies_only()
        };
        let doubled = OceanParams {
            diffusivity: 4000.0,
            eddies: Some(EddyParams {
                diffusivity: 4000.0,
                ..EddyParams::default()
            }),
            ..base
        };
        let unscaled = OceanParams {
            rotation_scaling: false,
            ..base
        };
        let earth = run(base, EARTH_ROTATION);
        assert_eq!(
            run(base, EARTH_ROTATION / 2.0),
            run(doubled, EARTH_ROTATION)
        );
        assert_eq!(run(unscaled, EARTH_ROTATION / 2.0), earth);
        assert_ne!(run(base, EARTH_ROTATION / 2.0), earth);
    }

    #[test]
    #[should_panic(expected = "eddy transport needs a deep layer")]
    fn eddy_transport_needs_a_deep_layer() {
        let grid = Grid::new(4);
        let params = OceanParams {
            eddies: Some(EddyParams::default()),
            ..OceanParams::slab(50.0)
        };
        let land = vec![0.0; grid.n()];
        Ocean::new(params, &grid.ops(), &land, vec![290.0; grid.n()], ROTATION);
    }
}
