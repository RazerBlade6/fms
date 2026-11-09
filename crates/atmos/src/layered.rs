//! Layered dry primitive-equation model in sigma coordinates: the dynamical
//! core of the GCM, with any number of layers ([`Levels`]) on the TRiSK
//! C-grid. With two layers it is the Phillips two-level model.
//!
//! # Vertical structure
//!
//! `σ = p / p_s` runs from 0 at the model top (p = 0) to 1 at the surface.
//! Layer `k` (0 at the top) spans `[σ_{k−½}, σ_{k+½}]`, with pressure
//! thickness `δpₖ = δσₖ p_s`. Interface `j` (between layers `j` and `j + 1`)
//! is at `σⱼ = σ_{j+½}`.
//!
//! # Prognostic variables
//!
//! - `p_s` surface pressure at cells (Pa)
//! - `uₖ` normal wind at edges in each layer (m s⁻¹)
//! - `Θₖ = p_s θₖ`, mass-weighted potential temperature at cells (Pa K)
//! - optionally `Qₖ = p_s qₖ`, mass-weighted specific humidity at cells (Pa),
//!   carried like `Θ` (see [`crate::transport`] for the semi-implicit
//!   stepper's positive-definite scheme)
//!
//! Mass, `Θ` and `Q` are linear invariants of the equations, so RK4 conserves
//! their totals to round-off.
//!
//! # Equations (per layer k)
//!
//! ```text
//! Fₖ = δpₖ,ₑ uₖ,   Dₖ = ∇·Fₖ
//! ∂p_s/∂t      = −Σₖ Dₖ
//! Wⱼ           = σⱼ Σ_{k>j} Dₖ − (1 − σⱼ) Σ_{k≤j} Dₖ     (p_s σ̇ at interface j, downward)
//! δσₖ ∂Θₖ/∂t   = −∇·(Fₖ θₖ,ₑ) + W_{k−1} θ̂_{k−1} − Wₖ θ̂ₖ
//! δσₖ ∂Qₖ/∂t   = −∇·(Fₖ qₖ,ₑ) + W_{k−1} q̂_{k−1} − Wₖ q̂ₖ   (q̂: see below)
//! ∂uₖ/∂t       = qₖ F⊥ₖ − ∇(Kₖ + Φₖ) − c_p θₖ,ₑ ∇Pₖ + Vₖ
//! Vₖ           = −[Wₖ,ₑ (u_{k+1} − uₖ) + W_{k−1},ₑ (uₖ − u_{k−1})] / (2 δpₖ,ₑ)
//! ```
//!
//! with `qₖ = (ζₖ + f)/δpₖ` the potential vorticity, `Kₖ` the TRiSK kinetic
//! energy, `θ̂ⱼ = (θⱼ + θ_{j+1})/2` the interface potential temperature, and
//! terms at the model top and the surface (`W_{−1}`, `W_{N−1}`) zero.
//!
//! The vertical discretization follows Arakawa & Suarez (1983). Layer Exner
//! functions `Pₖ` are the mass-weighted means of `Π = (p/p₀)^κ` over each
//! layer ([`Levels::layer_exner`]). The geopotentials satisfy the hydrostatic
//! relation across each interface, with the same `θ̂` as the vertical
//! advection of θ,
//!
//! ```text
//! Φⱼ − Φ_{j+1} = c_p θ̂ⱼ (P_{j+1} − Pⱼ)
//! ```
//!
//! and their mass-weighted mean is fixed by integrating `Φ dp` by parts,
//!
//! ```text
//! Σₖ δσₖ Φₖ = Φ_s + c_p κ Σₖ δσₖ θₖ Pₖ
//! ```
//!
//! These are exactly the choices that make the spatial scheme conserve total
//! energy `Σ Aᵢ [Σₖ δpₖ (Kₖ + c_p θₖ Pₖ) + p_s Φ_s] / g`; for constant θ they
//! reduce to the exact layer-mean geopotential. With them an isentropic
//! atmosphere at rest over any terrain stays exactly at rest.
//!
//! Humidity is concentrated low in each layer, so across an interface it is
//! carried at the interface values of the two layers' own profiles
//! ([`crate::radiation::interface_humidity_ratios`]): the mean of the two
//! here, the upwind one in the semi-implicit stepper's positive scheme.

use fms_grid::Operators;
use fms_grid::operators::MIN_CHUNK;
use rayon::prelude::*;

use crate::diffusion::{Hyperdiffusion, HyperdiffusionOptions};
use crate::dissipation::KineticHeating;
use crate::forcing::{EnergyTerms, Forcing};
use crate::helmholtz::dot;
pub use crate::levels::Levels;
use crate::levels::{Eigen, MAX_LAYERS};
use crate::radiation::{RadiationParams, interface_humidity_ratios};

/// Reference pressure for potential temperature (Pa).
pub const REFERENCE_PRESSURE: f64 = 1.0e5;

/// Thermodynamic constants of the (dry) air.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Gas {
    /// Specific gas constant `R` (J kg⁻¹ K⁻¹).
    pub gas_constant: f64,
    /// Specific heat at constant pressure `c_p` (J kg⁻¹ K⁻¹).
    pub heat_capacity: f64,
}

impl Gas {
    /// Earth's dry air.
    pub const EARTH_AIR: Self = Self {
        gas_constant: 287.04,
        heat_capacity: 1004.64,
    };

    /// `κ = R / c_p`.
    pub fn kappa(&self) -> f64 {
        self.gas_constant / self.heat_capacity
    }

    /// Exner function `Π = (p/p₀)^κ`.
    pub fn exner(&self, p: f64) -> f64 {
        fms_math::pow(p / REFERENCE_PRESSURE, self.kappa())
    }
}

/// The prognostic state. Layer 0 is the top.
#[derive(Clone, Debug, PartialEq)]
pub struct LayeredState {
    /// Surface pressure at cells (Pa).
    pub ps: Vec<f64>,
    /// Normal wind at edges, per layer (m s⁻¹).
    pub u: Vec<Vec<f64>>,
    /// Mass-weighted potential temperature `p_s θ` at cells, per layer (Pa K).
    pub theta_ps: Vec<Vec<f64>>,
    /// Mass-weighted specific humidity `p_s q` at cells, per layer (Pa);
    /// no layers for a dry atmosphere.
    pub q_ps: Vec<Vec<f64>>,
}

impl LayeredState {
    /// A dry state from surface pressure, winds and potential temperatures
    /// (one of each per layer).
    ///
    /// # Panics
    /// If the winds and temperatures have different numbers of layers.
    pub fn new(ps: Vec<f64>, u: Vec<Vec<f64>>, theta: Vec<Vec<f64>>) -> Self {
        assert_eq!(u.len(), theta.len(), "one wind and one θ per layer");
        let theta_ps = theta
            .into_iter()
            .map(|t| t.iter().zip(&ps).map(|(t, p)| t * p).collect())
            .collect();
        Self {
            ps,
            u,
            theta_ps,
            q_ps: Vec::new(),
        }
    }

    /// The state with specific humidity `q` per layer (kg kg⁻¹).
    ///
    /// # Panics
    /// If there isn't one humidity per layer and cell.
    pub fn with_humidity(mut self, q: Vec<Vec<f64>>) -> Self {
        assert_eq!(q.len(), self.layers(), "one humidity per layer");
        for q in &q {
            assert_eq!(q.len(), self.ps.len(), "one humidity per cell");
        }
        self.q_ps = q
            .into_iter()
            .map(|q| q.iter().zip(&self.ps).map(|(q, p)| q * p).collect())
            .collect();
        self
    }

    /// Number of layers.
    pub fn layers(&self) -> usize {
        self.u.len()
    }

    /// Whether the state carries humidity.
    pub fn is_moist(&self) -> bool {
        !self.q_ps.is_empty()
    }

    /// Specific humidity of layer `k` (kg kg⁻¹).
    ///
    /// # Panics
    /// If the state is dry.
    pub fn humidity(&self, k: usize) -> Vec<f64> {
        self.q_ps[k]
            .iter()
            .zip(&self.ps)
            .map(|(q, p)| q / p)
            .collect()
    }

    /// The same state without humidity.
    pub fn dry(&self) -> Self {
        Self {
            ps: self.ps.clone(),
            u: self.u.clone(),
            theta_ps: self.theta_ps.clone(),
            q_ps: Vec::new(),
        }
    }

    /// Potential temperature of layer `k` (K).
    pub fn theta(&self, k: usize) -> Vec<f64> {
        self.theta_ps[k]
            .iter()
            .zip(&self.ps)
            .map(|(t, p)| t / p)
            .collect()
    }

    /// A dry state of zeros.
    pub(crate) fn dry_zeros(layers: usize, n_cells: usize, n_edges: usize) -> Self {
        Self {
            ps: vec![0.0; n_cells],
            u: vec![vec![0.0; n_edges]; layers],
            theta_ps: vec![vec![0.0; n_cells]; layers],
            q_ps: Vec::new(),
        }
    }

    /// Zeros shaped like `like` (moist if it is).
    pub(crate) fn zeros_like(like: &Self) -> Self {
        let z = |v: &Vec<Vec<f64>>| v.iter().map(|x| vec![0.0; x.len()]).collect();
        Self {
            ps: vec![0.0; like.ps.len()],
            u: z(&like.u),
            theta_ps: z(&like.theta_ps),
            q_ps: z(&like.q_ps),
        }
    }

    /// Every field: `p_s`, then the winds, `Θ` and `Q` of each layer.
    pub(crate) fn parts(&self) -> Vec<&[f64]> {
        let mut parts: Vec<&[f64]> = vec![&self.ps];
        parts.extend(self.u.iter().map(Vec::as_slice));
        parts.extend(self.theta_ps.iter().map(Vec::as_slice));
        parts.extend(self.q_ps.iter().map(Vec::as_slice));
        parts
    }

    /// [`parts`](Self::parts), mutably.
    pub(crate) fn parts_mut(&mut self) -> Vec<&mut [f64]> {
        let mut parts: Vec<&mut [f64]> = vec![&mut self.ps];
        parts.extend(self.u.iter_mut().map(Vec::as_mut_slice));
        parts.extend(self.theta_ps.iter_mut().map(Vec::as_mut_slice));
        parts.extend(self.q_ps.iter_mut().map(Vec::as_mut_slice));
        parts
    }

    /// `self = base + a · tendency`.
    fn set_axpy(&mut self, base: &Self, a: f64, tendency: &Self) {
        for ((out, x), y) in self
            .parts_mut()
            .into_iter()
            .zip(base.parts())
            .zip(tendency.parts())
        {
            out.par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .zip(x)
                .zip(y)
                .for_each(|((o, &x), &y)| *o = x + a * y);
        }
    }
}

/// Linear gravity-wave terms about a resting, flat reference state, in the
/// variables `y = (p_s', θ₀', …, θ_{N−1}')` and layer divergences
/// `δ = (δ₀, …, δ_{N−1})`: `∂y/∂t = B δ` and `∂uₖ/∂t = −∇(J y)ₖ`. Their
/// product `A = J B` is the gravity-wave matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct Linearization {
    /// `J` (`N × (N+1)`): rows `Mₖ`, columns `(p_s', θ₀', …)`, where
    /// `Mₖ = Φₖ + c_p θ̄ₖ Pₖ` is the layer's pressure-gradient potential.
    pub j: Vec<Vec<f64>>,
    /// `B` (`(N+1) × N`): rows `(p_s', θ₀', …)`, columns `(δ₀, …)`.
    pub b: Vec<Vec<f64>>,
}

/// Global integrals for checking conservation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayeredDiagnostics {
    /// Total mass `Σ Aᵢ p_s / g` (kg).
    pub mass: f64,
    /// Total mass-weighted potential temperature `Σₖ Σ Aᵢ δpₖ θₖ / g` (kg K).
    pub theta: f64,
    /// Total kinetic energy (J).
    pub kinetic_energy: f64,
    /// Total energy: kinetic + enthalpy + surface potential (J). Conserved by
    /// the spatial scheme.
    pub energy: f64,
    /// Total water vapour `Σₖ Σ Aᵢ δσₖ Qₖ / g` (kg; 0 when dry).
    pub water: f64,
}

/// Diagnostic fields for one layer.
#[derive(Clone, Debug, Default)]
struct LayerWork {
    dp_edge: Vec<f64>,
    flux: Vec<f64>,
    divergence: Vec<f64>,
    theta: Vec<f64>,
    theta_edge: Vec<f64>,
    exner: Vec<f64>,
    geopotential: Vec<f64>,
    kinetic: Vec<f64>,
    bernoulli: Vec<f64>,
    vorticity: Vec<f64>,
    dp_vertex: Vec<f64>,
    q_vertex: Vec<f64>,
    q_edge: Vec<f64>,
    pv_flux: Vec<f64>,
    gradient: Vec<f64>,
    exner_gradient: Vec<f64>,
    theta_flux_divergence: Vec<f64>,
}

#[derive(Clone, Debug, Default)]
struct Work {
    layers: Vec<LayerWork>,
    /// `p_s σ̇` at each interface (cells and edges), and `θ̂` there.
    w: Vec<Vec<f64>>,
    w_edge: Vec<Vec<f64>>,
    theta_interface: Vec<Vec<f64>>,
    /// Per cell, every layer's Exner function then geopotential.
    column: Vec<f64>,
    cell_scratch: Vec<f64>,
    edge_scratch: Vec<f64>,
}

/// The layered dry dynamical core on one mesh.
pub struct LayeredModel<'a> {
    ops: Operators<'a>,
    gravity: f64,
    gas: Gas,
    levels: Levels,
    /// Coriolis parameter at vertices (s⁻¹).
    coriolis: Vec<f64>,
    /// Surface geopotential `g · z_s` at cells (m² s⁻²).
    surface_geopotential: Vec<f64>,
    work: Work,
    stages: [LayeredState; 4],
    stage_state: LayeredState,
    diffusion: Option<Hyperdiffusion>,
    forcing: Option<Box<dyn Forcing + 'a>>,
    /// Profile within the layers that sets humidity at the interfaces.
    humidity_profile: RadiationParams,
    /// Whether diffusion returns lost kinetic energy as heat.
    heat_dissipation: bool,
}

impl<'a> LayeredModel<'a> {
    /// A model with gravity (m s⁻²), gas constants, vertical `levels`,
    /// Coriolis parameter at vertices (s⁻¹) and surface geopotential at
    /// cells (m² s⁻²).
    ///
    /// # Panics
    /// If field lengths don't match the mesh.
    pub fn new(
        ops: Operators<'a>,
        gravity: f64,
        gas: Gas,
        levels: Levels,
        coriolis: Vec<f64>,
        surface_geopotential: Vec<f64>,
    ) -> Self {
        let m = ops.mesh();
        let (nc, ne, nv) = (m.n_cells(), m.n_edges(), m.n_vertices());
        assert_eq!(coriolis.len(), nv, "one Coriolis value per vertex");
        assert_eq!(surface_geopotential.len(), nc, "one geopotential per cell");
        let n = levels.count();
        let layer = || LayerWork {
            dp_edge: vec![0.0; ne],
            flux: vec![0.0; ne],
            divergence: vec![0.0; nc],
            theta: vec![0.0; nc],
            theta_edge: vec![0.0; ne],
            exner: vec![0.0; nc],
            geopotential: vec![0.0; nc],
            kinetic: vec![0.0; nc],
            bernoulli: vec![0.0; nc],
            vorticity: vec![0.0; nv],
            dp_vertex: vec![0.0; nv],
            q_vertex: vec![0.0; nv],
            q_edge: vec![0.0; ne],
            pv_flux: vec![0.0; ne],
            gradient: vec![0.0; ne],
            exner_gradient: vec![0.0; ne],
            theta_flux_divergence: vec![0.0; nc],
        };
        Self {
            ops,
            gravity,
            gas,
            levels,
            coriolis,
            surface_geopotential,
            work: Work {
                layers: (0..n).map(|_| layer()).collect(),
                w: vec![vec![0.0; nc]; n - 1],
                w_edge: vec![vec![0.0; ne]; n - 1],
                theta_interface: vec![vec![0.0; nc]; n - 1],
                column: vec![0.0; 2 * n * nc],
                cell_scratch: vec![0.0; nc],
                edge_scratch: vec![0.0; ne],
            },
            stages: std::array::from_fn(|_| LayeredState::dry_zeros(n, nc, ne)),
            stage_state: LayeredState::dry_zeros(n, nc, ne),
            diffusion: None,
            forcing: None,
            humidity_profile: RadiationParams::default(),
            heat_dissipation: false,
        }
    }

    /// The model with humidity profiles within the layers following the
    /// radiation's (its lapse rate and humidity shape): use the forcing's
    /// radiation parameters. Default: [`RadiationParams::default`].
    pub fn with_humidity_profile(mut self, radiation: RadiationParams) -> Self {
        self.humidity_profile = radiation;
        self
    }

    fn parameters(&self) -> Parameters<'_, 'a> {
        Parameters {
            ops: &self.ops,
            gas: self.gas,
            levels: &self.levels,
            gravity: self.gravity,
            coriolis: &self.coriolis,
            surface_geopotential: &self.surface_geopotential,
            humidity_profile: &self.humidity_profile,
        }
    }

    /// Per interface and cell: humidity at the interface relative to each
    /// adjacent layer's mean (`[upper layer at its bottom, lower layer at
    /// its top]`; [`crate::radiation::interface_humidity_ratios`]).
    pub fn interface_ratios(&self, state: &LayeredState) -> Vec<[Vec<f64>; 2]> {
        interface_ratios(&self.parameters(), state)
    }

    /// The model with hyperdiffusion ([`crate::diffusion`]). Without it the
    /// core is inviscid and collapses to the grid scale once eddies mature.
    pub fn with_hyperdiffusion(mut self, options: HyperdiffusionOptions) -> Self {
        let mut diffusion = Hyperdiffusion::new(&self.ops, options);
        if self.heat_dissipation {
            diffusion.set_heating(Some(self.kinetic_heating()));
        }
        self.diffusion = Some(diffusion);
        self
    }

    /// The model with the kinetic energy hyperdiffusion removes returned as
    /// heat where it is lost ([`crate::dissipation`]), so it conserves total
    /// energy. Off by default (Held–Suarez leaves it out).
    pub fn with_dissipative_heating(mut self) -> Self {
        self.heat_dissipation = true;
        let heating = self.kinetic_heating();
        if let Some(d) = &mut self.diffusion {
            d.set_heating(Some(heating));
        }
        self
    }

    /// A [`KineticHeating`] for this model's mesh, gas and layers.
    pub fn kinetic_heating(&self) -> KineticHeating {
        KineticHeating::new(&self.ops, self.gas, self.levels.clone())
    }

    /// The hyperdiffusion, if any.
    pub fn hyperdiffusion(&self) -> Option<&Hyperdiffusion> {
        self.diffusion.as_ref()
    }

    /// The model with external forcing ([`crate::forcing`]), e.g.
    /// [`crate::HeldSuarez`].
    pub fn with_forcing(mut self, forcing: Box<dyn Forcing + 'a>) -> Self {
        self.forcing = Some(forcing);
        self
    }

    /// Whether the model has external forcing.
    pub fn has_forcing(&self) -> bool {
        self.forcing.is_some()
    }

    /// The external forcing, if any.
    pub fn forcing(&self) -> Option<&(dyn Forcing + 'a)> {
        self.forcing.as_deref()
    }

    /// Tell the forcing a step of `dt` seconds has completed with `state`
    /// ([`Forcing::advance`]). The steppers call this; call it yourself only
    /// when stepping by other means.
    pub fn advance_forcing(&mut self, state: &LayeredState, dt: f64) {
        if let Some(f) = &mut self.forcing {
            f.advance(state, dt);
        }
    }

    /// The operators.
    pub fn operators(&self) -> &Operators<'a> {
        &self.ops
    }

    /// Gas constants.
    pub fn gas(&self) -> Gas {
        self.gas
    }

    /// Gravity (m s⁻²).
    pub fn gravity(&self) -> f64 {
        self.gravity
    }

    /// The vertical levels.
    pub fn levels(&self) -> &Levels {
        &self.levels
    }

    /// Number of layers.
    pub fn layers(&self) -> usize {
        self.levels.count()
    }

    /// Coriolis parameter at vertices (s⁻¹).
    pub fn coriolis(&self) -> &[f64] {
        &self.coriolis
    }

    /// Surface geopotential at cells (m² s⁻²).
    pub fn surface_geopotential(&self) -> &[f64] {
        &self.surface_geopotential
    }

    /// Layer pressure-thickness fractions `δσₖ`.
    pub fn dsigma(&self) -> &[f64] {
        self.levels.dsigma()
    }

    /// Layer Exner functions `Pₖ` for surface pressure `ps`.
    pub fn layer_exner(&self, ps: f64) -> Vec<f64> {
        self.levels.layer_exner(self.gas, ps)
    }

    /// Layer geopotentials `Φₖ` (m² s⁻²) for one column with surface
    /// pressure `ps`, layer potential temperatures `theta` and surface
    /// geopotential `phi_s`.
    pub fn geopotentials_at(&self, ps: f64, theta: &[f64], phi_s: f64) -> Vec<f64> {
        let exner = self.layer_exner(ps);
        let mut phi = vec![0.0; self.layers()];
        geopotentials(self.gas, &self.levels, theta, &exner, phi_s, &mut phi);
        phi
    }

    /// The vertical-structure matrix of linear gravity waves about a resting,
    /// flat state with surface pressure `ps` and layer potential temperatures
    /// `theta`: `∂²δ/∂t² = −∇²(A δ)` for the layer divergences `δ`.
    ///
    /// `A = J B`, where `B` gives the tendencies of `(p_s', θₖ')` from the
    /// divergences (linearized continuity, thermodynamics and interface
    /// fluxes) and `J` the change of each layer's pressure-gradient
    /// potential `Mₖ = Φₖ + c_p θ̄ₖ Pₖ`. Its eigenvalues are `−c²` of the
    /// external and internal modes ([`LayeredModel::gravity_wave_speeds`]);
    /// the semi-implicit scheme treats exactly these terms implicitly.
    pub fn gravity_wave_matrix(&self, ps: f64, theta: &[f64]) -> Vec<Vec<f64>> {
        let Linearization { j, b } = self.linearization(ps, theta);
        let n = self.layers();
        (0..n)
            .map(|r| {
                (0..n)
                    .map(|c| (0..=n).map(|i| j[r][i] * b[i][c]).sum())
                    .collect()
            })
            .collect()
    }

    /// The linear gravity-wave terms about a resting, flat state with surface
    /// pressure `ps` and layer potential temperatures `theta` (see
    /// [`Linearization`]).
    pub fn linearization(&self, ps: f64, theta: &[f64]) -> Linearization {
        let n = self.layers();
        let ds = self.dsigma();
        let sigma = self.levels.interior();
        // The interface mass flux per unit p̄ from each layer's divergence:
        // Wⱼ/p̄ = Σₗ m[j][l] δₗ.
        let m: Vec<Vec<f64>> = sigma
            .iter()
            .enumerate()
            .map(|(j, &s)| {
                (0..n)
                    .map(|l| if l > j { s * ds[l] } else { -(1.0 - s) * ds[l] })
                    .collect()
            })
            .collect();
        let th_hat: Vec<f64> = (0..n - 1)
            .map(|j| 0.5 * (theta[j] + theta[j + 1]))
            .collect();
        // Continuity: ∂p_s/∂t = −p̄ Σ δσₗ δₗ. Thermodynamics, in advective
        // form (continuity removed): ∂θₖ'/∂t = −[W_{k−1} (θₖ − θ̂_{k−1})
        // + Wₖ (θ̂ₖ − θₖ)] / (δσₖ p̄).
        let mut b = vec![vec![0.0; n]; n + 1];
        for l in 0..n {
            b[0][l] = -ps * ds[l];
        }
        for k in 0..n {
            for l in 0..n {
                let mut v = 0.0;
                if k > 0 {
                    v += (theta[k] - th_hat[k - 1]) * m[k - 1][l];
                }
                if k + 1 < n {
                    v += (th_hat[k] - theta[k]) * m[k][l];
                }
                b[k + 1][l] = -v / ds[k];
            }
        }
        // J by central differences of Mₖ = Φₖ + c_p θ̄ₖ Pₖ (θ̄ frozen in the
        // second term, as in the pressure-gradient force).
        let cp = self.gas.heat_capacity;
        let potential = |ps: f64, th: &[f64]| -> Vec<f64> {
            let phi = self.geopotentials_at(ps, th, 0.0);
            let exner = self.layer_exner(ps);
            (0..n).map(|k| phi[k] + cp * theta[k] * exner[k]).collect()
        };
        let (ep, et) = (1.0e-5 * ps.abs().max(1.0), 1.0e-3);
        let mut j = vec![vec![0.0; n + 1]; n];
        let (hi, lo) = (potential(ps + ep, theta), potential(ps - ep, theta));
        for k in 0..n {
            j[k][0] = (hi[k] - lo[k]) / (2.0 * ep);
        }
        for l in 0..n {
            let (mut up, mut down) = (theta.to_vec(), theta.to_vec());
            up[l] += et;
            down[l] -= et;
            let (hi, lo) = (potential(ps, &up), potential(ps, &down));
            for k in 0..n {
                j[k][1 + l] = (hi[k] - lo[k]) / (2.0 * et);
            }
        }
        Linearization { j, b }
    }

    /// Phase speeds (m s⁻¹) of linear gravity waves about a resting, flat
    /// state, fastest (external) first: `c² = −eigenvalues` of
    /// [`LayeredModel::gravity_wave_matrix`].
    pub fn gravity_wave_speeds(&self, ps: f64, theta: &[f64]) -> Vec<f64> {
        Eigen::of(&self.gravity_wave_matrix(ps, theta))
            .values
            .iter()
            .map(|l| (-l).max(0.0).sqrt())
            .collect()
    }

    /// Geopotential of each layer (m² s⁻²) at every cell.
    // Cells and layers index several arrays together.
    #[allow(clippy::needless_range_loop)]
    pub fn geopotential(&self, state: &LayeredState) -> Vec<Vec<f64>> {
        let n = self.layers();
        let mut phi = vec![vec![0.0; state.ps.len()]; n];
        let (mut exner, mut theta, mut out) =
            ([0.0; MAX_LAYERS], [0.0; MAX_LAYERS], [0.0; MAX_LAYERS]);
        for i in 0..state.ps.len() {
            let ps = state.ps[i];
            self.levels.layer_exner_into(self.gas, ps, &mut exner[..n]);
            for k in 0..n {
                theta[k] = state.theta_ps[k][i] / ps;
            }
            geopotentials(
                self.gas,
                &self.levels,
                &theta[..n],
                &exner[..n],
                self.surface_geopotential[i],
                &mut out[..n],
            );
            for k in 0..n {
                phi[k][i] = out[k];
            }
        }
        phi
    }

    /// Adiabatic, inviscid time tendencies of `state` into `out` (no
    /// diffusion: see [`diffusion_tendency`](Self::diffusion_tendency)).
    pub fn tendency(&mut self, state: &LayeredState, out: &mut LayeredState) {
        let params = Parameters {
            ops: &self.ops,
            gas: self.gas,
            levels: &self.levels,
            gravity: self.gravity,
            coriolis: &self.coriolis,
            surface_geopotential: &self.surface_geopotential,
            humidity_profile: &self.humidity_profile,
        };
        compute_tendency(params, &mut self.work, state, out);
    }

    /// Hyperdiffusion tendencies of `state` into `out` (zero without
    /// hyperdiffusion; `out.ps` is always zero).
    pub fn diffusion_tendency(&mut self, state: &LayeredState, out: &mut LayeredState) {
        for x in out.parts_mut() {
            x.fill(0.0);
        }
        if let Some(d) = &mut self.diffusion {
            d.add_tendency(&self.ops, state, out);
        }
    }

    /// Forcing tendencies of `state` into `out` (zero without forcing;
    /// `out.ps` is always zero).
    pub fn forcing_tendency(&mut self, state: &LayeredState, out: &mut LayeredState) {
        for x in out.parts_mut() {
            x.fill(0.0);
        }
        if let Some(f) = &mut self.forcing {
            f.add_tendency(state, out);
        }
    }

    /// Rate of change of total energy (W) caused by a tendency of the winds
    /// and `Θ` at fixed surface pressure, such as
    /// [`diffusion_tendency`](Self::diffusion_tendency): the exact derivative
    /// of the energy in [`diagnostics`](Self::diagnostics) along it, split
    /// into the part through the winds and the part through `Θ`.
    /// `tendency.ps` is ignored.
    pub fn energy_tendency(&self, state: &LayeredState, tendency: &LayeredState) -> EnergyTerms {
        let ops = &self.ops;
        let r2 = ops.radius() * ops.radius();
        let g = ops.geometry();
        let ds = self.dsigma();
        let n = self.layers();
        let edges = ops.mesh().edges();
        let edge_areas = ops.edge_areas();
        let mut rate = EnergyTerms::default();
        let gravity = self.gravity;
        let cp = self.gas.heat_capacity;
        // Layer Exner functions, once per cell.
        let mut exner = vec![0.0; n * state.ps.len()];
        exner
            .par_chunks_mut(n)
            .with_min_len(MIN_CHUNK / 16)
            .zip(state.ps.par_iter())
            .for_each(|(e, &ps)| self.levels.layer_exner_into(self.gas, ps, e));
        for (k, &dsk) in ds.iter().enumerate() {
            // ∂E/∂uₑ = 2 Aₑ δpₑ uₑ / g, with δpₑ the mean of the two cells.
            let kinetic: Vec<f64> = (0..edges.len())
                .into_par_iter()
                .with_min_len(MIN_CHUNK)
                .map(|e| {
                    let [a, b] = edges[e].cells;
                    let dp = dsk * (state.ps[a as usize] + state.ps[b as usize]);
                    edge_areas[e] * dp * state.u[k][e] / gravity
                })
                .collect();
            rate.kinetic += dot(&kinetic, &tendency.u[k]);
            // ∂E/∂Θₖ = Aᵢ δσₖ c_p Πₖ / g.
            let enthalpy: Vec<f64> = (0..state.ps.len())
                .into_par_iter()
                .with_min_len(MIN_CHUNK)
                .map(|i| g.cell_areas[i] * r2 * dsk * cp * exner[i * n + k] / gravity)
                .collect();
            rate.thermal += dot(&enthalpy, &tendency.theta_ps[k]);
        }
        rate
    }

    /// Advance `state` by `dt` seconds with classic fourth-order Runge–Kutta,
    /// including hyperdiffusion and forcing (and [`Forcing::advance`] at the
    /// end).
    pub fn step(&mut self, state: &mut LayeredState, dt: f64) {
        if self.stage_state.is_moist() != state.is_moist() {
            self.stages = std::array::from_fn(|_| LayeredState::zeros_like(state));
            self.stage_state = LayeredState::zeros_like(state);
        }
        let Self {
            ops,
            gas,
            levels,
            gravity,
            coriolis,
            surface_geopotential,
            work,
            stages,
            stage_state,
            diffusion,
            forcing,
            humidity_profile,
            heat_dissipation: _,
        } = self;
        let params = Parameters {
            ops,
            gas: *gas,
            levels,
            gravity: *gravity,
            coriolis,
            surface_geopotential,
            humidity_profile,
        };
        let mut full = |x: &LayeredState, out: &mut LayeredState| {
            compute_tendency(params, work, x, out);
            if let Some(d) = diffusion.as_mut() {
                d.add_tendency(ops, x, out);
            }
            if let Some(f) = forcing.as_mut() {
                f.add_tendency(x, out);
            }
        };
        let [k1, k2, k3, k4] = stages;
        full(state, k1);
        stage_state.set_axpy(state, 0.5 * dt, k1);
        full(stage_state, k2);
        stage_state.set_axpy(state, 0.5 * dt, k2);
        full(stage_state, k3);
        stage_state.set_axpy(state, dt, k3);
        full(stage_state, k4);
        let tendencies = [k1.parts(), k2.parts(), k3.parts(), k4.parts()];
        for (j, x) in state.parts_mut().into_iter().enumerate() {
            let [a, b, c, d] = [
                tendencies[0][j],
                tendencies[1][j],
                tendencies[2][j],
                tendencies[3][j],
            ];
            x.par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, x)| *x += dt / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]));
        }
        if let Some(f) = forcing.as_mut() {
            f.advance(state, dt);
        }
    }

    /// Relative vorticity of layer `k` at vertices (s⁻¹).
    pub fn vorticity(&self, state: &LayeredState, k: usize) -> Vec<f64> {
        let mut zeta = vec![0.0; self.ops.mesh().n_vertices()];
        self.ops.curl(&state.u[k], &mut zeta);
        zeta
    }

    /// Mass, θ, kinetic and total energy of `state`.
    pub fn diagnostics(&self, state: &LayeredState) -> LayeredDiagnostics {
        let ops = &self.ops;
        let (m, g) = (ops.mesh(), ops.geometry());
        let r2 = ops.radius() * ops.radius();
        let ds = self.dsigma();
        let n = self.layers();
        let theta: Vec<Vec<f64>> = (0..n).map(|k| state.theta(k)).collect();
        let mut kinetic = vec![vec![0.0; m.n_cells()]; n];
        for (u, ke) in state.u.iter().zip(&mut kinetic) {
            ops.kinetic_energy(u, ke);
        }
        let (mut mass, mut th, mut ke, mut energy) = (0.0, 0.0, 0.0, 0.0);
        let water = self.water(state);
        let mut exner = [0.0; MAX_LAYERS];
        for i in 0..m.n_cells() {
            let a = g.cell_areas[i] * r2 / self.gravity;
            let ps = state.ps[i];
            self.levels.layer_exner_into(self.gas, ps, &mut exner[..n]);
            mass += a * ps;
            energy += a * ps * self.surface_geopotential[i];
            for k in 0..n {
                let dp = ds[k] * ps;
                th += a * dp * theta[k][i];
                ke += a * dp * kinetic[k][i];
                energy +=
                    a * dp * (kinetic[k][i] + self.gas.heat_capacity * theta[k][i] * exner[k]);
            }
        }
        LayeredDiagnostics {
            mass,
            theta: th,
            kinetic_energy: ke,
            energy,
            water,
        }
    }

    /// Total water vapour of `state` (kg; 0 when dry).
    pub fn water(&self, state: &LayeredState) -> f64 {
        self.water_tendency(state)
    }

    /// Global integral `Σₖ Σ Aᵢ δσₖ xₖ / g` of a humidity field or tendency
    /// (`x.q_ps`): water in kg, or its rate of change in kg s⁻¹.
    pub fn water_tendency(&self, x: &LayeredState) -> f64 {
        if !x.is_moist() {
            return 0.0;
        }
        let ops = &self.ops;
        let a = &ops.geometry().cell_areas;
        let scale = ops.radius() * ops.radius() / self.gravity;
        self.dsigma()
            .iter()
            .zip(&x.q_ps)
            .map(|(ds, q)| ds * scale * dot(a, q))
            .sum()
    }

    /// An RK4 time step (s) at Courant number `courant` for the fastest
    /// signal: the external (Lamb) gravity wave plus the wind,
    /// `courant · min dₑ / (max |u| + √(γ R T_max))`.
    pub fn stable_dt(&self, state: &LayeredState, courant: f64) -> f64 {
        let g = self.ops.geometry();
        let min_d = g.dual_edge_lengths.iter().fold(f64::MAX, |a, &b| a.min(b)) * self.ops.radius();
        let max_u = state
            .u
            .iter()
            .flatten()
            .fold(0.0f64, |a, &u| a.max(u.abs()));
        let mut max_t = 0.0f64;
        for (k, theta_ps) in state.theta_ps.iter().enumerate() {
            for (tp, &ps) in theta_ps.iter().zip(&state.ps) {
                max_t = max_t.max(tp / ps * self.layer_exner(ps)[k]);
            }
        }
        let gamma = 1.0 / (1.0 - self.gas.kappa());
        courant * min_d / (max_u + (gamma * self.gas.gas_constant * max_t).sqrt())
    }
}

#[derive(Clone, Copy)]
struct Parameters<'p, 'a> {
    ops: &'p Operators<'a>,
    gas: Gas,
    levels: &'p Levels,
    gravity: f64,
    coriolis: &'p [f64],
    surface_geopotential: &'p [f64],
    humidity_profile: &'p RadiationParams,
}

fn interface_ratios(params: &Parameters, state: &LayeredState) -> Vec<[Vec<f64>; 2]> {
    let n = params.levels.count();
    let ratios: Vec<Vec<[f64; 2]>> = (0..state.ps.len())
        .into_par_iter()
        .with_min_len(MIN_CHUNK / 16)
        .map(|i| {
            let ps = state.ps[i];
            let mut exner = [0.0; MAX_LAYERS];
            params
                .levels
                .layer_exner_into(params.gas, ps, &mut exner[..n]);
            let t: Vec<f64> = (0..n)
                .map(|k| state.theta_ps[k][i] / ps * exner[k])
                .collect();
            interface_humidity_ratios(
                params.humidity_profile,
                params.gas,
                params.levels,
                ps,
                params.gravity,
                &t,
            )
        })
        .collect();
    (0..n - 1)
        .map(|j| {
            [
                ratios.iter().map(|r| r[j][0]).collect(),
                ratios.iter().map(|r| r[j][1]).collect(),
            ]
        })
        .collect()
}

/// Layer geopotentials `Φₖ` at one cell (see the module docs), from layer
/// potential temperatures `theta` and Exner functions `exner`, into `out`.
fn geopotentials(
    gas: Gas,
    levels: &Levels,
    theta: &[f64],
    exner: &[f64],
    phi_s: f64,
    out: &mut [f64],
) {
    let n = levels.count();
    let ds = levels.dsigma();
    let cp = gas.heat_capacity;
    // X = c_p κ Σ δσₖ θₖ Pₖ.
    let mut sum = 0.0;
    for k in 0..n {
        sum += ds[k] * theta[k] * exner[k];
    }
    let x = cp * gas.kappa() * sum;
    // Offsets above the lowest layer: cₖ = Φₖ − Φ_{N−1}, from the
    // hydrostatic jumps across the interfaces below k.
    let mut offset = [0.0; MAX_LAYERS];
    for j in (0..n - 1).rev() {
        let theta_hat = 0.5 * (theta[j] + theta[j + 1]);
        let jump = cp * theta_hat * (exner[j + 1] - exner[j]);
        offset[j] = offset[j + 1] + jump;
    }
    // The lowest layer from the mass-weighted mean: Σ δσₖ Φₖ = Φ_s + X.
    let mut mean_offset = 0.0;
    for k in 0..n {
        mean_offset += ds[k] * offset[k];
    }
    let lowest = phi_s + x - mean_offset;
    for k in 0..n {
        out[k] = lowest + offset[k];
    }
}

fn compute_tendency(
    params: Parameters,
    work: &mut Work,
    state: &LayeredState,
    out: &mut LayeredState,
) {
    let Parameters {
        ops,
        gas,
        levels,
        coriolis,
        surface_geopotential,
        ..
    } = params;
    let n = levels.count();
    let ds = levels.dsigma();
    let sigma = levels.interior();
    let cp = gas.heat_capacity;
    let Work {
        layers,
        w,
        w_edge,
        theta_interface,
        column,
        cell_scratch,
        edge_scratch,
    } = work;

    // Thickness fluxes, divergence and θ per layer.
    ops.cell_to_edge(&state.ps, edge_scratch);
    for (k, l) in layers.iter_mut().enumerate() {
        l.dp_edge
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(edge_scratch.par_iter())
            .for_each(|(d, &p)| *d = ds[k] * p);
        l.flux
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(e, f)| *f = l.dp_edge[e] * state.u[k][e]);
        ops.divergence(&l.flux, &mut l.divergence);
        l.theta
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, t)| *t = state.theta_ps[k][i] / state.ps[i]);
    }
    // Exner functions and geopotentials need every layer's θ.
    {
        let layers = &*layers;
        column
            .par_chunks_mut(2 * n)
            .with_min_len(MIN_CHUNK / 16)
            .enumerate()
            .for_each(|(i, c)| {
                let (exner, phi) = c.split_at_mut(n);
                levels.layer_exner_into(gas, state.ps[i], exner);
                let mut theta = [0.0; MAX_LAYERS];
                for (k, l) in layers.iter().enumerate() {
                    theta[k] = l.theta[i];
                }
                geopotentials(
                    gas,
                    levels,
                    &theta[..n],
                    exner,
                    surface_geopotential[i],
                    phi,
                );
            });
    }
    for (k, l) in layers.iter_mut().enumerate() {
        let column = &*column;
        l.exner
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(l.geopotential.par_iter_mut())
            .enumerate()
            .for_each(|(i, (e, g))| {
                *e = column[2 * n * i + k];
                *g = column[2 * n * i + n + k];
            });
    }
    // Interface mass fluxes Wⱼ = σⱼ Σ_{k>j} Dₖ − (1 − σⱼ) Σ_{k≤j} Dₖ, and θ̂ⱼ.
    for (j, &s) in sigma.iter().enumerate() {
        let layers = &*layers;
        w[j].par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(theta_interface[j].par_iter_mut())
            .enumerate()
            .for_each(|(i, (w, t))| {
                let mut above = layers[0].divergence[i];
                for l in &layers[1..=j] {
                    above += l.divergence[i];
                }
                let mut below = layers[j + 1].divergence[i];
                for l in &layers[j + 2..] {
                    below += l.divergence[i];
                }
                *w = s * below - (1.0 - s) * above;
                *t = 0.5 * (layers[j].theta[i] + layers[j + 1].theta[i]);
            });
        ops.cell_to_edge(&w[j], &mut w_edge[j]);
    }

    // Surface pressure.
    {
        let layers = &*layers;
        out.ps
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, d)| {
                let mut total = layers[0].divergence[i];
                for l in &layers[1..] {
                    total += l.divergence[i];
                }
                *d = -total;
            });
    }

    for k in 0..n {
        let l = &mut layers[k];
        // Thermodynamics: δσₖ ∂Θₖ/∂t = −∇·(Fₖ θₖ,ₑ) + W_{k−1} θ̂_{k−1} − Wₖ θ̂ₖ.
        ops.cell_to_edge(&l.theta, &mut l.theta_edge);
        edge_scratch
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(e, f)| *f = l.flux[e] * l.theta_edge[e]);
        ops.divergence(edge_scratch, &mut l.theta_flux_divergence);
        let (w, theta_interface) = (&*w, &*theta_interface);
        out.theta_ps[k]
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, d)| {
                let mut vertical = 0.0;
                if k > 0 {
                    vertical += w[k - 1][i] * theta_interface[k - 1][i];
                }
                if k + 1 < n {
                    vertical -= w[k][i] * theta_interface[k][i];
                }
                *d = (-l.theta_flux_divergence[i] + vertical) / ds[k];
            });

        // Momentum: PV flux, Bernoulli gradient, Exner-gradient term.
        ops.kinetic_energy(&state.u[k], &mut l.kinetic);
        l.bernoulli
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, b)| *b = l.kinetic[i] + l.geopotential[i]);
        ops.gradient(&l.bernoulli, &mut l.gradient);
        ops.gradient(&l.exner, &mut l.exner_gradient);

        ops.curl(&state.u[k], &mut l.vorticity);
        cell_scratch
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(state.ps.par_iter())
            .for_each(|(d, &p)| *d = ds[k] * p);
        ops.cell_to_vertex(cell_scratch, &mut l.dp_vertex);
        l.q_vertex
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(v, q)| *q = (l.vorticity[v] + coriolis[v]) / l.dp_vertex[v]);
        ops.vertex_to_edge(&l.q_vertex, &mut l.q_edge);
        ops.pv_flux(&l.q_edge, &l.flux, &mut l.pv_flux);
    }

    // Humidity, centred like θ (only where both the state and the output
    // carry it): δσₖ ∂Qₖ/∂t = −∇·(Fₖ qₖ,ₑ) + W_{k−1} q̂_{k−1} − Wₖ q̂ₖ.
    if state.is_moist() && out.is_moist() {
        let ratio = interface_ratios(&params, state);
        let q: Vec<Vec<f64>> = (0..n)
            .map(|k| {
                state.q_ps[k]
                    .par_iter()
                    .with_min_len(MIN_CHUNK)
                    .zip(state.ps.par_iter())
                    .map(|(q, p)| q / p)
                    .collect()
            })
            .collect();
        // q̂ⱼ: the mean of the two layers' profiles at the interface.
        let q_hat: Vec<Vec<f64>> = (0..n - 1)
            .map(|j| {
                let (q, r) = (&q, &ratio[j]);
                (0..state.ps.len())
                    .into_par_iter()
                    .with_min_len(MIN_CHUNK)
                    .map(|i| 0.5 * (r[0][i] * q[j][i] + r[1][i] * q[j + 1][i]))
                    .collect()
            })
            .collect();
        for k in 0..n {
            let l = &mut layers[k];
            ops.cell_to_edge(&q[k], &mut l.theta_edge);
            edge_scratch
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(e, f)| *f = l.flux[e] * l.theta_edge[e]);
            ops.divergence(edge_scratch, cell_scratch);
            let (w, q_hat, cells) = (&*w, &q_hat, &*cell_scratch);
            out.q_ps[k]
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, d)| {
                    let mut vertical = 0.0;
                    if k > 0 {
                        vertical += w[k - 1][i] * q_hat[k - 1][i];
                    }
                    if k + 1 < n {
                        vertical -= w[k][i] * q_hat[k][i];
                    }
                    *d = (-cells[i] + vertical) / ds[k];
                });
            // Restore θ at edges for the momentum equation below.
            ops.cell_to_edge(&l.theta, &mut l.theta_edge);
        }
    }

    // Assemble the wind tendencies, including vertical advection.
    let layers = &*layers;
    let w_edge = &*w_edge;
    for (k, out_u) in out.u.iter_mut().enumerate() {
        let l = &layers[k];
        let u = &state.u;
        out_u
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(e, d)| {
                let mut vertical = 0.0;
                if k + 1 < n {
                    vertical += w_edge[k][e] * (u[k + 1][e] - u[k][e]);
                }
                if k > 0 {
                    vertical += w_edge[k - 1][e] * (u[k][e] - u[k - 1][e]);
                }
                *d = l.pv_flux[e]
                    - l.gradient[e]
                    - cp * l.theta_edge[e] * l.exner_gradient[e]
                    - vertical / (2.0 * l.dp_edge[e]);
            });
    }
}
