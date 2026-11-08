//! Two-layer dry primitive-equation model in sigma coordinates: the
//! dynamical core of the GCM (the Phillips two-level model on the TRiSK
//! C-grid).
//!
//! # Vertical structure
//!
//! `σ = p / p_s` runs from 0 at the model top (p = 0) to 1 at the surface.
//! Layer 0 (upper) spans `[0, σᵢ]`, layer 1 (lower) `[σᵢ, 1]`, so the pressure
//! thicknesses are `δp₀ = σᵢ p_s` and `δp₁ = (1 − σᵢ) p_s`. The default
//! `σᵢ = ½` gives equal-mass layers.
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
//! ∂p_s/∂t      = −(D₀ + D₁)
//! W            = σᵢ D₁ − (1 − σᵢ) D₀           (p_s σ̇ at the interface, downward)
//! δσₖ ∂Θₖ/∂t   = −∇·(Fₖ θₖ,ₑ) ∓ W θ̂             (− for the upper layer)
//! δσₖ ∂Qₖ/∂t   = −∇·(Fₖ qₖ,ₑ) ∓ W q̂             (q̂: see below)
//! ∂uₖ/∂t       = qₖ F⊥ₖ − ∇(Kₖ + Φₖ) − c_p θₖ,ₑ ∇Pₖ + Vₖ
//! V₀ = −Wₑ (u₁ − u₀) / (2 δp₀,ₑ),   V₁ = Wₑ (u₀ − u₁) / (2 δp₁,ₑ)
//! ```
//!
//! with `qₖ = (ζₖ + f)/δpₖ` the potential vorticity, `Kₖ` the TRiSK kinetic
//! energy and `θ̂ = (θ₀ + θ₁)/2` the interface potential temperature.
//!
//! The vertical discretization follows Arakawa & Suarez (1983). Layer Exner
//! functions are the mass-weighted means of `Π = (p/p₀)^κ` over each layer:
//!
//! ```text
//! P₀ = Πᵢ / (1 + κ),   P₁ = (p_s Π_s − pᵢ Πᵢ) / ((1 + κ) δp₁)
//! ```
//!
//! Humidity is concentrated low in each layer, so across the interface it is
//! carried at the interface values of the layers' own profiles
//! ([`crate::radiation::interface_humidity_ratios`]): the mean of the two
//! here, the upwind one in the semi-implicit stepper's positive scheme.
//!
//! The geopotentials are
//!
//! ```text
//! X  = (c_p κ / p_s) Σₖ δpₖ θₖ Pₖ
//! Φ₁ = Φ_s + X − σᵢ c_p θ̂ (P₁ − P₀)
//! Φ₀ = Φ₁ + c_p θ̂ (P₁ − P₀)
//! ```
//!
//! These are exactly the choices that make the spatial scheme conserve total
//! energy `Σ Aᵢ [Σₖ δpₖ (Kₖ + c_p θₖ Pₖ) + p_s Φ_s] / g`; for constant θ they
//! reduce to the exact layer-mean geopotential. With them an isentropic
//! atmosphere at rest over any terrain stays exactly at rest.

use fms_geo::Operators;
use fms_geo::operators::MIN_CHUNK;
use rayon::prelude::*;

use crate::diffusion::{Hyperdiffusion, HyperdiffusionOptions};
use crate::forcing::{EnergyTerms, Forcing};
use crate::helmholtz::dot;
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

/// The prognostic state. Index 0 is the upper layer, 1 the lower.
#[derive(Clone, Debug, PartialEq)]
pub struct TwoLayerState {
    /// Surface pressure at cells (Pa).
    pub ps: Vec<f64>,
    /// Normal wind at edges, per layer (m s⁻¹).
    pub u: [Vec<f64>; 2],
    /// Mass-weighted potential temperature `p_s θ` at cells, per layer (Pa K).
    pub theta_ps: [Vec<f64>; 2],
    /// Mass-weighted specific humidity `p_s q` at cells, per layer (Pa);
    /// empty for a dry atmosphere.
    pub q_ps: [Vec<f64>; 2],
}

impl TwoLayerState {
    /// A dry state from surface pressure, winds and potential temperatures.
    pub fn new(ps: Vec<f64>, u: [Vec<f64>; 2], theta: [Vec<f64>; 2]) -> Self {
        let theta_ps = theta.map(|t| t.iter().zip(&ps).map(|(t, p)| t * p).collect());
        Self {
            ps,
            u,
            theta_ps,
            q_ps: [Vec::new(), Vec::new()],
        }
    }

    /// The state with specific humidity `q` per layer (kg kg⁻¹).
    ///
    /// # Panics
    /// If a layer's length doesn't match the cells.
    pub fn with_humidity(mut self, q: [Vec<f64>; 2]) -> Self {
        for q in &q {
            assert_eq!(q.len(), self.ps.len(), "one humidity per cell");
        }
        self.q_ps = q.map(|q| q.iter().zip(&self.ps).map(|(q, p)| q * p).collect());
        self
    }

    /// Whether the state carries humidity.
    pub fn is_moist(&self) -> bool {
        !self.q_ps[0].is_empty()
    }

    /// Specific humidity of layer `k` (kg kg⁻¹; empty when dry).
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
            q_ps: [Vec::new(), Vec::new()],
            ..self.clone()
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

    fn dry_zeros(n_cells: usize, n_edges: usize) -> Self {
        Self::new(
            vec![0.0; n_cells],
            [vec![0.0; n_edges], vec![0.0; n_edges]],
            [vec![0.0; n_cells], vec![0.0; n_cells]],
        )
    }

    /// Zeros shaped like `like` (moist if it is).
    pub(crate) fn zeros_like(like: &Self) -> Self {
        let z = |v: &Vec<f64>| vec![0.0; v.len()];
        Self {
            ps: z(&like.ps),
            u: [z(&like.u[0]), z(&like.u[1])],
            theta_ps: [z(&like.theta_ps[0]), z(&like.theta_ps[1])],
            q_ps: [z(&like.q_ps[0]), z(&like.q_ps[1])],
        }
    }

    pub(crate) fn parts(&self) -> [&[f64]; 7] {
        [
            &self.ps,
            &self.u[0],
            &self.u[1],
            &self.theta_ps[0],
            &self.theta_ps[1],
            &self.q_ps[0],
            &self.q_ps[1],
        ]
    }

    pub(crate) fn parts_mut(&mut self) -> [&mut [f64]; 7] {
        let [u0, u1] = &mut self.u;
        let [t0, t1] = &mut self.theta_ps;
        let [q0, q1] = &mut self.q_ps;
        [&mut self.ps, u0, u1, t0, t1, q0, q1]
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
/// variables `y = (p_s', θ₀', θ₁')` and layer divergences `δ = (δ₀, δ₁)`:
/// `∂y/∂t = B δ` and `∂uₖ/∂t = −∇(J y)ₖ`. Their product `A = J B` is the
/// gravity-wave matrix.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Linearization {
    /// `J`: rows `(M₀, M₁)`, columns `(p_s', θ₀', θ₁')`, where
    /// `Mₖ = Φₖ + c_p θ̄ₖ Pₖ` is the layer's pressure-gradient potential.
    pub j: [[f64; 3]; 2],
    /// `B`: rows `(p_s', θ₀', θ₁')`, columns `(δ₀, δ₁)`.
    pub b: [[f64; 2]; 3],
}

/// Global integrals for checking conservation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoLayerDiagnostics {
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
    layers: [LayerWork; 2],
    /// `p_s σ̇` at the interface (cells and edges).
    w: Vec<f64>,
    w_edge: Vec<f64>,
    theta_interface: Vec<f64>,
    cell_scratch: Vec<f64>,
    edge_scratch: Vec<f64>,
}

/// The two-layer dry dynamical core on one mesh.
pub struct TwoLayer<'a> {
    ops: Operators<'a>,
    gravity: f64,
    gas: Gas,
    /// Interface σ.
    sigma: f64,
    /// Coriolis parameter at vertices (s⁻¹).
    coriolis: Vec<f64>,
    /// Surface geopotential `g · z_s` at cells (m² s⁻²).
    surface_geopotential: Vec<f64>,
    work: Work,
    stages: [TwoLayerState; 4],
    stage_state: TwoLayerState,
    diffusion: Option<Hyperdiffusion>,
    forcing: Option<Box<dyn Forcing + 'a>>,
    /// Profile within the layers that sets humidity at the interface.
    humidity_profile: RadiationParams,
}

impl<'a> TwoLayer<'a> {
    /// A model with gravity (m s⁻²), gas constants, interface `sigma` (½ for
    /// equal-mass layers), Coriolis parameter at vertices (s⁻¹) and surface
    /// geopotential at cells (m² s⁻²).
    ///
    /// # Panics
    /// If `sigma` is not strictly between 0 and 1, or field lengths don't
    /// match the mesh.
    pub fn new(
        ops: Operators<'a>,
        gravity: f64,
        gas: Gas,
        sigma: f64,
        coriolis: Vec<f64>,
        surface_geopotential: Vec<f64>,
    ) -> Self {
        assert!(sigma > 0.0 && sigma < 1.0, "interface σ must be in (0, 1)");
        let m = ops.mesh();
        let (nc, ne, nv) = (m.n_cells(), m.n_edges(), m.n_vertices());
        assert_eq!(coriolis.len(), nv, "one Coriolis value per vertex");
        assert_eq!(surface_geopotential.len(), nc, "one geopotential per cell");
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
            sigma,
            coriolis,
            surface_geopotential,
            work: Work {
                layers: [layer(), layer()],
                w: vec![0.0; nc],
                w_edge: vec![0.0; ne],
                theta_interface: vec![0.0; nc],
                cell_scratch: vec![0.0; nc],
                edge_scratch: vec![0.0; ne],
            },
            stages: std::array::from_fn(|_| TwoLayerState::dry_zeros(nc, ne)),
            stage_state: TwoLayerState::dry_zeros(nc, ne),
            diffusion: None,
            forcing: None,
            humidity_profile: RadiationParams::default(),
        }
    }

    /// The model with humidity profiles within the layers following the
    /// radiation's (its lapse rate and humidity shape): use the forcing's
    /// radiation parameters. Default: [`RadiationParams::default`].
    pub fn with_humidity_profile(mut self, radiation: RadiationParams) -> Self {
        self.humidity_profile = radiation;
        self
    }

    /// Per cell: humidity at the interface relative to each layer's mean
    /// (`[upper at its bottom, lower at its top]`;
    /// [`crate::radiation::interface_humidity_ratios`]).
    pub fn interface_ratios(&self, state: &TwoLayerState) -> [Vec<f64>; 2] {
        interface_ratios(
            &Parameters {
                ops: &self.ops,
                gas: self.gas,
                sigma: self.sigma,
                gravity: self.gravity,
                coriolis: &self.coriolis,
                surface_geopotential: &self.surface_geopotential,
                humidity_profile: &self.humidity_profile,
            },
            state,
        )
    }

    /// The model with hyperdiffusion ([`crate::diffusion`]). Without it the
    /// core is inviscid and collapses to the grid scale once eddies mature.
    pub fn with_hyperdiffusion(mut self, options: HyperdiffusionOptions) -> Self {
        self.diffusion = Some(Hyperdiffusion::new(&self.ops, options));
        self
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
    pub fn advance_forcing(&mut self, state: &TwoLayerState, dt: f64) {
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

    /// Interface σ.
    pub fn sigma(&self) -> f64 {
        self.sigma
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
    pub fn dsigma(&self) -> [f64; 2] {
        [self.sigma, 1.0 - self.sigma]
    }

    /// Layer Exner functions `[P₀, P₁]` for surface pressure `ps`.
    pub fn layer_exner(&self, ps: f64) -> [f64; 2] {
        layer_exner(self.gas, self.sigma, ps)
    }

    /// Layer geopotentials `[Φ₀, Φ₁]` (m² s⁻²) for one column with surface
    /// pressure `ps`, layer potential temperatures `theta` and surface
    /// geopotential `phi_s`.
    pub fn geopotentials_at(&self, ps: f64, theta: [f64; 2], phi_s: f64) -> [f64; 2] {
        let (p0, p1) = geopotentials(self.gas, self.sigma, ps, theta, phi_s);
        [p0, p1]
    }

    /// The vertical-structure matrix of linear gravity waves about a resting,
    /// flat state with surface pressure `ps` and layer potential temperatures
    /// `theta`: `∂²δ/∂t² = −∇²(A δ)` for the layer divergences `δ = (δ₀, δ₁)`.
    ///
    /// `A = J B`, where `B` gives the tendencies of `(p_s', θ₀', θ₁')` from the
    /// divergences (linearized continuity, thermodynamics and interface flux)
    /// and `J` the change of each layer's pressure-gradient potential
    /// `Mₖ = Φₖ + c_p θ̄ₖ Pₖ`. Its eigenvalues are `c²` of the external and
    /// internal modes ([`TwoLayer::gravity_wave_speeds`]); the semi-implicit
    /// scheme treats exactly these terms implicitly.
    pub fn gravity_wave_matrix(&self, ps: f64, theta: [f64; 2]) -> [[f64; 2]; 2] {
        let Linearization { j, b } = self.linearization(ps, theta);
        let mut a = [[0.0; 2]; 2];
        for (r, row) in a.iter_mut().enumerate() {
            for (c, out) in row.iter_mut().enumerate() {
                *out = (0..3).map(|i| j[r][i] * b[i][c]).sum();
            }
        }
        a
    }

    /// The linear gravity-wave terms about a resting, flat state with surface
    /// pressure `ps` and layer potential temperatures `theta` (see
    /// [`Linearization`]).
    pub fn linearization(&self, ps: f64, theta: [f64; 2]) -> Linearization {
        let s = self.sigma;
        let th_hat = 0.5 * (theta[0] + theta[1]);
        let b = [
            [-ps * s, -ps * (1.0 - s)],
            [
                -(1.0 - s) * (theta[0] - th_hat),
                (1.0 - s) * (theta[0] - th_hat),
            ],
            [-s * (th_hat - theta[1]), s * (th_hat - theta[1])],
        ];
        // J by central differences of Mₖ = Φₖ + c_p θ̄ₖ Pₖ (θ̄ frozen in the
        // second term, as in the pressure-gradient force).
        let cp = self.gas.heat_capacity;
        let potential = |ps: f64, th: [f64; 2]| -> [f64; 2] {
            let phi = self.geopotentials_at(ps, th, 0.0);
            let exner = self.layer_exner(ps);
            [
                phi[0] + cp * theta[0] * exner[0],
                phi[1] + cp * theta[1] * exner[1],
            ]
        };
        let (ep, et) = (1.0e-5 * ps.abs().max(1.0), 1.0e-3);
        let mut j = [[0.0; 3]; 2];
        let (hi, lo) = (potential(ps + ep, theta), potential(ps - ep, theta));
        for k in 0..2 {
            j[k][0] = (hi[k] - lo[k]) / (2.0 * ep);
        }
        for l in 0..2 {
            let (mut up, mut down) = (theta, theta);
            up[l] += et;
            down[l] -= et;
            let (hi, lo) = (potential(ps, up), potential(ps, down));
            for k in 0..2 {
                j[k][1 + l] = (hi[k] - lo[k]) / (2.0 * et);
            }
        }
        Linearization { j, b }
    }

    /// Phase speeds `[external, internal]` (m s⁻¹) of linear gravity waves
    /// about a resting, flat state: `c² = −eigenvalues` of
    /// [`TwoLayer::gravity_wave_matrix`].
    pub fn gravity_wave_speeds(&self, ps: f64, theta: [f64; 2]) -> [f64; 2] {
        let a = self.gravity_wave_matrix(ps, theta);
        let (tr, det) = (a[0][0] + a[1][1], a[0][0] * a[1][1] - a[0][1] * a[1][0]);
        let disc = (0.25 * tr * tr - det).max(0.0).sqrt();
        let (l1, l2) = (0.5 * tr - disc, 0.5 * tr + disc);
        [(-l1).max(0.0).sqrt(), (-l2).max(0.0).sqrt()]
    }

    /// Geopotential of each layer (m² s⁻²) at every cell.
    pub fn geopotential(&self, state: &TwoLayerState) -> [Vec<f64>; 2] {
        let theta = [state.theta(0), state.theta(1)];
        let mut phi = [vec![0.0; state.ps.len()], vec![0.0; state.ps.len()]];
        for i in 0..state.ps.len() {
            let (p0, p1) = geopotentials(
                self.gas,
                self.sigma,
                state.ps[i],
                [theta[0][i], theta[1][i]],
                self.surface_geopotential[i],
            );
            phi[0][i] = p0;
            phi[1][i] = p1;
        }
        phi
    }

    /// Adiabatic, inviscid time tendencies of `state` into `out` (no
    /// diffusion: see [`diffusion_tendency`](Self::diffusion_tendency)).
    pub fn tendency(&mut self, state: &TwoLayerState, out: &mut TwoLayerState) {
        let params = Parameters {
            ops: &self.ops,
            gas: self.gas,
            sigma: self.sigma,
            gravity: self.gravity,
            coriolis: &self.coriolis,
            surface_geopotential: &self.surface_geopotential,
            humidity_profile: &self.humidity_profile,
        };
        compute_tendency(params, &mut self.work, state, out);
    }

    /// Hyperdiffusion tendencies of `state` into `out` (zero without
    /// hyperdiffusion; `out.ps` is always zero).
    pub fn diffusion_tendency(&mut self, state: &TwoLayerState, out: &mut TwoLayerState) {
        for x in out.parts_mut() {
            x.fill(0.0);
        }
        if let Some(d) = &mut self.diffusion {
            d.add_tendency(&self.ops, state, out);
        }
    }

    /// Forcing tendencies of `state` into `out` (zero without forcing;
    /// `out.ps` is always zero).
    pub fn forcing_tendency(&mut self, state: &TwoLayerState, out: &mut TwoLayerState) {
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
    pub fn energy_tendency(&self, state: &TwoLayerState, tendency: &TwoLayerState) -> EnergyTerms {
        let ops = &self.ops;
        let r2 = ops.radius() * ops.radius();
        let g = ops.geometry();
        let ds = self.dsigma();
        let edges = ops.mesh().edges();
        let edge_areas = ops.edge_areas();
        let mut rate = EnergyTerms::default();
        let gravity = self.gravity;
        let cp = self.gas.heat_capacity;
        // Layer Exner functions, once per cell.
        let exner: Vec<[f64; 2]> = state
            .ps
            .par_iter()
            .with_min_len(MIN_CHUNK)
            .map(|&ps| self.layer_exner(ps))
            .collect();
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
                .map(|i| g.cell_areas[i] * r2 * dsk * cp * exner[i][k] / gravity)
                .collect();
            rate.thermal += dot(&enthalpy, &tendency.theta_ps[k]);
        }
        rate
    }

    /// Advance `state` by `dt` seconds with classic fourth-order Runge–Kutta,
    /// including hyperdiffusion and forcing (and [`Forcing::advance`] at the
    /// end).
    pub fn step(&mut self, state: &mut TwoLayerState, dt: f64) {
        if self.stage_state.is_moist() != state.is_moist() {
            self.stages = std::array::from_fn(|_| TwoLayerState::zeros_like(state));
            self.stage_state = TwoLayerState::zeros_like(state);
        }
        let Self {
            ops,
            gas,
            sigma,
            gravity,
            coriolis,
            surface_geopotential,
            work,
            stages,
            stage_state,
            diffusion,
            forcing,
            humidity_profile,
        } = self;
        let params = Parameters {
            ops,
            gas: *gas,
            sigma: *sigma,
            gravity: *gravity,
            coriolis,
            surface_geopotential,
            humidity_profile,
        };
        let mut full = |x: &TwoLayerState, out: &mut TwoLayerState| {
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
    pub fn vorticity(&self, state: &TwoLayerState, k: usize) -> Vec<f64> {
        let mut zeta = vec![0.0; self.ops.mesh().n_vertices()];
        self.ops.curl(&state.u[k], &mut zeta);
        zeta
    }

    /// Mass, θ, kinetic and total energy of `state`.
    pub fn diagnostics(&self, state: &TwoLayerState) -> TwoLayerDiagnostics {
        let ops = &self.ops;
        let (m, g) = (ops.mesh(), ops.geometry());
        let r2 = ops.radius() * ops.radius();
        let ds = self.dsigma();
        let theta = [state.theta(0), state.theta(1)];
        let mut kinetic = [vec![0.0; m.n_cells()], vec![0.0; m.n_cells()]];
        for (u, ke) in state.u.iter().zip(&mut kinetic) {
            ops.kinetic_energy(u, ke);
        }
        let (mut mass, mut th, mut ke, mut energy) = (0.0, 0.0, 0.0, 0.0);
        let water = self.water(state);
        for i in 0..m.n_cells() {
            let a = g.cell_areas[i] * r2 / self.gravity;
            let ps = state.ps[i];
            let exner = self.layer_exner(ps);
            mass += a * ps;
            energy += a * ps * self.surface_geopotential[i];
            for k in 0..2 {
                let dp = ds[k] * ps;
                th += a * dp * theta[k][i];
                ke += a * dp * kinetic[k][i];
                energy +=
                    a * dp * (kinetic[k][i] + self.gas.heat_capacity * theta[k][i] * exner[k]);
            }
        }
        TwoLayerDiagnostics {
            mass,
            theta: th,
            kinetic_energy: ke,
            energy,
            water,
        }
    }

    /// Total water vapour of `state` (kg; 0 when dry).
    pub fn water(&self, state: &TwoLayerState) -> f64 {
        self.water_tendency(state)
    }

    /// Global integral `Σₖ Σ Aᵢ δσₖ xₖ / g` of a humidity field or tendency
    /// (`x.q_ps`): water in kg, or its rate of change in kg s⁻¹.
    pub fn water_tendency(&self, x: &TwoLayerState) -> f64 {
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
    pub fn stable_dt(&self, state: &TwoLayerState, courant: f64) -> f64 {
        let g = self.ops.geometry();
        let min_d = g.dual_edge_lengths.iter().fold(f64::MAX, |a, &b| a.min(b)) * self.ops.radius();
        let max_u = state
            .u
            .iter()
            .flatten()
            .fold(0.0f64, |a, &u| a.max(u.abs()));
        let mut max_t = 0.0f64;
        for k in 0..2 {
            for (tp, &ps) in state.theta_ps[k].iter().zip(&state.ps) {
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
    sigma: f64,
    gravity: f64,
    coriolis: &'p [f64],
    surface_geopotential: &'p [f64],
    humidity_profile: &'p RadiationParams,
}

fn interface_ratios(params: &Parameters, state: &TwoLayerState) -> [Vec<f64>; 2] {
    let ratios: Vec<[f64; 2]> = (0..state.ps.len())
        .into_par_iter()
        .with_min_len(MIN_CHUNK / 16)
        .map(|i| {
            let ps = state.ps[i];
            let exner = layer_exner(params.gas, params.sigma, ps);
            let t = [0, 1].map(|k| state.theta_ps[k][i] / ps * exner[k]);
            interface_humidity_ratios(
                params.humidity_profile,
                params.gas,
                ps,
                params.sigma,
                params.gravity,
                t,
            )
        })
        .collect();
    [
        ratios.iter().map(|r| r[0]).collect(),
        ratios.iter().map(|r| r[1]).collect(),
    ]
}

/// Layer Exner functions for surface pressure `ps`.
fn layer_exner(gas: Gas, sigma: f64, ps: f64) -> [f64; 2] {
    let kappa = gas.kappa();
    let pi = sigma * ps;
    let (exner_i, exner_s) = (gas.exner(pi), gas.exner(ps));
    [
        exner_i / (1.0 + kappa),
        (ps * exner_s - pi * exner_i) / ((1.0 + kappa) * (1.0 - sigma) * ps),
    ]
}

/// Layer geopotentials `(Φ₀, Φ₁)` at one cell (see the module docs).
fn geopotentials(gas: Gas, sigma: f64, ps: f64, theta: [f64; 2], phi_s: f64) -> (f64, f64) {
    let [p0, p1] = layer_exner(gas, sigma, ps);
    let cp = gas.heat_capacity;
    let theta_hat = 0.5 * (theta[0] + theta[1]);
    let x = cp * gas.kappa() * (sigma * theta[0] * p0 + (1.0 - sigma) * theta[1] * p1);
    let jump = cp * theta_hat * (p1 - p0);
    let phi1 = phi_s + x - sigma * jump;
    (phi1 + jump, phi1)
}

fn compute_tendency(
    params: Parameters,
    work: &mut Work,
    state: &TwoLayerState,
    out: &mut TwoLayerState,
) {
    let Parameters {
        ops,
        gas,
        sigma,
        coriolis,
        surface_geopotential,
        ..
    } = params;
    let ds = [sigma, 1.0 - sigma];
    let cp = gas.heat_capacity;
    let Work {
        layers,
        w,
        w_edge,
        theta_interface,
        cell_scratch,
        edge_scratch,
    } = work;

    // Thickness fluxes, divergence, θ, Exner and geopotential per layer.
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
    {
        let [l0, l1] = layers;
        // Exner functions and geopotentials need both layers' θ.
        l0.exner
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(l1.exner.par_iter_mut())
            .zip(l0.geopotential.par_iter_mut())
            .zip(l1.geopotential.par_iter_mut())
            .enumerate()
            .for_each(|(i, (((e0, e1), g0), g1))| {
                let ps = state.ps[i];
                [*e0, *e1] = layer_exner(gas, sigma, ps);
                (*g0, *g1) = geopotentials(
                    gas,
                    sigma,
                    ps,
                    [l0.theta[i], l1.theta[i]],
                    surface_geopotential[i],
                );
            });
        // Interface mass flux W = σᵢ D₁ − (1 − σᵢ) D₀, and θ̂.
        w.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, w)| *w = ds[0] * l1.divergence[i] - ds[1] * l0.divergence[i]);
        theta_interface
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, t)| *t = 0.5 * (l0.theta[i] + l1.theta[i]));
    }
    ops.cell_to_edge(w, w_edge);

    // Surface pressure.
    {
        let [l0, l1] = &*layers;
        out.ps
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, d)| *d = -(l0.divergence[i] + l1.divergence[i]));
    }

    for k in 0..2 {
        let l = &mut layers[k];
        // Thermodynamics: δσₖ ∂Θₖ/∂t = −∇·(Fₖ θₖ,ₑ) ∓ W θ̂.
        ops.cell_to_edge(&l.theta, &mut l.theta_edge);
        edge_scratch
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(e, f)| *f = l.flux[e] * l.theta_edge[e]);
        ops.divergence(edge_scratch, &mut l.theta_flux_divergence);
        let sign = if k == 0 { -1.0 } else { 1.0 };
        out.theta_ps[k]
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, d)| {
                *d = (-l.theta_flux_divergence[i] + sign * w[i] * theta_interface[i]) / ds[k];
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
    // carry it): δσₖ ∂Qₖ/∂t = −∇·(Fₖ qₖ,ₑ) ∓ W q̂.
    if state.is_moist() && out.is_moist() {
        let ratio = interface_ratios(&params, state);
        let q: [Vec<f64>; 2] = std::array::from_fn(|k| {
            state.q_ps[k]
                .par_iter()
                .with_min_len(MIN_CHUNK)
                .zip(state.ps.par_iter())
                .map(|(q, p)| q / p)
                .collect()
        });
        for k in 0..2 {
            let l = &mut layers[k];
            ops.cell_to_edge(&q[k], &mut l.theta_edge);
            edge_scratch
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(e, f)| *f = l.flux[e] * l.theta_edge[e]);
            ops.divergence(edge_scratch, cell_scratch);
            let sign = if k == 0 { -1.0 } else { 1.0 };
            let (q, ratio) = (&q, &ratio);
            out.q_ps[k]
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, d)| {
                    let q_hat = 0.5 * (ratio[0][i] * q[0][i] + ratio[1][i] * q[1][i]);
                    *d = (-cell_scratch[i] + sign * w[i] * q_hat) / ds[k];
                });
            // Restore θ at edges for the momentum equation below.
            ops.cell_to_edge(&l.theta, &mut l.theta_edge);
        }
    }

    // Assemble the wind tendencies, including vertical advection.
    let [l0, l1] = &*layers;
    let (u0, u1) = (&state.u[0], &state.u[1]);
    let [out0, out1] = &mut out.u;
    out0.par_iter_mut()
        .with_min_len(MIN_CHUNK)
        .zip(out1.par_iter_mut())
        .enumerate()
        .for_each(|(e, (d0, d1))| {
            let shear = u1[e] - u0[e];
            *d0 = l0.pv_flux[e]
                - l0.gradient[e]
                - cp * l0.theta_edge[e] * l0.exner_gradient[e]
                - w_edge[e] * shear / (2.0 * l0.dp_edge[e]);
            *d1 = l1.pv_flux[e]
                - l1.gradient[e]
                - cp * l1.theta_edge[e] * l1.exner_gradient[e]
                - w_edge[e] * shear / (2.0 * l1.dp_edge[e]);
        });
}
