//! Held & Suarez (1994) forcing: the standard benchmark for dry dynamical
//! cores. Temperature relaxes toward a zonally symmetric radiative-equilibrium
//! profile, and Rayleigh friction damps winds in the boundary layer. The
//! resulting climate (jets, surface easterlies and westerlies, baroclinic
//! eddies) is compared between models.
//!
//! # Forcing (per unit mass, at σ = p/p_s)
//!
//! ```text
//! ∂T/∂t = −k_T (T − T_eq)            ∂u/∂t = −k_v u
//! T_eq  = max(200 K, [315 K − ΔT_y sin²φ − Δθ_z ln(p/p₀) cos²φ] (p/p₀)^κ)
//! k_T   = k_a + (k_s − k_a) f(σ) cos⁴φ,   k_v = k_f f(σ),   f = max(0, (σ − σ_b)/(1 − σ_b))
//! ```
//!
//! with `k_f = 1/day`, `k_a = 1/40 day`, `k_s = 1/4 day`, `σ_b = 0.7`,
//! `ΔT_y = 60 K`, `Δθ_z = 10 K`, `p₀ = 1000 hPa`.
//!
//! # Layers
//!
//! The model's layer values are mass-weighted means over σ, so the forcing is
//! averaged the same way. Relaxing `T` at fixed `p` is relaxing θ toward
//! `θ_eq = T_eq (p₀/p)^κ`. For a layer with uniform θₖ:
//!
//! ```text
//! ∂θₖ/∂t = −⟨k_T⟩ θₖ + ⟨k_T θ_eq⟩,    ∂uₖ/∂t = −⟨k_v⟩ uₖ
//! ```
//!
//! where `⟨·⟩` is the mean over the layer's σ range. Friction acts only in
//! layers reaching below `σ_b` (for two equal-mass layers, the lower at
//! `⟨k_v⟩ = 0.3/day`); boundary-layer cooling likewise. The averages are computed by Gauss–Legendre quadrature
//! in `x = σ^(1−κ)`, which makes the stratospheric branch `∝ σ^(−κ)` exact
//! despite its singularity at the model top, with panel breaks at `σ_b`.
//! Per node the dependence on `p_s` separates (`ln p = ln σ + ln p_s`), so
//! evaluating the forcing costs no transcendental functions per node.

use fms_grid::operators::MIN_CHUNK;
use rayon::prelude::*;

use crate::forcing::Forcing;
use crate::layered::{LayeredModel, LayeredState, REFERENCE_PRESSURE};

const DAY: f64 = 86_400.0;

/// Held–Suarez parameters. [`Default`] is the published set.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeldSuarezParams {
    /// Boundary-layer friction rate `k_f` (s⁻¹).
    pub friction: f64,
    /// Free-atmosphere relaxation rate `k_a` (s⁻¹).
    pub relaxation_free: f64,
    /// Surface relaxation rate `k_s` (s⁻¹) at the equator.
    pub relaxation_surface: f64,
    /// Top of the boundary layer, `σ_b`.
    pub sigma_boundary: f64,
    /// Equator-to-pole temperature difference `ΔT_y` (K).
    pub delta_t_y: f64,
    /// Static-stability parameter `Δθ_z` (K).
    pub delta_theta_z: f64,
    /// Equatorial surface temperature (K).
    pub t_surface: f64,
    /// Stratospheric temperature floor (K).
    pub t_stratosphere: f64,
    /// Reference pressure `p₀` (Pa).
    pub p0: f64,
}

impl Default for HeldSuarezParams {
    fn default() -> Self {
        Self {
            friction: 1.0 / DAY,
            relaxation_free: 1.0 / (40.0 * DAY),
            relaxation_surface: 1.0 / (4.0 * DAY),
            sigma_boundary: 0.7,
            delta_t_y: 60.0,
            delta_theta_z: 10.0,
            t_surface: 315.0,
            t_stratosphere: 200.0,
            p0: REFERENCE_PRESSURE,
        }
    }
}

/// A quadrature node within a layer.
#[derive(Clone, Copy, Debug)]
struct Node {
    /// Weight: the node's share of the layer mean (they sum to 1).
    weight: f64,
    /// `ln σ`.
    ln_sigma: f64,
    /// `σ^(−κ)`.
    sigma_pow: f64,
    /// Boundary-layer factor `f(σ)`.
    boundary: f64,
}

/// Panels per stretch of σ between breaks. The kink where the stratospheric
/// floor takes over falls inside a panel; 16 panels put the error there
/// below 0.02 K.
const PANELS: usize = 16;
/// 4-point Gauss–Legendre on [−1, 1].
const GAUSS: [(f64, f64); 4] = [
    (-0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
    (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
];

/// Held–Suarez forcing for one model.
#[derive(Clone, Debug)]
pub struct HeldSuarez {
    params: HeldSuarezParams,
    kappa: f64,
    /// `sin²φ` per cell.
    sin2: Vec<f64>,
    nodes: Vec<Vec<Node>>,
    /// Layer-mean friction rate `⟨k_v⟩` (s⁻¹).
    friction: Vec<f64>,
}

impl HeldSuarez {
    /// Forcing for `model` (its mesh, gas and levels).
    pub fn new(model: &LayeredModel, params: HeldSuarezParams) -> Self {
        let kappa = model.gas().kappa();
        let levels = model.levels();
        let sin2 = model
            .operators()
            .mesh()
            .cell_centers()
            .iter()
            .map(|p| p.z * p.z)
            .collect();
        let nodes: Vec<Vec<Node>> = (0..levels.count())
            .map(|k| {
                let (a, b) = levels.bounds(k);
                layer_nodes(a, b, params.sigma_boundary, kappa)
            })
            .collect();
        let friction = nodes
            .iter()
            .map(|n| params.friction * n.iter().map(|n| n.weight * n.boundary).sum::<f64>())
            .collect();
        Self {
            params,
            kappa,
            sin2,
            nodes,
            friction,
        }
    }

    /// The parameters.
    pub fn params(&self) -> HeldSuarezParams {
        self.params
    }

    /// Layer-mean friction rates, top first (s⁻¹).
    pub fn friction(&self) -> &[f64] {
        &self.friction
    }

    /// Radiative-equilibrium potential temperature (K) at latitude with
    /// `sin²φ = sin2` and pressure `p` (Pa).
    pub fn theta_eq(&self, sin2: f64, p: f64) -> f64 {
        let hs = &self.params;
        let ratio = p / hs.p0;
        let troposphere = hs.t_surface
            - hs.delta_t_y * sin2
            - hs.delta_theta_z * fms_math::ln(ratio) * (1.0 - sin2);
        let stratosphere = hs.t_stratosphere * fms_math::pow(ratio, -self.kappa);
        troposphere.max(stratosphere)
    }

    /// For layer `k`, the relaxation rate `⟨k_T⟩` (s⁻¹) and `⟨k_T θ_eq⟩`
    /// (K s⁻¹) for a column at `sin²φ = sin2` with surface pressure `ps`.
    fn relaxation(&self, k: usize, sin2: f64, ln_ps: f64, ps_pow: f64) -> (f64, f64) {
        let hs = &self.params;
        let cos2 = 1.0 - sin2;
        let cos4 = cos2 * cos2;
        let troposphere = hs.t_surface - hs.delta_t_y * sin2;
        let strat = hs.t_stratosphere * ps_pow;
        let (mut rate, mut forcing) = (0.0, 0.0);
        for n in &self.nodes[k] {
            let theta_eq = (troposphere - hs.delta_theta_z * (n.ln_sigma + ln_ps) * cos2)
                .max(strat * n.sigma_pow);
            let k_t = hs.relaxation_free
                + (hs.relaxation_surface - hs.relaxation_free) * n.boundary * cos4;
            rate += n.weight * k_t;
            forcing += n.weight * k_t * theta_eq;
        }
        (rate, forcing)
    }

    /// Layer-mean relaxation rates `⟨k_T⟩` (s⁻¹) and equilibrium potential
    /// temperatures `⟨k_T θ_eq⟩ / ⟨k_T⟩` (K) at a cell for surface pressure
    /// `ps` (Pa), top layer first.
    pub fn equilibrium(&self, cell: usize, ps: f64) -> Vec<(f64, f64)> {
        let (ln_ps, ps_pow) = self.surface_terms(ps);
        (0..self.nodes.len())
            .map(|k| {
                let (rate, forcing) = self.relaxation(k, self.sin2[cell], ln_ps, ps_pow);
                (rate, forcing / rate)
            })
            .collect()
    }

    /// The standard starting point: at rest, surface pressure `p₀`
    /// everywhere, each layer's θ the global mean of its equilibrium θ, plus
    /// deterministic noise of amplitude `noise` (K) to break the symmetry.
    /// For a flat planet.
    pub fn initial_state(&self, model: &LayeredModel, noise: f64) -> LayeredState {
        let ops = model.operators();
        let areas = &ops.geometry().cell_areas;
        let (nc, ne) = (ops.mesh().n_cells(), ops.mesh().n_edges());
        let p0 = self.params.p0;
        let total: f64 = areas.iter().sum();
        let mean = |k: usize| -> f64 {
            (0..nc)
                .map(|i| areas[i] * self.equilibrium(i, p0)[k].1)
                .sum::<f64>()
                / total
        };
        let n = model.layers();
        let means: Vec<f64> = (0..n).map(mean).collect();
        let mut x = 0x2545_F491_4F6C_DD1Du64;
        let theta = means
            .iter()
            .map(|&m| {
                (0..nc)
                    .map(|_| {
                        x ^= x << 13;
                        x ^= x >> 7;
                        x ^= x << 17;
                        m + noise * ((x >> 11) as f64 / (1u64 << 52) as f64 - 1.0)
                    })
                    .collect()
            })
            .collect();
        LayeredState::new(vec![p0; nc], vec![vec![0.0; ne]; n], theta)
    }

    /// `ln(p_s/p₀)` and `(p₀/p_s)^κ`.
    fn surface_terms(&self, ps: f64) -> (f64, f64) {
        let ln_ps = fms_math::ln(ps / self.params.p0);
        (ln_ps, fms_math::exp(-self.kappa * ln_ps))
    }
}

impl Forcing for HeldSuarez {
    fn add_tendency(&mut self, state: &LayeredState, out: &mut LayeredState) {
        let this = &*self;
        for (k, out) in out.theta_ps.iter_mut().enumerate() {
            out.par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, d)| {
                    let ps = state.ps[i];
                    let (ln_ps, ps_pow) = this.surface_terms(ps);
                    let (r, f) = this.relaxation(k, this.sin2[i], ln_ps, ps_pow);
                    // ∂Θ/∂t = p_s ∂θ/∂t at fixed p_s.
                    *d += ps * f - r * state.theta_ps[k][i];
                });
        }
        for k in 0..self.friction.len() {
            let rate = self.friction[k];
            if rate > 0.0 {
                out.u[k]
                    .par_iter_mut()
                    .with_min_len(MIN_CHUNK)
                    .zip(state.u[k].par_iter())
                    .for_each(|(d, &u)| *d -= rate * u);
            }
        }
    }
}

/// Quadrature nodes for the mean over `σ ∈ [a, b]`.
fn layer_nodes(a: f64, b: f64, sigma_b: f64, kappa: f64) -> Vec<Node> {
    let e = 1.0 - kappa;
    let to_x = |s: f64| fms_math::pow(s, e);
    let mut breaks = vec![a];
    if sigma_b > a && sigma_b < b {
        breaks.push(sigma_b);
    }
    breaks.push(b);
    let mut nodes = Vec::new();
    for w in breaks.windows(2) {
        let (x0, x1) = (to_x(w[0]), to_x(w[1]));
        let h = (x1 - x0) / PANELS as f64;
        for p in 0..PANELS {
            let mid = x0 + (p as f64 + 0.5) * h;
            for (t, g) in GAUSS {
                let x = mid + 0.5 * h * t;
                let sigma = fms_math::pow(x, 1.0 / e);
                // dσ/dx = σ^κ / (1 − κ).
                let weight = g * 0.5 * h * fms_math::pow(sigma, kappa) / e;
                nodes.push(Node {
                    weight,
                    ln_sigma: fms_math::ln(sigma),
                    sigma_pow: fms_math::pow(sigma, -kappa),
                    boundary: ((sigma - sigma_b) / (1.0 - sigma_b)).max(0.0),
                });
            }
        }
    }
    let total: f64 = nodes.iter().map(|n| n.weight).sum();
    for n in &mut nodes {
        n.weight /= total;
    }
    nodes
}
