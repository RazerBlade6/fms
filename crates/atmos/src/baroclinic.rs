//! Balanced zonal flows with vertical shear for a two-layer model, and the
//! linear theory of their baroclinic instability (Phillips 1954; the
//! two-level model in Holton, *An Introduction to Dynamic Meteorology*, §8.2).
//!
//! The flow is solid-body rotation in each layer, `uₖ = Uₖ cos φ`. Steady
//! balance (gradient wind plus thermal wind) needs, for each layer,
//!
//! ```text
//! ∇Φₖ + c_p θₖ ∇Pₖ = −∇Gₖ,   Gₖ = (aΩUₖ + Uₖ²/2) sin²φ
//! ```
//!
//! (the same balance as Williamson case 2, with the two-layer pressure-gradient
//! force). With the static stability `θ₀ − θ₁` held fixed, this is a 2×2
//! linear ODE in latitude for the interface potential temperature `θ̂` and the
//! surface pressure, integrated here from the equator to the poles with RK4
//! and interpolated onto the cells.

use std::f64::consts::FRAC_PI_2;

use fms_geo::coords::to_lat_lon;

use crate::layered::{LayeredModel, LayeredState};

/// Latitude steps for integrating the balance equations from 0 to 90°.
const TABLE_STEPS: usize = 4000;

/// A balanced, vertically sheared zonal flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShearedFlow {
    /// Equatorial speed of each layer's solid-body flow, `[upper, lower]`
    /// (m s⁻¹, eastward).
    pub speed: [f64; 2],
    /// Static stability: `θ₀ − θ₁`, upper minus lower layer (K).
    pub stability: f64,
    /// Interface potential temperature `θ̂` at the equator (K).
    pub theta_equator: f64,
    /// Surface pressure at the equator (Pa).
    pub ps_equator: f64,
    /// Planetary rotation rate (s⁻¹).
    pub omega: f64,
}

/// `(θ̂, p_s)` against latitude, on a uniform table from 0 to 90°.
struct Profile {
    values: Vec<(f64, f64)>,
}

impl Profile {
    fn at(&self, lat: f64) -> (f64, f64) {
        let x = lat.abs() / FRAC_PI_2 * TABLE_STEPS as f64;
        let i = (x.floor() as usize).min(TABLE_STEPS - 1);
        let f = x - i as f64;
        let (a, b) = (self.values[i], self.values[i + 1]);
        (a.0 + f * (b.0 - a.0), a.1 + f * (b.1 - a.1))
    }
}

impl ShearedFlow {
    /// The balanced state on `model`'s mesh (flat terrain is assumed).
    ///
    /// # Panics
    /// If the model doesn't have two layers.
    pub fn state(&self, model: &LayeredModel) -> LayeredState {
        let ops = model.operators();
        let (m, g) = (ops.mesh(), ops.geometry());
        let a = ops.radius();
        let profile = self.profile(model);
        let mut ps = Vec::with_capacity(m.n_cells());
        let mut theta = vec![
            Vec::with_capacity(m.n_cells()),
            Vec::with_capacity(m.n_cells()),
        ];
        for &p in m.cell_centers() {
            let (th, p_s) = profile.at(to_lat_lon(p).0);
            ps.push(p_s);
            theta[0].push(th + 0.5 * self.stability);
            theta[1].push(th - 0.5 * self.stability);
        }
        // Winds from vertex stream functions ψ = −a U sin φ, so the initial
        // flow is exactly non-divergent.
        let u = self
            .speed
            .iter()
            .map(|&speed| {
                m.edges()
                    .iter()
                    .zip(&g.edge_lengths)
                    .map(|(e, &l)| {
                        let psi = |v: u32| -a * speed * m.vertices()[v as usize].z;
                        -(psi(e.vertices[1]) - psi(e.vertices[0])) / (l * a)
                    })
                    .collect()
            })
            .collect();
        LayeredState::new(ps, u, theta)
    }

    /// Integrate the balance equations from the equator to the pole.
    fn profile(&self, model: &LayeredModel) -> Profile {
        assert_eq!(model.layers(), 2, "the sheared flow is for two layers");
        let a = model.operators().radius();
        let cp = model.gas().heat_capacity;
        let s = self.stability;
        // d/dφ of Gₖ = (aΩUₖ + Uₖ²/2) sin²φ.
        let coeff = self.speed.map(|u| a * self.omega * u + 0.5 * u * u);
        let slope = |lat: f64, y: (f64, f64)| -> (f64, f64) {
            let (th, ps) = y;
            let (sin, cos) = fms_math::sin_cos(lat);
            let dg = coeff.map(|c| c * 2.0 * sin * cos);
            // Φₖ and Pₖ as functions of (θ̂, p_s), differentiated numerically.
            let phi =
                |th: f64, ps: f64| model.geopotentials_at(ps, &[th + 0.5 * s, th - 0.5 * s], 0.0);
            let (et, ep) = (1e-3, 1.0);
            let (pt_hi, pt_lo) = (phi(th + et, ps), phi(th - et, ps));
            let (pp_hi, pp_lo) = (phi(th, ps + ep), phi(th, ps - ep));
            let (ex_hi, ex_lo) = (model.layer_exner(ps + ep), model.layer_exner(ps - ep));
            let dphi_dth = [
                (pt_hi[0] - pt_lo[0]) / (2.0 * et),
                (pt_hi[1] - pt_lo[1]) / (2.0 * et),
            ];
            let theta = [th + 0.5 * s, th - 0.5 * s];
            let dforce_dps: [f64; 2] = std::array::from_fn(|k| {
                (pp_hi[k] - pp_lo[k]) / (2.0 * ep)
                    + cp * theta[k] * (ex_hi[k] - ex_lo[k]) / (2.0 * ep)
            });
            // [dΦ₀/dθ̂  dF₀/dp_s] [θ̂'] = −[G₀']
            // [dΦ₁/dθ̂  dF₁/dp_s] [p_s']    [G₁']
            let det = dphi_dth[0] * dforce_dps[1] - dforce_dps[0] * dphi_dth[1];
            let d_th = (-dg[0] * dforce_dps[1] + dforce_dps[0] * dg[1]) / det;
            let d_ps = (-dphi_dth[0] * dg[1] + dg[0] * dphi_dth[1]) / det;
            (d_th, d_ps)
        };
        let h = FRAC_PI_2 / TABLE_STEPS as f64;
        let mut y = (self.theta_equator, self.ps_equator);
        let mut values = Vec::with_capacity(TABLE_STEPS + 1);
        values.push(y);
        for i in 0..TABLE_STEPS {
            let lat = i as f64 * h;
            let add = |y: (f64, f64), k: (f64, f64), f: f64| (y.0 + f * k.0, y.1 + f * k.1);
            let k1 = slope(lat, y);
            let k2 = slope(lat + 0.5 * h, add(y, k1, 0.5 * h));
            let k3 = slope(lat + 0.5 * h, add(y, k2, 0.5 * h));
            let k4 = slope(lat + h, add(y, k3, h));
            y = (
                y.0 + h / 6.0 * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0),
                y.1 + h / 6.0 * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1),
            );
            values.push(y);
        }
        Profile { values }
    }

    /// Linear-theory growth rate (s⁻¹) of zonal wavenumber `m` at latitude
    /// `lat`, from the two-level quasi-geostrophic model on the local β-plane
    /// (Holton §8.2.2, with no meridional structure):
    ///
    /// ```text
    /// δ = β²λ⁴ / (k⁴(k² + 2λ²)²) − U_T² (2λ² − k²) / (k² + 2λ²)
    /// growth = k √(−δ) if δ < 0, else 0
    /// ```
    ///
    /// with `U_T = (U₀ − U₁) cos φ / 2`, `k = m / (a cos φ)`, and `2λ² = f²/c²`,
    /// `c` the model's internal gravity-wave speed
    /// ([`LayeredModel::gravity_wave_speeds`]). Returns `(growth, critical)`: the growth
    /// rate, and the smallest `speed[0] − speed[1]` (equatorial speeds, as in
    /// [`ShearedFlow::speed`]) for which any wavenumber is unstable at this
    /// latitude, from `U_T > β/(2λ²)`.
    pub fn phillips_growth_rate(&self, model: &LayeredModel, m: u32, lat: f64) -> (f64, f64) {
        let a = model.operators().radius();
        let (sin, cos) = fms_math::sin_cos(lat);
        let f = 2.0 * self.omega * sin;
        let beta = 2.0 * self.omega * cos / a;
        let k = f64::from(m) / (a * cos);
        let ut = 0.5 * (self.speed[0] - self.speed[1]) * cos;
        let profile = self.profile(model);
        let (th, ps) = profile.at(lat);
        // The deformation wavenumber from the model's own internal gravity-wave
        // speed: 2λ² = f² / c_int².
        let c_int = model
            .gravity_wave_speeds(ps, &[th + 0.5 * self.stability, th - 0.5 * self.stability])[1];
        let lambda2 = f * f / (2.0 * c_int * c_int);
        let k2 = k * k;
        let delta = beta * beta * lambda2 * lambda2 / (k2 * k2 * (k2 + 2.0 * lambda2).powi(2))
            - ut * ut * (2.0 * lambda2 - k2) / (k2 + 2.0 * lambda2);
        let growth = if delta < 0.0 {
            k * (-delta).sqrt()
        } else {
            0.0
        };
        (growth, beta / lambda2 / cos)
    }
}
