//! The surface layer: drag, heat and moisture exchange between the surface
//! and the bottom model layer, the 10 m wind, and entrainment at the top of
//! the boundary layer.
//!
//! The bottom layer is treated as a well-mixed boundary layer: its wind,
//! potential temperature and specific humidity hold from the top of a thin
//! surface layer up to its top. Below, the wind follows the neutral log
//! profile `V(z) ∝ ln(z/z₀)`, so the 10 m wind is the layer's wind times
//! `ln(10/z₀)/ln(z_b/z₀)`, with `z_b` the height of the layer's mean level.
//! Fluxes use bulk formulas at 10 m:
//!
//! - stress `τ = ρ C_D V V₁₀`, sensible heat `ρ c_p C_H V (θ_s − θ_b) Π_s` and
//!   evaporation `ρ C_E V (q*(T_s) − q_b)`, with `V` the 10 m wind plus a
//!   gust from the boundary layer's convective velocity `w*` (Beljaars 1995)
//!   and a calm-air minimum, both in quadrature;
//! - neutral coefficients from the roughness lengths,
//!   `C_D = k²/ln²(10/z₀)` and `C_H = k²/(ln(10/z₀) ln(10/z₀ₕ))`; over the
//!   ocean `z₀` follows Charnock's relation plus a smooth-flow term, and the
//!   heat and moisture lengths are ECMWF's (`0.40ν/u*`, `0.62ν/u*`);
//! - stability from the bulk Richardson number at 10 m: Louis (1979) in
//!   unstable air, the long-tailed `1/(1 + 10 Ri)` in stable air.
//!
//! The boundary layer exchanges air with the layer above at the entrainment
//! velocity of mixed-layer theory, `w_e = A w*³/(h Δb)` (Lilly 1968; Deardorff
//! 1979), with `Δb` the buoyancy jump at its top. Where surface air saturates
//! below the top and is still buoyant there, shallow cumulus also carry its
//! water and momentum up with a cloud-base mass flux `0.03 ρ w*` (Grant
//! 2001); they reach the inversion at neutral buoyancy, so they carry no
//! heat across it.

use crate::radiation::saturation_humidity;

/// Von Kármán's constant.
pub const KARMAN: f64 = 0.4;
/// Kinematic viscosity of air (m² s⁻¹).
const VISCOSITY: f64 = 1.5e-5;
/// Height of the surface wind (m).
pub const WIND_HEIGHT: f64 = 10.0;
/// Water's molar mass over dry air's, as `1/ε − 1` for virtual temperature.
const VIRTUAL: f64 = 0.608;

/// Surface-layer parameters. [`Default`] is Earth-like.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceLayerParams {
    /// Charnock's constant: ocean roughness `α u*²/g`.
    pub charnock: f64,
    /// Gust factor: `β w*` added in quadrature to the wind (Beljaars 1995).
    pub gust: f64,
    /// Calm-air wind (m s⁻¹) added in quadrature: mesoscale variability.
    pub calm: f64,
    /// Entrainment efficiency `A` of mixed-layer theory.
    pub entrainment: f64,
    /// Largest entrainment velocity (m s⁻¹): keeps the explicit exchange
    /// well below the time step's limit.
    pub max_entrainment: f64,
    /// Shallow cumulus mass flux at cloud base per unit `ρ w*` (Grant 2001).
    pub cumulus: f64,
    /// Factor on the drag coefficient, for sensitivity experiments (1).
    pub drag_factor: f64,
}

impl Default for SurfaceLayerParams {
    fn default() -> Self {
        Self {
            charnock: 0.018,
            gust: 1.2,
            calm: 1.0,
            entrainment: 0.2,
            max_entrainment: 0.05,
            cumulus: 0.03,
            drag_factor: 1.0,
        }
    }
}

/// What lies under the boundary layer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Surface {
    /// Open water: roughness from the wind (Charnock).
    Ocean,
    /// A fixed momentum roughness length (m); heat and moisture use a tenth.
    Rough(f64),
}

/// The air at the bottom model layer and the surface under it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceAir {
    /// Speed of the bottom layer's wind (m s⁻¹).
    pub wind: f64,
    /// Height of the bottom layer's mean level above the surface (m).
    pub height: f64,
    /// The bottom layer's potential temperature (K).
    pub theta: f64,
    /// The bottom layer's specific humidity.
    pub humidity: f64,
    /// Surface temperature (K).
    pub surface_temperature: f64,
    /// Surface pressure (Pa).
    pub ps: f64,
    /// Exner function at the surface, `(p_s/p₀)^κ`.
    pub surface_exner: f64,
    /// The boundary layer's convective velocity `w*` (m s⁻¹).
    pub convective_velocity: f64,
    /// Gravity (m s⁻²).
    pub gravity: f64,
}

/// Exchange between the surface and the bottom layer.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Exchange {
    /// 10 m wind over the bottom layer's wind.
    pub wind_ratio: f64,
    /// 10 m wind speed (m s⁻¹), without gusts.
    pub wind: f64,
    /// Effective speed for the fluxes, with gusts (m s⁻¹).
    pub speed: f64,
    /// Drag coefficient at 10 m, with stability.
    pub drag: f64,
    /// Heat exchange coefficient at 10 m, with stability.
    pub heat: f64,
    /// Moisture exchange coefficient at 10 m, with stability.
    pub moisture: f64,
    /// Momentum roughness length (m).
    pub roughness: f64,
    /// Bulk Richardson number at 10 m.
    pub richardson: f64,
    /// Friction velocity (m s⁻¹).
    pub friction_velocity: f64,
}

impl Exchange {
    /// Momentum exchange `ρ C_D V` (kg m⁻² s⁻¹) per unit of 10 m wind.
    pub fn momentum(&self, density: f64) -> f64 {
        density * self.drag * self.speed
    }
}

/// Neutral exchange coefficient between heights `z` over roughness lengths
/// `a` and `b`.
fn neutral(z: f64, a: f64, b: f64) -> f64 {
    KARMAN * KARMAN / (fms_math::ln(z / a) * fms_math::ln(z / b))
}

/// Louis (1979) stability factors `(f_m, f_h)` for bulk Richardson number
/// `ri` at height `z` over roughness `z0`, with neutral drag `cn`; in stable
/// air the long-tailed `1/(1 + 10 Ri)`.
pub fn stability(ri: f64, z: f64, z0: f64, cn: f64) -> (f64, f64) {
    if ri >= 0.0 {
        let f = 1.0 / (1.0 + 10.0 * ri);
        (f, f)
    } else {
        // b = c = 5.
        let c = 75.0 * cn * (-ri * z / z0).sqrt();
        (1.0 - 10.0 * ri / (1.0 + c), 1.0 - 15.0 * ri / (1.0 + c))
    }
}

/// Surface exchange for `air` over `surface`.
pub fn exchange(params: &SurfaceLayerParams, surface: Surface, air: &SurfaceAir) -> Exchange {
    let z = WIND_HEIGHT;
    let height = air.height.max(2.0 * z);
    let gust2 = (params.gust * air.convective_velocity).powi(2) + params.calm * params.calm;
    // Virtual potential temperatures of the air and of saturated air at the
    // surface.
    let theta_s = air.surface_temperature / air.surface_exner;
    let q_s = saturation_humidity(air.surface_temperature, air.ps);
    let tv_air = air.theta * (1.0 + VIRTUAL * air.humidity);
    let tv_surface = theta_s * (1.0 + VIRTUAL * q_s);
    let roughness = |u_star: f64| match surface {
        Surface::Ocean => {
            let u = u_star.max(1e-3);
            let z0 = params.charnock * u * u / air.gravity + 0.11 * VISCOSITY / u;
            (z0, 0.40 * VISCOSITY / u, 0.62 * VISCOSITY / u)
        }
        Surface::Rough(z0) => (z0, 0.1 * z0, 0.1 * z0),
    };
    // Roughness and drag depend on each other over water: a few fixed-point
    // passes from a moderate sea converge to well under a percent.
    let mut u_star = 0.25;
    let mut result = Exchange::default();
    for _ in 0..4 {
        let (z0, z0h, z0q) = roughness(u_star);
        let ratio = fms_math::ln(z / z0) / fms_math::ln(height / z0);
        let wind = ratio * air.wind;
        let speed = (wind * wind + gust2).sqrt();
        let ri =
            air.gravity * z * (tv_air - tv_surface) / (0.5 * (tv_air + tv_surface) * speed * speed);
        let cn = neutral(z, z0, z0);
        let (fm, fh) = stability(ri, z, z0, cn);
        let drag = params.drag_factor * cn * fm;
        u_star = drag.sqrt() * speed;
        result = Exchange {
            wind_ratio: ratio,
            wind,
            speed,
            drag,
            heat: neutral(z, z0, z0h) * fh,
            moisture: neutral(z, z0, z0q) * fh,
            roughness: z0,
            richardson: ri,
            friction_velocity: u_star,
        };
    }
    result
}

/// Convective velocity `w*` (m s⁻¹) of a boundary layer `depth` m deep with
/// surface sensible heat flux `sensible` (W m⁻²) and evaporation
/// `evaporation` (kg m⁻² s⁻¹) under air of potential temperature `theta`,
/// density `density` and heat capacity `cp`, plus the turbulence driven from
/// above by `cloud_top_cooling` (W m⁻²), the longwave cooling of clouds at
/// its top (Lock et al. 2000); zero in stable air.
#[allow(clippy::too_many_arguments)]
pub fn convective_velocity(
    gravity: f64,
    depth: f64,
    theta: f64,
    density: f64,
    cp: f64,
    sensible: f64,
    evaporation: f64,
    cloud_top_cooling: f64,
) -> f64 {
    let buoyancy = gravity / theta
        * ((sensible + cloud_top_cooling.max(0.0)) / (density * cp)
            + VIRTUAL * theta * evaporation / density);
    if buoyancy > 0.0 {
        fms_math::pow(buoyancy * depth, 1.0 / 3.0)
    } else {
        0.0
    }
}

/// Entrainment velocity (m s⁻¹) at the top of a boundary layer `depth` m
/// deep with convective velocity `w_star`, under a jump in virtual potential
/// temperature `jump` (K) above air of virtual potential temperature
/// `theta_v`.
pub fn entrainment_velocity(
    params: &SurfaceLayerParams,
    gravity: f64,
    depth: f64,
    w_star: f64,
    theta_v: f64,
    jump: f64,
) -> f64 {
    let db = gravity * jump.max(0.0) / theta_v;
    let rising = params.entrainment * w_star.powi(3);
    if rising <= 0.0 {
        0.0
    } else if db * depth <= rising / params.max_entrainment {
        params.max_entrainment
    } else {
        rising / (depth * db)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn air(wind: f64, theta: f64, ts: f64) -> SurfaceAir {
        let ps = 1.0e5;
        SurfaceAir {
            wind,
            height: 650.0,
            theta,
            humidity: 0.8 * saturation_humidity(theta, ps),
            surface_temperature: ts,
            ps,
            surface_exner: 1.0,
            convective_velocity: 0.0,
            gravity: 9.81,
        }
    }

    #[test]
    fn neutral_ocean_drag_matches_observations() {
        // Over the sea, C_D at 10 m rises with the wind: ~1.0e-3 at 5 m/s
        // to ~1.5e-3 at 15 m/s (Large & Pond 1981; Smith 1988), and C_H,
        // C_E ~1.1–1.2e-3. Neutral: air at the sea's θ, dry enough not to
        // matter.
        let p = SurfaceLayerParams {
            calm: 0.0,
            ..Default::default()
        };
        let neutral = |v10: f64| {
            // Find the layer wind that gives the 10 m wind.
            let mut a = air(v10 / 0.7, 290.0, 290.0);
            a.humidity = saturation_humidity(290.0, a.ps);
            let mut e = exchange(&p, Surface::Ocean, &a);
            a.wind = v10 / e.wind_ratio;
            e = exchange(&p, Surface::Ocean, &a);
            e
        };
        let (slow, fast) = (neutral(5.0), neutral(15.0));
        assert!((slow.wind - 5.0).abs() < 0.1 && (fast.wind - 15.0).abs() < 0.2);
        assert!(slow.drag > 0.9e-3 && slow.drag < 1.25e-3, "{slow:?}");
        assert!(fast.drag > 1.4e-3 && fast.drag < 1.9e-3, "{fast:?}");
        for e in [slow, fast] {
            assert!(e.richardson.abs() < 1e-3);
            assert!(e.heat > 0.9e-3 && e.heat < 1.3e-3, "{e:?}");
            assert!(e.moisture > e.heat);
            // The 10 m wind is ~70% of the boundary layer's over the sea.
            assert!(e.wind_ratio > 0.65 && e.wind_ratio < 0.8, "{e:?}");
        }
    }

    #[test]
    fn stability_damps_and_convection_enhances_the_exchange() {
        let p = SurfaceLayerParams::default();
        let neutral = exchange(&p, Surface::Ocean, &air(8.0, 290.0, 290.0));
        // Warm air over a cold sea: stable, weaker exchange.
        let stable = exchange(&p, Surface::Ocean, &air(8.0, 295.0, 285.0));
        // Cold air over a warm sea: unstable, stronger.
        let unstable = exchange(&p, Surface::Ocean, &air(8.0, 285.0, 295.0));
        assert!(stable.richardson > 0.0 && unstable.richardson < 0.0);
        assert!(stable.heat < 0.7 * neutral.heat, "{stable:?} {neutral:?}");
        assert!(
            unstable.heat > 1.1 * neutral.heat,
            "{unstable:?} {neutral:?}"
        );
        // Over rough land the drag is several times the sea's and the 10 m
        // wind a smaller share of the layer's.
        let land = exchange(&p, Surface::Rough(0.1), &air(8.0, 290.0, 290.0));
        assert!(
            land.drag > 3.0 * neutral.drag && land.wind_ratio < 0.6,
            "{land:?}"
        );
        // Gusts keep the exchange going in calm, convecting air.
        let mut calm = air(0.0, 290.0, 290.0);
        calm.convective_velocity = 1.5;
        let e = exchange(&p, Surface::Ocean, &calm);
        assert!((e.speed - (1.8f64.powi(2) + 1.0).sqrt()).abs() < 1e-9);
    }

    #[test]
    fn entrainment_follows_mixed_layer_theory() {
        let p = SurfaceLayerParams::default();
        // Subtropical trade-wind layer: w* ~1 m/s, 1 km deep, a 5 K jump:
        // w_e = 0.2 / (1000 · 9.81 · 5 / 300) ≈ 1.2 cm/s.
        let w = entrainment_velocity(&p, 9.81, 1000.0, 1.0, 300.0, 5.0);
        assert!((w - 0.2 / (1000.0 * 9.81 * 5.0 / 300.0)).abs() < 1e-12);
        // No convection, no entrainment; no jump, the cap.
        assert_eq!(entrainment_velocity(&p, 9.81, 1000.0, 0.0, 300.0, 5.0), 0.0);
        assert_eq!(
            entrainment_velocity(&p, 9.81, 1000.0, 1.0, 300.0, 0.0),
            p.max_entrainment
        );
        // w* from a typical ocean flux: 15 W/m² sensible, 120 W/m² latent.
        let w_star =
            convective_velocity(9.81, 1000.0, 300.0, 1.15, 1004.0, 15.0, 120.0 / 2.5e6, 0.0);
        assert!(w_star > 0.4 && w_star < 0.9, "{w_star}");
        assert_eq!(
            convective_velocity(9.81, 1000.0, 300.0, 1.15, 1004.0, -10.0, 0.0, 0.0),
            0.0
        );
        // A cloud top cooling by 40 W/m² stirs a layer the surface doesn't.
        let cloudy = convective_velocity(9.81, 1000.0, 300.0, 1.15, 1004.0, 0.0, 0.0, 40.0);
        assert!(cloudy > 0.7 && cloudy < 1.3, "{cloudy}");
    }
}
