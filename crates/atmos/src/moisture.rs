//! Moist physics of one two-layer column: convection, large-scale
//! condensation, re-evaporation of falling rain, rain or snow, and the
//! clouds the radiation sees. Pure functions of the column's state, so the
//! 3-D model ([`crate::ColumnPhysics`]) and tests share them.
//!
//! # Humidity within a layer
//!
//! Each layer carries its mean specific humidity `qₖ`. Within the layer,
//! relative humidity follows the Manabe–Wetherald shape (moistest at the
//! bottom) scaled to hold `qₖ`, and temperature follows the radiation's
//! in-layer lapse rate ([`crate::radiation::layer_saturation`]). The layer's
//! relative humidity is `qₖ / q̄*ₖ`, its water over its saturation capacity.
//!
//! # Convection: Betts–Miller relaxation
//!
//! A parcel of surface air (the lower layer extrapolated to the ground, with
//! the profile's surface humidity) is lifted: dry-adiabatically to its
//! condensation level, then along the pseudo-adiabat. Its temperature,
//! averaged over each layer, is the reference profile; in the upper layer it
//! is bounded below by the stratosphere's temperature, taken as the skin
//! temperature `(OLR/2σ)^{1/4}` of the column. If the upper layer is colder
//! than its reference, the column is unstable, and both layers relax
//! towards the reference temperature and a reference humidity over
//! [`MoistureParams::convection_time`] (Betts 1986; Frierson 2007). The
//! reference humidity has the in-layer profile's shape with
//! [`MoistureParams::convection_humidity`] at the surface, so the upper
//! layer, which includes the dry stratosphere, is referred to a far lower
//! mean relative humidity than the lower:
//!
//! - **Deep** (the relaxation removes water): it rains out, and the
//!   reference temperature shifts uniformly so the column's heating is
//!   exactly the latent heat released (enthalpy plus latent energy is
//!   conserved).
//! - **Shallow** (it would add water): no rain; humidity and temperature are
//!   only redistributed between the layers, each conserved.
//!
//! The moist adiabat's lapse rate falls as the climate warms, so the upper
//! layer warms more than the surface: the lapse-rate feedback.
//!
//! # Large-scale condensation
//!
//! A layer more humid than [`MoistureParams::condensation_humidity`]
//! condenses the excess over [`MoistureParams::condensation_time`], with its
//! latent heat warming the layer (allowing for the warming raising the
//! saturation, as an implicit step would). The upper layer's rain falls
//! through the lower layer, which re-evaporates a share proportional to its
//! dryness.
//!
//! # Rain or snow
//!
//! Precipitation reaching the ground is snow below freezing at the surface
//! air temperature (a linear ramp over a few kelvin). Snow releases its
//! latent heat of fusion in the lower layer; the surface pays it back when
//! the snow melts.
//!
//! # Clouds
//!
//! Two sheets: high cloud in the upper layer, low cloud in the lower.
//! Stratiform cover grows linearly with each layer's relative humidity
//! between two bounds, `c = (RH − RH₁)/(RH₂ − RH₁)` clamped to [0, 1]
//! (SPEEDY (Molteni 2003) uses the square, for the most humid level; layer
//! means vary less, and the square can't give Earth's ~0.5 low cloud at its
//! layer-mean humidity ~0.63 while still leaving a condensing layer short of
//! overcast). As in SPEEDY, deep convection adds anvil cloud to the upper
//! sheet, growing with the square root of its rain.

use crate::radiation::{LayerSaturation, RadiationParams, layer_saturation, saturation_humidity};
use crate::two_layer::Gas;

/// Latent heat of vaporization (J kg⁻¹).
pub const LATENT_HEAT: f64 = 2.5e6;
/// Latent heat of fusion (J kg⁻¹).
pub const FUSION_HEAT: f64 = 3.34e5;
/// Gas constant of water vapour (J kg⁻¹ K⁻¹).
pub const VAPOUR_GAS_CONSTANT: f64 = 461.5;
/// Seconds per day, for rain in mm/day.
const DAY: f64 = 86_400.0;

/// Cloud parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CloudParams {
    /// Level σ of the high and low cloud sheets.
    pub sigma: [f64; 2],
    /// Albedo of overcast high and low cloud.
    pub albedo: [f64; 2],
    /// Longwave emissivity of overcast high and low cloud.
    pub emissivity: [f64; 2],
    /// Layer relative humidities at which stratiform cloud starts, upper
    /// and lower layer.
    pub critical_humidity: [f64; 2],
    /// Layer relative humidities at which stratiform cover is complete.
    pub overcast_humidity: [f64; 2],
    /// Anvil cover per √(mm/day) of convective rain.
    pub convective: f64,
    /// Largest anvil cover.
    pub max_convective: f64,
}

impl Default for CloudParams {
    fn default() -> Self {
        Self {
            sigma: [0.25, 0.8],
            albedo: [0.4, 0.5],
            emissivity: [0.8, 1.0],
            critical_humidity: [0.15, 0.3],
            overcast_humidity: [0.7, 1.0],
            convective: 0.15,
            max_convective: 0.6,
        }
    }
}

/// Moist-physics parameters. [`Default`] is Earth-like.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MoistureParams {
    /// Layer relative humidity above which water condenses, upper and lower
    /// layer (below 1: a layer saturates in places before its mean does).
    pub condensation_humidity: [f64; 2],
    /// Large-scale condensation time (s).
    pub condensation_time: f64,
    /// Surface relative humidity of convection's reference profile (the
    /// Manabe–Wetherald shape above it).
    pub convection_humidity: f64,
    /// Convective relaxation time (s).
    pub convection_time: f64,
    /// Share of the upper layer's rain a bone-dry lower layer re-evaporates.
    pub reevaporation: f64,
    /// Surface air temperatures (K) below which all precipitation is snow,
    /// and above which all is rain.
    pub snow_temperature: (f64, f64),
    /// Clouds.
    pub clouds: CloudParams,
}

impl Default for MoistureParams {
    fn default() -> Self {
        Self {
            condensation_humidity: [0.45, 0.7],
            condensation_time: 4.0 * 3600.0,
            convection_humidity: 0.7,
            convection_time: 2.0 * 3600.0,
            reevaporation: 0.5,
            snow_temperature: (271.15, 275.15),
            clouds: CloudParams::default(),
        }
    }
}

/// One column's state, as the moist physics sees it.
#[derive(Clone, Copy, Debug)]
pub struct MoistColumn {
    /// Surface pressure (Pa).
    pub ps: f64,
    /// Layer temperatures, upper first (K).
    pub temperature: [f64; 2],
    /// Layer specific humidities (kg kg⁻¹).
    pub humidity: [f64; 2],
    /// Surface air temperature (K).
    pub surface_air: f64,
    /// Stratosphere temperature (K): the floor of the upper layer's
    /// convective reference.
    pub stratosphere: f64,
}

/// The moist physics' tendencies and diagnostics for one column.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MoistTendency {
    /// Temperature tendency per layer (K s⁻¹).
    pub temperature: [f64; 2],
    /// Specific-humidity tendency per layer (s⁻¹).
    pub humidity: [f64; 2],
    /// Rain and snow reaching the surface (kg m⁻² s⁻¹).
    pub rain: f64,
    /// See `rain`.
    pub snow: f64,
    /// The part of the precipitation from convection (kg m⁻² s⁻¹).
    pub convective: f64,
    /// Cloud cover of the high and low sheets.
    pub clouds: [f64; 2],
    /// Layer relative humidities.
    pub relative_humidity: [f64; 2],
}

/// Geometry and constants the column physics needs.
#[derive(Clone, Copy, Debug)]
pub struct Setting<'p> {
    /// Moisture parameters.
    pub params: &'p MoistureParams,
    /// Radiation parameters (for the in-layer profile).
    pub radiation: &'p RadiationParams,
    /// Gas constants.
    pub gas: Gas,
    /// Interface σ.
    pub sigma: f64,
    /// Gravity (m s⁻²).
    pub gravity: f64,
}

impl Setting<'_> {
    fn saturation(&self, ps: f64, t: [f64; 2]) -> [LayerSaturation; 2] {
        layer_saturation(self.radiation, self.gas, ps, self.sigma, self.gravity, t)
    }

    fn masses(&self, ps: f64) -> [f64; 2] {
        [
            self.sigma * ps / self.gravity,
            (1.0 - self.sigma) * ps / self.gravity,
        ]
    }
}

/// `dq*/dT` (K⁻¹) at temperature `t` with saturation humidity `q`
/// (Clausius–Clapeyron).
fn saturation_slope(t: f64, q: f64) -> f64 {
    LATENT_HEAT * q / (VAPOUR_GAS_CONSTANT * t * t)
}

/// Gauss–Legendre nodes on [−1, 1] (as the radiation's layer quadrature).
const GAUSS: [(f64, f64); 4] = [
    (-0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
    (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
];

/// Temperature of a parcel lifted from `(t, q)` at pressure `p` to pressure
/// `p_to < p`: dry adiabat (q conserved) until saturated, then the
/// pseudo-adiabat `dT/d ln p = (R T + L q*)/(c_p + L² q*/(R_v T²))`
/// (midpoint steps of at most 0.05 in ln p). Returns the temperature and
/// humidity at `p_to`.
pub fn lift(gas: Gas, t: f64, q: f64, p: f64, p_to: f64) -> (f64, f64) {
    let (r, cp) = (gas.gas_constant, gas.heat_capacity);
    let span = fms_math::ln(p / p_to);
    if span <= 0.0 {
        return (t, q);
    }
    let steps = (span / 0.05).ceil().max(1.0);
    let h = span / steps;
    let (mut t, mut q, mut lnp) = (t, q, fms_math::ln(p));
    let slope = |t: f64, lnp: f64, q: f64| -> (f64, bool) {
        let qs = saturation_humidity(t, fms_math::exp(lnp));
        if q < qs {
            // Dry: θ conserved, dT/d ln p = κ T.
            (r * t / cp, false)
        } else {
            (
                (r * t + LATENT_HEAT * qs)
                    / (cp + LATENT_HEAT * LATENT_HEAT * qs / (VAPOUR_GAS_CONSTANT * t * t)),
                true,
            )
        }
    };
    for _ in 0..steps as u64 {
        let (s1, _) = slope(t, lnp, q);
        let (s2, saturated) = slope(t - 0.5 * h * s1, lnp - 0.5 * h, q);
        t -= h * s2;
        lnp -= h;
        if saturated {
            q = q.min(saturation_humidity(t, fms_math::exp(lnp)));
        }
    }
    (t, q)
}

/// The convective reference: layer-mean temperatures of a surface-air
/// parcel lifted through the column (the upper layer's bounded below by
/// `column.stratosphere`).
pub fn convective_reference(setting: &Setting, column: &MoistColumn, surface_q: f64) -> [f64; 2] {
    let (ps, sigma) = (column.ps, setting.sigma);
    // Lift through the lower layer's nodes, then the upper's, bottom up.
    let mut nodes: Vec<(usize, f64, f64)> = Vec::with_capacity(8);
    for (k, (a, b)) in [(0.0, sigma), (sigma, 1.0)].into_iter().enumerate() {
        for (x, w) in GAUSS {
            nodes.push((k, 0.5 * (a + b) + 0.5 * (b - a) * x, 0.5 * w));
        }
    }
    nodes.sort_by(|a, b| b.1.total_cmp(&a.1));
    let (mut t, mut q, mut p) = (column.surface_air, surface_q, ps);
    let mut reference = [0.0; 2];
    for (k, s, w) in nodes {
        let p_to = s * ps;
        (t, q) = lift(setting.gas, t, q, p, p_to);
        p = p_to;
        let value = if k == 0 {
            t.max(column.stratosphere)
        } else {
            t
        };
        reference[k] += w * value;
    }
    reference
}

/// Stratiform cloud cover for layer relative humidity `rh`, from none at
/// `rh_c` to overcast at `rh_o`.
pub fn stratiform_cover(rh: f64, rh_c: f64, rh_o: f64) -> f64 {
    ((rh - rh_c) / (rh_o - rh_c)).clamp(0.0, 1.0)
}

/// The moist physics of one column.
pub fn moist_physics(setting: &Setting, column: &MoistColumn) -> MoistTendency {
    let p = setting.params;
    let (gas, ps) = (setting.gas, column.ps);
    let (cp, l) = (gas.heat_capacity, LATENT_HEAT);
    let m = setting.masses(ps);
    let t = column.temperature;
    let q = column.humidity.map(|q| q.max(0.0));
    let sat = setting.saturation(ps, t);
    let rh = [q[0] / sat[0].mean, q[1] / sat[1].mean];
    let mut out = MoistTendency {
        relative_humidity: rh,
        ..MoistTendency::default()
    };

    // Convection.
    let q_star_surface = saturation_humidity(column.surface_air, ps);
    let surface_q = (q[1] / sat[1].shaped * q_star_surface).min(q_star_surface);
    let reference = convective_reference(setting, column, surface_q);
    if reference[0] > t[0] {
        let tau = p.convection_time;
        let ref_sat = setting.saturation(ps, reference);
        let mut dt = [(reference[0] - t[0]) / tau, (reference[1] - t[1]) / tau];
        let mut dq = [
            (p.convection_humidity * ref_sat[0].shaped - q[0]) / tau,
            (p.convection_humidity * ref_sat[1].shaped - q[1]) / tau,
        ];
        let rain = -(dq[0] * m[0] + dq[1] * m[1]);
        let total = m[0] + m[1];
        if rain > 0.0 {
            // Deep: heating equals the latent heat released.
            let c = (l * rain - cp * (dt[0] * m[0] + dt[1] * m[1])) / (cp * total);
            dt = dt.map(|d| d + c);
            out.convective = rain;
        } else {
            // Shallow: redistribute heat and water, conserving both.
            let cq = -(dq[0] * m[0] + dq[1] * m[1]) / total;
            let ct = -(dt[0] * m[0] + dt[1] * m[1]) / total;
            dq = dq.map(|d| d + cq);
            dt = dt.map(|d| d + ct);
        }
        out.temperature = dt;
        out.humidity = dq;
    }

    // Large-scale condensation, per layer.
    let mut condensed = [0.0; 2];
    for k in 0..2 {
        let rh_c = p.condensation_humidity[k];
        let excess = q[k] - rh_c * sat[k].mean;
        if excess > 0.0 {
            let gain = 1.0 + rh_c * l / cp * saturation_slope(t[k], sat[k].mean);
            let rate = excess / gain / p.condensation_time;
            condensed[k] = rate * m[k];
            out.humidity[k] -= rate;
            out.temperature[k] += l * rate / cp;
        }
    }
    // The upper layer's rain re-evaporates partly in the lower layer, up to
    // the lower layer's condensation threshold.
    let room =
        ((p.condensation_humidity[1] * sat[1].mean - q[1]) * m[1] / p.condensation_time).max(0.0);
    let evaporated = (p.reevaporation * (1.0 - rh[1]).clamp(0.0, 1.0) * condensed[0]).min(room);
    out.humidity[1] += evaporated / m[1];
    out.temperature[1] -= l * evaporated / (cp * m[1]);

    // Rain or snow.
    let total = out.convective + condensed[0] + condensed[1] - evaporated;
    let (cold, warm) = p.snow_temperature;
    let snow_share = ((warm - column.surface_air) / (warm - cold)).clamp(0.0, 1.0);
    out.snow = snow_share * total;
    out.rain = total - out.snow;
    out.temperature[1] += FUSION_HEAT * out.snow / (cp * m[1]);

    // Clouds.
    let c = &p.clouds;
    let anvil = (c.convective * (out.convective * DAY).sqrt()).min(c.max_convective);
    let cover = |k: usize| stratiform_cover(rh[k], c.critical_humidity[k], c.overcast_humidity[k]);
    out.clouds = [1.0 - (1.0 - cover(0)) * (1.0 - anvil), cover(1)];
    out
}

/// The surface humidity (kg kg⁻¹) of the lower layer's profile for surface
/// air at temperature `surface_air`, capped at saturation.
pub fn surface_humidity(
    setting: &Setting,
    ps: f64,
    temperature: [f64; 2],
    q_lower: f64,
    surface_air: f64,
) -> f64 {
    let sat = setting.saturation(ps, temperature);
    let q_star = saturation_humidity(surface_air, ps);
    (q_lower.max(0.0) / sat[1].shaped * q_star).min(q_star)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setting<'p>(params: &'p MoistureParams, radiation: &'p RadiationParams) -> Setting<'p> {
        Setting {
            params,
            radiation,
            gas: Gas::EARTH_AIR,
            sigma: 0.5,
            gravity: 9.81,
        }
    }

    fn column(t: [f64; 2], rh: [f64; 2], surface_air: f64, setting: &Setting) -> MoistColumn {
        let sat = setting.saturation(1e5, t);
        MoistColumn {
            ps: 1e5,
            temperature: t,
            humidity: [rh[0] * sat[0].mean, rh[1] * sat[1].mean],
            surface_air,
            stratosphere: 210.0,
        }
    }

    /// Tropical, unstable; midlatitude, stable and saturated; polar, snowing.
    fn columns(s: &Setting) -> [MoistColumn; 3] {
        [
            column([228.0, 290.0], [0.3, 0.8], 300.0, s),
            column([232.0, 272.0], [0.6, 0.9], 285.0, s),
            column([215.0, 250.0], [0.5, 0.85], 262.0, s),
        ]
    }

    #[test]
    fn columns_conserve_energy_and_water() {
        let (p, r) = (MoistureParams::default(), RadiationParams::default());
        let s = setting(&p, &r);
        let cp = s.gas.heat_capacity;
        for c in columns(&s) {
            let out = moist_physics(&s, &c);
            let m = s.masses(c.ps);
            let heating = cp * (m[0] * out.temperature[0] + m[1] * out.temperature[1]);
            let water = m[0] * out.humidity[0] + m[1] * out.humidity[1];
            let precipitation = out.rain + out.snow;
            let latent = LATENT_HEAT * precipitation + FUSION_HEAT * out.snow;
            assert!(precipitation >= 0.0 && out.rain >= 0.0 && out.snow >= 0.0);
            // To round-off of the layers' (partly cancelling) terms.
            let heat_scale =
                cp * (m[0] * out.temperature[0].abs() + m[1] * out.temperature[1].abs());
            let water_scale = m[0] * out.humidity[0].abs() + m[1] * out.humidity[1].abs();
            assert!(
                (heating - latent).abs() < 1e-12 * heat_scale.max(1e-3),
                "{heating} vs {latent}"
            );
            assert!(
                (water + precipitation).abs() < 1e-12 * water_scale.max(1e-12),
                "{water} vs {precipitation}"
            );
        }
    }

    #[test]
    fn convection_rains_where_the_column_is_unstable() {
        let (p, r) = (MoistureParams::default(), RadiationParams::default());
        let s = setting(&p, &r);
        let [tropical, midlatitude, polar] = columns(&s).map(|c| moist_physics(&s, &c));
        // Warm, moist surface air under a cold upper layer: deep convection
        // dries and warms the column and rains.
        assert!(tropical.convective > 0.0 && tropical.temperature[0] > 0.0);
        assert!(tropical.clouds[0] > 0.0);
        // A stable column doesn't convect; saturated layers condense.
        assert_eq!(midlatitude.convective, 0.0);
        assert!(midlatitude.rain > 0.0 && midlatitude.clouds[1] > 0.5);
        // Below freezing at the surface, it snows.
        assert!(polar.snow > 0.0 && polar.rain == 0.0);
    }

    #[test]
    fn parcels_follow_the_dry_then_the_moist_adiabat() {
        let gas = Gas::EARTH_AIR;
        // Dry: θ conserved.
        let (t, q) = lift(gas, 300.0, 0.0, 1e5, 5e4);
        assert!((t / gas.exner(5e4) - 300.0).abs() < 0.01, "{t}");
        assert_eq!(q, 0.0);
        // Saturated from the start: warmer than the dry adiabat, and keeps
        // only what saturation allows.
        let q0 = saturation_humidity(300.0, 1e5);
        let (t, q) = lift(gas, 300.0, q0, 1e5, 5e4);
        assert!(t > 300.0 * gas.exner(5e4) + 15.0, "{t}");
        assert!((q - saturation_humidity(t, 5e4)).abs() < 1e-4 * q0);
        // Saturated at 300 K and 1000 hPa (θₑ ≈ 360 K): ~274 K at 500 hPa.
        assert!((270.0..280.0).contains(&t), "{t}");
    }

    #[test]
    fn cloud_cover_grows_with_humidity() {
        let c = CloudParams::default();
        let cover = |rh| stratiform_cover(rh, c.critical_humidity[1], c.overcast_humidity[1]);
        assert_eq!(cover(0.2), 0.0);
        assert!(cover(0.5) < cover(0.7) && cover(0.7) < 1.0);
        assert_eq!(cover(1.2), 1.0);
    }
}
