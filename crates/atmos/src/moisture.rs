//! Moist physics of one layered column: convection, large-scale
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
//! A parcel of surface air (the lowest layer extrapolated to the ground,
//! with the profile's surface humidity) is lifted: dry-adiabatically to its
//! condensation level, then along the pseudo-adiabat. Its temperature,
//! averaged over each layer, is the reference profile; in the top layer it
//! is bounded below by the stratosphere's temperature, taken as the skin
//! temperature `(OLR/2σ)^{1/4}` of the column. If a layer above the lowest
//! is colder than its reference, the column is unstable: every layer from
//! the highest such one (the cloud top) down to the surface relaxes towards
//! the reference temperature and a reference humidity over
//! [`MoistureParams::convection_time`] (Betts 1986; Frierson 2007). The
//! reference humidity has the in-layer profile's shape with
//! [`MoistureParams::convection_humidity`] at the surface, so a layer that
//! includes the dry stratosphere is referred to a far lower mean relative
//! humidity than one near the surface:
//!
//!
//! Temperature relaxes only in the cloud layers, all but the lowest: the
//! subcloud layer (the boundary layer) keeps its temperature, as in Betts
//! (1986), and only its humidity relaxes.
//!
//! - **Deep** (the relaxation removes water): it rains out, and the
//!   cloud layers' reference temperature shifts uniformly so their heating
//!   is exactly the latent heat released (enthalpy plus latent energy is
//!   conserved).
//! - **Shallow** (it would add water): no rain. Temperature is only
//!   redistributed between the cloud layers (enthalpy conserved), and
//!   water only moves from layers above their reference humidity to those
//!   below it, at the rate the givers supply: none is created, and no layer
//!   dries faster than its own relaxation, so none is drawn below zero.
//!
//! The moist adiabat's lapse rate falls as the climate warms, so the upper
//! layers warm more than the surface: the lapse-rate feedback.
//!
//! # Large-scale condensation
//!
//! A layer more humid than its condensation threshold
//! ([`MoistureParams::condensation_humidity`]) condenses the excess over
//! [`MoistureParams::condensation_time`], with its latent heat warming the
//! layer (allowing for the warming raising the saturation, as an implicit
//! step would). The rain falls through the layers below, each of which
//! re-evaporates a share proportional to its dryness.
//!
//! # Rain or snow
//!
//! Precipitation reaching the ground is snow below freezing at the surface
//! air temperature (a linear ramp over a few kelvin). Snow releases its
//! latent heat of fusion in the lowest layer; the surface pays it back when
//! the snow melts.
//!
//! # Clouds
//!
//! Two sheets, high and low ([`CloudParams::sigma`]), each diagnosed from
//! the relative humidity of the layer it lies in. Stratiform cover grows
//! linearly with that humidity between two bounds,
//! `c = (RH − RH₁)/(RH₂ − RH₁)` clamped to [0, 1] (SPEEDY (Molteni 2003)
//! uses the square, for the most humid level; layer means vary less, and
//! the square can't give Earth's ~0.5 low cloud at its layer-mean humidity
//! ~0.63 while still leaving a condensing layer short of overcast). As in
//! SPEEDY, deep convection adds anvil cloud to the high sheet, growing with
//! the square root of its rain.

use crate::layered::Gas;
use crate::levels::{Levels, MAX_LAYERS};
use crate::radiation::{
    LayerSaturation, RadiationParams, layer_saturation, saturation_humidity, saturation_pressure,
    specific_humidity,
};

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
    /// Level σ of the high and low cloud sheets; each takes its cover from
    /// the layer it lies in. Keep them off the layer interfaces: a sheet on
    /// one belongs to the layer above it, which at the climate's σ = 0.25 is
    /// the nearly dry stratospheric layer.
    pub sigma: [f64; 2],
    /// Albedo of overcast high and low cloud.
    pub albedo: [f64; 2],
    /// Longwave emissivity of overcast high and low cloud.
    pub emissivity: [f64; 2],
    /// Layer relative humidities at which stratiform cloud starts, high
    /// and low sheet.
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
            sigma: [0.3, 0.9],
            albedo: [0.4, 0.5],
            emissivity: [0.8, 1.0],
            critical_humidity: [0.45, 0.5],
            overcast_humidity: [0.85, 1.0],
            convective: 0.15,
            max_convective: 0.6,
        }
    }
}

/// σ at the centre of the climate's boundary layer (σ 0.85–1).
pub const BOUNDARY_LAYER_CENTRE: f64 = 0.925;

/// Moist-physics parameters. [`Default`] is Earth-like.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MoistureParams {
    /// Layer relative humidity above which water condenses (below 1: a
    /// layer saturates in places before its mean does), for layers centred
    /// at σ = ¼, ¾ and [`BOUNDARY_LAYER_CENTRE`] (the climate's boundary
    /// layer, which holds near 0.7 without raining); layers in between
    /// interpolate linearly in the σ of their centre, and layers beyond take
    /// the nearer value. Two equal-mass layers get exactly the first two.
    pub condensation_humidity: [f64; 3],
    /// Large-scale condensation time (s).
    pub condensation_time: f64,
    /// Bottom-layer relative humidity below which air climbing a slope
    /// doesn't rain.
    pub orographic_humidity: f64,
    /// Surface relative humidity of convection's reference profile (the
    /// Manabe–Wetherald shape above it).
    pub convection_humidity: f64,
    /// Convective relaxation time (s).
    pub convection_time: f64,
    /// Share of the rain falling into it that a bone-dry layer re-evaporates.
    pub reevaporation: f64,
    /// Surface air temperatures (K) below which all precipitation is snow,
    /// and above which all is rain.
    pub snow_temperature: (f64, f64),
    /// Clouds.
    pub clouds: CloudParams,
}

impl MoistureParams {
    /// The condensation threshold of a layer spanning `[top, bottom]` in σ
    /// ([`condensation_humidity`](Self::condensation_humidity)).
    pub fn condensation_threshold(&self, top: f64, bottom: f64) -> f64 {
        let centre = 0.5 * (top + bottom);
        let [upper, lower, boundary] = self.condensation_humidity;
        if centre <= 0.75 {
            let x = ((centre - 0.25) / 0.5).clamp(0.0, 1.0);
            upper * (1.0 - x) + lower * x
        } else {
            let x = ((centre - 0.75) / (BOUNDARY_LAYER_CENTRE - 0.75)).min(1.0);
            lower * (1.0 - x) + boundary * x
        }
    }
}

impl Default for MoistureParams {
    fn default() -> Self {
        Self {
            condensation_humidity: [0.45, 0.7, 0.9],
            condensation_time: 4.0 * 3600.0,
            orographic_humidity: 0.5,
            convection_humidity: 0.7,
            convection_time: 2.0 * 3600.0,
            reevaporation: 0.5,
            snow_temperature: (271.15, 275.15),
            clouds: CloudParams::default(),
        }
    }
}

/// One column's state, as the moist physics sees it.
#[derive(Clone, Debug)]
pub struct MoistColumn {
    /// Surface pressure (Pa).
    pub ps: f64,
    /// Layer temperatures, top first (K).
    pub temperature: Vec<f64>,
    /// Layer specific humidities (kg kg⁻¹).
    pub humidity: Vec<f64>,
    /// Surface air temperature (K).
    pub surface_air: f64,
    /// Surface air specific humidity (kg kg⁻¹): the convective parcel's.
    pub surface_humidity: f64,
    /// Stratosphere temperature (K): the floor of the top layer's
    /// convective reference.
    pub stratosphere: f64,
    /// Upslope ascent forcing rain out of the bottom layer's air (m s⁻¹; zero
    /// on a lee slope, or over smooth terrain).
    pub orographic_ascent: f64,
}

/// The moist physics' tendencies and diagnostics for one column.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MoistTendency {
    /// Temperature tendency per layer (K s⁻¹).
    pub temperature: Vec<f64>,
    /// Specific-humidity tendency per layer (s⁻¹).
    pub humidity: Vec<f64>,
    /// Rain and snow reaching the surface (kg m⁻² s⁻¹).
    pub rain: f64,
    /// See `rain`.
    pub snow: f64,
    /// The part of the precipitation from convection (kg m⁻² s⁻¹).
    pub convective: f64,
    /// The part of the precipitation from air forced over sub-grid hills
    /// (kg m⁻² s⁻¹).
    pub orographic: f64,
    /// Water the convection adds to the lowest layer (negative if it
    /// dries it; kg m⁻² s⁻¹).
    pub convective_boundary: f64,
    /// Large-scale condensation per layer (kg m⁻² s⁻¹).
    pub condensed: Vec<f64>,
    /// Re-evaporation of falling rain per layer (kg m⁻² s⁻¹).
    pub reevaporated: Vec<f64>,
    /// Cloud cover of the high and low sheets.
    pub clouds: [f64; 2],
    /// Layer relative humidities.
    pub relative_humidity: Vec<f64>,
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
    /// The layers.
    pub levels: &'p Levels,
    /// Gravity (m s⁻²).
    pub gravity: f64,
}

impl Setting<'_> {
    fn saturation(&self, ps: f64, t: &[f64]) -> Vec<LayerSaturation> {
        layer_saturation(self.radiation, self.gas, self.levels, ps, self.gravity, t)
    }

    fn masses(&self, ps: f64) -> Vec<f64> {
        self.levels
            .dsigma()
            .iter()
            .map(|ds| ds * ps / self.gravity)
            .collect()
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

/// Intervals of [`SaturationTable`], their width (K) and where they start:
/// 148.15 to 353.15 K, with 0 °C on a boundary, where the formula changes.
const TABLE_INTERVALS: usize = 1640;
const TABLE_STEP: f64 = 0.125;
/// Intervals below 0 °C.
const TABLE_ICE: usize = 1000;
const TABLE_START: f64 = 273.15 - TABLE_ICE as f64 * TABLE_STEP;

/// [`saturation_pressure`] tabulated for [`lift`], which spends most of its
/// time there: per interval, the values and (scaled) slopes at its ends, for
/// cubic Hermite interpolation. Agrees with the formula to a few parts in
/// 10⁹ (at the cold end; closer when warmer).
struct SaturationTable(Vec<[f64; 4]>);

impl SaturationTable {
    fn get() -> &'static Self {
        static TABLE: std::sync::OnceLock<SaturationTable> = std::sync::OnceLock::new();
        TABLE.get_or_init(|| {
            // The Magnus forms of `saturation_pressure` and their slopes,
            // each interval with the branch it lies in.
            let end = |i: usize, ice: bool| -> (f64, f64) {
                let c = TABLE_START + TABLE_STEP * i as f64 - 273.15;
                let (e0, a, b) = if ice {
                    (611.21, 22.587, 273.86)
                } else {
                    (610.94, 17.625, 243.04)
                };
                let e = e0 * fms_math::exp(a * c / (c + b));
                (e, TABLE_STEP * e * a * b / ((c + b) * (c + b)))
            };
            Self(
                (0..TABLE_INTERVALS)
                    .map(|i| {
                        let ice = i < TABLE_ICE;
                        let ((f0, d0), (f1, d1)) = (end(i, ice), end(i + 1, ice));
                        [f0, f1, d0, d1]
                    })
                    .collect(),
            )
        })
    }

    /// Saturation vapour pressure (Pa) at `t` (K); the formula itself
    /// outside the table.
    fn pressure(&self, t: f64) -> f64 {
        let x = (t - TABLE_START) / TABLE_STEP;
        if !(x >= 0.0 && x < TABLE_INTERVALS as f64) {
            return saturation_pressure(t);
        }
        let i = x as usize;
        let u = x - i as f64;
        let [f0, f1, d0, d1] = self.0[i];
        let (u2, u3) = (u * u, u * u * u);
        f0 * (2.0 * u3 - 3.0 * u2 + 1.0)
            + d0 * (u3 - 2.0 * u2 + u)
            + f1 * (3.0 * u2 - 2.0 * u3)
            + d1 * (u3 - u2)
    }
}

/// Temperature of a parcel lifted from `(t, q)` at pressure `p` to pressure
/// `p_to < p`: dry adiabat (q conserved) until saturated, then the
/// pseudo-adiabat `dT/d ln p = (R T + L q*)/(c_p + L² q*/(R_v T²))`
/// (midpoint steps of at most 0.05 in ln p). Returns the temperature and
/// humidity at `p_to`.
pub fn lift(gas: Gas, t: f64, q: f64, p: f64, p_to: f64) -> (f64, f64) {
    let (t, q, _) = lift_from(gas, t, q, p, p_to, None);
    (t, q)
}

/// [`lift`], given the saturation humidity at the start if it is known, and
/// returning the one at the end if the last step computed it.
fn lift_from(
    gas: Gas,
    t: f64,
    q: f64,
    p: f64,
    p_to: f64,
    known: Option<f64>,
) -> (f64, f64, Option<f64>) {
    let (r, cp) = (gas.gas_constant, gas.heat_capacity);
    let span = fms_math::ln(p / p_to);
    if span <= 0.0 {
        return (t, q, known);
    }
    let steps = (span / 0.05).ceil().max(1.0);
    let h = span / steps;
    // Pressure falls by this factor every half step.
    let half = fms_math::exp(-0.5 * h);
    let steps = steps as u64;
    let table = SaturationTable::get();
    let saturation_humidity = |t: f64, p: f64| specific_humidity(table.pressure(t), p);
    let (mut t, mut q, mut p, mut known) = (t, q, p, known);
    let slope = |t: f64, qs: f64, q: f64| -> (f64, bool) {
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
    for i in 0..steps {
        let qs = known.take().unwrap_or_else(|| saturation_humidity(t, p));
        let (s1, _) = slope(t, qs, q);
        let (t_mid, p_mid) = (t - 0.5 * h * s1, p * half);
        let (s2, saturated) = slope(t_mid, saturation_humidity(t_mid, p_mid), q);
        t -= h * s2;
        // The last step ends exactly at `p_to`.
        p = if i + 1 == steps { p_to } else { p_mid * half };
        if saturated {
            let qs = saturation_humidity(t, p);
            q = q.min(qs);
            known = Some(qs);
        }
    }
    (t, q, known)
}

/// The convective reference: layer-mean temperatures of a surface-air
/// parcel lifted through the column (the top layer's bounded below by
/// `column.stratosphere`).
pub fn convective_reference(setting: &Setting, column: &MoistColumn, surface_q: f64) -> Vec<f64> {
    let ps = column.ps;
    let n = setting.levels.count();
    // Every layer's nodes, lifted through bottom up.
    let mut nodes = [(0usize, 0.0f64, 0.0f64); 4 * MAX_LAYERS];
    let nodes = &mut nodes[..4 * n];
    for k in 0..n {
        let (a, b) = setting.levels.bounds(k);
        for (i, (x, w)) in GAUSS.into_iter().enumerate() {
            nodes[4 * k + i] = (k, 0.5 * (a + b) + 0.5 * (b - a) * x, 0.5 * w);
        }
    }
    nodes.sort_by(|a, b| b.1.total_cmp(&a.1));
    let (mut t, mut q, mut p) = (column.surface_air, surface_q, ps);
    let mut reference = vec![0.0; n];
    // The parcel only cools as it rises: once it is colder than the
    // stratosphere, the rest of the top layer is at the floor.
    let mut floored = false;
    let mut known = None;
    for &(k, s, w) in nodes.iter() {
        if floored {
            reference[k] += w * column.stratosphere;
            continue;
        }
        let p_to = s * ps;
        (t, q, known) = lift_from(setting.gas, t, q, p, p_to, known);
        p = p_to;
        let value = if k == 0 {
            floored = t <= column.stratosphere;
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
    let levels = setting.levels;
    let n = levels.count();
    let bottom = n - 1;
    let m = setting.masses(ps);
    let t = &column.temperature;
    let q: Vec<f64> = column.humidity.iter().map(|q| q.max(0.0)).collect();
    let sat = setting.saturation(ps, t);
    let rh: Vec<f64> = (0..n).map(|k| q[k] / sat[k].mean).collect();
    let mut out = MoistTendency {
        temperature: vec![0.0; n],
        humidity: vec![0.0; n],
        relative_humidity: rh.clone(),
        ..MoistTendency::default()
    };

    // Convection, from the highest layer colder than the parcel down.
    let surface_q = column
        .surface_humidity
        .clamp(0.0, saturation_humidity(column.surface_air, ps));
    let reference = convective_reference(setting, column, surface_q);
    if let Some(top) = (0..bottom).find(|&k| reference[k] > t[k]) {
        let tau = p.convection_time;
        let ref_sat = setting.saturation(ps, &reference);
        let layers = top..n;
        // The subcloud (lowest) layer keeps its temperature: its reference
        // includes moist-adiabatic warming above the condensation level
        // that the well-mixed layer doesn't have.
        let mut dt: Vec<f64> = layers
            .clone()
            .map(|k| {
                if k == bottom {
                    0.0
                } else {
                    (reference[k] - t[k]) / tau
                }
            })
            .collect();
        let mut dq: Vec<f64> = layers
            .clone()
            .map(|k| (p.convection_humidity * ref_sat[k].shaped - q[k]) / tau)
            .collect();
        let mc = &m[top..];
        let weighted = |x: &[f64]| -> f64 {
            let mut s = x[0] * mc[0];
            for (x, m) in x.iter().zip(mc).skip(1) {
                s += x * m;
            }
            s
        };
        let rain = -weighted(&dq);
        let mut total = mc[0];
        for m in &mc[1..] {
            total += m;
        }
        if rain > 0.0 {
            // Deep: heating equals the latent heat released, all of it in
            // the cloud layers above the lowest (the subcloud layer keeps
            // its own relaxation; Betts 1986).
            let cloud = total - mc[mc.len() - 1];
            let c = (l * rain - cp * weighted(&dt)) / (cp * cloud);
            let last = dt.len() - 1;
            dt[..last].iter_mut().for_each(|d| *d += c);
            out.convective = rain;
        } else {
            // Shallow: redistribute heat and water, conserving both. Water
            // moves only from layers above their reference to those below
            // it, at the rate the givers supply, so no layer dries faster
            // than its own relaxation (and none is drawn below zero).
            let (mut give, mut take) = (0.0, 0.0);
            for (d, m) in dq.iter().zip(mc) {
                if *d < 0.0 {
                    give -= d * m;
                } else {
                    take += d * m;
                }
            }
            let scale = if take > 0.0 { give / take } else { 0.0 };
            dq.iter_mut()
                .filter(|d| **d > 0.0)
                .for_each(|d| *d *= scale);
            let ct = -weighted(&dt) / (total - mc[mc.len() - 1]);
            let last = dt.len() - 1;
            dt[..last].iter_mut().for_each(|d| *d += ct);
        }
        out.convective_boundary = dq[dq.len() - 1] * mc[mc.len() - 1];
        out.temperature[top..].copy_from_slice(&dt);
        out.humidity[top..].copy_from_slice(&dq);
    }

    // Large-scale condensation, per layer.
    let threshold: Vec<f64> = (0..n)
        .map(|k| {
            let (a, b) = levels.bounds(k);
            p.condensation_threshold(a, b)
        })
        .collect();
    let mut condensed = vec![0.0; n];
    for k in 0..n {
        let rh_c = threshold[k];
        let excess = q[k] - rh_c * sat[k].mean;
        if excess > 0.0 {
            let gain = 1.0 + rh_c * l / cp * saturation_slope(t[k], sat[k].mean);
            let rate = excess / gain / p.condensation_time;
            condensed[k] = rate * m[k];
            out.humidity[k] -= rate;
            out.temperature[k] += l * rate / cp;
        }
    }
    // Air forced up a slope rains out in proportion to the moisture it
    // carries up, `ε ρ q w` (Smith 1979), once it's humid enough, and no
    // faster than the large-scale condensation.
    if column.orographic_ascent > 0.0 {
        let density = ps / (gas.gas_constant * column.surface_air);
        let humid =
            ((rh[bottom] - p.orographic_humidity) / (1.0 - p.orographic_humidity)).clamp(0.0, 1.0);
        let rate = q[bottom]
            * humid
            * (column.orographic_ascent * density / m[bottom]).min(1.0 / p.condensation_time);
        out.orographic = rate * m[bottom];
        out.humidity[bottom] -= rate;
        out.temperature[bottom] += l * rate / cp;
    }
    // Rain from above re-evaporates partly in each layer it falls through,
    // up to that layer's condensation threshold.
    let mut evaporated = vec![0.0; n];
    let mut falling = condensed[0];
    for k in 1..n {
        let room = ((threshold[k] * sat[k].mean - q[k]) * m[k] / p.condensation_time).max(0.0);
        evaporated[k] = (p.reevaporation * (1.0 - rh[k]).clamp(0.0, 1.0) * falling).min(room);
        out.humidity[k] += evaporated[k] / m[k];
        out.temperature[k] -= l * evaporated[k] / (cp * m[k]);
        falling += condensed[k] - evaporated[k];
    }

    // Rain or snow.
    let mut total = out.convective + out.orographic;
    for c in &condensed {
        total += c;
    }
    for e in &evaporated {
        total -= e;
    }
    let (cold, warm) = p.snow_temperature;
    let snow_share = ((warm - column.surface_air) / (warm - cold)).clamp(0.0, 1.0);
    out.snow = snow_share * total;
    out.rain = total - out.snow;
    out.temperature[bottom] += FUSION_HEAT * out.snow / (cp * m[bottom]);
    out.condensed = condensed;
    out.reevaporated = evaporated;

    // Clouds.
    let c = &p.clouds;
    let anvil = (c.convective * (out.convective * DAY).sqrt()).min(c.max_convective);
    let cover = |sheet: usize| {
        let k = levels.layer_of(c.sigma[sheet]);
        stratiform_cover(
            rh[k],
            c.critical_humidity[sheet],
            c.overcast_humidity[sheet],
        )
    };
    out.clouds = [1.0 - (1.0 - cover(0)) * (1.0 - anvil), cover(1)];
    out
}

/// The surface humidity (kg kg⁻¹) of the lowest layer's profile for surface
/// air at temperature `surface_air`, capped at saturation.
pub fn surface_humidity(
    setting: &Setting,
    ps: f64,
    temperature: &[f64],
    q_lowest: f64,
    surface_air: f64,
) -> f64 {
    let sat = setting.saturation(ps, temperature);
    let q_star = saturation_humidity(surface_air, ps);
    (q_lowest.max(0.0) / sat[sat.len() - 1].shaped * q_star).min(q_star)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two() -> &'static Levels {
        static L: std::sync::OnceLock<Levels> = std::sync::OnceLock::new();
        L.get_or_init(|| Levels::two(0.5))
    }

    fn setting<'p>(
        params: &'p MoistureParams,
        radiation: &'p RadiationParams,
        levels: &'p Levels,
    ) -> Setting<'p> {
        Setting {
            params,
            radiation,
            gas: Gas::EARTH_AIR,
            levels,
            gravity: 9.81,
        }
    }

    /// A column with layer temperatures `temperature` and relative
    /// humidities interpolated from `rh` (at σ = ¼ and ¾) at each layer's
    /// centre.
    fn column(temperature: Vec<f64>, rh: [f64; 2], surface_air: f64, s: &Setting) -> MoistColumn {
        let levels = s.levels;
        let sat = s.saturation(1e5, &temperature);
        let humidity = (0..levels.count())
            .map(|k| {
                let (a, b) = levels.bounds(k);
                let x = ((0.5 * (a + b) - 0.25) / 0.5).clamp(0.0, 1.0);
                (rh[0] * (1.0 - x) + rh[1] * x) * sat[k].mean
            })
            .collect::<Vec<f64>>();
        let bottom = levels.count() - 1;
        MoistColumn {
            ps: 1e5,
            surface_humidity: surface_humidity(s, 1e5, &temperature, humidity[bottom], surface_air),
            humidity,
            temperature,
            surface_air,
            stratosphere: 210.0,
            orographic_ascent: 0.0,
        }
    }

    /// Layer means of the temperature profile `T_s σ^β`, no colder than
    /// `floor`.
    fn profile(s: &Setting, surface: f64, beta: f64, floor: f64) -> Vec<f64> {
        (0..s.levels.count())
            .map(|k| {
                let (a, b) = s.levels.bounds(k);
                GAUSS
                    .iter()
                    .map(|&(x, w)| {
                        let sigma = 0.5 * (a + b) + 0.5 * (b - a) * x;
                        0.5 * w * (surface * fms_math::pow(sigma, beta)).max(floor)
                    })
                    .sum()
            })
            .collect()
    }

    /// Tropical, unstable (every layer 3 K colder than the lifted surface
    /// parcel); midlatitude, stable and saturated; polar, snowing.
    fn columns(s: &Setting) -> [MoistColumn; 3] {
        let tropical = {
            let (rh, surface_air) = ([0.3, 0.8], 300.0);
            let c = column(profile(s, 300.0, 0.19, 200.0), rh, surface_air, s);
            let reference = convective_reference(s, &c, c.surface_humidity);
            column(
                reference.iter().map(|t| t - 3.0).collect(),
                rh,
                surface_air,
                s,
            )
        };
        [
            tropical,
            column(profile(s, 285.0, 0.1, 215.0), [0.6, 0.9], 285.0, s),
            column(profile(s, 262.0, 0.08, 205.0), [0.5, 0.85], 262.0, s),
        ]
    }

    fn layer_sets() -> Vec<Levels> {
        vec![
            Levels::two(0.5),
            Levels::new(&[0.4, 0.85]),
            Levels::equal_mass(5),
        ]
    }

    #[test]
    // Layers index the masses and both tendencies together.
    #[allow(clippy::needless_range_loop)]
    fn columns_conserve_energy_and_water() {
        let (p, r) = (MoistureParams::default(), RadiationParams::default());
        for levels in layer_sets() {
            let s = setting(&p, &r, &levels);
            let cp = s.gas.heat_capacity;
            for c in columns(&s) {
                let out = moist_physics(&s, &c);
                let m = s.masses(c.ps);
                let (mut heating, mut water, mut heat_scale, mut water_scale) =
                    (0.0, 0.0, 0.0, 0.0);
                for k in 0..levels.count() {
                    heating += cp * m[k] * out.temperature[k];
                    water += m[k] * out.humidity[k];
                    heat_scale += cp * m[k] * out.temperature[k].abs();
                    water_scale += m[k] * out.humidity[k].abs();
                }
                let precipitation = out.rain + out.snow;
                let latent = LATENT_HEAT * precipitation + FUSION_HEAT * out.snow;
                assert!(precipitation >= 0.0 && out.rain >= 0.0 && out.snow >= 0.0);
                // To round-off of the layers' (partly cancelling) terms.
                assert!(
                    (heating - latent).abs() < 1e-12 * heat_scale.max(1e-3),
                    "{} layers: {heating} vs {latent}",
                    levels.count()
                );
                assert!(
                    (water + precipitation).abs() < 1e-12 * water_scale.max(1e-12),
                    "{} layers: {water} vs {precipitation}",
                    levels.count()
                );
            }
        }
    }

    #[test]
    fn orographic_rain_takes_its_water_and_heat_from_the_bottom_layer() {
        let (p, r) = (MoistureParams::default(), RadiationParams::default());
        for levels in layer_sets() {
            let s = setting(&p, &r, &levels);
            let [_, base, _] = columns(&s);
            let bottom = levels.count() - 1;
            let mass = levels.dsigma()[bottom] * base.ps / s.gravity;
            let lifted = MoistColumn {
                orographic_ascent: 0.02,
                ..base.clone()
            };
            let (with, without) = (moist_physics(&s, &lifted), moist_physics(&s, &base));
            assert_eq!(without.orographic, 0.0);
            assert!(with.orographic > 0.0);
            let total = |t: &MoistTendency| t.rain + t.snow;
            assert!((total(&with) - total(&without) - with.orographic).abs() < 1e-15);
            let dried = (without.humidity[bottom] - with.humidity[bottom]) * mass;
            assert!((dried - with.orographic).abs() < 1e-15 * mass);
            let warmed = (with.temperature[bottom] - without.temperature[bottom])
                * s.gas.heat_capacity
                * mass;
            let expected = (LATENT_HEAT
                + FUSION_HEAT * (with.snow - without.snow) / with.orographic)
                * with.orographic;
            assert!((warmed - expected).abs() < 1e-9 * expected.abs());
            // Too dry to rain, however hard it climbs.
            let t = base.temperature.clone();
            let air = base.surface_air;
            let dry = MoistColumn {
                orographic_ascent: 0.02,
                ..column(t, [0.1, 0.3], air, &s)
            };
            assert_eq!(moist_physics(&s, &dry).orographic, 0.0);
        }
    }

    #[test]
    fn convection_rains_where_the_column_is_unstable() {
        let (p, r) = (MoistureParams::default(), RadiationParams::default());
        for levels in layer_sets() {
            let s = setting(&p, &r, &levels);
            let [tropical, midlatitude, polar] = columns(&s).map(|c| moist_physics(&s, &c));
            let n = levels.count();
            // Warm, moist surface air under cold upper layers: deep
            // convection dries the column, warms the cloud layers and
            // rains; convection only dries the subcloud layer.
            assert!(tropical.convective > 0.0, "{n}");
            assert!(
                tropical.temperature[..n - 1].iter().all(|&t| t > 0.0),
                "{n}"
            );
            assert!(tropical.humidity[n - 1] < 0.0, "{n}");
            assert!(tropical.clouds[0] > 0.0);
            // A stable column doesn't convect; saturated layers condense.
            assert_eq!(midlatitude.convective, 0.0, "{n}");
            assert!(midlatitude.rain > 0.0 && midlatitude.clouds[1] > 0.5, "{n}");
            // Below freezing at the surface, it snows.
            assert!(polar.snow > 0.0 && polar.rain == 0.0, "{n}");
        }
    }

    #[test]
    fn two_equal_layers_get_the_condensation_thresholds_exactly() {
        let p = MoistureParams::default();
        let l = two();
        let [upper, lower, boundary] = p.condensation_humidity;
        assert_eq!(
            p.condensation_threshold(l.bounds(0).0, l.bounds(0).1),
            upper
        );
        assert_eq!(
            p.condensation_threshold(l.bounds(1).0, l.bounds(1).1),
            lower
        );
        // In between, linear in σ; beyond, the nearer value.
        let mid = p.condensation_threshold(0.4, 0.6);
        assert!((mid - 0.5 * (upper + lower)).abs() < 1e-15);
        assert_eq!(p.condensation_threshold(0.85, 1.0), boundary);
        assert_eq!(p.condensation_threshold(0.0, 0.1), upper);
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
    fn the_saturation_table_matches_the_formula() {
        let table = SaturationTable::get();
        let mut worst = 0.0f64;
        for i in 0..=250_000 {
            let t = 140.0 + 220.0 * f64::from(i) / 250_000.0;
            let exact = saturation_pressure(t);
            worst = worst.max(((table.pressure(t) - exact) / exact).abs());
        }
        assert!(worst < 1e-8, "{worst:e}");
        assert_eq!(table.pressure(273.15), saturation_pressure(273.15));
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
