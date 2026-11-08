//! Semi-grey radiation for one column: shortwave absorption and reflection,
//! and two-stream longwave transfer in three spectral intervals.
//!
//! # Longwave
//!
//! The thermal spectrum (wavenumber `ν`, cm⁻¹) is split into intervals, each
//! grey:
//!
//! - **CO₂ band**: `667 ± W/2`, opaque. Its width follows the
//!   exponential-wing band model. CO₂ line strength falls off as
//!   `exp(−|ν − ν₀|/l)` away from the band centre, so the band is opaque out to
//!   where `κ(ν) · u = 1`, with `u` the pressure-broadened absorber amount
//!   `∝ C p_s² / g`. Hence `W ≈ W₀ + 2l ln(C p_s² g₀ / (C₀ p₀² g))`
//!   ([`co2_band_width`]), and the forcing is logarithmic in CO₂ because the
//!   band edges move outward into regions that emit from the warm lower
//!   atmosphere or surface (Wilson & Gea-Banacloche 2012; Jeevanjee et al.
//!   2021). `l` is calibrated, not the spectroscopic 10.2 cm⁻¹: two layers
//!   can't resolve the stratosphere's cooling and readjustment, which adds
//!   about a third to the instantaneous forcing, so `l` is set to give the
//!   adjusted forcing (~3.7 W m⁻² per doubling) at the model's own
//!   equilibrium.
//! - **Window**: 800–1200 cm⁻¹ (less any CO₂ band), absorbed only by the
//!   water-vapour continuum (self-broadened, `∝ q e`).
//! - **Water vapour**: the rest, in two parts: the cores of the rotation and
//!   vibration bands (below 550 and above 1400 cm⁻¹), and their much weaker
//!   edges either side of the window (550–800, 1200–1400), where the CO₂
//!   band's wings also lie. Within each, water's lines span orders of magnitude in
//!   strength, so the band is split into [`WATER_TERMS`] equal parts whose
//!   absorption coefficients are a decade apart (a minimal k-distribution),
//!   each acting on the pressure-weighted water column `∫ q (p/p₀) dp/g`.
//!   As water increases the strong parts saturate and the weak ones take
//!   over, so outgoing longwave falls roughly with `ln w`, as observed,
//!   rather than collapsing once a single grey band turns opaque.
//!
//! Each layer's Planck emission `B` (the interval's share of `σT⁴`, from
//! the exact series for the Planck integral) is treated as linear in optical
//! depth between its edges. The lower layer's edges follow the convective
//! lapse rate through its mean temperature: the interface above, the surface
//! air below. The upper layer runs from the interface to its own mean
//! temperature at the top (an isothermal, stratosphere-like cap). So in an
//! opaque band a layer emits upward from its top and downward from its
//! bottom: the lower layer warms the surface from its warm base, and
//! neighbouring layers don't spuriously exchange heat at their mean
//! temperatures. The surface is black. Optical depths are diffuse (the 1.66
//! factor is absorbed into the coefficients).
//!
//! # Water vapour
//!
//! Relative humidity is fixed (Manabe & Wetherald 1967):
//! `RH(σ) = RH_s (σ − 0.02)/0.98`. Saturation vapour pressure is Magnus over
//! water above freezing and over ice below. Within a layer, temperature
//! follows the convective lapse rate Γ through the layer's mean temperature,
//! `T(p) = Tₖ (p/pₖ)^{RΓ/g}` (pₖ: the pressure whose Exner function is the
//! layer mean); columns are integrated by Gauss–Legendre quadrature.
//!
//! # Shortwave
//!
//! Absorbed on the way down (slant path `1/μ`) and again after reflection
//! (diffuse path): ozone in the upper layer, and water vapour in both, in
//! the near-infrared only (a fixed share of sunlight) and saturating as
//! `√w`. A cloud deck at the layer interface reflects a fixed share (a
//! stand-in until the model diagnoses clouds), and the surface reflects its
//! own albedo; reflections between the two are summed exactly. Putting the
//! clouds' reflection at the surface instead, below the water vapour, would
//! make water's shortwave absorption feedback several times too strong.
//!
//! Every flux is in W m⁻² and every flux budget closes exactly: absorbed
//! sunlight minus outgoing longwave equals the net heating of the layers and
//! the surface.

use crate::two_layer::{Gas, REFERENCE_PRESSURE};

/// Stefan–Boltzmann constant (W m⁻² K⁻⁴).
pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;
/// Second radiation constant `hc/k` (cm K).
const C2: f64 = 1.438_776_877;
/// CO₂ band centre (cm⁻¹).
const CO2_CENTRE: f64 = 667.0;
/// Window edges (cm⁻¹).
const WINDOW: (f64, f64) = (800.0, 1200.0);
/// Water's molar mass over dry air's.
const EPSILON: f64 = 0.622;
/// Reference gravity for the CO₂ absorber amount (m s⁻²).
const G0: f64 = 9.81;

/// Radiation parameters. [`Default`] is calibrated to Earth.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadiationParams {
    /// CO₂ band width `W₀` (cm⁻¹) at the reference CO₂ and pressure.
    pub co2_width: f64,
    /// Band-wing decay length `l` (cm⁻¹): each band edge moves out by `l`
    /// per e-fold of CO₂.
    pub co2_wing: f64,
    /// Reference CO₂ (ppm) for `W₀`.
    pub co2_reference: f64,
    /// Optical depth of the CO₂ band per unit σ (large: opaque).
    pub co2_opacity: f64,
    /// Water-vapour band absorption coefficient (m² kg⁻¹ of
    /// pressure-weighted water).
    pub water_opacity: f64,
    /// Water's opacity at the band edges (550–800 and 1200–1400 cm⁻¹)
    /// relative to the cores.
    pub water_edge: f64,
    /// Window continuum coefficient (m² kg⁻¹ of water, per 1 kPa of vapour
    /// pressure).
    pub continuum_opacity: f64,
    /// Surface relative humidity `RH_s` of the Manabe–Wetherald profile.
    pub surface_humidity: f64,
    /// Shortwave absorption by ozone (optical depth, upper layer).
    pub ozone_absorption: f64,
    /// Shortwave absorption by water vapour: near-infrared optical depth
    /// per `√(kg m⁻²)`.
    pub water_absorption: f64,
    /// Share of sunlight in the near-infrared water-vapour bands.
    pub near_infrared: f64,
    /// Lapse rate assumed within a layer for the water-vapour profile
    /// (K m⁻¹).
    pub lapse_rate: f64,
    /// Effective cosine of the solar zenith angle for the downward beam.
    pub mu: f64,
    /// Surface albedo.
    pub albedo: f64,
    /// Albedo of a cloud deck at the layer interface: a fixed stand-in for
    /// clouds until the model diagnoses them.
    pub cloud_albedo: f64,
}

impl Default for RadiationParams {
    /// Longwave calibrated in a single column (`examples/rcm.rs`): an
    /// observed global-mean Earth column (surface 288 K, layers 272 and
    /// 226 K) radiates ~244 W m⁻² to space and loses ~2 W m⁻² more per
    /// kelvin of warming at fixed relative humidity, and the model's own
    /// equilibrium gains ~3.6 W m⁻² per CO₂ doubling. The cloud albedo is
    /// calibrated in 3-D (`examples/climate.rs`): five years at 2.5k cells
    /// give a global mean of ~287 K and a planetary albedo of ~0.29.
    fn default() -> Self {
        Self {
            co2_width: 160.0,
            co2_wing: 30.0,
            co2_reference: 400.0,
            co2_opacity: 50.0,
            water_opacity: 60.0,
            water_edge: 0.02,
            continuum_opacity: 0.03,
            surface_humidity: 0.77,
            ozone_absorption: 0.03,
            water_absorption: 0.04,
            near_infrared: 0.35,
            lapse_rate: 6.5e-3,
            mu: 0.5,
            albedo: 0.08,
            cloud_albedo: 0.33,
        }
    }
}

/// One column's inputs.
#[derive(Clone, Copy, Debug)]
pub struct Column {
    /// Surface pressure (Pa).
    pub ps: f64,
    /// Surface temperature (K).
    pub surface_temperature: f64,
    /// Layer potential temperatures, upper first (K).
    pub theta: [f64; 2],
    /// Interface σ.
    pub sigma: f64,
    /// Gravity (m s⁻²).
    pub gravity: f64,
    /// Top-of-atmosphere insolation (W m⁻²).
    pub insolation: f64,
    /// CO₂ (ppm).
    pub co2: f64,
}

/// One column's fluxes (W m⁻²) and water (kg m⁻²).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Fluxes {
    /// Absorbed sunlight (insolation minus reflected).
    pub absorbed: f64,
    /// Outgoing longwave.
    pub outgoing: f64,
    /// Net radiative heating of each layer (shortwave + longwave).
    pub layers: [f64; 2],
    /// Shortwave absorbed by the surface.
    pub surface_shortwave: f64,
    /// Downward longwave at the surface.
    pub surface_down: f64,
    /// Upward longwave from the surface.
    pub surface_up: f64,
    /// Water vapour in each layer.
    pub water: [f64; 2],
    /// Width of the opaque CO₂ band (cm⁻¹).
    pub co2_width: f64,
}

impl Fluxes {
    /// Net radiative heating of the surface.
    pub fn surface(&self) -> f64 {
        self.surface_shortwave + self.surface_down - self.surface_up
    }
}

/// Saturation vapour pressure (Pa): Magnus over water at or above 0 °C, over
/// ice below.
pub fn saturation_pressure(t: f64) -> f64 {
    let c = t - 273.15;
    if c >= 0.0 {
        610.94 * fms_math::exp(17.625 * c / (c + 243.04))
    } else {
        611.21 * fms_math::exp(22.587 * c / (c + 273.86))
    }
}

/// Specific humidity (kg kg⁻¹) at vapour pressure `e` and pressure `p`. The
/// vapour pressure is capped at a quarter of `p`: beyond that the dry-air
/// model this rests on stops making sense.
pub fn specific_humidity(e: f64, p: f64) -> f64 {
    let e = e.min(0.25 * p);
    EPSILON * e / (p - (1.0 - EPSILON) * e)
}

/// `∫ₓ^∞ t³/(eᵗ − 1) dt`: exponential series above 1, Bernoulli series
/// below (both converge to round-off).
fn planck_tail(x: f64) -> f64 {
    const TOTAL: f64 =
        std::f64::consts::PI * std::f64::consts::PI * std::f64::consts::PI * std::f64::consts::PI
            / 15.0;
    if x <= 0.0 {
        return TOTAL;
    }
    if x >= 1.0 {
        let terms = (2.0 + 40.0 / x) as usize;
        (1..=terms)
            .map(|n| {
                let n = n as f64;
                fms_math::exp(-n * x)
                    * (x * x * x / n
                        + 3.0 * x * x / (n * n)
                        + 6.0 * x / (n * n * n)
                        + 6.0 / (n * n * n * n))
            })
            .sum()
    } else {
        // t/(eᵗ − 1) = Σ Bₖ tᵏ/k!, so ∫₀ˣ t² · that = Σ Bₖ x^{k+3}/(k! (k+3)).
        const B: [(f64, i32); 10] = [
            (1.0, 0),
            (-0.5, 1),
            (1.0 / 12.0, 2),
            (-1.0 / 720.0, 4),
            (1.0 / 30_240.0, 6),
            (-1.0 / 1_209_600.0, 8),
            (1.0 / 47_900_160.0, 10),
            (-691.0 / 1_307_674_368_000.0, 12),
            (1.338_253_653_068_467_9e-11, 14),
            (-3.389_680_296_322_583e-13, 16),
        ];
        let head: f64 = B
            .iter()
            .map(|&(b, k)| b * x.powi(k + 3) / f64::from(k + 3))
            .sum();
        TOTAL - head
    }
}

/// Fraction of black-body emission at temperature `t` between wavenumbers
/// `a` and `b` (cm⁻¹).
pub fn planck_fraction(a: f64, b: f64, t: f64) -> f64 {
    let g = |nu: f64| {
        if nu.is_infinite() {
            0.0
        } else {
            planck_tail(C2 * nu / t)
        }
    };
    let total = planck_tail(0.0);
    ((g(a) - g(b)) / total).max(0.0)
}

/// Which absorber controls an interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Co2,
    Window,
    /// The cores of water's rotation (< 550 cm⁻¹) and vibration
    /// (> 1400 cm⁻¹) bands.
    Water,
    /// The weak edges of those bands either side of the window, where the
    /// CO₂ band's wings lie.
    WaterEdge,
}

/// Where water's band cores end (cm⁻¹).
const WATER_CORES: (f64, f64) = (550.0, 1400.0);

/// The spectral intervals for a CO₂ band of width `w`.
fn intervals(w: f64) -> Vec<(f64, f64, Kind)> {
    let (lo, hi) = ((CO2_CENTRE - 0.5 * w).max(0.0), CO2_CENTRE + 0.5 * w);
    // Breakpoints, then classify each piece by its midpoint.
    let mut cuts = vec![
        0.0,
        WATER_CORES.0,
        WINDOW.0,
        WINDOW.1,
        WATER_CORES.1,
        f64::INFINITY,
    ];
    if w > 0.0 {
        cuts.extend([lo, hi]);
    }
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    cuts.windows(2)
        .filter(|c| c[1] > c[0])
        .map(|c| {
            let mid = if c[1].is_infinite() {
                c[0] + 1.0
            } else {
                0.5 * (c[0] + c[1])
            };
            let kind = if w > 0.0 && (lo..hi).contains(&mid) {
                Kind::Co2
            } else if (WINDOW.0..WINDOW.1).contains(&mid) {
                Kind::Window
            } else if (WATER_CORES.0..WATER_CORES.1).contains(&mid) {
                Kind::WaterEdge
            } else {
                Kind::Water
            };
            (c[0], c[1], kind)
        })
        .collect()
}

/// The water band's k-distribution: (share of the band, absorption
/// coefficient relative to [`RadiationParams::water_opacity`]).
pub const WATER_TERMS: [(f64, f64); 4] = [(0.25, 10.0), (0.25, 1.0), (0.25, 0.1), (0.25, 0.01)];

/// `(1 − e^{−τ}(1 + τ))/τ`: the extra emission from a source that varies
/// linearly across a layer of optical depth `τ` (transmissivity `t`), per
/// unit difference between its far and near edges.
fn linear_source_term(tau: f64, t: f64) -> f64 {
    if tau < 1e-4 {
        tau * (0.5 - tau / 3.0)
    } else {
        (1.0 - t * (1.0 + tau)) / tau
    }
}

/// Width of the opaque CO₂ band (cm⁻¹):
/// `W = 2l ln(1 + (e^{W₀/2l} − 1) u)`, with `u` the pressure-broadened CO₂
/// amount relative to the reference (`C p_s² g₀ / (C₀ p₀² g)`). For a thick
/// band this is the band-wing law `W₀ + 2l ln u`; for traces of CO₂ it
/// narrows in proportion to the amount (the weak-line limit) instead of
/// vanishing.
pub fn co2_band_width(params: &RadiationParams, column: &Column) -> f64 {
    let p0 = REFERENCE_PRESSURE;
    let u =
        column.co2.max(0.0) / params.co2_reference * (column.ps / p0).powi(2) * G0 / column.gravity;
    let l2 = 2.0 * params.co2_wing;
    let scale = fms_math::exp(params.co2_width / l2) - 1.0;
    l2 * fms_math::ln(1.0 + scale * u)
}

/// Gauss–Legendre nodes on [−1, 1].
const GAUSS: [(f64, f64); 4] = [
    (-0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
    (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
];

/// Per layer: water (kg m⁻²), pressure-weighted water and continuum
/// absorber amount.
fn water_columns(
    params: &RadiationParams,
    gas: Gas,
    column: &Column,
    exner: [f64; 2],
) -> [(f64, f64, f64); 2] {
    let kappa = gas.kappa();
    let p0 = REFERENCE_PRESSURE;
    let beta = gas.gas_constant * params.lapse_rate / column.gravity;
    let bounds = [(0.0, column.sigma), (column.sigma, 1.0)];
    std::array::from_fn(|k| {
        let (a, b) = bounds[k];
        let t_mean = column.theta[k] * exner[k];
        let p_mean = p0 * fms_math::pow(exner[k], 1.0 / kappa);
        let (mut w, mut wp, mut wc) = (0.0, 0.0, 0.0);
        for (x, g) in GAUSS {
            let sigma = 0.5 * (a + b) + 0.5 * (b - a) * x;
            let p = sigma * column.ps;
            let t = t_mean * fms_math::pow(p / p_mean, beta);
            let rh = (params.surface_humidity * (sigma - 0.02) / 0.98).max(0.0);
            let e = rh * saturation_pressure(t);
            let q = specific_humidity(e, p);
            let mass = 0.5 * g * (b - a) * column.ps / column.gravity;
            w += mass * q;
            wp += mass * q * p / p0;
            wc += mass * q * e.min(0.25 * p) / 1000.0;
        }
        (w, wp, wc)
    })
}

/// Radiative fluxes of one column.
pub fn column_fluxes(params: &RadiationParams, gas: Gas, column: &Column) -> Fluxes {
    let dsigma = [column.sigma, 1.0 - column.sigma];
    let kappa = gas.kappa();
    let p0 = REFERENCE_PRESSURE;
    // Layer temperatures: θ times the layer-mean Exner function.
    let exner = {
        let pi = column.sigma * column.ps;
        let (ei, es) = (gas.exner(pi), gas.exner(column.ps));
        [
            ei / (1.0 + kappa),
            (column.ps * es - pi * ei) / ((1.0 + kappa) * (1.0 - column.sigma) * column.ps),
        ]
    };
    let water = water_columns(params, gas, column, exner);
    let temperature = [column.theta[0] * exner[0], column.theta[1] * exner[1]];

    // Shortwave: transmissivity of each layer along a path of `m` air masses.
    let transmit = |k: usize, m: f64| {
        let ozone = if k == 0 { params.ozone_absorption } else { 0.0 };
        let vapour = params.near_infrared
            * (1.0 - fms_math::exp(-m * params.water_absorption * water[k].0.sqrt()));
        fms_math::exp(-m * ozone) * (1.0 - vapour)
    };
    let down = [transmit(0, 1.0 / params.mu), transmit(1, 1.0 / params.mu)];
    let up = [transmit(0, 1.66), transmit(1, 1.66)];
    let s = column.insolation.max(0.0);
    let mut layers = [0.0; 2];
    // Through the upper layer to the cloud deck at the interface.
    layers[0] += s * (1.0 - down[0]);
    let at_cloud = s * down[0];
    let (ac, as_) = (params.cloud_albedo, params.albedo);
    // Below the cloud: the direct beam through the lower layer, then
    // reflections between surface and cloud base (diffuse transmissivity
    // `up[1]` each pass), summed as a geometric series.
    let below = (1.0 - ac) * at_cloud;
    let first = below * down[1];
    let surface_in = first / (1.0 - as_ * ac * up[1] * up[1]);
    let rising = as_ * surface_in * up[1];
    layers[1] +=
        below * (1.0 - down[1]) + as_ * surface_in * (1.0 - up[1]) + ac * rising * (1.0 - up[1]);
    let surface_shortwave = (1.0 - as_) * surface_in;
    // Upward through the upper layer: the cloud's reflection of the beam
    // from above plus what gets back through the cloud from below.
    let leaving = ac * at_cloud + (1.0 - ac) * rising;
    layers[0] += leaving * (1.0 - up[0]);
    let escaped = leaving * up[0];

    // Longwave, interval by interval.
    let width = co2_band_width(params, column);
    let ts = column.surface_temperature;
    let black = |t: f64| STEFAN_BOLTZMANN * t * t * t * t;
    let (mut outgoing, mut surface_down) = (0.0, 0.0);
    let mut surface_up = 0.0;
    // Edge temperatures: the lower layer follows the convective lapse rate
    // through its mean, up to the interface and down to the surface air; the
    // upper layer's top is at its own mean temperature (an isothermal,
    // stratosphere-like cap).
    let lower_profile = |p: f64| {
        let p_mean = p0 * fms_math::pow(exner[1], 1.0 / kappa);
        let beta = gas.gas_constant * params.lapse_rate / column.gravity;
        temperature[1] * fms_math::pow(p / p_mean, beta)
    };
    let interface = lower_profile(column.sigma * column.ps);
    let surface_air = lower_profile(column.ps);
    for (a, b, kind) in intervals(width) {
        let terms: &[(f64, f64)] = match kind {
            Kind::Water | Kind::WaterEdge => &WATER_TERMS,
            _ => &[(1.0, 1.0)],
        };
        // Planck emission of the interval at each layer's edges.
        let planck = |t: f64| black(t) * planck_fraction(a, b, t);
        let middle = planck(interface);
        // (upper edge, lower edge) of each layer.
        let edges = [
            (planck(temperature[0]), middle),
            (middle, planck(surface_air)),
        ];
        let u_planck = planck(ts);
        for &(share, scale) in terms {
            let tau: [f64; 2] = std::array::from_fn(|k| match kind {
                Kind::Co2 => params.co2_opacity * dsigma[k],
                Kind::Window => params.continuum_opacity * water[k].2,
                Kind::Water => scale * params.water_opacity * water[k].1,
                Kind::WaterEdge => scale * params.water_edge * params.water_opacity * water[k].1,
            });
            let trans = tau.map(|t| fms_math::exp(-t));
            // Emission up (from the layer top) and down (from its bottom) for
            // a source linear in optical depth.
            let emit: [(f64, f64); 2] = std::array::from_fn(|k| {
                let (b_top, b_bottom) = edges[k];
                let slope = linear_source_term(tau[k], trans[k]) * (b_bottom - b_top);
                (
                    share * ((1.0 - trans[k]) * b_top + slope),
                    share * ((1.0 - trans[k]) * b_bottom - slope),
                )
            });
            let u_surface = share * u_planck;
            // Upward: surface → lower layer → upper layer → space.
            let u_mid = u_surface * trans[1] + emit[1].0;
            let u_top = u_mid * trans[0] + emit[0].0;
            // Downward: space → upper → lower → surface.
            let d_mid = emit[0].1;
            let d_surface = d_mid * trans[1] + emit[1].1;
            layers[0] += (1.0 - trans[0]) * u_mid - emit[0].0 - emit[0].1;
            layers[1] += (1.0 - trans[1]) * (u_surface + d_mid) - emit[1].0 - emit[1].1;
            outgoing += u_top;
            surface_down += d_surface;
            surface_up += u_surface;
        }
    }
    Fluxes {
        absorbed: s - escaped,
        outgoing,
        layers,
        surface_shortwave,
        surface_down,
        surface_up,
        water: [water[0].0, water[1].0],
        co2_width: width,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn earthlike() -> Column {
        Column {
            ps: 1.0e5,
            surface_temperature: 288.0,
            theta: [330.0, 290.0],
            sigma: 0.5,
            gravity: 9.81,
            insolation: 340.0,
            co2: 400.0,
        }
    }

    #[test]
    fn planck_fractions_are_exact() {
        // The whole spectrum, and the classic Wien-peak split: half of the
        // emission lies above x = hcν/kT ≈ 3.503.
        assert!((planck_fraction(0.0, f64::INFINITY, 288.0) - 1.0).abs() < 1e-14);
        let half = 3.503_018 * 288.0 / C2;
        assert!((planck_fraction(half, f64::INFINITY, 288.0) - 0.5).abs() < 1e-5);
        // Both series agree where they meet.
        let (a, b) = (planck_tail(1.0 - 1e-12), planck_tail(1.0));
        assert!((a - b).abs() < 1e-12, "{a} vs {b}");
        // Fractions over a partition sum to 1.
        let sum: f64 = intervals(160.0)
            .iter()
            .map(|&(a, b, _)| planck_fraction(a, b, 250.0))
            .sum();
        assert!((sum - 1.0).abs() < 1e-13);
    }

    #[test]
    fn spectral_intervals_follow_the_band() {
        let kinds = |w| intervals(w).iter().map(|i| i.2).collect::<Vec<_>>();
        use Kind::*;
        assert_eq!(kinds(0.0), [Water, WaterEdge, Window, WaterEdge, Water]);
        assert_eq!(
            kinds(160.0),
            [Water, WaterEdge, Co2, WaterEdge, Window, WaterEdge, Water]
        );
        // A band wide enough to eat into the window.
        let wide = intervals(400.0);
        assert!(wide.iter().any(|i| i.2 == Kind::Co2 && i.1 > WINDOW.0));
    }

    #[test]
    fn the_budget_closes() {
        let params = RadiationParams::default();
        let f = column_fluxes(&params, Gas::EARTH_AIR, &earthlike());
        let net = f.absorbed - f.outgoing;
        let heating = f.layers[0] + f.layers[1] + f.surface();
        assert!((net - heating).abs() < 1e-10, "{net} vs {heating}");
    }

    #[test]
    fn an_isothermal_column_radiates_as_a_black_body() {
        // Everything at one temperature, no sunlight: the planet radiates as
        // a black body, and that is exactly what the column loses (the top
        // cools to space).
        let gas = Gas::EARTH_AIR;
        let t = 260.0;
        let mut c = earthlike();
        c.insolation = 0.0;
        c.surface_temperature = t;
        let ps: f64 = c.ps;
        let pi = 0.5 * ps;
        let kappa = gas.kappa();
        let exner = [
            gas.exner(pi) / (1.0 + kappa),
            (ps * gas.exner(ps) - pi * gas.exner(pi)) / ((1.0 + kappa) * 0.5 * ps),
        ];
        c.theta = [t / exner[0], t / exner[1]];
        // No lapse rate within the layers either.
        let params = RadiationParams {
            lapse_rate: 0.0,
            ..RadiationParams::default()
        };
        let f = column_fluxes(&params, gas, &c);
        let scale = STEFAN_BOLTZMANN * t.powi(4);
        assert!((f.outgoing - scale).abs() < 1e-9 * scale);
        let net = f.layers[0] + f.layers[1] + f.surface();
        assert!((net + f.outgoing).abs() < 1e-9 * scale);
        // The surface, under an atmosphere at its own temperature, only
        // loses what escapes through both layers.
        assert!(f.surface() < 0.0 && f.surface() > -0.5 * scale);
    }

    #[test]
    fn more_co2_traps_more_heat() {
        let params = RadiationParams::default();
        let olr = |co2| {
            let mut c = earthlike();
            c.co2 = co2;
            column_fluxes(&params, Gas::EARTH_AIR, &c).outgoing
        };
        let values: Vec<f64> = [0.0, 1.0, 100.0, 400.0, 800.0, 1e4, 1e5].map(olr).to_vec();
        assert!(values.windows(2).all(|w| w[1] <= w[0]), "{values:?}");
        // Logarithmic while the band is thick and clear of the window:
        // nearly equal steps per doubling.
        let (a, b, c) = (olr(200.0), olr(400.0), olr(800.0));
        assert!(((a - b) / (b - c) - 1.0).abs() < 0.1, "{a} {b} {c}");
    }
}
