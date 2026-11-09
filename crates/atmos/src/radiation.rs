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
//! depth between its edges. An interface's temperature is that of the layer
//! below it, along the convective lapse rate through that layer's mean
//! temperature; the lowest layer's bottom edge is the surface air, and the
//! top layer's top is at its own mean temperature (an isothermal,
//! stratosphere-like cap). So in an opaque band a layer emits upward from its
//! top and downward from its bottom: the lowest layer warms the surface from
//! its warm base, and neighbouring layers don't spuriously exchange heat at
//! their mean temperatures. The surface is black. Optical depths are diffuse (the 1.66
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
//! Two bands, each absorbed on the way down (slant path `1/μ`) and again
//! after reflection (diffuse path): the visible, absorbed by ozone at the
//! top, and the near-infrared (a fixed share of sunlight), absorbed by water
//! vapour. Water's absorption saturates with the slant path `u` from the
//! top (transmissivity `exp(−a√u)`), so each slab absorbs what's left of the
//! band between its top and bottom paths and the column's absorption
//! doesn't depend on how it's layered. A cloud deck at the layer interface reflects a fixed share (a
//! stand-in until the model diagnoses clouds), and the surface reflects its
//! own albedo; reflections between the two are summed exactly. Putting the
//! clouds' reflection at the surface instead, below the water vapour, would
//! make water's shortwave absorption feedback several times too strong.
//!
//! Every flux is in W m⁻² and every flux budget closes exactly: absorbed
//! sunlight minus outgoing longwave equals the net heating of the layers and
//! the surface.

use crate::layered::{Gas, REFERENCE_PRESSURE};
use crate::levels::Levels;

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
    /// Shortwave absorption by ozone: optical depth of the visible band, at
    /// the top.
    pub ozone_absorption: f64,
    /// Near-infrared absorption by water vapour: optical depth per
    /// `√(kg m⁻²)` of slant water path from the top.
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
    /// Albedo of the fixed cloud deck at mid-level that the single column
    /// (fixed relative humidity) uses in place of clouds.
    pub cloud_albedo: f64,
    /// Reflectance of the clear air (Rayleigh scattering), applied as a
    /// non-absorbing sheet at the atmosphere's mass centre (σ = ½).
    pub rayleigh: f64,
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
            ozone_absorption: 0.046,
            water_absorption: 0.062,
            near_infrared: 0.35,
            lapse_rate: 6.5e-3,
            mu: 0.5,
            albedo: 0.08,
            cloud_albedo: 0.33,
            rayleigh: 0.06,
        }
    }
}

/// One column's inputs.
#[derive(Clone, Copy, Debug)]
pub struct Column<'c> {
    /// Surface pressure (Pa).
    pub ps: f64,
    /// Surface temperature (K).
    pub surface_temperature: f64,
    /// Layer potential temperatures, top first (K).
    pub theta: &'c [f64],
    /// The layers.
    pub levels: &'c Levels,
    /// Gravity (m s⁻²).
    pub gravity: f64,
    /// Top-of-atmosphere insolation (W m⁻²).
    pub insolation: f64,
    /// CO₂ (ppm).
    pub co2: f64,
    /// Water vapour.
    pub humidity: Humidity<'c>,
    /// Cloud sheets (unused ones have zero reflectance and emissivity).
    pub clouds: &'c [Cloud],
}

/// Where a column's water vapour comes from.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Humidity<'c> {
    /// Fixed relative humidity: the Manabe–Wetherald profile
    /// `RH(σ) = RH_s (σ − 0.02)/0.98` with
    /// [`RadiationParams::surface_humidity`].
    Fixed,
    /// Each layer's mean specific humidity (kg kg⁻¹), distributed within the
    /// layer with the same profile shape.
    Layers(&'c [f64]),
}

/// A cloud sheet: an infinitely thin, grey layer of partial cloud at one
/// level, reflecting sunlight (without absorbing it) and absorbing and
/// emitting longwave at the local air temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cloud {
    /// Level σ.
    pub sigma: f64,
    /// Share of sunlight reflected (cloud fraction × cloud albedo).
    pub reflectance: f64,
    /// Longwave emissivity (cloud fraction × cloud emissivity).
    pub emissivity: f64,
}

impl Cloud {
    /// No cloud.
    pub const NONE: Self = Self {
        sigma: 0.5,
        reflectance: 0.0,
        emissivity: 0.0,
    };

    /// A deck at `sigma` reflecting `albedo` of sunlight, transparent to
    /// longwave (the fixed stand-in for clouds of M2 step 8).
    pub fn deck(sigma: f64, albedo: f64) -> Self {
        Self {
            sigma,
            reflectance: albedo,
            emissivity: 0.0,
        }
    }

    fn is_empty(&self) -> bool {
        self.reflectance == 0.0 && self.emissivity == 0.0
    }
}

/// One column's fluxes (W m⁻²) and water (kg m⁻²).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Fluxes {
    /// Absorbed sunlight (insolation minus reflected).
    pub absorbed: f64,
    /// Outgoing longwave.
    pub outgoing: f64,
    /// Net radiative heating of each layer (shortwave + longwave).
    pub layers: Vec<f64>,
    /// Shortwave absorbed by the surface.
    pub surface_shortwave: f64,
    /// Downward longwave at the surface.
    pub surface_down: f64,
    /// Upward longwave from the surface.
    pub surface_up: f64,
    /// Water vapour in each layer (kg m⁻²).
    pub water: Vec<f64>,
    /// Width of the opaque CO₂ band (cm⁻¹).
    pub co2_width: f64,
}

impl Fluxes {
    /// Net radiative heating of the surface.
    pub fn surface(&self) -> f64 {
        self.surface_shortwave + self.surface_down - self.surface_up
    }

    /// Water vapour in the whole column (kg m⁻²).
    pub fn total_water(&self) -> f64 {
        self.water.iter().sum()
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
        // e^{−nx}, by repeated multiplication.
        let step = fms_math::exp(-x);
        let mut decay = 1.0;
        (1..=terms)
            .map(|n| {
                let n = n as f64;
                decay *= step;
                decay
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

/// Shape of relative humidity in σ (Manabe & Wetherald 1967): moistest at
/// the surface, dry at the top.
fn humidity_shape(sigma: f64) -> f64 {
    ((sigma - 0.02) / 0.98).max(0.0)
}

/// Constants of the in-layer profile for one set of levels, gas and
/// exponent `β = RΓ/g`: `(σ/σ̄ₖ)^β` where the model always needs it. `σ̄ₖ` is
/// the level whose Exner function is layer `k`'s mean, as a fraction of the
/// surface pressure (the same for every column).
#[derive(Clone, Debug)]
pub(crate) struct ProfileShape {
    kappa: f64,
    beta: f64,
    mean_sigma: Vec<f64>,
    /// At each layer's Gauss nodes.
    nodes: Vec<[f64; 4]>,
    /// At each interface, for the layer above it and the one below.
    interfaces: Vec<[f64; 2]>,
}

impl ProfileShape {
    fn new(levels: &Levels, gas: Gas, beta: f64) -> Self {
        let kappa = gas.kappa();
        let mean_sigma: Vec<f64> = levels
            .layer_exner(gas, REFERENCE_PRESSURE)
            .into_iter()
            .map(|e| fms_math::pow(e, 1.0 / kappa))
            .collect();
        let factor = |k: usize, s: f64| fms_math::pow(s / mean_sigma[k], beta);
        let nodes = (0..levels.count())
            .map(|k| {
                let (a, b) = levels.bounds(k);
                GAUSS.map(|(x, _)| factor(k, 0.5 * (a + b) + 0.5 * (b - a) * x))
            })
            .collect();
        let interfaces = levels
            .interior()
            .iter()
            .enumerate()
            .map(|(j, &s)| [factor(j, s), factor(j + 1, s)])
            .collect();
        Self {
            kappa,
            beta,
            mean_sigma,
            nodes,
            interfaces,
        }
    }

    /// The levels' memoized constants if they are for this gas and `β`,
    /// otherwise the same values computed afresh.
    fn of(levels: &Levels, gas: Gas, beta: f64) -> std::borrow::Cow<'_, Self> {
        let shape = levels.profile_shape(|| Self::new(levels, gas, beta));
        if shape.kappa == gas.kappa() && shape.beta == beta {
            std::borrow::Cow::Borrowed(shape)
        } else {
            std::borrow::Cow::Owned(Self::new(levels, gas, beta))
        }
    }
}

/// The temperature profile within the layers: each follows the lapse rate
/// `Γ` through its mean temperature `Tₖ` at `pₖ` (the pressure whose Exner
/// function is the layer mean): `T(p) = Tₖ (p/pₖ)^{RΓ/g}`.
#[derive(Clone, Debug)]
pub(crate) struct Profile<'l> {
    temperature: &'l [f64],
    shape: std::borrow::Cow<'l, ProfileShape>,
    beta: f64,
    ps: f64,
    levels: &'l Levels,
}

impl<'l> Profile<'l> {
    pub(crate) fn new(
        params: &RadiationParams,
        gas: Gas,
        levels: &'l Levels,
        ps: f64,
        gravity: f64,
        t: &'l [f64],
    ) -> Self {
        let beta = gas.gas_constant * params.lapse_rate / gravity;
        Self {
            temperature: t,
            shape: ProfileShape::of(levels, gas, beta),
            beta,
            ps,
            levels,
        }
    }

    pub(crate) fn at(&self, s: f64) -> f64 {
        self.in_layer(self.levels.layer_of(s), s)
    }

    /// Layer `k`'s profile at `s` (extended beyond the layer if need be).
    pub(crate) fn in_layer(&self, k: usize, s: f64) -> f64 {
        self.temperature[k] * fms_math::pow(s / self.shape.mean_sigma[k], self.beta)
    }

    /// Layer `k`'s profile at its Gauss node `g`.
    fn at_node(&self, k: usize, g: usize) -> f64 {
        self.temperature[k] * self.shape.nodes[k][g]
    }

    /// At interface `j`, the profile of the layer above it (`side` 0) or
    /// below it (1).
    fn at_interface(&self, j: usize, side: usize) -> f64 {
        self.temperature[j + side] * self.shape.interfaces[j][side]
    }

    /// Temperature at interface `j` (below layer `j`): the layer below's
    /// profile there.
    fn interface(&self, j: usize) -> f64 {
        self.at_interface(j, 1)
    }
}

/// Saturation specific humidity (kg kg⁻¹) at temperature `t` and pressure
/// `p`.
pub fn saturation_humidity(t: f64, p: f64) -> f64 {
    specific_humidity(saturation_pressure(t), p)
}

/// Each layer's saturation, for the moist physics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerSaturation {
    /// Mass-weighted mean saturation specific humidity `q̄*` (kg kg⁻¹).
    pub mean: f64,
    /// Mass-weighted mean of the humidity profile shape times saturation:
    /// a layer holding `qₖ` has relative humidity `(qₖ / shaped) ·
    /// shape(σ)` at each level.
    pub shaped: f64,
}

/// Saturation of each layer at layer-mean temperatures `t` (K, top first),
/// along the in-layer temperature profile (by Gauss–Legendre quadrature).
pub fn layer_saturation(
    params: &RadiationParams,
    gas: Gas,
    levels: &Levels,
    ps: f64,
    gravity: f64,
    t: &[f64],
) -> Vec<LayerSaturation> {
    let profile = Profile::new(params, gas, levels, ps, gravity, t);
    saturation_of(&profile)
}

fn saturation_of(profile: &Profile) -> Vec<LayerSaturation> {
    let ps = profile.ps;
    (0..profile.levels.count())
        .map(|k| {
            let (a, b) = profile.levels.bounds(k);
            let (mut mean, mut shaped) = (0.0, 0.0);
            for (g, (x, w)) in GAUSS.into_iter().enumerate() {
                let s = 0.5 * (a + b) + 0.5 * (b - a) * x;
                let q = saturation_humidity(profile.at_node(k, g), s * ps);
                mean += 0.5 * w * q;
                shaped += 0.5 * w * humidity_shape(s) * q;
            }
            LayerSaturation { mean, shaped }
        })
        .collect()
}

/// Per interface, specific humidity there relative to each adjacent layer's
/// mean, along each layer's own profile: `[upper layer at its bottom, lower
/// layer at its top]`. Water is concentrated low in each layer, so the lower
/// layer's air at the interface holds far less than its mean, and the upper
/// layer's far more. The dynamics carries humidity across the interfaces at
/// these values.
pub fn interface_humidity_ratios(
    params: &RadiationParams,
    gas: Gas,
    levels: &Levels,
    ps: f64,
    gravity: f64,
    t: &[f64],
) -> Vec<[f64; 2]> {
    let profile = Profile::new(params, gas, levels, ps, gravity, t);
    let sat = saturation_of(&profile);
    levels
        .interior()
        .iter()
        .enumerate()
        .map(|(j, &sigma)| {
            std::array::from_fn(|side| {
                let k = j + side;
                let at = humidity_shape(sigma)
                    * saturation_humidity(profile.at_interface(j, side), sigma * ps);
                if sat[k].shaped > 0.0 {
                    at / sat[k].shaped
                } else {
                    1.0
                }
            })
        })
        .collect()
}

/// A slab of gas between two levels, within one layer.
#[derive(Clone, Copy, Debug, Default)]
struct Slab {
    layer: usize,
    top: f64,
    bottom: f64,
    /// Water (kg m⁻²), pressure-weighted water and continuum absorber.
    water: f64,
    weighted: f64,
    continuum: f64,
}

/// One element of the column, top to bottom.
#[derive(Clone, Copy, Debug)]
enum Element {
    Slab(usize),
    Sheet(usize),
}

/// The column cut into slabs at the interfaces and the cloud levels, with
/// each slab's water: `q(p) = Aₖ · shape(σ) · q*(T(p), p)`, where `Aₖ` is
/// fixed ([`Humidity::Fixed`]) or set so the layer holds its humidity.
fn slabs(
    params: &RadiationParams,
    column: &Column,
    sheets: &[Cloud],
    profile: &Profile,
) -> (Vec<Slab>, Vec<Element>) {
    let levels = column.levels;
    let n = levels.count();
    let mut cuts = levels.interfaces().to_vec();
    for c in sheets {
        if !c.is_empty() && c.sigma > 0.0 && c.sigma < 1.0 {
            cuts.push(c.sigma);
        }
    }
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    let mut slabs: Vec<Slab> = cuts
        .windows(2)
        .map(|w| Slab {
            layer: levels.layer_of(0.5 * (w[0] + w[1])),
            top: w[0],
            bottom: w[1],
            ..Slab::default()
        })
        .collect();
    // Nodes: (slab, mass, σ, shape × q*).
    let p0 = REFERENCE_PRESSURE;
    let mut nodes = Vec::with_capacity(4 * slabs.len());
    let mut shaped = vec![0.0; n];
    for (j, slab) in slabs.iter().enumerate() {
        let (a, b) = (slab.top, slab.bottom);
        for (x, w) in GAUSS {
            let s = 0.5 * (a + b) + 0.5 * (b - a) * x;
            let mass = 0.5 * w * (b - a) * column.ps / column.gravity;
            let q = humidity_shape(s) * saturation_humidity(profile.at(s), s * column.ps);
            shaped[slab.layer] += mass * q;
            nodes.push((j, mass, s, q));
        }
    }
    let ds = levels.dsigma();
    let amplitude: Vec<f64> = match column.humidity {
        Humidity::Fixed => vec![params.surface_humidity; n],
        Humidity::Layers(q) => (0..n)
            .map(|k| {
                let water = q[k].max(0.0) * ds[k] * column.ps / column.gravity;
                if shaped[k] > 0.0 {
                    water / shaped[k]
                } else {
                    0.0
                }
            })
            .collect(),
    };
    for (j, mass, s, q) in nodes {
        let slab = &mut slabs[j];
        let q = amplitude[slab.layer] * q;
        let p = s * column.ps;
        let e = (q * p / (EPSILON + (1.0 - EPSILON) * q)).min(0.25 * p);
        slab.water += mass * q;
        slab.weighted += mass * q * p / p0;
        slab.continuum += mass * q * e / 1000.0;
    }
    // Elements: slabs, with each cloud sheet at its level.
    let mut elements = Vec::with_capacity(slabs.len() + sheets.len());
    let mut order: Vec<usize> = (0..sheets.len())
        .filter(|&c| !sheets[c].is_empty())
        .collect();
    order.sort_by(|&a, &b| sheets[a].sigma.total_cmp(&sheets[b].sigma));
    let mut next = order.into_iter().peekable();
    for (j, slab) in slabs.iter().enumerate() {
        while let Some(&c) = next.peek() {
            if sheets[c].sigma <= slab.top {
                elements.push(Element::Sheet(c));
                next.next();
            } else {
                break;
            }
        }
        elements.push(Element::Slab(j));
    }
    elements.extend(next.map(Element::Sheet));
    (slabs, elements)
}

/// Diffusivity factor: the mean slant path of diffuse light (air masses).
const DIFFUSIVITY: f64 = 1.66;

/// Sunlight in one band through the column's elements.
struct Band {
    /// Absorbed by each element (zero for sheets).
    absorbed: Vec<f64>,
    /// Absorbed by the surface.
    surface: f64,
    /// Reflected to space.
    escaped: f64,
}

/// The adding method for one band of sunlight `s` (W m⁻²): sheets reflect
/// and transmit without absorbing; `slab(j)` gives slab `j`'s direct and
/// diffuse transmissivities; reflections between everything are summed
/// exactly.
fn shortwave_band(
    elements: &[Element],
    sheets: &[Cloud],
    albedo: f64,
    s: f64,
    slab: impl Fn(usize) -> (f64, f64),
) -> Band {
    // Per element: reflectance and direct and diffuse transmissivities.
    let optics: Vec<(f64, f64, f64)> = elements
        .iter()
        .map(|e| match *e {
            Element::Slab(j) => {
                let (direct, diffuse) = slab(j);
                (0.0, direct, diffuse)
            }
            Element::Sheet(c) => {
                let r = sheets[c].reflectance;
                (r, 1.0 - r, 1.0 - r)
            }
        })
        .collect();
    let n = elements.len();
    // Direct beam down; diffuse reflectance R̃ of everything below each
    // interface, and the diffuse light U leaving each interface upward that
    // the direct beam generates below it.
    let mut direct = vec![0.0; n + 1];
    direct[0] = s;
    for (j, &(_, t_dir, _)) in optics.iter().enumerate() {
        direct[j + 1] = direct[j] * t_dir;
    }
    let mut reflect = vec![0.0; n + 1];
    let mut rising = vec![0.0; n + 1];
    reflect[n] = albedo;
    rising[n] = albedo * direct[n];
    for j in (0..n).rev() {
        let (r, _, t) = optics[j];
        let denom = 1.0 - r * reflect[j + 1];
        reflect[j] = r + t * t * reflect[j + 1] / denom;
        rising[j] = t * rising[j + 1] / denom + r * direct[j];
    }
    let mut down = vec![0.0; n + 1];
    for j in 0..n {
        let (r, _, t) = optics[j];
        down[j + 1] = (t * down[j] + r * rising[j + 1]) / (1.0 - r * reflect[j + 1]);
    }
    let up = |j: usize| reflect[j] * down[j] + rising[j];
    let absorbed = (0..n)
        .map(|j| {
            let (r, _, t) = optics[j];
            if r > 0.0 {
                0.0
            } else {
                (direct[j] - direct[j + 1]) + (1.0 - t) * (down[j] + up(j + 1))
            }
        })
        .collect();
    Band {
        absorbed,
        surface: (1.0 - albedo) * (down[n] + direct[n]),
        escaped: up(0),
    }
}

/// Radiative fluxes of one column.
pub fn column_fluxes(params: &RadiationParams, gas: Gas, column: &Column) -> Fluxes {
    let levels = column.levels;
    let n_layers = levels.count();
    let exner = levels.layer_exner(gas, column.ps);
    let temperature: Vec<f64> = (0..n_layers).map(|k| column.theta[k] * exner[k]).collect();
    let profile = Profile::new(params, gas, levels, column.ps, column.gravity, &temperature);
    // The clouds, then the clear air's Rayleigh scattering.
    let mut sheets = column.clouds.to_vec();
    sheets.push(Cloud {
        sigma: 0.5,
        reflectance: params.rayleigh,
        emissivity: 0.0,
    });
    let (slabs, elements) = slabs(params, column, &sheets, &profile);
    let mut water = vec![0.0; n_layers];
    for slab in &slabs {
        water[slab.layer] += slab.water;
    }
    let sheet_layer = |c: usize| levels.layer_of(sheets[c].sigma);
    let mut layers = vec![0.0; n_layers];

    // Shortwave, in two bands: the visible (ozone at the top) and the
    // near-infrared (water vapour). Each band runs the adding method
    // through the same sheets and surface. Water's absorption saturates
    // with the path from the top: transmissivity `exp(−a√(mW))` along `m`
    // air masses through `W` of water, so a slab between cumulative paths
    // `W₁` and `W₂` transmits `exp(−a(√(mW₂) − √(mW₁)))` and the column's
    // absorption doesn't depend on how it's cut into slabs.
    let mut path = vec![0.0; slabs.len() + 1];
    for (j, slab) in slabs.iter().enumerate() {
        path[j + 1] = path[j] + slab.water;
    }
    let s = column.insolation.max(0.0);
    let albedo = params.albedo;
    let visible = |j: usize, m: f64| {
        if j == 0 {
            fms_math::exp(-m * params.ozone_absorption)
        } else {
            1.0
        }
    };
    let infrared = |j: usize, m: f64| {
        let a = params.water_absorption;
        fms_math::exp(-a * ((m * path[j + 1]).sqrt() - (m * path[j]).sqrt()))
    };
    let mut surface_shortwave = 0.0;
    let mut escaped = 0.0;
    for (share, transmit) in [
        (
            1.0 - params.near_infrared,
            &visible as &dyn Fn(usize, f64) -> f64,
        ),
        (params.near_infrared, &infrared),
    ] {
        let band = shortwave_band(&elements, &sheets, albedo, share * s, |j| {
            (transmit(j, 1.0 / params.mu), transmit(j, DIFFUSIVITY))
        });
        for (j, e) in elements.iter().enumerate() {
            if let Element::Slab(i) = *e {
                layers[slabs[i].layer] += band.absorbed[j];
            }
        }
        surface_shortwave += band.surface;
        escaped += band.escaped;
    }

    // Longwave, interval by interval. Edge temperatures: each interface
    // that of the layer below it, along the convective lapse rate through
    // its mean; the lowest layer's bottom the surface air; the top layer's
    // top its own mean temperature (an isothermal, stratosphere-like cap).
    let width = co2_band_width(params, column);
    let ts = column.surface_temperature;
    let black = |t: f64| STEFAN_BOLTZMANN * t * t * t * t;
    let (mut outgoing, mut surface_down) = (0.0, 0.0);
    let mut surface_up = 0.0;
    let interface: Vec<f64> = (0..n_layers - 1).map(|j| profile.interface(j)).collect();
    let surface_air = profile.in_layer(n_layers - 1, 1.0);
    let cloud_temperature: Vec<f64> = sheets.iter().map(|c| profile.at(c.sigma)).collect();
    let ds = levels.dsigma();
    let n = elements.len();
    let mut edges = vec![(0.0, 0.0); n_layers];
    let mut tau = vec![0.0; slabs.len()];
    let mut clear = vec![0.0; slabs.len()];
    let mut sheet_planck = vec![0.0; cloud_temperature.len()];
    let mut layer_tau = vec![0.0; n_layers];
    let mut above = vec![0.0; n_layers];
    let mut emit = vec![(0.0, 0.0); n];
    let mut rise = vec![0.0; n + 1];
    let mut fall = vec![0.0; n + 1];
    for (a, b, kind) in intervals(width) {
        let terms: &[(f64, f64)] = match kind {
            Kind::Water | Kind::WaterEdge => &WATER_TERMS,
            _ => &[(1.0, 1.0)],
        };
        // Planck emission of the interval at each layer's edges.
        let planck = |t: f64| black(t) * planck_fraction(a, b, t);
        let mut top = planck(temperature[0]);
        for (k, edge) in edges.iter_mut().enumerate() {
            let bottom = if k + 1 < n_layers {
                planck(interface[k])
            } else {
                planck(surface_air)
            };
            *edge = (top, bottom);
            top = bottom;
        }
        for (s, &t) in sheet_planck.iter_mut().zip(&cloud_temperature) {
            *s = planck(t);
        }
        let u_planck = planck(ts);
        for &(share, scale) in terms {
            layer_tau.fill(0.0);
            for (t, slab) in tau.iter_mut().zip(&slabs) {
                *t = match kind {
                    Kind::Co2 => params.co2_opacity * (slab.bottom - slab.top),
                    Kind::Window => params.continuum_opacity * slab.continuum,
                    Kind::Water => scale * params.water_opacity * slab.weighted,
                    Kind::WaterEdge => {
                        scale * params.water_edge * params.water_opacity * slab.weighted
                    }
                };
                layer_tau[slab.layer] += *t;
            }
            for (c, t) in clear.iter_mut().zip(&tau) {
                *c = fms_math::exp(-t);
            }
            // Emission of each element up (from its top) and down (from its
            // bottom). A slab's source is linear in optical depth across its
            // layer, so its edges' values are interpolated in τ.
            above.fill(0.0);
            for (j, e) in elements.iter().enumerate() {
                emit[j] = match *e {
                    Element::Slab(i) => {
                        let slab = &slabs[i];
                        let k = slab.layer;
                        let (b_top, b_bottom) = edges[k];
                        // Position of the slab's edges within its layer, in
                        // optical depth (or σ, if the layer is transparent).
                        let (x0, x1) = if layer_tau[k] > 1e-12 {
                            (above[k] / layer_tau[k], (above[k] + tau[i]) / layer_tau[k])
                        } else {
                            let top = levels.bounds(k).0;
                            ((slab.top - top) / ds[k], (slab.bottom - top) / ds[k])
                        };
                        above[k] += tau[i];
                        let (top, bottom) = (
                            b_top + (b_bottom - b_top) * x0,
                            b_top + (b_bottom - b_top) * x1,
                        );
                        let t = clear[i];
                        let slope = linear_source_term(tau[i], t) * (bottom - top);
                        (
                            share * ((1.0 - t) * top + slope),
                            share * ((1.0 - t) * bottom - slope),
                        )
                    }
                    Element::Sheet(c) => {
                        let e = sheets[c].emissivity * share * sheet_planck[c];
                        (e, e)
                    }
                };
            }
            let transmit = |j: usize| match elements[j] {
                Element::Slab(i) => clear[i],
                Element::Sheet(c) => 1.0 - sheets[c].emissivity,
            };
            // Upward: surface → … → space; downward: space → … → surface.
            let u_surface = share * u_planck;
            rise[n] = u_surface;
            for j in (0..n).rev() {
                rise[j] = rise[j + 1] * transmit(j) + emit[j].0;
            }
            fall[0] = 0.0;
            for j in 0..n {
                fall[j + 1] = fall[j] * transmit(j) + emit[j].1;
            }
            for (j, e) in elements.iter().enumerate() {
                let k = match *e {
                    Element::Slab(i) => slabs[i].layer,
                    Element::Sheet(c) => sheet_layer(c),
                };
                layers[k] += (1.0 - transmit(j)) * (rise[j + 1] + fall[j]) - emit[j].0 - emit[j].1;
            }
            outgoing += rise[0];
            surface_down += fall[n];
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
        water,
        co2_width: width,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const THETA: [f64; 2] = [330.0, 290.0];
    const DECK: [Cloud; 2] = [
        Cloud {
            sigma: 0.5,
            reflectance: 0.33,
            emissivity: 0.0,
        },
        Cloud::NONE,
    ];

    fn two() -> &'static Levels {
        static L: std::sync::OnceLock<Levels> = std::sync::OnceLock::new();
        L.get_or_init(|| Levels::two(0.5))
    }

    fn earthlike() -> Column<'static> {
        Column {
            ps: 1.0e5,
            surface_temperature: 288.0,
            theta: &THETA,
            levels: two(),
            gravity: 9.81,
            insolation: 340.0,
            co2: 400.0,
            humidity: Humidity::Fixed,
            clouds: &DECK,
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
    fn more_layers_resolve_the_same_column() {
        // The same atmosphere (each layer at the two-layer profile's
        // temperature there, fixed relative humidity) cut into more
        // layers: the budget closes, the water is the same, and the
        // longwave and shortwave agree to within the resampled profile.
        // Measured against two layers: OLR +2.3 and +4.9 W/m², surface
        // longwave +4.0 and +4.5, water +0.5% and +0.3%, surface sunlight
        // −0.03 and −0.01 (the
        // atmosphere's +0.6 and +0.4, mostly the resampled water; O-27).
        let params = RadiationParams::default();
        let gas = Gas::EARTH_AIR;
        let base = earthlike();
        let reference = column_fluxes(&params, gas, &base);
        let t2: Vec<f64> = two()
            .layer_exner(gas, base.ps)
            .iter()
            .zip(THETA)
            .map(|(e, th)| e * th)
            .collect();
        let profile = Profile::new(&params, gas, two(), base.ps, base.gravity, &t2);
        for levels in [
            Levels::new(&[0.2, 0.5, 0.8]),
            Levels::new(&[0.1, 0.3, 0.5, 0.7, 0.85]),
        ] {
            // Layer means of the two-layer profile (an isothermal cap above
            // the upper layer's mean level).
            let p_mean = two().mean_pressures(gas, base.ps);
            let exner = levels.layer_exner(gas, base.ps);
            let theta: Vec<f64> = (0..levels.count())
                .map(|k| {
                    let (a, b) = levels.bounds(k);
                    let t: f64 = GAUSS
                        .iter()
                        .map(|&(x, w)| {
                            let s = 0.5 * (a + b) + 0.5 * (b - a) * x;
                            let t = if s * base.ps < p_mean[0] {
                                t2[0]
                            } else {
                                profile.at(s)
                            };
                            0.5 * w * t
                        })
                        .sum();
                    t / exner[k]
                })
                .collect();
            let column = Column {
                theta: &theta,
                levels: &levels,
                ..base
            };
            let f = column_fluxes(&params, gas, &column);
            let heating: f64 = f.layers.iter().sum::<f64>() + f.surface();
            assert!((f.absorbed - f.outgoing - heating).abs() < 1e-10);
            assert!((f.total_water() / reference.total_water() - 1.0).abs() < 0.01);
            assert!((f.outgoing - reference.outgoing).abs() < 6.0, "{f:?}");
            assert!(
                (f.surface_down - reference.surface_down).abs() < 6.0,
                "{f:?}"
            );
            assert!(
                (f.surface_shortwave - reference.surface_shortwave).abs() < 0.2,
                "{f:?}"
            );
            assert!((f.absorbed - reference.absorbed).abs() < 1.5, "{f:?}");
        }
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
        let exner = two().layer_exner(gas, c.ps);
        let theta = [t / exner[0], t / exner[1]];
        c.theta = &theta;
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

    #[test]
    fn an_invisible_sheet_changes_nothing() {
        // Cutting the layers at extra levels changes only the quadrature of
        // the water profile: the slabs' sources are linear in optical depth
        // across each layer either way.
        let params = RadiationParams::default();
        let mut c = earthlike();
        c.humidity = Humidity::Layers(&[3e-4, 6e-3]);
        c.clouds = &[Cloud::NONE; 2];
        let plain = column_fluxes(&params, Gas::EARTH_AIR, &c);
        let faint = |sigma| Cloud {
            sigma,
            reflectance: 1e-300,
            emissivity: 1e-300,
        };
        let sheets = [faint(0.25), faint(0.8)];
        c.clouds = &sheets;
        let cut = column_fluxes(&params, Gas::EARTH_AIR, &c);
        for (a, b) in [
            (plain.outgoing, cut.outgoing),
            (plain.absorbed, cut.absorbed),
            (plain.surface_down, cut.surface_down),
            (plain.layers[0], cut.layers[0]),
            (plain.layers[1], cut.layers[1]),
            (plain.water[1], cut.water[1]),
        ] {
            assert!((a - b).abs() < 0.05 + 1e-3 * a.abs(), "{a} vs {b}");
        }
        // The layer holds exactly the water it was given.
        assert!((cut.water[1] - 6e-3 * 0.5 * c.ps / c.gravity).abs() < 1e-9);
    }

    #[test]
    fn clouds_reflect_trap_and_keep_the_budget() {
        let params = RadiationParams::default();
        let mut c = earthlike();
        c.humidity = Humidity::Layers(&[3e-4, 6e-3]);
        c.clouds = &[Cloud::NONE; 2];
        let clear = column_fluxes(&params, Gas::EARTH_AIR, &c);
        let sheets = [
            Cloud {
                sigma: 0.25,
                reflectance: 0.2,
                emissivity: 0.6,
            },
            Cloud {
                sigma: 0.8,
                reflectance: 0.3,
                emissivity: 0.7,
            },
        ];
        c.clouds = &sheets;
        let cloudy = column_fluxes(&params, Gas::EARTH_AIR, &c);
        assert!(cloudy.absorbed < clear.absorbed - 50.0);
        assert!(cloudy.outgoing < clear.outgoing - 10.0);
        // Low cloud radiates down onto the surface.
        assert!(cloudy.surface_down > clear.surface_down);
        let net = cloudy.absorbed - cloudy.outgoing;
        let heating = cloudy.layers[0] + cloudy.layers[1] + cloudy.surface();
        assert!((net - heating).abs() < 1e-10, "{net} vs {heating}");
    }
}
