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
use crate::levels::{Levels, MAX_LAYERS};

/// Stefan–Boltzmann constant (W m⁻² K⁻⁴).
pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;
/// Second radiation constant `hc/k` (cm K).
const C2: f64 = 1.438_776_877;
/// Index of the visible band of sunlight.
pub const VISIBLE: usize = 0;
/// Index of the near-infrared band of sunlight.
pub const NEAR_INFRARED: usize = 1;

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
    /// Surface albedo in each band of sunlight ([`VISIBLE`],
    /// [`NEAR_INFRARED`]).
    pub albedo: [f64; 2],
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
            albedo: [0.08; 2],
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
    /// The visible band's part of `surface_shortwave`.
    pub surface_visible: f64,
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

    /// Values in a restart row of a column with `layers` layers.
    pub(crate) fn restart_width(layers: usize) -> usize {
        7 + 2 * layers
    }

    /// Append the fluxes to `row` as a restart row, in declaration order
    /// with each per-layer field written in full.
    ///
    /// # Panics
    /// If a per-layer field doesn't have `layers` values.
    pub(crate) fn push_restart_row(&self, layers: usize, row: &mut Vec<f64>) {
        let Self {
            absorbed,
            outgoing,
            layers: heating,
            surface_shortwave,
            surface_visible,
            surface_down,
            surface_up,
            water,
            co2_width,
        } = self;
        assert!(heating.len() == layers && water.len() == layers);
        row.extend([*absorbed, *outgoing]);
        row.extend(heating);
        row.extend([
            *surface_shortwave,
            *surface_visible,
            *surface_down,
            *surface_up,
        ]);
        row.extend(water);
        row.push(*co2_width);
    }

    /// The fluxes of a row written by [`push_restart_row`](Self::push_restart_row),
    /// which has [`restart_width`](Self::restart_width) values.
    pub(crate) fn from_restart_row(row: &[f64], layers: usize) -> Self {
        let (head, rest) = row.split_at(2);
        let (heating, rest) = rest.split_at(layers);
        let (surface, rest) = rest.split_at(4);
        let (water, rest) = rest.split_at(layers);
        Self {
            absorbed: head[0],
            outgoing: head[1],
            layers: heating.to_vec(),
            surface_shortwave: surface[0],
            surface_visible: surface[1],
            surface_down: surface[2],
            surface_up: surface[3],
            water: water.to_vec(),
            co2_width: rest[0],
        }
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
    let total = planck_tail(0.0);
    ((planck_tail_at(a, t) - planck_tail_at(b, t)) / total).max(0.0)
}

/// [`planck_tail`] above wavenumber `nu` (cm⁻¹) at temperature `t`.
fn planck_tail_at(nu: f64, t: f64) -> f64 {
    if nu.is_infinite() {
        0.0
    } else {
        planck_tail(C2 * nu / t)
    }
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

/// Most cloud sheets a column may have ([`Column::clouds`]).
pub const MAX_CLOUDS: usize = 4;
/// The clouds and the clear air's Rayleigh sheet.
const MAX_SHEETS: usize = MAX_CLOUDS + 1;
/// The layers, cut at each sheet.
const MAX_SLABS: usize = MAX_LAYERS + MAX_SHEETS;
/// The slabs and the sheets between them.
const MAX_ELEMENTS: usize = MAX_SLABS + MAX_SHEETS;
/// Spectral intervals: the pieces between at most eight breakpoints.
const MAX_INTERVALS: usize = 7;

/// A list of at most `N` values, kept on the stack.
#[derive(Clone, Copy, Debug)]
struct Fixed<T, const N: usize> {
    items: [T; N],
    len: usize,
}

impl<T: Copy, const N: usize> Fixed<T, N> {
    /// Empty, with every slot holding `value`.
    fn filled(value: T) -> Self {
        Self {
            items: [value; N],
            len: 0,
        }
    }

    /// Append `value`.
    ///
    /// # Panics
    /// If the list is full.
    fn push(&mut self, value: T) {
        self.items[self.len] = value;
        self.len += 1;
    }

    fn clear(&mut self) {
        self.len = 0;
    }
}

impl<T, const N: usize> std::ops::Deref for Fixed<T, N> {
    type Target = [T];

    fn deref(&self) -> &[T] {
        &self.items[..self.len]
    }
}

impl<T, const N: usize> std::ops::DerefMut for Fixed<T, N> {
    fn deref_mut(&mut self) -> &mut [T] {
        &mut self.items[..self.len]
    }
}

/// The spectral intervals for a CO₂ band of width `w`.
#[cfg(test)]
fn intervals(w: f64) -> Fixed<(f64, f64, Kind), MAX_INTERVALS> {
    spectrum(w).intervals
}

/// The spectral intervals for a CO₂ band of width `w`, and their
/// breakpoints: each interval lies between two of them.
#[derive(Clone, Copy, Debug)]
struct Spectrum {
    /// Sorted, without repeats.
    cuts: Fixed<f64, 8>,
    intervals: Fixed<(f64, f64, Kind), MAX_INTERVALS>,
    /// The cuts each interval lies between.
    bounds: Fixed<(usize, usize), MAX_INTERVALS>,
}

impl Spectrum {
    /// `black(t) · planck_fraction(a, b, t)` of each interval `(a, b)` at
    /// temperature `t`, into `out`; each breakpoint's tail is computed once
    /// for the two intervals it bounds.
    fn planck(&self, t: f64, out: &mut [f64]) {
        let mut tail = [0.0; 8];
        for (g, &nu) in tail.iter_mut().zip(self.cuts.iter()) {
            *g = planck_tail_at(nu, t);
        }
        let total = planck_tail(0.0);
        let emission = black(t);
        for (o, &(a, b)) in out.iter_mut().zip(self.bounds.iter()) {
            *o = emission * ((tail[a] - tail[b]) / total).max(0.0);
        }
    }
}

/// [`spectrum`] of each width, for the intervals alone.
fn spectrum(w: f64) -> Spectrum {
    let (lo, hi) = ((CO2_CENTRE - 0.5 * w).max(0.0), CO2_CENTRE + 0.5 * w);
    // Breakpoints, then classify each piece by its midpoint.
    let mut cuts = Fixed::<f64, 8>::filled(0.0);
    for cut in [
        0.0,
        WATER_CORES.0,
        WINDOW.0,
        WINDOW.1,
        WATER_CORES.1,
        f64::INFINITY,
    ] {
        cuts.push(cut);
    }
    if w > 0.0 {
        cuts.push(lo);
        cuts.push(hi);
    }
    cuts.sort_by(f64::total_cmp);
    let mut unique = Fixed::<f64, 8>::filled(0.0);
    for &cut in cuts.iter() {
        if unique.last() != Some(&cut) {
            unique.push(cut);
        }
    }
    let mut intervals = Fixed::filled((0.0, 0.0, Kind::Water));
    let mut bounds = Fixed::filled((0, 0));
    for (j, c) in unique.windows(2).enumerate().filter(|(_, c)| c[1] > c[0]) {
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
        intervals.push((c[0], c[1], kind));
        bounds.push((j, j + 1));
    }
    Spectrum {
        cuts: unique,
        intervals,
        bounds,
    }
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

    /// Layer `k`'s profile at its Gauss node `g`: [`at`](Self::at) that
    /// node, from the memoized constants.
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
#[derive(Clone, Copy, Debug, Default, PartialEq)]
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
    let mut out = vec![LayerSaturation::default(); levels.count()];
    layer_saturation_into(params, gas, levels, ps, gravity, t, &mut out);
    out
}

/// [`layer_saturation`] into `out` (one per layer).
pub fn layer_saturation_into(
    params: &RadiationParams,
    gas: Gas,
    levels: &Levels,
    ps: f64,
    gravity: f64,
    t: &[f64],
    out: &mut [LayerSaturation],
) {
    let profile = Profile::new(params, gas, levels, ps, gravity, t);
    saturation_into(&profile, out);
}

fn saturation_into(profile: &Profile, out: &mut [LayerSaturation]) {
    let ps = profile.ps;
    let layers = (0..profile.levels.count()).map(|k| {
        let (a, b) = profile.levels.bounds(k);
        let (mut mean, mut shaped) = (0.0, 0.0);
        for (g, (x, w)) in GAUSS.into_iter().enumerate() {
            let s = 0.5 * (a + b) + 0.5 * (b - a) * x;
            let q = saturation_humidity(profile.at_node(k, g), s * ps);
            mean += 0.5 * w * q;
            shaped += 0.5 * w * humidity_shape(s) * q;
        }
        LayerSaturation { mean, shaped }
    });
    for (out, sat) in out.iter_mut().zip(layers) {
        *out = sat;
    }
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
    let mut out = vec![[0.0; 2]; levels.count() - 1];
    interface_humidity_ratios_into(params, gas, levels, ps, gravity, t, &mut out);
    out
}

/// [`interface_humidity_ratios`] into `out` (one per interior interface).
pub fn interface_humidity_ratios_into(
    params: &RadiationParams,
    gas: Gas,
    levels: &Levels,
    ps: f64,
    gravity: f64,
    t: &[f64],
    out: &mut [[f64; 2]],
) {
    let profile = Profile::new(params, gas, levels, ps, gravity, t);
    let mut sat = [LayerSaturation::default(); MAX_LAYERS];
    saturation_into(&profile, &mut sat[..levels.count()]);
    for ((j, &sigma), ratios) in levels.interior().iter().enumerate().zip(out) {
        *ratios = std::array::from_fn(|side| {
            let k = j + side;
            let at = humidity_shape(sigma)
                * saturation_humidity(profile.at_interface(j, side), sigma * ps);
            if sat[k].shaped > 0.0 {
                at / sat[k].shaped
            } else {
                1.0
            }
        });
    }
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

impl Slab {
    /// Whether `other` has the same optical depths in every interval.
    fn same_optics(&self, other: &Self) -> bool {
        self.top.to_bits() == other.top.to_bits()
            && self.bottom.to_bits() == other.bottom.to_bits()
            && self.weighted.to_bits() == other.weighted.to_bits()
            && self.continuum.to_bits() == other.continuum.to_bits()
    }
}

/// One element of the column, top to bottom.
#[derive(Clone, Copy, Debug)]
enum Element {
    Slab(usize),
    Sheet(usize),
}

/// Diffusivity factor: the mean slant path of diffuse light (air masses).
const DIFFUSIVITY: f64 = 1.66;

/// What a column's fluxes share whatever its clouds: its temperatures, its
/// CO₂ band, and in each spectral interval the Planck emission of the
/// layers' edges, of the surface and of the Rayleigh sheet.
#[derive(Clone, Copy, Debug)]
struct Shared {
    width: f64,
    spectrum: Spectrum,
    /// Per interval: at the top layer's top, each interface and the
    /// surface air (`N + 1`).
    edges: [[f64; MAX_LAYERS + 1]; MAX_INTERVALS],
    /// Per interval: the surface's and the Rayleigh sheet's.
    surface: [f64; MAX_INTERVALS],
    rayleigh: [f64; MAX_INTERVALS],
    /// The top slab's ozone transmissivities: direct and diffuse.
    ozone: [f64; 2],
}

/// One sky's work through a column: its sheets, slabs and elements, and the
/// fluxes as they add up.
#[derive(Clone, Copy, Debug)]
struct Sky {
    /// The clouds, then the clear air's Rayleigh sheet.
    sheets: Fixed<Cloud, MAX_SHEETS>,
    sheet_layer: [usize; MAX_SHEETS],
    /// Each sheet's temperature (the profile's at its level; unset for
    /// empty sheets, which aren't elements, and the Rayleigh sheet, whose
    /// emission is shared), and its Planck emission in the current interval.
    sheet_temperature: [f64; MAX_SHEETS],
    sheet_planck: [[f64; MAX_INTERVALS]; MAX_SHEETS],
    slabs: Fixed<Slab, MAX_SLABS>,
    elements: Fixed<Element, MAX_ELEMENTS>,
    /// Per slab: a slab of the reference sky with the same optics (whose
    /// transmissivities are this one's), if any.
    twin: [Option<usize>; MAX_SLABS],
    /// Per Gauss node of each slab: mass, σ and `shape × q*`.
    nodes: Fixed<(f64, f64, f64), { 4 * MAX_SLABS }>,
    shaped: [f64; MAX_LAYERS],
    amplitude: [f64; MAX_LAYERS],
    // Work arrays.
    optics: [(f64, f64, f64); MAX_ELEMENTS],
    direct: [f64; MAX_ELEMENTS + 1],
    reflect: [f64; MAX_ELEMENTS + 1],
    rising: [f64; MAX_ELEMENTS + 1],
    down: [f64; MAX_ELEMENTS + 1],
    absorbed: [f64; MAX_ELEMENTS],
    tau: [f64; MAX_SLABS],
    clear: [f64; MAX_SLABS],
    layer_tau: [f64; MAX_LAYERS],
    above: [f64; MAX_LAYERS],
    emit: [(f64, f64); MAX_ELEMENTS],
    rise: [f64; MAX_ELEMENTS + 1],
    fall: [f64; MAX_ELEMENTS + 1],
    // The fluxes.
    water: [f64; MAX_LAYERS],
    layers: [f64; MAX_LAYERS],
    surface_shortwave: f64,
    surface_visible: f64,
    escaped: f64,
    outgoing: f64,
    surface_down: f64,
    surface_up: f64,
}

impl Sky {
    const EMPTY: Self = Self {
        sheets: Fixed {
            items: [Cloud::NONE; MAX_SHEETS],
            len: 0,
        },
        sheet_layer: [0; MAX_SHEETS],
        sheet_temperature: [0.0; MAX_SHEETS],
        sheet_planck: [[0.0; MAX_INTERVALS]; MAX_SHEETS],
        slabs: Fixed {
            items: [Slab {
                layer: 0,
                top: 0.0,
                bottom: 0.0,
                water: 0.0,
                weighted: 0.0,
                continuum: 0.0,
            }; MAX_SLABS],
            len: 0,
        },
        elements: Fixed {
            items: [Element::Slab(0); MAX_ELEMENTS],
            len: 0,
        },
        twin: [None; MAX_SLABS],
        nodes: Fixed {
            items: [(0.0, 0.0, 0.0); 4 * MAX_SLABS],
            len: 0,
        },
        shaped: [0.0; MAX_LAYERS],
        amplitude: [0.0; MAX_LAYERS],
        optics: [(0.0, 0.0, 0.0); MAX_ELEMENTS],
        direct: [0.0; MAX_ELEMENTS + 1],
        reflect: [0.0; MAX_ELEMENTS + 1],
        rising: [0.0; MAX_ELEMENTS + 1],
        down: [0.0; MAX_ELEMENTS + 1],
        absorbed: [0.0; MAX_ELEMENTS],
        tau: [0.0; MAX_SLABS],
        clear: [0.0; MAX_SLABS],
        layer_tau: [0.0; MAX_LAYERS],
        above: [0.0; MAX_LAYERS],
        emit: [(0.0, 0.0); MAX_ELEMENTS],
        rise: [0.0; MAX_ELEMENTS + 1],
        fall: [0.0; MAX_ELEMENTS + 1],
        water: [0.0; MAX_LAYERS],
        layers: [0.0; MAX_LAYERS],
        surface_shortwave: 0.0,
        surface_visible: 0.0,
        escaped: 0.0,
        outgoing: 0.0,
        surface_down: 0.0,
        surface_up: 0.0,
    };
}

/// Working storage for the radiation of one column at a time, reused from
/// column to column (each call sets everything it reads).
#[derive(Clone, Debug)]
pub(crate) struct RadiationScratch {
    shared: Shared,
    /// The clear sky (the reference), and the sky with the column's clouds.
    skies: [Sky; 2],
}

impl RadiationScratch {
    pub(crate) fn new() -> Self {
        Self {
            shared: Shared {
                width: 0.0,
                spectrum: spectrum(0.0),
                edges: [[0.0; MAX_LAYERS + 1]; MAX_INTERVALS],
                surface: [0.0; MAX_INTERVALS],
                rayleigh: [0.0; MAX_INTERVALS],
                ozone: [0.0; 2],
            },
            skies: [Sky::EMPTY; 2],
        }
    }
}

/// Radiative fluxes of one column.
///
/// # Panics
/// If the column has more than [`MAX_CLOUDS`] cloud sheets.
pub fn column_fluxes(params: &RadiationParams, gas: Gas, column: &Column) -> Fluxes {
    let mut out = Fluxes::default();
    let mut scratch = RadiationScratch::new();
    let [_, sky] = &mut scratch.skies;
    let shared = &mut scratch.shared;
    let profile_t = shared_values(params, gas, column, shared);
    let profile = Profile::new(
        params,
        gas,
        column.levels,
        column.ps,
        column.gravity,
        &profile_t,
    );
    sky_setup(params, column, column.clouds, &profile, None, sky);
    sky_shortwave(params, column, shared, sky);
    sky_planck(shared, sky);
    let levels = column.levels;
    for (index, &(_, _, kind)) in shared.spectrum.intervals.iter().enumerate() {
        for &(share, scale) in terms(kind) {
            sky_longwave(params, levels, shared, index, kind, share, scale, None, sky);
        }
    }
    sky_fluxes(levels.count(), column, shared, sky, &mut out);
    out
}

/// The fluxes of `column` with its clouds, into `all`, and with none
/// (clear sky), into `clear`: the same values as [`column_fluxes`] of each,
/// sharing what doesn't depend on the clouds.
///
/// # Panics
/// If the column has more than [`MAX_CLOUDS`] cloud sheets.
pub(crate) fn sky_fluxes_into(
    params: &RadiationParams,
    gas: Gas,
    column: &Column,
    scratch: &mut RadiationScratch,
    all: &mut Fluxes,
    clear: &mut Fluxes,
) {
    let RadiationScratch { shared, skies } = scratch;
    let [reference, sky] = skies;
    let profile_t = shared_values(params, gas, column, shared);
    let profile = Profile::new(
        params,
        gas,
        column.levels,
        column.ps,
        column.gravity,
        &profile_t,
    );
    let none = [Cloud::NONE; MAX_CLOUDS];
    assert!(
        column.clouds.len() <= MAX_CLOUDS,
        "at most {MAX_CLOUDS} clouds"
    );
    sky_setup(
        params,
        column,
        &none[..column.clouds.len()],
        &profile,
        None,
        reference,
    );
    sky_setup(
        params,
        column,
        column.clouds,
        &profile,
        Some(reference),
        sky,
    );
    sky_shortwave(params, column, shared, reference);
    sky_shortwave(params, column, shared, sky);
    sky_planck(shared, reference);
    sky_planck(shared, sky);
    let levels = column.levels;
    for (index, &(_, _, kind)) in shared.spectrum.intervals.iter().enumerate() {
        for &(share, scale) in terms(kind) {
            sky_longwave(
                params, levels, shared, index, kind, share, scale, None, reference,
            );
            sky_longwave(
                params,
                levels,
                shared,
                index,
                kind,
                share,
                scale,
                Some(&reference.clear),
                sky,
            );
        }
    }
    let n = levels.count();
    sky_fluxes(n, column, shared, sky, all);
    sky_fluxes(n, column, shared, reference, clear);
}

/// The parts of an interval's band and their absorption coefficients.
fn terms(kind: Kind) -> &'static [(f64, f64)] {
    match kind {
        Kind::Water | Kind::WaterEdge => &WATER_TERMS,
        _ => &[(1.0, 1.0)],
    }
}

/// Fill `shared` for `column`, and return its layer temperatures (for the
/// profile).
fn shared_values(
    params: &RadiationParams,
    gas: Gas,
    column: &Column,
    shared: &mut Shared,
) -> [f64; MAX_LAYERS] {
    let levels = column.levels;
    let n_layers = levels.count();
    let mut exner = [0.0; MAX_LAYERS];
    levels.layer_exner_into(gas, column.ps, &mut exner[..n_layers]);
    let mut temperature = [0.0; MAX_LAYERS];
    for k in 0..n_layers {
        temperature[k] = column.theta[k] * exner[k];
    }
    let profile = Profile::new(params, gas, levels, column.ps, column.gravity, &temperature);
    // Edge temperatures: each interface that of the layer below it, along
    // the convective lapse rate through its mean; the lowest layer's bottom
    // the surface air; the top layer's top its own mean temperature (an
    // isothermal, stratosphere-like cap).
    let mut edge_temperature = [0.0; MAX_LAYERS + 1];
    edge_temperature[0] = temperature[0];
    for j in 0..n_layers - 1 {
        edge_temperature[j + 1] = profile.interface(j);
    }
    edge_temperature[n_layers] = profile.in_layer(n_layers - 1, 1.0);
    let rayleigh = profile.at(RAYLEIGH_SIGMA);
    let ts = column.surface_temperature;
    shared.width = co2_band_width(params, column);
    shared.spectrum = spectrum(shared.width);
    let spectrum = &shared.spectrum;
    let mut planck = [0.0; MAX_INTERVALS];
    for (k, &t) in edge_temperature[..=n_layers].iter().enumerate() {
        spectrum.planck(t, &mut planck);
        for (edges, &p) in shared
            .edges
            .iter_mut()
            .zip(&planck[..spectrum.intervals.len()])
        {
            edges[k] = p;
        }
    }
    spectrum.planck(rayleigh, &mut shared.rayleigh);
    spectrum.planck(ts, &mut shared.surface);
    shared.ozone =
        [1.0 / params.mu, DIFFUSIVITY].map(|m| fms_math::exp(-m * params.ozone_absorption));
    temperature
}

/// Black-body emission (W m⁻²) at `t`.
fn black(t: f64) -> f64 {
    STEFAN_BOLTZMANN * t * t * t * t
}

/// Level of the clear air's Rayleigh sheet: the atmosphere's mass centre.
const RAYLEIGH_SIGMA: f64 = 0.5;

/// Set up `sky` with cloud sheets `clouds`: its sheets, its slabs (the
/// layers cut at the sheets) with each slab's water, `q(p) = Aₖ · shape(σ) ·
/// q*(T(p), p)` where `Aₖ` is fixed ([`Humidity::Fixed`]) or set so the
/// layer holds its humidity, and its elements top to bottom. A slab with
/// the same edges as one of the `reference` sky has the same water at its
/// nodes.
fn sky_setup(
    params: &RadiationParams,
    column: &Column,
    clouds: &[Cloud],
    profile: &Profile,
    reference: Option<&Sky>,
    sky: &mut Sky,
) {
    let levels = column.levels;
    let n = levels.count();
    assert!(clouds.len() <= MAX_CLOUDS, "at most {MAX_CLOUDS} clouds");
    // The clouds, then the clear air's Rayleigh scattering.
    sky.sheets.clear();
    for &c in clouds {
        sky.sheets.push(c);
    }
    sky.sheets.push(Cloud {
        sigma: RAYLEIGH_SIGMA,
        reflectance: params.rayleigh,
        emissivity: 0.0,
    });
    let rayleigh = sky.sheets.len() - 1;
    for (c, sheet) in sky.sheets.iter().enumerate() {
        sky.sheet_layer[c] = levels.layer_of(sheet.sigma);
        sky.sheet_temperature[c] = if c == rayleigh || sheet.is_empty() {
            0.0
        } else {
            profile.at(sheet.sigma)
        };
    }
    let sheets = &sky.sheets;
    let mut cuts = Fixed::<f64, { MAX_LAYERS + 1 + MAX_SHEETS }>::filled(0.0);
    for &s in levels.interfaces() {
        cuts.push(s);
    }
    for c in sheets.iter() {
        if !c.is_empty() && c.sigma > 0.0 && c.sigma < 1.0 {
            cuts.push(c.sigma);
        }
    }
    cuts.sort_by(f64::total_cmp);
    sky.slabs.clear();
    let mut last: Option<f64> = None;
    for &cut in cuts.iter() {
        if let Some(top) = last {
            if top == cut {
                continue;
            }
            sky.slabs.push(Slab {
                layer: levels.layer_of(0.5 * (top + cut)),
                top,
                bottom: cut,
                ..Slab::default()
            });
        }
        last = Some(cut);
    }
    // Nodes: (mass, σ, shape × q*).
    sky.nodes.clear();
    sky.shaped[..n].fill(0.0);
    let interfaces = levels.interfaces();
    for slab in sky.slabs.iter() {
        let (a, b) = (slab.top, slab.bottom);
        // The same nodes in the reference sky (the same values), or a
        // whole layer's (whose profile there is memoized).
        let twin = reference.and_then(|r| {
            r.slabs
                .iter()
                .position(|s| s.top.to_bits() == a.to_bits() && s.bottom.to_bits() == b.to_bits())
        });
        let whole = a == interfaces[slab.layer] && b == interfaces[slab.layer + 1];
        for (g, (x, w)) in GAUSS.into_iter().enumerate() {
            let s = 0.5 * (a + b) + 0.5 * (b - a) * x;
            let mass = 0.5 * w * (b - a) * column.ps / column.gravity;
            let q = match (twin, reference) {
                (Some(j), Some(r)) => r.nodes[4 * j + g].2,
                _ => {
                    let t = if whole {
                        profile.at_node(slab.layer, g)
                    } else {
                        profile.at(s)
                    };
                    humidity_shape(s) * saturation_humidity(t, s * column.ps)
                }
            };
            sky.shaped[slab.layer] += mass * q;
            sky.nodes.push((mass, s, q));
        }
    }
    let ds = levels.dsigma();
    for k in 0..n {
        sky.amplitude[k] = match column.humidity {
            Humidity::Fixed => params.surface_humidity,
            Humidity::Layers(q) => {
                let water = q[k].max(0.0) * ds[k] * column.ps / column.gravity;
                if sky.shaped[k] > 0.0 {
                    water / sky.shaped[k]
                } else {
                    0.0
                }
            }
        };
    }
    let p0 = REFERENCE_PRESSURE;
    for (j, slab) in sky.slabs.iter_mut().enumerate() {
        for &(mass, s, q) in &sky.nodes[4 * j..4 * j + 4] {
            let q = sky.amplitude[slab.layer] * q;
            let p = s * column.ps;
            let e = (q * p / (EPSILON + (1.0 - EPSILON) * q)).min(0.25 * p);
            slab.water += mass * q;
            slab.weighted += mass * q * p / p0;
            slab.continuum += mass * q * e / 1000.0;
        }
    }
    for (j, slab) in sky.slabs.iter().enumerate() {
        sky.twin[j] = reference.and_then(|r| r.slabs.iter().position(|s| s.same_optics(slab)));
    }
    // Elements: slabs, with each cloud sheet at its level (sorted by level,
    // in sheet order where levels are equal).
    let mut order = Fixed::<usize, MAX_SHEETS>::filled(0);
    for c in 0..sheets.len() {
        if !sheets[c].is_empty() {
            order.push(c);
        }
    }
    order.sort_by(|&a, &b| sheets[a].sigma.total_cmp(&sheets[b].sigma));
    sky.elements.clear();
    let mut next = order.iter().copied().peekable();
    for (j, slab) in sky.slabs.iter().enumerate() {
        while let Some(&c) = next.peek() {
            if sheets[c].sigma <= slab.top {
                sky.elements.push(Element::Sheet(c));
                next.next();
            } else {
                break;
            }
        }
        sky.elements.push(Element::Slab(j));
    }
    for c in next {
        sky.elements.push(Element::Sheet(c));
    }
    sky.water[..n].fill(0.0);
    for slab in sky.slabs.iter() {
        sky.water[slab.layer] += slab.water;
    }
    sky.layers[..n].fill(0.0);
    sky.outgoing = 0.0;
    sky.surface_down = 0.0;
    sky.surface_up = 0.0;
}

/// Shortwave, in two bands: the visible (ozone at the top) and the
/// near-infrared (water vapour). Each band runs the adding method through
/// the same sheets and surface. Water's absorption saturates with the path
/// from the top: transmissivity `exp(−a√(mW))` along `m` air masses through
/// `W` of water, so a slab between cumulative paths `W₁` and `W₂` transmits
/// `exp(−a(√(mW₂) − √(mW₁)))` and the column's absorption doesn't depend on
/// how it's cut into slabs.
fn sky_shortwave(params: &RadiationParams, column: &Column, shared: &Shared, sky: &mut Sky) {
    let mut path = [0.0; MAX_SLABS + 1];
    for (j, slab) in sky.slabs.iter().enumerate() {
        path[j + 1] = path[j] + slab.water;
    }
    let s = column.insolation.max(0.0);
    sky.surface_shortwave = 0.0;
    sky.surface_visible = 0.0;
    sky.escaped = 0.0;
    for (b, share) in [
        (VISIBLE, 1.0 - params.near_infrared),
        (NEAR_INFRARED, params.near_infrared),
    ] {
        let transmit = |j: usize, m: f64, diffuse: usize| {
            if b == VISIBLE {
                if j == 0 { shared.ozone[diffuse] } else { 1.0 }
            } else {
                let a = params.water_absorption;
                fms_math::exp(-a * ((m * path[j + 1]).sqrt() - (m * path[j]).sqrt()))
            }
        };
        let (surface, escaped) = shortwave_band(sky, params.albedo[b], share * s, |j| {
            (transmit(j, 1.0 / params.mu, 0), transmit(j, DIFFUSIVITY, 1))
        });
        if b == VISIBLE {
            sky.surface_visible = surface;
        }
        for (j, e) in sky.elements.iter().enumerate() {
            if let Element::Slab(i) = *e {
                sky.layers[sky.slabs[i].layer] += sky.absorbed[j];
            }
        }
        sky.surface_shortwave += surface;
        sky.escaped += escaped;
    }
}

/// The adding method for one band of sunlight `s` (W m⁻²) through `sky`'s
/// elements: sheets reflect and transmit without absorbing; `slab(j)` gives
/// slab `j`'s direct and diffuse transmissivities; reflections between
/// everything are summed exactly. Each element's absorption goes into
/// `sky.absorbed` (zero for sheets); returns what the surface absorbs and
/// what escapes to space.
fn shortwave_band(
    sky: &mut Sky,
    albedo: f64,
    s: f64,
    slab: impl Fn(usize) -> (f64, f64),
) -> (f64, f64) {
    let n = sky.elements.len();
    // Per element: reflectance and direct and diffuse transmissivities.
    for (o, e) in sky.optics.iter_mut().zip(sky.elements.iter()) {
        *o = match *e {
            Element::Slab(j) => {
                let (direct, diffuse) = slab(j);
                (0.0, direct, diffuse)
            }
            Element::Sheet(c) => {
                let r = sky.sheets[c].reflectance;
                (r, 1.0 - r, 1.0 - r)
            }
        };
    }
    let optics = &sky.optics[..n];
    // Direct beam down; diffuse reflectance R̃ of everything below each
    // interface, and the diffuse light U leaving each interface upward that
    // the direct beam generates below it.
    let direct = &mut sky.direct;
    direct[0] = s;
    for (j, &(_, t_dir, _)) in optics.iter().enumerate() {
        direct[j + 1] = direct[j] * t_dir;
    }
    let (reflect, rising) = (&mut sky.reflect, &mut sky.rising);
    reflect[n] = albedo;
    rising[n] = albedo * direct[n];
    for j in (0..n).rev() {
        let (r, _, t) = optics[j];
        let denom = 1.0 - r * reflect[j + 1];
        reflect[j] = r + t * t * reflect[j + 1] / denom;
        rising[j] = t * rising[j + 1] / denom + r * direct[j];
    }
    let down = &mut sky.down;
    down[0] = 0.0;
    for j in 0..n {
        let (r, _, t) = optics[j];
        down[j + 1] = (t * down[j] + r * rising[j + 1]) / (1.0 - r * reflect[j + 1]);
    }
    let up = |j: usize| reflect[j] * down[j] + rising[j];
    for (j, absorbed) in sky.absorbed[..n].iter_mut().enumerate() {
        let (r, _, t) = optics[j];
        *absorbed = if r > 0.0 {
            0.0
        } else {
            (direct[j] - direct[j + 1]) + (1.0 - t) * (down[j] + up(j + 1))
        };
    }
    ((1.0 - albedo) * (down[n] + direct[n]), up(0))
}

/// The Planck emission of `sky`'s cloud sheets in each interval (the
/// Rayleigh sheet's is shared; empty sheets aren't elements).
fn sky_planck(shared: &Shared, sky: &mut Sky) {
    let rayleigh = sky.sheets.len() - 1;
    for (c, sheet) in sky.sheets.iter().enumerate() {
        if c != rayleigh && !sheet.is_empty() {
            shared
                .spectrum
                .planck(sky.sheet_temperature[c], &mut sky.sheet_planck[c]);
        }
    }
}

/// One part (`share` of the band, absorption coefficient `scale`) of
/// longwave interval `index` (of absorber `kind`) through `sky`, top to
/// bottom and back. Slabs with a twin take their transmissivities from
/// `reference` (the reference sky's, in this part).
#[allow(clippy::too_many_arguments)]
fn sky_longwave(
    params: &RadiationParams,
    levels: &Levels,
    shared: &Shared,
    index: usize,
    kind: Kind,
    share: f64,
    scale: f64,
    reference: Option<&[f64; MAX_SLABS]>,
    sky: &mut Sky,
) {
    let n_layers = levels.count();
    let ds = levels.dsigma();
    let edges = &shared.edges[index];
    let Sky {
        sheets,
        sheet_layer,
        sheet_planck,
        slabs,
        elements,
        twin,
        tau,
        clear,
        layer_tau,
        above,
        emit,
        rise,
        fall,
        layers,
        outgoing,
        surface_down,
        surface_up,
        ..
    } = sky;
    layer_tau[..n_layers].fill(0.0);
    for (t, slab) in tau.iter_mut().zip(slabs.iter()) {
        *t = match kind {
            Kind::Co2 => params.co2_opacity * (slab.bottom - slab.top),
            Kind::Window => params.continuum_opacity * slab.continuum,
            Kind::Water => scale * params.water_opacity * slab.weighted,
            Kind::WaterEdge => scale * params.water_edge * params.water_opacity * slab.weighted,
        };
        layer_tau[slab.layer] += *t;
    }
    for (i, (c, t)) in clear.iter_mut().zip(&tau[..slabs.len()]).enumerate() {
        *c = match (twin[i], reference) {
            (Some(j), Some(r)) => r[j],
            _ => fms_math::exp(-t),
        };
    }
    // Emission of each element up (from its top) and down (from its
    // bottom). A slab's source is linear in optical depth across its layer,
    // so its edges' values are interpolated in τ.
    above[..n_layers].fill(0.0);
    for (j, e) in elements.iter().enumerate() {
        emit[j] = match *e {
            Element::Slab(i) => {
                let slab = &slabs[i];
                let k = slab.layer;
                let (b_top, b_bottom) = (edges[k], edges[k + 1]);
                // Position of the slab's edges within its layer, in optical
                // depth (or σ, if the layer is transparent).
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
                let planck = if c == sheets.len() - 1 {
                    shared.rayleigh[index]
                } else {
                    sheet_planck[c][index]
                };
                let e = sheets[c].emissivity * share * planck;
                (e, e)
            }
        };
    }
    let transmit = |j: usize| match elements[j] {
        Element::Slab(i) => clear[i],
        Element::Sheet(c) => 1.0 - sheets[c].emissivity,
    };
    // Upward: surface → … → space; downward: space → … → surface.
    let n = elements.len();
    let u_surface = share * shared.surface[index];
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
            Element::Sheet(c) => sheet_layer[c],
        };
        layers[k] += (1.0 - transmit(j)) * (rise[j + 1] + fall[j]) - emit[j].0 - emit[j].1;
    }
    *outgoing += rise[0];
    *surface_down += fall[n];
    *surface_up += u_surface;
}

/// `sky`'s fluxes into `out`.
fn sky_fluxes(n: usize, column: &Column, shared: &Shared, sky: &Sky, out: &mut Fluxes) {
    let s = column.insolation.max(0.0);
    out.absorbed = s - sky.escaped;
    out.outgoing = sky.outgoing;
    out.layers.clear();
    out.layers.extend_from_slice(&sky.layers[..n]);
    out.surface_shortwave = sky.surface_shortwave;
    out.surface_visible = sky.surface_visible;
    out.surface_down = sky.surface_down;
    out.surface_up = sky.surface_up;
    out.water.clear();
    out.water.extend_from_slice(&sky.water[..n]);
    out.co2_width = shared.width;
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
    fn each_band_of_sunlight_reflects_off_its_own_albedo() {
        let params = RadiationParams {
            ozone_absorption: 0.0,
            rayleigh: 0.0,
            albedo: [0.5, 0.1],
            ..RadiationParams::default()
        };
        let column = Column {
            clouds: &[Cloud::NONE; 2],
            ..earthlike()
        };
        let base = column_fluxes(&params, Gas::EARTH_AIR, &column);
        let expected = (1.0 - 0.5) * (1.0 - params.near_infrared) * column.insolation;
        assert!((base.surface_visible - expected).abs() < 1e-9);
        let bright = RadiationParams {
            albedo: [0.5, 0.6],
            ..params
        };
        let more = column_fluxes(&bright, Gas::EARTH_AIR, &column);
        assert_eq!(more.surface_visible, base.surface_visible);
        assert!(more.surface_shortwave < base.surface_shortwave - 1.0);
        assert!(more.absorbed < base.absorbed - 1.0);
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
