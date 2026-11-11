//! Column physics for a radiatively driven planet: [`crate::radiation`],
//! [`crate::moisture`], and a surface of ocean and land tiles with bulk
//! sensible and latent heat fluxes. [`ColumnPhysics`] is the 3-D forcing; [`SingleColumn`]
//! is a 1-D radiative–convective model with fixed relative humidity, kept
//! to check the radiation against Earth's global-mean column.
//!
//! # Pieces (3-D)
//!
//! - **Surface**: each cell has an ocean, a sea-ice and a land tile, weighted
//!   by its land fraction and ice cover, each with its own temperature,
//!   albedo, roughness and fluxes. The ocean tile is the mixed layer of
//!   [`fms_ocean::Ocean`], which moves heat sideways and down (with the
//!   open water's wind stress, if it has Ekman transport), and the ice
//!   tile is its ice; the land tile is a skin with no heat
//!   capacity over a soil heat store, with a bucket of soil water and a pack
//!   of snow. Each absorbs sunlight and longwave, emits as a black
//!   body, and pays the latent heat of what evaporates; open water also pays
//!   the fusion of the snow that falls on it and melts. The radiation sees one surface per
//!   cell, with the tiles' mean albedo and mean σT⁴, and splits its fluxes
//!   between them so nothing is lost.
//! - **Skin and soil**: each step the skin temperature is solved so that
//!   what the skin absorbs balances emission, the turbulent fluxes and
//!   conduction into the soil, implicitly with the soil's own update
//!   (`solve_skin`). The soil then takes exactly the net flux the
//!   atmosphere was given, so a skin solved inexactly loses no energy.
//!   Radiation runs less often than the skin changes: the land tile emits at
//!   its current skin temperature, and the difference from what the last
//!   radiation call counted goes into the bottom layer, which is where most
//!   of the surface's longwave is absorbed.
//! - **Soil water**: Manabe's (1969) bucket. Rain and snow fill it, and what
//!   overflows its capacity runs off; evaporation falls in proportion to its
//!   water below a critical level. Runoff leaves the planet's water until
//!   rivers (M3).
//! - **Snow**: what falls on land piles up, brightening the tile as it
//!   covers it. A skin that would rise above freezing under snow stays at
//!   freezing, and the surplus melts the snow into the bucket. Evaporation
//!   takes snow first, paying its fusion too. The pack owes the heat of
//!   fusion the air was given when it froze, so the planet's stored heat
//!   counts it as a deficit, `−L_f × snow`. Land fully covered for a whole
//!   year is an ice sheet, brighter and smoother than seasonal snow; land
//!   starting cold enough starts as one, under a thick pack.
//! - **Boundary layer**: the bottom layer is a well-mixed boundary layer
//!   ([`crate::surface_layer`]): surface air has its potential temperature
//!   and specific humidity. Bulk formulas with stability-dependent
//!   coefficients and an ocean roughness give the surface's sensible heat,
//!   evaporation and drag on the bottom layer's wind; the layer above is
//!   entrained at mixed-layer theory's rate, exchanging dry static energy,
//!   water and momentum, and shallow cumulus carry its water and momentum
//!   up where its air saturates below its top.
//! - **Moisture**: convection, large-scale condensation, re-evaporation,
//!   rain or snow, and diagnosed clouds ([`crate::moisture`]).
//! - **Sub-grid hills**: wind blowing up the terrain's slope rains out a
//!   share of the moisture it carries (Smith 1979), rising with the spread
//!   of the cell's sub-grid terrain, so rain falls on windward slopes and not
//!   on lee ones. The spread also slows the wind by a drag whose lost kinetic
//!   energy heats the boundary layer like the rest of the drag.
//! - **Dissipation**: the kinetic energy drag and entrainment remove heats
//!   the boundary layer ([`crate::dissipation`]).
//!
//! # Time stepping
//!
//! Radiation runs every [`PhysicsParams::radiation_interval`] with sunlight
//! averaged over the coming interval and the clouds of the latest moist
//! physics; surface fluxes and the surface tiles update every step
//! ([`Forcing::advance`]). The atmosphere and the tiles use the same stored
//! fluxes and the same precipitation, so the planet's energy and water
//! budgets close exactly: what the atmosphere and surface gain is what was
//! absorbed from the Sun minus what was emitted to space (friction's lost
//! kinetic energy comes back as heat), and the atmosphere's water changes
//! by evaporation minus precipitation.

use fms_grid::Operators;
use fms_grid::operators::{MIN_CHUNK, MIN_COLUMNS};
use fms_grid::zonal::east_north;
use fms_ocean::{IceForcing, JumpLayers, Ocean, OceanParams, OceanSnapshot};
use fms_planet::Planet;
use glam::DVec3;
use rayon::prelude::*;

use crate::dissipation::KineticHeating;
use crate::forcing::{CellField, Forcing, ForcingBudget};
use crate::layered::{Gas, LayeredModel, LayeredState, REFERENCE_PRESSURE};
use crate::levels::{Levels, MAX_LAYERS};
use crate::moisture::{
    FUSION_HEAT, MoistColumn, MoistTendency, MoistureParams, Setting, VAPOUR_GAS_CONSTANT, lift,
    moist_physics,
};
use crate::radiation::{
    Cloud, Column, Fluxes, Humidity, NEAR_INFRARED, Profile, RadiationParams, STEFAN_BOLTZMANN,
    VISIBLE, column_fluxes, saturation_humidity, saturation_pressure, specific_humidity,
};
use crate::surface_layer::{
    Surface, SurfaceAir, SurfaceLayerParams, convective_velocity, entrainment_velocity, exchange,
};

pub use crate::moisture::LATENT_HEAT;

const DAY: f64 = 86_400.0;

/// The climate's interior layer interfaces (σ): four layers. The bottom one
/// (~1.4 km) is the boundary layer. The top is split at 250 hPa so that the
/// Hadley outflow (σ 0.25–0.5) is not averaged with the stratosphere, which
/// inflates its moist static energy and suppresses the subtropical eddy
/// momentum flux.
pub const CLIMATE_LEVELS: [f64; 3] = [0.25, 0.5, 0.85];
/// Index of each cell's ocean tile.
const OCEAN: usize = 0;
/// Index of each cell's land tile.
const LAND: usize = 1;
/// Index of each cell's sea-ice tile.
const ICE: usize = 2;
/// The temperature snow melts at (K).
const MELTING_POINT: f64 = 273.15;

/// The land tile. [`Default`] is a mean Earth-like land surface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LandParams {
    /// Surface albedo.
    pub albedo: f64,
    /// Momentum roughness length (m).
    pub roughness: f64,
    /// Depth of the soil heat store (m).
    pub soil_depth: f64,
    /// Volumetric heat capacity of the soil (J m⁻³ K⁻¹).
    pub soil_heat_capacity: f64,
    /// Thermal conductance between the skin and the soil store (W m⁻² K⁻¹):
    /// soil's ~1 W m⁻¹ K⁻¹ over the ~15 cm that the day's heating reaches.
    pub skin_conductance: f64,
    /// Most water the soil holds before it runs off (kg m⁻², or mm).
    pub soil_water_capacity: f64,
    /// Fraction of the capacity above which the soil evaporates freely.
    pub critical_wetness: f64,
    /// Albedo of land fully covered by snow, a mean over fresh and old snow
    /// and what vegetation shows through.
    pub snow_albedo: f64,
    /// Snow (kg m⁻²) that covers the land fully.
    pub snow_covering: f64,
    /// Albedo of an ice sheet: old, clean snow with nothing showing through.
    pub ice_sheet_albedo: f64,
    /// Momentum roughness length of an ice sheet (m).
    pub ice_sheet_roughness: f64,
    /// Years of unbroken full snow cover that make land an ice sheet;
    /// infinite for no ice sheets.
    pub ice_sheet_years: f64,
    /// Land starting colder than this (K) starts as an ice sheet. A sheet
    /// needs its bright surface to outlast the summer, so one can't grow
    /// from bare ground under seasonal snow's albedo.
    pub ice_sheet_start_temperature: f64,
    /// Snow (kg m⁻²) on land that starts as an ice sheet.
    pub ice_sheet_start_pack: f64,
}

impl Default for LandParams {
    fn default() -> Self {
        Self {
            albedo: 0.2,
            roughness: 0.1,
            soil_depth: 1.0,
            soil_heat_capacity: 2.5e6,
            skin_conductance: 7.0,
            soil_water_capacity: 150.0,
            critical_wetness: 0.75,
            snow_albedo: 0.65,
            snow_covering: 10.0,
            ice_sheet_albedo: 0.82,
            ice_sheet_roughness: 1.0e-3,
            ice_sheet_years: 1.0,
            ice_sheet_start_temperature: 267.0,
            ice_sheet_start_pack: 1000.0,
        }
    }
}

impl LandParams {
    /// Heat capacity of the soil store (J m⁻² K⁻¹).
    pub fn heat_capacity(&self) -> f64 {
        self.soil_depth * self.soil_heat_capacity
    }

    /// Soil water (kg m⁻²) at which evaporation starts to be limited.
    fn critical_water(&self) -> f64 {
        self.critical_wetness * self.soil_water_capacity
    }

    /// The albedo of land under `snow` (kg m⁻²).
    fn albedo_under(&self, snow: f64) -> f64 {
        let cover = (snow / self.snow_covering).min(1.0);
        self.albedo + cover * (self.snow_albedo - self.albedo)
    }

    /// The albedo of land under `snow` (kg m⁻²), which an ice sheet overrides.
    fn surface_albedo(&self, snow: f64, ice_sheet: bool) -> f64 {
        if ice_sheet {
            self.ice_sheet_albedo
        } else {
            self.albedo_under(snow)
        }
    }

    /// The momentum roughness length of the land (m).
    fn surface_roughness(&self, ice_sheet: bool) -> f64 {
        if ice_sheet {
            self.ice_sheet_roughness
        } else {
            self.roughness
        }
    }

    /// How long land has been fully covered (s), `age` before a step of `dt`
    /// that left `snow` (kg m⁻²) on it. Losing full cover restarts the clock.
    fn cover_age(&self, age: f64, snow: f64, dt: f64) -> f64 {
        if snow >= self.snow_covering {
            age + dt
        } else {
            0.0
        }
    }

    /// The share (0–1) of the potential evaporation that soil holding
    /// `water` (kg m⁻², with any snow) supplies.
    fn wetness(&self, water: f64) -> f64 {
        (water / self.critical_water()).min(1.0)
    }

    /// Soil `water` (kg m⁻²) after gaining `gain` (kg m⁻²; a loss if
    /// negative): the new amount and what overflowed. Water can't run out:
    /// evaporation is proportional to soil water and snow below the critical
    /// level, and a step evaporates far less than that (~2 of 112 mm at
    /// most), taking snow first.
    fn fill(&self, water: f64, gain: f64) -> (f64, f64) {
        let new = water + gain;
        let overflow = (new - self.soil_water_capacity).max(0.0);
        (new - overflow, overflow)
    }

    /// The conductance the skin feels over a step of `dt` seconds. The soil
    /// warms during the step by what it conducts, which is solved together
    /// with the skin: `Λ (T_skin − T_soil⁺)` with
    /// `T_soil⁺ = T_soil + Δt G / C` is `Λ/(1 + Λ Δt/C) (T_skin − T_soil)`.
    fn step_conductance(&self, dt: f64) -> f64 {
        self.skin_conductance / (1.0 + self.skin_conductance * dt / self.heat_capacity())
    }
}

/// Sub-grid hills. [`Default`] is Earth-like.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrographyParams {
    /// Share of the moisture carried up a slope that rains out, over hills
    /// at least [`full_relief`](Self::full_relief) high (Smith 1979's `ε`).
    pub rain_efficiency: f64,
    /// Height spread (m) at which the hills lift the air fully.
    pub full_relief: f64,
    /// Horizontal spacing of the hills, for their steepness (m).
    pub hill_spacing: f64,
    /// Drag coefficient of the hills per unit of their steepness, height
    /// over spacing; a simplified Lott & Miller (1997).
    pub drag: f64,
    /// Strongest rate of the hills' drag (s⁻¹), for the explicit step.
    pub max_drag_rate: f64,
}

impl Default for OrographyParams {
    fn default() -> Self {
        Self {
            rain_efficiency: 0.5,
            full_relief: 500.0,
            hill_spacing: 100e3,
            drag: 1.0,
            max_drag_rate: 1.0 / (3.0 * 3600.0),
        }
    }
}

/// Column-physics parameters. [`Default`] is Earth-like.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhysicsParams {
    /// Radiation.
    pub radiation: RadiationParams,
    /// Moisture (3-D only).
    pub moisture: MoistureParams,
    /// CO₂ (ppm).
    pub co2: f64,
    /// The ocean under the ocean tile.
    pub ocean: OceanParams,
    /// Shortening the ocean's spin-up.
    pub spin_up: OceanSpinUp,
    /// The land tile.
    pub land: LandParams,
    /// Sub-grid hills.
    pub orography: OrographyParams,
    /// The surface layer and boundary-layer top (3-D).
    pub surface_layer: SurfaceLayerParams,
    /// The single column's bulk transfer coefficient for heat and moisture.
    pub transfer: f64,
    /// The single column's gust speed, added in quadrature to the wind
    /// (m s⁻¹).
    pub gust: f64,
    /// The single column's lapse rate (K m⁻¹) from its lower layer's mean
    /// down to the surface air, and its critical lapse rate for convection
    /// between its two layer means. Earth's observed layer means (272 and
    /// 226 K, 10.5 km apart) sit at 4.4 K/km: less than the tropospheric
    /// 6.5 K/km because the upper layer includes the stratosphere.
    pub lapse_rate: f64,
    /// The single column's convective relaxation time (s).
    pub convection_time: f64,
    /// Time between radiation calls (s).
    pub radiation_interval: f64,
}

/// Asynchronous spin-up of the ocean: `jumps` years, the first ending
/// `first_year` years into the run and each next one `every` years after the
/// last, end in a jump, which repeats the change of the `layers` over that
/// whole year for `extra_years` more. The years between jumps run in true
/// time. Runs that report a climate need unaccelerated years after the last
/// jump.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OceanSpinUp {
    /// How many years end in a jump; 0 for none.
    pub jumps: u32,
    /// Years of change added by each jump.
    pub extra_years: f64,
    /// Most any layer of any cell may change in one jump (K).
    pub limit: f64,
    /// Which layers jump.
    pub layers: JumpLayers,
    /// The first year that ends in a jump, counting the run's first year
    /// as 1; at least 1.
    pub first_year: u32,
    /// Years from one jump to the next; at least 1.
    pub every: u32,
}

impl Default for OceanSpinUp {
    fn default() -> Self {
        Self {
            jumps: 0,
            extra_years: 3.0,
            limit: 3.0,
            layers: JumpLayers::Both,
            first_year: 2,
            every: 1,
        }
    }
}

impl OceanSpinUp {
    /// The spin-up for an ocean with currents: the deep layer, the slow
    /// store, jumps 20 years at the ends of years 4 and 7, and the mixed
    /// layer follows it. Year 4 is after the first years' transient and the
    /// first ice has melted; two true-time years between the jumps let the
    /// currents settle after the first, which a jump would otherwise
    /// amplify.
    pub fn with_currents() -> Self {
        Self {
            jumps: 2,
            extra_years: 20.0,
            limit: 3.0,
            layers: JumpLayers::Deep,
            first_year: 4,
            every: 3,
        }
    }

    /// The spin-up for an ocean without currents: both layers jump 3 years
    /// at the ends of years 2, 3 and 4, mostly to hurry the mixed layer and
    /// the ice.
    pub fn still() -> Self {
        Self {
            jumps: 3,
            ..Self::default()
        }
    }

    /// Whether the run's `year`th year (from 1) ends in a jump.
    pub fn ends_in_jump(&self, year: u32) -> bool {
        let Some(since) = year.checked_sub(self.first_year) else {
            return false;
        };
        self.every > 0 && since % self.every == 0 && since / self.every < self.jumps
    }

    /// The year that ends in the last jump, or `None` without jumps.
    pub fn last_jump_year(&self) -> Option<u32> {
        (self.jumps > 0).then(|| {
            (self.jumps - 1)
                .saturating_mul(self.every)
                .saturating_add(self.first_year)
        })
    }
}

impl Default for PhysicsParams {
    fn default() -> Self {
        Self {
            radiation: RadiationParams::default(),
            moisture: MoistureParams::default(),
            co2: 420.0,
            ocean: OceanParams::default(),
            spin_up: OceanSpinUp::default(),
            land: LandParams::default(),
            orography: OrographyParams::default(),
            surface_layer: SurfaceLayerParams::default(),
            transfer: 1.3e-3,
            gust: 3.0,
            lapse_rate: 4.4e-3,
            convection_time: 2.0 * 3600.0,
            radiation_interval: 3.0 * 3600.0,
        }
    }
}

impl PhysicsParams {
    /// Heat capacity of the ocean tile's mixed layer (J m⁻² K⁻¹).
    pub fn heat_capacity(&self) -> f64 {
        self.ocean.mixed_capacity()
    }

    /// Heat capacity of the land tile's soil store (J m⁻² K⁻¹).
    pub fn land_heat_capacity(&self) -> f64 {
        self.land.heat_capacity()
    }
}

/// Sea-surface temperature profiles of the Aqua-Planet Experiment (Neale &
/// Hoskins 2000; Blackburn et al. 2013), zonally and hemispherically
/// symmetric, 0 °C poleward of 60°.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ApeSst {
    /// `27 (1 − sin²(3φ/2))` °C: APE's main experiment.
    Control,
    /// The mean of CONTROL and FLAT (`27 (1 − sin⁴(3φ/2))` °C): closest to
    /// the observed Pacific profile.
    Qobs,
}

impl ApeSst {
    /// Sea-surface temperature (K) at latitude `lat` (radians).
    pub fn at(self, lat: f64) -> f64 {
        let celsius = if lat.abs() < std::f64::consts::FRAC_PI_3 {
            let s2 = fms_math::sin(1.5 * lat).powi(2);
            match self {
                Self::Control => 27.0 * (1.0 - s2),
                Self::Qobs => 27.0 * (1.0 - 0.5 * (s2 + s2 * s2)),
            }
        } else {
            0.0
        };
        273.15 + celsius
    }
}

/// Air temperature at the surface (K): the lowest layer's mean temperature
/// extrapolated down the convective lapse rate, given that layer's Exner
/// function `exner` and potential temperature `theta_lowest`.
pub fn surface_air_temperature(
    params: &PhysicsParams,
    gas: Gas,
    gravity: f64,
    ps: f64,
    exner: f64,
    theta_lowest: f64,
) -> f64 {
    let p_mean = REFERENCE_PRESSURE * fms_math::pow(exner, 1.0 / gas.kappa());
    let beta = gas.gas_constant * params.lapse_rate / gravity;
    theta_lowest * exner * fms_math::pow(ps / p_mean, beta)
}

/// Sensible and latent heat flux from the surface (W m⁻², positive upward)
/// for surface air at temperature `air`.
pub fn surface_fluxes(
    params: &PhysicsParams,
    gas: Gas,
    ps: f64,
    surface_temperature: f64,
    air: f64,
    wind: f64,
) -> (f64, f64) {
    let density = ps / (gas.gas_constant * air);
    let speed = (wind * wind + params.gust * params.gust).sqrt();
    let k = density * params.transfer * speed;
    let sensible = k * gas.heat_capacity * (surface_temperature - air);
    let q_surface = specific_humidity(saturation_pressure(surface_temperature), ps);
    let q_air = specific_humidity(
        params.radiation.surface_humidity * saturation_pressure(air),
        ps,
    );
    (sensible, k * LATENT_HEAT * (q_surface - q_air))
}

/// Potential temperatures a two-layer column (interface `sigma`) relaxes to
/// if it is less stable than the critical lapse rate, else `None`. The
/// adjusted column keeps its enthalpy `Σ δσₖ Πₖ θₖ`.
pub fn convective_target(
    params: &PhysicsParams,
    gas: Gas,
    sigma: f64,
    gravity: f64,
    ps: f64,
    theta: [f64; 2],
) -> Option<[f64; 2]> {
    let exner = Levels::two(sigma).layer_exner(gas, ps);
    let exner = [exner[0], exner[1]];
    let t = [theta[0] * exner[0], theta[1] * exner[1]];
    // Thickness between the levels whose Exner functions are the layer means.
    let p = exner.map(|e| REFERENCE_PRESSURE * fms_math::pow(e, 1.0 / gas.kappa()));
    let dz = gas.gas_constant * 0.5 * (t[0] + t[1]) / gravity * fms_math::ln(p[1] / p[0]);
    let drop = params.lapse_rate * dz;
    if t[1] - t[0] <= drop {
        return None;
    }
    let ds = [sigma, 1.0 - sigma];
    let enthalpy = ds[0] * t[0] + ds[1] * t[1];
    // T₀ = T₁ − drop with the same enthalpy.
    let t1 = enthalpy + ds[0] * drop;
    Some([(t1 - drop) / exner[0], t1 / exner[1]])
}

/// Heating rate of `Θₖ = p_s θₖ` (Pa K s⁻¹) from `q` W m⁻² absorbed by a
/// layer.
pub(crate) fn theta_tendency(gas: Gas, gravity: f64, dsigma: f64, exner: f64, q: f64) -> f64 {
    gravity * q / (gas.heat_capacity * dsigma * exner)
}

/// One tile's exchange with the bottom layer.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct TileFlux {
    /// Sensible heat (W m⁻², upward).
    sensible: f64,
    /// Evaporation (kg m⁻² s⁻¹).
    evaporation: f64,
    /// Drag on the bottom layer's wind (s⁻¹).
    drag: f64,
    /// 10 m wind speed (m s⁻¹) and surface stress (N m⁻²).
    wind: f64,
    stress: f64,
    /// Eastward and northward stress of the surface on the air (N m⁻²).
    stress_east: f64,
    stress_north: f64,
}

impl TileFlux {
    /// The mean of an ocean and a land tile, weighted by land `fraction`.
    fn blend(ocean: Self, land: Self, fraction: f64) -> Self {
        if fraction == 0.0 {
            return ocean;
        }
        if fraction == 1.0 {
            return land;
        }
        let mean = |a: f64, b: f64| (1.0 - fraction) * a + fraction * b;
        Self {
            sensible: mean(ocean.sensible, land.sensible),
            evaporation: mean(ocean.evaporation, land.evaporation),
            drag: mean(ocean.drag, land.drag),
            wind: mean(ocean.wind, land.wind),
            stress: mean(ocean.stress, land.stress),
            stress_east: mean(ocean.stress_east, land.stress_east),
            stress_north: mean(ocean.stress_north, land.stress_north),
        }
    }
}

/// Most Newton iterations of the skin balance, and the change (K) at which
/// it has converged.
const SKIN_ITERATIONS: usize = 8;
const SKIN_TOLERANCE: f64 = 1e-3;
/// Largest change of the skin temperature in one iteration (K): the exchange
/// coefficients jump as the surface layer changes stability.
const SKIN_MAX_CHANGE: f64 = 10.0;

/// The snow (kg m⁻², at most `snow`) that `surplus` heat (W m⁻²) melts in
/// `dt` seconds.
fn snow_melted(snow: f64, surplus: f64, dt: f64) -> f64 {
    (surplus.max(0.0) * dt / FUSION_HEAT).min(snow)
}

/// Black-body emission (W m⁻²) at temperature `t`.
fn black_body(t: f64) -> f64 {
    STEFAN_BOLTZMANN * t * t * t * t
}

/// The skin temperature (K) where `absorbed` radiation (W m⁻²) balances
/// black-body emission, conduction `conductance (T − soil)` and
/// `turbulent(T)`, which returns the upward turbulent flux (W m⁻²) and an
/// estimate of its slope (W m⁻² K⁻¹). A secant iteration from `start`, safe
/// against overshooting: the exchange coefficients change with stability, so
/// a slope that holds them fixed is out by a factor of two or more in light
/// winds and at night. Whatever imbalance is left goes to the soil.
fn solve_skin(
    start: f64,
    absorbed: f64,
    soil: f64,
    conductance: f64,
    turbulent: impl Fn(f64) -> (f64, f64),
) -> f64 {
    let mut t = start;
    // The root lies above any skin that gains heat and below any that loses
    // it. A step outside those bounds bisects them instead.
    let (mut low, mut high) = (f64::NEG_INFINITY, f64::INFINITY);
    let mut previous = (f64::NAN, f64::NAN);
    for _ in 0..SKIN_ITERATIONS {
        let (flux, slope) = turbulent(t);
        let residual = absorbed - black_body(t) - flux - conductance * (t - soil);
        if residual >= 0.0 {
            low = t;
        } else {
            high = t;
        }
        // After the first pass, the slope actually seen since the last one.
        let frozen = 4.0 * STEFAN_BOLTZMANN * t * t * t + slope + conductance;
        let secant = (previous.0 - residual) / (t - previous.1);
        let slope = if secant > 0.0 && secant.is_finite() {
            secant
        } else {
            frozen
        };
        previous = (residual, t);
        let newton = t + (residual / slope).clamp(-SKIN_MAX_CHANGE, SKIN_MAX_CHANGE);
        let next = if (low..=high).contains(&newton) {
            newton
        } else {
            0.5 * (low + high)
        };
        let change = next - t;
        t = next;
        if change.abs() < SKIN_TOLERANCE {
            break;
        }
    }
    t
}

/// The surface radiation of one cell's tiles (W m⁻²): sunlight absorbed and
/// longwave emitted, per unit area of each tile.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct TileRadiation {
    shortwave: [f64; 3],
    longwave_up: [f64; 3],
}

/// What the radiation sees of a cell's tiles.
#[derive(Clone, Copy, Debug, PartialEq)]
struct TiledSurface {
    /// Land fraction.
    fraction: f64,
    /// Share of the rest under sea ice.
    ice: f64,
    /// Each tile's albedo in each band of sunlight.
    albedo: [[f64; 2]; 3],
    /// Each tile's temperature (K).
    temperature: [f64; 3],
}

impl TiledSurface {
    /// Each tile's share of the cell.
    fn weights(&self) -> [f64; 3] {
        let sea = 1.0 - self.fraction;
        let mut w = [0.0; 3];
        w[OCEAN] = sea * (1.0 - self.ice);
        w[LAND] = self.fraction;
        w[ICE] = sea * self.ice;
        w
    }

    /// Area-weighted albedo in each band.
    fn mean_albedo(&self) -> [f64; 2] {
        let w = self.weights();
        std::array::from_fn(|b| {
            w[OCEAN] * self.albedo[OCEAN][b]
                + w[LAND] * self.albedo[LAND][b]
                + w[ICE] * self.albedo[ICE][b]
        })
    }

    /// `T⁴` of each tile, and their area-weighted mean.
    fn fourth_powers(&self) -> ([f64; 3], f64) {
        let fourth = self.temperature.map(|t| t * t * t * t);
        let w = self.weights();
        (
            fourth,
            w[OCEAN] * fourth[OCEAN] + w[LAND] * fourth[LAND] + w[ICE] * fourth[ICE],
        )
    }

    /// The single temperature whose black-body emission is the tiles' mean.
    fn radiating_temperature(&self) -> f64 {
        let w = self.weights();
        if w[OCEAN] == 1.0 {
            self.temperature[OCEAN]
        } else if w[LAND] == 1.0 {
            self.temperature[LAND]
        } else {
            self.fourth_powers().1.sqrt().sqrt()
        }
    }

    /// Divide the cell's surface `shortwave`, of which `visible` is in the
    /// visible band, and `longwave_up` between the tiles, conserving both:
    /// in each band a tile absorbs in proportion to `1 − α`, and it emits in
    /// proportion to `T⁴`.
    fn split(&self, shortwave: f64, visible: f64, longwave_up: f64) -> TileRadiation {
        let mean_albedo = self.mean_albedo();
        let share = |a: &[f64; 2], b: usize| (1.0 - a[b]) / (1.0 - mean_albedo[b]);
        let (fourth, mean) = self.fourth_powers();
        TileRadiation {
            // Written as a correction so that equal bands cost no round-off.
            shortwave: self.albedo.map(|a| {
                let infrared = share(&a, NEAR_INFRARED);
                shortwave * infrared + visible * (share(&a, VISIBLE) - infrared)
            }),
            longwave_up: fourth.map(|t| longwave_up * (t / mean)),
        }
    }
}

/// What the top of one column's boundary layer depends on.
pub(crate) struct BoundaryColumn {
    /// Gas, gravity (m s⁻²), and the lapse rate (K m⁻¹) along which the layer
    /// above reaches the top.
    pub gas: Gas,
    /// See `gas`.
    pub gravity: f64,
    /// See `gas`.
    pub lapse_rate: f64,
    /// σ of the layer's top, and its depth (m).
    pub top: f64,
    /// See `top`.
    pub depth: f64,
    /// Surface pressure (Pa).
    pub ps: f64,
    /// The layer's potential temperature (K), its specific humidity floored
    /// at zero, and that humidity as the state holds it.
    pub theta: f64,
    /// See `theta`.
    pub humidity: f64,
    /// See `theta`.
    pub layer_humidity: f64,
    /// Surface air temperature (K) and density (kg m⁻³).
    pub air: f64,
    /// See `air`.
    pub density: f64,
    /// The layer above's potential temperature (K), specific humidity and
    /// layer Exner function.
    pub theta_above: f64,
    /// See `theta_above`.
    pub humidity_above: f64,
    /// See `theta_above`.
    pub exner_above: f64,
    /// Surface sensible heat (W m⁻²), evaporation (kg m⁻² s⁻¹), and clouds'
    /// longwave cooling of the layer (W m⁻²).
    pub sensible: f64,
    /// See `sensible`.
    pub evaporation: f64,
    /// See `sensible`.
    pub cloud_top_cooling: f64,
}

/// The exchange at the top of one column's boundary layer.
pub(crate) struct BoundaryTop {
    /// The convective velocity `w*` and entrainment velocity (m s⁻¹).
    pub convective_velocity: f64,
    /// See `convective_velocity`.
    pub entrainment_velocity: f64,
    /// Whether shallow cumulus rise through the top.
    pub cumulus: bool,
    /// Density of the air at the top (kg m⁻³).
    pub density: f64,
}

/// Entrainment and shallow cumulus at the top of the boundary layer, as a
/// function of one column.
pub(crate) fn boundary_layer_top(c: &BoundaryColumn, sl: &SurfaceLayerParams) -> BoundaryTop {
    let (gas, g) = (c.gas, c.gravity);
    let (rd, cp) = (gas.gas_constant, gas.heat_capacity);
    let (ps, top) = (c.ps, c.top);
    // Clouds' longwave cooling of the layer drives turbulence from its top.
    let w_star = convective_velocity(
        g,
        c.depth,
        c.theta,
        c.density,
        cp,
        c.sensible,
        c.evaporation,
        c.cloud_top_cooling,
    );
    // Entrainment of the layer above, by the jump in virtual θ at the
    // boundary layer's top: the layer above extrapolated down to it along
    // the profile's lapse rate.
    let p_above = REFERENCE_PRESSURE * fms_math::pow(c.exner_above, 1.0 / gas.kappa());
    let beta = rd * sl.top_lapse_factor * c.lapse_rate / g;
    let t_top = c.theta_above * c.exner_above * fms_math::pow(top * ps / p_above, beta);
    let virtual_factor = |q: f64| 1.0 + 0.608 * q.max(0.0);
    let (below_factor, above_factor) = (
        virtual_factor(c.layer_humidity),
        virtual_factor(c.humidity_above),
    );
    let below = c.theta * below_factor;
    let jump = t_top / gas.exner(top * ps) * above_factor - below;
    let w_e = entrainment_velocity(sl, g, c.depth, w_star, below, jump);
    // Shallow cumulus: surface air saturated by the top and still buoyant
    // there.
    let p_top = top * ps;
    let (t_parcel, q_parcel) = lift(gas, c.air, c.humidity, ps, p_top);
    let cumulus =
        q_parcel < c.humidity && t_parcel * (1.0 + 0.608 * q_parcel) > t_top * above_factor;
    BoundaryTop {
        convective_velocity: w_star,
        entrainment_velocity: w_e,
        cumulus,
        density: p_top / (rd * t_top),
    }
}

/// Exchange between the surface, the boundary layer and the layer above in
/// one cell: the area-weighted mean of its tiles' fluxes, and each tile's.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct SurfaceFlux {
    /// Each tile's flux; unset for a tile with no area.
    tiles: [TileFlux; 3],
    /// The land tile's skin temperature these fluxes were computed at (K).
    skin: f64,
    /// The ice tile's surface temperature they were computed at (K).
    ice_skin: f64,
    /// The bottom layer's air brought dry-adiabatically to the surface (K).
    air: f64,
    /// Sensible heat (W m⁻², upward).
    sensible: f64,
    /// Evaporation (kg m⁻² s⁻¹).
    evaporation: f64,
    /// Drag on the bottom layer's wind (s⁻¹).
    drag: f64,
    /// Ascent of the bottom layer's wind up the terrain's slope, in the
    /// share the sub-grid hills lift (m s⁻¹; zero on a lee slope).
    upslope: f64,
    /// Air entrained across the boundary layer's top (kg m⁻² s⁻¹).
    entrainment: f64,
    /// Shallow cumulus mass flux out of the boundary layer (kg m⁻² s⁻¹).
    cumulus: f64,
    /// Entrainment velocity (m s⁻¹).
    entrainment_velocity: f64,
    /// The boundary layer's convective velocity `w*` (m s⁻¹).
    convective_velocity: f64,
    /// 10 m wind speed (m s⁻¹) and surface stress (N m⁻²).
    wind: f64,
    stress: f64,
    /// Water carried out of the boundary layer by entrainment and shallow
    /// cumulus (kg m⁻² s⁻¹).
    export: f64,
    /// Eastward and northward stress of the surface on the air (N m⁻²).
    stress_east: f64,
    stress_north: f64,
}

/// The radiative–convective forcing of the layered model: radiation, moist
/// physics, the surface tiles, surface fluxes and boundary-layer friction. Needs
/// a moist state.
pub struct ColumnPhysics<'a> {
    params: PhysicsParams,
    gas: Gas,
    levels: Levels,
    gravity: f64,
    ops: Operators<'a>,
    planet: &'a Planet,
    /// Model time (s since t = 0).
    time: f64,
    next_radiation: f64,
    /// Land fraction of each cell.
    land_fraction: Vec<f64>,
    /// Sub-grid spread of each cell's surface height (m).
    relief: Vec<f64>,
    /// Gradient of each cell's surface height (m per m).
    relief_slope: Vec<DVec3>,
    /// The ocean under each cell's ocean tile.
    ocean: Ocean,
    /// When the current model year ends (s).
    next_year: f64,
    /// The ocean at the start of the year being measured for a jump.
    year_start: Option<OceanSnapshot>,
    /// Model years that have ended, while jumps remain.
    years_ended: u32,
    jumps_done: u32,
    /// Heat the jumps have added to the ocean (J); not part of the budget.
    jumped: f64,
    /// Skin temperature of each cell's land tile (K).
    land_temperature: Vec<f64>,
    /// Temperature of each cell's soil store (K).
    soil_temperature: Vec<f64>,
    /// Water in each cell's soil (kg per m² of land).
    soil_water: Vec<f64>,
    /// Snow on each cell's land (kg per m² of land).
    snow: Vec<f64>,
    /// How long each cell's land has been fully covered by snow (s).
    cover_age: Vec<f64>,
    /// The unbroken cover that makes land an ice sheet (s).
    ice_sheet_age: f64,
    /// Last radiation call, all-sky and clear-sky, and its division between
    /// the tiles.
    radiation: Vec<Fluxes>,
    clear: Vec<Fluxes>,
    tile_radiation: Vec<TileRadiation>,
    /// Surface fluxes, from the last step.
    surface: Vec<SurfaceFlux>,
    /// The moist physics of the last tendency evaluation.
    moist: Vec<MoistTendency>,
    heating: KineticHeating,
    /// Whether the surface temperature is held fixed (an infinite heat
    /// reservoir) instead of evolving.
    prescribed: bool,
    /// Heat the surface has taken up since the start (J).
    uptake: f64,
    /// Kinetic energy turned to heat by the last tendency evaluation (W),
    /// and since the start (J).
    dissipation: f64,
    dissipated: f64,
    absorbed: f64,
    outgoing: f64,
    evaporated: f64,
    precipitated: f64,
    runoff: f64,
    ocean_exchange: f64,
}

impl<'a> ColumnPhysics<'a> {
    /// Physics for `model` on `planet` with no land, starting at time `time`
    /// (s) from `state` with surface temperatures `surface_temperature` (K).
    ///
    /// # Panics
    /// If `state` is dry, or the spin-up has jumps with a first year or a
    /// spacing of 0.
    pub fn new(
        model: &LayeredModel<'a>,
        planet: &'a Planet,
        params: PhysicsParams,
        time: f64,
        state: &LayeredState,
        surface_temperature: Vec<f64>,
    ) -> Self {
        let no_land = vec![0.0; surface_temperature.len()];
        Self::with_land(
            model,
            planet,
            params,
            time,
            state,
            surface_temperature,
            no_land,
        )
    }

    /// As [`new`](Self::new), with land covering `land_fraction` (0–1) of
    /// each cell; both tiles, and the soil, start at `surface_temperature`.
    ///
    /// # Panics
    /// If `state` is dry, a land fraction is outside 0–1, or the spin-up has
    /// jumps with a first year or a spacing of 0.
    pub fn with_land(
        model: &LayeredModel<'a>,
        planet: &'a Planet,
        params: PhysicsParams,
        time: f64,
        state: &LayeredState,
        surface_temperature: Vec<f64>,
        land_fraction: Vec<f64>,
    ) -> Self {
        assert!(state.is_moist(), "column physics needs a moist state");
        let n = state.ps.len();
        assert_eq!(surface_temperature.len(), n);
        assert_eq!(land_fraction.len(), n);
        assert!(
            land_fraction.iter().all(|f| (0.0..=1.0).contains(f)),
            "land fractions must be within 0–1"
        );
        let spin_up = params.spin_up;
        assert!(
            spin_up.jumps == 0 || (spin_up.first_year >= 1 && spin_up.every >= 1),
            "ocean jumps need a first year and a spacing of at least 1"
        );
        let ocean = Ocean::new(
            params.ocean,
            model.operators(),
            &land_fraction,
            surface_temperature.clone(),
            planet.rotation_rate(),
        );
        let mut physics = Self {
            params,
            gas: model.gas(),
            levels: model.levels().clone(),
            gravity: model.gravity(),
            ops: *model.operators(),
            planet,
            time,
            next_radiation: time,
            land_fraction,
            relief: vec![0.0; n],
            relief_slope: vec![DVec3::ZERO; n],
            land_temperature: surface_temperature.clone(),
            soil_temperature: surface_temperature,
            soil_water: vec![params.land.critical_water(); n],
            snow: vec![0.0; n],
            cover_age: vec![0.0; n],
            ice_sheet_age: params.land.ice_sheet_years * planet.year_seconds(),
            ocean,
            next_year: time + planet.year_seconds(),
            year_start: None,
            years_ended: 0,
            jumps_done: 0,
            jumped: 0.0,
            radiation: vec![Fluxes::default(); n],
            clear: vec![Fluxes::default(); n],
            tile_radiation: vec![TileRadiation::default(); n],
            surface: vec![SurfaceFlux::default(); n],
            moist: vec![MoistTendency::default(); n],
            heating: model.kinetic_heating(),
            prescribed: false,
            uptake: 0.0,
            dissipation: 0.0,
            dissipated: 0.0,
            absorbed: 0.0,
            outgoing: 0.0,
            evaporated: 0.0,
            precipitated: 0.0,
            runoff: 0.0,
            ocean_exchange: 0.0,
        };
        if physics.ice_sheet_age.is_finite() {
            let land = params.land;
            for i in 0..n {
                if physics.land_fraction[i] > 0.0
                    && physics.land_temperature[i] < land.ice_sheet_start_temperature
                {
                    physics.snow[i] = land.ice_sheet_start_pack;
                    physics.cover_age[i] = physics.ice_sheet_age;
                }
            }
        }
        // Water that starts below freezing is ice from the start.
        physics.ocean.rearrange_ice();
        // The first year is measured from here when it ends in a jump.
        if params.spin_up.ends_in_jump(1) {
            physics.year_start = Some(physics.ocean.snapshot());
        }
        // Clouds for the first radiation call; radiation (clear sky) for the
        // first convection's stratosphere temperature.
        physics.radiate(state);
        physics.update_moisture(state);
        physics.radiate(state);
        physics.update_surface_fluxes(state, None);
        physics
    }

    /// The physics with each cell's mean surface `height` (m), whose slope
    /// the wind climbs, and sub-grid spread of it `relief` (m), which rains
    /// and drags on the air crossing it.
    ///
    /// # Panics
    /// If either doesn't have one value per cell.
    pub fn with_relief(mut self, relief: Vec<f64>, height: &[f64]) -> Self {
        assert_eq!(relief.len(), self.relief.len());
        assert_eq!(height.len(), self.relief.len());
        let (mesh, geometry) = (self.ops.mesh(), self.ops.geometry());
        let mut gradient = vec![0.0; mesh.edges().len()];
        self.ops.gradient(height, &mut gradient);
        self.relief_slope = geometry.reconstruct_vectors(mesh, &gradient);
        self.relief = relief;
        self
    }

    /// The physics with the surface temperature held at its starting values
    /// (fixed sea-surface temperature, as models are tuned; e.g. [`ApeSst`]):
    /// the surface becomes an infinite heat reservoir whose uptake the
    /// budget counts as stored heat. The atmosphere settles in about a month
    /// instead of the surface's years.
    ///
    /// # Panics
    /// If there is land, which always evolves.
    pub fn with_prescribed_surface(mut self) -> Self {
        assert!(
            self.land_fraction.iter().all(|f| *f == 0.0),
            "only the ocean's temperature can be prescribed"
        );
        self.prescribed = true;
        self
    }

    /// Total heat of the ocean (J).
    fn ocean_heat(&self) -> f64 {
        let areas = &self.ops.geometry().cell_areas;
        let r2 = self.ops.radius() * self.ops.radius();
        (0..areas.len())
            .map(|i| (1.0 - self.land_fraction[i]) * (areas[i] * r2 * self.ocean.heat(i)))
            .sum()
    }

    /// A model year has ended: jump if it was measured, and start measuring
    /// the next if that ends in a jump.
    fn end_spin_up_year(&mut self) {
        let spin_up = self.params.spin_up;
        self.years_ended += 1;
        if let Some(start) = self.year_start.take() {
            debug_assert!(spin_up.ends_in_jump(self.years_ended));
            let before = self.ocean_heat();
            self.ocean
                .jump(&start, spin_up.extra_years, spin_up.limit, spin_up.layers);
            self.jumped += self.ocean_heat() - before;
            self.jumps_done += 1;
        }
        if spin_up.ends_in_jump(self.years_ended + 1) {
            self.year_start = Some(self.ocean.snapshot());
        }
    }

    /// The parameters.
    pub fn params(&self) -> &PhysicsParams {
        &self.params
    }

    /// Area-weighted mean surface temperature per cell (K).
    pub fn surface_temperature(&self) -> Vec<f64> {
        (0..self.land_fraction.len())
            .map(|i| {
                let water = self.ocean.mixed_temperature()[i];
                let sea = match self.ocean.ice_fraction(i) {
                    0.0 => water,
                    ice => (1.0 - ice) * water + ice * self.ocean.ice_temperature(i),
                };
                match self.land_fraction[i] {
                    0.0 => sea,
                    1.0 => self.land_temperature[i],
                    f => (1.0 - f) * sea + f * self.land_temperature[i],
                }
            })
            .collect()
    }

    /// Skin temperature of each cell's land tile (K); meaningless where the
    /// cell has no land.
    pub fn land_temperature(&self) -> &[f64] {
        &self.land_temperature
    }

    /// Temperature of each cell's soil store (K); meaningless where the cell
    /// has no land.
    pub fn soil_temperature(&self) -> &[f64] {
        &self.soil_temperature
    }

    /// Water in each cell's soil (kg per m² of land); meaningless where the
    /// cell has no land.
    pub fn soil_water(&self) -> &[f64] {
        &self.soil_water
    }

    /// Snow on each cell's land (kg per m² of land); meaningless where the
    /// cell has no land.
    pub fn snow(&self) -> &[f64] {
        &self.snow
    }

    /// Whether cell `i`'s land is an ice sheet: snow has covered it through
    /// a whole year, so no seasonal snow counts.
    fn is_ice_sheet(&self, i: usize) -> bool {
        self.cover_age[i] >= self.ice_sheet_age
    }

    /// What the radiation sees of cell `i`'s tiles.
    fn tiled_surface(&self, i: usize) -> TiledSurface {
        TiledSurface {
            fraction: self.land_fraction[i],
            ice: self.ocean.ice_fraction(i),
            albedo: [
                self.params.radiation.albedo,
                [self.params
                    .land
                    .surface_albedo(self.snow[i], self.is_ice_sheet(i)); 2],
                self.ocean.ice_albedo(i, self.params.radiation.albedo),
            ],
            temperature: [
                self.ocean.mixed_temperature()[i],
                self.land_temperature[i],
                self.ocean.ice_temperature(i),
            ],
        }
    }

    /// The last radiation call's fluxes per cell.
    pub fn fluxes(&self) -> &[Fluxes] {
        &self.radiation
    }

    /// The last moist physics per cell.
    pub fn moisture(&self) -> &[MoistTendency] {
        &self.moist
    }

    /// Model time (s).
    pub fn time(&self) -> f64 {
        self.time
    }

    fn setting(&self) -> Setting<'_> {
        Setting {
            params: &self.params.moisture,
            radiation: &self.params.radiation,
            gas: self.gas,
            levels: &self.levels,
            gravity: self.gravity,
        }
    }

    /// One column's view for the moist physics.
    fn moist_column(&self, state: &LayeredState, i: usize) -> MoistColumn {
        let ps = state.ps[i];
        let n = self.levels.count();
        let mut exner = [0.0; MAX_LAYERS];
        self.levels.layer_exner_into(self.gas, ps, &mut exner[..n]);
        let theta: Vec<f64> = (0..n).map(|k| state.theta_ps[k][i] / ps).collect();
        // The stratosphere sits near the skin temperature of the column's
        // outgoing longwave.
        let olr = self.radiation[i].outgoing.max(1.0);
        let humidity: Vec<f64> = (0..n).map(|k| state.q_ps[k][i] / ps).collect();
        // The well-mixed boundary layer reaches the surface.
        let surface_air = theta[n - 1] * self.gas.exner(ps);
        MoistColumn {
            ps,
            temperature: (0..n).map(|k| theta[k] * exner[k]).collect(),
            surface_humidity: humidity[n - 1].clamp(0.0, saturation_humidity(surface_air, ps)),
            humidity,
            surface_air,
            stratosphere: fms_math::pow(olr / (2.0 * STEFAN_BOLTZMANN), 0.25),
            orographic_ascent: self.surface[i].upslope,
        }
    }

    fn update_moisture(&mut self, state: &LayeredState) {
        let moist: Vec<MoistTendency> = (0..state.ps.len())
            .into_par_iter()
            .with_min_len(MIN_COLUMNS)
            .map(|i| moist_physics(&self.setting(), &self.moist_column(state, i)))
            .collect();
        self.moist = moist;
    }

    /// The surface fluxes for the coming step of `dt` seconds, first solving
    /// the land skin's balance; with no `dt` the skin stays where it is.
    fn update_surface_fluxes(&mut self, state: &LayeredState, dt: Option<f64>) {
        let (mesh, geometry) = (self.ops.mesh(), self.ops.geometry());
        let n = self.levels.count();
        let bottom = n - 1;
        let winds = geometry.reconstruct_vectors(mesh, &state.u[bottom]);
        let centers = mesh.cell_centers();
        let (gas, g) = (self.gas, self.gravity);
        let (rd, cp) = (gas.gas_constant, gas.heat_capacity);
        let sl = &self.params.surface_layer;
        let ds = self.levels.dsigma();
        let top = self.levels.interfaces()[bottom];
        let old = &self.surface;
        let surface: Vec<SurfaceFlux> = (0..state.ps.len())
            .into_par_iter()
            .with_min_len(MIN_COLUMNS)
            .map(|i| {
                let (e, nn) = east_north(centers[i], winds[i]);
                let ps = state.ps[i];
                let mut exner = [0.0; MAX_LAYERS];
                self.levels.layer_exner_into(gas, ps, &mut exner[..n]);
                let theta = state.theta_ps[bottom][i] / ps;
                let q = (state.q_ps[bottom][i] / ps).max(0.0);
                let surface_exner = gas.exner(ps);
                let air = theta * surface_exner;
                // Heights of the layer's mean level and top (hypsometric).
                let t = theta * exner[bottom];
                let scale = rd * t / g;
                let p_mean = REFERENCE_PRESSURE * fms_math::pow(exner[bottom], 1.0 / gas.kappa());
                let height = scale * fms_math::ln(ps / p_mean);
                let depth = -scale * fms_math::ln(top);
                let density = ps / (rd * air);
                let wind = (e * e + nn * nn).sqrt();
                let air_over = |ts: f64| SurfaceAir {
                    wind,
                    height,
                    theta,
                    humidity: q,
                    surface_temperature: ts,
                    ps,
                    surface_exner,
                    convective_velocity: old[i].convective_velocity,
                    gravity: g,
                };
                // Dew condenses whatever the soil holds.
                let supplied =
                    |deficit: f64, wetness: f64| if deficit > 0.0 { wetness } else { 1.0 };
                let tile = |ts: f64, surface: Surface, wetness: f64| {
                    let x = exchange(sl, surface, &air_over(ts));
                    let deficit = saturation_humidity(ts, ps) - q;
                    let momentum = x.momentum(density);
                    TileFlux {
                        sensible: density * cp * x.heat * x.speed * (ts - air),
                        evaporation: density
                            * x.moisture
                            * x.speed
                            * deficit
                            * supplied(deficit, wetness),
                        drag: g * momentum * x.wind_ratio / (ps * ds[bottom]),
                        wind: x.wind,
                        stress: momentum * x.wind,
                        stress_east: -momentum * x.wind_ratio * e,
                        stress_north: -momentum * x.wind_ratio * nn,
                    }
                };
                // A surface's turbulent heat loss at `t`, and its slope.
                // Exchange coefficients are frozen in the slope only.
                let losses = |t: f64, surface: Surface, wetness: f64, latent: f64| {
                    let x = exchange(sl, surface, &air_over(t));
                    let heat = density * cp * x.heat * x.speed;
                    let water = latent * density * x.moisture * x.speed;
                    let qs = saturation_humidity(t, ps);
                    let dqs = qs * LATENT_HEAT / (VAPOUR_GAS_CONSTANT * t * t);
                    let water = water * supplied(qs - q, wetness);
                    (heat * (t - air) + water * (qs - q), heat + water * dqs)
                };
                let fraction = self.land_fraction[i];
                let mut tiles = [TileFlux::default(); 3];
                if fraction < 1.0 {
                    tiles[OCEAN] = tile(self.ocean.mixed_temperature()[i], Surface::Ocean, 1.0);
                }
                let mut skin = self.land_temperature[i];
                if fraction > 0.0 {
                    let rough =
                        Surface::Rough(self.params.land.surface_roughness(self.is_ice_sheet(i)));
                    let snow = self.snow[i];
                    let wetness = self.params.land.wetness(self.soil_water[i] + snow);
                    if let Some(dt) = dt {
                        let absorbed =
                            self.tile_radiation[i].shortwave[LAND] + self.radiation[i].surface_down;
                        // Evaporating snow also pays for its melting.
                        let latent = LATENT_HEAT + if snow > 0.0 { FUSION_HEAT } else { 0.0 };
                        skin = solve_skin(
                            skin,
                            absorbed,
                            self.soil_temperature[i],
                            self.params.land.step_conductance(dt),
                            |t| losses(t, rough, wetness, latent),
                        );
                        if snow > 0.0 {
                            skin = skin.min(MELTING_POINT);
                        }
                    }
                    tiles[LAND] = tile(skin, rough, wetness);
                }
                let ice = self.ocean.ice_fraction(i);
                let mut ice_skin = self.ocean.ice_temperature(i);
                if ice > 0.0 {
                    let p = self.ocean.ice_params().expect("ice needs parameters");
                    let rough = Surface::Rough(p.roughness);
                    if let Some(dt) = dt {
                        let (conductance, interior) = self.ocean.ice_balance(i, dt);
                        let absorbed =
                            self.tile_radiation[i].shortwave[ICE] + self.radiation[i].surface_down;
                        let latent = LATENT_HEAT + FUSION_HEAT;
                        ice_skin = solve_skin(ice_skin, absorbed, interior, conductance, |t| {
                            losses(t, rough, 1.0, latent)
                        })
                        .min(p.melting);
                    }
                    tiles[ICE] = tile(ice_skin, rough, 1.0);
                }
                let sea = TileFlux::blend(tiles[OCEAN], tiles[ICE], ice);
                let mean = TileFlux::blend(sea, tiles[LAND], fraction);
                let hills = &self.params.orography;
                let hill_drag = (g * density * hills.drag * self.relief[i] / hills.hill_spacing
                    * wind
                    / (ps * ds[bottom]))
                    .min(hills.max_drag_rate);
                let (sensible, evaporation) = (mean.sensible, mean.evaporation);
                // Clouds' longwave cooling of the layer drives turbulence
                // from its top.
                let cloud_top_cooling =
                    self.clear[i].layers[bottom] - self.radiation[i].layers[bottom];
                let above = bottom - 1;
                let layer_top = boundary_layer_top(
                    &BoundaryColumn {
                        gas,
                        gravity: g,
                        lapse_rate: self.params.radiation.lapse_rate,
                        top,
                        depth,
                        ps,
                        theta,
                        humidity: q,
                        layer_humidity: state.q_ps[bottom][i] / ps,
                        air,
                        density,
                        theta_above: state.theta_ps[above][i] / ps,
                        humidity_above: state.q_ps[above][i] / ps,
                        exner_above: exner[above],
                        sensible,
                        evaporation,
                        cloud_top_cooling,
                    },
                    sl,
                );
                let (w_star, w_e, cumulus, density_top) = (
                    layer_top.convective_velocity,
                    layer_top.entrainment_velocity,
                    layer_top.cumulus,
                    layer_top.density,
                );
                SurfaceFlux {
                    tiles,
                    skin,
                    ice_skin,
                    air,
                    sensible,
                    evaporation,
                    drag: mean.drag + hill_drag,
                    upslope: hills.rain_efficiency
                        * (self.relief[i] / hills.full_relief).min(1.0)
                        * winds[i].dot(self.relief_slope[i]).max(0.0),
                    entrainment: density_top * w_e,
                    cumulus: if cumulus {
                        sl.cumulus * density_top * w_star
                    } else {
                        0.0
                    },
                    entrainment_velocity: w_e,
                    convective_velocity: w_star,
                    wind: mean.wind,
                    stress: mean.stress,
                    export: (density_top * w_e
                        + if cumulus {
                            sl.cumulus * density_top * w_star
                        } else {
                            0.0
                        })
                        * (state.q_ps[bottom][i] - state.q_ps[above][i])
                        / ps,
                    stress_east: mean.stress_east,
                    stress_north: mean.stress_north,
                }
            })
            .collect();
        self.surface = surface;
        for (i, flux) in self.surface.iter().enumerate() {
            if self.land_fraction[i] > 0.0 {
                self.land_temperature[i] = flux.skin;
            }
        }
    }

    /// Radiation for the interval starting now, with sunlight averaged over
    /// it, all-sky and clear-sky.
    fn radiate(&mut self, state: &LayeredState) {
        const SAMPLES: usize = 6;
        let centers = self.ops.mesh().cell_centers();
        let n = centers.len();
        let mut insolation = vec![0.0; n];
        let mut sample = vec![0.0; n];
        let h = self.params.radiation_interval / SAMPLES as f64;
        for j in 0..SAMPLES {
            self.planet
                .insolation_into(centers, self.time + (j as f64 + 0.5) * h, &mut sample);
            for (a, s) in insolation.iter_mut().zip(&sample) {
                *a += s / SAMPLES as f64;
            }
        }
        let (params, gas, gravity) = (&self.params, self.gas, self.gravity);
        let levels = &self.levels;
        let layers = levels.count();
        let clouds = &params.moisture.clouds;
        let moist = &self.moist;
        let this = &*self;
        let results: Vec<(Fluxes, Fluxes, TileRadiation)> = (0..n)
            .into_par_iter()
            .with_min_len(MIN_COLUMNS)
            .map(|i| {
                let ps = state.ps[i];
                let surface = this.tiled_surface(i);
                let radiation = RadiationParams {
                    albedo: surface.mean_albedo(),
                    ..params.radiation
                };
                let cover = moist[i].clouds;
                let theta: Vec<f64> = (0..layers).map(|k| state.theta_ps[k][i] / ps).collect();
                let q: Vec<f64> = (0..layers).map(|k| state.q_ps[k][i] / ps).collect();
                let sheets: [Cloud; 2] = std::array::from_fn(|c| Cloud {
                    sigma: clouds.sigma[c],
                    reflectance: cover[c] * clouds.albedo[c],
                    emissivity: cover[c] * clouds.emissivity[c],
                });
                let mut column = Column {
                    ps,
                    surface_temperature: surface.radiating_temperature(),
                    theta: &theta,
                    levels,
                    gravity,
                    insolation: insolation[i],
                    co2: params.co2,
                    humidity: Humidity::Layers(&q),
                    clouds: &sheets,
                };
                let all = column_fluxes(&radiation, gas, &column);
                column.clouds = &[Cloud::NONE; 2];
                let clear = column_fluxes(&radiation, gas, &column);
                let tiles =
                    surface.split(all.surface_shortwave, all.surface_visible, all.surface_up);
                (all, clear, tiles)
            })
            .collect();
        let n_results = results.len();
        let (mut all, mut clear, mut tiles) = (
            Vec::with_capacity(n_results),
            Vec::with_capacity(n_results),
            Vec::with_capacity(n_results),
        );
        for (a, c, t) in results {
            all.push(a);
            clear.push(c);
            tiles.push(t);
        }
        self.radiation = all;
        self.clear = clear;
        self.tile_radiation = tiles;
        self.next_radiation = self.time + self.params.radiation_interval;
    }
}

impl Forcing for ColumnPhysics<'_> {
    fn add_tendency(&mut self, state: &LayeredState, out: &mut LayeredState) {
        assert!(
            state.is_moist() && out.is_moist(),
            "column physics needs a moist state"
        );
        self.update_moisture(state);
        let (gas, gravity) = (self.gas, self.gravity);
        let levels = &self.levels;
        let n = levels.count();
        let lowest = n - 1;
        let ds = levels.dsigma();
        let mut exner = vec![0.0; n * state.ps.len()];
        exner
            .par_chunks_mut(n)
            .with_min_len(MIN_COLUMNS)
            .zip(state.ps.par_iter())
            .for_each(|(e, &ps)| levels.layer_exner_into(gas, ps, e));
        let (radiation, surface, moist) = (&self.radiation, &self.surface, &self.moist);
        let (tiles, land, skin) = (
            &self.tile_radiation,
            &self.land_fraction,
            &self.land_temperature,
        );
        let ocean = &self.ocean;
        let above = lowest - 1;
        let cp = gas.heat_capacity;
        let top = levels.interfaces()[lowest];
        out.theta_ps
            .par_iter_mut()
            .zip(out.q_ps.par_iter_mut())
            .enumerate()
            .for_each(|(k, (dt, dq))| {
                let exner = &exner;
                dt.par_iter_mut()
                    .with_min_len(MIN_CHUNK)
                    .zip(dq.par_iter_mut())
                    .enumerate()
                    .for_each(|(i, (dt, dq))| {
                        let ps = state.ps[i];
                        let e = exner[i * n + k];
                        let m = &moist[i];
                        let s = &surface[i];
                        let mut heating = radiation[i].layers[k];
                        let mut water = ps * m.humidity[k];
                        if k == lowest || k == above {
                            // Entrainment: dry static energy `c_p Π θ + gz`
                            // and water mixed across the boundary layer's top.
                            let theta = |k: usize| state.theta_ps[k][i] / ps;
                            let q = |k: usize| state.q_ps[k][i] / ps;
                            let exner_top = gas.exner(top * ps);
                            let heat =
                                s.entrainment * cp * exner_top * (theta(above) - theta(lowest));
                            let moisture =
                                gravity * (s.entrainment + s.cumulus) * (q(above) - q(lowest));
                            if k == lowest {
                                heating += s.sensible + heat;
                                // What the land emits beyond what the last
                                // radiation call counted.
                                if land[i] > 0.0 {
                                    heating += land[i]
                                        * (black_body(skin[i]) - tiles[i].longwave_up[LAND]);
                                }
                                let ice = (1.0 - land[i]) * ocean.ice_fraction(i);
                                if ice > 0.0 {
                                    heating +=
                                        ice * (black_body(s.ice_skin) - tiles[i].longwave_up[ICE]);
                                }
                                water += gravity * s.evaporation / ds[k] + moisture / ds[k];
                            } else {
                                heating -= heat;
                                water -= moisture / ds[k];
                            }
                        }
                        *dt += theta_tendency(gas, gravity, ds[k], e, heating)
                            + ps * m.temperature[k] / e;
                        *dq += water;
                    });
            });
        // Surface drag on the boundary layer's wind and entrainment of the
        // layer above's momentum; the kinetic energy they remove heats the
        // boundary layer.
        let edges = self.ops.mesh().edges();
        let (u_low, u_up) = (&state.u[lowest], &state.u[above]);
        let mut du_low = vec![0.0; edges.len()];
        let mut du_up = vec![0.0; edges.len()];
        du_low
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(du_up.par_iter_mut())
            .enumerate()
            .for_each(|(e, (dl, du))| {
                let [a, b] = edges[e].cells.map(|c| c as usize);
                let drag = 0.5 * (surface[a].drag + surface[b].drag);
                let exchange = |s: &SurfaceFlux| s.entrainment + s.cumulus;
                let mixing = 0.5 * (exchange(&surface[a]) + exchange(&surface[b])) * gravity
                    / (0.5 * (state.ps[a] + state.ps[b]));
                let shear = u_up[e] - u_low[e];
                *dl = -drag * u_low[e] + mixing * shear / ds[lowest];
                *du = -mixing * shear / ds[above];
            });
        let heating = &mut self.heating;
        heating.prepare(&state.ps);
        let ops = &self.ops;
        heating.add(
            ops,
            &state.ps,
            lowest,
            lowest,
            u_low,
            &du_low,
            &mut out.theta_ps[lowest],
        );
        heating.add(
            ops,
            &state.ps,
            above,
            lowest,
            u_up,
            &du_up,
            &mut out.theta_ps[lowest],
        );
        self.dissipation = -(heating.rate(ops, &state.ps, lowest, gravity, u_low, &du_low)
            + heating.rate(ops, &state.ps, above, gravity, u_up, &du_up));
        for (k, du) in [(lowest, &du_low), (above, &du_up)] {
            out.u[k]
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .zip(du.par_iter())
                .for_each(|(o, &d)| *o += d);
        }
    }

    fn advance(&mut self, state: &LayeredState, dt: f64) {
        // The surface and the budgets use the same stored fluxes and
        // precipitation the atmosphere was given.
        let areas = &self.ops.geometry().cell_areas;
        let r2 = self.ops.radius() * self.ops.radius();
        let land_capacity = self.params.land_heat_capacity();
        let (radiation, surface, moist) = (&self.radiation, &self.surface, &self.moist);
        let (tiles, land) = (&self.tile_radiation, &self.land_fraction);
        let skin = &self.land_temperature;
        let net = |i: usize, tile: usize| {
            let flux = surface[i].tiles[tile];
            let emitted = match tile {
                LAND => black_body(skin[i]),
                ICE => black_body(surface[i].ice_skin),
                _ => tiles[i].longwave_up[tile],
            };
            // Snow landing on open water melts at once; on land and ice
            // it's a pack. Ice that evaporates pays for its thawing first.
            let fusion = match tile {
                OCEAN => FUSION_HEAT * moist[i].snow,
                ICE => FUSION_HEAT * flux.evaporation.max(0.0),
                _ => 0.0,
            };
            tiles[i].shortwave[tile] + radiation[i].surface_down
                - emitted
                - flux.sensible
                - LATENT_HEAT * flux.evaporation
                - fusion
        };
        if self.prescribed {
            let uptake: f64 = (0..areas.len()).map(|i| areas[i] * net(i, OCEAN)).sum();
            self.uptake += dt * r2 * uptake;
        } else {
            self.ocean.add_surface_heat(dt, |i| net(i, OCEAN));
            self.ocean.step_ice(dt, |i| IceForcing {
                net: net(i, ICE),
                skin: surface[i].ice_skin,
                stress: surface[i].tiles[OCEAN].stress,
                sublimation: FUSION_HEAT * surface[i].tiles[ICE].evaporation.max(0.0),
                snowfall: moist[i].snow,
            });
            // The open water's stress on the air, turned round.
            self.ocean.set_wind_stress(|i| {
                let water = &surface[i].tiles[OCEAN];
                [-water.stress_east, -water.stress_north]
            });
            self.ocean.transport(dt);
        }
        let (mut absorbed, mut outgoing, mut evaporated, mut precipitated) = (0.0, 0.0, 0.0, 0.0);
        let (mut runoff, mut ocean_exchange) = (0.0, 0.0);
        for i in 0..areas.len() {
            let a = areas[i];
            let (f, falling) = (land[i], moist[i].rain + moist[i].snow);
            if f < 1.0 {
                let ice = self.ocean.ice_fraction(i);
                let evaporation = if ice > 0.0 {
                    (1.0 - ice) * surface[i].tiles[OCEAN].evaporation
                        + ice * surface[i].tiles[ICE].evaporation
                } else {
                    surface[i].tiles[OCEAN].evaporation
                };
                ocean_exchange += a * (1.0 - f) * (evaporation - falling);
            }
            if f > 0.0 {
                let evaporation = surface[i].tiles[LAND].evaporation;
                let pack = self.snow[i] + dt * moist[i].snow;
                // Evaporation takes the snow first, which costs its fusion.
                let sublimated = (dt * evaporation).clamp(0.0, pack);
                let mut heat = net(i, LAND) - FUSION_HEAT * sublimated / dt;
                // A skin held at freezing: what the balance leaves after
                // conduction into the soil melts the snow.
                let melted = if skin[i] == MELTING_POINT {
                    let conducted = self.params.land.step_conductance(dt)
                        * (skin[i] - self.soil_temperature[i]);
                    snow_melted(pack - sublimated, heat - conducted, dt)
                } else {
                    0.0
                };
                heat -= FUSION_HEAT * melted / dt;
                self.snow[i] = pack - sublimated - melted;
                self.cover_age[i] = self
                    .params
                    .land
                    .cover_age(self.cover_age[i], self.snow[i], dt);
                self.soil_temperature[i] += dt * heat / land_capacity;
                let gain = dt * (moist[i].rain - evaporation) + sublimated + melted;
                let (water, overflow) = self.params.land.fill(self.soil_water[i], gain);
                self.soil_water[i] = water;
                runoff += a * f * overflow;
            }
            absorbed += a * radiation[i].absorbed;
            outgoing += a * radiation[i].outgoing;
            evaporated += a * surface[i].evaporation;
            precipitated += a * (moist[i].rain + moist[i].snow);
        }
        self.absorbed += dt * r2 * absorbed;
        self.outgoing += dt * r2 * outgoing;
        self.evaporated += dt * r2 * evaporated;
        self.precipitated += dt * r2 * precipitated;
        self.runoff += r2 * runoff;
        self.ocean_exchange += dt * r2 * ocean_exchange;
        self.dissipated += dt * self.dissipation;
        self.time += dt;
        if self.jumps_done < self.params.spin_up.jumps
            && !self.prescribed
            && self.time >= self.next_year - 1e-6
        {
            self.next_year += self.planet.year_seconds();
            self.end_spin_up_year();
        }
        let radiating = self.time >= self.next_radiation - 1e-6;
        if radiating && !self.prescribed {
            // Before the fluxes, so they and the next sunlight see the same
            // ice.
            self.ocean.rearrange_ice();
        }
        let solve = (!self.prescribed).then_some(dt);
        self.update_surface_fluxes(state, solve);
        if radiating {
            self.radiate(state);
        }
    }

    fn budget(&self) -> Option<ForcingBudget> {
        let areas = &self.ops.geometry().cell_areas;
        let r2 = self.ops.radius() * self.ops.radius();
        let stored: f64 = if self.prescribed {
            self.uptake
        } else {
            let land = self.params.land_heat_capacity();
            (0..areas.len())
                .map(|i| {
                    let (a, f) = (areas[i], self.land_fraction[i]);
                    f * (a * r2 * (land * self.soil_temperature[i] - FUSION_HEAT * self.snow[i]))
                })
                .sum::<f64>()
                + self.ocean_heat()
                - self.jumped
        };
        Some(ForcingBudget {
            absorbed: self.absorbed,
            outgoing: self.outgoing,
            stored,
            evaporated: self.evaporated,
            precipitated: self.precipitated,
            dissipated: self.dissipated,
            soil_water: (0..areas.len())
                .map(|i| areas[i] * r2 * self.land_fraction[i] * self.soil_water[i])
                .sum(),
            snow: (0..areas.len())
                .map(|i| areas[i] * r2 * self.land_fraction[i] * self.snow[i])
                .sum(),
            runoff: self.runoff,
            ocean_exchange: self.ocean_exchange,
        })
    }

    fn cell_fields(&self) -> Vec<CellField> {
        let rad = |f: fn(&Fluxes) -> f64| self.radiation.iter().map(f).collect();
        let lowest = self.levels.count() - 1;
        let moist = |f: fn(&MoistTendency) -> f64| self.moist.iter().map(f).collect();
        let field = |name, description, units, values| CellField {
            name,
            description,
            units,
            values,
        };
        vec![
            field(
                "land.soil_water",
                "Soil water, per unit area of the whole cell",
                "kg m-2",
                self.soil_water
                    .iter()
                    .zip(&self.land_fraction)
                    .map(|(w, f)| w * f)
                    .collect(),
            ),
            field(
                "land.snow",
                "Snow on land, per unit area of the whole cell",
                "kg m-2",
                self.snow
                    .iter()
                    .zip(&self.land_fraction)
                    .map(|(s, f)| s * f)
                    .collect(),
            ),
            field(
                "land.ice_sheet",
                "Share of the whole cell that is ice sheet",
                "1",
                (0..self.land_fraction.len())
                    .map(|i| {
                        if self.is_ice_sheet(i) {
                            self.land_fraction[i]
                        } else {
                            0.0
                        }
                    })
                    .collect(),
            ),
            field(
                "ocean.mixed_temperature",
                "Mixed-layer temperature (its start where a cell has no ocean)",
                "K",
                self.ocean.mixed_temperature().to_vec(),
            ),
            field(
                "ocean.ice_fraction",
                "Share of the ocean under sea ice",
                "1",
                (0..self.land_fraction.len())
                    .map(|i| self.ocean.ice_fraction(i))
                    .collect(),
            ),
            field(
                "ocean.ice_thickness",
                "Mean thickness of the sea ice where there is any",
                "m",
                (0..self.land_fraction.len())
                    .map(|i| self.ocean.ice_thickness(i))
                    .collect(),
            ),
            field(
                "ocean.snow_on_ice",
                "Snow on the sea ice, per unit area of ice",
                "kg m-2",
                (0..self.land_fraction.len())
                    .map(|i| self.ocean.snow_on_ice(i))
                    .collect(),
            ),
            field(
                "ocean.ice_temperature",
                "Surface temperature of the sea ice or its snow (freezing where there is none)",
                "K",
                (0..self.land_fraction.len())
                    .map(|i| self.ocean.ice_temperature(i))
                    .collect(),
            ),
            field(
                "ocean.deep_temperature",
                "Deep-layer temperature (its start where a cell has no ocean)",
                "K",
                self.ocean.deep_temperature().to_vec(),
            ),
            field(
                "ocean.ekman_heating",
                "Heating of the ocean's two layers by Ekman transport, per unit area of ocean",
                "W m-2",
                (0..self.land_fraction.len())
                    .map(|i| self.ocean.ekman_heating(i))
                    .collect(),
            ),
            field(
                "ocean.ekman_mixed_heating",
                "Heating of the mixed layer by Ekman transport, per unit area of ocean",
                "W m-2",
                (0..self.land_fraction.len())
                    .map(|i| self.ocean.ekman_mixed_heating(i))
                    .collect(),
            ),
            field(
                "ocean.eddy_heating",
                "Heating of the ocean's two layers by eddy transport, per unit area of ocean",
                "W m-2",
                (0..self.land_fraction.len())
                    .map(|i| self.ocean.eddy_heating(i))
                    .collect(),
            ),
            field(
                "ocean.eddy_mixed_heating",
                "Heating of the mixed layer by eddy transport, per unit area of ocean",
                "W m-2",
                (0..self.land_fraction.len())
                    .map(|i| self.ocean.eddy_mixed_heating(i))
                    .collect(),
            ),
            field(
                "ocean.diffusion_heating",
                "Heating of the ocean's two layers by horizontal diffusion, per unit area of ocean",
                "W m-2",
                (0..self.land_fraction.len())
                    .map(|i| self.ocean.diffusion_heating(i))
                    .collect(),
            ),
            field(
                "ocean.ekman_heat_gained",
                "Heat Ekman transport has brought into the ocean's two layers since the start, \
                 per unit area of ocean; its change over a span is the mean heating, from every step",
                "J m-2",
                (0..self.land_fraction.len())
                    .map(|i| self.ocean.transported_heat(i)[0])
                    .collect(),
            ),
            field(
                "ocean.eddy_heat_gained",
                "Heat eddy transport has brought into the ocean's two layers since the start, \
                 per unit area of ocean",
                "J m-2",
                (0..self.land_fraction.len())
                    .map(|i| self.ocean.transported_heat(i)[1])
                    .collect(),
            ),
            field(
                "ocean.diffusion_heat_gained",
                "Heat horizontal diffusion has brought into the ocean's two layers since the start, \
                 per unit area of ocean",
                "J m-2",
                (0..self.land_fraction.len())
                    .map(|i| self.ocean.transported_heat(i)[2])
                    .collect(),
            ),
            field(
                "ocean.upwelling",
                "Deep water welling up into the mixed layer (negative where it sinks), per unit area of ocean",
                "kg m-2 s-1",
                (0..self.land_fraction.len())
                    .map(|i| self.ocean.upwelling(i))
                    .collect(),
            ),
            field(
                "ocean.basin",
                "Ocean basin of the wind-driven gyres, numbered from 1 by area (0 where the cell is \
                 not at least half ocean or has no vertex inside a basin)",
                "1",
                (0..self.land_fraction.len())
                    .map(|i| {
                        self.ocean
                            .basins()
                            .cell_basin(i as u32)
                            .map_or(0.0, |b| b as f64 + 1.0)
                    })
                    .collect(),
            ),
            field(
                "rad.surface_temperature",
                "Surface temperature",
                "K",
                self.surface_temperature(),
            ),
            field(
                "rad.outgoing_longwave",
                "Outgoing longwave radiation",
                "W m-2",
                rad(|f| f.outgoing),
            ),
            field(
                "rad.absorbed_shortwave",
                "Absorbed sunlight",
                "W m-2",
                rad(|f| f.absorbed),
            ),
            field(
                "rad.surface_shortwave",
                "Sunlight absorbed by the surface",
                "W m-2",
                rad(|f| f.surface_shortwave),
            ),
            field(
                "rad.surface_longwave",
                "Net longwave loss of the surface",
                "W m-2",
                rad(|f| f.surface_up - f.surface_down),
            ),
            field(
                "rad.clear_surface_longwave",
                "Net longwave loss of the surface under clear skies",
                "W m-2",
                self.clear
                    .iter()
                    .map(|f| f.surface_up - f.surface_down)
                    .collect(),
            ),
            field(
                "rad.shortwave_cloud_effect",
                "Shortwave cloud radiative effect",
                "W m-2",
                self.radiation
                    .iter()
                    .zip(&self.clear)
                    .map(|(a, c)| a.absorbed - c.absorbed)
                    .collect(),
            ),
            field(
                "rad.longwave_cloud_effect",
                "Longwave cloud radiative effect",
                "W m-2",
                self.radiation
                    .iter()
                    .zip(&self.clear)
                    .map(|(a, c)| c.outgoing - a.outgoing)
                    .collect(),
            ),
            field(
                "atm.water_vapour",
                "Water-vapour column",
                "kg m-2",
                self.radiation.iter().map(Fluxes::total_water).collect(),
            ),
            field(
                "atm.precipitation",
                "Precipitation",
                "mm day-1",
                moist(|m| (m.rain + m.snow) * DAY),
            ),
            field(
                "atm.orographic_precipitation",
                "Rain and snow from air forced over sub-grid hills",
                "mm day-1",
                moist(|m| m.orographic * DAY),
            ),
            field(
                "atm.snowfall",
                "Snowfall (water equivalent)",
                "mm day-1",
                moist(|m| m.snow * DAY),
            ),
            field(
                "atm.convective_precipitation",
                "Convective precipitation",
                "mm day-1",
                moist(|m| m.convective * DAY),
            ),
            field(
                "atm.evaporation",
                "Evaporation",
                "mm day-1",
                self.surface.iter().map(|s| s.evaporation * DAY).collect(),
            ),
            field(
                "atm.latent_heat_flux",
                "Surface latent heat flux",
                "W m-2",
                self.surface
                    .iter()
                    .map(|s| s.evaporation * LATENT_HEAT)
                    .collect(),
            ),
            field(
                "atm.sensible_heat_flux",
                "Surface sensible heat flux",
                "W m-2",
                self.surface.iter().map(|s| s.sensible).collect(),
            ),
            field(
                "atm.boundary_layer_moist_heating",
                "Moist-physics heating of the boundary layer",
                "K day-1",
                self.moist
                    .iter()
                    .map(|m| m.temperature[lowest] * DAY)
                    .collect(),
            ),
            field(
                "atm.surface_wind",
                "10 m wind speed",
                "m s-1",
                self.surface.iter().map(|s| s.wind).collect(),
            ),
            field(
                "atm.surface_stress",
                "Surface wind stress",
                "N m-2",
                self.surface.iter().map(|s| s.stress).collect(),
            ),
            field(
                "atm.shallow_cumulus",
                "Shallow cumulus mass flux out of the boundary layer",
                "g m-2 s-1",
                self.surface.iter().map(|s| 1e3 * s.cumulus).collect(),
            ),
            field(
                "moist.condensation_boundary",
                "Large-scale condensation in the boundary layer",
                "mm day-1",
                self.moist
                    .iter()
                    .map(|m| m.condensed[lowest] * DAY)
                    .collect(),
            ),
            field(
                "moist.condensation_above",
                "Large-scale condensation in the layer above the boundary layer",
                "mm day-1",
                self.moist
                    .iter()
                    .map(|m| m.condensed[lowest - 1] * DAY)
                    .collect(),
            ),
            field(
                "moist.condensation_upper",
                "Large-scale condensation in the higher layers",
                "mm day-1",
                self.moist
                    .iter()
                    .map(|m| m.condensed[..lowest - 1].iter().sum::<f64>() * DAY)
                    .collect(),
            ),
            field(
                "moist.reevaporation_boundary",
                "Rain re-evaporated in the boundary layer",
                "mm day-1",
                self.moist
                    .iter()
                    .map(|m| m.reevaporated[lowest] * DAY)
                    .collect(),
            ),
            field(
                "moist.reevaporation_above",
                "Rain re-evaporated above the boundary layer",
                "mm day-1",
                self.moist
                    .iter()
                    .map(|m| m.reevaporated[..lowest].iter().sum::<f64>() * DAY)
                    .collect(),
            ),
            field(
                "moist.convective_boundary",
                "Water convection adds to the boundary layer",
                "mm day-1",
                self.moist
                    .iter()
                    .map(|m| m.convective_boundary * DAY)
                    .collect(),
            ),
            field(
                "atm.surface_air_temperature",
                "The bottom layer's air brought dry-adiabatically to the surface",
                "K",
                self.surface.iter().map(|s| s.air).collect(),
            ),
            field(
                "atm.boundary_layer_export",
                "Water carried out of the boundary layer by entrainment and cumulus",
                "mm day-1",
                self.surface.iter().map(|s| s.export * DAY).collect(),
            ),
            field(
                "atm.surface_stress_east",
                "Eastward surface stress on the air",
                "N m-2",
                self.surface.iter().map(|s| s.stress_east).collect(),
            ),
            field(
                "atm.surface_stress_north",
                "Northward surface stress on the air",
                "N m-2",
                self.surface.iter().map(|s| s.stress_north).collect(),
            ),
            field(
                "rad.free_troposphere_heating",
                "Radiative heating of the layers above the boundary layer",
                "W m-2",
                self.radiation
                    .iter()
                    .map(|f| f.layers[..lowest].iter().sum())
                    .collect(),
            ),
            field(
                "rad.clear_free_troposphere_heating",
                "Clear-sky radiative heating of the layers above the boundary layer",
                "W m-2",
                self.clear
                    .iter()
                    .map(|f| f.layers[..lowest].iter().sum())
                    .collect(),
            ),
            field(
                "atm.convective_velocity",
                "Boundary-layer convective velocity w*",
                "m s-1",
                self.surface.iter().map(|s| s.convective_velocity).collect(),
            ),
            field(
                "atm.entrainment",
                "Entrainment velocity at the boundary-layer top",
                "mm s-1",
                self.surface
                    .iter()
                    .map(|s| 1e3 * s.entrainment_velocity)
                    .collect(),
            ),
            field(
                "atm.high_cloud",
                "High cloud cover",
                "1",
                moist(|m| m.clouds[0]),
            ),
            field(
                "atm.low_cloud",
                "Low cloud cover",
                "1",
                moist(|m| m.clouds[1]),
            ),
            field(
                "atm.cloud_cover",
                "Total cloud cover (random overlap)",
                "1",
                moist(|m| 1.0 - (1.0 - m.clouds[0]) * (1.0 - m.clouds[1])),
            ),
            field(
                "atm.relative_humidity_upper",
                "Top-layer relative humidity",
                "1",
                moist(|m| m.relative_humidity[0]),
            ),
            field(
                "atm.relative_humidity_lower",
                "Lowest-layer relative humidity",
                "1",
                self.moist
                    .iter()
                    .map(|m| m.relative_humidity[lowest])
                    .collect(),
            ),
        ]
    }
}

/// A moist start for the 3-D model on `planet`: everywhere the single
/// column's equilibrium (with `params`' CO₂) at rest, with its fixed-profile
/// humidity, plus 0.1 K of seeded noise in θ to break the symmetry. Returns
/// the state and the surface temperature (K).
///
/// The single column has two equal-mass layers; other layers take the mean
/// temperature of its profile over their range (each layer's lapse rate
/// through its mean, isothermal above the upper layer's mean level).
///
/// Over terrain (the model's surface geopotential) the same temperatures
/// hold on pressure surfaces, so the column is hydrostatic and at rest: the
/// planet's surface pressure is the sea-level value, and each cell's follows
/// from its height.
pub fn initial_state(
    model: &LayeredModel,
    planet: &Planet,
    params: &PhysicsParams,
) -> (LayeredState, f64) {
    let mut column = SingleColumn::earth(planet.mean_flux());
    column.params.co2 = params.co2;
    column.ps = planet.params().surface_pressure;
    column.gravity = model.gravity();
    let eq = column.equilibrium();
    let (gas, levels, ps) = (model.gas(), model.levels(), column.ps);
    let exner = levels.layer_exner(gas, ps);
    let two = Levels::two(column.sigma);
    let temperature: Vec<f64> = if *levels == two {
        eq.temperature.to_vec()
    } else {
        let profile = Profile::new(
            &params.radiation,
            gas,
            &two,
            ps,
            column.gravity,
            &eq.temperature,
        );
        let cap = two.mean_pressures(gas, ps)[0];
        (0..levels.count())
            .map(|k| {
                let (a, b) = levels.bounds(k);
                GAUSS
                    .iter()
                    .map(|&(x, w)| {
                        let s = 0.5 * (a + b) + 0.5 * (b - a) * x;
                        let t = if s * ps < cap {
                            eq.temperature[0]
                        } else {
                            profile.at(s)
                        };
                        0.5 * w * t
                    })
                    .sum()
            })
            .collect()
    };
    let theta0: Vec<f64> = temperature.iter().zip(&exner).map(|(t, e)| t / e).collect();
    let q = fixed_humidity(params, gas, levels, model.gravity(), ps, &theta0);
    let mesh = model.operators().mesh();
    let (nc, ne) = (mesh.n_cells(), mesh.n_edges());
    let n = levels.count();
    let reference = RestProfile::new(
        gas,
        ps,
        theta0[n - 1] * gas.exner(ps),
        &levels.mean_pressures(gas, ps),
        &temperature,
    );
    // Surface pressure, θ and humidity of each column at rest.
    let columns: Vec<(f64, Vec<f64>, Vec<f64>)> = model
        .surface_geopotential()
        .par_iter()
        .with_min_len(MIN_COLUMNS)
        .map(|&phi| {
            if phi == 0.0 {
                return (ps, theta0.clone(), q.clone());
            }
            let p = reference.pressure_at(phi);
            let theta: Vec<f64> = levels
                .mean_pressures(gas, p)
                .iter()
                .zip(levels.layer_exner(gas, p))
                .map(|(&pk, e)| reference.temperature(pk) / e)
                .collect();
            let q = fixed_humidity(params, gas, levels, model.gravity(), p, &theta);
            (p, theta, q)
        })
        .collect();
    let mut x = 0x2545_F491_4F6C_DD1Du64;
    let theta: Vec<Vec<f64>> = (0..n)
        .map(|k| {
            columns
                .iter()
                .map(|c| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    c.1[k] + 0.1 * ((x >> 11) as f64 / (1u64 << 52) as f64 - 1.0)
                })
                .collect()
        })
        .collect();
    let humidity = (0..n)
        .map(|k| columns.iter().map(|c| c.2[k]).collect())
        .collect();
    let surface_pressure = columns.iter().map(|c| c.0).collect();
    debug_assert_eq!(columns.len(), nc);
    let state =
        LayeredState::new(surface_pressure, vec![vec![0.0; ne]; n], theta).with_humidity(humidity);
    (state, eq.surface_temperature)
}

/// Warmest minus coldest starting temperature per unit of the annual-mean
/// insolation's relative contrast (K): Earth's seas, 29 K across 0.72.
const START_CONTRAST: f64 = 40.0;

/// [`initial_surface_temperature`], warmer where the year's mean sunlight is
/// above the global mean and colder where it is below, with the global mean
/// unchanged, so the ocean starts nearer its balance than from one value.
pub fn initial_surface_temperature_by_latitude(
    model: &LayeredModel,
    planet: &Planet,
    params: &PhysicsParams,
    sea_level: f64,
) -> Vec<f64> {
    const SAMPLES: usize = 48;
    let year = planet.year_seconds();
    let ops = model.operators();
    let sunlight: Vec<f64> = ops
        .mesh()
        .cell_centers()
        .iter()
        .map(|&p| {
            let lat = fms_grid::coords::to_lat_lon(p).0;
            (0..SAMPLES)
                .map(|n| planet.daily_mean_insolation(lat, year * n as f64 / SAMPLES as f64))
                .sum::<f64>()
                / SAMPLES as f64
        })
        .collect();
    let areas = &ops.geometry().cell_areas;
    let mean =
        sunlight.iter().zip(areas).map(|(s, a)| s * a).sum::<f64>() / areas.iter().sum::<f64>();
    initial_surface_temperature(model, params, sea_level)
        .into_iter()
        .zip(&sunlight)
        .map(|(t, s)| t + START_CONTRAST * (s / mean - 1.0))
        .collect()
}

/// Starting surface temperature per cell (K): `sea_level` (K) at sea level,
/// falling with height at the radiation's lapse rate.
pub fn initial_surface_temperature(
    model: &LayeredModel,
    params: &PhysicsParams,
    sea_level: f64,
) -> Vec<f64> {
    let lapse = params.radiation.lapse_rate / model.gravity();
    model
        .surface_geopotential()
        .iter()
        .map(|phi| sea_level - lapse * phi)
        .collect()
}

/// A horizontally uniform temperature profile: linear in ln p between the
/// surface at sea level and each layer's mean level, isothermal above the
/// top one, and the hydrostatic pressure at a given height in it.
struct RestProfile {
    gas_constant: f64,
    /// `(ln p, T)` from sea level upward.
    nodes: Vec<(f64, f64)>,
}

impl RestProfile {
    /// The profile with temperature `surface` (K) at sea-level pressure `ps`
    /// (Pa), and `temperature[k]` at pressures `pressure[k]`, layers ordered
    /// top down.
    fn new(gas: Gas, ps: f64, surface: f64, pressure: &[f64], temperature: &[f64]) -> Self {
        let mut nodes = vec![(fms_math::ln(ps), surface)];
        nodes.extend(
            pressure
                .iter()
                .zip(temperature)
                .rev()
                .map(|(&p, &t)| (fms_math::ln(p), t)),
        );
        Self {
            gas_constant: gas.gas_constant,
            nodes,
        }
    }

    /// Temperature (K) at pressure `p` (Pa).
    fn temperature(&self, p: f64) -> f64 {
        let x = fms_math::ln(p);
        for pair in self.nodes.windows(2) {
            let ((xa, ta), (xb, tb)) = (pair[0], pair[1]);
            if x >= xb {
                return ta + (tb - ta) * (xa - x) / (xa - xb);
            }
        }
        self.nodes[self.nodes.len() - 1].1
    }

    /// Pressure (Pa) at geopotential `phi` (m² s⁻², at or above sea level).
    ///
    /// Within a segment `T = T_a + s d` with `d = ln p_a − ln p`, so the
    /// hypsometric equation `c = (Φ − Φ_a)/R = T_a d + s d²/2` has the
    /// root `d = 2c / (T_a + √(T_a² + 2 s c))`, regular at `s = 0`.
    fn pressure_at(&self, phi: f64) -> f64 {
        let mut c = phi / self.gas_constant;
        for pair in self.nodes.windows(2) {
            let ((xa, ta), (xb, tb)) = (pair[0], pair[1]);
            let span = xa - xb;
            let thickness = 0.5 * (ta + tb) * span;
            if c <= thickness {
                let s = (tb - ta) / span;
                let d = 2.0 * c / (ta + (ta * ta + 2.0 * s * c).sqrt());
                return fms_math::exp(xa - d);
            }
            c -= thickness;
        }
        let (x, t) = self.nodes[self.nodes.len() - 1];
        fms_math::exp(x - c / t)
    }
}

/// Gauss–Legendre nodes on [−1, 1].
const GAUSS: [(f64, f64); 4] = [
    (-0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
    (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
];

/// Initial humidity (kg kg⁻¹, per layer) for a column at fixed relative
/// humidity: the Manabe–Wetherald profile of the single column.
pub fn fixed_humidity(
    params: &PhysicsParams,
    gas: Gas,
    levels: &Levels,
    gravity: f64,
    ps: f64,
    theta: &[f64],
) -> Vec<f64> {
    let f = column_fluxes(
        &params.radiation,
        gas,
        &Column {
            ps,
            surface_temperature: 288.0,
            theta,
            levels,
            gravity,
            insolation: 0.0,
            co2: params.co2,
            humidity: Humidity::Fixed,
            clouds: &[Cloud::NONE; 2],
        },
    );
    let ds = levels.dsigma();
    (0..levels.count())
        .map(|k| f.water[k] * gravity / (ds[k] * ps))
        .collect()
}

/// A single-column radiative–convective model: the same physics as
/// [`ColumnPhysics`] without dynamics, under a fixed insolation. For
/// calibration and quick experiments (CO₂, pressure, gravity, sunlight).
#[derive(Clone, Copy, Debug)]
pub struct SingleColumn {
    /// Parameters (the slab depth only sets how fast it converges).
    pub params: PhysicsParams,
    /// Gas constants.
    pub gas: Gas,
    /// Interface σ.
    pub sigma: f64,
    /// Gravity (m s⁻²).
    pub gravity: f64,
    /// Surface pressure (Pa).
    pub ps: f64,
    /// Mean insolation (W m⁻²).
    pub insolation: f64,
    /// Fixed wind for the surface fluxes (m s⁻¹).
    pub wind: f64,
    /// Layer potential temperatures (K).
    pub theta: [f64; 2],
    /// Surface temperature (K).
    pub surface_temperature: f64,
}

/// A single column's equilibrium.
#[derive(Clone, Debug, PartialEq)]
pub struct Equilibrium {
    /// Radiative fluxes.
    pub fluxes: Fluxes,
    /// Surface temperature (K).
    pub surface_temperature: f64,
    /// Layer temperatures, upper first (K).
    pub temperature: [f64; 2],
    /// Sensible and latent heat flux (W m⁻²).
    pub surface_fluxes: (f64, f64),
    /// Steps taken.
    pub steps: u64,
}

impl SingleColumn {
    /// Albedo of the single column's fixed cloud deck, calibrated so the
    /// column reaches Earth's global-mean climate (with the clear air's
    /// Rayleigh scattering, [`RadiationParams::rayleigh`], on top). One
    /// column under the planet's mean sunlight isn't the planet's mean (the
    /// physics is nonlinear in temperature), so this is not the 3-D model's
    /// cloud albedo; the 3-D model diagnoses its clouds
    /// ([`crate::moisture`]).
    pub const SINGLE_COLUMN_CLOUD_ALBEDO: f64 = 0.27;

    /// An Earth-like column under the global-mean sunlight of a planet with
    /// mean stellar flux `s0`, with [`Self::SINGLE_COLUMN_CLOUD_ALBEDO`].
    pub fn earth(s0: f64) -> Self {
        let mut params = PhysicsParams {
            ocean: OceanParams::slab(1.0),
            ..PhysicsParams::default()
        };
        params.radiation.cloud_albedo = Self::SINGLE_COLUMN_CLOUD_ALBEDO;
        Self {
            params,
            gas: Gas::EARTH_AIR,
            sigma: 0.5,
            gravity: 9.81,
            ps: 101_325.0,
            insolation: s0 / 4.0,
            wind: 5.0,
            theta: [330.0, 290.0],
            surface_temperature: 288.0,
        }
    }

    /// The column's two layers.
    pub fn levels(&self) -> Levels {
        Levels::two(self.sigma)
    }

    /// Fluxes of the current state.
    pub fn fluxes(&self) -> Fluxes {
        column_fluxes(
            &self.params.radiation,
            self.gas,
            &Column {
                ps: self.ps,
                surface_temperature: self.surface_temperature,
                theta: &self.theta,
                levels: &self.levels(),
                gravity: self.gravity,
                insolation: self.insolation,
                co2: self.params.co2,
                humidity: Humidity::Fixed,
                clouds: &[
                    Cloud::deck(self.sigma, self.params.radiation.cloud_albedo),
                    Cloud::NONE,
                ],
            },
        )
    }

    fn layer_exner(&self) -> [f64; 2] {
        let e = self.levels().layer_exner(self.gas, self.ps);
        [e[0], e[1]]
    }

    /// Sensible and latent heat flux (W m⁻²).
    pub fn surface_fluxes(&self) -> (f64, f64) {
        let air = surface_air_temperature(
            &self.params,
            self.gas,
            self.gravity,
            self.ps,
            self.layer_exner()[1],
            self.theta[1],
        );
        surface_fluxes(
            &self.params,
            self.gas,
            self.ps,
            self.surface_temperature,
            air,
            self.wind,
        )
    }

    /// Set the layer temperatures (K, upper first).
    pub fn set_temperatures(&mut self, t: [f64; 2]) {
        let exner = self.layer_exner();
        self.theta = [t[0] / exner[0], t[1] / exner[1]];
    }

    /// One step of `dt` seconds (forward Euler, radiation every step).
    pub fn step(&mut self, dt: f64) {
        let f = self.fluxes();
        let (sensible, latent) = self.surface_fluxes();
        let exner = self.layer_exner();
        let ds = [self.sigma, 1.0 - self.sigma];
        let heating = [
            f.layers[0] + ds[0] * latent,
            f.layers[1] + sensible + ds[1] * latent,
        ];
        for k in 0..2 {
            self.theta[k] +=
                dt * theta_tendency(self.gas, self.gravity, ds[k], exner[k], heating[k]) / self.ps;
        }
        if let Some(target) = convective_target(
            &self.params,
            self.gas,
            self.sigma,
            self.gravity,
            self.ps,
            self.theta,
        ) {
            let a = (dt / self.params.convection_time).min(1.0);
            for (theta, target) in self.theta.iter_mut().zip(target) {
                *theta += a * (target - *theta);
            }
        }
        self.surface_temperature +=
            dt * (f.surface() - sensible - latent) / self.params.heat_capacity();
    }

    /// Step until the top-of-atmosphere imbalance and the surface
    /// temperature settle (or 50 model years pass).
    pub fn equilibrium(mut self) -> Equilibrium {
        let dt = 3600.0;
        let mut steps = 0;
        let mut last = f64::INFINITY;
        while steps < (50.0 * 365.0 * 24.0) as u64 {
            for _ in 0..240 {
                self.step(dt);
            }
            steps += 240;
            let f = self.fluxes();
            let change = (self.surface_temperature - last).abs();
            last = self.surface_temperature;
            if !self.surface_temperature.is_finite() {
                break;
            }
            if change < 1e-6 && (f.absorbed - f.outgoing).abs() < 1e-3 {
                break;
            }
        }
        let exner = self.layer_exner();
        Equilibrium {
            fluxes: self.fluxes(),
            surface_temperature: self.surface_temperature,
            temperature: [self.theta[0] * exner[0], self.theta[1] * exner[1]],
            surface_fluxes: self.surface_fluxes(),
            steps,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GAS: Gas = Gas::EARTH_AIR;

    /// The years among the first 30 that end in a jump.
    fn jump_years(jumps: u32, first_year: u32, every: u32) -> Vec<u32> {
        let spin_up = OceanSpinUp {
            jumps,
            first_year,
            every,
            ..OceanSpinUp::default()
        };
        (1..=30).filter(|&y| spin_up.ends_in_jump(y)).collect()
    }

    #[test]
    fn the_spin_up_presets_jump_at_their_years() {
        let jumps = |spin_up: OceanSpinUp| -> Vec<u32> {
            (1..=30).filter(|&y| spin_up.ends_in_jump(y)).collect()
        };
        let currents = OceanSpinUp::with_currents();
        assert_eq!(jumps(currents), [4, 7]);
        assert_eq!(currents.layers, JumpLayers::Deep);
        assert_eq!((currents.extra_years, currents.limit), (20.0, 3.0));
        let still = OceanSpinUp::still();
        assert_eq!(jumps(still), [2, 3, 4]);
        assert_eq!(still.layers, JumpLayers::Both);
        assert_eq!((still.extra_years, still.limit), (3.0, 3.0));
        assert_eq!(OceanSpinUp::default().jumps, 0, "no jumps by default");
    }

    #[test]
    fn ocean_jumps_end_the_scheduled_years() {
        let default = OceanSpinUp::default();
        assert_eq!((default.first_year, default.every), (2, 1));
        assert_eq!(default.layers, JumpLayers::Both);
        assert_eq!(jump_years(3, 2, 1), [2, 3, 4], "the default schedule");
        assert_eq!(jump_years(0, 2, 1), [] as [u32; 0], "no jumps");
        assert_eq!(jump_years(3, 4, 1), [4, 5, 6]);
        assert_eq!(jump_years(3, 3, 2), [3, 5, 7]);
        assert_eq!(jump_years(2, 1, 3), [1, 4]);
        assert_eq!(jump_years(1, 5, 4), [5]);
        let last = |jumps, first_year, every| {
            OceanSpinUp {
                jumps,
                first_year,
                every,
                ..default
            }
            .last_jump_year()
        };
        assert_eq!(last(3, 2, 1), Some(4));
        assert_eq!(last(3, 3, 2), Some(7));
        assert_eq!(last(0, 2, 1), None);
    }

    fn mixed_surface() -> TiledSurface {
        TiledSurface {
            fraction: 0.3,
            ice: 0.0,
            albedo: [[0.08; 2], [0.2; 2], [0.5; 2]],
            temperature: [290.0, 305.0, 260.0],
        }
    }

    #[test]
    fn soil_supplies_evaporation_in_proportion_to_its_water_below_the_critical_level() {
        let land = LandParams::default();
        let critical = land.critical_wetness * land.soil_water_capacity;
        assert_eq!(land.wetness(0.0), 0.0);
        assert_eq!(land.wetness(0.5 * critical), 0.5);
        assert_eq!(land.wetness(critical), 1.0);
        assert_eq!(land.wetness(land.soil_water_capacity), 1.0);
    }

    #[test]
    fn the_bucket_overflows_only_what_it_cannot_hold() {
        let land = LandParams::default();
        let capacity = land.soil_water_capacity;
        assert_eq!(land.fill(50.0, 10.0), (60.0, 0.0));
        assert_eq!(land.fill(50.0, -10.0), (40.0, 0.0));
        let (water, overflow) = land.fill(capacity - 1.0, 5.0);
        assert_eq!(water, capacity);
        assert!((overflow - 4.0).abs() < 1e-12);
    }

    #[test]
    fn snow_brightens_the_land_as_it_covers_it() {
        let land = LandParams::default();
        assert_eq!(land.albedo_under(0.0), land.albedo);
        assert_eq!(land.albedo_under(land.snow_covering), land.snow_albedo);
        assert_eq!(
            land.albedo_under(5.0 * land.snow_covering),
            land.snow_albedo
        );
        let half = land.albedo_under(0.5 * land.snow_covering);
        assert!((half - 0.5 * (land.albedo + land.snow_albedo)).abs() < 1e-12);
    }

    #[test]
    fn only_an_ice_sheet_takes_the_ice_sheet_surface() {
        let land = LandParams::default();
        for snow in [0.0, 5.0, land.snow_covering, 500.0] {
            assert_eq!(land.surface_albedo(snow, false), land.albedo_under(snow));
        }
        assert_eq!(land.surface_roughness(false), land.roughness);
        assert_eq!(land.surface_albedo(500.0, true), land.ice_sheet_albedo);
        assert_eq!(land.surface_roughness(true), land.ice_sheet_roughness);
    }

    #[test]
    fn losing_full_cover_restarts_the_clock() {
        let land = LandParams::default();
        let full = land.snow_covering;
        assert_eq!(land.cover_age(100.0, full, 10.0), 110.0);
        assert_eq!(land.cover_age(0.0, 5.0 * full, 10.0), 10.0);
        assert_eq!(land.cover_age(1.0e9, 0.99 * full, 10.0), 0.0);
    }

    #[test]
    fn snow_melts_by_the_heat_available_but_no_more_than_there_is() {
        let dt = 1000.0;
        let surplus = 334.0;
        let melted = snow_melted(10.0, surplus, dt);
        assert!((melted * FUSION_HEAT - surplus * dt).abs() < 1e-9);
        assert_eq!(snow_melted(0.25, surplus, dt), 0.25);
        assert_eq!(snow_melted(10.0, -surplus, dt), 0.0);
    }

    #[test]
    fn the_skin_balance_finds_its_exact_root() {
        let (target, air, soil, drag, conductance) = (300.0, 290.0, 285.0, 40.0, 7.0);
        let absorbed = black_body(target) + drag * (target - air) + conductance * (target - soil);
        let turbulent = |t: f64| (drag * (t - air), drag);
        let t = solve_skin(280.0, absorbed, soil, conductance, turbulent);
        assert!((t - target).abs() < 1e-3, "{t}");
    }

    #[test]
    fn a_skin_with_only_radiation_sits_at_its_black_body_temperature() {
        let t = solve_skin(260.0, 400.0, 0.0, 0.0, |_| (0.0, 0.0));
        let exact = fms_math::pow(400.0 / STEFAN_BOLTZMANN, 0.25);
        assert!((t - exact).abs() < 1e-3, "{t} vs {exact}");
    }

    #[test]
    fn skin_and_soil_follow_the_analytic_solution_at_any_step() {
        // With a linear air flux in place of the emission (cancelled in the
        // closure), the soil relaxes to `T_eq` at a rate `κ/C`, with
        // `κ = Λ B/(Λ + B)`; implicit steps shrink its distance by exactly
        // `1/(1 + κ Δt/C)`, and are stable however long.
        let land = LandParams::default();
        let (lambda, capacity) = (land.skin_conductance, land.heat_capacity());
        let (b, air, absorbed, start) = (40.0, 290.0, 200.0, 280.0);
        let equilibrium = air + absorbed / b;
        let kappa = lambda * b / (lambda + b);
        let run = |dt: f64, steps: i32| {
            let turbulent = |t: f64| {
                let cubed = 4.0 * STEFAN_BOLTZMANN * t * t * t;
                (b * (t - air) - black_body(t), b - cubed)
            };
            let (mut skin, mut soil) = (start, start);
            let mut previous = (equilibrium - soil).abs();
            for _ in 0..steps {
                let g = land.step_conductance(dt);
                skin = solve_skin(skin, absorbed, soil, g, turbulent);
                soil += dt * (absorbed - b * (skin - air)) / capacity;
                let distance = (equilibrium - soil).abs();
                assert!(distance <= previous, "overshot at Δt {dt}");
                previous = distance;
            }
            soil
        };
        let distance = start - equilibrium;
        for (dt, steps) in [(900.0, 960), (86_400.0, 10), (864_000.0, 10)] {
            let exact = equilibrium + distance / (1.0 + kappa * dt / capacity).powi(steps);
            let soil = run(dt, steps);
            assert!((soil - exact).abs() < 1e-6, "Δt {dt}: {soil} vs {exact}");
        }
        // Short steps converge on the continuous solution.
        let time = 900.0 * 960.0;
        let continuous = equilibrium + distance * fms_math::exp(-kappa * time / capacity);
        assert!((run(900.0, 960) - continuous).abs() < 0.03);
    }

    #[test]
    fn splitting_radiation_between_tiles_conserves_it() {
        let s = mixed_surface();
        let t = s.split(180.0, 110.0, 400.0);
        let mean = |x: [f64; 3]| (1.0 - s.fraction) * x[OCEAN] + s.fraction * x[LAND];
        assert!((mean(t.shortwave) - 180.0).abs() < 1e-10);
        assert!((mean(t.longwave_up) - 400.0).abs() < 1e-10);
        // The brighter tile absorbs less, the warmer tile emits more.
        assert!(t.shortwave[LAND] < t.shortwave[OCEAN]);
        assert!(t.longwave_up[LAND] > t.longwave_up[OCEAN]);
    }

    #[test]
    fn the_radiating_temperature_has_the_tiles_mean_emission() {
        let s = mixed_surface();
        let t = s.radiating_temperature();
        assert!((t * t * t * t / s.fourth_powers().1 - 1.0).abs() < 1e-14);
        assert!(t > 290.0 && t < 305.0);
    }

    #[test]
    fn equal_band_albedos_split_sunlight_as_one_band_does_bit_for_bit() {
        let s = mixed_surface();
        let co_albedo = 1.0 - s.mean_albedo()[VISIBLE];
        for visible in [0.0, 65.3, 180.0] {
            let t = s.split(180.0, visible, 400.0);
            for tile in [OCEAN, LAND] {
                let one_band = 180.0 * ((1.0 - s.albedo[tile][VISIBLE]) / co_albedo);
                assert_eq!(t.shortwave[tile], one_band);
            }
        }
    }

    #[test]
    fn tiles_split_sunlight_band_by_band() {
        let s = TiledSurface {
            albedo: [[0.06, 0.1], [0.5, 0.3], [0.5, 0.3]],
            ..mixed_surface()
        };
        let (visible, infrared) = (110.0, 70.0);
        let t = s.split(visible + infrared, visible, 400.0);
        let mean = s.mean_albedo();
        for tile in [OCEAN, LAND] {
            let a = s.albedo[tile];
            let expected = visible / (1.0 - mean[VISIBLE]) * (1.0 - a[VISIBLE])
                + infrared / (1.0 - mean[NEAR_INFRARED]) * (1.0 - a[NEAR_INFRARED]);
            assert!((t.shortwave[tile] - expected).abs() < 1e-12);
        }
        let total = (1.0 - s.fraction) * t.shortwave[OCEAN] + s.fraction * t.shortwave[LAND];
        assert!((total - 180.0).abs() < 1e-12, "{total}");
    }

    #[test]
    fn an_ice_tile_shares_the_cells_radiation_conservatively() {
        let s = TiledSurface {
            ice: 0.4,
            albedo: [[0.06, 0.1], [0.5, 0.3], [0.57, 0.31]],
            temperature: [272.0, 280.0, 255.0],
            ..mixed_surface()
        };
        let w = s.weights();
        assert!((w[OCEAN] + w[LAND] + w[ICE] - 1.0).abs() < 1e-15);
        let t = s.split(180.0, 110.0, 300.0);
        let mean = |x: [f64; 3]| w[OCEAN] * x[OCEAN] + w[LAND] * x[LAND] + w[ICE] * x[ICE];
        assert!((mean(t.shortwave) - 180.0).abs() < 1e-10);
        assert!((mean(t.longwave_up) - 300.0).abs() < 1e-10);
        assert!(t.shortwave[ICE] < t.shortwave[OCEAN]);
    }

    #[test]
    fn a_tile_covering_the_whole_cell_is_used_exactly() {
        for (fraction, tile) in [(0.0, OCEAN), (1.0, LAND)] {
            let s = TiledSurface {
                fraction,
                ..mixed_surface()
            };
            let t = s.split(180.0, 110.0, 400.0);
            assert_eq!(s.radiating_temperature(), s.temperature[tile]);
            assert_eq!(s.mean_albedo(), s.albedo[tile]);
            assert_eq!(t.shortwave[tile], 180.0);
            assert_eq!(t.longwave_up[tile], 400.0);
        }
    }

    #[test]
    fn blending_tiles_weights_them_by_land_fraction() {
        let ocean = TileFlux {
            sensible: 10.0,
            evaporation: 1e-5,
            drag: 1e-5,
            ..TileFlux::default()
        };
        let land = TileFlux {
            sensible: 50.0,
            evaporation: 3e-5,
            drag: 5e-5,
            ..TileFlux::default()
        };
        assert_eq!(TileFlux::blend(ocean, land, 0.0), ocean);
        assert_eq!(TileFlux::blend(ocean, land, 1.0), land);
        let mean = TileFlux::blend(ocean, land, 0.25);
        assert!((mean.sensible - 20.0).abs() < 1e-12);
        assert!((mean.drag - 2e-5).abs() < 1e-17);
    }

    #[test]
    fn an_isothermal_rest_profile_has_the_barometric_pressure() {
        let t = 260.0;
        let profile = RestProfile::new(GAS, 1.0e5, t, &[3.0e4, 7.0e4, 9.0e4], &[t; 3]);
        for z in [100.0, 1500.0, 6000.0, 12_000.0] {
            let phi = 9.81 * z;
            let expected = 1.0e5 * fms_math::exp(-phi / (GAS.gas_constant * t));
            let p = profile.pressure_at(phi);
            assert!(
                (p / expected - 1.0).abs() < 1e-12,
                "z {z}: {p} vs {expected}"
            );
        }
    }

    #[test]
    fn pressure_at_inverts_the_hypsometric_integral() {
        let (p, t) = ([3.0e4, 7.0e4, 9.0e4], [225.0, 265.0, 282.0]);
        let profile = RestProfile::new(GAS, 1.0e5, 290.0, &p, &t);
        // Layer thicknesses by the trapezoid in ln p, in segments of 1 Pa.
        for phi in [500.0, 8_000.0, 30_000.0, 60_000.0] {
            let target = profile.pressure_at(phi);
            let mut integral = 0.0;
            let mut a = 1.0e5;
            while a > target {
                let b = (a - 1.0_f64).max(target);
                let mid = 0.5 * (profile.temperature(a) + profile.temperature(b));
                integral += GAS.gas_constant * mid * fms_math::ln(a / b);
                a = b;
            }
            assert!((integral / phi - 1.0).abs() < 1e-4, "{integral} vs {phi}");
        }
    }

    #[test]
    fn the_profile_passes_through_its_nodes() {
        let (p, t) = ([3.0e4, 7.0e4, 9.0e4], [225.0, 265.0, 282.0]);
        let profile = RestProfile::new(GAS, 1.0e5, 290.0, &p, &t);
        for (p, t) in p.iter().zip(t) {
            assert!((profile.temperature(*p) - t).abs() < 1e-9);
        }
        assert_eq!(profile.temperature(1.0e3), 225.0);
    }
}
