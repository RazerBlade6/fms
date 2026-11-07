//! User-facing planet parameters, validation and clamping.

use std::fmt;

/// Mean solar flux at 1 AU from a 1 L☉ star (W m⁻²).
pub const SOLAR_CONSTANT_1AU: f64 = 1361.0;

/// Everything the user can set about a planet. `Default` is Earth.
///
/// Values outside the supported ranges are clamped (with a warning) by
/// [`Planet::new`](crate::Planet::new). The ranges keep the climate model
/// numerically stable; see `docs/PLAN.md` §2.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanetParams {
    /// Mean radius (m).
    pub radius: f64,
    /// Surface gravity (m s⁻²).
    pub gravity: f64,
    /// Length of the solar day, one day/night cycle (s). Rotation is prograde.
    pub solar_day: f64,
    /// Year length, in solar days.
    pub year_days: f64,
    /// Axial tilt / obliquity (degrees).
    pub obliquity_deg: f64,
    /// Orbital eccentricity.
    pub eccentricity: f64,
    /// Solar longitude at perihelion (degrees): where in the seasonal cycle the
    /// planet is closest to its star. 0° = northern vernal equinox, 90° =
    /// northern summer solstice. Earth: ~282.9° (early January).
    pub perihelion_longitude_deg: f64,
    /// Stellar luminosity (L☉).
    pub luminosity: f64,
    /// Orbital semi-major axis (AU).
    pub semi_major_axis_au: f64,
    /// Mean surface pressure (Pa).
    pub surface_pressure: f64,
    /// CO₂ concentration (ppm by volume).
    pub co2_ppm: f64,
    /// Fraction of the surface covered by ocean, 0–1.
    pub ocean_fraction: f64,
}

impl Default for PlanetParams {
    fn default() -> Self {
        Self::earth()
    }
}

impl PlanetParams {
    /// Present-day Earth.
    pub const fn earth() -> Self {
        Self {
            radius: 6.371e6,
            gravity: 9.81,
            solar_day: 86_400.0,
            year_days: 365.2422,
            obliquity_deg: 23.44,
            eccentricity: 0.0167,
            perihelion_longitude_deg: 282.9,
            luminosity: 1.0,
            semi_major_axis_au: 1.0,
            surface_pressure: 101_325.0,
            co2_ppm: 420.0,
            ocean_fraction: 0.71,
        }
    }
}

/// Supported range of each parameter.
pub(crate) struct Limit {
    pub name: &'static str,
    pub min: f64,
    pub max: f64,
    pub get: fn(&PlanetParams) -> f64,
    pub set: fn(&mut PlanetParams, f64),
}

macro_rules! limit {
    ($field:ident, $min:expr, $max:expr) => {
        Limit {
            name: stringify!($field),
            min: $min,
            max: $max,
            get: |p| p.$field,
            set: |p, v| p.$field = v,
        }
    };
}

/// Supported ranges. Wide enough for most fantasy worlds; extremes such as
/// tidal locking or 90° tilt are out of scope for now (`docs/PLAN.md` M7).
pub(crate) const LIMITS: &[Limit] = &[
    limit!(radius, 1.0e6, 2.0e7),
    limit!(gravity, 1.0, 30.0),
    limit!(solar_day, 4.0 * 3600.0, 240.0 * 3600.0),
    limit!(year_days, 20.0, 5000.0),
    limit!(obliquity_deg, 0.0, 60.0),
    limit!(eccentricity, 0.0, 0.3),
    limit!(perihelion_longitude_deg, 0.0, 360.0),
    limit!(luminosity, 0.05, 20.0),
    limit!(semi_major_axis_au, 0.05, 50.0),
    limit!(surface_pressure, 1.0e4, 1.0e6),
    limit!(co2_ppm, 0.0, 100_000.0),
    limit!(ocean_fraction, 0.0, 1.0),
];

/// Stellar flux outside which the climate is probably runaway-hot or frozen;
/// such values are allowed but flagged.
pub(crate) const FLUX_WARNING_RANGE: (f64, f64) = (500.0, 3000.0);

/// A parameter that was changed or flagged during validation.
#[derive(Clone, Debug, PartialEq)]
pub enum Adjustment {
    /// The value was outside the supported range and was clamped.
    Clamped {
        /// Parameter name.
        parameter: &'static str,
        /// Value supplied.
        from: f64,
        /// Value used.
        to: f64,
    },
    /// The value was NaN or infinite and was replaced by Earth's value.
    Replaced {
        /// Parameter name.
        parameter: &'static str,
        /// Value used.
        to: f64,
    },
    /// The perihelion longitude was wrapped into `[0, 360)`.
    Wrapped {
        /// Value supplied.
        from: f64,
        /// Value used.
        to: f64,
    },
    /// Allowed, but the resulting stellar flux is far from habitable.
    ExtremeFlux {
        /// Mean stellar flux at the planet (W m⁻²).
        flux: f64,
    },
}

impl fmt::Display for Adjustment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Adjustment::Clamped {
                parameter,
                from,
                to,
            } => {
                write!(
                    f,
                    "{parameter} {from} is outside the supported range; using {to}"
                )
            }
            Adjustment::Replaced { parameter, to } => {
                write!(f, "{parameter} was not a finite number; using {to}")
            }
            Adjustment::Wrapped { from, to } => {
                write!(f, "perihelion longitude {from}° wrapped to {to}°")
            }
            Adjustment::ExtremeFlux { flux } => write!(
                f,
                "stellar flux {flux:.0} W/m² is far from Earth's 1361 W/m²; expect a frozen or \
                 runaway climate"
            ),
        }
    }
}

/// Validate and clamp `params`, returning the usable values and what changed.
pub(crate) fn sanitize(mut params: PlanetParams) -> (PlanetParams, Vec<Adjustment>) {
    let earth = PlanetParams::earth();
    let mut adjustments = Vec::new();
    for limit in LIMITS {
        let value = (limit.get)(&params);
        if !value.is_finite() {
            let to = (limit.get)(&earth);
            (limit.set)(&mut params, to);
            adjustments.push(Adjustment::Replaced {
                parameter: limit.name,
                to,
            });
        } else if limit.name == "perihelion_longitude_deg" {
            let to = value.rem_euclid(360.0);
            if to != value {
                (limit.set)(&mut params, to);
                adjustments.push(Adjustment::Wrapped { from: value, to });
            }
        } else if value < limit.min || value > limit.max {
            let to = value.clamp(limit.min, limit.max);
            (limit.set)(&mut params, to);
            adjustments.push(Adjustment::Clamped {
                parameter: limit.name,
                from: value,
                to,
            });
        }
    }
    let flux = SOLAR_CONSTANT_1AU * params.luminosity / params.semi_major_axis_au.powi(2);
    if flux < FLUX_WARNING_RANGE.0 || flux > FLUX_WARNING_RANGE.1 {
        adjustments.push(Adjustment::ExtremeFlux { flux });
    }
    (params, adjustments)
}
