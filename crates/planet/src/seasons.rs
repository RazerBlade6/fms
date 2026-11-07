//! Astronomical seasons from the planet's orbit.

use std::f64::consts::{FRAC_PI_2, TAU};

use crate::Planet;
use crate::orbit::eccentric_from_true;

/// Hemisphere, for naming seasons.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Hemisphere {
    /// Northern hemisphere.
    North,
    /// Southern hemisphere.
    South,
}

/// Astronomical season, bounded by equinoxes and solstices.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Season {
    /// From the vernal (spring) equinox.
    Spring,
    /// From the summer solstice.
    Summer,
    /// From the autumnal equinox.
    Autumn,
    /// From the winter solstice.
    Winter,
}

impl Season {
    /// All seasons, in order from the vernal equinox.
    pub const ALL: [Season; 4] = [
        Season::Spring,
        Season::Summer,
        Season::Autumn,
        Season::Winter,
    ];
}

impl Planet {
    /// Season in `hemisphere` at time `t`, from the solar longitude: northern
    /// spring is 0–90°, summer 90–180°, autumn 180–270°, winter 270–360°. The
    /// southern hemisphere is offset by half a year.
    pub fn season(&self, t: f64, hemisphere: Hemisphere) -> Season {
        let lon = self.orbit_state(t).solar_longitude;
        let quarter = ((lon / FRAC_PI_2).floor() as usize) % 4;
        let quarter = match hemisphere {
            Hemisphere::North => quarter,
            Hemisphere::South => (quarter + 2) % 4,
        };
        Season::ALL[quarter]
    }

    /// Time in `[0, year)` at which the solar longitude equals `lon` (radians).
    /// Exact: inverts true anomaly → eccentric anomaly → mean anomaly.
    pub fn time_of_solar_longitude(&self, lon: f64) -> f64 {
        let orbit = self.orbit();
        let e = orbit.eccentricity;
        let nu = lon - orbit.perihelion_longitude;
        let ecc = eccentric_from_true(nu, e);
        let mean = ecc - e * fms_math::sin(ecc);
        let equinox = self.orbit_state(0.0).mean_anomaly;
        ((mean - equinox).rem_euclid(TAU) / TAU) * orbit.period
    }

    /// Start times in `[0, year)` of northern spring, summer, autumn and winter:
    /// the vernal equinox, summer solstice, autumnal equinox and winter solstice.
    pub fn season_starts(&self) -> [f64; 4] {
        [0.0, 1.0, 2.0, 3.0].map(|q| self.time_of_solar_longitude(q * FRAC_PI_2))
    }

    /// Lengths of northern spring, summer, autumn and winter, in solar days.
    /// Equal on a circular orbit; on an eccentric one, seasons near
    /// perihelion are shorter (Kepler's second law).
    pub fn season_lengths_days(&self) -> [f64; 4] {
        let starts = self.season_starts();
        let year = self.year_seconds();
        [0, 1, 2, 3].map(|i| {
            let end = if i == 3 { year } else { starts[i + 1] };
            (end - starts[i]) / self.solar_day()
        })
    }
}
