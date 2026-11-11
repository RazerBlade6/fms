//! Planet parameters, orbit and insolation.
//!
//! [`PlanetParams`] holds what the user sets. [`Planet::new`] validates and clamps
//! it, then derives everything else: rotation rate, Coriolis parameter, orbital
//! period, stellar flux, and the Sun's position at any time.
//!
//! **Time convention:** `t = 0` is the northern vernal equinox at midnight at
//! longitude 0, matching the scheduler's calendar (year fraction 0 is the
//! equinox). The Sun's position uses *mean* solar time; the equation of time (a
//! few minutes on Earth) is ignored.
//!
//! See `docs/climate/PLAN.md` §4.1 for the design.

mod orbit;
mod params;
mod seasons;

use std::f64::consts::{PI, TAU};

use fms_grid::coords::from_lat_lon;
use glam::DVec3;
use rayon::prelude::*;

pub use orbit::{Orbit, OrbitState, eccentric_anomaly, eccentric_from_true, true_anomaly};
pub use params::{Adjustment, PlanetParams, SOLAR_CONSTANT_1AU};
pub use seasons::{Hemisphere, Season};

/// A validated planet with derived quantities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Planet {
    params: PlanetParams,
    orbit: Orbit,
    rotation_rate: f64,
    mean_flux: f64,
}

impl Planet {
    /// Validate `params` (clamping out-of-range values) and derive the rest.
    /// Returns the planet and a list of every change or warning.
    pub fn new(params: PlanetParams) -> (Self, Vec<Adjustment>) {
        let (params, adjustments) = params::sanitize(params);
        let period = params.year_days * params.solar_day;
        // Prograde rotation: the solar day is longer than the sidereal day because
        // the planet also moves along its orbit: 1/T_sid = 1/T_sol + 1/T_orbit.
        let sidereal_day = 1.0 / (1.0 / params.solar_day + 1.0 / period);
        let orbit = Orbit::new(
            period,
            params.eccentricity,
            params.obliquity_deg.to_radians(),
            params.perihelion_longitude_deg.to_radians(),
        );
        let planet = Self {
            params,
            orbit,
            rotation_rate: TAU / sidereal_day,
            mean_flux: SOLAR_CONSTANT_1AU * params.luminosity / params.semi_major_axis_au.powi(2),
        };
        (planet, adjustments)
    }

    /// Present-day Earth.
    pub fn earth() -> Self {
        Self::new(PlanetParams::earth()).0
    }

    /// The validated parameters.
    pub fn params(&self) -> &PlanetParams {
        &self.params
    }

    /// The orbit.
    pub fn orbit(&self) -> &Orbit {
        &self.orbit
    }

    /// Rotation rate Ω relative to the stars (rad s⁻¹).
    pub fn rotation_rate(&self) -> f64 {
        self.rotation_rate
    }

    /// Sidereal day: one rotation relative to the stars (s).
    pub fn sidereal_day(&self) -> f64 {
        TAU / self.rotation_rate
    }

    /// Solar day (s).
    pub fn solar_day(&self) -> f64 {
        self.params.solar_day
    }

    /// Orbital period (s).
    pub fn year_seconds(&self) -> f64 {
        self.orbit.period
    }

    /// Coriolis parameter `f = 2Ω sin φ` at latitude `lat` (radians).
    pub fn coriolis(&self, lat: f64) -> f64 {
        2.0 * self.rotation_rate * fms_math::sin(lat)
    }

    /// Coriolis parameter at a unit vector (`2Ω z`).
    pub fn coriolis_at(&self, p: DVec3) -> f64 {
        2.0 * self.rotation_rate * p.z
    }

    /// Mean stellar flux at the semi-major axis, `S₀` (W m⁻²).
    pub fn mean_flux(&self) -> f64 {
        self.mean_flux
    }

    /// Orbital state at time `t`.
    pub fn orbit_state(&self, t: f64) -> OrbitState {
        self.orbit.state(t)
    }

    /// Unit vector towards the Sun at time `t`, in planet-fixed coordinates
    /// (+z north, +x at longitude 0). The subsolar point is at latitude δ and at
    /// the longitude where it is local noon: noon at longitude 0 at half a day.
    pub fn sun_direction(&self, t: f64) -> DVec3 {
        let state = self.orbit.state(t);
        let day_fraction = (t / self.params.solar_day).rem_euclid(1.0);
        let subsolar_lon = PI - TAU * day_fraction;
        from_lat_lon(state.declination, subsolar_lon)
    }

    /// Top-of-atmosphere flux at the unit vector `p` at time `t` (W m⁻²).
    pub fn insolation(&self, p: DVec3, t: f64) -> f64 {
        let state = self.orbit.state(t);
        let sun = self.sun_direction(t);
        self.mean_flux / state.distance_ratio.powi(2) * p.dot(sun).max(0.0)
    }

    /// Top-of-atmosphere flux at many points at time `t`, in parallel.
    pub fn insolation_into(&self, points: &[DVec3], t: f64, out: &mut [f64]) {
        assert_eq!(points.len(), out.len());
        let state = self.orbit.state(t);
        let sun = self.sun_direction(t);
        let flux = self.mean_flux / state.distance_ratio.powi(2);
        out.par_iter_mut()
            .zip(points)
            .for_each(|(o, p)| *o = flux * p.dot(sun).max(0.0));
    }

    /// Top-of-atmosphere flux averaged over one day at latitude `lat` (radians),
    /// using the orbital state at time `t` (assumed constant over the day).
    ///
    /// `Q = (S₀/π)(a/r)² (h₀ sin φ sin δ + cos φ cos δ sin h₀)` with the sunset hour
    /// angle `cos h₀ = −tan φ tan δ`, clamped for polar day and night.
    pub fn daily_mean_insolation(&self, lat: f64, t: f64) -> f64 {
        let state = self.orbit.state(t);
        let (sin_lat, cos_lat) = fms_math::sin_cos(lat);
        let (sin_dec, cos_dec) = fms_math::sin_cos(state.declination);
        let cos_h0 = if cos_lat * cos_dec == 0.0 {
            // At a pole: day if the Sun is above that pole's horizon.
            if sin_lat * sin_dec > 0.0 { -1.0 } else { 1.0 }
        } else {
            (-(sin_lat * sin_dec) / (cos_lat * cos_dec)).clamp(-1.0, 1.0)
        };
        let h0 = fms_math::acos(cos_h0);
        self.mean_flux / PI / state.distance_ratio.powi(2)
            * (h0 * sin_lat * sin_dec + cos_lat * cos_dec * fms_math::sin(h0))
    }
}
