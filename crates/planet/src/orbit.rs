//! Keplerian orbit and the Sun's position in the sky.

use std::f64::consts::TAU;

/// Solve Kepler's equation `E − e·sin E = M` for the eccentric anomaly `E`.
///
/// Newton iteration from `E₀ = M + e·sin M`; converges to machine precision in a
/// few iterations for `e < 0.9`.
pub fn eccentric_anomaly(mean_anomaly: f64, e: f64) -> f64 {
    let m = mean_anomaly.rem_euclid(TAU);
    let mut ecc = m + e * m.sin();
    for _ in 0..50 {
        let delta = (ecc - e * ecc.sin() - m) / (1.0 - e * ecc.cos());
        ecc -= delta;
        if delta.abs() < 1e-15 {
            break;
        }
    }
    ecc
}

/// True anomaly from the eccentric anomaly.
pub fn true_anomaly(eccentric_anomaly: f64, e: f64) -> f64 {
    let half = eccentric_anomaly / 2.0;
    2.0 * ((1.0 + e).sqrt() * half.sin()).atan2((1.0 - e).sqrt() * half.cos())
}

/// Eccentric anomaly from the true anomaly (inverse of [`true_anomaly`]).
pub fn eccentric_from_true(true_anomaly: f64, e: f64) -> f64 {
    let half = true_anomaly / 2.0;
    2.0 * ((1.0 - e).sqrt() * half.sin()).atan2((1.0 + e).sqrt() * half.cos())
}

/// Where the planet is on its orbit at one instant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrbitState {
    /// Mean anomaly (radians, `[0, 2π)`).
    pub mean_anomaly: f64,
    /// True anomaly: angle from perihelion (radians).
    pub true_anomaly: f64,
    /// Star–planet distance divided by the semi-major axis.
    pub distance_ratio: f64,
    /// Solar longitude: the Sun's position along the ecliptic as seen from the
    /// planet (radians, `[0, 2π)`). 0 = northern vernal equinox, π/2 = northern
    /// summer solstice.
    pub solar_longitude: f64,
    /// Solar declination (radians): latitude of the subsolar point.
    pub declination: f64,
}

/// A Keplerian orbit with a fixed spin axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Orbit {
    /// Orbital period (s).
    pub period: f64,
    /// Eccentricity.
    pub eccentricity: f64,
    /// Obliquity (radians).
    pub obliquity: f64,
    /// Solar longitude at perihelion (radians).
    pub perihelion_longitude: f64,
    /// Mean anomaly at `t = 0`, which is the northern vernal equinox.
    mean_anomaly_at_equinox: f64,
}

impl Orbit {
    /// An orbit whose `t = 0` is the northern vernal equinox.
    pub fn new(period: f64, eccentricity: f64, obliquity: f64, perihelion_longitude: f64) -> Self {
        // At the equinox the solar longitude is 0, so the true anomaly is −ϖ.
        let nu0 = -perihelion_longitude;
        let e0 = eccentric_from_true(nu0, eccentricity);
        let m0 = (e0 - eccentricity * e0.sin()).rem_euclid(TAU);
        Self {
            period,
            eccentricity,
            obliquity,
            perihelion_longitude,
            mean_anomaly_at_equinox: m0,
        }
    }

    /// Orbital state `t` seconds after a northern vernal equinox.
    pub fn state(&self, t: f64) -> OrbitState {
        let e = self.eccentricity;
        let mean_anomaly = (self.mean_anomaly_at_equinox + TAU * t / self.period).rem_euclid(TAU);
        let ecc = eccentric_anomaly(mean_anomaly, e);
        let nu = true_anomaly(ecc, e);
        let solar_longitude = (nu + self.perihelion_longitude).rem_euclid(TAU);
        let declination = (self.obliquity.sin() * solar_longitude.sin()).asin();
        OrbitState {
            mean_anomaly,
            true_anomaly: nu,
            distance_ratio: 1.0 - e * ecc.cos(),
            solar_longitude,
            declination,
        }
    }
}
