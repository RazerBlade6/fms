//! Latitude/longitude conversions.
//!
//! Convention: `+z` is the north pole, `+x` is (lat 0, lon 0), `+y` is
//! (lat 0, lon 90° E). Angles are in radians; latitude is in `[-π/2, π/2]`,
//! longitude in `(-π, π]`.

use glam::DVec3;

/// Unit vector for a latitude/longitude pair (radians).
#[inline]
pub fn from_lat_lon(lat: f64, lon: f64) -> DVec3 {
    let (sin_lat, cos_lat) = fms_math::sin_cos(lat);
    let (sin_lon, cos_lon) = fms_math::sin_cos(lon);
    DVec3::new(cos_lat * cos_lon, cos_lat * sin_lon, sin_lat)
}

/// Latitude and longitude (radians) of a unit vector.
#[inline]
pub fn to_lat_lon(p: DVec3) -> (f64, f64) {
    let lat = fms_math::atan2(p.z, fms_math::hypot(p.x, p.y));
    let lon = fms_math::atan2(p.y, p.x);
    (lat, lon)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn round_trip() {
        for &(lat, lon) in &[
            (0.0, 0.0),
            (0.5, -2.0),
            (-1.2, 3.0),
            (FRAC_PI_2 - 1e-9, 1.0),
        ] {
            let (la, lo) = to_lat_lon(from_lat_lon(lat, lon));
            assert!(
                (la - lat).abs() < 1e-12 && (lo - lon).abs() < 1e-9,
                "{lat},{lon} -> {la},{lo}"
            );
        }
    }

    #[test]
    fn axes() {
        assert!((from_lat_lon(FRAC_PI_2, 0.0) - DVec3::Z).length() < 1e-15);
        assert!((from_lat_lon(0.0, FRAC_PI_2) - DVec3::Y).length() < 1e-15);
        assert!((from_lat_lon(0.0, 0.0) - DVec3::X).length() < 1e-15);
    }
}
