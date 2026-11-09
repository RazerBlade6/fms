//! An idealized continent, for developing the land surface before running on
//! real terrain: one rectangular landmass in an ocean, with a low plain and
//! a north–south mountain ridge near its west coast (where the mid-latitude
//! westerlies meet it first).

use fms_grid::coords::to_lat_lon;
use glam::DVec3;
use rayon::prelude::*;

/// One rectangular continent. [`Default`] spans 35° S to 60° N.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Continent {
    /// Western and eastern coasts (degrees east).
    pub west_deg: f64,
    /// See `west_deg`.
    pub east_deg: f64,
    /// Southern and northern coasts (degrees north).
    pub south_deg: f64,
    /// See `south_deg`.
    pub north_deg: f64,
    /// Height of the plain (m).
    pub plain: f64,
    /// Height of the ridge above the plain (m).
    pub ridge: f64,
    /// Longitude of the ridge's crest (degrees east).
    pub ridge_lon_deg: f64,
    /// The ridge's Gaussian width (degrees of longitude).
    pub ridge_width_deg: f64,
    /// Latitude over which the ridge tapers to nothing at the continent's
    /// northern and southern coasts (degrees).
    pub ridge_taper_deg: f64,
    /// Depth of the surrounding ocean (m, positive).
    pub ocean_depth: f64,
}

impl Default for Continent {
    fn default() -> Self {
        Self {
            west_deg: -30.0,
            east_deg: 30.0,
            south_deg: -35.0,
            north_deg: 60.0,
            plain: 200.0,
            ridge: 2500.0,
            ridge_lon_deg: -18.0,
            ridge_width_deg: 4.0,
            ridge_taper_deg: 10.0,
            ocean_depth: 4000.0,
        }
    }
}

impl Continent {
    /// Elevation (m above sea level) at unit vector `p`.
    pub fn elevation_at(&self, p: DVec3) -> f64 {
        let (lat, lon) = to_lat_lon(p);
        let (lat, lon) = (lat.to_degrees(), lon.to_degrees());
        if lon < self.west_deg
            || lon > self.east_deg
            || lat < self.south_deg
            || lat > self.north_deg
        {
            return -self.ocean_depth;
        }
        let from_coast = (lat - self.south_deg).min(self.north_deg - lat);
        let t = (from_coast / self.ridge_taper_deg).clamp(0.0, 1.0);
        let taper = t * t * (3.0 - 2.0 * t);
        let x = (lon - self.ridge_lon_deg) / self.ridge_width_deg;
        self.plain + self.ridge * taper * fms_math::exp(-0.5 * x * x)
    }

    /// Elevation (m) at each of `points`.
    pub fn elevation(&self, points: &[DVec3]) -> Vec<f64> {
        points.par_iter().map(|&p| self.elevation_at(p)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn point(lat_deg: f64, lon_deg: f64) -> DVec3 {
        let (lat, lon) = (lat_deg.to_radians(), lon_deg.to_radians());
        DVec3::new(
            fms_math::cos(lat) * fms_math::cos(lon),
            fms_math::cos(lat) * fms_math::sin(lon),
            fms_math::sin(lat),
        )
    }

    #[test]
    fn land_is_inside_the_coasts_and_the_ridge_peaks_at_its_crest() {
        let c = Continent::default();
        assert_eq!(c.elevation_at(point(10.0, 120.0)), -c.ocean_depth);
        assert_eq!(c.elevation_at(point(70.0, 0.0)), -c.ocean_depth);
        let crest = c.elevation_at(point(15.0, c.ridge_lon_deg));
        assert!((crest - (c.plain + c.ridge)).abs() < 1e-6, "crest {crest}");
        let east = c.elevation_at(point(15.0, 25.0));
        assert!((east - c.plain).abs() < 1.0, "plain {east}");
    }

    #[test]
    fn the_ridge_tapers_away_at_the_northern_and_southern_coasts() {
        let c = Continent::default();
        let near_north = c.elevation_at(point(c.north_deg - 0.01, c.ridge_lon_deg));
        assert!((near_north - c.plain).abs() < 1.0, "{near_north}");
    }
}
