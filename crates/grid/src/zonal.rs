//! Zonal means: averages over latitude bands, the standard way to summarize
//! a general circulation (jets, Hadley cells, temperature contrasts).

use glam::DVec3;

use crate::coords::to_lat_lon;
use crate::geometry::MeshGeometry;
use crate::mesh::SphereMesh;

/// Equal-width latitude bands, each cell assigned whole to the band holding
/// its centre.
#[derive(Clone, Debug)]
pub struct LatitudeBands {
    band_of_cell: Vec<u32>,
    /// Cell area in each band (unit sphere).
    areas: Vec<f64>,
}

impl LatitudeBands {
    /// `count` bands from the south pole to the north pole.
    ///
    /// # Panics
    /// If `count` is zero or some band holds no cell centre (too many bands
    /// for the mesh).
    pub fn new(mesh: &SphereMesh, geometry: &MeshGeometry, count: usize) -> Self {
        assert!(count > 0, "need at least one band");
        let width = std::f64::consts::PI / count as f64;
        let band_of_cell: Vec<u32> = mesh
            .cell_centers()
            .iter()
            .map(|&p| {
                let (lat, _) = to_lat_lon(p);
                (((lat + std::f64::consts::FRAC_PI_2) / width) as usize).min(count - 1) as u32
            })
            .collect();
        let mut areas = vec![0.0; count];
        for (&b, &a) in band_of_cell.iter().zip(&geometry.cell_areas) {
            areas[b as usize] += a;
        }
        assert!(
            areas.iter().all(|&a| a > 0.0),
            "{count} bands is too many for {} cells",
            mesh.n_cells()
        );
        Self {
            band_of_cell,
            areas,
        }
    }

    /// Number of bands.
    pub fn len(&self) -> usize {
        self.areas.len()
    }

    /// Whether there are no bands (never true).
    pub fn is_empty(&self) -> bool {
        self.areas.is_empty()
    }

    /// Cell area in each band (unit sphere), south to north.
    pub fn areas(&self) -> &[f64] {
        &self.areas
    }

    /// Latitude (radians) of each band's centre, south to north.
    pub fn latitudes(&self) -> Vec<f64> {
        let width = std::f64::consts::PI / self.len() as f64;
        (0..self.len())
            .map(|b| -std::f64::consts::FRAC_PI_2 + (b as f64 + 0.5) * width)
            .collect()
    }

    /// The band holding each cell.
    pub fn band_of_cell(&self) -> &[u32] {
        &self.band_of_cell
    }

    /// Area-weighted mean of a cell field over each band. Sums in cell
    /// order, so the result is deterministic.
    pub fn mean(&self, geometry: &MeshGeometry, values: &[f64]) -> Vec<f64> {
        assert_eq!(values.len(), self.band_of_cell.len());
        let mut sums = vec![0.0; self.len()];
        for ((&b, &v), &a) in self
            .band_of_cell
            .iter()
            .zip(values)
            .zip(&geometry.cell_areas)
        {
            sums[b as usize] += a * v;
        }
        sums.iter().zip(&self.areas).map(|(s, a)| s / a).collect()
    }
}

/// Eastward and northward components of tangent vectors at points `p`
/// (zero at the poles, where east is undefined).
pub fn east_north(p: DVec3, v: DVec3) -> (f64, f64) {
    let east = DVec3::Z.cross(p);
    let len = east.length();
    if len < 1e-12 {
        return (0.0, 0.0);
    }
    let east = east / len;
    let north = p.cross(east);
    (v.dot(east), v.dot(north))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coords::from_lat_lon;
    use crate::mesh::LloydOptions;

    #[test]
    fn bands_cover_the_sphere_and_average_exactly() {
        let (m, _) = SphereMesh::centroidal(8, LloydOptions::default());
        let g = MeshGeometry::new(&m);
        let bands = LatitudeBands::new(&m, &g, 12);
        let total: f64 = bands.areas.iter().sum();
        assert!((total - 4.0 * std::f64::consts::PI).abs() < 1e-9);
        // A constant averages to itself; z = sin(lat) is antisymmetric.
        assert!(
            bands
                .mean(&g, &vec![2.5; m.n_cells()])
                .iter()
                .all(|v| (v - 2.5).abs() < 1e-12)
        );
        let z: Vec<f64> = m.cell_centers().iter().map(|p| p.z).collect();
        let means = bands.mean(&g, &z);
        let lats = bands.latitudes();
        for (m, lat) in means.iter().zip(&lats) {
            assert!((m - fms_math::sin(*lat)).abs() < 0.15, "{m} at {lat}");
        }
    }

    #[test]
    fn east_and_north_components() {
        let p = from_lat_lon(0.3, 1.1);
        let east = DVec3::Z.cross(p).normalize();
        let north = p.cross(east);
        let (e, n) = east_north(p, 3.0 * east - 2.0 * north);
        assert!((e - 3.0).abs() < 1e-12 && (n + 2.0).abs() < 1e-12);
        assert_eq!(east_north(DVec3::Z, DVec3::X), (0.0, 0.0));
    }
}
