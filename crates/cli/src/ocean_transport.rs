//! Northward ocean heat transport from the heating each transport leaves in
//! the ocean's columns: the paper's Eq. 8 (Bhatnagar et al. 2026) as a sum
//! over cells.

use fms_grid::LatitudeBands;

/// Northward heat transport (W) across the southern edge of each band and
/// the north pole, south to north: `Q(φ) = −Σ heating × area` over the cells
/// whose centres lie in the bands south of `φ`. `heating` is per unit of
/// ocean area (W m⁻²) and `ocean_area` each cell's ocean (m²). The first
/// value, at the south pole, is 0; the last is minus the total heating,
/// which is zero for a transport that only moves heat. Sums in cell order,
/// then band order, so the result is deterministic.
pub fn northward_transport(bands: &LatitudeBands, heating: &[f64], ocean_area: &[f64]) -> Vec<f64> {
    assert_eq!(heating.len(), bands.band_of_cell().len());
    assert_eq!(ocean_area.len(), heating.len());
    let mut gained = vec![0.0; bands.len()];
    for ((&b, &q), &a) in bands.band_of_cell().iter().zip(heating).zip(ocean_area) {
        gained[b as usize] += q * a;
    }
    let mut transport = Vec::with_capacity(bands.len() + 1);
    let mut q = 0.0;
    transport.push(q);
    for g in gained {
        q -= g;
        transport.push(q);
    }
    transport
}

/// Latitudes (degrees) of the band edges [`northward_transport`] reports at,
/// from −90 to 90.
pub fn edge_latitudes(bands: &LatitudeBands) -> Vec<f64> {
    let width = 180.0 / bands.len() as f64;
    (0..=bands.len())
        .map(|k| -90.0 + k as f64 * width)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use fms_grid::{LloydOptions, MeshGeometry, SphereMesh};

    const RADIUS: f64 = 6.371e6;

    fn grid(k: u32) -> (SphereMesh, MeshGeometry) {
        let (mesh, _) = SphereMesh::centroidal(k, LloydOptions::default());
        let geometry = MeshGeometry::new(&mesh);
        (mesh, geometry)
    }

    #[test]
    fn a_source_and_a_sink_carry_their_heat_between_them() {
        // A cell near 50° S heated and one near 30° N cooled by the same
        // power: the ocean carries it north across every latitude between
        // them and nowhere else.
        let (mesh, geometry) = grid(8);
        let bands = LatitudeBands::new(&mesh, &geometry, 18);
        let area: Vec<f64> = geometry
            .cell_areas
            .iter()
            .map(|a| a * RADIUS * RADIUS)
            .collect();
        let nearest = |z: f64| {
            (0..mesh.n_cells())
                .min_by(|&a, &b| {
                    let d = |i: usize| (mesh.cell_centers()[i].z - z).abs();
                    d(a).total_cmp(&d(b))
                })
                .unwrap()
        };
        let (source, sink) = (nearest(-0.766), nearest(0.5));
        let power = 1e15;
        let mut heating = vec![0.0; mesh.n_cells()];
        heating[source] = -power / area[source];
        heating[sink] = power / area[sink];
        let q = northward_transport(&bands, &heating, &area);
        let (from, to) = (bands.band_of_cell()[source], bands.band_of_cell()[sink]);
        for (k, (q, lat)) in q.iter().zip(edge_latitudes(&bands)).enumerate() {
            let expected = if k as u32 > from && k as u32 <= to {
                power
            } else {
                0.0
            };
            assert!((q - expected).abs() < 1e-6 * power, "{lat}°: {q}");
        }
    }

    #[test]
    fn a_heating_in_sin_latitude_gives_the_analytic_transport() {
        // Heating c·sin φ (W m⁻², all ocean) gives Q(φ) = −∫ c sin φ 2πR²
        // cos φ dφ from the south pole = π R² c cos² φ. The bands take whole
        // cells, so a band edge sits up to half a cell from where the sum
        // splits.
        let (mesh, geometry) = grid(16);
        let bands = LatitudeBands::new(&mesh, &geometry, 18);
        let area: Vec<f64> = geometry
            .cell_areas
            .iter()
            .map(|a| a * RADIUS * RADIUS)
            .collect();
        let c = 50.0;
        let heating: Vec<f64> = mesh.cell_centers().iter().map(|p| c * p.z).collect();
        let q = northward_transport(&bands, &heating, &area);
        let peak = std::f64::consts::PI * RADIUS * RADIUS * c;
        let mut worst = 0.0f64;
        for (q, lat) in q.iter().zip(edge_latitudes(&bands)) {
            let cos = fms_math::cos(lat.to_radians());
            worst = worst.max((q - peak * cos * cos).abs() / peak);
        }
        // Measured 0.8% of the peak; the closure at the north pole is exact
        // to round-off.
        eprintln!("largest error {worst:.4} of the peak {:.3} PW", peak / 1e15);
        assert!(worst < 0.02, "{worst}");
        let closure = q.last().unwrap() / peak;
        assert!(closure.abs() < 1e-12, "{closure}");
    }
}
