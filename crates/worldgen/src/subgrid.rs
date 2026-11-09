//! Terrain statistics the climate grid needs about terrain it can't resolve.
//!
//! The atmosphere sits on the *surface height*: the land surface, or sea level
//! over ocean (`max(elevation, 0)`). Each climate cell gets the area-weighted
//! mean of that over its terrain children, its spread (sub-grid orography, for
//! drag and orographic rain), the mean slope, and the land fraction.

use fms_geo::{MeshGeometry, SphereMesh, TwoLevelMap};
use rayon::prelude::*;

/// Surface height (m): the land surface, or sea level (0) over ocean.
pub fn surface_height(elevation: &[f64]) -> Vec<f64> {
    elevation.par_iter().map(|&z| z.max(0.0)).collect()
}

/// Whether a terrain cell is land: above sea level.
pub fn is_land(elevation: f64) -> bool {
    elevation > 0.0
}

/// Magnitude of the surface-height gradient in each cell (m per m), from the
/// TRiSK edge gradients reconstructed to a cell vector (Perot). `radius` is
/// the planet radius in metres.
pub fn slope(surface: &[f64], mesh: &SphereMesh, geometry: &MeshGeometry, radius: f64) -> Vec<f64> {
    let gradient: Vec<f64> = mesh
        .edges()
        .par_iter()
        .zip(&geometry.dual_edge_lengths)
        .map(|(e, &d)| (surface[e.cells[1] as usize] - surface[e.cells[0] as usize]) / (d * radius))
        .collect();
    geometry
        .reconstruct_vectors(mesh, &gradient)
        .par_iter()
        .map(|v| v.length())
        .collect()
}

/// A few grid-scale diffusion passes over a climate-grid `height` (any
/// units), conserving its area mean: the terrain the σ-coordinate dynamics
/// can carry. Steep cell-to-cell steps leave a persistent error in the
/// pressure gradient (see `docs/OBSERVATIONS.md` O-30). Each pass diffuses
/// by a fixed fraction of the smallest cell's area, so the filter's width is
/// a fixed number of cells at any resolution, and it is monotone: heights
/// stay within their original range, and a flat surface stays exactly flat.
pub fn smooth_height(
    height: &[f64],
    mesh: &SphereMesh,
    geometry: &MeshGeometry,
    passes: usize,
) -> Vec<f64> {
    const DIFFUSION: f64 = 0.08;
    let smallest = geometry
        .cell_areas
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let mut h = height.to_vec();
    let mut net = vec![0.0; h.len()];
    for _ in 0..passes {
        net.fill(0.0);
        for (e, edge) in mesh.edges().iter().enumerate() {
            let [a, b] = edge.cells.map(|c| c as usize);
            let flux = geometry.edge_lengths[e] / geometry.dual_edge_lengths[e] * (h[b] - h[a]);
            net[a] += flux;
            net[b] -= flux;
        }
        for ((h, n), area) in h.iter_mut().zip(&net).zip(&geometry.cell_areas) {
            *h += DIFFUSION * smallest / area * n;
        }
    }
    h
}

/// Sub-grid terrain statistics, one value per climate cell.
#[derive(Clone, Debug, PartialEq)]
pub struct TerrainStats {
    /// Fraction of the cell's area that is land (0–1).
    pub land_fraction: Vec<f64>,
    /// Mean surface height (m).
    pub surface_height: Vec<f64>,
    /// Standard deviation of surface height within the cell (m).
    pub surface_height_stddev: Vec<f64>,
    /// Mean slope magnitude of the terrain cells within the cell (m per m).
    pub mean_slope: Vec<f64>,
}

impl TerrainStats {
    /// Aggregate terrain-grid `elevation` and its `slope` (from [`slope`]) to
    /// the climate grid.
    pub fn compute(elevation: &[f64], slope: &[f64], map: &TwoLevelMap) -> Self {
        let n = map.n_coarse();
        let restrict = |fine: &[f64]| {
            let mut coarse = vec![0.0; n];
            map.restrict_mean(fine, &mut coarse);
            coarse
        };
        let surface = surface_height(elevation);
        let land: Vec<f64> = elevation
            .par_iter()
            .map(|&z| if is_land(z) { 1.0 } else { 0.0 })
            .collect();

        let surface_height = restrict(&surface);
        let mut variance = vec![0.0; n];
        map.restrict_variance(&surface, &surface_height, &mut variance);
        Self {
            // Clamp away round-off in the weighted sum of 0s and 1s.
            land_fraction: restrict(&land).iter().map(|f| f.clamp(0.0, 1.0)).collect(),
            surface_height_stddev: variance.iter().map(|v| v.max(0.0).sqrt()).collect(),
            surface_height,
            mean_slope: restrict(slope),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fms_geo::LloydOptions;

    fn grid() -> (SphereMesh, MeshGeometry) {
        let (mesh, _) = SphereMesh::centroidal(8, LloydOptions::default());
        let geometry = MeshGeometry::new(&mesh);
        (mesh, geometry)
    }

    #[test]
    fn smoothing_conserves_the_mean_and_stays_in_range() {
        let (mesh, geometry) = grid();
        let height: Vec<f64> = (0..mesh.n_cells())
            .map(|i| if i % 7 == 0 { 3000.0 } else { 0.0 })
            .collect();
        let smooth = smooth_height(&height, &mesh, &geometry, 10);
        let mean = |h: &[f64]| -> f64 {
            h.iter()
                .zip(&geometry.cell_areas)
                .map(|(x, a)| x * a)
                .sum::<f64>()
        };
        assert!((mean(&smooth) / mean(&height) - 1.0).abs() < 1e-13);
        let peak = smooth.iter().copied().fold(f64::MIN, f64::max);
        let lowest = smooth.iter().copied().fold(f64::MAX, f64::min);
        assert!(peak < 3000.0, "peak {peak}");
        assert!(lowest >= 0.0, "lowest {lowest}");
    }

    #[test]
    fn a_flat_surface_stays_exactly_flat() {
        let (mesh, geometry) = grid();
        let flat = vec![0.0; mesh.n_cells()];
        assert_eq!(smooth_height(&flat, &mesh, &geometry, 10), flat);
    }
}
