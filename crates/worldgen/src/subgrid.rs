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
