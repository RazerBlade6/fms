//! Resampling a raster onto a mesh: the area-weighted mean of the pixels in
//! each cell.
//!
//! A pixel belongs to the cell whose centre is nearest its own centre, which is
//! exactly the Voronoi cell containing it. Each cell tests only the pixels in
//! the bounding box of a cap around it, against its own neighbours, so every
//! pixel is counted in exactly one cell (ties go to the lower cell index) and
//! the work runs in parallel with a fixed summation order per cell.
//!
//! When the raster is coarser than the mesh some cells contain no pixel centre;
//! those take a bilinear sample at their centre instead.

use std::f64::consts::FRAC_PI_2;

use fms_geo::coords::to_lat_lon;
use fms_geo::spherical::arc_length;
use fms_geo::{CellId, SphereMesh};
use glam::DVec3;
use rayon::prelude::*;

use crate::raster::Raster;

/// The area-weighted mean of `raster` over each cell of `mesh`.
pub fn resample(raster: &Raster, mesh: &SphereMesh) -> Vec<f64> {
    let centers = mesh.cell_centers();
    cell_sums(raster, mesh)
        .into_iter()
        .zip(centers)
        .map(|((sum, weight), &p)| {
            if weight > 0.0 {
                sum / weight
            } else {
                raster.sample(p)
            }
        })
        .collect()
}

/// Per cell: `Σ w·z` and `Σ w` over the pixels it contains, with `w` the
/// pixel's area divided by the column spacing.
pub(crate) fn cell_sums(raster: &Raster, mesh: &SphereMesh) -> Vec<(f64, f64)> {
    let rows: Vec<(f64, f64, f64)> = (0..raster.height())
        .map(|r| {
            let (sin, cos) = fms_math::sin_cos(raster.lat(r));
            (sin, cos, raster.row_weight(r))
        })
        .collect();
    let cols: Vec<(f64, f64)> = (0..raster.width())
        .map(|c| fms_math::sin_cos(raster.lon(c)))
        .collect();
    let centers = mesh.cell_centers();
    let vertices = mesh.vertices();

    (0..mesh.n_cells() as CellId)
        .into_par_iter()
        .map(|cell| {
            let center = centers[cell as usize];
            let neighbors: Vec<(CellId, DVec3)> = mesh
                .cell_neighbors(cell)
                .iter()
                .map(|&n| (n, centers[n as usize]))
                .collect();
            // The cell lies inside the cap through its farthest vertex.
            let radius = mesh
                .cell_vertices(cell)
                .iter()
                .map(|&v| arc_length(center, vertices[v as usize]))
                .fold(0.0, f64::max)
                * (1.0 + 1e-9)
                + 1e-12;

            let (mut sum, mut weight) = (0.0, 0.0);
            for_each_pixel_in_cap(raster, center, radius, |row, col| {
                let (sin_lat, cos_lat, w) = rows[row];
                let (sin_lon, cos_lon) = cols[col];
                let p = DVec3::new(cos_lat * cos_lon, cos_lat * sin_lon, sin_lat);
                let d = p.dot(center);
                let inside = neighbors.iter().all(|&(n, q)| {
                    let dn = p.dot(q);
                    d > dn || (d == dn && cell < n)
                });
                if inside {
                    sum += w * raster.get(row, col);
                    weight += w;
                }
            });
            (sum, weight)
        })
        .collect()
}

/// Call `f(row, col)` for every pixel whose centre is in the latitude/longitude
/// bounding box of the spherical cap of angular `radius` around `center`.
fn for_each_pixel_in_cap(
    raster: &Raster,
    center: DVec3,
    radius: f64,
    mut f: impl FnMut(usize, usize),
) {
    let (lat, lon) = to_lat_lon(center);
    // Rows count southward, so the cap's north edge gives the first row.
    let first = (raster.row_of(lat + radius).ceil() as i64).max(0);
    let last = (raster.row_of(lat - radius).floor() as i64).min(raster.height() as i64 - 1);
    if first > last {
        return; // The cap fits between two rows of pixel centres.
    }
    let (first, last) = (first as usize, last as usize);

    let width = raster.width();
    let covers_pole = lat + radius >= FRAC_PI_2 || lat - radius <= -FRAC_PI_2;
    let ratio = fms_math::sin(radius) / fms_math::cos(lat);
    let all_columns = covers_pole || ratio >= 1.0;
    let (col_start, n_cols) = if all_columns {
        (0, width)
    } else {
        let half = fms_math::asin(ratio);
        // Unwrapped column indices; `col_of` is periodic so work relative to it.
        let x_center = raster.col_of(lon);
        let span = half / raster.spacing();
        let lo = (x_center - span).ceil() as i64;
        let hi = (x_center + span).floor() as i64;
        let n = (hi - lo + 1).max(0) as usize;
        (lo, n.min(width))
    };

    for row in first..=last {
        for i in 0..n_cols {
            let col = (col_start + i as i64).rem_euclid(width as i64) as usize;
            f(row, col);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fms_geo::LloydOptions;

    #[test]
    fn every_pixel_is_counted_exactly_once() {
        let (mesh, _) = SphereMesh::centroidal(12, LloydOptions::default());
        for (w, h) in [(360, 180), (97, 41), (1440, 720)] {
            let raster = Raster::from_fn(w, h, |_, _| 1.0);
            let total: f64 = (0..h).map(|r| raster.row_weight(r) * w as f64).sum();
            let counted: f64 = cell_sums(&raster, &mesh).iter().map(|s| s.1).sum();
            assert!(
                ((counted - total) / total).abs() < 1e-12,
                "{w}×{h}: counted {counted}, total {total}"
            );
        }
    }

    #[test]
    fn pixel_weight_in_each_cell_matches_its_area() {
        let (mesh, _) = SphereMesh::centroidal(12, LloydOptions::default());
        let raster = Raster::from_fn(1440, 720, |_, _| 1.0);
        let areas = mesh.cell_areas();
        let dlon = raster.spacing();
        let worst = cell_sums(&raster, &mesh)
            .iter()
            .zip(&areas)
            .map(|(&(_, w), &a)| (w * dlon / a - 1.0).abs())
            .fold(0.0, f64::max);
        // ~1000 pixels per cell; quantization error ~ perimeter / area.
        assert!(worst < 0.05, "worst relative area error {worst}");
    }
}
