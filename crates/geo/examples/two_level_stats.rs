//! Time the fine ↔ coarse mapping and point location at realistic sizes.
//!
//! `cargo run --release -p fms-geo --example two_level_stats -- [coarse] [fine]`

use std::time::Instant;

use fms_geo::coords::from_lat_lon;
use fms_geo::{CellLocator, LloydOptions, MeshGeometry, SphereMesh, TwoLevelMap};
use rayon::prelude::*;

fn main() {
    let mut args = std::env::args()
        .skip(1)
        .filter_map(|s| s.parse::<u32>().ok());
    let (kc, kf) = (args.next().unwrap_or(45), args.next().unwrap_or(400));

    let (coarse, _) = SphereMesh::centroidal(kc, LloydOptions::default());
    let (fine, _) = SphereMesh::centroidal(kf, LloydOptions::default());
    let fine_areas = MeshGeometry::new(&fine).cell_areas;
    println!(
        "coarse {} cells, fine {} cells",
        coarse.n_cells(),
        fine.n_cells()
    );

    let t = Instant::now();
    let map = TwoLevelMap::new(&coarse, &fine, &fine_areas);
    println!("TwoLevelMap::new: {:.2?}", t.elapsed());
    let kids: Vec<usize> = (0..coarse.n_cells() as u32)
        .map(|c| map.children(c).len())
        .collect();
    println!(
        "children per coarse cell: min {} max {} mean {:.1}",
        kids.iter().min().unwrap(),
        kids.iter().max().unwrap(),
        fine.n_cells() as f64 / coarse.n_cells() as f64
    );

    let field: Vec<f64> = fine.cell_centers().iter().map(|p| p.x).collect();
    let mut out = vec![0.0; coarse.n_cells()];
    let t = Instant::now();
    map.restrict_mean(&field, &mut out);
    println!("restrict_mean: {:.2?}", t.elapsed());
    let mut back = vec![0.0; fine.n_cells()];
    let t = Instant::now();
    map.prolong(&out, &mut back);
    println!("prolong: {:.2?}", t.elapsed());

    let locator = CellLocator::new(&fine);
    let (w, h) = (4096usize, 2048usize);
    let t = Instant::now();
    let raster: Vec<u32> = (0..w * h)
        .into_par_iter()
        .map(|i| {
            let (x, y) = (i % w, i / w);
            let lat =
                std::f64::consts::FRAC_PI_2 - std::f64::consts::PI * (y as f64 + 0.5) / h as f64;
            let lon = -std::f64::consts::PI + std::f64::consts::TAU * (x as f64 + 0.5) / w as f64;
            locator.locate(&fine, from_lat_lon(lat, lon))
        })
        .collect();
    println!(
        "rasterize {w}x{h} on fine mesh: {:.2?} ({} px)",
        t.elapsed(),
        raster.len()
    );
}
