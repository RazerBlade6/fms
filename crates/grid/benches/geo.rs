//! Grid construction, geometry, fine/coarse transfers and point location.
//!
//! `cargo bench -p fms-grid` (compare against a saved baseline with
//! `-- --save-baseline before` / `-- --baseline before`).

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use fms_grid::coords::from_lat_lon;
use fms_grid::{CellLocator, LloydOptions, MeshGeometry, SphereMesh, TwoLevelMap};

const COARSE: u32 = 45; // 20k cells, the default climate grid
const FINE: u32 = 160; // 256k cells, the CLI's default terrain grid

fn build(c: &mut Criterion) {
    let mut g = c.benchmark_group("build");
    g.sample_size(10);
    g.bench_function("geodesic_f160", |b| {
        b.iter(|| SphereMesh::geodesic(black_box(FINE)))
    });
    g.bench_function("centroidal_f45", |b| {
        b.iter(|| SphereMesh::centroidal(black_box(COARSE), LloydOptions::default()))
    });
    let (mesh, _) = SphereMesh::centroidal(FINE, LloydOptions::default());
    g.bench_function("geometry_f160", |b| {
        b.iter(|| MeshGeometry::new(black_box(&mesh)))
    });
    g.bench_function("hilbert_order_f160", |b| {
        b.iter(|| black_box(&mesh).hilbert_ordered())
    });
    g.finish();
}

fn transfers(c: &mut Criterion) {
    let (coarse, _) = SphereMesh::centroidal(COARSE, LloydOptions::default());
    let (fine, _) = SphereMesh::centroidal(FINE, LloydOptions::default());
    let fine_geo = MeshGeometry::new(&fine);
    let coarse_geo = MeshGeometry::new(&coarse);

    let mut g = c.benchmark_group("two_level_f160_to_f45");
    g.sample_size(20);
    g.bench_function("build_map", |b| {
        b.iter(|| TwoLevelMap::new(&coarse, &fine, black_box(&fine_geo.cell_areas)))
    });
    let map = TwoLevelMap::new(&coarse, &fine, &fine_geo.cell_areas);
    let fine_field: Vec<f64> = fine.cell_centers().iter().map(|p| p.x).collect();
    let mut coarse_field = vec![0.0; coarse.n_cells()];
    let mut back = vec![0.0; fine.n_cells()];
    g.bench_function("restrict_mean", |b| {
        b.iter(|| map.restrict_mean(black_box(&fine_field), &mut coarse_field))
    });
    g.bench_function("restrict_sum", |b| {
        b.iter(|| map.restrict_sum(black_box(&fine_field), &mut coarse_field))
    });
    g.bench_function("prolong", |b| {
        b.iter(|| map.prolong(black_box(&coarse_field), &mut back))
    });
    g.finish();

    let mut g = c.benchmark_group("operators_f45");
    let u: Vec<f64> = (0..coarse.n_edges())
        .map(|e| (e as f64 * 0.1).fract())
        .collect();
    g.bench_function("perot_reconstruction", |b| {
        b.iter(|| coarse_geo.reconstruct_vectors(&coarse, black_box(&u)))
    });
    g.finish();
}

fn locate(c: &mut Criterion) {
    let (fine, _) = SphereMesh::centroidal(FINE, LloydOptions::default());
    let locator = CellLocator::new(&fine);
    let points: Vec<_> = (0..10_000)
        .map(|i| {
            let f = f64::from(i);
            from_lat_lon(
                (f * 0.618_034).fract() * 3.0 - 1.5,
                (f * 0.414_214).fract() * 6.2 - 3.1,
            )
        })
        .collect();
    let mut g = c.benchmark_group("locate_f160");
    g.bench_function("cells_10k_points", |b| {
        b.iter(|| {
            points
                .iter()
                .map(|&p| locator.locate(&fine, p))
                .fold(0u64, |a, c| a + u64::from(c))
        })
    });
    g.bench_function("triangles_10k_points", |b| {
        b.iter(|| {
            points
                .iter()
                .map(|&p| locator.locate_triangle(&fine, p).0)
                .fold(0u64, |a, c| a + u64::from(c))
        })
    });
    g.finish();
}

criterion_group!(benches, build, transfers, locate);
criterion_main!(benches);
