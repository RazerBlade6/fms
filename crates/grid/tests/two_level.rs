//! Point location and fine ↔ coarse grid mapping.

use std::f64::consts::PI;
use std::sync::OnceLock;

use fms_grid::{CellLocator, LloydOptions, MeshGeometry, SphereMesh, TwoLevelMap};
use glam::DVec3;

struct Grids {
    coarse: SphereMesh,
    fine: SphereMesh,
    coarse_areas: Vec<f64>,
    fine_areas: Vec<f64>,
    map: TwoLevelMap,
}

/// Coarse frequency 10 (1002 cells) and fine frequency 40 (16002 cells).
fn grids() -> &'static Grids {
    static GRIDS: OnceLock<Grids> = OnceLock::new();
    GRIDS.get_or_init(|| {
        let (coarse, _) = SphereMesh::centroidal(10, LloydOptions::default());
        let (fine, _) = SphereMesh::centroidal(40, LloydOptions::default());
        let coarse_areas = MeshGeometry::new(&coarse).cell_areas;
        let fine_areas = MeshGeometry::new(&fine).cell_areas;
        let map = TwoLevelMap::new(&coarse, &fine, &fine_areas);
        Grids {
            coarse,
            fine,
            coarse_areas,
            fine_areas,
            map,
        }
    })
}

fn random_points(n: usize, seed: u64) -> Vec<DVec3> {
    let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        (x >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..n)
        .map(|_| {
            // Uniform on the sphere.
            let z = 2.0 * next() - 1.0;
            let phi = 2.0 * PI * next();
            let r = (1.0 - z * z).sqrt();
            DVec3::new(r * fms_math::cos(phi), r * fms_math::sin(phi), z)
        })
        .collect()
}

fn brute_force_nearest(mesh: &SphereMesh, p: DVec3) -> u32 {
    let centers = mesh.cell_centers();
    (0..mesh.n_cells())
        .max_by(|&a, &b| p.dot(centers[a]).total_cmp(&p.dot(centers[b])))
        .unwrap() as u32
}

#[test]
fn locate_matches_brute_force() {
    let g = grids();
    for mesh in [&g.coarse, &g.fine] {
        let locator = CellLocator::new(mesh);
        for p in random_points(2000, 7) {
            assert_eq!(locator.locate(mesh, p), brute_force_nearest(mesh, p));
        }
        // Poles and the date line are bucket-table edge cases.
        for p in [
            DVec3::Z,
            -DVec3::Z,
            DVec3::new(-1.0, 1e-12, 0.0).normalize(),
            DVec3::new(-1.0, -1e-12, 0.0).normalize(),
        ] {
            assert_eq!(locator.locate(mesh, p), brute_force_nearest(mesh, p));
        }
    }
}

#[test]
fn locate_cell_centers_returns_the_cell() {
    let mesh = &grids().fine;
    let locator = CellLocator::new(mesh);
    for (c, &p) in mesh.cell_centers().iter().enumerate() {
        assert_eq!(locator.locate(mesh, p), c as u32);
    }
}

#[test]
fn locate_triangle_contains_point() {
    let mesh = &grids().coarse;
    let locator = CellLocator::new(mesh);
    for p in random_points(2000, 11) {
        let (t, w) = locator.locate_triangle(mesh, p);
        assert!(w.iter().all(|&x| x >= 0.0), "{w:?}");
        assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        // Weights reproduce the (gnomonically projected) point.
        let [a, b, c] = mesh.triangles()[t as usize].map(|i| mesh.cell_centers()[i as usize]);
        let q = a * w[0] + b * w[1] + c * w[2];
        assert!(q.normalize().dot(p) > 1.0 - 1e-12);
    }
}

#[test]
fn children_partition_fine_cells() {
    let g = grids();
    let mut seen = vec![false; g.fine.n_cells()];
    for c in 0..g.coarse.n_cells() as u32 {
        let kids = g.map.children(c);
        assert!(!kids.is_empty());
        for &f in kids {
            assert!(!seen[f as usize]);
            seen[f as usize] = true;
            assert_eq!(g.map.parents()[f as usize], c);
        }
        // Child areas approximate the parent's area (cells are assigned whole).
        let child_area: f64 = kids.iter().map(|&f| g.fine_areas[f as usize]).sum();
        let ratio = child_area / g.coarse_areas[c as usize];
        assert!(
            (0.75..1.25).contains(&ratio),
            "coarse cell {c}: area ratio {ratio}"
        );
    }
    assert!(seen.iter().all(|&s| s));
}

#[test]
fn restriction_conserves_and_preserves_constants() {
    let g = grids();
    let mut coarse = vec![0.0; g.coarse.n_cells()];

    g.map
        .restrict_mean(&vec![3.5; g.fine.n_cells()], &mut coarse);
    assert!(coarse.iter().all(|&x| (x - 3.5).abs() < 1e-12));

    let water: Vec<f64> = g
        .fine
        .cell_centers()
        .iter()
        .map(|p| 1.0 + p.x * p.y)
        .collect();
    g.map.restrict_sum(&water, &mut coarse);
    let (fine_total, coarse_total) = (water.iter().sum::<f64>(), coarse.iter().sum::<f64>());
    assert!((fine_total - coarse_total).abs() < 1e-9 * fine_total);

    let mut variance = vec![0.0; g.coarse.n_cells()];
    g.map
        .restrict_mean(&vec![2.0; g.fine.n_cells()], &mut coarse);
    g.map
        .restrict_variance(&vec![2.0; g.fine.n_cells()], &coarse, &mut variance);
    assert!(variance.iter().all(|&v| v.abs() < 1e-24));
}

#[test]
fn restrict_mean_approximates_smooth_fields() {
    let g = grids();
    let field = |p: DVec3| fms_math::sin(2.0 * p.x) + p.y * p.z;
    let fine: Vec<f64> = g.fine.cell_centers().iter().map(|&p| field(p)).collect();
    let mut coarse = vec![0.0; g.coarse.n_cells()];
    g.map.restrict_mean(&fine, &mut coarse);
    let max_err = g
        .coarse
        .cell_centers()
        .iter()
        .zip(&coarse)
        .map(|(&p, &v)| (field(p) - v).abs())
        .fold(0.0, f64::max);
    // Coarse spacing is ~0.12 rad; field gradients are O(2).
    assert!(max_err < 0.05, "{max_err}");
}

#[test]
fn prolong_preserves_constants_and_is_accurate() {
    let g = grids();
    let mut fine = vec![0.0; g.fine.n_cells()];
    g.map.prolong(&vec![-1.25; g.coarse.n_cells()], &mut fine);
    assert!(fine.iter().all(|&x| (x + 1.25).abs() < 1e-12));

    let field = |p: DVec3| fms_math::sin(2.0 * p.x) + p.y * p.z;
    let coarse: Vec<f64> = g.coarse.cell_centers().iter().map(|&p| field(p)).collect();
    g.map.prolong(&coarse, &mut fine);
    let max_err = g
        .fine
        .cell_centers()
        .iter()
        .zip(&fine)
        .map(|(&p, &v)| (field(p) - v).abs())
        .fold(0.0, f64::max);
    // Second-order in coarse spacing (~0.12 rad).
    assert!(max_err < 0.02, "{max_err}");
}
