//! Topology and geometry invariants of the icosahedral Voronoi mesh.

use std::f64::consts::PI;

use fms_geo::{LloydOptions, SphereMesh};

const FREQUENCIES: [u32; 5] = [1, 2, 5, 8, 16];

#[test]
fn exactly_twelve_pentagons_rest_hexagons() {
    for k in FREQUENCIES {
        let mesh = SphereMesh::geodesic(k);
        let pentagons = (0..mesh.n_cells())
            .filter(|&c| mesh.cell_neighbors(c as u32).len() == 5)
            .count();
        let hexagons = (0..mesh.n_cells())
            .filter(|&c| mesh.cell_neighbors(c as u32).len() == 6)
            .count();
        assert_eq!(pentagons, 12, "frequency {k}");
        assert_eq!(pentagons + hexagons, mesh.n_cells(), "frequency {k}");
    }
}

#[test]
fn euler_characteristic_is_two() {
    for k in FREQUENCIES {
        let mesh = SphereMesh::geodesic(k);
        let chi = mesh.n_vertices() as i64 - mesh.n_edges() as i64 + mesh.n_cells() as i64;
        assert_eq!(chi, 2, "frequency {k}");
        assert_eq!(mesh.n_cells(), 10 * (k as usize).pow(2) + 2);
    }
}

#[test]
fn cell_areas_sum_to_sphere_area() {
    for k in FREQUENCIES {
        let mesh = SphereMesh::geodesic(k);
        let total: f64 = mesh.cell_areas().iter().sum();
        assert!((total - 4.0 * PI).abs() < 1e-11, "frequency {k}: {total}");
    }
}

#[test]
fn adjacency_is_symmetric_and_consistent_with_edges() {
    let mesh = SphereMesh::geodesic(8);
    for a in 0..mesh.n_cells() as u32 {
        for &b in mesh.cell_neighbors(a) {
            assert!(
                mesh.cell_neighbors(b).contains(&a),
                "{a} -> {b} not symmetric"
            );
        }
    }
    let degree_sum: usize = (0..mesh.n_cells() as u32)
        .map(|c| mesh.cell_neighbors(c).len())
        .sum();
    assert_eq!(degree_sum, 2 * mesh.n_edges());
    for e in mesh.edges() {
        let [a, b] = e.cells;
        assert!(a < b);
        // Both end vertices are shared by the two cells.
        for v in e.vertices {
            assert!(mesh.cell_vertices(a).contains(&v) && mesh.cell_vertices(b).contains(&v));
        }
    }
}

#[test]
fn edges_have_cell_zero_on_the_left() {
    let mesh = SphereMesh::geodesic(5);
    let c = mesh.cell_centers();
    let v = mesh.vertices();
    for e in mesh.edges() {
        let (from, to) = (v[e.vertices[0] as usize], v[e.vertices[1] as usize]);
        let left = c[e.cells[0] as usize];
        assert!(from.cross(to).dot(left) > 0.0);
    }
}

#[test]
fn geodesic_mesh_is_delaunay() {
    for k in FREQUENCIES {
        assert!(SphereMesh::geodesic(k).is_delaunay(), "frequency {k}");
    }
}

#[test]
fn lloyd_converges_to_centroidal_and_stays_valid() {
    let geodesic = SphereMesh::geodesic(16);
    let before = geodesic.quality();
    let (mesh, report) = SphereMesh::centroidal(
        16,
        LloydOptions {
            max_iterations: 2000,
            tolerance: 1e-6,
        },
    );
    let after = mesh.quality();

    assert!(report.converged, "{report:?}");
    assert!(after.max_centroid_offset < 1e-5, "{after:?}");
    assert!(after.max_centroid_offset < before.max_centroid_offset);
    assert!(
        after.area_ratio() < before.area_ratio(),
        "before {before:?}, after {after:?}"
    );
    assert!(mesh.is_delaunay());

    let total: f64 = mesh.cell_areas().iter().sum();
    assert!((total - 4.0 * PI).abs() < 1e-11);
    for p in mesh.cell_centers() {
        assert!((p.length() - 1.0).abs() < 1e-14);
    }
}

#[test]
fn multilevel_lloyd_converges_and_stays_valid() {
    // 40 > COARSEST_LEVEL, so this optimizes frequency 20 and then refines to 40.
    let (mesh, report) = SphereMesh::centroidal(40, LloydOptions::default());
    assert_eq!(mesh.n_cells(), 10 * 40 * 40 + 2);
    assert!(report.converged, "{report:?}");
    assert!(mesh.quality().max_centroid_offset < 1e-5);
    assert!(mesh.is_delaunay());
    let total: f64 = mesh.cell_areas().iter().sum();
    assert!((total - 4.0 * PI).abs() < 1e-10);
}
