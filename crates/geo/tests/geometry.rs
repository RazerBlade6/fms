//! Metric terms, mimetic operator identities and TRiSK properties.

use std::f64::consts::PI;
use std::sync::OnceLock;

use fms_geo::{LloydOptions, MeshGeometry, SphereMesh};
use glam::DVec3;

/// A centroidal, Hilbert-ordered mesh shared by the tests (built once).
fn mesh(k: u32) -> &'static (SphereMesh, MeshGeometry) {
    static M16: OnceLock<(SphereMesh, MeshGeometry)> = OnceLock::new();
    static M32: OnceLock<(SphereMesh, MeshGeometry)> = OnceLock::new();
    let cell = match k {
        16 => &M16,
        32 => &M32,
        _ => panic!("unsupported test frequency"),
    };
    cell.get_or_init(|| {
        let (m, _) = SphereMesh::centroidal(k, LloydOptions::default());
        let g = MeshGeometry::new(&m);
        (m, g)
    })
}

/// Deterministic pseudo-random values in [-1, 1] (no RNG dependency needed).
fn noise(n: usize, seed: u64) -> Vec<f64> {
    let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    (0..n)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            (x >> 11) as f64 / (1u64 << 52) as f64 - 1.0
        })
        .collect()
}

#[test]
fn edge_frames_are_orthonormal_and_oriented() {
    let (m, g) = mesh(16);
    for (e, edge) in m.edges().iter().enumerate() {
        let (x, n, t) = (g.edge_points[e], g.edge_normals[e], g.edge_tangents[e]);
        assert!((x.length() - 1.0).abs() < 1e-14);
        assert!((n.length() - 1.0).abs() < 1e-14 && (t.length() - 1.0).abs() < 1e-14);
        assert!(x.dot(n).abs() < 1e-14 && x.dot(t).abs() < 1e-14 && n.dot(t).abs() < 1e-14);
        let to_cell1 =
            m.cell_centers()[edge.cells[1] as usize] - m.cell_centers()[edge.cells[0] as usize];
        assert!(n.dot(to_cell1) > 0.0);
        let along =
            m.vertices()[edge.vertices[1] as usize] - m.vertices()[edge.vertices[0] as usize];
        assert!(
            t.dot(along) > 0.0,
            "tangent must point from vertices[0] to vertices[1]"
        );
    }
}

#[test]
fn kites_partition_cells_and_sphere() {
    let (m, g) = mesh(16);
    for c in 0..m.n_cells() as u32 {
        let kites: f64 = g.cell_kites(c).iter().sum();
        let area = g.cell_areas[c as usize];
        assert!(
            (kites - area).abs() < 1e-12 * area.max(1.0),
            "cell {c}: {kites} vs {area}"
        );
    }
    let dual_total: f64 = g.vertex_areas.iter().sum();
    assert!((dual_total - 4.0 * PI).abs() < 1e-10);
}

#[test]
fn lengths_are_consistent() {
    let (m, g) = mesh(16);
    // For an orthogonal Voronoi–Delaunay pair, Σ ½·lₑ·dₑ approximates the sphere
    // area (exact in the plane).
    let total: f64 = (0..m.n_edges())
        .map(|e| 0.5 * g.edge_lengths[e] * g.dual_edge_lengths[e])
        .sum();
    assert!((total / (4.0 * PI) - 1.0).abs() < 1e-3, "{total}");
    assert!(g.edge_lengths.iter().all(|&l| l > 0.0));
}

#[test]
fn curl_of_gradient_is_zero() {
    let (m, g) = mesh(16);
    let h = noise(m.n_cells(), 1);
    let grad: Vec<f64> = m
        .edges()
        .iter()
        .enumerate()
        .map(|(e, edge)| {
            (h[edge.cells[1] as usize] - h[edge.cells[0] as usize]) / g.dual_edge_lengths[e]
        })
        .collect();
    for v in 0..m.n_vertices() as u32 {
        let signs = g.vertex_edge_signs(v);
        let circulation: f64 = m
            .vertex_edges(v)
            .iter()
            .zip(signs)
            .map(|(&e, s)| s * g.dual_edge_lengths[e as usize] * grad[e as usize])
            .sum();
        assert!(circulation.abs() < 1e-13, "vertex {v}: {circulation}");
    }
}

#[test]
fn divergence_of_stream_function_flow_is_zero() {
    let (m, g) = mesh(16);
    let psi = noise(m.n_vertices(), 2);
    // Normal velocity from a vertex stream function: u = −(ψ₁ − ψ₀)/lₑ.
    let u: Vec<f64> = m
        .edges()
        .iter()
        .enumerate()
        .map(|(e, edge)| {
            -(psi[edge.vertices[1] as usize] - psi[edge.vertices[0] as usize]) / g.edge_lengths[e]
        })
        .collect();
    for c in 0..m.n_cells() as u32 {
        let flux: f64 = m
            .cell_edges(c)
            .iter()
            .zip(g.cell_edge_signs(c))
            .map(|(&e, s)| s * g.edge_lengths[e as usize] * u[e as usize])
            .sum();
        assert!(flux.abs() < 1e-13, "cell {c}: {flux}");
    }
}

#[test]
fn trisk_weights_are_antisymmetric() {
    let (m, g) = mesh(16);
    for e in 0..m.n_edges() as u32 {
        let (neighbors, weights) = g.trisk_stencil(e);
        for (&f, &w) in neighbors.iter().zip(weights) {
            let (back, back_w) = g.trisk_stencil(f);
            let pos = back
                .iter()
                .position(|&x| x == e)
                .expect("stencils are symmetric");
            assert!(
                (w + back_w[pos]).abs() < 1e-14,
                "w[{e},{f}] = {w}, w[{f},{e}] = {}",
                back_w[pos]
            );
        }
    }
}

fn reconstruct_tangential(m: &SphereMesh, g: &MeshGeometry, u: &[f64]) -> Vec<f64> {
    (0..m.n_edges() as u32)
        .map(|e| {
            let (neighbors, weights) = g.trisk_stencil(e);
            let sum: f64 = neighbors
                .iter()
                .zip(weights)
                .map(|(&f, &w)| w * g.edge_lengths[f as usize] * u[f as usize])
                .sum();
            sum / g.dual_edge_lengths[e as usize]
        })
        .collect()
}

/// (max, rms) error of the TRiSK tangential reconstruction of solid-body
/// rotation with unit maximum speed.
fn solid_body_error(k: u32) -> (f64, f64) {
    let (m, g) = mesh(k);
    let omega = DVec3::new(0.3, -0.2, 1.0).normalize();
    let velocity = |x: DVec3| omega.cross(x);
    let u: Vec<f64> = (0..m.n_edges())
        .map(|e| velocity(g.edge_points[e]).dot(g.edge_normals[e]))
        .collect();
    let v = reconstruct_tangential(m, g, &u);
    let errors: Vec<f64> = (0..m.n_edges())
        .map(|e| v[e] - velocity(g.edge_points[e]).dot(g.edge_tangents[e]))
        .collect();
    let max = errors.iter().fold(0.0f64, |a, e| a.max(e.abs()));
    let rms = (errors.iter().map(|e| e * e).sum::<f64>() / errors.len() as f64).sqrt();
    (max, rms)
}

#[test]
fn trisk_reconstructs_solid_body_rotation() {
    // TRiSK's tangential reconstruction is first-order in RMS but does not
    // converge in max-norm on geodesic grids (Peixoto 2016): local errors near
    // grid irregularities stay at a fixed level. We check both behaviours, plus
    // the sign (a sign error would give errors near 2).
    let ((max16, rms16), (max32, rms32)) = (solid_body_error(16), solid_body_error(32));
    println!("k=16 max {max16:.3e} rms {rms16:.3e}; k=32 max {max32:.3e} rms {rms32:.3e}");
    assert!(max16 < 0.01 && max32 < 0.01, "max error too large");
    assert!(
        rms32 < rms16 / 1.7,
        "RMS error should be ~first order: {rms16} -> {rms32}"
    );
}

#[test]
fn trisk_maps_divergence_to_dual_cells() {
    // Thuburn et al. (2009): the divergence of the reconstructed perpendicular
    // flux on each dual triangle equals the kite-weighted average of the primal
    // divergences of its three cells. This is what keeps geostrophic modes steady.
    let (m, g) = mesh(16);
    let f = noise(m.n_edges(), 3);
    let perp = reconstruct_tangential(m, g, &f);
    let div: Vec<f64> = (0..m.n_cells() as u32)
        .map(|c| {
            let flux: f64 = m
                .cell_edges(c)
                .iter()
                .zip(g.cell_edge_signs(c))
                .map(|(&e, s)| s * g.edge_lengths[e as usize] * f[e as usize])
                .sum();
            flux / g.cell_areas[c as usize]
        })
        .collect();
    for v in 0..m.n_vertices() as u32 {
        // Outward across dual edge e from triangle v is +tₑ at vertices[0], −tₑ at vertices[1].
        let dual_div: f64 = m
            .vertex_edges(v)
            .iter()
            .map(|&e| {
                let sign = if m.edges()[e as usize].vertices[0] == v {
                    1.0
                } else {
                    -1.0
                };
                sign * g.dual_edge_lengths[e as usize] * perp[e as usize]
            })
            .sum();
        let kites = g.vertex_kites(v);
        let averaged: f64 = m.triangles()[v as usize]
            .iter()
            .zip(kites)
            .map(|(&c, kite)| {
                kite / g.cell_areas[c as usize] * div[c as usize] * g.cell_areas[c as usize]
            })
            .sum();
        assert!(
            (dual_div - averaged).abs() < 1e-12,
            "vertex {v}: {dual_div} vs {averaged}"
        );
    }
}

#[test]
fn hilbert_order_improves_locality_and_keeps_topology() {
    let (geodesic, _) = (SphereMesh::geodesic(32), ());
    let ordered = geodesic.hilbert_ordered();
    let spread = |m: &SphereMesh| -> f64 {
        let total: u64 = (0..m.n_cells() as u32)
            .flat_map(|c| {
                m.cell_neighbors(c)
                    .iter()
                    .map(move |&n| u64::from(c.abs_diff(n)))
            })
            .sum();
        total as f64 / (2 * m.n_edges()) as f64
    };
    let (before, after) = (spread(&geodesic), spread(&ordered));
    println!("mean neighbour index distance: {before:.1} -> {after:.1}");
    assert!(after < before / 2.0);
    assert_eq!(ordered.n_cells(), geodesic.n_cells());
    assert_eq!(ordered.n_edges(), geodesic.n_edges());
    assert!(ordered.is_delaunay());
    let total: f64 = ordered.cell_areas().iter().sum();
    assert!((total - 4.0 * PI).abs() < 1e-11);
}

/// (max, rms) error of Perot reconstruction of solid-body rotation, relative to
/// the maximum speed (1).
fn perot_error(k: u32) -> (f64, f64) {
    let (m, g) = mesh(k);
    let omega = DVec3::new(0.3, -0.2, 1.0).normalize();
    let u: Vec<f64> = (0..m.n_edges())
        .map(|e| omega.cross(g.edge_points[e]).dot(g.edge_normals[e]))
        .collect();
    let v = g.reconstruct_vectors(m, &u);
    let errors: Vec<f64> = m
        .cell_centers()
        .iter()
        .zip(&v)
        .map(|(&x, &r)| (r - omega.cross(x)).length())
        .collect();
    let max = errors.iter().fold(0.0f64, |a, &e| a.max(e));
    let rms = (errors.iter().map(|e| e * e).sum::<f64>() / errors.len() as f64).sqrt();
    (max, rms)
}

#[test]
fn perot_reconstructs_tangent_vectors() {
    let (m, g) = mesh(16);
    // Results are tangent to the sphere.
    let u = noise(m.n_edges(), 5);
    for (x, v) in m.cell_centers().iter().zip(g.reconstruct_vectors(m, &u)) {
        assert!(x.dot(v).abs() < 1e-12);
    }
    let ((max16, rms16), (max32, rms32)) = (perot_error(16), perot_error(32));
    println!("perot: k=16 max {max16:.3e} rms {rms16:.3e}; k=32 max {max32:.3e} rms {rms32:.3e}");
    // Max error does not shrink (largest at pentagons), as with TRiSK.
    assert!(max16 < 0.03 && max32 < 0.03);
    assert!(rms32 < rms16, "should converge: {rms16} -> {rms32}");
}

#[test]
fn averaging_onto_cells_preserves_constants() {
    let (m, g) = mesh(16);
    assert!(
        g.vertices_to_cells(m, &vec![2.5; m.n_vertices()])
            .iter()
            .all(|&v| (v - 2.5).abs() < 1e-12)
    );
    assert!(
        g.edges_to_cells(m, &vec![-1.0; m.n_edges()])
            .iter()
            .all(|&v| (v + 1.0).abs() < 1e-12)
    );
}
