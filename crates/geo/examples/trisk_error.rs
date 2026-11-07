//! TRiSK tangential-velocity reconstruction error for solid-body rotation.
//!
//! `cargo run --release -p fms-geo --example trisk_error -- 16 32 64`

use fms_geo::{LloydOptions, MeshGeometry, SphereMesh};
use glam::DVec3;

fn main() {
    let omega = DVec3::new(0.3, -0.2, 1.0).normalize();
    for k in std::env::args()
        .skip(1)
        .filter_map(|s| s.parse::<u32>().ok())
    {
        let (m, _) = SphereMesh::centroidal(k, LloydOptions::default());
        let g = MeshGeometry::new(&m);
        let velocity = |x: DVec3| omega.cross(x);
        let u: Vec<f64> = (0..m.n_edges())
            .map(|e| velocity(g.edge_points[e]).dot(g.edge_normals[e]))
            .collect();
        let (mut max, mut sq) = (0.0f64, 0.0);
        for e in 0..m.n_edges() as u32 {
            let (nb, w) = g.trisk_stencil(e);
            let v: f64 = nb
                .iter()
                .zip(w)
                .map(|(&f, &w)| w * g.edge_lengths[f as usize] * u[f as usize])
                .sum::<f64>()
                / g.dual_edge_lengths[e as usize];
            let err = v - velocity(g.edge_points[e as usize]).dot(g.edge_tangents[e as usize]);
            max = max.max(err.abs());
            sq += err * err;
        }
        let rms = (sq / m.n_edges() as f64).sqrt();
        println!(
            "k={k:4} cells={:8} max={max:.3e} rms={rms:.3e}",
            m.n_cells()
        );
    }
}
