//! Where the TRiSK divergence errors are largest: the hexagons next to the 12
//! pentagons, where the maximum error stays fixed as resolution increases
//! (OBSERVATIONS.md O-1).
//!
//! ```sh
//! cargo run --release -p fms-geo --example divergence_error
//! ```
use fms_geo::{LloydOptions, MeshGeometry, Operators, SphereMesh};
use glam::DVec3;
fn t(p: DVec3, g: DVec3) -> DVec3 {
    g - p * p.dot(g)
}
fn main() {
    for k in [16u32, 32, 64] {
        let (m, _) = SphereMesh::centroidal(k, LloydOptions::default());
        let g = MeshGeometry::new(&m);
        let ops = Operators::new(&m, &g, 1.0);
        let vel = |p: DVec3| t(p, DVec3::new(2.0 * p.x, -2.0 * p.y, 0.3));
        let u: Vec<f64> = g
            .edge_points
            .iter()
            .zip(&g.edge_normals)
            .map(|(&x, &n)| vel(x).dot(n))
            .collect();
        let mut div = vec![0.0; m.n_cells()];
        ops.divergence(&u, &mut div);
        let pent: Vec<DVec3> = (0..m.n_cells() as u32)
            .filter(|&c| m.cell_neighbors(c).len() == 5)
            .map(|c| m.cell_centers()[c as usize])
            .collect();
        let mut errs: Vec<(f64, f64)> = m
            .cell_centers()
            .iter()
            .zip(&div)
            .map(|(&p, &d)| {
                let exact = -6.0 * (p.x * p.x - p.y * p.y) - 0.6 * p.z;
                let dist = pent
                    .iter()
                    .map(|&q| fms_math::acos(p.dot(q).clamp(-1.0, 1.0)))
                    .fold(9.0, f64::min)
                    / m.mean_spacing();
                ((d - exact).abs() / 6.0, dist)
            })
            .collect();
        errs.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        let far: f64 = errs
            .iter()
            .filter(|e| e.1 > 3.0)
            .map(|e| e.0)
            .fold(0.0, f64::max);
        println!(
            "k={k}: top-5 errors (rel, cells from pentagon): {:?}; max beyond 3 cells: {far:.2e}",
            errs[..5]
                .iter()
                .map(|e| (format!("{:.2e}", e.0), format!("{:.1}", e.1)))
                .collect::<Vec<_>>()
        );
    }
}
