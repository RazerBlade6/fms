//! Print mesh size, build time and quality for geodesic vs. centroidal meshes.
//!
//! `cargo run --release -p fms-geo --example mesh_stats -- <frequency> [tolerance]`

use std::time::Instant;

use fms_geo::{LloydOptions, MeshGeometry, SphereMesh};

fn main() {
    let mut args = std::env::args().skip(1);
    let k: u32 = args.next().and_then(|s| s.parse().ok()).unwrap_or(45);
    let tolerance: f64 = args.next().and_then(|s| s.parse().ok()).unwrap_or(1e-6);

    let geodesic = SphereMesh::geodesic(k);
    let q0 = geodesic.quality();
    println!("frequency {k}: {} cells", geodesic.n_cells());
    println!(
        "  geodesic:   area ratio {:.4}, max centroid offset {:.3e}",
        q0.area_ratio(),
        q0.max_centroid_offset
    );

    let t = Instant::now();
    let options = LloydOptions {
        tolerance,
        ..LloydOptions::default()
    };
    let (mesh, report) = SphereMesh::centroidal(k, options);
    let elapsed = t.elapsed();
    let q1 = mesh.quality();
    println!(
        "  centroidal: area ratio {:.4}, max centroid offset {:.3e}",
        q1.area_ratio(),
        q1.max_centroid_offset
    );
    println!("  {report:?}, total {elapsed:.2?}");
    println!("  delaunay: {}", mesh.is_delaunay());

    let t = Instant::now();
    let geometry = MeshGeometry::new(&mesh);
    println!(
        "  geometry + TRiSK weights: {:.2?} ({} edges)",
        t.elapsed(),
        geometry.edge_lengths.len()
    );
}
