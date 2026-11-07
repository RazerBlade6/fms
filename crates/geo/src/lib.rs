//! Icosahedral geodesic grids, two-level mapping, projections and picking.
//!
//! The core type is [`SphereMesh`]: a spherical centroidal Voronoi tessellation
//! (12 pentagons, the rest hexagons) built by subdividing an icosahedron and
//! optimizing it with Lloyd's algorithm. All geometry is on the unit sphere.
//!
//! See `docs/PLAN.md` for the design.

pub mod geometry;
pub mod hilbert;
pub mod icosahedron;
pub mod mesh;
pub mod spherical;

pub use geometry::MeshGeometry;
pub use mesh::{
    CellId, Edge, EdgeId, LloydOptions, LloydReport, MeshQuality, SphereMesh, VertexId,
};
