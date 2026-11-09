//! Icosahedral geodesic grids, two-level mapping, projections and picking.
//!
//! The core type is [`SphereMesh`]: a spherical centroidal Voronoi tessellation
//! (12 pentagons, the rest hexagons) built by subdividing an icosahedron and
//! optimizing it with Lloyd's algorithm. All geometry is on the unit sphere.
//!
//! See `docs/PLAN.md` for the design.

pub mod coords;
pub mod geometry;
pub mod hilbert;
pub mod icosahedron;
pub mod locate;
pub mod mesh;
pub mod operators;
pub mod spherical;
pub mod two_level;
pub mod zonal;

pub use geometry::MeshGeometry;
pub use locate::CellLocator;
pub use mesh::{
    CellId, Edge, EdgeId, LloydOptions, LloydReport, MeshQuality, SphereMesh, VertexId,
};
pub use operators::Operators;
pub use two_level::TwoLevelMap;
pub use zonal::LatitudeBands;
