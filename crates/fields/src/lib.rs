//! Field registry and struct-of-arrays per-element storage.
//!
//! A **field** is a named, typed array with one value per cell, edge or vertex of
//! one of the two grids (terrain or climate): `terrain.elevation`,
//! `atmos.surface_temperature`, `atmos.wind_normal_lower`, and later
//! `soil.fertility` or `civ.population`. Names are `namespace.field`, with the
//! owning system or crate as the namespace, so independent crates cannot collide.
//! Every layer of world state is a field, so map modes, saving and new systems all
//! work the same way for every layer.
//!
//! ```
//! use fms_fields::{FieldRegistry, FieldSpec, GridSizes, Location, MeshCounts};
//!
//! let counts = MeshCounts { cells: 12, edges: 30, vertices: 20 };
//! let mut fields = FieldRegistry::new(GridSizes { terrain: counts, climate: counts });
//! let elevation = fields
//!     .register(FieldSpec::new("terrain.elevation", "m", Location::TERRAIN_CELLS), 0.0f32)
//!     .unwrap();
//!
//! fields.write(elevation)[3] = 1200.0;
//! assert_eq!(fields.read(elevation)[3], 1200.0);
//! // Other crates find the field by name, with the type checked once.
//! assert_eq!(fields.id::<f32>("terrain.elevation").unwrap(), elevation);
//! ```
//!
//! See `docs/DESIGN.md` for the design.

mod column;
mod meta;
mod registry;

pub use column::{Column, DType, FieldType};
pub use meta::{
    Buffering, ColorScale, Element, FieldMeta, FieldSpec, GridLevel, GridSizes, Location,
    MeshCounts, Quantity,
};
pub use registry::{
    ColumnRef, FieldError, FieldId, FieldKey, FieldMut, FieldRef, FieldRegistry, validate_name,
};
