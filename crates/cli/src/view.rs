//! Turning any field into something drawable: a scalar per cell, plus vectors
//! for fields that are normal components of a vector field.

use fms_fields::{Element, FieldKey, Quantity};
use fms_sim::World;
use glam::DVec3;

/// A field resampled to cell centres.
pub struct View {
    /// One scalar per cell (vector fields: the speed).
    pub values: Vec<f64>,
    /// One tangent vector per cell, for vector fields.
    pub vectors: Option<Vec<DVec3>>,
}

/// Resample `key` to the cells of its grid.
///
/// - cell fields are used as they are;
/// - vertex fields are kite-area averaged onto cells;
/// - edge scalars are length-weighted averages over each cell's edges;
/// - edge normal components are reconstructed into vectors (Perot), shown as
///   speed.
pub fn view(world: &World, key: FieldKey) -> View {
    let spec = &world.fields.meta(key).spec;
    let grid = world.grid(spec.location.grid);
    let raw = world.fields.read_column(key).to_f64_vec();
    match (spec.location.element, spec.quantity) {
        (Element::Cell, _) => View {
            values: raw,
            vectors: None,
        },
        (Element::Vertex, _) => View {
            values: grid.geometry.vertices_to_cells(&grid.mesh, &raw),
            vectors: None,
        },
        (Element::Edge, Quantity::Scalar) => View {
            values: grid.geometry.edges_to_cells(&grid.mesh, &raw),
            vectors: None,
        },
        (Element::Edge, Quantity::NormalComponent) => {
            let vectors = grid.geometry.reconstruct_vectors(&grid.mesh, &raw);
            View {
                values: vectors.iter().map(|v| v.length()).collect(),
                vectors: Some(vectors),
            }
        }
    }
}

/// Whether a field can provide vectors for an arrow overlay.
pub fn is_vector(world: &World, key: FieldKey) -> bool {
    let spec = &world.fields.meta(key).spec;
    spec.location.element == Element::Edge && spec.quantity == Quantity::NormalComponent
}

/// Legend title for a field, noting how non-cell fields were resampled.
pub fn title(world: &World, key: FieldKey) -> String {
    let spec = &world.fields.meta(key).spec;
    let how = match (spec.location.element, spec.quantity) {
        (Element::Cell, _) => "",
        (Element::Vertex, _) => " AT CELLS",
        (Element::Edge, Quantity::Scalar) => " AT CELLS",
        (Element::Edge, Quantity::NormalComponent) => " SPEED",
    };
    if spec.units.is_empty() {
        format!("{}{how}", spec.name)
    } else {
        format!("{}{how} ({})", spec.name, spec.units)
    }
}
