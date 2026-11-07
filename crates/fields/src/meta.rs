//! Where a field lives and how to describe it.

use fms_geo::SphereMesh;

use crate::column::DType;

/// Which of the two grids a field lives on.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum GridLevel {
    /// The fine terrain grid (editing, hydrology, rendering).
    Terrain,
    /// The coarse climate grid (atmosphere and ocean).
    Climate,
}

/// Which mesh element a field has one value per.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Element {
    /// Voronoi cells: scalars such as temperature or elevation.
    Cell,
    /// Voronoi edges: normal velocities and fluxes.
    Edge,
    /// Voronoi vertices (dual triangles): vorticity.
    Vertex,
}

/// Grid level plus element kind: determines a field's length.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Location {
    /// The grid.
    pub grid: GridLevel,
    /// The element kind.
    pub element: Element,
}

impl Location {
    /// One value per terrain cell.
    pub const TERRAIN_CELLS: Self = Self {
        grid: GridLevel::Terrain,
        element: Element::Cell,
    };
    /// One value per climate cell.
    pub const CLIMATE_CELLS: Self = Self {
        grid: GridLevel::Climate,
        element: Element::Cell,
    };
    /// One value per climate edge.
    pub const CLIMATE_EDGES: Self = Self {
        grid: GridLevel::Climate,
        element: Element::Edge,
    };
    /// One value per climate vertex.
    pub const CLIMATE_VERTICES: Self = Self {
        grid: GridLevel::Climate,
        element: Element::Vertex,
    };
}

/// Element counts of one mesh.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MeshCounts {
    /// Number of cells.
    pub cells: usize,
    /// Number of edges.
    pub edges: usize,
    /// Number of vertices.
    pub vertices: usize,
}

impl MeshCounts {
    /// Counts of `mesh`.
    pub fn of(mesh: &SphereMesh) -> Self {
        Self {
            cells: mesh.n_cells(),
            edges: mesh.n_edges(),
            vertices: mesh.n_vertices(),
        }
    }

    /// Count for one element kind.
    pub fn get(&self, element: Element) -> usize {
        match element {
            Element::Cell => self.cells,
            Element::Edge => self.edges,
            Element::Vertex => self.vertices,
        }
    }
}

/// Element counts of both grids.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GridSizes {
    /// The terrain grid.
    pub terrain: MeshCounts,
    /// The climate grid.
    pub climate: MeshCounts,
}

impl GridSizes {
    /// Number of values a field at `location` holds.
    pub fn len(&self, location: Location) -> usize {
        match location.grid {
            GridLevel::Terrain => self.terrain.get(location.element),
            GridLevel::Climate => self.climate.get(location.element),
        }
    }
}

/// Suggested colour scale for displaying a field. Renderers interpret it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ColorScale {
    /// Monotonic ramp for magnitudes (default).
    Sequential,
    /// Two-sided ramp around a midpoint (anomalies, signed velocities).
    Diverging,
    /// Distinct colours per integer value (classes, ids).
    Categorical,
    /// Hypsometric tint: ocean blues below zero, land colours above.
    Elevation,
    /// Cold-to-hot ramp.
    Temperature,
    /// Dry-to-wet ramp.
    Moisture,
}

/// Whether a field has one buffer or a front/back pair.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Buffering {
    /// One buffer, read and written in place.
    Single,
    /// Front (read) and back (write) buffers, exchanged with
    /// [`FieldRegistry::swap`](crate::FieldRegistry::swap). Lets a system write the
    /// next state while others (or a UI snapshot) read the current one.
    Double,
}

/// Description of a field, supplied at registration.
#[derive(Clone, Debug, PartialEq)]
pub struct FieldSpec {
    /// Unique name of the form `namespace.field`, both parts `snake_case`, where
    /// the namespace is the owning system or crate, e.g. `"atmos.surface_temperature"`.
    pub name: &'static str,
    /// One-line human description.
    pub description: &'static str,
    /// SI unit symbol, e.g. `"K"`, `"m"`, `"m s-1"`; empty for dimensionless.
    pub units: &'static str,
    /// Where the field lives.
    pub location: Location,
    /// Physically valid range, if any. Used for display limits and to detect
    /// numerical blow-ups.
    pub valid_range: Option<(f64, f64)>,
    /// Suggested colour scale.
    pub color_scale: ColorScale,
    /// Single or double buffered.
    pub buffering: Buffering,
    /// Whether the field is part of the saved state. Derived fields that can be
    /// recomputed from others are not.
    pub persistent: bool,
}

impl FieldSpec {
    /// A single-buffered, persistent field with no range and a sequential scale.
    /// Adjust the rest with struct-update syntax.
    pub const fn new(name: &'static str, units: &'static str, location: Location) -> Self {
        Self {
            name,
            description: "",
            units,
            location,
            valid_range: None,
            color_scale: ColorScale::Sequential,
            buffering: Buffering::Single,
            persistent: true,
        }
    }
}

/// A registered field's full metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct FieldMeta {
    /// The registration spec.
    pub spec: FieldSpec,
    /// Element type.
    pub dtype: DType,
    /// Number of values.
    pub len: usize,
}
