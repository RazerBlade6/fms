//! Point probes: read fields at one location, for time series.
//!
//! A probe picks the cell nearest a latitude/longitude on each field's grid
//! (it does not interpolate, so it reports what the model actually holds there)
//! and reads one value per field each time it is sampled. Non-cell fields are
//! resampled to that cell the same way the map renderer does, so a probe and a
//! map always agree:
//!
//! - vertex fields: kite-area average;
//! - edge scalars: edge-length average;
//! - edge normal components (winds, currents): Perot reconstruction, reported as
//!   two columns, the eastward and northward components.

use std::fmt;

use fms_fields::{Element, FieldKey, GridLevel, Quantity};
use fms_grid::CellId;
use fms_grid::coords::{from_lat_lon, to_lat_lon};
use glam::DVec3;

use crate::World;

/// What a probe column reports from its field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Component {
    /// The field's value, resampled to the cell if needed.
    Value,
    /// Eastward component of a vector field.
    East,
    /// Northward component of a vector field.
    North,
}

/// The cell a probe reads on one grid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Site {
    /// Which grid.
    pub grid: GridLevel,
    /// The cell nearest the requested point.
    pub cell: CellId,
    /// Latitude of the cell centre (degrees).
    pub lat_deg: f64,
    /// Longitude of the cell centre (degrees, −180 to 180).
    pub lon_deg: f64,
}

/// One output column of a probe.
#[derive(Clone, Debug)]
pub struct ProbeColumn {
    /// The field read.
    pub key: FieldKey,
    /// What is read from it.
    pub component: Component,
    /// Column name: the field name, with `.east` / `.north` for vector components.
    pub label: String,
    /// Units of the field.
    pub units: &'static str,
    /// Where it is read.
    pub site: Site,
}

/// Why a probe could not be made.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ProbeError {
    /// Latitude outside −90..90, or a non-finite coordinate.
    InvalidLocation {
        /// Requested latitude (degrees).
        lat_deg: f64,
        /// Requested longitude (degrees).
        lon_deg: f64,
    },
}

impl fmt::Display for ProbeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::InvalidLocation { lat_deg, lon_deg } => write!(
                f,
                "invalid probe location ({lat_deg}, {lon_deg}): latitude must be in \
                 -90..90 and both must be finite"
            ),
        }
    }
}

impl std::error::Error for ProbeError {}

/// Reads a set of fields at one location.
#[derive(Clone, Debug)]
pub struct Probe {
    lat_deg: f64,
    lon_deg: f64,
    columns: Vec<ProbeColumn>,
}

impl Probe {
    /// A probe at (`lat_deg`, `lon_deg`) reading `keys`, in order. Vector fields
    /// give two columns (east, then north); every other field gives one.
    pub fn new(
        world: &World,
        lat_deg: f64,
        lon_deg: f64,
        keys: &[FieldKey],
    ) -> Result<Self, ProbeError> {
        if !(lat_deg.is_finite() && lon_deg.is_finite() && (-90.0..=90.0).contains(&lat_deg)) {
            return Err(ProbeError::InvalidLocation { lat_deg, lon_deg });
        }
        let point = from_lat_lon(lat_deg.to_radians(), lon_deg.to_radians());
        let site = |grid_level| {
            let grid = world.grid(grid_level);
            let cell = grid.locator.locate(&grid.mesh, point);
            let (lat, lon) = to_lat_lon(grid.mesh.cell_centers()[cell as usize]);
            Site {
                grid: grid_level,
                cell,
                lat_deg: lat.to_degrees(),
                lon_deg: lon.to_degrees(),
            }
        };
        let mut columns = Vec::new();
        for &key in keys {
            let spec = &world.fields.meta(key).spec;
            let site = site(spec.location.grid);
            let column = |component, label| ProbeColumn {
                key,
                component,
                label,
                units: spec.units,
                site,
            };
            if spec.location.element == Element::Edge && spec.quantity == Quantity::NormalComponent
            {
                columns.push(column(Component::East, format!("{}.east", spec.name)));
                columns.push(column(Component::North, format!("{}.north", spec.name)));
            } else {
                columns.push(column(Component::Value, spec.name.to_owned()));
            }
        }
        Ok(Self {
            lat_deg,
            lon_deg,
            columns,
        })
    }

    /// The requested latitude (degrees).
    pub fn lat_deg(&self) -> f64 {
        self.lat_deg
    }

    /// The requested longitude (degrees).
    pub fn lon_deg(&self) -> f64 {
        self.lon_deg
    }

    /// The output columns, in order.
    pub fn columns(&self) -> &[ProbeColumn] {
        &self.columns
    }

    /// Read every column now, appending one value per column to `out`.
    pub fn sample_into(&self, world: &World, out: &mut Vec<f64>) {
        for c in &self.columns {
            let field = world.fields.read_column(c.key);
            let get = |i: u32| field.get_f64(i as usize);
            let grid = world.grid(c.site.grid);
            let (mesh, geometry) = (&grid.mesh, &grid.geometry);
            let spec = &world.fields.meta(c.key).spec;
            let value = match (spec.location.element, c.component) {
                (Element::Cell, _) => get(c.site.cell),
                (Element::Vertex, _) => geometry.vertex_average(mesh, c.site.cell, get),
                (Element::Edge, Component::Value) => geometry.edge_average(mesh, c.site.cell, get),
                (Element::Edge, component) => {
                    let v = geometry.reconstruct_vector(mesh, c.site.cell, get);
                    let (east, north) = local_frame(mesh.cell_centers()[c.site.cell as usize]);
                    match component {
                        Component::East => v.dot(east),
                        _ => v.dot(north),
                    }
                }
            };
            out.push(value);
        }
    }

    /// Read every column now.
    pub fn sample(&self, world: &World) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.columns.len());
        self.sample_into(world, &mut out);
        out
    }
}

/// Unit east and north vectors at `p`. At the poles, where east is undefined,
/// "east" is taken as the direction of longitude 90° E.
pub fn local_frame(p: DVec3) -> (DVec3, DVec3) {
    let k = DVec3::Z.cross(p);
    let east = if k.length_squared() > 1e-24 {
        k.normalize()
    } else {
        DVec3::Y
    };
    (east, p.cross(east).normalize())
}
