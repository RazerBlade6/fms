//! Terrain: elevation on the terrain grid, and the fields derived from it on
//! both grids (slope, and the climate grid's sub-grid statistics).
//!
//! Elevation is the saved state; everything else is recomputed by
//! [`World::set_elevation`] whenever it changes.

use std::path::Path;
use std::time::Instant;

use fms_fields::{ColorScale, FieldError, FieldId, FieldRegistry, FieldSpec, Location};
use fms_worldgen::subgrid::{slope, smooth_height, surface_height};
use fms_worldgen::{Continent, HeightmapError, LoadOptions, TerrainStats};

use crate::World;

/// Elevation of a world with no terrain yet: a uniform ocean floor (m), about
/// Earth's mean ocean depth.
pub const BLANK_OCEAN_ELEVATION: f64 = -4000.0;

/// Pixels per terrain-cell spacing to keep when reading a heightmap. Finer
/// files are box-averaged down to this while loading; ~16 pixels per cell is
/// plenty for a good cell mean.
const PIXELS_PER_CELL_SPACING: f64 = 4.0;

/// Smoothing passes for the terrain the atmosphere sees. On Earth at 10k
/// cells, ten bring the persistent σ-coordinate wind error from 5.7 to
/// 0.87 m/s, near the 0.58 of forty, for 9% of the highest cell's height
/// ([`smooth_height`]; `docs/OBSERVATIONS.md` O-30).
const SURFACE_SMOOTHING_PASSES: usize = 10;

/// Handles to the terrain fields.
#[derive(Clone, Copy, Debug)]
pub struct TerrainFields {
    /// Height above sea level of each terrain cell (m); negative is sea floor.
    pub elevation: FieldId<f64>,
    /// Surface slope of each terrain cell (m per m).
    pub slope: FieldId<f64>,
    /// Land fraction of each climate cell.
    pub land_fraction: FieldId<f64>,
    /// Mean surface height of each climate cell (m; sea level over ocean).
    pub surface_height: FieldId<f64>,
    /// Sub-grid standard deviation of surface height in each climate cell (m).
    pub surface_height_stddev: FieldId<f64>,
    /// Mean slope within each climate cell (m per m).
    pub mean_slope: FieldId<f64>,
}

/// What [`World::load_heightmap`] read, and how long it took.
#[derive(Clone, Copy, Debug)]
pub struct HeightmapImport {
    /// Raster size after any averaging while reading (columns, rows).
    pub size: (usize, usize),
    /// Reading and decoding the file (s).
    pub read_seconds: f64,
    /// Resampling onto the terrain grid and recomputing derived fields (s).
    pub resample_seconds: f64,
}

pub(crate) fn register(fields: &mut FieldRegistry) -> Result<TerrainFields, FieldError> {
    let derived = |name, description, units, location| FieldSpec {
        description,
        persistent: false,
        ..FieldSpec::new(name, units, location)
    };
    let slope_spec = |name, description, location| FieldSpec {
        valid_range: Some((0.0, 10.0)),
        ..derived(name, description, "m m-1", location)
    };
    Ok(TerrainFields {
        elevation: fields.register(
            FieldSpec {
                description: "Height above sea level; negative is sea floor",
                valid_range: Some((-12_000.0, 10_000.0)),
                color_scale: ColorScale::Elevation,
                ..FieldSpec::new("terrain.elevation", "m", Location::TERRAIN_CELLS)
            },
            BLANK_OCEAN_ELEVATION,
        )?,
        slope: fields.register(
            slope_spec(
                "terrain.slope",
                "Surface slope (sea surface over ocean)",
                Location::TERRAIN_CELLS,
            ),
            0.0,
        )?,
        land_fraction: fields.register(
            FieldSpec {
                valid_range: Some((0.0, 1.0)),
                ..derived(
                    "terrain.land_fraction",
                    "Fraction of the cell that is land",
                    "",
                    Location::CLIMATE_CELLS,
                )
            },
            0.0,
        )?,
        surface_height: fields.register(
            FieldSpec {
                valid_range: Some((0.0, 10_000.0)),
                color_scale: ColorScale::Elevation,
                ..derived(
                    "terrain.surface_height",
                    "Mean surface height (sea level over ocean)",
                    "m",
                    Location::CLIMATE_CELLS,
                )
            },
            0.0,
        )?,
        surface_height_stddev: fields.register(
            FieldSpec {
                valid_range: Some((0.0, 5_000.0)),
                ..derived(
                    "terrain.surface_height_stddev",
                    "Sub-grid spread of surface height",
                    "m",
                    Location::CLIMATE_CELLS,
                )
            },
            0.0,
        )?,
        mean_slope: fields.register(
            slope_spec(
                "terrain.mean_slope",
                "Mean surface slope of the terrain cells inside",
                Location::CLIMATE_CELLS,
            ),
            0.0,
        )?,
    })
}

impl World {
    /// Replace the terrain-grid elevation (m) and recompute every field
    /// derived from it.
    ///
    /// # Panics
    /// If `elevation` doesn't have one value per terrain cell.
    pub fn set_elevation(&mut self, elevation: Vec<f64>) {
        assert_eq!(
            elevation.len(),
            self.terrain.mesh.n_cells(),
            "one elevation per terrain cell"
        );
        self.fields
            .write(self.core.terrain.elevation)
            .copy_from_slice(&elevation);
        self.update_terrain();
    }

    /// Read a global heightmap (GeoTIFF, or greyscale PNG with `png_range` the
    /// elevations of black and white) and make it the terrain.
    pub fn load_heightmap(
        &mut self,
        path: &Path,
        png_range: Option<(f64, f64)>,
    ) -> Result<HeightmapImport, HeightmapError> {
        let t = Instant::now();
        let options = LoadOptions {
            png_range,
            target_spacing: Some(self.terrain.mesh.mean_spacing() / PIXELS_PER_CELL_SPACING),
        };
        let raster = fms_worldgen::load(path, &options)?;
        let read_seconds = t.elapsed().as_secs_f64();
        let t = Instant::now();
        self.set_elevation(fms_worldgen::resample(&raster, &self.terrain.mesh));
        Ok(HeightmapImport {
            size: (raster.width(), raster.height()),
            read_seconds,
            resample_seconds: t.elapsed().as_secs_f64(),
        })
    }

    /// Make the terrain the default idealized [`Continent`].
    pub fn set_idealized_continent(&mut self) {
        let elevation = Continent::default().elevation(self.terrain.mesh.cell_centers());
        self.set_elevation(elevation);
    }

    /// Surface geopotential of each climate cell (m² s⁻²): gravity times the
    /// mean surface height, smoothed so the σ-coordinate dynamics can carry
    /// it. Exactly zero on a world with no land.
    pub fn surface_geopotential(&self) -> Vec<f64> {
        let gravity = self.planet.params().gravity;
        smooth_height(
            &self.fields.read(self.core.terrain.surface_height),
            &self.climate.mesh,
            &self.climate.geometry,
            SURFACE_SMOOTHING_PASSES,
        )
        .iter()
        .map(|z| gravity * z)
        .collect()
    }

    /// Recompute slope and the climate grid's terrain statistics from the
    /// current elevation.
    pub(crate) fn update_terrain(&mut self) {
        let f = self.core.terrain;
        let elevation = self.fields.read(f.elevation).to_vec();
        let radius = self.planet.params().radius;
        let s = slope(
            &surface_height(&elevation),
            &self.terrain.mesh,
            &self.terrain.geometry,
            radius,
        );
        let stats = TerrainStats::compute(&elevation, &s, &self.two_level);
        self.fields.write(f.slope).copy_from_slice(&s);
        for (id, values) in [
            (f.land_fraction, &stats.land_fraction),
            (f.surface_height, &stats.surface_height),
            (f.surface_height_stddev, &stats.surface_height_stddev),
            (f.mean_slope, &stats.mean_slope),
        ] {
            self.fields.write(id).copy_from_slice(values);
        }
    }
}
