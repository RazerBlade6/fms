//! Terrain: elevation on the terrain grid, and the fields derived from it on
//! both grids (slope and the sea mask, and the climate grid's sub-grid
//! statistics).
//!
//! Elevation is the saved state; everything else is recomputed by
//! [`World::set_elevation`] whenever it changes.

use std::path::Path;
use std::time::Instant;

use fms_fields::{ColorScale, FieldError, FieldId, FieldRegistry, FieldSpec, Location};
use fms_worldgen::subgrid::{slope, smooth_height, surface_height};
use fms_worldgen::{Continent, HeightmapError, LoadOptions, SeaMask, TerrainStats};

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
pub const SURFACE_SMOOTHING_PASSES: usize = 10;

/// Handles to the terrain fields.
#[derive(Clone, Copy, Debug)]
pub struct TerrainFields {
    /// Height above sea level of each terrain cell (m); negative is sea floor.
    pub elevation: FieldId<f64>,
    /// Surface slope of each terrain cell (m per m).
    pub slope: FieldId<f64>,
    /// Whether each terrain cell is sea (1) or land (0), from [`SeaMask`]:
    /// land below sea level away from the ocean is land.
    pub sea: FieldId<f64>,
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
        sea: fields.register(
            FieldSpec {
                valid_range: Some((0.0, 1.0)),
                color_scale: ColorScale::Categorical,
                ..derived(
                    "terrain.sea",
                    "1 where the cell is sea, 0 where it is land (including land below sea level)",
                    "",
                    Location::TERRAIN_CELLS,
                )
            },
            1.0,
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
        self.surface_geopotential_smoothed(SURFACE_SMOOTHING_PASSES)
    }

    /// [`surface_geopotential`](Self::surface_geopotential) with `passes`
    /// smoothing passes instead of the default, for experiments; 0 is the
    /// plain cell means.
    pub fn surface_geopotential_smoothed(&self, passes: usize) -> Vec<f64> {
        let gravity = self.planet.params().gravity;
        smooth_height(
            &self.fields.read(self.core.terrain.surface_height),
            &self.climate.mesh,
            &self.climate.geometry,
            passes,
        )
        .iter()
        .map(|z| gravity * z)
        .collect()
    }

    /// Recompute slope, the sea mask and the climate grid's terrain
    /// statistics from the current elevation.
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
        let sea = SeaMask::new(
            &elevation,
            &self.terrain.mesh,
            &self.terrain.geometry.cell_areas,
        )
        .sea;
        let stats = TerrainStats::compute(&elevation, &sea, &s, &self.two_level);
        self.fields.write(f.slope).copy_from_slice(&s);
        for (out, &sea) in self.fields.write(f.sea).iter_mut().zip(&sea) {
            *out = if sea { 1.0 } else { 0.0 };
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::WorldConfig;
    use fms_grid::coords::from_lat_lon;
    use fms_worldgen::MIN_SEA_FRACTION;
    use glam::DVec3;

    fn world(terrain_frequency: u32, climate_frequency: u32, radius: f64) -> World {
        let mut config = WorldConfig {
            terrain_frequency,
            climate_frequency,
            ..WorldConfig::default()
        };
        config.planet.radius = radius;
        World::new(config)
    }

    fn hash(values: &[f64]) -> u64 {
        values.iter().fold(0xcbf2_9ce4_8422_2325u64, |h, v| {
            (h ^ v.to_bits()).wrapping_mul(0x0100_0000_01b3)
        })
    }

    #[test]
    fn the_blank_ocean_and_the_continent_are_as_before_the_sea_mask_bit_for_bit() {
        // Hashes of `terrain.land_fraction` (blank ocean, then the idealized
        // continent) and of the continent's other climate-grid statistics,
        // taken before the sea mask existed, when land was `elevation > 0`.
        const BEFORE: [(u32, u32, [u64; 5]); 2] = [
            (
                20,
                8,
                [
                    0xb877_cbdc_903d_31ed,
                    0x2518_e739_7557_d01c,
                    0xdc91_5ce8_fa65_39a8,
                    0xe9a2_232a_b0c9_0e50,
                    0x013c_64ef_8368_c7ea,
                ],
            ),
            (
                64,
                16,
                [
                    0xdc0e_1f65_c1aa_77ed,
                    0xc7ae_3e7b_e31b_71ae,
                    0x489c_6ced_db80_b77c,
                    0xa377_8449_8eb2_d0d7,
                    0x6ae5_2bdf_e449_0a2b,
                ],
            ),
        ];
        for (terrain, climate, before) in BEFORE {
            let mut w = world(terrain, climate, 6.371e6);
            let t = w.core.terrain;
            assert!(w.fields.read(t.sea).iter().all(|&s| s == 1.0));
            let blank = hash(&w.fields.read(t.land_fraction));
            w.set_idealized_continent();
            let mut now = vec![blank];
            for id in [
                t.land_fraction,
                t.surface_height,
                t.surface_height_stddev,
                t.mean_slope,
            ] {
                now.push(hash(&w.fields.read(id)));
            }
            assert_eq!(now, before, "{terrain}/{climate}: {now:#018x?}");
            let elevation = w.fields.read(t.elevation);
            for (s, z) in w.fields.read(t.sea).iter().zip(elevation.iter()) {
                assert_eq!(*s == 1.0, *z <= 0.0);
            }
        }
    }

    #[test]
    fn a_dry_depression_is_land_and_an_enclosed_sea_is_sea_on_any_planet() {
        // At frequency 128 a terrain cell is about 6e-6 of the sphere: a
        // depression of 0.4 of the threshold has a few cells and a sea of 4
        // times it a few dozen. Both lie inside a continent covering the
        // northern hemisphere.
        let cap = |fraction: f64| fms_math::acos(1.0 - 2.0 * fraction);
        let angle = |a: DVec3, b: DVec3| fms_math::atan2(a.cross(b).length(), a.dot(b));
        let dry = from_lat_lon(0.6, 0.0);
        let sea = from_lat_lon(0.6, 1.5);
        let holes = [
            (dry, cap(0.4 * MIN_SEA_FRACTION), -80.0),
            (sea, cap(4.0 * MIN_SEA_FRACTION), -300.0),
        ];
        let mut masks = Vec::new();
        for radius in [6.371e6, 3.39e6] {
            let mut w = world(128, 8, radius);
            let elevation = w
                .terrain
                .mesh
                .cell_centers()
                .iter()
                .map(|&p| {
                    if p.z < 0.0 {
                        return -4000.0;
                    }
                    holes
                        .iter()
                        .find(|&&(c, r, _)| angle(p, c) <= r)
                        .map_or(500.0, |h| h.2)
                })
                .collect();
            w.set_elevation(elevation);
            let t = w.core.terrain;
            let elevation = w.fields.read(t.elevation).to_vec();
            let mask = w.fields.read(t.sea).to_vec();
            let cells = |depth: f64| {
                (0..elevation.len())
                    .filter(|&c| elevation[c] == depth)
                    .collect::<Vec<_>>()
            };
            let (dry_cells, sea_cells) = (cells(-80.0), cells(-300.0));
            assert!(dry_cells.len() >= 2 && sea_cells.len() >= 20);
            assert!(dry_cells.iter().all(|&c| mask[c] == 0.0));
            assert!(sea_cells.iter().all(|&c| mask[c] == 1.0));
            // The depression's climate cell is all land.
            let parent = w.two_level.parents()[dry_cells[0]] as usize;
            assert!(w.fields.read(t.land_fraction)[parent] > 1.0 - 1e-12);
            masks.push(mask);
        }
        assert_eq!(masks[0], masks[1]);
    }
}
