//! Simulation orchestration: builds a world (grids, fields, planet, systems) and
//! steps it.
//!
//! For now a world holds the two grids, the planet, a handful of geometry
//! fields, the terrain (elevation plus its derived statistics), and the
//! insolation systems. Later milestones add the GCM, ocean,
//! hydrology and climate systems here.
//!
//! See `docs/PLAN.md` §3 for the design.

mod insolation;
pub mod probe;
mod terrain;

use std::time::Instant;

use fms_fields::{
    ColorScale, FieldError, FieldId, FieldRegistry, FieldSpec, GridLevel, GridSizes, Location,
    MeshCounts,
};
use fms_geo::coords::to_lat_lon;
use fms_geo::{CellLocator, LloydOptions, MeshGeometry, SphereMesh, TwoLevelMap};
use fms_planet::{Adjustment, Planet, PlanetParams};
use fms_sched::{Calendar, Scheduler, StepReport};

pub use insolation::{DailyInsolation, Insolation};
pub use probe::{Probe, ProbeError};
pub use terrain::{BLANK_OCEAN_ELEVATION, HeightmapImport, TerrainFields};

/// One of the two grids, with its precomputed geometry and point locator.
pub struct Grid {
    /// The mesh.
    pub mesh: SphereMesh,
    /// Metric terms and operator weights.
    pub geometry: MeshGeometry,
    /// Point → cell lookup.
    pub locator: CellLocator,
}

impl Grid {
    /// Build a centroidal grid of the given frequency.
    pub fn new(frequency: u32) -> Self {
        let (mesh, _) = SphereMesh::centroidal(frequency, LloydOptions::default());
        let geometry = MeshGeometry::new(&mesh);
        let locator = CellLocator::new(&mesh);
        Self {
            mesh,
            geometry,
            locator,
        }
    }
}

/// How to build a world.
#[derive(Clone, Copy, Debug)]
pub struct WorldConfig {
    /// Terrain grid frequency (`10·f² + 2` cells). 400 gives 1.6M cells.
    pub terrain_frequency: u32,
    /// Climate grid frequency. 45 gives 20k cells.
    pub climate_frequency: u32,
    /// Base steps per solar day. Until `dt` is derived from the CFL limit
    /// (roadmap M2 step 13), this is set directly.
    pub steps_per_day: u32,
    /// The planet.
    pub planet: PlanetParams,
}

impl Default for WorldConfig {
    fn default() -> Self {
        Self {
            terrain_frequency: 400,
            climate_frequency: 45,
            steps_per_day: 48,
            planet: PlanetParams::earth(),
        }
    }
}

/// Handles to the fields every world has.
#[derive(Clone, Copy, Debug)]
pub struct CoreFields {
    /// Latitude of each terrain cell centre (degrees).
    pub terrain_latitude: FieldId<f64>,
    /// Area of each terrain cell (m²).
    pub terrain_cell_area: FieldId<f64>,
    /// Latitude of each climate cell centre (degrees).
    pub climate_latitude: FieldId<f64>,
    /// Area of each climate cell (m²).
    pub climate_cell_area: FieldId<f64>,
    /// Instantaneous top-of-atmosphere insolation on climate cells (W m⁻²).
    pub toa_insolation: FieldId<f64>,
    /// Daily-mean top-of-atmosphere insolation on climate cells (W m⁻²).
    pub daily_insolation: FieldId<f64>,
    /// Elevation and the fields derived from it.
    pub terrain: TerrainFields,
}

/// Timings from building a world, for diagnostics.
#[derive(Clone, Copy, Debug)]
pub struct BuildTimings {
    /// Building and optimizing the terrain grid (s).
    pub terrain_grid: f64,
    /// Building and optimizing the climate grid (s).
    pub climate_grid: f64,
    /// Building the fine ↔ coarse map (s).
    pub two_level: f64,
}

/// A complete simulated world.
pub struct World {
    /// The planet (validated).
    pub planet: Planet,
    /// Changes made while validating the planet parameters.
    pub adjustments: Vec<Adjustment>,
    /// Fine grid.
    pub terrain: Grid,
    /// Coarse grid.
    pub climate: Grid,
    /// Fine ↔ coarse mapping.
    pub two_level: TwoLevelMap,
    /// All fields.
    pub fields: FieldRegistry,
    /// Handles to the built-in fields.
    pub core: CoreFields,
    /// The scheduler and clock.
    pub scheduler: Scheduler,
    /// How long construction took.
    pub timings: BuildTimings,
}

impl World {
    /// Build a world: grids, mapping, planet, built-in fields and systems. The
    /// terrain starts as a blank ocean ([`BLANK_OCEAN_ELEVATION`]); use
    /// [`World::load_heightmap`] or [`World::set_elevation`] to change it.
    pub fn new(config: WorldConfig) -> Self {
        let (planet, adjustments) = Planet::new(config.planet);

        let t = Instant::now();
        let terrain = Grid::new(config.terrain_frequency);
        let terrain_grid = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let climate = Grid::new(config.climate_frequency);
        let climate_grid = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let two_level =
            TwoLevelMap::new(&climate.mesh, &terrain.mesh, &terrain.geometry.cell_areas);
        let two_level_time = t.elapsed().as_secs_f64();

        let mut fields = FieldRegistry::new(GridSizes {
            terrain: MeshCounts::of(&terrain.mesh),
            climate: MeshCounts::of(&climate.mesh),
        });
        let core = register_core_fields(&mut fields, &planet, &terrain, &climate)
            .expect("built-in field names are valid and unique");

        let calendar = Calendar::new(
            planet.solar_day(),
            config.steps_per_day,
            planet.params().year_days,
        )
        .expect("validated planet gives a valid calendar");
        let mut scheduler = Scheduler::new(calendar);
        let centers = climate.mesh.cell_centers().to_vec();
        scheduler
            .add(Insolation::new(
                planet,
                centers.clone(),
                core.toa_insolation,
            ))
            .expect("built-in system is valid");
        scheduler
            .add(DailyInsolation::new(planet, centers, core.daily_insolation))
            .expect("built-in system is valid");

        let mut world = Self {
            planet,
            adjustments,
            terrain,
            climate,
            two_level,
            fields,
            core,
            scheduler,
            timings: BuildTimings {
                terrain_grid,
                climate_grid,
                two_level: two_level_time,
            },
        };
        world.update_terrain();
        world
    }

    /// The grid a field lives on.
    pub fn grid(&self, level: GridLevel) -> &Grid {
        match level {
            GridLevel::Terrain => &self.terrain,
            GridLevel::Climate => &self.climate,
        }
    }

    /// Simulated seconds since the start (`t = 0` is the northern vernal equinox).
    pub fn time(&self) -> f64 {
        self.scheduler.current_step() as f64 * self.scheduler.calendar().dt()
    }

    /// Advance one base step.
    pub fn step(&mut self) -> StepReport {
        self.scheduler.step(&self.fields)
    }

    /// Advance `n` base steps.
    pub fn run_steps(&mut self, n: u64) {
        self.scheduler.run_steps(&self.fields, n);
    }
}

fn register_core_fields(
    fields: &mut FieldRegistry,
    planet: &Planet,
    terrain: &Grid,
    climate: &Grid,
) -> Result<CoreFields, FieldError> {
    let r2 = planet.params().radius.powi(2);
    let latitudes = |g: &Grid| {
        g.mesh
            .cell_centers()
            .iter()
            .map(|&p| to_lat_lon(p).0.to_degrees())
            .collect()
    };
    let areas = |g: &Grid| g.geometry.cell_areas.iter().map(|a| a * r2).collect();
    let derived = |name, description, units, location| FieldSpec {
        description,
        persistent: false,
        ..FieldSpec::new(name, units, location)
    };
    let latitude_spec = |name, location| FieldSpec {
        valid_range: Some((-90.0, 90.0)),
        color_scale: ColorScale::Diverging,
        ..derived(name, "Latitude of the cell centre", "degrees", location)
    };
    let insolation_spec = |name, description| FieldSpec {
        valid_range: Some((0.0, 2.0 * planet.mean_flux())),
        color_scale: ColorScale::Temperature,
        ..derived(name, description, "W m-2", Location::CLIMATE_CELLS)
    };

    Ok(CoreFields {
        terrain_latitude: fields.register_with(
            latitude_spec("geo.terrain_latitude", Location::TERRAIN_CELLS),
            latitudes(terrain),
        )?,
        terrain_cell_area: fields.register_with(
            derived(
                "geo.terrain_cell_area",
                "Cell area",
                "m2",
                Location::TERRAIN_CELLS,
            ),
            areas(terrain),
        )?,
        climate_latitude: fields.register_with(
            latitude_spec("geo.climate_latitude", Location::CLIMATE_CELLS),
            latitudes(climate),
        )?,
        climate_cell_area: fields.register_with(
            derived(
                "geo.climate_cell_area",
                "Cell area",
                "m2",
                Location::CLIMATE_CELLS,
            ),
            areas(climate),
        )?,
        toa_insolation: fields.register(
            insolation_spec(
                "planet.toa_insolation",
                "Top-of-atmosphere sunlight, instantaneous",
            ),
            0.0,
        )?,
        daily_insolation: fields.register(
            insolation_spec(
                "planet.daily_insolation",
                "Top-of-atmosphere sunlight, daily mean",
            ),
            0.0,
        )?,
        terrain: terrain::register(fields)?,
    })
}
