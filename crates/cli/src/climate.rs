//! `fms-climate run`: the layered model (four layers by default, the bottom one
//! the boundary layer) driven by radiation and moist physics over slab
//! ocean and land tiles and the world's terrain, on the world's planet
//! (orbit, CO₂, surface pressure). Spins up, then averages the final year: zonal means, a global
//! summary with the planet's energy and water budgets, and annual-mean maps,
//! in `output/climate/`. On request, the sea temperatures are prescribed
//! (`--fixed-sea-temperature`) and the run is checked for a pattern locked
//! to the grid (`--imprint-block`, [`crate::imprint`]).

use std::fmt::Write as _;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::{Args, ValueEnum};
use fms_atmos::column_physics::{
    ApeSst, CLIMATE_LEVELS, LATENT_HEAT, LandParams, OceanSpinUp, initial_state,
    initial_surface_temperature_by_latitude,
};
use fms_atmos::moisture::{CloudParams, MoistureParams};
use fms_atmos::{
    ColumnPhysics, Forcing, GridImprint, LayeredSemiImplicit, PhysicsParams, ZonalClimatology,
};
use fms_fields::{ColorScale, FieldId, FieldSpec, Location};
use fms_grid::LatitudeBands;
use fms_ocean::{
    EddyParams, EkmanParams, GyreParams, IceParams, JumpLayers, OceanParams, WET_OCEAN_FRACTION,
};
use glam::DVec3;

use crate::basins::write_gyre_map;
use crate::chart::{self, Series, TimeAxis};
use crate::layered_run::{
    DynamicsArgs, Fields, Output, Plan, SAMPLES_PER_DAY, build_model, describe, integrate,
    sample_every, stepper, write_map,
};
use crate::ocean_transport::{edge_latitudes, northward_transport};
use crate::render::write_png;
use crate::restart;
use crate::{ImageArgs, WorldArgs, build_world};

/// Design wind for the step: radiatively driven jets and storms reach
/// ~100 m s⁻¹ (stronger than Held–Suarez's).
const MAX_WIND: f64 = 120.0;

#[derive(Args)]
pub struct ClimateArgs {
    #[command(flatten)]
    world: WorldArgs,
    // The spin-up options' defaults depend on whether the ocean has
    // currents (`preset`), so clap leaves them unset.
    /// Simulated years; statistics come from the last one [default: 9 with
    /// currents, 6 without].
    #[arg(long)]
    years: Option<f64>,
    /// Years that end with the ocean jumping forward, from
    /// `--jump-first-year` on, one every `--jump-every`; 0 for a run in true
    /// time. Needs `--years` long enough for two years in true time after
    /// the last jump [default: 2 with currents, 3 without].
    #[arg(long)]
    ocean_jumps: Option<u32>,
    /// Years of the ocean's change added by each jump [default: 20 with
    /// currents, 3 without].
    #[arg(long)]
    jump_years: Option<f64>,
    /// Most any ocean temperature may change in one jump (K) [default: 3].
    #[arg(long)]
    jump_limit: Option<f64>,
    /// Which ocean layers jump: `both`, or `deep` alone, the slow store, which
    /// the mixed layer then follows [default: deep with currents, both
    /// without].
    #[arg(long, value_enum)]
    jump_layers: Option<JumpLayersArg>,
    /// The first year that ends with an ocean jump, counting the run's first
    /// year as 1 [default: 4 with currents, 2 without].
    #[arg(long, value_parser = clap::value_parser!(u32).range(1..))]
    jump_first_year: Option<u32>,
    /// Years from one ocean jump to the next. The years between run in true
    /// time, and each jump repeats only the change of the year it ends
    /// [default: 3 with currents, 1 without].
    #[arg(long, value_parser = clap::value_parser!(u32).range(1..))]
    jump_every: Option<u32>,
    /// Horizontal heat diffusivity of the ocean (m²/s); 0 for none. The
    /// paper's code uses 6250 when eddy transport is off.
    #[arg(long, default_value_t = OceanParams::default().diffusivity)]
    ocean_diffusivity: f64,
    /// Let a mixed layer colder than the deep layer stay on top.
    #[arg(long)]
    no_ocean_convection: bool,
    /// Turn off Ekman transport, the wind-driven flow in the mixed layer that
    /// the deep layer returns, with upwelling and downwelling between.
    #[arg(long)]
    no_ocean_ekman: bool,
    /// Keep the paper's Eq. 4 at the equator instead of blending in the
    /// Sverdrup transport of the wind's curl; ignored with
    /// `--no-ocean-ekman`.
    #[arg(long)]
    no_ocean_sverdrup: bool,
    /// Turn off Gent–McWilliams eddy transport, the flow in the mixed layer
    /// down the slope of the isotherms that the deep layer returns.
    #[arg(long)]
    no_ocean_eddies: bool,
    /// Turn off the wind-driven gyres, Stommel's depth-integrated flow with
    /// its western boundary currents, which both layers carry their share of.
    #[arg(long)]
    no_ocean_gyres: bool,
    /// Depth the gyres' flow is spread over (m); the mixed and deep layers
    /// carry their depths' share of it, the rest flows below them. Ignored
    /// with `--no-ocean-gyres`.
    #[arg(long, default_value_t = GyreParams::default().depth)]
    ocean_gyre_depth: f64,
    /// The gyres' linear bottom friction (1/s), which sets how wide the
    /// western boundary currents are and how much flows round Antarctica.
    /// Ignored with `--no-ocean-gyres`.
    #[arg(long, default_value_t = GyreParams::default().friction)]
    ocean_gyre_friction: f64,
    /// Hours between the ocean's transports (currents, eddies, diffusion,
    /// convective adjustment), rounded to whole steps; 0 for every step.
    /// The sea's surface heat and ice are always stepped every step.
    /// Default: the radiation's interval.
    #[arg(long)]
    ocean_interval: Option<f64>,
    /// Without currents (`--no-ocean-ekman`, `--no-ocean-eddies` and
    /// `--no-ocean-gyres`), the deep layer starts at the mixed layer's
    /// temperature but no warmer than this (K); `inf` for no limit. Ignored
    /// with currents, where the deep layer always starts at the mixed layer's
    /// temperature.
    #[arg(long, default_value_t = OceanParams::default().deep_start_limit)]
    ocean_deep_start_limit: f64,
    /// Give sea ice and its snow the albedo of open water, which turns off
    /// the ice-albedo feedback.
    #[arg(long)]
    no_ice_albedo: bool,
    /// Keep the seasonal snow surface on land that stays snow-covered.
    #[arg(long)]
    no_ice_sheet: bool,
    /// Land starting colder than this (K) starts as an ice sheet.
    #[arg(long, default_value_t = LandParams::default().ice_sheet_start_temperature)]
    ice_sheet_start_temperature: f64,
    /// Smoothing passes over the terrain the atmosphere sees; 0 for the
    /// plain cell means.
    #[arg(long, default_value_t = fms_sim::SURFACE_SMOOTHING_PASSES)]
    terrain_smoothing: usize,
    /// Specific humidity (kg/kg) below which low cloud thins; 0 for none.
    #[arg(long, default_value_t = CloudParams::default().freeze_dry_humidity)]
    freeze_dry_humidity: f64,
    /// Latitude bands for the zonal means (default: about one per grid
    /// spacing).
    #[arg(long)]
    bands: Option<usize>,
    /// The GIF covers this many final days.
    #[arg(long, default_value_t = 30.0)]
    movie_days: f64,
    /// GIF frames.
    #[arg(long, default_value_t = 121)]
    frames: u32,
    /// GIF frames per second.
    #[arg(long, default_value_t = 12)]
    fps: u16,
    /// Write the run's final state to FILE when it has run to the end, to
    /// continue from with `--restart`. A run that saves its state may end
    /// before two years in true time follow the last ocean jump.
    #[arg(long, value_name = "FILE")]
    save_restart: Option<PathBuf>,
    /// Start from the state saved in FILE by `--save-restart`. Every setting
    /// still comes from this command line, never from the file; the run says
    /// whether it continues the saved one exactly or is a warm start, and
    /// how the state was adapted to changed settings. `--years` counts from
    /// the restart and may be 1; steps, days, years and months go on from
    /// the file's.
    #[arg(long, value_name = "FILE")]
    restart: Option<PathBuf>,
    /// Hold the sea temperatures at the Aqua-Planet Experiment's CONTROL
    /// profile, a function of latitude alone (27 °C at the equator, 0 °C
    /// poleward of 60°), instead of running the ocean: no currents, sea ice
    /// or ocean jumps (the jump options are ignored), and `--years` may be
    /// any number above 1. The world must be all ocean.
    #[arg(long, conflicts_with_all = ["restart", "save_restart"])]
    fixed_sea_temperature: bool,
    /// Check the run for a pattern locked to the grid, in blocks of DAYS:
    /// from `--imprint-spinup` days into this run to its end, the samples
    /// go to two halves in alternating blocks, and the time mean of the
    /// dynamics' fields, precipitation, evaporation and surface temperature
    /// is split into what the halves share and what they don't. Writes
    /// `imprint.csv`, maps and a block of the summary. The check isn't
    /// saved in restart files: after `--restart` it covers this run only.
    #[arg(long, value_name = "DAYS")]
    imprint_block: Option<f64>,
    /// Days from the start of this run before the grid-imprint check
    /// samples; ignored without `--imprint-block`.
    #[arg(long, value_name = "DAYS", default_value_t = 200.0)]
    imprint_spinup: f64,
    #[command(flatten)]
    dynamics: DynamicsArgs,
    #[command(flatten)]
    image: ImageArgs,
}

/// `--jump-layers`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum JumpLayersArg {
    /// The mixed layer and the deep layer.
    Both,
    /// The deep layer alone.
    Deep,
}

impl From<JumpLayersArg> for JumpLayers {
    fn from(layers: JumpLayersArg) -> Self {
        match layers {
            JumpLayersArg::Both => JumpLayers::Both,
            JumpLayersArg::Deep => JumpLayers::Deep,
        }
    }
}

/// The default spin-up and run length (years) for an ocean with or without
/// currents.
fn preset(currents: bool) -> (OceanSpinUp, f64) {
    if currents {
        (OceanSpinUp::with_currents(), 9.0)
    } else {
        (OceanSpinUp::still(), 6.0)
    }
}

impl ClimateArgs {
    /// The spin-up and run length (years): the preset for an ocean with or
    /// without `currents`, with each option given replacing its value.
    fn spin_up(&self, currents: bool) -> (OceanSpinUp, f64) {
        let (preset, years) = preset(currents);
        let spin_up = OceanSpinUp {
            jumps: self.ocean_jumps.unwrap_or(preset.jumps),
            extra_years: self.jump_years.unwrap_or(preset.extra_years),
            limit: self.jump_limit.unwrap_or(preset.limit),
            layers: self.jump_layers.map_or(preset.layers, JumpLayers::from),
            first_year: self.jump_first_year.unwrap_or(preset.first_year),
            every: self.jump_every.unwrap_or(preset.every),
        };
        (spin_up, self.years.unwrap_or(years))
    }

    /// The ocean the options ask for.
    fn ocean(&self) -> OceanParams {
        let radiation = PhysicsParams::default().radiation;
        OceanParams {
            diffusivity: self.ocean_diffusivity,
            convective_adjustment: !self.no_ocean_convection,
            ekman: (!self.no_ocean_ekman).then(|| EkmanParams {
                sverdrup: !self.no_ocean_sverdrup,
                ..EkmanParams::default()
            }),
            eddies: (!self.no_ocean_eddies).then(EddyParams::default),
            gyres: (!self.no_ocean_gyres).then(|| GyreParams {
                depth: self.ocean_gyre_depth,
                friction: self.ocean_gyre_friction,
                ..GyreParams::default()
            }),
            deep_start_limit: self.ocean_deep_start_limit,
            ice: OceanParams::default().ice.map(|ice| {
                if self.no_ice_albedo {
                    IceParams {
                        albedo_limit: radiation.albedo,
                        snow_albedo: radiation.albedo[0],
                        ..ice
                    }
                } else {
                    ice
                }
            }),
            ..OceanParams::default()
        }
    }
}

/// The spin-up, for the run's description.
fn describe_spin_up(spin_up: &OceanSpinUp) -> String {
    let Some(last) = spin_up.last_jump_year() else {
        return "ocean spin-up: none (true time)".into();
    };
    let years: Vec<String> = (spin_up.first_year..=last)
        .filter(|&y| spin_up.ends_in_jump(y))
        .map(|y| y.to_string())
        .collect();
    let years = match years.split_last() {
        Some((last, [])) => format!("the end of year {last}"),
        Some((last, rest)) => format!("the ends of years {} and {last}", rest.join(", ")),
        None => unreachable!("a spin-up with a last jump has jumps"),
    };
    let layers = match spin_up.layers {
        JumpLayers::Both => "two-layer",
        JumpLayers::Deep => "deep-layer",
    };
    format!(
        "ocean spin-up: {} {layers} jump{} of {} years at {years}, at most {} K",
        spin_up.jumps,
        if spin_up.jumps == 1 { "" } else { "s" },
        spin_up.extra_years,
        spin_up.limit
    )
}

/// Checks that a run of `years` leaves at least two years in true time
/// after the last of the `spin_up`'s jumps.
fn check_years(years: f64, spin_up: &OceanSpinUp) -> Result<(), String> {
    let Some(last) = spin_up.last_jump_year() else {
        return Ok(());
    };
    let least = u64::from(last) + 2;
    if years < least as f64 {
        return Err(format!(
            "--years must be at least {least} with --ocean-jumps {} --jump-first-year {} --jump-every {}: \
             the last jump ends year {last}, and two years in true time must follow it",
            spin_up.jumps, spin_up.first_year, spin_up.every
        ));
    }
    Ok(())
}

/// [`check_years`] for a run whose last year is `years`, counted from the
/// spin-up's start; one that saves a restart (`saving`) may end sooner, to
/// be continued, with a warning.
fn check_end_year(years: f64, spin_up: &OceanSpinUp, saving: bool) -> Result<(), String> {
    match check_years(years, spin_up) {
        Err(e) if saving => {
            eprintln!(
                "warning: {e}; running anyway, as --save-restart saves the state to continue from"
            );
            Ok(())
        }
        checked => checked,
    }
}

/// Writes the final year's mean of each field for every cell, with the
/// height of the terrain the atmosphere sees (m), to compare runs over a part
/// of the map such as an ice sheet.
fn write_annual_cells(
    path: &Path,
    centers: &[DVec3],
    land: &[f64],
    height: &[f64],
    names: &[&str],
    means: &[Vec<f64>],
) -> Result<(), String> {
    let mut rows = format!(
        "latitude_deg,longitude_deg,land_fraction,surface_height_m,{}\n",
        names.join(",")
    );
    for (i, p) in centers.iter().enumerate() {
        let _ = write!(
            rows,
            "{:.3},{:.3},{:.6},{:.1}",
            fms_math::asin(p.z / p.length()).to_degrees(),
            fms_math::atan2(p.y, p.x).to_degrees(),
            land[i],
            height[i],
        );
        for mean in means {
            let _ = write!(rows, ",{:.6}", mean[i]);
        }
        rows.push('\n');
    }
    std::fs::write(path, rows).map_err(|e| format!("cannot write {}: {e}", path.display()))
}

/// The extra fields of the grid-imprint check, as `(name, units, mirror)`:
/// the physics' cell fields of those names.
const IMPRINT_FIELDS: [(&str, &str, f64); 3] = [
    ("atm.precipitation", "mm/day", 1.0),
    ("atm.evaporation", "mm/day", 1.0),
    ("rad.surface_temperature", "K", 1.0),
];

/// The maps of the ocean's state and currents, which a prescribed sea
/// leaves at their starting values: not drawn.
const OCEAN_MAPS: [&str; 3] = [
    "ocean.deep_temperature",
    "ocean.gyre_streamfunction",
    "ocean.gyre_heating",
];

pub fn run(args: ClimateArgs) -> Result<(), String> {
    let run_started = Instant::now();
    let ocean = args.ocean();
    let fixed = args.fixed_sea_temperature;
    let (mut spin_up, years) = args.spin_up(ocean.has_currents());
    if fixed {
        // A prescribed sea doesn't jump.
        spin_up.jumps = 0;
    }
    if let Some(block) = args.imprint_block
        && (block.is_nan() || block <= 0.0)
    {
        return Err("--imprint-block must be positive".into());
    }
    if args.imprint_spinup.is_nan() || args.imprint_spinup < 0.0 {
        return Err("--imprint-spinup must be 0 or more".into());
    }
    let saving = args.save_restart.is_some();
    if args.restart.is_some() {
        // Checked against the spin-up once the restart says where it is.
        if years.is_nan() || years < 1.0 {
            return Err(
                "--years must be at least 1 with --restart (the last year is averaged)".into(),
            );
        }
    } else {
        if years <= 1.0 {
            return Err("--years must be more than 1 (the last year is averaged)".into());
        }
        check_end_year(years, &spin_up, saving)?;
    }
    let mut world = build_world(&args.world)?;
    let fields = Fields::register(&mut world)?;
    let planet_params = *world.planet.params();
    let params = PhysicsParams {
        ocean,
        land: LandParams {
            ice_sheet_years: if args.no_ice_sheet {
                f64::INFINITY
            } else {
                LandParams::default().ice_sheet_years
            },
            ice_sheet_start_temperature: args.ice_sheet_start_temperature,
            ..LandParams::default()
        },
        moisture: MoistureParams {
            clouds: CloudParams {
                freeze_dry_humidity: args.freeze_dry_humidity,
                ..CloudParams::default()
            },
            ..MoistureParams::default()
        },
        co2: planet_params.co2_ppm,
        spin_up,
        ..PhysicsParams::default()
    };
    // Annual-mean maps of these of the physics' cell fields
    // (`ColumnPhysics::cell_fields`).
    let maps = [
        (
            "rad.surface_temperature",
            "Surface temperature",
            "K",
            ColorScale::Temperature,
        ),
        (
            "rad.outgoing_longwave",
            "Outgoing longwave radiation",
            "W m-2",
            ColorScale::Sequential,
        ),
        (
            "rad.absorbed_shortwave",
            "Absorbed sunlight",
            "W m-2",
            ColorScale::Sequential,
        ),
        (
            "atm.water_vapour",
            "Water-vapour column",
            "kg m-2",
            ColorScale::Moisture,
        ),
        (
            "atm.precipitation",
            "Precipitation",
            "mm day-1",
            ColorScale::Moisture,
        ),
        (
            "atm.orographic_precipitation",
            "Orographic precipitation",
            "mm day-1",
            ColorScale::Moisture,
        ),
        (
            "atm.evaporation",
            "Evaporation",
            "mm day-1",
            ColorScale::Moisture,
        ),
        (
            "atm.cloud_cover",
            "Total cloud cover",
            "1",
            ColorScale::Sequential,
        ),
        (
            "atm.relative_humidity_lower",
            "Lower-layer relative humidity",
            "1",
            ColorScale::Moisture,
        ),
        (
            "land.soil_water",
            "Soil water (per unit of cell area)",
            "kg m-2",
            ColorScale::Moisture,
        ),
        (
            "land.snow",
            "Snow on land (per unit of cell area)",
            "kg m-2",
            ColorScale::Moisture,
        ),
        (
            "land.ice_sheet",
            "Ice sheet (per unit of cell area)",
            "1",
            ColorScale::Sequential,
        ),
        (
            "ocean.ice_fraction",
            "Sea-ice cover of the ocean",
            "1",
            ColorScale::Sequential,
        ),
        (
            "ocean.snow_on_ice",
            "Snow on sea ice, per unit area of ice",
            "kg m-2",
            ColorScale::Sequential,
        ),
        (
            "ocean.ice_thickness",
            "Mean sea-ice thickness",
            "m",
            ColorScale::Sequential,
        ),
        (
            "ocean.deep_temperature",
            "Deep-ocean temperature",
            "K",
            ColorScale::Temperature,
        ),
        (
            "ocean.gyre_streamfunction",
            "Streamfunction of the wind-driven gyres (clockwise round a maximum)",
            "Sv",
            ColorScale::Diverging,
        ),
        (
            "ocean.gyre_heating",
            "Heating of the ocean by the wind-driven gyres, per unit area of ocean",
            "W m-2",
            ColorScale::Diverging,
        ),
        (
            "rad.shortwave_cloud_effect",
            "Shortwave cloud radiative effect",
            "W m-2",
            ColorScale::Sequential,
        ),
        (
            "rad.longwave_cloud_effect",
            "Longwave cloud radiative effect",
            "W m-2",
            ColorScale::Sequential,
        ),
    ];
    // Further global means for the summary.
    let extra = [
        "atm.convective_precipitation",
        "atm.snowfall",
        "atm.latent_heat_flux",
        "atm.sensible_heat_flux",
        "atm.high_cloud",
        "atm.low_cloud",
        "atm.relative_humidity_upper",
        "rad.surface_shortwave",
        "rad.surface_longwave",
        "rad.clear_surface_longwave",
        "atm.surface_air_temperature",
        "atm.surface_stress_east",
        "ocean.mixed_temperature",
    ];
    let names: Vec<&str> = maps.iter().map(|m| m.0).chain(extra).collect();
    // Monthly zonal means over the whole run, for the seasons and the drift.
    let monthly = [
        "rad.surface_temperature",
        "ocean.mixed_temperature",
        "ocean.deep_temperature",
        "rad.absorbed_shortwave",
        "rad.outgoing_longwave",
        "atm.precipitation",
        "atm.cloud_cover",
        "ocean.snow_on_ice",
        "atm.surface_air_temperature",
    ];
    let mut ids: Vec<FieldId<f64>> = Vec::new();
    for (name, description, units, color_scale) in maps {
        let spec = FieldSpec {
            description,
            color_scale,
            persistent: false,
            ..FieldSpec::new(name, units, Location::CLIMATE_CELLS)
        };
        ids.push(
            world
                .fields
                .register(spec, 0.0)
                .map_err(|e| e.to_string())?,
        );
    }
    // The gyres' ocean basins, which are fixed: drawn once, at the start.
    let basin_id = world
        .fields
        .register(
            FieldSpec {
                description: "Ocean basin of the wind-driven gyres, numbered from 1 by area \
                              (0 where there is none)",
                color_scale: ColorScale::Categorical,
                persistent: false,
                ..FieldSpec::new("ocean.basin", "", Location::CLIMATE_CELLS)
            },
            0.0,
        )
        .map_err(|e| e.to_string())?;
    let world = world;

    // Start from the single column's equilibrium, at rest, with its
    // humidity.
    let geopotential = world.surface_geopotential_smoothed(args.terrain_smoothing);
    let height: Vec<f64> = geopotential
        .iter()
        .map(|phi| phi / planet_params.gravity)
        .collect();
    let setup = build_model(
        &world,
        &args.dynamics,
        MAX_WIND,
        &CLIMATE_LEVELS,
        geopotential,
    )?;
    // Radiation every whole number of steps, within 3 h and an eighth of the day.
    let mut params = params.with_step(setup.dt, setup.day);
    if let Some(hours) = args.ocean_interval {
        if hours.is_nan() || hours < 0.0 {
            return Err("--ocean-interval must be 0 or more".into());
        }
        params.transport_interval = (hours * 3600.0 / setup.dt).round().max(1.0) * setup.dt;
    }
    let model = setup
        .model
        .with_humidity_profile(params.radiation)
        .with_dissipative_heating();
    let (mesh, geometry) = (&world.climate.mesh, &world.climate.geometry);
    let nc = mesh.n_cells();
    let land = world.fields.read(world.core.terrain.land_fraction).to_vec();
    if fixed && land.iter().any(|f| *f > 0.0) {
        return Err(
            "--fixed-sea-temperature: only the ocean's temperature can be prescribed, and this \
             world has land"
                .into(),
        );
    }
    let relief = world
        .fields
        .read(world.core.terrain.surface_height_stddev)
        .to_vec();
    let surface_height = world
        .fields
        .read(world.core.terrain.surface_height)
        .to_vec();
    let settings = restart::Settings {
        planet: &planet_params,
        physics: &params,
        sigma: model.levels().interfaces(),
        dt: setup.dt,
        diffusion: model.hyperdiffusion().map(|d| d.options()),
        terrain: [
            model.surface_geopotential(),
            &relief,
            &surface_height,
            &land,
        ],
    }
    .fingerprint();
    let (state, surface_temperature) = initial_state(&model, &world.planet, &params);
    eprintln!("climate: single-column start at {surface_temperature:.1} K");
    let start_surface = if fixed {
        mesh.cell_centers()
            .iter()
            .map(|&p| ApeSst::Control.at(fms_grid::coords::to_lat_lon(p).0))
            .collect()
    } else {
        initial_surface_temperature_by_latitude(&model, &world.planet, &params, surface_temperature)
    };
    let physics = ColumnPhysics::with_land(
        &model,
        &world.planet,
        params,
        0.0,
        &state,
        start_surface,
        land.clone(),
    )
    .with_relief(relief, &surface_height);
    let physics = if fixed {
        physics.with_prescribed_surface()
    } else {
        physics
    };
    // Where each wanted field sits in the physics' list.
    let fields_now = physics.cell_fields();
    let imprint_index: Vec<usize> = IMPRINT_FIELDS
        .iter()
        .map(|(n, _, _)| {
            fields_now
                .iter()
                .position(|f| f.name == *n)
                .ok_or_else(|| format!("the physics has no field {n}"))
        })
        .collect::<Result<_, _>>()?;
    let index: Vec<usize> = names
        .iter()
        .map(|n| {
            fields_now
                .iter()
                .position(|f| f.name == *n)
                .ok_or_else(|| format!("the physics has no field {n}"))
        })
        .collect::<Result<_, _>>()?;
    let monthly_index: Vec<usize> = monthly
        .iter()
        .map(|n| {
            fields_now
                .iter()
                .position(|f| f.name == *n)
                .ok_or_else(|| format!("the physics has no field {n}"))
        })
        .collect::<Result<_, _>>()?;
    let basins = fields_now
        .iter()
        .find(|f| f.name == "ocean.basin")
        .expect("the physics reports the ocean's basins")
        .values
        .clone();
    let ice_index = ["ocean.ice_fraction", "ocean.ice_thickness"].map(|n| {
        fields_now
            .iter()
            .position(|f| f.name == n)
            .expect("the physics reports its ice")
    });
    // The heat each ocean transport has brought into each column since the
    // start, whose change over the final year is its mean heating.
    let gained_index = [
        "ocean.ekman_heat_gained",
        "ocean.eddy_heat_gained",
        "ocean.gyre_heat_gained",
        "ocean.diffusion_heat_gained",
    ]
    .map(|n| {
        fields_now
            .iter()
            .position(|f| f.name == n)
            .expect("the physics reports the ocean's transport")
    });
    let model = model.with_forcing(Box::new(physics));
    let bands = args
        .bands
        .unwrap_or_else(|| 2 * ((nc as f64 / 10.0).sqrt().round() as usize).max(8));
    let mut climate = ZonalClimatology::new(&model, bands);
    let lat_bands = LatitudeBands::new(mesh, geometry, bands);
    // The grid-imprint check, in blocks of whole samples.
    let mut imprint = args.imprint_block.map(|days| {
        let samples_per_day = setup.steps_per_day as f64 / sample_every(setup.steps_per_day) as f64;
        let block = ((days * samples_per_day).round() as u64).max(1);
        let (upper, lower) = crate::layered_run::upper_lower(&model);
        GridImprint::with_fields(&model, upper, lower, block, bands, &IMPRINT_FIELDS)
    });
    let surface = if fixed {
        format!(
            "; {} steps to the {:.4} h solar day, radiation every {} ({:.2} h); radiation and moisture over sea \
             temperatures prescribed as the Aqua-Planet Experiment's CONTROL profile (no ocean currents, sea ice or \
             ocean spin-up), CO₂ {} ppm",
            setup.steps_per_day,
            setup.day / 3600.0,
            fms_atmos::radiation_steps(setup.dt, setup.day),
            params.radiation_interval / 3600.0,
            params.co2,
        )
    } else {
        format!(
            "; {} steps to the {:.4} h solar day, radiation every {} ({:.2} h), ocean transport every {} ({:.2} h); radiation and moisture over a {} m mixed layer on a {} m deep ocean {}, and a {} m soil store, CO₂ {} ppm; {}",
            setup.steps_per_day,
            setup.day / 3600.0,
            fms_atmos::radiation_steps(setup.dt, setup.day),
            params.radiation_interval / 3600.0,
            (params.transport_interval / setup.dt).round(),
            params.transport_interval / 3600.0,
            params.ocean.mixed_depth,
            params.ocean.deep_depth,
            currents(&params.ocean),
            params.land.soil_depth,
            params.co2,
            describe_spin_up(&params.spin_up)
        )
    };
    let description = describe(
        &model,
        &state,
        setup.dt,
        (years * world.planet.year_seconds() / setup.day).round(),
        &surface,
    );
    eprintln!("climate: {description}");

    let dt = setup.dt;
    let per_year = (world.planet.year_seconds() / dt).round() as u64;
    let mut si = stepper(model, state, dt);
    // Steps, days, years and months count on from a restart's.
    let mut description = description;
    let start_step = match &args.restart {
        Some(path) => {
            let notes = restart::load(path, &mut si, settings)?;
            let report = restart::report(path, &si, &notes, per_year, setup.day);
            eprintln!("climate: {report}");
            description.push('\n');
            description.push_str(&report);
            let years_run = si.steps() as f64 / per_year as f64;
            check_end_year(years_run + years, &spin_up, saving)?;
            si.steps()
        }
        None => 0,
    };
    let total = start_step + (years * per_year as f64).round() as u64;
    let stats_from = total - per_year;
    let imprint_from =
        start_step + (args.imprint_spinup * setup.steps_per_day as f64).round() as u64;
    let movie_from = total
        .saturating_sub((args.movie_days * setup.steps_per_day as f64).round() as u64)
        .max(start_step);
    let out = Output {
        command: "run",
        dir: args.image.out_dir.join("climate"),
        image: &args.image,
        fps: args.fps,
    };
    let plan = Plan {
        total,
        frames: u64::from(args.frames),
        frames_from: movie_from,
    };
    let mut sums: Vec<Vec<f64>> = vec![vec![0.0; nc]; names.len()];
    let mut samples = 0.0;
    // The budgets and the ocean's transported heat at the start of the
    // statistics year.
    let snapshot = |si: &LayeredSemiImplicit| {
        let forcing = si.model().forcing().expect("climate has forcing");
        let fields = forcing.cell_fields();
        let start = (
            si.model().diagnostics(si.state()),
            si.energy_budget(),
            forcing.budget().expect("physics keeps a budget"),
            si.water(),
            si.water_budget(),
        );
        let gained: Vec<Vec<f64>> = gained_index
            .iter()
            .map(|&k| fields[k].values.clone())
            .collect();
        (start, (si.steps(), gained))
    };
    // A restart may begin the statistics year.
    let (mut start, mut gained_start) = if stats_from <= si.steps() {
        let (start, gained) = snapshot(&si);
        (Some(start), Some(gained))
    } else {
        (None, None)
    };
    let mut saved = Ok(());
    std::fs::create_dir_all(&out.dir)
        .map_err(|e| format!("cannot create {}: {e}", out.dir.display()))?;
    world.fields.write(basin_id).copy_from_slice(&basins);
    write_map(&world, basin_id.key(), None, &out.dir, "FIXED", &args.image)?;
    let path = out.dir.join("monthly_zonal.csv");
    let file = std::fs::File::create(&path)
        .map_err(|e| format!("cannot create {}: {e}", path.display()))?;
    let mut monthly_file = std::io::BufWriter::new(file);
    let mut monthly_result = writeln!(
        monthly_file,
        "year,month,latitude_deg,area_share,ice_cover,ice_volume_m,{}",
        monthly.join(",")
    );
    // Each month's rain belt, from the same zonal means.
    let path = out.dir.join("rain_belt.csv");
    let file = std::fs::File::create(&path)
        .map_err(|e| format!("cannot create {}: {e}", path.display()))?;
    let mut rain_file = std::io::BufWriter::new(file);
    let mut rain_result = writeln!(
        rain_file,
        "year,month,centroid_latitude_deg,tropical_rain_mm_day,subtropical_rain_mm_day"
    );
    let rain_index = 2 + monthly
        .iter()
        .position(|n| *n == "atm.precipitation")
        .expect("precipitation is a monthly field");
    let temperature_index = 2 + monthly
        .iter()
        .position(|n| *n == "rad.surface_temperature")
        .expect("surface temperature is a monthly field");
    let band_latitudes = lat_bands.latitudes();
    let band_degrees: Vec<f64> = band_latitudes.iter().map(|l| l.to_degrees()).collect();
    // The last 12 completed months' zonal-mean precipitation, by month number.
    let mut last_months: Vec<(u64, Vec<f64>)> = Vec::new();
    // And their zonal-mean surface temperature.
    let mut last_temperature: Vec<Vec<f64>> = Vec::new();
    let band_shares: Vec<f64> = {
        let total: f64 = lat_bands.areas().iter().sum();
        lat_bands.areas().iter().map(|a| a / total).collect()
    };
    let mut month_sums: Vec<Vec<f64>> = vec![vec![0.0; nc]; monthly.len() + 2];
    let mut month_samples = 0.0;
    let steps_in_month = |month: u64| (month as f64 * per_year as f64 / 12.0).round() as u64;
    // Months already complete when a restart begins.
    let mut month_done = 0_u64;
    while steps_in_month(month_done + 1) <= start_step {
        month_done += 1;
    }
    let first_month = month_done;
    let log = integrate(&world, &fields, &mut si, &plan, &out, &description, |si| {
        if si.steps() == total
            && let Some(path) = &args.save_restart
        {
            saved = restart::save(path, si, settings);
        }
        let forcing = si.model().forcing().expect("climate has forcing");
        {
            let fields = forcing.cell_fields();
            let (cover, thickness) = (&fields[ice_index[0]].values, &fields[ice_index[1]].values);
            for (i, land) in land.iter().enumerate() {
                let c = cover[i] * (1.0 - land);
                month_sums[0][i] += c;
                month_sums[1][i] += c * thickness[i];
            }
            for (sum, &f) in month_sums[2..].iter_mut().zip(&monthly_index) {
                for (s, v) in sum.iter_mut().zip(&fields[f].values) {
                    *s += v;
                }
            }
            month_samples += 1.0;
            if let Some(imprint) = &mut imprint
                && si.steps() > imprint_from
            {
                let extra: Vec<&[f64]> = imprint_index
                    .iter()
                    .map(|&f| fields[f].values.as_slice())
                    .collect();
                imprint.sample_with(si.model(), si.state(), &extra);
            }
        }
        if si.steps() >= steps_in_month(month_done + 1) {
            let zonal: Vec<Vec<f64>> = month_sums
                .iter()
                .map(|m| {
                    let mean: Vec<f64> = m.iter().map(|v| v / month_samples).collect();
                    lat_bands.mean(geometry, &mean)
                })
                .collect();
            let (year, month) = (month_done / 12 + 1, month_done % 12 + 1);
            let mut rows = String::new();
            for b in 0..band_latitudes.len() {
                let _ = write!(
                    rows,
                    "{year},{month},{:.3},{:.6}",
                    band_latitudes[b].to_degrees(),
                    band_shares[b],
                );
                for z in &zonal {
                    let _ = write!(rows, ",{:.6}", z[b]);
                }
                rows.push('\n');
            }
            if monthly_result.is_ok() {
                monthly_result = monthly_file
                    .write_all(rows.as_bytes())
                    .and_then(|()| monthly_file.flush());
            }
            let rain = &zonal[rain_index];
            let value = |v: Option<f64>| v.map_or(String::new(), |v| format!("{v:.6}"));
            if rain_result.is_ok() {
                rain_result = writeln!(
                    rain_file,
                    "{year},{month},{},{},{}",
                    value(fms_atmos::rain_belt(&band_degrees, rain)),
                    value(fms_atmos::belt_mean(&band_degrees, rain, 0.0, 10.0)),
                    value(fms_atmos::belt_mean(&band_degrees, rain, 20.0, 30.0)),
                )
                .and_then(|()| rain_file.flush());
            }
            if last_months.len() == crate::circulation::RAIN_BELT_MONTHS {
                last_months.remove(0);
                last_temperature.remove(0);
            }
            last_months.push((month, rain.clone()));
            last_temperature.push(zonal[temperature_index].clone());
            for m in &mut month_sums {
                m.fill(0.0);
            }
            month_samples = 0.0;
            month_done += 1;
        }
        if start.is_none() && si.steps() >= stats_from {
            let (now, gained) = snapshot(si);
            start = Some(now);
            gained_start = Some(gained);
        }
        if si.steps() > stats_from {
            climate.sample(si.model(), si.state());
            let fields = forcing.cell_fields();
            for (sum, &i) in sums.iter_mut().zip(&index) {
                for (s, v) in sum.iter_mut().zip(&fields[i].values) {
                    *s += v;
                }
            }
            samples += 1.0;
        }
    })?;
    saved?;
    monthly_result.map_err(|e| format!("cannot write the monthly means: {e}"))?;
    rain_result.map_err(|e| format!("cannot write the monthly rain belt: {e}"))?;
    let rain_belt =
        crate::circulation::RainBelt::of(&band_degrees, &last_months, month_done - first_month);
    let seasonal_range = crate::circulation::SeasonalRange::of(
        &band_degrees,
        &last_temperature,
        month_done - first_month,
    );
    let means: Vec<Vec<f64>> = sums
        .iter()
        .map(|s| s.iter().map(|v| v / samples).collect())
        .collect();
    write_annual_cells(
        &out.dir.join("annual_cells.csv"),
        mesh.cell_centers(),
        &land,
        &height,
        &names,
        &means,
    )?;
    // A prescribed sea has no transport.
    let transport = gained_start
        .filter(|_| !fixed)
        .map(|(from, gained)| {
            let fields = si.model().forcing().expect("forcing").cell_fields();
            let seconds = (si.steps() - from) as f64 * dt;
            let heating: [Vec<f64>; 4] = std::array::from_fn(|j| {
                fields[gained_index[j]]
                    .values
                    .iter()
                    .zip(&gained[j])
                    .map(|(now, then)| (now - then) / seconds)
                    .collect()
            });
            ocean_heat_transport(
                si.model(),
                &lat_bands,
                &land,
                &heating,
                &out.dir,
                args.image.width,
            )
        })
        .transpose()?;
    let overturning = climate.overturning();
    let (upper, lower) = crate::layered_run::upper_lower(si.model());
    let all_layers = climate.means();
    let circulation = fms_atmos::Circulation::of(&all_layers, &overturning, upper, lower);
    let contrast = fms_atmos::circulation::equator_pole_contrast(
        &all_layers.latitude,
        &all_layers.temperature[lower],
    );
    let zonal = crate::layered_run::upper_and_lower(si.model(), all_layers);
    let zonal_physics: Vec<Vec<f64>> = means.iter().map(|m| lat_bands.mean(geometry, m)).collect();

    // Annual-mean maps. The gyres' fields are drawn on the cells wet for the
    // gyres alone, centred on zero.
    let gyre_maps = ["ocean.gyre_streamfunction", "ocean.gyre_heating"];
    for ((id, mean), name) in ids.iter().zip(&means).zip(&names) {
        if fixed && OCEAN_MAPS.contains(name) {
            continue;
        }
        if gyre_maps.contains(name) {
            let wet_only: Vec<f64> = mean
                .iter()
                .zip(&land)
                .map(|(&v, &f)| {
                    if 1.0 - f >= WET_OCEAN_FRACTION {
                        v
                    } else {
                        f64::NAN
                    }
                })
                .collect();
            let units = world.fields.meta(id.key()).spec.units;
            write_gyre_map(
                &world.climate,
                &wet_only,
                name,
                units,
                "FINAL-YEAR MEAN",
                &args.image,
                &out.dir,
            )?;
            continue;
        }
        world.fields.write(*id).copy_from_slice(mean);
        write_map(
            &world,
            id.key(),
            None,
            &out.dir,
            "FINAL-YEAR MEAN",
            &args.image,
        )?;
    }

    // Zonal means: CSV and chart.
    let mut csv = String::from(
        "latitude_deg,wind_upper_m_s,wind_lower_m_s,temperature_upper_k,temperature_lower_k,surface_pressure_pa",
    );
    // Physics fields by name (units: `ColumnPhysics::cell_fields`).
    for name in &names {
        let _ = write!(csv, ",{name}");
    }
    csv.push('\n');
    for (i, latitude) in zonal.latitude.iter().enumerate() {
        let _ = write!(
            csv,
            "{},{},{},{},{},{}",
            latitude,
            zonal.wind[0][i],
            zonal.wind[1][i],
            zonal.temperature[0][i],
            zonal.temperature[1][i],
            zonal.surface_pressure[i],
        );
        for z in &zonal_physics {
            let _ = write!(csv, ",{}", z[i]);
        }
        csv.push('\n');
    }
    let path = out.dir.join("zonal_mean.csv");
    std::fs::write(&path, csv).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!("zonal means -> {}", path.display());
    let net: Vec<f64> = zonal_physics[2]
        .iter()
        .zip(&zonal_physics[1])
        .map(|(a, o)| a - o)
        .collect();
    let series = [
        Series {
            title: "SURFACE TEMPERATURE (K)".into(),
            values: &zonal_physics[0],
        },
        Series {
            title: "UPPER-LAYER EASTWARD WIND (M/S)".into(),
            values: &zonal.wind[0],
        },
        Series {
            title: "LOWER-LAYER EASTWARD WIND (M/S)".into(),
            values: &zonal.wind[1],
        },
        Series {
            title: "ABSORBED SUNLIGHT MINUS OUTGOING LONGWAVE (W/M2)".into(),
            values: &net,
        },
        Series {
            title: "PRECIPITATION (MM/DAY)".into(),
            values: &zonal_physics[4],
        },
        Series {
            title: "EVAPORATION (MM/DAY)".into(),
            values: &zonal_physics[5],
        },
        Series {
            title: "CLOUD COVER".into(),
            values: &zonal_physics[6],
        },
    ];
    let axis = TimeAxis {
        times: &zonal.latitude,
        caption: "LATITUDE (DEGREES)",
        hours: false,
    };
    let header = format!("CLIMATE: FINAL-YEAR ZONAL MEANS, {nc} CELLS");
    let (pixels, height) = chart::render_chart(args.image.width, &header, &axis, &series);
    let path = out.dir.join("zonal_mean.png");
    write_png(&path, args.image.width, height, &pixels, &chart::palette())
        .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!("zonal-mean chart -> {}", path.display());
    crate::circulation::write_overturning(
        &overturning,
        &out.dir,
        args.image.width,
        &format!("CLIMATE: FINAL-YEAR OVERTURNING, {nc} CELLS"),
    )?;
    let measures: Vec<(&str, Option<f64>)> = [("equator_pole_contrast_k", Some(contrast))]
        .into_iter()
        .chain(rain_belt.measures())
        .chain(seasonal_range.measures())
        .collect();
    crate::circulation::write_measures(&circulation, &measures, &out.dir)?;
    rain_belt.write_chart(
        &out.dir,
        args.image.width,
        &format!("CLIMATE: RAIN BELT, FINAL 12 MONTHS, {nc} CELLS"),
    )?;

    // Summary.
    let area: f64 = geometry.cell_areas.iter().sum();
    let global = |v: &[f64]| -> f64 {
        v.iter()
            .zip(&geometry.cell_areas)
            .map(|(x, a)| x * a)
            .sum::<f64>()
            / area
    };
    let mut summary = format!("fms-climate run: {description}\n\n");
    let g =
        |name: &str| global(&means[names.iter().position(|n| *n == name).expect("known field")]);
    let (asr, olr) = (g("rad.absorbed_shortwave"), g("rad.outgoing_longwave"));
    let _ = writeln!(
        summary,
        "Final-year global means: surface temperature {:.2} K, absorbed sunlight {asr:.1} W/m², \
         outgoing longwave {olr:.1} W/m², planetary albedo {:.3}",
        g("rad.surface_temperature"),
        1.0 - asr / (world.planet.mean_flux() / 4.0),
    );
    let _ = writeln!(
        summary,
        "Clouds: cover {:.2} (high {:.2}, low {:.2}), radiative effect {:+.1} W/m² shortwave, {:+.1} longwave",
        g("atm.cloud_cover"),
        g("atm.high_cloud"),
        g("atm.low_cloud"),
        g("rad.shortwave_cloud_effect"),
        g("rad.longwave_cloud_effect"),
    );
    let _ = writeln!(
        summary,
        "Water: vapour {:.1} kg/m², relative humidity {:.2} upper / {:.2} lower layer; precipitation {:.2} mm/day \
         ({:.2} convective, {:.2} snow), evaporation {:.2} mm/day",
        g("atm.water_vapour"),
        g("atm.relative_humidity_upper"),
        g("atm.relative_humidity_lower"),
        g("atm.precipitation"),
        g("atm.convective_precipitation"),
        g("atm.snowfall"),
        g("atm.evaporation"),
    );
    let _ = writeln!(
        summary,
        "Surface (W/m²): sunlight absorbed {:.1}, net longwave loss {:.1} ({:.1} under clear skies), \
         latent heat {:.1}, sensible heat {:.1}",
        g("rad.surface_shortwave"),
        g("rad.surface_longwave"),
        g("rad.clear_surface_longwave"),
        g("atm.latent_heat_flux"),
        g("atm.sensible_heat_flux"),
    );
    if let Some(line) = transport {
        summary.push_str(&line);
    }
    if let Some((d0, e0, b0, w0, wb0)) = start {
        let d1 = si.model().diagnostics(si.state());
        let (e1, b1, w1, wb1) = (
            si.energy_budget(),
            si.model()
                .forcing()
                .and_then(|f| f.budget())
                .expect("budget"),
            si.water(),
            si.water_budget(),
        );
        let seconds = (si.steps() - stats_from) as f64 * dt;
        let r2 = si.model().operators().radius().powi(2);
        let wm2 = |j: f64| j / seconds / (area * r2);
        let toa = (b1.absorbed - b0.absorbed) - (b1.outgoing - b0.outgoing);
        // Latent energy of the water the budget closes on (the mean of the two
        // leapfrog levels).
        let atmosphere = (d1.energy - d0.energy) + LATENT_HEAT * (w1 - w0);
        let slab = b1.stored - b0.stored;
        // Friction's and diffusion's lost kinetic energy comes back as heat;
        // what's left of them is round-off.
        let friction = e1.forcing.kinetic - e0.forcing.kinetic;
        let returned = b1.dissipated - b0.dissipated;
        let diffusion = e1.diffusion.total() - e0.diffusion.total();
        let _ = writeln!(
            summary,
            "Final-year energy budget (W/m²): net in at the top of the atmosphere {:+.3} = atmosphere uptake {:+.3} \
             (enthalpy and latent) + {} {:+.3} + residual {:+.4}; friction dissipated {:.3} \
             (returned as heat, net {:+.1e}), diffusion's net {:+.1e}",
            wm2(toa),
            wm2(atmosphere),
            // A prescribed sea is a reservoir that takes up whatever the
            // surface's fluxes bring it.
            if fixed {
                "the prescribed sea's uptake (its temperature held, so any amount)"
            } else {
                "surface uptake"
            },
            wm2(slab),
            wm2(toa - atmosphere - slab + friction + returned + diffusion),
            -wm2(friction),
            wm2(friction + returned),
            wm2(diffusion),
        );
        let exchanged = (b1.evaporated - b0.evaporated) - (b1.precipitated - b0.precipitated);
        let (filter, fixer) = (wb1.filter - wb0.filter, wb1.fixer - wb0.fixer);
        let _ = writeln!(
            summary,
            "Final-year water budget (kg): change {:+.4e} = evaporation − precipitation {:+.4e} + time filter {:+.2e} \
             + fixer {:+.2e} + residual {:+.1e}",
            w1 - w0,
            exchanged,
            filter,
            fixer,
            (w1 - w0) - exchanged - filter - fixer,
        );
        let (soil, runoff) = (b1.soil_water - b0.soil_water, b1.runoff - b0.runoff);
        let snow = b1.snow - b0.snow;
        let ocean = b1.ocean_exchange - b0.ocean_exchange;
        let _ = writeln!(
            summary,
            "Final-year land water (kg): soil change {:+.4e}, snow change {:+.4e}, runoff {:.4e}; atmosphere, land and runoff \
             change {:+.4e} = ocean evaporation − precipitation {:+.4e} + filter and fixer, residual {:+.1e}",
            soil,
            snow,
            runoff,
            (w1 - w0) + soil + snow + runoff,
            ocean,
            (w1 - w0) + soil + snow + runoff - ocean - filter - fixer,
        );
    }
    summary.push_str(&crate::circulation::summary(
        &circulation,
        &overturning,
        Some(contrast),
        Some((&rain_belt, &seasonal_range)),
    ));
    let _ = writeln!(
        summary,
        "Means of {samples} samples ({SAMPLES_PER_DAY} per day) over the final year."
    );
    if let Some(imprint) = &imprint {
        let block_days = args.imprint_block.expect("an imprint check has blocks");
        let results = if imprint.complete_blocks() >= 2 {
            let days = |step: u64| step as f64 * dt / setup.day;
            let results = imprint.results();
            crate::imprint::write(
                imprint,
                &results,
                &world.climate,
                &out.dir,
                &args.image,
                &format!("DAYS {:.0} TO {:.0}", days(imprint_from), days(si.steps())),
            )?;
            results
        } else {
            Vec::new()
        };
        summary.push_str(&crate::imprint::summary(imprint, &results, block_days));
        if !results.is_empty() && land.iter().any(|f| *f > 0.0) {
            summary.push_str(
                "  With land, the time mean is not expected to be zonally symmetric: the residuals \
                 hold the land's pattern as well as any of the grid's.\n",
            );
        }
        if !results.is_empty() && args.restart.is_some() {
            summary.push_str(
                "  The check isn't saved in restart files: it covers this run only, from its \
                 spin-up on.\n",
            );
        }
    }
    let path = out.dir.join("summary.txt");
    std::fs::write(&path, &summary).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    print!("{summary}");
    // Setup runs to the loop's start; output is everything after its end,
    // the time loop's own maps and GIF included.
    let whole = run_started.elapsed();
    let setup = log
        .loop_started
        .map_or(whole, |started| started.duration_since(run_started));
    let output = whole.saturating_sub(setup + log.loop_time);
    eprintln!(
        "timing: whole run {:.1} s: setup {:.1} s, time loop {:.1} s, output {:.1} s",
        whole.as_secs_f64(),
        setup.as_secs_f64(),
        log.loop_time.as_secs_f64(),
        output.as_secs_f64()
    );
    Ok(())
}

/// Writes the final year's northward ocean heat transport by each transport
/// and their total, from each one's mean heating per unit of ocean (W m⁻²,
/// Ekman, eddies, gyres, diffusion), as `ocean_heat_transport.csv` and a
/// chart, and returns the summary's line.
fn ocean_heat_transport(
    model: &fms_atmos::LayeredModel<'_>,
    bands: &LatitudeBands,
    land: &[f64],
    heating: &[Vec<f64>; 4],
    dir: &Path,
    width: usize,
) -> Result<String, String> {
    let ops = model.operators();
    let r2 = ops.radius() * ops.radius();
    let area: Vec<f64> = ops
        .geometry()
        .cell_areas
        .iter()
        .zip(land)
        .map(|(a, f)| (1.0 - f) * a * r2)
        .collect();
    let mut q: Vec<Vec<f64>> = heating
        .iter()
        .map(|h| {
            northward_transport(bands, h, &area)
                .iter()
                .map(|w| w / 1e15)
                .collect()
        })
        .collect();
    let total: Vec<f64> = (0..q[0].len())
        .map(|k| q.iter().map(|c| c[k]).sum())
        .collect();
    q.push(total);
    let latitudes = edge_latitudes(bands);
    let names = ["Ekman", "eddies", "gyres", "diffusion", "total"];
    let mut csv = String::from("latitude_deg,ekman_pw,eddies_pw,gyres_pw,diffusion_pw,total_pw\n");
    for (k, lat) in latitudes.iter().enumerate() {
        let _ = write!(csv, "{lat:.3}");
        for c in &q {
            let _ = write!(csv, ",{:.6}", c[k]);
        }
        csv.push('\n');
    }
    let path = dir.join("ocean_heat_transport.csv");
    std::fs::write(&path, csv).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!("ocean heat transport -> {}", path.display());
    let series: Vec<Series> = names
        .iter()
        .zip(&q)
        .map(|(name, values)| Series {
            title: format!(
                "NORTHWARD OCEAN HEAT TRANSPORT, {} (PW)",
                name.to_uppercase()
            ),
            values,
        })
        .collect();
    let axis = TimeAxis {
        times: &latitudes,
        caption: "LATITUDE (DEGREES)",
        hours: false,
    };
    let (pixels, height) = chart::render_chart(
        width,
        "CLIMATE: FINAL-YEAR NORTHWARD OCEAN HEAT TRANSPORT",
        &axis,
        &series,
    );
    let path = dir.join("ocean_heat_transport.png");
    write_png(&path, width, height, &pixels, &chart::palette())
        .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!("ocean heat transport chart -> {}", path.display());
    let extreme = |values: &[f64], north: bool| {
        let pick = |a: &(usize, &f64), b: &(usize, &f64)| a.1.total_cmp(b.1);
        let (k, v) = if north {
            values.iter().enumerate().max_by(pick)
        } else {
            values.iter().enumerate().min_by(pick)
        }
        .expect("bands");
        format!("{v:+.3} at {:.0}°", latitudes[k])
    };
    let mut line = String::from(
        "Ocean heat transport (PW, final-year mean of every step; most northward / most southward):",
    );
    for (name, values) in names.iter().zip(&q) {
        let _ = write!(
            line,
            " {name} {} / {};",
            extreme(values, true),
            extreme(values, false)
        );
    }
    let equator = latitudes.iter().position(|l| l.abs() < 1e-9);
    if let Some(k) = equator {
        let _ = write!(line, " total at the equator {:+.3};", q[4][k]);
    }
    let _ = writeln!(
        line,
        " at the north pole (closure) {:+.1e}",
        q[4].last().expect("bands")
    );
    Ok(line)
}

/// How the ocean's currents read in the run's description.
fn currents(ocean: &OceanParams) -> String {
    let mut parts = Vec::new();
    match ocean.ekman {
        Some(EkmanParams { sverdrup: true, .. }) => {
            parts.push("Ekman transport (Sverdrup at the equator)".to_string());
        }
        Some(_) => parts.push("Ekman transport".into()),
        None => {}
    }
    if ocean.eddies.is_some() {
        parts.push("eddy transport".into());
    }
    if let Some(gyres) = ocean.gyres {
        parts.push(format!("wind-driven gyres over {} m", gyres.depth));
    }
    match parts.split_last() {
        None => "without currents".into(),
        Some((last, [])) => format!("with {last}"),
        Some((last, rest)) => format!("with {} and {last}", rest.join(", ")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_description_names_the_currents() {
        let default = OceanParams::default();
        assert_eq!(
            currents(&default),
            "with Ekman transport (Sverdrup at the equator), eddy transport and wind-driven gyres over 700 m"
        );
        assert_eq!(currents(&OceanParams::still()), "without currents");
        let plain = OceanParams {
            ekman: Some(EkmanParams {
                sverdrup: false,
                ..EkmanParams::default()
            }),
            gyres: None,
            ..default
        };
        assert_eq!(currents(&plain), "with Ekman transport and eddy transport");
        assert_eq!(
            currents(&OceanParams {
                ekman: None,
                gyres: None,
                ..default
            }),
            "with eddy transport"
        );
        assert_eq!(
            currents(&OceanParams {
                gyres: Some(GyreParams::default()),
                ..OceanParams::still()
            }),
            "with wind-driven gyres over 700 m"
        );
    }

    #[test]
    fn the_gyre_options_set_the_gyres() {
        let ocean = |options: &[&str]| parse(options).ocean().gyres;
        assert_eq!(ocean(&[]), Some(GyreParams::default()));
        assert_eq!(ocean(&["--no-ocean-gyres"]), None);
        assert_eq!(
            ocean(&["--ocean-gyre-depth", "500", "--ocean-gyre-friction", "1e-6"]),
            Some(GyreParams {
                depth: 500.0,
                friction: 1e-6,
                ..GyreParams::default()
            })
        );
    }

    #[test]
    fn two_years_in_true_time_follow_the_last_jump() {
        let spin_up = |jumps, first_year, every| OceanSpinUp {
            jumps,
            first_year,
            every,
            ..OceanSpinUp::default()
        };
        // The default schedule: jumps end years 2–4, so 6 years.
        assert!(check_years(6.0, &spin_up(3, 2, 1)).is_ok());
        let short = check_years(5.0, &spin_up(3, 2, 1)).unwrap_err();
        assert!(short.contains("at least 6"), "{short}");
        assert!(short.contains("ends year 4"), "{short}");
        // As before: N + 3 years for N jumps from the second year on.
        for n in 1..6 {
            let years = f64::from(n + 3);
            assert!(check_years(years, &spin_up(n, 2, 1)).is_ok());
            assert!(check_years(years - 0.5, &spin_up(n, 2, 1)).is_err());
        }
        // Jumps ending years 3, 5 and 7 need 9 years; 4, 5 and 6 need 8.
        assert!(check_years(9.0, &spin_up(3, 3, 2)).is_ok());
        assert!(check_years(8.0, &spin_up(3, 3, 2)).is_err());
        assert!(check_years(8.0, &spin_up(3, 4, 1)).is_ok());
        assert!(check_years(7.0, &spin_up(3, 4, 1)).is_err());
        // No jumps, no limit.
        assert!(check_years(2.0, &spin_up(0, 2, 1)).is_ok());
    }

    /// The climate options parsed from `options`.
    fn try_parse(options: &[&str]) -> Result<ClimateArgs, clap::Error> {
        use clap::Parser;
        #[derive(Parser)]
        struct Cli {
            #[command(flatten)]
            climate: ClimateArgs,
        }
        let argv = std::iter::once("run").chain(options.iter().copied());
        Cli::try_parse_from(argv).map(|cli| cli.climate)
    }

    fn parse(options: &[&str]) -> ClimateArgs {
        try_parse(options).unwrap()
    }

    #[test]
    fn a_prescribed_sea_runs_without_restarts() {
        assert!(parse(&["--fixed-sea-temperature"]).fixed_sea_temperature);
        assert!(!parse(&[]).fixed_sea_temperature);
        for other in ["--restart", "--save-restart"] {
            let error = try_parse(&["--fixed-sea-temperature", other, "x.fmsrst"])
                .err()
                .expect("conflicts");
            assert_eq!(error.kind(), clap::error::ErrorKind::ArgumentConflict);
        }
        let check = parse(&["--imprint-block", "50"]);
        assert_eq!(
            (check.imprint_block, check.imprint_spinup),
            (Some(50.0), 200.0)
        );
        assert_eq!(parse(&[]).imprint_block, None);
    }

    #[test]
    fn the_spin_up_defaults_to_the_preset_for_the_ocean() {
        let none = parse(&[]);
        assert_eq!(none.spin_up(true), (OceanSpinUp::with_currents(), 9.0));
        assert_eq!(none.spin_up(false), (OceanSpinUp::still(), 6.0));
        // The presets themselves.
        let currents = OceanSpinUp::with_currents();
        assert_eq!(
            (currents.jumps, currents.first_year, currents.every),
            (2, 4, 3)
        );
        let still = OceanSpinUp::still();
        assert_eq!((still.jumps, still.first_year, still.every), (3, 2, 1));
    }

    #[test]
    fn a_spin_up_option_replaces_only_its_own_value() {
        let every = parse(&["--jump-every", "2"]);
        let (spin_up, years) = every.spin_up(true);
        assert_eq!(
            spin_up,
            OceanSpinUp {
                every: 2,
                ..OceanSpinUp::with_currents()
            }
        );
        assert_eq!(years, 9.0);
        let (spin_up, years) = every.spin_up(false);
        assert_eq!(
            spin_up,
            OceanSpinUp {
                every: 2,
                ..OceanSpinUp::still()
            }
        );
        assert_eq!(years, 6.0);
        let both = parse(&["--jump-layers", "both", "--years", "12"]);
        assert_eq!(
            both.spin_up(true),
            (
                OceanSpinUp {
                    layers: JumpLayers::Both,
                    ..OceanSpinUp::with_currents()
                },
                12.0
            )
        );
        let rest = parse(&[
            "--jump-years",
            "5",
            "--jump-limit",
            "1",
            "--jump-first-year",
            "3",
        ]);
        let (spin_up, _) = rest.spin_up(false);
        assert_eq!(
            (spin_up.extra_years, spin_up.limit, spin_up.first_year),
            (5.0, 1.0, 3)
        );
        assert_eq!((spin_up.jumps, spin_up.every), (3, 1));
        assert!(try_parse(&["--jump-every", "0"]).is_err());
        assert!(try_parse(&["--jump-first-year", "0"]).is_err());
    }

    #[test]
    fn no_ocean_jumps_is_true_time_for_the_years_given() {
        let true_time = parse(&["--ocean-jumps", "0"]);
        for currents in [true, false] {
            let (spin_up, years) = true_time.spin_up(currents);
            assert_eq!(spin_up.jumps, 0);
            assert_eq!(spin_up.last_jump_year(), None);
            assert_eq!(years, preset(currents).1);
            assert!(check_years(2.0, &spin_up).is_ok());
            assert_eq!(
                describe_spin_up(&spin_up),
                "ocean spin-up: none (true time)"
            );
        }
        let (_, years) = parse(&["--ocean-jumps", "0", "--years", "20"]).spin_up(true);
        assert_eq!(years, 20.0);
    }

    #[test]
    fn both_presets_leave_two_true_time_years() {
        for currents in [true, false] {
            let (spin_up, years) = preset(currents);
            assert!(check_years(years, &spin_up).is_ok());
            assert!(check_years(years - 1.0, &spin_up).is_err());
        }
        // The currents' preset with a third jump, at year 10, needs 12.
        let (spin_up, years) = parse(&["--ocean-jumps", "3"]).spin_up(true);
        let error = check_years(years, &spin_up).unwrap_err();
        assert!(error.contains("at least 12"), "{error}");
    }

    #[test]
    fn the_description_names_the_spin_up() {
        assert_eq!(
            describe_spin_up(&OceanSpinUp::with_currents()),
            "ocean spin-up: 2 deep-layer jumps of 20 years at the ends of years 4 and 7, at most 3 K"
        );
        assert_eq!(
            describe_spin_up(&OceanSpinUp::still()),
            "ocean spin-up: 3 two-layer jumps of 3 years at the ends of years 2, 3 and 4, at most 3 K"
        );
        let one = OceanSpinUp {
            jumps: 1,
            extra_years: 2.5,
            ..OceanSpinUp::with_currents()
        };
        assert_eq!(
            describe_spin_up(&one),
            "ocean spin-up: 1 deep-layer jump of 2.5 years at the end of year 4, at most 3 K"
        );
    }
}
