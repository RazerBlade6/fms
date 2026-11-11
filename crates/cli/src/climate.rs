//! `fms-climate climate`: the layered model (four layers by default, the bottom one
//! the boundary layer) driven by radiation and moist physics over slab
//! ocean and land tiles and the world's terrain, on the world's planet
//! (orbit, CO₂, surface pressure). Spins up, then averages the final year: zonal means, a global
//! summary with the planet's energy and water budgets, and annual-mean maps,
//! in `output/climate/`.

use std::fmt::Write as _;
use std::io::Write as _;
use std::path::Path;

use clap::Args;
use fms_atmos::column_physics::{
    CLIMATE_LEVELS, LATENT_HEAT, LandParams, OceanSpinUp, initial_state,
    initial_surface_temperature_by_latitude,
};
use fms_atmos::moisture::{CloudParams, MoistureParams};
use fms_atmos::{ColumnPhysics, Forcing, PhysicsParams, ZonalClimatology};
use fms_fields::{ColorScale, FieldId, FieldSpec, Location};
use fms_grid::LatitudeBands;
use fms_ocean::{EddyParams, EkmanParams, IceParams, OceanParams};
use glam::DVec3;

use crate::chart::{self, Series, TimeAxis};
use crate::layered_run::{
    DynamicsArgs, Fields, Output, Plan, SAMPLES_PER_DAY, build_model, describe, integrate, stepper,
    write_map,
};
use crate::ocean_transport::{edge_latitudes, northward_transport};
use crate::render::write_png;
use crate::{ImageArgs, WorldArgs, build_world};

/// Design wind for the step: radiatively driven jets and storms reach
/// ~100 m s⁻¹ (stronger than Held–Suarez's).
const MAX_WIND: f64 = 120.0;

#[derive(Args)]
pub struct ClimateArgs {
    #[command(flatten)]
    world: WorldArgs,
    /// Simulated years; statistics come from the last one.
    #[arg(long, default_value_t = 6.0)]
    years: f64,
    /// Years, from the second on, that end with the ocean jumping forward; 0
    /// for a run in true time. Needs `--years` of at least this plus 3.
    #[arg(long, default_value_t = 3)]
    ocean_jumps: u32,
    /// Years of the ocean's change added by each jump.
    #[arg(long, default_value_t = OceanSpinUp::default().extra_years)]
    jump_years: f64,
    /// Most any ocean temperature may change in one jump (K).
    #[arg(long, default_value_t = OceanSpinUp::default().limit)]
    jump_limit: f64,
    /// Horizontal heat diffusivity of the ocean (m²/s); 0 for none.
    #[arg(long, default_value_t = OceanParams::default().diffusivity)]
    ocean_diffusivity: f64,
    /// Let a mixed layer colder than the deep layer stay on top.
    #[arg(long)]
    no_ocean_convection: bool,
    /// Move the ocean with the wind: Ekman transport in the mixed layer,
    /// returned by the deep layer, with upwelling and downwelling between.
    #[arg(long)]
    ocean_ekman: bool,
    /// With `--ocean-ekman`, keep the paper's Eq. 4 at the equator instead
    /// of blending in the Sverdrup transport of the wind's curl; ignored
    /// without it.
    #[arg(long)]
    no_ocean_sverdrup: bool,
    /// Gent–McWilliams eddy transport: a flow in the mixed layer down the
    /// slope of the isotherms, returned by the deep layer.
    #[arg(long)]
    ocean_eddies: bool,
    /// The deep layer starts at the mixed layer's temperature, but no warmer
    /// than this (K); `inf` for no limit.
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
    #[command(flatten)]
    dynamics: DynamicsArgs,
    #[command(flatten)]
    image: ImageArgs,
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

pub fn run(args: ClimateArgs) -> Result<(), String> {
    if args.years <= 1.0 {
        return Err("--years must be more than 1 (the last year is averaged)".into());
    }
    if args.ocean_jumps > 0 && args.years < f64::from(args.ocean_jumps) + 3.0 {
        return Err(format!(
            "--years must be at least {} with --ocean-jumps {}: one year to start and two after the last jump",
            args.ocean_jumps + 3,
            args.ocean_jumps
        ));
    }
    let mut world = build_world(&args.world)?;
    let fields = Fields::register(&mut world)?;
    let planet_params = *world.planet.params();
    let radiation = PhysicsParams::default().radiation;
    let ocean = OceanParams {
        diffusivity: args.ocean_diffusivity,
        convective_adjustment: !args.no_ocean_convection,
        ekman: args.ocean_ekman.then(|| EkmanParams {
            sverdrup: !args.no_ocean_sverdrup,
            ..EkmanParams::default()
        }),
        eddies: args.ocean_eddies.then(EddyParams::default),
        deep_start_limit: args.ocean_deep_start_limit,
        ice: OceanParams::default().ice.map(|ice| {
            if args.no_ice_albedo {
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
    };
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
        spin_up: OceanSpinUp {
            jumps: args.ocean_jumps,
            extra_years: args.jump_years,
            limit: args.jump_limit,
        },
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
    let model = setup
        .model
        .with_humidity_profile(params.radiation)
        .with_dissipative_heating();
    let (mesh, geometry) = (&world.climate.mesh, &world.climate.geometry);
    let nc = mesh.n_cells();
    let (state, surface_temperature) = initial_state(&model, &world.planet, &params);
    eprintln!("climate: single-column start at {surface_temperature:.1} K");
    let physics = ColumnPhysics::with_land(
        &model,
        &world.planet,
        params,
        0.0,
        &state,
        initial_surface_temperature_by_latitude(
            &model,
            &world.planet,
            &params,
            surface_temperature,
        ),
        world.fields.read(world.core.terrain.land_fraction).to_vec(),
    )
    .with_relief(
        world
            .fields
            .read(world.core.terrain.surface_height_stddev)
            .to_vec(),
        &world.fields.read(world.core.terrain.surface_height),
    );
    // Where each wanted field sits in the physics' list.
    let fields_now = physics.cell_fields();
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
        "ocean.diffusion_heat_gained",
    ]
    .map(|n| {
        fields_now
            .iter()
            .position(|f| f.name == n)
            .expect("the physics reports the ocean's transport")
    });
    let land = world.fields.read(world.core.terrain.land_fraction).to_vec();
    let model = model.with_forcing(Box::new(physics));
    let bands = args
        .bands
        .unwrap_or_else(|| 2 * ((nc as f64 / 10.0).sqrt().round() as usize).max(8));
    let mut climate = ZonalClimatology::new(&model, bands);
    let lat_bands = LatitudeBands::new(mesh, geometry, bands);
    let description = describe(
        &model,
        &state,
        setup.dt,
        (args.years * world.planet.year_seconds() / 86_400.0).round(),
        &format!(
            "; radiation and moisture over a {} m mixed layer on a {} m deep ocean{}{} and a {} m soil store, CO₂ {} ppm",
            params.ocean.mixed_depth,
            params.ocean.deep_depth,
            match params.ocean.ekman {
                Some(EkmanParams { sverdrup: true, .. }) => {
                    " with Ekman transport, Sverdrup at the equator"
                }
                Some(_) => " with Ekman transport",
                None => "",
            },
            if params.ocean.eddies.is_some() {
                " with eddy transport"
            } else {
                ""
            },
            params.land.soil_depth,
            params.co2
        ),
    );
    eprintln!("climate: {description}");

    let dt = setup.dt;
    let per_year = (world.planet.year_seconds() / dt).round() as u64;
    let total = (args.years * per_year as f64).round() as u64;
    let stats_from = total - per_year;
    let movie_from =
        total.saturating_sub((args.movie_days * setup.steps_per_day as f64).round() as u64);
    let out = Output {
        command: "climate",
        dir: args.image.out_dir.join("climate"),
        image: &args.image,
        fps: args.fps,
    };
    let plan = Plan {
        total,
        frames: u64::from(args.frames),
        frames_from: movie_from,
    };
    let mut si = stepper(model, state, dt);
    let mut sums: Vec<Vec<f64>> = vec![vec![0.0; nc]; names.len()];
    let mut samples = 0.0;
    let mut start = None;
    let mut gained_start: Option<(u64, Vec<Vec<f64>>)> = None;
    std::fs::create_dir_all(&out.dir)
        .map_err(|e| format!("cannot create {}: {e}", out.dir.display()))?;
    let path = out.dir.join("monthly_zonal.csv");
    let file = std::fs::File::create(&path)
        .map_err(|e| format!("cannot create {}: {e}", path.display()))?;
    let mut monthly_file = std::io::BufWriter::new(file);
    let mut monthly_result = writeln!(
        monthly_file,
        "year,month,latitude_deg,area_share,ice_cover,ice_volume_m,{}",
        monthly.join(",")
    );
    let band_latitudes = lat_bands.latitudes();
    let band_shares: Vec<f64> = {
        let total: f64 = lat_bands.areas().iter().sum();
        lat_bands.areas().iter().map(|a| a / total).collect()
    };
    let mut month_sums: Vec<Vec<f64>> = vec![vec![0.0; nc]; monthly.len() + 2];
    let mut month_samples = 0.0;
    let mut month_done = 0_u64;
    let steps_in_month = |month: u64| (month as f64 * per_year as f64 / 12.0).round() as u64;
    integrate(&world, &fields, &mut si, &plan, &out, &description, |si| {
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
            for m in &mut month_sums {
                m.fill(0.0);
            }
            month_samples = 0.0;
            month_done += 1;
        }
        if start.is_none() && si.steps() >= stats_from {
            start = Some((
                si.model().diagnostics(si.state()),
                si.energy_budget(),
                forcing.budget().expect("physics keeps a budget"),
                si.water(),
                si.water_budget(),
            ));
            let fields = forcing.cell_fields();
            gained_start = Some((
                si.steps(),
                gained_index
                    .iter()
                    .map(|&k| fields[k].values.clone())
                    .collect(),
            ));
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
    monthly_result.map_err(|e| format!("cannot write the monthly means: {e}"))?;
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
    let transport = gained_start
        .map(|(from, gained)| {
            let fields = si.model().forcing().expect("forcing").cell_fields();
            let seconds = (si.steps() - from) as f64 * dt;
            let heating: [Vec<f64>; 3] = std::array::from_fn(|j| {
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
    let zonal = crate::layered_run::upper_and_lower(si.model(), climate.means());
    let zonal_physics: Vec<Vec<f64>> = means.iter().map(|m| lat_bands.mean(geometry, m)).collect();

    // Annual-mean maps.
    for (id, mean) in ids.iter().zip(&means) {
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

    // Summary.
    let area: f64 = geometry.cell_areas.iter().sum();
    let global = |v: &[f64]| -> f64 {
        v.iter()
            .zip(&geometry.cell_areas)
            .map(|(x, a)| x * a)
            .sum::<f64>()
            / area
    };
    let mut summary = format!("fms-climate climate: {description}\n\n");
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
             (enthalpy and latent) + surface uptake {:+.3} + residual {:+.4}; friction dissipated {:.3} \
             (returned as heat, net {:+.1e}), diffusion's net {:+.1e}",
            wm2(toa),
            wm2(atmosphere),
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
    let _ = writeln!(
        summary,
        "Means of {samples} samples ({SAMPLES_PER_DAY} per day) over the final year."
    );
    let path = out.dir.join("summary.txt");
    std::fs::write(&path, &summary).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    print!("{summary}");
    Ok(())
}

/// Writes the final year's northward ocean heat transport by each transport
/// and their total, from each one's mean heating per unit of ocean (W m⁻²,
/// Ekman, eddies, diffusion), as `ocean_heat_transport.csv` and a chart, and
/// returns the summary's line.
fn ocean_heat_transport(
    model: &fms_atmos::LayeredModel<'_>,
    bands: &LatitudeBands,
    land: &[f64],
    heating: &[Vec<f64>; 3],
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
    let names = ["Ekman", "eddies", "diffusion", "total"];
    let mut csv = String::from("latitude_deg,ekman_pw,eddies_pw,diffusion_pw,total_pw\n");
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
        let _ = write!(line, " total at the equator {:+.3};", q[3][k]);
    }
    let _ = writeln!(
        line,
        " at the north pole (closure) {:+.1e}",
        q[3].last().expect("bands")
    );
    Ok(line)
}
