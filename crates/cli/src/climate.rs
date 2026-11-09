//! `fms climate`: the layered model (four layers by default, the bottom one
//! the boundary layer) driven by radiation and moist physics over slab
//! ocean and land tiles and the world's terrain, on the world's planet
//! (orbit, CO₂, surface pressure). Spins up, then averages the final year: zonal means, a global
//! summary with the planet's energy and water budgets, and annual-mean maps,
//! in `output/climate/`.

use std::fmt::Write as _;

use clap::Args;
use fms_atmos::column_physics::{
    CLIMATE_LEVELS, LATENT_HEAT, initial_state, initial_surface_temperature,
};
use fms_atmos::{ColumnPhysics, Forcing, PhysicsParams, ZonalClimatology};
use fms_fields::{ColorScale, FieldId, FieldSpec, Location};
use fms_geo::LatitudeBands;

use crate::chart::{self, Series, TimeAxis};
use crate::layered_run::{
    DynamicsArgs, Fields, Output, Plan, SAMPLES_PER_DAY, build_model, describe, integrate, stepper,
    write_map,
};
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
    #[arg(long, default_value_t = 4.0)]
    years: f64,
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

pub fn run(args: ClimateArgs) -> Result<(), String> {
    if args.years <= 1.0 {
        return Err("--years must be more than 1 (the last year is averaged)".into());
    }
    let mut world = build_world(&args.world)?;
    let fields = Fields::register(&mut world)?;
    let planet_params = *world.planet.params();
    let params = PhysicsParams {
        co2: planet_params.co2_ppm,
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
    ];
    let names: Vec<&str> = maps.iter().map(|m| m.0).chain(extra).collect();
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
    let setup = build_model(
        &world,
        &args.dynamics,
        MAX_WIND,
        &CLIMATE_LEVELS,
        world.surface_geopotential(),
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
        initial_surface_temperature(&model, &params, surface_temperature),
        world.fields.read(world.core.terrain.land_fraction).to_vec(),
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
            "; radiation and moisture over a {} m ocean slab and a {} m land slab, CO₂ {} ppm",
            params.slab_depth, params.land.slab_depth, params.co2
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
    integrate(&world, &fields, &mut si, &plan, &out, &description, |si| {
        let forcing = si.model().forcing().expect("climate has forcing");
        if start.is_none() && si.steps() >= stats_from {
            start = Some((
                si.model().diagnostics(si.state()),
                si.energy_budget(),
                forcing.budget().expect("physics keeps a budget"),
                si.water(),
                si.water_budget(),
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
    let means: Vec<Vec<f64>> = sums
        .iter()
        .map(|s| s.iter().map(|v| v / samples).collect())
        .collect();
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
    let mut summary = format!("fms climate: {description}\n\n");
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
             (enthalpy and latent) + slab uptake {:+.3} + residual {:+.4}; friction dissipated {:.3} \
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
