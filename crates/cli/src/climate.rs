//! `fms climate`: the two-layer model driven by radiation over a slab
//! aquaplanet, on the world's planet (orbit, CO₂, surface pressure). Spins
//! up, then averages the final year: zonal means, a global summary with the
//! planet's energy budget, and annual-mean maps, in `output/climate/`.

use std::fmt::Write as _;

use clap::Args;
use fms_atmos::{
    ColumnPhysics, PhysicsParams, SingleColumn, TwoLayerState, ZonalClimatology,
    column_physics::layer_exner,
};
use fms_fields::{ColorScale, FieldId, FieldSpec, Location};
use fms_geo::LatitudeBands;

use crate::chart::{self, Series, TimeAxis};
use crate::render::write_png;
use crate::two_layer_run::{
    DynamicsArgs, Fields, Output, Plan, SAMPLES_PER_DAY, build_model, describe, integrate, stepper,
    write_map,
};
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
    /// Albedo of the stand-in cloud deck (default: Earth-calibrated).
    #[arg(long)]
    cloud_albedo: Option<f64>,
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
    let mut params = PhysicsParams {
        co2: planet_params.co2_ppm,
        ..PhysicsParams::default()
    };
    if let Some(a) = args.cloud_albedo {
        params.radiation.cloud_albedo = a;
    }
    // Annual-mean maps of the physics' cell fields (in the order
    // `ColumnPhysics::cell_fields` lists them).
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
            "rad.water_vapour",
            "Water-vapour column",
            "kg m-2",
            ColorScale::Sequential,
        ),
        (
            "rad.latent_heat_flux",
            "Surface latent heat flux",
            "W m-2",
            ColorScale::Sequential,
        ),
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

    // Start from the single column's equilibrium, at rest.
    let mut column = SingleColumn::earth(world.planet.mean_flux());
    column.params = PhysicsParams {
        slab_depth: 1.0,
        ..params
    };
    column.ps = planet_params.surface_pressure;
    column.gravity = planet_params.gravity;
    let eq = column.equilibrium();
    eprintln!(
        "climate: single-column start at {:.1} K (layers {:.1}/{:.1} K)",
        eq.surface_temperature, eq.temperature[0], eq.temperature[1]
    );
    let setup = build_model(&world, &args.dynamics, MAX_WIND)?;
    let model = setup.model;
    let (mesh, geometry) = (&world.climate.mesh, &world.climate.geometry);
    let (nc, ne) = (mesh.n_cells(), mesh.n_edges());
    let ps = planet_params.surface_pressure;
    let exner = layer_exner(model.gas(), model.sigma(), ps);
    let mut x = 0x2545_F491_4F6C_DD1Du64;
    let theta: [Vec<f64>; 2] = std::array::from_fn(|k| {
        (0..nc)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                eq.temperature[k] / exner[k] + 0.1 * ((x >> 11) as f64 / (1u64 << 52) as f64 - 1.0)
            })
            .collect()
    });
    let state = TwoLayerState::new(vec![ps; nc], [vec![0.0; ne], vec![0.0; ne]], theta);
    let physics = ColumnPhysics::new(
        &model,
        &world.planet,
        params,
        0.0,
        &state,
        vec![eq.surface_temperature; nc],
    );
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
            "; radiation over a {} m slab, CO₂ {} ppm, cloud albedo {}",
            params.slab_depth, params.co2, params.radiation.cloud_albedo
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
    let mut sums: Vec<Vec<f64>> = vec![vec![0.0; nc]; ids.len()];
    let mut samples = 0.0;
    let mut start = None;
    integrate(&world, &fields, &mut si, &plan, &out, &description, |si| {
        let forcing = si.model().forcing().expect("climate has forcing");
        if start.is_none() && si.steps() >= stats_from {
            start = Some((
                si.model().diagnostics(si.state()),
                si.energy_budget(),
                forcing.budget().expect("physics keeps a budget"),
            ));
        }
        if si.steps() > stats_from {
            climate.sample(si.model(), si.state());
            for (sum, field) in sums.iter_mut().zip(forcing.cell_fields()) {
                for (s, v) in sum.iter_mut().zip(&field.values) {
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
    let zonal = climate.means();
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
        "latitude_deg,wind_upper_m_s,wind_lower_m_s,temperature_upper_k,temperature_lower_k,surface_pressure_pa,surface_temperature_k,outgoing_longwave_w_m2,absorbed_shortwave_w_m2\n",
    );
    for (i, latitude) in zonal.latitude.iter().enumerate() {
        let _ = writeln!(
            csv,
            "{},{},{},{},{},{},{},{},{}",
            latitude,
            zonal.wind[0][i],
            zonal.wind[1][i],
            zonal.temperature[0][i],
            zonal.temperature[1][i],
            zonal.surface_pressure[i],
            zonal_physics[0][i],
            zonal_physics[1][i],
            zonal_physics[2][i]
        );
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
            title: "WATER VAPOUR (KG/M2)".into(),
            values: &zonal_physics[3],
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
    let (asr, olr) = (global(&means[2]), global(&means[1]));
    let _ = writeln!(
        summary,
        "Final-year global means: surface temperature {:.2} K, absorbed sunlight {asr:.1} W/m², \
         outgoing longwave {olr:.1} W/m², planetary albedo {:.3}, water vapour {:.1} kg/m², \
         latent heat flux {:.1} W/m²",
        global(&means[0]),
        1.0 - asr / (world.planet.mean_flux() / 4.0),
        global(&means[3]),
        global(&means[4]),
    );
    if let Some((d0, e0, b0)) = start {
        let d1 = si.model().diagnostics(si.state());
        let (e1, b1) = (
            si.energy_budget(),
            si.model()
                .forcing()
                .and_then(|f| f.budget())
                .expect("budget"),
        );
        let seconds = (si.steps() - stats_from) as f64 * dt;
        let r2 = si.model().operators().radius().powi(2);
        let wm2 = |j: f64| j / seconds / (area * r2);
        let toa = (b1.absorbed - b0.absorbed) - (b1.outgoing - b0.outgoing);
        let atmosphere = d1.energy - d0.energy;
        let slab = b1.stored - b0.stored;
        let friction = e1.forcing.kinetic - e0.forcing.kinetic;
        let diffusion = e1.diffusion.total() - e0.diffusion.total();
        let _ = writeln!(
            summary,
            "Final-year energy budget (W/m²): net in at the top of the atmosphere {:+.3} = atmosphere uptake {:+.3} \
             + slab uptake {:+.3} + kinetic energy removed by friction {:.3} and diffusion {:.3} + residual {:+.4}",
            wm2(toa),
            wm2(atmosphere),
            wm2(slab),
            -wm2(friction),
            -wm2(diffusion),
            wm2(toa - atmosphere - slab + friction + diffusion),
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
