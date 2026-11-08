//! `fms baroclinic`: run the two-layer dry core on a balanced, vertically
//! sheared zonal flow seeded with small waves, and render the baroclinic
//! storms that grow from it.

use std::time::Instant;

use clap::Args;
use fms_atmos::{
    Gas, ShearedFlow, TwoLayer, TwoLayerSemiImplicit, TwoLayerSemiImplicitOptions, TwoLayerState,
};
use fms_fields::{ColorScale, FieldId, FieldKey, FieldSpec, Location, Quantity};
use fms_geo::Operators;
use fms_geo::coords::to_lat_lon;
use fms_sim::World;

use crate::render::{format_value, write_gif, write_png};
use crate::{ImageArgs, WorldArgs, build_world, capture, ensure_dir, paint};

const DAY: f64 = 86_400.0;

#[derive(Args)]
pub struct BaroclinicArgs {
    #[command(flatten)]
    world: WorldArgs,
    /// Upper-layer minus lower-layer wind at the equator (m s⁻¹). The flow is
    /// unstable above about 11 m s⁻¹ for the default stability.
    #[arg(long, default_value_t = 40.0)]
    shear: f64,
    /// Static stability: upper minus lower layer potential temperature (K).
    #[arg(long, default_value_t = 40.0)]
    stability: f64,
    /// Seed amplitude (K) of each of zonal waves 4–8 at 45° N and S.
    #[arg(long, default_value_t = 0.5)]
    seed: f64,
    /// Simulated days. The core has no dissipation yet (ROADMAP M2 step 7),
    /// so once the storms mature (about day 7) fronts collapse to the grid
    /// scale and the run blows up.
    #[arg(long, default_value_t = 7.0)]
    days: f64,
    /// GIF frames (spread evenly over the run, including the start).
    #[arg(long, default_value_t = 57)]
    frames: u32,
    /// GIF frames per second.
    #[arg(long, default_value_t = 12)]
    fps: u16,
    /// Design maximum wind for the time step (m s⁻¹).
    #[arg(long, default_value_t = 80.0)]
    max_wind: f64,
    #[command(flatten)]
    image: ImageArgs,
}

struct Fields {
    surface_pressure: FieldId<f64>,
    theta_lower: FieldId<f64>,
    vorticity_lower: FieldId<f64>,
    wind_upper: FieldId<f64>,
    wind_lower: FieldId<f64>,
}

fn register(world: &mut World) -> Result<Fields, String> {
    let err = |e: fms_fields::FieldError| e.to_string();
    let spec = |name, description, units, location, color_scale| FieldSpec {
        description,
        color_scale,
        persistent: false,
        ..FieldSpec::new(name, units, location)
    };
    let wind = |name, description| FieldSpec {
        quantity: Quantity::NormalComponent,
        ..spec(
            name,
            description,
            "m s-1",
            Location::CLIMATE_EDGES,
            ColorScale::Sequential,
        )
    };
    let f = &mut world.fields;
    Ok(Fields {
        surface_pressure: f
            .register(
                spec(
                    "dyn.surface_pressure",
                    "Surface pressure",
                    "hPa",
                    Location::CLIMATE_CELLS,
                    ColorScale::Diverging,
                ),
                0.0,
            )
            .map_err(err)?,
        theta_lower: f
            .register(
                spec(
                    "dyn.theta_lower",
                    "Lower-layer potential temperature",
                    "K",
                    Location::CLIMATE_CELLS,
                    ColorScale::Temperature,
                ),
                0.0,
            )
            .map_err(err)?,
        vorticity_lower: f
            .register(
                spec(
                    "dyn.vorticity_lower",
                    "Lower-layer relative vorticity",
                    "s-1",
                    Location::CLIMATE_VERTICES,
                    ColorScale::Diverging,
                ),
                0.0,
            )
            .map_err(err)?,
        wind_upper: f
            .register(wind("dyn.wind_upper", "Upper-layer wind"), 0.0)
            .map_err(err)?,
        wind_lower: f
            .register(wind("dyn.wind_lower", "Lower-layer wind"), 0.0)
            .map_err(err)?,
    })
}

pub fn run(mut args: BaroclinicArgs) -> Result<(), String> {
    let mut world = build_world(&args.world)?;
    let fields = register(&mut world)?;
    let world = world;
    let params = *world.planet.params();
    let (mesh, geometry) = (&world.climate.mesh, &world.climate.geometry);
    let omega = world.planet.rotation_rate();
    let radius = params.radius;

    let model = TwoLayer::new(
        Operators::new(mesh, geometry, radius),
        params.gravity,
        Gas::EARTH_AIR,
        0.5,
        mesh.vertices().iter().map(|p| 2.0 * omega * p.z).collect(),
        vec![0.0; mesh.n_cells()],
    );
    let flow = ShearedFlow {
        speed: [10.0 + args.shear, 10.0],
        stability: args.stability,
        theta_equator: 315.0,
        ps_equator: params.surface_pressure,
        omega,
    };
    let mut state = flow.state(&model);
    // Seed waves 4–8 with fixed phases, in both hemispheres.
    for k in 0..2 {
        for (i, &p) in mesh.cell_centers().iter().enumerate() {
            let (lat, lon) = to_lat_lon(p);
            let band = |center: f64| {
                fms_math::exp(-((lat - center.to_radians()) / 15f64.to_radians()).powi(2))
            };
            let waves: f64 = (4..=8)
                .map(|m| fms_math::cos(f64::from(m) * lon + 0.7 * f64::from(m)))
                .sum();
            state.theta_ps[k][i] += state.ps[i] * args.seed * (band(45.0) + band(-45.0)) * waves;
        }
    }
    let explicit_dt = model.stable_dt(&state, 0.5);
    let dt = TwoLayerSemiImplicit::stable_dt(&model, args.max_wind, 0.6);
    let steps_per_day = (DAY / dt).ceil() as u64;
    let dt = DAY / steps_per_day as f64;
    let total = (args.days * steps_per_day as f64).round() as u64;
    let frame_count = u64::from(args.frames.max(2));
    eprintln!(
        "baroclinic: {} cells, semi-implicit dt = {dt:.0} s ({:.1}× the explicit limit), {total} steps over {} days",
        mesh.n_cells(),
        dt / explicit_dt,
        args.days
    );

    let publish = |state: &TwoLayerState, model: &TwoLayer| {
        world
            .fields
            .write(fields.surface_pressure)
            .iter_mut()
            .zip(&state.ps)
            .for_each(|(o, p)| *o = p / 100.0);
        world
            .fields
            .write(fields.theta_lower)
            .copy_from_slice(&state.theta(1));
        world
            .fields
            .write(fields.vorticity_lower)
            .copy_from_slice(&model.vorticity(state, 1));
        world
            .fields
            .write(fields.wind_upper)
            .copy_from_slice(&state.u[0]);
        world
            .fields
            .write(fields.wind_lower)
            .copy_from_slice(&state.u[1]);
    };

    let key = fields.surface_pressure.key();
    let overlay = (!args.image.no_arrows).then_some(fields.wind_lower.key());
    let mut frames = Vec::with_capacity(frame_count as usize);
    publish(&state, &model);
    frames.push({
        let mut f = capture(&world, key, overlay);
        f.stamp = "DAY 0.0".into();
        f
    });
    let start = model.diagnostics(&state);
    let t = Instant::now();
    let mut si = TwoLayerSemiImplicit::new(
        model,
        state,
        dt,
        TwoLayerSemiImplicitOptions::default(),
        0.5,
    );
    for k in 1..frame_count {
        let target = total * k / (frame_count - 1);
        while si.steps() < target {
            si.step();
        }
        let s = si.state();
        let max_u = s.u.iter().flatten().fold(0.0f64, |a, &u| a.max(u.abs()));
        if std::env::var("FMS_VERBOSE").is_ok() {
            eprintln!(
                "day {:5.2}: max |u| {max_u:.1} m/s",
                si.steps() as f64 * dt / DAY
            );
        }
        if !max_u.is_finite() || s.ps.iter().any(|p| !p.is_finite()) {
            return Err(format!(
                "the model blew up by day {:.1} (no dissipation yet; see ROADMAP M2 step 7)",
                si.steps() as f64 * dt / DAY
            ));
        }
        publish(s, si.model());
        let mut frame = capture(&world, key, overlay);
        frame.stamp = format!("DAY {:.1}", si.steps() as f64 * dt / DAY);
        frames.push(frame);
    }
    let elapsed = t.elapsed();
    let end = si.model().diagnostics(si.state());
    let stats = si.solver_stats();
    eprintln!(
        "ran in {elapsed:.1?} ({:.0} steps/s, {:.1} CG iterations per solve); mass change {:+.1e}, energy change {:+.1e}",
        total as f64 / elapsed.as_secs_f64(),
        stats.iterations as f64 / stats.solves.max(1) as f64,
        end.mass / start.mass - 1.0,
        end.energy / start.energy - 1.0,
    );

    args.image.out_dir = args.image.out_dir.join("baroclinic");
    ensure_dir(&args.image.out_dir)?;
    let (images, colormap, (w, h)) = paint(&world, key, overlay, &frames, &args.image);
    let path = args.image.out_dir.join("dyn.surface_pressure.gif");
    write_gif(&path, w, h, &images, colormap.palette(), args.fps)
        .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!(
        "dyn.surface_pressure, {frame_count} frames -> {}",
        path.display()
    );

    let keys: [FieldKey; 5] = [
        key,
        fields.theta_lower.key(),
        fields.vorticity_lower.key(),
        fields.wind_upper.key(),
        fields.wind_lower.key(),
    ];
    for key in keys {
        let is_wind = key == fields.wind_upper.key() || key == fields.wind_lower.key();
        let overlay = (is_wind && !args.image.no_arrows).then_some(key);
        let mut frame = capture(&world, key, overlay);
        frame.stamp = format!("DAY {:.1}", args.days);
        let (images, colormap, (w, h)) = paint(
            &world,
            key,
            overlay,
            std::slice::from_ref(&frame),
            &args.image,
        );
        let name = world.fields.meta(key).spec.name;
        let path = args.image.out_dir.join(format!("{name}.png"));
        write_png(&path, w, h, &images[0], colormap.palette())
            .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
        let (min, max) = colormap.range();
        println!(
            "{name} [{} .. {}] -> {}",
            format_value(min),
            format_value(max),
            path.display()
        );
    }
    Ok(())
}
