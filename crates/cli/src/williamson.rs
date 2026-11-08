//! `fms williamson`: run a Williamson shallow-water test case on the climate
//! grid and render it, to check the dynamical core by eye.

use std::time::Instant;

use clap::Args;
use fms_atmos::williamson::{Case, DAY, GRAVITY, RADIUS};
use fms_atmos::{Diagnostics, SemiImplicit, SemiImplicitOptions, ShallowWater, State};
use fms_fields::{ColorScale, FieldId, FieldKey, FieldSpec, Location, Quantity};
use fms_geo::Operators;
use fms_sim::World;

use crate::render::{format_value, write_gif, write_png};
use crate::{ImageArgs, WorldArgs, build_world, capture, ensure_dir, paint};

#[derive(Args)]
pub struct WilliamsonArgs {
    #[command(flatten)]
    world: WorldArgs,
    /// Test case: 2 (steady zonal flow), 5 (flow over a mountain) or 6
    /// (Rossby–Haurwitz wave).
    #[arg(long, value_parser = clap::value_parser!(u8).range(2..=6))]
    case: u8,
    /// Case 2 only: angle between the flow and the rotation axis (degrees).
    #[arg(long, default_value_t = 0.0)]
    alpha: f64,
    /// Simulated days (default: the case's standard length, 5, 15 or 14).
    #[arg(long)]
    days: Option<f64>,
    /// GIF frames (spread evenly over the run, including the start).
    #[arg(long, default_value_t = 61)]
    frames: u32,
    /// GIF frames per second.
    #[arg(long, default_value_t = 12)]
    fps: u16,
    /// Step semi-implicitly (leapfrog, implicit gravity waves) instead of with
    /// explicit RK4. Allows a step limited by the wind, not gravity waves.
    #[arg(long)]
    semi_implicit: bool,
    /// Courant number for the time step. Default: 0.5 for RK4 (gravity waves
    /// plus wind), 0.6 semi-implicit (wind and Coriolis only).
    #[arg(long)]
    courant: Option<f64>,
    #[command(flatten)]
    image: ImageArgs,
}

struct Fields {
    surface: FieldId<f64>,
    vorticity: FieldId<f64>,
    wind: FieldId<f64>,
    error: Option<FieldId<f64>>,
}

fn register(world: &mut World, with_error: bool) -> Result<Fields, String> {
    let err = |e: fms_fields::FieldError| e.to_string();
    let spec = |name, description, units, location, color_scale| FieldSpec {
        description,
        color_scale,
        persistent: false,
        ..FieldSpec::new(name, units, location)
    };
    let f = &mut world.fields;
    Ok(Fields {
        surface: f
            .register(
                spec(
                    "sw.surface_height",
                    "Free-surface height h + b",
                    "m",
                    Location::CLIMATE_CELLS,
                    ColorScale::Sequential,
                ),
                0.0,
            )
            .map_err(err)?,
        vorticity: f
            .register(
                spec(
                    "sw.vorticity",
                    "Relative vorticity",
                    "s-1",
                    Location::CLIMATE_VERTICES,
                    ColorScale::Diverging,
                ),
                0.0,
            )
            .map_err(err)?,
        wind: f
            .register(
                FieldSpec {
                    quantity: Quantity::NormalComponent,
                    ..spec(
                        "sw.wind",
                        "Velocity",
                        "m s-1",
                        Location::CLIMATE_EDGES,
                        ColorScale::Sequential,
                    )
                },
                0.0,
            )
            .map_err(err)?,
        error: if with_error {
            Some(
                f.register(
                    spec(
                        "sw.height_error",
                        "Surface height minus the exact solution",
                        "m",
                        Location::CLIMATE_CELLS,
                        ColorScale::Diverging,
                    ),
                    0.0,
                )
                .map_err(err)?,
            )
        } else {
            None
        },
    })
}

pub fn run(mut args: WilliamsonArgs) -> Result<(), String> {
    let case = match args.case {
        2 => Case::SteadyZonal {
            alpha: args.alpha.to_radians(),
        },
        5 => Case::Mountain,
        6 => Case::RossbyHaurwitz,
        n => return Err(format!("no test case {n}; choose 2, 5 or 6")),
    };
    let mut world = build_world(&args.world)?;
    let fields = register(&mut world, matches!(case, Case::SteadyZonal { .. }))?;
    let world = world;
    let (mesh, geometry) = (&world.climate.mesh, &world.climate.geometry);

    let setup = case.setup(mesh, geometry);
    let topography = setup.topography.clone();
    let model = ShallowWater::new(
        Operators::new(mesh, geometry, RADIUS),
        GRAVITY,
        setup.coriolis,
        setup.topography,
    );
    let state = setup.state;
    let days = args.days.unwrap_or(case.days());
    let explicit_dt = model.stable_dt(&state, 0.5);
    let dt_max = if args.semi_implicit {
        SemiImplicit::stable_dt(&model, case.max_wind(), args.courant.unwrap_or(0.6))
    } else {
        model.stable_dt(&state, args.courant.unwrap_or(0.5))
    };
    let steps_per_day = (DAY / dt_max).ceil() as u64;
    let dt = DAY / steps_per_day as f64;
    let total = (days * steps_per_day as f64).round() as u64;
    let frame_count = u64::from(args.frames.max(2));
    eprintln!(
        "{}: {} cells, {}, dt = {dt:.1} s ({:.1}× the explicit limit), {total} steps over {days} days",
        case.name(),
        mesh.n_cells(),
        if args.semi_implicit {
            "semi-implicit"
        } else {
            "RK4"
        },
        dt / explicit_dt,
    );
    let initial = state.clone();
    let mut stepper = if args.semi_implicit {
        let options = SemiImplicitOptions::for_state(&state);
        Stepper::SemiImplicit(Box::new(SemiImplicit::new(model, state, dt, options, 0.5)))
    } else {
        Stepper::Explicit {
            model: Box::new(model),
            state,
        }
    };

    let publish = |state: &State| {
        let mut surface = world.fields.write(fields.surface);
        for ((s, h), b) in surface.iter_mut().zip(&state.h).zip(&topography) {
            *s = h + b;
        }
        if let Some(id) = fields.error {
            let mut error = world.fields.write(id);
            for ((e, s), &p) in error
                .iter_mut()
                .zip(surface.iter())
                .zip(mesh.cell_centers())
            {
                *e = s - case.surface_height(p);
            }
        }
        world
            .fields
            .write(fields.vorticity)
            .copy_from_slice(&model_vorticity(mesh, geometry, state));
        world.fields.write(fields.wind).copy_from_slice(&state.u);
    };

    let key = fields.surface.key();
    let overlay = (!args.image.no_arrows).then_some(fields.wind.key());
    let start = stepper.diagnostics();
    let t = Instant::now();
    let mut frames = Vec::with_capacity(frame_count as usize);
    let mut done = 0;
    publish(&initial);
    frames.push({
        let mut frame = capture(&world, key, overlay);
        frame.stamp = "DAY 0.0".into();
        frame
    });
    // The semi-implicit stepper takes its first step on creation.
    if args.semi_implicit {
        done = 1;
    }
    for k in 1..frame_count {
        let target = total * k / (frame_count - 1);
        while done < target {
            stepper.step(dt);
            done += 1;
        }
        publish(stepper.state());
        let mut frame = capture(&world, key, overlay);
        frame.stamp = format!("DAY {:.1}", done as f64 * dt / DAY);
        frames.push(frame);
    }
    let elapsed = t.elapsed();
    let end = stepper.diagnostics();
    eprintln!(
        "ran in {elapsed:.1?} ({:.0} steps/s); mass change {:+.1e}, energy change {:+.1e}",
        total as f64 / elapsed.as_secs_f64(),
        end.mass / start.mass - 1.0,
        end.energy / start.energy - 1.0,
    );
    if let Stepper::SemiImplicit(si) = &stepper {
        let stats = si.solver_stats();
        eprintln!(
            "Helmholtz: {:.1} conjugate-gradient iterations per step, worst residual {:.1e}",
            stats.iterations as f64 / stats.solves.max(1) as f64,
            stats.worst_residual
        );
    }

    args.image.out_dir = args.image.out_dir.join(case.name());
    ensure_dir(&args.image.out_dir)?;
    let (images, colormap, (w, h)) = paint(&world, key, overlay, &frames, &args.image);
    let path = args.image.out_dir.join("sw.surface_height.gif");
    write_gif(&path, w, h, &images, colormap.palette(), args.fps)
        .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!(
        "sw.surface_height, {frame_count} frames -> {}",
        path.display()
    );

    // Final state of every field, as PNGs.
    let mut keys: Vec<FieldKey> = vec![key, fields.vorticity.key(), fields.wind.key()];
    keys.extend(fields.error.map(|id| id.key()));
    for key in keys {
        let overlay = (key == fields.wind.key() && !args.image.no_arrows).then_some(key);
        let mut frame = capture(&world, key, overlay);
        frame.stamp = format!("DAY {days:.1}");
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

/// Either time-stepping scheme behind one interface.
enum Stepper<'a> {
    Explicit {
        model: Box<ShallowWater<'a>>,
        state: State,
    },
    SemiImplicit(Box<SemiImplicit<'a>>),
}

impl Stepper<'_> {
    fn step(&mut self, dt: f64) {
        match self {
            Self::Explicit { model, state } => model.step(state, dt),
            Self::SemiImplicit(si) => {
                si.step();
            }
        }
    }

    fn state(&self) -> &State {
        match self {
            Self::Explicit { state, .. } => state,
            Self::SemiImplicit(si) => si.state(),
        }
    }

    fn diagnostics(&self) -> Diagnostics {
        match self {
            Self::Explicit { model, state } => model.diagnostics(state),
            Self::SemiImplicit(si) => si.model().diagnostics(si.state()),
        }
    }
}

fn model_vorticity(
    mesh: &fms_geo::SphereMesh,
    geometry: &fms_geo::MeshGeometry,
    state: &State,
) -> Vec<f64> {
    let mut zeta = vec![0.0; mesh.n_vertices()];
    Operators::new(mesh, geometry, RADIUS).curl(&state.u, &mut zeta);
    zeta
}
