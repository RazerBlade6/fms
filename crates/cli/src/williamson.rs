//! `fms williamson`: run a Williamson shallow-water test case on the climate
//! grid and render it, to check the dynamical core by eye.

use std::time::Instant;

use clap::Args;
use fms_atmos::williamson::{Case, DAY, GRAVITY, RADIUS};
use fms_atmos::{ShallowWater, State};
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
    /// Courant number for the time step (RK4 is stable to about 0.7).
    #[arg(long, default_value_t = 0.5)]
    courant: f64,
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
    let mut model = ShallowWater::new(
        Operators::new(mesh, geometry, RADIUS),
        GRAVITY,
        setup.coriolis,
        setup.topography,
    );
    let mut state = setup.state;
    let days = args.days.unwrap_or(case.days());
    let dt_max = model.stable_dt(&state, args.courant);
    let steps_per_day = (DAY / dt_max).ceil() as u64;
    let dt = DAY / steps_per_day as f64;
    let total = (days * steps_per_day as f64).round() as u64;
    let frame_count = u64::from(args.frames.max(2));
    eprintln!(
        "{}: {} cells, dt = {dt:.1} s, {total} steps over {days} days",
        case.name(),
        mesh.n_cells()
    );

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
    let start = model.diagnostics(&state);
    let t = Instant::now();
    let mut frames = Vec::with_capacity(frame_count as usize);
    let mut done = 0;
    for k in 0..frame_count {
        let target = total * k / (frame_count - 1);
        for _ in done..target {
            model.step(&mut state, dt);
        }
        done = target;
        publish(&state);
        let mut frame = capture(&world, key, overlay);
        frame.stamp = format!("DAY {:.1}", done as f64 * dt / DAY);
        frames.push(frame);
    }
    let elapsed = t.elapsed();
    let end = model.diagnostics(&state);
    eprintln!(
        "ran in {elapsed:.1?} ({:.0} steps/s); mass change {:+.1e}, energy change {:+.1e}",
        total as f64 / elapsed.as_secs_f64(),
        end.mass / start.mass - 1.0,
        end.energy / start.energy - 1.0,
    );

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

fn model_vorticity(
    mesh: &fms_geo::SphereMesh,
    geometry: &fms_geo::MeshGeometry,
    state: &State,
) -> Vec<f64> {
    let mut zeta = vec![0.0; mesh.n_vertices()];
    Operators::new(mesh, geometry, RADIUS).curl(&state.u, &mut zeta);
    zeta
}
