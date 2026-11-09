//! Shared machinery for commands that run the layered core (`fms
//! baroclinic`, `fms held-suarez`, `fms climate`): model setup from the
//! command line, the checked time loop, the diagnostics log, the blow-up
//! dump, and the maps and GIF of the dynamics fields. "Upper" means the layer
//! holding 250 hPa (the jets), "lower" the lowest layer.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::Args;
use fms_atmos::{
    BlowUp, Gas, HyperdiffusionOptions, LayeredModel, LayeredSemiImplicit,
    LayeredSemiImplicitOptions, LayeredState, Levels,
};
use fms_fields::{ColorScale, FieldId, FieldKey, FieldSpec, Location, Quantity};
use fms_geo::Operators;
use fms_geo::coords::to_lat_lon;
use fms_sim::World;

use crate::render::{format_value, write_gif, write_png};
use crate::run_log::{Record, RunLog};
use crate::{Frame, ImageArgs, capture, coord_label, ensure_dir, paint};

pub const DAY: f64 = 86_400.0;
/// Diagnostics samples per day.
pub const SAMPLES_PER_DAY: u64 = 4;

/// Layers, time step and dissipation options.
#[derive(Args)]
pub struct DynamicsArgs {
    /// Interior layer interfaces in σ = p/p_s, top down, comma-separated
    /// (e.g. `0.4,0.85` for three layers with a thin bottom one). Default:
    /// `0.5` (two layers), or the climate's four layers for `fms climate`.
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub levels: Option<Vec<f64>>,
    /// Design maximum wind for the time step (m s⁻¹). The run stops if a
    /// wind outruns the step (its Courant limit is this / 0.6).
    #[arg(long)]
    pub max_wind: Option<f64>,
    /// Hyperdiffusion order (4 = ∇⁴, 6 = ∇⁶, 8 = ∇⁸).
    #[arg(long, default_value_t = 4)]
    pub diffusion_order: u32,
    /// Hyperdiffusion e-folding time of the shortest resolved wave (hours).
    /// Default: the strongest the time step handles with a 2× margin.
    #[arg(long)]
    pub e_folding: Option<f64>,
    /// Run without hyperdiffusion (inviscid; collapses to the grid scale).
    #[arg(long)]
    pub inviscid: bool,
}

/// A flat layered model on the world's climate grid, with its step.
pub struct Setup<'w> {
    pub model: LayeredModel<'w>,
    pub dt: f64,
    pub steps_per_day: u64,
}

/// The model and step for `args`, with `max_wind` as the design wind and
/// `levels` as the interior interfaces unless the user gave others.
pub fn build_model<'w>(
    world: &'w World,
    args: &DynamicsArgs,
    max_wind: f64,
    levels: &[f64],
    surface_geopotential: Vec<f64>,
) -> Result<Setup<'w>, String> {
    let params = *world.planet.params();
    let (mesh, geometry) = (&world.climate.mesh, &world.climate.geometry);
    let omega = world.planet.rotation_rate();
    let interior = args.levels.as_deref().unwrap_or(levels);
    if interior.len() + 1 > fms_atmos::levels::MAX_LAYERS
        || !interior.windows(2).all(|w| w[0] < w[1])
        || interior.iter().any(|&s| s <= 0.0 || s >= 1.0)
    {
        return Err(format!(
            "--levels must increase strictly within (0, 1), at most {} layers",
            fms_atmos::levels::MAX_LAYERS
        ));
    }
    let mut model = LayeredModel::new(
        Operators::new(mesh, geometry, params.radius),
        params.gravity,
        Gas::EARTH_AIR,
        Levels::new(interior),
        mesh.vertices().iter().map(|p| 2.0 * omega * p.z).collect(),
        surface_geopotential,
    );
    let max_wind = args.max_wind.unwrap_or(max_wind);
    if max_wind <= 0.0 {
        return Err("--max-wind must be positive".into());
    }
    let dt = LayeredSemiImplicit::stable_dt(&model, max_wind, 0.6);
    let steps_per_day = (DAY / dt).ceil() as u64;
    let dt = DAY / steps_per_day as f64;
    if !args.diffusion_order.is_multiple_of(2) || args.diffusion_order < 2 {
        return Err("--diffusion-order must be even and at least 2".into());
    }
    if !args.inviscid {
        let options = match args.e_folding {
            Some(hours) if hours > 0.0 => HyperdiffusionOptions {
                order: args.diffusion_order,
                e_folding: hours * 3600.0,
            },
            Some(_) => return Err("--e-folding must be positive".into()),
            None => HyperdiffusionOptions::strongest_for_step(
                model.operators(),
                args.diffusion_order,
                dt,
            ),
        };
        model = model.with_hyperdiffusion(options);
        let d = model.hyperdiffusion().expect("just set");
        if dt * d.max_rate() >= 1.0 {
            return Err(format!(
                "hyperdiffusion too strong for the {dt:.0} s step: lengthen --e-folding to at least {:.1} h",
                d.options().e_folding / 3600.0 * dt * d.max_rate() * 1.25
            ));
        }
    }
    Ok(Setup {
        model,
        dt,
        steps_per_day,
    })
}

/// One line describing the run, for logs and reports.
pub fn describe(
    model: &LayeredModel,
    state: &LayeredState,
    dt: f64,
    days: f64,
    extra: &str,
) -> String {
    let diffusion = match model.hyperdiffusion() {
        Some(d) => format!(
            "∇{} hyperdiffusion, ν = {:.2e} m^{} s⁻¹ (grid-scale e-folding {:.1} h)",
            d.options().order,
            d.coefficient(),
            d.options().order,
            d.options().e_folding / 3600.0
        ),
        None => "inviscid".into(),
    };
    format!(
        "{} cells, {} layers (interfaces σ = {:?}), semi-implicit dt = {dt:.0} s ({:.1}× the explicit limit), {days} days; {diffusion}{extra}",
        model.operators().mesh().n_cells(),
        model.layers(),
        model.levels().interior(),
        dt / model.stable_dt(state, 0.5),
    )
}

/// The dynamics fields shown on maps.
pub struct Fields {
    surface_pressure: FieldId<f64>,
    theta_lower: FieldId<f64>,
    vorticity_lower: FieldId<f64>,
    wind_upper: FieldId<f64>,
    wind_lower: FieldId<f64>,
}

impl Fields {
    pub fn register(world: &mut World) -> Result<Self, String> {
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
        Ok(Self {
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
                        "Lowest-layer potential temperature",
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
                        "Lowest-layer relative vorticity",
                        "s-1",
                        Location::CLIMATE_VERTICES,
                        ColorScale::Diverging,
                    ),
                    0.0,
                )
                .map_err(err)?,
            wind_upper: f
                .register(wind("dyn.wind_upper", "Wind of the layer at 250 hPa"), 0.0)
                .map_err(err)?,
            wind_lower: f
                .register(wind("dyn.wind_lower", "Lowest-layer wind"), 0.0)
                .map_err(err)?,
        })
    }

    fn keys(&self) -> [FieldKey; 5] {
        [
            self.surface_pressure.key(),
            self.theta_lower.key(),
            self.vorticity_lower.key(),
            self.wind_upper.key(),
            self.wind_lower.key(),
        ]
    }

    /// Copy a state into the fields.
    pub fn publish(&self, world: &World, state: &LayeredState, model: &LayeredModel) {
        let f = &world.fields;
        f.write(self.surface_pressure)
            .iter_mut()
            .zip(&state.ps)
            .for_each(|(o, p)| *o = p / 100.0);
        let (upper, lower) = upper_lower(model);
        f.write(self.theta_lower)
            .copy_from_slice(&state.theta(lower));
        f.write(self.vorticity_lower)
            .copy_from_slice(&model.vorticity(state, lower));
        f.write(self.wind_upper).copy_from_slice(&state.u[upper]);
        f.write(self.wind_lower).copy_from_slice(&state.u[lower]);
    }
}

/// The layers reported as "upper" (the one holding 250 hPa) and "lower"
/// (the lowest).
pub fn upper_lower(model: &LayeredModel) -> (usize, usize) {
    (model.levels().layer_of(0.25), model.layers() - 1)
}

/// Zonal means of just the "upper" and "lower" layers ([`upper_lower`]), as
/// a two-layer table.
pub fn upper_and_lower(
    model: &LayeredModel,
    means: fms_atmos::ZonalMeans,
) -> fms_atmos::ZonalMeans {
    let (upper, lower) = upper_lower(model);
    fms_atmos::ZonalMeans {
        wind: vec![means.wind[upper].clone(), means.wind[lower].clone()],
        temperature: vec![
            means.temperature[upper].clone(),
            means.temperature[lower].clone(),
        ],
        ..means
    }
}

/// Where and how to write images.
pub struct Output<'a> {
    pub command: &'a str,
    pub dir: PathBuf,
    pub image: &'a ImageArgs,
    pub fps: u16,
}

/// What to record while running.
pub struct Plan {
    /// Steps to run (the stepper starts at step 1).
    pub total: u64,
    /// GIF frames, spread evenly over steps `frames_from..=total`.
    pub frames: u64,
    pub frames_from: u64,
}

/// Run `si` to `plan.total`, checking every step. Samples the diagnostics
/// [`SAMPLES_PER_DAY`] times a day (calling `on_sample` each time), records
/// the GIF frames, and at the end writes the GIF, `diagnostics.csv` and maps
/// of the final state to `out.dir`. On a blow-up, writes the dump to
/// `out.dir/blowup` and returns an error.
pub fn integrate(
    world: &World,
    fields: &Fields,
    si: &mut LayeredSemiImplicit,
    plan: &Plan,
    out: &Output,
    setup: &str,
    mut on_sample: impl FnMut(&LayeredSemiImplicit),
) -> Result<RunLog, String> {
    ensure_dir(&out.dir)?;
    let dt = si.dt();
    let per_day = (DAY / dt).round() as u64;
    let sample_every = (per_day / SAMPLES_PER_DAY).max(1);
    let key = fields.surface_pressure.key();
    let overlay = (!out.image.no_arrows).then_some(fields.wind_lower.key());
    let frame_count = plan.frames.max(2);
    let frame_step = |j: u64| {
        plan.frames_from + (plan.total - plan.frames_from.min(plan.total)) * j / (frame_count - 1)
    };
    let mut frames: Vec<Frame> = Vec::with_capacity(frame_count as usize);
    let mut next_frame = 0;
    let capture_due = |si: &LayeredSemiImplicit, frames: &mut Vec<Frame>, next: &mut u64| {
        while *next < frame_count && si.steps() >= frame_step(*next) {
            fields.publish(world, si.state(), si.model());
            let mut frame = capture(world, key, overlay);
            frame.stamp = format!("DAY {:.1}", si.steps() as f64 * dt / DAY);
            frames.push(frame);
            *next += 1;
        }
    };

    let mut log = RunLog::default();
    let health = si
        .check()
        .map_err(|b| format!("initial state is unhealthy: {b}"))?;
    log.push(Record {
        step: si.steps(),
        day: si.steps() as f64 * dt / DAY,
        diagnostics: si.model().diagnostics(si.state()),
        budget: si.energy_budget(),
        health,
    });
    capture_due(si, &mut frames, &mut next_frame);
    let t = Instant::now();
    let mut last_good = si.state().clone();
    while si.steps() < plan.total {
        last_good.clone_from(si.state());
        si.step();
        let day = si.steps() as f64 * dt / DAY;
        let health = match si.check() {
            Ok(h) => h,
            Err(blow_up) => {
                let dump = out.dir.join("blowup");
                write_dump(
                    &dump, world, fields, si, &last_good, &log, &blow_up, setup, out,
                )?;
                write_animation(world, key, overlay, &frames, out)?;
                return Err(format!(
                    "the model blew up on day {day:.2}: {blow_up}. Diagnostic dump in {}",
                    dump.display()
                ));
            }
        };
        if si.steps() % sample_every == 0 || si.steps() == plan.total {
            log.push(Record {
                step: si.steps(),
                day,
                diagnostics: si.model().diagnostics(si.state()),
                budget: si.energy_budget(),
                health,
            });
            on_sample(si);
            if std::env::var("FMS_VERBOSE").is_ok() {
                eprintln!("day {day:6.2}: max |u| {:.1} m/s", health.max_wind);
            }
        }
        capture_due(si, &mut frames, &mut next_frame);
    }
    let elapsed = t.elapsed();
    let stats = si.solver_stats();
    let (first, last) = (log.first().copied(), log.last().copied());
    if let (Some(first), Some(last)) = (first, last) {
        let (dm, _, de) = log.drift(&last);
        eprintln!(
            "ran in {elapsed:.1?} ({:.0} steps/s, {:.1} CG iterations per solve); mass change {dm:+.1e}, \
             energy {:+.1e} (budget residual {de:+.1e})",
            (plan.total - first.step) as f64 / elapsed.as_secs_f64(),
            stats.iterations as f64 / stats.solves.max(1) as f64,
            last.diagnostics.energy / first.diagnostics.energy - 1.0,
        );
    }
    write_animation(world, key, overlay, &frames, out)?;
    let csv = out.dir.join("diagnostics.csv");
    log.write_csv(&csv)?;
    println!("diagnostics -> {}", csv.display());
    fields.publish(world, si.state(), si.model());
    write_maps(
        world,
        fields,
        &out.dir,
        &format!("DAY {:.1}", si.steps() as f64 * dt / DAY),
        out.image,
    )?;
    Ok(log)
}

fn write_animation(
    world: &World,
    key: FieldKey,
    overlay: Option<FieldKey>,
    frames: &[Frame],
    out: &Output,
) -> Result<(), String> {
    if frames.is_empty() {
        return Ok(());
    }
    let (images, colormap, (w, h)) = paint(world, key, overlay, frames, out.image);
    let path = out.dir.join("dyn.surface_pressure.gif");
    write_gif(&path, w, h, &images, colormap.palette(), out.fps)
        .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!(
        "dyn.surface_pressure, {} frames -> {}",
        frames.len(),
        path.display()
    );
    Ok(())
}

/// A PNG of each field, as currently published.
fn write_maps(
    world: &World,
    fields: &Fields,
    dir: &Path,
    stamp: &str,
    image: &ImageArgs,
) -> Result<(), String> {
    for key in fields.keys() {
        let is_wind = key == fields.wind_upper.key() || key == fields.wind_lower.key();
        let overlay = (is_wind && !image.no_arrows).then_some(key);
        write_map(world, key, overlay, dir, stamp, image)?;
    }
    Ok(())
}

/// A PNG of one field, as currently published, named after it.
pub fn write_map(
    world: &World,
    key: FieldKey,
    overlay: Option<FieldKey>,
    dir: &Path,
    stamp: &str,
    image: &ImageArgs,
) -> Result<(), String> {
    let mut frame = capture(world, key, overlay);
    frame.stamp = stamp.into();
    let (images, colormap, (w, h)) =
        paint(world, key, overlay, std::slice::from_ref(&frame), image);
    let name = world.fields.meta(key).spec.name;
    let path = dir.join(format!("{name}.png"));
    write_png(&path, w, h, &images[0], colormap.palette())
        .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    let (min, max) = colormap.range();
    println!(
        "{name} [{} .. {}] -> {}",
        format_value(min),
        format_value(max),
        path.display()
    );
    Ok(())
}

/// The blow-up dump: `report.txt`, `diagnostics.csv`, and maps of the last
/// healthy state.
#[allow(clippy::too_many_arguments)]
fn write_dump(
    dir: &Path,
    world: &World,
    fields: &Fields,
    si: &LayeredSemiImplicit,
    last_good: &LayeredState,
    log: &RunLog,
    blow_up: &BlowUp,
    setup: &str,
    out: &Output,
) -> Result<(), String> {
    ensure_dir(dir)?;
    let geometry = &world.climate.geometry;
    let point = if blow_up.variable.at_edges() {
        geometry.edge_points[blow_up.index]
    } else {
        world.climate.mesh.cell_centers()[blow_up.index]
    };
    let (lat, lon) = to_lat_lon(point);
    let dt = si.dt();
    let day = si.steps() as f64 * dt / DAY;
    let limits = si.limits();
    let mut report = String::new();
    let _ = writeln!(
        report,
        "fms {} blew up at step {} (day {day:.3})\n",
        out.command,
        si.steps()
    );
    let _ = writeln!(report, "Problem:  {blow_up}");
    let _ = writeln!(
        report,
        "Where:    {}",
        coord_label(lat.to_degrees(), lon.to_degrees(), " ")
    );
    let _ = writeln!(report, "Run:      {setup}");
    let _ = writeln!(
        report,
        "Limits:   |u| ≤ {:.1} m/s (Courant 1), p_s {:.0}–{:.0} hPa, θ {:.0}–{:.0} K\n",
        limits.max_wind,
        limits.surface_pressure.0 / 100.0,
        limits.surface_pressure.1 / 100.0,
        limits.theta.0,
        limits.theta.1
    );
    if let Ok(h) = fms_atmos::safety::check(last_good, &limits) {
        let (e, layer) = h.max_wind_at;
        let (wlat, wlon) = to_lat_lon(geometry.edge_points[e]);
        let _ = writeln!(
            report,
            "Last healthy state (step {}): max |u| {:.1} m/s ({} layer, at {}), p_s {:.0}–{:.0} hPa, θ {:.0}–{:.0} K\n",
            si.steps() - 1,
            h.max_wind,
            if layer == 0 { "upper" } else { "lower" },
            coord_label(wlat.to_degrees(), wlon.to_degrees(), " "),
            h.surface_pressure.0 / 100.0,
            h.surface_pressure.1 / 100.0,
            h.theta.0,
            h.theta.1
        );
    }
    let _ = writeln!(
        report,
        "Recent diagnostics (residual = energy change not accounted for by forcing and dissipation):"
    );
    report.push_str(&log.tail_table(12));
    let _ = writeln!(
        report,
        "\nMaps of the last healthy state are in this directory."
    );
    let path = dir.join("report.txt");
    std::fs::write(&path, &report).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    eprintln!("{report}");
    log.write_csv(&dir.join("diagnostics.csv"))?;
    fields.publish(world, last_good, si.model());
    write_maps(
        world,
        fields,
        dir,
        &format!("LAST HEALTHY DAY {:.2}", day - dt / DAY),
        out.image,
    )
}

/// The default semi-implicit options.
pub fn stepper<'w>(
    model: LayeredModel<'w>,
    state: LayeredState,
    dt: f64,
) -> LayeredSemiImplicit<'w> {
    LayeredSemiImplicit::new(model, state, dt, LayeredSemiImplicitOptions::default(), 0.5)
}
