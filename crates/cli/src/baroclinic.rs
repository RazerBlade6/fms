//! `fms baroclinic`: run the two-layer dry core on a balanced, vertically
//! sheared zonal flow seeded with small waves, and render the baroclinic
//! storms that grow from it.
//!
//! The state is checked after every step ([`fms_atmos::safety`]). If it blows
//! up, the run stops and writes a diagnostic dump to `baroclinic/blowup/`: a
//! report, the diagnostics history, and maps of the last healthy state.

use std::fmt::Write as _;
use std::path::Path;
use std::time::Instant;

use clap::Args;
use fms_atmos::{
    BlowUp, Gas, HyperdiffusionOptions, ShearedFlow, TwoLayer, TwoLayerSemiImplicit,
    TwoLayerSemiImplicitOptions, TwoLayerState,
};
use fms_fields::{ColorScale, FieldId, FieldKey, FieldSpec, Location, Quantity};
use fms_geo::Operators;
use fms_geo::coords::to_lat_lon;
use fms_sim::World;

use crate::render::{format_value, write_gif, write_png};
use crate::run_log::{Record, RunLog};
use crate::{Frame, ImageArgs, WorldArgs, build_world, capture, coord_label, ensure_dir, paint};

const DAY: f64 = 86_400.0;
/// Diagnostics samples per day.
const SAMPLES_PER_DAY: u64 = 4;

#[derive(Args)]
pub struct BaroclinicArgs {
    #[command(flatten)]
    world: WorldArgs,
    /// Upper-layer minus lower-layer wind at the equator (m s⁻¹). The flow is
    /// unstable above about 11 m s⁻¹ for the default stability.
    #[arg(long, default_value_t = 30.0)]
    shear: f64,
    /// Static stability: upper minus lower layer potential temperature (K).
    #[arg(long, default_value_t = 40.0)]
    stability: f64,
    /// Seed amplitude (K) of each of zonal waves 4–8 at 45° N and S.
    #[arg(long, default_value_t = 0.5)]
    seed: f64,
    /// Simulated days.
    #[arg(long, default_value_t = 30.0)]
    days: f64,
    /// GIF frames (spread evenly over the run, including the start).
    #[arg(long, default_value_t = 121)]
    frames: u32,
    /// GIF frames per second.
    #[arg(long, default_value_t = 12)]
    fps: u16,
    /// Design maximum wind for the time step (m s⁻¹). Storms on the default
    /// jet peak near 100 m s⁻¹; the run stops if a wind outruns the step.
    #[arg(long, default_value_t = 120.0)]
    max_wind: f64,
    /// Hyperdiffusion order (4 = ∇⁴, 6 = ∇⁶, 8 = ∇⁸).
    #[arg(long, default_value_t = 4)]
    diffusion_order: u32,
    /// Hyperdiffusion e-folding time of the shortest resolved wave (hours).
    /// Default: the strongest the time step handles with a 2× margin.
    #[arg(long)]
    e_folding: Option<f64>,
    /// Run without hyperdiffusion (inviscid; collapses to the grid scale).
    #[arg(long)]
    inviscid: bool,
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

impl Fields {
    fn keys(&self) -> [FieldKey; 5] {
        [
            self.surface_pressure.key(),
            self.theta_lower.key(),
            self.vorticity_lower.key(),
            self.wind_upper.key(),
            self.wind_lower.key(),
        ]
    }
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

fn publish(world: &World, fields: &Fields, state: &TwoLayerState, model: &TwoLayer) {
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
}

pub fn run(mut args: BaroclinicArgs) -> Result<(), String> {
    let mut world = build_world(&args.world)?;
    let fields = register(&mut world)?;
    let world = world;
    let params = *world.planet.params();
    let (mesh, geometry) = (&world.climate.mesh, &world.climate.geometry);
    let omega = world.planet.rotation_rate();
    let radius = params.radius;

    let mut model = TwoLayer::new(
        Operators::new(mesh, geometry, radius),
        params.gravity,
        Gas::EARTH_AIR,
        0.5,
        mesh.vertices().iter().map(|p| 2.0 * omega * p.z).collect(),
        vec![0.0; mesh.n_cells()],
    );
    let dt = TwoLayerSemiImplicit::stable_dt(&model, args.max_wind, 0.6);
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
    }
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
    if let Some(d) = model.hyperdiffusion()
        && dt * d.max_rate() >= 1.0
    {
        return Err(format!(
            "hyperdiffusion too strong for the {dt:.0} s step: lengthen --e-folding to at least {:.1} h",
            d.options().e_folding / 3600.0 * dt * d.max_rate() * 1.25
        ));
    }
    let total = (args.days * steps_per_day as f64).round() as u64;
    let frame_count = u64::from(args.frames.max(2));
    let sample_every = (steps_per_day / SAMPLES_PER_DAY).max(1);
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
    let setup = format!(
        "{} cells, semi-implicit dt = {dt:.0} s ({:.1}× the explicit limit), {total} steps over {} days; {diffusion}",
        mesh.n_cells(),
        dt / explicit_dt,
        args.days
    );
    eprintln!("baroclinic: {setup}");

    let key = fields.surface_pressure.key();
    let overlay = (!args.image.no_arrows).then_some(fields.wind_lower.key());
    let mut frames = Vec::with_capacity(frame_count as usize);
    publish(&world, &fields, &state, &model);
    frames.push({
        let mut f = capture(&world, key, overlay);
        f.stamp = "DAY 0.0".into();
        f
    });
    let t = Instant::now();
    let mut log = RunLog::default();
    let start_diagnostics = model.diagnostics(&state);
    let mut si = TwoLayerSemiImplicit::new(
        model,
        state.clone(),
        dt,
        TwoLayerSemiImplicitOptions::default(),
        0.5,
    );
    let health = si
        .check()
        .map_err(|b| format!("initial state is unhealthy: {b}"))?;
    log.push(Record {
        step: 0,
        day: 0.0,
        diagnostics: start_diagnostics,
        dissipated: 0.0,
        health,
    });
    args.image.out_dir = args.image.out_dir.join("baroclinic");
    ensure_dir(&args.image.out_dir)?;

    // The last healthy state, for the dump.
    let mut last_good = state;
    let mut next_frame = 1;
    while si.steps() < total {
        last_good.clone_from(si.state());
        si.step();
        let day = si.steps() as f64 * dt / DAY;
        let health = match si.check() {
            Ok(h) => h,
            Err(blow_up) => {
                let dump = args.image.out_dir.join("blowup");
                write_dump(
                    &dump,
                    &world,
                    &fields,
                    &si,
                    &last_good,
                    &log,
                    &blow_up,
                    &setup,
                    &args.image,
                )?;
                write_animation(&world, key, overlay, &frames, &args)?;
                return Err(format!(
                    "the model blew up on day {day:.2}: {blow_up}. Diagnostic dump in {}",
                    dump.display()
                ));
            }
        };
        if si.steps() % sample_every == 0 || si.steps() == total {
            log.push(Record {
                step: si.steps(),
                day,
                diagnostics: si.model().diagnostics(si.state()),
                dissipated: si.dissipated_energy(),
                health,
            });
            if std::env::var("FMS_VERBOSE").is_ok() {
                eprintln!("day {day:6.2}: max |u| {:.1} m/s", health.max_wind);
            }
        }
        while next_frame < frame_count && si.steps() >= total * next_frame / (frame_count - 1) {
            publish(&world, &fields, si.state(), si.model());
            let mut frame = capture(&world, key, overlay);
            frame.stamp = format!("DAY {day:.1}");
            frames.push(frame);
            next_frame += 1;
        }
    }
    let elapsed = t.elapsed();
    let stats = si.solver_stats();
    let last = log.last().copied().expect("the log has the start");
    let (dm, _, de) = log.drift(&last);
    eprintln!(
        "ran in {elapsed:.1?} ({:.0} steps/s, {:.1} CG iterations per solve); mass change {dm:+.1e}, \
         energy {:+.1e} of which dissipated {:+.1e} (residual {de:+.1e})",
        total as f64 / elapsed.as_secs_f64(),
        stats.iterations as f64 / stats.solves.max(1) as f64,
        last.diagnostics.energy / start_diagnostics.energy - 1.0,
        si.dissipated_energy() / start_diagnostics.energy,
    );

    write_animation(&world, key, overlay, &frames, &args)?;
    let csv = args.image.out_dir.join("diagnostics.csv");
    log.write_csv(&csv)?;
    println!("diagnostics -> {}", csv.display());
    publish(&world, &fields, si.state(), si.model());
    write_maps(
        &world,
        &fields,
        &args.image.out_dir,
        &format!("DAY {:.1}", args.days),
        &args.image,
    )
}

fn write_animation(
    world: &World,
    key: FieldKey,
    overlay: Option<FieldKey>,
    frames: &[Frame],
    args: &BaroclinicArgs,
) -> Result<(), String> {
    let (images, colormap, (w, h)) = paint(world, key, overlay, frames, &args.image);
    let path = args.image.out_dir.join("dyn.surface_pressure.gif");
    write_gif(&path, w, h, &images, colormap.palette(), args.fps)
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
    }
    Ok(())
}

/// The blow-up dump: `report.txt`, `diagnostics.csv`, and maps of the last
/// healthy state.
#[allow(clippy::too_many_arguments)]
fn write_dump(
    dir: &Path,
    world: &World,
    fields: &Fields,
    si: &TwoLayerSemiImplicit,
    last_good: &TwoLayerState,
    log: &RunLog,
    blow_up: &BlowUp,
    setup: &str,
    image: &ImageArgs,
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
        "fms baroclinic blew up at step {} (day {day:.3})\n",
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
        "Recent diagnostics (residual = energy change not accounted for by dissipation):"
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
    publish(world, fields, last_good, si.model());
    write_maps(
        world,
        fields,
        dir,
        &format!("LAST HEALTHY DAY {:.2}", day - dt / DAY),
        image,
    )
}

#[cfg(test)]
mod tests {
    use clap::Parser;

    use super::*;

    #[derive(Parser)]
    struct Cli {
        #[command(flatten)]
        args: BaroclinicArgs,
    }

    #[test]
    fn a_blow_up_stops_the_run_and_writes_a_dump() {
        // Inviscid, with the step sized for 45 m/s winds: the storms outrun
        // its Courant limit (100 m/s at this resolution) within days.
        let dir = std::env::temp_dir().join(format!("fms-blowup-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let cli = Cli::parse_from([
            "fms",
            "--terrain-frequency",
            "16",
            "--climate-frequency",
            "16",
            "--shear",
            "40",
            "--max-wind",
            "45",
            "--inviscid",
            "--days",
            "30",
            "--frames",
            "4",
            "--width",
            "128",
            "--out-dir",
            dir.to_str().unwrap(),
        ]);
        let err = run(cli.args).unwrap_err();
        assert!(err.contains("blew up"), "{err}");
        let dump = dir.join("baroclinic/blowup");
        let report = std::fs::read_to_string(dump.join("report.txt")).unwrap();
        assert!(report.contains("Courant limit"), "{report}");
        assert!(report.contains("Recent diagnostics"), "{report}");
        for file in [
            "diagnostics.csv",
            "dyn.surface_pressure.png",
            "dyn.wind_upper.png",
            "dyn.vorticity_lower.png",
        ] {
            assert!(dump.join(file).exists(), "missing {file}");
        }
        assert!(dir.join("baroclinic/dyn.surface_pressure.gif").exists());
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
