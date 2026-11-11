//! `fms-climate`: headless runner for the climate simulator.
//!
//! Builds a world, runs it, and renders fields to equirectangular PNGs and GIFs
//! (default output directory: `output/`).
//!
//! ```sh
//! fms-climate fields                                   # list fields
//! fms-climate dump --all                               # PNG of every cell field at t = 0
//! fms-climate dump -f planet.toa_insolation --hours 12 # one field after half a day
//! fms-climate animate                                  # a year of daily-mean insolation
//! fms-climate animate -f planet.toa_insolation --days 1 --frames 48
//! fms-climate patterns                                 # analytic test fields (edges, vertices, vectors)
//! fms-climate probe --lat 0 --lon 0 -f planet.toa_insolation --days 2  # time series at a point
//! fms-climate williamson --case 5                     # shallow-water flow over a mountain
//! fms-climate baroclinic                              # two-layer core: storms on a sheared jet
//! fms-climate held-suarez --climate-frequency 32      # dynamical-core benchmark, 1200 days
//! fms-climate climate --climate-frequency 16 --co2 840  # radiative climate, final-year means
//! fms-climate basins --heightmap earth.tif             # sea mask and ocean basins
//! ```

mod arrows;
mod baroclinic;
mod basins;
mod chart;
mod climate;
mod colormap;
mod font;
mod held_suarez;
mod layered_run;
mod ocean_transport;
mod patterns;
mod render;
mod run_log;
mod view;
mod williamson;

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use clap::{Args, Parser, Subcommand};
use fms_fields::{ColorScale, FieldKey};
use fms_planet::Hemisphere;
use fms_sim::{Probe, World, WorldConfig};
use glam::DVec3;

use crate::arrows::ArrowLayout;
use crate::colormap::Colormap;
use crate::render::{
    CellRaster, Legend, Overlay, format_value, image_height, render, write_gif, write_png,
};
use crate::view::{View, is_vector, title, view};

#[derive(Parser)]
#[command(
    name = "fms-climate",
    version,
    about = "fms-climate: headless climate simulator"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// List every field with its location, units and description.
    Fields(WorldArgs),
    /// Run the world for a while, then write a PNG of each requested field.
    Dump(DumpArgs),
    /// Run the world and write a GIF of one field over time.
    Animate(AnimateArgs),
    /// Write analytic test fields (edge, vertex and vector fields) to
    /// `output/patterns/`, to check rendering before real dynamics exist.
    Patterns(PatternArgs),
    /// Record fields at one location over time: a CSV plus a line chart.
    Probe(ProbeArgs),
    /// Run a Williamson shallow-water test case (2, 5 or 6) on the climate
    /// grid; write a GIF and final-state PNGs to `output/williamsonN/`.
    Williamson(williamson::WilliamsonArgs),
    /// Run the two-layer dry core on a sheared jet seeded with small waves
    /// and render the growing baroclinic storms to `output/baroclinic/`.
    Baroclinic(baroclinic::BaroclinicArgs),
    /// Run the Held–Suarez dynamical-core benchmark: spin up from rest under
    /// idealized heating and friction, then average the zonal-mean climate;
    /// write profiles, a summary, maps and a GIF to `output/held_suarez/`.
    HeldSuarez(held_suarez::HeldSuarezArgs),
    /// Run the two-layer model driven by radiation over a slab aquaplanet
    /// on the world's planet; average the final year and write zonal means,
    /// a summary with the energy budget, and annual-mean maps to
    /// `output/climate/`.
    Climate(climate::ClimateArgs),
    /// Write the world's sea mask, land fraction and ocean basins (the
    /// gyres' domains) as maps, with tables of the bodies below sea level
    /// and of the basins, without running anything.
    Basins(basins::BasinsArgs),
}

#[derive(Args, Clone)]
struct WorldArgs {
    /// Terrain grid frequency (10·f² + 2 cells; 160 ≈ 256k cells, 400 ≈ 1.6M).
    #[arg(long, default_value_t = 160)]
    terrain_frequency: u32,
    /// Climate grid frequency (45 ≈ 20k cells).
    #[arg(long, default_value_t = 45)]
    climate_frequency: u32,
    /// Base steps per solar day.
    #[arg(long, default_value_t = 48)]
    steps_per_day: u32,
    /// Global heightmap to use as the terrain: a GeoTIFF such as ETOPO 2022
    /// (`scripts/fetch-etopo.sh`), or a greyscale PNG with --heightmap-range.
    /// Without one the world is a blank ocean.
    #[arg(long)]
    heightmap: Option<PathBuf>,
    /// Elevations (m) of black and white in a PNG heightmap, `MIN,MAX`.
    #[arg(long, value_parser = parse_range, allow_hyphen_values = true)]
    heightmap_range: Option<(f64, f64)>,
    /// Use an idealized continent with a mountain ridge as the terrain.
    #[arg(long, conflicts_with = "heightmap")]
    continent: bool,
    /// Atmospheric CO₂ (ppm, 0–100,000). Default: Earth's 420.
    #[arg(long)]
    co2: Option<f64>,
    /// Surface pressure (bar, 0.1–10). Default: Earth's 1.013.
    #[arg(long)]
    surface_pressure: Option<f64>,
    /// Axial tilt (degrees). Default: Earth's 23.44.
    #[arg(long)]
    obliquity: Option<f64>,
    /// Orbital eccentricity. Default: Earth's 0.0167.
    #[arg(long)]
    eccentricity: Option<f64>,
    /// Solar day: noon to noon (hours, 4–240). Default: Earth's 24.
    #[arg(long)]
    solar_day: Option<f64>,
    /// Year (solar days, 20–5,000). Default: Earth's 365.2422.
    #[arg(long)]
    year_days: Option<f64>,
}

#[derive(Args)]
struct DumpArgs {
    #[command(flatten)]
    world: WorldArgs,
    /// Field to render (repeatable).
    #[arg(short, long = "field")]
    fields: Vec<String>,
    /// Render every field.
    #[arg(long)]
    all: bool,
    /// Simulated days to run before rendering.
    #[arg(long, default_value_t = 0.0)]
    days: f64,
    /// Extra simulated hours to run before rendering.
    #[arg(long, default_value_t = 0.0)]
    hours: f64,
    #[command(flatten)]
    image: ImageArgs,
}

#[derive(Args)]
struct AnimateArgs {
    #[command(flatten)]
    world: WorldArgs,
    /// Field to animate.
    #[arg(short, long, default_value = "planet.daily_insolation")]
    field: String,
    /// Simulated days the animation spans (default: one planet year).
    #[arg(long)]
    days: Option<f64>,
    /// Number of frames.
    #[arg(long, default_value_t = 73)]
    frames: u32,
    /// Frames per second.
    #[arg(long, default_value_t = 12)]
    fps: u16,
    #[command(flatten)]
    image: ImageArgs,
}

#[derive(Args)]
struct ProbeArgs {
    #[command(flatten)]
    world: WorldArgs,
    /// Latitude in degrees (-90 to 90; south is negative).
    #[arg(long, allow_hyphen_values = true)]
    lat: f64,
    /// Longitude in degrees (west is negative).
    #[arg(long, allow_hyphen_values = true)]
    lon: f64,
    /// Field to record (repeatable). Vector fields give east and north columns.
    #[arg(short, long = "field", required = true)]
    fields: Vec<String>,
    /// Simulated days to record.
    #[arg(long, default_value_t = 1.0)]
    days: f64,
    /// Extra simulated hours to record.
    #[arg(long, default_value_t = 0.0)]
    hours: f64,
    /// Record every N base steps, starting with the first.
    #[arg(long, default_value_t = 1, value_parser = clap::value_parser!(u64).range(1..))]
    every: u64,
    /// Chart width in pixels.
    #[arg(long, default_value_t = 1024)]
    width: usize,
    /// Output directory.
    #[arg(long, default_value = "output")]
    out_dir: PathBuf,
}

#[derive(Args)]
struct PatternArgs {
    #[command(flatten)]
    world: WorldArgs,
    #[command(flatten)]
    image: ImageArgs,
}

#[derive(Args)]
struct ImageArgs {
    /// Image width in pixels (height is half).
    #[arg(long, default_value_t = 1024)]
    width: usize,
    /// Fixed colour range `MIN,MAX` (default: the data range).
    #[arg(long, value_parser = parse_range, allow_hyphen_values = true)]
    range: Option<(f64, f64)>,
    /// Output directory.
    #[arg(long, default_value = "output")]
    out_dir: PathBuf,
    /// Blend between cell centres instead of drawing each cell flat. Prettier,
    /// but hides the true resolution. Ignored for categorical fields.
    #[arg(long)]
    smooth: bool,
    /// Overlay arrows from this vector field (an edge normal-component field,
    /// e.g. a wind). Vector fields show their own arrows by default.
    #[arg(long)]
    vectors: Option<String>,
    /// Draw no arrows, even for vector fields.
    #[arg(long, conflicts_with = "vectors")]
    no_arrows: bool,
    /// Approximate spacing between arrows (pixels).
    #[arg(long, default_value_t = 28)]
    arrow_spacing: usize,
}

fn parse_range(s: &str) -> Result<(f64, f64), String> {
    let (a, b) = s.split_once(',').ok_or("expected MIN,MAX")?;
    let parse = |x: &str| x.trim().parse::<f64>().map_err(|e| e.to_string());
    let (min, max) = (parse(a)?, parse(b)?);
    if max > min {
        Ok((min, max))
    } else {
        Err("MAX must be greater than MIN".into())
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Command::Fields(args) => build_world(&args).and_then(|w| list_fields(&w)),
        Command::Dump(args) => dump(args),
        Command::Animate(args) => animate(args),
        Command::Patterns(args) => render_patterns(args),
        Command::Probe(args) => probe(args),
        Command::Williamson(args) => williamson::run(args),
        Command::Baroclinic(args) => baroclinic::run(args),
        Command::HeldSuarez(args) => held_suarez::run(args),
        Command::Climate(args) => climate::run(args),
        Command::Basins(args) => basins::run(args),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("error: {message}");
            ExitCode::FAILURE
        }
    }
}

fn build_world(args: &WorldArgs) -> Result<World, String> {
    let t = Instant::now();
    let mut config = WorldConfig {
        terrain_frequency: args.terrain_frequency,
        climate_frequency: args.climate_frequency,
        steps_per_day: args.steps_per_day,
        ..WorldConfig::default()
    };
    if let Some(co2) = args.co2 {
        config.planet.co2_ppm = co2;
    }
    if let Some(bar) = args.surface_pressure {
        config.planet.surface_pressure = bar * 1.0e5;
    }
    if let Some(degrees) = args.obliquity {
        config.planet.obliquity_deg = degrees;
    }
    if let Some(eccentricity) = args.eccentricity {
        config.planet.eccentricity = eccentricity;
    }
    if let Some(hours) = args.solar_day {
        config.planet.solar_day = hours * 3600.0;
    }
    if let Some(days) = args.year_days {
        config.planet.year_days = days;
    }
    let mut world = World::new(config);
    eprintln!(
        "world: {} terrain cells, {} climate cells, dt = {} s (built in {:.1?}: terrain grid {:.1}s, \
         climate grid {:.1}s, mapping {:.2}s)",
        world.terrain.mesh.n_cells(),
        world.climate.mesh.n_cells(),
        world.scheduler.calendar().dt(),
        t.elapsed(),
        world.timings.terrain_grid,
        world.timings.climate_grid,
        world.timings.two_level,
    );
    for a in &world.adjustments {
        eprintln!("planet: {a}");
    }
    if args.continent {
        world.set_idealized_continent();
    }
    if let Some(path) = &args.heightmap {
        let import = world
            .load_heightmap(path, args.heightmap_range)
            .map_err(|e| e.to_string())?;
        eprintln!(
            "heightmap: {} read at {}×{} in {:.1}s, resampled in {:.2}s",
            path.display(),
            import.size.0,
            import.size.1,
            import.read_seconds,
            import.resample_seconds,
        );
    }
    Ok(world)
}

fn list_fields(world: &World) -> Result<(), String> {
    for key in world.fields.keys() {
        let meta = world.fields.meta(key);
        let s = &meta.spec;
        println!(
            "{:<30} {:?} {:<6} {:<8} {:>9} values  {}",
            s.name,
            s.location.grid,
            format!("{:?}", s.location.element),
            format!("[{}]", s.units),
            meta.len,
            s.description
        );
    }
    Ok(())
}

fn ensure_dir(dir: &Path) -> Result<(), String> {
    std::fs::create_dir_all(dir).map_err(|e| format!("cannot create {}: {e}", dir.display()))
}

fn steps_for(world: &World, days: f64) -> u64 {
    (days * f64::from(world.scheduler.calendar().steps_per_day()))
        .round()
        .max(0.0) as u64
}

/// A human-readable date stamp for legends: the time of the last completed
/// step, which is when the fields were last computed.
fn date_label(world: &World) -> String {
    date_at(world, world.scheduler.current_step().saturating_sub(1))
}

/// Date stamp of the start of `step`.
fn date_at(world: &World, step: u64) -> String {
    let calendar = *world.scheduler.calendar();
    let date = calendar.date(step);
    let hour = date.day_fraction(&calendar) * 24.0;
    let season = world
        .planet
        .season(step as f64 * calendar.dt(), Hemisphere::North);
    format!(
        "YEAR {} DAY {} {:02}:{:02} N {:?}",
        date.year,
        date.day_of_year + 1,
        hour.floor() as u32,
        ((hour.fract()) * 60.0).round() as u32 % 60,
        season,
    )
}

/// Print where simulation time went: throughput and a per-system breakdown.
fn print_profile(world: &World, wall: std::time::Duration) {
    let steps = world.scheduler.current_step();
    let calendar = world.scheduler.calendar();
    let days = steps as f64 / f64::from(calendar.steps_per_day());
    let secs = wall.as_secs_f64().max(1e-9);
    eprintln!(
        "simulated {days:.2} days ({steps} steps) in {wall:.2?}: {:.0} steps/s, {:.2} days/s",
        steps as f64 / secs,
        days / secs
    );
    let profile = world.scheduler.profile();
    let busy: f64 = profile.iter().map(|p| p.total.as_secs_f64()).sum();
    eprintln!(
        "  {:<28} {:>8} {:>11} {:>11} {:>6}",
        "system", "runs", "total", "mean", "share"
    );
    for p in profile {
        let total = p.total.as_secs_f64();
        let mean = if p.runs > 0 {
            total / p.runs as f64
        } else {
            0.0
        };
        eprintln!(
            "  {:<28} {:>8} {:>9.2}ms {:>9.1}µs {:>5.1}%",
            p.name,
            p.runs,
            total * 1e3,
            mean * 1e6,
            if busy > 0.0 {
                100.0 * total / busy
            } else {
                0.0
            }
        );
    }
    eprintln!(
        "  {:<28} {:>8} {:>9.2}ms  (systems in one stage run in parallel, so this can exceed wall time)",
        "all systems",
        "",
        busy * 1e3
    );
}

fn resolve(world: &World, name: &str) -> Result<FieldKey, String> {
    world.fields.key(name).map_err(|e| e.to_string())
}

/// The field whose arrows to overlay on `key`, if any.
fn overlay_for(
    world: &World,
    key: FieldKey,
    image: &ImageArgs,
) -> Result<Option<FieldKey>, String> {
    if image.no_arrows {
        return Ok(None);
    }
    match &image.vectors {
        Some(name) => {
            let k = resolve(world, name)?;
            if is_vector(world, k) {
                Ok(Some(k))
            } else {
                Err(format!(
                    "`{name}` is not a vector (edge normal-component) field"
                ))
            }
        }
        None => Ok(is_vector(world, key).then_some(key)),
    }
}

fn data_range(values: &[f64]) -> (f64, f64) {
    let (min, max) = values
        .iter()
        .filter(|v| v.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
            (lo.min(v), hi.max(v))
        });
    if !min.is_finite() {
        (0.0, 1.0)
    } else if max > min {
        (min, max)
    } else {
        (min - 0.5, max + 0.5)
    }
}

/// One moment of a field, resampled for drawing.
struct Frame {
    view: View,
    arrows: Option<Vec<DVec3>>,
    stamp: String,
}

fn capture(world: &World, key: FieldKey, overlay: Option<FieldKey>) -> Frame {
    let arrows = overlay.map(|k| view(world, k).vectors.expect("overlay is a vector field"));
    Frame {
        view: view(world, key),
        arrows,
        stamp: date_label(world),
    }
}

/// Render frames of one field with a shared colour range and arrow scale.
/// Returns the images, the colour map and the image size.
fn paint(
    world: &World,
    key: FieldKey,
    overlay: Option<FieldKey>,
    frames: &[Frame],
    image: &ImageArgs,
) -> (Vec<Vec<u8>>, Colormap, (usize, usize)) {
    let spec = &world.fields.meta(key).spec;
    let (min, max) = image.range.unwrap_or_else(|| {
        let all: Vec<f64> = frames
            .iter()
            .flat_map(|f| f.view.values.iter().copied())
            .collect();
        data_range(&all)
    });
    let colormap = Colormap::new(spec.color_scale, min, max);
    let smooth = image.smooth && spec.color_scale != ColorScale::Categorical;
    let raster = CellRaster::new(world.grid(spec.location.grid), image.width, smooth);
    let layout = overlay.map(|k| {
        let grid = world.grid(world.fields.meta(k).spec.location.grid);
        ArrowLayout::new(grid, image.width, image.arrow_spacing)
    });
    let max_speed = frames
        .iter()
        .flat_map(|f| f.arrows.iter().flatten())
        .map(|v| v.length())
        .fold(0.0, f64::max);
    let mut heading = title(world, key);
    if let Some(k) = overlay.filter(|&k| k != key) {
        heading.push_str(&format!(" + {}", world.fields.meta(k).spec.name));
    }
    if max_speed > 0.0 {
        heading.push_str(&format!(" ARROW MAX {}", format_value(max_speed)));
    }
    let images = frames
        .iter()
        .map(|f| {
            let overlay = match (&layout, &f.arrows) {
                (Some(layout), Some(vectors)) => Some(Overlay {
                    layout,
                    vectors,
                    max_speed,
                }),
                _ => None,
            };
            let legend = Legend {
                title: &heading,
                subtitle: &f.stamp,
            };
            render(
                &raster,
                &f.view.values,
                &colormap,
                &legend,
                overlay.as_ref(),
            )
        })
        .collect();
    (images, colormap, (raster.width, image_height(&raster)))
}

fn write_one(
    world: &World,
    key: FieldKey,
    overlay: Option<FieldKey>,
    image: &ImageArgs,
) -> Result<(), String> {
    let frame = capture(world, key, overlay);
    let (images, colormap, (w, h)) =
        paint(world, key, overlay, std::slice::from_ref(&frame), image);
    let name = world.fields.meta(key).spec.name;
    let path = image.out_dir.join(format!("{name}.png"));
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

fn dump(args: DumpArgs) -> Result<(), String> {
    let mut world = build_world(&args.world)?;
    let keys: Vec<FieldKey> = if args.all {
        world.fields.keys().collect()
    } else if args.fields.is_empty() {
        return Err("give at least one --field NAME, or --all (see `fms-climate fields`)".into());
    } else {
        args.fields
            .iter()
            .map(|n| resolve(&world, n))
            .collect::<Result<_, _>>()?
    };
    let overlays = keys
        .iter()
        .map(|&k| overlay_for(&world, k, &args.image))
        .collect::<Result<Vec<_>, _>>()?;

    // Always run at least one step so every system has produced its output.
    let steps = steps_for(&world, args.days + args.hours / 24.0);
    let t = Instant::now();
    world.run_steps(steps.max(1));
    print_profile(&world, t.elapsed());
    ensure_dir(&args.image.out_dir)?;
    for (key, overlay) in keys.into_iter().zip(overlays) {
        write_one(&world, key, overlay, &args.image)?;
    }
    Ok(())
}

fn animate(args: AnimateArgs) -> Result<(), String> {
    let mut world = build_world(&args.world)?;
    let key = resolve(&world, &args.field)?;
    let overlay = overlay_for(&world, key, &args.image)?;
    let days = args.days.unwrap_or(world.planet.params().year_days);
    let frame_count = args.frames.max(1);
    let total_steps = steps_for(&world, days);
    let frame_step = |k: u32| total_steps * u64::from(k) / u64::from(frame_count);

    // Simulate first, then pick one colour range and arrow scale for all frames.
    let t = Instant::now();
    let mut frames = Vec::with_capacity(frame_count as usize);
    world.run_steps(1);
    for k in 0..frame_count {
        let target = frame_step(k).max(1);
        let now = world.scheduler.current_step();
        world.run_steps(target.saturating_sub(now));
        frames.push(capture(&world, key, overlay));
    }
    print_profile(&world, t.elapsed());

    let (images, colormap, (w, h)) = paint(&world, key, overlay, &frames, &args.image);
    ensure_dir(&args.image.out_dir)?;
    let name = world.fields.meta(key).spec.name;
    let path = args.image.out_dir.join(format!("{name}.gif"));
    write_gif(&path, w, h, &images, colormap.palette(), args.fps)
        .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    let (min, max) = colormap.range();
    println!(
        "{name} [{} .. {}], {frame_count} frames -> {}",
        format_value(min),
        format_value(max),
        path.display()
    );
    Ok(())
}

fn render_patterns(mut args: PatternArgs) -> Result<(), String> {
    let mut world = build_world(&args.world)?;
    let fields = patterns::register(&mut world)?;
    world.run_steps(1);
    args.image.out_dir = args.image.out_dir.join("patterns");
    ensure_dir(&args.image.out_dir)?;
    for (key, default_overlay) in fields {
        let overlay = if args.image.no_arrows {
            None
        } else if args.image.vectors.is_some() {
            overlay_for(&world, key, &args.image)?
        } else {
            default_overlay
        };
        write_one(&world, key, overlay, &args.image)?;
    }
    Ok(())
}

/// Longitude wrapped to (-180, 180].
fn wrap_lon(lon: f64) -> f64 {
    let l = (lon + 180.0).rem_euclid(360.0) - 180.0;
    if l == -180.0 { 180.0 } else { l }
}

/// `33.9S 18.4E`-style label, joined with `sep`.
fn coord_label(lat: f64, lon: f64, sep: &str) -> String {
    let part = |v: f64, pos: char, neg: char| {
        let s = format!("{:.2}", v.abs());
        let s = s.trim_end_matches('0').trim_end_matches('.');
        let hemisphere = if v < 0.0 && s != "0" { neg } else { pos };
        format!("{s}{hemisphere}")
    };
    format!(
        "{}{sep}{}",
        part(lat, 'N', 'S'),
        part(wrap_lon(lon), 'E', 'W')
    )
}

/// Run the world, recording the probe's columns at the steps given by
/// `--every`, then write `probe_<location>.csv` and `.png`.
fn probe(args: ProbeArgs) -> Result<(), String> {
    let mut world = build_world(&args.world)?;
    let keys = args
        .fields
        .iter()
        .map(|n| resolve(&world, n))
        .collect::<Result<Vec<_>, _>>()?;
    let probe = Probe::new(&world, args.lat, args.lon, &keys).map_err(|e| e.to_string())?;
    let n_cols = probe.columns().len();

    // A sample after step `s` holds the fields computed for time `s·dt`.
    let total = steps_for(&world, args.days + args.hours / 24.0).max(1);
    let mut steps = Vec::new();
    let mut values = Vec::new(); // row-major: one row of `n_cols` per sample
    let t = Instant::now();
    for step in 0..total {
        world.step();
        if step % args.every == 0 {
            steps.push(step);
            probe.sample_into(&world, &mut values);
        }
    }
    print_profile(&world, t.elapsed());

    ensure_dir(&args.out_dir)?;
    let tag = coord_label(args.lat, args.lon, "_");
    let csv_path = args.out_dir.join(format!("probe_{tag}.csv"));
    write_probe_csv(&csv_path, &world, &probe, &steps, &values)
        .map_err(|e| format!("cannot write {}: {e}", csv_path.display()))?;

    let calendar = *world.scheduler.calendar();
    let hours = (total as f64) / f64::from(calendar.steps_per_day()) <= 3.0;
    let per_day = if hours { 24.0 } else { 1.0 };
    let times: Vec<f64> = steps
        .iter()
        .map(|&s| s as f64 * calendar.dt() / calendar.day_seconds() * per_day)
        .collect();
    let columns: Vec<Vec<f64>> = (0..n_cols)
        .map(|j| values.iter().skip(j).step_by(n_cols).copied().collect())
        .collect();
    let series: Vec<chart::Series<'_>> = probe
        .columns()
        .iter()
        .zip(&columns)
        .map(|(c, v)| chart::Series {
            title: format!(
                "{}{}  {:?} CELL {} AT {}",
                c.label,
                units_suffix(c.units),
                c.site.grid,
                c.site.cell,
                coord_label(c.site.lat_deg, c.site.lon_deg, " ")
            ),
            values: v,
        })
        .collect();
    let header = format!(
        "PROBE {}   {} TO {}",
        coord_label(args.lat, args.lon, " "),
        date_at(&world, steps[0]),
        date_at(&world, *steps.last().expect("at least one sample")),
    );
    let axis = chart::TimeAxis {
        times: &times,
        caption: if hours { "HOURS" } else { "DAYS" },
        hours,
    };
    let (pixels, height) = chart::render_chart(args.width, &header, &axis, &series);
    let png_path = args.out_dir.join(format!("probe_{tag}.png"));
    write_png(&png_path, args.width, height, &pixels, &chart::palette())
        .map_err(|e| format!("cannot write {}: {e}", png_path.display()))?;

    for (c, v) in probe.columns().iter().zip(&columns) {
        let finite: Vec<f64> = v.iter().copied().filter(|x| x.is_finite()).collect();
        let (min, max) = finite
            .iter()
            .fold((f64::NAN, f64::NAN), |(lo, hi), &x| (x.min(lo), x.max(hi)));
        let mean = finite.iter().sum::<f64>() / finite.len().max(1) as f64;
        println!(
            "{}{}: min {}, max {}, mean {}{}",
            c.label,
            units_suffix(c.units),
            format_value(min),
            format_value(max),
            format_value(mean),
            if finite.len() < v.len() {
                format!(" ({} non-finite)", v.len() - finite.len())
            } else {
                String::new()
            }
        );
    }
    println!(
        "{} samples -> {}, {}",
        steps.len(),
        csv_path.display(),
        png_path.display()
    );
    Ok(())
}

fn units_suffix(units: &str) -> String {
    if units.is_empty() {
        String::new()
    } else {
        format!(" ({units})")
    }
}

/// One row per sample: step, time, date, then one column per probe column.
/// Floats are written in Rust's shortest round-trip form, so values are exact.
fn write_probe_csv(
    path: &Path,
    world: &World,
    probe: &Probe,
    steps: &[u64],
    values: &[f64],
) -> std::io::Result<()> {
    use std::io::Write;
    let calendar = *world.scheduler.calendar();
    let mut out = std::io::BufWriter::new(std::fs::File::create(path)?);
    write!(out, "step,time_s,day,year,day_of_year,hour")?;
    for c in probe.columns() {
        if c.units.is_empty() {
            write!(out, ",{}", c.label)?;
        } else {
            write!(out, ",{} [{}]", c.label, c.units)?;
        }
    }
    writeln!(out)?;
    let n = probe.columns().len();
    for (&step, row) in steps.iter().zip(values.chunks(n)) {
        let time = step as f64 * calendar.dt();
        let date = calendar.date(step);
        write!(
            out,
            "{step},{time},{},{},{},{}",
            time / calendar.day_seconds(),
            date.year,
            date.day_of_year + 1,
            date.day_fraction(&calendar) * 24.0
        )?;
        for v in row {
            write!(out, ",{v}")?;
        }
        writeln!(out)?;
    }
    out.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coordinate_labels() {
        assert_eq!(coord_label(33.9, 18.4, " "), "33.9N 18.4E");
        assert_eq!(coord_label(-33.925, -70.0, "_"), "33.92S_70W");
        assert_eq!(coord_label(0.0, 180.0, " "), "0N 180E");
        assert_eq!(coord_label(-0.001, 190.0, " "), "0N 170W");
        assert_eq!(coord_label(90.0, -180.0, " "), "90N 180E");
        assert_eq!(wrap_lon(540.0), 180.0);
        assert_eq!(wrap_lon(-190.0), 170.0);
    }

    #[test]
    fn ranges_parse() {
        assert_eq!(parse_range("-10,25.5"), Ok((-10.0, 25.5)));
        assert_eq!(parse_range(" 0 , 1 "), Ok((0.0, 1.0)));
        assert!(parse_range("5,5").is_err());
        assert!(parse_range("1").is_err());
        assert!(parse_range("a,b").is_err());
    }

    #[test]
    fn data_range_handles_edge_cases() {
        assert_eq!(data_range(&[1.0, f64::NAN, 3.0]), (1.0, 3.0));
        assert_eq!(data_range(&[2.0, 2.0]), (1.5, 2.5));
        assert_eq!(data_range(&[f64::NAN]), (0.0, 1.0));
    }

    #[test]
    fn cli_definition_is_valid() {
        use clap::CommandFactory;
        Cli::command().debug_assert();
    }
}
