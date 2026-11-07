//! `fms`: headless runner for the Fantasy Map Simulator.
//!
//! Builds a world, runs it, and renders fields to equirectangular PNGs and GIFs
//! (default output directory: `output/`).
//!
//! ```sh
//! fms fields                                   # list fields
//! fms dump --all                               # PNG of every cell field at t = 0
//! fms dump -f planet.toa_insolation --hours 12 # one field after half a day
//! fms animate                                  # a year of daily-mean insolation
//! fms animate -f planet.toa_insolation --days 1 --frames 48
//! fms patterns                                 # analytic test fields (edges, vertices, vectors)
//! ```

mod arrows;
mod colormap;
mod font;
mod patterns;
mod render;
mod view;

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use clap::{Args, Parser, Subcommand};
use fms_fields::{ColorScale, FieldKey};
use fms_planet::Hemisphere;
use fms_sim::{World, WorldConfig};
use glam::DVec3;

use crate::arrows::ArrowLayout;
use crate::colormap::Colormap;
use crate::render::{
    CellRaster, Legend, Overlay, format_value, image_height, render, write_gif, write_png,
};
use crate::view::{View, is_vector, title, view};

#[derive(Parser)]
#[command(
    name = "fms",
    version,
    about = "Fantasy Map Simulator: headless runner"
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
        Command::Fields(args) => list_fields(&build_world(&args)),
        Command::Dump(args) => dump(args),
        Command::Animate(args) => animate(args),
        Command::Patterns(args) => render_patterns(args),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("error: {message}");
            ExitCode::FAILURE
        }
    }
}

fn build_world(args: &WorldArgs) -> World {
    let t = Instant::now();
    let world = World::new(WorldConfig {
        terrain_frequency: args.terrain_frequency,
        climate_frequency: args.climate_frequency,
        steps_per_day: args.steps_per_day,
        ..WorldConfig::default()
    });
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
    world
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
    let calendar = *world.scheduler.calendar();
    let step = world.scheduler.current_step().saturating_sub(1);
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
    let mut world = build_world(&args.world);
    let keys: Vec<FieldKey> = if args.all {
        world.fields.keys().collect()
    } else if args.fields.is_empty() {
        return Err("give at least one --field NAME, or --all (see `fms fields`)".into());
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
    let mut world = build_world(&args.world);
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
    let mut world = build_world(&args.world);
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

#[cfg(test)]
mod tests {
    use super::*;

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
