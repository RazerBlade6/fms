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
//! ```

mod colormap;
mod font;
mod render;

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use clap::{Args, Parser, Subcommand};
use fms_fields::{ColorScale, Element, FieldKey};
use fms_sched::Hemisphere;
use fms_sim::{World, WorldConfig};

use crate::colormap::Colormap;
use crate::render::{CellRaster, Legend, format_value, image_height, render, write_gif, write_png};

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
    /// Render every cell field.
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
            "{:<28} {:?} {:?} {:<8} {:>9} values  {}",
            s.name,
            s.location.grid,
            s.location.element,
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
    let date = calendar.date(world.scheduler.current_step().saturating_sub(1));
    let hour = date.day_fraction(&calendar) * 24.0;
    format!(
        "YEAR {} DAY {} {:02}:{:02} {:?}",
        date.year,
        date.day_of_year + 1,
        hour.floor() as u32,
        ((hour.fract()) * 60.0).round() as u32 % 60,
        date.season(&calendar, Hemisphere::North),
    )
}

fn resolve(world: &World, name: &str) -> Result<FieldKey, String> {
    let key = world.fields.key(name).map_err(|e| e.to_string())?;
    let element = world.fields.meta(key).spec.location.element;
    if element != Element::Cell {
        return Err(format!(
            "`{name}` lives on {element:?}s; only cell fields can be rendered so far"
        ));
    }
    Ok(key)
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

fn title(world: &World, key: FieldKey) -> String {
    let spec = &world.fields.meta(key).spec;
    if spec.units.is_empty() {
        spec.name.to_owned()
    } else {
        format!("{} ({})", spec.name, spec.units)
    }
}

fn dump(args: DumpArgs) -> Result<(), String> {
    let mut world = build_world(&args.world);
    let names: Vec<String> = if args.all {
        world
            .fields
            .keys()
            .filter(|&k| world.fields.meta(k).spec.location.element == Element::Cell)
            .map(|k| world.fields.meta(k).spec.name.to_owned())
            .collect()
    } else if args.fields.is_empty() {
        return Err("give at least one --field NAME, or --all (see `fms fields`)".into());
    } else {
        args.fields.clone()
    };
    let keys = names
        .iter()
        .map(|n| resolve(&world, n))
        .collect::<Result<Vec<_>, _>>()?;

    // Always run at least one step so every system has produced its output.
    let steps = steps_for(&world, args.days + args.hours / 24.0);
    world.run_steps(steps.max(1));
    ensure_dir(&args.image.out_dir)?;
    let stamp = date_label(&world);

    for key in keys {
        let meta = world.fields.meta(key);
        let grid = world.grid(meta.spec.location.grid);
        let smooth = args.image.smooth && meta.spec.color_scale != ColorScale::Categorical;
        let raster = CellRaster::new(grid, args.image.width, smooth);
        let values = world.fields.read_column(key).to_f64_vec();
        let (min, max) = args.image.range.unwrap_or_else(|| data_range(&values));
        let colormap = Colormap::new(meta.spec.color_scale, min, max);
        let pixels = render(
            &raster,
            &values,
            &colormap,
            &Legend {
                title: &title(&world, key),
                subtitle: &stamp,
            },
        );
        let path = args.image.out_dir.join(format!("{}.png", meta.spec.name));
        write_png(
            &path,
            raster.width,
            image_height(&raster),
            &pixels,
            colormap.palette(),
        )
        .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
        println!(
            "{} [{} .. {}] -> {}",
            meta.spec.name,
            format_value(min),
            format_value(max),
            path.display()
        );
    }
    Ok(())
}

fn animate(args: AnimateArgs) -> Result<(), String> {
    let mut world = build_world(&args.world);
    let key = resolve(&world, &args.field)?;
    let days = args.days.unwrap_or(world.planet.params().year_days);
    let frames = args.frames.max(1);
    let total_steps = steps_for(&world, days);
    let frame_step = |k: u32| total_steps * u64::from(k) / u64::from(frames);

    // Simulate first, then pick one colour range for every frame.
    let t = Instant::now();
    let mut snapshots = Vec::with_capacity(frames as usize);
    world.run_steps(1);
    for k in 0..frames {
        let target = frame_step(k).max(1);
        let now = world.scheduler.current_step();
        world.run_steps(target.saturating_sub(now));
        snapshots.push((
            world.fields.read_column(key).to_f64_vec(),
            date_label(&world),
        ));
    }
    eprintln!(
        "simulated {days} days ({total_steps} steps) in {:.1?}",
        t.elapsed()
    );

    let (min, max) = args.image.range.unwrap_or_else(|| {
        let all: Vec<f64> = snapshots
            .iter()
            .flat_map(|(v, _)| v.iter().copied())
            .collect();
        data_range(&all)
    });
    let meta = world.fields.meta(key);
    let colormap = Colormap::new(meta.spec.color_scale, min, max);
    let smooth = args.image.smooth && meta.spec.color_scale != ColorScale::Categorical;
    let raster = CellRaster::new(
        world.grid(meta.spec.location.grid),
        args.image.width,
        smooth,
    );
    let title = title(&world, key);
    let images: Vec<Vec<u8>> = snapshots
        .iter()
        .map(|(values, stamp)| {
            render(
                &raster,
                values,
                &colormap,
                &Legend {
                    title: &title,
                    subtitle: stamp,
                },
            )
        })
        .collect();

    ensure_dir(&args.image.out_dir)?;
    let path = args.image.out_dir.join(format!("{}.gif", meta.spec.name));
    write_gif(
        &path,
        raster.width,
        image_height(&raster),
        &images,
        colormap.palette(),
        args.fps,
    )
    .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!(
        "{} [{} .. {}], {} frames -> {}",
        meta.spec.name,
        format_value(min),
        format_value(max),
        frames,
        path.display()
    );
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
