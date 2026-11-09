//! `fms held-suarez`: the Held & Suarez (1994) dynamical-core benchmark on
//! the two-layer model. Spins up from rest, then averages the zonal-mean
//! climate; writes profiles, a summary, and the usual maps and GIF (of the
//! last days) to `output/held_suarez/`.

use std::fmt::Write as _;

use clap::Args;
use fms_atmos::{HeldSuarez, HeldSuarezParams, ZonalClimatology, ZonalMeans};

use crate::chart::{self, Series, TimeAxis};
use crate::layered_run::{
    DynamicsArgs, Fields, Output, Plan, build_model, describe, integrate, stepper,
};
use crate::render::write_png;
use crate::{ImageArgs, WorldArgs, build_world};

/// Design wind: Held–Suarez winds peak near 65 m s⁻¹.
const MAX_WIND: f64 = 80.0;

#[derive(Args)]
pub struct HeldSuarezArgs {
    #[command(flatten)]
    world: WorldArgs,
    /// Simulated days, including spin-up.
    #[arg(long, default_value_t = 1200.0)]
    days: f64,
    /// Spin-up days excluded from the statistics.
    #[arg(long, default_value_t = 200.0)]
    spinup: f64,
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

pub fn run(args: HeldSuarezArgs) -> Result<(), String> {
    if !(0.0..args.days).contains(&args.spinup) {
        return Err("--spinup must be shorter than --days".into());
    }
    let mut world = build_world(&args.world)?;
    let fields = Fields::register(&mut world)?;
    let world = world;
    let setup = build_model(
        &world,
        &args.dynamics,
        MAX_WIND,
        &[0.5],
        vec![0.0; world.climate.mesh.n_cells()],
    )?;
    let hs = HeldSuarez::new(&setup.model, HeldSuarezParams::default());
    let state = hs.initial_state(&setup.model, 0.1);
    let model = setup.model.with_forcing(Box::new(hs));
    let n_cells = world.climate.mesh.n_cells();
    let bands = args
        .bands
        .unwrap_or_else(|| 2 * ((n_cells as f64 / 10.0).sqrt().round() as usize).max(8));
    let mut climate = ZonalClimatology::new(&model, bands);
    let description = describe(
        &model,
        &state,
        setup.dt,
        args.days,
        &format!(
            "; Held–Suarez forcing, statistics after day {}",
            args.spinup
        ),
    );
    eprintln!("held-suarez: {description}");

    let per_day = setup.steps_per_day;
    let total = (args.days * per_day as f64).round() as u64;
    let spinup = (args.spinup * per_day as f64).round() as u64;
    let movie_from = total.saturating_sub((args.movie_days * per_day as f64).round() as u64);
    let out = Output {
        command: "held-suarez",
        dir: args.image.out_dir.join("held_suarez"),
        image: &args.image,
        fps: args.fps,
    };
    let plan = Plan {
        total,
        frames: u64::from(args.frames),
        frames_from: movie_from,
    };
    let mut si = stepper(model, state, setup.dt);
    let mut spun_up = None;
    let log = integrate(&world, &fields, &mut si, &plan, &out, &description, |si| {
        if si.steps() > spinup {
            climate.sample(si.model(), si.state());
        }
        if spun_up.is_none() && si.steps() >= spinup {
            spun_up = Some((si.model().diagnostics(si.state()), si.energy_budget()));
        }
    })?;
    let means = crate::layered_run::upper_and_lower(si.model(), climate.means());

    // Profiles: CSV and chart.
    let mut csv = String::from(
        "latitude_deg,wind_upper_m_s,wind_lower_m_s,temperature_upper_k,temperature_lower_k,surface_pressure_pa\n",
    );
    for i in 0..means.latitude.len() {
        let _ = writeln!(
            csv,
            "{},{},{},{},{},{}",
            means.latitude[i],
            means.wind[0][i],
            means.wind[1][i],
            means.temperature[0][i],
            means.temperature[1][i],
            means.surface_pressure[i]
        );
    }
    let path = out.dir.join("zonal_mean.csv");
    std::fs::write(&path, csv).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!("zonal means -> {}", path.display());
    let hpa: Vec<f64> = means.surface_pressure.iter().map(|p| p / 100.0).collect();
    let series = [
        Series {
            title: "UPPER-LAYER EASTWARD WIND (M/S)".into(),
            values: &means.wind[0],
        },
        Series {
            title: "LOWER-LAYER EASTWARD WIND (M/S)".into(),
            values: &means.wind[1],
        },
        Series {
            title: "UPPER-LAYER TEMPERATURE (K)".into(),
            values: &means.temperature[0],
        },
        Series {
            title: "LOWER-LAYER TEMPERATURE (K)".into(),
            values: &means.temperature[1],
        },
        Series {
            title: "SURFACE PRESSURE (HPA)".into(),
            values: &hpa,
        },
    ];
    let header = format!(
        "HELD-SUAREZ ZONAL MEANS, DAYS {:.0} TO {:.0}, {} CELLS",
        args.spinup, args.days, n_cells
    );
    let axis = TimeAxis {
        times: &means.latitude,
        caption: "LATITUDE (DEGREES)",
        hours: false,
    };
    let width = args.image.width;
    let (pixels, height) = chart::render_chart(width, &header, &axis, &series);
    let path = out.dir.join("zonal_mean.png");
    write_png(&path, width, height, &pixels, &chart::palette())
        .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!("zonal-mean chart -> {}", path.display());

    // Summary.
    let mut summary = format!("fms held-suarez: {description}\n\n");
    summary.push_str(&jets(&means));
    if let (Some((d0, b0)), Some(end)) = (spun_up, log.last()) {
        let seconds = (end.step - spinup) as f64 * si.dt();
        let ops = si.model().operators();
        let area = 4.0 * std::f64::consts::PI * ops.radius() * ops.radius();
        let wm2 = |j: f64| j / seconds / area;
        let (d1, b1) = (end.diagnostics, end.budget);
        let change = d1.energy - d0.energy;
        let heating = b1.forcing.thermal - b0.forcing.thermal;
        let friction = b1.forcing.kinetic - b0.forcing.kinetic;
        let diffusion = b1.diffusion.total() - b0.diffusion.total();
        let _ = writeln!(
            summary,
            "\nEnergy budget after spin-up (W/m²): dE/dt {:+.3} = Newtonian heating {:+.3} \
             + friction {:+.3} + diffusion {:+.3} + residual {:+.4}",
            wm2(change),
            wm2(heating),
            wm2(friction),
            wm2(diffusion),
            wm2(change - heating - friction - diffusion)
        );
    }
    let _ = writeln!(
        summary,
        "Time means of {} samples ({} per day).",
        climate.samples(),
        crate::layered_run::SAMPLES_PER_DAY
    );
    let path = out.dir.join("summary.txt");
    std::fs::write(&path, &summary).map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    print!("{summary}");
    Ok(())
}

/// The jets and their hemispheric agreement, as text.
fn jets(m: &ZonalMeans) -> String {
    let mut out = String::new();
    let peak = |v: &[f64], north: bool| {
        (0..v.len())
            .filter(|&i| (m.latitude[i] > 0.0) == north && m.latitude[i].abs() < 75.0)
            .fold((0, f64::NEG_INFINITY), |a, i| {
                if v[i] > a.1 { (i, v[i]) } else { a }
            })
    };
    for (name, north) in [("North", true), ("South", false)] {
        let (i, up) = peak(&m.wind[0], north);
        let (j, low) = peak(&m.wind[1], north);
        let _ = writeln!(
            out,
            "{name}: upper-layer jet {up:.1} m/s at {:.1}°; lower-layer westerlies {low:.1} m/s at {:.1}°",
            m.latitude[i], m.latitude[j]
        );
    }
    let tropics: Vec<f64> = (0..m.latitude.len())
        .filter(|&i| m.latitude[i].abs() <= 30.0)
        .map(|i| m.wind[1][i])
        .collect();
    let mean = |v: &[f64], lo: f64, hi: f64| {
        let x: Vec<f64> = (0..m.latitude.len())
            .filter(|&i| (lo..=hi).contains(&m.latitude[i].abs()))
            .map(|i| v[i])
            .collect();
        x.iter().sum::<f64>() / x.len().max(1) as f64
    };
    let _ = writeln!(
        out,
        "Lower-layer wind within 30° of the equator: {:.2} m/s; lower-layer temperature equator − pole: {:.1} K",
        tropics.iter().sum::<f64>() / tropics.len().max(1) as f64,
        mean(&m.temperature[1], 0.0, 10.0) - mean(&m.temperature[1], 80.0, 90.0)
    );
    let mirror = m.mirrored();
    let worst = |a: &[f64], b: &[f64]| {
        a.iter()
            .zip(b)
            .fold(0.0f64, |w, (x, y)| w.max((x - y).abs()))
    };
    let _ = writeln!(
        out,
        "North vs south, largest difference: upper wind {:.2} m/s, lower wind {:.2} m/s, lower temperature {:.2} K",
        worst(&m.wind[0], &mirror.wind[0]),
        worst(&m.wind[1], &mirror.wind[1]),
        worst(&m.temperature[1], &mirror.temperature[1])
    );
    out
}
