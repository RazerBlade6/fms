//! `fms baroclinic`: run the two-layer dry core on a balanced, vertically
//! sheared zonal flow seeded with small waves, and render the baroclinic
//! storms that grow from it.
//!
//! The state is checked after every step ([`fms_atmos::safety`]). If it blows
//! up, the run stops and writes a diagnostic dump to `baroclinic/blowup/`: a
//! report, the diagnostics history, and maps of the last healthy state.

use clap::Args;
use fms_atmos::ShearedFlow;
use fms_geo::coords::to_lat_lon;

use crate::layered_run::{
    DynamicsArgs, Fields, Output, Plan, build_model, describe, integrate, stepper,
};
use crate::{ImageArgs, WorldArgs, build_world};

/// Design wind: storms on the default jet peak near 100 m s⁻¹.
const MAX_WIND: f64 = 120.0;

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
    #[command(flatten)]
    dynamics: DynamicsArgs,
    #[command(flatten)]
    image: ImageArgs,
}

pub fn run(args: BaroclinicArgs) -> Result<(), String> {
    let mut world = build_world(&args.world)?;
    let fields = Fields::register(&mut world)?;
    let world = world;
    let setup = build_model(&world, &args.dynamics, MAX_WIND)?;
    if setup.model.layers() != 2 {
        return Err("fms baroclinic is a two-layer experiment (Phillips's model): leave --levels at one interface".into());
    }
    let model = setup.model;
    let flow = ShearedFlow {
        speed: [10.0 + args.shear, 10.0],
        stability: args.stability,
        theta_equator: 315.0,
        ps_equator: world.planet.params().surface_pressure,
        omega: world.planet.rotation_rate(),
    };
    let mut state = flow.state(&model);
    // Seed waves 4–8 with fixed phases, in both hemispheres.
    for k in 0..2 {
        for (i, &p) in world.climate.mesh.cell_centers().iter().enumerate() {
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
    let description = describe(&model, &state, setup.dt, args.days, "");
    eprintln!("baroclinic: {description}");
    let total = (args.days * setup.steps_per_day as f64).round() as u64;
    let mut si = stepper(model, state, setup.dt);
    let out = Output {
        command: "baroclinic",
        dir: args.image.out_dir.join("baroclinic"),
        image: &args.image,
        fps: args.fps,
    };
    let plan = Plan {
        total,
        frames: u64::from(args.frames),
        frames_from: 0,
    };
    integrate(&world, &fields, &mut si, &plan, &out, &description, |_| {})?;
    Ok(())
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
