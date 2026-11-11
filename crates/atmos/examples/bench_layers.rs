//! Wall-clock cost of the 3-D model by number of layers: the moist
//! aquaplanet (radiation, moisture, slab) or dry Held–Suarez, run for a few
//! days and repeated, reporting each run and the median. Run one at a time on
//! an idle machine and log the conditions with the numbers in
//! `docs/BENCHMARKS.md`.
//!
//! ```sh
//! cargo run --release -p fms-atmos --example bench_layers -- \
//!     [--mode moist|hs] [--levels 0.5] [--days 30] [--repeats 3] [--k 16] [--reference means]
//! ```
//!
//! `--reference means` uses the start's global means as the semi-implicit
//! reference (step 9's choice) instead of the isothermal default.

use std::time::Instant;

use fms_atmos::column_physics::initial_state;
use fms_atmos::layered_si::global_means;
use fms_atmos::williamson::{GRAVITY, RADIUS};
use fms_atmos::{
    ColumnPhysics, Gas, HeldSuarez, HeldSuarezParams, HyperdiffusionOptions, LayeredModel,
    LayeredSemiImplicit, LayeredSemiImplicitOptions, Levels, PhysicsParams,
};
use fms_grid::{LloydOptions, MeshGeometry, Operators, SphereMesh};
use fms_planet::Planet;

const DAY: f64 = 86_400.0;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let arg = |name: &str, default: &str| -> String {
        args.iter()
            .position(|a| a == name)
            .map_or(default.to_string(), |i| args[i + 1].clone())
    };
    let mode = arg("--mode", "moist");
    let interior: Vec<f64> = arg("--levels", "0.5")
        .split(',')
        .map(|s| s.parse().expect("--levels: comma-separated σ"))
        .collect();
    let days: f64 = arg("--days", "30").parse().unwrap();
    let repeats: usize = arg("--repeats", "3").parse().unwrap();
    let k: u32 = arg("--k", "16").parse().unwrap();
    let means = arg("--reference", "isothermal") == "means";

    let (mesh, _) = SphereMesh::centroidal(k, LloydOptions::default());
    let geometry = MeshGeometry::new(&mesh);
    let planet = Planet::earth();
    let levels = Levels::new(&interior);
    println!(
        "{mode}, {} cells, {} layers {:?}, {days} days, reference {}, {} rayon threads",
        mesh.n_cells(),
        levels.count(),
        interior,
        if means { "global means" } else { "isothermal" },
        rayon::current_num_threads()
    );
    let mut runs = Vec::new();
    for run in 0..repeats {
        let model = LayeredModel::new(
            Operators::new(&mesh, &geometry, RADIUS),
            GRAVITY,
            Gas::EARTH_AIR,
            levels.clone(),
            mesh.vertices()
                .iter()
                .map(|p| 2.0 * planet.rotation_rate() * p.z)
                .collect(),
            vec![0.0; mesh.n_cells()],
        );
        let max_wind = if mode == "hs" { 80.0 } else { 120.0 };
        let dt = LayeredSemiImplicit::stable_dt(&model, max_wind, 0.6);
        let dt = DAY / (DAY / dt).ceil();
        let diffusion = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
        let setup = Instant::now();
        let (model, state) = if mode == "hs" {
            let hs = HeldSuarez::new(&model, HeldSuarezParams::default());
            let state = hs.initial_state(&model, 0.1);
            (
                model
                    .with_hyperdiffusion(diffusion)
                    .with_forcing(Box::new(hs)),
                state,
            )
        } else {
            let params = PhysicsParams::default();
            let model = model.with_humidity_profile(params.radiation);
            let (state, ts) = initial_state(&model, &planet, &params);
            let physics = ColumnPhysics::new(
                &model,
                &planet,
                params.with_step(dt, planet.solar_day()),
                0.0,
                &state,
                vec![ts; mesh.n_cells()],
            );
            (
                model
                    .with_hyperdiffusion(diffusion)
                    .with_forcing(Box::new(physics)),
                state,
            )
        };
        let options = LayeredSemiImplicitOptions {
            reference: means.then(|| global_means(model.operators(), &state)),
            ..LayeredSemiImplicitOptions::default()
        };
        let mut si = LayeredSemiImplicit::new(model, state, dt, options, 0.5);
        let setup = setup.elapsed().as_secs_f64();
        let steps = (days * DAY / dt).round() as u64;
        let start = Instant::now();
        while si.steps() < steps {
            si.step();
        }
        let wall = start.elapsed().as_secs_f64();
        si.check().expect("the benchmark run blew up");
        let stats = si.solver_stats();
        println!(
            "  run {}: setup {setup:.2} s, {steps} steps of {dt:.0} s in {wall:.2} s = {:.2} ms/step, {:.1} model days/s, CG {:.1}/solve",
            run + 1,
            1e3 * wall / steps as f64,
            days / wall,
            stats.iterations as f64 / stats.solves as f64
        );
        runs.push(wall / steps as f64);
    }
    runs.sort_by(f64::total_cmp);
    println!("  median {:.2} ms/step", 1e3 * runs[runs.len() / 2]);
}
