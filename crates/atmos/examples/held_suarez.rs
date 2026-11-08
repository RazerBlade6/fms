//! The Held–Suarez benchmark on the two-layer core: spin up, then print the
//! time-mean zonal-mean climate, hemispheric asymmetry and energy budget.
//!
//! ```sh
//! cargo run --release -p fms-atmos --example held_suarez [frequency] [days] [spin-up days]
//! ```
//!
//! Environment: `FMS_MAX_WIND` (design wind for the step, default 80 m/s),
//! `FMS_BANDS` (latitude bands, default 36).

use std::time::Instant;

use fms_atmos::williamson::{GRAVITY, OMEGA, RADIUS};
use fms_atmos::{
    Gas, HeldSuarez, HeldSuarezParams, HyperdiffusionOptions, TwoLayer, TwoLayerSemiImplicit,
    TwoLayerSemiImplicitOptions, ZonalClimatology,
};
use fms_geo::{LloydOptions, MeshGeometry, Operators, SphereMesh};

const DAY: f64 = 86_400.0;

fn env<T: std::str::FromStr>(name: &str, default: T) -> T {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

fn main() {
    let mut args = std::env::args().skip(1);
    let k: u32 = args.next().map_or(16, |s| s.parse().unwrap());
    let days: f64 = args.next().map_or(500.0, |s| s.parse().unwrap());
    let spinup: f64 = args.next().map_or(200.0, |s| s.parse().unwrap());
    let (mesh, _) = SphereMesh::centroidal(k, LloydOptions::default());
    let geometry = MeshGeometry::new(&mesh);
    let model = TwoLayer::new(
        Operators::new(&mesh, &geometry, RADIUS),
        GRAVITY,
        Gas::EARTH_AIR,
        0.5,
        mesh.vertices().iter().map(|p| 2.0 * OMEGA * p.z).collect(),
        vec![0.0; mesh.n_cells()],
    );
    let dt = TwoLayerSemiImplicit::stable_dt(&model, env("FMS_MAX_WIND", 80.0), 0.6);
    let dt = DAY / (DAY / dt).ceil();
    let diffusion = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
    let hs = HeldSuarez::new(&model, HeldSuarezParams::default());
    println!(
        "frequency {k}: {} cells, dt {dt:.0} s, ∇⁴ τ = {:.1} h; friction {:.3}/day (lower)",
        mesh.n_cells(),
        diffusion.e_folding / 3600.0,
        hs.friction()[1] * DAY
    );
    let state = hs.initial_state(&model, 0.1);
    let model = model
        .with_hyperdiffusion(diffusion)
        .with_forcing(Box::new(hs));
    let start = model.diagnostics(&state);
    let mut climate = ZonalClimatology::new(&model, env("FMS_BANDS", 36));
    let mut si = TwoLayerSemiImplicit::new(
        model,
        state,
        dt,
        TwoLayerSemiImplicitOptions::default(),
        0.5,
    );
    let per_day = (DAY / dt).round() as u64;
    let t = Instant::now();
    let mut peak: f64 = 0.0;
    let mut at_spinup = None;
    while si.steps() < (days * per_day as f64) as u64 {
        si.step();
        let h = match si.check() {
            Ok(h) => h,
            Err(b) => {
                println!("BLOW-UP day {:.2}: {b}", si.steps() as f64 / per_day as f64);
                return;
            }
        };
        peak = peak.max(h.max_wind);
        let day = si.steps() as f64 / per_day as f64;
        if si.steps() % (per_day / 4).max(1) == 0 && day > spinup {
            climate.sample(si.model(), si.state());
        }
        if si.steps() == (spinup * per_day as f64) as u64 {
            at_spinup = Some((si.model().diagnostics(si.state()), si.energy_budget()));
        }
        if si.steps() % (50 * per_day) == 0 {
            let d = si.model().diagnostics(si.state());
            let b = si.energy_budget();
            println!(
                "day {day:>5.0}: max |u| {:>5.1}  KE {:.3e}  E−E₀ {:+.2e}  friction {:+.2e}  newtonian {:+.2e}  diffusion {:+.2e}  residual {:+.2e}",
                h.max_wind,
                d.kinetic_energy,
                d.energy - start.energy,
                b.forcing.kinetic,
                b.forcing.thermal,
                b.diffusion.total(),
                d.energy - start.energy - b.forcing.total() - b.diffusion.total(),
            );
        }
    }
    let elapsed = t.elapsed().as_secs_f64();
    println!(
        "{days} days in {elapsed:.0} s ({:.0} steps/s); peak |u| {peak:.1} m/s; {} samples",
        si.steps() as f64 / elapsed,
        climate.samples()
    );
    if let Some((d0, b0)) = at_spinup {
        let d1 = si.model().diagnostics(si.state());
        let b1 = si.energy_budget();
        let years = (days - spinup) * DAY;
        let area = 4.0 * std::f64::consts::PI * RADIUS * RADIUS;
        let wm2 = |j: f64| j / years / area;
        println!(
            "after spin-up (W/m²): dE/dt {:+.3}  friction {:+.3}  newtonian {:+.3}  diffusion kinetic {:+.3} thermal {:+.3}  residual {:+.4}",
            wm2(d1.energy - d0.energy),
            wm2(b1.forcing.kinetic - b0.forcing.kinetic),
            wm2(b1.forcing.thermal - b0.forcing.thermal),
            wm2(b1.diffusion.kinetic - b0.diffusion.kinetic),
            wm2(b1.diffusion.thermal - b0.diffusion.thermal),
            wm2(d1.energy
                - d0.energy
                - (b1.forcing.total() - b0.forcing.total())
                - (b1.diffusion.total() - b0.diffusion.total())),
        );
    }
    let m = climate.means();
    let mirror = m.mirrored();
    println!(" lat   u_up  u_low   T_up  T_low   p_s    | SH mirror u_up u_low T_low");
    for i in 0..m.latitude.len() {
        println!(
            "{:>5.1} {:>6.1} {:>6.1} {:>6.1} {:>6.1} {:>7.1} | {:>6.1} {:>6.1} {:>6.1}",
            m.latitude[i],
            m.wind[0][i],
            m.wind[1][i],
            m.temperature[0][i],
            m.temperature[1][i],
            m.surface_pressure[i] / 100.0,
            mirror.wind[0][i],
            mirror.wind[1][i],
            mirror.temperature[1][i],
        );
    }
}
