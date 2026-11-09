//! The Held–Suarez benchmark on the layered core (two layers unless
//! `FMS_LEVELS` gives other interior interfaces, comma-separated): spin up,
//! then print the time-mean zonal-mean climate, hemispheric asymmetry, energy
//! budget and eddy momentum fluxes.
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
    Gas, HeldSuarez, HeldSuarezParams, HyperdiffusionOptions, LayeredModel, LayeredSemiImplicit,
    LayeredSemiImplicitOptions, Levels, ZonalClimatology,
};
use fms_grid::{LatitudeBands, LloydOptions, MeshGeometry, Operators, SphereMesh};

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
    let levels = std::env::var("FMS_LEVELS").map_or(Levels::two(0.5), |v| {
        let interior: Vec<f64> = v.split(',').map(|s| s.parse().unwrap()).collect();
        Levels::new(&interior)
    });
    println!("levels {:?}", levels.interfaces());
    let n_layers = levels.count();
    let dsigma = levels.dsigma().to_vec();
    let model = LayeredModel::new(
        Operators::new(&mesh, &geometry, RADIUS),
        GRAVITY,
        Gas::EARTH_AIR,
        levels,
        mesh.vertices().iter().map(|p| 2.0 * OMEGA * p.z).collect(),
        vec![0.0; mesh.n_cells()],
    );
    let dt = LayeredSemiImplicit::stable_dt(&model, env("FMS_MAX_WIND", 80.0), 0.6);
    let dt = DAY / (DAY / dt).ceil();
    let diffusion = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
    let hs = HeldSuarez::new(&model, HeldSuarezParams::default());
    println!(
        "frequency {k}: {} cells, dt {dt:.0} s, ∇⁴ τ = {:.1} h; friction {:.3}/day (lower)",
        mesh.n_cells(),
        diffusion.e_folding / 3600.0,
        hs.friction()[n_layers - 1] * DAY
    );
    let friction = hs.friction().to_vec();
    let state = hs.initial_state(&model, 0.1);
    let model = model
        .with_hyperdiffusion(diffusion)
        .with_forcing(Box::new(hs));
    let start = model.diagnostics(&state);
    let nb: usize = env("FMS_BANDS", 36);
    let mut climate = ZonalClimatology::new(&model, nb);
    let bands = LatitudeBands::new(&mesh, &geometry, nb);
    // Per band and layer: zonal means of u, v and u·v.
    let mut fluxes = vec![vec![[0.0f64; 3]; n_layers]; nb];
    let mut si =
        LayeredSemiImplicit::new(model, state, dt, LayeredSemiImplicitOptions::default(), 0.5);
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
            for (k, u) in si.state().u.iter().enumerate() {
                let v = geometry.reconstruct_vectors(&mesh, u);
                let (east, north): (Vec<f64>, Vec<f64>) = mesh
                    .cell_centers()
                    .iter()
                    .zip(&v)
                    .map(|(&p, &v)| fms_grid::zonal::east_north(p, v))
                    .unzip();
                let uv: Vec<f64> = east.iter().zip(&north).map(|(e, n)| e * n).collect();
                for (j, f) in [&east, &north, &uv].into_iter().enumerate() {
                    for (b, x) in bands.mean(&geometry, f).into_iter().enumerate() {
                        fluxes[b][k][j] += x;
                    }
                }
            }
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
    let lowest = n_layers - 1;
    println!(" lat   u_up  u_low   T_up  T_low   p_s    | SH mirror u_up u_low T_low");
    for i in 0..m.latitude.len() {
        println!(
            "{:>5.1} {:>6.1} {:>6.1} {:>6.1} {:>6.1} {:>7.1} | {:>6.1} {:>6.1} {:>6.1}",
            m.latitude[i],
            m.wind[0][i],
            m.wind[lowest][i],
            m.temperature[0][i],
            m.temperature[lowest][i],
            m.surface_pressure[i] / 100.0,
            mirror.wind[0][i],
            mirror.wind[lowest][i],
            mirror.temperature[lowest][i],
        );
    }
    print_eddy_fluxes(
        &m.latitude,
        &fluxes,
        climate.samples() as f64,
        &dsigma,
        &friction,
    );
}

/// The eddy momentum flux `[u'v']` per layer and the column momentum budget:
/// in a steady state the friction on the column equals the convergence of the
/// momentum flux, `−Σ m ⟨k_v⟩ [u] = (1/(a cos²φ)) ∂(cos²φ Σ m [uv])/∂φ`.
fn print_eddy_fluxes(
    latitude: &[f64],
    fluxes: &[Vec<[f64; 3]>],
    samples: f64,
    dsigma: &[f64],
    friction: &[f64],
) {
    let n = dsigma.len();
    let mass: Vec<f64> = dsigma.iter().map(|d| d * 1e5 / GRAVITY).collect();
    let mean = |b: usize, k: usize, j: usize| fluxes[b][k][j] / samples;
    let eddy = |b: usize, k: usize| mean(b, k, 2) - mean(b, k, 0) * mean(b, k, 1);
    let column = |b: usize| (0..n).map(|k| mass[k] * mean(b, k, 2)).sum::<f64>();
    let rad: Vec<f64> = latitude.iter().map(|l| l.to_radians()).collect();
    let header: Vec<String> = (0..n).map(|k| format!("uv'_{k}")).collect();
    println!(
        "\n  Eddy momentum flux [u'v'] (m²/s²; positive is poleward in the north) and the column budget (mN/m²)"
    );
    println!(
        "  lat   {}  column_uv'  friction  flux_div",
        header
            .iter()
            .map(|h| format!("{h:>7}"))
            .collect::<Vec<_>>()
            .join(" ")
    );
    for b in 0..latitude.len() {
        let layers: Vec<String> = (0..n).map(|k| format!("{:>7.1}", eddy(b, k))).collect();
        let column_eddy = (0..n).map(|k| dsigma[k] * eddy(b, k)).sum::<f64>();
        let drag = -(0..n)
            .map(|k| mass[k] * friction[k] * mean(b, k, 0))
            .sum::<f64>();
        let divergence = if b == 0 || b + 1 == latitude.len() {
            f64::NAN
        } else {
            let g = |c: usize| fms_math::cos(rad[c]).powi(2) * column(c);
            (g(b + 1) - g(b - 1))
                / (rad[b + 1] - rad[b - 1])
                / (RADIUS * fms_math::cos(rad[b]).powi(2))
        };
        println!(
            "  {:>5.1} {}  {:>9.1}  {:>8.1}  {:>8.1}",
            latitude[b],
            layers.join(" "),
            column_eddy,
            1e3 * drag,
            1e3 * divergence,
        );
    }
}
