//! Baroclinic instability in the two-layer model: growth rates of small
//! wavenumber-m perturbations on a balanced sheared flow, against two-level
//! quasi-geostrophic theory.
//!
//! ```sh
//! cargo run --release -p fms-atmos --example baroclinic [frequency]
//! ```

use fms_atmos::williamson::{GRAVITY, OMEGA, RADIUS};
use fms_atmos::{Gas, ShearedFlow, TwoLayer, TwoLayerState};
use fms_geo::coords::to_lat_lon;
use fms_geo::{LloydOptions, MeshGeometry, Operators, SphereMesh};

const DAY: f64 = 86_400.0;

fn main() {
    let k: u32 = std::env::args().nth(1).map_or(32, |s| s.parse().unwrap());
    let (mesh, _) = SphereMesh::centroidal(k, LloydOptions::default());
    let geometry = MeshGeometry::new(&mesh);
    println!("frequency {k}: {} cells", mesh.n_cells());
    let coriolis: Vec<f64> = mesh.vertices().iter().map(|p| 2.0 * OMEGA * p.z).collect();
    let model = || {
        TwoLayer::new(
            Operators::new(&mesh, &geometry, RADIUS),
            GRAVITY,
            Gas::EARTH_AIR,
            0.5,
            coriolis.clone(),
            vec![0.0; mesh.n_cells()],
        )
    };

    // 1. Balance: the unperturbed flow should stay put.
    let flow = |shear: f64| ShearedFlow {
        speed: [10.0 + shear, 10.0],
        stability: 40.0,
        theta_equator: 315.0,
        ps_equator: 1.0e5,
        omega: OMEGA,
    };
    {
        let mut m = model();
        let s0 = flow(30.0).state(&m);
        let mut s = s0.clone();
        let dt = step(&m, &s);
        for _ in 0..(DAY / dt) as u32 {
            m.step(&mut s, dt);
        }
        let du =
            s.u.iter()
                .zip(&s0.u)
                .flat_map(|(a, b)| a.iter().zip(b))
                .fold(0.0f64, |w, (a, b)| w.max((a - b).abs()));
        let dps =
            s.ps.iter()
                .zip(&s0.ps)
                .fold(0.0f64, |w, (a, b)| w.max((a - b).abs()));
        println!("balance after 1 day: max |Δu| {du:.3} m/s, max |Δp_s| {dps:.1} Pa");
    }

    // 2. Growth rates.
    let m0 = model();
    let shears: Vec<f64> = std::env::var("FMS_SHEARS").map_or(vec![10.0, 20.0, 30.0, 45.0], |v| {
        v.split(',').map(|x| x.parse().unwrap()).collect()
    });
    let waves: Vec<u32> = std::env::var("FMS_WAVES").map_or(vec![4, 5, 6, 7, 8, 10], |v| {
        v.split(',').map(|x| x.parse().unwrap()).collect()
    });
    for shear in shears {
        let f = flow(shear);
        let (_, critical) = f.phillips_growth_rate(&m0, 6, 45f64.to_radians());
        println!("shear {shear} m/s (critical at 45°: {critical:.1} m/s)");
        let control = run(&mesh, model(), f, None);
        for &wave in &waves {
            let (theory, _) = f.phillips_growth_rate(&m0, wave, 45f64.to_radians());
            let perturbed = run(&mesh, model(), f, Some(wave));
            let series: Vec<f64> = perturbed
                .iter()
                .zip(&control)
                .map(|(p, c)| amplitude(&mesh, &geometry, p, c, wave))
                .collect();
            let measured = fit_growth(&series);
            if std::env::var("FMS_TRACE").is_ok() {
                let line: Vec<String> = series
                    .iter()
                    .step_by(4)
                    .map(|a| format!("{a:.0e}"))
                    .collect();
                println!("    m={wave} amplitude (K) daily: {}", line.join(" "));
            }
            println!(
                "  m = {wave:>2}: measured {:>6.3} /day   theory {:>6.3} /day   ratio {}",
                measured * DAY,
                theory * DAY,
                if theory > 0.0 {
                    format!("{:.2}", measured / theory)
                } else {
                    "-".into()
                }
            );
        }
    }
}

/// Samples per day, and run length (days).
const PER_DAY: u32 = 4;
const DAYS: u32 = 10;
/// Seed amplitude (K).
const SEED: f64 = 0.01;

fn step(model: &TwoLayer, state: &TwoLayerState) -> f64 {
    let dt = model.stable_dt(state, 0.5);
    DAY / (DAY / dt).ceil()
}

/// Run the sheared flow, seeded with a wavenumber-`wave` θ bump at 45° N if
/// given, and return θ̂ every 6 hours.
fn run(
    mesh: &SphereMesh,
    mut model: TwoLayer,
    flow: ShearedFlow,
    wave: Option<u32>,
) -> Vec<Vec<f64>> {
    let mut state = flow.state(&model);
    if let Some(wave) = wave {
        for k in 0..2 {
            for (i, &p) in mesh.cell_centers().iter().enumerate() {
                let (lat, lon) = to_lat_lon(p);
                let bump =
                    fms_math::exp(-((lat - 45f64.to_radians()) / 15f64.to_radians()).powi(2));
                state.theta_ps[k][i] +=
                    state.ps[i] * SEED * bump * fms_math::cos(f64::from(wave) * lon);
            }
        }
    }
    let dt = step(&model, &flow.state(&model));
    let per_sample = (DAY / f64::from(PER_DAY) / dt).round() as u32;
    let theta_hat = |s: &TwoLayerState| -> Vec<f64> {
        (0..s.ps.len())
            .map(|i| 0.5 * (s.theta_ps[0][i] + s.theta_ps[1][i]) / s.ps[i])
            .collect()
    };
    let mut out = vec![theta_hat(&state)];
    for _ in 0..PER_DAY * DAYS {
        for _ in 0..per_sample {
            model.step(&mut state, dt);
        }
        out.push(theta_hat(&state));
    }
    out
}

/// Least-squares exponential growth rate (s⁻¹), fitted from day 2 (once the
/// growing mode dominates the seed) until the wave reaches 0.5 K (before it
/// turns nonlinear).
fn fit_growth(series: &[f64]) -> f64 {
    let dt = DAY / f64::from(PER_DAY);
    let fit: Vec<(f64, f64)> = series
        .iter()
        .enumerate()
        .skip(2 * PER_DAY as usize)
        .take_while(|(_, a)| **a < 0.5)
        .map(|(n, a)| (n as f64 * dt, fms_math::ln(*a)))
        .collect();
    if fit.len() < 3 {
        return f64::NAN;
    }
    let n = fit.len() as f64;
    let (mt, my) = fit
        .iter()
        .fold((0.0, 0.0), |(a, b), (t, y)| (a + t / n, b + y / n));
    let (num, den) = fit.iter().fold((0.0, 0.0), |(a, b), (t, y)| {
        (a + (t - mt) * (y - my), b + (t - mt).powi(2))
    });
    num / den
}

/// Amplitude (K) of zonal wavenumber `wave` in `perturbed − control` over
/// 30–60° N, area-weighted.
fn amplitude(
    mesh: &SphereMesh,
    geometry: &MeshGeometry,
    perturbed: &[f64],
    control: &[f64],
    wave: u32,
) -> f64 {
    let (mut re, mut im, mut area) = (0.0, 0.0, 0.0);
    for (i, &p) in mesh.cell_centers().iter().enumerate() {
        let (lat, lon) = to_lat_lon(p);
        if (30f64.to_radians()..60f64.to_radians()).contains(&lat) {
            let a = geometry.cell_areas[i];
            let d = perturbed[i] - control[i];
            let (s, c) = fms_math::sin_cos(f64::from(wave) * lon);
            re += a * d * c;
            im += a * d * s;
            area += a;
        }
    }
    2.0 * (re * re + im * im).sqrt() / area
}
