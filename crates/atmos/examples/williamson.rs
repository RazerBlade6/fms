//! Run the Williamson shallow-water test cases and print error norms,
//! conservation and timing.
//!
//! ```sh
//! cargo run --release -p fms-atmos --example williamson
//! ```

use std::time::Instant;

use fms_atmos::williamson::{Case, DAY, GRAVITY, RADIUS, error_norms};
use fms_atmos::{ShallowWater, State};
use fms_geo::coords::to_lat_lon;
use fms_geo::{LloydOptions, MeshGeometry, Operators, SphereMesh};

/// Courant number for the time step.
const COURANT: f64 = 0.5;

fn main() {
    for k in [16u32, 32, 64] {
        let (mesh, _) = SphereMesh::centroidal(k, LloydOptions::default());
        let geometry = MeshGeometry::new(&mesh);
        println!("--- frequency {k}: {} cells", mesh.n_cells());
        for case in [
            Case::SteadyZonal { alpha: 0.0 },
            Case::SteadyZonal {
                alpha: std::f64::consts::FRAC_PI_4,
            },
            Case::Mountain,
            Case::RossbyHaurwitz,
        ] {
            run(&mesh, &geometry, case);
        }
    }
}

fn run(mesh: &SphereMesh, geometry: &MeshGeometry, case: Case) {
    let ops = Operators::new(mesh, geometry, RADIUS);
    let setup = case.setup(mesh, geometry);
    let mut model = ShallowWater::new(ops, GRAVITY, setup.coriolis, setup.topography.clone());
    let mut state = setup.state;
    let dt = model.stable_dt(&state, COURANT);
    let steps_per_day = (DAY / dt).ceil() as u64;
    let dt = DAY / steps_per_day as f64;
    let start = model.diagnostics(&state);
    let phase0 = wave4_phase(mesh, &state);
    let mut shift = 0.0;
    let mut last = phase0;
    let t = Instant::now();
    for _ in 0..case.days() as u64 {
        for _ in 0..steps_per_day {
            model.step(&mut state, dt);
        }
        // Unwrap the wave-4 phase day by day (it moves ~50° of phase per day).
        let p = wave4_phase(mesh, &state);
        let mut d = p - last;
        d = (d + std::f64::consts::PI).rem_euclid(2.0 * std::f64::consts::PI)
            - std::f64::consts::PI;
        shift += d;
        last = p;
    }
    let elapsed = t.elapsed().as_secs_f64();
    let end = model.diagnostics(&state);
    let surface: Vec<f64> = state
        .h
        .iter()
        .zip(&setup.topography)
        .map(|(h, b)| h + b)
        .collect();
    let (lo, hi) = surface
        .iter()
        .fold((f64::MAX, f64::MIN), |(a, b), &h| (a.min(h), b.max(h)));
    print!(
        "{:<12} dt {dt:>6.1}s  {:>6} steps in {elapsed:>6.2}s  mass {:+.1e}  energy {:+.1e}  enstrophy {:+.1e}  h [{lo:.0}, {hi:.0}]",
        format!("{}{}", case.name(), alpha_label(case)),
        steps_per_day * case.days() as u64,
        end.mass / start.mass - 1.0,
        end.energy / start.energy - 1.0,
        end.potential_enstrophy / start.potential_enstrophy - 1.0,
    );
    match case {
        Case::SteadyZonal { .. } => {
            let exact: Vec<f64> = mesh
                .cell_centers()
                .iter()
                .map(|&p| case.surface_height(p))
                .collect();
            let areas = &geometry.cell_areas;
            let (l1, l2, li) = error_norms(&surface, &exact, areas);
            println!("  l1 {l1:.2e} l2 {l2:.2e} linf {li:.2e}");
        }
        Case::RossbyHaurwitz => {
            // The wave-4 phase moves at −4·(angular speed).
            let speed = -shift / 4.0 / (case.days() * DAY);
            let theory = Case::rossby_haurwitz_phase_speed();
            println!(
                "  phase speed {:.3}°/day ({:.0}% of non-divergent theory)",
                speed.to_degrees() * DAY,
                100.0 * speed / theory
            );
        }
        Case::Mountain => println!(),
    }
    let _ = phase0;
}

fn alpha_label(case: Case) -> String {
    match case {
        Case::SteadyZonal { alpha } => format!(" α={:.0}°", alpha.to_degrees()),
        _ => String::new(),
    }
}

/// Phase (radians) of the zonal wavenumber-4 component of `h` in the band
/// 30°–50° N.
fn wave4_phase(mesh: &SphereMesh, state: &State) -> f64 {
    let (mut re, mut im) = (0.0, 0.0);
    for (&p, &h) in mesh.cell_centers().iter().zip(&state.h) {
        let (lat, lon) = to_lat_lon(p);
        if (30f64.to_radians()..50f64.to_radians()).contains(&lat) {
            let (s, c) = fms_math::sin_cos(4.0 * lon);
            re += h * c;
            im -= h * s;
        }
    }
    fms_math::atan2(im, re)
}
