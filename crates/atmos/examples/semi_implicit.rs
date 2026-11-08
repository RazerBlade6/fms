//! Semi-implicit stepping of the Williamson cases at multiples of the explicit
//! (RK4) time step: errors, conservation, solver work and speed.
//!
//! ```sh
//! cargo run --release -p fms-atmos --example semi_implicit
//! ```

use std::time::Instant;

use fms_atmos::williamson::{Case, DAY, GRAVITY, RADIUS, error_norms};
use fms_atmos::{SemiImplicit, SemiImplicitOptions, ShallowWater, State};
use fms_geo::coords::to_lat_lon;
use fms_geo::{LloydOptions, MeshGeometry, Operators, SphereMesh};

fn main() {
    let k: u32 = std::env::args().nth(1).map_or(32, |s| s.parse().unwrap());
    let (mesh, _) = SphereMesh::centroidal(k, LloydOptions::default());
    let geometry = MeshGeometry::new(&mesh);
    println!("frequency {k}: {} cells", mesh.n_cells());
    for case in [
        Case::SteadyZonal { alpha: 0.0 },
        Case::SteadyZonal {
            alpha: std::f64::consts::FRAC_PI_4,
        },
        Case::Mountain,
        Case::RossbyHaurwitz,
    ] {
        for multiple in [1.0, 2.0, 4.0, 6.0, 8.0] {
            run(&mesh, &geometry, case, multiple);
        }
    }
}

fn run(mesh: &SphereMesh, geometry: &MeshGeometry, case: Case, multiple: f64) {
    let setup = case.setup(mesh, geometry);
    let topography = setup.topography.clone();
    let model = ShallowWater::new(
        Operators::new(mesh, geometry, RADIUS),
        GRAVITY,
        setup.coriolis,
        setup.topography,
    );
    let explicit_dt = model.stable_dt(&setup.state, 0.5);
    let steps_per_day = (DAY / (multiple * explicit_dt)).ceil() as u64;
    let dt = DAY / steps_per_day as f64;
    let start = model.diagnostics(&setup.state);
    let phase0 = wave4_phase(mesh, &setup.state);
    let options = SemiImplicitOptions::for_state(&setup.state);
    let t = Instant::now();
    let mut si = SemiImplicit::new(model, setup.state, dt, options, 0.5);
    let total = steps_per_day * case.days() as u64;
    let (mut last, mut shift, mut blown) = (phase0, 0.0, false);
    while si.steps() < total {
        si.step();
        if si.steps() % steps_per_day == 0 {
            let p = wave4_phase(mesh, si.state());
            shift += (p - last + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
                - std::f64::consts::PI;
            last = p;
            if !si.state().h.iter().all(|h| h.is_finite() && *h > 0.0) {
                blown = true;
                break;
            }
        }
    }
    let elapsed = t.elapsed().as_secs_f64();
    let label = format!(
        "{}{} dt×{multiple}",
        case.name(),
        match case {
            Case::SteadyZonal { alpha } => format!(" α={:.0}°", alpha.to_degrees()),
            _ => String::new(),
        }
    );
    if blown {
        println!("{label:<26} dt {dt:>6.0}s  BLEW UP at step {}", si.steps());
        return;
    }
    let end = si.model().diagnostics(si.state());
    let stats = si.solver_stats();
    let surface: Vec<f64> = si
        .state()
        .h
        .iter()
        .zip(&topography)
        .map(|(h, b)| h + b)
        .collect();
    let (lo, hi) = surface
        .iter()
        .fold((f64::MAX, f64::MIN), |(a, b), &h| (a.min(h), b.max(h)));
    print!(
        "{label:<26} dt {dt:>6.0}s {total:>5} steps {elapsed:>5.2}s  CG {:.1} it/solve  mass {:+.1e} energy {:+.1e}  h [{lo:.0}, {hi:.0}]",
        stats.iterations as f64 / stats.solves as f64,
        end.mass / start.mass - 1.0,
        end.energy / start.energy - 1.0,
    );
    match case {
        Case::SteadyZonal { .. } => {
            let exact: Vec<f64> = mesh
                .cell_centers()
                .iter()
                .map(|&p| case.surface_height(p))
                .collect();
            let (_, l2, li) = error_norms(&surface, &exact, &geometry.cell_areas);
            println!("  l2 {l2:.2e} linf {li:.2e}");
        }
        Case::RossbyHaurwitz => {
            let speed = -shift / 4.0 / (case.days() * DAY);
            println!(
                "  phase {:.0}%",
                100.0 * speed / Case::rossby_haurwitz_phase_speed()
            );
        }
        Case::Mountain => println!(),
    }
}

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
