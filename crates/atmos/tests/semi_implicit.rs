//! Semi-implicit time stepping: the Helmholtz solver, the Williamson cases at
//! steps several times the explicit limit, and thread-count determinism.

use std::f64::consts::{FRAC_PI_4, PI};
use std::sync::OnceLock;

use fms_atmos::helmholtz::dot;
use fms_atmos::williamson::{Case, DAY, GRAVITY, RADIUS, error_norms};
use fms_atmos::{Helmholtz, SemiImplicit, SemiImplicitOptions, ShallowWater, State};
use fms_geo::coords::to_lat_lon;
use fms_geo::{LloydOptions, MeshGeometry, Operators, SphereMesh};

/// Courant number for the semi-implicit step: measured stable to ~0.75 on
/// case 6 (the fastest winds), so this leaves a margin.
const COURANT: f64 = 0.6;

fn mesh() -> &'static (SphereMesh, MeshGeometry) {
    static M: OnceLock<(SphereMesh, MeshGeometry)> = OnceLock::new();
    M.get_or_init(|| {
        let (m, _) = SphereMesh::centroidal(16, LloydOptions::default());
        let g = MeshGeometry::new(&m);
        (m, g)
    })
}

fn noise(n: usize, seed: u64) -> Vec<f64> {
    let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    (0..n)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            (x >> 11) as f64 / (1u64 << 52) as f64 - 1.0
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Helmholtz solver
// ---------------------------------------------------------------------------

#[test]
fn helmholtz_recovers_a_manufactured_solution() {
    let (m, g) = mesh();
    let ops = Operators::new(m, g, RADIUS);
    // c = Δt² g H̄ for a 3000 s step over 10 km of fluid.
    let c = 3000.0f64.powi(2) * GRAVITY * 1.0e4;
    let truth = noise(m.n_cells(), 1);
    // b = (I − c∇²) x.
    let mut grad = vec![0.0; m.n_edges()];
    let mut lap = vec![0.0; m.n_cells()];
    ops.gradient(&truth, &mut grad);
    ops.divergence(&grad, &mut lap);
    let b: Vec<f64> = truth.iter().zip(&lap).map(|(x, l)| x - c * l).collect();

    let mut solver = Helmholtz::new(ops, c);
    solver.tolerance = 1e-12;
    let mut x = vec![0.0; m.n_cells()];
    let report = solver.solve(&b, &mut x);
    let err = x
        .iter()
        .zip(&truth)
        .fold(0.0f64, |a, (x, t)| a.max((x - t).abs()));
    eprintln!("{report:?}, max error {err:.1e}");
    assert!(report.relative_residual <= 1e-12 && report.iterations < 200);
    assert!(err < 1e-9, "max error {err}");
}

#[test]
fn helmholtz_with_no_coupling_is_the_identity() {
    let (m, g) = mesh();
    let mut solver = Helmholtz::new(Operators::new(m, g, RADIUS), 0.0);
    let b = noise(m.n_cells(), 2);
    let mut x = vec![0.0; m.n_cells()];
    let report = solver.solve(&b, &mut x);
    assert!(report.iterations <= 1);
    assert!(x.iter().zip(&b).all(|(x, b)| (x - b).abs() < 1e-15));
}

#[test]
fn dot_products_are_order_independent_of_threads() {
    let a = noise(100_003, 3);
    let b = noise(100_003, 4);
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let four = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    assert_eq!(
        one.install(|| dot(&a, &b)).to_bits(),
        four.install(|| dot(&a, &b)).to_bits()
    );
}

// ---------------------------------------------------------------------------
// Williamson cases at the semi-implicit step
// ---------------------------------------------------------------------------

struct Run {
    state: State,
    topography: Vec<f64>,
    mass: f64,
    energy: f64,
    /// The step as a multiple of the explicit RK4 limit.
    speedup: f64,
    wave4_shift: f64,
}

fn run(case: Case, days: u32) -> Run {
    let (m, g) = mesh();
    let setup = case.setup(m, g);
    let topography = setup.topography.clone();
    let model = ShallowWater::new(
        Operators::new(m, g, RADIUS),
        GRAVITY,
        setup.coriolis,
        setup.topography,
    );
    let explicit = model.stable_dt(&setup.state, 0.5);
    let steps = (DAY / SemiImplicit::stable_dt(&model, case.max_wind(), COURANT)).ceil() as u64;
    let dt = DAY / steps as f64;
    let start = model.diagnostics(&setup.state);
    let options = SemiImplicitOptions::for_state(&setup.state);
    let mut si = SemiImplicit::new(model, setup.state, dt, options, 0.5);
    let (mut last, mut shift) = (wave4_phase(m, si.state()), 0.0);
    for day in 1..=u64::from(days) {
        while si.steps() < day * steps {
            si.step();
        }
        let p = wave4_phase(m, si.state());
        shift += (p - last + PI).rem_euclid(2.0 * PI) - PI;
        last = p;
    }
    let end = si.model().diagnostics(si.state());
    Run {
        state: si.state().clone(),
        topography,
        mass: end.mass / start.mass - 1.0,
        energy: end.energy / start.energy - 1.0,
        speedup: dt / explicit,
        wave4_shift: -shift / 4.0,
    }
}

fn wave4_phase(m: &SphereMesh, state: &State) -> f64 {
    let (mut re, mut im) = (0.0, 0.0);
    for (&p, &h) in m.cell_centers().iter().zip(&state.h) {
        let (lat, lon) = to_lat_lon(p);
        if (30f64.to_radians()..50f64.to_radians()).contains(&lat) {
            let (s, c) = fms_math::sin_cos(4.0 * lon);
            re += h * c;
            im -= h * s;
        }
    }
    fms_math::atan2(im, re)
}

fn check_conservation(r: &Run) {
    assert!(r.mass.abs() < 1e-13, "mass drift {}", r.mass);
    // Leapfrog with the RAW filter isn't energy-conserving, but drifts slowly.
    // Measured: 1e-7 to 5e-6 over the standard runs.
    assert!(r.energy.abs() < 2e-5, "energy drift {}", r.energy);
}

#[test]
fn case_2_at_the_semi_implicit_step() {
    // Explicit RK4 at k = 16 gives l2 = 5.4e-4 (α = 0) and 3.7e-4 (α = 45°).
    for (alpha, rk4_l2) in [(0.0, 5.44e-4), (FRAC_PI_4, 3.65e-4)] {
        let case = Case::SteadyZonal { alpha };
        let r = run(case, 5);
        check_conservation(&r);
        let (m, g) = mesh();
        let surface = &r.state.h;
        let exact: Vec<f64> = m
            .cell_centers()
            .iter()
            .map(|&p| case.surface_height(p))
            .collect();
        let (_, l2, linf) = error_norms(surface, &exact, &g.cell_areas);
        eprintln!(
            "α = {:.0}°: step {:.1}× explicit, l2 {l2:.2e}, l∞ {linf:.2e}",
            alpha.to_degrees(),
            r.speedup
        );
        assert!(r.speedup > 4.0, "step only {:.1}× explicit", r.speedup);
        assert!(l2 < 1.1 * rk4_l2 && linf < 2e-3, "l2 {l2}, l∞ {linf}");
    }
}

#[test]
fn case_5_at_the_semi_implicit_step() {
    let r = run(Case::Mountain, 15);
    check_conservation(&r);
    let surface: Vec<f64> = r
        .state
        .h
        .iter()
        .zip(&r.topography)
        .map(|(h, b)| h + b)
        .collect();
    let (lo, hi) = surface
        .iter()
        .fold((f64::MAX, f64::MIN), |(a, b), &h| (a.min(h), b.max(h)));
    eprintln!(
        "step {:.1}× explicit, day 15 h [{lo:.0}, {hi:.0}]",
        r.speedup
    );
    assert!(r.speedup > 4.0);
    assert!((4980.0..5100.0).contains(&lo) && (5900.0..5990.0).contains(&hi));
}

#[test]
fn case_6_at_the_semi_implicit_step() {
    let r = run(Case::RossbyHaurwitz, 14);
    check_conservation(&r);
    let ratio = r.wave4_shift / (14.0 * DAY) / Case::rossby_haurwitz_phase_speed();
    eprintln!(
        "step {:.1}× explicit, phase speed {:.1}% of theory",
        r.speedup,
        100.0 * ratio
    );
    assert!(r.speedup > 4.0);
    assert!((0.9..1.0).contains(&ratio), "phase speed ratio {ratio}");
}

#[test]
fn results_do_not_depend_on_the_thread_count() {
    let (m, g) = mesh();
    let simulate = || {
        let setup = Case::Mountain.setup(m, g);
        let model = ShallowWater::new(
            Operators::new(m, g, RADIUS),
            GRAVITY,
            setup.coriolis,
            setup.topography,
        );
        let options = SemiImplicitOptions::for_state(&setup.state);
        let mut si = SemiImplicit::new(model, setup.state, 3600.0, options, 0.5);
        for _ in 0..24 {
            si.step();
        }
        si.state().clone()
    };
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let four = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    let (a, b) = (one.install(simulate), four.install(simulate));
    let bits = |s: &State| -> Vec<u64> { s.h.iter().chain(&s.u).map(|v| v.to_bits()).collect() };
    assert_eq!(bits(&a), bits(&b));
}
