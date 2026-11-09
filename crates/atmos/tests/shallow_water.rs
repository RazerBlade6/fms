//! Shallow-water model: exact discrete properties, and the Williamson et al.
//! (1992) test cases 2, 5 and 6.

use std::f64::consts::{FRAC_PI_4, PI};
use std::sync::OnceLock;

use fms_atmos::williamson::{Case, DAY, GRAVITY, RADIUS, error_norms};
use fms_atmos::{ShallowWater, State};
use fms_grid::coords::to_lat_lon;
use fms_grid::{LloydOptions, MeshGeometry, Operators, SphereMesh};

fn mesh(k: u32) -> &'static (SphereMesh, MeshGeometry) {
    static M: [OnceLock<(SphereMesh, MeshGeometry)>; 2] = [OnceLock::new(), OnceLock::new()];
    let slot = match k {
        16 => &M[0],
        32 => &M[1],
        _ => panic!("unsupported test frequency"),
    };
    slot.get_or_init(|| {
        let (m, _) = SphereMesh::centroidal(k, LloydOptions::default());
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

/// Result of running a case for its standard length.
struct Run {
    state: State,
    topography: Vec<f64>,
    mass: f64,
    energy: f64,
    /// Unwrapped eastward shift of the wave-4 pattern (radians).
    wave4_shift: f64,
}

/// Run `case` for `days` at Courant number 0.5, with a whole number of steps
/// per day.
fn run(k: u32, case: Case, days: u32) -> Run {
    let (m, g) = mesh(k);
    let setup = case.setup(m, g);
    let topography = setup.topography.clone();
    let mut model = ShallowWater::new(
        Operators::new(m, g, RADIUS),
        GRAVITY,
        setup.coriolis,
        setup.topography,
    );
    let mut state = setup.state;
    let steps = (DAY / model.stable_dt(&state, 0.5)).ceil() as u32;
    let dt = DAY / f64::from(steps);
    let start = model.diagnostics(&state);
    let (mut last, mut shift) = (wave4_phase(m, &state), 0.0);
    for _ in 0..days {
        for _ in 0..steps {
            model.step(&mut state, dt);
        }
        let p = wave4_phase(m, &state);
        shift += (p - last + PI).rem_euclid(2.0 * PI) - PI;
        last = p;
    }
    let end = model.diagnostics(&state);
    Run {
        state,
        topography,
        mass: end.mass / start.mass - 1.0,
        energy: end.energy / start.energy - 1.0,
        wave4_shift: -shift / 4.0,
    }
}

/// Phase of the zonal wavenumber-4 component of `h` between 30° and 50° N.
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

// ---------------------------------------------------------------------------
// Exact discrete properties
// ---------------------------------------------------------------------------

#[test]
fn lake_at_rest_stays_at_rest() {
    // Flat free surface over the case-5 mountain, no flow: every tendency is
    // zero, so topography alone doesn't drive spurious flow.
    let (m, g) = mesh(16);
    let topography: Vec<f64> = m
        .cell_centers()
        .iter()
        .map(|&p| Case::Mountain.topography(p))
        .collect();
    let state = State {
        h: topography.iter().map(|b| 5000.0 - b).collect(),
        u: vec![0.0; m.n_edges()],
    };
    let coriolis = m
        .vertices()
        .iter()
        .map(|&p| Case::Mountain.coriolis(p))
        .collect();
    let mut model = ShallowWater::new(Operators::new(m, g, RADIUS), GRAVITY, coriolis, topography);
    let mut tendency = state.clone();
    model.tendency(&state, &mut tendency);
    assert!(tendency.h.iter().all(|&d| d == 0.0));
    // h + b may round to 5000 ± 1 ulp, so allow a round-off-sized gradient.
    let worst = tendency.u.iter().fold(0.0f64, |a, &d| a.max(d.abs()));
    assert!(worst < 1e-15, "spurious acceleration {worst} m s-2");
}

#[test]
fn semi_discrete_energy_is_conserved() {
    // For any state, the spatial scheme's tendencies change the total energy
    // E = Σᵢ Aᵢ [hᵢKᵢ + g hᵢ(hᵢ/2 + bᵢ)] by exactly zero:
    // dE/dt = Σₑ lₑdₑ hₑ uₑ ∂uₑ + Σᵢ Aᵢ (Kᵢ + g(hᵢ + bᵢ)) ∂hᵢ.
    let (m, g) = mesh(16);
    let ops = Operators::new(m, g, RADIUS);
    let h: Vec<f64> = noise(m.n_cells(), 1)
        .iter()
        .map(|x| 5000.0 + 500.0 * x)
        .collect();
    let u: Vec<f64> = noise(m.n_edges(), 2).iter().map(|x| 30.0 * x).collect();
    let b: Vec<f64> = noise(m.n_cells(), 3)
        .iter()
        .map(|x| 1000.0 + 1000.0 * x)
        .collect();
    let coriolis = m.vertices().iter().map(|p| 2.0 * 7.292e-5 * p.z).collect();
    let state = State { h, u };
    let mut model = ShallowWater::new(ops, GRAVITY, coriolis, b.clone());
    let mut d = state.clone();
    model.tendency(&state, &mut d);

    let r2 = RADIUS * RADIUS;
    let mut h_edge = vec![0.0; m.n_edges()];
    let mut kinetic = vec![0.0; m.n_cells()];
    ops.cell_to_edge(&state.h, &mut h_edge);
    ops.kinetic_energy(&state.u, &mut kinetic);
    let edge_terms: Vec<f64> = (0..m.n_edges())
        .map(|e| g.edge_lengths[e] * g.dual_edge_lengths[e] * r2 * h_edge[e] * state.u[e] * d.u[e])
        .collect();
    let cell_terms: Vec<f64> = (0..m.n_cells())
        .map(|i| g.cell_areas[i] * r2 * (kinetic[i] + GRAVITY * (state.h[i] + b[i])) * d.h[i])
        .collect();
    let total: f64 = edge_terms.iter().chain(&cell_terms).sum();
    let scale: f64 = edge_terms.iter().chain(&cell_terms).map(|x| x.abs()).sum();
    assert!(total.abs() < 1e-13 * scale, "dE/dt = {total} of {scale}");
}

// ---------------------------------------------------------------------------
// Williamson test cases
// ---------------------------------------------------------------------------

fn tc2_errors(k: u32, alpha: f64) -> (f64, f64, f64) {
    let case = Case::SteadyZonal { alpha };
    let r = run(k, case, case.days() as u32);
    assert!(r.mass.abs() < 1e-13, "mass drift {}", r.mass);
    assert!(r.energy.abs() < 1e-8, "energy drift {}", r.energy);
    let (m, g) = mesh(k);
    let exact: Vec<f64> = m
        .cell_centers()
        .iter()
        .map(|&p| case.surface_height(p))
        .collect();
    error_norms(&r.state.h, &exact, &g.cell_areas)
}

#[test]
fn case_2_steady_zonal_flow() {
    // Measured (l1, l2, l∞) after 5 days:
    //   k = 16: α = 0: 4.2e-4, 5.4e-4, 1.4e-3; α = 45°: 3.0e-4, 3.7e-4, 9.2e-4
    //   k = 32: α = 0: 1.1e-4, 1.4e-4, 5.9e-4; α = 45°: 7.6e-5, 9.6e-5, 5.0e-4
    // l1 and l2 converge at second order; l∞ more slowly, limited by the
    // pentagons (OBSERVATIONS.md O-1).
    for alpha in [0.0, FRAC_PI_4] {
        let coarse = tc2_errors(16, alpha);
        let fine = tc2_errors(32, alpha);
        eprintln!(
            "α = {:.0}°: k=16 {coarse:?}  k=32 {fine:?}",
            alpha.to_degrees()
        );
        assert!(fine.1 < 2e-4 && fine.2 < 1e-3, "{fine:?}");
        assert!(coarse.1 / fine.1 > 3.5, "l2 should converge at ~2nd order");
        assert!(coarse.2 / fine.2 > 1.5, "l∞ should decrease");
    }
}

#[test]
fn case_5_flow_over_a_mountain() {
    let r = run(16, Case::Mountain, 15);
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
    // Reference solutions at day 15 span about 5050–5950 m. Measured:
    // [5043, 5958] at k = 16, [5036, 5951] at k = 32.
    eprintln!(
        "day 15 surface height [{lo:.0}, {hi:.0}] m, energy drift {:.1e}",
        r.energy
    );
    assert!((4980.0..5100.0).contains(&lo) && (5900.0..5990.0).contains(&hi));
    assert!(r.mass.abs() < 1e-13, "mass drift {}", r.mass);
    assert!(r.energy.abs() < 1e-7, "energy drift {}", r.energy);
    assert!(r.state.h.iter().all(|&h| h > 0.0));
}

#[test]
fn case_6_rossby_haurwitz_wave() {
    let r = run(16, Case::RossbyHaurwitz, 14);
    let speed = r.wave4_shift / (14.0 * DAY);
    let ratio = speed / Case::rossby_haurwitz_phase_speed();
    // Shallow-water Rossby–Haurwitz waves move a little slower than the
    // non-divergent theory. Measured: 97% at k = 16, 95% at k = 32 and 64.
    eprintln!(
        "phase speed {:.2}°/day ({:.1}% of theory), energy drift {:.1e}",
        speed.to_degrees() * DAY,
        100.0 * ratio,
        r.energy
    );
    assert!((0.9..1.0).contains(&ratio), "phase speed ratio {ratio}");
    assert!(r.mass.abs() < 1e-13, "mass drift {}", r.mass);
    assert!(r.energy.abs() < 1e-7, "energy drift {}", r.energy);
    // The wave keeps its shape: the wave-4 amplitude band stays in range.
    let (lo, hi) = r
        .state
        .h
        .iter()
        .fold((f64::MAX, f64::MIN), |(a, b), &h| (a.min(h), b.max(h)));
    assert!(
        (8000.0..8700.0).contains(&lo) && (10300.0..10700.0).contains(&hi),
        "[{lo}, {hi}]"
    );
}
