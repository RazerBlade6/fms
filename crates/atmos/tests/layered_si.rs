//! Semi-implicit layered stepping: the RK4 acceptance tests (conservation,
//! rest over a mountain, baroclinic instability) rerun at a step several
//! times the explicit limit, on two layers and on more, plus thread-count
//! determinism.

use std::sync::OnceLock;

use fms_atmos::williamson::{Case, GRAVITY, OMEGA, RADIUS};
use fms_atmos::{
    Gas, HyperdiffusionOptions, LayeredModel, LayeredSemiImplicit, LayeredSemiImplicitOptions,
    LayeredState, Levels, ShearedFlow,
};
use fms_grid::coords::to_lat_lon;
use fms_grid::{LloydOptions, MeshGeometry, Operators, SphereMesh};

const GAS: Gas = Gas::EARTH_AIR;
const DAY: f64 = 86_400.0;
/// Design maximum wind and Courant number for the semi-implicit step.
const MAX_WIND: f64 = 60.0;
const COURANT: f64 = 0.6;

fn mesh() -> &'static (SphereMesh, MeshGeometry) {
    static M: OnceLock<(SphereMesh, MeshGeometry)> = OnceLock::new();
    M.get_or_init(|| {
        let (m, _) = SphereMesh::centroidal(16, LloydOptions::default());
        let g = MeshGeometry::new(&m);
        (m, g)
    })
}

fn model(mountain: bool) -> LayeredModel<'static> {
    model_with(mountain, Levels::two(0.5))
}

fn model_with(mountain: bool, levels: Levels) -> LayeredModel<'static> {
    let (m, g) = mesh();
    let phi_s = if mountain {
        m.cell_centers()
            .iter()
            .map(|&p| GRAVITY * Case::Mountain.topography(p))
            .collect()
    } else {
        vec![0.0; m.n_cells()]
    };
    LayeredModel::new(
        Operators::new(m, g, RADIUS),
        GRAVITY,
        GAS,
        levels,
        m.vertices().iter().map(|p| 2.0 * OMEGA * p.z).collect(),
        phi_s,
    )
}

fn sheared(shear: f64) -> ShearedFlow {
    ShearedFlow {
        speed: [10.0 + shear, 10.0],
        stability: 40.0,
        theta_equator: 315.0,
        ps_equator: 1.0e5,
        omega: OMEGA,
    }
}

/// A semi-implicit stepper at the advective/inertial limit, with a whole
/// number of steps per day.
fn stepper(model: LayeredModel<'static>, state: LayeredState) -> LayeredSemiImplicit<'static> {
    let dt = LayeredSemiImplicit::stable_dt(&model, MAX_WIND, COURANT);
    let dt = DAY / (DAY / dt).ceil();
    LayeredSemiImplicit::new(model, state, dt, LayeredSemiImplicitOptions::default(), 0.5)
}

fn run_days(si: &mut LayeredSemiImplicit, days: f64) {
    let target = (days * DAY / si.dt()).round() as u64;
    while si.steps() < target {
        si.step();
    }
}

fn max_wind(s: &LayeredState) -> f64 {
    s.u.iter().flatten().fold(0.0f64, |a, &u| a.max(u.abs()))
}

#[test]
fn conserves_mass_and_theta_at_a_long_step() {
    let model = model(true);
    let state = sheared(30.0).state(&model);
    let explicit = model.stable_dt(&state, 0.5);
    let start = model.diagnostics(&state);
    let mut si = stepper(model, state);
    run_days(&mut si, 10.0);
    let end = si.model().diagnostics(si.state());
    let rel = |a: f64, b: f64| b / a - 1.0;
    let stats = si.solver_stats();
    eprintln!(
        "step {:.1}× explicit, CG {:.1}/solve: mass {:+.1e} θ {:+.1e} energy {:+.1e}",
        si.dt() / explicit,
        stats.iterations as f64 / stats.solves as f64,
        rel(start.mass, end.mass),
        rel(start.theta, end.theta),
        rel(start.energy, end.energy),
    );
    assert!(
        si.dt() > 5.0 * explicit,
        "step only {:.1}× explicit",
        si.dt() / explicit
    );
    // The new p_s and Θ come from the divergence of the new winds, which has
    // no net source: exact whatever the solver tolerance.
    assert!(rel(start.mass, end.mass).abs() < 1e-13);
    assert!(rel(start.theta, end.theta).abs() < 1e-13);
    // Leapfrog with the RAW filter isn't energy-conserving. Measured: −3.1e-5
    // over 10 days with the default isothermal reference (−1.6e-5 with the
    // state's own mean profile; RK4: −3e-7).
    assert!(rel(start.energy, end.energy).abs() < 5e-5);
}

#[test]
fn more_layers_conserve_mass_and_theta_at_a_long_step() {
    // A stratified flow, sheared from 10 m/s at the bottom to 40 m/s at the
    // top, over the mountain for 10 days, on two, three (a thin bottom
    // layer) and five layers. It starts out of balance, so gravity waves
    // radiate and the RAW filter damps them: energy falls by 1.0e-4,
    // 1.6e-4 and 1.4e-4 (the balanced two-layer flow above: 3.1e-5).
    let (m, g) = mesh();
    for levels in [
        Levels::two(0.5),
        Levels::new(&[0.4, 0.85]),
        Levels::equal_mass(5),
    ] {
        let n = levels.count();
        let model = model_with(true, levels);
        let exner = model.layer_exner(1.0e5);
        // T = 288 K (p/p₀)^{RΓ/g} for Γ = 6.5 K/km, no colder than 215 K.
        let beta = GAS.gas_constant * 6.5e-3 / GRAVITY;
        let theta: Vec<Vec<f64>> = (0..n)
            .map(|k| {
                let t = (288.0 * fms_math::pow(exner[k], beta / GAS.kappa())).max(215.0);
                vec![t / exner[k]; m.n_cells()]
            })
            .collect();
        let u = (0..n)
            .map(|k| {
                let speed = 10.0 + 30.0 * (n - 1 - k) as f64 / (n - 1) as f64;
                m.edges()
                    .iter()
                    .zip(&g.edge_lengths)
                    .map(|(e, &l)| {
                        let psi = |v: u32| -RADIUS * speed * m.vertices()[v as usize].z;
                        -(psi(e.vertices[1]) - psi(e.vertices[0])) / (l * RADIUS)
                    })
                    .collect()
            })
            .collect();
        let state = LayeredState::new(vec![1.0e5; m.n_cells()], u, theta);
        let explicit = model.stable_dt(&state, 0.5);
        let start = model.diagnostics(&state);
        let mut si = stepper(model, state);
        run_days(&mut si, 10.0);
        si.check().expect("the layered run blew up");
        let end = si.model().diagnostics(si.state());
        let rel = |a: f64, b: f64| b / a - 1.0;
        eprintln!(
            "{n} layers, step {:.1}× explicit, waves {:.0?} m/s: mass {:+.1e} θ {:+.1e} energy {:+.1e}",
            si.dt() / explicit,
            si.wave_speeds(),
            rel(start.mass, end.mass),
            rel(start.theta, end.theta),
            rel(start.energy, end.energy),
        );
        assert!(rel(start.mass, end.mass).abs() < 1e-13);
        assert!(rel(start.theta, end.theta).abs() < 1e-13);
        assert!(rel(start.energy, end.energy).abs() < 5e-4);
    }
}

#[test]
fn isentropic_rest_over_a_mountain_stays_at_rest() {
    // Balanced exactly by the spatial scheme; the semi-implicit split moves
    // the linear part of that balance into the Helmholtz solve, so rest is
    // kept to the solver tolerance.
    for levels in [
        Levels::two(0.5),
        Levels::new(&[0.4, 0.85]),
        Levels::equal_mass(5),
    ] {
        isentropic_rest(levels);
    }
}

fn isentropic_rest(levels: Levels) {
    let n = levels.count();
    let model = model_with(true, levels);
    let theta0 = 300.0;
    let exner_ref = GAS.exner(1.0e5);
    let ps: Vec<f64> = model
        .surface_geopotential()
        .iter()
        .map(|phi| {
            1.0e5
                * fms_math::pow(
                    exner_ref - phi / (GAS.heat_capacity * theta0),
                    1.0 / GAS.kappa(),
                )
        })
        .collect();
    let (nc, ne) = (ps.len(), mesh().0.n_edges());
    let state = LayeredState::new(ps, vec![vec![0.0; ne]; n], vec![vec![theta0; nc]; n]);
    let mut si = stepper(model, state);
    run_days(&mut si, 5.0);
    let wind = max_wind(si.state());
    eprintln!("{n} layers, isentropic rest, day 5: max |u| {wind:.1e} m/s");
    assert!(wind < 1e-3, "{n} layers: spurious wind {wind} m/s");
}

#[test]
fn isothermal_rest_over_a_mountain_stays_nearly_at_rest() {
    // RK4 gives 0.105 m/s at this resolution (0.027 m/s at k = 32), for
    // two layers or more.
    for levels in [
        Levels::two(0.5),
        Levels::new(&[0.4, 0.85]),
        Levels::equal_mass(5),
    ] {
        isothermal_rest(levels);
    }
}

fn isothermal_rest(levels: Levels) {
    let n = levels.count();
    let model = model_with(true, levels);
    let t0 = 260.0;
    let ps: Vec<f64> = model
        .surface_geopotential()
        .iter()
        .map(|phi| 1.0e5 * fms_math::exp(-phi / (GAS.gas_constant * t0)))
        .collect();
    let mut theta = vec![vec![0.0; ps.len()]; n];
    for (i, &p) in ps.iter().enumerate() {
        let e = model.layer_exner(p);
        for k in 0..n {
            theta[k][i] = t0 / e[k];
        }
    }
    let ne = mesh().0.n_edges();
    let state = LayeredState::new(ps, vec![vec![0.0; ne]; n], theta);
    let mut si = stepper(model, state);
    run_days(&mut si, 5.0);
    let wind = max_wind(si.state());
    eprintln!("{n} layers, isothermal rest, day 5: max |u| {wind:.3} m/s");
    assert!(wind < 0.2, "{n} layers: spurious wind {wind} m/s");
}

/// Growth rate (per day) of a 0.01 K wavenumber-6 θ seed at 45° N, against
/// an unperturbed control, fitted from day 2 to day 8; optionally with the
/// default (strongest stable) ∇⁴ hyperdiffusion.
fn growth(shear: f64, diffuse: bool) -> f64 {
    let (m, g) = mesh();
    let flow = sheared(shear);
    let seeded = |seed: f64| {
        let mut model = model(false);
        if diffuse {
            let dt = LayeredSemiImplicit::stable_dt(&model, MAX_WIND, COURANT);
            let dt = DAY / (DAY / dt).ceil();
            let options = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
            model = model.with_hyperdiffusion(options);
        }
        let mut state = flow.state(&model);
        for k in 0..2 {
            for (i, &p) in m.cell_centers().iter().enumerate() {
                let (lat, lon) = to_lat_lon(p);
                let bump =
                    fms_math::exp(-((lat - 45f64.to_radians()) / 15f64.to_radians()).powi(2));
                state.theta_ps[k][i] += state.ps[i] * seed * bump * fms_math::cos(6.0 * lon);
            }
        }
        stepper(model, state)
    };
    let (mut control, mut perturbed) = (seeded(0.0), seeded(0.01));
    let theta_hat =
        |s: &LayeredState, i: usize| 0.5 * (s.theta_ps[0][i] + s.theta_ps[1][i]) / s.ps[i];
    let mut samples = Vec::new();
    for day in 0..=8 {
        run_days(&mut control, f64::from(day));
        run_days(&mut perturbed, f64::from(day));
        let (mut re, mut im, mut area) = (0.0, 0.0, 0.0);
        for (i, &p) in m.cell_centers().iter().enumerate() {
            let (lat, lon) = to_lat_lon(p);
            if (30f64.to_radians()..60f64.to_radians()).contains(&lat) {
                let d = g.cell_areas[i]
                    * (theta_hat(perturbed.state(), i) - theta_hat(control.state(), i));
                let (s, c) = fms_math::sin_cos(6.0 * lon);
                re += d * c;
                im += d * s;
                area += g.cell_areas[i];
            }
        }
        samples.push((
            f64::from(day),
            fms_math::ln(2.0 * (re * re + im * im).sqrt() / area),
        ));
    }
    let fit = &samples[2..];
    let n = fit.len() as f64;
    let (mt, my) = fit
        .iter()
        .fold((0.0, 0.0), |(a, b), (t, y)| (a + t / n, b + y / n));
    let (num, den) = fit.iter().fold((0.0, 0.0), |(a, b), (t, y)| {
        (a + (t - mt) * (y - my), b + (t - mt).powi(2))
    });
    num / den
}

#[test]
fn baroclinic_instability_at_a_long_step() {
    // Same experiment as the RK4 test in layered.rs. RK4 measures 0.033 and
    // 0.827 /day; the semi-implicit scheme −0.021 and 0.827 /day.
    let reference = model(false);
    let lat = 45f64.to_radians();
    let stable = growth(10.0, false);
    let (theory, _) = sheared(45.0).phillips_growth_rate(&reference, 6, lat);
    let unstable = growth(45.0, false);
    let ratio = unstable / (theory * DAY);
    eprintln!(
        "growth: {stable:.3} /day at 10 m/s, {unstable:.3} /day at 45 m/s (ratio {ratio:.2})"
    );
    assert!(stable < 0.05, "subcritical wave grew at {stable} /day");
    assert!((0.6..1.0).contains(&ratio), "growth ratio {ratio}");
}

#[test]
fn hyperdiffusion_leaves_baroclinic_growth_alone() {
    // Scale-selective: wavenumber 6 is far from the grid scale, so its growth
    // barely changes.
    let (inviscid, diffused) = (growth(45.0, false), growth(45.0, true));
    eprintln!("growth: {inviscid:.3} /day inviscid, {diffused:.3} /day with ∇⁴");
    assert!(
        (diffused / inviscid - 1.0).abs() < 0.03,
        "{inviscid} vs {diffused}"
    );
}

#[test]
fn results_do_not_depend_on_the_thread_count() {
    let simulate = || {
        let model = model(true);
        let state = sheared(30.0).state(&model);
        let mut si = stepper(model, state);
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
    let bits = |s: &LayeredState| -> Vec<u64> {
        s.ps.iter()
            .chain(s.u.iter().flatten())
            .chain(s.theta_ps.iter().flatten())
            .map(|v| v.to_bits())
            .collect()
    };
    assert_eq!(bits(&a), bits(&b));
}
