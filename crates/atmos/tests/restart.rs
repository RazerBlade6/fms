//! Restarts of the semi-implicit stepper (M2.13.5): a run saved part way and
//! resumed in a fresh stepper ends bit for bit where the straight run does,
//! dry, under Held–Suarez forcing, moist, and across thread counts; and a
//! file that doesn't fit is refused without touching the stepper.

use std::sync::OnceLock;

use fms_atmos::williamson::{Case, GRAVITY, OMEGA, RADIUS};
use fms_atmos::{
    Gas, HeldSuarez, HeldSuarezParams, HyperdiffusionOptions, LayeredModel, LayeredSemiImplicit,
    LayeredSemiImplicitOptions, LayeredState, Levels, Restart, RestartError, ShearedFlow,
};
use fms_grid::coords::to_lat_lon;
use fms_grid::{LloydOptions, MeshGeometry, Operators, SphereMesh};

const GAS: Gas = Gas::EARTH_AIR;
const DAY: f64 = 86_400.0;
const MAX_WIND: f64 = 60.0;
const COURANT: f64 = 0.6;
/// Steps of the straight run, and of the run that is saved.
const STEPS: u64 = 30;
const SAVED_AT: u64 = 13;

type Mesh = (SphereMesh, MeshGeometry);

fn mesh_of(cell: &'static OnceLock<Mesh>, frequency: u32) -> &'static Mesh {
    cell.get_or_init(|| {
        let (m, _) = SphereMesh::centroidal(frequency, LloydOptions::default());
        let g = MeshGeometry::new(&m);
        (m, g)
    })
}

fn mesh() -> &'static Mesh {
    static M: OnceLock<Mesh> = OnceLock::new();
    mesh_of(&M, 16)
}

fn coarse_mesh() -> &'static Mesh {
    static M: OnceLock<Mesh> = OnceLock::new();
    mesh_of(&M, 8)
}

fn model_on(mesh: &'static Mesh, mountain: bool, levels: Levels) -> LayeredModel<'static> {
    let (m, g) = mesh;
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

/// The step for winds up to `max_wind`, a whole number per day.
fn step_for(model: &LayeredModel, max_wind: f64) -> f64 {
    let dt = LayeredSemiImplicit::stable_dt(model, max_wind, COURANT);
    DAY / (DAY / dt).ceil()
}

fn diffused(model: LayeredModel<'static>, dt: f64) -> LayeredModel<'static> {
    let options = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
    model.with_hyperdiffusion(options)
}

fn stepper(
    model: LayeredModel<'static>,
    state: LayeredState,
    dt: f64,
) -> LayeredSemiImplicit<'static> {
    LayeredSemiImplicit::new(model, state, dt, LayeredSemiImplicitOptions::default(), 0.5)
}

/// A stratified atmosphere at 1000 hPa, its wind sheared from 10 m/s at the
/// bottom to 40 m/s at the top (as in `layered_si.rs`).
fn stratified(model: &LayeredModel) -> LayeredState {
    let (m, g) = (model.operators().mesh(), model.operators().geometry());
    let n = model.layers();
    let exner = model.layer_exner(1.0e5);
    let beta = GAS.gas_constant * 6.5e-3 / GRAVITY;
    let theta = (0..n)
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
    LayeredState::new(vec![1.0e5; m.n_cells()], u, theta)
}

/// Case 1: four layers over the mountain, sheared, with hyperdiffusion and
/// dissipative heating.
fn dry() -> LayeredSemiImplicit<'static> {
    let model = model_on(mesh(), true, Levels::new(&[0.25, 0.5, 0.85]));
    let dt = step_for(&model, MAX_WIND);
    let model = diffused(model.with_dissipative_heating(), dt);
    let state = stratified(&model);
    stepper(model, state, dt)
}

/// Case 2: Held–Suarez, with hyperdiffusion on three layers, or alone (the
/// stepper's forcing-only path) on two.
fn held_suarez(levels: Levels, diffuse: bool) -> LayeredSemiImplicit<'static> {
    let model = model_on(mesh(), false, levels);
    let dt = step_for(&model, 80.0);
    let hs = HeldSuarez::new(&model, HeldSuarezParams::default());
    let state = hs.initial_state(&model, 0.1);
    let model = if diffuse { diffused(model, dt) } else { model };
    stepper(model.with_forcing(Box::new(hs)), state, dt)
}

/// Case 3: the moisture test's storm carrying a square plume of humidity.
fn moist() -> LayeredSemiImplicit<'static> {
    let model = model_on(mesh(), false, Levels::two(0.5));
    let dt = step_for(&model, 120.0);
    let model = diffused(model, dt);
    let flow = ShearedFlow {
        speed: [40.0, 10.0],
        stability: 40.0,
        theta_equator: 315.0,
        ps_equator: 1.0e5,
        omega: OMEGA,
    };
    let mut state = flow.state(&model);
    let centers = mesh().0.cell_centers();
    for k in 0..2 {
        for (i, &p) in centers.iter().enumerate() {
            let (lat, lon) = to_lat_lon(p);
            let band =
                |c: f64| fms_math::exp(-((lat - c.to_radians()) / 15f64.to_radians()).powi(2));
            let waves: f64 = (4..=8)
                .map(|m| fms_math::cos(f64::from(m) * lon + 0.7 * f64::from(m)))
                .sum();
            state.theta_ps[k][i] += state.ps[i] * 0.5 * (band(45.0) + band(-45.0)) * waves;
        }
    }
    let q: Vec<f64> = centers
        .iter()
        .map(|&p| {
            let (lat, lon) = to_lat_lon(p);
            let inside = (lat.to_degrees() - 45.0).abs() < 12.0 && lon.to_degrees().abs() < 30.0;
            if inside { 0.015 } else { 0.0 }
        })
        .collect();
    let state = state.with_humidity(vec![q.iter().map(|q| 0.1 * q).collect(), q]);
    stepper(model, state, dt)
}

fn run_to(si: &mut LayeredSemiImplicit, steps: u64) {
    while si.steps() < steps {
        si.step();
    }
}

/// Through the file format and back.
fn through_bytes(r: &Restart) -> Restart {
    Restart::read_from(&mut r.to_bytes().as_slice()).expect("a restart reads back")
}

/// Saved at `SAVED_AT` by one stepper.
fn saved(build: impl Fn() -> LayeredSemiImplicit<'static>) -> Restart {
    let mut b = build();
    run_to(&mut b, SAVED_AT);
    through_bytes(&b.save_restart())
}

/// A fresh stepper with `restart` loaded, run to `STEPS`.
fn resumed(
    build: impl Fn() -> LayeredSemiImplicit<'static>,
    restart: &Restart,
) -> LayeredSemiImplicit<'static> {
    let mut c = build();
    c.load_restart(restart).expect("the restart loads");
    assert_eq!(c.steps(), SAVED_AT);
    run_to(&mut c, STEPS);
    c
}

/// The two restarts are the same bit for bit, naming the first array that
/// differs if not.
fn assert_identical(straight: &Restart, resumed: &Restart) {
    assert_eq!(straight.step(), resumed.step());
    let names = |r: &Restart| r.names().map(str::to_owned).collect::<Vec<_>>();
    assert_eq!(names(straight), names(resumed));
    for ((name, a), (_, b)) in straight.arrays().zip(resumed.arrays()) {
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert!(bits(a) == bits(b), "`{name}` differs after resuming");
    }
    assert!(straight.to_bytes() == resumed.to_bytes());
}

/// Run straight to `STEPS`, and saved at `SAVED_AT` then resumed in a fresh
/// stepper; both end bit for bit the same. Returns the two steppers.
fn resume(
    build: impl Fn() -> LayeredSemiImplicit<'static>,
) -> (LayeredSemiImplicit<'static>, LayeredSemiImplicit<'static>) {
    let mut a = build();
    run_to(&mut a, STEPS);
    let c = resumed(&build, &saved(&build));
    assert_identical(&a.save_restart(), &c.save_restart());
    (a, c)
}

#[test]
fn a_dry_run_resumes_bit_for_bit() {
    let (a, _) = resume(dry);
    assert!(a.energy_budget().diffusion.kinetic < 0.0);
    assert!(
        a.energy_budget().diffusion.thermal > 0.0,
        "dissipative heating is on"
    );
}

#[test]
fn a_held_suarez_run_resumes_bit_for_bit() {
    let (a, _) = resume(|| held_suarez(Levels::new(&[0.4, 0.85]), true));
    assert!(a.energy_budget().forcing.thermal != 0.0);
    let (a, _) = resume(|| held_suarez(Levels::two(0.5), false));
    assert!(a.energy_budget().forcing.kinetic < 0.0);
    assert_eq!(a.energy_budget().diffusion.total(), 0.0);
}

#[test]
fn a_moist_run_resumes_bit_for_bit() {
    let (a, c) = resume(moist);
    let (wa, wc) = (a.water_budget(), c.water_budget());
    assert!(wa.filter != 0.0 && wa.max_substeps > 0, "{wa:?}");
    let bits =
        |w: fms_atmos::WaterBudget| [w.sources, w.filter, w.fixer, w.borrowed].map(f64::to_bits);
    assert_eq!(bits(wa), bits(wc));
    assert_eq!(wa.max_substeps, wc.max_substeps);
    assert_eq!(a.water().to_bits(), c.water().to_bits());
}

#[test]
fn a_run_saved_on_one_thread_resumes_on_four() {
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let four = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    let mut a = dry();
    run_to(&mut a, STEPS);
    let restart = one.install(|| saved(dry));
    let c = four.install(|| resumed(dry, &restart).save_restart());
    assert_identical(&a.save_restart(), &c);
}

/// Loading `from` into `si` fails, and leaves `si` as it was.
fn refused(si: &mut LayeredSemiImplicit, from: &Restart) -> RestartError {
    let before = si.save_restart().to_bytes();
    let err = si
        .load_restart(from)
        .expect_err("the restart should be refused");
    assert!(
        si.save_restart().to_bytes() == before,
        "a refused load changed the stepper ({err})"
    );
    err
}

/// A copy of `r` with `edit` applied to each array (`None` drops it).
fn edited(r: &Restart, edit: impl Fn(&str, &[f64]) -> Option<Vec<f64>>) -> Restart {
    let mut out = Restart::new(r.grid(), r.step());
    for (name, values) in r.arrays() {
        if let Some(values) = edit(name, values) {
            out.put(name, &values);
        }
    }
    out
}

#[test]
fn restarts_that_do_not_fit_are_refused() {
    let levels = || Levels::new(&[0.25, 0.5, 0.85]);
    let build = || stepper_at_rest(mesh(), levels(), false);
    let good = build().save_restart();
    let mut si = build();
    si.step();

    // Another mesh.
    let other = stepper_at_rest(coarse_mesh(), levels(), false).save_restart();
    let err = refused(&mut si, &other);
    assert!(matches!(err, RestartError::DifferentGrid { .. }), "{err}");
    eprintln!("{err}");

    // Other levels: fewer, or as many in other places.
    for other in [Levels::new(&[0.4, 0.85]), Levels::new(&[0.3, 0.5, 0.85])] {
        let file = stepper_at_rest(mesh(), other, false).save_restart();
        let err = refused(&mut si, &file);
        assert!(matches!(err, RestartError::Mismatch(_)), "{err}");
        eprintln!("{err}");
    }

    // Moist into dry, and dry into moist.
    let moist_file = stepper_at_rest(mesh(), levels(), true).save_restart();
    let err = refused(&mut si, &moist_file);
    assert!(matches!(err, RestartError::Mismatch(_)), "{err}");
    eprintln!("{err}");
    let err = refused(&mut stepper_at_rest(mesh(), levels(), true), &good);
    assert!(matches!(err, RestartError::Mismatch(_)), "{err}");
    eprintln!("{err}");

    // An array missing, or the wrong length; the last array missing shows
    // nothing was copied before the check failed.
    for dropped in [
        "dynamics.modal.0",
        "dynamics.previous.u.3",
        "dynamics.solver_stats",
    ] {
        let file = edited(&good, |name, v| (name != dropped).then(|| v.to_vec()));
        let err = refused(&mut si, &file);
        assert!(
            matches!(&err, RestartError::Missing { name } if name == dropped),
            "{err}"
        );
    }
    let file = edited(&good, |name, v| {
        let mut v = v.to_vec();
        if name == "dynamics.previous.theta_ps.1" {
            v.pop();
        }
        Some(v)
    });
    let err = refused(&mut si, &file);
    assert!(matches!(err, RestartError::WrongLength { .. }), "{err}");
    eprintln!("{err}");

    // And the good file loads.
    si.load_restart(&good).unwrap();
    assert!(si.save_restart().to_bytes() == good.to_bytes());
}

/// A quick stepper on `mesh`: isothermal air at rest, dry or moist.
fn stepper_at_rest(
    mesh: &'static Mesh,
    levels: Levels,
    moist: bool,
) -> LayeredSemiImplicit<'static> {
    let n = levels.count();
    let model = model_on(mesh, false, levels);
    let (nc, ne) = (mesh.0.n_cells(), mesh.0.n_edges());
    let exner = model.layer_exner(1.0e5);
    let theta = (0..n).map(|k| vec![260.0 / exner[k]; nc]).collect();
    let mut state = LayeredState::new(vec![1.0e5; nc], vec![vec![0.0; ne]; n], theta);
    if moist {
        state = state.with_humidity(vec![vec![1e-3; nc]; n]);
    }
    let dt = step_for(&model, MAX_WIND);
    stepper(model, state, dt)
}
