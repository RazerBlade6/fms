//! Restarts of the semi-implicit stepper (M2.13.5): a run saved part way and
//! resumed in a fresh stepper ends bit for bit where the straight run does,
//! dry, under Held–Suarez forcing, moist, coupled to the column physics with
//! its ocean, sea ice and land, and across thread counts; and a file that
//! doesn't fit is refused without touching the stepper.

use std::sync::OnceLock;

use fms_atmos::column_physics::{
    CLIMATE_LEVELS, OceanSpinUp, initial_state, initial_surface_temperature_by_latitude,
};
use fms_atmos::williamson::{Case, GRAVITY, OMEGA, RADIUS};
use fms_atmos::{
    CellField, ColumnPhysics, ForcingBudget, Gas, HeldSuarez, HeldSuarezParams,
    HyperdiffusionOptions, LayeredModel, LayeredSemiImplicit, LayeredSemiImplicitOptions,
    LayeredState, Levels, OceanParams, PhysicsParams, Restart, RestartError, ShearedFlow,
};
use fms_grid::coords::to_lat_lon;
use fms_grid::{LloydOptions, MeshGeometry, Operators, SphereMesh};
use fms_planet::{Planet, PlanetParams};
use glam::DVec3;

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
    let notes = c.load_restart(restart).expect("the restart loads");
    assert!(notes.is_empty(), "the same stepper adapted: {notes:?}");
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

// Coupled runs: the column physics, with its ocean, ice and land, at 642
// cells on the climate's four layers.

type Stepper = LayeredSemiImplicit<'static>;

fn earth() -> &'static Planet {
    static P: OnceLock<Planet> = OnceLock::new();
    P.get_or_init(Planet::earth)
}

/// Earth with a 20-day year, the shortest a planet may have, so that a year
/// ends within a short run.
fn short_year() -> &'static Planet {
    static P: OnceLock<Planet> = OnceLock::new();
    P.get_or_init(|| {
        Planet::new(PlanetParams {
            year_days: 20.0,
            ..PlanetParams::earth()
        })
        .0
    })
}

/// No land.
fn sea(_: DVec3) -> f64 {
    0.0
}

/// A continent within 45° of the prime meridian, with partly land cells
/// along its coasts (as in `climate.rs`).
fn continent(p: DVec3) -> f64 {
    let lon = fms_math::atan2(p.y, p.x).abs().to_degrees();
    ((45.0 - lon) / 10.0 + 0.5).clamp(0.0, 1.0)
}

/// No terrain: a cell's sub-grid height spread and mean height (m).
fn flat(_: DVec3) -> (f64, f64) {
    (0.0, 0.0)
}

/// A north–south ridge along the prime meridian, inside the continent:
/// its sub-grid spread and mean height (m).
fn ridge(p: DVec3) -> (f64, f64) {
    let lon = fms_math::atan2(p.y, p.x).to_degrees() / 8.0;
    let shape = fms_math::exp(-0.5 * lon * lon);
    (800.0 * shape, 2500.0 * shape)
}

/// How a coupled run is built: everything a restarted run is built from
/// again.
#[derive(Clone, Copy)]
struct World {
    planet: fn() -> &'static Planet,
    params: PhysicsParams,
    land: fn(DVec3) -> f64,
    relief: fn(DVec3) -> (f64, f64),
    /// The sea-level temperature of a start warmer at the equator than at
    /// the poles (K); `None` for the single column's everywhere.
    sea_level: Option<f64>,
    /// Sea temperatures held fixed.
    prescribed: bool,
    /// The design wind, which sets the step (m s⁻¹).
    max_wind: f64,
    /// Steps between the ocean's transports, instead of the radiation's.
    ocean_steps: Option<f64>,
}

impl World {
    /// An aquaplanet over the default ocean, with currents and ice.
    fn aquaplanet() -> Self {
        Self {
            planet: earth,
            params: PhysicsParams::default(),
            land: sea,
            relief: flat,
            sea_level: None,
            prescribed: false,
            max_wind: 120.0,
            ocean_steps: None,
        }
    }

    fn build(&self) -> Stepper {
        let (m, g) = coarse_mesh();
        let planet = (self.planet)();
        let centers = m.cell_centers();
        let model = LayeredModel::new(
            Operators::new(m, g, RADIUS),
            GRAVITY,
            GAS,
            Levels::new(&CLIMATE_LEVELS),
            m.vertices()
                .iter()
                .map(|p| 2.0 * planet.rotation_rate() * p.z)
                .collect(),
            centers
                .iter()
                .map(|&p| GRAVITY * (self.relief)(p).1)
                .collect(),
        );
        let (state, ts) = initial_state(&model, planet, &self.params);
        let dt = step_for(&model, self.max_wind);
        let diffusion = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
        let surface = match self.sea_level {
            Some(t) => initial_surface_temperature_by_latitude(&model, planet, &self.params, t),
            None => vec![ts; m.n_cells()],
        };
        let mut params = self.params.with_step(dt, planet.solar_day());
        if let Some(steps) = self.ocean_steps {
            params.transport_interval = steps * dt;
        }
        let physics = ColumnPhysics::with_land(
            &model,
            planet,
            params,
            0.0,
            &state,
            surface,
            centers.iter().map(|&p| (self.land)(p)).collect(),
        )
        .with_relief(
            centers.iter().map(|&p| (self.relief)(p).0).collect(),
            &centers
                .iter()
                .map(|&p| (self.relief)(p).1)
                .collect::<Vec<_>>(),
        );
        let physics = if self.prescribed {
            physics.with_prescribed_surface()
        } else {
            physics
        };
        let model = model
            .with_hyperdiffusion(diffusion)
            .with_dissipative_heating()
            .with_forcing(Box::new(physics));
        stepper(model, state, dt)
    }
}

/// The forcing's budget, as bits.
fn budget_bits(si: &Stepper) -> [u64; 10] {
    let ForcingBudget {
        absorbed,
        outgoing,
        stored,
        evaporated,
        precipitated,
        dissipated,
        soil_water,
        snow,
        runoff,
        ocean_exchange,
    } = si.model().forcing().unwrap().budget().unwrap();
    [
        absorbed,
        outgoing,
        stored,
        evaporated,
        precipitated,
        dissipated,
        soil_water,
        snow,
        runoff,
        ocean_exchange,
    ]
    .map(f64::to_bits)
}

fn fields(si: &Stepper) -> Vec<CellField> {
    si.model().forcing().unwrap().cell_fields()
}

/// The forcing's outputs, its budget and every cell field, are the same
/// bit for bit in both steppers.
fn assert_same_outputs(a: &Stepper, b: &Stepper, when: &str) {
    assert_eq!(budget_bits(a), budget_bits(b), "the budget differs {when}");
    let (fa, fb) = (fields(a), fields(b));
    assert_eq!(fa.len(), fb.len());
    for (x, y) in fa.iter().zip(&fb) {
        assert_eq!(x.name, y.name);
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert!(
            bits(&x.values) == bits(&y.values),
            "the cell field `{}` differs {when}",
            x.name
        );
    }
}

/// What [`resume_coupled`] leaves.
struct Coupled {
    /// The file saved part way.
    file: Restart,
    /// The cell fields of the saved run when it was saved.
    saved_fields: Vec<CellField>,
    /// The straight run, at its end.
    straight: Stepper,
}

/// Run `build`'s stepper straight to `steps`; and to `saved_at`, saved,
/// through bytes, loaded into a fresh stepper and run on to `steps`. Right
/// after loading, the fresh stepper saves the file again and reports what
/// the saved one did; at the end it is the straight run bit for bit, its
/// outputs included.
fn resume_coupled(build: impl Fn() -> Stepper, saved_at: u64, steps: u64) -> Coupled {
    let mut a = build();
    run_to(&mut a, steps);
    let mut b = build();
    run_to(&mut b, saved_at);
    let file = through_bytes(&b.save_restart());
    let mut c = build();
    let notes = c.load_restart(&file).expect("the restart loads");
    assert!(notes.is_empty(), "the same stepper adapted: {notes:?}");
    assert_eq!(c.steps(), saved_at);
    assert!(
        c.save_restart().to_bytes() == file.to_bytes(),
        "loaded and saved at once, the file changed"
    );
    assert_same_outputs(&b, &c, "right after loading");
    let saved_fields = fields(&b);
    drop(b);
    run_to(&mut c, steps);
    assert_identical(&a.save_restart(), &c.save_restart());
    assert_same_outputs(&a, &c, "at the end");
    Coupled {
        file,
        saved_fields,
        straight: a,
    }
}

/// The array `name` of `r`, whatever its length.
fn array<'r>(r: &'r Restart, name: &str) -> &'r [f64] {
    r.arrays()
        .find(|(n, _)| *n == name)
        .unwrap_or_else(|| panic!("no array `{name}`"))
        .1
}

/// `physics.clock`: time, end of year, steps since radiation and since the
/// ocean's transport, years ended, jumps done.
fn clock(r: &Restart) -> [f64; 6] {
    r.get_array("physics.clock").unwrap()
}

fn field<'f>(fields: &'f [CellField], name: &str) -> &'f [f64] {
    &fields.iter().find(|f| f.name == name).unwrap().values
}

#[test]
fn an_aquaplanet_with_sea_ice_resumes_bit_for_bit() {
    // Cold poles: ice from the start. Saved on a step that is neither a
    // radiation step nor one of the ocean's transports.
    let world = World {
        sea_level: Some(280.0),
        ..World::aquaplanet()
    };
    let run = resume_coupled(|| world.build(), 13, 30);
    let [_, _, since_radiation, since_transport, ..] = clock(&run.file);
    assert!(
        since_radiation > 0.0 && since_transport > 0.0,
        "saved between calls: {since_radiation} {since_transport}"
    );
    assert!(array(&run.file, "ocean.wind.samples")[0] > 0.0);
    assert!(
        array(&run.file, "ice.fraction").iter().any(|&f| f > 0.0),
        "no sea ice at the save point"
    );
    assert!(
        array(&run.file, "ocean.mass_flux.ekman")
            .iter()
            .any(|&f| f != 0.0),
        "the currents never ran"
    );
    assert!(!run.file.contains("physics.year_start.mixed"));
}

#[test]
fn a_world_with_land_resumes_bit_for_bit_across_a_jump() {
    // A continent with a ridge under a 20-day year, its first two years
    // ending in a jump: saved half way through the first, which is being
    // measured for its jump, after the gyres have been solved, with snow on
    // the cold land; the jump comes after the resume.
    let world = World {
        planet: short_year,
        params: PhysicsParams {
            spin_up: OceanSpinUp {
                jumps: 2,
                first_year: 1,
                every: 1,
                ..OceanSpinUp::default()
            },
            ..PhysicsParams::default()
        },
        land: continent,
        relief: ridge,
        sea_level: Some(282.0),
        ..World::aquaplanet()
    };
    let per_day = (DAY / world.build().dt()).round() as u64;
    let (saved_at, steps) = (10 * per_day + 1, 20 * per_day + 3);
    let run = resume_coupled(|| world.build(), saved_at, steps);
    let file = &run.file;
    assert!(
        file.contains("physics.year_start.mixed"),
        "no year measured"
    );
    assert_eq!(clock(file)[5], 0.0, "jumped before the save");
    assert_eq!(clock(&run.straight.save_restart())[5], 1.0, "no jump after");
    let [_, solves] = file.get_array("ocean.gyres.clock").unwrap();
    assert!(solves >= 1.0, "the gyres were never solved");
    assert!(array(file, "ocean.gyres.psi").iter().any(|&p| p != 0.0));
    let land = array(file, "physics.land_fraction");
    let on_land = |name: &str| -> Vec<f64> {
        array(file, name)
            .iter()
            .zip(land)
            .filter(|(_, f)| **f > 0.0)
            .map(|(x, _)| *x)
            .collect()
    };
    assert!(on_land("physics.snow").iter().any(|&s| s > 0.0), "no snow");
    let water = on_land("physics.soil_water");
    let spread =
        water.iter().fold(0.0f64, |m, &w| m.max(w)) - water.iter().fold(f64::MAX, |m, &w| m.min(w));
    assert!(spread > 1.0, "the soil water never changed: {spread}");
    assert!(
        field(&run.saved_fields, "atm.orographic_precipitation")
            .iter()
            .any(|&r| r > 0.0),
        "no rain from the hills"
    );
}

#[test]
fn a_coupled_run_saved_on_one_thread_resumes_on_four() {
    // A continent with a ridge: the gyres are first solved after the resume.
    let world = World {
        land: continent,
        relief: ridge,
        ..World::aquaplanet()
    };
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let four = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    let (saved_at, steps) = (13, 30);
    let mut a = world.build();
    run_to(&mut a, steps);
    let file = one.install(|| {
        let mut b = world.build();
        run_to(&mut b, saved_at);
        through_bytes(&b.save_restart())
    });
    let c = four.install(|| {
        let mut c = world.build();
        assert!(c.load_restart(&file).unwrap().is_empty());
        run_to(&mut c, steps);
        c
    });
    assert_identical(&a.save_restart(), &c.save_restart());
    assert_same_outputs(&a, &c, "at the end");
    let [_, solves] = file.get_array("ocean.gyres.clock").unwrap();
    assert_eq!(solves, 0.0);
    assert!(a.model().forcing().unwrap().cell_fields().len() > 10);
    let [_, solves] = a.save_restart().get_array("ocean.gyres.clock").unwrap();
    assert!(solves >= 1.0);
}

#[test]
fn a_run_over_prescribed_sea_temperatures_resumes_bit_for_bit() {
    let world = World {
        prescribed: true,
        ..World::aquaplanet()
    };
    let run = resume_coupled(|| world.build(), 13, 30);
    assert!(array(&run.file, "physics.budget")[0] != 0.0, "no uptake");
}

/// Loading `from` into `si` fails and leaves it as it was, outputs
/// included.
fn refused_coupled(si: &mut Stepper, from: &Restart) -> RestartError {
    let before = (si.save_restart().to_bytes(), budget_bits(si), fields(si));
    let err = refused(si, from);
    assert!(
        (si.save_restart().to_bytes(), budget_bits(si), fields(si)) == before,
        "a refused load changed the stepper ({err})"
    );
    err
}

#[test]
fn coupled_restarts_that_do_not_fit_are_refused() {
    let world = World::aquaplanet();
    let good = {
        let mut b = world.build();
        run_to(&mut b, 3);
        b.save_restart()
    };
    let mut si = world.build();
    si.step();

    // Other land.
    let file = World {
        land: continent,
        ..world
    }
    .build()
    .save_restart();
    let err = refused_coupled(&mut si, &file);
    assert!(
        matches!(&err, RestartError::Mismatch(m) if m.contains("land")),
        "{err}"
    );
    eprintln!("{err}");

    // The dynamics alone.
    let file = edited(&good, |name, v| {
        name.starts_with("dynamics.").then(|| v.to_vec())
    });
    let err = refused_coupled(&mut si, &file);
    assert!(matches!(err, RestartError::Mismatch(_)), "{err}");
    eprintln!("{err}");

    // Ice in the file and not in the run, and the reverse: warm starts.
    let no_ice = World {
        params: PhysicsParams {
            ocean: OceanParams {
                ice: None,
                ..OceanParams::default()
            },
            ..PhysicsParams::default()
        },
        ..world
    };
    let notes = no_ice.build().load_restart(&good).unwrap();
    assert!(
        notes.len() == 1 && notes[0].contains("sea ice is off"),
        "{notes:?}"
    );
    let notes = world
        .build()
        .load_restart(&no_ice.build().save_restart())
        .unwrap();
    assert!(
        notes.len() == 1 && notes[0].contains("sea ice is on"),
        "{notes:?}"
    );

    // An array of the physics, the ocean or the ice missing, or one value
    // short. The ice's last array is the last checked: the dynamics and
    // the physics are untouched when it fails.
    for name in [
        "dynamics.dt",
        "physics.intervals",
        "physics.clock",
        "physics.soil_water",
        "physics.surface",
        "physics.moist",
        "ocean.components",
        "ocean.mixed",
        "ocean.mass_flux.gyre",
        "ocean.gyres.stats",
        "ice.thickness",
        "ice.limited",
    ] {
        let file = edited(&good, |n, v| (n != name).then(|| v.to_vec()));
        let err = refused_coupled(&mut si, &file);
        assert!(
            matches!(&err, RestartError::Missing { name: n } if n == name),
            "{name}: {err}"
        );
        let file = edited(&good, |n, v| {
            let mut v = v.to_vec();
            if n == name {
                v.pop();
            }
            Some(v)
        });
        let err = refused_coupled(&mut si, &file);
        assert!(
            matches!(&err, RestartError::WrongLength { name: n, .. } if n == name),
            "{name}: {err}"
        );
    }

    // And the good file loads, ignoring arrays it doesn't know.
    let mut file = good.clone();
    file.put("cli.something", &[1.0]);
    assert!(si.load_restart(&file).unwrap().is_empty());
    assert!(si.save_restart().to_bytes() == good.to_bytes());
}

// Warm starts: a file loaded into a run built with other settings is
// adapted to them, and says how.

/// The stepper is healthy at the end of each of `days` more days.
fn runs_on(si: &mut Stepper, days: u64) {
    let per_day = (DAY / si.dt()).round() as u64;
    for _ in 0..days {
        run_to(si, si.steps() + per_day);
        si.check().expect("the warm start stays healthy");
    }
}

#[test]
fn a_changed_step_rebases_the_step_count() {
    let world = World::aquaplanet();
    let mut b = world.build();
    run_to(&mut b, 40);
    let file = through_bytes(&b.save_restart());
    let slower = World {
        max_wind: 240.0,
        ..world
    };
    let mut si = slower.build();
    assert!(si.dt() < b.dt());
    let notes = si.load_restart(&file).unwrap();
    eprintln!("{notes:#?}");
    let expected = (40.0 * b.dt() / si.dt()).round() as u64;
    assert_eq!(si.steps(), expected);
    assert!(
        notes[0].contains("time step changed")
            && notes[0].contains(&format!("counted as {expected} ")),
        "{notes:?}"
    );
    // The radiation and the ocean's transport follow the step.
    assert!(
        notes.iter().any(|n| n.contains("intervals changed")),
        "{notes:?}"
    );
    runs_on(&mut si, 3);
}

#[test]
fn changed_physics_intervals_are_noted() {
    let world = World::aquaplanet();
    let mut b = world.build();
    run_to(&mut b, 13);
    let file = through_bytes(&b.save_restart());
    let every_step = World {
        ocean_steps: Some(1.0),
        ..world
    };
    let mut si = every_step.build();
    let notes = si.load_restart(&file).unwrap();
    eprintln!("{notes:#?}");
    assert!(
        notes.len() == 1 && notes[0].contains("intervals changed"),
        "{notes:?}"
    );
    assert_eq!(si.steps(), 13);
    runs_on(&mut si, 2);
}

/// A 20-day year with one ocean jump, ending year `first_year`.
fn one_jump(first_year: u32) -> World {
    World {
        planet: short_year,
        params: PhysicsParams {
            spin_up: OceanSpinUp {
                jumps: 1,
                first_year,
                every: 1,
                ..OceanSpinUp::default()
            },
            ..PhysicsParams::default()
        },
        ..World::aquaplanet()
    }
}

/// Jumps done, and the heat they added.
fn jumps(si: &Stepper) -> (f64, f64) {
    let r = si.save_restart();
    (clock(&r)[5], array(&r, "physics.budget")[9])
}

#[test]
fn a_changed_jump_schedule_is_reconciled() {
    // Saved with the jump at the end of year 3: half way through year 1,
    // which isn't measured, and at the start of year 2. Then with it at the
    // end of year 1: half way through year 1, being measured.
    let late = one_jump(3);
    let mut b = late.build();
    let per_day = (DAY / b.dt()).round() as u64;
    run_to(&mut b, 10 * per_day);
    let mid_year = through_bytes(&b.save_restart());
    run_to(&mut b, 20 * per_day);
    let year_start = through_bytes(&b.save_restart());
    let [time, next_year, ..] = clock(&year_start);
    assert!(
        (next_year - time - 20.0 * DAY).abs() < 1e-6,
        "not at a year's start: {time} {next_year}"
    );
    assert!(!mid_year.contains("physics.year_start.mixed"));
    assert!(!year_start.contains("physics.year_start.mixed"));
    let mut b = one_jump(1).build();
    run_to(&mut b, 10 * per_day);
    let measured = through_bytes(&b.save_restart());
    assert!(measured.contains("physics.year_start.mixed"));

    // No jump wanted where the file measures one: dropped, and the year
    // ends without one.
    let mut si = late.build();
    let notes = si.load_restart(&measured).unwrap();
    eprintln!("{notes:#?}");
    assert!(
        notes.len() == 1 && notes[0].contains("dropped"),
        "{notes:?}"
    );
    run_to(&mut si, 21 * per_day);
    assert_eq!(jumps(&si), (0.0, 0.0));

    // A jump wanted at the end of year 2, loaded at its start: measured from
    // the loaded ocean, and it jumps.
    let mut si = one_jump(2).build();
    let notes = si.load_restart(&year_start).unwrap();
    eprintln!("{notes:#?}");
    assert!(
        notes.len() == 1 && notes[0].contains("measured from here"),
        "{notes:?}"
    );
    run_to(&mut si, 41 * per_day);
    let (done, heat) = jumps(&si);
    assert!(done == 1.0 && heat != 0.0, "no jump: {done} {heat}");

    // A jump wanted at the end of year 1, loaded half way through it: that
    // jump is skipped.
    let mut si = one_jump(1).build();
    let notes = si.load_restart(&mid_year).unwrap();
    eprintln!("{notes:#?}");
    assert!(
        notes.len() == 1 && notes[0].contains("skipped"),
        "{notes:?}"
    );
    run_to(&mut si, 21 * per_day);
    assert_eq!(jumps(&si), (0.0, 0.0));
    si.check().expect("healthy");
}

#[test]
fn years_a_saved_run_did_not_count_are_counted_on_loading() {
    // Over prescribed sea temperatures years aren't counted: saved a day
    // into year 2, its year still ends at day 20. Loaded where a jump ends
    // year 2, year 1 is counted as ended, once, and year 2's jump skipped.
    let prescribed = World {
        prescribed: true,
        ..one_jump(1)
    };
    let mut b = prescribed.build();
    let per_day = (DAY / b.dt()).round() as u64;
    run_to(&mut b, 21 * per_day);
    let file = through_bytes(&b.save_restart());
    let [time, next_year, _, _, years_ended, _] = clock(&file);
    assert!(time > next_year && years_ended == 0.0, "{time} {next_year}");
    // Measured from the start, for the jump at the end of year 1 that a
    // prescribed surface never takes.
    assert!(file.contains("physics.year_start.mixed"));
    let mut si = one_jump(2).build();
    let notes = si.load_restart(&file).unwrap();
    eprintln!("{notes:#?}");
    assert!(
        notes.len() == 2
            && notes[0].contains("1 more are counted as ended, so year 2 is in progress, and")
            && notes[1].contains("skipped"),
        "{notes:?}"
    );
    run_to(&mut si, 23 * per_day);
    let [_, next_year, _, _, years_ended, jumps_done] = clock(&si.save_restart());
    assert_eq!((years_ended, jumps_done), (1.0, 0.0));
    assert!((next_year - 40.0 * DAY).abs() < 1e-6, "{next_year}");
}

#[test]
fn ocean_components_switched_on_or_off_start_warm_in_a_coupled_run() {
    // A continent, so the gyres are solved; the file has every component.
    let world = World {
        land: continent,
        ..World::aquaplanet()
    };
    let mut b = world.build();
    run_to(&mut b, 13);
    let full = through_bytes(&b.save_restart());
    let without = |ocean: OceanParams| World {
        params: PhysicsParams {
            ocean,
            ..PhysicsParams::default()
        },
        ..world
    };
    let base = OceanParams::default();
    for other in [
        OceanParams {
            ekman: None,
            ..base
        },
        OceanParams {
            eddies: None,
            ..base
        },
        OceanParams {
            gyres: None,
            ..base
        },
        OceanParams { ice: None, ..base },
    ] {
        let mut si = without(other).build();
        let notes = si.load_restart(&full).unwrap();
        eprintln!("{notes:?}");
        assert!(notes.len() == 1 && notes[0].contains(" off "), "{notes:?}");
        runs_on(&mut si, 1);
        let mut si = world.build();
        let mut b = without(other).build();
        run_to(&mut b, 13);
        let notes = si.load_restart(&b.save_restart()).unwrap();
        assert!(notes.len() == 1 && notes[0].contains(" on "), "{notes:?}");
        runs_on(&mut si, 1);
    }
}
