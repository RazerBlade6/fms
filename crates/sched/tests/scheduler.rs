//! Scheduler behaviour: rates, stage planning, determinism, access checks.

use std::sync::{Arc, Mutex};

use fms_fields::{FieldId, FieldRegistry, FieldSpec, GridSizes, Location, MeshCounts};
use fms_sched::{
    Access, Calendar, Rate, SchedError, Scheduler, StepContext, System, SystemDesc, Tier, Tiers,
};

fn registry() -> FieldRegistry {
    let counts = MeshCounts {
        cells: 64,
        edges: 186,
        vertices: 124,
    };
    FieldRegistry::new(GridSizes {
        terrain: counts,
        climate: counts,
    })
}

fn field(r: &mut FieldRegistry, name: &'static str, value: f64) -> FieldId<f64> {
    r.register(FieldSpec::new(name, "", Location::CLIMATE_CELLS), value)
        .unwrap()
}

/// A closure-backed system for tests.
struct TestSystem<F> {
    desc: SystemDesc,
    body: F,
}

impl<F: FnMut(&StepContext<'_>) + Send> System for TestSystem<F> {
    fn describe(&self) -> SystemDesc {
        self.desc.clone()
    }
    fn run(&mut self, ctx: &StepContext<'_>) {
        (self.body)(ctx)
    }
}

fn system<F: FnMut(&StepContext<'_>) + Send>(
    name: &'static str,
    access: Access,
    rate: Rate,
    body: F,
) -> TestSystem<F> {
    TestSystem {
        desc: SystemDesc {
            name,
            access,
            rate,
            tiers: Tiers::ALL,
        },
        body,
    }
}

/// Records which steps a system ran on.
fn recorder(
    name: &'static str,
    rate: Rate,
    log: &Arc<Mutex<Vec<(&'static str, u64)>>>,
) -> impl System + use<> {
    let log = Arc::clone(log);
    system(name, Access::new(), rate, move |ctx| {
        log.lock().unwrap().push((name, ctx.date().step))
    })
}

#[test]
fn rates_fire_on_the_right_steps() {
    let calendar = Calendar::new(10.0, 4, 24.0).unwrap(); // 4 steps/day, 24-day year, 2-day months
    let log = Arc::new(Mutex::new(Vec::new()));
    let mut s = Scheduler::new(calendar);
    s.add(recorder("t.every3", Rate::Every(3), &log)).unwrap();
    s.add(recorder("t.daily", Rate::Daily, &log)).unwrap();
    s.add(recorder("t.monthly", Rate::Monthly, &log)).unwrap();
    s.add(recorder("t.yearly", Rate::Yearly, &log)).unwrap();
    s.add(recorder("t.demand", Rate::OnDemand, &log)).unwrap();
    let fields = registry();
    s.run_steps(&fields, 4 * 24 * 2); // two years

    let steps = |name: &str| -> Vec<u64> {
        log.lock()
            .unwrap()
            .iter()
            .filter(|(n, _)| *n == name)
            .map(|&(_, s)| s)
            .collect()
    };
    assert_eq!(steps("t.every3"), (0..192).step_by(3).collect::<Vec<_>>());
    assert_eq!(steps("t.daily"), (0..192).step_by(4).collect::<Vec<_>>());
    assert_eq!(steps("t.monthly"), (0..192).step_by(8).collect::<Vec<_>>());
    assert_eq!(steps("t.yearly"), [0, 96]);
    assert!(steps("t.demand").is_empty());
}

#[test]
fn on_demand_runs_once_per_request() {
    let log = Arc::new(Mutex::new(Vec::new()));
    let mut s = Scheduler::new(Calendar::earth());
    let id = s.add(recorder("t.demand", Rate::OnDemand, &log)).unwrap();
    let fields = registry();
    s.run_steps(&fields, 3);
    s.request(id);
    s.run_steps(&fields, 3);
    assert_eq!(*log.lock().unwrap(), [("t.demand", 3)]);
}

#[test]
fn tiers_filter_systems() {
    let log = Arc::new(Mutex::new(Vec::new()));
    let mut s = Scheduler::new(Calendar::earth());
    let dynamics = recorder("t.dynamics", Rate::Every(1), &log);
    struct Tiered<S>(S, Tiers);
    impl<S: System> System for Tiered<S> {
        fn describe(&self) -> SystemDesc {
            SystemDesc {
                tiers: self.1,
                ..self.0.describe()
            }
        }
        fn run(&mut self, ctx: &StepContext<'_>) {
            self.0.run(ctx)
        }
    }
    s.add(Tiered(dynamics, Tiers::DYNAMICS)).unwrap();
    s.add(Tiered(
        recorder("t.climate", Rate::Every(1), &log),
        Tiers::CLIMATE_MODE,
    ))
    .unwrap();
    let fields = registry();
    s.run_steps(&fields, 2);
    s.set_tier(Tier::ClimateMode);
    s.run_steps(&fields, 2);
    assert_eq!(
        *log.lock().unwrap(),
        [
            ("t.dynamics", 0),
            ("t.dynamics", 1),
            ("t.climate", 2),
            ("t.climate", 3)
        ]
    );
}

#[test]
fn stages_respect_conflicts_and_registration_order() {
    let mut r = registry();
    let (x, y, z) = (
        field(&mut r, "t.x", 0.0),
        field(&mut r, "t.y", 0.0),
        field(&mut r, "t.z", 0.0),
    );
    let mut s = Scheduler::new(Calendar::earth());
    let noop = |_: &StepContext<'_>| {};
    let a = s
        .add(system("t.a", Access::new().write(x), Rate::Every(1), noop))
        .unwrap();
    let b = s
        .add(system("t.b", Access::new().read(x), Rate::Every(1), noop))
        .unwrap();
    let c = s
        .add(system(
            "t.c",
            Access::new().read(y).write(z),
            Rate::Every(1),
            noop,
        ))
        .unwrap();
    let d = s
        .add(system("t.d", Access::new().write(x), Rate::Every(1), noop))
        .unwrap();
    let e = s
        .add(system("t.e", Access::new().read(y), Rate::Every(1), noop))
        .unwrap();
    // a: stage 0. b reads x after a writes it: 1. c is independent: 0.
    // d writes x, conflicting with a and b: 2. e only reads y (c reads it too): 0.
    assert_eq!(s.plan(), vec![vec![a, c, e], vec![b], vec![d]]);
}

/// Run a fixed workload; return final field values. Mixes order-sensitive
/// writes to a shared field with parallel independent writers.
fn workload() -> Vec<Vec<f64>> {
    let mut r = registry();
    let x = field(&mut r, "t.x", 1.0);
    let outs: Vec<FieldId<f64>> = ["t.o0", "t.o1", "t.o2", "t.o3"]
        .iter()
        .map(|&n| field(&mut r, n, 0.0))
        .collect();
    let mut s = Scheduler::new(Calendar::earth());
    s.add(system(
        "t.double",
        Access::new().write(x),
        Rate::Every(1),
        move |ctx| {
            for v in ctx.write(x).iter_mut() {
                *v = *v * 1.000_1 + 0.3;
            }
        },
    ))
    .unwrap();
    for (i, &o) in outs.iter().enumerate() {
        let name: &'static str = Box::leak(format!("t.out{i}").into_boxed_str());
        s.add(system(
            name,
            Access::new().read(x).write(o),
            Rate::Every(1 + i as u32),
            move |ctx| {
                let src = ctx.read(x);
                for (k, v) in ctx.write(o).iter_mut().enumerate() {
                    *v = fms_math::sin(*v + src[k] * (i + 1) as f64);
                }
            },
        ))
        .unwrap();
    }
    s.add(system(
        "t.shift",
        Access::new().write(x),
        Rate::Every(2),
        move |ctx| {
            for (k, v) in ctx.write(x).iter_mut().enumerate() {
                *v = fms_math::cos(*v - k as f64) * 10.0;
            }
        },
    ))
    .unwrap();
    s.run_steps(&r, 500);
    std::iter::once(x)
        .chain(outs)
        .map(|id| r.read(id).to_vec())
        .collect()
}

#[test]
fn results_are_identical_across_thread_counts() {
    let run = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(workload)
    };
    let one = run(1);
    assert_eq!(one, run(8));
    assert_eq!(one, run(3));
    // Bitwise equality, not just approximate.
    let bits = |v: &Vec<Vec<f64>>| v.iter().flatten().map(|f| f.to_bits()).collect::<Vec<_>>();
    assert_eq!(bits(&one), bits(&run(8)));
}

#[test]
#[should_panic(expected = "system `t.sneaky` tried to write field `t.x` without declaring it")]
fn undeclared_access_panics() {
    let mut r = registry();
    let x = field(&mut r, "t.x", 0.0);
    let mut s = Scheduler::new(Calendar::earth());
    s.add(system(
        "t.sneaky",
        Access::new().read(x),
        Rate::Every(1),
        move |ctx| {
            ctx.write(x)[0] = 1.0;
        },
    ))
    .unwrap();
    s.step(&r);
}

#[test]
fn interval_seconds_tracks_time_since_last_run() {
    let seen = Arc::new(Mutex::new(Vec::new()));
    let mut s = Scheduler::new(Calendar::new(100.0, 10, 50.0).unwrap()); // dt = 10 s
    let log = Arc::clone(&seen);
    let id = s
        .add(system(
            "t.slow",
            Access::new(),
            Rate::Every(5),
            move |ctx| log.lock().unwrap().push(ctx.interval_seconds()),
        ))
        .unwrap();
    let fields = registry();
    s.run_steps(&fields, 12); // runs at 0, 5, 10
    s.request(id);
    s.step(&fields); // step 12: requested, 2 steps after 10
    assert_eq!(*seen.lock().unwrap(), [50.0, 50.0, 50.0, 20.0]);
}

#[test]
fn add_validates_systems() {
    let mut s = Scheduler::new(Calendar::earth());
    let noop = |_: &StepContext<'_>| {};
    let a = s
        .add(system("t.a", Access::new(), Rate::Daily, noop))
        .unwrap();
    assert_eq!(
        s.add(system("t.a", Access::new(), Rate::Daily, noop))
            .unwrap_err(),
        SchedError::DuplicateName("t.a")
    );
    assert_eq!(
        s.add(system("noname", Access::new(), Rate::Daily, noop))
            .unwrap_err(),
        SchedError::InvalidName("noname")
    );
    assert_eq!(
        s.add(system("t.zero", Access::new(), Rate::Every(0), noop))
            .unwrap_err(),
        SchedError::ZeroRate("t.zero")
    );
    assert_eq!(s.describe(a).name, "t.a");
    assert_eq!(a.index(), 0);
}

#[test]
fn timings_are_reported_and_accumulated() {
    let mut s = Scheduler::new(Calendar::earth());
    let slow = s
        .add(system("t.slow", Access::new(), Rate::Every(1), |_| {
            std::thread::sleep(std::time::Duration::from_millis(2))
        }))
        .unwrap();
    let rare = s
        .add(system("t.rare", Access::new(), Rate::Every(3), |_| {}))
        .unwrap();
    let fields = registry();
    let report = s.step(&fields);
    assert_eq!(report.timings.len(), 2);
    let slow_time = report.timings.iter().find(|(id, _)| *id == slow).unwrap().1;
    assert!(slow_time >= std::time::Duration::from_millis(2));
    assert!(report.elapsed >= slow_time);

    s.run_steps(&fields, 5); // steps 1..=5: slow runs 5 more, rare runs at 3
    let profile = s.profile();
    assert_eq!((profile[0].name, profile[0].runs), ("t.slow", 6));
    assert_eq!((profile[1].id, profile[1].runs), (rare, 2));
    assert!(profile[0].total >= std::time::Duration::from_millis(12));
    s.reset_profile();
    assert!(s.profile().iter().all(|p| p.runs == 0 && p.total.is_zero()));
}
