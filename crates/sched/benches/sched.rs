//! Scheduler overhead per step, with trivial systems.
//!
//! `cargo bench -p fms-sched`

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use fms_fields::{FieldId, FieldRegistry, FieldSpec, GridSizes, Location, MeshCounts};
use fms_sched::{Access, Calendar, Rate, Scheduler, StepContext, System, SystemDesc, Tiers};

struct Touch {
    name: &'static str,
    access: Access,
    out: Option<FieldId<f64>>,
}

impl System for Touch {
    fn describe(&self) -> SystemDesc {
        SystemDesc {
            name: self.name,
            access: self.access.clone(),
            rate: Rate::Every(1),
            tiers: Tiers::ALL,
        }
    }
    fn run(&mut self, ctx: &StepContext<'_>) {
        if let Some(out) = self.out {
            ctx.write(out)[0] += 1.0;
        }
    }
}

fn scheduler(c: &mut Criterion) {
    let counts = MeshCounts {
        cells: 1000,
        edges: 2994,
        vertices: 1996,
    };
    let mut fields = FieldRegistry::new(GridSizes {
        terrain: counts,
        climate: counts,
    });
    let mut s = Scheduler::new(Calendar::earth());
    // 20 systems: pairs writing a shared field (forcing stages) plus independent ones.
    for i in 0..20 {
        let name: &'static str = Box::leak(format!("bench.system_{i}").into_boxed_str());
        let field_name: &'static str = Box::leak(format!("bench.field_{}", i / 2).into_boxed_str());
        let id = fields
            .id::<f64>(field_name)
            .or_else(|_| {
                fields.register(FieldSpec::new(field_name, "", Location::CLIMATE_CELLS), 0.0)
            })
            .unwrap();
        s.add(Touch {
            name,
            access: Access::new().write(id),
            out: Some(id),
        })
        .unwrap();
    }
    c.bench_function("scheduler/step_20_systems", |b| {
        b.iter(|| black_box(s.step(&fields)))
    });
}

criterion_group!(benches, scheduler);
criterion_main!(benches);
