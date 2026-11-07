//! Insolation over a climate-sized set of points.
//!
//! `cargo bench -p fms-planet`

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use fms_geo::coords::from_lat_lon;
use fms_planet::Planet;

fn insolation(c: &mut Criterion) {
    let earth = Planet::earth();
    let points: Vec<_> = (0..20_252)
        .map(|i| {
            let f = f64::from(i);
            from_lat_lon(
                (f * 0.618_034).fract() * 3.0 - 1.5,
                (f * 0.414_214).fract() * 6.2 - 3.1,
            )
        })
        .collect();
    let mut out = vec![0.0; points.len()];
    c.bench_function("planet/insolation_20k", |b| {
        b.iter(|| earth.insolation_into(black_box(&points), 1.0e7, &mut out))
    });
    c.bench_function("planet/daily_mean_20k", |b| {
        b.iter(|| {
            points
                .iter()
                .map(|p| earth.daily_mean_insolation(fms_math::asin(p.z), 1.0e7))
                .sum::<f64>()
        })
    });
    c.bench_function("planet/orbit_state", |b| {
        b.iter(|| earth.orbit_state(black_box(1.234e7)))
    });
}

criterion_group!(benches, insolation);
criterion_main!(benches);
