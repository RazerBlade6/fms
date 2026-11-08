//! Moisture (M2 step 9), dynamics side: humidity carried through a
//! baroclinic storm by the semi-implicit stepper's flux-corrected transport
//! stays positive and within its range, and water is exactly accounted for.
//! (Consistency with the air's mass: `transport.rs`'s unit tests.)

use std::sync::OnceLock;

use fms_atmos::williamson::{GRAVITY, OMEGA, RADIUS};
use fms_atmos::{
    Gas, HyperdiffusionOptions, LayeredModel, LayeredSemiImplicit, LayeredSemiImplicitOptions,
    LayeredState, Levels, ShearedFlow,
};
use fms_geo::coords::to_lat_lon;
use fms_geo::{LloydOptions, MeshGeometry, Operators, SphereMesh};

const DAY: f64 = 86_400.0;
const MAX_WIND: f64 = 120.0;

fn mesh() -> &'static (SphereMesh, MeshGeometry) {
    static M: OnceLock<(SphereMesh, MeshGeometry)> = OnceLock::new();
    M.get_or_init(|| {
        let (m, _) = SphereMesh::centroidal(16, LloydOptions::default());
        let g = MeshGeometry::new(&m);
        (m, g)
    })
}

/// A diffused model and its step (as in `dissipation.rs`).
fn model() -> (LayeredModel<'static>, f64) {
    let (m, g) = mesh();
    let model = LayeredModel::new(
        Operators::new(m, g, RADIUS),
        GRAVITY,
        Gas::EARTH_AIR,
        Levels::two(0.5),
        m.vertices().iter().map(|p| 2.0 * OMEGA * p.z).collect(),
        vec![0.0; m.n_cells()],
    );
    let dt = LayeredSemiImplicit::stable_dt(&model, MAX_WIND, 0.6);
    let dt = DAY / (DAY / dt).ceil();
    let options = HyperdiffusionOptions::strongest_for_step(model.operators(), 4, dt);
    (model.with_hyperdiffusion(options), dt)
}

/// The `fms baroclinic` storm: a 30 m/s-shear jet seeded with waves 4–8.
fn storm(model: &LayeredModel) -> LayeredState {
    let flow = ShearedFlow {
        speed: [40.0, 10.0],
        stability: 40.0,
        theta_equator: 315.0,
        ps_equator: 1.0e5,
        omega: OMEGA,
    };
    let mut state = flow.state(model);
    for k in 0..2 {
        for (i, &p) in mesh().0.cell_centers().iter().enumerate() {
            let (lat, lon) = to_lat_lon(p);
            let band =
                |c: f64| fms_math::exp(-((lat - c.to_radians()) / 15f64.to_radians()).powi(2));
            let waves: f64 = (4..=8)
                .map(|m| fms_math::cos(f64::from(m) * lon + 0.7 * f64::from(m)))
                .sum();
            state.theta_ps[k][i] += state.ps[i] * 0.5 * (band(45.0) + band(-45.0)) * waves;
        }
    }
    state
}

fn run(state: LayeredState, days: f64) -> LayeredSemiImplicit<'static> {
    let (model, dt) = model();
    let mut si =
        LayeredSemiImplicit::new(model, state, dt, LayeredSemiImplicitOptions::default(), 0.5);
    let steps = (days * DAY / dt).round() as u64;
    while si.steps() < steps {
        si.step();
        si.check().expect("storm stays healthy");
    }
    si
}

#[test]
fn a_sharp_plume_stays_positive_and_water_is_conserved() {
    // A square plume of moist air in the storm's path, in both layers.
    let (model, _) = model();
    let state = storm(&model);
    let centers = mesh().0.cell_centers();
    let q: Vec<f64> = centers
        .iter()
        .map(|&p| {
            let (lat, lon) = to_lat_lon(p);
            let inside = (lat.to_degrees() - 45.0).abs() < 12.0 && lon.to_degrees().abs() < 30.0;
            if inside { 0.015 } else { 0.0 }
        })
        .collect();
    let state = state.with_humidity(vec![q.iter().map(|q| 0.1 * q).collect(), q.clone()]);
    let mut si = run(state, 0.0);
    let start = si.water();
    let days = 12.0;
    let steps = (days * DAY / si.dt()).round() as u64;
    let mut most_negative = 0.0f64;
    let mut highest = 0.0f64;
    while si.steps() < steps {
        si.step();
        for k in 0..2 {
            for q in si.state().humidity(k) {
                most_negative = most_negative.min(q);
                highest = highest.max(q);
            }
        }
    }
    let b = si.water_budget();
    let residual = si.water() - start - b.sources - b.filter - b.fixer;
    eprintln!(
        "water {start:.6e} kg; filter {:+.2e}, fixer {:+.2e}, borrowed {:.2e}, residual {residual:+.1e}; max substeps {}; q in [{most_negative:.1e}, {highest:.5}]",
        b.filter, b.fixer, b.borrowed, b.max_substeps
    );
    // Positive. The maximum may rise above the plume's: where air
    // converges and rises, only the dry top of the lower layer's profile
    // leaves across the interface, so the layer's mean humidity builds up (in
    // the full model it condenses). Measured 0.025.
    assert!(most_negative >= 0.0);
    assert!(highest <= 0.04, "humidity reached {highest}");
    // The budget closes to round-off. The fixer is a trace; the time
    // filter's clamping is small even for this square plume (measured
    // +6e-4 of the water over 12 days; smooth fields need far less).
    assert!(residual.abs() < 1e-12 * start, "residual {residual:e}");
    assert!(b.filter.abs() < 2e-3 * start);
    assert!(b.fixer < 1e-9 * start);
    // The plume has been stirred: water has spread well beyond the box.
    let spread = si
        .state()
        .humidity(1)
        .iter()
        .filter(|&&q| q > 1e-4 && q < 0.0149)
        .count();
    assert!(spread > 100, "only {spread} cells in between");
}
