//! The overturning diagnostic on a real mesh (M2.14.1): the sampled mass
//! flux across each band boundary is exactly the model's flux divergence
//! summed over the cells south of it, and an analytic overturning comes back
//! with its cells, edges and strength.

use fms_atmos::williamson::{GRAVITY, OMEGA, RADIUS};
use fms_atmos::{Circulation, Gas, LayeredModel, LayeredState, Levels, ZonalClimatology};
use fms_grid::coords::to_lat_lon;
use fms_grid::zonal::east_north;
use fms_grid::{LatitudeBands, LloydOptions, MeshGeometry, Operators, SphereMesh};

const BANDS: usize = 36;
const PS: f64 = 1e5;

/// The upper layer's northward wind (m/s); the lower layer has its negative.
fn v(lat: f64) -> f64 {
    10.0 * fms_math::sin(6.0 * lat)
}

/// The analytic streamfunction at σ 0.5 (kg/s): the upper layer's
/// northward mass flux around the latitude circle.
fn analytic(lat: f64) -> f64 {
    2.0 * std::f64::consts::PI * RADIUS * fms_math::cos(lat) * v(lat) * PS * 0.5 / GRAVITY
}

#[test]
fn the_overturning_closes_exactly_and_matches_an_analytic_flow() {
    let (mesh, _) = SphereMesh::centroidal(16, LloydOptions::default());
    let geometry = MeshGeometry::new(&mesh);
    let model = LayeredModel::new(
        Operators::new(&mesh, &geometry, RADIUS),
        GRAVITY,
        Gas::EARTH_AIR,
        Levels::two(0.5),
        mesh.vertices().iter().map(|p| 2.0 * OMEGA * p.z).collect(),
        vec![0.0; mesh.n_cells()],
    );
    // Northward wind V(φ) projected on each edge normal.
    let upper: Vec<f64> = geometry
        .edge_points
        .iter()
        .zip(&geometry.edge_normals)
        .map(|(&p, &n)| v(to_lat_lon(p).0) * east_north(p, n).1)
        .collect();
    let lower: Vec<f64> = upper.iter().map(|u| -u).collect();
    let state = LayeredState::new(
        vec![PS; mesh.n_cells()],
        vec![upper, lower],
        vec![vec![300.0; mesh.n_cells()]; 2],
    );
    let mut climate = ZonalClimatology::new(&model, BANDS);
    climate.sample(&model, &state);
    let o = climate.overturning();
    assert_eq!(o.latitude.len(), BANDS + 1);
    assert_eq!(o.sigma, [0.5, 1.0]);

    // Exact bookkeeping: each layer's flux across boundary b is the mass
    // leaving the cells south of it, from the model's own divergence.
    let ops = model.operators();
    let bands = LatitudeBands::new(&mesh, &geometry, BANDS);
    let mut ps_edge = vec![0.0; mesh.n_edges()];
    ops.cell_to_edge(&state.ps, &mut ps_edge);
    let largest = o.streamfunction[0]
        .iter()
        .fold(0.0f64, |m, f| m.max(f.abs()));
    let mut worst = 0.0f64;
    for k in 0..2 {
        let flux: Vec<f64> = ps_edge
            .iter()
            .zip(&state.u[k])
            .map(|(p, u)| model.dsigma()[k] * p * u)
            .collect();
        let mut divergence = vec![0.0; mesh.n_cells()];
        ops.divergence(&flux, &mut divergence);
        for b in 0..=BANDS {
            let leaving: f64 = (0..mesh.n_cells())
                .filter(|&i| (bands.band_of_cell()[i] as usize) < b)
                .map(|i| divergence[i] * geometry.cell_areas[i] * RADIUS * RADIUS / GRAVITY)
                .sum();
            let sampled = o.streamfunction[k][b]
                - if k > 0 {
                    o.streamfunction[k - 1][b]
                } else {
                    0.0
                };
            worst = worst.max((sampled - leaving).abs() / largest);
        }
    }
    eprintln!("bookkeeping: largest difference {worst:.1e} of the largest flux");
    assert!(worst < 1e-10, "{worst:e}");

    // Analytic: the streamfunction at σ 0.5 follows 2π a cos φ V ps δσ / g,
    // to within the jagged boundaries' error.
    let peak = (0..=9000)
        .map(|j| analytic((j as f64 / 100.0).to_radians()).abs())
        .fold(0.0f64, f64::max);
    let psi = o.at(0.5);
    let error = o.latitude.iter().zip(psi).fold(0.0f64, |m, (&l, &p)| {
        m.max((p - analytic(l.to_radians())).abs())
    }) / peak;
    eprintln!(
        "analytic: peak {:.1} × 10⁹ kg/s, largest error {:.2}% of it",
        peak / 1e9,
        100.0 * error
    );
    // Measured: 4.9%, at ±5° where ψ is steepest: the discrete latitude
    // circle (cells whose centre is south of it) zigzags about the true one
    // by an effective 0.5° there; 0.1° at the equator (ψ = −22 × 10⁹ kg/s).
    assert!(error < 0.08, "{error}");
    // The layers cancel exactly, so the surface closes to round-off.
    eprintln!("residual {:e} kg/s", o.residual());
    assert!(o.residual() <= 1e-12 * peak, "{}", o.residual());

    // The cells: Hadley, Ferrel and polar in each hemisphere (measured:
    // Hadley edges at 30.4° N and 30.3° S, 0.7% weaker than the analytic
    // peak, read at the boundary at 15°).
    let c = Circulation::of(&climate.means(), &o, 0, 1);
    eprintln!("cells {:?}", c.cells);
    assert_eq!(c.cells_in(true), 3);
    assert_eq!(c.cells_in(false), 3);
    for north in [true, false] {
        let edge = c.hadley_edge(north).unwrap();
        let strength = c.hadley(north).unwrap().strength;
        assert!((edge - 30.0).abs() < 2.0, "Hadley edge {edge}°");
        assert!(
            (strength.abs() / peak - 1.0).abs() < 0.05,
            "Hadley strength {strength:e} against {peak:e}"
        );
        assert_eq!(strength > 0.0, north);
    }
}
