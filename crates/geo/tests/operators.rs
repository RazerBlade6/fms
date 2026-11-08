//! TRiSK operators: convergence against analytic spherical-harmonic fields, and
//! the discrete identities that make a model built on them conserve mass and
//! energy.

use std::sync::OnceLock;

use fms_geo::{LloydOptions, MeshGeometry, Operators, SphereMesh};
use glam::DVec3;

/// Centroidal meshes shared by the tests (each built once).
fn mesh(k: u32) -> &'static (SphereMesh, MeshGeometry) {
    static M: [OnceLock<(SphereMesh, MeshGeometry)>; 3] =
        [OnceLock::new(), OnceLock::new(), OnceLock::new()];
    let slot = match k {
        16 => &M[0],
        32 => &M[1],
        64 => &M[2],
        _ => panic!("unsupported test frequency"),
    };
    slot.get_or_init(|| {
        let (m, _) = SphereMesh::centroidal(k, LloydOptions::default());
        let g = MeshGeometry::new(&m);
        (m, g)
    })
}

const FREQUENCIES: [u32; 3] = [16, 32, 64];

fn log2(x: f64) -> f64 {
    fms_math::ln(x) / std::f64::consts::LN_2
}

/// Deterministic pseudo-random values in [-1, 1].
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
// Analytic fields. Homogeneous harmonic polynomials restricted to the unit
// sphere are spherical harmonics: degree l has ∇²Y = −l(l+1) Y.
// ---------------------------------------------------------------------------

/// Project a 3-D gradient onto the tangent plane at `p`.
fn tangent(p: DVec3, g: DVec3) -> DVec3 {
    g - p * p.dot(g)
}

/// Scalar test field: `xyz` (l = 3) plus `½(x² − y²)` (l = 2).
fn h(p: DVec3) -> f64 {
    p.x * p.y * p.z + 0.5 * (p.x * p.x - p.y * p.y)
}
fn grad_h(p: DVec3) -> DVec3 {
    tangent(p, DVec3::new(p.y * p.z + p.x, p.x * p.z - p.y, p.x * p.y))
}
fn laplacian_h(p: DVec3) -> f64 {
    -12.0 * p.x * p.y * p.z - 6.0 * 0.5 * (p.x * p.x - p.y * p.y)
}

/// Velocity potential `χ = x² − y² + 0.3 z` and stream function `ψ = xyz`.
fn chi_laplacian(p: DVec3) -> f64 {
    -6.0 * (p.x * p.x - p.y * p.y) - 2.0 * 0.3 * p.z
}
fn psi_laplacian(p: DVec3) -> f64 {
    -12.0 * p.x * p.y * p.z
}
/// `u = ∇χ + k × ∇ψ`: divergence ∇²χ, vorticity ∇²ψ.
fn velocity(p: DVec3) -> DVec3 {
    let grad_chi = tangent(p, DVec3::new(2.0 * p.x, -2.0 * p.y, 0.3));
    let grad_psi = tangent(p, DVec3::new(p.y * p.z, p.x * p.z, p.x * p.y));
    grad_chi + p.cross(grad_psi)
}

/// (max, rms) of `got − exact`, relative to the largest |exact|.
fn errors(got: &[f64], exact: &[f64]) -> (f64, f64) {
    let scale = exact.iter().fold(0.0f64, |a, &b| a.max(b.abs()));
    let max = got
        .iter()
        .zip(exact)
        .fold(0.0f64, |a, (g, e)| a.max((g - e).abs()));
    let rms = (got
        .iter()
        .zip(exact)
        .map(|(g, e)| (g - e).powi(2))
        .sum::<f64>()
        / got.len() as f64)
        .sqrt();
    (max / scale, rms / scale)
}

/// Errors at each test frequency, and the observed orders between them.
struct Convergence {
    name: &'static str,
    errors: Vec<(f64, f64)>,
}

impl Convergence {
    fn measure(name: &'static str, f: impl Fn(&SphereMesh, &MeshGeometry) -> (f64, f64)) -> Self {
        let errors = FREQUENCIES
            .iter()
            .map(|&k| {
                let (m, g) = mesh(k);
                f(m, g)
            })
            .collect();
        let c = Self { name, errors };
        eprintln!(
            "{:<16} max {:?} (orders {:?})  rms {:?} (orders {:?})",
            c.name,
            c.errors
                .iter()
                .map(|e| format!("{:.2e}", e.0))
                .collect::<Vec<_>>(),
            c.orders(|e| e.0),
            c.errors
                .iter()
                .map(|e| format!("{:.2e}", e.1))
                .collect::<Vec<_>>(),
            c.orders(|e| e.1),
        );
        c
    }

    /// log₂(error ratio) between successive doublings of resolution.
    fn orders(&self, pick: impl Fn(&(f64, f64)) -> f64) -> Vec<String> {
        self.errors
            .windows(2)
            .map(|w| format!("{:.2}", log2(pick(&w[0]) / pick(&w[1]))))
            .collect()
    }

    /// Assert the order at the finest doubling is at least `max_order` and
    /// `rms_order` (use 0 to only require the error to stay bounded).
    fn assert_orders(&self, max_order: f64, rms_order: f64) {
        let [.., a, b] = self.errors[..] else {
            unreachable!()
        };
        let max = log2(a.0 / b.0);
        let rms = log2(a.1 / b.1);
        assert!(
            max >= max_order && rms >= rms_order,
            "{}: orders max {max:.2}, rms {rms:.2}; wanted ≥ {max_order}, ≥ {rms_order}",
            self.name
        );
    }
}

fn normal_components(g: &MeshGeometry, f: impl Fn(DVec3) -> DVec3) -> Vec<f64> {
    g.edge_points
        .iter()
        .zip(&g.edge_normals)
        .map(|(&x, &n)| f(x).dot(n))
        .collect()
}

// ---------------------------------------------------------------------------
// Convergence
// ---------------------------------------------------------------------------

#[test]
fn gradient_converges_at_second_order() {
    let c = Convergence::measure("gradient", |m, g| {
        let ops = Operators::new(m, g, 1.0);
        let hc: Vec<f64> = m.cell_centers().iter().map(|&p| h(p)).collect();
        let mut got = vec![0.0; m.n_edges()];
        ops.gradient(&hc, &mut got);
        errors(&got, &normal_components(g, grad_h))
    });
    c.assert_orders(1.8, 1.8);
}

#[test]
fn divergence_converges() {
    let c = Convergence::measure("divergence", |m, g| {
        let ops = Operators::new(m, g, 1.0);
        let u = normal_components(g, velocity);
        let mut got = vec![0.0; m.n_cells()];
        ops.divergence(&u, &mut got);
        let exact: Vec<f64> = m.cell_centers().iter().map(|&p| chi_laplacian(p)).collect();
        errors(&got, &exact)
    });
    // RMS converges at first order. The maximum stays at ~1.1%, in the
    // hexagons next to the pentagons (OBSERVATIONS.md O-1; see the
    // `divergence_error` example). Measured max 1.24e-2, 1.15e-2, 1.13e-2.
    c.assert_orders(-0.2, 0.9);
    assert!(c.errors.iter().all(|e| e.0 < 0.015), "{:?}", c.errors);
}

#[test]
fn laplacian_converges() {
    // div(grad h): the operator the Helmholtz solver and hyperdiffusion use.
    let c = Convergence::measure("laplacian", |m, g| {
        let ops = Operators::new(m, g, 1.0);
        let hc: Vec<f64> = m.cell_centers().iter().map(|&p| h(p)).collect();
        let mut grad = vec![0.0; m.n_edges()];
        let mut got = vec![0.0; m.n_cells()];
        ops.gradient(&hc, &mut grad);
        ops.divergence(&grad, &mut got);
        let exact: Vec<f64> = m.cell_centers().iter().map(|&p| laplacian_h(p)).collect();
        errors(&got, &exact)
    });
    // As for divergence: RMS converges, the maximum (at the pentagons) stays
    // at ~1.2%.
    c.assert_orders(-0.2, 0.9);
    assert!(c.errors.iter().all(|e| e.0 < 0.015), "{:?}", c.errors);
}

#[test]
fn curl_converges_at_first_order() {
    let c = Convergence::measure("curl", |m, g| {
        let ops = Operators::new(m, g, 1.0);
        let u = normal_components(g, velocity);
        let mut got = vec![0.0; m.n_vertices()];
        ops.curl(&u, &mut got);
        let exact: Vec<f64> = m.vertices().iter().map(|&p| psi_laplacian(p)).collect();
        errors(&got, &exact)
    });
    // The circulation is the mean vorticity over the dual triangle, which is
    // second-order accurate at the triangle's centroid but only first order at
    // the Voronoi vertex (its circumcentre), where it is stored.
    c.assert_orders(0.9, 0.9);
}

#[test]
fn tangential_velocity_converges_in_rms() {
    let c = Convergence::measure("tangential", |m, g| {
        let ops = Operators::new(m, g, 1.0);
        let u = normal_components(g, velocity);
        let mut got = vec![0.0; m.n_edges()];
        ops.tangential(&u, &mut got);
        let exact: Vec<f64> = g
            .edge_points
            .iter()
            .zip(&g.edge_tangents)
            .map(|(&x, &t)| velocity(x).dot(t))
            .collect();
        errors(&got, &exact)
    });
    // Max error stays fixed near the pentagons (OBSERVATIONS.md O-1).
    c.assert_orders(-0.3, 0.8);
    assert!(c.errors.iter().all(|e| e.0 < 0.02), "{:?}", c.errors);
}

#[test]
fn kinetic_energy_converges_in_rms() {
    let c = Convergence::measure("kinetic energy", |m, g| {
        let ops = Operators::new(m, g, 1.0);
        let u = normal_components(g, velocity);
        let mut got = vec![0.0; m.n_cells()];
        ops.kinetic_energy(&u, &mut got);
        let exact: Vec<f64> = m
            .cell_centers()
            .iter()
            .map(|&p| 0.5 * velocity(p).length_squared())
            .collect();
        errors(&got, &exact)
    });
    c.assert_orders(-0.3, 0.8);
    assert!(c.errors.iter().all(|e| e.0 < 0.05), "{:?}", c.errors);
}

#[test]
fn averages_converge() {
    let c = Convergence::measure("cell → vertex", |m, g| {
        let ops = Operators::new(m, g, 1.0);
        let hc: Vec<f64> = m.cell_centers().iter().map(|&p| h(p)).collect();
        let mut got = vec![0.0; m.n_vertices()];
        ops.cell_to_vertex(&hc, &mut got);
        let exact: Vec<f64> = m.vertices().iter().map(|&p| h(p)).collect();
        errors(&got, &exact)
    });
    // First order, like curl: a kite-weighted mean belongs at the dual
    // triangle's centroid, not its circumcentre.
    c.assert_orders(0.9, 0.9);
    // cell → edge is centred on the dual edge: second order. vertex → edge is
    // evaluated where the dual edge crosses the Voronoi edge, which is off the
    // Voronoi edge's midpoint where cells are distorted: first order in the
    // maximum, second in RMS.
    for (name, to_edge, max_order) in [("cell → edge", true, 1.8), ("vertex → edge", false, 0.9)]
    {
        let c = Convergence::measure(name, |m, g| {
            let ops = Operators::new(m, g, 1.0);
            let mut got = vec![0.0; m.n_edges()];
            if to_edge {
                let hc: Vec<f64> = m.cell_centers().iter().map(|&p| h(p)).collect();
                ops.cell_to_edge(&hc, &mut got);
            } else {
                let hv: Vec<f64> = m.vertices().iter().map(|&p| h(p)).collect();
                ops.vertex_to_edge(&hv, &mut got);
            }
            let exact: Vec<f64> = g.edge_points.iter().map(|&p| h(p)).collect();
            errors(&got, &exact)
        });
        c.assert_orders(max_order, 1.8);
    }
}

#[test]
fn operators_scale_with_the_radius() {
    let (m, g) = mesh(16);
    let (unit, earth) = (Operators::new(m, g, 1.0), Operators::new(m, g, 6.371e6));
    let hc = noise(m.n_cells(), 7);
    let u = noise(m.n_edges(), 8);
    let mut a = vec![0.0; m.n_edges()];
    let mut b = vec![0.0; m.n_edges()];
    unit.gradient(&hc, &mut a);
    earth.gradient(&hc, &mut b);
    assert!(
        a.iter()
            .zip(&b)
            .all(|(x, y)| (x / 6.371e6 - y).abs() < 1e-18)
    );
    let mut a = vec![0.0; m.n_vertices()];
    let mut b = vec![0.0; m.n_vertices()];
    unit.curl(&u, &mut a);
    earth.curl(&u, &mut b);
    assert!(
        a.iter()
            .zip(&b)
            .all(|(x, y)| (x / 6.371e6 - y).abs() < 1e-18)
    );
    // Dimensionless: tangential velocity and kinetic energy don't scale.
    let mut a = vec![0.0; m.n_cells()];
    let mut b = vec![0.0; m.n_cells()];
    unit.kinetic_energy(&u, &mut a);
    earth.kinetic_energy(&u, &mut b);
    assert_eq!(a, b);
}

// ---------------------------------------------------------------------------
// Discrete identities (exact to round-off, on random fields)
// ---------------------------------------------------------------------------

const R: f64 = 6.371e6;

fn sum_weighted(a: &[f64], w: &[f64]) -> f64 {
    a.iter().zip(w).map(|(x, y)| x * y).sum()
}

fn cell_areas(g: &MeshGeometry) -> Vec<f64> {
    g.cell_areas.iter().map(|a| a * R * R).collect()
}

#[test]
fn divergence_has_no_net_source() {
    // Gauss: Σ Aᵢ div(F)ᵢ = 0, so a flux-form continuity equation conserves mass.
    let (m, g) = mesh(32);
    let ops = Operators::new(m, g, R);
    let f = noise(m.n_edges(), 1);
    let mut div = vec![0.0; m.n_cells()];
    ops.divergence(&f, &mut div);
    let areas = cell_areas(g);
    let scale: f64 = div.iter().zip(&areas).map(|(d, a)| (d * a).abs()).sum();
    assert!(sum_weighted(&div, &areas).abs() < 1e-13 * scale);
}

#[test]
fn total_circulation_is_zero() {
    // Stokes on a closed surface: Σ A_v ζ_v = 0.
    let (m, g) = mesh(32);
    let ops = Operators::new(m, g, R);
    let u = noise(m.n_edges(), 2);
    let mut zeta = vec![0.0; m.n_vertices()];
    ops.curl(&u, &mut zeta);
    let areas: Vec<f64> = g.vertex_areas.iter().map(|a| a * R * R).collect();
    let scale: f64 = zeta.iter().zip(&areas).map(|(z, a)| (z * a).abs()).sum();
    assert!(sum_weighted(&zeta, &areas).abs() < 1e-13 * scale);
}

#[test]
fn curl_of_gradient_vanishes() {
    let (m, g) = mesh(32);
    let ops = Operators::new(m, g, R);
    let hc = noise(m.n_cells(), 3);
    let mut grad = vec![0.0; m.n_edges()];
    let mut zeta = vec![0.0; m.n_vertices()];
    ops.gradient(&hc, &mut grad);
    ops.curl(&grad, &mut zeta);
    let scale = grad.iter().fold(0.0f64, |a, b| a.max(b.abs())) / (m.mean_spacing() * R);
    assert!(zeta.iter().all(|z| z.abs() < 1e-12 * scale));
}

#[test]
fn gradient_and_divergence_are_adjoint() {
    // Σₑ 2Aₑ grad(h)ₑ Fₑ = −Σᵢ Aᵢ hᵢ div(F)ᵢ (2Aₑ = lₑdₑ): the work done by
    // pressure gradients equals the potential energy released by convergence.
    let (m, g) = mesh(32);
    let ops = Operators::new(m, g, R);
    let hc = noise(m.n_cells(), 4);
    let f = noise(m.n_edges(), 5);
    let mut grad = vec![0.0; m.n_edges()];
    let mut div = vec![0.0; m.n_cells()];
    ops.gradient(&hc, &mut grad);
    ops.divergence(&f, &mut div);
    let edge_areas = ops.edge_areas();
    let lhs: f64 = (0..m.n_edges())
        .map(|e| 2.0 * edge_areas[e] * grad[e] * f[e])
        .sum();
    let rhs: f64 = -sum_weighted(
        &hc.iter().zip(&div).map(|(a, b)| a * b).collect::<Vec<_>>(),
        &cell_areas(g),
    );
    assert!((lhs - rhs).abs() < 1e-12 * lhs.abs(), "{lhs} vs {rhs}");
}

#[test]
fn kinetic_energy_sums_to_the_edge_energy() {
    // Σᵢ Aᵢ Kᵢ = Σₑ Aₑ uₑ²: cell kinetic energy redistributes edge energy.
    let (m, g) = mesh(32);
    let ops = Operators::new(m, g, R);
    let u = noise(m.n_edges(), 6);
    let mut ke = vec![0.0; m.n_cells()];
    ops.kinetic_energy(&u, &mut ke);
    let cells = sum_weighted(&ke, &cell_areas(g));
    let edges: f64 = ops
        .edge_areas()
        .iter()
        .zip(&u)
        .map(|(a, u)| a * u * u)
        .sum();
    assert!((cells - edges).abs() < 1e-12 * edges);
}

#[test]
fn cell_to_vertex_conserves_the_area_integral() {
    // Kites partition both cells and dual triangles: Σ A_v h_v = Σ Aᵢ hᵢ.
    let (m, g) = mesh(32);
    let ops = Operators::new(m, g, R);
    let hc = noise(m.n_cells(), 9);
    let mut hv = vec![0.0; m.n_vertices()];
    ops.cell_to_vertex(&hc, &mut hv);
    let a = sum_weighted(&hc, &g.cell_areas);
    let b = sum_weighted(&hv, &g.vertex_areas);
    assert!((a - b).abs() < 1e-12, "{a} vs {b}");
}

#[test]
fn linear_shallow_water_conserves_energy() {
    // ∂h/∂t = −H div(u), ∂u/∂t = −g grad(h) + f v(u). Total energy
    // E = Σₑ Aₑ H uₑ² + Σᵢ Aᵢ g hᵢ²/2, whose kinetic part is H Σᵢ Aᵢ Kᵢ (the
    // mean of uₙ² over edge directions is ½|u|²). Its tendency is zero for any
    // state: the pressure terms cancel by adjointness and the Coriolis term
    // (antisymmetric TRiSK weights) does no work.
    let (m, g) = mesh(32);
    let ops = Operators::new(m, g, R);
    let (depth, gravity, f) = (1000.0, 9.81, 1.0e-4);
    let h = noise(m.n_cells(), 10);
    let u = noise(m.n_edges(), 11);

    let mut div = vec![0.0; m.n_cells()];
    let mut grad = vec![0.0; m.n_edges()];
    let mut v = vec![0.0; m.n_edges()];
    ops.divergence(&u, &mut div);
    ops.gradient(&h, &mut grad);
    ops.tangential(&u, &mut v);
    let dh: Vec<f64> = div.iter().map(|d| -depth * d).collect();
    let du: Vec<f64> = grad
        .iter()
        .zip(&v)
        .map(|(gr, v)| -gravity * gr + f * v)
        .collect();

    let edge_areas = ops.edge_areas();
    let kinetic: f64 = (0..m.n_edges())
        .map(|e| 2.0 * edge_areas[e] * depth * u[e] * du[e])
        .sum();
    let potential: f64 = (0..m.n_cells())
        .map(|i| g.cell_areas[i] * R * R * gravity * h[i] * dh[i])
        .sum();
    let scale = kinetic.abs() + potential.abs();
    assert!(
        (kinetic + potential).abs() < 1e-12 * scale,
        "dE/dt = {} (terms {kinetic}, {potential})",
        kinetic + potential
    );
}
