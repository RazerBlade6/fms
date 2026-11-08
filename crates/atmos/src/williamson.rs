//! Williamson et al. (1992) shallow-water test cases 2, 5 and 6: the standard
//! checks for a dynamical core on the sphere.
//!
//! Velocities are initialized from a stream function at vertices
//! (`uₑ = −(ψ(v₁) − ψ(v₀)) / lₑ`), so the initial flow is exactly
//! non-divergent on the discrete mesh, as in MPAS.

use std::f64::consts::PI;

use fms_geo::coords::to_lat_lon;
use fms_geo::{MeshGeometry, SphereMesh};
use glam::DVec3;

use crate::shallow_water::State;

/// Earth radius used by the test suite (m).
pub const RADIUS: f64 = 6.37122e6;
/// Earth rotation rate used by the test suite (s⁻¹).
pub const OMEGA: f64 = 7.292e-5;
/// Gravity used by the test suite (m s⁻²).
pub const GRAVITY: f64 = 9.80616;
/// Seconds per day.
pub const DAY: f64 = 86_400.0;

/// One of the test cases.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Case {
    /// Test case 2: steady, geostrophically balanced zonal flow, rotated by
    /// `alpha` (radians) relative to the poles. The exact solution is the
    /// initial state.
    SteadyZonal {
        /// Angle between the flow's axis and the rotation axis (0 to π/2).
        alpha: f64,
    },
    /// Test case 5: zonal flow impinging on an isolated conical mountain.
    /// No analytic solution; checked by conservation and convergence.
    Mountain,
    /// Test case 6: Rossby–Haurwitz wave of wavenumber 4, which moves
    /// eastward almost without changing shape.
    RossbyHaurwitz,
}

/// Initial conditions and forcing for a case on one mesh.
pub struct Setup {
    /// Initial state.
    pub state: State,
    /// Coriolis parameter at vertices (s⁻¹).
    pub coriolis: Vec<f64>,
    /// Bottom topography at cells (m).
    pub topography: Vec<f64>,
}

// Case 2 and 5 flow speed and depth.
const TC2_U0: f64 = 2.0 * PI * RADIUS / (12.0 * DAY);
const TC2_GH0: f64 = 2.94e4;
const TC5_U0: f64 = 20.0;
const TC5_H0: f64 = 5960.0;
const TC5_MOUNTAIN_HEIGHT: f64 = 2000.0;
const TC5_MOUNTAIN_RADIUS: f64 = PI / 9.0;
const TC5_MOUNTAIN_LAT: f64 = PI / 6.0;
const TC5_MOUNTAIN_LON: f64 = -PI / 2.0; // 3π/2 in Williamson's [0, 2π)

/// Rossby–Haurwitz wave parameters.
const RH_OMEGA: f64 = 7.848e-6;
const RH_K: f64 = 7.848e-6;
const RH_R: i32 = 4;
const RH_H0: f64 = 8000.0;

impl Case {
    /// The standard run length (days): 5, 15 and 14 for cases 2, 5 and 6.
    pub fn days(&self) -> f64 {
        match self {
            Self::SteadyZonal { .. } => 5.0,
            Self::Mountain => 15.0,
            Self::RossbyHaurwitz => 14.0,
        }
    }

    /// The fastest wind (m s⁻¹) over the standard run, with a little margin,
    /// for choosing a semi-implicit time step: case 2 keeps its initial 38.6;
    /// case 5 accelerates from 20 to about 40; case 6 peaks near 100.
    pub fn max_wind(&self) -> f64 {
        match self {
            Self::SteadyZonal { .. } => TC2_U0,
            Self::Mountain => 42.0,
            Self::RossbyHaurwitz => 105.0,
        }
    }

    /// Short name for files and messages.
    pub fn name(&self) -> &'static str {
        match self {
            Self::SteadyZonal { .. } => "williamson2",
            Self::Mountain => "williamson5",
            Self::RossbyHaurwitz => "williamson6",
        }
    }

    /// Axis of the solid-body flow (and of the Coriolis parameter in case 2).
    fn axis(&self) -> DVec3 {
        match *self {
            Self::SteadyZonal { alpha } => {
                let (s, c) = fms_math::sin_cos(alpha);
                DVec3::new(-s, 0.0, c)
            }
            _ => DVec3::Z,
        }
    }

    /// Coriolis parameter at a point. Case 2 tilts the rotation axis with the
    /// flow, so the balanced state is the same at any `alpha`.
    pub fn coriolis(&self, p: DVec3) -> f64 {
        2.0 * OMEGA * p.dot(self.axis())
    }

    /// Bottom topography at a point (m).
    pub fn topography(&self, p: DVec3) -> f64 {
        match self {
            Self::Mountain => {
                let (lat, lon) = to_lat_lon(p);
                let dlon = (lon - TC5_MOUNTAIN_LON + PI).rem_euclid(2.0 * PI) - PI;
                let dlat = lat - TC5_MOUNTAIN_LAT;
                let r = (dlon * dlon + dlat * dlat).sqrt().min(TC5_MOUNTAIN_RADIUS);
                TC5_MOUNTAIN_HEIGHT * (1.0 - r / TC5_MOUNTAIN_RADIUS)
            }
            _ => 0.0,
        }
    }

    /// Stream function at a point (m² s⁻¹).
    pub fn stream_function(&self, p: DVec3) -> f64 {
        match self {
            Self::SteadyZonal { .. } => -RADIUS * TC2_U0 * p.dot(self.axis()),
            Self::Mountain => -RADIUS * TC5_U0 * p.z,
            Self::RossbyHaurwitz => {
                let (lat, lon) = to_lat_lon(p);
                let (sin_lat, cos_lat) = fms_math::sin_cos(lat);
                let a2 = RADIUS * RADIUS;
                -a2 * RH_OMEGA * sin_lat
                    + a2 * RH_K
                        * cos_lat.powi(RH_R)
                        * sin_lat
                        * fms_math::cos(f64::from(RH_R) * lon)
            }
        }
    }

    /// Initial free-surface height `h + b` at a point (m). For case 2 this is
    /// also the exact solution at all times.
    pub fn surface_height(&self, p: DVec3) -> f64 {
        match self {
            Self::SteadyZonal { .. } => {
                let s = p.dot(self.axis());
                (TC2_GH0 - (RADIUS * OMEGA * TC2_U0 + 0.5 * TC2_U0 * TC2_U0) * s * s) / GRAVITY
            }
            Self::Mountain => {
                TC5_H0 - (RADIUS * OMEGA * TC5_U0 + 0.5 * TC5_U0 * TC5_U0) * p.z * p.z / GRAVITY
            }
            Self::RossbyHaurwitz => {
                let (lat, lon) = to_lat_lon(p);
                let c = fms_math::cos(lat);
                let (w, k, r) = (RH_OMEGA, RH_K, f64::from(RH_R));
                let cr = c.powi(RH_R);
                let a = 0.5 * w * (2.0 * OMEGA + w) * c * c
                    + 0.25
                        * k
                        * k
                        * ((r + 1.0) * cr * cr * c * c + (2.0 * r * r - r - 2.0) * cr * cr
                            - 2.0 * r * r * c.powi(2 * RH_R - 2));
                let b = 2.0 * (OMEGA + w) * k / ((r + 1.0) * (r + 2.0))
                    * cr
                    * ((r * r + 2.0 * r + 2.0) - (r + 1.0) * (r + 1.0) * c * c);
                let cc = 0.25 * k * k * cr * cr * ((r + 1.0) * c * c - (r + 2.0));
                let gh = GRAVITY * RH_H0
                    + RADIUS
                        * RADIUS
                        * (a + b * fms_math::cos(r * lon) + cc * fms_math::cos(2.0 * r * lon));
                gh / GRAVITY
            }
        }
    }

    /// Angular phase speed of the Rossby–Haurwitz wave for non-divergent flow
    /// (rad s⁻¹, eastward): `(R(3 + R)ω − 2Ω) / ((1 + R)(2 + R))`. The
    /// divergent shallow-water wave moves slightly slower.
    pub fn rossby_haurwitz_phase_speed() -> f64 {
        let r = f64::from(RH_R);
        (r * (3.0 + r) * RH_OMEGA - 2.0 * OMEGA) / ((1.0 + r) * (2.0 + r))
    }

    /// Initial state, Coriolis parameter and topography on `mesh`.
    pub fn setup(&self, mesh: &SphereMesh, geometry: &MeshGeometry) -> Setup {
        let topography: Vec<f64> = mesh
            .cell_centers()
            .iter()
            .map(|&p| self.topography(p))
            .collect();
        let h: Vec<f64> = mesh
            .cell_centers()
            .iter()
            .zip(&topography)
            .map(|(&p, b)| self.surface_height(p) - b)
            .collect();
        let psi: Vec<f64> = mesh
            .vertices()
            .iter()
            .map(|&p| self.stream_function(p))
            .collect();
        let u = mesh
            .edges()
            .iter()
            .zip(&geometry.edge_lengths)
            .map(|(e, &l)| {
                -(psi[e.vertices[1] as usize] - psi[e.vertices[0] as usize]) / (l * RADIUS)
            })
            .collect();
        Setup {
            state: State { h, u },
            coriolis: mesh.vertices().iter().map(|&p| self.coriolis(p)).collect(),
            topography,
        }
    }
}

/// Williamson's normalized error norms of `got` against `exact` over cells
/// with areas `areas`: `(l1, l2, l∞)`.
pub fn error_norms(got: &[f64], exact: &[f64], areas: &[f64]) -> (f64, f64, f64) {
    let (mut e1, mut x1, mut e2, mut x2, mut emax, mut xmax) = (0.0, 0.0, 0.0, 0.0, 0.0f64, 0.0f64);
    for ((&g, &x), &a) in got.iter().zip(exact).zip(areas) {
        let d = (g - x).abs();
        e1 += a * d;
        x1 += a * x.abs();
        e2 += a * d * d;
        x2 += a * x * x;
        emax = emax.max(d);
        xmax = xmax.max(x.abs());
    }
    (e1 / x1, (e2 / x2).sqrt(), emax / xmax)
}
