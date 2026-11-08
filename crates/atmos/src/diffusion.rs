//! Scale-selective hyperdiffusion `∇²ⁿ` for the layered model: the
//! dissipation that stops the enstrophy cascade piling up at the grid scale
//! (docs/OBSERVATIONS.md O-14).
//!
//! # Operators
//!
//! With `L_c = divergence ∘ gradient` (cells) and
//! `L_v = −curl ∘ tangential_gradient` (vertices), the TRiSK vector Laplacian
//! is `L u = gradient(δ) − tangential_gradient(ζ)`. Divergence and curl pass
//! through it exactly (`div L u = L_c δ`, `curl L u = L_v ζ`), so
//!
//! ```text
//! Lⁿ u = gradient(L_cⁿ⁻¹ δ) − tangential_gradient(L_vⁿ⁻¹ ζ)
//! ```
//!
//! and the tendencies are
//!
//! ```text
//! ∂uₖ/∂t = (−1)ⁿ⁺¹ ν Lⁿ uₖ
//! ∂Θₖ/∂t = (−1)ⁿ⁺¹ ν divergence(p̂_s gradient(L_cⁿ⁻¹ θₖ))
//! ```
//!
//! with `p̂_s` the surface pressure at edges. Potential temperature is
//! diffused (not `Θ = p_s θ`), in flux form: total `Θ` is conserved exactly,
//! and a uniform θ is left alone, so an isentropic atmosphere over mountains
//! stays exactly at rest. Surface pressure is not diffused.
//!
//! # Strength
//!
//! One coefficient `ν` serves winds and θ alike. It is set from an e-folding
//! time `τ` for the shortest mode the cells resolve: `ν = 1/(τ λ_cⁿ)`, where
//! `λ_c` is the largest eigenvalue of `−L_c` on the actual mesh (power
//! iteration). So the damping scales with resolution by itself, and a wave of
//! spherical-harmonic degree `l` on a planet of radius `a` e-folds in
//! `τ (λ_c a² / l(l+1))ⁿ`.
//!
//! The vertex Laplacian reaches higher: on a hexagonal C-grid there are about
//! twice as many vertices as cells, and the extra rotational modes zig-zag
//! between neighbouring vertices with `λ_v ≈ 3.7 λ_c`. They are damped
//! `(λ_v/λ_c)ⁿ` times faster than `τ` (14× for `∇⁴`, 190× for `∇⁸`), which is
//! what limits the time step of explicit diffusion:
//! [`Hyperdiffusion::max_rate`].
//!
//! The kinetic energy removed is not returned as heat; see
//! [`crate::LayeredModel::energy_tendency`] for the budget.

use fms_geo::Operators;
use fms_geo::operators::MIN_CHUNK;
use rayon::prelude::*;

use crate::helmholtz::dot;
use crate::layered::LayeredState;

/// Settings for [`Hyperdiffusion`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HyperdiffusionOptions {
    /// Order of the operator: 4 for `∇⁴`, 6 for `∇⁶`, 8 for `∇⁸`. Even, ≥ 2.
    pub order: u32,
    /// E-folding time (s) of the shortest resolvable mode.
    pub e_folding: f64,
}

impl Default for HyperdiffusionOptions {
    /// `∇⁴` with a grid-scale e-folding time of 3 hours.
    fn default() -> Self {
        Self {
            order: 4,
            e_folding: 3.0 * 3600.0,
        }
    }
}

impl HyperdiffusionOptions {
    /// The strongest `order` hyperdiffusion that lagged leapfrog stepping at
    /// `dt` handles with a 2× margin: `dt · max_rate = ½`. Since stable steps
    /// scale with the grid spacing, so does this e-folding time — the
    /// grid-scale advective time, the natural timescale for grid-scale
    /// dissipation.
    pub fn strongest_for_step(ops: &Operators, order: u32, dt: f64) -> Self {
        let unit = Hyperdiffusion::new(
            ops,
            Self {
                order,
                e_folding: 1.0,
            },
        );
        // max_rate scales as 1/τ.
        Self {
            order,
            e_folding: 2.0 * dt * unit.max_rate(),
        }
    }
}

/// Power iterations for the largest Laplacian eigenvalues.
const POWER_ITERATIONS: usize = 200;

/// `∇²ⁿ` hyperdiffusion on one mesh.
#[derive(Clone, Debug)]
pub struct Hyperdiffusion {
    options: HyperdiffusionOptions,
    /// Largest eigenvalues of `−L_c` and `−L_v` (m⁻²).
    lambda: [f64; 2],
    /// Coefficient `ν` (m²ⁿ s⁻¹).
    nu: f64,
    radius: f64,
    work: Work,
}

#[derive(Clone, Debug)]
struct Work {
    cells: Vec<f64>,
    vertices: Vec<f64>,
    edges: Vec<f64>,
    edges2: Vec<f64>,
    ps_edge: Vec<f64>,
}

impl Hyperdiffusion {
    /// Hyperdiffusion on the mesh of `ops`.
    ///
    /// # Panics
    /// If the order is odd or below 2, or the e-folding time isn't positive.
    pub fn new(ops: &Operators, options: HyperdiffusionOptions) -> Self {
        assert!(
            options.order >= 2 && options.order.is_multiple_of(2),
            "hyperdiffusion order must be even and at least 2"
        );
        assert!(options.e_folding > 0.0, "e-folding time must be positive");
        let m = ops.mesh();
        let (nc, ne, nv) = (m.n_cells(), m.n_edges(), m.n_vertices());
        let mut work = Work {
            cells: vec![0.0; nc],
            vertices: vec![0.0; nv],
            edges: vec![0.0; ne],
            edges2: vec![0.0; ne],
            ps_edge: vec![0.0; ne],
        };
        let lambda = [
            largest_eigenvalue(ops, Grid::Cells, &mut work),
            largest_eigenvalue(ops, Grid::Vertices, &mut work),
        ];
        let n = options.order / 2;
        let nu = 1.0 / (options.e_folding * lambda[0].powi(n as i32));
        Self {
            options,
            lambda,
            nu,
            radius: ops.radius(),
            work,
        }
    }

    /// The settings.
    pub fn options(&self) -> HyperdiffusionOptions {
        self.options
    }

    /// Largest eigenvalues of the cell and vertex Laplacians (m⁻²).
    pub fn largest_eigenvalues(&self) -> [f64; 2] {
        self.lambda
    }

    /// The coefficient `ν` (m²ⁿ s⁻¹).
    pub fn coefficient(&self) -> f64 {
        self.nu
    }

    /// Damping rate (s⁻¹) of a wave of spherical-harmonic degree `l`.
    pub fn rate_at_degree(&self, l: u32) -> f64 {
        let k2 = f64::from(l) * f64::from(l + 1) / (self.radius * self.radius);
        self.nu * k2.powi((self.options.order / 2) as i32)
    }

    /// The fastest damping rate of any mode (s⁻¹): the zig-zag vertex modes.
    /// Explicit diffusion lagged over a leapfrog step is stable for
    /// `Δt · max_rate < 1`; RK4 for `Δt · max_rate < 2.78`.
    pub fn max_rate(&self) -> f64 {
        let n = (self.options.order / 2) as i32;
        self.nu * self.lambda[0].max(self.lambda[1]).powi(n)
    }

    /// Adds the diffusion tendencies of `state` to `out` (every layer's winds
    /// and `Θ`; `out.ps` is untouched).
    pub fn add_tendency(&mut self, ops: &Operators, state: &LayeredState, out: &mut LayeredState) {
        let n = self.options.order / 2;
        let sign = if n % 2 == 1 { 1.0 } else { -1.0 };
        let nu = self.nu;
        let Work {
            cells,
            vertices,
            edges,
            edges2,
            ps_edge,
        } = &mut self.work;
        ops.cell_to_edge(&state.ps, ps_edge);
        for k in 0..state.layers() {
            // Wind: divergent part into `edges`, rotational part into `edges2`.
            ops.divergence(&state.u[k], cells);
            for _ in 1..n {
                cell_laplacian(ops, cells, edges);
            }
            ops.gradient(cells, edges);
            ops.curl(&state.u[k], vertices);
            for _ in 1..n {
                vertex_laplacian(ops, vertices, edges2);
            }
            ops.tangential_gradient(vertices, edges2);
            let (d, r) = (&*edges, &*edges2);
            out.u[k]
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(e, o)| *o += sign * nu * (d[e] - r[e]));

            // θ, in flux form.
            cells
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, c)| *c = state.theta_ps[k][i] / state.ps[i]);
            for _ in 1..n {
                cell_laplacian(ops, cells, edges);
            }
            ops.gradient(cells, edges);
            edges
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .zip(ps_edge.par_iter())
                .for_each(|(f, &p)| *f *= p);
            ops.divergence(edges, cells);
            out.theta_ps[k]
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .zip(cells.par_iter())
                .for_each(|(o, &c)| *o += sign * nu * c);
        }
    }
}

#[derive(Clone, Copy)]
enum Grid {
    Cells,
    Vertices,
}

/// `x ← L_c x`, using `edges` as scratch.
fn cell_laplacian(ops: &Operators, x: &mut [f64], edges: &mut [f64]) {
    ops.gradient(x, edges);
    ops.divergence(edges, x);
}

/// `x ← L_v x = −curl(tangential_gradient(x))`, using `edges` as scratch.
fn vertex_laplacian(ops: &Operators, x: &mut [f64], edges: &mut [f64]) {
    ops.tangential_gradient(x, edges);
    ops.curl(edges, x);
    x.par_iter_mut()
        .with_min_len(MIN_CHUNK)
        .for_each(|v| *v = -*v);
}

/// Largest eigenvalue of `−L` on cells or vertices, by power iteration with
/// the area-weighted Rayleigh quotient (`−L` is self-adjoint and positive
/// semi-definite in that inner product). Deterministic: a fixed start vector
/// and fixed-order sums.
fn largest_eigenvalue(ops: &Operators, grid: Grid, work: &mut Work) -> f64 {
    let g = ops.geometry();
    let r2 = ops.radius() * ops.radius();
    let (areas, x): (Vec<f64>, &mut Vec<f64>) = match grid {
        Grid::Cells => (
            g.cell_areas.iter().map(|a| a * r2).collect(),
            &mut work.cells,
        ),
        Grid::Vertices => (
            g.vertex_areas.iter().map(|a| a * r2).collect(),
            &mut work.vertices,
        ),
    };
    // A fixed pseudo-random start (xorshift), so every grid-scale mode is
    // present.
    let mut s = 0x9E37_79B9_7F4A_7C15u64;
    for v in x.iter_mut() {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        *v = (s >> 11) as f64 / (1u64 << 52) as f64 - 1.0;
    }
    let mut y = x.clone();
    let mut lambda = 0.0;
    for _ in 0..POWER_ITERATIONS {
        y.copy_from_slice(x);
        match grid {
            Grid::Cells => cell_laplacian(ops, &mut y, &mut work.edges),
            Grid::Vertices => vertex_laplacian(ops, &mut y, &mut work.edges),
        }
        // y = L x; Rayleigh quotient of −L.
        let weighted: Vec<f64> = x.iter().zip(&areas).map(|(v, a)| v * a).collect();
        let norm2 = dot(&weighted, x);
        lambda = -dot(&weighted, &y) / norm2;
        let scale = 1.0 / dot(&areas, &y.iter().map(|v| v * v).collect::<Vec<_>>()).sqrt();
        x.iter_mut().zip(&y).for_each(|(x, y)| *x = -y * scale);
    }
    lambda
}
