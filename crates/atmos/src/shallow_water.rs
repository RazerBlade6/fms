//! Single-layer shallow-water model on the sphere: the TRiSK scheme of
//! Ringler et al. (2010), stepped explicitly with fourth-order Runge–Kutta.
//!
//! Prognostic variables are the fluid thickness `h` at cells and the normal
//! velocity `u` at edges. In vector-invariant form:
//!
//! ```text
//! ∂h/∂t = −∇·(h u)
//! ∂u/∂t = −(ζ + f) k × u − ∇(g (h + b) + K)
//! ```
//!
//! with `b` the bottom topography, `ζ` relative vorticity, `f` the Coriolis
//! parameter and `K = |u|²/2`. Discretely: the thickness flux uses the
//! two-cell mean of `h` at edges, the Coriolis term is the energy-conserving
//! potential-vorticity flux ([`Operators::pv_flux`]) with `q = (ζ + f)/h` at
//! vertices averaged to edges, and `K` is the Ringler kinetic energy. With
//! these choices the semi-discrete equations conserve mass exactly and total
//! energy exactly (time stepping adds a small dissipation).

use fms_grid::Operators;
use fms_grid::operators::MIN_CHUNK;
use rayon::prelude::*;

/// The prognostic state.
#[derive(Clone, Debug, PartialEq)]
pub struct State {
    /// Fluid thickness at cells (m).
    pub h: Vec<f64>,
    /// Normal velocity at edges, along `nₑ` (m s⁻¹).
    pub u: Vec<f64>,
}

impl State {
    fn zeros(n_cells: usize, n_edges: usize) -> Self {
        Self {
            h: vec![0.0; n_cells],
            u: vec![0.0; n_edges],
        }
    }

    /// `self = base + a · tendency`.
    fn set_axpy(&mut self, base: &State, a: f64, tendency: &State) {
        let axpy = |out: &mut [f64], x: &[f64], y: &[f64]| {
            out.par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .zip(x)
                .zip(y)
                .for_each(|((o, &x), &y)| *o = x + a * y);
        };
        axpy(&mut self.h, &base.h, &tendency.h);
        axpy(&mut self.u, &base.u, &tendency.u);
    }
}

/// Global integrals for checking conservation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Diagnostics {
    /// Total fluid volume `Σ Aᵢ hᵢ` (m³). Conserved exactly.
    pub mass: f64,
    /// Total energy per unit density, kinetic plus potential (m⁵ s⁻²).
    /// Conserved by the spatial scheme; RK4 dissipates it very slightly.
    pub energy: f64,
    /// Total potential enstrophy `Σ A_v h_v q_v² / 2` (m s⁻²). Not exactly
    /// conserved by this (energy-conserving) scheme, but should drift slowly.
    pub potential_enstrophy: f64,
}

/// Diagnostic fields from the last tendency evaluation.
#[derive(Clone, Debug, Default)]
struct Work {
    h_edge: Vec<f64>,
    flux: Vec<f64>,
    divergence: Vec<f64>,
    kinetic: Vec<f64>,
    bernoulli: Vec<f64>,
    gradient: Vec<f64>,
    vorticity: Vec<f64>,
    h_vertex: Vec<f64>,
    q_vertex: Vec<f64>,
    q_edge: Vec<f64>,
    pv_flux: Vec<f64>,
}

/// The shallow-water equations on one mesh, with their scratch buffers.
pub struct ShallowWater<'a> {
    ops: Operators<'a>,
    gravity: f64,
    /// Coriolis parameter at vertices (s⁻¹).
    coriolis: Vec<f64>,
    /// Bottom topography at cells (m).
    topography: Vec<f64>,
    work: Work,
    /// RK4 stage tendencies and the intermediate state.
    stages: [State; 4],
    stage_state: State,
}

impl<'a> ShallowWater<'a> {
    /// A model with gravity `gravity` (m s⁻²), Coriolis parameter `coriolis`
    /// at vertices (s⁻¹) and bottom topography `topography` at cells (m).
    ///
    /// # Panics
    /// If the field lengths don't match the mesh.
    pub fn new(ops: Operators<'a>, gravity: f64, coriolis: Vec<f64>, topography: Vec<f64>) -> Self {
        let m = ops.mesh();
        let (nc, ne, nv) = (m.n_cells(), m.n_edges(), m.n_vertices());
        assert_eq!(coriolis.len(), nv, "one Coriolis value per vertex");
        assert_eq!(topography.len(), nc, "one topography value per cell");
        let work = Work {
            h_edge: vec![0.0; ne],
            flux: vec![0.0; ne],
            divergence: vec![0.0; nc],
            kinetic: vec![0.0; nc],
            bernoulli: vec![0.0; nc],
            gradient: vec![0.0; ne],
            vorticity: vec![0.0; nv],
            h_vertex: vec![0.0; nv],
            q_vertex: vec![0.0; nv],
            q_edge: vec![0.0; ne],
            pv_flux: vec![0.0; ne],
        };
        Self {
            ops,
            gravity,
            coriolis,
            topography,
            work,
            stages: std::array::from_fn(|_| State::zeros(nc, ne)),
            stage_state: State::zeros(nc, ne),
        }
    }

    /// The operators.
    pub fn operators(&self) -> &Operators<'a> {
        &self.ops
    }

    /// Coriolis parameter at vertices (s⁻¹).
    pub fn coriolis(&self) -> &[f64] {
        &self.coriolis
    }

    /// Bottom topography at cells (m).
    pub fn topography(&self) -> &[f64] {
        &self.topography
    }

    /// Gravity (m s⁻²).
    pub fn gravity(&self) -> f64 {
        self.gravity
    }

    /// Time tendencies of `state` into `out`.
    pub fn tendency(&mut self, state: &State, out: &mut State) {
        compute_tendency(
            &self.ops,
            self.gravity,
            &self.coriolis,
            &self.topography,
            &mut self.work,
            state,
            out,
        );
    }

    /// Advance `state` by `dt` seconds with classic fourth-order Runge–Kutta.
    pub fn step(&mut self, state: &mut State, dt: f64) {
        let Self {
            ops,
            gravity,
            coriolis,
            topography,
            work,
            stages,
            stage_state,
        } = self;
        let mut tendency = |s: &State, out: &mut State| {
            compute_tendency(ops, *gravity, coriolis, topography, work, s, out);
        };
        let [k1, k2, k3, k4] = stages;
        tendency(state, k1);
        stage_state.set_axpy(state, 0.5 * dt, k1);
        tendency(stage_state, k2);
        stage_state.set_axpy(state, 0.5 * dt, k2);
        tendency(stage_state, k3);
        stage_state.set_axpy(state, dt, k3);
        tendency(stage_state, k4);
        let combine = |x: &mut [f64], a: &[f64], b: &[f64], c: &[f64], d: &[f64]| {
            x.par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, x)| {
                    *x += dt / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]);
                });
        };
        combine(&mut state.h, &k1.h, &k2.h, &k3.h, &k4.h);
        combine(&mut state.u, &k1.u, &k2.u, &k3.u, &k4.u);
    }

    /// Relative vorticity of `state` at vertices (s⁻¹).
    pub fn vorticity(&self, state: &State) -> Vec<f64> {
        let mut zeta = vec![0.0; self.ops.mesh().n_vertices()];
        self.ops.curl(&state.u, &mut zeta);
        zeta
    }

    /// Mass, energy and potential enstrophy of `state`.
    pub fn diagnostics(&self, state: &State) -> Diagnostics {
        let ops = &self.ops;
        let (m, g) = (ops.mesh(), ops.geometry());
        let r2 = ops.radius() * ops.radius();
        let mut kinetic = vec![0.0; m.n_cells()];
        ops.kinetic_energy(&state.u, &mut kinetic);
        let mut zeta = vec![0.0; m.n_vertices()];
        ops.curl(&state.u, &mut zeta);
        let mut h_vertex = vec![0.0; m.n_vertices()];
        ops.cell_to_vertex(&state.h, &mut h_vertex);

        let (mut mass, mut energy) = (0.0, 0.0);
        for (((&area, &h), &k), &b) in g
            .cell_areas
            .iter()
            .zip(&state.h)
            .zip(&kinetic)
            .zip(&self.topography)
        {
            let a = area * r2;
            mass += a * h;
            energy += a * (h * k + self.gravity * h * (0.5 * h + b));
        }
        let potential_enstrophy = (0..m.n_vertices())
            .map(|v| {
                let q = (zeta[v] + self.coriolis[v]) / h_vertex[v];
                g.vertex_areas[v] * r2 * h_vertex[v] * q * q / 2.0
            })
            .sum();
        Diagnostics {
            mass,
            energy,
            potential_enstrophy,
        }
    }

    /// A time step (s) at Courant number `courant` for the fastest signal,
    /// gravity waves plus advection: `courant · min dₑ / max(|u| + √(g h))`.
    /// About 0.5 is safely stable with RK4.
    pub fn stable_dt(&self, state: &State, courant: f64) -> f64 {
        let g = self.ops.geometry();
        let min_d = g.dual_edge_lengths.iter().fold(f64::MAX, |a, &b| a.min(b));
        let max_c = state.h.iter().fold(0.0f64, |a, &h| a.max(h)) * self.gravity;
        let max_u = state.u.iter().fold(0.0f64, |a, &u| a.max(u.abs()));
        courant * min_d * self.ops.radius() / (max_u + max_c.sqrt())
    }
}

fn compute_tendency(
    ops: &Operators,
    gravity: f64,
    coriolis: &[f64],
    topography: &[f64],
    w: &mut Work,
    state: &State,
    out: &mut State,
) {
    // Continuity: thickness flux with the two-cell mean thickness.
    ops.cell_to_edge(&state.h, &mut w.h_edge);
    w.flux
        .par_iter_mut()
        .with_min_len(MIN_CHUNK)
        .zip(&w.h_edge)
        .zip(&state.u)
        .for_each(|((f, &h), &u)| *f = h * u);
    ops.divergence(&w.flux, &mut w.divergence);
    out.h
        .par_iter_mut()
        .with_min_len(MIN_CHUNK)
        .zip(&w.divergence)
        .for_each(|(dh, &d)| *dh = -d);

    // Bernoulli function g (h + b) + K and its gradient.
    ops.kinetic_energy(&state.u, &mut w.kinetic);
    w.bernoulli
        .par_iter_mut()
        .with_min_len(MIN_CHUNK)
        .enumerate()
        .for_each(|(i, b)| *b = gravity * (state.h[i] + topography[i]) + w.kinetic[i]);
    ops.gradient(&w.bernoulli, &mut w.gradient);

    // Potential vorticity q = (ζ + f) / h at vertices, then at edges.
    ops.curl(&state.u, &mut w.vorticity);
    ops.cell_to_vertex(&state.h, &mut w.h_vertex);
    w.q_vertex
        .par_iter_mut()
        .with_min_len(MIN_CHUNK)
        .enumerate()
        .for_each(|(v, q)| *q = (w.vorticity[v] + coriolis[v]) / w.h_vertex[v]);
    ops.vertex_to_edge(&w.q_vertex, &mut w.q_edge);
    ops.pv_flux(&w.q_edge, &w.flux, &mut w.pv_flux);

    out.u
        .par_iter_mut()
        .with_min_len(MIN_CHUNK)
        .zip(&w.pv_flux)
        .zip(&w.gradient)
        .for_each(|((du, &q), &g)| *du = q - g);
}
