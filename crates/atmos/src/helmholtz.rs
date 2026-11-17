//! Solver for the Helmholtz equation `(I − c ∇²) x = b` on cells, which the
//! semi-implicit time step reduces to.
//!
//! `∇²` is the TRiSK Laplacian `divergence(gradient(·))`. Multiplied by the
//! cell areas the system is symmetric positive definite (gradient and
//! divergence are adjoint), so it is solved by conjugate gradient with a
//! diagonal (Jacobi) preconditioner.
//!
//! Every dot product sums fixed-size chunks in a fixed order, so the solution
//! and the iteration count are identical on any thread count.

use fms_grid::Operators;
use fms_grid::operators::{for_each_mut, for_each_mut2, splits};
use rayon::prelude::*;

/// Convergence report of one solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolveReport {
    /// Iterations taken.
    pub iterations: u32,
    /// Final residual norm relative to the right-hand side's.
    pub relative_residual: f64,
}

/// A Helmholtz solver for one mesh and coefficient `c`.
pub struct Helmholtz<'a> {
    ops: Operators<'a>,
    c: f64,
    /// Relative residual at which to stop.
    pub tolerance: f64,
    /// Iteration limit.
    pub max_iterations: u32,
    /// Unit-sphere cell areas (the symmetrizing weights).
    areas: Vec<f64>,
    /// Inverse of the matrix diagonal.
    inverse_diagonal: Vec<f64>,
    // Scratch.
    r: Vec<f64>,
    z: Vec<f64>,
    p: Vec<f64>,
    ap: Vec<f64>,
    edge: Vec<f64>,
}

impl<'a> Helmholtz<'a> {
    /// A solver for `(I − c ∇²) x = b`, with `c ≥ 0` in m².
    pub fn new(ops: Operators<'a>, c: f64) -> Self {
        assert!(c >= 0.0, "c must be non-negative");
        let (m, g) = (ops.mesh(), ops.geometry());
        let r2 = ops.radius() * ops.radius();
        let areas = g.cell_areas.clone();
        // Diagonal of A(I − c∇²): Aᵢ + (c/r²) Σₑ lₑ/dₑ over the cell's edges.
        let inverse_diagonal = (0..m.n_cells() as u32)
            .map(|i| {
                let coupling: f64 = m
                    .cell_edges(i)
                    .iter()
                    .map(|&e| g.edge_lengths[e as usize] / g.dual_edge_lengths[e as usize])
                    .sum();
                1.0 / (areas[i as usize] + c / r2 * coupling)
            })
            .collect();
        let (nc, ne) = (m.n_cells(), m.n_edges());
        Self {
            ops,
            c,
            tolerance: 1e-12,
            max_iterations: 500,
            areas,
            inverse_diagonal,
            r: vec![0.0; nc],
            z: vec![0.0; nc],
            p: vec![0.0; nc],
            ap: vec![0.0; nc],
            edge: vec![0.0; ne],
        }
    }

    /// The coefficient `c` (m²).
    pub fn coefficient(&self) -> f64 {
        self.c
    }

    /// `out = A (I − c∇²) x`, the symmetric form of the operator.
    fn apply(ops: &Operators, c: f64, areas: &[f64], edge: &mut [f64], x: &[f64], out: &mut [f64]) {
        ops.gradient(x, edge);
        ops.divergence(edge, out);
        for_each_mut(out, |i, o| *o = areas[i] * (x[i] - c * *o));
    }

    /// Solve `(I − c∇²) x = b`, starting from the guess already in `x`.
    pub fn solve(&mut self, b: &[f64], x: &mut [f64]) -> SolveReport {
        let Self {
            ops,
            c,
            tolerance,
            max_iterations,
            areas,
            inverse_diagonal,
            r,
            z,
            p,
            ap,
            edge,
        } = self;
        // Residual of the symmetric system A(I − c∇²)x = A b.
        Self::apply(ops, *c, areas, edge, x, ap);
        for_each_mut(r, |i, r| *r = areas[i] * b[i] - ap[i]);
        // |A b|, as dot(A b, A b) without storing A b.
        let b_norm = sum_of(b.len(), |i| {
            let ab = areas[i] * b[i];
            ab * ab
        })
        .sqrt();
        if b_norm == 0.0 {
            x.fill(0.0);
            return SolveReport {
                iterations: 0,
                relative_residual: 0.0,
            };
        }
        precondition(inverse_diagonal, r, z);
        p.copy_from_slice(z);
        let mut rz = dot(r, z);
        let mut residual = dot(r, r).sqrt() / b_norm;
        let mut iterations = 0;
        while residual > *tolerance && iterations < *max_iterations {
            Self::apply(ops, *c, areas, edge, p, ap);
            let alpha = rz / dot(p, ap);
            for_each_mut2(x, r, |i, x, r| {
                *x += alpha * p[i];
                *r -= alpha * ap[i];
            });
            precondition(inverse_diagonal, r, z);
            let rz_next = dot(r, z);
            let beta = rz_next / rz;
            rz = rz_next;
            for_each_mut(p, |i, p| *p = z[i] + beta * *p);
            residual = dot(r, r).sqrt() / b_norm;
            iterations += 1;
        }
        SolveReport {
            iterations,
            relative_residual: residual,
        }
    }
}

fn precondition(inverse_diagonal: &[f64], r: &[f64], z: &mut [f64]) {
    for_each_mut(z, |i, z| *z = inverse_diagonal[i] * r[i]);
}

/// Elements per partial sum of [`dot`]. Results depend on it: changing it
/// changes every solve in the last bits.
const DOT_CHUNK: usize = 4096;

/// Most chunk sums [`dot`] keeps on the stack (arrays up to 256k
/// elements); longer ones collect them in a `Vec`.
const STACK_CHUNKS: usize = 64;

/// Dot product with a fixed summation order: chunk sums (in parallel when
/// the arrays are long enough to split), then the chunk sums in order. A
/// short array is summed on the calling thread without allocating.
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    chunked(a.len(), |lo, hi| {
        a[lo..hi]
            .iter()
            .zip(&b[lo..hi])
            .map(|(x, y)| x * y)
            .sum::<f64>()
    })
}

/// `Σᵢ term(i)` for `i` in `0..len`, in [`dot`]'s order: each chunk of
/// [`DOT_CHUNK`] terms summed left to right, then the chunk sums in order.
/// So `sum_of(n, |i| a[i] * b[i])` is bit for bit `dot(a, b)`, without
/// storing the factors.
pub fn sum_of(len: usize, term: impl Fn(usize) -> f64 + Sync + Send) -> f64 {
    chunked(len, |lo, hi| (lo..hi).map(&term).sum::<f64>())
}

/// The sum of `chunk(lo, hi)` over consecutive [`DOT_CHUNK`]-element
/// ranges of `0..len`, in order; the chunks computed side by side only when
/// the length [`splits`].
fn chunked(len: usize, chunk: impl Fn(usize, usize) -> f64 + Sync + Send) -> f64 {
    let chunk = |c: usize| chunk(c * DOT_CHUNK, ((c + 1) * DOT_CHUNK).min(len));
    let chunks = len.div_ceil(DOT_CHUNK);
    if !splits(len) {
        (0..chunks).map(chunk).sum()
    } else if chunks <= STACK_CHUNKS {
        let mut sums = [0.0; STACK_CHUNKS];
        sums[..chunks]
            .par_iter_mut()
            .enumerate()
            .for_each(|(c, s)| *s = chunk(c));
        sums[..chunks].iter().sum()
    } else {
        (0..chunks)
            .into_par_iter()
            .map(chunk)
            .collect::<Vec<f64>>()
            .iter()
            .sum()
    }
}
