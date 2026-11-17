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
use fms_grid::operators::{for_each_mut, splits};
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

    /// Solve `(I − c∇²) x = b`, starting from the guess already in `x`.
    ///
    /// Each pass over the cells does all of an iteration's work that needs
    /// no global sum first, and accumulates the sums it feeds, each in
    /// [`dot`]'s order: so it gives what separate passes and dots give.
    // The fused loops index several arrays by one cell index.
    #[allow(clippy::needless_range_loop)]
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
        let (ops, c, areas, inverse_diagonal) = (*ops, *c, &*areas, &*inverse_diagonal);
        // Residual of the symmetric system A(I − c∇²)x = A b, preconditioned,
        // the first search direction, and |A b|², r·z and r·r.
        ops.gradient(x, edge);
        let [bb, mut rz, rr] = {
            let (x, edge) = (&*x, &*edge);
            fused3(r, z, p, |o, r, z, p| {
                let mut s = [neutral(); 3];
                for j in 0..r.len() {
                    let i = o + j;
                    let ax = areas[i] * (x[i] - c * ops.divergence_at(i, |e| edge[e]));
                    let ab = areas[i] * b[i];
                    r[j] = ab - ax;
                    z[j] = inverse_diagonal[i] * r[j];
                    p[j] = z[j];
                    s[0] += ab * ab;
                    s[1] += r[j] * z[j];
                    s[2] += r[j] * r[j];
                }
                s
            })
        };
        let b_norm = bb.sqrt();
        if b_norm == 0.0 {
            x.fill(0.0);
            return SolveReport {
                iterations: 0,
                relative_residual: 0.0,
            };
        }
        let mut residual = rr.sqrt() / b_norm;
        let mut iterations = 0;
        while residual > *tolerance && iterations < *max_iterations {
            // ap = A(I − c∇²) p, and p·ap.
            ops.gradient(p, edge);
            let [pap] = {
                let (p, edge) = (&*p, &*edge);
                fused1(ap, |o, ap| {
                    let mut s = [neutral()];
                    for j in 0..ap.len() {
                        let i = o + j;
                        ap[j] = areas[i] * (p[i] - c * ops.divergence_at(i, |e| edge[e]));
                        s[0] += p[i] * ap[j];
                    }
                    s
                })
            };
            let alpha = rz / pap;
            // The step, the new residual preconditioned, r·z and r·r.
            let [rz_next, rr] = {
                let (p, ap) = (&*p, &*ap);
                fused3(x, r, z, |o, x, r, z| {
                    let mut s = [neutral(); 2];
                    for j in 0..x.len() {
                        let i = o + j;
                        x[j] += alpha * p[i];
                        r[j] -= alpha * ap[i];
                        z[j] = inverse_diagonal[i] * r[j];
                        s[0] += r[j] * z[j];
                        s[1] += r[j] * r[j];
                    }
                    s
                })
            };
            let beta = rz_next / rz;
            rz = rz_next;
            let z = &*z;
            for_each_mut(p, |i, p| *p = z[i] + beta * *p);
            residual = rr.sqrt() / b_norm;
            iterations += 1;
        }
        SolveReport {
            iterations,
            relative_residual: residual,
        }
    }
}

/// What a sum of no terms is: where [`dot`]'s sums (`Iterator::sum`) start.
#[inline]
fn neutral() -> f64 {
    std::iter::empty::<f64>().sum()
}

/// Over `a`'s consecutive [`DOT_CHUNK`]-element chunks (side by side when
/// its length [`splits`]), `body(offset, chunk)`, which returns the chunk's
/// `K` partial sums; returns each summed over the chunks in order. With
/// partial sums taken left to right from [`neutral`], each total is summed
/// exactly as [`dot`] sums.
fn fused1<const K: usize>(
    a: &mut [f64],
    body: impl Fn(usize, &mut [f64]) -> [f64; K] + Sync + Send,
) -> [f64; K] {
    let chunks = a.len().div_ceil(DOT_CHUNK);
    let mut totals = [neutral(); K];
    let mut add = |s: [f64; K]| {
        for (t, s) in totals.iter_mut().zip(s) {
            *t += s;
        }
    };
    if !splits(a.len()) {
        for (c, a) in a.chunks_mut(DOT_CHUNK).enumerate() {
            add(body(c * DOT_CHUNK, a));
        }
    } else if chunks <= STACK_CHUNKS {
        let mut sums = [[0.0; K]; STACK_CHUNKS];
        sums[..chunks]
            .par_iter_mut()
            .zip(a.par_chunks_mut(DOT_CHUNK))
            .enumerate()
            .for_each(|(c, (s, a))| *s = body(c * DOT_CHUNK, a));
        sums[..chunks].iter().for_each(|&s| add(s));
    } else {
        let sums: Vec<[f64; K]> = a
            .par_chunks_mut(DOT_CHUNK)
            .enumerate()
            .map(|(c, a)| body(c * DOT_CHUNK, a))
            .collect();
        sums.into_iter().for_each(add);
    }
    totals
}

/// [`fused1`] over three slices of the same length at once.
fn fused3<const K: usize>(
    a: &mut [f64],
    b: &mut [f64],
    c: &mut [f64],
    body: impl Fn(usize, &mut [f64], &mut [f64], &mut [f64]) -> [f64; K] + Sync + Send,
) -> [f64; K] {
    assert!(a.len() == b.len() && a.len() == c.len());
    let chunks = a.len().div_ceil(DOT_CHUNK);
    let mut totals = [neutral(); K];
    let mut add = |s: [f64; K]| {
        for (t, s) in totals.iter_mut().zip(s) {
            *t += s;
        }
    };
    if !splits(a.len()) {
        for (i, ((a, b), c)) in a
            .chunks_mut(DOT_CHUNK)
            .zip(b.chunks_mut(DOT_CHUNK))
            .zip(c.chunks_mut(DOT_CHUNK))
            .enumerate()
        {
            add(body(i * DOT_CHUNK, a, b, c));
        }
    } else if chunks <= STACK_CHUNKS {
        let mut sums = [[0.0; K]; STACK_CHUNKS];
        sums[..chunks]
            .par_iter_mut()
            .zip(
                a.par_chunks_mut(DOT_CHUNK)
                    .zip(b.par_chunks_mut(DOT_CHUNK))
                    .zip(c.par_chunks_mut(DOT_CHUNK)),
            )
            .enumerate()
            .for_each(|(i, (s, ((a, b), c)))| *s = body(i * DOT_CHUNK, a, b, c));
        sums[..chunks].iter().for_each(|&s| add(s));
    } else {
        let sums: Vec<[f64; K]> = a
            .par_chunks_mut(DOT_CHUNK)
            .zip(b.par_chunks_mut(DOT_CHUNK))
            .zip(c.par_chunks_mut(DOT_CHUNK))
            .enumerate()
            .map(|(i, ((a, b), c))| body(i * DOT_CHUNK, a, b, c))
            .collect();
        sums.into_iter().for_each(add);
    }
    totals
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
