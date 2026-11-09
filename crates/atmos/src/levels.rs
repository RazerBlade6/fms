//! The vertical grid: `N` layers in `σ = p/p_s`, top first, and the small
//! dense linear algebra of the vertical modes.
//!
//! Layer `k` spans `[σ_{k−½}, σ_{k+½}]`, from `σ_{½} = 0` at the model top
//! to `σ_{N+½} = 1` at the surface. Its pressure thickness is
//! `δpₖ = δσₖ p_s`, and its Exner function is the mass-weighted mean of
//! `Π = (p/p₀)^κ` over the layer (Arakawa & Suarez 1983):
//!
//! ```text
//! Pₖ = (p_{k+½} Π_{k+½} − p_{k−½} Π_{k−½}) / ((1 + κ) δpₖ)
//! ```

// Small dense matrices and per-layer arrays are indexed together.
#![allow(clippy::needless_range_loop)]

use std::sync::OnceLock;

use crate::layered::Gas;
use crate::radiation::ProfileShape;

/// Most layers a model may have (bounds the per-column scratch arrays).
pub const MAX_LAYERS: usize = 16;

/// The σ levels of the layer interfaces.
#[derive(Clone, Debug, PartialEq)]
pub struct Levels {
    /// All interfaces, `0` at the top to `1` at the surface (`N + 1` values).
    interfaces: Vec<f64>,
    /// `δσₖ` per layer.
    dsigma: Vec<f64>,
    profile: ProfileCache,
}

/// Memoized constants of the in-layer temperature profile
/// ([`crate::radiation`]). Not part of the levels' identity.
#[derive(Clone, Debug, Default)]
struct ProfileCache(OnceLock<ProfileShape>);

impl PartialEq for ProfileCache {
    fn eq(&self, _: &Self) -> bool {
        true
    }
}

impl Levels {
    /// Layers with interior interfaces `interior` (strictly increasing, in
    /// (0, 1)); `&[]` is a single layer.
    ///
    /// # Panics
    /// If the interfaces aren't strictly increasing within (0, 1), or there
    /// would be more than [`MAX_LAYERS`] layers.
    pub fn new(interior: &[f64]) -> Self {
        assert!(interior.len() < MAX_LAYERS, "at most {MAX_LAYERS} layers");
        let mut interfaces = Vec::with_capacity(interior.len() + 2);
        interfaces.push(0.0);
        interfaces.extend_from_slice(interior);
        interfaces.push(1.0);
        assert!(
            interfaces.windows(2).all(|w| w[0] < w[1]),
            "interfaces must increase strictly within (0, 1): {interior:?}"
        );
        let dsigma = interfaces.windows(2).map(|w| w[1] - w[0]).collect();
        Self {
            interfaces,
            dsigma,
            profile: ProfileCache::default(),
        }
    }

    /// Two layers split at `sigma`.
    pub fn two(sigma: f64) -> Self {
        Self::new(&[sigma])
    }

    /// `n` layers of equal mass.
    ///
    /// # Panics
    /// If `n` is 0 or more than [`MAX_LAYERS`].
    pub fn equal_mass(n: usize) -> Self {
        assert!(n >= 1, "at least one layer");
        let interior: Vec<f64> = (1..n).map(|k| k as f64 / n as f64).collect();
        Self::new(&interior)
    }

    /// The memoized profile constants, computed by `init` on first use.
    pub(crate) fn profile_shape(&self, init: impl FnOnce() -> ProfileShape) -> &ProfileShape {
        self.profile.0.get_or_init(init)
    }

    /// Number of layers `N`.
    pub fn count(&self) -> usize {
        self.dsigma.len()
    }

    /// All `N + 1` interfaces, top (0) to surface (1).
    pub fn interfaces(&self) -> &[f64] {
        &self.interfaces
    }

    /// The `N − 1` interior interfaces; interface `j` lies between layers
    /// `j` and `j + 1`.
    pub fn interior(&self) -> &[f64] {
        &self.interfaces[1..self.interfaces.len() - 1]
    }

    /// Layer thicknesses `δσₖ`.
    pub fn dsigma(&self) -> &[f64] {
        &self.dsigma
    }

    /// Top and bottom σ of layer `k`.
    pub fn bounds(&self, k: usize) -> (f64, f64) {
        (self.interfaces[k], self.interfaces[k + 1])
    }

    /// The layer containing level `sigma` (an interface belongs to the layer
    /// above it).
    pub fn layer_of(&self, sigma: f64) -> usize {
        (0..self.count())
            .find(|&k| sigma <= self.interfaces[k + 1])
            .unwrap_or(self.count() - 1)
    }

    /// Layer Exner functions `Pₖ` for surface pressure `ps`, into `out`
    /// (one per layer).
    pub fn layer_exner_into(&self, gas: Gas, ps: f64, out: &mut [f64]) {
        let kappa = gas.kappa();
        let n = self.count();
        // p Π at the top of the current layer (0 at the model top).
        let mut top = 0.0;
        for k in 0..n {
            let s = self.interfaces[k + 1];
            let p = s * ps;
            let exner = gas.exner(p);
            let bottom = p * exner;
            out[k] = if k == 0 {
                // p_{½} = 0: P₀ = Π_{3/2} / (1 + κ).
                exner / (1.0 + kappa)
            } else {
                (bottom - top) / ((1.0 + kappa) * self.dsigma[k] * ps)
            };
            top = bottom;
        }
    }

    /// Layer Exner functions `Pₖ` for surface pressure `ps`.
    pub fn layer_exner(&self, gas: Gas, ps: f64) -> Vec<f64> {
        let mut out = vec![0.0; self.count()];
        self.layer_exner_into(gas, ps, &mut out);
        out
    }

    /// The pressure (Pa) whose Exner function is each layer's mean.
    pub fn mean_pressures(&self, gas: Gas, ps: f64) -> Vec<f64> {
        self.layer_exner(gas, ps)
            .into_iter()
            .map(|e| crate::layered::REFERENCE_PRESSURE * fms_math::pow(e, 1.0 / gas.kappa()))
            .collect()
    }
}

/// Eigen-decomposition `A = V diag(λ) V⁻¹` of a small real matrix with real
/// eigenvalues (the gravity-wave matrix), the eigenvalues in increasing
/// order.
///
/// With more than one interface the gravity-wave matrix isn't self-adjoint
/// (the arithmetic-mean interface θ of the vertical scheme makes it slightly
/// non-normal), so this is the general decomposition. Eigenvalues may repeat:
/// the internal modes of an unstratified atmosphere all have zero speed.
#[derive(Clone, Debug, PartialEq)]
pub struct Eigen {
    /// Eigenvalues, increasing.
    pub values: Vec<f64>,
    /// Columns are the eigenvectors (unit length).
    pub vectors: Vec<Vec<f64>>,
    /// `V⁻¹`.
    pub inverse: Vec<Vec<f64>>,
}

impl Eigen {
    /// Decompose `a` (square, rows). One- and two-layer matrices use closed
    /// forms; larger ones the shifted QR algorithm for the eigenvalues and
    /// (subspace) inverse iteration for the vectors.
    ///
    /// # Panics
    /// If the eigenvalues are complex or `a` isn't diagonalizable.
    pub fn of(a: &[Vec<f64>]) -> Self {
        let n = a.len();
        if n == 1 {
            return Self {
                values: vec![a[0][0]],
                vectors: vec![vec![1.0]],
                inverse: vec![vec![1.0]],
            };
        }
        if n == 2 {
            return Self::of_2x2([[a[0][0], a[0][1]], [a[1][0], a[1][1]]]);
        }
        let values = eigenvalues(a);
        let scale = a
            .iter()
            .flatten()
            .fold(0.0f64, |m, &x| m.max(x.abs()))
            .max(1e-100);
        // Groups of equal eigenvalues, each with a basis of its eigenspace.
        let mut columns: Vec<Vec<f64>> = Vec::with_capacity(n);
        let mut start = 0;
        while start < n {
            let mut end = start + 1;
            while end < n && values[end] - values[start] <= 1e-10 * scale {
                end += 1;
            }
            let mean = values[start..end].iter().sum::<f64>() / (end - start) as f64;
            columns.extend(eigenspace(a, mean, end - start, start, scale));
            start = end;
        }
        let vectors: Vec<Vec<f64>> = (0..n)
            .map(|r| columns.iter().map(|c| c[r]).collect())
            .collect();
        let inverse = invert(&vectors).expect("the gravity-wave matrix must be diagonalizable");
        Self {
            values,
            vectors,
            inverse,
        }
    }

    fn of_2x2(a: [[f64; 2]; 2]) -> Self {
        let (tr, det) = (a[0][0] + a[1][1], a[0][0] * a[1][1] - a[0][1] * a[1][0]);
        let disc = (0.25 * tr * tr - det).max(0.0).sqrt();
        let lambda = [0.5 * tr - disc, 0.5 * tr + disc];
        let vector = |l: f64, fallback: [f64; 2]| -> [f64; 2] {
            let (x, y) = if a[0][1].abs() >= a[1][0].abs() && a[0][1] != 0.0 {
                (a[0][1], l - a[0][0])
            } else if a[1][0] != 0.0 {
                (l - a[1][1], a[1][0])
            } else {
                (fallback[0], fallback[1])
            };
            let n = (x * x + y * y).sqrt();
            [x / n, y / n]
        };
        let (e0, e1) = (vector(lambda[0], [1.0, 0.0]), vector(lambda[1], [0.0, 1.0]));
        let v = [[e0[0], e1[0]], [e0[1], e1[1]]];
        let d = v[0][0] * v[1][1] - v[0][1] * v[1][0];
        let v_inv = [[v[1][1] / d, -v[0][1] / d], [-v[1][0] / d, v[0][0] / d]];
        Self {
            values: lambda.to_vec(),
            vectors: v.iter().map(|r| r.to_vec()).collect(),
            inverse: v_inv.iter().map(|r| r.to_vec()).collect(),
        }
    }
}

/// Real eigenvalues of `a` by the shifted QR algorithm with deflation,
/// increasing.
///
/// # Panics
/// If they don't converge (complex eigenvalues).
fn eigenvalues(a: &[Vec<f64>]) -> Vec<f64> {
    let mut h: Vec<Vec<f64>> = a.to_vec();
    let mut n = h.len();
    let scale = a
        .iter()
        .flatten()
        .fold(0.0f64, |m, &x| m.max(x.abs()))
        .max(f64::MIN_POSITIVE);
    let mut values = Vec::with_capacity(n);
    while n > 1 {
        let mut converged = false;
        for _ in 0..1000 {
            // The last row's off-diagonal part (the matrix isn't Hessenberg).
            let off: f64 = (0..n - 1).map(|c| h[n - 1][c].abs()).sum();
            if off <= 1e-15 * scale {
                converged = true;
                break;
            }
            // Wilkinson shift: the trailing 2×2's eigenvalue nearer its
            // last diagonal entry (its real part if complex).
            let (p, q, r, s) = (
                h[n - 2][n - 2],
                h[n - 2][n - 1],
                h[n - 1][n - 2],
                h[n - 1][n - 1],
            );
            let (mid, half) = (0.5 * (p + s), 0.5 * (p - s));
            let disc = half * half + q * r;
            let mu = if disc >= 0.0 {
                let root = disc.sqrt();
                let (l1, l2) = (mid + root, mid - root);
                if (l1 - s).abs() < (l2 - s).abs() {
                    l1
                } else {
                    l2
                }
            } else {
                mid
            };
            for (i, row) in h.iter_mut().enumerate().take(n) {
                row[i] -= mu;
            }
            let (q_m, r_m) = qr(&h, n);
            // h ← R Q + μ I on the leading n × n block.
            for i in 0..n {
                for j in 0..n {
                    h[i][j] = (0..n).map(|k| r_m[i][k] * q_m[k][j]).sum::<f64>();
                }
                h[i][i] += mu;
            }
        }
        assert!(converged, "the gravity-wave matrix has complex eigenvalues");
        values.push(h[n - 1][n - 1]);
        n -= 1;
    }
    values.push(h[0][0]);
    values.sort_by(f64::total_cmp);
    values
}

/// Householder QR of the leading `n × n` block of `a`: `(Q, R)`.
fn qr(a: &[Vec<f64>], n: usize) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let mut r: Vec<Vec<f64>> = (0..n).map(|i| a[i][..n].to_vec()).collect();
    let mut q: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    for c in 0..n - 1 {
        let norm = (c..n).map(|i| r[i][c] * r[i][c]).sum::<f64>().sqrt();
        if norm == 0.0 {
            continue;
        }
        let alpha = if r[c][c] > 0.0 { -norm } else { norm };
        let mut v: Vec<f64> = (0..n).map(|i| if i < c { 0.0 } else { r[i][c] }).collect();
        v[c] -= alpha;
        let vv: f64 = v.iter().map(|x| x * x).sum();
        if vv == 0.0 {
            continue;
        }
        // R ← (I − 2vvᵀ/vᵀv) R;  Q ← Q (I − 2vvᵀ/vᵀv).
        for j in 0..n {
            let d: f64 = (c..n).map(|i| v[i] * r[i][j]).sum::<f64>() * 2.0 / vv;
            for i in c..n {
                r[i][j] -= d * v[i];
            }
        }
        for row in q.iter_mut() {
            let d: f64 = (c..n).map(|i| row[i] * v[i]).sum::<f64>() * 2.0 / vv;
            for i in c..n {
                row[i] -= d * v[i];
            }
        }
    }
    (q, r)
}

/// A basis (unit vectors) of the eigenspace of `a` for an eigenvalue
/// `lambda` of multiplicity `m`, by subspace inverse iteration from start
/// vectors chosen by `seed`.
fn eigenspace(a: &[Vec<f64>], lambda: f64, m: usize, seed: usize, scale: f64) -> Vec<Vec<f64>> {
    let n = a.len();
    // Shifted just off the eigenvalue, so the system is solvable but
    // amplifies the eigenspace by ~10¹⁰ per iteration.
    let shift = lambda + 1e-10 * scale;
    let shifted: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| a[i][j] - if i == j { shift } else { 0.0 })
                .collect()
        })
        .collect();
    let mut basis: Vec<Vec<f64>> = (0..m)
        .map(|c| {
            (0..n)
                .map(|i| 1.0 + 0.37 * ((i * 7 + (seed + c) * 13) % 11) as f64)
                .collect()
        })
        .collect();
    for _ in 0..4 {
        for v in &mut basis {
            *v = solve(&shifted, v).expect("inverse iteration needs a solvable system");
        }
        orthonormalize(&mut basis);
    }
    // A fixed sign for a single vector: the largest component positive.
    if m == 1 {
        let v = &mut basis[0];
        let big = v
            .iter()
            .copied()
            .fold(0.0f64, |b, x| if x.abs() > b.abs() { x } else { b });
        if big < 0.0 {
            v.iter_mut().for_each(|x| *x = -*x);
        }
    }
    basis
}

/// Modified Gram–Schmidt.
fn orthonormalize(vectors: &mut [Vec<f64>]) {
    for i in 0..vectors.len() {
        for j in 0..i {
            let d: f64 = vectors[i].iter().zip(&vectors[j]).map(|(a, b)| a * b).sum();
            let vj = vectors[j].clone();
            for (x, y) in vectors[i].iter_mut().zip(vj) {
                *x -= d * y;
            }
        }
        let norm = vectors[i].iter().map(|x| x * x).sum::<f64>().sqrt();
        vectors[i].iter_mut().for_each(|x| *x /= norm);
    }
}

/// Solve `m x = b` by Gaussian elimination with partial pivoting; `None` if
/// singular.
fn solve(m: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = m.len();
    let mut a: Vec<Vec<f64>> = m
        .iter()
        .zip(b)
        .map(|(row, &b)| {
            let mut r = row.clone();
            r.push(b);
            r
        })
        .collect();
    for c in 0..n {
        let p = (c..n).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs()))?;
        if a[p][c] == 0.0 {
            return None;
        }
        a.swap(c, p);
        for i in c + 1..n {
            let f = a[i][c] / a[c][c];
            for j in c..=n {
                a[i][j] -= f * a[c][j];
            }
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let s: f64 = (i + 1..n).map(|j| a[i][j] * x[j]).sum();
        x[i] = (a[i][n] - s) / a[i][i];
    }
    Some(x)
}

/// Inverse of a square matrix; `None` if singular.
fn invert(m: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = m.len();
    let columns: Option<Vec<Vec<f64>>> = (0..n)
        .map(|c| {
            let e: Vec<f64> = (0..n).map(|i| if i == c { 1.0 } else { 0.0 }).collect();
            solve(m, &e)
        })
        .collect();
    let columns = columns?;
    Some(
        (0..n)
            .map(|r| columns.iter().map(|c| c[r]).collect())
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn levels_and_exner() {
        let l = Levels::equal_mass(4);
        assert_eq!(l.count(), 4);
        assert_eq!(l.interior(), &[0.25, 0.5, 0.75]);
        assert_eq!(l.layer_of(0.1), 0);
        assert_eq!(l.layer_of(0.25), 0);
        assert_eq!(l.layer_of(0.26), 1);
        assert_eq!(l.layer_of(1.0), 3);
        // The mass-weighted Exner functions average to the column's:
        // Σ δσₖ Pₖ = Π_s / (1 + κ).
        let gas = Gas::EARTH_AIR;
        let ps = 98_000.0;
        let p = l.layer_exner(gas, ps);
        let mean: f64 = p.iter().zip(l.dsigma()).map(|(p, d)| p * d).sum();
        assert!((mean - gas.exner(ps) / (1.0 + gas.kappa())).abs() < 1e-15);
        assert!(p.windows(2).all(|w| w[0] < w[1]));
        // Two layers: the closed forms.
        let two = Levels::two(0.4);
        let p = two.layer_exner(gas, ps);
        let (pi, kappa) = (0.4 * ps, gas.kappa());
        assert_eq!(p[0], gas.exner(pi) / (1.0 + kappa));
        assert_eq!(
            p[1],
            (ps * gas.exner(ps) - pi * gas.exner(pi)) / ((1.0 + kappa) * (1.0 - 0.4) * ps)
        );
    }

    #[test]
    fn eigen_decomposition() {
        // A matrix with real eigenvalues (self-adjoint for weights w), and
        // the zero matrix (an eigenvalue of multiplicity 3).
        let w = [0.2, 0.3, 0.35, 0.15];
        let sym = [
            [-4.0, 1.0, 0.5, 0.1],
            [1.0, -3.0, 0.2, 0.3],
            [0.5, 0.2, -2.0, 0.2],
            [0.1, 0.3, 0.2, -1.0],
        ];
        let a: Vec<Vec<f64>> = (0..4)
            .map(|r| (0..4).map(|c| sym[r][c] / w[r]).collect())
            .collect();
        for a in [a, vec![vec![0.0; 3]; 3]] {
            let n = a.len();
            let e = Eigen::of(&a);
            assert!(e.values.windows(2).all(|p| p[0] <= p[1]));
            // A V = V Λ and V V⁻¹ = I, repeated eigenvalues included.
            for r in 0..n {
                for c in 0..n {
                    let av: f64 = (0..n).map(|k| a[r][k] * e.vectors[k][c]).sum();
                    assert!(
                        (av - e.vectors[r][c] * e.values[c]).abs() < 1e-12,
                        "{r} {c}"
                    );
                    let vi: f64 = (0..n).map(|k| e.vectors[r][k] * e.inverse[k][c]).sum();
                    assert!((vi - if r == c { 1.0 } else { 0.0 }).abs() < 1e-12);
                }
            }
        }
    }
}
