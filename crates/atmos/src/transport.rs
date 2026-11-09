//! Positive-definite humidity transport for the semi-implicit layered
//! stepper: flux-corrected transport (Zalesak 1979) around the leapfrog
//! step.
//!
//! # The scheme
//!
//! The leapfrog step takes `Qₖ = p_s qₖ` from level `n−1` to `n+1` over
//! `T = 2Δt`. Centred humidity fluxes (the "high-order" scheme, the same as
//! θ's) make spurious negative values and overshoots near sharp gradients;
//! upwind fluxes (the "low-order" scheme) are positive and free of new
//! extrema but smear everything. FCT takes the upwind solution and adds back
//! as much of the difference between the two fluxes, edge by edge, as it can
//! without leaving the range of the nearby values:
//!
//! 1. **Low order**: upwind from level `n−1` with the step's mass fluxes,
//!    in as many forward substeps as positivity needs (each cell may lose at
//!    most 90% of its water per substep), plus the physics sources.
//! 2. **High order**: centred at level `n` (leapfrog).
//!
//! Across each layer interface, humidity is carried at the interface values
//! of the two layers' own profiles (`rₖ qₖ`, with ratios `rₖ` from
//! [`crate::radiation::interface_humidity_ratios`]): the upwind layer's in
//! the low-order flux, the mean of the two in the high-order flux.
//! 3. **Limit** the antidiffusive flux (high minus low) on each edge, and
//!    across each interface, so that every cell stays within the range of
//!    `q` at `n−1`, `n` and in the low-order solution, over itself, its
//!    neighbours and the adjacent layers.
//!
//! The mass fluxes are those the semi-implicit step actually used: per layer
//! `F̄ₖ = Fₖⁿ + p̄ δσₖ (uₖⁿ⁻¹ + uₖⁿ⁺¹ − 2uₖⁿ)/2` (the nonlinear flux at `n`
//! with its linear part averaged over the old and new levels), and the
//! interface fluxes `W̄ⱼ` they imply (as in [`crate::layered`]). These change
//! each layer's mass from `n−1` to `n+1` exactly as the surface-pressure
//! equation did, so humidity and mass are consistent: a uniform `q` stays
//! uniform. Everything is in flux form, so total water changes only through
//! the sources.

// Loops over the layers index several per-layer arrays together.
#![allow(clippy::needless_range_loop)]

use fms_grid::Operators;
use fms_grid::operators::MIN_CHUNK;
use rayon::prelude::*;

/// Largest share of a cell's water an upwind substep may move out.
const MAX_OUTFLOW: f64 = 0.9;

/// Inputs of one transport step (`N` layers, `N − 1` interfaces).
pub struct TransportStep<'s> {
    /// Layer thickness fractions `δσₖ`.
    pub dsigma: &'s [f64],
    /// Length of the step, `2Δt` for leapfrog (s).
    pub interval: f64,
    /// Surface pressure at the start (`n−1`) and end (`n+1`) (Pa).
    pub ps_old: &'s [f64],
    /// See `ps_old`.
    pub ps_new: &'s [f64],
    /// Mass flux `F̄ₖ` per layer at edges (Pa m s⁻¹), positive from
    /// `cells[0]` to `cells[1]`.
    pub flux: Vec<&'s [f64]>,
    /// Mass flux `W̄ⱼ` across each interface (`p_s σ̇`, Pa s⁻¹), positive
    /// downward.
    pub interface: Vec<&'s [f64]>,
    /// `Qₖ` at the start (Pa).
    pub q_old: Vec<&'s [f64]>,
    /// Specific humidity `qₖ` at the middle level (kg kg⁻¹).
    pub q_mid: Vec<&'s [f64]>,
    /// Sources of `Qₖ` (Pa s⁻¹): physics, constant over the step.
    pub source: Vec<&'s [f64]>,
    /// Per interface, humidity there relative to each adjacent layer's mean:
    /// the upper layer at its bottom, the lower at its top.
    pub interface_ratio: Vec<[&'s [f64]; 2]>,
}

/// What one transport step did.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TransportReport {
    /// Upwind substeps taken.
    pub substeps: u32,
}

/// Scratch space for [`Transport::step`].
pub struct Transport {
    low_flux: Vec<Vec<f64>>,
    anti: Vec<Vec<f64>>,
    low_interface: Vec<Vec<f64>>,
    anti_interface: Vec<Vec<f64>>,
    content: Vec<Vec<f64>>,
    q_low: Vec<Vec<f64>>,
    q_max: Vec<Vec<f64>>,
    q_min: Vec<Vec<f64>>,
    r_plus: Vec<Vec<f64>>,
    r_minus: Vec<Vec<f64>>,
    edge_amount: Vec<f64>,
    divergence: Vec<f64>,
    q: Vec<Vec<f64>>,
}

impl Transport {
    /// Scratch for `layers` layers on a mesh with `n_cells` cells and
    /// `n_edges` edges.
    pub fn new(layers: usize, n_cells: usize, n_edges: usize) -> Self {
        let cells = |n: usize| vec![vec![0.0; n_cells]; n];
        let edges = || vec![vec![0.0; n_edges]; layers];
        let interfaces = layers.saturating_sub(1);
        Self {
            low_flux: edges(),
            anti: edges(),
            low_interface: cells(interfaces),
            anti_interface: cells(interfaces),
            content: cells(layers),
            q_low: cells(layers),
            q_max: cells(layers),
            q_min: cells(layers),
            r_plus: cells(layers),
            r_minus: cells(layers),
            edge_amount: vec![0.0; n_edges],
            divergence: vec![0.0; n_cells],
            q: cells(layers),
        }
    }

    /// Transport humidity over one step, writing `Qₖ` at the end into `out`.
    pub fn step(
        &mut self,
        ops: &Operators,
        input: &TransportStep,
        out: Vec<&mut [f64]>,
    ) -> TransportReport {
        let mesh = ops.mesh();
        let g = ops.geometry();
        let r = ops.radius();
        let edges = mesh.edges();
        let ds = input.dsigma;
        let n = ds.len();
        let t = input.interval;
        let w = &input.interface;
        let ratio = &input.interface_ratio;

        // Substeps: each cell may move out at most MAX_OUTFLOW of its water
        // per substep, at its smallest mass over the step.
        let outflow_rate = |c: usize, k: usize| -> f64 {
            let cid = c as u32;
            let horizontal: f64 = mesh
                .cell_edges(cid)
                .iter()
                .zip(g.cell_edge_signs(cid))
                .map(|(&e, &s)| {
                    (s * input.flux[k][e as usize]).max(0.0) * g.edge_lengths[e as usize]
                })
                .sum::<f64>()
                / (g.cell_areas[c] * r);
            // Down through the bottom, up through the top.
            let mut vertical = 0.0;
            if k + 1 < n {
                vertical += ratio[k][0][c] * w[k][c].max(0.0);
            }
            if k > 0 {
                vertical += ratio[k - 1][1][c] * (-w[k - 1][c]).max(0.0);
            }
            (horizontal + vertical) / (ds[k] * input.ps_old[c].min(input.ps_new[c]))
        };
        let max_rate = (0..mesh.n_cells())
            .into_par_iter()
            .with_min_len(MIN_CHUNK)
            .map(|c| {
                let mut m = outflow_rate(c, 0);
                for k in 1..n {
                    m = m.max(outflow_rate(c, k));
                }
                m
            })
            .reduce(|| 0.0, f64::max);
        let substeps = ((t * max_rate / MAX_OUTFLOW).ceil() as u32).max(1);
        let h = t / f64::from(substeps);

        // 1. Low order: upwind substeps from n−1, with the sources.
        for k in 0..n {
            self.content[k]
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .zip(input.q_old[k].par_iter())
                .for_each(|(c, &q)| *c = ds[k] * q);
            self.low_flux[k].fill(0.0);
        }
        for x in &mut self.low_interface {
            x.fill(0.0);
        }
        for m in 0..substeps {
            let s = f64::from(m) / f64::from(substeps);
            for k in 0..n {
                let content = &self.content[k];
                self.q[k]
                    .par_iter_mut()
                    .with_min_len(MIN_CHUNK)
                    .enumerate()
                    .for_each(|(i, q)| {
                        let ps = input.ps_old[i] + s * (input.ps_new[i] - input.ps_old[i]);
                        *q = content[i] / (ds[k] * ps);
                    });
            }
            for k in 0..n {
                let (q, flux) = (&self.q[k], input.flux[k]);
                self.edge_amount
                    .par_iter_mut()
                    .with_min_len(MIN_CHUNK)
                    .zip(self.low_flux[k].par_iter_mut())
                    .enumerate()
                    .for_each(|(e, (x, total))| {
                        let [a, b] = edges[e].cells;
                        let upwind = if flux[e] >= 0.0 {
                            q[a as usize]
                        } else {
                            q[b as usize]
                        };
                        *x = h * flux[e] * upwind;
                        *total += *x;
                    });
                ops.divergence(&self.edge_amount, &mut self.divergence);
                let (div, source) = (&self.divergence, input.source[k]);
                self.content[k]
                    .par_iter_mut()
                    .with_min_len(MIN_CHUNK)
                    .enumerate()
                    .for_each(|(i, c)| *c += -div[i] + h * ds[k] * source[i]);
            }
            // Across each interface (downward positive: from layer j to j+1).
            for j in 0..n - 1 {
                let (upper, lower) = self.content.split_at_mut(j + 1);
                let (c0, c1) = (&mut upper[j], &mut lower[0]);
                let (q, w, ratio) = (&self.q, w[j], ratio[j]);
                c0.par_iter_mut()
                    .with_min_len(MIN_CHUNK)
                    .zip(c1.par_iter_mut())
                    .zip(self.low_interface[j].par_iter_mut())
                    .enumerate()
                    .for_each(|(i, ((c0, c1), total))| {
                        let upwind = if w[i] >= 0.0 {
                            ratio[0][i] * q[j][i]
                        } else {
                            ratio[1][i] * q[j + 1][i]
                        };
                        let x = h * w[i] * upwind;
                        *c0 -= x;
                        *c1 += x;
                        *total += x;
                    });
            }
        }
        for k in 0..n {
            let content = &self.content[k];
            self.q_low[k]
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, q)| *q = content[i] / (ds[k] * input.ps_new[i]));
        }

        // 2. Antidiffusive fluxes: centred at n minus upwind.
        for k in 0..n {
            let (q, flux, low) = (input.q_mid[k], input.flux[k], &self.low_flux[k]);
            self.anti[k]
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(e, a)| {
                    let [i, j] = edges[e].cells;
                    let centred = 0.5 * (q[i as usize] + q[j as usize]);
                    *a = t * flux[e] * centred - low[e];
                });
        }
        for j in 0..n - 1 {
            let (q, w, ratio) = (&input.q_mid, w[j], ratio[j]);
            let low = &self.low_interface[j];
            self.anti_interface[j]
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, a)| {
                    let centred = 0.5 * (ratio[0][i] * q[j][i] + ratio[1][i] * q[j + 1][i]);
                    *a = t * w[i] * centred - low[i];
                });
        }

        // 3. Bounds: q at n−1, n and in the low-order solution, over the cell,
        // its neighbours and the adjacent layers.
        let local = |k: usize, i: usize| -> (f64, f64) {
            let values = [
                input.q_old[k][i] / input.ps_old[i],
                input.q_mid[k][i],
                self.q_low[k][i],
            ];
            values
                .iter()
                .fold((f64::MAX, f64::MIN), |(lo, hi), &v| (lo.min(v), hi.max(v)))
        };
        for k in 0..n {
            let [q_min, q_max] = [&mut self.q_min[k], &mut self.q_max[k]];
            q_min
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .zip(q_max.par_iter_mut())
                .enumerate()
                .for_each(|(i, (lo, hi))| {
                    let (mut a, mut b) = local(k, i);
                    for l in [k.checked_sub(1), (k + 1 < n).then_some(k + 1)]
                        .into_iter()
                        .flatten()
                    {
                        let (c, d) = local(l, i);
                        (a, b) = (a.min(c), b.max(d));
                    }
                    for &c in mesh.cell_neighbors(i as u32) {
                        let (c, d) = local(k, c as usize);
                        (a, b) = (a.min(c), b.max(d));
                    }
                    *lo = a.max(0.0);
                    *hi = b;
                });
        }

        // Zalesak's ratios: the share of the incoming (outgoing)
        // antidiffusive water each cell can take (give) within its bounds.
        for k in 0..n {
            let (anti, anti_v) = (&self.anti[k], &self.anti_interface);
            let (content, q_min, q_max) = (&self.content[k], &self.q_min[k], &self.q_max[k]);
            let [r_plus, r_minus] = [&mut self.r_plus[k], &mut self.r_minus[k]];
            r_plus
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .zip(r_minus.par_iter_mut())
                .enumerate()
                .for_each(|(i, (rp, rm))| {
                    let cid = i as u32;
                    let scale = 1.0 / (g.cell_areas[i] * r);
                    let (mut incoming, mut outgoing) = (0.0, 0.0);
                    for (&e, &s) in mesh.cell_edges(cid).iter().zip(g.cell_edge_signs(cid)) {
                        let x = s * anti[e as usize] * g.edge_lengths[e as usize] * scale;
                        if x > 0.0 {
                            outgoing += x;
                        } else {
                            incoming -= x;
                        }
                    }
                    // Downward positive: out through the bottom, in through
                    // the top.
                    let mut vertical = |v: f64| {
                        if v > 0.0 {
                            outgoing += v;
                        } else {
                            incoming -= v;
                        }
                    };
                    if k + 1 < n {
                        vertical(anti_v[k][i]);
                    }
                    if k > 0 {
                        vertical(-anti_v[k - 1][i]);
                    }
                    let mass = ds[k] * input.ps_new[i];
                    let room = mass * q_max[i] - content[i];
                    let spare = content[i] - mass * q_min[i];
                    *rp = if incoming > 0.0 {
                        (room / incoming).clamp(0.0, 1.0)
                    } else {
                        0.0
                    };
                    *rm = if outgoing > 0.0 {
                        (spare / outgoing).clamp(0.0, 1.0)
                    } else {
                        0.0
                    };
                });
        }

        // Limit and apply.
        for (k, out) in out.into_iter().enumerate() {
            let (anti, rp, rm) = (&self.anti[k], &self.r_plus[k], &self.r_minus[k]);
            self.edge_amount
                .par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(e, x)| {
                    let [a, b] = edges[e].cells;
                    let (a, b) = (a as usize, b as usize);
                    let c = if anti[e] >= 0.0 {
                        rm[a].min(rp[b])
                    } else {
                        rp[a].min(rm[b])
                    };
                    *x = c * anti[e];
                });
            ops.divergence(&self.edge_amount, &mut self.divergence);
            let (content, div) = (&self.content[k], &self.divergence);
            let (anti_v, r_plus, r_minus) = (&self.anti_interface, &self.r_plus, &self.r_minus);
            // The limited flux across interface j (downward positive).
            let limited = |j: usize, i: usize| -> f64 {
                let a = anti_v[j][i];
                let c = if a >= 0.0 {
                    r_minus[j][i].min(r_plus[j + 1][i])
                } else {
                    r_plus[j][i].min(r_minus[j + 1][i])
                };
                c * a
            };
            out.par_iter_mut()
                .with_min_len(MIN_CHUNK)
                .enumerate()
                .for_each(|(i, q)| {
                    let mut vertical = 0.0;
                    if k > 0 {
                        vertical += limited(k - 1, i);
                    }
                    if k + 1 < n {
                        vertical -= limited(k, i);
                    }
                    *q = (content[i] - div[i] + vertical) / ds[k];
                });
        }
        TransportReport { substeps }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::levels::Levels;
    use fms_grid::{LloydOptions, MeshGeometry, SphereMesh};

    /// A random step on a small mesh: edge mass fluxes per layer, the
    /// interface fluxes and end surface pressure they imply.
    struct Setup {
        ps_old: Vec<f64>,
        ps_new: Vec<f64>,
        flux: Vec<Vec<f64>>,
        interface: Vec<Vec<f64>>,
    }

    fn setup(ops: &Operators, levels: &Levels) -> Setup {
        let mesh = ops.mesh();
        let (nc, ne) = (mesh.n_cells(), mesh.n_edges());
        let n = levels.count();
        let ds = levels.dsigma();
        let mut x = 0x9E37_79B9_7F4A_7C15u64;
        let mut rand = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            (x >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        // Winds up to 100 m/s, independent per layer: strong enough to need
        // upwind substeps.
        let flux: Vec<Vec<f64>> = (0..n)
            .map(|k| (0..ne).map(|_| ds[k] * 1e5 * 200.0 * rand()).collect())
            .collect();
        let mut d = vec![vec![0.0; nc]; n];
        for k in 0..n {
            ops.divergence(&flux[k], &mut d[k]);
        }
        let interface = levels
            .interior()
            .iter()
            .enumerate()
            .map(|(j, &s)| {
                (0..nc)
                    .map(|i| {
                        let above: f64 = (0..=j).map(|k| d[k][i]).sum();
                        let below: f64 = (j + 1..n).map(|k| d[k][i]).sum();
                        s * below - (1.0 - s) * above
                    })
                    .collect()
            })
            .collect();
        let ps_old: Vec<f64> = (0..nc).map(|_| 1e5 + 1e3 * rand()).collect();
        let t = 4000.0;
        let ps_new = (0..nc)
            .map(|i| ps_old[i] - t * (0..n).map(|k| d[k][i]).sum::<f64>())
            .collect();
        Setup {
            ps_old,
            ps_new,
            flux,
            interface,
        }
    }

    #[test]
    fn uniform_humidity_stays_uniform_and_patches_stay_positive() {
        let (mesh, _) = SphereMesh::centroidal(8, LloydOptions::default());
        let geometry = MeshGeometry::new(&mesh);
        let ops = Operators::new(&mesh, &geometry, 6.371e6);
        let nc = mesh.n_cells();
        for levels in [
            Levels::two(0.5),
            Levels::new(&[0.4, 0.85]),
            Levels::equal_mass(5),
        ] {
            let n = levels.count();
            let setup = setup(&ops, &levels);
            let ones = vec![1.0; nc];
            let zeros = vec![0.0; nc];
            let mut transport = Transport::new(n, nc, mesh.n_edges());
            let mut step = |q: &[Vec<f64>]| {
                let q_old: Vec<Vec<f64>> = q
                    .iter()
                    .map(|q| q.iter().zip(&setup.ps_old).map(|(q, p)| q * p).collect())
                    .collect();
                let mut out = vec![vec![0.0; nc]; n];
                let report = transport.step(
                    &ops,
                    &TransportStep {
                        dsigma: levels.dsigma(),
                        interval: 4000.0,
                        ps_old: &setup.ps_old,
                        ps_new: &setup.ps_new,
                        flux: setup.flux.iter().map(Vec::as_slice).collect(),
                        interface: setup.interface.iter().map(Vec::as_slice).collect(),
                        q_old: q_old.iter().map(Vec::as_slice).collect(),
                        q_mid: q.iter().map(Vec::as_slice).collect(),
                        source: vec![zeros.as_slice(); n],
                        interface_ratio: vec![[ones.as_slice(), ones.as_slice()]; n - 1],
                    },
                    out.iter_mut().map(Vec::as_mut_slice).collect(),
                );
                assert!(report.substeps > 1, "the test should need substeps");
                (q_old, out)
            };
            // Mass and humidity move together: uniform q stays uniform.
            let (_, out) = step(&vec![vec![0.01; nc]; n]);
            for k in 0..n {
                for i in 0..nc {
                    assert!((out[k][i] / setup.ps_new[i] / 0.01 - 1.0).abs() < 1e-12);
                }
            }
            // Patches stay within their bounds, and water is conserved.
            let patch: Vec<f64> = (0..nc)
                .map(|i| if i % 7 < 2 { 0.02 } else { 0.0 })
                .collect();
            let (q_old, out) = step(&vec![patch; n]);
            let total = |q: &[Vec<f64>]| -> f64 {
                q.iter()
                    .zip(levels.dsigma())
                    .map(|(q, ds)| {
                        ds * q
                            .iter()
                            .zip(&geometry.cell_areas)
                            .map(|(q, a)| q * a)
                            .sum::<f64>()
                    })
                    .sum()
            };
            assert!(
                (total(&out) / total(&q_old) - 1.0).abs() < 1e-13,
                "{n} layers"
            );
            for k in 0..n {
                for i in 0..nc {
                    let q = out[k][i] / setup.ps_new[i];
                    // To round-off (the stepper fills any trace below zero).
                    assert!(
                        (-1e-15..=0.02 * (1.0 + 1e-12)).contains(&q),
                        "{n} layers: q {q}"
                    );
                }
            }
        }
    }
}
