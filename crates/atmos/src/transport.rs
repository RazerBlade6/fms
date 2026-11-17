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
use fms_grid::operators::{MIN_CHUNK, for_each_mut, for_each_mut2, for_each_mut3};
use rayon::prelude::*;

use crate::levels::MAX_LAYERS;

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
    pub flux: &'s [Vec<f64>],
    /// Mass flux `W̄ⱼ` across each interface (`p_s σ̇`, Pa s⁻¹), positive
    /// downward.
    pub interface: &'s [Vec<f64>],
    /// `Qₖ` at the start (Pa).
    pub q_old: &'s [Vec<f64>],
    /// Specific humidity `qₖ` at the middle level (kg kg⁻¹).
    pub q_mid: &'s [Vec<f64>],
    /// Sources of `Qₖ` (Pa s⁻¹): physics, constant over the step; `None`
    /// for none (zero).
    pub source: Option<&'s [Vec<f64>]>,
    /// Per interface, humidity there relative to each adjacent layer's mean:
    /// the upper layer at its bottom, the lower at its top.
    pub interface_ratio: &'s [[Vec<f64>; 2]],
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
    /// Per cell and layer, the range of q at n−1, n and in the low-order
    /// solution.
    local_min: Vec<Vec<f64>>,
    local_max: Vec<Vec<f64>>,
    r_plus: Vec<Vec<f64>>,
    r_minus: Vec<Vec<f64>>,
    edge_amount: Vec<Vec<f64>>,
    divergence: Vec<Vec<f64>>,
    q: Vec<Vec<f64>>,
    /// The sources when there are none.
    zeros: Vec<f64>,
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
            local_min: cells(layers),
            local_max: cells(layers),
            r_plus: cells(layers),
            r_minus: cells(layers),
            edge_amount: edges(),
            divergence: cells(layers),
            q: cells(layers),
            zeros: vec![0.0; n_cells],
        }
    }

    /// Transport humidity over one step, writing `Qₖ` at the end into `out`.
    pub fn step(
        &mut self,
        ops: &Operators,
        input: &TransportStep,
        out: &mut [Vec<f64>],
    ) -> TransportReport {
        let mut nothing = [(); MAX_LAYERS];
        let n = input.dsigma.len();
        self.step_then(ops, input, out, &mut nothing[..n], |_, _, _| {})
    }

    /// [`step`](Self::step), then `then(k, Qₖ, &mut payload[k])` for each
    /// layer as soon as its `Qₖ` is written, in the same task: per-layer
    /// work that would otherwise need a parallel region of its own.
    pub fn step_then<T: Send>(
        &mut self,
        ops: &Operators,
        input: &TransportStep,
        out: &mut [Vec<f64>],
        payload: &mut [T],
        then: impl Fn(usize, &mut [f64], &mut T) + Sync + Send,
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
        // Per layer, then over the layers: a maximum is exact in any order.
        let mut layer_max = [0.0; MAX_LAYERS];
        layer_max[..n]
            .par_iter_mut()
            .enumerate()
            .for_each(|(k, m)| {
                *m = (0..mesh.n_cells())
                    .into_par_iter()
                    .with_min_len(MIN_CHUNK)
                    .map(|c| outflow_rate(c, k))
                    .reduce(|| 0.0, f64::max);
            });
        let max_rate = layer_max[..n].iter().fold(0.0, |a: f64, &b| a.max(b));
        let substeps = ((t * max_rate / MAX_OUTFLOW).ceil() as u32).max(1);
        let h = t / f64::from(substeps);

        // 1. Low order: upwind substeps from n−1, with the sources. Each
        // substep's horizontal part per layer (the first also starting the
        // layer's water from n−1), then the exchange across the interfaces.
        for x in &mut self.low_interface {
            x.fill(0.0);
        }
        let zeros = &self.zeros;
        for m in 0..substeps {
            let s = f64::from(m) / f64::from(substeps);
            self.content
                .par_iter_mut()
                .zip(self.q.par_iter_mut())
                .zip(self.low_flux.par_iter_mut())
                .zip(self.edge_amount.par_iter_mut())
                .zip(self.divergence.par_iter_mut())
                .enumerate()
                .for_each(
                    |(k, ((((content, q), low_flux), edge_amount), divergence))| {
                        if m == 0 {
                            let q_old = &input.q_old[k];
                            for_each_mut(content, |i, c| *c = ds[k] * q_old[i]);
                            low_flux.fill(0.0);
                        }
                        {
                            let content = &*content;
                            for_each_mut(q, |i, q| {
                                let ps = input.ps_old[i] + s * (input.ps_new[i] - input.ps_old[i]);
                                *q = content[i] / (ds[k] * ps);
                            });
                        }
                        let (q, flux) = (&*q, &input.flux[k]);
                        for_each_mut2(edge_amount, low_flux, |e, x, total| {
                            let [a, b] = edges[e].cells;
                            let upwind = if flux[e] >= 0.0 {
                                q[a as usize]
                            } else {
                                q[b as usize]
                            };
                            *x = h * flux[e] * upwind;
                            *total += *x;
                        });
                        ops.divergence(edge_amount, divergence);
                        let source = input.source.map_or(zeros, |s| &s[k]);
                        let div = &*divergence;
                        for_each_mut(content, |i, c| *c += -div[i] + h * ds[k] * source[i]);
                    },
                );
            // Across each interface (downward positive: from layer j to j+1):
            // cheap, and it couples the layers, so on this thread.
            for j in 0..n - 1 {
                let (upper, lower) = self.content.split_at_mut(j + 1);
                let (c0, c1) = (&mut upper[j], &mut lower[0]);
                let (q, w, ratio) = (&self.q, &w[j], &ratio[j]);
                for_each_mut3(c0, c1, &mut self.low_interface[j], |i, c0, c1, total| {
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

        // Per layer: the low-order solution, the range of q at n−1, n and in
        // it, and the antidiffusive fluxes (centred at n minus upwind) along
        // the edges and across the interface below.
        {
            let (content_all, low_flux_all, low_interface) =
                (&self.content, &self.low_flux, &self.low_interface);
            let mut anti_below: [Option<&mut Vec<f64>>; MAX_LAYERS] = {
                let mut a = self.anti_interface.iter_mut();
                std::array::from_fn(|_| a.next())
            };
            self.q_low
                .par_iter_mut()
                .zip(self.local_min.par_iter_mut())
                .zip(self.local_max.par_iter_mut())
                .zip(self.anti.par_iter_mut())
                .zip(anti_below[..n].par_iter_mut())
                .enumerate()
                .for_each(
                    |(k, ((((q_low, local_min), local_max), anti), anti_below))| {
                        let content = &content_all[k];
                        for_each_mut(q_low, |i, q| *q = content[i] / (ds[k] * input.ps_new[i]));
                        let q_low = &*q_low;
                        for_each_mut2(local_min, local_max, |i, lo, hi| {
                            let values = [
                                input.q_old[k][i] / input.ps_old[i],
                                input.q_mid[k][i],
                                q_low[i],
                            ];
                            (*lo, *hi) = values
                                .iter()
                                .fold((f64::MAX, f64::MIN), |(lo, hi), &v| (lo.min(v), hi.max(v)));
                        });
                        let (q, flux, low) = (&input.q_mid[k], &input.flux[k], &low_flux_all[k]);
                        for_each_mut(anti, |e, a| {
                            let [i, j] = edges[e].cells;
                            let centred = 0.5 * (q[i as usize] + q[j as usize]);
                            *a = t * flux[e] * centred - low[e];
                        });
                        if let Some(anti_below) = anti_below.as_deref_mut() {
                            let j = k;
                            let (q, w, ratio) = (input.q_mid, &w[j], &ratio[j]);
                            let low = &low_interface[j];
                            for_each_mut(anti_below, |i, a| {
                                let centred =
                                    0.5 * (ratio[0][i] * q[j][i] + ratio[1][i] * q[j + 1][i]);
                                *a = t * w[i] * centred - low[i];
                            });
                        }
                    },
                );
        }

        // Per layer: the bounds (q at n−1, n and in the low-order solution,
        // over the cell, its neighbours and the adjacent layers), then
        // Zalesak's ratios: the share of the incoming (outgoing)
        // antidiffusive water each cell can take (give) within them.
        let (local_min, local_max) = (&self.local_min, &self.local_max);
        let local = |k: usize, i: usize| (local_min[k][i], local_max[k][i]);
        let (anti_all, anti_v) = (&self.anti, &self.anti_interface);
        let content_all = &self.content;
        self.q_min
            .par_iter_mut()
            .zip(self.q_max.par_iter_mut())
            .zip(self.r_plus.par_iter_mut())
            .zip(self.r_minus.par_iter_mut())
            .enumerate()
            .for_each(|(k, (((q_min, q_max), r_plus), r_minus))| {
                for_each_mut2(q_min, q_max, |i, lo, hi| {
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
                let anti = &anti_all[k];
                let (content, q_min, q_max) = (&content_all[k], &*q_min, &*q_max);
                for_each_mut2(r_plus, r_minus, |i, rp, rm| {
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
            });

        // Limit and apply, per layer.
        let (anti_all, anti_v) = (&self.anti, &self.anti_interface);
        let (content_all, r_plus, r_minus) = (&self.content, &self.r_plus, &self.r_minus);
        out.par_iter_mut()
            .zip(self.edge_amount.par_iter_mut())
            .zip(self.divergence.par_iter_mut())
            .zip(payload.par_iter_mut())
            .enumerate()
            .for_each(|(k, (((out, edge_amount), divergence), payload))| {
                let (anti, rp, rm) = (&anti_all[k], &r_plus[k], &r_minus[k]);
                for_each_mut(edge_amount, |e, x| {
                    let [a, b] = edges[e].cells;
                    let (a, b) = (a as usize, b as usize);
                    let c = if anti[e] >= 0.0 {
                        rm[a].min(rp[b])
                    } else {
                        rp[a].min(rm[b])
                    };
                    *x = c * anti[e];
                });
                ops.divergence(edge_amount, divergence);
                let (content, div) = (&content_all[k], &*divergence);
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
                for_each_mut(out, |i, q| {
                    let mut vertical = 0.0;
                    if k > 0 {
                        vertical += limited(k - 1, i);
                    }
                    if k + 1 < n {
                        vertical -= limited(k, i);
                    }
                    *q = (content[i] - div[i] + vertical) / ds[k];
                });
                then(k, out, payload);
            });
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
            let ratios = vec![[vec![1.0; nc], vec![1.0; nc]]; n - 1];
            let zeros = vec![vec![0.0; nc]; n];
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
                        flux: &setup.flux,
                        interface: &setup.interface,
                        q_old: &q_old,
                        q_mid: q,
                        source: Some(&zeros),
                        interface_ratio: &ratios,
                    },
                    &mut out,
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
