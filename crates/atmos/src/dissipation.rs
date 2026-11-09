//! Dissipative heating: kinetic energy that friction, mixing or diffusion
//! removes from the winds, returned as heat, so total energy is conserved
//! ([`crate::LayeredModel::energy_tendency`]) instead of leaking out of the
//! planet ([OBSERVATIONS.md](../../../docs/OBSERVATIONS.md) O-21).
//!
//! A wind tendency `δuₑ` of layer `k` changes the kinetic energy at rate
//! `Aₑ δσₖ (p_a + p_b) uₑ δuₑ / g` (the exact derivative of the model's
//! kinetic energy, `Aₑ` the edge's area and `p_a`, `p_b` its cells' surface
//! pressures). Half of each edge's rate goes to each of its cells, as
//! heating of a chosen layer there; so the heating's enthalpy rate cancels
//! the winds' kinetic-energy rate exactly.

use fms_grid::Operators;
use fms_grid::operators::MIN_CHUNK;
use rayon::prelude::*;

use crate::layered::Gas;
use crate::levels::Levels;

/// Turns kinetic-energy tendencies into heating.
#[derive(Clone, Debug)]
pub struct KineticHeating {
    gas: Gas,
    levels: Levels,
    /// Edge areas `Aₑ` (m²).
    edge_areas: Vec<f64>,
    /// `1/(Aᵢ c_p)` per cell (m⁻² per J kg⁻¹ K⁻¹).
    cell_weights: Vec<f64>,
    /// Layer Exner functions per cell (`n` per cell), from [`prepare`](Self::prepare).
    exner: Vec<f64>,
}

impl KineticHeating {
    /// For the mesh of `ops` and the model's gas and layers.
    pub fn new(ops: &Operators, gas: Gas, levels: Levels) -> Self {
        let r2 = ops.radius() * ops.radius();
        let cell_weights = ops
            .geometry()
            .cell_areas
            .iter()
            .map(|a| 1.0 / (a * r2 * gas.heat_capacity))
            .collect();
        Self {
            gas,
            exner: vec![0.0; levels.count() * ops.mesh().n_cells()],
            levels,
            edge_areas: ops.edge_areas(),
            cell_weights,
        }
    }

    /// Layer Exner functions for surface pressure `ps`; call before
    /// [`add`](Self::add) whenever `ps` changes.
    pub fn prepare(&mut self, ps: &[f64]) {
        let (gas, levels) = (self.gas, &self.levels);
        self.exner
            .par_chunks_mut(levels.count())
            .with_min_len(MIN_CHUNK / 128)
            .zip(ps.par_iter())
            .for_each(|(e, &p)| levels.layer_exner_into(gas, p, e));
    }

    /// Rate of change of kinetic energy (W) from the wind tendency `du` of
    /// layer `k` (winds `u`), with gravity `gravity`.
    pub fn rate(
        &self,
        ops: &Operators,
        ps: &[f64],
        k: usize,
        gravity: f64,
        u: &[f64],
        du: &[f64],
    ) -> f64 {
        let edges = ops.mesh().edges();
        let ds = self.levels.dsigma()[k];
        (0..edges.len())
            .into_par_iter()
            .with_min_len(MIN_CHUNK)
            .map(|e| {
                let [a, b] = edges[e].cells;
                self.edge_areas[e] * ds * (ps[a as usize] + ps[b as usize]) * u[e] * du[e] / gravity
            })
            .sum()
    }

    /// Add to `theta_ps` (the `Θ` tendency of layer `into`, Pa K s⁻¹) the
    /// heating that returns the kinetic energy the wind tendency `du` of
    /// layer `from` (winds `u`) removes; negative where `du` adds energy.
    #[allow(clippy::too_many_arguments)]
    pub fn add(
        &self,
        ops: &Operators,
        ps: &[f64],
        from: usize,
        into: usize,
        u: &[f64],
        du: &[f64],
        theta_ps: &mut [f64],
    ) {
        let mesh = ops.mesh();
        let edges = mesh.edges();
        let n = self.levels.count();
        let ds = self.levels.dsigma();
        let ratio = ds[from] / ds[into];
        let areas = &self.edge_areas;
        theta_ps
            .par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, t)| {
                let mut gain = 0.0;
                for &e in mesh.cell_edges(i as u32) {
                    let e = e as usize;
                    let [a, b] = edges[e].cells;
                    gain += 0.5 * areas[e] * (ps[a as usize] + ps[b as usize]) * u[e] * du[e];
                }
                // Θ rate = g q / (c_p δσ Π) with q = −gain δσ_from / (g Aᵢ).
                *t -= gain * ratio * self.cell_weights[i] / self.exner[i * n + into];
            });
    }
}
