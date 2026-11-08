//! Discrete TRiSK operators on the C-grid (Thuburn et al. 2009; Ringler et al.
//! 2010), the scheme used by MPAS.
//!
//! Scalars live at cells, normal velocities `uₑ` (along `nₑ`) at edges, and
//! vorticity at vertices; see [`MeshGeometry`] for the sign conventions. Every
//! operator works in physical units: lengths are scaled by the planet radius
//! given to [`Operators::new`], so a gradient of metres comes out in m m⁻¹ and
//! a curl of m s⁻¹ in s⁻¹.
//!
//! Each operator writes into a caller-supplied slice, allocates nothing, and
//! computes every output element independently in a fixed order, so results
//! don't depend on the thread count.
//!
//! The operators are *mimetic*: their discrete forms keep the identities the
//! continuous ones have, exactly (to round-off), on any mesh:
//!
//! - `curl(gradient(h)) = 0` and `divergence` of a stream-function flow is 0
//! - `Σᵢ Aᵢ divergence(F)ᵢ = 0` (Gauss: no net source, so mass is conserved)
//! - `Σₑ 2Aₑ gradient(h)ₑ Fₑ = −Σᵢ Aᵢ hᵢ divergence(F)ᵢ`, with `Aₑ = lₑdₑ/2`
//!   (gradient and divergence are adjoint, so pressure work exchanges energy
//!   without creating it)
//! - `Σᵢ Aᵢ kinetic_energy(u)ᵢ = Σₑ Aₑ uₑ²`
//! - `Σ_v A_v cell_to_vertex(h)_v = Σᵢ Aᵢ hᵢ`
//! - TRiSK weights are antisymmetric, so the Coriolis term does no work
//!
//! Accuracy on centroidal grids, measured against spherical harmonics:
//!
//! | Operator | Max error | RMS error |
//! |---|---|---|
//! | gradient, cell → edge | 2nd order | 2nd order |
//! | vertex → edge | 1st order | 2nd order |
//! | curl, cell → vertex | 1st order | 1st order |
//! | divergence, Laplacian | fixed, ~1.1% at the pentagons | 1st order |
//! | tangential velocity | fixed, ~0.6% at the pentagons | ~1.5 order |
//! | kinetic energy | fixed, ~1.5% at the pentagons | ~1.3 order |
//!
//! The fixed local errors at the pentagons are tracked in [OBSERVATIONS.md] O-1.
//!
//! [OBSERVATIONS.md]: https://github.com/RazerBlade6/fantasy-map-simulator/blob/main/docs/OBSERVATIONS.md

use rayon::prelude::*;

use crate::geometry::MeshGeometry;
use crate::mesh::{CellId, SphereMesh, VertexId};

/// Smallest number of elements handed to one thread. Each operator does
/// little work per element, so splitting finer costs more in scheduling than
/// it saves (measured: 7× slower on a 2.5k-cell mesh with 15 threads).
pub const MIN_CHUNK: usize = 4096;

/// The TRiSK operators for one mesh at one planet radius.
#[derive(Clone, Copy, Debug)]
pub struct Operators<'a> {
    mesh: &'a SphereMesh,
    geometry: &'a MeshGeometry,
    radius: f64,
}

impl<'a> Operators<'a> {
    /// Operators on `mesh` scaled to a sphere of `radius` (m).
    ///
    /// # Panics
    /// If `radius` is not positive.
    pub fn new(mesh: &'a SphereMesh, geometry: &'a MeshGeometry, radius: f64) -> Self {
        assert!(radius > 0.0, "radius must be positive");
        Self {
            mesh,
            geometry,
            radius,
        }
    }

    /// The mesh.
    pub fn mesh(&self) -> &'a SphereMesh {
        self.mesh
    }

    /// The mesh geometry.
    pub fn geometry(&self) -> &'a MeshGeometry {
        self.geometry
    }

    /// The planet radius (m).
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Area `Aₑ = lₑdₑ/2` of the rhombus around each edge (m²): the weight of
    /// edge quantities in global sums. Total kinetic energy per unit mass is
    /// `Σₑ Aₑ uₑ²` (the mean of `uₙ²` over edge directions is `½|u|²`).
    pub fn edge_areas(&self) -> Vec<f64> {
        let r2 = self.radius * self.radius;
        let g = self.geometry;
        g.edge_lengths
            .iter()
            .zip(&g.dual_edge_lengths)
            .map(|(l, d)| 0.5 * l * d * r2)
            .collect()
    }

    /// Gradient of a cell scalar, normal to each edge:
    /// `(h[cells[1]] − h[cells[0]]) / dₑ`.
    pub fn gradient(&self, h: &[f64], out: &mut [f64]) {
        self.check(h.len(), self.mesh.n_cells(), out.len(), self.mesh.n_edges());
        let edges = self.mesh.edges();
        let d = &self.geometry.dual_edge_lengths;
        let r = self.radius;
        out.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(e, out)| {
                let [a, b] = edges[e].cells;
                *out = (h[b as usize] - h[a as usize]) / (d[e] * r);
            });
    }

    /// Divergence of an edge flux (normal component per edge), per cell:
    /// `(1/Aᵢ) Σₑ n_{e,i} lₑ Fₑ`.
    pub fn divergence(&self, flux: &[f64], out: &mut [f64]) {
        self.check(
            flux.len(),
            self.mesh.n_edges(),
            out.len(),
            self.mesh.n_cells(),
        );
        let g = self.geometry;
        let r = self.radius;
        out.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(c, out)| {
                let c = c as CellId;
                let sum: f64 = self
                    .mesh
                    .cell_edges(c)
                    .iter()
                    .zip(g.cell_edge_signs(c))
                    .map(|(&e, &s)| s * g.edge_lengths[e as usize] * flux[e as usize])
                    .sum();
                *out = sum / (g.cell_areas[c as usize] * r);
            });
    }

    /// Relative vorticity of an edge velocity, per vertex:
    /// `(1/A_v) Σₑ t_{e,v} dₑ uₑ` (circulation counter-clockwise, seen from
    /// outside the sphere).
    pub fn curl(&self, u: &[f64], out: &mut [f64]) {
        self.check(
            u.len(),
            self.mesh.n_edges(),
            out.len(),
            self.mesh.n_vertices(),
        );
        let g = self.geometry;
        let r = self.radius;
        out.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(v, out)| {
                let v = v as VertexId;
                let sum: f64 = self
                    .mesh
                    .vertex_edges(v)
                    .iter()
                    .zip(g.vertex_edge_signs(v))
                    .map(|(&e, s)| s * g.dual_edge_lengths[e as usize] * u[e as usize])
                    .sum();
                *out = sum / (g.vertex_areas[v as usize] * r);
            });
    }

    /// Tangential component (along `tₑ`) of the velocity at each edge, from
    /// the normal components of the edges around it (TRiSK):
    /// `vₑ = (1/dₑ) Σ_{e'} w_{e,e'} l_{e'} u_{e'}`.
    pub fn tangential(&self, u: &[f64], out: &mut [f64]) {
        self.check(u.len(), self.mesh.n_edges(), out.len(), self.mesh.n_edges());
        let g = self.geometry;
        out.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(e, out)| {
                let (stencil, weights) = g.trisk_stencil(e as u32);
                let sum: f64 = stencil
                    .iter()
                    .zip(weights)
                    .map(|(&f, &w)| w * g.edge_lengths[f as usize] * u[f as usize])
                    .sum();
                *out = sum / g.dual_edge_lengths[e];
            });
    }

    /// The potential-vorticity flux `q F⊥` at each edge, in the
    /// energy-conserving TRiSK form (Ringler et al. 2010, eq. 49):
    /// `(1/dₑ) Σ_{e'} w_{e,e'} l_{e'} F_{e'} (qₑ + q_{e'})/2`.
    ///
    /// `flux` is the thickness flux `Fₑ = hₑuₑ` and `q` the potential vorticity
    /// at edges. In the vector-invariant momentum equation this is the
    /// `+q·(F along tₑ)` term, i.e. minus the normal component of `(ζ + f) k × u`.
    /// With constant `q` it equals `q · tangential(flux)`. Because the weights
    /// are antisymmetric and the `q` average symmetric, it does no work:
    /// `Σₑ 2Aₑ Fₑ · pv_flux(q, F)ₑ = 0`.
    pub fn pv_flux(&self, q: &[f64], flux: &[f64], out: &mut [f64]) {
        let n = self.mesh.n_edges();
        self.check(flux.len(), n, out.len(), n);
        assert_eq!(q.len(), n, "q has the wrong length");
        let g = self.geometry;
        out.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(e, out)| {
                let (stencil, weights) = g.trisk_stencil(e as u32);
                let sum: f64 = stencil
                    .iter()
                    .zip(weights)
                    .map(|(&f, &w)| {
                        let f = f as usize;
                        w * g.edge_lengths[f] * flux[f] * 0.5 * (q[e] + q[f])
                    })
                    .sum();
                *out = sum / g.dual_edge_lengths[e];
            });
    }

    /// Kinetic energy per unit mass at each cell (Ringler et al. 2010):
    /// `Kᵢ = (1/Aᵢ) Σₑ (Aₑ/2) uₑ²`, with `Aₑ = lₑdₑ/2`.
    pub fn kinetic_energy(&self, u: &[f64], out: &mut [f64]) {
        self.check(u.len(), self.mesh.n_edges(), out.len(), self.mesh.n_cells());
        let g = self.geometry;
        out.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(c, out)| {
                let c = c as CellId;
                let sum: f64 = self
                    .mesh
                    .cell_edges(c)
                    .iter()
                    .map(|&e| {
                        let e = e as usize;
                        0.25 * g.edge_lengths[e] * g.dual_edge_lengths[e] * u[e] * u[e]
                    })
                    .sum();
                *out = sum / g.cell_areas[c as usize];
            });
    }

    /// Cell scalar at edges: the mean of the two cells either side.
    pub fn cell_to_edge(&self, h: &[f64], out: &mut [f64]) {
        self.check(h.len(), self.mesh.n_cells(), out.len(), self.mesh.n_edges());
        let edges = self.mesh.edges();
        out.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(e, out)| {
                let [a, b] = edges[e].cells;
                *out = 0.5 * (h[a as usize] + h[b as usize]);
            });
    }

    /// Cell scalar at vertices: the mean of the three cells around each,
    /// weighted by the kite area each shares with the dual triangle.
    pub fn cell_to_vertex(&self, h: &[f64], out: &mut [f64]) {
        self.check(
            h.len(),
            self.mesh.n_cells(),
            out.len(),
            self.mesh.n_vertices(),
        );
        let g = self.geometry;
        let triangles = self.mesh.triangles();
        out.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(v, out)| {
                let kites = g.vertex_kites(v as VertexId);
                let sum: f64 = triangles[v]
                    .iter()
                    .zip(kites)
                    .map(|(&c, k)| k * h[c as usize])
                    .sum();
                *out = sum / g.vertex_areas[v];
            });
    }

    /// Vertex scalar at edges: the mean of the edge's two end vertices.
    pub fn vertex_to_edge(&self, q: &[f64], out: &mut [f64]) {
        self.check(
            q.len(),
            self.mesh.n_vertices(),
            out.len(),
            self.mesh.n_edges(),
        );
        let edges = self.mesh.edges();
        out.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(e, out)| {
                let [a, b] = edges[e].vertices;
                *out = 0.5 * (q[a as usize] + q[b as usize]);
            });
    }

    fn check(&self, input: usize, expected_in: usize, output: usize, expected_out: usize) {
        assert_eq!(input, expected_in, "input has the wrong length");
        assert_eq!(output, expected_out, "output has the wrong length");
    }
}
