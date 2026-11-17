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
//! - `Σₑ 2Aₑ uₑ tangential_gradient(ψ)ₑ = Σ_v A_v ψ_v curl(u)_v`, and
//!   `divergence(tangential_gradient(ψ)) = 0`
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
//! | vector Laplacian | not convergent | fixed: ~1% divergent, ~12% rotational part |
//!
//! The fixed local errors at the pentagons are tracked in [OBSERVATIONS.md] O-1.
//!
//! [OBSERVATIONS.md]: https://github.com/RazerBlade6/fms/blob/main/docs/climate/OBSERVATIONS.md

use rayon::prelude::*;

use crate::geometry::MeshGeometry;
use crate::mesh::{CellId, SphereMesh, VertexId};

/// Smallest number of elements handed to one thread. Each operator does
/// little work per element, so splitting finer costs more in scheduling than
/// it saves (measured: 7× slower on a 2.5k-cell mesh with 15 threads).
pub const MIN_CHUNK: usize = 4096;

/// Smallest number of columns handed to one thread where each is expensive
/// and their costs differ (the column physics).
pub const MIN_COLUMNS: usize = 32;

/// Whether a loop over `len` elements in chunks of at least [`MIN_CHUNK`] is
/// split between threads: rayon splits a range only while both halves keep
/// `MIN_CHUNK` elements, so below `2 × MIN_CHUNK` it is one chunk.
#[inline]
pub fn splits(len: usize) -> bool {
    len >= 2 * MIN_CHUNK
}

/// `f(i, &mut out[i])` for every element: in chunks of at least
/// [`MIN_CHUNK`] on the pool's threads when that [`splits`], otherwise a
/// plain loop on the calling thread, which pushes no job and wakes no
/// worker. Each element is computed the same way either way.
#[inline]
pub fn for_each_mut<T: Send>(out: &mut [T], f: impl Fn(usize, &mut T) + Sync + Send) {
    if splits(out.len()) {
        out.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .enumerate()
            .for_each(|(i, o)| f(i, o));
    } else {
        for (i, o) in out.iter_mut().enumerate() {
            f(i, o);
        }
    }
}

/// [`for_each_mut`] over two slices of the same length at once.
///
/// # Panics
/// If the lengths differ.
#[inline]
pub fn for_each_mut2<A: Send, B: Send>(
    a: &mut [A],
    b: &mut [B],
    f: impl Fn(usize, &mut A, &mut B) + Sync + Send,
) {
    assert_eq!(a.len(), b.len(), "the slices must have the same length");
    if splits(a.len()) {
        a.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(b.par_iter_mut())
            .enumerate()
            .for_each(|(i, (a, b))| f(i, a, b));
    } else {
        for (i, (a, b)) in a.iter_mut().zip(b.iter_mut()).enumerate() {
            f(i, a, b);
        }
    }
}

/// [`for_each_mut`] over three slices of the same length at once.
///
/// # Panics
/// If the lengths differ.
#[inline]
pub fn for_each_mut3<A: Send, B: Send, C: Send>(
    a: &mut [A],
    b: &mut [B],
    c: &mut [C],
    f: impl Fn(usize, &mut A, &mut B, &mut C) + Sync + Send,
) {
    assert!(
        a.len() == b.len() && a.len() == c.len(),
        "the slices must have the same length"
    );
    if splits(a.len()) {
        a.par_iter_mut()
            .with_min_len(MIN_CHUNK)
            .zip(b.par_iter_mut())
            .zip(c.par_iter_mut())
            .enumerate()
            .for_each(|(i, ((a, b), c))| f(i, a, b, c));
    } else {
        for (i, ((a, b), c)) in a.iter_mut().zip(b.iter_mut()).zip(c.iter_mut()).enumerate() {
            f(i, a, b, c);
        }
    }
}

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
        for_each_mut(out, |e, out| *out = self.gradient_at(e, |c| h[c]));
    }

    /// [`gradient`](Self::gradient) at edge `e` of the cell scalar `h(c)`.
    /// The operators compute each element with these `…_at` functions, so a
    /// loop that fuses them with other work gets the same values bit for bit.
    #[inline]
    pub fn gradient_at(&self, e: usize, h: impl Fn(usize) -> f64) -> f64 {
        let [a, b] = self.mesh.edges()[e].cells;
        (h(b as usize) - h(a as usize)) / (self.geometry.dual_edge_lengths[e] * self.radius)
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
        for_each_mut(out, |c, out| *out = self.divergence_at(c, |e| flux[e]));
    }

    /// [`divergence`](Self::divergence) at cell `c` of the edge flux `flux(e)`.
    #[inline]
    pub fn divergence_at(&self, c: usize, flux: impl Fn(usize) -> f64) -> f64 {
        let g = self.geometry;
        let id = c as CellId;
        let sum: f64 = self
            .mesh
            .cell_edges(id)
            .iter()
            .zip(g.cell_edge_signs(id))
            .map(|(&e, &s)| s * g.edge_lengths[e as usize] * flux(e as usize))
            .sum();
        sum / (g.cell_areas[c] * self.radius)
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
        for_each_mut(out, |v, out| *out = self.curl_at(v, |e| u[e]));
    }

    /// [`curl`](Self::curl) at vertex `v` of the edge velocity `u(e)`.
    #[inline]
    pub fn curl_at(&self, v: usize, u: impl Fn(usize) -> f64) -> f64 {
        let g = self.geometry;
        let id = v as VertexId;
        let sum: f64 = self
            .mesh
            .vertex_edges(id)
            .iter()
            .zip(g.vertex_edge_signs(id))
            .map(|(&e, s)| s * g.dual_edge_lengths[e as usize] * u(e as usize))
            .sum();
        sum / (g.vertex_areas[v] * self.radius)
    }

    /// Derivative of a vertex scalar along each edge tangent `tₑ`:
    /// `(ψ[vertices[1]] − ψ[vertices[0]]) / lₑ`.
    ///
    /// The adjoint of [`curl`](Self::curl):
    /// `Σₑ 2Aₑ uₑ tangential_gradient(ψ)ₑ = Σ_v A_v ψ_v curl(u)_v`. Its divergence
    /// is zero (a closed loop of differences around each cell), so
    /// `−tangential_gradient(ψ)` is the flow of stream function `ψ`, and the
    /// vector Laplacian is `∇²u = gradient(δ) − tangential_gradient(ζ)`.
    pub fn tangential_gradient(&self, psi: &[f64], out: &mut [f64]) {
        self.check(
            psi.len(),
            self.mesh.n_vertices(),
            out.len(),
            self.mesh.n_edges(),
        );
        for_each_mut(out, |e, out| {
            *out = self.tangential_gradient_at(e, |v| psi[v])
        });
    }

    /// [`tangential_gradient`](Self::tangential_gradient) at edge `e` of the
    /// vertex scalar `psi(v)`.
    #[inline]
    pub fn tangential_gradient_at(&self, e: usize, psi: impl Fn(usize) -> f64) -> f64 {
        let [a, b] = self.mesh.edges()[e].vertices;
        (psi(b as usize) - psi(a as usize)) / (self.geometry.edge_lengths[e] * self.radius)
    }

    /// Tangential component (along `tₑ`) of the velocity at each edge, from
    /// the normal components of the edges around it (TRiSK):
    /// `vₑ = (1/dₑ) Σ_{e'} w_{e,e'} l_{e'} u_{e'}`.
    pub fn tangential(&self, u: &[f64], out: &mut [f64]) {
        self.check(u.len(), self.mesh.n_edges(), out.len(), self.mesh.n_edges());
        let g = self.geometry;
        for_each_mut(out, |e, out| {
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
        for_each_mut(out, |e, out| *out = self.pv_flux_at(e, q, flux));
    }

    /// [`pv_flux`](Self::pv_flux) at edge `e`.
    #[inline]
    pub fn pv_flux_at(&self, e: usize, q: &[f64], flux: &[f64]) -> f64 {
        let g = self.geometry;
        let (stencil, weights) = g.trisk_stencil(e as u32);
        let sum: f64 = stencil
            .iter()
            .zip(weights)
            .map(|(&f, &w)| {
                let f = f as usize;
                w * g.edge_lengths[f] * flux[f] * 0.5 * (q[e] + q[f])
            })
            .sum();
        sum / g.dual_edge_lengths[e]
    }

    /// Kinetic energy per unit mass at each cell (Ringler et al. 2010):
    /// `Kᵢ = (1/Aᵢ) Σₑ (Aₑ/2) uₑ²`, with `Aₑ = lₑdₑ/2`.
    pub fn kinetic_energy(&self, u: &[f64], out: &mut [f64]) {
        self.check(u.len(), self.mesh.n_edges(), out.len(), self.mesh.n_cells());
        for_each_mut(out, |c, out| *out = self.kinetic_energy_at(c, u));
    }

    /// [`kinetic_energy`](Self::kinetic_energy) at cell `c`.
    #[inline]
    pub fn kinetic_energy_at(&self, c: usize, u: &[f64]) -> f64 {
        let g = self.geometry;
        let sum: f64 = self
            .mesh
            .cell_edges(c as CellId)
            .iter()
            .map(|&e| {
                let e = e as usize;
                0.25 * g.edge_lengths[e] * g.dual_edge_lengths[e] * u[e] * u[e]
            })
            .sum();
        sum / g.cell_areas[c]
    }

    /// Cell scalar at edges: the mean of the two cells either side.
    pub fn cell_to_edge(&self, h: &[f64], out: &mut [f64]) {
        self.check(h.len(), self.mesh.n_cells(), out.len(), self.mesh.n_edges());
        for_each_mut(out, |e, out| *out = self.cell_to_edge_at(e, |c| h[c]));
    }

    /// [`cell_to_edge`](Self::cell_to_edge) at edge `e` of the cell scalar `h(c)`.
    #[inline]
    pub fn cell_to_edge_at(&self, e: usize, h: impl Fn(usize) -> f64) -> f64 {
        let [a, b] = self.mesh.edges()[e].cells;
        0.5 * (h(a as usize) + h(b as usize))
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
        for_each_mut(out, |v, out| *out = self.cell_to_vertex_at(v, |c| h[c]));
    }

    /// [`cell_to_vertex`](Self::cell_to_vertex) at vertex `v` of the cell
    /// scalar `h(c)`.
    #[inline]
    pub fn cell_to_vertex_at(&self, v: usize, h: impl Fn(usize) -> f64) -> f64 {
        let g = self.geometry;
        let kites = g.vertex_kites(v as VertexId);
        let sum: f64 = self.mesh.triangles()[v]
            .iter()
            .zip(kites)
            .map(|(&c, k)| k * h(c as usize))
            .sum();
        sum / g.vertex_areas[v]
    }

    /// Vertex scalar at edges: the mean of the edge's two end vertices.
    pub fn vertex_to_edge(&self, q: &[f64], out: &mut [f64]) {
        self.check(
            q.len(),
            self.mesh.n_vertices(),
            out.len(),
            self.mesh.n_edges(),
        );
        for_each_mut(out, |e, out| *out = self.vertex_to_edge_at(e, |v| q[v]));
    }

    /// [`vertex_to_edge`](Self::vertex_to_edge) at edge `e` of the vertex
    /// scalar `q(v)`.
    #[inline]
    pub fn vertex_to_edge_at(&self, e: usize, q: impl Fn(usize) -> f64) -> f64 {
        let [a, b] = self.mesh.edges()[e].vertices;
        0.5 * (q(a as usize) + q(b as usize))
    }

    fn check(&self, input: usize, expected_in: usize, output: usize, expected_out: usize) {
        assert_eq!(input, expected_in, "input has the wrong length");
        assert_eq!(output, expected_out, "output has the wrong length");
    }
}
