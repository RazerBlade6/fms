//! Metric quantities and TRiSK weights for finite-volume operators on a
//! [`SphereMesh`].
//!
//! Conventions follow Thuburn et al. (2009) and Ringler et al. (2010), the scheme
//! used by MPAS:
//!
//! - Scalars (mass, temperature, …) live at **cells**.
//! - Each **edge** `e` carries a normal velocity `uₑ` along `nₑ`, the unit vector
//!   at the edge point pointing from `cells[0]` to `cells[1]`.
//! - `tₑ = k × nₑ` (with `k` the local vertical) is the edge tangent. It points
//!   from `vertices[0]` to `vertices[1]`.
//! - `lₑ` is the length of the Voronoi edge (between its two vertices); `dₑ` is
//!   the distance between the two cell centres.
//! - Vorticity lives at **vertices** (the dual triangles).
//!
//! With these, the basic operators are:
//!
//! - gradient at an edge: `(h[cells[1]] − h[cells[0]]) / dₑ`
//! - divergence at a cell: `(1/Aᵢ) Σₑ n_{e,i} lₑ uₑ`
//! - curl at a vertex: `(1/A_v) Σₑ t_{e,v} dₑ uₑ`
//! - tangential velocity at an edge (TRiSK): `vₑ = (1/dₑ) Σ_{e'} w_{e,e'} l_{e'} u_{e'}`
//!
//! All quantities are on the unit sphere: scale lengths by the radius `r` and
//! areas by `r²`.

use glam::DVec3;
use rayon::prelude::*;

use crate::mesh::{CellId, EdgeId, SphereMesh, VertexId};
use crate::spherical::{arc_length, triangle_area};

/// Precomputed metric terms and operator weights for a mesh.
#[derive(Clone, Debug)]
pub struct MeshGeometry {
    /// Voronoi cell areas `Aᵢ`.
    pub cell_areas: Vec<f64>,
    /// Dual (triangle) areas `A_v`: the sum of the three kites at each vertex.
    pub vertex_areas: Vec<f64>,
    /// Edge points: midpoint of the arc between the two cell centres.
    pub edge_points: Vec<DVec3>,
    /// Unit normals `nₑ` at edge points, from `cells[0]` to `cells[1]`.
    pub edge_normals: Vec<DVec3>,
    /// Unit tangents `tₑ = k × nₑ` at edge points, from `vertices[0]` to `vertices[1]`.
    pub edge_tangents: Vec<DVec3>,
    /// Voronoi edge lengths `lₑ`.
    pub edge_lengths: Vec<f64>,
    /// Distances between cell centres `dₑ`.
    pub dual_edge_lengths: Vec<f64>,
    /// Kite areas, aligned with [`SphereMesh::cell_vertices`]: the area shared by a
    /// cell and the dual triangle of each of its vertices.
    kite_areas: Vec<f64>,
    /// Same kites indexed by vertex, aligned with that triangle's cell order.
    vertex_kites: Vec<[f64; 3]>,
    /// `n_{e,i}` for each cell edge, aligned with [`SphereMesh::cell_edges`]: `+1`
    /// if `nₑ` points out of the cell, `−1` otherwise.
    cell_edge_signs: Vec<f64>,
    /// `t_{e,v}` for each vertex edge, aligned with [`SphereMesh::vertex_edges`]:
    /// `+1` if `nₑ` points counter-clockwise around the vertex.
    vertex_edge_signs: Vec<[f64; 3]>,
    /// CSR offsets into `trisk_edges` / `trisk_weights`; length `n_edges + 1`.
    trisk_offsets: Vec<u32>,
    trisk_edges: Vec<EdgeId>,
    trisk_weights: Vec<f64>,
    /// Offsets into `kite_areas` / `cell_edge_signs` (copied from the mesh rings).
    ring_offsets: Vec<u32>,
}

impl MeshGeometry {
    /// Compute all metric terms and TRiSK weights for `mesh`.
    pub fn new(mesh: &SphereMesh) -> Self {
        let centers = mesh.cell_centers();
        let vertices = mesh.vertices();
        let edges = mesh.edges();

        let mut ring_offsets = Vec::with_capacity(mesh.n_cells() + 1);
        ring_offsets.push(0u32);
        for c in 0..mesh.n_cells() as CellId {
            let next = ring_offsets[c as usize] + mesh.cell_vertices(c).len() as u32;
            ring_offsets.push(next);
        }

        let edge_points: Vec<DVec3> = edges
            .par_iter()
            .map(|e| (centers[e.cells[0] as usize] + centers[e.cells[1] as usize]).normalize())
            .collect();
        let edge_normals: Vec<DVec3> = edges
            .par_iter()
            .zip(&edge_points)
            .map(|(e, &x)| {
                let d = centers[e.cells[1] as usize] - centers[e.cells[0] as usize];
                (d - x * d.dot(x)).normalize()
            })
            .collect();
        let edge_tangents: Vec<DVec3> = edge_points
            .iter()
            .zip(&edge_normals)
            .map(|(&x, &n)| x.cross(n))
            .collect();
        let edge_lengths: Vec<f64> = edges
            .par_iter()
            .map(|e| {
                arc_length(
                    vertices[e.vertices[0] as usize],
                    vertices[e.vertices[1] as usize],
                )
            })
            .collect();
        let dual_edge_lengths: Vec<f64> = edges
            .par_iter()
            .map(|e| arc_length(centers[e.cells[0] as usize], centers[e.cells[1] as usize]))
            .collect();

        let cell_areas = mesh.cell_areas();

        // Kites: the quadrilateral (cell centre, edge point before the vertex,
        // vertex, edge point after the vertex). Vertex `k` of a cell lies between
        // its edges `k` and `k + 1`.
        let mut kite_areas = vec![0.0; *ring_offsets.last().unwrap() as usize];
        let mut cell_edge_signs = vec![0.0; kite_areas.len()];
        for c in 0..mesh.n_cells() as CellId {
            let start = ring_offsets[c as usize] as usize;
            let ring = mesh.cell_vertices(c);
            let cell_edges = mesh.cell_edges(c);
            let m = ring.len();
            let xc = centers[c as usize];
            for k in 0..m {
                let xv = vertices[ring[k] as usize];
                let before = edge_points[cell_edges[k] as usize];
                let after = edge_points[cell_edges[(k + 1) % m] as usize];
                kite_areas[start + k] =
                    triangle_area(xc, before, xv) + triangle_area(xc, xv, after);
                cell_edge_signs[start + k] = if edges[cell_edges[k] as usize].cells[0] == c {
                    1.0
                } else {
                    -1.0
                };
            }
        }

        let mut vertex_kites = vec![[0.0; 3]; mesh.n_vertices()];
        for c in 0..mesh.n_cells() as CellId {
            let start = ring_offsets[c as usize] as usize;
            for (k, &v) in mesh.cell_vertices(c).iter().enumerate() {
                let slot = mesh.triangles()[v as usize]
                    .iter()
                    .position(|&x| x == c)
                    .expect("cell is a corner of its vertex's triangle");
                vertex_kites[v as usize][slot] = kite_areas[start + k];
            }
        }
        let vertex_areas: Vec<f64> = vertex_kites.iter().map(|k| k.iter().sum()).collect();

        let vertex_edge_signs: Vec<[f64; 3]> = (0..mesh.n_vertices() as VertexId)
            .map(|v| {
                let xv = vertices[v as usize];
                mesh.vertex_edges(v).map(|e| {
                    let ccw = xv.cross(edge_points[e as usize]);
                    if edge_normals[e as usize].dot(ccw) > 0.0 {
                        1.0
                    } else {
                        -1.0
                    }
                })
            })
            .collect();

        let mut geometry = Self {
            cell_areas,
            vertex_areas,
            edge_points,
            edge_normals,
            edge_tangents,
            edge_lengths,
            dual_edge_lengths,
            kite_areas,
            vertex_kites,
            cell_edge_signs,
            vertex_edge_signs,
            trisk_offsets: Vec::new(),
            trisk_edges: Vec::new(),
            trisk_weights: Vec::new(),
            ring_offsets,
        };
        geometry.compute_trisk_weights(mesh);
        geometry
    }

    /// TRiSK weights (Thuburn et al. 2009, eq. 33; Ringler et al. 2010).
    ///
    /// For edge `e` and each of its cells `i`, walk counter-clockwise around `i`
    /// from `e`. For each other edge `e'` reached, with `S` the sum of the kite
    /// fractions `R_{i,v}` of the vertices passed so far:
    ///
    /// `w_{e,e'} = n_{e,i} · n_{e',i} · (½ − S)`
    ///
    /// Kite fractions are normalized to sum to exactly 1 per cell, which makes the
    /// weights exactly antisymmetric (`w_{e,e'} = −w_{e',e}`), and that is what
    /// makes the Coriolis term energy-neutral.
    fn compute_trisk_weights(&mut self, mesh: &SphereMesh) {
        let edges = mesh.edges();
        let mut offsets = Vec::with_capacity(edges.len() + 1);
        offsets.push(0u32);
        let mut list_edges = Vec::new();
        let mut list_weights = Vec::new();

        for (e, edge) in edges.iter().enumerate() {
            for &c in &edge.cells {
                let start = self.ring_offsets[c as usize] as usize;
                let cell_edges = mesh.cell_edges(c);
                let m = cell_edges.len();
                let kites = &self.kite_areas[start..start + m];
                let kite_total: f64 = kites.iter().sum();
                let signs = &self.cell_edge_signs[start..start + m];
                let j = cell_edges
                    .iter()
                    .position(|&x| x as usize == e)
                    .expect("edge on cell");
                let mut s = 0.0;
                for step in 1..m {
                    // Vertex between edge (j + step − 1) and edge (j + step).
                    s += kites[(j + step - 1) % m] / kite_total;
                    let k = (j + step) % m;
                    list_edges.push(cell_edges[k]);
                    list_weights.push(signs[j] * signs[k] * (0.5 - s));
                }
            }
            offsets.push(list_edges.len() as u32);
        }

        self.trisk_offsets = offsets;
        self.trisk_edges = list_edges;
        self.trisk_weights = list_weights;
    }

    /// Kite areas of `cell`, aligned with [`SphereMesh::cell_vertices`].
    pub fn cell_kites(&self, cell: CellId) -> &[f64] {
        &self.kite_areas[self.ring_range(cell)]
    }

    /// Kite areas at `vertex`, aligned with that triangle's cells
    /// ([`SphereMesh::triangles`]).
    pub fn vertex_kites(&self, vertex: VertexId) -> [f64; 3] {
        self.vertex_kites[vertex as usize]
    }

    /// `n_{e,i}` for each edge of `cell`, aligned with [`SphereMesh::cell_edges`].
    pub fn cell_edge_signs(&self, cell: CellId) -> &[f64] {
        &self.cell_edge_signs[self.ring_range(cell)]
    }

    /// `t_{e,v}` for each edge of `vertex`, aligned with [`SphereMesh::vertex_edges`].
    pub fn vertex_edge_signs(&self, vertex: VertexId) -> [f64; 3] {
        self.vertex_edge_signs[vertex as usize]
    }

    /// TRiSK stencil of `edge`: the edges `e'` and weights `w_{e,e'}`.
    pub fn trisk_stencil(&self, edge: EdgeId) -> (&[EdgeId], &[f64]) {
        let r = self.trisk_offsets[edge as usize] as usize
            ..self.trisk_offsets[edge as usize + 1] as usize;
        (&self.trisk_edges[r.clone()], &self.trisk_weights[r])
    }

    fn ring_range(&self, cell: CellId) -> std::ops::Range<usize> {
        self.ring_offsets[cell as usize] as usize..self.ring_offsets[cell as usize + 1] as usize
    }
}
