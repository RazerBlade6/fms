//! Icosahedral Voronoi mesh on the unit sphere.
//!
//! The *primal* mesh is a Delaunay triangulation whose points are the cell
//! generators (cell centres). Its *dual* is the Voronoi tessellation: one cell per
//! point (12 pentagons, the rest hexagons), one Voronoi vertex per triangle (the
//! triangle's circumcentre), and one Voronoi edge per triangulation edge.
//!
//! All geometry is on the unit sphere; callers scale by the planet radius.

use glam::DVec3;
use rayon::prelude::*;

use crate::icosahedron;
use crate::spherical::{arc_length, circumcenter, polygon_moment, triangle_area};

/// Index of a Voronoi cell (equivalently, a triangulation point).
pub type CellId = u32;
/// Index of a Voronoi vertex (equivalently, a triangulation triangle).
pub type VertexId = u32;
/// Index of a Voronoi edge (equivalently, a triangulation edge).
pub type EdgeId = u32;

/// A Voronoi edge: the boundary between two adjacent cells.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Edge {
    /// The two cells this edge separates, with `cells[0] < cells[1]`.
    pub cells: [CellId; 2],
    /// The two Voronoi vertices at the ends of this edge. Walking from
    /// `vertices[0]` to `vertices[1]`, `cells[0]` is on the left (seen from
    /// outside the sphere).
    pub vertices: [VertexId; 2],
}

/// Options for Lloyd's algorithm (spherical centroidal Voronoi tessellation).
#[derive(Clone, Copy, Debug)]
pub struct LloydOptions {
    /// Upper bound on the number of iterations.
    pub max_iterations: u32,
    /// Stop when no generator moves further than this fraction of the mean cell
    /// spacing in one iteration. Relative, so it means the same at every resolution.
    pub tolerance: f64,
}

impl Default for LloydOptions {
    fn default() -> Self {
        Self {
            max_iterations: 5000,
            tolerance: 1e-6,
        }
    }
}

/// Outcome of a Lloyd optimization run.
#[derive(Clone, Copy, Debug)]
pub struct LloydReport {
    /// Iterations performed (at the finest level, for multilevel runs).
    pub iterations: u32,
    /// Largest generator displacement in the final iteration, as a fraction of the
    /// mean cell spacing.
    pub final_max_displacement: f64,
    /// Whether `final_max_displacement` fell below the tolerance.
    pub converged: bool,
}

/// Summary statistics of mesh quality.
#[derive(Clone, Copy, Debug)]
pub struct MeshQuality {
    /// Smallest cell area (unit sphere).
    pub min_area: f64,
    /// Largest cell area (unit sphere).
    pub max_area: f64,
    /// Largest distance between a generator and its cell's centroid, as a
    /// fraction of the mean cell spacing. Zero for a perfect centroidal Voronoi
    /// tessellation.
    pub max_centroid_offset: f64,
}

impl MeshQuality {
    /// Ratio of largest to smallest cell area (1.0 is perfectly uniform).
    pub fn area_ratio(&self) -> f64 {
        self.max_area / self.min_area
    }
}

/// Icosahedral Voronoi mesh on the unit sphere.
#[derive(Clone, Debug)]
pub struct SphereMesh {
    frequency: u32,
    /// Cell generators (unit vectors).
    cell_centers: Vec<DVec3>,
    /// Delaunay triangles (CCW from outside); triangle `t` is dual to vertex `t`.
    triangles: Vec<[CellId; 3]>,
    /// Voronoi vertices (triangle circumcentres, unit vectors).
    vertices: Vec<DVec3>,
    /// CSR offsets into `cell_vertices` / `cell_neighbors`; length `n_cells + 1`.
    ring_offsets: Vec<u32>,
    /// Voronoi vertices of each cell, CCW from outside.
    cell_vertices: Vec<VertexId>,
    /// Neighbours of each cell, CCW from outside. Neighbour `i` lies across the
    /// Voronoi edge from vertex `i - 1` to vertex `i` (cyclically).
    cell_neighbors: Vec<CellId>,
    /// Edges of each cell, aligned with `cell_neighbors`: edge `i` separates the
    /// cell from neighbour `i`.
    cell_edges: Vec<EdgeId>,
    /// Voronoi edges.
    edges: Vec<Edge>,
    /// The three edges meeting at each Voronoi vertex.
    vertex_edges: Vec<[EdgeId; 3]>,
}

impl SphereMesh {
    /// Build the geodesic mesh of the given frequency (`10·f² + 2` cells) without
    /// optimization. Use [`SphereMesh::optimize`] to make it centroidal.
    ///
    /// # Panics
    /// If `frequency` is zero or too large for `u32` indices.
    pub fn geodesic(frequency: u32) -> Self {
        Self::from_triangulation(frequency, icosahedron::subdivide(frequency))
    }

    /// Build a spherical centroidal Voronoi tessellation of the given frequency.
    ///
    /// Plain Lloyd iteration converges very slowly on fine meshes, because
    /// large-scale distortion only diffuses out a cell at a time. So this works
    /// multilevel: it optimizes a coarse mesh (frequency ≤ [`COARSEST_LEVEL`])
    /// to convergence, refines it by a small factor, re-optimizes, and repeats.
    /// Each finer level starts with the large scales already correct.
    ///
    /// The chain of levels comes from the small prime factors (≤ 7) of `frequency`,
    /// so frequencies like `400 = 25·2⁴` are fast. A frequency with a large prime
    /// factor falls back to plain Lloyd from that factor up, which is much slower.
    pub fn centroidal(frequency: u32, options: LloydOptions) -> (Self, LloydReport) {
        let levels = refinement_levels(frequency);
        let mut mesh = Self::geodesic(levels[0]);
        let mut report = mesh.optimize(options);
        for pair in levels.windows(2) {
            let factor = pair[1] / pair[0];
            let tri = icosahedron::refine(&mesh.triangulation(), factor);
            mesh = Self::from_triangulation(pair[1], tri);
            report = mesh.optimize(options);
        }
        (mesh.hilbert_ordered(), report)
    }

    /// Renumber cells and vertices along a Hilbert space-filling curve, so cells
    /// that are close on the sphere are close in memory. Geometry and topology
    /// are unchanged; only indices move. Edges follow the cell order.
    pub fn hilbert_ordered(&self) -> Self {
        let key = |p: DVec3| crate::hilbert::key(p);
        let mut order: Vec<CellId> = (0..self.n_cells() as CellId).collect();
        order.sort_by_key(|&c| key(self.cell_centers[c as usize]));
        let mut new_index = vec![0 as CellId; self.n_cells()];
        for (new, &old) in order.iter().enumerate() {
            new_index[old as usize] = new as CellId;
        }
        let points = order
            .iter()
            .map(|&c| self.cell_centers[c as usize])
            .collect();
        let mut triangles: Vec<(u64, [CellId; 3])> = self
            .triangles
            .iter()
            .zip(&self.vertices)
            .map(|(t, &v)| (key(v), t.map(|c| new_index[c as usize])))
            .collect();
        triangles.sort_by_key(|&(k, _)| k);
        let tri = icosahedron::Triangulation {
            points,
            triangles: triangles.into_iter().map(|(_, t)| t).collect(),
        };
        Self::from_triangulation(self.frequency, tri)
    }

    /// Build the Voronoi dual of a triangulation whose topology is that of the
    /// icosahedral lattice of the given `frequency`.
    fn from_triangulation(frequency: u32, tri: icosahedron::Triangulation) -> Self {
        let icosahedron::Triangulation { points, triangles } = tri;
        let n_cells = points.len();

        // Triangles incident to each point, as CSR.
        let mut counts = vec![0u32; n_cells + 1];
        for t in &triangles {
            for &p in t {
                counts[p as usize + 1] += 1;
            }
        }
        let mut ring_offsets = counts;
        for i in 0..n_cells {
            ring_offsets[i + 1] += ring_offsets[i];
        }
        let mut incident = vec![0 as VertexId; ring_offsets[n_cells] as usize];
        let mut fill = ring_offsets.clone();
        for (t, tri) in triangles.iter().enumerate() {
            for &p in tri {
                incident[fill[p as usize] as usize] = t as VertexId;
                fill[p as usize] += 1;
            }
        }

        // Order each ring CCW. For triangle (a, b, c) rotated so the centre `a` is
        // first, the next triangle CCW around `a` is the one whose `b` equals this
        // triangle's `c`. The neighbour shared by consecutive triangles is that `b`.
        let rotated = |t: VertexId, a: CellId| -> (CellId, CellId) {
            let [x, y, z] = triangles[t as usize];
            if x == a {
                (y, z)
            } else if y == a {
                (z, x)
            } else {
                (x, y)
            }
        };
        let mut cell_vertices = vec![0 as VertexId; incident.len()];
        let mut cell_neighbors = vec![0 as CellId; incident.len()];
        for a in 0..n_cells {
            let (start, end) = (ring_offsets[a] as usize, ring_offsets[a + 1] as usize);
            let ring = &incident[start..end];
            let centre = a as CellId;
            let mut current = ring[0];
            for slot in 0..ring.len() {
                cell_vertices[start + slot] = current;
                let (b, c) = rotated(current, centre);
                // Neighbour `slot + 1` lies between vertex `slot` and `slot + 1`.
                cell_neighbors[start + (slot + 1) % ring.len()] = c;
                if slot == 0 {
                    // The neighbour preceding vertex 0 is its `b`.
                    cell_neighbors[start] = b;
                }
                current = *ring
                    .iter()
                    .find(|&&t| rotated(t, centre).0 == c)
                    .expect("triangulation is a closed manifold");
            }
            debug_assert_eq!(current, ring[0], "ring around cell {a} does not close");
        }

        // One Voronoi edge per neighbour pair, recorded from the lower-indexed cell.
        let mut edges = Vec::with_capacity(incident.len() / 2);
        for a in 0..n_cells {
            let (start, end) = (ring_offsets[a] as usize, ring_offsets[a + 1] as usize);
            let len = end - start;
            for slot in 0..len {
                let b = cell_neighbors[start + slot];
                if (a as CellId) < b {
                    let from = cell_vertices[start + (slot + len - 1) % len];
                    let to = cell_vertices[start + slot];
                    edges.push(Edge {
                        cells: [a as CellId, b],
                        vertices: [from, to],
                    });
                }
            }
        }

        // Edge ids on each cell, aligned with its neighbour list.
        let mut cell_edges = vec![0 as EdgeId; incident.len()];
        let slot_of = |cell: CellId, neighbor: CellId| -> usize {
            let start = ring_offsets[cell as usize] as usize;
            let end = ring_offsets[cell as usize + 1] as usize;
            start
                + cell_neighbors[start..end]
                    .iter()
                    .position(|&n| n == neighbor)
                    .expect("adjacency is symmetric")
        };
        let mut vertex_edges = vec![[EdgeId::MAX; 3]; triangles.len()];
        for (id, e) in edges.iter().enumerate() {
            let [a, b] = e.cells;
            cell_edges[slot_of(a, b)] = id as EdgeId;
            cell_edges[slot_of(b, a)] = id as EdgeId;
            for v in e.vertices {
                let slots = &mut vertex_edges[v as usize];
                let free = slots
                    .iter()
                    .position(|&s| s == EdgeId::MAX)
                    .expect("3 edges per vertex");
                slots[free] = id as EdgeId;
            }
        }

        let mut mesh = Self {
            frequency,
            cell_centers: points,
            triangles,
            vertices: Vec::new(),
            ring_offsets,
            cell_vertices,
            cell_neighbors,
            cell_edges,
            edges,
            vertex_edges,
        };
        mesh.update_vertices();
        mesh
    }

    /// Recompute Voronoi vertices from the current generators.
    fn update_vertices(&mut self) {
        let centers = &self.cell_centers;
        self.vertices = self
            .triangles
            .par_iter()
            .map(|&[a, b, c]| {
                circumcenter(
                    centers[a as usize],
                    centers[b as usize],
                    centers[c as usize],
                )
            })
            .collect();
    }

    /// Run Lloyd's algorithm: repeatedly move each generator to the centroid of
    /// its Voronoi cell. Topology is kept fixed; check [`SphereMesh::is_delaunay`]
    /// afterwards to confirm it is still valid.
    pub fn optimize(&mut self, options: LloydOptions) -> LloydReport {
        let mut report = LloydReport {
            iterations: 0,
            final_max_displacement: f64::INFINITY,
            converged: false,
        };
        while report.iterations < options.max_iterations {
            let new_centers: Vec<DVec3> = (0..self.n_cells())
                .into_par_iter()
                .map(|c| self.centroid(c as CellId))
                .collect();
            let max_move = new_centers
                .par_iter()
                .zip(&self.cell_centers)
                .map(|(&n, &o)| arc_length(n, o))
                .reduce(|| 0.0, f64::max)
                / self.mean_spacing();
            self.cell_centers = new_centers;
            self.update_vertices();
            report.iterations += 1;
            report.final_max_displacement = max_move;
            if max_move < options.tolerance {
                report.converged = true;
                break;
            }
        }
        report
    }

    /// The primal triangulation (generators and Delaunay triangles).
    pub fn triangulation(&self) -> icosahedron::Triangulation {
        icosahedron::Triangulation {
            points: self.cell_centers.clone(),
            triangles: self.triangles.clone(),
        }
    }

    /// Mean distance between neighbouring generators on the unit sphere,
    /// estimated from the area of an average hexagon.
    pub fn mean_spacing(&self) -> f64 {
        // A regular hexagon with centre spacing d has area (√3 / 2)·d².
        let mean_area = 4.0 * std::f64::consts::PI / self.n_cells() as f64;
        (mean_area * 2.0 / 3f64.sqrt()).sqrt()
    }

    /// Subdivision frequency this mesh was built from.
    pub fn frequency(&self) -> u32 {
        self.frequency
    }

    /// Number of Voronoi cells.
    pub fn n_cells(&self) -> usize {
        self.cell_centers.len()
    }

    /// Number of Voronoi vertices (= Delaunay triangles).
    pub fn n_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Number of Voronoi edges.
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }

    /// Cell generators (unit vectors).
    pub fn cell_centers(&self) -> &[DVec3] {
        &self.cell_centers
    }

    /// Voronoi vertex positions (unit vectors).
    pub fn vertices(&self) -> &[DVec3] {
        &self.vertices
    }

    /// Delaunay triangles; triangle `t` is dual to Voronoi vertex `t`.
    pub fn triangles(&self) -> &[[CellId; 3]] {
        &self.triangles
    }

    /// Voronoi edges.
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }

    fn ring_range(&self, cell: CellId) -> std::ops::Range<usize> {
        self.ring_offsets[cell as usize] as usize..self.ring_offsets[cell as usize + 1] as usize
    }

    /// Voronoi vertices of `cell`, counter-clockwise from outside.
    pub fn cell_vertices(&self, cell: CellId) -> &[VertexId] {
        &self.cell_vertices[self.ring_range(cell)]
    }

    /// Neighbours of `cell`, counter-clockwise from outside. Neighbour `i` lies
    /// across the edge from vertex `i - 1` to vertex `i` (cyclically).
    pub fn cell_neighbors(&self, cell: CellId) -> &[CellId] {
        &self.cell_neighbors[self.ring_range(cell)]
    }

    /// Edges of `cell`, aligned with [`SphereMesh::cell_neighbors`]: edge `i`
    /// separates the cell from neighbour `i`.
    pub fn cell_edges(&self, cell: CellId) -> &[EdgeId] {
        &self.cell_edges[self.ring_range(cell)]
    }

    /// The three edges meeting at Voronoi vertex `vertex` (in no particular order).
    pub fn vertex_edges(&self, vertex: VertexId) -> [EdgeId; 3] {
        self.vertex_edges[vertex as usize]
    }

    /// Area of `cell` on the unit sphere.
    pub fn cell_area(&self, cell: CellId) -> f64 {
        let c = self.cell_centers[cell as usize];
        let ring = self.cell_vertices(cell);
        (0..ring.len())
            .map(|i| {
                let p = self.vertices[ring[i] as usize];
                let q = self.vertices[ring[(i + 1) % ring.len()] as usize];
                triangle_area(c, p, q)
            })
            .sum()
    }

    /// Areas of all cells on the unit sphere.
    pub fn cell_areas(&self) -> Vec<f64> {
        (0..self.n_cells())
            .into_par_iter()
            .map(|c| self.cell_area(c as CellId))
            .collect()
    }

    /// Centroid of `cell` (area-weighted, projected onto the sphere).
    pub fn centroid(&self, cell: CellId) -> DVec3 {
        let ring = self.cell_vertices(cell);
        polygon_moment(ring.iter().map(|&v| self.vertices[v as usize])).normalize()
    }

    /// Whether the triangulation satisfies the Delaunay property: no point lies
    /// strictly inside the circumcircle of any triangle it is adjacent to. This
    /// guarantees the stored topology really is the Voronoi diagram of the
    /// generators.
    pub fn is_delaunay(&self) -> bool {
        const EPS: f64 = 1e-12;
        (0..self.n_cells()).into_par_iter().all(|a| {
            let ring = self.cell_vertices(a as CellId);
            let neighbors = self.cell_neighbors(a as CellId);
            // Every neighbour of `a` must lie outside the circumcircles of all
            // triangles around `a`.
            ring.iter().all(|&t| {
                let n = self.vertices[t as usize];
                let threshold = n.dot(self.cell_centers[a]);
                let tri = self.triangles[t as usize];
                neighbors
                    .iter()
                    .filter(|b| !tri.contains(b))
                    .all(|&b| n.dot(self.cell_centers[b as usize]) <= threshold + EPS)
            })
        })
    }

    /// Area uniformity and centroidality statistics.
    pub fn quality(&self) -> MeshQuality {
        let areas = self.cell_areas();
        let max_centroid_offset = (0..self.n_cells())
            .into_par_iter()
            .map(|c| arc_length(self.centroid(c as CellId), self.cell_centers[c]))
            .reduce(|| 0.0, f64::max)
            / self.mean_spacing();
        MeshQuality {
            min_area: areas.iter().copied().fold(f64::INFINITY, f64::min),
            max_area: areas.iter().copied().fold(0.0, f64::max),
            max_centroid_offset,
        }
    }
}

/// Frequencies at or below this are optimized directly by Lloyd's algorithm.
pub const COARSEST_LEVEL: u32 = 32;

/// Chain of frequencies from coarsest to `frequency`, each dividing the next by a
/// prime factor ≤ 7.
fn refinement_levels(frequency: u32) -> Vec<u32> {
    let mut levels = vec![frequency];
    let mut k = frequency;
    while k > COARSEST_LEVEL {
        match [2, 3, 5, 7].into_iter().find(|&p| k.is_multiple_of(p)) {
            Some(p) => {
                k /= p;
                levels.push(k);
            }
            None => break,
        }
    }
    levels.reverse();
    levels
}

#[cfg(test)]
mod tests {
    use super::refinement_levels;

    #[test]
    fn refinement_chains() {
        assert_eq!(refinement_levels(16), vec![16]);
        assert_eq!(refinement_levels(400), vec![25, 50, 100, 200, 400]);
        assert_eq!(refinement_levels(45), vec![15, 45]);
        assert_eq!(refinement_levels(37), vec![37]);
        assert_eq!(refinement_levels(2 * 37), vec![37, 74]);
    }
}
