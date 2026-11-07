//! Point location: which cell (or Delaunay triangle) contains a point.

use glam::DVec3;

use crate::coords::from_lat_lon;
use crate::mesh::{CellId, SphereMesh, VertexId};

/// Fast point → cell lookup for a [`SphereMesh`].
///
/// A Voronoi cell is the set of points nearest its generator, so locating a point
/// means finding its nearest generator. This walks greedily over the Delaunay
/// graph (always moving to a neighbour closer to the point), which is guaranteed
/// to reach the nearest generator on a Delaunay triangulation. A coarse
/// latitude/longitude table of starting cells keeps walks to a few steps.
///
/// The locator does not borrow the mesh; pass the same mesh it was built from.
#[derive(Clone, Debug)]
pub struct CellLocator {
    n_lat: usize,
    n_lon: usize,
    start: Vec<CellId>,
}

impl CellLocator {
    /// Build the start table for `mesh`.
    pub fn new(mesh: &SphereMesh) -> Self {
        // Roughly 16 cells per bucket keeps walks short and the table small.
        let n_lat = ((mesh.n_cells() as f64 / 32.0).sqrt().ceil() as usize).max(4);
        let n_lon = 2 * n_lat;
        let mut start = Vec::with_capacity(n_lat * n_lon);
        let mut current = 0;
        for i in 0..n_lat {
            let lat = bucket_center(i, n_lat, -std::f64::consts::FRAC_PI_2, std::f64::consts::PI);
            // Snake through rows so consecutive buckets are adjacent.
            for jj in 0..n_lon {
                let j = if i % 2 == 0 { jj } else { n_lon - 1 - jj };
                let lon = bucket_center(j, n_lon, -std::f64::consts::PI, std::f64::consts::TAU);
                current = walk(mesh, from_lat_lon(lat, lon), current);
                start.push((i, j, current));
            }
        }
        start.sort_unstable_by_key(|&(i, j, _)| (i, j));
        let start = start.into_iter().map(|(_, _, c)| c).collect();
        Self {
            n_lat,
            n_lon,
            start,
        }
    }

    /// The cell containing unit vector `p`.
    pub fn locate(&self, mesh: &SphereMesh, p: DVec3) -> CellId {
        walk(mesh, p, self.start_cell(p))
    }

    fn start_cell(&self, p: DVec3) -> CellId {
        let lat = p.z.clamp(-1.0, 1.0).asin();
        let lon = p.y.atan2(p.x);
        let i = bucket_index(
            lat,
            self.n_lat,
            -std::f64::consts::FRAC_PI_2,
            std::f64::consts::PI,
        );
        let j = bucket_index(
            lon,
            self.n_lon,
            -std::f64::consts::PI,
            std::f64::consts::TAU,
        );
        self.start[i * self.n_lon + j]
    }

    /// The Delaunay triangle (dual vertex) containing `p`, and the barycentric
    /// weights of its three cells (aligned with [`SphereMesh::triangles`]).
    ///
    /// Weights are planar barycentric coordinates of `p` projected onto the
    /// triangle's plane: non-negative, summing to 1, and exact for functions that
    /// are linear in 3D position.
    pub fn locate_triangle(&self, mesh: &SphereMesh, p: DVec3) -> (VertexId, [f64; 3]) {
        let nearest = self.locate(mesh, p);
        // The containing triangle almost always touches the nearest cell; if not,
        // it touches one of its neighbours.
        let candidates =
            std::iter::once(nearest).chain(mesh.cell_neighbors(nearest).iter().copied());
        let mut best = (VertexId::MAX, [0.0; 3], f64::NEG_INFINITY);
        for c in candidates {
            for &t in mesh.cell_vertices(c) {
                let w = barycentric(mesh, t, p);
                let worst = w[0].min(w[1]).min(w[2]);
                if worst >= 0.0 {
                    return (t, w);
                }
                if worst > best.2 {
                    best = (t, w, worst);
                }
            }
        }
        // Only reachable through round-off on an edge: clamp the closest candidate.
        let w = best.1.map(|x| x.max(0.0));
        let sum = w[0] + w[1] + w[2];
        (best.0, w.map(|x| x / sum))
    }
}

/// Planar barycentric coordinates of `p` with respect to triangle `t`.
fn barycentric(mesh: &SphereMesh, t: VertexId, p: DVec3) -> [f64; 3] {
    let [a, b, c] = mesh.triangles()[t as usize].map(|i| mesh.cell_centers()[i as usize]);
    let wa = p.dot(b.cross(c));
    let wb = p.dot(c.cross(a));
    let wc = p.dot(a.cross(b));
    let sum = wa + wb + wc;
    [wa / sum, wb / sum, wc / sum]
}

/// Greedy walk from `start` to the generator nearest to `p`.
pub fn walk(mesh: &SphereMesh, p: DVec3, start: CellId) -> CellId {
    let centers = mesh.cell_centers();
    let mut current = start;
    let mut best = p.dot(centers[current as usize]);
    loop {
        let mut next = current;
        for &n in mesh.cell_neighbors(current) {
            let d = p.dot(centers[n as usize]);
            if d > best {
                best = d;
                next = n;
            }
        }
        if next == current {
            return current;
        }
        current = next;
    }
}

fn bucket_center(i: usize, n: usize, min: f64, span: f64) -> f64 {
    min + span * (i as f64 + 0.5) / n as f64
}

fn bucket_index(x: f64, n: usize, min: f64, span: f64) -> usize {
    (((x - min) / span * n as f64) as usize).min(n - 1)
}
