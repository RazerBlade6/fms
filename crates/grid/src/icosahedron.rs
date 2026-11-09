//! The base icosahedron and its geodesic subdivision into a triangulation of the
//! unit sphere.

use std::collections::HashMap;

use glam::DVec3;

/// Number of vertices of an icosahedron.
pub const VERTICES: usize = 12;
/// Number of edges of an icosahedron.
pub const EDGES: usize = 30;
/// Number of faces of an icosahedron.
pub const FACES: usize = 20;

/// Vertices of the regular icosahedron, normalized to the unit sphere.
pub fn vertices() -> [DVec3; VERTICES] {
    let phi = (1.0 + 5f64.sqrt()) / 2.0;
    [
        DVec3::new(-1.0, phi, 0.0),
        DVec3::new(1.0, phi, 0.0),
        DVec3::new(-1.0, -phi, 0.0),
        DVec3::new(1.0, -phi, 0.0),
        DVec3::new(0.0, -1.0, phi),
        DVec3::new(0.0, 1.0, phi),
        DVec3::new(0.0, -1.0, -phi),
        DVec3::new(0.0, 1.0, -phi),
        DVec3::new(phi, 0.0, -1.0),
        DVec3::new(phi, 0.0, 1.0),
        DVec3::new(-phi, 0.0, -1.0),
        DVec3::new(-phi, 0.0, 1.0),
    ]
    .map(DVec3::normalize)
}

/// Faces of the icosahedron, each counter-clockwise when seen from outside.
pub fn faces() -> [[u32; 3]; FACES] {
    let v = vertices();
    let mut faces = [
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    // Enforce outward (CCW) orientation rather than trusting the table.
    for f in &mut faces {
        let [a, b, c] = f.map(|i| v[i as usize]);
        if (b - a).cross(c - a).dot(a + b + c) < 0.0 {
            f.swap(1, 2);
        }
    }
    faces
}

/// A triangulation of the unit sphere.
#[derive(Clone, Debug)]
pub struct Triangulation {
    /// Unit-vector positions of the points (the future Voronoi cell generators).
    pub points: Vec<DVec3>,
    /// Triangles as point indices, counter-clockwise seen from outside.
    pub triangles: Vec<[u32; 3]>,
}

impl Triangulation {
    /// The icosahedron itself (frequency 1).
    pub fn icosahedron() -> Self {
        Self {
            points: vertices().to_vec(),
            triangles: faces().to_vec(),
        }
    }
}

/// Number of points produced by a geodesic subdivision of the given frequency.
pub fn point_count(frequency: u32) -> usize {
    let k = frequency as usize;
    10 * k * k + 2
}

/// Subdivide the icosahedron with the given `frequency` (≥ 1).
///
/// Each face is split into a triangular lattice with `frequency` segments per edge,
/// giving `10·frequency² + 2` points and `20·frequency²` triangles.
///
/// # Panics
/// If `frequency` is zero or so large that point indices overflow `u32`.
pub fn subdivide(frequency: u32) -> Triangulation {
    refine(&Triangulation::icosahedron(), frequency)
}

/// Split every triangle of `base` into a lattice with `factor` segments per edge,
/// projecting new points onto the sphere. Points shared between triangles (base
/// points and points on base edges) are emitted exactly once, and base points keep
/// their indices.
///
/// Refining an icosahedral mesh of frequency `k` by `factor` gives a mesh with the
/// same topology as frequency `k·factor`.
///
/// # Panics
/// If `factor` is zero or the result would overflow `u32` indices.
pub fn refine(base: &Triangulation, factor: u32) -> Triangulation {
    assert!(factor >= 1, "refinement factor must be at least 1");
    let m = factor;

    // Canonical numbering of base edges, keyed by (min, max).
    let mut edge_ids: HashMap<(u32, u32), u32> = HashMap::new();
    for t in &base.triangles {
        for (u, v) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let next = edge_ids.len() as u32;
            edge_ids.entry((u.min(v), u.max(v))).or_insert(next);
        }
    }

    let (m_us, n_base) = (m as usize, base.points.len());
    let n_points = n_base
        + edge_ids.len() * (m_us - 1)
        + base.triangles.len() * (m_us - 1) * m_us.saturating_sub(2) / 2;
    assert!(
        u32::try_from(n_points).is_ok(),
        "refinement produces too many points"
    );

    let mf = f64::from(m);
    let lerp = |p: DVec3, q: DVec3, t: f64| p + (q - p) * t;

    let mut points = vec![DVec3::ZERO; n_points];
    points[..n_base].copy_from_slice(&base.points);

    // Point index for step `s` (0 < s < m) along the base edge u → v.
    let edge_base = n_base as u32;
    let edge_point = |u: u32, v: u32, s: u32| -> u32 {
        let along = if u < v { s } else { m - s };
        edge_base + edge_ids[&(u.min(v), u.max(v))] * (m - 1) + (along - 1)
    };
    for (&(u, v), &id) in &edge_ids {
        for s in 1..m {
            let p = lerp(
                base.points[u as usize],
                base.points[v as usize],
                f64::from(s) / mf,
            );
            points[(edge_base + id * (m - 1) + (s - 1)) as usize] = p.normalize();
        }
    }

    let mut next_interior = edge_base + edge_ids.len() as u32 * (m - 1);
    let mut triangles = Vec::with_capacity(base.triangles.len() * m_us * m_us);
    // Per-triangle lattice index map, reused: lattice[i][j] for i + j ≤ m.
    let mut lattice = vec![0u32; (m_us + 1) * (m_us + 1)];
    let at = |i: u32, j: u32| (i * (m + 1) + j) as usize;

    for &[a, b, c] in &base.triangles {
        let (pa, pb, pc) = (
            base.points[a as usize],
            base.points[b as usize],
            base.points[c as usize],
        );
        for i in 0..=m {
            for j in 0..=(m - i) {
                let index = if i == 0 && j == 0 {
                    a
                } else if i == m {
                    b
                } else if j == m {
                    c
                } else if j == 0 {
                    edge_point(a, b, i)
                } else if i == 0 {
                    edge_point(a, c, j)
                } else if i + j == m {
                    edge_point(b, c, j)
                } else {
                    let p = pa + (pb - pa) * (f64::from(i) / mf) + (pc - pa) * (f64::from(j) / mf);
                    let id = next_interior;
                    next_interior += 1;
                    points[id as usize] = p.normalize();
                    id
                };
                lattice[at(i, j)] = index;
            }
        }
        for i in 0..m {
            for j in 0..(m - i) {
                triangles.push([
                    lattice[at(i, j)],
                    lattice[at(i + 1, j)],
                    lattice[at(i, j + 1)],
                ]);
                if i + j + 1 < m {
                    triangles.push([
                        lattice[at(i + 1, j)],
                        lattice[at(i + 1, j + 1)],
                        lattice[at(i, j + 1)],
                    ]);
                }
            }
        }
    }
    debug_assert_eq!(next_interior as usize, n_points);

    Triangulation { points, triangles }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn faces_are_outward() {
        let v = vertices();
        for [a, b, c] in faces().map(|f| f.map(|i| v[i as usize])) {
            assert!((b - a).cross(c - a).dot(a + b + c) > 0.0);
        }
    }

    #[test]
    fn base_has_thirty_unique_edges() {
        let mut edges = HashSet::new();
        for f in faces() {
            for (u, v) in [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])] {
                edges.insert((u.min(v), u.max(v)));
            }
        }
        assert_eq!(edges.len(), EDGES);
    }

    #[test]
    fn counts_match_formula() {
        for k in [1, 2, 3, 7, 16] {
            let t = subdivide(k);
            assert_eq!(t.points.len(), point_count(k));
            assert_eq!(t.triangles.len(), 20 * (k as usize).pow(2));
        }
    }

    #[test]
    fn points_are_unique_and_on_sphere() {
        let t = subdivide(9);
        for p in &t.points {
            assert!((p.length() - 1.0).abs() < 1e-14);
        }
        let mut sorted: Vec<[i64; 3]> = t
            .points
            .iter()
            .map(|p| (*p * 1e9).round().as_i64vec3().to_array())
            .collect();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), t.points.len(), "duplicate points");
    }

    #[test]
    fn refining_matches_direct_subdivision_counts() {
        let refined = refine(&subdivide(3), 4);
        let direct = subdivide(12);
        assert_eq!(refined.points.len(), direct.points.len());
        assert_eq!(refined.triangles.len(), direct.triangles.len());
    }

    #[test]
    fn triangles_are_outward() {
        let t = subdivide(6);
        for tri in &t.triangles {
            let [a, b, c] = tri.map(|i| t.points[i as usize]);
            assert!((b - a).cross(c - a).dot(a + b + c) > 0.0);
        }
    }
}
