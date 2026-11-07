//! Geometry primitives on the unit sphere.

use glam::DVec3;

/// Circumcentre of the spherical triangle `a, b, c` (counter-clockwise when seen
/// from outside the sphere), as a unit vector.
///
/// This is the point equidistant from all three vertices, i.e. the Voronoi vertex
/// dual to the triangle.
#[inline]
pub fn circumcenter(a: DVec3, b: DVec3, c: DVec3) -> DVec3 {
    (b - a).cross(c - a).normalize()
}

/// Area (spherical excess) of the spherical triangle `a, b, c` on the unit sphere.
///
/// Uses the Van Oosterom–Strackee formula, which is numerically stable for the
/// small triangles found on fine meshes. Always non-negative.
#[inline]
pub fn triangle_area(a: DVec3, b: DVec3, c: DVec3) -> f64 {
    let numerator = a.dot(b.cross(c)).abs();
    let denominator = 1.0 + a.dot(b) + b.dot(c) + c.dot(a);
    2.0 * numerator.atan2(denominator)
}

/// Great-circle distance (central angle, radians) between two unit vectors.
#[inline]
pub fn arc_length(a: DVec3, b: DVec3) -> f64 {
    // atan2 form is accurate for both tiny and near-antipodal angles.
    a.cross(b).length().atan2(a.dot(b))
}

/// First moment `∫ x dA` of a spherical polygon whose vertices `ring` are listed
/// counter-clockwise (seen from outside) and joined by great-circle arcs.
///
/// Exact for geodesic polygons: by Stokes' theorem the integral reduces to
/// `½ Σ θᵢ n̂ᵢ`, where `θᵢ` is the arc length of edge `i` and `n̂ᵢ` its unit normal.
/// Normalizing the result gives the polygon's centroid projected onto the sphere.
pub fn polygon_moment(ring: impl IntoIterator<Item = DVec3> + Clone) -> DVec3 {
    let first = ring.clone().into_iter().next();
    let mut sum = DVec3::ZERO;
    let mut iter = ring.into_iter().peekable();
    while let Some(p) = iter.next() {
        let q = match iter.peek() {
            Some(&q) => q,
            None => match first {
                Some(f) => f,
                None => break,
            },
        };
        let n = p.cross(q);
        let len = n.length();
        if len > 0.0 {
            sum += n * (arc_length(p, q) / len);
        }
    }
    sum * 0.5
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn octant_triangle_area_is_eighth_of_sphere() {
        let area = triangle_area(DVec3::X, DVec3::Y, DVec3::Z);
        assert!((area - 4.0 * PI / 8.0).abs() < 1e-14);
    }

    #[test]
    fn octant_circumcenter_is_diagonal() {
        let c = circumcenter(DVec3::X, DVec3::Y, DVec3::Z);
        assert!((c - DVec3::ONE.normalize()).length() < 1e-14);
    }

    #[test]
    fn arc_length_quarter_turn() {
        assert!((arc_length(DVec3::X, DVec3::Y) - FRAC_PI_2).abs() < 1e-15);
    }

    #[test]
    fn polygon_moment_matches_cap_integral() {
        // Regular polygon approximating a cap of colatitude α around +z. The
        // exact moment of a cap is π sin²α ẑ; a 4096-gon is close to the cap.
        let alpha: f64 = 0.3;
        let n = 4096;
        let ring: Vec<DVec3> = (0..n)
            .map(|i| {
                let phi = 2.0 * PI * f64::from(i) / f64::from(n);
                DVec3::new(
                    alpha.sin() * phi.cos(),
                    alpha.sin() * phi.sin(),
                    alpha.cos(),
                )
            })
            .collect();
        let m = polygon_moment(ring.iter().copied());
        assert!(m.x.abs() < 1e-12 && m.y.abs() < 1e-12);
        assert!((m.z - PI * alpha.sin().powi(2)).abs() < 1e-5);
    }

    #[test]
    fn polygon_moment_of_octant_is_symmetric() {
        let m = polygon_moment([DVec3::X, DVec3::Y, DVec3::Z]);
        // Exact value: each component equals π/4.
        assert!((m - DVec3::splat(PI / 4.0)).length() < 1e-14);
    }
}
