//! Which terrain cells are sea.
//!
//! A cell at or below sea level is sea only if its body of such cells
//! (joined through shared edges) is the planet's largest, the world ocean,
//! or covers at least [`MIN_SEA_FRACTION`] of the sphere. Smaller bodies are
//! dry land below sea level, like the Qattara Depression or Death Valley
//! (`docs/OBSERVATIONS.md` O-11); hydrology may fill them later.
//!
//! The size rule, rather than connection to the ocean alone, is because real
//! straits (the Bosporus, Gibraltar) are narrower than a terrain cell: the
//! seas behind them come out disconnected, and must stay sea.

use fms_grid::{CellId, SphereMesh};
use glam::DVec3;

/// Smallest body below sea level that is sea although it isn't the world
/// ocean, as a fraction of the sphere's surface: 5e-5 is about 25,500 km² on
/// Earth, between the largest dry depressions and the smallest enclosed seas
/// the terrain grid resolves.
pub const MIN_SEA_FRACTION: f64 = 5e-5;

/// One connected body of terrain cells at or below sea level.
#[derive(Clone, Debug, PartialEq)]
pub struct SeaBody {
    /// Its lowest-numbered cell, where the search found it.
    pub first_cell: CellId,
    /// Number of cells.
    pub cells: usize,
    /// Area as a fraction of the sphere's surface.
    pub area_fraction: f64,
    /// Area-weighted mean position of its cells, on the unit sphere.
    pub centre: DVec3,
    /// Its lowest elevation (m).
    pub lowest: f64,
    /// Whether this is the world ocean: the largest body (the lowest-numbered
    /// on a tie).
    pub world_ocean: bool,
    /// Whether its cells are sea: the world ocean, or at least the minimum
    /// size.
    pub sea: bool,
}

/// The sea mask of a terrain grid and the bodies below sea level it came
/// from.
#[derive(Clone, Debug, PartialEq)]
pub struct SeaMask {
    /// Whether each cell is sea.
    pub sea: Vec<bool>,
    /// Every body of cells at or below sea level, in order of their first
    /// cell.
    pub bodies: Vec<SeaBody>,
}

impl SeaMask {
    /// The sea mask of `elevation` (m, one per cell of `mesh`) with the
    /// default [`MIN_SEA_FRACTION`]. `cell_areas` are the cells' areas in any
    /// unit; their sum is the sphere.
    pub fn new(elevation: &[f64], mesh: &SphereMesh, cell_areas: &[f64]) -> Self {
        Self::with_min_fraction(elevation, mesh, cell_areas, MIN_SEA_FRACTION)
    }

    /// [`new`](Self::new) with bodies of at least `min_fraction` of the
    /// sphere kept as sea.
    ///
    /// A flood fill from each unvisited cell in index order, so the result is
    /// fixed by the mesh and the elevations alone, in time linear in the
    /// number of cells.
    ///
    /// # Panics
    /// If the lengths differ from the mesh's cells.
    pub fn with_min_fraction(
        elevation: &[f64],
        mesh: &SphereMesh,
        cell_areas: &[f64],
        min_fraction: f64,
    ) -> Self {
        let n = mesh.n_cells();
        assert_eq!(elevation.len(), n, "one elevation per cell");
        assert_eq!(cell_areas.len(), n, "one area per cell");
        let sphere: f64 = cell_areas.iter().sum();
        let below = |c: usize| elevation[c] <= 0.0;
        let centres = mesh.cell_centers();

        const NONE: u32 = u32::MAX;
        let mut body = vec![NONE; n];
        let mut bodies = Vec::new();
        let mut stack: Vec<CellId> = Vec::new();
        for start in 0..n {
            if !below(start) || body[start] != NONE {
                continue;
            }
            let id = bodies.len() as u32;
            body[start] = id;
            stack.push(start as CellId);
            let (mut cells, mut area, mut centre) = (0, 0.0, DVec3::ZERO);
            let mut lowest = f64::INFINITY;
            while let Some(c) = stack.pop() {
                cells += 1;
                lowest = lowest.min(elevation[c as usize]);
                area += cell_areas[c as usize];
                centre += cell_areas[c as usize] * centres[c as usize];
                for &next in mesh.cell_neighbors(c) {
                    if below(next as usize) && body[next as usize] == NONE {
                        body[next as usize] = id;
                        stack.push(next);
                    }
                }
            }
            bodies.push(SeaBody {
                first_cell: start as CellId,
                cells,
                area_fraction: area / sphere,
                centre: centre.normalize_or_zero(),
                lowest,
                world_ocean: false,
                sea: false,
            });
        }
        // The first largest, so a tie goes to the lowest-numbered.
        let largest = bodies
            .iter()
            .enumerate()
            .fold(None, |best: Option<(usize, f64)>, (i, b)| match best {
                Some((_, a)) if a >= b.area_fraction => best,
                _ => Some((i, b.area_fraction)),
            });
        if let Some((i, _)) = largest {
            bodies[i].world_ocean = true;
        }
        for b in &mut bodies {
            b.sea = b.world_ocean || b.area_fraction >= min_fraction;
        }
        let sea = body
            .iter()
            .map(|&b| b != NONE && bodies[b as usize].sea)
            .collect();
        Self { sea, bodies }
    }

    /// Whether every cell at or below sea level is sea, as it is when the
    /// world ocean is the only body.
    pub fn is_whole_ocean(&self) -> bool {
        self.bodies.iter().all(|b| b.sea)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fms_grid::coords::from_lat_lon;
    use fms_grid::{LloydOptions, MeshGeometry};

    struct Grid {
        mesh: SphereMesh,
        geometry: MeshGeometry,
    }

    fn grid(frequency: u32) -> Grid {
        let (mesh, _) = SphereMesh::centroidal(frequency, LloydOptions::default());
        let geometry = MeshGeometry::new(&mesh);
        Grid { mesh, geometry }
    }

    /// Angular radius (radians) of a cap covering `fraction` of the sphere.
    fn cap_radius(fraction: f64) -> f64 {
        fms_math::acos(1.0 - 2.0 * fraction)
    }

    fn angle(a: DVec3, b: DVec3) -> f64 {
        fms_math::atan2(a.cross(b).length(), a.dot(b))
    }

    fn mask(g: &Grid, elevation: &[f64], min_fraction: f64) -> SeaMask {
        SeaMask::with_min_fraction(elevation, &g.mesh, &g.geometry.cell_areas, min_fraction)
    }

    /// The sea mask's verdict on the cell nearest `p`.
    fn sea_at(g: &Grid, m: &SeaMask, p: DVec3) -> bool {
        let nearest = (0..g.mesh.n_cells())
            .min_by(|&a, &b| {
                angle(g.mesh.cell_centers()[a], p).total_cmp(&angle(g.mesh.cell_centers()[b], p))
            })
            .unwrap();
        m.sea[nearest]
    }

    /// Ocean outside a continental cap of `land` radians around (0°, 0°),
    /// with depressions (centre, radius) inside it.
    fn continent_with(g: &Grid, land: f64, holes: &[(DVec3, f64)]) -> Vec<f64> {
        let centre = from_lat_lon(0.0, 0.0);
        g.mesh
            .cell_centers()
            .iter()
            .map(|&p| {
                if angle(p, centre) > land {
                    -4000.0
                } else if holes.iter().any(|&(c, r)| angle(p, c) <= r) {
                    -50.0
                } else {
                    300.0
                }
            })
            .collect()
    }

    #[test]
    fn all_below_sea_level_is_all_sea() {
        let g = grid(16);
        let m = mask(&g, &vec![-10.0; g.mesh.n_cells()], MIN_SEA_FRACTION);
        assert!(m.sea.iter().all(|&s| s));
        assert_eq!(m.bodies.len(), 1);
        assert!(m.bodies[0].world_ocean && m.is_whole_ocean());
        assert!((m.bodies[0].area_fraction - 1.0).abs() < 1e-12);
        // Exactly zero is sea level, and counts as below.
        let m = mask(&g, &vec![0.0; g.mesh.n_cells()], MIN_SEA_FRACTION);
        assert!(m.sea.iter().all(|&s| s));
        let m = mask(&g, &vec![1.0; g.mesh.n_cells()], MIN_SEA_FRACTION);
        assert!(m.sea.iter().all(|&s| !s) && m.bodies.is_empty());
    }

    #[test]
    fn the_idealized_continent_is_the_old_rule() {
        let g = grid(32);
        let elevation = crate::Continent::default().elevation(g.mesh.cell_centers());
        let m = mask(&g, &elevation, MIN_SEA_FRACTION);
        assert_eq!(m.bodies.len(), 1);
        for (s, z) in m.sea.iter().zip(&elevation) {
            assert_eq!(*s, *z <= 0.0);
        }
    }

    #[test]
    fn a_small_depression_is_land_and_a_large_enclosed_sea_is_sea() {
        // At frequency 32 a cell is about 1e-4 of the sphere, so the
        // threshold here is 1e-2: a small depression of 0.3e-2 and a large
        // enclosed sea of 3e-2 inside a continent.
        let g = grid(32);
        let min = 1e-2;
        let small = from_lat_lon(0.0, 0.5);
        let large = from_lat_lon(0.0, -0.4);
        let elevation = continent_with(
            &g,
            1.2,
            &[(small, cap_radius(0.3e-2)), (large, cap_radius(3e-2))],
        );
        let m = mask(&g, &elevation, min);
        assert_eq!(m.bodies.len(), 3);
        assert!(!sea_at(&g, &m, small));
        assert!(sea_at(&g, &m, large));
        assert!(sea_at(&g, &m, from_lat_lon(0.0, std::f64::consts::PI)));
        assert!(!m.is_whole_ocean());
        // Every land cell by the old rule is still land.
        for (s, z) in m.sea.iter().zip(&elevation) {
            assert!(!*s || *z <= 0.0);
        }
    }

    #[test]
    fn a_depression_joined_to_the_ocean_by_a_one_cell_channel_is_sea() {
        let g = grid(32);
        let min = 1e-2;
        let hole = from_lat_lon(0.0, 0.3);
        let mut elevation = continent_with(&g, 1.0, &[(hole, cap_radius(0.3e-2))]);
        assert!(!sea_at(&g, &mask(&g, &elevation, min), hole));
        // Cut a channel one cell wide: the shortest path through the cell
        // graph from the depression to the ocean.
        let n = g.mesh.n_cells();
        let mut from = vec![u32::MAX; n];
        let mut queue: std::collections::VecDeque<usize> = (0..n)
            .filter(|&c| elevation[c] == -50.0)
            .inspect(|&c| from[c] = c as u32)
            .collect();
        let end = loop {
            let c = queue.pop_front().expect("the ocean is reachable");
            if elevation[c] == -4000.0 {
                break c;
            }
            for &next in g.mesh.cell_neighbors(c as CellId) {
                if from[next as usize] == u32::MAX {
                    from[next as usize] = c as u32;
                    queue.push_back(next as usize);
                }
            }
        };
        let mut c = from[end] as usize;
        let mut length = 0;
        while elevation[c] != -50.0 {
            elevation[c] = -20.0;
            c = from[c] as usize;
            length += 1;
        }
        assert!(length >= 2, "channel of {length} cells");
        let m = mask(&g, &elevation, min);
        assert!(sea_at(&g, &m, hole));
        assert!(m.is_whole_ocean() && m.bodies.len() == 1);
    }

    #[test]
    fn two_oceans_separated_by_land_are_both_sea() {
        // A land belt around the equator; the smaller ocean, south of it, is
        // still far above the threshold.
        let g = grid(16);
        let elevation: Vec<f64> = g
            .mesh
            .cell_centers()
            .iter()
            .map(|p| {
                if (-0.3..0.1).contains(&p.z) {
                    100.0
                } else {
                    -3000.0
                }
            })
            .collect();
        let m = mask(&g, &elevation, MIN_SEA_FRACTION);
        assert_eq!(m.bodies.len(), 2);
        assert_eq!(m.bodies.iter().filter(|b| b.world_ocean).count(), 1);
        let north = m.bodies.iter().find(|b| b.centre.z > 0.0).unwrap();
        assert!(north.world_ocean);
        assert!(m.is_whole_ocean());
    }

    #[test]
    fn the_threshold_is_on_area_whatever_the_resolution() {
        // The same shapes at two resolutions: depressions of 0.6 and 1.6
        // times the threshold, which needs a few dozen cells at the coarser.
        let min = 1e-2;
        let below = from_lat_lon(0.1, 0.4);
        let above = from_lat_lon(-0.1, -0.4);
        let holes = [
            (below, cap_radius(0.6 * min)),
            (above, cap_radius(1.6 * min)),
        ];
        for frequency in [32, 64] {
            let g = grid(frequency);
            let m = mask(&g, &continent_with(&g, 1.2, &holes), min);
            assert!(!sea_at(&g, &m, below), "frequency {frequency}");
            assert!(sea_at(&g, &m, above), "frequency {frequency}");
            for b in m.bodies.iter().filter(|b| !b.world_ocean) {
                let expected = if b.centre.dot(below) > b.centre.dot(above) {
                    0.6
                } else {
                    1.6
                };
                // The cells' areas approximate the cap's.
                assert!(
                    (b.area_fraction / min / expected - 1.0).abs() < 0.2,
                    "frequency {frequency}: {} of the threshold, expected {expected}",
                    b.area_fraction / min
                );
            }
        }
    }

    #[test]
    fn the_cell_areas_unit_does_not_matter() {
        // Areas in m² on a planet of any radius give the same mask.
        let g = grid(32);
        let min = 1e-2;
        let elevation = continent_with(
            &g,
            1.2,
            &[
                (from_lat_lon(0.0, 0.5), cap_radius(0.8 * min)),
                (from_lat_lon(0.0, -0.4), cap_radius(1.2 * min)),
            ],
        );
        let unit = mask(&g, &elevation, min);
        for radius in [6.371e6, 3.39e6] {
            let areas: Vec<f64> = g
                .geometry
                .cell_areas
                .iter()
                .map(|a| a * radius * radius)
                .collect();
            let scaled = SeaMask::with_min_fraction(&elevation, &g.mesh, &areas, min);
            assert_eq!(scaled.sea, unit.sea);
        }
    }
}
