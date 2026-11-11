//! Ocean basins on the climate mesh, for the wind-driven gyres.
//!
//! The gyres' streamfunction lives on the mesh's vertices, is held at zero on
//! coasts and is solved in each basin separately; the flow across an edge is
//! the difference of the streamfunction between its two vertices. Here:
//!
//! - a cell is **wet** if at least [`WET_OCEAN_FRACTION`] of it is ocean.
//!   This binary mask is the gyres' alone: Ekman and eddy transport and
//!   diffusion keep their fractional openings;
//! - a vertex is **interior** if its three cells are wet, and **coast**
//!   otherwise, so the flow across any edge with a dry cell is zero;
//! - a **basin** is a connected set of interior vertices, joined through mesh
//!   edges. Basins with fewer than [`MIN_BASIN_VERTICES`] are too small to
//!   hold a gyre and are dropped: their vertices count as coast.
//!
//! Everything is a function of the mesh and the land fraction, so it's fixed
//! for a run.

use fms_grid::{CellId, Operators, VertexId};

/// Smallest ocean fraction of a wet cell.
pub const WET_OCEAN_FRACTION: f64 = 0.5;

/// Fewest interior vertices a basin needs to be kept.
pub const MIN_BASIN_VERTICES: usize = 12;

/// One ocean basin: a connected set of interior vertices.
#[derive(Clone, Debug, PartialEq)]
pub struct Basin {
    /// Its interior vertices, in ascending order.
    pub vertices: Vec<VertexId>,
    /// Its area, the sum of its vertices' dual areas (m²).
    pub area: f64,
    /// Whether any of its vertices has a coast vertex as a neighbour. A basin
    /// covering a whole planet without land (an aquaplanet) has none.
    pub touches_coast: bool,
    /// The number of separate landmasses bordering it: connected sets of
    /// cells that aren't wet (joined through shared edges) with a cell at a
    /// coast vertex next to the basin. A landmass is counted once however
    /// long its coast, and a landmass reaching two basins counts in both.
    pub landmasses: usize,
}

/// The ocean basins of one mesh and land fraction, largest first.
#[derive(Clone, Debug, PartialEq)]
pub struct Basins {
    basins: Vec<Basin>,
    dropped: Vec<Basin>,
    wet: Vec<bool>,
    vertex_basin: Vec<Option<u32>>,
    cell_basin: Vec<Option<u32>>,
}

impl Basins {
    /// The basins on the mesh of `ops`, with land covering `land_fraction`
    /// of each cell, dropping those with fewer than [`MIN_BASIN_VERTICES`].
    pub fn new(ops: &Operators, land_fraction: &[f64]) -> Self {
        Self::with_min_vertices(ops, land_fraction, MIN_BASIN_VERTICES)
    }

    /// [`new`](Self::new), dropping basins with fewer than `min_vertices`
    /// interior vertices.
    ///
    /// Basins are numbered by area, largest first, a tie going to the one
    /// with the lowest-numbered vertex. Every search runs in index order, so
    /// the result depends on nothing but the mesh and the land fraction.
    ///
    /// # Panics
    /// If `land_fraction` doesn't have one value per cell.
    pub fn with_min_vertices(ops: &Operators, land_fraction: &[f64], min_vertices: usize) -> Self {
        let (mesh, geometry) = (ops.mesh(), ops.geometry());
        let (n_cells, n_vertices) = (mesh.n_cells(), mesh.n_vertices());
        assert_eq!(land_fraction.len(), n_cells, "one land fraction per cell");
        let wet: Vec<bool> = land_fraction
            .iter()
            .map(|f| 1.0 - f >= WET_OCEAN_FRACTION)
            .collect();
        let interior: Vec<bool> = mesh
            .triangles()
            .iter()
            .map(|cells| cells.iter().all(|&c| wet[c as usize]))
            .collect();
        // The vertex at the other end of each of a vertex's edges.
        let neighbours = |v: VertexId| {
            mesh.vertex_edges(v).map(|e| {
                let ends = mesh.edges()[e as usize].vertices;
                if ends[0] == v { ends[1] } else { ends[0] }
            })
        };

        // Landmasses: connected sets of cells that aren't wet.
        const NONE: u32 = u32::MAX;
        let mut landmass = vec![NONE; n_cells];
        let mut stack: Vec<u32> = Vec::new();
        let mut count = 0;
        for start in 0..n_cells {
            if wet[start] || landmass[start] != NONE {
                continue;
            }
            landmass[start] = count;
            stack.push(start as CellId);
            while let Some(c) = stack.pop() {
                for &next in mesh.cell_neighbors(c) {
                    if !wet[next as usize] && landmass[next as usize] == NONE {
                        landmass[next as usize] = count;
                        stack.push(next);
                    }
                }
            }
            count += 1;
        }

        // Connected sets of interior vertices.
        let mut component = vec![NONE; n_vertices];
        let mut found: Vec<Basin> = Vec::new();
        let r2 = ops.radius() * ops.radius();
        for start in 0..n_vertices {
            if !interior[start] || component[start] != NONE {
                continue;
            }
            let id = found.len() as u32;
            component[start] = id;
            stack.push(start as VertexId);
            let mut vertices = Vec::new();
            while let Some(v) = stack.pop() {
                vertices.push(v);
                for w in neighbours(v) {
                    if interior[w as usize] && component[w as usize] == NONE {
                        // The edge's two cells are cells of both its ends,
                        // so an edge between interior vertices is between
                        // wet cells.
                        debug_assert!(mesh.vertex_edges(v).iter().any(|&e| {
                            let edge = mesh.edges()[e as usize];
                            edge.vertices.contains(&w)
                                && edge.cells.iter().all(|&c| wet[c as usize])
                        }));
                        component[w as usize] = id;
                        stack.push(w);
                    }
                }
            }
            vertices.sort_unstable();
            let area = vertices
                .iter()
                .map(|&v| geometry.vertex_areas[v as usize] * r2)
                .sum();
            // Coast vertices next to the basin, and the land at them: an
            // interior neighbour is in this basin, so every other neighbour
            // is coast, and its cell off the shared edge is dry.
            let mut lands: Vec<u32> = Vec::new();
            let mut touches_coast = false;
            for &v in &vertices {
                for w in neighbours(v).into_iter().filter(|&w| !interior[w as usize]) {
                    touches_coast = true;
                    lands.extend(
                        mesh.triangles()[w as usize]
                            .iter()
                            .filter(|&&c| !wet[c as usize])
                            .map(|&c| landmass[c as usize]),
                    );
                }
            }
            lands.sort_unstable();
            lands.dedup();
            found.push(Basin {
                vertices,
                area,
                touches_coast,
                landmasses: lands.len(),
            });
        }

        found.sort_by(|a, b| {
            b.area
                .total_cmp(&a.area)
                .then(a.vertices[0].cmp(&b.vertices[0]))
        });
        let (basins, dropped): (Vec<Basin>, Vec<Basin>) = found
            .into_iter()
            .partition(|b| b.vertices.len() >= min_vertices);
        let mut vertex_basin = vec![None; n_vertices];
        for (k, basin) in basins.iter().enumerate() {
            for &v in &basin.vertices {
                vertex_basin[v as usize] = Some(k as u32);
            }
        }
        // A wet cell takes its interior vertices' basin, the largest if they
        // differ (the lowest number).
        let cell_basin = (0..n_cells)
            .map(|c| {
                if !wet[c] {
                    return None;
                }
                mesh.cell_vertices(c as CellId)
                    .iter()
                    .filter_map(|&v| vertex_basin[v as usize])
                    .min()
            })
            .collect();
        Self {
            basins,
            dropped,
            wet,
            vertex_basin,
            cell_basin,
        }
    }

    /// The kept basins, largest first; a basin's number is its index here.
    pub fn basins(&self) -> &[Basin] {
        &self.basins
    }

    /// The basins dropped as too small, largest first. Their vertices are
    /// coast.
    pub fn dropped(&self) -> &[Basin] {
        &self.dropped
    }

    /// Whether `cell` is wet: at least [`WET_OCEAN_FRACTION`] ocean.
    pub fn is_wet(&self, cell: CellId) -> bool {
        self.wet[cell as usize]
    }

    /// The basin of `vertex`, or `None` where it is coast (including the
    /// vertices of dropped basins).
    pub fn vertex_basin(&self, vertex: VertexId) -> Option<usize> {
        self.vertex_basin[vertex as usize].map(|b| b as usize)
    }

    /// The basin of `cell`, for display: for a wet cell the basin of its
    /// interior vertices, the largest if they differ; `None` for a cell that
    /// isn't wet or has no interior vertex in a kept basin.
    pub fn cell_basin(&self, cell: CellId) -> Option<usize> {
        self.cell_basin[cell as usize].map(|b| b as usize)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fms_grid::coords::{from_lat_lon, to_lat_lon};
    use fms_grid::{LloydOptions, MeshGeometry, SphereMesh};
    use glam::DVec3;

    const RADIUS: f64 = 6.371e6;

    struct Grid {
        mesh: SphereMesh,
        geometry: MeshGeometry,
    }

    impl Grid {
        fn new(k: u32) -> Self {
            let (mesh, _) = SphereMesh::centroidal(k, LloydOptions::default());
            let geometry = MeshGeometry::new(&mesh);
            Self { mesh, geometry }
        }

        fn ops(&self) -> Operators<'_> {
            Operators::new(&self.mesh, &self.geometry, RADIUS)
        }

        fn basins(&self, land: &[f64]) -> Basins {
            Basins::new(&self.ops(), land)
        }

        /// Land fraction 1 where `land` holds at the cell centre, else 0.
        fn land(&self, land: impl Fn(DVec3) -> bool) -> Vec<f64> {
            self.mesh
                .cell_centers()
                .iter()
                .map(|&p| if land(p) { 1.0 } else { 0.0 })
                .collect()
        }

        /// The basin of the cell nearest `p`.
        fn basin_at(&self, b: &Basins, p: DVec3) -> Option<usize> {
            b.cell_basin(self.nearest(p))
        }

        fn nearest(&self, p: DVec3) -> CellId {
            (0..self.mesh.n_cells())
                .max_by(|&a, &c| {
                    let centres = self.mesh.cell_centers();
                    centres[a].dot(p).total_cmp(&centres[c].dot(p))
                })
                .unwrap() as CellId
        }
    }

    fn angle(a: DVec3, b: DVec3) -> f64 {
        fms_math::atan2(a.cross(b).length(), a.dot(b))
    }

    fn assert_ordered(b: &Basins) {
        for list in [b.basins(), b.dropped()] {
            for pair in list.windows(2) {
                assert!(
                    pair[0].area > pair[1].area
                        || (pair[0].area == pair[1].area
                            && pair[0].vertices[0] < pair[1].vertices[0])
                );
            }
            for basin in list {
                assert!(basin.vertices.windows(2).all(|v| v[0] < v[1]));
            }
        }
    }

    #[test]
    fn a_planet_without_land_is_one_basin_without_coast() {
        let g = Grid::new(8);
        let b = g.basins(&vec![0.0; g.mesh.n_cells()]);
        assert_eq!(b.basins().len(), 1);
        assert!(b.dropped().is_empty());
        let ocean = &b.basins()[0];
        assert!(!ocean.touches_coast);
        assert_eq!(ocean.landmasses, 0);
        assert_eq!(ocean.vertices.len(), g.mesh.n_vertices());
        let sphere = 4.0 * std::f64::consts::PI * RADIUS * RADIUS;
        assert!((ocean.area / sphere - 1.0).abs() < 1e-12);
        assert!((0..g.mesh.n_cells()).all(|c| b.cell_basin(c as CellId) == Some(0)));
    }

    #[test]
    fn the_idealized_continent_borders_one_basin() {
        let g = Grid::new(16);
        let land = g.land(|p| {
            let (lat, lon) = to_lat_lon(p);
            (-35.0..=60.0).contains(&lat.to_degrees()) && (-30.0..=30.0).contains(&lon.to_degrees())
        });
        let b = g.basins(&land);
        assert_eq!(b.basins().len(), 1);
        let ocean = &b.basins()[0];
        assert!(ocean.touches_coast);
        assert_eq!(ocean.landmasses, 1);
        for c in 0..g.mesh.n_cells() as CellId {
            assert_eq!(b.is_wet(c), land[c as usize] == 0.0);
            if !b.is_wet(c) {
                assert_eq!(b.cell_basin(c), None);
            }
        }
        // Every vertex next to the continent is coast, and every other one
        // is in the ocean.
        for (v, cells) in g.mesh.triangles().iter().enumerate() {
            let coast = cells.iter().any(|&c| land[c as usize] > 0.0);
            assert_eq!(b.vertex_basin(v as VertexId).is_none(), coast);
        }
    }

    #[test]
    fn two_barriers_from_pole_to_pole_make_two_basins() {
        // Land along the meridians 0° and 180°, joined at the poles.
        let g = Grid::new(16);
        let b = g.basins(&g.land(|p| p.y.abs() < 0.08));
        assert_eq!(b.basins().len(), 2, "{:?}", b.dropped());
        assert_ordered(&b);
        let east = g.basin_at(&b, from_lat_lon(0.0, 1.5)).unwrap();
        let west = g.basin_at(&b, from_lat_lon(0.0, -1.5)).unwrap();
        assert_ne!(east, west);
        for basin in b.basins() {
            assert!(basin.touches_coast);
            assert_eq!(basin.landmasses, 1);
        }
    }

    #[test]
    fn an_enclosed_sea_is_its_own_basin() {
        // A ring of land around a sea.
        let g = Grid::new(16);
        let centre = from_lat_lon(0.5, 1.0);
        let b = g.basins(&g.land(|p| (0.3..0.6).contains(&angle(p, centre))));
        assert_eq!(b.basins().len(), 2);
        assert_ordered(&b);
        assert_eq!(g.basin_at(&b, centre), Some(1));
        assert_eq!(g.basin_at(&b, -centre), Some(0));
        assert!(b.basins().iter().all(|s| s.landmasses == 1));
    }

    #[test]
    fn a_channel_one_cell_wide_joins_nothing() {
        // Two seas in a world of land, then a channel one cell wide between
        // them: the shortest path of cells from one to the other.
        let g = Grid::new(16);
        let (a, z) = (from_lat_lon(0.0, -0.6), from_lat_lon(0.0, 0.6));
        let mut land = g.land(|p| angle(p, a) > 0.3 && angle(p, z) > 0.3);
        let seas = g.basins(&land);
        assert_eq!(seas.basins().len(), 2);

        let n = g.mesh.n_cells();
        let mut from = vec![u32::MAX; n];
        let start = g.nearest(a) as usize;
        from[start] = start as u32;
        let mut queue = std::collections::VecDeque::from([start]);
        let end = loop {
            let c = queue.pop_front().unwrap();
            if angle(g.mesh.cell_centers()[c], z) <= 0.3 {
                break c;
            }
            for &next in g.mesh.cell_neighbors(c as CellId) {
                if from[next as usize] == u32::MAX {
                    from[next as usize] = c as u32;
                    queue.push_back(next as usize);
                }
            }
        };
        let mut channel = Vec::new();
        let mut c = from[end] as usize;
        while land[c] == 1.0 {
            land[c] = 0.0;
            channel.push(c);
            c = from[c] as usize;
        }
        assert!(channel.len() >= 3, "channel of {}", channel.len());

        let b = g.basins(&land);
        assert_eq!(b.basins().len(), 2);
        assert_ne!(g.basin_at(&b, a), g.basin_at(&b, z));
        // No vertex has all three of its cells in the channel, so the
        // channel has no interior vertex of its own: its vertices are coast,
        // or (at its ends) in the sea they touch.
        let in_channel = |c: CellId| channel.contains(&(c as usize));
        for &c in &channel {
            for &v in g.mesh.cell_vertices(c as CellId) {
                let cells = g.mesh.triangles()[v as usize];
                assert!(!cells.iter().all(|&c| in_channel(c)));
                if let Some(k) = b.vertex_basin(v) {
                    let sea = cells.iter().find(|&&c| !in_channel(c)).unwrap();
                    assert_eq!(b.cell_basin(*sea), Some(k));
                }
            }
        }
    }

    #[test]
    fn a_tiny_sea_is_dropped() {
        // One cell and its neighbours: the six vertices of the middle cell
        // are interior, too few for a gyre.
        let g = Grid::new(16);
        let middle = g.nearest(from_lat_lon(0.4, 2.0));
        let mut land = g.land(|p| p.z > -0.3);
        land[middle as usize] = 0.0;
        for &c in g.mesh.cell_neighbors(middle) {
            land[c as usize] = 0.0;
        }
        let b = g.basins(&land);
        assert_eq!(b.basins().len(), 1);
        assert_eq!(b.dropped().len(), 1);
        let tiny = &b.dropped()[0];
        assert_eq!(tiny.vertices.len(), g.mesh.cell_vertices(middle).len());
        assert_eq!(tiny.landmasses, 1);
        assert!(tiny.vertices.iter().all(|&v| b.vertex_basin(v).is_none()));
        assert_eq!(b.cell_basin(middle), None);
        assert!(b.is_wet(middle));
        // Kept with a lower bar.
        let kept = Basins::with_min_vertices(&g.ops(), &land, 6);
        assert_eq!(kept.basins().len(), 2);
        assert_eq!(kept.cell_basin(middle), Some(1));
    }

    #[test]
    fn half_ocean_is_wet_and_less_is_not() {
        let g = Grid::new(8);
        let cell = g.nearest(from_lat_lon(0.2, 0.3));
        let vertices = g.mesh.cell_vertices(cell).to_vec();
        for (land, wet) in [
            (0.5, true),
            (0.5 - 1e-9, true),
            (0.5 + 1e-9, false),
            (1.0, false),
        ] {
            let mut fraction = vec![0.0; g.mesh.n_cells()];
            fraction[cell as usize] = land;
            let b = g.basins(&fraction);
            assert_eq!(b.is_wet(cell), wet, "land {land}");
            for &v in &vertices {
                assert_eq!(b.vertex_basin(v).is_some(), wet, "land {land}");
            }
            assert_eq!(b.basins()[0].touches_coast, !wet);
            assert_eq!(b.basins()[0].landmasses, usize::from(!wet));
        }
    }

    #[test]
    fn the_basins_are_the_same_on_a_mesh_built_twice() {
        let land = |g: &Grid| {
            g.land(|p| {
                p.y.abs() < 0.08
                    || angle(p, from_lat_lon(0.5, 1.0)) < 0.4
                    || angle(p, from_lat_lon(-0.6, -2.0)) < 0.2
            })
        };
        let (g1, g2) = (Grid::new(16), Grid::new(16));
        let (b1, b2) = (g1.basins(&land(&g1)), g2.basins(&land(&g2)));
        assert_eq!(b1, b2);
        assert!(b1.basins().len() >= 2);
        assert_ordered(&b1);
    }
}
