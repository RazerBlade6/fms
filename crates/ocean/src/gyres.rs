//! The wind-driven gyres: Stommel's steady balance for the depth-integrated
//! mass transport, solved on the mesh's vertices, basin by basin.
//!
//! The transport is `M = k × ∇Ψ` (kg m⁻¹ s⁻¹) for a streamfunction `Ψ`
//! (kg s⁻¹), so northward transport is `∂Ψ/∂x` and water goes clockwise
//! round a maximum of `Ψ`, seen from outside the sphere. With linear friction
//! `r` (s⁻¹) and the wind's stress `τ` (N m⁻²) it obeys
//!
//! `r ∇²Ψ + β ∂Ψ/∂x = k·∇×τ`,
//!
//! with `β = ∂f/∂y` the northward gradient of the Coriolis parameter
//! `f = 2Ω sin φ`. The `β` term is the divergence of `Ψ c`, with
//! `c = −k × ∇f` pointing east at `β` (m⁻¹ s⁻¹) and having no divergence, and
//! it is in that form, as a flux, that it is taken here.
//!
//! `Ψ` is zero on every coast vertex ([`Basins`]), so no water crosses a coast,
//! and the flow across a mesh edge is the difference of `Ψ` between its two
//! vertices ([`Gyres::edge_flux`]), so it has no divergence in any cell.
//!
//! **Discretisation.** An interior vertex's control volume is the triangle
//! joining its three cell centres (area `A_v R²`). Each side of it is the dual
//! edge of one of the vertex's three mesh edges, which leads to a neighbour
//! vertex `w`. Integrated over the triangle and multiplied by −1, the equation
//! at vertex `v` is
//!
//! `Σ_sides [ r (dₑ/lₑ)(Ψ_v − Ψ_w) + Gₑ Ψ_upwind ] = −A_v R² curl_v`   (kg s⁻²)
//!
//! - `r (dₑ/lₑ)(Ψ_v − Ψ_w)` is the finite-volume Laplacian: the gradient of
//!   `Ψ` along the mesh edge (length `lₑR`) times the side it crosses (length
//!   `dₑR`).
//! - `Gₑ` is the flux of `−c` out of the triangle across the side (s⁻¹). The
//!   flux of `k × ∇f` across a curve is the difference of `f` between its
//!   ends, so `Gₑ = t_{e,v} (f[cells[0]] − f[cells[1]])`, with `t_{e,v}` the
//!   sign of [`MeshGeometry::vertex_edge_signs`](fms_grid::MeshGeometry::vertex_edge_signs)
//!   and `f` at the cell centres. The three sum to zero, as `c` has no
//!   divergence.
//! - `Ψ` on the side is taken upwind of `−c`: `Ψ_v` where `Gₑ` leaves the
//!   triangle, `Ψ_w` where it enters. `−c` is the way Rossby waves carry `Ψ`,
//!   westward, and the boundary current forms where they arrive, at the
//!   western coast. So the matrix is an M-matrix (positive diagonal, no
//!   positive entry off it, each row summing to at least zero) and the
//!   solution is monotone: the western boundary current is held to about one
//!   vertex spacing without wiggles, though no mesh here resolves it. The
//!   price is a numerical friction of about `β Δx / 2`, with `Δx` the vertex
//!   spacing: 2.7 × 10⁻⁶ s⁻¹ at 30° on a mesh of 2,562 cells, as much as
//!   the default [`GyreParams::friction`], and half that at 10,242 cells.
//!   On the mesh's triangles, which point alternately one way and the
//!   other, that friction is first order only on average: at one vertex the
//!   upwind `β` term is off by a third or more, with opposite signs on
//!   neighbouring vertices. The solution still converges at first order.
//!
//! `f` takes the signed rotation rate, so on a planet turning backwards `c`
//! points west and the boundary currents are on the eastern coasts.
//!
//! **Basins.** Each basin with a coast ([`Basin::touches_coast`]) is its own
//! linear system, its unknowns its interior vertices in ascending order. A
//! basin without one, an aquaplanet's single ocean, is skipped and keeps
//! `Ψ = 0`: it would have no boundary condition, and a zonal wind's curl has
//! no Sverdrup balance round a periodic ocean. Every island is a coast with
//! `Ψ = 0` of its own, so no water goes round one: there is no circumpolar
//! current.
//!
//! **Solver.** The matrices depend only on the mesh, the land and the
//! rotation, so they are built once, with an incomplete LU factorisation
//! (ILU(0), on the matrix's own pattern) as a preconditioner. Each solve is
//! preconditioned BiCGSTAB from the last solution, to a relative residual
//! `‖b − AΨ‖/‖b‖` of [`TOLERANCE`]. Everything runs in serial, in index
//! order: the same input bits give the same output bits.

use fms_grid::{EdgeId, Operators, VertexId};

use crate::basins::{Basin, Basins};

/// Relative residual `‖b − AΨ‖/‖b‖` a solve stops at.
pub const TOLERANCE: f64 = 1e-8;

/// Most BiCGSTAB iterations of one basin's solve.
pub const MAX_ITERATIONS: usize = 1000;

/// The gyres' one parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GyreParams {
    /// Stommel's linear bottom friction `r` (s⁻¹). The western boundary
    /// current is `r/β` wide, 100 km at 30° with the default; the mesh widens
    /// it to a vertex spacing or two (see the [module](self) documentation).
    pub friction: f64,
}

impl Default for GyreParams {
    fn default() -> Self {
        Self { friction: 2e-6 }
    }
}

/// How one basin's last solve went.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolveStats {
    /// BiCGSTAB iterations; 0 when the starting `Ψ` was already good enough
    /// or the curl was zero.
    pub iterations: usize,
    /// The final relative residual `‖b − AΨ‖/‖b‖`, computed afresh, not from
    /// the iteration's recurrence; 0 for a zero curl.
    pub residual: f64,
    /// Whether the residual is under [`TOLERANCE`]. Otherwise the solve
    /// stopped at [`MAX_ITERATIONS`].
    pub converged: bool,
}

/// One side of a vertex's control volume, as the gyre equation sees it.
#[derive(Clone, Copy, Debug)]
struct Side {
    /// The vertex at the other end of the mesh edge the side crosses.
    neighbour: VertexId,
    /// `dₑ/lₑ`: the side's length over the mesh edge's.
    conductance: f64,
    /// The flux of `c = −k × ∇f` out of the vertex's triangle across the side
    /// (s⁻¹): `t_{e,v} (f[cells[1]] − f[cells[0]])`.
    flux: f64,
}

/// The three sides of `v`'s control volume, given `f` at the cell centres.
fn sides(ops: &Operators, coriolis: &[f64], v: VertexId) -> [Side; 3] {
    let (mesh, geometry) = (ops.mesh(), ops.geometry());
    let edges = mesh.vertex_edges(v);
    let signs = geometry.vertex_edge_signs(v);
    std::array::from_fn(|k| {
        let e = edges[k] as usize;
        let edge = mesh.edges()[e];
        let [a, b] = edge.cells;
        Side {
            neighbour: if edge.vertices[0] == v {
                edge.vertices[1]
            } else {
                edge.vertices[0]
            },
            conductance: geometry.dual_edge_lengths[e] / geometry.edge_lengths[e],
            flux: signs[k] * (coriolis[b as usize] - coriolis[a as usize]),
        }
    })
}

/// A sparse matrix in compressed rows, each row's columns ascending.
#[derive(Clone, Debug)]
struct Csr {
    offsets: Vec<u32>,
    columns: Vec<u32>,
    values: Vec<f64>,
    /// Where each row's diagonal entry is.
    diagonal: Vec<u32>,
}

impl Csr {
    fn rows(&self) -> usize {
        self.diagonal.len()
    }

    fn row(&self, i: usize) -> std::ops::Range<usize> {
        self.offsets[i] as usize..self.offsets[i + 1] as usize
    }

    /// `out = A x`.
    fn apply(&self, x: &[f64], out: &mut [f64]) {
        for (i, out) in out.iter_mut().enumerate() {
            *out = self
                .row(i)
                .map(|k| self.values[k] * x[self.columns[k] as usize])
                .sum();
        }
    }

    /// `out = b − A x`.
    fn residual(&self, b: &[f64], x: &[f64], out: &mut [f64]) {
        self.apply(x, out);
        for (out, b) in out.iter_mut().zip(b) {
            *out = b - *out;
        }
    }

    /// The incomplete LU factorisation on the matrix's own pattern: `L` (unit
    /// diagonal, not stored) below the diagonal and `U` on and above it.
    fn ilu0(&self) -> Vec<f64> {
        let mut lu = self.values.clone();
        for i in 0..self.rows() {
            let row = self.row(i);
            for k in row.start..self.diagonal[i] as usize {
                let pivot = self.columns[k] as usize;
                lu[k] /= lu[self.diagonal[pivot] as usize];
                let l = lu[k];
                // Subtract l times the pivot row's U from this row, where the
                // pattern has room.
                let mut j = k + 1;
                for p in self.diagonal[pivot] as usize + 1..self.row(pivot).end {
                    let column = self.columns[p];
                    while j < row.end && self.columns[j] < column {
                        j += 1;
                    }
                    if j < row.end && self.columns[j] == column {
                        lu[j] -= l * lu[p];
                    }
                }
            }
        }
        lu
    }

    /// `out = (LU)⁻¹ r` for the factors `lu` from [`ilu0`](Self::ilu0).
    fn precondition(&self, lu: &[f64], r: &[f64], out: &mut [f64]) {
        for i in 0..self.rows() {
            let below = self.offsets[i] as usize..self.diagonal[i] as usize;
            out[i] = r[i]
                - below
                    .map(|k| lu[k] * out[self.columns[k] as usize])
                    .sum::<f64>();
        }
        for i in (0..self.rows()).rev() {
            let d = self.diagonal[i] as usize;
            let above = d + 1..self.row(i).end;
            let sum: f64 = above.map(|k| lu[k] * out[self.columns[k] as usize]).sum();
            out[i] = (out[i] - sum) / lu[d];
        }
    }
}

/// One basin's linear system.
#[derive(Clone, Debug)]
struct System {
    /// Its interior vertices, ascending: unknown `i` is `Ψ` at `vertices[i]`.
    vertices: Vec<VertexId>,
    matrix: Csr,
    /// The ILU(0) factors, on `matrix`'s pattern.
    factors: Vec<f64>,
    /// `−A_v R²` of each unknown's vertex (m²), turning the curl into the
    /// right-hand side.
    scale: Vec<f64>,
}

impl System {
    fn new(ops: &Operators, coriolis: &[f64], basin: &Basin, friction: f64) -> Self {
        let n_vertices = ops.mesh().n_vertices();
        let mut unknown = vec![u32::MAX; n_vertices];
        for (i, &v) in basin.vertices.iter().enumerate() {
            unknown[v as usize] = i as u32;
        }
        let r2 = ops.radius() * ops.radius();
        let n = basin.vertices.len();
        let mut offsets = Vec::with_capacity(n + 1);
        offsets.push(0);
        let (mut columns, mut values) = (Vec::with_capacity(4 * n), Vec::with_capacity(4 * n));
        let mut diagonal = Vec::with_capacity(n);
        let mut entries: Vec<(u32, f64)> = Vec::with_capacity(4);
        for (i, &v) in basin.vertices.iter().enumerate() {
            entries.clear();
            let mut centre = 0.0;
            for side in sides(ops, coriolis, v) {
                let diffusion = friction * side.conductance;
                let outflow = -side.flux;
                centre += diffusion + outflow.max(0.0);
                // A neighbour that isn't an unknown is coast, where Ψ = 0.
                let w = unknown[side.neighbour as usize];
                if w != u32::MAX {
                    entries.push((w, -diffusion + outflow.min(0.0)));
                }
            }
            entries.push((i as u32, centre));
            entries.sort_unstable_by_key(|&(column, _)| column);
            for &(column, value) in &entries {
                if column == i as u32 {
                    diagonal.push(columns.len() as u32);
                }
                columns.push(column);
                values.push(value);
            }
            offsets.push(columns.len() as u32);
        }
        let matrix = Csr {
            offsets,
            columns,
            values,
            diagonal,
        };
        let factors = matrix.ilu0();
        let scale = basin
            .vertices
            .iter()
            .map(|&v| -ops.geometry().vertex_areas[v as usize] * r2)
            .collect();
        Self {
            vertices: basin.vertices.clone(),
            matrix,
            factors,
            scale,
        }
    }

    /// Solve for `Ψ` with the curl `curl` (N m⁻³, per vertex), starting from
    /// and writing to `psi` (kg s⁻¹, per vertex).
    fn solve(&self, curl: &[f64], psi: &mut [f64]) -> SolveStats {
        let b: Vec<f64> = self
            .vertices
            .iter()
            .zip(&self.scale)
            .map(|(&v, s)| s * curl[v as usize])
            .collect();
        let mut x: Vec<f64> = self.vertices.iter().map(|&v| psi[v as usize]).collect();
        let stats = bicgstab(&self.matrix, &self.factors, &b, &mut x);
        for (&v, x) in self.vertices.iter().zip(x) {
            psi[v as usize] = x;
        }
        stats
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

/// Solve `A x = b` by BiCGSTAB (van der Vorst 1992) with the ILU(0) factors
/// `lu` as a right preconditioner, so the residual it tracks is the true one.
/// Starts from `x`; restarts from the current `x` if the recurrence breaks
/// down or its residual drifts from the true one.
fn bicgstab(a: &Csr, lu: &[f64], b: &[f64], x: &mut [f64]) -> SolveStats {
    let n = b.len();
    let norm_b = norm(b);
    if norm_b == 0.0 {
        x.fill(0.0);
        return SolveStats {
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
    }
    let target = TOLERANCE * norm_b;
    let mut r = vec![0.0; n];
    a.residual(b, x, &mut r);
    let mut residual = norm(&r);
    let [mut p, mut v, mut s, mut t, mut p_hat, mut s_hat] = std::array::from_fn(|_| vec![0.0; n]);
    let mut iterations = 0;
    while residual > target && iterations < MAX_ITERATIONS {
        let shadow = r.clone();
        let (mut rho_old, mut alpha, mut omega) = (1.0, 1.0, 1.0);
        p.fill(0.0);
        v.fill(0.0);
        while iterations < MAX_ITERATIONS {
            iterations += 1;
            let rho = dot(&shadow, &r);
            if rho == 0.0 || !rho.is_finite() {
                break;
            }
            let beta = (rho / rho_old) * (alpha / omega);
            for k in 0..n {
                p[k] = r[k] + beta * (p[k] - omega * v[k]);
            }
            a.precondition(lu, &p, &mut p_hat);
            a.apply(&p_hat, &mut v);
            let shadow_v = dot(&shadow, &v);
            if shadow_v == 0.0 {
                break;
            }
            alpha = rho / shadow_v;
            for k in 0..n {
                s[k] = r[k] - alpha * v[k];
            }
            if norm(&s) <= target {
                for k in 0..n {
                    x[k] += alpha * p_hat[k];
                }
                break;
            }
            a.precondition(lu, &s, &mut s_hat);
            a.apply(&s_hat, &mut t);
            let tt = dot(&t, &t);
            omega = if tt > 0.0 { dot(&t, &s) / tt } else { 0.0 };
            for k in 0..n {
                x[k] += alpha * p_hat[k] + omega * s_hat[k];
                r[k] = s[k] - omega * t[k];
            }
            rho_old = rho;
            if omega == 0.0 || norm(&r) <= target {
                break;
            }
        }
        a.residual(b, x, &mut r);
        residual = norm(&r);
    }
    SolveStats {
        iterations,
        residual: residual / norm_b,
        converged: residual <= target,
    }
}

/// The gyres' streamfunction on one mesh, land and rotation.
#[derive(Clone, Debug)]
pub struct Gyres {
    params: GyreParams,
    /// One per basin, in [`Basins::basins`] order; `None` for a basin without
    /// a coast, which isn't solved.
    systems: Vec<Option<System>>,
    stats: Vec<Option<SolveStats>>,
    /// Per vertex (kg s⁻¹).
    psi: Vec<f64>,
    /// Each mesh edge's two vertices.
    ends: Vec<[VertexId; 2]>,
}

impl Gyres {
    /// The gyres of the `basins` on the mesh of `ops`, on a planet turning at
    /// `rotation` (rad s⁻¹, negative for retrograde), with `Ψ = 0` until the
    /// first [`solve`](Self::solve). Builds and factorises each basin's
    /// matrix.
    ///
    /// # Panics
    /// If `params.friction` isn't positive, or `basins` isn't of this mesh.
    pub fn new(ops: &Operators, basins: &Basins, rotation: f64, params: GyreParams) -> Self {
        assert!(params.friction > 0.0, "the gyres need friction");
        let mesh = ops.mesh();
        let coriolis: Vec<f64> = mesh
            .cell_centers()
            .iter()
            .map(|p| 2.0 * rotation * p.z)
            .collect();
        let systems: Vec<Option<System>> = basins
            .basins()
            .iter()
            .map(|basin| {
                basin
                    .touches_coast
                    .then(|| System::new(ops, &coriolis, basin, params.friction))
            })
            .collect();
        let stats = systems
            .iter()
            .map(|s| {
                s.as_ref().map(|_| SolveStats {
                    iterations: 0,
                    residual: 0.0,
                    converged: true,
                })
            })
            .collect();
        Self {
            params,
            systems,
            stats,
            psi: vec![0.0; mesh.n_vertices()],
            ends: mesh.edges().iter().map(|e| e.vertices).collect(),
        }
    }

    /// The parameters.
    pub fn params(&self) -> &GyreParams {
        &self.params
    }

    /// Solve every basin with a coast for the stress curl `curl` (N m⁻³, one
    /// per vertex), each from its last solution. Values at vertices that
    /// aren't in a solved basin are ignored.
    ///
    /// # Panics
    /// If `curl` doesn't have one value per vertex.
    pub fn solve(&mut self, curl: &[f64]) {
        assert_eq!(curl.len(), self.psi.len(), "one curl per vertex");
        for (system, stats) in self.systems.iter().zip(&mut self.stats) {
            if let Some(system) = system {
                *stats = Some(system.solve(curl, &mut self.psi));
            }
        }
    }

    /// The streamfunction at each vertex (kg s⁻¹): zero on coasts and in
    /// basins that aren't solved.
    pub fn streamfunction(&self) -> &[f64] {
        &self.psi
    }

    /// The mass flowing across `edge` from its `cells[0]` to its `cells[1]`
    /// (kg s⁻¹): `Ψ[vertices[0]] − Ψ[vertices[1]]`.
    ///
    /// The edge's tangent `tₑ = k × nₑ` runs from `vertices[0]` to
    /// `vertices[1]`, so the transport along the normal is
    /// `(k × ∇Ψ)·nₑ = −∇Ψ·tₑ`, and times the edge's length that is the
    /// difference. Summed round a cell, these cancel: the flow has no
    /// divergence. Across an edge with a dry cell both vertices are coast, so
    /// it is exactly zero.
    pub fn edge_flux(&self, edge: EdgeId) -> f64 {
        let [a, b] = self.ends[edge as usize];
        self.psi[a as usize] - self.psi[b as usize]
    }

    /// How each basin's last solve went, in [`Basins::basins`] order; `None`
    /// for a basin without a coast, which isn't solved.
    pub fn stats(&self) -> &[Option<SolveStats>] {
        &self.stats
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::EARTH_ROTATION;
    use fms_grid::coords::to_lat_lon;
    use fms_grid::{CellId, LloydOptions, MeshGeometry, SphereMesh};
    use glam::DVec3;
    use std::f64::consts::PI;

    const RADIUS: f64 = 6.371e6;
    const OMEGA: f64 = EARTH_ROTATION;
    /// Peak stress of the test wind (N m⁻²).
    const TAU0: f64 = 0.1;

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

        /// Latitude and longitude of each vertex (radians).
        fn vertex_lat_lon(&self) -> Vec<(f64, f64)> {
            self.mesh
                .vertices()
                .iter()
                .map(|&p| to_lat_lon(p))
                .collect()
        }

        /// The mean length of the mesh's edges (m): the vertex spacing.
        fn spacing(&self) -> f64 {
            let l = &self.geometry.edge_lengths;
            l.iter().sum::<f64>() / l.len() as f64 * RADIUS
        }
    }

    /// A box of ocean in a world of land, in degrees.
    #[derive(Clone, Copy)]
    struct Rectangle {
        south: f64,
        north: f64,
        west: f64,
        east: f64,
    }

    /// The test basin, about the size of the North Atlantic's subtropics.
    /// Its sides are a degree east of the meridians 0° and 60°, on which
    /// cell centres lie, so round-off can't put them on either side.
    const NORTH: Rectangle = Rectangle {
        south: 10.0,
        north: 50.0,
        west: 1.0,
        east: 61.0,
    };

    /// Its mirror image south of the equator.
    const SOUTH: Rectangle = Rectangle {
        south: -50.0,
        north: -10.0,
        ..NORTH
    };

    impl Rectangle {
        /// Land fraction 0 at cell centres in the box, else 1.
        fn land(&self, g: &Grid) -> Vec<f64> {
            g.mesh
                .cell_centers()
                .iter()
                .map(|&p| {
                    let (lat, lon) = to_lat_lon(p);
                    let inside = (self.south..=self.north).contains(&lat.to_degrees())
                        && (self.west..=self.east).contains(&lon.to_degrees());
                    if inside { 0.0 } else { 1.0 }
                })
                .collect()
        }

        /// The test wind's eastward stress `τ_x(φ) = −τ₀ cos(π(φ−φ_S)/(φ_N−φ_S))`
        /// and its derivative in latitude (N m⁻² rad⁻¹), at latitude `lat`
        /// (radians), for the box north of the equator; for the one south of
        /// it, the mirror image.
        fn stress(&self, lat: f64) -> (f64, f64) {
            let (south, north) = (NORTH.south.to_radians(), NORTH.north.to_radians());
            let mirrored = self.north < 0.0;
            let phi = if mirrored { -lat } else { lat };
            let width = north - south;
            let theta = PI * (phi - south) / width;
            let tau = -TAU0 * fms_math::cos(theta);
            let slope = TAU0 * fms_math::sin(theta) * PI / width;
            (tau, if mirrored { -slope } else { slope })
        }

        /// The wind's curl at each vertex (N m⁻³):
        /// `−(1/(R cos φ)) ∂(τ_x cos φ)/∂φ = (τ_x tan φ − ∂τ_x/∂φ)/R`.
        fn curl(&self, g: &Grid) -> Vec<f64> {
            g.vertex_lat_lon()
                .iter()
                .map(|&(lat, _)| {
                    let (tau, slope) = self.stress(lat);
                    (tau * fms_math::tan(lat) - slope) / RADIUS
                })
                .collect()
        }

        /// Whether `(lat, lon)` (radians) is at least `margin` degrees inside
        /// the box from the north and south and in the given longitudes.
        fn holds(&self, (lat, lon): (f64, f64), margin: f64, lons: (f64, f64)) -> bool {
            let (lat, lon) = (lat.to_degrees(), lon.to_degrees());
            (self.south + margin..=self.north - margin).contains(&lat)
                && (lons.0..=lons.1).contains(&lon)
        }
    }

    /// The solved gyres of a box.
    fn solve(g: &Grid, b: Rectangle, rotation: f64) -> (Basins, Gyres) {
        let basins = Basins::new(&g.ops(), &b.land(g));
        assert_eq!(basins.basins().len(), 1);
        let mut gyres = Gyres::new(&g.ops(), &basins, rotation, GyreParams::default());
        gyres.solve(&b.curl(g));
        (basins, gyres)
    }

    /// The vertex where `|Ψ|` is largest.
    fn extremum(gyres: &Gyres) -> usize {
        let psi = gyres.streamfunction();
        (0..psi.len())
            .max_by(|&a, &b| psi[a].abs().total_cmp(&psi[b].abs()))
            .unwrap()
    }

    /// Sverdrup's interior: `Ψ_S = −(1/β) ∫ curl dx` from the vertex east to
    /// the coast at `b.east`, which for a zonal wind is
    /// `−curl R² (λ_E − λ) / (2Ω)` (kg s⁻¹).
    fn sverdrup(b: Rectangle, curl: f64, lon: f64) -> f64 {
        -curl * RADIUS * RADIUS * (b.east.to_radians() - lon) / (2.0 * OMEGA)
    }

    /// The largest `|a − b|` over the vertices where `which` holds, and the
    /// largest `|b|` there.
    fn worst(a: &[f64], b: &[f64], which: impl Fn(usize) -> bool) -> (f64, f64) {
        (0..a.len())
            .filter(|&v| which(v))
            .fold((0.0, 0.0), |(e, m): (f64, f64), v| {
                (e.max((a[v] - b[v]).abs()), m.max(b[v].abs()))
            })
    }

    /// Northward transport per metre (kg m⁻¹ s⁻¹) across each edge facing
    /// within 30° of north or south, with the edge point's longitude
    /// (degrees); edges facing more east or west are left out.
    fn northward(g: &Grid, gyres: &Gyres) -> Vec<(f64, f64)> {
        (0..g.mesh.n_edges())
            .filter_map(|e| {
                let point = g.geometry.edge_points[e];
                let north = DVec3::Z - point * point.z;
                let facing = g.geometry.edge_normals[e].dot(north.normalize());
                (facing.abs() > 0.866).then(|| {
                    let width = g.geometry.edge_lengths[e] * RADIUS;
                    let transport = gyres.edge_flux(e as EdgeId) * facing.signum() / width;
                    (transport, to_lat_lon(point).1.to_degrees())
                })
            })
            .collect()
    }

    #[test]
    fn the_boundary_current_is_on_the_western_coast() {
        let g = Grid::new(24);
        let (_, gyres) = solve(&g, NORTH, OMEGA);
        let psi = gyres.streamfunction();
        let ll = g.vertex_lat_lon();
        let m = extremum(&gyres);
        // Anticyclonic: Ψ has a maximum, round which water goes clockwise.
        // The wind's curl is positive within 4° of the northern coast, where
        // a weak cyclonic gyre is all that has the other sign.
        let lowest = psi.iter().fold(0.0, |a: f64, &p| a.min(p));
        assert!(
            psi[m] > 0.0 && lowest > -0.02 * psi[m],
            "{} {lowest}",
            psi[m]
        );
        assert!(ll[m].1.to_degrees() < NORTH.west + 15.0, "{:?}", ll[m]);
        // The strongest northward flow is within a vertex spacing or two of
        // the western coast, and the strongest southward flow is weaker
        // than a fifth of it: the return flow is spread across the basin.
        let flow = northward(&g, &gyres);
        let (most, lon) = flow
            .iter()
            .fold((0.0, 0.0), |a, &b| if b.0 > a.0 { b } else { a });
        let least = flow.iter().fold(0.0, |a: f64, b| a.min(b.0));
        assert!(most > 0.0);
        let spacing = g.spacing() / (RADIUS * fms_math::cos(30f64.to_radians())) * 180.0 / PI;
        assert!(
            lon < NORTH.west + 2.0 * spacing,
            "{lon}° with spacing {spacing}°"
        );
        assert!(-least < 0.2 * most, "{least} {most}");
    }

    #[test]
    fn the_interior_is_in_sverdrup_balance() {
        // The largest error in the eastern two-thirds, 8° from the northern
        // and southern coasts and 10° from the eastern one, against the
        // largest Ψ_S there. With the default friction it is mostly
        // friction's, which smooths Ψ across latitudes: 0.21, 0.17, 0.15 and
        // 0.13 at 2,562, 5,762, 10,242 and 40,962 cells. With a tenth of the
        // friction the scheme's own error shows, and is first order: 0.12,
        // 0.079, 0.065 and 0.040.
        let error = |g: &Grid, friction: f64| {
            let basins = Basins::new(&g.ops(), &NORTH.land(g));
            let mut gyres = Gyres::new(&g.ops(), &basins, OMEGA, GyreParams { friction });
            let curl = NORTH.curl(g);
            gyres.solve(&curl);
            let ll = g.vertex_lat_lon();
            let expected: Vec<f64> = (0..ll.len())
                .map(|v| sverdrup(NORTH, curl[v], ll[v].1))
                .collect();
            let (e, m) = worst(gyres.streamfunction(), &expected, |v| {
                NORTH.holds(ll[v], 8.0, (NORTH.west + 20.0, NORTH.east - 10.0))
            });
            assert!(m > 10e9, "{m}");
            e / m
        };
        let (coarse, fine) = (Grid::new(16), Grid::new(32));
        let medium = Grid::new(24);
        let friction = GyreParams::default().friction;
        assert!(error(&coarse, friction) < 0.25);
        assert!(error(&medium, friction) < 0.2);
        let (coarse, fine) = (
            error(&coarse, friction / 10.0),
            error(&fine, friction / 10.0),
        );
        assert!(
            coarse < 0.14 && fine < 0.09 && fine < 0.7 * coarse,
            "{coarse} {fine}"
        );
    }

    #[test]
    fn a_retrograde_planet_has_eastern_boundary_currents() {
        let g = Grid::new(24);
        let (_, gyres) = solve(&g, NORTH, -OMEGA);
        let m = extremum(&gyres);
        let (_, lon) = to_lat_lon(g.mesh.vertices()[m]);
        assert!(lon.to_degrees() > NORTH.east - 15.0, "{}", lon.to_degrees());
        assert!(gyres.streamfunction()[m] > 0.0);
    }

    #[test]
    fn the_southern_basin_is_the_mirror_image() {
        let g = Grid::new(24);
        let (north, south) = (solve(&g, NORTH, OMEGA).1, solve(&g, SOUTH, OMEGA).1);
        let vertices = g.mesh.vertices();
        let mirror = |v: usize| {
            let p = vertices[v] * DVec3::new(1.0, 1.0, -1.0);
            (0..vertices.len())
                .max_by(|&a, &b| vertices[a].dot(p).total_cmp(&vertices[b].dot(p)))
                .unwrap()
        };
        let (psi, mirrored) = (north.streamfunction(), south.streamfunction());
        let (mut error, mut largest, mut offset) = (0.0_f64, 0.0_f64, 0.0_f64);
        for v in (0..psi.len()).filter(|&v| psi[v] != 0.0) {
            let w = mirror(v);
            // The icosahedron is symmetric about the equator, and so is the
            // mesh, to round-off.
            offset = offset.max((vertices[w] * DVec3::new(1.0, 1.0, -1.0) - vertices[v]).length());
            error = error.max((mirrored[w] + psi[v]).abs());
            largest = largest.max(psi[v].abs());
        }
        assert!(offset < 1e-12, "{offset}");
        // So the gyres differ by no more than the solver's tolerance allows
        // (8 × 10⁻¹⁰ of the largest, on 5,762 cells).
        assert!(error < 1e-8 * largest, "{error} of {largest}");
    }

    #[test]
    fn the_flow_has_no_divergence_and_stays_off_land() {
        let g = Grid::new(24);
        let (basins, gyres) = solve(&g, NORTH, OMEGA);
        let largest = (0..g.mesh.n_edges())
            .map(|e| gyres.edge_flux(e as EdgeId).abs())
            .fold(0.0, f64::max);
        assert!(largest > 1e9);
        for c in 0..g.mesh.n_cells() as CellId {
            let out: f64 = g
                .mesh
                .cell_edges(c)
                .iter()
                .zip(g.geometry.cell_edge_signs(c))
                .map(|(&e, s)| s * gyres.edge_flux(e))
                .sum();
            assert!(out.abs() < 1e-14 * largest, "cell {c}: {out}");
        }
        for (e, edge) in g.mesh.edges().iter().enumerate() {
            if !edge.cells.iter().all(|&c| basins.is_wet(c)) {
                assert_eq!(gyres.edge_flux(e as EdgeId), 0.0);
            }
        }
    }

    #[test]
    fn the_edge_flux_is_the_flow_of_k_cross_grad_psi() {
        // Ψ = Ψ₀ sin φ: ∇Ψ = (Ψ₀ cos φ / R) north, so M = k × ∇Ψ is
        // −(Ψ₀ cos φ / R) east, a westward flow, and across an edge of
        // length lR it carries −Ψ₀ cos φ (east · nₑ) l.
        let g = Grid::new(16);
        let basins = Basins::new(&g.ops(), &vec![0.0; g.mesh.n_cells()]);
        let mut gyres = Gyres::new(&g.ops(), &basins, OMEGA, GyreParams::default());
        let psi0 = 1e10;
        gyres.psi = g.mesh.vertices().iter().map(|p| psi0 * p.z).collect();
        let (mut error, mut largest) = (0.0_f64, 0.0_f64);
        for e in 0..g.mesh.n_edges() {
            let point = g.geometry.edge_points[e];
            let east = DVec3::Z.cross(point);
            // |east| is cos φ.
            let expected =
                -psi0 * east.dot(g.geometry.edge_normals[e]) * g.geometry.edge_lengths[e];
            error = error.max((gyres.edge_flux(e as EdgeId) - expected).abs());
            largest = largest.max(expected.abs());
        }
        assert!(error < 0.01 * largest, "{error} of {largest}");
    }

    #[test]
    fn zero_wind_gives_no_gyre_and_double_the_wind_twice_the_gyre() {
        let g = Grid::new(16);
        let basins = Basins::new(&g.ops(), &NORTH.land(&g));
        let mut gyres = Gyres::new(&g.ops(), &basins, OMEGA, GyreParams::default());
        let curl = NORTH.curl(&g);
        gyres.solve(&curl);
        let once = gyres.streamfunction().to_vec();
        gyres.solve(&vec![0.0; curl.len()]);
        assert!(gyres.streamfunction().iter().all(|&p| p == 0.0));
        let stats = gyres.stats()[0].unwrap();
        assert_eq!((stats.iterations, stats.residual), (0, 0.0));
        // From Ψ = 0, every step of the solve scales by 2 exactly.
        let doubled: Vec<f64> = curl.iter().map(|c| 2.0 * c).collect();
        gyres.solve(&doubled);
        let twice = gyres.streamfunction();
        assert!(once.iter().any(|&p| p != 0.0));
        assert!(once.iter().zip(twice).all(|(a, b)| 2.0 * a == *b));
    }

    #[test]
    fn the_operator_is_consistent_and_monotone() {
        let g = Grid::new(24);
        let ops = g.ops();
        let coriolis: Vec<f64> = g
            .mesh
            .cell_centers()
            .iter()
            .map(|p| 2.0 * OMEGA * p.z)
            .collect();
        // c has no divergence: the three fluxes of each vertex cancel.
        for v in 0..g.mesh.n_vertices() as VertexId {
            let s = sides(&ops, &coriolis, v);
            let total: f64 = s.iter().map(|s| s.flux).sum();
            let size: f64 = s.iter().map(|s| s.flux.abs()).sum();
            assert!(
                total.abs() <= 1e-14 * size.max(1e-20),
                "{v}: {total} of {size}"
            );
        }
        // For Ψ = cos φ sin λ (the vertex's y), β ∂Ψ/∂x = 2Ω cos φ cos λ / R²
        // (its x times 2Ω/R²).
        let vertices = g.mesh.vertices();
        let r2 = RADIUS * RADIUS;
        let areas = &g.geometry.vertex_areas;
        let expected: Vec<f64> = vertices.iter().map(|p| 2.0 * OMEGA * p.x / r2).collect();
        let advect = |upwind: bool| -> Vec<f64> {
            (0..vertices.len())
                .map(|v| {
                    let own = vertices[v].y;
                    let sum: f64 = sides(&ops, &coriolis, v as VertexId)
                        .iter()
                        .map(|s| {
                            let other = vertices[s.neighbour as usize].y;
                            let face = match upwind {
                                false => 0.5 * (own + other),
                                // Upwind of −c: its own where −c leaves.
                                true if s.flux < 0.0 => own,
                                true => other,
                            };
                            s.flux * face
                        })
                        .sum();
                    sum / (areas[v] * r2)
                })
                .collect()
        };
        // Centred, the flux form gets it to 1.4% RMS and 7.8% at worst, next
        // to the pentagons.
        let centred = advect(false);
        let (e, m) = worst(&centred, &expected, |_| true);
        let squares = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
        let difference: Vec<f64> = centred.iter().zip(&expected).map(|(a, b)| a - b).collect();
        let rms = (squares(&difference, &difference) / squares(&expected, &expected)).sqrt();
        assert!(e < 0.1 * m && rms < 0.02, "centred {} {rms}", e / m);
        // Upwind adds a diffusion of `|c·n|`-weighted differences. On the
        // mesh's triangles, which point alternately up and down, that is
        // first order only on average: at a single vertex the error is a
        // third of the answer or more, of one sign on the triangles pointing
        // one way and the other on the rest. Against a smooth weight it
        // cancels to 3 × 10⁻⁴.
        let upwind = advect(true);
        let (e, m) = worst(&upwind, &expected, |_| true);
        assert!(e > 0.3 * m);
        let weak = |a: &[f64]| {
            (0..a.len())
                .map(|v| areas[v] * expected[v] * a[v])
                .sum::<f64>()
        };
        let error = (weak(&upwind) - weak(&expected)) / weak(&expected);
        assert!(error.abs() < 1e-3, "upwind {error}");
        // An M-matrix: positive diagonal, nothing positive off it, every row
        // summing to at least zero, and more than zero next to the coast.
        let (_, gyres) = solve(&g, NORTH, OMEGA);
        let a = &gyres.systems[0].as_ref().unwrap().matrix;
        let mut strict = 0;
        for i in 0..a.rows() {
            let d = a.values[a.diagonal[i] as usize];
            let row = a.row(i);
            let off: f64 = row.clone().map(|k| a.values[k]).sum::<f64>() - d;
            assert!(d > 0.0);
            assert!(
                row.clone()
                    .all(|k| k == a.diagonal[i] as usize || a.values[k] <= 0.0)
            );
            assert!(d + off >= -1e-12 * d, "row {i}: {d} {off}");
            if d + off > 1e-6 * d {
                strict += 1;
            }
        }
        assert!(strict > 0 && strict < a.rows());
    }

    #[test]
    fn an_ocean_without_coast_has_no_gyre() {
        let g = Grid::new(16);
        let basins = Basins::new(&g.ops(), &vec![0.0; g.mesh.n_cells()]);
        assert!(!basins.basins()[0].touches_coast);
        let mut gyres = Gyres::new(&g.ops(), &basins, OMEGA, GyreParams::default());
        gyres.solve(&NORTH.curl(&g));
        assert_eq!(gyres.stats(), &[None]);
        assert!(gyres.streamfunction().iter().all(|&p| p == 0.0));
    }

    #[test]
    fn solves_converge_and_are_the_same_on_a_mesh_built_twice() {
        let (g1, g2) = (Grid::new(24), Grid::new(24));
        let (_, mut gyres) = solve(&g1, NORTH, OMEGA);
        assert_eq!(
            gyres.streamfunction(),
            solve(&g2, NORTH, OMEGA).1.streamfunction()
        );
        let stats = gyres.stats()[0].unwrap();
        assert!(stats.converged && stats.residual < TOLERANCE, "{stats:?}");
        assert!(stats.iterations > 1);
        // From the last solution, the same wind needs at most one more step.
        let before = gyres.streamfunction().to_vec();
        gyres.solve(&NORTH.curl(&g1));
        let again = gyres.stats()[0].unwrap();
        assert!(
            again.iterations <= 1 && again.residual < TOLERANCE,
            "{again:?}"
        );
        let (e, m) = worst(gyres.streamfunction(), &before, |_| true);
        assert!(e <= 1e-6 * m);
    }
}
