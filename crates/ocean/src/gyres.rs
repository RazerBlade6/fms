//! The wind-driven gyres: Stommel's steady balance for the depth-integrated
//! mass transport, solved on the mesh's vertices, with each landmass's coast
//! at one value of the streamfunction found by the island rule.
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
//! The flow across a mesh edge is the difference of `Ψ` between its two
//! vertices ([`Gyres::edge_flux`]), so it has no divergence in any cell.
//! `Ψ` is one constant along each coast group
//! ([`CoastGroup`](crate::CoastGroup)): every vertex touching a landmass,
//! with the dropped basins joining landmasses. So no water crosses a coast,
//! but it can flow between two landmasses, `Ψ_A − Ψ_B` of it, round an
//! island or through a strait.
//!
//! **Discretisation.** An interior vertex's control volume is the triangle
//! joining its three cell centres (area `A_v R²`). Each side of it is the dual
//! edge of one of the vertex's three mesh edges, which leads to a neighbour
//! vertex `w`. Integrated over the triangle and multiplied by −1, the equation
//! at vertex `v` is
//!
//! `Σ_sides [ r (dₑ/lₑ)(Ψ_v − Ψ_w) + Gₑ Ψ_upwind ] = −Σ_sides t_{e,v} dₑ R τₑ`   (kg s⁻²)
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
//! - The right-hand side is the circulation of the stress round the
//!   triangle (times −1), as [`Operators::curl`] takes it: `τₑ` is the
//!   stress along the edge's normal, the mean of its two cells'
//!   ([`edge_stress`]), and `A_v R² curl_v = Σ t_{e,v} dₑ R τₑ`.
//!
//! `f` takes the signed rotation rate, so on a planet turning backwards `c`
//! points west and the boundary currents are on the eastern coasts.
//!
//! **The island rule.** A coast group's one unknown has one equation, the sum
//! of the equations of all its vertices: the union of their triangles is its
//! control volume. Within it, an edge joining two of its vertices adds
//! nothing, as both ends hold the same `Ψ`: the friction across it is zero,
//! and the `β` flux and stress leave one triangle and enter the other. What
//! is left are the sides crossing edges from a vertex of the group to one
//! outside it: the friction and `β` flux across the group's outer boundary
//! balance the circulation of the wind round it (Godfrey's island rule, here
//! with Stommel's friction). Each such edge has two wet cells, as a dry cell
//! of it would put both its ends in that cell's group, so the stress there
//! is the ocean's. The row is still one of an M-matrix.
//!
//! **Components.** Basins and the groups bordering them, joined where they
//! touch, make connected components, each one linear system: its unknowns
//! are its basins' interior vertices, ascending, then its groups but one,
//! ascending. The one left out, the reference, is the group with the largest
//! landmass area, and holds `Ψ = 0`. A component without a group, an
//! aquaplanet's single ocean, is skipped and keeps `Ψ = 0`: it would have no
//! boundary condition, and a zonal wind's curl has no Sverdrup balance round
//! a periodic ocean.
//!
//! **Solver.** The matrices depend only on the mesh, the land and the
//! rotation, so they are built once, with an incomplete LU factorisation
//! (ILU(0), on the matrix's own pattern) as a preconditioner. Each solve is
//! preconditioned BiCGSTAB from the last solution, to a relative residual
//! `‖b − AΨ‖/‖b‖` of [`TOLERANCE`]. Everything runs in serial, in index
//! order: the same input bits give the same output bits.

use fms_grid::zonal::east_north;
use fms_grid::{EdgeId, Operators, VertexId};

use crate::basins::Basins;

/// Relative residual `‖b − AΨ‖/‖b‖` a solve stops at.
pub const TOLERANCE: f64 = 1e-8;

/// Most BiCGSTAB iterations of one component's solve.
pub const MAX_ITERATIONS: usize = 1000;

/// The gyres' parameters: the solver's friction, and how the ocean
/// ([`Ocean`](crate::Ocean)) feeds it the wind and spreads its flow over its
/// layers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GyreParams {
    /// Stommel's linear bottom friction `r` (s⁻¹). The western boundary
    /// current is `r/β` wide, 100 km at 30° with the default; the mesh widens
    /// it to a vertex spacing or two (see the [module](self) documentation).
    pub friction: f64,
    /// Depth the gyre's flow is spread over, at one velocity (m). The
    /// ocean's mixed layer carries `mixed_depth / depth` of each edge's flow
    /// and its deep layer `deep_depth / depth`, both the same way; the rest
    /// flows below the model's ocean and carries no heat. The default, 700 m,
    /// is about the depth of the wind-driven thermocline.
    pub depth: f64,
    /// E-folding time of the running mean of the wind's stress that the
    /// solver sees (s): the gyres answer the wind of the last month or so,
    /// not of the day.
    pub averaging: f64,
    /// Model time between two solves (s), each from the last solution.
    pub interval: f64,
}

impl Default for GyreParams {
    fn default() -> Self {
        Self {
            friction: 2e-6,
            depth: 700.0,
            averaging: 30.0 * 86_400.0,
            interval: 86_400.0,
        }
    }
}

/// How one component's last solve went.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolveStats {
    /// BiCGSTAB iterations; 0 when the starting `Ψ` was already good enough
    /// or the wind was zero.
    pub iterations: usize,
    /// The final relative residual `‖b − AΨ‖/‖b‖`, computed afresh, not from
    /// the iteration's recurrence; 0 for a zero wind.
    pub residual: f64,
    /// Whether the residual is under [`TOLERANCE`]. Otherwise the solve
    /// stopped at [`MAX_ITERATIONS`].
    pub converged: bool,
}

/// A connected set of basins and the coast groups bordering them, solved
/// together.
#[derive(Clone, Debug, PartialEq)]
pub struct Component {
    /// Its kept basins, ascending, numbered as in [`Basins::basins`].
    pub basins: Vec<usize>,
    /// Its coast groups, ascending, numbered as in [`Basins::groups`].
    pub groups: Vec<usize>,
    /// The group held at `Ψ = 0`, the one with the largest landmass area
    /// (the first of `groups`); `None` for a component without coast, which
    /// isn't solved.
    pub reference: Option<usize>,
    /// Unknowns at interior vertices of its basins: all of them, or none
    /// when it isn't solved.
    pub vertex_unknowns: usize,
    /// Unknowns for its coast groups: all but the reference.
    pub group_unknowns: usize,
}

/// Each mesh edge's stress along its normal (N m⁻²), from each cell's
/// eastward and northward stress (N m⁻²): the mean of the two cells' stress
/// vectors, dotted with the edge's normal, as the ocean's Ekman transport
/// takes it. This is what [`Gyres::solve`] takes. A cell on a pole, where
/// east is undefined, adds nothing.
///
/// # Panics
/// If `stress` doesn't have one value per cell.
pub fn edge_stress(ops: &Operators, stress: &[[f64; 2]]) -> Vec<f64> {
    let (mesh, geometry) = (ops.mesh(), ops.geometry());
    assert_eq!(stress.len(), mesh.n_cells(), "one stress per cell");
    let centres = mesh.cell_centers();
    mesh.edges()
        .iter()
        .zip(&geometry.edge_normals)
        .map(|(edge, &normal)| {
            let along = |c: u32| {
                // The normal's east and north components at the cell centre.
                let (east, north) = east_north(centres[c as usize], normal);
                let [e, n] = stress[c as usize];
                e * east + n * north
            };
            0.5 * (along(edge.cells[0]) + along(edge.cells[1]))
        })
        .collect()
}

/// One side of a vertex's control volume, as the gyre equation sees it.
#[derive(Clone, Copy, Debug)]
struct Side {
    /// The mesh edge the side crosses.
    edge: EdgeId,
    /// The vertex at the other end of it.
    neighbour: VertexId,
    /// `dₑ/lₑ`: the side's length over the mesh edge's.
    conductance: f64,
    /// The flux of `c = −k × ∇f` out of the vertex's triangle across the side
    /// (s⁻¹): `t_{e,v} (f[cells[1]] − f[cells[0]])`.
    flux: f64,
    /// `−t_{e,v} dₑ R` (m): times the edge's stress along its normal, the
    /// side's share of the right-hand side.
    stress_weight: f64,
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
            edge: edges[k],
            neighbour: if edge.vertices[0] == v {
                edge.vertices[1]
            } else {
                edge.vertices[0]
            },
            conductance: geometry.dual_edge_lengths[e] / geometry.edge_lengths[e],
            flux: signs[k] * (coriolis[b as usize] - coriolis[a as usize]),
            stress_weight: -signs[k] * geometry.dual_edge_lengths[e] * ops.radius(),
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

/// One component's linear system.
#[derive(Clone, Debug)]
struct System {
    /// The interior vertices that are unknowns, ascending: unknown `i` is `Ψ`
    /// at `vertices[i]`.
    vertices: Vec<VertexId>,
    /// The groups that are unknowns, ascending: unknown `vertices.len() + k`
    /// is `Ψ` of group `groups[k]`.
    groups: Vec<usize>,
    matrix: Csr,
    /// The ILU(0) factors, on `matrix`'s pattern.
    factors: Vec<f64>,
    /// The right-hand side: unknown `i`'s is the sum over
    /// `stencil[stencil_offsets[i]..stencil_offsets[i + 1]]` of each edge's
    /// stress along its normal times the weight (m).
    stencil_offsets: Vec<u32>,
    stencil: Vec<(EdgeId, f64)>,
}

/// Sort a row's `(column, value)` entries by column and sum those sharing
/// one, in the order they came.
fn merge(entries: &mut Vec<(u32, f64)>) {
    entries.sort_by_key(|&(column, _)| column);
    let mut merged = 0;
    for k in 0..entries.len() {
        if merged > 0 && entries[merged - 1].0 == entries[k].0 {
            entries[merged - 1].1 += entries[k].1;
        } else {
            entries[merged] = entries[k];
            merged += 1;
        }
    }
    entries.truncate(merged);
}

impl System {
    /// The system of `component`, with `f` at the cell centres.
    fn new(
        ops: &Operators,
        coriolis: &[f64],
        basins: &Basins,
        component: &Component,
        friction: f64,
    ) -> Self {
        const NONE: u32 = u32::MAX;
        let mut vertices: Vec<VertexId> = component
            .basins
            .iter()
            .flat_map(|&b| basins.basins()[b].vertices.iter().copied())
            .collect();
        vertices.sort_unstable();
        let groups: Vec<usize> = component.groups.iter().skip(1).copied().collect();
        // The unknown holding each vertex's Ψ; NONE where it is 0.
        let mut unknown = vec![NONE; ops.mesh().n_vertices()];
        for (i, &v) in vertices.iter().enumerate() {
            unknown[v as usize] = i as u32;
        }
        for (k, &g) in groups.iter().enumerate() {
            for &v in &basins.groups()[g].vertices {
                unknown[v as usize] = (vertices.len() + k) as u32;
            }
        }
        let n = vertices.len() + groups.len();
        let mut offsets = Vec::with_capacity(n + 1);
        offsets.push(0);
        let (mut columns, mut values) = (Vec::with_capacity(4 * n), Vec::with_capacity(4 * n));
        let mut diagonal = Vec::with_capacity(n);
        let mut stencil_offsets = Vec::with_capacity(n + 1);
        stencil_offsets.push(0);
        let mut stencil = Vec::with_capacity(3 * n);
        let mut entries: Vec<(u32, f64)> = Vec::new();
        // Each unknown's control volume: one vertex's triangle, or the
        // triangles of all of a group's vertices.
        let interior = vertices.iter().map(|v| (std::slice::from_ref(v), None));
        let coasts = groups
            .iter()
            .map(|&g| (&basins.groups()[g].vertices[..], Some(g)));
        for (i, (members, group)) in interior.chain(coasts).enumerate() {
            entries.clear();
            let mut centre = 0.0;
            for &v in members {
                for side in sides(ops, coriolis, v) {
                    // Inside a group both ends hold its Ψ: nothing crosses.
                    if group.is_some() && basins.vertex_group(side.neighbour) == group {
                        continue;
                    }
                    let diffusion = friction * side.conductance;
                    let outflow = -side.flux;
                    centre += diffusion + outflow.max(0.0);
                    let w = unknown[side.neighbour as usize];
                    if w != NONE {
                        entries.push((w, -diffusion + outflow.min(0.0)));
                    }
                    stencil.push((side.edge, side.stress_weight));
                }
            }
            entries.push((i as u32, centre));
            merge(&mut entries);
            for &(column, value) in &entries {
                if column == i as u32 {
                    diagonal.push(columns.len() as u32);
                }
                columns.push(column);
                values.push(value);
            }
            offsets.push(columns.len() as u32);
            stencil_offsets.push(stencil.len() as u32);
        }
        let matrix = Csr {
            offsets,
            columns,
            values,
            diagonal,
        };
        debug_assert_eq!(matrix.rows(), n, "every row has its diagonal");
        let factors = matrix.ilu0();
        Self {
            vertices,
            groups,
            matrix,
            factors,
            stencil_offsets,
            stencil,
        }
    }

    /// Solve for `Ψ` with `stress` along each edge's normal (N m⁻²),
    /// starting from and writing to `psi` (kg s⁻¹, per vertex) and
    /// `group_psi` (kg s⁻¹, per group). The groups' vertices are left for the
    /// caller.
    fn solve(&self, stress: &[f64], psi: &mut [f64], group_psi: &mut [f64]) -> SolveStats {
        let b: Vec<f64> = self
            .stencil_offsets
            .windows(2)
            .map(|w| {
                self.stencil[w[0] as usize..w[1] as usize]
                    .iter()
                    .map(|&(e, weight)| weight * stress[e as usize])
                    .sum()
            })
            .collect();
        let mut x: Vec<f64> = self
            .vertices
            .iter()
            .map(|&v| psi[v as usize])
            .chain(self.groups.iter().map(|&g| group_psi[g]))
            .collect();
        let stats = bicgstab(&self.matrix, &self.factors, &b, &mut x);
        let (at_vertices, at_groups) = x.split_at(self.vertices.len());
        for (&v, &x) in self.vertices.iter().zip(at_vertices) {
            psi[v as usize] = x;
        }
        for (&g, &x) in self.groups.iter().zip(at_groups) {
            group_psi[g] = x;
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

/// The basins and coast groups of `basins`, gathered into connected
/// components: a basin and a group are joined where one of the basin's
/// vertices is a neighbour of one of the group's, two groups where two of
/// their vertices are neighbours (a strait one edge wide). Numbered in the
/// order of their first basin, then of their first group.
fn components(ops: &Operators, basins: &Basins) -> Vec<Component> {
    let mesh = ops.mesh();
    let n_basins = basins.basins().len();
    let n = n_basins + basins.groups().len();
    // Union–find over the nodes: basins first, then groups.
    let mut parent: Vec<usize> = (0..n).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    let node = |v: VertexId| {
        basins
            .vertex_basin(v)
            .or_else(|| basins.vertex_group(v).map(|g| n_basins + g))
    };
    for v in 0..mesh.n_vertices() as VertexId {
        let Some(a) = node(v) else { continue };
        for e in mesh.vertex_edges(v) {
            let ends = mesh.edges()[e as usize].vertices;
            let w = if ends[0] == v { ends[1] } else { ends[0] };
            if let Some(b) = node(w) {
                let (a, b) = (root(&mut parent, a), root(&mut parent, b));
                parent[a.max(b)] = a.min(b);
            }
        }
    }
    let mut id = vec![usize::MAX; n];
    let mut found: Vec<Component> = Vec::new();
    for i in 0..n {
        let r = root(&mut parent, i);
        if id[r] == usize::MAX {
            id[r] = found.len();
            found.push(Component {
                basins: Vec::new(),
                groups: Vec::new(),
                reference: None,
                vertex_unknowns: 0,
                group_unknowns: 0,
            });
        }
        let c = &mut found[id[r]];
        if i < n_basins {
            c.basins.push(i);
            c.vertex_unknowns += basins.basins()[i].vertices.len();
        } else {
            c.groups.push(i - n_basins);
        }
    }
    for c in &mut found {
        // Groups are numbered largest first.
        c.reference = c.groups.first().copied();
        c.group_unknowns = c.groups.len().saturating_sub(1);
        if c.reference.is_none() {
            c.vertex_unknowns = 0;
        }
    }
    found
}

/// The gyres' streamfunction on one mesh, land and rotation.
#[derive(Clone, Debug)]
pub struct Gyres {
    params: GyreParams,
    components: Vec<Component>,
    /// One per component; `None` for one without coast, which isn't solved.
    systems: Vec<Option<System>>,
    stats: Vec<Option<SolveStats>>,
    /// Per vertex (kg s⁻¹).
    psi: Vec<f64>,
    /// Per coast group (kg s⁻¹).
    group_psi: Vec<f64>,
    /// Each coast group's vertices.
    group_vertices: Vec<Vec<VertexId>>,
    /// Each coast group's component.
    group_component: Vec<usize>,
    /// Each mesh edge's two vertices.
    ends: Vec<[VertexId; 2]>,
}

impl Gyres {
    /// The gyres of the `basins` on the mesh of `ops`, on a planet turning at
    /// `rotation` (rad s⁻¹, negative for retrograde), with `Ψ = 0` until the
    /// first [`solve`](Self::solve). Finds the components and builds and
    /// factorises each one's matrix.
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
        let components = components(ops, basins);
        let systems: Vec<Option<System>> = components
            .iter()
            .map(|c| {
                c.reference
                    .map(|_| System::new(ops, &coriolis, basins, c, params.friction))
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
        let mut group_component = vec![0; basins.groups().len()];
        for (k, c) in components.iter().enumerate() {
            for &g in &c.groups {
                group_component[g] = k;
            }
        }
        Self {
            params,
            components,
            systems,
            stats,
            psi: vec![0.0; mesh.n_vertices()],
            group_psi: vec![0.0; basins.groups().len()],
            group_vertices: basins.groups().iter().map(|g| g.vertices.clone()).collect(),
            group_component,
            ends: mesh.edges().iter().map(|e| e.vertices).collect(),
        }
    }

    /// The parameters.
    pub fn params(&self) -> &GyreParams {
        &self.params
    }

    /// The connected components, in the order of their first basin, then of
    /// their first group.
    pub fn components(&self) -> &[Component] {
        &self.components
    }

    /// Solve every component with a coast for the wind's `stress` along each
    /// mesh edge's normal (N m⁻², one per edge, from its `cells[0]` to its
    /// `cells[1]`; see [`edge_stress`]), each from its last solution. Values
    /// on edges with a cell that isn't wet, and on edges no solved
    /// component's equation reaches, are ignored.
    ///
    /// # Panics
    /// If `stress` doesn't have one value per edge.
    pub fn solve(&mut self, stress: &[f64]) {
        assert_eq!(stress.len(), self.ends.len(), "one stress per edge");
        for (system, stats) in self.systems.iter().zip(&mut self.stats) {
            if let Some(system) = system {
                *stats = Some(system.solve(stress, &mut self.psi, &mut self.group_psi));
                for &g in &system.groups {
                    for &v in &self.group_vertices[g] {
                        self.psi[v as usize] = self.group_psi[g];
                    }
                }
            }
        }
    }

    /// The streamfunction at each vertex (kg s⁻¹): its group's on coasts,
    /// zero on the reference group's and in components that aren't solved.
    pub fn streamfunction(&self) -> &[f64] {
        &self.psi
    }

    /// The streamfunction of each coast group, numbered as in
    /// [`Basins::groups`] (kg s⁻¹): zero for each component's reference.
    pub fn group_streamfunction(&self) -> &[f64] {
        &self.group_psi
    }

    /// The component of coast group `group`.
    pub fn group_component(&self, group: usize) -> usize {
        self.group_component[group]
    }

    /// The mass flowing across `edge` from its `cells[0]` to its `cells[1]`
    /// (kg s⁻¹): `Ψ[vertices[0]] − Ψ[vertices[1]]`.
    ///
    /// The edge's tangent `tₑ = k × nₑ` runs from `vertices[0]` to
    /// `vertices[1]`, so the transport along the normal is
    /// `(k × ∇Ψ)·nₑ = −∇Ψ·tₑ`, and times the edge's length that is the
    /// difference. Summed round a cell, these cancel: the flow has no
    /// divergence. Across an edge with a dry cell both vertices are in that
    /// cell's group, so it is exactly zero.
    pub fn edge_flux(&self, edge: EdgeId) -> f64 {
        let [a, b] = self.ends[edge as usize];
        self.psi[a as usize] - self.psi[b as usize]
    }

    /// How each component's last solve went, in [`components`](Self::components)
    /// order; `None` for a component without coast, which isn't solved.
    pub fn stats(&self) -> &[Option<SolveStats>] {
        &self.stats
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::EARTH_ROTATION;
    use fms_grid::coords::{from_lat_lon, to_lat_lon};
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

        /// The test wind's stress along each edge's normal (N m⁻²): `τ_x`
        /// at the edge point times its eastward unit vector's component
        /// along the normal.
        fn wind(&self, g: &Grid) -> Vec<f64> {
            zonal_wind(g, |lat| self.stress(lat).0)
        }

        /// Whether `(lat, lon)` (radians) is at least `margin` degrees inside
        /// the box from the north and south and in the given longitudes.
        fn holds(&self, (lat, lon): (f64, f64), margin: f64, lons: (f64, f64)) -> bool {
            let (lat, lon) = (lat.to_degrees(), lon.to_degrees());
            (self.south + margin..=self.north - margin).contains(&lat)
                && (lons.0..=lons.1).contains(&lon)
        }
    }

    /// The stress along each edge's normal (N m⁻²) of the eastward stress
    /// `tau(φ)` at the edge points.
    fn zonal_wind(g: &Grid, tau: impl Fn(f64) -> f64) -> Vec<f64> {
        g.geometry
            .edge_points
            .iter()
            .zip(&g.geometry.edge_normals)
            .map(|(&p, &n)| {
                let east = DVec3::Z.cross(p).normalize_or_zero();
                tau(to_lat_lon(p).0) * east.dot(n)
            })
            .collect()
    }

    /// The solved gyres of a box.
    fn solve(g: &Grid, b: Rectangle, rotation: f64) -> (Basins, Gyres) {
        let basins = Basins::new(&g.ops(), &b.land(g));
        assert_eq!(basins.basins().len(), 1);
        let mut gyres = Gyres::new(&g.ops(), &basins, rotation, GyreParams::default());
        gyres.solve(&b.wind(g));
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
            let params = GyreParams {
                friction,
                ..GyreParams::default()
            };
            let mut gyres = Gyres::new(&g.ops(), &basins, OMEGA, params);
            let curl = NORTH.curl(g);
            gyres.solve(&NORTH.wind(g));
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
        let wind = NORTH.wind(&g);
        gyres.solve(&wind);
        let once = gyres.streamfunction().to_vec();
        gyres.solve(&vec![0.0; wind.len()]);
        assert!(gyres.streamfunction().iter().all(|&p| p == 0.0));
        let stats = gyres.stats()[0].unwrap();
        assert_eq!((stats.iterations, stats.residual), (0, 0.0));
        // From Ψ = 0, every step of the solve scales by 2 exactly.
        let doubled: Vec<f64> = wind.iter().map(|c| 2.0 * c).collect();
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
        gyres.solve(&NORTH.wind(&g));
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
        gyres.solve(&NORTH.wind(&g1));
        let again = gyres.stats()[0].unwrap();
        assert!(
            again.iterations <= 1 && again.residual < TOLERANCE,
            "{again:?}"
        );
        let (e, m) = worst(gyres.streamfunction(), &before, |_| true);
        assert!(e <= 1e-6 * m);
    }

    /// The angle between two unit vectors (radians).
    fn angle(a: DVec3, b: DVec3) -> f64 {
        fms_math::atan2(a.cross(b).length(), a.dot(b))
    }

    /// Land fraction 1 where `land(lat, lon)` (degrees) holds at the cell
    /// centre, else 0.
    fn land_where(g: &Grid, land: impl Fn(f64, f64) -> bool) -> Vec<f64> {
        g.mesh
            .cell_centers()
            .iter()
            .map(|&p| {
                let (lat, lon) = to_lat_lon(p);
                if land(lat.to_degrees(), lon.to_degrees()) {
                    1.0
                } else {
                    0.0
                }
            })
            .collect()
    }

    /// Each component's matrix is an M-matrix: positive diagonal, nothing
    /// positive off it, every row summing to at least zero. Returns the
    /// number of rows summing to more than zero.
    fn assert_m_matrix(a: &Csr) -> usize {
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
        strict
    }

    /// The box with an island: its cells within `radius` degrees of 30° N
    /// 31° E, in the middle of the gyre's interior.
    fn island(g: &Grid, radius: f64) -> Vec<f64> {
        let centre = from_lat_lon(30f64.to_radians(), 31f64.to_radians());
        let mut land = NORTH.land(g);
        for (c, &p) in g.mesh.cell_centers().iter().enumerate() {
            if angle(p, centre).to_degrees() < radius {
                land[c] = 1.0;
            }
        }
        land
    }

    /// The solved gyres for `land` and the edge stress `wind`.
    fn solve_land(g: &Grid, land: &[f64], wind: &[f64]) -> (Basins, Gyres) {
        let basins = Basins::new(&g.ops(), land);
        let mut gyres = Gyres::new(&g.ops(), &basins, OMEGA, GyreParams::default());
        gyres.solve(wind);
        (basins, gyres)
    }

    #[test]
    fn an_island_in_the_gyre_has_flow_on_both_sides() {
        let g = Grid::new(24);
        let (basins, gyres) = solve_land(&g, &island(&g, 4.0), &NORTH.wind(&g));
        assert_eq!(basins.basins().len(), 1);
        assert_eq!(basins.groups().len(), 2);
        let c = &gyres.components()[0];
        assert_eq!(gyres.components().len(), 1);
        assert_eq!((c.reference, c.group_unknowns), (Some(0), 1));
        let island = &basins.groups()[1];
        assert!(island.landmasses.len() == 1 && island.vertices.len() >= 12);
        let psi_island = gyres.group_streamfunction()[1];
        let psi = gyres.streamfunction();
        assert!(
            island
                .vertices
                .iter()
                .all(|&v| psi[v as usize] == psi_island)
        );
        // Without the island, Ψ over its vertices spans `plain`.
        let (_, without) = solve(&g, NORTH, OMEGA);
        let plain = island
            .vertices
            .iter()
            .map(|&v| without.streamfunction()[v as usize])
            .fold((f64::MAX, f64::MIN), |(a, b), p| (a.min(p), b.max(p)));
        // The ring of vertices round the island.
        let ring = island
            .vertices
            .iter()
            .flat_map(|&v| g.mesh.vertex_edges(v))
            .flat_map(|e| g.mesh.edges()[e as usize].vertices)
            .filter(|&w| basins.vertex_group(w).is_none())
            .map(|w| psi[w as usize])
            .fold((f64::MAX, f64::MIN), |(a, b), p| (a.min(p), b.max(p)));
        // Ψ of the island (6 cells, 22 vertices) lies between its
        // neighbours' and within the range the gyre without it has there:
        // 9.29 Sv, against 7.11 to 9.92 round it and 7.27 to 10.41 (mean
        // 8.70) on its vertices without it.
        assert!(ring.0 < psi_island && psi_island < ring.1, "{ring:?}");
        assert!(plain.0 < psi_island && psi_island < plain.1, "{plain:?}");
        // So Ψ_I flows south between the island and the eastern coast,
        // and as much north between it and the western one, where without
        // the island rule nothing would pass on the eastern side.
        assert!(psi_island > 5e9, "{psi_island}");
    }

    /// Land poleward of 65° S and equatorward of 45° S, leaving a channel
    /// round the planet; with `barrier`, also a ridge across it at the
    /// meridian 0°, 12° wide, from the north side down to 55° S, leaving a
    /// gap 10° wide next to the cap.
    fn channel(g: &Grid, barrier: bool) -> Vec<f64> {
        land_where(g, |lat, lon| {
            lat < -65.0 || lat > -45.0 || (barrier && lon.abs() < 6.0 && lat > -55.0)
        })
    }

    /// The channel with two barriers, 12° wide, between which no latitude
    /// circle is open: one from the north side at the meridian 0° down to
    /// 60° S, one from the cap at 180° up to 58° S.
    fn staggered(g: &Grid) -> Vec<f64> {
        land_where(g, |lat, lon| {
            lat < -65.0
                || lat > -45.0
                || (lon.abs() < 6.0 && lat > -60.0)
                || (lon.abs() > 174.0 && lat < -58.0)
        })
    }

    /// The stress of an eastward wind `τ₀` along each edge's normal.
    fn eastward(g: &Grid) -> Vec<f64> {
        zonal_wind(g, |_| TAU0)
    }

    /// The transport round a channel, eastward (kg s⁻¹): Ψ of the cap,
    /// against the continent's 0.
    fn transport(basins: &Basins, gyres: &Gyres) -> f64 {
        assert_eq!(basins.groups().len(), 2);
        assert_eq!(gyres.components().len(), 1);
        assert_eq!(gyres.components()[0].reference, Some(0));
        // The continent is the larger, holding the north pole's latitudes.
        assert!(basins.groups()[0].area > basins.groups()[1].area);
        gyres.group_streamfunction()[1]
    }

    #[test]
    fn an_open_channel_carries_the_wind_over_the_friction() {
        // Zonally symmetric, r M_x = τ_x, so round a channel 20° wide the
        // flow is τ₀ R (20°) / r = 111 Sv. The scheme's numerical friction
        // takes some, as the mesh's sides don't follow the latitudes: the
        // ratio is 0.66 at 2,562 cells, 0.74 at 5,762, 0.80 at 10,242 and
        // 0.89 at 40,962, the shortfall halving with the spacing.
        let expected = TAU0 * RADIUS * 20f64.to_radians() / GyreParams::default().friction;
        let ratio = |k: u32, low: f64, high: f64| {
            let g = Grid::new(k);
            let (basins, gyres) = solve_land(&g, &channel(&g, false), &eastward(&g));
            let ratio = transport(&basins, &gyres) / expected;
            assert!((low..high).contains(&ratio), "{k}: {ratio}");
            ratio
        };
        let (coarse, fine) = (ratio(16, 0.6, 0.72), ratio(32, 0.75, 0.86));
        assert!(1.0 - fine < 0.65 * (1.0 - coarse), "{coarse} {fine}");
    }

    #[test]
    fn a_barrier_with_a_gap_lets_most_of_the_flow_through() {
        // The gap's latitudes, 55° S to 65° S, are open all round, and
        // along them only friction holds the flow back, as round the open
        // channel; the water from the rest goes round the barrier. So the
        // barrier takes little: 69.1 Sv through the gap against 73.2 round
        // the open channel at 2,562 cells, 80.3 against 89.3 at 10,242.
        // With no latitude open (two staggered barriers) `β` holds back
        // more, but not most: 51.2 and 64.4 Sv. This is the flat-bottomed
        // ocean's circumpolar current, which only friction limits.
        for k in [16, 32] {
            let g = Grid::new(k);
            let wind = eastward(&g);
            let flow = |land: &[f64]| {
                let (basins, gyres) = solve_land(&g, land, &wind);
                transport(&basins, &gyres)
            };
            let (open, gap, closed) = (
                flow(&channel(&g, false)),
                flow(&channel(&g, true)),
                flow(&staggered(&g)),
            );
            assert!(
                0.85 * open < gap && gap < 0.97 * open,
                "{k}: {gap} of {open}"
            );
            assert!(
                0.6 * open < closed && closed < 0.8 * open,
                "{k}: {closed} of {open}"
            );
        }
    }

    #[test]
    fn group_rows_are_the_sums_of_their_vertices_rows() {
        for (k, land) in [(24, true), (16, false)] {
            let g = Grid::new(k);
            let ops = g.ops();
            let land = if land {
                island(&g, 4.0)
            } else {
                channel(&g, true)
            };
            let (basins, gyres) = solve_land(&g, &land, &NORTH.wind(&g));
            let coriolis: Vec<f64> = g
                .mesh
                .cell_centers()
                .iter()
                .map(|p| 2.0 * OMEGA * p.z)
                .collect();
            let system = gyres.systems[0].as_ref().unwrap();
            let a = &system.matrix;
            let strict = assert_m_matrix(a);
            assert!(strict > 0 && strict < a.rows());
            // Each vertex's unknown, as the system numbers them.
            let mut unknown = vec![u32::MAX; g.mesh.n_vertices()];
            for (i, &v) in system.vertices.iter().enumerate() {
                unknown[v as usize] = i as u32;
            }
            for (k, &group) in system.groups.iter().enumerate() {
                for &v in &basins.groups()[group].vertices {
                    unknown[v as usize] = (system.vertices.len() + k) as u32;
                }
            }
            let stress = NORTH.wind(&g);
            let friction = GyreParams::default().friction;
            for (k, &group) in system.groups.iter().enumerate() {
                let i = system.vertices.len() + k;
                // The sum of the full rows of the group's vertices, inner
                // edges and all, with Ψ at each vertex its unknown's.
                let mut row = vec![0.0; a.rows()];
                let (mut b, mut size) = (0.0, 0.0);
                for &v in &basins.groups()[group].vertices {
                    for side in sides(&ops, &coriolis, v) {
                        let (diffusion, outflow) = (friction * side.conductance, -side.flux);
                        row[i] += diffusion + outflow.max(0.0);
                        let w = unknown[side.neighbour as usize];
                        if w != u32::MAX {
                            row[w as usize] += -diffusion + outflow.min(0.0);
                        }
                        let term = side.stress_weight * stress[side.edge as usize];
                        b += term;
                        size += term.abs();
                    }
                }
                let scale = row[i];
                for p in a.row(i) {
                    let column = a.columns[p] as usize;
                    assert!((a.values[p] - row[column]).abs() <= 1e-12 * scale);
                    row[column] = 0.0;
                }
                assert!(row.iter().all(|x| x.abs() <= 1e-12 * scale));
                let s = &system.stencil
                    [system.stencil_offsets[i] as usize..system.stencil_offsets[i + 1] as usize];
                let ours: f64 = s.iter().map(|&(e, w)| w * stress[e as usize]).sum();
                assert!((ours - b).abs() <= 1e-12 * size, "{ours} {b}");
            }
            // Inside a group, an edge's two sides cancel exactly, and no
            // water crosses it.
            for (e, edge) in g.mesh.edges().iter().enumerate() {
                let [v, w] = edge.vertices;
                let group = basins.vertex_group(v);
                if group.is_none() || group != basins.vertex_group(w) {
                    continue;
                }
                let side = |v: VertexId| {
                    sides(&ops, &coriolis, v)
                        .into_iter()
                        .find(|s| s.edge == e as EdgeId)
                        .unwrap()
                };
                let (a, b) = (side(v), side(w));
                assert_eq!(a.conductance, b.conductance);
                assert_eq!(a.flux, -b.flux);
                assert_eq!(a.stress_weight, -b.stress_weight);
                assert_eq!(gyres.edge_flux(e as EdgeId), 0.0);
            }
            // And the flow still has no divergence and stays off land.
            let largest = (0..g.mesh.n_edges())
                .map(|e| gyres.edge_flux(e as EdgeId).abs())
                .fold(0.0, f64::max);
            for c in 0..g.mesh.n_cells() as CellId {
                let out: f64 = g
                    .mesh
                    .cell_edges(c)
                    .iter()
                    .zip(g.geometry.cell_edge_signs(c))
                    .map(|(&e, s)| s * gyres.edge_flux(e))
                    .sum();
                assert!(out.abs() < 1e-14 * largest, "cell {c}: {out}");
                if !basins.is_wet(c) {
                    for &e in g.mesh.cell_edges(c) {
                        assert_eq!(gyres.edge_flux(e), 0.0);
                    }
                }
            }
        }
    }

    #[test]
    fn a_strait_one_edge_wide_carries_the_difference() {
        // Two one-cell islands across the ends of one edge between wet
        // cells: its two vertices are in different groups.
        let g = Grid::new(24);
        let centre = from_lat_lon(30f64.to_radians(), 40f64.to_radians());
        let e = (0..g.mesh.n_edges())
            .max_by(|&a, &b| {
                let p = &g.geometry.edge_points;
                p[a].dot(centre).total_cmp(&p[b].dot(centre))
            })
            .unwrap();
        let edge = g.mesh.edges()[e];
        let mut land = NORTH.land(&g);
        let islands = edge.vertices.map(|v| {
            let third = g.mesh.triangles()[v as usize]
                .into_iter()
                .find(|c| !edge.cells.contains(c))
                .unwrap();
            land[third as usize] = 1.0;
            third
        });
        assert!(!g.mesh.cell_neighbors(islands[0]).contains(&islands[1]));
        let (basins, gyres) = solve_land(&g, &land, &NORTH.wind(&g));
        assert!(edge.cells.iter().all(|&c| basins.is_wet(c)));
        let [a, b] = edge.vertices.map(|v| basins.vertex_group(v).unwrap());
        assert_ne!(a, b);
        assert_eq!(basins.groups().len(), 3);
        assert_eq!(gyres.components().len(), 1);
        let psi = gyres.group_streamfunction();
        let flow = gyres.edge_flux(e as EdgeId);
        assert_eq!(flow, psi[a] - psi[b]);
        assert!(flow.abs() > 1e6, "{flow}");
        // The two groups' rows are coupled through the edge.
        let system = gyres.systems[0].as_ref().unwrap();
        let n = system.vertices.len() as u32;
        let row = |g: usize| n + system.groups.iter().position(|&x| x == g).unwrap() as u32;
        let a_row = system.matrix.row(row(a) as usize);
        assert!(
            a_row
                .clone()
                .any(|k| system.matrix.columns[k] == row(b) && system.matrix.values[k] < 0.0)
        );
        assert_m_matrix(&system.matrix);
    }

    #[test]
    fn the_reference_is_the_largest_landmass() {
        // Land south of 20° S and north of 50° N: the south is the larger.
        let g = Grid::new(16);
        let land = land_where(&g, |lat, _| !(-20.0..=50.0).contains(&lat));
        let (basins, gyres) = solve_land(&g, &land, &eastward(&g));
        assert_eq!(basins.groups().len(), 2);
        let south = (0..g.mesh.n_cells() as CellId)
            .find(|&c| g.mesh.cell_centers()[c as usize].z < -0.99)
            .unwrap();
        assert_eq!(basins.cell_group(south), Some(0));
        assert_eq!(gyres.components()[0].reference, Some(0));
        let psi = gyres.streamfunction();
        assert!(
            basins.groups()[0]
                .vertices
                .iter()
                .all(|&v| psi[v as usize] == 0.0)
        );
        assert_ne!(gyres.group_streamfunction()[1], 0.0);
        assert_eq!(gyres.group_streamfunction()[0], 0.0);
    }

    #[test]
    fn groups_solve_the_same_twice_and_linearly() {
        let (g1, g2) = (Grid::new(16), Grid::new(16));
        let solved = |g: &Grid, scale: f64| {
            let wind: Vec<f64> = eastward(g).iter().map(|w| scale * w).collect();
            solve_land(g, &channel(g, true), &wind).1
        };
        let (once, again) = (solved(&g1, 1.0), solved(&g2, 1.0));
        assert_eq!(once.streamfunction(), again.streamfunction());
        assert_eq!(once.group_streamfunction(), again.group_streamfunction());
        let stats = once.stats()[0].unwrap();
        assert!(stats.converged && stats.residual < TOLERANCE, "{stats:?}");
        let twice = solved(&g1, 2.0);
        assert!(
            once.streamfunction()
                .iter()
                .zip(twice.streamfunction())
                .all(|(a, b)| 2.0 * a == *b)
        );
    }

    #[test]
    fn the_edge_stress_of_the_cells_stress_is_the_winds() {
        // The box's wind given at the cell centres, eastward and northward:
        // its mean along each edge's normal is the wind at the edge point
        // to second order, and the gyre barely changes.
        let g = Grid::new(24);
        let ops = g.ops();
        let cells: Vec<[f64; 2]> = g
            .mesh
            .cell_centers()
            .iter()
            .map(|&p| [NORTH.stress(to_lat_lon(p).0).0, 0.0])
            .collect();
        let from_cells = edge_stress(&ops, &cells);
        let at_edges = NORTH.wind(&g);
        // On the box's edges: 0.24% of the largest. (Near the poles, whose
        // cells' east swings round, it is worse, but no gyre is there.)
        let land = NORTH.land(&g);
        let wet = |e: usize| {
            g.mesh.edges()[e]
                .cells
                .iter()
                .all(|&c| land[c as usize] == 0.0)
        };
        let (e, m) = worst(&from_cells, &at_edges, wet);
        assert!(e < 0.004 * m, "{e} of {m}");
        let (a, b) = (
            solve_land(&g, &land, &from_cells).1,
            solve_land(&g, &land, &at_edges).1,
        );
        // And Ψ by 0.19% of the largest.
        let (e, m) = worst(a.streamfunction(), b.streamfunction(), |_| true);
        assert!(e < 0.004 * m, "{e} of {m}");
    }
}
