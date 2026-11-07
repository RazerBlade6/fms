//! Mapping between the fine terrain grid and the coarse climate grid.
//!
//! The two meshes are optimized independently, so they are not exactly nested.
//! Each fine cell is assigned wholly to the coarse cell containing its centre
//! (its *parent*). Aggregation (fine → coarse) uses these parent/child sets;
//! interpolation (coarse → fine) uses barycentric weights on the coarse Delaunay
//! triangle containing each fine cell centre.

use rayon::prelude::*;

use crate::locate::CellLocator;
use crate::mesh::{CellId, SphereMesh};

/// Barycentric interpolation stencil: three coarse cells and their weights.
pub type Stencil = ([CellId; 3], [f64; 3]);

/// Precomputed fine ↔ coarse mapping.
#[derive(Clone, Debug)]
pub struct TwoLevelMap {
    /// Coarse parent of each fine cell.
    parent: Vec<CellId>,
    /// CSR offsets into `children`; length `n_coarse + 1`.
    child_offsets: Vec<u32>,
    /// Fine children of each coarse cell, sorted by index.
    children: Vec<CellId>,
    /// Area of each child as a fraction of its parent's total child area,
    /// aligned with `children`.
    child_weights: Vec<f64>,
    /// For each fine cell: the three coarse cells around it and their
    /// barycentric weights.
    interpolation: Vec<Stencil>,
}

impl TwoLevelMap {
    /// Build the mapping. `fine_areas` are the fine cells' areas (any consistent
    /// unit), used to weight aggregation.
    ///
    /// # Panics
    /// If some coarse cell receives no fine cell. That cannot happen when the fine
    /// mesh is at least about twice the frequency of the coarse one.
    pub fn new(coarse: &SphereMesh, fine: &SphereMesh, fine_areas: &[f64]) -> Self {
        assert_eq!(fine_areas.len(), fine.n_cells());
        let locator = CellLocator::new(coarse);

        let located: Vec<(CellId, Stencil)> = fine
            .cell_centers()
            .par_iter()
            .map(|&p| {
                let parent = locator.locate(coarse, p);
                let (t, w) = locator.locate_triangle(coarse, p);
                (parent, (coarse.triangles()[t as usize], w))
            })
            .collect();
        let (parent, interpolation): (Vec<CellId>, Vec<_>) = located.into_iter().unzip();

        let n_coarse = coarse.n_cells();
        let mut child_offsets = vec![0u32; n_coarse + 1];
        for &p in &parent {
            child_offsets[p as usize + 1] += 1;
        }
        for i in 0..n_coarse {
            child_offsets[i + 1] += child_offsets[i];
        }
        let mut fill = child_offsets.clone();
        let mut children = vec![0 as CellId; fine.n_cells()];
        for (f, &p) in parent.iter().enumerate() {
            children[fill[p as usize] as usize] = f as CellId;
            fill[p as usize] += 1;
        }
        if let Some(empty) = (0..n_coarse).find(|&c| child_offsets[c] == child_offsets[c + 1]) {
            panic!(
                "coarse cell {empty} has no fine children; the fine mesh is too coarse \
                 (frequency {} vs {})",
                fine.frequency(),
                coarse.frequency()
            );
        }

        let mut child_weights = vec![0.0; children.len()];
        for c in 0..n_coarse {
            let range = child_offsets[c] as usize..child_offsets[c + 1] as usize;
            let total: f64 = children[range.clone()]
                .iter()
                .map(|&f| fine_areas[f as usize])
                .sum();
            for i in range {
                child_weights[i] = fine_areas[children[i] as usize] / total;
            }
        }

        Self {
            parent,
            child_offsets,
            children,
            child_weights,
            interpolation,
        }
    }

    /// Number of coarse cells.
    pub fn n_coarse(&self) -> usize {
        self.child_offsets.len() - 1
    }

    /// Number of fine cells.
    pub fn n_fine(&self) -> usize {
        self.parent.len()
    }

    /// Coarse parent of each fine cell.
    pub fn parents(&self) -> &[CellId] {
        &self.parent
    }

    /// Fine children of `coarse_cell`, sorted by index.
    pub fn children(&self, coarse_cell: CellId) -> &[CellId] {
        let c = coarse_cell as usize;
        &self.children[self.child_offsets[c] as usize..self.child_offsets[c + 1] as usize]
    }

    /// Interpolation stencil of `fine_cell`: three coarse cells and weights.
    pub fn interpolation(&self, fine_cell: CellId) -> Stencil {
        self.interpolation[fine_cell as usize]
    }

    /// Area-weighted mean of an *intensive* fine field (temperature, elevation, …)
    /// over each coarse cell's children.
    pub fn restrict_mean(&self, fine: &[f64], coarse: &mut [f64]) {
        assert_eq!(fine.len(), self.n_fine());
        assert_eq!(coarse.len(), self.n_coarse());
        coarse.par_iter_mut().enumerate().for_each(|(c, out)| {
            let range = self.child_offsets[c] as usize..self.child_offsets[c + 1] as usize;
            *out = self.children[range.clone()]
                .iter()
                .zip(&self.child_weights[range])
                .map(|(&f, &w)| w * fine[f as usize])
                .sum();
        });
    }

    /// Sum of an *extensive* fine field (water volume, mass, …) over each coarse
    /// cell's children. Conserves the global total exactly (up to round-off).
    pub fn restrict_sum(&self, fine: &[f64], coarse: &mut [f64]) {
        assert_eq!(fine.len(), self.n_fine());
        assert_eq!(coarse.len(), self.n_coarse());
        coarse.par_iter_mut().enumerate().for_each(|(c, out)| {
            *out = self
                .children(c as CellId)
                .iter()
                .map(|&f| fine[f as usize])
                .sum();
        });
    }

    /// Variance of a fine field about its area-weighted mean within each coarse
    /// cell (e.g. sub-grid orography). `mean` must come from
    /// [`TwoLevelMap::restrict_mean`] of the same field.
    pub fn restrict_variance(&self, fine: &[f64], mean: &[f64], coarse: &mut [f64]) {
        assert_eq!(fine.len(), self.n_fine());
        assert_eq!(mean.len(), self.n_coarse());
        assert_eq!(coarse.len(), self.n_coarse());
        coarse.par_iter_mut().enumerate().for_each(|(c, out)| {
            let range = self.child_offsets[c] as usize..self.child_offsets[c + 1] as usize;
            *out = self.children[range.clone()]
                .iter()
                .zip(&self.child_weights[range])
                .map(|(&f, &w)| w * (fine[f as usize] - mean[c]).powi(2))
                .sum();
        });
    }

    /// Barycentric interpolation of a coarse field onto fine cell centres.
    pub fn prolong(&self, coarse: &[f64], fine: &mut [f64]) {
        assert_eq!(coarse.len(), self.n_coarse());
        assert_eq!(fine.len(), self.n_fine());
        fine.par_iter_mut()
            .zip(&self.interpolation)
            .for_each(|(out, (cells, w))| {
                *out = w[0] * coarse[cells[0] as usize]
                    + w[1] * coarse[cells[1] as usize]
                    + w[2] * coarse[cells[2] as usize];
            });
    }
}
