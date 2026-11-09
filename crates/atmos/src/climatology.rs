//! Time-mean zonal means: the climate of a run, summarized by latitude.

use fms_grid::LatitudeBands;
use fms_grid::zonal::east_north;

use crate::layered::{LayeredModel, LayeredState};

/// Time-mean zonal means, per latitude band (south to north).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ZonalMeans {
    /// Band-centre latitudes (degrees).
    pub latitude: Vec<f64>,
    /// Eastward wind per layer, top first (m s⁻¹).
    pub wind: Vec<Vec<f64>>,
    /// Temperature per layer, top first (K).
    pub temperature: Vec<Vec<f64>>,
    /// Surface pressure (Pa).
    pub surface_pressure: Vec<f64>,
}

impl ZonalMeans {
    /// The north–south mirror image: band `i` swapped with band `n − 1 − i`
    /// (latitudes keep their order).
    pub fn mirrored(&self) -> Self {
        let rev = |v: &Vec<f64>| v.iter().rev().copied().collect::<Vec<_>>();
        Self {
            latitude: self.latitude.clone(),
            wind: self.wind.iter().map(rev).collect(),
            temperature: self.temperature.iter().map(rev).collect(),
            surface_pressure: rev(&self.surface_pressure),
        }
    }
}

/// Accumulates [`ZonalMeans`] over samples of a run.
#[derive(Clone, Debug)]
pub struct ZonalClimatology {
    bands: LatitudeBands,
    samples: u64,
    sums: ZonalMeans,
}

impl ZonalClimatology {
    /// Averages over `bands` equal-width latitude bands.
    pub fn new(model: &LayeredModel, bands: usize) -> Self {
        let ops = model.operators();
        let bands = LatitudeBands::new(ops.mesh(), ops.geometry(), bands);
        let n = bands.len();
        let sums = ZonalMeans {
            latitude: bands.latitudes().iter().map(|l| l.to_degrees()).collect(),
            wind: vec![vec![0.0; n]; model.layers()],
            temperature: vec![vec![0.0; n]; model.layers()],
            surface_pressure: vec![0.0; n],
        };
        Self {
            bands,
            samples: 0,
            sums,
        }
    }

    /// Add one sample.
    pub fn sample(&mut self, model: &LayeredModel, state: &LayeredState) {
        let ops = model.operators();
        let (mesh, geometry) = (ops.mesh(), ops.geometry());
        let centers = mesh.cell_centers();
        let add = |sum: &mut Vec<f64>, values: &[f64]| {
            for (s, m) in sum.iter_mut().zip(self.bands.mean(geometry, values)) {
                *s += m;
            }
        };
        for k in 0..model.layers() {
            let vectors = geometry.reconstruct_vectors(mesh, &state.u[k]);
            let east: Vec<f64> = centers
                .iter()
                .zip(&vectors)
                .map(|(&p, &v)| east_north(p, v).0)
                .collect();
            add(&mut self.sums.wind[k], &east);
            let temperature: Vec<f64> = state
                .ps
                .iter()
                .zip(&state.theta_ps[k])
                .map(|(&ps, &t)| t / ps * model.layer_exner(ps)[k])
                .collect();
            add(&mut self.sums.temperature[k], &temperature);
        }
        add(&mut self.sums.surface_pressure, &state.ps);
        self.samples += 1;
    }

    /// Samples taken.
    pub fn samples(&self) -> u64 {
        self.samples
    }

    /// The means so far.
    pub fn means(&self) -> ZonalMeans {
        let n = self.samples.max(1) as f64;
        let scale = |v: &Vec<f64>| v.iter().map(|x| x / n).collect::<Vec<_>>();
        ZonalMeans {
            latitude: self.sums.latitude.clone(),
            wind: self.sums.wind.iter().map(scale).collect(),
            temperature: self.sums.temperature.iter().map(scale).collect(),
            surface_pressure: scale(&self.sums.surface_pressure),
        }
    }
}
