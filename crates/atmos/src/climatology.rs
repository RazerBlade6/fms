//! Time-mean zonal means: the climate of a run, summarized by latitude, and
//! the meridional overturning ([`Overturning`]) from the northward mass flux
//! across the band boundaries.

use fms_grid::LatitudeBands;
use fms_grid::zonal::east_north;

use crate::circulation::Overturning;
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

/// An edge on one or more discrete latitude circles: it separates a cell in
/// a band south of each boundary `first..=last` from one north of it.
#[derive(Clone, Debug)]
struct Crossing {
    edge: usize,
    cells: [usize; 2],
    /// `±lₑ a / g`: the edge's length (m) over gravity, signed `+` when the
    /// edge normal (from `cells[0]` to `cells[1]`) points north across the
    /// boundaries.
    weight: f64,
    first: usize,
    last: usize,
}

/// Accumulates [`ZonalMeans`] and the northward mass flux across the band
/// boundaries ([`Overturning`]) over samples of a run.
#[derive(Clone, Debug)]
pub struct ZonalClimatology {
    bands: LatitudeBands,
    samples: u64,
    sums: ZonalMeans,
    /// The edges on each band boundary's discrete latitude circle, in edge
    /// order.
    crossings: Vec<Crossing>,
    /// `δσₖ` per layer.
    dsigma: Vec<f64>,
    /// σ at the bottom of each layer.
    sigma: Vec<f64>,
    /// Summed northward mass flux (kg s⁻¹) per layer, top first, and band
    /// boundary, south pole to north pole.
    flux: Vec<Vec<f64>>,
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
        // An edge between bands b0 < b1 lies on the latitude circles of
        // boundaries b0 + 1 ..= b1; its flux is northward there when its
        // normal points from the southern cell to the northern one.
        let band = bands.band_of_cell();
        let scale = ops.radius() / model.gravity();
        let lengths = &ops.geometry().edge_lengths;
        let crossings = ops
            .mesh()
            .edges()
            .iter()
            .enumerate()
            .filter_map(|(e, edge)| {
                let cells = edge.cells.map(|c| c as usize);
                let (b0, b1) = (band[cells[0]] as usize, band[cells[1]] as usize);
                let (sign, first, last) = match b0.cmp(&b1) {
                    std::cmp::Ordering::Equal => return None,
                    std::cmp::Ordering::Less => (1.0, b0 + 1, b1),
                    std::cmp::Ordering::Greater => (-1.0, b1 + 1, b0),
                };
                Some(Crossing {
                    edge: e,
                    cells,
                    weight: sign * lengths[e] * scale,
                    first,
                    last,
                })
            })
            .collect();
        let levels = model.levels();
        Self {
            bands,
            samples: 0,
            sums,
            crossings,
            dsigma: levels.dsigma().to_vec(),
            sigma: levels.interfaces()[1..].to_vec(),
            flux: vec![vec![0.0; n + 1]; model.layers()],
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
        // The model's own layer mass flux, δσₖ p̄ₛ uₑ (Pa m s⁻¹), times
        // lₑ a / g: kg s⁻¹ across each edge.
        for (k, flux) in self.flux.iter_mut().enumerate() {
            let (ds, u) = (self.dsigma[k], &state.u[k]);
            for c in &self.crossings {
                let ps_edge = 0.5 * (state.ps[c.cells[0]] + state.ps[c.cells[1]]);
                let f = ds * ps_edge * u[c.edge] * c.weight;
                for sum in &mut flux[c.first..=c.last] {
                    *sum += f;
                }
            }
        }
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

    /// The time-mean overturning so far: the northward mass flux above each
    /// layer's bottom, summed from the top down.
    pub fn overturning(&self) -> Overturning {
        let n = self.samples.max(1) as f64;
        let boundaries = self.bands.len() + 1;
        let width = 180.0 / self.bands.len() as f64;
        let mut above = vec![0.0; boundaries];
        let streamfunction = self
            .flux
            .iter()
            .map(|layer| {
                for (a, f) in above.iter_mut().zip(layer) {
                    *a += f / n;
                }
                above.clone()
            })
            .collect();
        Overturning {
            latitude: (0..boundaries).map(|b| -90.0 + b as f64 * width).collect(),
            sigma: self.sigma.clone(),
            streamfunction,
        }
    }
}
