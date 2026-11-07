//! Global equirectangular rasters: one sample per latitude/longitude pixel.

use std::f64::consts::{FRAC_PI_2, PI, TAU};
use std::fmt;

use fms_geo::coords::to_lat_lon;
use glam::DVec3;

/// Relative tolerance (in pixels) when checking that a raster spans the globe.
const COVERAGE_TOLERANCE: f64 = 0.01;

/// A global equirectangular raster.
///
/// Rows run north to south and columns west to east. Each sample is the value at
/// its pixel centre, and stands for the pixel's area. The columns cover all 360°
/// of longitude exactly once; the rows reach both poles.
#[derive(Clone, Debug)]
pub struct Raster {
    width: usize,
    height: usize,
    /// Row-major; row 0 is the northernmost.
    data: Vec<f32>,
    geometry: Geometry,
}

/// Where the pixels are: centre of pixel (0, 0) and the spacing, in radians.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Geometry {
    /// Longitude of column 0's centre.
    pub lon0: f64,
    /// Latitude of row 0's centre.
    pub lat0: f64,
    /// Column spacing (eastward). `width · dlon = 2π`.
    pub dlon: f64,
    /// Row spacing (southward).
    pub dlat: f64,
}

/// Why a raster could not be built.
#[derive(Clone, Debug, PartialEq)]
pub enum RasterError {
    /// The data length doesn't match `width · height`.
    WrongLength {
        /// `width · height`.
        expected: usize,
        /// The data supplied.
        actual: usize,
    },
    /// The raster doesn't span the whole globe.
    NotGlobal(String),
    /// Some pixels have no value (NaN or the file's no-data marker).
    MissingData {
        /// How many.
        count: usize,
    },
}

impl fmt::Display for RasterError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongLength { expected, actual } => {
                write!(f, "expected {expected} samples, got {actual}")
            }
            Self::NotGlobal(why) => write!(f, "not a global raster: {why}"),
            Self::MissingData { count } => {
                write!(f, "{count} pixels have no data; fill them before importing")
            }
        }
    }
}

impl std::error::Error for RasterError {}

impl Raster {
    /// A global raster whose pixels tile the sphere edge to edge: column 0 starts
    /// at 180° W and row 0 at the north pole. This is how a plain image
    /// (e.g. a PNG heightmap) is read.
    pub fn global(width: usize, height: usize, data: Vec<f32>) -> Result<Self, RasterError> {
        let dlon = TAU / width as f64;
        let dlat = PI / height as f64;
        Self::with_geometry(
            width,
            height,
            data,
            Geometry {
                lon0: -PI + 0.5 * dlon,
                lat0: FRAC_PI_2 - 0.5 * dlat,
                dlon,
                dlat,
            },
        )
    }

    /// A global raster sampled from `f(latitude, longitude)` (radians) at each
    /// pixel centre, laid out as in [`Raster::global`]. For synthetic heightmaps.
    pub fn from_fn(width: usize, height: usize, f: impl Fn(f64, f64) -> f64) -> Self {
        let mut raster =
            Self::global(width, height, vec![0.0; width * height]).expect("length matches");
        for row in 0..height {
            let lat = raster.lat(row);
            for col in 0..width {
                raster.data[row * width + col] = f(lat, raster.lon(col)) as f32;
            }
        }
        raster
    }

    /// Build with explicit geometry, checking that it covers the globe and that
    /// every pixel has a value.
    pub(crate) fn with_geometry(
        width: usize,
        height: usize,
        data: Vec<f32>,
        geometry: Geometry,
    ) -> Result<Self, RasterError> {
        if data.len() != width * height {
            return Err(RasterError::WrongLength {
                expected: width * height,
                actual: data.len(),
            });
        }
        check_global(width, height, geometry)?;
        let missing = data.iter().filter(|v| v.is_nan()).count();
        if missing > 0 {
            return Err(RasterError::MissingData { count: missing });
        }
        Ok(Self {
            width,
            height,
            data,
            geometry,
        })
    }

    /// Number of columns.
    pub fn width(&self) -> usize {
        self.width
    }

    /// Number of rows.
    pub fn height(&self) -> usize {
        self.height
    }

    /// Samples, row-major, row 0 northernmost.
    pub fn data(&self) -> &[f32] {
        &self.data
    }

    /// Value of one pixel.
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.data[row * self.width + col] as f64
    }

    /// Column spacing (radians).
    pub fn spacing(&self) -> f64 {
        self.geometry.dlon
    }

    /// Latitude of a row's pixel centres (radians).
    pub fn lat(&self, row: usize) -> f64 {
        self.geometry.lat0 - row as f64 * self.geometry.dlat
    }

    /// Longitude of a column's pixel centres (radians), in `[-π, π)` or a little
    /// outside it, depending on the source's origin.
    pub fn lon(&self, col: usize) -> f64 {
        self.geometry.lon0 + col as f64 * self.geometry.dlon
    }

    /// Area of one pixel in `row` divided by the column spacing:
    /// `sin(north edge) − sin(south edge)`, with the edges clipped to the poles.
    /// Pixels in a row all have the same area, so this is the row's weight in
    /// area-weighted means.
    pub fn row_weight(&self, row: usize) -> f64 {
        self.geometry.row_weight(row)
    }

    /// Fractional row of a latitude: 0 at row 0's centre, increasing southward.
    pub(crate) fn row_of(&self, lat: f64) -> f64 {
        (self.geometry.lat0 - lat) / self.geometry.dlat
    }

    /// Fractional column of a longitude, in `[0, width)`.
    pub(crate) fn col_of(&self, lon: f64) -> f64 {
        let x = (lon - self.geometry.lon0) / self.geometry.dlon;
        x.rem_euclid(self.width as f64)
    }

    /// Bilinear interpolation at a point on the unit sphere. Wraps in longitude;
    /// clamps to the first and last rows near the poles.
    pub fn sample(&self, p: DVec3) -> f64 {
        let (lat, lon) = to_lat_lon(p);
        let y = self.row_of(lat).clamp(0.0, (self.height - 1) as f64);
        let x = self.col_of(lon);
        let (r0, c0) = (y.floor() as usize, x.floor() as usize);
        let r1 = (r0 + 1).min(self.height - 1);
        let c1 = (c0 + 1) % self.width;
        let (fy, fx) = (y - r0 as f64, x - c0 as f64);
        let top = (1.0 - fx) * self.get(r0, c0) + fx * self.get(r0, c1);
        let bottom = (1.0 - fx) * self.get(r1, c0) + fx * self.get(r1, c1);
        (1.0 - fy) * top + fy * bottom
    }
}

impl Geometry {
    pub(crate) fn row_weight(&self, row: usize) -> f64 {
        let lat = self.lat0 - row as f64 * self.dlat;
        let north = (lat + 0.5 * self.dlat).min(FRAC_PI_2);
        let south = (lat - 0.5 * self.dlat).max(-FRAC_PI_2);
        (fms_math::sin(north) - fms_math::sin(south)).max(0.0)
    }

    /// The geometry after averaging `k × k` blocks of pixels.
    pub(crate) fn decimated(&self, k: usize) -> Self {
        let offset = 0.5 * (k as f64 - 1.0);
        Self {
            lon0: self.lon0 + offset * self.dlon,
            lat0: self.lat0 - offset * self.dlat,
            dlon: self.dlon * k as f64,
            dlat: self.dlat * k as f64,
        }
    }
}

fn check_global(width: usize, height: usize, g: Geometry) -> Result<(), RasterError> {
    if width == 0 || height == 0 || !(g.dlon > 0.0 && g.dlat > 0.0) {
        return Err(RasterError::NotGlobal(format!(
            "{width}×{height} pixels with spacing {:.4}° × {:.4}°",
            g.dlon.to_degrees(),
            g.dlat.to_degrees()
        )));
    }
    let span = width as f64 * g.dlon;
    if (span - TAU).abs() > COVERAGE_TOLERANCE * g.dlon {
        return Err(RasterError::NotGlobal(format!(
            "columns span {:.3}° of longitude, not 360°",
            span.to_degrees()
        )));
    }
    // The rows must reach both poles: the outer pixels' far edges at or past them.
    let north_edge = g.lat0 + 0.5 * g.dlat;
    let south_edge = g.lat0 - (height as f64 - 0.5) * g.dlat;
    let tolerance = COVERAGE_TOLERANCE * g.dlat;
    if north_edge < FRAC_PI_2 - tolerance || south_edge > -FRAC_PI_2 + tolerance {
        return Err(RasterError::NotGlobal(format!(
            "rows span {:.3}° to {:.3}° of latitude, not pole to pole",
            south_edge.to_degrees(),
            north_edge.to_degrees()
        )));
    }
    Ok(())
}

/// Builds a raster from rows as they're decoded, averaging `k × k` blocks of
/// pixels on the way (weighted by pixel area). Lets a large file be read at
/// reduced resolution without holding it all in memory.
pub(crate) struct RowSink {
    /// Columns to keep from each source row (drops a duplicated wrap column).
    src_width: usize,
    src_geometry: Geometry,
    k: usize,
    src_rows_seen: usize,
    out_width: usize,
    out: Vec<f32>,
    /// Running sums for the current block row: `Σ w·z` and `Σ w`.
    sums: Vec<f64>,
    weights: Vec<f64>,
}

impl RowSink {
    /// `src_width` columns are kept from each pushed row; `k` must divide it and
    /// the number of rows.
    pub(crate) fn new(src_width: usize, src_height: usize, geometry: Geometry, k: usize) -> Self {
        assert!(k >= 1 && src_width.is_multiple_of(k) && src_height.is_multiple_of(k));
        let out_width = src_width / k;
        Self {
            src_width,
            src_geometry: geometry,
            k,
            src_rows_seen: 0,
            out_width,
            out: Vec::with_capacity(out_width * (src_height / k)),
            sums: vec![0.0; out_width],
            weights: vec![0.0; out_width],
        }
    }

    /// Add the next source row. NaN marks a missing value; a block with no valid
    /// pixel stays NaN.
    pub(crate) fn push(&mut self, row: &[f32]) {
        let row = &row[..self.src_width];
        if self.k == 1 {
            self.out.extend_from_slice(row);
        } else {
            let w = self.src_geometry.row_weight(self.src_rows_seen);
            for (block, chunk) in row.chunks_exact(self.k).enumerate() {
                for &z in chunk {
                    if !z.is_nan() {
                        self.sums[block] += w * z as f64;
                        self.weights[block] += w;
                    }
                }
            }
            if (self.src_rows_seen + 1).is_multiple_of(self.k) {
                for (s, w) in self.sums.iter_mut().zip(&mut self.weights) {
                    self.out
                        .push(if *w > 0.0 { (*s / *w) as f32 } else { f32::NAN });
                    (*s, *w) = (0.0, 0.0);
                }
            }
        }
        self.src_rows_seen += 1;
    }

    /// The finished raster.
    pub(crate) fn finish(self) -> Result<Raster, RasterError> {
        let height = self.out.len() / self.out_width;
        Raster::with_geometry(
            self.out_width,
            height,
            self.out,
            self.src_geometry.decimated(self.k),
        )
    }
}

/// Largest block size `k ≤ max_k` that divides both dimensions (at least 1).
pub(crate) fn block_size(width: usize, height: usize, max_k: usize) -> usize {
    (1..=max_k.max(1))
        .rev()
        .find(|k| width.is_multiple_of(*k) && height.is_multiple_of(*k))
        .unwrap_or(1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn row_weights_sum_to_two() {
        // ∫ cos φ dφ over [-π/2, π/2] = 2.
        for (w, h) in [(8, 4), (360, 180), (7, 5)] {
            let r = Raster::global(w, h, vec![0.0; w * h]).unwrap();
            let total: f64 = (0..h).map(|row| r.row_weight(row)).sum();
            assert!((total - 2.0).abs() < 1e-12, "{w}×{h}: {total}");
        }
    }

    #[test]
    fn rejects_regional_rasters() {
        let g = Geometry {
            lon0: 0.0,
            lat0: 0.0,
            dlon: 0.01,
            dlat: 0.01,
        };
        assert!(matches!(
            Raster::with_geometry(10, 10, vec![0.0; 100], g),
            Err(RasterError::NotGlobal(_))
        ));
        assert!(matches!(
            Raster::global(4, 2, vec![0.0; 7]),
            Err(RasterError::WrongLength { .. })
        ));
        let mut data = vec![0.0; 8];
        data[3] = f32::NAN;
        assert_eq!(
            Raster::global(4, 2, data).unwrap_err(),
            RasterError::MissingData { count: 1 }
        );
    }

    #[test]
    fn decimation_is_an_area_weighted_block_mean() {
        let fine = Raster::from_fn(12, 6, |lat, lon| 100.0 * lat + 10.0 * lon);
        let mut sink = RowSink::new(12, 6, fine.geometry, 3);
        for row in 0..6 {
            sink.push(&fine.data[row * 12..(row + 1) * 12]);
        }
        let coarse = sink.finish().unwrap();
        assert_eq!((coarse.width(), coarse.height()), (4, 2));
        // The coarse pixel grid is exactly the 4×2 global layout.
        let expected = Raster::global(4, 2, vec![0.0; 8]).unwrap();
        assert!((coarse.lat(0) - expected.lat(0)).abs() < 1e-12);
        assert!((coarse.lon(0) - expected.lon(0)).abs() < 1e-12);
        for row in 0..2 {
            for col in 0..4 {
                let (mut s, mut w) = (0.0, 0.0);
                for r in 3 * row..3 * row + 3 {
                    for c in 3 * col..3 * col + 3 {
                        s += fine.row_weight(r) * fine.get(r, c);
                        w += fine.row_weight(r);
                    }
                }
                assert!((coarse.get(row, col) - s / w).abs() < 1e-4);
            }
        }
    }

    #[test]
    fn bilinear_sampling_wraps_in_longitude() {
        let r = Raster::from_fn(36, 18, |_, lon| fms_math::cos(lon));
        for lon in [-PI, -3.1, 3.1, PI - 1e-9] {
            let p = fms_geo::coords::from_lat_lon(0.1, lon);
            assert!((r.sample(p) - fms_math::cos(lon)).abs() < 0.01, "{lon}");
        }
    }

    #[test]
    fn block_size_divides_both_dimensions() {
        assert_eq!(block_size(21600, 10800, 7), 6);
        assert_eq!(block_size(21600, 10800, 0), 1);
        assert_eq!(block_size(21601, 10801, 4), 1);
    }
}
