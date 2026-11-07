//! Reading heightmap files: GeoTIFF (e.g. ETOPO 2022) and greyscale PNG.
//!
//! Both must be global equirectangular images (plate carrée), north up. Values
//! are elevations in metres relative to sea level, negative below it.

use std::f64::consts::{FRAC_PI_2, PI, TAU};
use std::fmt;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use tiff::decoder::{ChunkType, Decoder, DecodingResult, Limits};
use tiff::tags::Tag;

use crate::raster::{Geometry, Raster, RasterError, RowSink, block_size};

/// How to read a heightmap.
#[derive(Clone, Copy, Debug, Default)]
pub struct LoadOptions {
    /// For PNG: the elevations (m) of black and white. Required for PNG, since
    /// an image has no units; ignored for GeoTIFF, which stores metres.
    pub png_range: Option<(f64, f64)>,
    /// Average blocks of pixels while reading so the result is no finer than
    /// this spacing (radians). Saves memory on files much finer than the grid.
    /// `None` keeps full resolution.
    pub target_spacing: Option<f64>,
}

/// Why a heightmap could not be read.
#[derive(Debug)]
pub enum HeightmapError {
    /// The file could not be opened or read.
    Io(PathBuf, std::io::Error),
    /// The file is not a valid TIFF or PNG.
    Decode(PathBuf, String),
    /// The file is valid but not usable as a heightmap.
    Unsupported(PathBuf, String),
    /// The pixels don't form a usable global raster.
    Raster(PathBuf, RasterError),
}

impl fmt::Display for HeightmapError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(p, e) => write!(f, "{}: {e}", p.display()),
            Self::Decode(p, e) => write!(f, "{}: cannot decode: {e}", p.display()),
            Self::Unsupported(p, e) => write!(f, "{}: {e}", p.display()),
            Self::Raster(p, e) => write!(f, "{}: {e}", p.display()),
        }
    }
}

impl std::error::Error for HeightmapError {}

/// Read a heightmap, choosing the format from the extension (`.tif`/`.tiff` or
/// `.png`).
pub fn load(path: &Path, options: &LoadOptions) -> Result<Raster, HeightmapError> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(str::to_ascii_lowercase);
    match ext.as_deref() {
        Some("tif" | "tiff") => read_geotiff(path, options),
        Some("png") => read_png(path, options),
        _ => Err(unsupported(
            path,
            "unknown format; expected a .tif, .tiff or .png file",
        )),
    }
}

fn unsupported(path: &Path, why: impl Into<String>) -> HeightmapError {
    HeightmapError::Unsupported(path.to_path_buf(), why.into())
}

fn open(path: &Path) -> Result<BufReader<File>, HeightmapError> {
    File::open(path)
        .map(BufReader::new)
        .map_err(|e| HeightmapError::Io(path.to_path_buf(), e))
}

/// Block size for averaging, given the source spacing.
fn choose_block(width: usize, height: usize, spacing: f64, options: &LoadOptions) -> usize {
    let max_k = options
        .target_spacing
        .map_or(1, |t| (t / spacing).floor() as usize);
    block_size(width, height, max_k)
}

// ---------------------------------------------------------------------------
// GeoTIFF
// ---------------------------------------------------------------------------

/// GeoTIFF key: model type (1 projected, 2 geographic, 3 geocentric).
const GT_MODEL_TYPE: u16 = 1024;
/// GeoTIFF key: raster type (1 pixel is area, 2 pixel is point).
const GT_RASTER_TYPE: u16 = 1025;
const MODEL_GEOGRAPHIC: u16 = 2;
const RASTER_PIXEL_IS_POINT: u16 = 2;

fn read_geotiff(path: &Path, options: &LoadOptions) -> Result<Raster, HeightmapError> {
    let decode = |e: tiff::TiffError| HeightmapError::Decode(path.to_path_buf(), e.to_string());
    let mut dec = Decoder::new(open(path)?)
        .map_err(decode)?
        .with_limits(Limits::unlimited());
    let (w, h) = dec.dimensions().map_err(decode)?;
    let (w, h) = (w as usize, h as usize);
    match dec.colortype().map_err(decode)? {
        tiff::ColorType::Gray(_) => {}
        other => {
            return Err(unsupported(
                path,
                format!("expected one band of elevations, found {other:?}"),
            ));
        }
    }

    let scale = dec
        .find_tag(Tag::ModelPixelScaleTag)
        .map_err(decode)?
        .map(|v| v.into_f64_vec())
        .transpose()
        .map_err(decode)?;
    let tiepoint = dec
        .find_tag(Tag::ModelTiepointTag)
        .map_err(decode)?
        .map(|v| v.into_f64_vec())
        .transpose()
        .map_err(decode)?;
    let keys = dec
        .find_tag(Tag::GeoKeyDirectoryTag)
        .map_err(decode)?
        .map(|v| v.into_u16_vec())
        .transpose()
        .map_err(decode)?
        .unwrap_or_default();
    let nodata: Option<f64> = dec
        .find_tag(Tag::GdalNodata)
        .map_err(decode)?
        .and_then(|v| v.into_string().ok())
        .and_then(|s| s.trim_matches(char::from(0)).trim().parse().ok());

    if geo_key(&keys, GT_MODEL_TYPE).is_some_and(|t| t != MODEL_GEOGRAPHIC) {
        return Err(unsupported(
            path,
            "the raster is in a projected coordinate system; \
             only latitude/longitude (plate carrée) is supported",
        ));
    }
    let pixel_is_point = geo_key(&keys, GT_RASTER_TYPE) == Some(RASTER_PIXEL_IS_POINT);

    let geometry = match (scale, tiepoint) {
        (Some(s), Some(t)) if s.len() >= 2 && t.len() >= 6 => {
            let (sx, sy) = (s[0], s[1]);
            let (i, j, x, y) = (t[0], t[1], t[3], t[4]);
            // Pixel-is-area: the tiepoint is a pixel's corner; is-point: its centre.
            let half = if pixel_is_point { 0.0 } else { 0.5 };
            Geometry {
                lon0: (x + (half - i) * sx).to_radians(),
                lat0: (y - (half - j) * sy).to_radians(),
                dlon: sx.to_radians(),
                dlat: sy.to_radians(),
            }
        }
        // A plain TIFF: assume it tiles the globe like an image.
        _ => Geometry {
            lon0: -PI + 0.5 * TAU / w as f64,
            lat0: FRAC_PI_2 - 0.5 * PI / h as f64,
            dlon: TAU / w as f64,
            dlat: PI / h as f64,
        },
    };

    // Grid-registered rasters repeat the first column at +180°; drop it.
    let wraps_twice = ((w as f64 - 1.0) * geometry.dlon - TAU).abs() < 0.01 * geometry.dlon;
    let keep = if wraps_twice { w - 1 } else { w };
    let k = choose_block(keep, h, geometry.dlon, options);
    let mut sink = RowSink::new(keep, h, geometry, k);

    let to_f32 = |v: f64| match nodata {
        Some(nd) if v == nd => f32::NAN,
        _ => v as f32,
    };
    let n_chunks = match dec.get_chunk_type() {
        ChunkType::Strip => dec.strip_count().map_err(decode)?,
        ChunkType::Tile => dec.tile_count().map_err(decode)?,
    };
    match dec.get_chunk_type() {
        ChunkType::Strip => {
            for i in 0..n_chunks {
                let (cw, ch) = dec.chunk_data_dimensions(i);
                let values = samples(dec.read_chunk(i).map_err(decode)?, to_f32);
                for r in 0..ch as usize {
                    sink.push(&values[r * cw as usize..(r + 1) * cw as usize]);
                }
            }
        }
        ChunkType::Tile => {
            let (tw, th) = dec.chunk_dimensions();
            let (tw, th) = (tw as usize, th as usize);
            let across = w.div_ceil(tw);
            let mut band = vec![f32::NAN; w * th];
            for ty in 0..h.div_ceil(th) {
                let mut rows = 0;
                for tx in 0..across {
                    let index = (ty * across + tx) as u32;
                    let (cw, ch) = dec.chunk_data_dimensions(index);
                    let (cw, ch) = (cw as usize, ch as usize);
                    let values = samples(dec.read_chunk(index).map_err(decode)?, to_f32);
                    for r in 0..ch {
                        band[r * w + tx * tw..r * w + tx * tw + cw]
                            .copy_from_slice(&values[r * cw..(r + 1) * cw]);
                    }
                    rows = ch;
                }
                for r in 0..rows {
                    sink.push(&band[r * w..(r + 1) * w]);
                }
            }
        }
    }
    sink.finish()
        .map_err(|e| HeightmapError::Raster(path.to_path_buf(), e))
}

/// Look up a short-valued key in a GeoKeyDirectory: a 4-value header, then
/// `[key, location, count, value]` entries. Location 0 means the value is inline.
fn geo_key(directory: &[u16], key: u16) -> Option<u16> {
    let n = *directory.get(3)? as usize;
    directory
        .get(4..4 + 4 * n)?
        .as_chunks::<4>()
        .0
        .iter()
        .find(|e| e[0] == key && e[1] == 0)
        .map(|e| e[3])
}

fn samples(result: DecodingResult, f: impl Fn(f64) -> f32) -> Vec<f32> {
    fn map<T: Copy>(v: Vec<T>, g: impl Fn(T) -> f64, f: impl Fn(f64) -> f32) -> Vec<f32> {
        v.into_iter().map(|x| f(g(x))).collect()
    }
    match result {
        DecodingResult::U8(v) => map(v, f64::from, f),
        DecodingResult::U16(v) => map(v, f64::from, f),
        DecodingResult::U32(v) => map(v, f64::from, f),
        DecodingResult::U64(v) => map(v, |x| x as f64, f),
        DecodingResult::I8(v) => map(v, f64::from, f),
        DecodingResult::I16(v) => map(v, f64::from, f),
        DecodingResult::I32(v) => map(v, f64::from, f),
        DecodingResult::I64(v) => map(v, |x| x as f64, f),
        DecodingResult::F16(v) => map(v, f64::from, f),
        DecodingResult::F32(v) => map(v, f64::from, f),
        DecodingResult::F64(v) => map(v, |x| x, f),
    }
}

// ---------------------------------------------------------------------------
// PNG
// ---------------------------------------------------------------------------

fn read_png(path: &Path, options: &LoadOptions) -> Result<Raster, HeightmapError> {
    let Some((min, max)) = options.png_range else {
        return Err(unsupported(
            path,
            "a PNG heightmap needs the elevations of black and white (e.g. -8000,6000 m)",
        ));
    };
    let decode = |e: png::DecodingError| HeightmapError::Decode(path.to_path_buf(), e.to_string());
    let mut decoder = png::Decoder::new(open(path)?);
    // Expand low bit depths and palettes, but keep 16 bits.
    decoder.set_transformations(png::Transformations::EXPAND);
    let mut reader = decoder.read_info().map_err(decode)?;
    let (w, h) = {
        let info = reader.info();
        (info.width as usize, info.height as usize)
    };
    let (color, depth) = reader.output_color_type();
    if color != png::ColorType::Grayscale {
        return Err(unsupported(
            path,
            format!("expected a greyscale image, found {color:?}"),
        ));
    }
    let (bytes_per_sample, full_scale) = match depth {
        png::BitDepth::Sixteen => (2, 65535.0),
        _ => (1, 255.0),
    };
    let to_metres = |v: f64| (min + (max - min) * v / full_scale) as f32;

    let geometry = Geometry {
        lon0: -PI + 0.5 * TAU / w as f64,
        lat0: FRAC_PI_2 - 0.5 * PI / h as f64,
        dlon: TAU / w as f64,
        dlat: PI / h as f64,
    };
    let k = choose_block(w, h, geometry.dlon, options);
    let mut sink = RowSink::new(w, h, geometry, k);
    let mut row = vec![0.0f32; w];
    while let Some(raw) = reader.next_row().map_err(decode)? {
        let data = raw.data();
        for (out, s) in row.iter_mut().zip(data.chunks_exact(bytes_per_sample)) {
            let v = if bytes_per_sample == 2 {
                u16::from_be_bytes([s[0], s[1]]) as f64
            } else {
                s[0] as f64
            };
            *out = to_metres(v);
        }
        sink.push(&row);
    }
    sink.finish()
        .map_err(|e| HeightmapError::Raster(path.to_path_buf(), e))
}
