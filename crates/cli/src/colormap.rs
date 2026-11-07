//! Colour scales as 256-entry palettes, so the same indexed image works for PNG
//! and GIF.

use fms_fields::ColorScale;

/// Number of palette entries used for data values.
pub const LEVELS: usize = 253;
/// Palette index for missing / non-finite values.
pub const MISSING: u8 = 253;
/// Palette index for white (text, background).
pub const WHITE: u8 = 254;
/// Palette index for black.
pub const BLACK: u8 = 255;

/// A colour scale applied to a value range.
pub struct Colormap {
    scale: ColorScale,
    min: f64,
    max: f64,
    palette: [[u8; 3]; 256],
}

type Stops = &'static [(f64, [u8; 3])];

const VIRIDIS: Stops = &[
    (0.0, [68, 1, 84]),
    (0.25, [59, 82, 139]),
    (0.5, [33, 145, 140]),
    (0.75, [94, 201, 98]),
    (1.0, [253, 231, 37]),
];
const DIVERGING: Stops = &[
    (0.0, [33, 102, 172]),
    (0.25, [103, 169, 207]),
    (0.5, [247, 247, 247]),
    (0.75, [239, 138, 98]),
    (1.0, [178, 24, 43]),
];
const TEMPERATURE: Stops = &[
    (0.0, [48, 18, 59]),
    (0.15, [70, 107, 227]),
    (0.35, [40, 187, 236]),
    (0.5, [50, 241, 152]),
    (0.65, [164, 252, 60]),
    (0.8, [251, 185, 56]),
    (0.92, [228, 70, 10]),
    (1.0, [122, 4, 3]),
];
const MOISTURE: Stops = &[
    (0.0, [140, 81, 10]),
    (0.25, [216, 179, 101]),
    (0.5, [246, 232, 195]),
    (0.75, [90, 180, 172]),
    (1.0, [1, 102, 94]),
];
/// Hypsometric tint, ocean part: deep to shallow.
const OCEAN: Stops = &[
    (0.0, [8, 29, 88]),
    (0.7, [34, 94, 168]),
    (1.0, [161, 218, 180]),
];
/// Hypsometric tint, land part: lowland to peaks.
const LAND: Stops = &[
    (0.0, [26, 120, 60]),
    (0.3, [166, 189, 96]),
    (0.6, [140, 100, 60]),
    (0.86, [200, 190, 180]),
    (1.0, [255, 255, 255]),
];
/// Palette entries `0..OCEAN_LEVELS` are ocean, the rest land, so sea level is a
/// hard colour boundary.
const OCEAN_LEVELS: usize = LEVELS / 2;

const CATEGORICAL: &[[u8; 3]] = &[
    [31, 119, 180],
    [255, 127, 14],
    [44, 160, 44],
    [214, 39, 40],
    [148, 103, 189],
    [140, 86, 75],
    [227, 119, 194],
    [127, 127, 127],
    [188, 189, 34],
    [23, 190, 207],
    [174, 199, 232],
    [255, 187, 120],
    [152, 223, 138],
    [255, 152, 150],
    [197, 176, 213],
    [196, 156, 148],
    [247, 182, 210],
    [199, 199, 199],
    [219, 219, 141],
    [158, 218, 229],
];

fn interpolate(stops: Stops, t: f64) -> [u8; 3] {
    let t = t.clamp(0.0, 1.0);
    let i = stops
        .iter()
        .position(|&(s, _)| s >= t)
        .unwrap_or(stops.len() - 1)
        .max(1);
    let (t0, c0) = stops[i - 1];
    let (t1, c1) = stops[i];
    let f = if t1 > t0 { (t - t0) / (t1 - t0) } else { 0.0 };
    [0, 1, 2].map(|k| (f64::from(c0[k]) + f * (f64::from(c1[k]) - f64::from(c0[k]))).round() as u8)
}

impl Colormap {
    /// A colour scale over `[min, max]`.
    pub fn new(scale: ColorScale, min: f64, max: f64) -> Self {
        let mut palette = [[0u8; 3]; 256];
        for (i, entry) in palette.iter_mut().take(LEVELS).enumerate() {
            let t = i as f64 / (LEVELS - 1) as f64;
            *entry = match scale {
                ColorScale::Sequential => interpolate(VIRIDIS, t),
                ColorScale::Diverging => interpolate(DIVERGING, t),
                ColorScale::Temperature => interpolate(TEMPERATURE, t),
                ColorScale::Moisture => interpolate(MOISTURE, t),
                ColorScale::Elevation if i < OCEAN_LEVELS => {
                    interpolate(OCEAN, i as f64 / (OCEAN_LEVELS - 1) as f64)
                }
                ColorScale::Elevation => interpolate(
                    LAND,
                    (i - OCEAN_LEVELS) as f64 / (LEVELS - OCEAN_LEVELS - 1) as f64,
                ),
                ColorScale::Categorical => CATEGORICAL[i % CATEGORICAL.len()],
            };
        }
        palette[MISSING as usize] = [255, 0, 255];
        palette[WHITE as usize] = [255, 255, 255];
        palette[BLACK as usize] = [0, 0, 0];
        Self {
            scale,
            min,
            max,
            palette,
        }
    }

    /// The 256-colour palette.
    pub fn palette(&self) -> &[[u8; 3]; 256] {
        &self.palette
    }

    /// Value range `(min, max)`.
    pub fn range(&self) -> (f64, f64) {
        (self.min, self.max)
    }

    /// Palette index of a value.
    pub fn index(&self, v: f64) -> u8 {
        if !v.is_finite() {
            return MISSING;
        }
        let t = match self.scale {
            // Integer classes cycle through the palette.
            ColorScale::Categorical => {
                return (v.round().rem_euclid(CATEGORICAL.len() as f64)) as u8;
            }
            // Below sea level uses the ocean half, at or above it the land half.
            ColorScale::Elevation => return self.elevation_index(v),
            _ if self.max > self.min => (v - self.min) / (self.max - self.min),
            _ => 0.5,
        };
        (t.clamp(0.0, 1.0) * (LEVELS - 1) as f64).round() as u8
    }

    fn elevation_index(&self, v: f64) -> u8 {
        let land_levels = LEVELS - OCEAN_LEVELS;
        if v < 0.0 {
            // Deepest at `min` (or −1 m if the range has no ocean).
            let depth = if self.min < 0.0 { self.min } else { -1.0 };
            let t = (1.0 - v / depth).clamp(0.0, 1.0);
            ((t * OCEAN_LEVELS as f64) as usize).min(OCEAN_LEVELS - 1) as u8
        } else {
            let height = if self.max > 0.0 { self.max } else { 1.0 };
            let t = (v / height).clamp(0.0, 1.0);
            (OCEAN_LEVELS + (t * (land_levels - 1) as f64).round() as usize) as u8
        }
    }

    /// Palette index for position `t ∈ [0, 1]` along the colour bar.
    pub fn bar_index(&self, t: f64) -> u8 {
        match self.scale {
            ColorScale::Elevation => self.elevation_index(self.min + t * (self.max - self.min)),
            ColorScale::Categorical => {
                ((t * CATEGORICAL.len() as f64) as usize).min(CATEGORICAL.len() - 1) as u8
            }
            _ => (t.clamp(0.0, 1.0) * (LEVELS - 1) as f64).round() as u8,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn indices_span_the_palette_and_clamp() {
        let c = Colormap::new(ColorScale::Sequential, 10.0, 20.0);
        assert_eq!(c.index(10.0), 0);
        assert_eq!(c.index(20.0), (LEVELS - 1) as u8);
        assert_eq!(c.index(-1e9), 0);
        assert_eq!(c.index(1e9), (LEVELS - 1) as u8);
        assert_eq!(c.index(f64::NAN), MISSING);
        assert_eq!(c.palette()[0], [68, 1, 84]);
    }

    #[test]
    fn elevation_puts_sea_level_at_the_coast() {
        let c = Colormap::new(ColorScale::Elevation, -8000.0, 4000.0);
        let coast = OCEAN_LEVELS as u8;
        assert_eq!(c.index(0.0), coast);
        assert_eq!(c.index(-1.0), coast - 1);
        assert_eq!(c.index(-1e-9), coast - 1);
        // The two sides really are different colours (ocean is bluish).
        let [r, _, b] = c.palette()[coast as usize - 1];
        assert!(b > r);
        let [r, g, b] = c.palette()[coast as usize];
        assert!(g > b && g > r);
        assert_eq!(c.index(-8000.0), 0);
        assert_eq!(c.index(4000.0), (LEVELS - 1) as u8);
    }

    #[test]
    fn categorical_cycles() {
        let c = Colormap::new(ColorScale::Categorical, 0.0, 30.0);
        assert_eq!(c.index(3.0), 3);
        assert_eq!(c.index(23.0), 3);
    }

    #[test]
    fn degenerate_range_is_safe() {
        let c = Colormap::new(ColorScale::Diverging, 5.0, 5.0);
        assert_eq!(c.index(5.0), ((LEVELS - 1) as f64 * 0.5).round() as u8);
    }
}
