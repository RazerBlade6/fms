//! Equirectangular rendering of cell fields to indexed images, PNG and GIF.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use fms_geo::coords::from_lat_lon;
use fms_sim::Grid;
use rayon::prelude::*;

use glam::DVec3;

use crate::arrows::ArrowLayout;
use crate::colormap::{BLACK, Colormap, WHITE};
use crate::font::{draw_text, text_width};

/// Height of the legend strip below the map (pixels).
const FOOTER: usize = 44;
const TEXT_SCALE: usize = 2;

/// For each pixel of a `width × width/2` equirectangular map, which grid cells
/// it shows and with what weights. Computed once and reused for every frame.
pub struct CellRaster {
    /// Image width.
    pub width: usize,
    /// Map height (excluding the legend).
    pub height: usize,
    stencils: Vec<([u32; 3], [f32; 3])>,
}

impl CellRaster {
    /// Map every pixel centre to grid cells.
    ///
    /// With `smooth = false` each pixel shows the single cell containing it,
    /// which shows the true resolution. With `smooth = true` it blends the three
    /// cells around it (barycentric on the Delaunay triangle), which looks
    /// better but invents detail the model does not have.
    pub fn new(grid: &Grid, width: usize, smooth: bool) -> Self {
        let height = width / 2;
        let stencils = (0..width * height)
            .into_par_iter()
            .map(|i| {
                let (x, y) = (i % width, i / width);
                let lat = std::f64::consts::FRAC_PI_2
                    - std::f64::consts::PI * (y as f64 + 0.5) / height as f64;
                let lon =
                    -std::f64::consts::PI + std::f64::consts::TAU * (x as f64 + 0.5) / width as f64;
                let p = from_lat_lon(lat, lon);
                if smooth {
                    let (t, w) = grid.locator.locate_triangle(&grid.mesh, p);
                    (grid.mesh.triangles()[t as usize], w.map(|x| x as f32))
                } else {
                    let c = grid.locator.locate(&grid.mesh, p);
                    ([c, c, c], [1.0, 0.0, 0.0])
                }
            })
            .collect();
        Self {
            width,
            height,
            stencils,
        }
    }
}

/// Labels drawn in the legend.
pub struct Legend<'a> {
    /// Title, left-aligned on the first line (e.g. field name and units).
    pub title: &'a str,
    /// Right-aligned on the first line (e.g. the date).
    pub subtitle: &'a str,
}

/// Format a legend number compactly.
pub fn format_value(v: f64) -> String {
    if v == 0.0 {
        "0".into()
    } else if v.abs() >= 1e5 || v.abs() < 1e-2 {
        format!("{v:.2e}")
    } else if v.abs() >= 100.0 {
        format!("{v:.0}")
    } else {
        format!("{v:.2}")
    }
}

/// Vectors to draw as arrows over the map.
pub struct Overlay<'a> {
    /// Where to draw arrows.
    pub layout: &'a ArrowLayout,
    /// One vector per cell of the overlay field's grid.
    pub vectors: &'a [DVec3],
    /// Speed that maps to a full-length arrow.
    pub max_speed: f64,
}

/// Render `values` (one per cell) as an indexed image with a legend.
pub fn render(
    raster: &CellRaster,
    values: &[f64],
    colormap: &Colormap,
    legend: &Legend<'_>,
    overlay: Option<&Overlay<'_>>,
) -> Vec<u8> {
    let (w, h) = (raster.width, raster.height);
    let mut pixels = vec![WHITE; w * (h + FOOTER)];
    pixels[..w * h]
        .par_iter_mut()
        .zip(&raster.stencils)
        .for_each(|(p, (cells, weights))| {
            let v: f64 = cells
                .iter()
                .zip(weights)
                .filter(|&(_, &w)| w != 0.0)
                .map(|(&c, &w)| f64::from(w) * values[c as usize])
                .sum();
            *p = colormap.index(v);
        });
    if let Some(o) = overlay {
        o.layout
            .draw(&mut pixels[..w * h], w, h, o.vectors, o.max_speed);
    }

    // Legend: title line, colour bar, min/max labels.
    let margin = 8;
    let line = 5 * TEXT_SCALE;
    draw_text(
        &mut pixels,
        w,
        margin,
        h + 4,
        legend.title,
        TEXT_SCALE,
        BLACK,
    );
    let sub_w = text_width(legend.subtitle, TEXT_SCALE);
    if sub_w + margin < w {
        draw_text(
            &mut pixels,
            w,
            w - margin - sub_w,
            h + 4,
            legend.subtitle,
            TEXT_SCALE,
            BLACK,
        );
    }
    let bar_top = h + 4 + line + 4;
    let bar_h = 8;
    let (min, max) = colormap.range();
    let (min_s, max_s) = (format_value(min), format_value(max));
    let label_w = text_width(&min_s, TEXT_SCALE).max(text_width(&max_s, TEXT_SCALE));
    let bar_x0 = margin + label_w + 6;
    let bar_x1 = w
        .saturating_sub(margin + label_w + 6)
        .max(bar_x0 + 1)
        .min(w);
    let bar_x0 = bar_x0.min(bar_x1);
    for x in bar_x0..bar_x1 {
        let index = colormap.bar_index((x - bar_x0) as f64 / (bar_x1 - bar_x0 - 1).max(1) as f64);
        for y in bar_top..bar_top + bar_h {
            pixels[y * w + x] = index;
        }
    }
    let label_y = bar_top + bar_h / 2 - line / 2;
    draw_text(&mut pixels, w, margin, label_y, &min_s, TEXT_SCALE, BLACK);
    draw_text(
        &mut pixels,
        w,
        bar_x1 + 6,
        label_y,
        &max_s,
        TEXT_SCALE,
        BLACK,
    );
    pixels
}

/// Height of a rendered image including the legend.
pub fn image_height(raster: &CellRaster) -> usize {
    raster.height + FOOTER
}

/// Write an indexed image as a palette PNG.
pub fn write_png(
    path: &Path,
    width: usize,
    height: usize,
    pixels: &[u8],
    palette: &[[u8; 3]; 256],
) -> std::io::Result<()> {
    let file = BufWriter::new(File::create(path)?);
    let mut encoder = png::Encoder::new(file, width as u32, height as u32);
    encoder.set_color(png::ColorType::Indexed);
    encoder.set_depth(png::BitDepth::Eight);
    encoder.set_palette(palette.concat());
    let mut writer = encoder.write_header().map_err(std::io::Error::other)?;
    writer
        .write_image_data(pixels)
        .map_err(std::io::Error::other)?;
    writer.finish().map_err(std::io::Error::other)
}

/// Write indexed frames as a looping GIF sharing one palette.
pub fn write_gif(
    path: &Path,
    width: usize,
    height: usize,
    frames: &[Vec<u8>],
    palette: &[[u8; 3]; 256],
    fps: u16,
) -> std::io::Result<()> {
    let file = BufWriter::new(File::create(path)?);
    let flat = palette.concat();
    let mut encoder = gif::Encoder::new(file, width as u16, height as u16, &flat)
        .map_err(std::io::Error::other)?;
    encoder
        .set_repeat(gif::Repeat::Infinite)
        .map_err(std::io::Error::other)?;
    let delay = (100 / fps.max(1)).max(2);
    for pixels in frames {
        let mut frame =
            gif::Frame::from_indexed_pixels(width as u16, height as u16, pixels.clone(), None);
        frame.delay = delay;
        encoder.write_frame(&frame).map_err(std::io::Error::other)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use fms_fields::ColorScale;

    #[test]
    fn narrow_images_do_not_panic() {
        let grid = Grid::new(4);
        let values: Vec<f64> = (0..grid.mesh.n_cells()).map(|i| i as f64).collect();
        let colormap = Colormap::new(ColorScale::Sequential, 0.0, 100.0);
        for width in [2, 16, 40, 64, 200] {
            for smooth in [false, true] {
                let raster = CellRaster::new(&grid, width, smooth);
                let legend = Legend {
                    title: "TEST.FIELD (M)",
                    subtitle: "YEAR 0 DAY 1",
                };
                let pixels = render(&raster, &values, &colormap, &legend, None);
                assert_eq!(pixels.len(), width * image_height(&raster));
            }
        }
    }

    #[test]
    fn smooth_and_flat_agree_on_constant_fields() {
        let grid = Grid::new(6);
        let values = vec![42.0; grid.mesh.n_cells()];
        let colormap = Colormap::new(ColorScale::Sequential, 0.0, 100.0);
        let legend = Legend {
            title: "",
            subtitle: "",
        };
        let flat = render(
            &CellRaster::new(&grid, 128, false),
            &values,
            &colormap,
            &legend,
            None,
        );
        let smooth = render(
            &CellRaster::new(&grid, 128, true),
            &values,
            &colormap,
            &legend,
            None,
        );
        assert_eq!(flat, smooth);
    }

    #[test]
    fn value_formatting() {
        assert_eq!(format_value(0.0), "0");
        assert_eq!(format_value(1372.4), "1372");
        assert_eq!(format_value(3.391), "3.39");
        assert_eq!(format_value(2.06e9), "2.06e9");
        assert_eq!(format_value(-0.001), "-1.00e-3");
    }
}
