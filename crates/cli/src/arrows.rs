//! Arrow overlays for vector fields on equirectangular maps.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use fms_grid::coords::from_lat_lon;
use fms_sim::Grid;
use glam::DVec3;

use crate::colormap::{BLACK, WHITE};

/// Where arrows go, and the local east/north directions there.
pub struct ArrowLayout {
    anchors: Vec<Anchor>,
    spacing: f64,
}

struct Anchor {
    x: f64,
    y: f64,
    cell: u32,
    east: DVec3,
    north: DVec3,
}

impl ArrowLayout {
    /// Anchors roughly `spacing` pixels apart on the map, thinned towards the
    /// poles so arrows stay evenly spaced on the sphere.
    pub fn new(grid: &Grid, width: usize, spacing: usize) -> Self {
        let height = width / 2;
        let spacing = spacing.max(4) as f64;
        let rows = (height as f64 / spacing).floor().max(1.0) as usize;
        let mut anchors = Vec::new();
        for r in 0..rows {
            let y = (r as f64 + 0.5) * height as f64 / rows as f64;
            let lat = FRAC_PI_2 - PI * y / height as f64;
            let cos_lat = fms_math::cos(lat);
            let count = ((width as f64 * cos_lat / spacing).round() as usize).max(1);
            for k in 0..count {
                let x = (k as f64 + 0.5) * width as f64 / count as f64;
                let lon = -PI + TAU * x / width as f64;
                let p = from_lat_lon(lat, lon);
                let (sin_lon, cos_lon) = fms_math::sin_cos(lon);
                let sin_lat = fms_math::sin(lat);
                anchors.push(Anchor {
                    x,
                    y,
                    cell: grid.locator.locate(&grid.mesh, p),
                    east: DVec3::new(-sin_lon, cos_lon, 0.0),
                    north: DVec3::new(-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat),
                });
            }
        }
        Self { anchors, spacing }
    }

    /// Draw one arrow per anchor. Lengths scale with speed, so the fastest
    /// `max_speed` gets an arrow ~0.9 × spacing long.
    pub fn draw(
        &self,
        pixels: &mut [u8],
        width: usize,
        height: usize,
        vectors: &[DVec3],
        max_speed: f64,
    ) {
        if max_speed <= 0.0 {
            return;
        }
        let scale = 0.9 * self.spacing / max_speed;
        let mut lines = Vec::new();
        for a in &self.anchors {
            let v = vectors[a.cell as usize];
            // Screen y points down, so north is −y.
            let (dx, dy) = (v.dot(a.east) * scale, -v.dot(a.north) * scale);
            let len = fms_math::hypot(dx, dy);
            if len < 1.5 {
                continue;
            }
            let (x0, y0) = (a.x - dx / 2.0, a.y - dy / 2.0);
            let (x1, y1) = (a.x + dx / 2.0, a.y + dy / 2.0);
            lines.push((x0, y0, x1, y1));
            // Arrowhead: two strokes 25° either side of the shaft.
            let head = (len * 0.35).clamp(2.0, 6.0);
            let (ux, uy) = (dx / len, dy / len);
            let (c, s) = (0.906, 0.423);
            for side in [1.0, -1.0] {
                let (hx, hy) = (ux * c - side * uy * s, side * ux * s + uy * c);
                lines.push((x1, y1, x1 - hx * head, y1 - hy * head));
            }
        }
        // White halo first, then black core, so arrows read on any colour.
        for &(x0, y0, x1, y1) in &lines {
            draw_line(pixels, width, height, x0, y0, x1, y1, 1, WHITE);
        }
        for &(x0, y0, x1, y1) in &lines {
            draw_line(pixels, width, height, x0, y0, x1, y1, 0, BLACK);
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn draw_line(
    pixels: &mut [u8],
    width: usize,
    height: usize,
    x0: f64,
    y0: f64,
    x1: f64,
    y1: f64,
    radius: i64,
    color: u8,
) {
    let steps = (x1 - x0).abs().max((y1 - y0).abs()).ceil().max(1.0) as usize;
    for i in 0..=steps {
        let t = i as f64 / steps as f64;
        let (x, y) = (
            (x0 + t * (x1 - x0)).round() as i64,
            (y0 + t * (y1 - y0)).round() as i64,
        );
        for dy in -radius..=radius {
            for dx in -radius..=radius {
                let (px, py) = (x + dx, y + dy);
                if px >= 0 && py >= 0 && (px as usize) < width && (py as usize) < height {
                    pixels[py as usize * width + px as usize] = color;
                }
            }
        }
    }
}
