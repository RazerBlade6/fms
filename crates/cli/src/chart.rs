//! Line charts of time series, one stacked panel per series sharing the time
//! axis. Drawn as indexed images with the bitmap font, like the maps.

use crate::font::{draw_text, text_width};
use crate::render::format_value;

const BACKGROUND: u8 = 0;
const TEXT: u8 = 1;
const GRID: u8 = 2;
const AXIS: u8 = 3;
const LINE: u8 = 4;

const SCALE: usize = 2;
const LINE_H: usize = 5 * SCALE;
const HEADER: usize = 12 + LINE_H;
const PANEL_TITLE: usize = 8 + LINE_H;
const PLOT_H: usize = 140;
const PANEL_GAP: usize = 14;
const FOOTER: usize = 14 + 2 * LINE_H;
const MARGIN: usize = 10;

/// The chart palette: background, text, grid, axes and the data line.
pub fn palette() -> [[u8; 3]; 256] {
    let mut p = [[255, 255, 255]; 256];
    p[TEXT as usize] = [0, 0, 0];
    p[GRID as usize] = [225, 225, 225];
    p[AXIS as usize] = [120, 120, 120];
    p[LINE as usize] = [31, 119, 180];
    p
}

/// One panel: a title and one value per time sample.
pub struct Series<'a> {
    /// Panel title (field, units, where it was read).
    pub title: String,
    /// Values, aligned with the chart's times. Non-finite values leave gaps.
    pub values: &'a [f64],
}

/// The time axis.
pub struct TimeAxis<'a> {
    /// Sample times, in the axis unit, increasing.
    pub times: &'a [f64],
    /// Axis caption, e.g. `HOURS` or `DAYS`.
    pub caption: &'a str,
    /// Whether the unit is hours, which gets ticks at 1, 3, 6, 12, 24 … h.
    pub hours: bool,
}

/// Render the chart; returns the pixels and the image height.
pub fn render_chart(
    width: usize,
    header: &str,
    axis: &TimeAxis<'_>,
    series: &[Series<'_>],
) -> (Vec<u8>, usize) {
    let height = HEADER + series.len() * (PANEL_TITLE + PLOT_H + PANEL_GAP) + FOOTER;
    let mut px = vec![BACKGROUND; width * height];
    draw_text(&mut px, width, MARGIN, 6, header, SCALE, TEXT);

    // Every panel shares the left margin, wide enough for the longest y label.
    let y_axes: Vec<(f64, f64, Vec<f64>)> = series.iter().map(|s| y_axis(s.values)).collect();
    let label_w = y_axes
        .iter()
        .flat_map(|(_, _, ticks)| tick_labels(ticks))
        .map(|label| text_width(&label, SCALE))
        .max()
        .unwrap_or(0);
    let x0 = MARGIN + label_w + 6;
    let x1 = width.saturating_sub(MARGIN + 8).max(x0 + 2);

    let (t_min, t_max) = match (axis.times.first(), axis.times.last()) {
        (Some(&a), Some(&b)) if b > a => (a, b),
        (Some(&a), _) => (a - 0.5, a + 0.5),
        _ => (0.0, 1.0),
    };
    let x_of = |t: f64| x0 as f64 + (t - t_min) / (t_max - t_min) * (x1 - 1 - x0) as f64;
    let x_ticks = ticks(t_min, t_max, 8, axis.hours);

    let mut top = HEADER;
    for (s, (y_min, y_max, y_ticks)) in series.iter().zip(&y_axes) {
        draw_text(&mut px, width, x0, top + 2, &s.title, SCALE, TEXT);
        let (p0, p1) = (top + PANEL_TITLE, top + PANEL_TITLE + PLOT_H);
        let y_of = |v: f64| (p1 - 1) as f64 - (v - y_min) / (y_max - y_min) * (p1 - 1 - p0) as f64;

        for &t in &x_ticks {
            let x = x_of(t).round() as usize;
            vline(&mut px, width, x, p0, p1, GRID);
        }
        for (&v, label) in y_ticks.iter().zip(tick_labels(y_ticks)) {
            let y = y_of(v).round() as usize;
            hline(&mut px, width, x0, x1, y, GRID);
            let lx = x0 - 6 - text_width(&label, SCALE);
            draw_text(
                &mut px,
                width,
                lx,
                y.saturating_sub(LINE_H / 2),
                &label,
                SCALE,
                TEXT,
            );
        }
        hline(&mut px, width, x0, x1, p0, AXIS);
        hline(&mut px, width, x0, x1, p1 - 1, AXIS);
        vline(&mut px, width, x0, p0, p1, AXIS);
        vline(&mut px, width, x1 - 1, p0, p1, AXIS);

        let points: Vec<Option<(f64, f64)>> = axis
            .times
            .iter()
            .zip(s.values)
            .map(|(&t, &v)| v.is_finite().then(|| (x_of(t), y_of(v))))
            .collect();
        for pair in points.windows(2) {
            if let [Some(a), Some(b)] = pair {
                segment(&mut px, width, *a, *b, LINE);
            }
        }
        // Points with no finite neighbour would otherwise be invisible.
        for (i, p) in points.iter().enumerate() {
            let lonely = i.checked_sub(1).is_none_or(|j| points[j].is_none())
                && points.get(i + 1).is_none_or(|q| q.is_none());
            if let (Some(p), true) = (p, lonely) {
                segment(&mut px, width, *p, *p, LINE);
            }
        }
        top = p1 + PANEL_GAP;
    }

    // Time labels under the last panel.
    let label_y = top - PANEL_GAP + 4;
    for (&t, label) in x_ticks.iter().zip(tick_labels(&x_ticks)) {
        let w = text_width(&label, SCALE);
        let x = (x_of(t).round() as usize)
            .saturating_sub(w / 2)
            .min(width.saturating_sub(w));
        draw_text(&mut px, width, x, label_y, &label, SCALE, TEXT);
    }
    let cw = text_width(axis.caption, SCALE);
    let cx = (x0 + x1).saturating_sub(cw) / 2;
    draw_text(
        &mut px,
        width,
        cx,
        label_y + LINE_H + 4,
        axis.caption,
        SCALE,
        TEXT,
    );
    (px, height)
}

/// Padded value range and tick positions for one series.
fn y_axis(values: &[f64]) -> (f64, f64, Vec<f64>) {
    let (lo, hi) = values
        .iter()
        .filter(|v| v.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
            (lo.min(v), hi.max(v))
        });
    let (lo, hi) = if !lo.is_finite() {
        (0.0, 1.0)
    } else if hi > lo {
        let pad = 0.05 * (hi - lo);
        (lo - pad, hi + pad)
    } else {
        let pad = if lo == 0.0 { 0.5 } else { 0.05 * lo.abs() };
        (lo - pad, hi + pad)
    };
    (lo, hi, ticks(lo, hi, 4, false))
}

/// Labels for evenly spaced ticks, with as many decimals as the spacing
/// needs (`0.2 0.4`, `6 12`), so neighbours never print the same.
fn tick_labels(ticks: &[f64]) -> Vec<String> {
    let step = ticks
        .windows(2)
        .map(|w| w[1] - w[0])
        .fold(f64::INFINITY, f64::min);
    let huge = ticks.iter().any(|t| t.abs() >= 1e5);
    if huge || step < 1e-3 {
        return ticks.iter().map(|&t| format_value(t)).collect();
    }
    let decimals = if step.is_finite() {
        (-fms_math::log10(step) - 1e-9).ceil().max(0.0) as usize
    } else {
        // One tick: show it as a whole number if it is one.
        usize::from(ticks.iter().any(|t| t.fract() != 0.0)) * 2
    };
    ticks.iter().map(|t| format!("{t:.decimals$}")).collect()
}

/// About `target` evenly spaced round values in `[lo, hi]`.
pub fn ticks(lo: f64, hi: f64, target: usize, hours: bool) -> Vec<f64> {
    let step = nice_step((hi - lo) / target.max(1) as f64, hours);
    let first = (lo / step).ceil() as i64;
    let last = (hi / step).floor() as i64;
    // `+ 0.0` turns -0 into 0 so it prints as "0".
    (first..=last).map(|k| k as f64 * step + 0.0).collect()
}

/// The "round" step nearest `raw` (by ratio): 1, 2 or 5 × 10ⁿ, or for hours
/// 1, 2, 3, 6 or 12 h, then whole days (24 h × a round number of days).
fn nice_step(raw: f64, hours: bool) -> f64 {
    if !raw.is_finite() || raw <= 0.0 {
        return 1.0;
    }
    let nearest = |candidates: &mut dyn Iterator<Item = f64>| {
        candidates
            .map(|c| (c, (c / raw).max(raw / c)))
            .fold(
                (1.0, f64::INFINITY),
                |best, x| if x.1 < best.1 { x } else { best },
            )
            .0
    };
    if hours && raw < 18.0 {
        return nearest(&mut [1.0, 2.0, 3.0, 6.0, 12.0, 24.0].into_iter());
    }
    if hours {
        return 24.0 * nice_step(raw / 24.0, false);
    }
    let magnitude = 10f64.powi(fms_math::log10(raw).floor() as i32);
    nearest(&mut [1.0, 2.0, 5.0, 10.0].into_iter().map(|m| m * magnitude))
}

fn put(px: &mut [u8], width: usize, x: usize, y: usize, c: u8) {
    if x < width && y < px.len() / width {
        px[y * width + x] = c;
    }
}

fn hline(px: &mut [u8], width: usize, x0: usize, x1: usize, y: usize, c: u8) {
    for x in x0..x1 {
        put(px, width, x, y, c);
    }
}

fn vline(px: &mut [u8], width: usize, x: usize, y0: usize, y1: usize, c: u8) {
    for y in y0..y1 {
        put(px, width, x, y, c);
    }
}

/// A 2-pixel-wide line segment (Bresenham on rounded endpoints).
fn segment(px: &mut [u8], width: usize, a: (f64, f64), b: (f64, f64), c: u8) {
    let (mut x, mut y) = (a.0.round() as i64, a.1.round() as i64);
    let (xb, yb) = (b.0.round() as i64, b.1.round() as i64);
    let (dx, dy) = ((xb - x).abs(), -(yb - y).abs());
    let (sx, sy) = (if x < xb { 1 } else { -1 }, if y < yb { 1 } else { -1 });
    let mut err = dx + dy;
    loop {
        for (ox, oy) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
            let (qx, qy) = (x + ox, y + oy);
            if qx >= 0 && qy >= 0 {
                put(px, width, qx as usize, qy as usize, c);
            }
        }
        if x == xb && y == yb {
            break;
        }
        let e2 = 2 * err;
        if e2 >= dy {
            err += dy;
            x += sx;
        }
        if e2 <= dx {
            err += dx;
            y += sy;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ticks_are_round_and_inside_the_range() {
        let tenths = ticks(0.0, 1.0, 4, false);
        assert_eq!(tenths.len(), 6);
        assert!(
            tenths
                .iter()
                .enumerate()
                .all(|(k, &t)| (t - 0.2 * k as f64).abs() < 1e-12)
        );
        assert_eq!(ticks(-3.0, 12.0, 4, false), [0.0, 5.0, 10.0]);
        assert_eq!(
            ticks(0.0, 24.0, 8, true),
            [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0]
        );
        assert_eq!(
            ticks(0.0, 72.0, 8, true),
            [0.0, 12.0, 24.0, 36.0, 48.0, 60.0, 72.0]
        );
        assert_eq!(
            ticks(0.0, 365.0, 8, false),
            [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0]
        );
        assert_eq!(ticks(-1.0, -0.0, 1, false), [-1.0, 0.0]);
        assert!(ticks(0.0, 0.0, 4, false).len() == 1);
    }

    #[test]
    fn tick_labels_use_the_spacing_precision() {
        assert_eq!(tick_labels(&[0.0, 6.0, 12.0]), ["0", "6", "12"]);
        assert_eq!(
            tick_labels(&[436.2, 436.4, 436.6]),
            ["436.2", "436.4", "436.6"]
        );
        assert_eq!(tick_labels(&[-0.05, 0.0, 0.05]), ["-0.05", "0.00", "0.05"]);
        assert_eq!(tick_labels(&[0.0, 2e5]), ["0", "2.00e5"]);
        assert_eq!(tick_labels(&[3.0]), ["3"]);
        assert_eq!(tick_labels(&[0.5]), ["0.50"]);
        assert!(tick_labels(&[]).is_empty());
    }

    #[test]
    fn charts_have_the_expected_size_and_draw_data() {
        let times: Vec<f64> = (0..200).map(f64::from).collect();
        let a: Vec<f64> = times.iter().map(|t| fms_math::sin(t / 10.0)).collect();
        let b = vec![f64::NAN; 200];
        let c = vec![5.0; 200];
        let axis = TimeAxis {
            times: &times,
            caption: "HOURS",
            hours: true,
        };
        for width in [40, 300, 1024] {
            let series = [
                Series {
                    title: "A".into(),
                    values: &a,
                },
                Series {
                    title: "B".into(),
                    values: &b,
                },
                Series {
                    title: "C".into(),
                    values: &c,
                },
            ];
            let (px, h) = render_chart(width, "TEST", &axis, &series);
            assert_eq!(px.len(), width * h);
            if width >= 300 {
                assert!(px.contains(&LINE));
            }
        }
    }

    #[test]
    fn a_single_sample_still_shows() {
        let axis = TimeAxis {
            times: &[0.0],
            caption: "HOURS",
            hours: true,
        };
        let series = [Series {
            title: "ONE".into(),
            values: &[1.0],
        }];
        let (px, _) = render_chart(400, "", &axis, &series);
        assert!(px.contains(&LINE));
    }
}
