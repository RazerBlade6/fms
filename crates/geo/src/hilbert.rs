//! 3D Hilbert curve keys, used to give mesh elements a cache-friendly order.

use glam::DVec3;

/// Bits per axis; 3 × 21 = 63 bits fit in a `u64`.
const BITS: u32 = 21;

/// Hilbert-curve key of a point in the cube `[-1, 1]³` (e.g. a unit vector).
///
/// Points that are close in space get nearby keys, so sorting by key gives an
/// order with good memory locality.
pub fn key(p: DVec3) -> u64 {
    let max = f64::from((1u32 << BITS) - 1);
    let quantize = |x: f64| (((x.clamp(-1.0, 1.0) + 1.0) * 0.5 * max).round()) as u32;
    index([quantize(p.x), quantize(p.y), quantize(p.z)])
}

/// Hilbert index of integer coordinates with `BITS` bits each, using Skilling's
/// transpose algorithm ("Programming the Hilbert curve", AIP Conf. Proc. 707, 2004).
fn index(mut x: [u32; 3]) -> u64 {
    let m = 1u32 << (BITS - 1);
    // Inverse undo of the excess work.
    let mut q = m;
    while q > 1 {
        let p = q - 1;
        for i in 0..3 {
            if x[i] & q != 0 {
                x[0] ^= p;
            } else {
                let t = (x[0] ^ x[i]) & p;
                x[0] ^= t;
                x[i] ^= t;
            }
        }
        q >>= 1;
    }
    // Gray encode.
    for i in 1..3 {
        x[i] ^= x[i - 1];
    }
    let mut t = 0;
    let mut q = m;
    while q > 1 {
        if x[2] & q != 0 {
            t ^= q - 1;
        }
        q >>= 1;
    }
    for v in &mut x {
        *v ^= t;
    }
    // Interleave the transposed bits into a single index, most significant first.
    let mut h = 0u64;
    for b in (0..BITS).rev() {
        for v in x {
            h = (h << 1) | u64::from((v >> b) & 1);
        }
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;

    /// On a 2-bit-per-axis grid (embedded in the full-precision grid), consecutive
    /// Hilbert indices must be unit steps apart: the defining property of the curve.
    #[test]
    fn consecutive_keys_are_adjacent_cells() {
        let shift = BITS - 2;
        let mut cells: Vec<(u64, [u32; 3])> = Vec::new();
        for x in 0..4u32 {
            for y in 0..4u32 {
                for z in 0..4u32 {
                    let c = [x, y, z];
                    cells.push((index(c.map(|v| v << shift)), c));
                }
            }
        }
        cells.sort_unstable();
        for w in cells.windows(2) {
            let d: u32 = (0..3).map(|i| w[0].1[i].abs_diff(w[1].1[i])).sum();
            assert_eq!(d, 1, "{:?} -> {:?}", w[0].1, w[1].1);
        }
        let mut keys: Vec<u64> = cells.iter().map(|c| c.0).collect();
        keys.dedup();
        assert_eq!(keys.len(), 64);
    }
}
