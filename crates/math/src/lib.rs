//! Portable, bit-reproducible math functions.
//!
//! `f64::sin`, `cos`, `exp` and friends call the platform's C math library,
//! which differs between macOS, Windows and Linux in the last bit or so. Over
//! millions of time steps those differences grow, and the same save would
//! diverge on different machines. These wrappers use the pure-Rust
//! [`libm`](https://crates.io/crates/libm) port of musl's math library instead,
//! which gives identical results everywhere.
//!
//! Basic arithmetic, `sqrt`, `abs`, `floor`, `round`, `min`/`max` and `powi`
//! are already exactly specified (IEEE 754 or compiler-builtins), so they stay
//! as std methods. Everything else must come from here: the workspace
//! `clippy.toml` forbids the std versions in simulation code.

#![allow(clippy::disallowed_methods)]

/// Sine (radians).
#[inline]
pub fn sin(x: f64) -> f64 {
    libm::sin(x)
}

/// Cosine (radians).
#[inline]
pub fn cos(x: f64) -> f64 {
    libm::cos(x)
}

/// `(sin x, cos x)`.
#[inline]
pub fn sin_cos(x: f64) -> (f64, f64) {
    libm::sincos(x)
}

/// Tangent (radians).
#[inline]
pub fn tan(x: f64) -> f64 {
    libm::tan(x)
}

/// Arcsine, in `[-π/2, π/2]`.
#[inline]
pub fn asin(x: f64) -> f64 {
    libm::asin(x)
}

/// Arccosine, in `[0, π]`.
#[inline]
pub fn acos(x: f64) -> f64 {
    libm::acos(x)
}

/// Arctangent, in `(-π/2, π/2)`.
#[inline]
pub fn atan(x: f64) -> f64 {
    libm::atan(x)
}

/// Four-quadrant arctangent of `y / x`, in `(-π, π]`.
#[inline]
pub fn atan2(y: f64, x: f64) -> f64 {
    libm::atan2(y, x)
}

/// `eˣ`.
#[inline]
pub fn exp(x: f64) -> f64 {
    libm::exp(x)
}

/// `eˣ − 1`, accurate for small `x`.
#[inline]
pub fn exp_m1(x: f64) -> f64 {
    libm::expm1(x)
}

/// Natural logarithm.
#[inline]
pub fn ln(x: f64) -> f64 {
    libm::log(x)
}

/// `ln(1 + x)`, accurate for small `x`.
#[inline]
pub fn ln_1p(x: f64) -> f64 {
    libm::log1p(x)
}

/// Base-10 logarithm.
#[inline]
pub fn log10(x: f64) -> f64 {
    libm::log10(x)
}

/// `xʸ` for real `y` (use `f64::powi` for integer exponents).
#[inline]
pub fn pow(x: f64, y: f64) -> f64 {
    libm::pow(x, y)
}

/// `√(x² + y²)` without overflow.
#[inline]
pub fn hypot(x: f64, y: f64) -> f64 {
    libm::hypot(x, y)
}

/// Cube root.
#[inline]
pub fn cbrt(x: f64) -> f64 {
    libm::cbrt(x)
}

/// Hyperbolic tangent.
#[inline]
pub fn tanh(x: f64) -> f64 {
    libm::tanh(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_6, PI};

    #[test]
    fn values() {
        assert!((sin(FRAC_PI_6) - 0.5).abs() < 1e-15);
        assert!((cos(PI) + 1.0).abs() < 1e-15);
        assert_eq!(sin_cos(0.3), (sin(0.3), cos(0.3)));
        assert!((atan2(1.0, 0.0) - FRAC_PI_2).abs() < 1e-15);
        assert!((asin(1.0) - FRAC_PI_2).abs() < 1e-15);
        assert!((exp(ln(7.5)) - 7.5).abs() < 1e-14);
        assert_eq!(hypot(3.0, 4.0), 5.0);
        assert!((pow(2.0, 0.5) - 2f64.sqrt()).abs() < 1e-15);
    }

    /// Golden bit patterns: these must be identical on every platform. CI runs
    /// this on macOS, Windows and Linux.
    #[test]
    fn bit_exact_golden_values() {
        let checksum = (1..=1000)
            .flat_map(|i| {
                let x = f64::from(i) * 0.012_345_678_9;
                [
                    sin(x),
                    cos(x),
                    atan2(x, 1.0 - x),
                    exp(-x),
                    ln(x),
                    asin(x.fract()),
                    tan(x),
                ]
            })
            .fold(0u64, |h, v| h.rotate_left(5) ^ v.to_bits());
        assert_eq!(checksum, GOLDEN, "libm results changed: {checksum:#018x}");
    }

    const GOLDEN: u64 = 0x85f8_fe99_99b3_8789;
}
