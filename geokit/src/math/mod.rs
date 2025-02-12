/// TODO: Make config dependent
pub type Float = f64;
pub const ZERO: Float = 0.0;
pub const PI_2: Float = std::f64::consts::FRAC_PI_2;
pub const PI_4: Float = std::f64::consts::FRAC_PI_4;
pub const PI: Float = std::f64::consts::PI;
pub const TAU: Float = std::f64::consts::TAU;

/// Error-free transformation of the sum of 2 floating point numbers.
/// The return pair contains the floating-point sum of `a` and `b` and
/// the floating-point error
pub fn sum(a: Float, b: Float) -> (Float, Float) {
    let x = a + b;
    let z = x - a;
    let y = a - (x - z) + (b - z);
    (x, y)
}

pub mod arithmetic;
pub mod complex;
pub mod polynomial;

pub mod utils;
