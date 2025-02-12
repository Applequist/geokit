use super::Float;

/// Error-free transformation of the sum of 2 floating point numbers.
/// The return pair contains the floating-point sum of `a` and `b` and
/// the floating-point error
pub fn sum(a: Float, b: Float) -> (Float, Float) {
    let x = a + b;
    let z = x - a;
    let y = a - (x - z) + (b - z);
    (x, y)
}
