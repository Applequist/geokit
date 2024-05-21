/// A polynomial of a single variable of degree less than `D`.
/// All coefficients are stored on the stack.
pub struct Polynomial<const D: usize>([f64; D]);

impl<const D: usize> Polynomial<D> {
    /// Create a new polynomial of degree less than `D`.
    /// The n-th coefficient in the array is the coefficient of degree n
    /// of the polynomial.
    pub fn new(coef: [f64; D]) -> Self {
        Self(coef)
    }

    /// Evaluate this polynomial at the given real value `x` using Horner's method.
    pub fn eval_at(&self, x: f64) -> f64 {
        self.0
            .into_iter()
            .rev()
            .fold(0.0, |acc, coef_i| x * acc + coef_i)
    }
}

#[cfg(test)]
mod tests {
    use crate::math::polynomial::Polynomial;

    #[test]
    fn eval_at() {
        assert_eq!(Polynomial::new([0.0; 3]).eval_at(1.0), 0.0);
        assert_eq!(Polynomial::new([1., -1.0, 1.0]).eval_at(1.0), 1.0);
        assert_eq!(Polynomial::new([1.0, 0.0, 2.0]).eval_at(1.0), 3.0);
        assert_eq!(Polynomial::new([1.0, 0.0, 2.0]).eval_at(2.0), 9.0);
    }
}
