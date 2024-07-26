use std::ops::{Add, Mul};

use num::{Num, Zero};

/// A polynomial of a single variable of degree less than `D`.
/// All coefficients are stored on the stack.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Polynomial<const D: usize, T = f64>([T; D]);

impl<const D: usize, T> Polynomial<D, T>
where
    T: Num + Copy,
{
    /// Create a new polynomial of degree less than `D`.
    /// The n-th coefficient in the array is the coefficient of degree n
    /// of the polynomial.
    pub const fn new(coef: [T; D]) -> Self {
        Self(coef)
    }

    /// Evaluate this polynomial at the given value `x` using Horner's method.
    pub fn eval_at<V>(&self, x: V) -> V
    where
        V: Copy + Zero + Mul<Output = V> + Add<T, Output = V>,
    {
        self.0
            .into_iter()
            .rev()
            .fold(V::zero(), |acc, coef_i| x * acc + coef_i)
    }

    pub fn solve<V>(&self) -> Vec<V> {
        // TODO
        vec![]
    }
}

impl<const D: usize> Polynomial<D, f64> {
    pub fn fast_eval_at(&self, x: f64) -> f64 {
        self.0
            .into_iter()
            .rev()
            .fold(0.0f64, |acc, coef_i| acc.mul_add(x, coef_i))
    }
}

#[cfg(test)]
mod tests {
    use num::Zero;

    use crate::math::complex::Complex;
    use crate::math::polynomial::Polynomial;

    #[test]
    fn eval_at() {
        // eval_at(T)
        assert_eq!(Polynomial::new([0.0; 3]).eval_at(1.0), 0.0);
        assert_eq!(
            Polynomial::<3, f64>::new([0.0; 3]).fast_eval_at(1.0f64),
            0.0
        );
        assert_eq!(Polynomial::new([1., -1.0, 1.0]).eval_at(1.0), 1.0);
        assert_eq!(
            Polynomial::<3, f64>::new([1., -1.0, 1.0]).fast_eval_at(1.0f64),
            1.0
        );
        assert_eq!(Polynomial::new([1.0, 0.0, 2.0]).eval_at(1.0), 3.0);
        assert_eq!(
            Polynomial::<3, f64>::new([1.0, 0.0, 2.0]).fast_eval_at(1.0f64),
            3.0
        );
        assert_eq!(Polynomial::new([1.0, 0.0, 2.0]).eval_at(2.0), 9.0);
        assert_eq!(
            Polynomial::<3, f64>::new([1.0, 0.0, 2.0]).fast_eval_at(2.0f64),
            9.0
        );

        // eval_at(V)
        assert_eq!(
            Polynomial::new([1., 0., 1.]).eval_at(Complex::new(0., 1.)),
            Complex::zero()
        );
    }

    #[test]
    fn eq() {
        assert_eq!(
            Polynomial::new([1., -2., 1.]),
            Polynomial::new([1., -2., 1.])
        );
    }
}
