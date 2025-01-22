use num::Zero;

use crate::math::Float;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Complex {
    pub re: Float,
    pub im: Float,
}

impl Complex {
    pub const fn zero() -> Self {
        Self { re: 0.0, im: 0.0 }
    }

    pub fn new(re: Float, im: Float) -> Self {
        Self { re, im }
    }

    /// Return the squared absolute value of this complex.
    pub fn abs_sq(&self) -> Float {
        self.re * self.re + self.im * self.im
    }

    /// Return the absolute value of this complex.
    pub fn abs(&self) -> Float {
        self.re.hypot(self.im)
    }

    /// Return the argument **in radians** of this complex.
    pub fn arg(&self) -> Float {
        self.im.atan2(self.re)
    }

    /// Return the multiplicative inverse of this complex.
    pub fn inv(&self) -> Self {
        let m = self.abs_sq();
        if m == 0.0 {
            panic!("Zero is an invalid denominator!");
        }
        Self {
            re: self.re / m,
            im: -self.im / m,
        }
    }
}

impl Zero for Complex {
    fn zero() -> Self {
        Self::new(0.0, 0.0)
    }

    fn is_zero(&self) -> bool {
        self.re == 0.0 && self.im == 0.0
    }
}

impl Add for Complex {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl Add<Float> for Complex {
    type Output = Self;

    fn add(self, rhs: Float) -> Self::Output {
        Self {
            re: self.re + rhs,
            im: self.im,
        }
    }
}

impl Sub for Complex {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
}

impl Sub<Float> for Complex {
    type Output = Complex;

    fn sub(self, rhs: Float) -> Self::Output {
        Self {
            re: self.re - rhs,
            im: self.im,
        }
    }
}

impl Mul for Complex {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl Mul<Float> for Complex {
    type Output = Complex;

    fn mul(self, rhs: Float) -> Self::Output {
        Self {
            re: self.re * rhs,
            im: self.im * rhs,
        }
    }
}

impl Div for Complex {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inv()
    }
}

impl Div<Float> for Complex {
    type Output = Complex;

    fn div(self, rhs: Float) -> Self::Output {
        Self {
            re: self.re / rhs,
            im: self.im / rhs,
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::math::complex::Complex;
    use crate::math::PI;
    use crate::math::PI_4;

    #[test]
    fn test_eq() {
        assert_eq!(Complex::new(1., 0.), Complex::new(1., 0.));
        assert_ne!(Complex::new(1., 0.1), Complex::new(1., 0.));
    }

    #[test]
    fn test_arg() {
        assert_eq!(Complex::new(1., 1.).arg(), PI_4);
        assert_eq!(Complex::new(0.5, 3f64.sqrt() / 2.).arg(), PI / 3.);
        assert_eq!(Complex::new(-1., 0.).arg(), PI);
    }

    #[test]
    fn test_abs() {
        assert_eq!(Complex::new(1., 1.).abs(), 2f64.sqrt());
        assert_eq!(Complex::new(0.5, 3f64.sqrt() / 2.).abs(), 1.);
        assert_eq!(Complex::new(-1., 0.).abs(), 1.);
    }

    #[test]
    fn test_inv() {
        let z = Complex::new(1., 2.);
        let inv_z = z.inv();
        assert_eq!(z * inv_z, Complex::new(1., 0.));
    }

    #[test]
    fn test_ops() {
        // Addition
        assert_eq!(
            Complex::new(1., 0.) + Complex::new(0., 1.),
            Complex::new(1., 1.)
        );

        // Subtraction
        assert_eq!(
            Complex::new(1., 0.) - Complex::new(0., 1.),
            Complex::new(1., -1.)
        );

        // Multiplication
        assert_eq!(
            Complex::new(1., 1.) * Complex::new(1., -1.),
            Complex::new(2., 0.)
        );
        assert_eq!(Complex::new(1., 1.) * 2., Complex::new(2., 2.));

        // Division
        assert_eq!(
            Complex::new(1., 1.) / Complex::new(1., -1.),
            Complex::new(0., 1.)
        );
        assert_eq!(Complex::new(1., 1.) / 2., Complex::new(1. / 2., 1. / 2.));
    }
}
