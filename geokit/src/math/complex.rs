use crate::{math::fp::Float, quantities::angle::Angle, units::angle::RAD};
use num::Zero;
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

    pub fn conj(self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }

    /// Return the squared absolute value of this complex.
    pub fn abs_sq(self) -> Float {
        self.re * self.re + self.im * self.im
    }

    /// Return the module or absolute value of this complex.
    pub fn abs(self) -> Float {
        self.re.hypot(self.im)
    }

    /// Return the argument **in (-pi..pi] radians** of this complex.
    pub fn arg(self) -> Angle {
        self.im.atan2(self.re) * RAD
    }

    /// Return the multiplicative inverse of this complex.
    pub fn inv(self) -> Result<Self, &'static str> {
        if self.is_zero() {
            Err("Complex not invertible: |z| = 0.0")
        } else {
            Ok(self.conj() / self.abs_sq())
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
    type Output = Result<Self, &'static str>;

    fn div(self, rhs: Self) -> Self::Output {
        let inv = rhs.inv()?;
        Ok(self * inv)
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
    use crate::quantities::angle::Angle;

    #[test]
    fn test_eq() {
        assert_eq!(Complex::new(1., 0.), Complex::new(1., 0.));
        assert_ne!(Complex::new(1., 0.1), Complex::new(1., 0.));
    }

    #[test]
    fn test_arg() {
        assert_eq!(Complex::new(1., 1.).arg(), Angle::PI / 4.);
        assert_eq!(Complex::new(0.5, 3f64.sqrt() / 2.).arg(), Angle::PI / 3.);
        assert_eq!(Complex::new(-1., 0.).arg(), Angle::PI);
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
        let inv_z = z.inv().unwrap();
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
        let div = Complex::new(1., 1.) / Complex::new(1., -1.);
        assert!(div.is_ok());
        assert_eq!(div.unwrap(), Complex::new(0., 1.));

        assert_eq!(Complex::new(1., 1.) / 2., Complex::new(1. / 2., 1. / 2.));
    }
}
