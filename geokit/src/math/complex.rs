use std::ops::{Add, Div, Mul, Sub};

use num::{Float, One, Zero};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Complex<T = f64> {
    pub re: T,
    pub im: T,
}

impl<T: Float> Complex<T> {
    pub fn new(re: T, im: T) -> Self {
        Self { re, im }
    }

    /// Return the squared absolute value of this complex.
    pub fn abs_sq(&self) -> T {
        self.re.clone() * self.re.clone() + self.im.clone() * self.im.clone()
    }

    /// Return the absolute value of this complex.
    pub fn abs(&self) -> T {
        self.re.hypot(self.im)
    }

    /// Return the argument **in radians** of this complex.
    pub fn arg(&self) -> T {
        self.im.atan2(self.re)
    }

    /// Return the multiplicative inverse of this complex.
    pub fn inv(&self) -> Self {
        let m = self.abs_sq();
        if m == T::zero() {
            panic!("Zero is an invalid denominator!");
        }
        Self {
            re: self.re / m,
            im: -self.im / m,
        }
    }
}

impl<T: Float> Zero for Complex<T> {
    fn zero() -> Self {
        Self::new(T::zero(), T::zero())
    }

    fn is_zero(&self) -> bool {
        self.re.is_zero() && self.im.is_zero()
    }
}

impl<T: Float> Add for Complex<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl<T: Float> Add<T> for Complex<T> {
    type Output = Self;

    fn add(self, rhs: T) -> Self::Output {
        Self {
            re: self.re + rhs,
            im: self.im,
        }
    }
}

impl<T: Float> Sub for Complex<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
}

impl<T: Float> Sub<T> for Complex<T> {
    type Output = Complex<T>;

    fn sub(self, rhs: T) -> Self::Output {
        Self {
            re: self.re - rhs,
            im: self.im,
        }
    }
}

impl<T: Float> One for Complex<T> {
    fn one() -> Self {
        Self::new(T::one(), T::zero())
    }
}

impl<T: Float> Mul for Complex<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl<T: Float> Mul<T> for Complex<T> {
    type Output = Complex<T>;

    fn mul(self, rhs: T) -> Self::Output {
        Self {
            re: self.re * rhs,
            im: self.im * rhs,
        }
    }
}

impl<T: Float> Div for Complex<T> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inv()
    }
}

impl<T: Float> Div<T> for Complex<T> {
    type Output = Complex<T>;

    fn div(self, rhs: T) -> Self::Output {
        Self {
            re: self.re / rhs,
            im: self.im / rhs,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_PI_4, PI};

    use crate::math::complex::Complex;

    #[test]
    fn test_eq() {
        assert_eq!(Complex::new(1., 0.), Complex::new(1., 0.));
        assert_ne!(Complex::new(1., 0.1), Complex::new(1., 0.));
    }

    #[test]
    fn test_arg() {
        assert_eq!(Complex::new(1., 1.).arg(), FRAC_PI_4);
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
            Complex::new(0., 2. / 2f64.sqrt())
        );
        assert_eq!(Complex::new(1., 1.) / 2., Complex::new(1. / 2., 1. / 2.));
    }
}
