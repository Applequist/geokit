use crate::units::length::LengthUnit;
use approx::AbsDiffEq;
use derive_more::derive::{Add, AddAssign, Display, Sub, SubAssign};
use std::ops::{Div, DivAssign, Mul, MulAssign};

/// [Length] is a generic length value type.
///
/// # Creation
///
/// There are 2 ways to create a [Length] value:
/// 1. By multiplying a [f64] quantity by a [LengthUnit] value:
/// ```
/// let l: Length = 100. * US_FT;
/// ```
///
/// 2. Or by using [Length::new], passing a quantity and a [LengthUnit]:
/// ```
/// let l = Length::new(100., US_FT);
/// ```
/// The 2nd way can also be used in `const` context.
///
/// # Operations
///
/// [Length] supports the following operations:
/// - Addition, subtraction,
/// - negation (yep [Length] can be negative:)),
/// - multiplication by a scalar (left and right)
/// - division by a scalor (right)
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd, Add, AddAssign, Sub, SubAssign, Display)]
#[display("{} m", _0)]
pub struct Length(f64);

impl Length {
    /// Create a length value whose `qty` is given in `unit`.
    #[inline]
    pub const fn new(qty: f64, unit: LengthUnit) -> Self {
        Self(qty * unit.m_per_unit())
    }

    /// Create a zero valued length.
    #[inline]
    pub fn zero() -> Self {
        Self(0.0)
    }

    /// Return this length value in the given unit.
    ///
    /// # Example
    /// ```
    /// let l: Length = 10. * M;
    /// let l_ft: f64 = l_m.val(US_FT);
    /// assert_eq!(l, l_ft * US_FT);
    /// ```
    pub fn val(self, unit: LengthUnit) -> f64 {
        self.0 / unit.m_per_unit()
    }

    /// Return this length value in meters.
    pub fn m(self) -> f64 {
        self.0
    }
}

impl Mul<Length> for f64 {
    type Output = Length;

    fn mul(self, rhs: Length) -> Self::Output {
        Length(self * rhs.0)
    }
}

impl Mul<f64> for Length {
    type Output = Length;

    fn mul(self, rhs: f64) -> Self::Output {
        rhs * self
    }
}

impl MulAssign<f64> for Length {
    fn mul_assign(&mut self, rhs: f64) {
        self.0 *= rhs;
    }
}

impl Div<f64> for Length {
    type Output = Length;

    fn div(self, rhs: f64) -> Self::Output {
        Length(self.0 / rhs)
    }
}

impl DivAssign<f64> for Length {
    fn div_assign(&mut self, rhs: f64) {
        self.0 /= rhs;
    }
}

impl AbsDiffEq for Length {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

impl Mul<LengthUnit> for f64 {
    type Output = Length;

    fn mul(self, rhs: LengthUnit) -> Self::Output {
        Length::new(self, rhs)
    }
}

impl Mul<f64> for LengthUnit {
    type Output = Length;

    fn mul(self, rhs: f64) -> Self::Output {
        rhs * self
    }
}

#[cfg(test)]
mod test {
    use crate::units::length::M;

    #[test]
    fn length_display() {
        assert_eq!(format!("{}", 1.2 * M), "1.2 m");
    }
}
