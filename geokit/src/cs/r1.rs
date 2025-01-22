//! The R1 *abstract* coordinate space is used to represents points on the oriented line.
//!
//! Each point on the line is represented by the distance of the point from an origin on the line.
//!
//! This module provides the following value types:
//! - [Length] to represent the coordinate of a single R1 point coordinate.
//!

use crate::{math::Float, units::length::LengthUnit};
use approx::AbsDiffEq;
use derive_more::derive::{Add, AddAssign, Display, Sub, SubAssign};
use std::ops::{Div, DivAssign, Mul, MulAssign};

/// [Length] is a generic length value type used to expressed coordinates
/// in cartesion CS.
///
/// # Creation
///
/// There are 2 ways to create a [Length] value:
/// 1. By multiplying a [Float] quantity by a [LengthUnit] value:
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
pub struct Length(Float);

impl Length {
    pub const ZERO: Length = Length(0.0);

    /// Create a length value whose `qty` is given in `unit`.
    #[inline]
    pub const fn new(qty: Float, unit: LengthUnit) -> Self {
        Self(qty * unit.m_per_unit())
    }

    /// Return this length value in the given unit.
    ///
    /// # Example
    /// ```
    /// let l: Length = 10. * M;
    /// let l_ft: Float = l_m.val(US_FT);
    /// assert_eq!(l, l_ft * US_FT);
    /// ```
    pub fn val(self, unit: LengthUnit) -> Float {
        self.0 / unit.m_per_unit()
    }

    /// Return this length value in meters.
    pub fn m(self) -> Float {
        self.0
    }
}

impl Mul<Length> for Float {
    type Output = Length;

    fn mul(self, rhs: Length) -> Self::Output {
        Length(self * rhs.0)
    }
}

impl Mul<Float> for Length {
    type Output = Length;

    fn mul(self, rhs: Float) -> Self::Output {
        rhs * self
    }
}

impl MulAssign<Float> for Length {
    fn mul_assign(&mut self, rhs: Float) {
        self.0 *= rhs;
    }
}

impl Div<Float> for Length {
    type Output = Length;

    fn div(self, rhs: Float) -> Self::Output {
        Length(self.0 / rhs)
    }
}

impl DivAssign<Float> for Length {
    fn div_assign(&mut self, rhs: Float) {
        self.0 /= rhs;
    }
}

impl Div for Length {
    type Output = Float;

    fn div(self, rhs: Self) -> Self::Output {
        self.0 / rhs.0
    }
}

impl AbsDiffEq for Length {
    type Epsilon = Float;

    fn default_epsilon() -> Self::Epsilon {
        Float::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

impl Mul<LengthUnit> for Float {
    type Output = Length;

    fn mul(self, rhs: LengthUnit) -> Self::Output {
        Length::new(self, rhs)
    }
}

impl Mul<Float> for LengthUnit {
    type Output = Length;

    fn mul(self, rhs: Float) -> Self::Output {
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
