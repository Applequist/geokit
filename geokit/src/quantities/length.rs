use crate::{
    math::fp::Float,
    units::{angle::RAD, length::LengthUnit},
};
use approx::AbsDiffEq;
use derive_more::derive::{Add, AddAssign, Display, Sub, SubAssign};
use std::ops::{Div, DivAssign, Mul, MulAssign};

use super::angle::Angle;

/// [Length] is a generic length value type.
/// The internal representation is a [Float] value in meters.
///
/// # Creation
///
/// There are 2 ways to create a [Length] value:
/// 1. By multiplying a [Float] quantity by a [LengthUnit] value.
///    The quantity expressed in the given unit is then converted to meters:
/// ```
/// use geokit::quantities::length::Length;
/// use geokit::units::length::US_FT;
/// let l: Length = 100. * US_FT;
/// ```
///
/// 2. Or by using [Length::new], passing a quantity and a [LengthUnit].
///    The quantity expressed in the given unit is then converted to meters:
/// ```
/// use geokit::quantities::length::Length;
/// use geokit::units::length::US_FT;
/// let l = Length::new(100., US_FT);
/// ```
///
/// The 2nd way can also be used in `const` context.
///
/// # Operations
///
/// [Length] supports the following operations:
/// - Addition, subtraction,
/// - negation,
/// - division
/// - multiplication by a scalar (left and right)
/// - division by a scalor (right)
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd, Add, AddAssign, Sub, SubAssign, Display)]
#[display("{} m", _0)]
pub struct Length(Float);

impl Length {
    pub const ZERO: Length = Length(0.0);

    /// 1e-4 m is considered *tiny*
    #[inline]
    pub const fn tiny() -> Length {
        Length(1e-4)
    }

    /// 1e-2 m is considered *small*
    #[inline]
    pub const fn small() -> Length {
        Length(1e-2)
    }

    pub const fn default_epsilon() -> Length {
        Self::tiny()
    }

    /// Create a length value whose `qty` is given in `unit`.
    #[inline]
    pub const fn new(qty: Float, unit: LengthUnit) -> Self {
        Self(qty * unit.m_per_unit())
    }

    /// Return this length value in the given unit.
    ///
    /// # Example
    /// ```
    /// use geokit::quantities::length::Length;
    /// use geokit::units::length::{M, US_FT};
    /// let l_m: Length = 10. * M;
    /// let l_ft = l_m.val(US_FT);
    /// assert_eq!(l_m, l_ft * US_FT);
    /// ```
    pub fn val(self, unit: LengthUnit) -> Float {
        self.0 / unit.m_per_unit()
    }
    /// Return this length value in meters.
    pub fn m(self) -> Float {
        self.0
    }

    pub fn abs(self) -> Length {
        Self(self.0.abs())
    }

    pub fn atan2(self, other: Self) -> Angle {
        Angle::new(self.m().atan2(other.m()), RAD)
    }

    pub fn hypot(self, other: Self) -> Length {
        Self(self.m().hypot(other.m()))
    }
}

impl Mul<LengthUnit> for Float {
    type Output = Length;

    /// Multiplying a [Float] value by a [LengthUnit] returns a [Length]
    fn mul(self, rhs: LengthUnit) -> Self::Output {
        Length::new(self, rhs)
    }
}

impl Mul<Float> for LengthUnit {
    type Output = Length;

    /// Multiplying a [Float] value by a [LengthUnit] returns a [Length]
    fn mul(self, rhs: Float) -> Self::Output {
        rhs * self
    }
}

impl Mul<Length> for Float {
    type Output = Length;

    /// Mutliplying [Length] by a scalar returns a [Length].
    fn mul(self, rhs: Length) -> Self::Output {
        Length(self * rhs.0)
    }
}

impl Mul<Float> for Length {
    type Output = Length;

    /// Mutliplying [Length] by a scalar returns a [Length].
    fn mul(self, rhs: Float) -> Self::Output {
        rhs * self
    }
}

impl MulAssign<Float> for Length {
    /// Mutliplying [Length] by a scalar returns a [Length].
    fn mul_assign(&mut self, rhs: Float) {
        self.0 *= rhs;
    }
}

impl Div<Float> for Length {
    type Output = Length;

    /// Dividing a [Length] by a [Float] scalar returns a [Length].
    fn div(self, rhs: Float) -> Self::Output {
        Length(self.0 / rhs)
    }
}

impl DivAssign<Float> for Length {
    /// Divide this [Length] by a [Float] scalar.
    fn div_assign(&mut self, rhs: Float) {
        self.0 /= rhs;
    }
}

impl Div for Length {
    type Output = Float;

    /// Compute the ratio of two [Length] as a [Float] scalar.
    fn div(self, rhs: Self) -> Self::Output {
        self.0 / rhs.0
    }
}

impl AbsDiffEq for Length {
    type Epsilon = Length;

    fn default_epsilon() -> Self::Epsilon {
        Length::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        (*self - *other).abs() <= epsilon
    }
}

/// [ArcLength] represents the length of an arc of circle.
/// It is convenient to represent the following type of operations
/// used in geometric computations:
/// [ArcLength] = [Length] * [Angle] or [Angle] * [Length]
/// [Angle] = [ArcLength] / [Length]
pub struct ArcLength(pub Length);

impl ArcLength {
    #[inline]
    pub fn length(self) -> Length {
        self.0
    }
}

impl Mul<Angle> for Length {
    type Output = ArcLength;

    #[inline]
    fn mul(self, rhs: Angle) -> Self::Output {
        ArcLength(self * rhs.rad())
    }
}

impl Mul<Length> for Angle {
    type Output = ArcLength;

    fn mul(self, rhs: Length) -> Self::Output {
        rhs * self
    }
}

impl Div<Length> for ArcLength {
    type Output = Angle;

    fn div(self, rhs: Length) -> Self::Output {
        (self.0 / rhs) * RAD
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
