//! Provide systems of axes for various cartesian CS:
//! - [GeocentricAxes] for geocentric cartesian CS.
//! - [ProjectedAxes] for projected cartesian CS.
//!
//! Coordinates in cartesian CS are defined using distance from planes along axes.
//  For that this module defines a [Length] type.

use crate::units::length::LengthUnit;

/// The [GeocentricAxes] enum defines the *coordinates system* part of **geocentric** CS.
/// Geocentric CS uses [meter][crate::units::length::M] unit by default for all axes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeocentricAxes {
    /// Coordinates are x, y, and z in meters.
    XYZ,
}

/// A [ProjectedAxes] value defines the **coordinates system** part of a [ProjectedCrs].
/// That is:
/// - the ordering and direction of the axes,
/// - the horizontal length unit used for easting and northing,
/// - the height length unit used for the ellipsoidal height.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProjectedAxes {
    /// The axes for a 3D projected Crs. Coordinates are, in order:
    /// - easting positive eastward,
    /// - northing positive northward,
    /// - ellipsoidal height positive up.
    EastNorthUp {
        /// the length unit for easting and northing, eg
        /// `M` or `US_FT`
        horiz_unit: LengthUnit,
        /// the length unit for easting and northing, eg
        /// `M` or `US_FT`
        height_unit: LengthUnit,
    },
    /// The axes for a 2D projected Crs. Coordinates are, in order:
    /// - easting positive eastward,
    /// - northing positive northward,
    EastNorth {
        /// the length unit for easting and northing, eg
        /// `M` or `US_FT`
        horiz_unit: LengthUnit,
    },
}

impl ProjectedAxes {
    pub fn dim(&self) -> usize {
        match self {
            ProjectedAxes::EastNorthUp {
                horiz_unit: _,
                height_unit: _,
            } => 3,
            ProjectedAxes::EastNorth { horiz_unit: _ } => 2,
        }
    }
}

use approx::AbsDiffEq;
use derive_more::derive::{Add, AddAssign, Display, Sub, SubAssign};
use std::ops::{Div, DivAssign, Mul, MulAssign};

use super::r1::Real;

/// [Length] is a generic length value type used to expressed coordinates
/// in cartesion CS.
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
pub struct Length(Real);

impl Length {
    pub const ZERO: Length = Length(0.0);

    /// Create a length value whose `qty` is given in `unit`.
    #[inline]
    pub const fn new(qty: f64, unit: LengthUnit) -> Self {
        Self(qty * unit.m_per_unit())
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
