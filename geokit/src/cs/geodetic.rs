//! Geodetic coordinates spaces are used to represent points on or near the surface of an ellipsoid
//! of revolution.
//!
//! Points in those spaces are represented using:
//! - longitude,
//! - latitude,
//! - and in the 3D case, ellipsoidal height.
//!
//! A special set of axes is said to be **normalized**.
//!
//! This module a [system of axes][GeodeticAxes] for geodetic CS and
//! a set of value types to represent **normalized** geodetic coordinates.

use crate::math::Float;
use crate::units::angle::{AngleUnit, RAD};
use crate::units::length::{LengthUnit, M};
use approx::AbsDiffEq;
use derive_more::derive::{Display, Neg};
use std::ops::{Add, AddAssign, Sub, SubAssign};

use super::s1::{Angle, Interval};

/// A [GeodeticAxes] defines the possible set of axes used in **geodetic** CS.
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum GeodeticAxes {
    /// Coordinates are given in the following order:
    /// - longitude positive east of prime meridian, using `angle_unit` [AngleUnit]
    /// - latitude positive north of equatorial plane, using `angle_unit` [AngleUnit],
    /// - ellipsoidal height positive upward, using `height_unit` [LengthUnit]
    EastNorthUp {
        /// The angle unit used for longitude and latitude, eg
        /// `DEG` or `GRAD`
        angle_unit: AngleUnit,
        /// The length unit used for ellipsoidal height, eg
        /// 'M' or 'US_FOOT'
        height_unit: LengthUnit,
    },
    /// Coordinates are given in the following order:
    /// - latitude positive north of equatorial plane, using `angle_unit` [AngleUnit],
    /// - longitude positive east of prime meridian, using `angle_unit` [AngleUnit]
    /// - ellipsoidal height positive upward, using `height_unit` [LengthUnit]
    NorthEastUp {
        /// The angle unit used for longitude and latitude, eg
        /// `DEG` or `GRAD`
        angle_unit: AngleUnit,
        /// The length unit used for ellipsoidal height, eg
        /// 'M' or 'US_FOOT'
        height_unit: LengthUnit,
    },
    /// Coordinates are, in the given order:
    /// - longitude positive east of prime meridian, using `angle_unit` [AngleUnit]
    /// - latitude positive north of equatorial plane, using `angle_unit` [AngleUnit],
    EastNorth {
        /// The angle unit used for longitude and latitude, eg
        /// `DEG` or `GRAD`
        angle_unit: AngleUnit,
    },
    /// Coordinates are, in the given order:
    /// - latitude positive north of equatorial plane, using `angle_unit` [AngleUnit],
    /// - longitude positive east of prime meridian, using `angle_unit` [AngleUnit]
    NorthEast {
        /// The angle unit used for longitude and latitude, eg
        /// `DEG` or `GRAD`
        angle_unit: AngleUnit,
    },
    /// Coordinates are, in the given order:
    /// - latitude positive north of equatorial plane, using `angle_unit` [AngleUnit],
    /// - longitude positive west of prime meridian, using `angle_unit` [AngleUnit]
    NorthWest {
        /// The angle unit used for longitude and latitude, eg
        /// `DEG` or `GRAD`
        angle_unit: AngleUnit,
    },
}

impl GeodeticAxes {
    /// Return the dimension (2D or 3D) of the coordinate system.
    pub fn dim(&self) -> usize {
        match self {
            GeodeticAxes::EastNorthUp {
                angle_unit: _,
                height_unit: _,
            }
            | GeodeticAxes::NorthEastUp {
                angle_unit: _,
                height_unit: _,
            } => 3,
            GeodeticAxes::EastNorth { angle_unit: _ }
            | GeodeticAxes::NorthEast { angle_unit: _ }
            | GeodeticAxes::NorthWest { angle_unit: _ } => 2,
        }
    }
}

impl Default for GeodeticAxes {
    /// Return the [`GeodeticAxes`] used in **normalized geodetic coordinates**.
    fn default() -> Self {
        Self::EastNorthUp {
            angle_unit: RAD,
            height_unit: M,
        }
    }
}

/// A longitude coordinate in [-pi..pi] radians.
/// You can add, subtract an [Angle] from [Lon],
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd, Neg, Display)]
#[display("{}", self.0.to_dms())]
pub struct Lon(Angle);

impl Lon {
    pub const MIN: Lon = Lon(Angle::M_PI);
    pub const ZERO: Lon = Lon(Angle::ZERO);
    pub const MAX: Lon = Lon(Angle::PI);

    /// Create a new longitude value from a given angle.
    /// The angle is wrapped into [-pi..pi].
    pub fn new(val: Angle) -> Self {
        Self(val.wrapped())
    }

    /// Create a new longitude value from an [Angle] ***ALREADY*** in [-pi..pi]
    pub const fn const_new(val: Angle) -> Self {
        Self(val)
    }

    /// Create a new longitue value from a *dms* angle value.
    /// The angle value is converted to radians and wrapped into [-pi, pi].
    ///
    /// # Parameters
    ///
    /// - `d`: the number of degrees. Determine the sign of the returned angle.
    /// - `m`: the number of minutes. Must be >= 0.
    /// - `s`: the number of seconds with fractional part. Must be >= 0.
    ///
    /// # Examples
    ///
    /// ```
    /// let a: Lon = Lon::dms(2., 20., 14.02500);
    /// assert!(Lon::dms(-12., 45., 59.1234).rad() < 0.0);
    /// ```
    pub fn dms(d: Float, m: Float, s: Float) -> Self {
        Self::new(Angle::dms(d, m, s))
    }

    /// Normalize the longitude into (-pi..pi] such that any point on a parallel
    /// has a unique *normalized* longitude.
    pub fn normalize(self) -> Self {
        if self <= Self::MIN {
            Self::MAX
        } else {
            self
        }
    }

    /// Return the angle of the longitude.
    #[inline]
    pub fn angle(self) -> Angle {
        self.0
    }

    /// Return the longitude as a raw angle value **in radians**.
    #[inline]
    pub fn rad(self) -> Float {
        self.0.rad()
    }
}

impl Add<Angle> for Lon {
    type Output = Self;

    fn add(self, rhs: Angle) -> Self::Output {
        Self::new(self.0 + rhs)
    }
}

impl Add<Lon> for Angle {
    type Output = Lon;

    fn add(self, rhs: Lon) -> Self::Output {
        rhs + self
    }
}

impl AddAssign<Angle> for Lon {
    fn add_assign(&mut self, rhs: Angle) {
        self.0 = (self.0 + rhs).normalized();
    }
}

impl Sub<Angle> for Lon {
    type Output = Self;

    fn sub(self, rhs: Angle) -> Self::Output {
        // Watch out for infinite recursion with self - rhs
        Self::new(self.0 - rhs)
    }
}

impl SubAssign<Angle> for Lon {
    fn sub_assign(&mut self, rhs: Angle) {
        self.0 = (self.0 - rhs).normalized();
    }
}

impl AbsDiffEq for Lon {
    type Epsilon = Float;

    fn default_epsilon() -> Self::Epsilon {
        Float::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

impl Sub for Lon {
    type Output = LonInterval;

    fn sub(self, rhs: Self) -> Self::Output {
        LonInterval::new(rhs, self)
    }
}

/// A **closed** longitude interval represented by its lower and upper bounds (inclusive).
/// Note that the lower may be greater that the higher bound, in this case the interval contains
/// the point on the ante-meridian.
/// The lower and higher bound representations may also be inverted in special cases to allow the
/// empty and the full intervals as well:
/// - the empty interval is represented as the inverted interval `[pi..-pi]`
/// - the full interval is represented as the `[-pi..pi]`
/// - other intervals are represented as:
///     - `[lo..hi]` for lo != -pi
///     - `[pi..hi]` for lo == -pi and hi != pi,
///     - `[lo..hi]` otherwise
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct LonInterval(Interval);

impl LonInterval {
    /// Return an empty interval.
    pub fn empty() -> Self {
        Self(Interval::empty())
    }

    pub fn singleton(lon: Lon) -> Self {
        Self(Interval::singleton(lon.angle()))
    }

    /// Return a full interval.
    pub fn full() -> Self {
        Self(Interval::full())
    }

    /// Create a new interval.
    pub fn new(lo: Lon, hi: Lon) -> Self {
        Self(Interval::new(lo.angle(), hi.angle()))
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn is_full(&self) -> bool {
        self.0.is_full()
    }

    fn contains_antemeridian(&self) -> bool {
        self.0.is_inverted()
    }

    /// Return the **positive** length of this interval.
    /// The length is 0 if the interval is empty or a singleton.
    pub fn length(&self) -> Angle {
        self.0.length()
    }
}

/// [Lat] represents a latitude coordinate in [-pi/2..pi/2] radians.
///
/// # Creation
///
/// There are 2 ways to create a [Lat]:
/// - by using [Lat::new] passing an [Angle] value. The
/// value is **clamped** to [-pi/2 .. pi/2]:
/// ```
/// let lat = Lat::new(45. * DEG);
/// ```
///
/// - by using [Lat::dms] passing degrees/minutes/seconds values.
/// The value is also clamped to [-pi/2..pi/2]:
/// ```
/// lat n = Lat::dms(44.0, 43., 27.183);
/// ```
///
/// # Operations
///
/// [Lat] supports the following operations:
/// - Negation
/// - saturating addition and subtraction of [Angle].
/// - subtraction resulting in an [Angle].
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq, Neg, Display)]
#[display("{}", self.0.to_dms())]
pub struct Lat(Angle);

impl Lat {
    pub const MIN: Lat = Lat(Angle::M_PI_2);
    pub const ZERO: Lat = Lat(Angle::ZERO);
    pub const MAX: Lat = Lat(Angle::PI_2);

    /// Create a new latitude value.
    /// The angle is clamped in [-pi/2..pi/2]
    pub fn new(val: Angle) -> Self {
        Lat(val.clamped(Angle::M_PI_2, Angle::PI_2))
    }

    /// Create a new latitude value from a *dms* angle value.
    /// The angle value is converted to radians and clamped into [-pi/2, pi/2].
    ///
    /// # Parameters
    ///
    /// - `d`: the number of degrees. Determine the sign of the returned angle.
    /// - `m`: the number of minutes. Must be >= 0.
    /// - `s`: the number of seconds with fractional part. Must be >= 0.
    ///
    /// # Examples
    ///
    /// ```
    /// let a: Lat = Lat::dms(2., 20., 14.02500);
    /// assert!(Lat::dms(-12., 45., 59.1234).rad() < 0.0);
    /// ```
    pub fn dms(d: Float, m: Float, s: Float) -> Self {
        Self::new(Angle::dms(d, m, s))
    }

    /// Return the latitude angle.
    #[inline]
    pub fn angle(self) -> Angle {
        self.0
    }

    /// Return this latitude as a raw angle value in radians.
    #[inline]
    pub fn rad(self) -> Float {
        self.0.rad()
    }
}

impl Add<Angle> for Lat {
    type Output = Self;

    fn add(self, rhs: Angle) -> Self::Output {
        Self::new(self.0 + rhs)
    }
}

impl Add<Lat> for Angle {
    type Output = Lat;

    fn add(self, rhs: Lat) -> Self::Output {
        rhs + self
    }
}

impl AddAssign<Angle> for Lat {
    fn add_assign(&mut self, rhs: Angle) {
        self.0 = (self.0 + rhs).clamped(Angle::M_PI_2, Angle::PI_2)
    }
}
impl Sub<Angle> for Lat {
    type Output = Self;

    fn sub(self, rhs: Angle) -> Self::Output {
        Self::new(self.0 - rhs)
    }
}

impl SubAssign<Angle> for Lat {
    fn sub_assign(&mut self, rhs: Angle) {
        self.0 = (self.0 - rhs).clamped(Angle::M_PI_2, Angle::PI_2);
    }
}

impl AbsDiffEq for Lat {
    type Epsilon = Float;

    fn default_epsilon() -> Self::Epsilon {
        Float::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

// TODO: Add the following:
// - a LatInterval.
// - a Height similar to Length
// - a LonLat and LonLatHeight
// - a LonLatRect
// - a LonLatHeightRect

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

    use approx::assert_abs_diff_eq;

    use crate::cs::geodetic::{Lat, Lon};
    use crate::units::angle::{DEG, RAD};

    #[test]
    fn test_lon() {
        // wrapping
        assert_eq!(Lon::new(2.0 * PI * RAD), Lon::new(0.0 * RAD));
        assert_abs_diff_eq!(
            Lon::new(185.0 * DEG),
            Lon::new(-175.0 * DEG),
            epsilon = 1e-15
        );
        // normalization
        assert_eq!(Lon::new(-PI * RAD).normalize(), Lon::new(PI * RAD));
        // equality
        assert_eq!(Lon::new(FRAC_PI_2 * RAD), Lon::new(90.0 * DEG));
        assert_ne!(Lon::new(90. * DEG), Lon::new(91.0 * DEG));
        // ops
        assert_eq!(Lon::new(90. * DEG) + 45.0 * DEG, Lon::new(135.0 * DEG));
        assert_eq!(Lon::new(90. * DEG) - 45.0 * DEG, Lon::new(45.0 * DEG));
        // display
        assert_eq!(format!("{}", Lon::new(90. * DEG)), "  90° 00′ 00.00000000″");
    }

    #[test]
    fn test_lat() {
        // clamping
        assert_eq!(Lat::new(91. * DEG), Lat::MAX);
        assert_eq!(Lat::new(-90.01 * DEG), Lat::MIN);
        // equality
        assert_eq!(Lat::new(45. * DEG), Lat::new(FRAC_PI_4 * RAD));
        assert_ne!(Lat::new(0. * DEG), Lat::new(0.001 * DEG));
        // ops
        assert_eq!(Lat::new(45. * DEG + 10. * DEG), Lat::new(55. * DEG));
        assert_eq!(Lat::new(45. * DEG) + 50. * DEG, Lat::MAX);
        assert_eq!(Lat::new(45. * DEG) - 90. * DEG, Lat::new(-45. * DEG));
        assert_eq!(Lat::new(-45. * DEG) - 50. * DEG, Lat::MIN);
        // Display
        assert_eq!(format!("{}", Lat::new(45. * DEG)), "  45° 00′ 00.00000000″");
    }
}
