//! Provide a value type to work with [Azimuth].

use crate::math::Float;
use crate::quantities::angle::Angle;
use approx::AbsDiffEq;
use derive_more::derive::Display;
use std::ops::{Add, Sub};

/// An [Azimuth] value represents a direction **in (-pi..pi] radians**, positive **clockwise** from North.
///
/// # Creation
///
/// There are 2 ways to create an [Azimuth] value:
/// 1. By using [Azimuth::new], passing an [Angle] value,
/// ```
/// use geokit::units::angle::DEG;
/// use geokit::cs::azimuth::Azimuth;
/// let az = Azimuth::new(33. * DEG);
/// ```
/// 2. Or by specifying DMS values:
/// ```
/// use geokit::cs::azimuth::Azimuth;
/// let az = Azimuth::dms(-12., 34., 56.123);
/// ```
/// In both cases, the angle value is wrapped into (-pi, pi] so that
/// any direction has unique [Azimuth] value.
///
/// # Operations
///
/// [Azimuth] supports the following operations:
/// - Subtraction: return the oriented angle between twa azimuth (same orientation
/// as [Azimuth], eg positive clockwise.
/// - Addition and subtraction of [Angle] values. Return an [Azimuth]
#[derive(Debug, Copy, Clone, PartialEq, Display)]
#[display("{}", self.0.to_dms())]
pub struct Azimuth(Angle);

impl Azimuth {
    pub const NORTH: Azimuth = Azimuth(Angle::ZERO);
    pub const EAST: Azimuth = Azimuth(Angle::PI_2);
    pub const SOUTH: Azimuth = Azimuth(Angle::PI);
    pub const WEST: Azimuth = Azimuth(Angle::M_PI_2);

    /// Create a new azimuth value.
    ///
    /// The given angle `val` is wrapped into (-pi, pi].
    pub fn new(val: Angle) -> Self {
        Self(val.normalized())
    }

    /// Create a new azimuth value from a *dms* angle value.
    /// The angle value is converted to radians and wrapped into (-pi, pi].
    ///
    /// # Parameters
    ///
    /// - `d`: the number of degrees. Determine the sign of the returned angle.
    /// - `m`: the number of minutes. Must be >= 0 and <= 59.
    /// - `s`: the number of seconds with fractional part. Must be >= 0 and <= 60.
    pub fn dms(d: Float, m: Float, s: Float) -> Self {
        Self::new(Angle::dms(d, m, s))
    }

    /// Return the azimtuh value as an [Angle]
    #[inline]
    pub fn angle(self) -> Angle {
        self.0
    }

    /// Return this azimuth as a raw angle value (-pi..pi] radians.
    #[inline]
    pub fn rad(self) -> Float {
        self.0.rad()
    }

    #[inline]
    pub fn sin(self) -> Float {
        self.0.sin()
    }

    #[inline]
    pub fn cos(self) -> Float {
        self.0.cos()
    }

    #[inline]
    pub fn sin_cos(self) -> (Float, Float) {
        self.0.sin_cos()
    }

    #[inline]
    pub fn tan(self) -> Float {
        self.0.tan()
    }
}

impl Add<Angle> for Azimuth {
    type Output = Azimuth;

    fn add(self, rhs: Angle) -> Self::Output {
        Azimuth::new(self.0 + rhs)
    }
}

impl Add<Azimuth> for Angle {
    type Output = Azimuth;

    fn add(self, rhs: Azimuth) -> Self::Output {
        rhs + self
    }
}

impl Sub for Azimuth {
    type Output = Angle;

    fn sub(self, rhs: Self) -> Self::Output {
        rhs.0.diff_to(self.0).normalized()
    }
}

impl AbsDiffEq for Azimuth {
    type Epsilon = Float;

    fn default_epsilon() -> Self::Epsilon {
        Float::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

#[cfg(test)]
mod test {

    use super::Azimuth;
    use crate::units::angle::DEG;
    use approx::assert_abs_diff_eq;

    #[test]
    fn subtraction() {
        let ne = Azimuth::new(10. * DEG);
        let se = Azimuth::new(150. * DEG);
        let sw = Azimuth::new(-170. * DEG);
        let nw = Azimuth::new(-30. * DEG);

        assert_abs_diff_eq!(se - ne, 140. * DEG, epsilon = 1e-15);
        assert_abs_diff_eq!(sw - ne, 180. * DEG, epsilon = 1e-15);
        assert_abs_diff_eq!(nw - ne, -40. * DEG, epsilon = 1e-15);

        assert_abs_diff_eq!(sw - se, 40. * DEG, epsilon = 1e-15);
        assert_abs_diff_eq!(nw - se, 180. * DEG, epsilon = 1e-15);
        assert_abs_diff_eq!(ne - se, -140. * DEG, epsilon = 1e-15);

        assert_abs_diff_eq!(nw - sw, 140. * DEG, epsilon = 1e-15);
        assert_abs_diff_eq!(ne - sw, 180. * DEG, epsilon = 1e-15);
        assert_abs_diff_eq!(se - sw, -40. * DEG, epsilon = 1e-15);

        assert_abs_diff_eq!(ne - nw, 40. * DEG, epsilon = 1e-15);
        assert_abs_diff_eq!(se - nw, 180. * DEG, epsilon = 1e-15);
        assert_abs_diff_eq!(sw - nw, -140. * DEG, epsilon = 1e-15);
    }

    #[test]
    fn turn_to() {}
}
