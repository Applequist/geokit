//! Provide a value type to work with [Azimuth].

use crate::quantities::angle::Angle;
use crate::units::angle::DEG;
use crate::{math::fp::Float, units::angle::RAD};
use approx::AbsDiffEq;
use derive_more::derive::Display;
use std::ops::{Add, Sub};

/// An [Azimuth] value represents a direction **in (-pi..pi] radians**, positive **clockwise** from North.
///
/// [Azimuth] supports the following operations:
/// - Subtraction: return the oriented angle between two azimuth (same orientation
/// as [Azimuth], eg positive clockwise.
/// - Addition and subtraction of [Angle] values. Return an [Azimuth]
#[derive(Debug, Copy, Clone, PartialEq, Display)]
#[display("{}", self.0.dms())]
pub struct Azimuth(Angle);

impl Azimuth {
    /// The azimuth pointing North.
    pub const NORTH: Azimuth = Azimuth(Angle::ZERO);
    /// The azimuth pointing East.
    pub const EAST: Azimuth = Azimuth(Angle::PI_2);
    /// The azimuth pointing South.
    pub const SOUTH: Azimuth = Azimuth(Angle::PI);
    /// The azimuth pointing West.
    pub const WEST: Azimuth = Azimuth(Angle::M_PI_2);

    /// Creates a new azimuth value.
    ///
    /// The given angle `val` is normalized into (-pi, pi].
    pub fn new(val: Angle) -> Self {
        Self(val.normalized())
    }

    /// Creates a [Azimuth] with the given angle **in radians**.
    ///
    /// The given angle is normalized into (-pi, pi].
    pub fn rad(val_rad: Float) -> Self {
        Self::new(val_rad * RAD)
    }

    /// Creates a [Azimuth] with the given angle **in degrees**.
    ///
    /// The given angle is converted to radians and normalized in (-pi, pi]
    pub fn deg(val_deg: Float) -> Self {
        Self::new(val_deg * DEG)
    }

    /// Creates a new azimuth value from a *dms* angle value.
    ///
    /// The angle value is converted to radians and wrapped into (-pi, pi].
    ///
    /// # Parameters
    ///
    /// - `d`: the number of degrees. Determine the sign of the returned angle.
    /// - `m`: the number of minutes. Must be >= 0 and <= 59.
    /// - `s`: the number of seconds with fractional part. Must be >= 0 and <= 60.
    pub fn dms(d: Float, m: Float, s: Float) -> Self {
        debug_assert!(d > -180. && m <= 180., "degrees must be in (-180..180]");
        debug_assert!(m >= 0. && m <= 59.0, "minutes must be in [0..59]");
        debug_assert!(s >= 0. && s < 60., "seconds must be in [0..60)");
        let f = d.signum();
        let deg = d.abs() + m / 60. + s / 3600.;
        Self::new(f * deg * DEG)
    }

    /// Return the azimtuh value as an [Angle]
    #[inline]
    pub fn angle(self) -> Angle {
        self.0
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

    /// Returns an [Azimuth] whose angle is the **normalized** angle `self.angle() + rhs`.
    ///
    /// # Example
    ///
    /// ```
    /// # use geokit::cs::azimuth::Azimuth;
    /// # use geokit::quantities::angle::Angle;
    /// # use geokit::units::angle::DEG;
    ///
    /// assert_eq!(Azimuth::EAST + 180. * DEG, Azimuth::WEST);
    /// ```
    fn add(self, rhs: Angle) -> Self::Output {
        Azimuth::new(self.0 + rhs)
    }
}

impl Add<Azimuth> for Angle {
    type Output = Azimuth;

    /// Returns an [Azimuth] whose angle is the **normalized** angle `self + rhs.angle()`.
    fn add(self, rhs: Azimuth) -> Self::Output {
        rhs + self
    }
}

impl Sub for Azimuth {
    type Output = Angle;

    /// Returns the **oriented** [Angle] in (-pi, pi] radians from `rhs` to `self`.
    ///
    /// The returned angle has the same orientation as [Azimuth].
    ///
    /// # Examples
    ///
    /// ```
    /// # use geokit::cs::azimuth::Azimuth;
    /// # use geokit::quantities::angle::Angle;
    /// # use geokit::units::angle::DEG;
    ///
    /// assert_eq!(Azimuth::EAST - Azimuth::NORTH, 90. * DEG);
    /// assert_eq!(Azimuth::WEST - Azimuth::NORTH, -90. * DEG);
    ///
    /// ```
    fn sub(self, rhs: Self) -> Self::Output {
        rhs.0.diff_to(self.0).normalized()
    }
}

impl AbsDiffEq for Azimuth {
    type Epsilon = Angle;

    fn default_epsilon() -> Self::Epsilon {
        Angle::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        (*self - *other).abs() <= epsilon
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

        assert_abs_diff_eq!(se - ne, 140. * DEG);
        assert_abs_diff_eq!(sw - ne, 180. * DEG);
        assert_abs_diff_eq!(nw - ne, -40. * DEG);

        assert_abs_diff_eq!(sw - se, 40. * DEG);
        assert_abs_diff_eq!(nw - se, 180. * DEG);
        assert_abs_diff_eq!(ne - se, -140. * DEG);

        assert_abs_diff_eq!(nw - sw, 140. * DEG);
        assert_abs_diff_eq!(ne - sw, 180. * DEG);
        assert_abs_diff_eq!(se - sw, -40. * DEG);

        assert_abs_diff_eq!(ne - nw, 40. * DEG);
        assert_abs_diff_eq!(se - nw, 180. * DEG);
        assert_abs_diff_eq!(sw - nw, -140. * DEG);
    }
}
