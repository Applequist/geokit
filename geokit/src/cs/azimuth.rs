//! Provide a value type to work with [Azimuth].

use crate::math::Float;

use super::s1::Angle;
use approx::AbsDiffEq;
use derive_more::derive::Display;
use std::ops::Add;

/// An [Azimuth] value represents a direction **in (-pi..pi] radians**, positive **clockwise** from North.
///
/// # Creation
///
/// There are 2 ways to create an [Azimuth] value:
/// 1. By using [Azimuth::new], passing an [Angle] value,
/// ```
/// let az = Azimuth::new(33. * DEG);
/// ```
/// 2. Or by specifying DMS values:
/// ```
/// let az = Azimuth::dms(-12., 34., 56.123);
/// ```
/// In both cases, the angle value is wrapped into (-pi, pi] so that
/// any direction has unique [Azimuth] value.
///
/// # Operations
///
/// [Azimuth] supports the following operations:
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

    /// Returns the smallest turn from this azimuth to the `other` azimuth.
    fn turn_to(self, other: Self) -> ToAz {
        // convert from (-pi..pi] to [0 2pi) clockwise
        let start = if self.0 < Angle::ZERO {
            self.0 + Angle::TWO_PI
        } else {
            self.0
        };

        let end = if other.0 < Angle::ZERO {
            other.0 + Angle::TWO_PI
        } else {
            other.0
        };

        let mut delta = end - start;
        if delta > Angle::PI {
            delta -= Angle::TWO_PI;
        } else if delta < Angle::M_PI {
            delta += Angle::TWO_PI;
        }

        if (delta.abs() - Angle::PI).abs().rad() < 1e-15 {
            ToAz::Antipodal
        } else {
            // NOTE: Azimuth is positive CW hence the inversion
            // of case
            if delta < Angle::ZERO {
                ToAz::Ccw(delta.abs())
            } else {
                ToAz::Cw(delta.abs())
            }
        }
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

impl AbsDiffEq for Azimuth {
    type Epsilon = Float;

    fn default_epsilon() -> Self::Epsilon {
        Float::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Display)]
pub enum ToAz {
    Cw(Angle),
    Ccw(Angle),
    Antipodal,
}

impl AbsDiffEq for ToAz {
    type Epsilon = Float;

    fn default_epsilon() -> Self::Epsilon {
        Float::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        match (self, other) {
            (ToAz::Antipodal, ToAz::Antipodal) => true,
            (ToAz::Cw(l), ToAz::Cw(r)) => l.abs_diff_eq(r, epsilon),
            (ToAz::Ccw(l), ToAz::Ccw(r)) => l.abs_diff_eq(r, epsilon),
            _ => false,
        }
    }
}

#[cfg(test)]
mod test {

    use approx::assert_abs_diff_eq;

    use crate::{cs::azimuth::ToAz, units::angle::DEG};

    use super::Azimuth;

    #[test]
    fn subtraction() {
        let ne = Azimuth::new(10. * DEG);
        let se = Azimuth::new(150. * DEG);
        let sw = Azimuth::new(-170. * DEG);
        let nw = Azimuth::new(-30. * DEG);

        assert_abs_diff_eq!(ne.turn_to(se), ToAz::Cw(140. * DEG), epsilon = 1e-15);
        assert_abs_diff_eq!(ne.turn_to(sw), ToAz::Antipodal, epsilon = 1e-15);
        assert_abs_diff_eq!(ne.turn_to(nw), ToAz::Ccw(40. * DEG));

        assert_abs_diff_eq!(se.turn_to(sw), ToAz::Cw(40. * DEG), epsilon = 1e-15);
        assert_abs_diff_eq!(se.turn_to(nw), ToAz::Antipodal, epsilon = 1e-15);
        assert_abs_diff_eq!(se.turn_to(ne), ToAz::Ccw(140. * DEG), epsilon = 1e-15);

        assert_abs_diff_eq!(sw.turn_to(nw), ToAz::Cw(140. * DEG), epsilon = 1e-15);
        assert_abs_diff_eq!(sw.turn_to(ne), ToAz::Antipodal, epsilon = 1e-15);
        assert_abs_diff_eq!(sw.turn_to(se), ToAz::Ccw(40. * DEG), epsilon = 1e-15);

        assert_abs_diff_eq!(nw.turn_to(ne), ToAz::Cw(40. * DEG), epsilon = 1e-15);
        assert_abs_diff_eq!(nw.turn_to(se), ToAz::Antipodal, epsilon = 1e-15);
        assert_abs_diff_eq!(nw.turn_to(sw), ToAz::Ccw(140. * DEG), epsilon = 1e-15);
    }
}
