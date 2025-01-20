use crate::quantity::angle::Angle;
use approx::AbsDiffEq;
use derive_more::derive::Display;

/// An azimuth direction **in (-pi..pi] radians**, positive **clockwise** from North.
#[derive(Debug, Copy, Clone, PartialEq, Default, Display)]
#[display("{}", self.0.to_dms())]
pub struct Azimuth(Angle);

impl Azimuth {
    pub const NORTH: Azimuth = Azimuth(Angle::zero());
    pub const EAST: Azimuth = Azimuth(Angle::PI_2);
    pub const SOUTH: Azimuth = Azimuth(Angle::PI);
    pub const WEST: Azimuth = Azimuth(Angle::M_PI_2);

    pub fn new(val: Angle) -> Self {
        let mut az = val.wrap();
        if az <= Angle::M_PI {
            az = Angle::PI;
        }
        Self(az)
    }

    /// Create a new azimuth value from a *dms* angle value.
    /// The angle value is converted to radians and wrapped into (-pi, pi].
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
    pub fn dms(d: f64, m: f64, s: f64) -> Self {
        Self::new(Angle::dms(d, m, s))
    }

    #[inline]
    pub const fn zero() -> Self {
        Azimuth(Angle::zero())
    }

    #[inline]
    pub fn angle(self) -> Angle {
        self.0
    }

    /// Return this azimuth as a raw angle value [-pi..pi] radians.
    #[inline]
    pub fn rad(self) -> f64 {
        self.0.rad()
    }
}

impl AbsDiffEq for Azimuth {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

pub mod geodetic;
