use approx::AbsDiffEq;
use derive_more::derive::Display;
use num::Zero;
use std::f64::consts::PI;

use crate::quantity::angle::units::{DEG, RAD};
use crate::quantity::angle::Angle;

/// An azimuth direction **in (-pi..pi] radians**, positive **clockwise** from North.
#[derive(Debug, Copy, Clone, PartialEq, Default, Display)]
#[display("{}", self.angle().to_dms())]
pub struct Azimuth(f64);

impl Azimuth {
    pub const NORTH: Azimuth = Azimuth(0.0);
    pub const EAST: Azimuth = Azimuth(90.0 * DEG.rad_per_unit());
    pub const SOUTH: Azimuth = Azimuth(180. * DEG.rad_per_unit());
    pub const WEST: Azimuth = Azimuth(-90.0 * DEG.rad_per_unit());

    pub fn new(val: Angle) -> Self {
        let mut az = val.wrap().rad();
        if az <= -PI {
            az = PI;
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

    pub fn zero() -> Self {
        Azimuth(0.0)
    }

    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    pub fn angle(self) -> Angle {
        self.0 * RAD
    }

    /// Return this azimuth as a raw angle value [-pi..pi] radians.
    #[inline]
    pub fn rad(self) -> f64 {
        self.0
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
