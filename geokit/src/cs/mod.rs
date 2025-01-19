use approx::AbsDiffEq;
use num::Zero;
use std::f64::consts::PI;
use std::fmt::{Display, Formatter};

use crate::quantity::angle::units::Deg;
use crate::quantity::angle::{formatters::DMS, Angle};

/// An azimuth direction **in (-pi..pi] radians**, positive **clockwise** from North.
#[derive(Debug, Copy, Clone, PartialEq, Default)]
pub struct Azimuth(f64);

impl Azimuth {
    pub const NORTH: Azimuth = Azimuth(0.0);
    pub const EAST: Azimuth = Azimuth(90.0 * Deg::UNIT.rad_per_unit());
    pub const SOUTH: Azimuth = Azimuth(180. * Deg::UNIT.rad_per_unit());
    pub const WEST: Azimuth = Azimuth(-90.0 * Deg::UNIT.rad_per_unit());

    pub fn new(val: Angle) -> Self {
        let mut az = val.wrap().rad();
        if az <= -PI {
            az = PI;
        }
        Self(az)
    }

    pub fn zero() -> Self {
        Azimuth(0.0)
    }

    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
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

impl Display for Azimuth {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", DMS::from_rad(self.0))
    }
}

pub mod geodetic;
