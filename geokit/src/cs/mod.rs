use approx::AbsDiffEq;

use crate::units::angle::Angle;

/// An azimuth direction starting at 0 pointing in the geographic North direction and
/// positive **clockwise**.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Azimuth(f64);

impl Azimuth {
    pub fn new<A: Angle>(val: A) -> Self {
        Self(val.to_radians().wrap(0.))
    }

    /// Return this azimuth as a raw angle value [-pi..pi] radians.
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
