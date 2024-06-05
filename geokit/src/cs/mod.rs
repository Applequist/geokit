use approx::AbsDiffEq;
use std::fmt::{Display, Formatter};

use crate::quantity::angle::{wrap, DMS};

/// An azimuth direction in radians, positive **clockwise** from north.
#[derive(Debug, Copy, Clone, PartialEq, Default)]
pub struct Azimuth(f64);

impl Azimuth {
    pub fn new(val: f64) -> Self {
        Self(wrap(val, 0.))
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
