use approx::AbsDiffEq;
use num::traits::WrappingAdd;
use num::Bounded;
use std::f64::consts::{FRAC_PI_2, PI};
use std::fmt::{Display, Formatter};
use std::ops::{Add, Neg, Sub};

use crate::units::angle::{Angle, Radians};

/// A longitude coordinate in [-pi..pi] radians.
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct Lon(f64);

impl Lon {
    pub const MIN: Lon = Lon(-PI);
    pub const MAX: Lon = Lon(PI);

    /// Create a new longitude value with a given angle.
    /// The angle is converted to radians and wrapped into [-pi..pi].
    pub fn new<U: Angle>(val: U) -> Self {
        Self(val.to_radians().wrap(-PI))
    }

    /// Normalize the longitude into (-pi..pi].
    pub fn normalize(self) -> Self {
        if self <= Self::MIN {
            Self::MAX
        } else {
            self
        }
    }

    /// Return the longitude as a raw angle value **in radians**.
    #[inline]
    pub fn rad(self) -> f64 {
        self.0
    }
}

impl<T> From<T> for Lon
where
    T: Angle,
{
    fn from(value: T) -> Self {
        Self::new(value)
    }
}

impl<U> Add<U> for Lon
where
    U: Angle,
{
    type Output = Self;

    fn add(self, rhs: U) -> Self::Output {
        // Watch out for infinite recursion with self + rhs
        Self::new(Radians(self.0) + rhs.to_radians())
    }
}

impl<U> Sub<U> for Lon
where
    U: Angle,
{
    type Output = Self;

    fn sub(self, rhs: U) -> Self::Output {
        // Watch out for infinite recursion with self - rhs
        Self::new(Radians(self.0) - rhs.to_radians())
    }
}

impl AbsDiffEq for Lon {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

impl Display for Lon {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "lon = {} rad", self.0)
    }
}

/// A latitude coordinate in [-pi/2..pi/2]  radians.
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct Lat(f64);

impl Lat {
    pub const MIN: Lat = Lat(-FRAC_PI_2);
    pub const MAX: Lat = Lat(FRAC_PI_2);

    /// Create a new latitude value with the given angle.
    /// The angle value is converted to radians and clamped into [-pi/2..pi/2].
    pub fn new<A: Angle>(val: A) -> Self {
        Lat(val.to_radians().clamp(-FRAC_PI_2, FRAC_PI_2))
    }

    /// Return this latitude as a raw angle value in radians.
    #[inline]
    pub fn rad(self) -> f64 {
        self.0
    }
}

impl<U> From<U> for Lat
where
    U: Angle,
{
    fn from(value: U) -> Self {
        Self::new(value)
    }
}

impl<U> Add<U> for Lat
where
    U: Angle,
{
    type Output = Self;

    fn add(self, rhs: U) -> Self::Output {
        Self::new(Radians(self.0) + rhs.to_radians())
    }
}

impl<U> Sub<U> for Lat
where
    U: Angle,
{
    type Output = Self;

    fn sub(self, rhs: U) -> Self::Output {
        Self::new(Radians(self.0) - rhs.to_radians())
    }
}

impl AbsDiffEq for Lat {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

impl Display for Lat {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "lat = {} rad", self.0)
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts;
    use std::f64::consts::{FRAC_PI_4, PI};

    use crate::cs::geodetic::{Lat, Lon};
    use crate::units::angle::{Degrees, Radians, DEG, RAD};

    #[test]
    fn test_lon() {
        // wrapping
        assert_eq!(Lon::new(Radians(2.0 * consts::PI)), Lon::new(Radians(0.0)));
        assert_eq!(Lon::new(185.0 * DEG), Lon::new(Degrees(-175.0)));
        // normalization
        assert_eq!(Lon::new(-PI * RAD).normalize(), Lon::new(PI * RAD));
        // equality
        assert_eq!(Lon::new(consts::FRAC_PI_2 * RAD), Lon::new(90.0 * DEG));
        assert_ne!(Lon::new(90. * DEG), Lon::new(91.0 * DEG));
        // ops
        assert_eq!(Lon::new(90. * DEG) + 45.0 * DEG, Lon::new(135.0 * DEG));
        assert_eq!(Lon::new(90. * DEG) - 45.0 * DEG, Lon::new(45.0 * DEG));
        // display
        assert_eq!(
            format!("{}", Lon::new(90. * DEG)),
            "lon = 1.5707963267948966 rad"
        );
    }

    #[test]
    fn test_lat() {
        // clamping
        assert_eq!(Lat::new(Degrees(91.)), Lat::MAX);
        assert_eq!(Lat::new(Degrees(-90.01)), Lat::MIN);
        // equality
        assert_eq!(Lat::new(Degrees(45.)), Lat::new(Radians(FRAC_PI_4)));
        assert_ne!(Lat::new(Degrees(0.)), Lat::new(Degrees(0.001)));
        // ops
        assert_eq!(
            Lat::new(Degrees(45.)) + Degrees(10.),
            Lat::new(Degrees(55.))
        );
        assert_eq!(Lat::new(Degrees(45.)) + Degrees(50.), Lat::MAX);
        assert_eq!(
            Lat::new(Degrees(45.)) - Degrees(90.),
            Lat::new(Degrees(-45.))
        );
        assert_eq!(Lat::new(Degrees(-45.)) - Degrees(50.), Lat::MIN);
        // Display
        assert_eq!(
            format!("{}", Lat::new(Degrees(45.))),
            "lat = 0.7853981633974483 rad"
        );
    }
}
