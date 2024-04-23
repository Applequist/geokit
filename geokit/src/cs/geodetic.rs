use crate::units::angle::{Angle, Radians};
use std::f64::consts::PI;
use std::fmt::{Display, Formatter};
use std::ops::{Add, Sub};

/// A longitude value in [-pi..pi] radians.
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct Lon(f64);

impl Lon {
    /// Create a new longitude value with a given angle.
    /// The angle is wrapped into [-pi..pi]
    pub fn new<U: Angle>(val: U) -> Self {
        Self(Self::rem_two_pi(val.to_radians().0))
    }

    /// Normalize the longitude into (-pi..pi]
    pub fn normalize(self) -> Self {
        let mut lon = self.0;
        if lon <= -PI {
            lon = PI;
        }
        Self(lon)
    }

    /// Warning: hacky code!
    /// This function aims at reproducing the C/C++ remainder function.
    /// It is mainly used to normalize longitude in [-pi, pi].
    fn rem_two_pi(val: f64) -> f64 {
        let mut na = val;
        let q = (na / (2.0 * PI)).trunc();
        na -= q * 2.0 * PI;
        while na < -PI {
            na += 2.0 * PI
        }
        while na > PI {
            na -= 2.0 * PI
        }
        na
    }

    pub fn rad(self) -> f64 {
        self.0
    }
}

impl Angle for Lon {

    fn to_radians(self) -> Radians {
        Radians(self.0)
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

impl Display for Lon {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} rad", self.0)
    }
}

#[cfg(test)]
mod tests {
    use crate::cs::geodetic::Lon;
    use crate::units::angle::{Degrees, Radians, DEG, RAD};
    use std::f64::consts;
    use std::f64::consts::PI;

    #[test]
    fn test_wrapping() {
        assert_eq!(Lon::new(Radians(2.0 * consts::PI)), Lon::new(Radians(0.0)));
        assert_eq!(Lon::new(185.0 * DEG), Lon::new(Degrees(-175.0)));
    }

    #[test]
    fn test_normalize() {
        assert_eq!(Lon::new(-PI * RAD).normalize(), Lon::new(PI * RAD));
    }

    #[test]
    fn test_ops() {
        assert_eq!(Lon::new(consts::FRAC_PI_2 * RAD), Lon::new(90.0 * DEG));
        assert_ne!(Lon::new(90. * DEG), Lon::new(91.0 * DEG));
        assert_eq!(Lon::new(90. * DEG) + 45.0 * DEG, Lon::new(135.0 * DEG));
        assert_eq!(Lon::new(90. * DEG) - 45.0 * DEG, Lon::new(45.0 * DEG));
    }
}
