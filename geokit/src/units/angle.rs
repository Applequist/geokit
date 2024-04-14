//! This module defines some types used to express angles in various units.
use std::f64::consts;
use std::ops::{Add, Sub, Mul, Div};
use std::fmt::{Display, Formatter};

/// Trait implemented by types whose values measuring some angle
/// can be converted to radians.
pub trait Angle {

    /// Convert this angle value into a raw angle value expressed in radians.
    fn to_radians(&self) -> f64;
}

macro_rules! angle_unit {
    ( $name:ident, $unit:ident, $abbr:literal, $to_radian:expr) => {
        #[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
        pub struct $name(pub f64);

        impl $name {
            const ABBR: &'static str = $abbr;
        }

        impl Angle for $name {
            fn to_radians(&self) -> f64 {
                self.0 * ($to_radian as f64)
            }
        }

        impl Add for $name {
            type Output = $name;
            fn add(self, rhs: Self) -> Self::Output {
                Self(self.0 + rhs.0)
            }
        }

        impl Sub for $name {
            type Output = $name;
            fn sub(self, rhs: Self) -> Self::Output {
                Self(self.0 - rhs.0)
            }
        }

        impl Mul<f64> for $name {
            type Output = $name;
            fn mul(self, rhs: f64) -> Self::Output {
                Self(self.0 * rhs)
            }
        }

        impl Mul<$name> for f64 {
            type Output = $name;
            fn mul(self, rhs: $name) -> Self::Output {
                $name(self * rhs.0)
            }
        }

        impl Div<f64> for $name {
            type Output = $name;
            fn div(self, rhs: f64) -> Self::Output {
                $name(self.0 / rhs)
            }
        }

        impl Display for $name {
            fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
                write!(f, "{} {}", self.0, Self::ABBR)
            }
        }

        pub static $unit: $name = $name(1.0);
    };
}

angle_unit!(Radian, RAD, "rad", 1.0);
angle_unit!(Degree, DEG, "deg", consts::PI / 180.0);
angle_unit!(Gradian, GRAD, "grad", consts::PI / 200.0);
angle_unit!(Arcsec, SEC, "sec", consts::PI / 648_000.0);

#[cfg(test)]
mod tests {
    use std::f64::consts;
    use crate::units::angle::{Angle, Degree, Radian};

    #[test]
    fn test_to_radians() {
        assert_eq!(Radian(consts::FRAC_PI_2).to_radians(), consts::FRAC_PI_2);
        assert_eq!(Degree(90.0).to_radians(), consts::FRAC_PI_2);
    }

    #[test]
    fn test_ops() {
        assert_eq!(Degree(90.0), Degree(90.0));
        assert_ne!(Degree(90.0), Degree(91.0));
        assert_eq!(Degree(90.0) + Degree(90.0), Degree(180.0));
        assert_eq!(Degree(90.0) - Degree(45.0), Degree(45.0));
        assert_eq!(2.0 * Degree(90.0), Degree(180.0));
        assert_eq!(Degree(90.0) * 2.0, Degree(180.0));
        assert_eq!(Degree(90.0) / 2.0, Degree(45.0));
        assert!(Degree(90.0) < Degree(90.1));
    }
}