//! This module defines some types used to express angles in various units.
use std::f64::consts;
use std::fmt::{Display, Formatter};
use std::ops::{Add, Div, Mul, Sub};

macro_rules! angle_unit {
    ( $name:ident, $unit:ident, $abbr:literal, $to_radians:expr) => {
        #[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
        pub struct $name(pub f64);

        impl $name {
            pub const ABBR: &'static str = $abbr;
            pub const RAD_PER_UNIT: f64 = $to_radians;

            /// Returns the unit value, aka how many radians per unit
            #[inline]
            pub fn unit() -> f64 {
                Self::RAD_PER_UNIT
            }

            /// Returns the raw angle value **in radians**.
            #[inline]
            pub fn rad(self) -> f64 {
                self.0 * Self::RAD_PER_UNIT
            }
        }

        impl Add for $name {
            type Output = $name;
            #[inline]
            fn add(self, rhs: Self) -> Self::Output {
                Self(self.0 + rhs.0)
            }
        }

        impl Sub for $name {
            type Output = $name;
            #[inline]
            fn sub(self, rhs: Self) -> Self::Output {
                Self(self.0 - rhs.0)
            }
        }

        impl Mul<f64> for $name {
            type Output = $name;
            #[inline]
            fn mul(self, rhs: f64) -> Self::Output {
                Self(self.0 * rhs)
            }
        }

        impl Mul<$name> for f64 {
            type Output = $name;
            #[inline]
            fn mul(self, rhs: $name) -> Self::Output {
                $name(self * rhs.0)
            }
        }

        impl Div<f64> for $name {
            type Output = $name;
            #[inline]
            fn div(self, rhs: f64) -> Self::Output {
                $name(self.0 / rhs)
            }
        }

        impl Display for $name {
            fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
                write!(f, "{} {}", self.0, Self::ABBR)
            }
        }

        pub const $unit: $name = $name(1.0);
    };
}

angle_unit!(Radians, RAD, "rad", 1.0);
angle_unit!(Degrees, DEG, "deg", consts::PI / 180.0);
angle_unit!(Gradians, GRAD, "grad", consts::PI / 200.0);
angle_unit!(Arcsecs, SEC, "sec", consts::PI / 648_000.0);

/// Trait implemented by types whose values measuring some angle
/// can be converted to radians.
pub trait Angle {
    /// Convert this angle value into a raw angle value expressed in radians.
    fn to_radians(self) -> Radians;
}

macro_rules! angle_impl {
    ( $($name:ident),* ) => {
        $(impl Angle for $name {
            fn to_radians(self) -> Radians {
                Radians(self.rad())
            }
        })*
    };
}

angle_impl!(Radians, Degrees, Gradians, Arcsecs);

#[cfg(test)]
mod tests {
    use crate::units::angle::{Angle, Degrees, Radians};
    use std::f64::consts;

    #[test]
    fn test_to_radians() {
        assert_eq!(Radians(consts::FRAC_PI_2).rad(), consts::FRAC_PI_2);
        assert_eq!(Degrees(90.0).to_radians(), Radians(consts::FRAC_PI_2));
    }

    #[test]
    fn test_ops() {
        assert_eq!(Degrees(90.0), Degrees(90.0));
        assert_ne!(Degrees(90.0), Degrees(91.0));
        assert_eq!(Degrees(90.0) + Degrees(90.0), Degrees(180.0));
        assert_eq!(Degrees(90.0) - Degrees(45.0), Degrees(45.0));
        assert_eq!(2.0 * Degrees(90.0), Degrees(180.0));
        assert_eq!(Degrees(90.0) * 2.0, Degrees(180.0));
        assert_eq!(Degrees(90.0) / 2.0, Degrees(45.0));
        assert!(Degrees(90.0) < Degrees(90.1));
    }
}
