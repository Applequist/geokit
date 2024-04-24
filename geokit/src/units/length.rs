//! This module defines some types used to express length in various units.
use std::fmt::{Display, Formatter};
use std::ops::{Add, Div, Mul, Sub};

macro_rules! length_unit {
    ( $name:ident, $unit:ident, $abbr:literal, $to_meters:expr) => {
        #[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
        pub struct $name(pub f64);

        impl $name {
            pub const ABBR: &'static str = $abbr;
            pub const M_PER_UNIT: f64 = $to_meters;

            /// Returns the unit value, aka how many meters per unit
            #[inline]
            pub fn unit() -> f64 {
                Self::M_PER_UNIT
            }

            #[inline]
            pub fn m(self) -> f64 {
                self.0 * Self::M_PER_UNIT
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

        impl Div<$name> for $name {
            type Output = f64;
            #[inline]
            fn div(self, rhs: $name) -> Self::Output {
                self.0 / rhs.0
            }
        }

        impl Display for $name {
            fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
                write!(f, "{} {}", self.0, Self::ABBR)
            }
        }

        pub const $unit: $name = $name(1.);
    };
}

length_unit!(Meters, M, "m", 1.0);

/// A trait implemented by types whose values measuring some length
/// can be converted to a *raw* value in meter.
pub trait Length {
    /// Convert this length value into a raw length value expressed in meters.
    fn to_meters(self) -> Meters;
}

macro_rules! impl_length {
    ($($name:ident),*) => {
        $(impl Length for $name {
            fn to_meters(self) -> Meters {
                Meters(self.m())
            }

        })*
    };
}

impl_length!(Meters);

#[cfg(test)]
mod tests {
    use crate::units::length::Meters;

    #[test]
    fn test_to_meters() {
        assert_eq!(Meters(1.5).m(), 1.5);
    }
    #[test]
    fn test_ops() {
        assert_eq!(Meters(1.), Meters(1.));
        assert_ne!(Meters(1.), Meters(1.1));
        assert_eq!(Meters(1.) + Meters(2.), Meters(3.));
        assert_eq!(Meters(2.) - Meters(1.), Meters(1.));
        assert_eq!(2. * Meters(1.), Meters(2.));
        assert_eq!(Meters(1.5) * 2.0, Meters(3.));
        assert_eq!(Meters(4.0) / 2.0, Meters(2.));
        assert_eq!(Meters(4.0) / Meters(2.), 2.0);
        assert!(Meters(1.) < Meters(1.1));
    }
}
