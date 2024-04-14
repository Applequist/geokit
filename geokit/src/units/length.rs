//! This module defines some types used to express length in various units.
use std::fmt::{Display, Formatter};
use std::ops::{Add, Sub, Mul, Div};

/// A trait implemented by types whose values measuring some length
/// can be converted to a *raw* value in meter.
pub trait Length<T = f64> {
    /// Convert this length value into a raw length value expressed in meters.
    fn to_meters(&self) -> T;
}

macro_rules! length_unit {
    ( $name:ident, $unit:ident, $abbr:literal, $to_meter:expr) => {
        #[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
        pub struct $name(pub f64);

        impl $name {
            const ABBR: &'static str = $abbr;
        }

        impl Length for $name {
            fn to_meters(&self) -> f64 {
                self.0 * ($to_meter as f64)
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

        impl Div<$name> for $name {
            type Output = f64;
            fn div(self, rhs: $name) -> Self::Output {
                self.0 / rhs.0
            }
        }

        impl Display for $name {
            fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
                write!(f, "{} {}", self.0, Self::ABBR)
            }
        }

        pub static $unit: $name = $name(1.);
    };
}

length_unit!(Meter, M, "m", 1.0);

#[cfg(test)]
mod tests {
    use crate::units::length::{Meter, Length};

    #[test]
    fn test_to_meters() {
        assert_eq!(Meter(1.5).to_meters(), 1.5);
    }
    #[test]
    fn test_ops() {
        assert_eq!(Meter(1.), Meter(1.));
        assert_ne!(Meter(1.), Meter(1.1));
        assert_eq!(Meter(1.) + Meter(2.), Meter(3.));
        assert_eq!(Meter(2.) - Meter(1.), Meter(1.));
        assert_eq!(2. * Meter(1.), Meter(2.));
        assert_eq!(Meter(1.5) * 2.0, Meter(3.));
        assert_eq!(Meter(4.0) / 2.0, Meter(2.));
        assert_eq!(Meter(4.0) / Meter(2.), 2.0);
        assert!(Meter(1.) < Meter(1.1));
    }

}