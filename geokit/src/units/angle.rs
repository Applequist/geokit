//! This module defines some types used to express angles in various units.
use std::f64::consts;
use std::f64::consts::{FRAC_PI_2, PI};
use std::fmt::{Display, Formatter};
use std::ops::{Add, Div, Mul, Neg, Sub};

use approx::AbsDiffEq;

/// Trait implemented by types whose values are measuring some angle
/// and can be converted to the same angle expressed in [Radians](radians).
pub trait Angle {
    /// Convert this angle value into [Radians].
    fn to_radians(self) -> Radians
    where
        Self: Sized,
    {
        Radians(self.rad())
    }

    /// Convert this angle value into a raw angle value **in radians**.
    fn rad(self) -> f64;
}

macro_rules! angle_unit {
    ( $(#[$outer:meta])* ($name:ident, $unit:ident, $abbr:literal, $to_radians:expr) ) => {
        $(#[$outer])*
        #[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
        pub struct $name(pub f64);

        impl $name {
            /// The unit abbreviation
            pub const ABBR: &'static str = $abbr;
            /// How many radians is 1 unit
            pub const RAD_PER_UNIT: f64 = $to_radians;

            /// Create a new angle from the given one effectively converting
            /// the angle value from the given angle's unit to this unit.
            pub fn from_angle<T: Angle>(angle: T) -> Self {
                Self::from_rad(angle.rad())
            }

            /// Create a new angle measured in this unit from a raw angle value **in radians**.
            pub fn from_rad(rad: f64) -> Self {
                Self(rad / Self::RAD_PER_UNIT)
            }

            /// Returns the unit value, aka how many radians per unit
            #[inline]
            pub fn unit() -> Radians {
                Radians(Self::RAD_PER_UNIT)
            }
        }

        impl Angle for $name {
            #[inline]
            fn rad(self) -> f64 {
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

        impl Neg for $name {
            type Output = $name;
            #[inline]
            fn neg(self) -> Self::Output {
                Self(-self.0)
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

        impl AbsDiffEq for $name {
            type Epsilon = f64;

            fn default_epsilon() -> Self::Epsilon {
                f64::EPSILON
            }

            fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
                self.0.abs_diff_eq(&other.0, epsilon)
            }
        }

        /// A constant angle value of 1 unit.
        pub const $unit: $name = $name(1.0);
    };
}

angle_unit!(
    /// An angle measured in radians.
    /// 1 radian is the angle at the center of arc of radius 1 and arc length 1.
    (Radians, RAD, "rad", 1.0)
);

impl Radians {
    pub const PI_2: Radians = Radians(FRAC_PI_2);
    pub const PI: Radians = Radians(PI);
    pub const TWO_PI: Radians = Radians(2. * PI);

    /// Wrap the angle value in `[min, min + 2 * PI]`
    /// Typical values of `min` are `-PI` and `2. * PI`.
    pub fn wrap(self, min: f64) -> f64 {
        let mut na = self.0;
        na = na % (2. * PI);
        while na < min {
            na += 2. * PI;
        }
        while na > min + (2. * PI) {
            na -= 2. * PI;
        }
        na
    }

    /// Clamp the angle value in [min, max].
    pub fn clamp(self, min: f64, max: f64) -> f64 {
        self.0.clamp(min, max)
    }
}

angle_unit!(
    /// An angle measured in degrees (abbr 'deg') with:
    /// 1 deg = PI / 180 rad = 0.017453292519943295 rad.
    (Degrees, DEG, "deg", consts::PI / 180.0)
);

angle_unit!(
    /// An angle measured in gradians (abbr 'grad') with:
    /// 1 grad = PI / 200 rad = 0.015707963267948967 rad.
    (Gradians, GRAD, "grad", consts::PI / 200.0)
);

angle_unit!(
    /// An angle measured in arc seconds (abbr 'sec') with:
    /// 1 sec = PI / 648 000 rad = 0.00000484813681109536 rad.
    (Arcsecs, SEC, "sec", consts::PI / 648_000.0)
);

impl Degrees {
    pub fn dms(deg: f64, min: f64, sec: f64) -> Self {
        let f = deg.signum();
        let deg = deg as f64 + f * min as f64 / 60. + f * sec / 3600.;
        Self(deg)
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::FRAC_PI_2;

    use crate::units::angle::{Angle, Degrees, Radians};

    #[test]
    fn test_to_radians() {
        assert_eq!(Radians(FRAC_PI_2).rad(), FRAC_PI_2);
        assert_eq!(Degrees(90.0).to_radians(), Radians(FRAC_PI_2));
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

    #[test]
    fn test_dms() {
        let deg = Degrees::dms(26., 30., 0.0);
        assert_eq!(deg.0, 26.5);
        let deg = Degrees::dms(-26., 30., 0.0);
        assert_eq!(deg.0, -26.5);
        let deg = Degrees::dms(-0., 30., 0.0);
        assert_eq!(deg.0, -0.5);
    }
}
