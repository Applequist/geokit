use approx::AbsDiffEq;
use derive_more::derive::{Add, AddAssign, Neg, Sub, SubAssign};

use crate::math::utils::wrap;
use std::{
    f64::consts::PI,
    ops::{Div, DivAssign, Mul, MulAssign},
};

/// An angle value.
/// Angles are created by multiplying a real valued quantity by a unit, eg:
/// ```
/// let a: Angle = 1.0 * Deg;
/// ```
/// Or using the [dms] function.
#[derive(Debug, Copy, Clone, PartialEq, Default, Add, AddAssign, Sub, SubAssign, Neg)]
pub struct Angle(f64);

impl Angle {
    pub const PI: Angle = Angle(PI);
    pub const M_PI: Angle = Angle(-PI);

    /// Return the angle value in radians.
    pub fn rad(&self) -> f64 {
        self.0
    }

    /// Wrap this angle into [-PI, PI] radians.
    pub fn wrap(self) -> Angle {
        Angle(wrap(self.0, PI))
    }
}

impl Mul<Angle> for f64 {
    type Output = Angle;

    fn mul(self, rhs: Angle) -> Self::Output {
        Angle(self * rhs.0)
    }
}

impl Mul<f64> for Angle {
    type Output = Angle;

    fn mul(self, rhs: f64) -> Self::Output {
        Angle(self.0 * rhs)
    }
}

impl MulAssign<f64> for Angle {
    fn mul_assign(&mut self, rhs: f64) {
        self.0 *= rhs;
    }
}

impl Div<f64> for Angle {
    type Output = Angle;

    fn div(self, rhs: f64) -> Self::Output {
        Angle(self.0 / rhs)
    }
}

impl DivAssign<f64> for Angle {
    fn div_assign(&mut self, rhs: f64) {
        self.0 /= rhs;
    }
}

impl AbsDiffEq for Angle {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

/// Create an 'Angle' value from a degree/minute/second value.
/// The returned angle has the same sign as the degrees parameter `d`.
///
/// # Parameters
///
/// - `d`: the number of degrees. Determine the sign of the returned angle.
/// - `m`: the number of minutes. Must be >= 0.
/// - `s`: the number of seconds with fractional part. Must be >= 0.
///
/// # Examples
///
/// ```
/// use geokit::quantity::angle::dms;
/// let a: Angle = dms(2., 20., 14.02500);
/// assert!(dms(-12., 45., 59.1234).rad() < 0.0);
/// ```
pub fn dms(d: f64, m: f64, s: f64) -> Angle {
    debug_assert!(m >= 0. && s >= 0., "minutes and seconds must be >= 0.0");
    let f = d.signum();
    let deg = d.abs() + m / 60. + s / 3600.;
    f * deg * units::DEG
}

pub mod units {
    use std::{f64::consts::PI, ops::Mul};

    use super::Angle;

    /// A to-radians angle converter.
    #[derive(Debug, Copy, Clone, PartialEq)]
    pub struct AngleUnit(f64, f64);

    impl AngleUnit {
        #[inline]
        pub const fn rad_per_unit(&self) -> f64 {
            self.0 / self.1
        }
    }

    impl Mul<AngleUnit> for f64 {
        type Output = Angle;

        fn mul(self, rhs: AngleUnit) -> Self::Output {
            Angle(self * rhs.0 / rhs.1)
        }
    }

    impl Mul<AngleUnit> for i32 {
        type Output = Angle;

        fn mul(self, rhs: AngleUnit) -> Self::Output {
            Angle((self as f64) * rhs.0 / rhs.1)
        }
    }

    pub const RAD: AngleUnit = AngleUnit(1.0, 1.0);
    pub const DEG: AngleUnit = AngleUnit(PI, 180.0);
    pub const GRAD: AngleUnit = AngleUnit(PI, 200.0);
    pub const SEC: AngleUnit = AngleUnit(PI, 648_000.0);
}

pub mod formatters {
    use std::fmt::{Display, Formatter};

    use super::Angle;

    /// Parse a raw angle value given in radians into degrees, minutes and seconds.
    ///
    /// Use to display angle in DMS format:
    /// ```
    /// use geokit::quantity::angle::{dms, DMS};
    /// assert_eq!(format!("{}", DMS::from_rad(dms(2., 20., 14.02500))), "   2° 20′ 14.02500000″");
    /// ```
    #[derive(Copy, Clone, Debug)]
    pub struct DMS(f64, f64, f64);

    impl DMS {
        pub fn from_rad(rad: f64) -> Self {
            let sgn = rad.signum();
            let mut deg = rad.abs().to_degrees();
            let d = deg.floor();
            deg = 60. * deg.fract();
            let m;
            let s;
            if (deg - deg.round()).abs() < 1e-8 {
                m = deg.round();
                s = 0.0;
            } else {
                m = deg.floor();
                s = 60. * deg.fract();
            }
            DMS(sgn * d, m, s)
        }

        pub fn deg(&self) -> f64 {
            self.0
        }

        pub fn min(&self) -> f64 {
            self.1
        }

        pub fn sec(&self) -> f64 {
            self.2
        }
    }

    impl From<Angle> for DMS {
        fn from(value: Angle) -> Self {
            DMS::from_rad(value.rad())
        }
    }

    impl Display for DMS {
        fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
            write!(f, "{:4}° {:02}′ {:011.8}″", self.0, self.1, self.2)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::quantity::angle::{
        dms,
        formatters::DMS,
        units::{DEG, GRAD, RAD, SEC},
    };

    use std::f64::consts::PI;

    #[test]
    fn units() {
        assert_eq!(1. * DEG, (PI / 180.0) * RAD);
        assert_eq!(1. * GRAD, (PI / 200.0) * RAD);
        assert_eq!(1. * SEC, (PI / 648_000.0) * RAD);
    }

    #[test]
    fn formatters() {
        assert_eq!(
            format!("{}", DMS::from(dms(-33., 6., 22.01545))),
            " -33° 06′ 22.01545000″"
        );
        assert_eq!(
            format!("{}", DMS::from(dms(45., 0., 0.))),
            "  45° 00′ 00.00000000″"
        );
        // assert_eq!(format!("{}", DMS::fromom_rad(FRAC_PI_4)), "  45° 00′ 00.00000000″");
        assert_eq!(
            format!("{}", DMS::from(dms(-179., 59., 59.1234))),
            "-179° 59′ 59.12340000″"
        );
        assert_eq!(
            format!("{}", DMS::from(dms(-33., 26., 0.))),
            " -33° 26′ 00.00000000″"
        );
        // Check rounding
        assert_eq!(
            format!("{}", DMS::from(dms(37., 19., 54.95367))),
            "  37° 19′ 54.95367000″"
        );
    }
}
