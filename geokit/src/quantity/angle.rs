use crate::quantity::angle::units::DEG;
use num::{FromPrimitive, Num};
use std::f64::consts::PI;

/// Wrap a raw angle value expressed in radians in `[min, min + 2 * PI]`.
/// Typical values of `min` are `-PI` and `2. * PI`.
pub fn wrap(rad: f64, min: f64) -> f64 {
    let mut na = rad;
    na = na % (2. * PI);
    while na < min {
        na += 2. * PI;
    }
    while na > min + (2. * PI) {
        na -= 2. * PI;
    }
    na
}

/// Return the raw angle value in radians of a degree/minute/second value.
///
/// ```
/// use geokit::cs::geodetic::Lon;
/// use geokit::quantity::angle::dms;
/// let lon = Lon::new(dms(2., 20., 14.02500));
/// ```
pub fn dms(d: f64, m: f64, s: f64) -> f64 {
    let f = d.signum();
    let deg = d + f * m / 60. + f * s / 3600.;
    deg * DEG
}

/// Defines conversion rate from angle unit to radians.
pub mod units {
    use std::f64::consts::PI;

    /// 1 radian. The reference angle unit.
    pub const RAD: f64 = 1.0;

    /// 1 degree: `1 deg = PI / 180 rad`
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use geokit::quantity::angle::units::DEG;
    /// assert_eq!(1.0 * DEG, PI / 180.);
    /// ```
    pub const DEG: f64 = PI / 180.0;

    /// 1 gradian: `1 grad = PI / 200. rad`
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use geokit::quantity::angle::units::GRAD;
    /// assert_eq!(1.0 * GRAD, PI / 200.);
    /// ```
    pub const GRAD: f64 = PI / 200.0;

    /// 1 arcsec: `1 sec = PI / 648_000. rad`
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use geokit::quantity::angle::units::SEC;
    /// assert_eq!(1. * SEC, PI / 648_000.);
    /// ```
    pub const SEC: f64 = PI / 648_000.0;
}
