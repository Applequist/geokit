use crate::quantity::angle::units::DEG;
use std::f64::consts::PI;
use std::fmt::{Display, Formatter};

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
    pub fn from_rad(rad: f64) -> DMS {
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

    /// Return the sign of the angle represented by this value.
    pub fn sign(&self) -> f64 {
        self.0.signum()
    }

    /// Return the absolute value of the degree part.
    pub fn deg(&self) -> f64 {
        self.0.abs()
    }

    /// Return the absolute value of the minute part.
    pub fn min(&self) -> f64 {
        self.1
    }

    /// Return the absolute value of the second part.
    pub fn sec(&self) -> f64 {
        self.2
    }
}

impl Display for DMS {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:4}° {:02}′ {:011.8}″", self.0, self.1, self.2)
    }
}

/// Return the raw angle value in radians of a degree/minute/second value.
/// The returned angle has the same size as the degrees parameter `d`.
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
/// use geokit::cs::geodetic::Lon;
/// use geokit::quantity::angle::dms;
/// let lon = Lon::new(dms(2., 20., 14.02500));
/// assert!(dms(-12., 45., 59.1234) < 0.0);
/// ```
pub fn dms(d: f64, m: f64, s: f64) -> f64 {
    debug_assert!(m >= 0. && s >= 0., "minutes and seconds must be >= 0.0");
    let f = d.signum();
    let deg = d.abs() + m / 60. + s / 3600.;
    f * deg * DEG
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

#[cfg(test)]
mod tests {
    use crate::quantity::angle::{dms, DMS};

    #[test]
    fn dms_display() {
        assert_eq!(
            format!("{}", DMS(-33., 6., 22.01545)),
            " -33° 06′ 22.01545000″"
        );
        assert_eq!(format!("{}", DMS(45., 0., 0.)), "  45° 00′ 00.00000000″");
        // assert_eq!(format!("{}", DMS::fromom_rad(FRAC_PI_4)), "  45° 00′ 00.00000000″");
        assert_eq!(
            format!("{}", DMS(-179., 59., 59.1234)),
            "-179° 59′ 59.12340000″"
        );
        assert_eq!(
            format!("{}", DMS::from_rad(dms(-33., 26., 0.))),
            " -33° 26′ 00.00000000″"
        );
    }

    #[test]
    fn dms_rounding() {
        assert_eq!(
            format!("{}", DMS::from_rad(dms(37., 19., 54.95367))),
            "  37° 19′ 54.95367000″"
        );
    }
}
