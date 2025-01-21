use crate::{
    math::utils::remainder,
    units::angle::{AngleUnit, DEG},
};
use approx::AbsDiffEq;
use derive_more::derive::{Add, AddAssign, Display, Neg, Sub, SubAssign};
use std::{
    f64::consts::{FRAC_PI_2, PI},
    ops::{Div, DivAssign, Mul, MulAssign},
};

/// [Angle] is a generic 1-dimensional angle value type.
///
/// # Creation
///
/// There are several ways to create an [Angle] value:
/// - by using [Angle::new], passing a quantity and an [AngleUnit]:
/// ```
/// use crate::quantities::angle::{Angle, units::GRAD};
/// let a = Angle::new(100., GRAD);
/// ```
/// - by using one of the convenience constructors:
/// ```
/// use crate::quantities::angle::{Angle};
/// let dms = Angle::dms(-110., 45., 53.234);
/// ```
/// - by multiplying a [f64] quantity by an [AngleUnit] value:
/// ```
/// use crate::quantities::angle::{Angle, units::GRAD};
/// let a: Angle = 100. * GRAD;
/// ```
///
/// # Operations
///
/// [Angle] supports the following operations:
/// - Addition, subtraction,
/// - negation
/// - multiplication by a scalar (left and right)
/// - division by a scalar (right)
///
/// To make sure that all points on the unit circle can be represented by a unique
/// angle value, [Angle] also has a ***normalization*** operation that wrap its value
/// in the (-pi, pi] raddians range.
#[derive(
    Debug, Copy, Clone, PartialEq, PartialOrd, Add, AddAssign, Sub, SubAssign, Neg, Display,
)]
#[display("{} rad", _0)]
pub struct Angle(f64);

impl Angle {
    pub const ZERO: Angle = Angle(0.0);
    pub const PI_2: Angle = Angle(FRAC_PI_2);
    pub const PI: Angle = Angle(PI);
    pub const TWO_PI: Angle = Angle(2. * PI);
    pub const M_PI_2: Angle = Angle(-FRAC_PI_2);
    pub const M_PI: Angle = Angle(-PI);

    /// Create an angle value whose `qty` is given in `unit`.
    #[inline]
    pub const fn new(qty: f64, unit: AngleUnit) -> Self {
        Angle(qty * unit.rad_per_unit())
    }

    /// Create an angle value from a degree/minute/second values.
    ///
    /// # Parameters
    ///
    /// - `d`: the number of degrees. Determine the sign of the returned angle.
    /// - `m`: the number of minutes. Must be >= 0. and <= 59.
    /// - `s`: the number of seconds with fractional part. Must be >= 0 and < 60.
    ///
    pub fn dms(d: f64, m: f64, s: f64) -> Angle {
        debug_assert!(m >= 0. && m <= 59.0, "minutes must be in [0..59]");
        debug_assert!(s >= 0. && s < 60., "seconds must be in [0..60)");
        let f = d.signum();
        let deg = d.abs() + m / 60. + s / 3600.;
        f * deg * DEG
    }

    /// Return the angle value clamped in [-pi/2..pi/2]
    pub(crate) fn clamped(self, min: Angle, max: Angle) -> Self {
        debug_assert!(
            min < max,
            "Expected min < max. Got min = {} and max = {}",
            min,
            max
        );
        Self(self.0.clamp(min.rad(), max.rad()))
    }

    /// Return the angle value wrapped in [-pi, pi].
    pub(crate) fn wrapped(self) -> Self {
        let mut a = remainder(self.0, 2. * PI);
        if a < -PI {
            a = PI;
        }
        Self(a)
    }

    /// Return this angle value wrapped into (-PI, PI] radians.
    ///
    /// # Example:
    /// ```
    /// let a = 185. * DEG;
    /// let normalized  = a.normalized();
    /// assert_eq!(normalized, -175. * DEG);
    /// ```
    pub fn normalized(self) -> Self {
        let mut a = remainder(self.0, 2. * PI);
        if a <= -PI {
            a = PI;
        }
        Self(a)
    }

    pub fn normalize(&mut self) {
        let mut a = remainder(self.0, 2. * PI);
        if a <= -PI {
            a = PI;
        }
        self.0 = a;
    }

    pub fn abs(self) -> Angle {
        Angle(self.0.abs())
    }

    /// Return this angle value in the given unit.
    ///
    /// # Example
    /// ```
    /// let a = Angle::PI_2;
    /// let a_deg = a.val(DEG);
    /// assert_eq!(a, a_deg * DEG);
    /// ```
    pub fn val(self, unit: AngleUnit) -> f64 {
        self.0 * unit.1 / unit.0
    }

    /// Return this angle value in radians.
    #[inline]
    pub fn rad(self) -> f64 {
        self.0
    }

    /// Convert this angle into a [Dms] value for formatting.
    pub fn to_dms(self) -> Dms {
        let rad = self.0;
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
        Dms(sgn * d, m, s)
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

impl Mul<AngleUnit> for f64 {
    type Output = Angle;

    fn mul(self, rhs: AngleUnit) -> Self::Output {
        Angle::new(self, rhs)
    }
}

impl Mul<f64> for AngleUnit {
    type Output = Angle;

    fn mul(self, rhs: f64) -> Self::Output {
        rhs * self
    }
}

impl Mul<AngleUnit> for i32 {
    type Output = Angle;

    fn mul(self, rhs: AngleUnit) -> Self::Output {
        Angle::new(self as f64, rhs)
    }
}

impl Mul<i32> for AngleUnit {
    type Output = Angle;

    fn mul(self, rhs: i32) -> Self::Output {
        (rhs as f64) * self
    }
}

/// An angle value expressed in degrees, minutes and seconds.
///
/// Use to display angle in DMS format:
/// ```
/// use geokit::quantity::angle::{dms, DMS};
/// assert_eq!(format!("{}", DMS::from_rad(dms(2., 20., 14.02500))), "   2° 20′ 14.02500000″");
/// ```
#[derive(Copy, Clone, Debug, Display)]
#[display("{:4}° {:02}′ {:011.8}″", _0, _1, _2)]
pub struct Dms(f64, f64, f64);

impl Dms {
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

/// An [Interval] is a **closed** interval on the unit circle
/// represented by its lower and upper bounds (inclusive).
///
/// Note that the lower may be greater that the higher bound,
/// in this case the interval contains the point on the ante-meridian.
///
/// The lower and higher bound representations may also be inverted
/// in special cases to allow the empty and the full intervals as well:
/// - the empty interval is represented as the inverted interval `[pi..-pi]`
/// - the full interval is represented as the `[-pi..pi]`
/// - other intervals are represented as:
///     - `[lo..hi]` for lo != -pi
///     - `[-pi..hi]` for lo == -pi and hi != pi,
///     - `[lo..hi]` otherwise
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Interval {
    lo: Angle,
    hi: Angle,
}

impl Interval {
    /// Return an empty interval.
    pub fn empty() -> Self {
        Self {
            lo: Angle::PI,
            hi: Angle::M_PI,
        }
    }

    /// Return a full interval.
    pub fn full() -> Self {
        Self {
            lo: Angle::M_PI,
            hi: Angle::PI,
        }
    }

    /// Create a new interval.
    pub fn new(lo: Angle, hi: Angle) -> Self {
        let mut low = lo.wrapped();
        let hig = hi.wrapped();
        if low == Angle::M_PI && hig != Angle::PI {
            low = Angle::PI;
        }
        Self { lo: low, hi: hig }
    }

    pub fn is_empty(&self) -> bool {
        self.lo == Angle::PI && self.hi == Angle::M_PI
    }

    pub fn is_full(&self) -> bool {
        self.lo == Angle::M_PI && self.hi == Angle::PI
    }

    fn is_inverted(&self) -> bool {
        self.lo > self.hi
    }

    /// Return the **positive** length of this interval.
    /// The length is 0 if the interval is empty or a singleton.
    pub fn length(&self) -> Angle {
        let mut length = self.hi - self.lo;
        if length < Angle::ZERO {
            length += 2. * Angle::PI;
        }
        debug_assert!(length >= Angle::ZERO);
        length
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_PI_4, PI};

    use approx::assert_abs_diff_eq;

    use super::Dms;
    use crate::{
        cs::s1::{Angle, Interval},
        units::angle::{DEG, GRAD, RAD, SEC},
    };

    #[test]
    fn conversions() {
        assert_eq!(1. * DEG, (PI / 180.0) * RAD);
        assert_eq!(1. * GRAD, (PI / 200.0) * RAD);
        assert_eq!(1. * SEC, (PI / 648_000.0) * RAD);
    }

    #[test]
    fn normalization() {
        assert_eq!((-180. * DEG).normalized(), 180. * DEG);
        assert_eq!(Angle::M_PI.normalized(), 180. * DEG);
    }

    #[test]
    fn dms_display() {
        assert_eq!(
            format!("{}", Dms(-33., 6., 22.01545)),
            " -33° 06′ 22.01545000″"
        );
        assert_eq!(format!("{}", Dms(45., 0., 0.)), "  45° 00′ 00.00000000″");
        assert_eq!(
            format!("{}", (FRAC_PI_4 * RAD).to_dms()),
            "  45° 00′ 00.00000000″"
        );
        assert_eq!(
            format!("{}", Dms(-179., 59., 59.1234)),
            "-179° 59′ 59.12340000″"
        );
        assert_eq!(format!("{}", Dms(-33., 26., 0.)), " -33° 26′ 00.00000000″");
        // Check rounding
        assert_eq!(
            format!("{}", Dms(37., 19., 54.95367)),
            "  37° 19′ 54.95367000″"
        );
    }

    #[test]
    fn test_interval() {
        assert!(Interval::empty().is_empty());
        assert!(Interval::empty().is_inverted());
        assert!(Interval::full().is_full());

        // Length
        assert_eq!(Interval::empty().length(), Angle::ZERO);
        assert_eq!(Interval::new(20. * DEG, 20. * DEG).length(), Angle::ZERO);
        assert_eq!(Interval::full().length(), 2. * PI * RAD);
        assert_eq!(Interval::new(0. * DEG, 160. * DEG).length(), 160. * DEG);
        assert_abs_diff_eq!(
            Interval::new(170. * DEG, -170. * DEG).length(),
            20. * DEG,
            epsilon = 1e-15
        );
        assert_abs_diff_eq!(
            Interval::new(170. * DEG, 130. * DEG).length(),
            320. * DEG,
            epsilon = 1e-15
        );
    }
}
