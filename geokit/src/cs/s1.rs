//! The S1 *abstract* coordinate space is used to represent points from the unit circle.
//! Each point on the circle is represented by the angle at the center of the circle between
//! an origin and the point, in the range [-pi..pi].
//!
//! This module provides the following value types:
//! - [Angle] to represent a single S1 point coordinate
//! - and [interval][Interval] to represent a segment on S1.
//!
//! These types are used to define other value types in more concrete CS like 2D- and 3D-
//! geodetic CS.

use crate::{
    math::{utils::remainder, Float, PI, PI_2, TAU},
    units::angle::{AngleUnit, DEG},
};
use approx::AbsDiffEq;
use derive_more::derive::{Add, AddAssign, Display, Neg, Sub, SubAssign};
use std::ops::{Div, DivAssign, Mul, MulAssign};

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
/// - by multiplying a [Float] quantity by an [AngleUnit] value:
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
pub struct Angle(Float);

impl Angle {
    pub const ZERO: Angle = Angle(0.0);
    pub const PI_2: Angle = Angle(PI_2);
    pub const PI: Angle = Angle(PI);
    pub const TWO_PI: Angle = Angle(TAU);
    pub const M_PI_2: Angle = Angle(-PI_2);
    pub const M_PI: Angle = Angle(-PI);

    /// Create an angle value whose `qty` is given in `unit`.
    #[inline]
    pub const fn new(qty: Float, unit: AngleUnit) -> Self {
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
    pub fn dms(d: Float, m: Float, s: Float) -> Angle {
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
        let mut a = remainder(self.0, TAU);
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
        let mut a = remainder(self.0, TAU);
        if a <= -PI {
            a = PI;
        }
        Self(a)
    }

    pub fn normalize(&mut self) {
        let mut a = remainder(self.0, TAU);
        if a <= -PI {
            a = PI;
        }
        self.0 = a;
    }

    /// Return the absolute value of this angle.
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
    pub fn val(self, unit: AngleUnit) -> Float {
        self.0 * unit.1 / unit.0
    }

    /// Return this angle value in radians.
    #[inline]
    pub fn rad(self) -> Float {
        self.0
    }

    #[inline]
    pub fn sin(self) -> Float {
        self.0.sin()
    }

    #[inline]
    pub fn cos(self) -> Float {
        self.0.cos()
    }

    #[inline]
    pub fn sin_cos(self) -> (Float, Float) {
        self.0.sin_cos()
    }

    #[inline]
    pub fn tan(self) -> Float {
        self.0.tan()
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

impl Mul<Angle> for Float {
    type Output = Angle;

    fn mul(self, rhs: Angle) -> Self::Output {
        Angle(self * rhs.0)
    }
}

impl Mul<Float> for Angle {
    type Output = Angle;

    fn mul(self, rhs: Float) -> Self::Output {
        Angle(self.0 * rhs)
    }
}

impl MulAssign<Float> for Angle {
    fn mul_assign(&mut self, rhs: Float) {
        self.0 *= rhs;
    }
}

impl Div<Float> for Angle {
    type Output = Angle;

    fn div(self, rhs: Float) -> Self::Output {
        Angle(self.0 / rhs)
    }
}

impl DivAssign<Float> for Angle {
    fn div_assign(&mut self, rhs: Float) {
        self.0 /= rhs;
    }
}

impl AbsDiffEq for Angle {
    type Epsilon = Float;

    fn default_epsilon() -> Self::Epsilon {
        Float::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

impl Mul<AngleUnit> for Float {
    type Output = Angle;

    fn mul(self, rhs: AngleUnit) -> Self::Output {
        Angle::new(self, rhs)
    }
}

impl Mul<Float> for AngleUnit {
    type Output = Angle;

    fn mul(self, rhs: Float) -> Self::Output {
        rhs * self
    }
}

impl Mul<AngleUnit> for i32 {
    type Output = Angle;

    fn mul(self, rhs: AngleUnit) -> Self::Output {
        Angle::new(self as Float, rhs)
    }
}

impl Mul<i32> for AngleUnit {
    type Output = Angle;

    fn mul(self, rhs: i32) -> Self::Output {
        (rhs as Float) * self
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
pub struct Dms(Float, Float, Float);

impl Dms {
    pub fn deg(&self) -> Float {
        self.0
    }

    pub fn min(&self) -> Float {
        self.1
    }

    pub fn sec(&self) -> Float {
        self.2
    }
}

/// An [Interval] is a **closed** interval on the unit circle
/// represented by its lower and upper bounds (inclusive).
///
/// Note that the lower may be greater that the higher bound,
/// in this case the interval contains the point (-1, 0) of the unit circle.
///
/// The lower and higher bound representations may *also* be inverted
/// in special cases to allow the empty and the full intervals as well:
/// - the empty interval is represented as the inverted interval `[pi..-pi]`
/// - the full interval is represented as the `[-pi..pi]`
/// - other intervals are represented as:
///     - `[lo..hi]` for lo != -pi
///     - `[-pi..hi]` for lo == -pi and hi != pi,
///     - `[lo..hi]` otherwise
#[derive(Debug, Copy, Clone, PartialEq, Display)]
#[display("lo = {}, hi = {}", lo, hi)]
pub struct Interval {
    lo: Angle,
    hi: Angle,
}

impl Interval {
    /// Create an empty interval.
    pub fn empty() -> Self {
        Self {
            lo: Angle::PI,
            hi: Angle::M_PI,
        }
    }

    /// Create an interval only containing the single value `p`.
    pub fn singleton(p: Angle) -> Self {
        let p_n = p.normalized();
        Self { lo: p_n, hi: p_n }
    }

    /// Create a full interval [-pi..pi].
    pub fn full() -> Self {
        Self {
            lo: Angle::M_PI,
            hi: Angle::PI,
        }
    }

    /// Create a new interval.
    ///
    /// The actual interval
    /// - the empty interval is represented as the inverted interval `[pi..-pi]`
    /// - the full interval is represented as the `[-pi..pi]`
    ///     - `[lo..hi]` for lo != -pi
    ///     - `[-pi..hi]` for lo == -pi and hi != pi,
    ///     - `[lo..hi]` otherwise
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

    pub(crate) fn is_inverted(&self) -> bool {
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

    use approx::assert_abs_diff_eq;

    use super::Dms;
    use crate::{
        cs::s1::{Angle, Interval},
        math::{PI, PI_4},
        units::angle::{DEG, GRAD, RAD, SEC},
    };

    #[test]
    fn angle_unit_equivalence() {
        assert_eq!(1. * DEG, (PI / 180.0) * RAD);
        assert_eq!(1. * GRAD, (PI / 200.0) * RAD);
        assert_eq!(1. * SEC, (PI / 648_000.0) * RAD);
    }

    #[test]
    fn angle_normalization() {
        assert_eq!((-180. * DEG).normalized(), 180. * DEG);
        assert_eq!(Angle::M_PI.normalized(), 180. * DEG);
    }

    #[test]
    fn angle_dms_display() {
        assert_eq!(
            format!("{}", Dms(-33., 6., 22.01545)),
            " -33° 06′ 22.01545000″"
        );
        assert_eq!(format!("{}", Dms(45., 0., 0.)), "  45° 00′ 00.00000000″");
        assert_eq!(
            format!("{}", (PI_4 * RAD).to_dms()),
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

    fn check_interval(
        tag: &str,
        ab: Interval,
        empty: bool,
        full: bool,
        inverted: bool,
        length: Angle,
    ) {
        dbg!(tag);
        assert_eq!(ab.is_empty(), empty, "Expected {} empty", ab);
        assert_eq!(ab.is_full(), full, "Expected {} full", ab);
        assert_eq!(ab.is_inverted(), inverted, "Expected {} inverted", ab);
        assert_abs_diff_eq!(ab.length(), length, epsilon = 1e-15);
    }

    #[test]
    fn empty_interval() {
        check_interval("empty", Interval::empty(), true, false, true, Angle::ZERO);
    }

    #[test]
    fn singleton_interval() {
        check_interval(
            "singleton",
            Interval::singleton(10. * DEG),
            false,
            false,
            false,
            Angle::ZERO,
        );
    }

    #[test]
    fn full_interval() {
        check_interval("full", Interval::full(), false, true, false, Angle::TWO_PI);
    }

    #[test]
    fn interval() {
        // Non-inverted intervals
        let ab = Interval::new(5. * DEG, 15 * DEG);
        check_interval("[5..10]", ab, false, false, false, 10. * DEG);

        let cd = Interval::new(-110. * DEG, 90. * DEG);
        check_interval("[-110..90]", cd, false, false, false, 200. * DEG);

        let ef = Interval::new(-170. * DEG, -30. * DEG);
        check_interval("[-170..-30]", ef, false, false, false, 140. * DEG);

        // Inverted intervals
        let ab_inv = Interval::new(15. * DEG, 5. * DEG);
        check_interval("[15..5]", ab_inv, false, false, true, 350. * DEG);

        let cd_inv = Interval::new(90. * DEG, -110. * DEG);
        check_interval("[90..-110]", cd_inv, false, false, true, 160. * DEG);

        let ef_inv = Interval::new(-30. * DEG, -170. * DEG);
        check_interval("[-30..-170]", ef_inv, false, false, true, 220. * DEG);
    }
}
