use crate::{
    math::fp::{remainder, Float, PI, PI_2, TAU},
    units::angle::{AngleUnit, DEG},
};
use approx::AbsDiffEq;
use derive_more::derive::{Add, AddAssign, Display, Neg, Sub, SubAssign};
use std::{
    fmt::LowerExp,
    ops::{Div, DivAssign, Mul, MulAssign},
};

/// [Angle] represents a 1-dimensional angle value.
/// The internal representation is a [Float] value in radians.
///
/// # Creation
///
/// There are several ways to create an [Angle] value:
/// - by using [Angle::new], passing a quantity and an [AngleUnit].
///   The quantity expressed in the given unit is then converted to radians.
/// ```
/// use geokit::quantities::angle::Angle;
/// use geokit::units::angle::GRAD;
/// let a = Angle::new(100., GRAD);
/// ```
/// - by multiplying a [Float] quantity by an [AngleUnit] value:
///   The quantity expressed in the given unit is then converted to radians.
/// ```
/// use geokit::quantities::angle::Angle;
/// use geokit::units::angle::GRAD;
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
/// - [Wrapping][Self::wrapped()] in [-pi, pi]
/// - [Normalization][Self::normalized()] in (-pi, pi]
/// - [Conversion][Self::val()] to other [AngleUnit]
/// - trigonometric operations like cos, sin, tan...
/// - formatting using [DMS][Self::to_dms].
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

    /// 1e-12 rad.
    /// Less than 7e-6 m at the equator
    #[inline]
    pub const fn tiny() -> Angle {
        Angle(1e-12)
    }

    /// 1e-9 rad.
    /// Less than 7e-3 m at the equator
    #[inline]
    pub const fn small() -> Angle {
        Angle(1e-9)
    }

    /// Default to `Self::tiny`.
    pub const fn default_epsilon() -> Angle {
        Self::tiny()
    }

    /// Create an angle value whose `qty` is given in `unit`.
    #[inline]
    pub const fn new(qty: Float, unit: AngleUnit) -> Self {
        Angle(qty * unit.rad_per_unit())
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
    /// use geokit::units::angle::DEG;
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

    /// [Normalize][Self::normalized] in place.
    pub fn normalize(&mut self) {
        let mut a = remainder(self.0, TAU);
        if a <= -PI {
            a = PI;
        }
        self.0 = a;
    }

    /// Return this angle value in the given unit.
    ///
    /// # Example
    /// ```
    /// use geokit::units::angle::DEG;
    /// use geokit::quantities::angle::Angle;
    /// let a = Angle::PI_2;
    /// let a_deg = a.val(DEG);
    /// assert_eq!(a, a_deg * DEG);
    /// ```
    pub fn val(self, unit: AngleUnit) -> Float {
        self.0 * unit.1 / unit.0
    }

    /// Return the absolute value of this angle.
    pub fn abs(self) -> Angle {
        Angle(self.0.abs())
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

    /// Returns the smallest *signed* angle from this angle to the `other` angle.
    pub fn diff_to(self, other: Self) -> Angle {
        (other - self).wrapped()
    }

    pub fn deg(self) -> Deg {
        Deg(self.val(DEG))
    }

    /// Convert this angle into a [Dms] value for formatting.
    pub fn dms(self) -> Dms {
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

    pub fn signum(&self) -> f64 {
        self.0.signum()
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
    type Epsilon = Angle;

    fn default_epsilon() -> Self::Epsilon {
        Angle::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        (*self - *other).abs() <= epsilon
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

impl LowerExp for Angle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        LowerExp::fmt(&self.0, f)?;
        f.write_str(" rad")
    }
}

/// Use to format angle using their values in degrees.
#[derive(Copy, Clone, PartialEq, Debug, Display)]
#[display("{} deg", _0)]
pub struct Deg(Float);

impl LowerExp for Deg {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        LowerExp::fmt(&self.0, f)?;
        f.write_str(" deg")
    }
}

/// An angle value expressed in degrees, minutes and seconds.
///
/// Use to display angle in DMS format:
/// ```
/// use geokit::quantities::angle::Angle;
/// use geokit::units::angle::DEG;
/// assert_eq!(format!("{}", Angle::new(2.33432, DEG).dms()), "   2° 20′ 03.55200000″");
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

#[cfg(test)]
mod tests {

    use super::Dms;
    use crate::{
        math::fp::{Float, PI, PI_4},
        quantities::angle::Angle,
        units::angle::{DEG, GRAD, RAD, SEC},
    };
    use approx::assert_abs_diff_eq;

    #[test]
    fn angle_unit_equivalence() {
        assert_eq!(1. * DEG, (PI / 180.0) * RAD);
        assert_eq!(1. * GRAD, (PI / 200.0) * RAD);
        assert_eq!(1. * SEC, (PI / 648_000.0) * RAD);
    }

    #[test]
    fn angle_wrapping() {
        let input_expected = [
            ("< -180", -210. * DEG, 150. * DEG),
            // FIX: Robustness issue
            // -180. - n * eps < -pi for n >= 65 on my test
            //(
            //    "-180 - eps",
            //    (-180. - Float::EPSILON) * DEG,
            //    (180. - Float::EPSILON) * DEG,
            //),
            (
                "-180 + eps",
                (-180. + Float::EPSILON) * DEG,
                (-180. + Float::EPSILON) * DEG,
            ),
            ("-180", -180. * DEG, -180. * DEG),
            ("-90", -90. * DEG, -90. * DEG),
            ("0", 0. * DEG, 0. * DEG),
            ("90", 90. * DEG, 90. * DEG),
            (
                "180 - eps",
                (180. - Float::EPSILON) * DEG,
                (180. - Float::EPSILON) * DEG,
            ),
            ("180", 180. * DEG, 180. * DEG),
            // FIX: Robustness issue
            //(
            //    "180 + eps",
            //    (180. + Float::EPSILON) * DEG,
            //    (-180. + Float::EPSILON) * DEG,
            //),
        ];

        for (t, a, e) in input_expected {
            let wrapped = a.wrapped();
            assert!(
                (wrapped - e).rad().abs() < 1e-15,
                "case {}: expected {}, got {}",
                t,
                e,
                wrapped
            );
        }
    }

    #[test]
    fn angle_normalization() {
        assert_eq!((180. * DEG).normalized(), 180. * DEG);
        assert_eq!((-180. * DEG).normalized(), 180. * DEG);
        assert_eq!(Angle::M_PI.normalized(), 180. * DEG);
    }

    #[test]
    fn angle_diff_to() {
        let a1 = 10. * DEG;
        let a2 = 80. * DEG;
        let a3 = 110. * DEG;
        let a4 = 300. * DEG;
        let a5 = -10. * DEG;
        let a6 = -90. * DEG;
        let a7 = 90. * DEG;

        assert_abs_diff_eq!(a1.diff_to(a2), 70. * DEG, epsilon = Angle(1e-15));
        assert_abs_diff_eq!(a1.diff_to(a3), 100. * DEG, epsilon = Angle(1e-15));
        assert_abs_diff_eq!(a1.diff_to(a4), -70. * DEG, epsilon = Angle(1e-15));
        assert_abs_diff_eq!(a1.diff_to(a5), -20. * DEG, epsilon = Angle(1e-15));
        assert_abs_diff_eq!(a7.diff_to(a6), -180. * DEG, epsilon = Angle(1e-15));
        assert_abs_diff_eq!(a6.diff_to(a7), 180. * DEG, epsilon = Angle(1e-15));
    }

    #[test]
    fn angle_dms_display() {
        assert_eq!(
            format!("{}", Dms(-33., 6., 22.01545)),
            " -33° 06′ 22.01545000″"
        );
        assert_eq!(format!("{}", Dms(45., 0., 0.)), "  45° 00′ 00.00000000″");
        assert_eq!(format!("{}", (PI_4 * RAD).dms()), "  45° 00′ 00.00000000″");
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
}
