use std::f64::consts::{FRAC_PI_2, PI};
use std::fmt::{Display, Formatter};
use std::ops::{Add, AddAssign, Sub, SubAssign};

use crate::quantity::angle::formatters::DMS;
use crate::quantity::angle::units::Rad;
use crate::quantity::angle::Angle;
use approx::AbsDiffEq;
use derive_more::derive::Neg;
use num::Zero;

/// A longitude coordinate in [-pi..pi] radians.
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq, Default, Neg)]
pub struct Lon(f64);

impl Lon {
    pub const MIN: Lon = Lon(-PI);
    pub const MAX: Lon = Lon(PI);

    /// Create a new longitude value with a given raw angle value in radians.
    /// The angle is wrapped into [-pi..pi].
    pub fn new(val: Angle) -> Self {
        Self(val.wrap().rad())
    }

    pub fn zero() -> Self {
        Lon(0.0)
    }

    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    /// Normalize the longitude into (-pi..pi].
    pub fn normalize(self) -> Self {
        if self <= Self::MIN {
            Self::MAX
        } else {
            self
        }
    }

    /// Return the longitude as a raw angle value **in radians**.
    #[inline]
    pub fn rad(self) -> f64 {
        self.0
    }
}

impl Add<Angle> for Lon {
    type Output = Self;

    fn add(self, rhs: Angle) -> Self::Output {
        Self::new(self.0 * Rad + rhs)
    }
}

impl Add<Lon> for Angle {
    type Output = Lon;

    fn add(self, rhs: Lon) -> Self::Output {
        rhs + self
    }
}

impl AddAssign<Angle> for Lon {
    fn add_assign(&mut self, rhs: Angle) {
        self.0 = (self.0 * Rad + rhs).wrap().rad();
    }
}

impl Sub<Angle> for Lon {
    type Output = Self;

    fn sub(self, rhs: Angle) -> Self::Output {
        // Watch out for infinite recursion with self - rhs
        Self::new(self.0 * Rad - rhs)
    }
}

impl AbsDiffEq for Lon {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

impl Sub for Lon {
    type Output = LonInterval;

    fn sub(self, rhs: Self) -> Self::Output {
        LonInterval::new(rhs, self)
    }
}

impl Display for Lon {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", DMS::from_rad(self.0))
    }
}

/// A **closed** longitude interval represented by its lower and upper bounds (inclusive).
/// Note that the lower may be greater that the higher bound, in this case the interval contains
/// the point on the ante-meridian.
/// The lower and higher bound representations may also be inverted in special cases to allow the
/// empty and the full intervals as well:
/// - the empty interval is represented as the inverted interval `[pi..-pi]`
/// - the full interval is represented as the `[-pi..pi]`
/// - other intervals are represented as:
///     - `[lo..hi]` for lo != -pi
///     - `[pi..hi]` for lo == -pi and hi != pi,
///     - `[lo..hi]` otherwise
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct LonInterval {
    lo: f64,
    hi: f64,
}

impl LonInterval {
    /// Return an empty interval.
    pub fn empty() -> Self {
        Self { lo: PI, hi: -PI }
    }

    /// Return a full interval.
    pub fn full() -> Self {
        Self { lo: -PI, hi: PI }
    }

    /// Create a new interval.
    pub fn new(lo: Lon, hi: Lon) -> Self {
        let mut low = lo;
        if lo == Lon::MIN && hi != Lon::MAX {
            low = Lon::MAX;
        }
        Self {
            lo: low.rad(),
            hi: hi.rad(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.lo == PI && self.hi == -PI
    }

    pub fn is_full(&self) -> bool {
        self.lo == -PI && self.hi == PI
    }

    fn is_inverted(&self) -> bool {
        self.lo > self.hi
    }

    /// Return the **positive** length of this interval.
    /// The length is 0 if the interval is empty or a singleton.
    pub fn length(&self) -> Angle {
        let mut length = self.hi - self.lo;
        if length < 0. {
            length += 2. * PI;
        }
        debug_assert!(length >= 0.);
        length * Rad
    }
}

/// A latitude coordinate in [-pi/2..pi/2]  radians.
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq, Default, Neg)]
pub struct Lat(f64);

impl Lat {
    pub const MIN: Lat = Lat(-FRAC_PI_2);
    pub const MAX: Lat = Lat(FRAC_PI_2);

    /// Create a new latitude value with the given raw angle value in radians.
    /// The angle value is clamped into [-pi/2..pi/2].
    pub fn new(val: Angle) -> Self {
        Lat(val.rad().clamp(-FRAC_PI_2, FRAC_PI_2))
    }

    pub fn zero() -> Self {
        Lat(0.0)
    }

    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    /// Return this latitude as a raw angle value in radians.
    #[inline]
    pub fn rad(self) -> f64 {
        self.0
    }
}

impl Add<Angle> for Lat {
    type Output = Self;

    fn add(self, rhs: Angle) -> Self::Output {
        Self::new(self.0 * Rad + rhs)
    }
}

impl Add<Lat> for Angle {
    type Output = Lat;

    fn add(self, rhs: Lat) -> Self::Output {
        rhs + self
    }
}

impl AddAssign<Angle> for Lat {
    fn add_assign(&mut self, rhs: Angle) {
        self.0 = (self.0 * Rad + rhs).rad().clamp(-PI, PI)
    }
}
impl Sub<Angle> for Lat {
    type Output = Self;

    fn sub(self, rhs: Angle) -> Self::Output {
        Self::new(self.0 * Rad - rhs)
    }
}

impl SubAssign<Angle> for Lat {
    fn sub_assign(&mut self, rhs: Angle) {
        self.0 = (self.0 * Rad - rhs).rad().clamp(-PI, PI);
    }
}

impl AbsDiffEq for Lat {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

impl Display for Lat {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", DMS::from_rad(self.0))
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

    use approx::assert_abs_diff_eq;

    use crate::cs::geodetic::{Lat, Lon, LonInterval};
    use crate::quantity::angle::units::{Deg, Rad};

    #[test]
    fn test_lon() {
        // wrapping
        assert_eq!(Lon::new(2.0 * PI * Rad), Lon::new(0.0 * Rad));
        assert_eq!(Lon::new(185.0 * Deg), Lon::new(-175.0 * Deg));
        // normalization
        assert_eq!(Lon::new(-PI * Rad).normalize(), Lon::new(PI * Rad));
        // equality
        assert_eq!(Lon::new(FRAC_PI_2 * Rad), Lon::new(90.0 * Deg));
        assert_ne!(Lon::new(90. * Deg), Lon::new(91.0 * Deg));
        // ops
        assert_eq!(Lon::new(90. * Deg) + 45.0 * Deg, Lon::new(135.0 * Deg));
        assert_eq!(Lon::new(90. * Deg) - 45.0 * Deg, Lon::new(45.0 * Deg));
        // display
        assert_eq!(format!("{}", Lon::new(90. * Deg)), "  90° 00′ 00.00000000″");
    }

    #[test]
    fn test_lat() {
        // clamping
        assert_eq!(Lat::new(91. * Deg), Lat::MAX);
        assert_eq!(Lat::new(-90.01 * Deg), Lat::MIN);
        // equality
        assert_eq!(Lat::new(45. * Deg), Lat::new(FRAC_PI_4 * Rad));
        assert_ne!(Lat::new(0. * Deg), Lat::new(0.001 * Deg));
        // ops
        assert_eq!(Lat::new(45. * Deg + 10. * Deg), Lat::new(55. * Deg));
        assert_eq!(Lat::new(45. * Deg) + 50. * Deg, Lat::MAX);
        assert_eq!(Lat::new(45. * Deg) - 90. * Deg, Lat::new(-45. * Deg));
        assert_eq!(Lat::new(-45. * Deg) - 50. * Deg, Lat::MIN);
        // Display
        assert_eq!(format!("{}", Lat::new(45. * Deg)), "  45° 00′ 00.00000000″");
    }

    #[test]
    fn test_lon_interval() {
        assert!(LonInterval::empty().is_empty());
        assert!(LonInterval::empty().is_inverted());
        assert!(LonInterval::full().is_full());

        // Length
        assert_eq!(LonInterval::empty().length(), 0. * Rad);
        assert_eq!(LonInterval::full().length(), 2. * PI * Rad);
        assert_eq!(
            LonInterval::new(Lon::new(0. * Deg), Lon::new(160. * Deg)).length(),
            160. * Deg
        );
        assert_abs_diff_eq!(
            LonInterval::new(Lon::new(170. * Deg), Lon::new(-170. * Deg)).length(),
            20. * Deg,
            epsilon = 1e-15
        );
        assert_abs_diff_eq!(
            LonInterval::new(Lon::new(170. * Deg), Lon::new(130. * Deg)).length(),
            320. * Deg,
            epsilon = 1e-15
        );
    }
}
