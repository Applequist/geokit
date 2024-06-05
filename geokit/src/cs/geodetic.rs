use std::f64::consts::{FRAC_PI_2, PI};
use std::fmt::{Display, Formatter};
use std::ops::{Add, Neg, Sub};

use crate::quantity::angle::{wrap, DMS};
use approx::AbsDiffEq;

/// A longitude coordinate in [-pi..pi] radians.
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq, Default)]
pub struct Lon(f64);

impl Lon {
    pub const MIN: Lon = Lon(-PI);
    pub const MAX: Lon = Lon(PI);

    /// Create a new longitude value with a given raw angle value in radians.
    /// The angle is wrapped into [-pi..pi].
    pub fn new(val: f64) -> Self {
        Self(wrap(val, -PI))
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

impl Add<f64> for Lon {
    type Output = Self;

    fn add(self, rhs: f64) -> Self::Output {
        Self::new(self.0 + rhs)
    }
}

impl Add<Lon> for f64 {
    type Output = Lon;

    fn add(self, rhs: Lon) -> Self::Output {
        rhs + self
    }
}

impl Neg for Lon {
    type Output = Lon;

    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl Sub<f64> for Lon {
    type Output = Self;

    fn sub(self, rhs: f64) -> Self::Output {
        // Watch out for infinite recursion with self - rhs
        Self::new(self.0 - rhs)
    }
}

impl Sub<Lon> for f64 {
    type Output = Lon;

    fn sub(self, rhs: Lon) -> Self::Output {
        Lon::new(self - rhs.0)
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

    pub fn length(&self) -> Option<f64> {
        let mut length = self.hi - self.lo;
        if length < 0. {
            length += 2. * PI;
        }
        if length > 0. {
            Some(length)
        } else {
            None
        }
    }
}

/// A latitude coordinate in [-pi/2..pi/2]  radians.
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq, Default)]
pub struct Lat(f64);

impl Lat {
    pub const MIN: Lat = Lat(-FRAC_PI_2);
    pub const MAX: Lat = Lat(FRAC_PI_2);

    /// Create a new latitude value with the given raw angle value in radians.
    /// The angle value is clamped into [-pi/2..pi/2].
    pub fn new(val: f64) -> Self {
        Lat(val.clamp(-FRAC_PI_2, FRAC_PI_2))
    }

    /// Return this latitude as a raw angle value in radians.
    #[inline]
    pub fn rad(self) -> f64 {
        self.0
    }
}

impl Add<f64> for Lat {
    type Output = Self;

    fn add(self, rhs: f64) -> Self::Output {
        Self::new(self.0 + rhs)
    }
}

impl Add<Lat> for f64 {
    type Output = Lat;

    fn add(self, rhs: Lat) -> Self::Output {
        Lat::new(self + rhs.0)
    }
}

impl Neg for Lat {
    type Output = Lat;

    fn neg(self) -> Self::Output {
        Lat::new(-self.0)
    }
}

impl Sub<f64> for Lat {
    type Output = Self;

    fn sub(self, rhs: f64) -> Self::Output {
        Self::new(self.0 - rhs)
    }
}

impl Sub<Lat> for f64 {
    type Output = Lat;

    fn sub(self, rhs: Lat) -> Self::Output {
        Lat::new(self - rhs.0)
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
    use crate::quantity::angle::units::DEG;

    #[test]
    fn test_lon() {
        // wrapping
        assert_eq!(Lon::new(2.0 * PI), Lon::new(0.0));
        assert_eq!(Lon::new(185.0 * DEG), Lon::new(-175.0 * DEG));
        // normalization
        assert_eq!(Lon::new(-PI).normalize(), Lon::new(PI));
        // equality
        assert_eq!(Lon::new(FRAC_PI_2), Lon::new(90.0 * DEG));
        assert_ne!(Lon::new(90. * DEG), Lon::new(91.0 * DEG));
        // ops
        assert_eq!(Lon::new(90. * DEG) + 45.0 * DEG, Lon::new(135.0 * DEG));
        assert_eq!(Lon::new(90. * DEG) - 45.0 * DEG, Lon::new(45.0 * DEG));
        // display
        assert_eq!(format!("{}", Lon::new(90. * DEG)), "  90° 00′ 00.00000000″");
    }

    #[test]
    fn test_lat() {
        // clamping
        assert_eq!(Lat::new(91. * DEG), Lat::MAX);
        assert_eq!(Lat::new(-90.01 * DEG), Lat::MIN);
        // equality
        assert_eq!(Lat::new(45. * DEG), Lat::new(FRAC_PI_4));
        assert_ne!(Lat::new(0. * DEG), Lat::new(0.001 * DEG));
        // ops
        assert_eq!(Lat::new(45. * DEG + 10. * DEG), Lat::new(55. * DEG));
        assert_eq!(Lat::new(45. * DEG) + 50. * DEG, Lat::MAX);
        assert_eq!(Lat::new(45. * DEG) - 90. * DEG, Lat::new(-45. * DEG));
        assert_eq!(Lat::new(-45. * DEG) - 50. * DEG, Lat::MIN);
        // Display
        assert_eq!(format!("{}", Lat::new(45. * DEG)), "  45° 00′ 00.00000000″");
    }

    #[test]
    fn test_lon_interval() {
        assert!(LonInterval::empty().is_empty());
        assert!(LonInterval::empty().is_inverted());
        assert!(LonInterval::full().is_full());

        // Length
        assert_eq!(LonInterval::empty().length(), None);
        assert_eq!(LonInterval::full().length(), Some(2. * PI));
        assert_eq!(
            LonInterval::new(Lon::new(0. * DEG), Lon::new(160. * DEG)).length(),
            Some(160. * DEG)
        );
        assert_abs_diff_eq!(
            LonInterval::new(Lon::new(170. * DEG), Lon::new(-170. * DEG))
                .length()
                .unwrap(),
            20. * DEG,
            epsilon = 1e-15
        );
    }
}
