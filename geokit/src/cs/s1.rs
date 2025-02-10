//! The S1 space is used to represent points on the unit circle.
//! Each point on the circle is represented by the angle at the center of the circle between
//! an origin and the point, in the range [-pi..pi].
//!
//! This module provides the following value types:
//! - [Angle] to represent a single S1 point coordinate
//! - and [interval][Interval] to represent a segment on S1.
//!
//! These types are used to define other value types in more concrete CS like 2D- and 3D-
//! geodetic CS.

use crate::quantities::angle::Angle;
use derive_more::derive::Display;

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
    /// FIX: Both singleton and empty interval have zero length
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

    use crate::{cs::s1::Interval, quantities::angle::Angle, units::angle::DEG};

    fn check_interval(
        tag: &str,
        ab: Interval,
        empty: bool,
        full: bool,
        inverted: bool,
        length: Angle,
    ) {
        assert_eq!(ab.is_empty(), empty, "Expected {} empty", ab);
        assert_eq!(ab.is_full(), full, "Expected {} full", ab);
        assert_eq!(ab.is_inverted(), inverted, "Expected {} inverted", ab);
        assert_abs_diff_eq!(ab.length(), length);
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
        let ab = Interval::new(5. * DEG, 15. * DEG);
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
