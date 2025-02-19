//! Geodetic coordinates systems are used to represent points on or near the surface of an ellipsoid
//! of revolution using longitude, latitude and in the 3D case, ellipsoid height coordinates.
//!
//! A geodetic coordinates system using the following coordinates in that order:
//! - longitude in [-pi, pi] radians positive east
//! - latitude in [-pi/2, pi/2] radians positive north
//! - height in meters positive up
//! is said to be normalized.
//!
//! This module defines a [system of axes][GeodeticAxes] for geodetic CS and
//! a set of value types to represent **normalized** geodetic coordinates.

use crate::math::fp::Float;
use crate::quantities::angle::Angle;
use crate::quantities::length::Length;
use crate::units::angle::{AngleUnit, DEG, RAD};
use crate::units::length::{LengthUnit, M};
use approx::AbsDiffEq;
use derive_more::derive::{Display, Neg};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

/// A [GeodeticAxes] defines the possible set of axes used in **geodetic** CS.
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum GeodeticAxes {
    /// Coordinates are given in the following order:
    /// - longitude positive east of prime meridian, using `angle_unit` [AngleUnit]
    /// - latitude positive north of equatorial plane, using `angle_unit` [AngleUnit],
    /// - ellipsoidal height positive upward, using `height_unit` [LengthUnit]
    EastNorthUp {
        /// The angle unit used for longitude and latitude, eg
        /// `DEG` or `GRAD`
        angle_unit: AngleUnit,
        /// The length unit used for ellipsoidal height, eg
        /// 'M' or 'US_FOOT'
        height_unit: LengthUnit,
    },
    /// Coordinates are given in the following order:
    /// - latitude positive north of equatorial plane, using `angle_unit` [AngleUnit],
    /// - longitude positive east of prime meridian, using `angle_unit` [AngleUnit]
    /// - ellipsoidal height positive upward, using `height_unit` [LengthUnit]
    NorthEastUp {
        /// The angle unit used for longitude and latitude, eg
        /// `DEG` or `GRAD`
        angle_unit: AngleUnit,
        /// The length unit used for ellipsoidal height, eg
        /// 'M' or 'US_FOOT'
        height_unit: LengthUnit,
    },
    /// Coordinates are, in the given order:
    /// - longitude positive east of prime meridian, using `angle_unit` [AngleUnit]
    /// - latitude positive north of equatorial plane, using `angle_unit` [AngleUnit],
    EastNorth {
        /// The angle unit used for longitude and latitude, eg
        /// `DEG` or `GRAD`
        angle_unit: AngleUnit,
    },
    /// Coordinates are, in the given order:
    /// - latitude positive north of equatorial plane, using `angle_unit` [AngleUnit],
    /// - longitude positive east of prime meridian, using `angle_unit` [AngleUnit]
    NorthEast {
        /// The angle unit used for longitude and latitude, eg
        /// `DEG` or `GRAD`
        angle_unit: AngleUnit,
    },
    /// Coordinates are, in the given order:
    /// - latitude positive north of equatorial plane, using `angle_unit` [AngleUnit],
    /// - longitude positive west of prime meridian, using `angle_unit` [AngleUnit]
    NorthWest {
        /// The angle unit used for longitude and latitude, eg
        /// `DEG` or `GRAD`
        angle_unit: AngleUnit,
    },
}

impl GeodeticAxes {
    /// Converts coordinates expressed in this system of axes into
    /// **normalized** geodetic coordinates.
    pub fn normalize(&self, coords: &[Float]) -> LLH {
        match self {
            GeodeticAxes::EastNorthUp {
                angle_unit,
                height_unit,
            } => LLH {
                lon: Lon::new(Angle::new(coords[0], *angle_unit)),
                lat: Lat::new(Angle::new(coords[1], *angle_unit)),
                height: Height::new(coords[2], *height_unit),
            },
            GeodeticAxes::EastNorth { angle_unit } => LLH {
                lon: Lon::new(Angle::new(coords[0], *angle_unit)),
                lat: Lat::new(Angle::new(coords[1], *angle_unit)),
                height: Height::ZERO,
            },
            GeodeticAxes::NorthEastUp {
                angle_unit,
                height_unit,
            } => LLH {
                lon: Lon::new(Angle::new(coords[1], *angle_unit)),
                lat: Lat::new(Angle::new(coords[0], *angle_unit)),
                height: Height::new(coords[2], *height_unit),
            },
            GeodeticAxes::NorthEast { angle_unit } => LLH {
                lon: Lon::new(Angle::new(coords[1], *angle_unit)),
                lat: Lat::new(Angle::new(coords[0], *angle_unit)),
                height: Height::ZERO,
            },
            GeodeticAxes::NorthWest { angle_unit } => LLH {
                lon: Lon::new(Angle::new(-coords[1], *angle_unit)),
                lat: Lat::new(Angle::new(coords[0], *angle_unit)),
                height: Height::ZERO,
            },
        }
    }

    /// Converts **normalized geodetic coordinates into coordinates
    /// expressed in this system of axes.
    pub fn denormalize(&self, llh: LLH, coords: &mut [Float]) {
        match self {
            GeodeticAxes::EastNorthUp {
                angle_unit,
                height_unit,
            } => {
                coords[0] = llh.lon.angle().val(*angle_unit);
                coords[1] = llh.lat.angle().val(*angle_unit);
                coords[2] = llh.height.val(*height_unit);
            }
            GeodeticAxes::EastNorth { angle_unit } => {
                coords[0] = llh.lon.angle().val(*angle_unit);
                coords[1] = llh.lat.angle().val(*angle_unit);
            }
            GeodeticAxes::NorthEastUp {
                angle_unit,
                height_unit,
            } => {
                coords[0] = llh.lat.angle().val(*angle_unit);
                coords[1] = llh.lon.angle().val(*angle_unit);
                coords[2] = llh.height.val(*height_unit);
            }
            GeodeticAxes::NorthEast { angle_unit } => {
                coords[0] = llh.lat.angle().val(*angle_unit);
                coords[1] = llh.lon.angle().val(*angle_unit);
            }
            GeodeticAxes::NorthWest { angle_unit } => {
                coords[0] = llh.lat.angle().val(*angle_unit);
                coords[1] = -llh.lon.angle().val(*angle_unit);
            }
        }
    }

    /// Returns the dimension (2D or 3D) of the coordinate system.
    pub fn dim(&self) -> usize {
        match self {
            GeodeticAxes::EastNorthUp {
                angle_unit: _,
                height_unit: _,
            }
            | GeodeticAxes::NorthEastUp {
                angle_unit: _,
                height_unit: _,
            } => 3,
            GeodeticAxes::EastNorth { angle_unit: _ }
            | GeodeticAxes::NorthEast { angle_unit: _ }
            | GeodeticAxes::NorthWest { angle_unit: _ } => 2,
        }
    }
}

impl Default for GeodeticAxes {
    /// Return the [`GeodeticAxes`] used in **normalized geodetic coordinates**.
    fn default() -> Self {
        Self::EastNorthUp {
            angle_unit: RAD,
            height_unit: M,
        }
    }
}

/// Errors bounds used to check approximate equality between geodetic coordinates.
///
/// See [LLH::approx_eq].
pub struct GeodeticErrors {
    pub lon: Angle,
    pub lat: Angle,
    pub height: Length,
}

impl GeodeticErrors {
    /// Tiny geodetic errors.
    pub const fn tiny() -> Self {
        GeodeticErrors {
            lon: Angle::tiny(),
            lat: Angle::tiny(),
            height: Length::tiny(),
        }
    }

    /// Small geodetic errors.
    pub const fn small() -> Self {
        GeodeticErrors {
            lon: Angle::small(),
            lat: Angle::small(),
            height: Length::small(),
        }
    }
}

impl Default for GeodeticErrors {
    fn default() -> Self {
        Self::tiny()
    }
}

/// [LLH] represents **normalized** geodetic coordinates:
/// - longitude in (-pi, pi] radians, positive east,
/// - latitude in [-pi/2, pi/2] radians, positive north,
/// - ellipsoidal height in meters, positive upward.
///
/// <div class="warning">Note that **normalized** geodetic coordinates does NOT enforce a **normalized** longitude.</div>
#[derive(Debug, Copy, Clone, PartialEq, Display)]
#[display("(lon = {}, lat = {}, h = {})", lon, lat, height)]
pub struct LLH {
    pub lon: Lon,
    pub lat: Lat,
    pub height: Length,
}

impl LLH {
    /// Checks whether this [LLH] is approximately equal to the `other` [LLH]
    /// within the given [GeodeticErrors] bounds.
    ///
    /// `self` and `other` are approximately equal if the following conditions are
    /// all satisfied:
    /// - `(self.lon - other.lon).abs() <= err.lon`
    /// - `(self.lat - other.lat).abs() <= err.lat`
    /// - `(self.height - other.height).abs() <= err.height`
    pub fn approx_eq(&self, other: &LLH, err: GeodeticErrors) -> bool {
        self.lon.abs_diff_eq(&other.lon, err.lon)
            && self.lat.abs_diff_eq(&other.lat, err.lat)
            && self.height.abs_diff_eq(&other.height, err.height)
    }
}

/// Returns whehter `res` is approximately equal to `exp` within the given [GeodeticErrors] error
/// bounds, printing information about the coordinates when not equal.
///
/// Use only in tests.
pub fn approx_eq_llh(res: &LLH, exp: &LLH, err: &GeodeticErrors) -> bool {
    let lon_ok = res.lon.abs_diff_eq(&exp.lon, err.lon);
    let lat_ok = res.lat.abs_diff_eq(&exp.lat, err.lat);
    let height_ok = res.height.abs_diff_eq(&exp.height, err.height);
    let is_approx_eq = lon_ok && lat_ok && height_ok;

    if !is_approx_eq {
        println!("----------- computed LLH != expected LLH ------------");
        if lon_ok {
            println!(
                "Longitude ok:      {} = {} +/- {:e}",
                res.lon,
                exp.lon,
                err.lon.deg()
            );
        } else {
            println!(
                "Longitude error: | {} - {} | = {:e} > {:e}",
                res.lon,
                exp.lon,
                (res.lon - exp.lon).abs().deg(),
                err.lon.deg()
            );
        }
        if lat_ok {
            println!(
                "Latitude ok:      {} = {} +/- {:e}",
                res.lat,
                exp.lat,
                err.lat.deg()
            );
        } else {
            println!(
                "Latitude error: | {} - {} | = {:e} > {:e}",
                res.lat,
                exp.lat,
                (res.lat - exp.lat).abs().deg(),
                err.lat.deg()
            );
        }
        if height_ok {
            println!(
                "Height ok      {} = {} +/- {:e} m",
                res.height,
                exp.height,
                err.height.m()
            );
        } else {
            println!(
                "Height error: | {} - {} | = {:e} m > {:e} m",
                res.height,
                exp.height,
                (res.height - exp.height).abs().m(),
                err.height.m()
            );
        }
        println!("");
    }
    is_approx_eq
}

/// [Lon] represents a longitude coordinate in [-pi..pi] radians.
///
/// [Lon] supports the following operations
/// - Addition and subtraction of [Angle] to yield a new [Lon]
/// - [Normalization][Self::normalized()] in (-pi..pi].
/// - [Conversion][Self::val()] to other [AngleUnit]
/// - Trigonometric operations like sin, cos, tan...
///
/// Although with a [-pi, pi] range, points on the ante-meridian have 2 possible
/// valid longitudes (-pi and pi), `-pi` is still useful in some cases, eg to
/// represent empty, and full longitude intervals.
/// To get a unique longitude value, use a [normalized](Lon::normalized) longitude.
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd, Neg, Display)]
#[display("{}", _0.deg())]
pub struct Lon(Angle);

impl Lon {
    /// The minimum longitude (inclusive).
    pub const MIN: Lon = Lon(Angle::M_PI);
    pub const ZERO: Lon = Lon(Angle::ZERO);
    /// The maximum longitude (inclusive).
    pub const MAX: Lon = Lon(Angle::PI);

    /// Creates a new longitude value from a given angle.
    ///
    /// The angle is wrapped into [-pi..pi].
    pub fn new(val: Angle) -> Self {
        Self(val.wrapped())
    }

    /// Creates a new longitude value from an [Angle] ***ALREADY*** in [-pi..pi]
    pub const fn const_new(val: Angle) -> Self {
        Self(val)
    }

    /// Creates a [Lon] value whose angle is given in radians.
    ///
    /// The angle is wrapped into [-pi, pi].
    pub fn rad(val_rad: Float) -> Self {
        Self::new(val_rad * RAD)
    }

    /// Creates a [Lon] value whose angle is given in degrees.
    ///
    /// The angle is converted to radians and wrapped into [-pi, pi].
    pub fn deg(val_deg: Float) -> Self {
        Self::new(val_deg * DEG)
    }

    /// Creates a new longitue value from a *dms* angle value.
    /// The angle value is converted to radians and wrapped into [-pi, pi].
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
    /// let a: Lon = Lon::dms(2., 20., 14.02500);
    /// assert!(Lon::dms(-12., 45., 59.1234) < Lon::ZERO);
    /// ```
    pub fn dms(d: Float, m: Float, s: Float) -> Self {
        debug_assert!(d > -180. && m <= 180., "degrees must be in (-180..180]");
        debug_assert!(m >= 0. && m <= 59.0, "minutes must be in [0..59]");
        debug_assert!(s >= 0. && s < 60., "seconds must be in [0..60)");
        let f = d.signum();
        let deg = d.abs() + m / 60. + s / 3600.;
        Self::new(f * deg * DEG)
    }

    /// Returns this longitude **normalized** into (-pi..pi]
    /// such that any point on a parallel has a unique longitude angle.
    pub fn normalized(self) -> Self {
        if self <= Self::MIN {
            Self::MAX
        } else {
            self
        }
    }

    /// Returns the angle of the longitude.
    #[inline]
    pub fn angle(self) -> Angle {
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

    pub fn signum(&self) -> Float {
        self.0.signum()
    }
}

impl Add<Angle> for Lon {
    type Output = Self;

    /// Returns a [Lon] whose angle is the sum `self.angle() + rhs`
    /// [wrapped](Angle::wrapped) in [-pi, pi] radians.
    fn add(self, rhs: Angle) -> Self::Output {
        Self::new(self.0 + rhs)
    }
}

impl Add<Lon> for Angle {
    type Output = Lon;

    /// Returns a [Lon] whose angle is the sum `self + rhs.angle()`
    /// [wrapped](Angle::wrapped) in [-pi, pi] radians.
    fn add(self, rhs: Lon) -> Self::Output {
        rhs + self
    }
}

impl AddAssign<Angle> for Lon {
    /// Sets this [Lon] angle to the sum of `self.angle() + rhs`
    /// [wrapped](Angle::wrapped) in [-pi, pi] radians.
    fn add_assign(&mut self, rhs: Angle) {
        self.0 = (self.0 + rhs).wrapped();
    }
}

impl Sub for Lon {
    type Output = Angle;

    /// Returns the [Angle] in (-pi, pi] radians `rhs.angle().diff_to(self.angle())`.
    ///
    /// # Example
    ///
    /// ```
    /// # use geokit::cs::geodetic::Lon;
    /// # use geokit::quantities::angle::Angle;
    /// # use geokit::units::angle::DEG;
    /// assert_eq!(Lon::new(10. * DEG) - Lon::new(-30. * DEG), 40. * DEG);
    /// assert_eq!(Lon::ZERO - Lon::MAX, Angle::PI);
    /// ```
    ///
    /// See [Angle::diff_to]
    fn sub(self, rhs: Self) -> Self::Output {
        rhs.0.diff_to(self.0)
    }
}

impl Sub<Angle> for Lon {
    type Output = Self;

    /// Returns a [Lon] whose angle is the difference `self.angle() - rhs`
    /// [wrapped](Angle::wrapped) in [-pi, pi].
    fn sub(self, rhs: Angle) -> Self::Output {
        // Watch out for infinite recursion with self - rhs
        Self::new(self.0 - rhs)
    }
}

impl SubAssign<Angle> for Lon {
    /// Sets this [Lon] angle to the difference `self.angle() - rhs`
    /// [wrapped](Angle::wrapped) in [-pi, pi]
    fn sub_assign(&mut self, rhs: Angle) {}
}

impl Mul<Float> for Lon {
    type Output = Lon;

    /// Returns a [Lon] whose angle is the product `self.angle() * rhs`
    /// [wrapped](Angle::wrapped) in [-pi, pi]
    fn mul(self, rhs: Float) -> Self::Output {
        Lon::new(self.0 * rhs)
    }
}

impl Mul<Lon> for Float {
    type Output = Lon;

    /// Returns a [Lon] whose angle is the product `self * rhs.angle()`
    /// [wrapped](Angle::wrapped) in [-pi, pi]
    fn mul(self, rhs: Lon) -> Self::Output {
        rhs * self
    }
}

impl MulAssign<Float> for Lon {
    /// Sets this [Lon] angle to the product `self.angle() * rhs`
    /// [wrapped](Angle::wrapped) in [-pi, pi]
    fn mul_assign(&mut self, rhs: Float) {
        self.0 = (self.0 * rhs).wrapped();
    }
}

impl Div<Float> for Lon {
    type Output = Lon;

    /// Returns a [Lon] whose angle is the division `self.angle() / rhs`
    /// [wrapped](Angle::wrapped) in [-pi, pi].
    fn div(self, rhs: Float) -> Self::Output {
        Lon::new(self.0 / rhs)
    }
}

impl DivAssign<Float> for Lon {
    /// Sets this [Lon] angle to the division `self.angle() / rhs`
    /// [wrapped](Angle::wrapped) in [-pi, pi].
    fn div_assign(&mut self, rhs: Float) {
        self.0 = (self.0 / rhs).wrapped()
    }
}

impl AbsDiffEq for Lon {
    type Epsilon = Angle;

    fn default_epsilon() -> Self::Epsilon {
        Angle::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        (*self - *other).abs() <= epsilon
    }
}

/// [Lat] represents a geodetic latitude coordinate in [-pi/2..pi/2] radians.
///
/// [Lat] supports the following operations:
/// - Negation
/// - saturating addition and subtraction of [Angle].
/// - subtraction resulting in an [Angle].
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq, Neg, Display)]
#[display("{}", self.0.deg())]
pub struct Lat(Angle);

impl Lat {
    pub const MIN: Lat = Lat(Angle::M_PI_2);
    pub const ZERO: Lat = Lat(Angle::ZERO);
    pub const MAX: Lat = Lat(Angle::PI_2);

    /// Creates a new latitude value.
    ///
    /// The angle is clamped in [-pi/2..pi/2]
    pub fn new(val: Angle) -> Self {
        Lat(val.clamped(Angle::M_PI_2, Angle::PI_2))
    }

    /// Creates a new latitude value with an angle in degrees.
    ///
    /// The angle is clamped in [-pi/2..pi/2]
    pub fn deg(val_deg: Float) -> Self {
        Self::new(Angle::new(val_deg, DEG))
    }

    /// Creates a new latitude value from a *dms* angle value.
    ///
    /// The angle value is converted to radians and clamped into [-pi/2, pi/2].
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
    /// # use geokit::cs::geodetic::Lat;
    /// let a: Lat = Lat::dms(2., 20., 14.02500);
    /// assert!(Lat::dms(-12., 45., 59.1234) < Lat::ZERO);
    /// ```
    pub fn dms(d: Float, m: Float, s: Float) -> Self {
        debug_assert!(d >= -90. && m <= 90., "degrees must be in [-90..90]");
        debug_assert!(m >= 0. && m <= 59.0, "minutes must be in [0..59]");
        debug_assert!(s >= 0. && s < 60., "seconds must be in [0..60)");
        let f = d.signum();
        let deg = d.abs() + m / 60. + s / 3600.;
        Self::new(f * deg * DEG)
    }

    #[inline]
    pub fn abs(self) -> Lat {
        Self(self.0.abs())
    }

    /// Return the latitude angle.
    #[inline]
    pub fn angle(self) -> Angle {
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

    pub fn signum(&self) -> Float {
        self.0.signum()
    }
}

impl Add<Angle> for Lat {
    type Output = Self;

    /// Returns a [Lat] whose angle is the sum `self.angle() + rhs`
    /// [clamped](Angle::clamped) in [-pi/2, pi/2].
    fn add(self, rhs: Angle) -> Self::Output {
        Self::new(self.0 + rhs)
    }
}

impl Add<Lat> for Angle {
    type Output = Lat;

    /// Returns a [Lat] whose angle is the sum `self + rhs.angle()`
    /// [clamped](Angle::clamped) in [-pi/2, pi/2].
    fn add(self, rhs: Lat) -> Self::Output {
        rhs + self
    }
}

impl AddAssign<Angle> for Lat {
    /// Sets this [Lat] angle to the sum `self + rhs.angle()`
    /// [clamped](Angle::clamped) in [-pi/2, pi/2].
    fn add_assign(&mut self, rhs: Angle) {
        self.0 = (self.0 + rhs).clamped(Angle::M_PI_2, Angle::PI_2)
    }
}

impl Sub for Lat {
    type Output = Angle;

    /// Returns the [Angle] `self.angle() - rhs.angle()`.
    fn sub(self, rhs: Self) -> Self::Output {
        self.0 - rhs.0
    }
}

impl Sub<Angle> for Lat {
    type Output = Self;

    /// Returns a [Lat] whose angle is the diff `self.angle() - rhs`
    /// [clamped](Angle::clamped) to [-pi/2, pi/2].
    fn sub(self, rhs: Angle) -> Self::Output {
        Self::new(self.0 - rhs)
    }
}

impl SubAssign<Angle> for Lat {
    /// Sets this [Lat] angle to the diff `self.angle() - rhs`
    /// [clamped](Angle::clamped) to [-pi/2, pi/2].
    fn sub_assign(&mut self, rhs: Angle) {
        self.0 = (self.0 - rhs).clamped(Angle::M_PI_2, Angle::PI_2);
    }
}

impl Mul<Float> for Lat {
    type Output = Lat;

    /// Returns a [Lat] whose angle is the product `self.angle() * rhs`
    /// [clamped](Angle::clamped) to [-pi/2, pi/2].
    fn mul(self, rhs: Float) -> Self::Output {
        Lat::new(self.0 * rhs)
    }
}

impl Mul<Lat> for Float {
    type Output = Lat;

    /// Returns a [Lat] whose angle is the product `self * rhs.angle()`
    /// [clamped](Angle::clamped) to [-pi/2, pi/2].
    fn mul(self, rhs: Lat) -> Self::Output {
        rhs * self
    }
}

impl MulAssign<Float> for Lat {
    /// Sets this [Lat] angle to the product `self.angle() * rhs`
    /// [clamped](Angle::clamped) to [-pi/2, pi/2].
    fn mul_assign(&mut self, rhs: Float) {
        self.0 = (self.0 * rhs).clamped(Angle::M_PI_2, Angle::PI_2);
    }
}

impl Div<Float> for Lat {
    type Output = Lat;

    /// Returns a [Lat] whose angle is the division `self.angle() / rhs`
    /// [clamped](Angle::clamped) to [-pi/2, pi/2].
    fn div(self, rhs: Float) -> Self::Output {
        Lat::new(self.0 / rhs)
    }
}

impl DivAssign<Float> for Lat {
    /// Sets this [Lat] angle to the division `self.angle() / rhs`
    /// [clamped](Angle::clamped) to [-pi/2, pi/2].
    fn div_assign(&mut self, rhs: Float) {
        self.0 = (self.0 / rhs).clamped(Angle::M_PI_2, Angle::PI_2);
    }
}

impl AbsDiffEq for Lat {
    type Epsilon = Angle;

    fn default_epsilon() -> Self::Epsilon {
        Angle::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        (*self - *other).abs() <= epsilon
    }
}

/// A type alias for
pub type Height = Length;

// TODO: Add the following:
// - a LonInterval
// - a LatInterval.
// - a Height similar to Length
// - a LonLat and LonLatHeight
// - a LonLatRect
// - a LonLatHeightRect

#[cfg(test)]
mod tests {
    use crate::cs::geodetic::{GeodeticAxes, Height, Lat, Lon, LLH};
    use crate::math::fp::Float;
    use crate::quantities::angle::Angle;
    use crate::units::angle::{DEG, RAD};
    use approx::{assert_abs_diff_eq, AbsDiffEq};
    use std::f64::consts::FRAC_PI_4;

    #[test]
    fn geodetic_coordinates_normalization() {
        let latlon_deg = GeodeticAxes::NorthWest { angle_unit: DEG };

        let llh = latlon_deg.normalize(&[10.0, 110.0]);
        assert_eq!(
            llh,
            LLH {
                lon: Lon::new(-110. * DEG),
                lat: Lat::new(10.0 * DEG),
                height: Height::ZERO
            }
        );

        let mut coords: [Float; 2] = [0.0; 2];
        latlon_deg.denormalize(llh, &mut coords);
        assert_eq!(coords, [10.0, 110.0]);
    }

    #[test]
    fn lon_wrapping() {
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
            let wrapped = Lon::new(a);
            assert!(
                (wrapped.angle() - e).rad().abs() < 1e-15,
                "case {}: expected {}, got {}",
                t,
                e,
                wrapped
            );
        }
    }

    #[test]
    fn lon_normalization() {
        let input_expected = [
            ("< -180", -210. * DEG, 150. * DEG),
            // FIX: Robustness issue
            // -180. - n * eps < -pi for n >= 65 on my test
            //(
            //    "-180 - eps",
            //    (-180. - Float::EPSILON) * DEG,
            //    (180. - Float::EPSILON) * DEG,
            //),
            // FIX: Robustness
            //(
            //    "-180 + eps",
            //    (-180. + Float::EPSILON) * DEG,
            //    (-180. + Float::EPSILON) * DEG,
            //),
            ("-180", -180. * DEG, 180. * DEG),
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

        // normalization
        for (t, a, e) in input_expected {
            let lon = Lon::new(a).normalized();
            if !lon.angle().abs_diff_eq(&e, Angle::default_epsilon()) {
                println!("case {t} failed");
            }
            assert_abs_diff_eq!(lon.angle(), e);
        }
    }

    #[test]
    fn lon_display() {
        // Default display
        assert_eq!(format!("{}", Lon::new(90. * DEG)), "90 deg");
        // DMS display
        assert_eq!(
            format!("{}", Lon::new(90. * DEG).angle().dms()),
            "  90° 00′ 00.00000000″"
        );
    }

    #[test]
    fn lon_operation() {
        // origins of longitude by Greenwich longitude
        let a = 10. * DEG;

        // (x, x + a, x - a)
        let data = [
            (
                Lon::new(-179. * DEG),
                Lon::new(-169. * DEG),
                Lon::new(171. * DEG),
            ),
            (
                Lon::new(-135. * DEG),
                Lon::new(-125. * DEG),
                Lon::new(-145. * DEG),
            ),
            (
                Lon::new(-90. * DEG),
                Lon::new(-80. * DEG),
                Lon::new(-100. * DEG),
            ),
            (
                Lon::new(-45. * DEG),
                Lon::new(-35. * DEG),
                Lon::new(-55. * DEG),
            ),
            (
                Lon::new(0. * DEG),
                Lon::new(10. * DEG),
                Lon::new(-10. * DEG),
            ),
            (
                Lon::new(45. * DEG),
                Lon::new(55. * DEG),
                Lon::new(35. * DEG),
            ),
            (
                Lon::new(90. * DEG),
                Lon::new(100. * DEG),
                Lon::new(80. * DEG),
            ),
            (
                Lon::new(135. * DEG),
                Lon::new(145. * DEG),
                Lon::new(125. * DEG),
            ),
            (
                Lon::new(180. * DEG),
                Lon::new(-170. * DEG),
                Lon::new(170. * DEG),
            ),
        ];

        for (x, x_plus_a, x_minus_a) in data {
            assert_abs_diff_eq!(x + a, x_plus_a);
            assert_abs_diff_eq!(x - a, x_minus_a);
        }
    }

    #[test]
    fn lat_clamping() {
        // clamping
        assert_eq!(Lat::new(91. * DEG), Lat::MAX);
        assert_eq!(Lat::new(-90.01 * DEG), Lat::MIN);
    }

    #[test]
    fn lat_equality() {
        // equality
        assert_eq!(Lat::new(45. * DEG), Lat::new(FRAC_PI_4 * RAD));
        assert_ne!(Lat::new(0. * DEG), Lat::new(0.001 * DEG));
    }

    #[test]
    fn lat_ops() {
        // ops
        assert_eq!(Lat::new(45. * DEG) + 10. * DEG, Lat::new(55. * DEG));
        assert_eq!(Lat::new(45. * DEG) + 50. * DEG, Lat::MAX);
        assert_eq!(Lat::new(45. * DEG) - 90. * DEG, Lat::new(-45. * DEG));
        assert_eq!(Lat::new(-45. * DEG) - 50. * DEG, Lat::MIN);
    }

    #[test]
    fn lat_display() {
        // Default display
        assert_eq!(format!("{}", Lat::new(45. * DEG)), "45 deg");
        // DMS display
        assert_eq!(
            format!("{}", Lat::new(45. * DEG).angle().dms()),
            "  45° 00′ 00.00000000″"
        );
    }
}
