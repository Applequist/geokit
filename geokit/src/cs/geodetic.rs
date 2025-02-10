//! aeodetic coordinates spaces are used to represent points on or near the surface of an ellipsoid
//! of revolution.
//!
//! Points in those spaces are represented using:
//! - longitude,
//! - latitude,
//! - and in the 3D case, ellipsoidal height.
//!
//! A special set of axes is said to be **normalized**.
//!
//! This module a [system of axes][GeodeticAxes] for geodetic CS and
//! a set of value types to represent **normalized** geodetic coordinates.

use crate::math::Float;
use crate::quantities::angle::Angle;
use crate::quantities::length::Length;
use crate::quantities::Convertible;
use crate::units::angle::{AngleUnit, DEG, RAD};
use crate::units::length::{LengthUnit, M};
use approx::AbsDiffEq;
use derive_more::derive::{Display, Neg};
use std::ops::{Add, AddAssign, Div, Sub, SubAssign};

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
    /// Convert coordinates expressed in this system of axes into
    /// [**normalized** geodetic coordinates][LLH].
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

    /// Convert [**normalized geodetic coordinates][LLH] into coordinates
    /// expressed in this system of axes.
    pub fn denormalize(&self, llh: LLH, coords: &mut [Float]) {
        match self {
            GeodeticAxes::EastNorthUp {
                angle_unit,
                height_unit,
            } => {
                coords[0] = llh.lon.val(*angle_unit);
                coords[1] = llh.lat.val(*angle_unit);
                coords[2] = llh.height.val(*height_unit);
            }
            GeodeticAxes::EastNorth { angle_unit } => {
                coords[0] = llh.lon.val(*angle_unit);
                coords[1] = llh.lat.val(*angle_unit);
            }
            GeodeticAxes::NorthEastUp {
                angle_unit,
                height_unit,
            } => {
                coords[0] = llh.lat.val(*angle_unit);
                coords[1] = llh.lon.val(*angle_unit);
                coords[2] = llh.height.val(*height_unit);
            }
            GeodeticAxes::NorthEast { angle_unit } => {
                coords[0] = llh.lat.val(*angle_unit);
                coords[1] = llh.lon.val(*angle_unit);
            }
            GeodeticAxes::NorthWest { angle_unit } => {
                coords[0] = llh.lat.val(*angle_unit);
                coords[1] = -llh.lon.val(*angle_unit);
            }
        }
    }

    /// Return the dimension (2D or 3D) of the coordinate system.
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

pub struct GeodeticErrors {
    pub lon_err: Angle,
    pub lat_err: Angle,
    pub height_err: Length,
}

impl GeodeticErrors {
    pub const fn super_tiny() -> Self {
        GeodeticErrors {
            lon_err: Angle::super_tiny(),
            lat_err: Angle::super_tiny(),
            height_err: Length::super_tiny(),
        }
    }

    pub const fn tiny() -> Self {
        GeodeticErrors {
            lon_err: Angle::tiny(),
            lat_err: Angle::tiny(),
            height_err: Length::tiny(),
        }
    }

    pub const fn small() -> Self {
        GeodeticErrors {
            lon_err: Angle::small(),
            lat_err: Angle::small(),
            height_err: Length::small(),
        }
    }
}

impl Default for GeodeticErrors {
    fn default() -> Self {
        GeodeticErrors::super_tiny()
    }
}

/// [LLH] represents **normalized** geodeetic coordinates:
/// - longitude in radians, positive east,
/// - latitude in radians, positive north,
/// - ellipsoidal height in meters, positive upward.
#[derive(Debug, Copy, Clone, PartialEq, Display)]
#[display("(lon = {}, lat = {}, h = {})", lon, lat, height)]
pub struct LLH {
    pub lon: Lon,
    pub lat: Lat,
    pub height: Length,
}

impl LLH {
    /// Check whether this [LLH] is equal to the `other` [LLH] within the given
    /// [GeodeticErrors] bounds.
    pub fn approx_eq(&self, other: &LLH, err: GeodeticErrors) -> bool {
        self.lon.abs_diff_eq(&other.lon, err.lon_err)
            && self.lat.abs_diff_eq(&other.lat, err.lat_err)
            && self.height.abs_diff_eq(&other.height, err.height_err)
    }
}

/// [Lon] represents a longitude coordinate in [-pi..pi] radians.
///
/// # Operations
///
/// [Lon] supports the following operations
/// - Addition and subtraction of [Angle] to yield a new [Lon]
/// - Addition and subtraction of [Lon] to yield a new [Lon]: used for change or origin of
/// longitude.
/// - [Normalization][Self::normalized()] in (-pi..pi].
/// - [Conversion][Self::val()] to other [AngleUnit]
/// - Trigonometric operations like sin, cos, tan...
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd, Neg, Display)]
#[display("{}", _0.deg())]
pub struct Lon(Angle);

impl Lon {
    pub const MIN: Lon = Lon(Angle::M_PI);
    pub const ZERO: Lon = Lon(Angle::ZERO);
    pub const MAX: Lon = Lon(Angle::PI);

    /// Create a new longitude value from a given angle.
    /// The angle is wrapped into [-pi..pi].
    pub fn new(val: Angle) -> Self {
        Self(val.wrapped())
    }

    /// Create a new longitude value from an [Angle] ***ALREADY*** in [-pi..pi]
    pub const fn const_new(val: Angle) -> Self {
        Self(val)
    }

    pub fn deg(val_deg: Float) -> Self {
        Self::new(Angle::new(val_deg, DEG))
    }

    /// Create a new longitue value from a *dms* angle value.
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
    /// assert!(Lon::dms(-12., 45., 59.1234).rad() < 0.0);
    /// ```
    pub fn dms(d: Float, m: Float, s: Float) -> Self {
        debug_assert!(d > -180. && m <= 180., "degrees must be in (-180..180]");
        debug_assert!(m >= 0. && m <= 59.0, "minutes must be in [0..59]");
        debug_assert!(s >= 0. && s < 60., "seconds must be in [0..60)");
        let f = d.signum();
        let deg = d.abs() + m / 60. + s / 3600.;
        Self::new(f * deg * DEG)
    }

    /// Normalize the longitude into (-pi..pi] such that any point on a parallel
    /// has a unique *normalized* longitude.
    pub fn normalize(self) -> Self {
        if self <= Self::MIN {
            Self::MAX
        } else {
            self
        }
    }

    /// Return the angle of the longitude.
    #[inline]
    pub fn angle(self) -> Angle {
        self.0
    }

    /// Return the longitude as a raw angle value **in radians**.
    #[inline]
    pub fn rad(self) -> Float {
        self.0.rad()
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
}

impl Convertible for Lon {
    type Unit = AngleUnit;

    #[inline]
    fn val(self, unit: AngleUnit) -> Float {
        self.0.val(unit)
    }
}

impl Add for Lon {
    type Output = Lon;

    /// Add two longitudes.
    /// This is useful to change reference of longitude.
    ///
    /// # Example
    ///
    /// To compute the longitude wrt Greenwich prime meridian
    /// from a longitude `lon` wrt a prime meridian of Greenwich longitude `lon_pm`:
    /// ```
    /// use geokit::units::angle::DEG;
    /// use approx::assert_abs_diff_eq;
    /// let lon_pm = -30. * DEG; // Greenwich longitude of new prime meridian
    /// let lon = 10. * DEG; // longitude wrt prime meridian a Greenwich longitude `lon_pm`
    /// let converted_lon = lon + lon_pm; // Greenwich longitude
    /// assert_abs_diff_eq!(converted_lon, -20. * DEG, epsilon = 1e-15);
    /// ```
    fn add(self, rhs: Self) -> Self::Output {
        Lon::new(self.0 + rhs.0)
    }
}

impl Add<Angle> for Lon {
    type Output = Self;

    /// Return a new longitude `rhs` radians east (rhs > 0) or
    /// west (rhs < 0) of this longitude.
    /// The resulting longitude is [wrapped][Angle::wrapped] in [-pi..pi]
    fn add(self, rhs: Angle) -> Self::Output {
        Self::new(self.0 + rhs)
    }
}

impl Add<Lon> for Angle {
    type Output = Lon;

    /// Return a new longitude by adding `self` radians east (self > 0) or i
    /// west (self < 0) to this longitude and [wrapping][Angle::wrapped]
    /// the result in [-pi..pi].
    fn add(self, rhs: Lon) -> Self::Output {
        rhs + self
    }
}

impl AddAssign<Angle> for Lon {
    /// Modify this longitude by adding `rhs` radians to it and
    /// [wrapping][Angle::wrapped] the result in [-pi..pi]
    fn add_assign(&mut self, rhs: Angle) {
        self.0 = (self.0 + rhs).wrapped();
    }
}

impl Sub for Lon {
    type Output = Lon;

    /// Subtract two longitudes.
    /// This is useful to change reference of longitude.
    ///
    /// # Example
    ///
    /// To compute the longitude wrt a prime meridian of Greenwich longitude `lon_pm`
    /// from a Greenwich longitude `lon`:
    /// ```
    /// use geokit::units::angle::DEG;
    /// let lon = 10. * DEG; // Greenwich longitude
    /// let lon_pm = -30. * DEG; // Greenwich longitude of new prime meridian
    /// let converted_lon = lon - lon_pm;
    /// assert_eq!(converted_lon, 40. * DEG);
    /// ```
    fn sub(self, rhs: Self) -> Self::Output {
        Lon::new(self.0 - rhs.0)
    }
}

impl Sub<Angle> for Lon {
    type Output = Self;

    /// Return a new longitude by subtracting `rhs` radians from this
    /// longitude and [wrapping][Self::wrapped] the result in [-pi, pi].
    fn sub(self, rhs: Angle) -> Self::Output {
        // Watch out for infinite recursion with self - rhs
        Self::new(self.0 - rhs)
    }
}

impl SubAssign<Angle> for Lon {
    /// Modify this longitude by subtracting `rhs` radians from it
    /// and wrapping the result in [-pi, pi]
    fn sub_assign(&mut self, rhs: Angle) {
        self.0 = (self.0 - rhs).normalized();
    }
}

impl AbsDiffEq for Lon {
    type Epsilon = Angle;

    fn default_epsilon() -> Self::Epsilon {
        Angle::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        (*self - *other).0.abs() <= epsilon
    }
}

/// [Lat] represents a latitude coordinate (geodetic or other auxiliary ones)
/// in [-pi/2..pi/2] radians.
///
/// # Creation
///
/// There are 2 ways to create a [Lat]:
/// - by using [Lat::new] passing an [Angle] value. The
/// value is **clamped** to [-pi/2 .. pi/2]:
/// ```
/// use geokit::cs::geodetic::Lat;
/// use geokit::units::angle::DEG;
/// let lat = Lat::new(45. * DEG);
/// ```
///
/// - by using [Lat::dms] passing degrees/minutes/seconds values.
/// The value is also clamped to [-pi/2..pi/2]:
/// ```
/// use geokit::cs::geodetic::Lat;
/// let n = Lat::dms(44.0, 43., 27.183);
/// ```
///
/// # Operations
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

    /// Create a new latitude value.
    /// The angle is clamped in [-pi/2..pi/2]
    pub fn new(val: Angle) -> Self {
        Lat(val.clamped(Angle::M_PI_2, Angle::PI_2))
    }

    /// Create a new latitude value with an angle in degrees.
    /// The angle is clamped in [-pi/2..pi/2]
    pub fn deg(val_deg: Float) -> Self {
        Self::new(Angle::new(val_deg, DEG))
    }

    /// Create a new latitude value from a *dms* angle value.
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
    /// use geokit::cs::geodetic::Lat;
    /// let a: Lat = Lat::dms(2., 20., 14.02500);
    /// assert!(Lat::dms(-12., 45., 59.1234).rad() < 0.0);
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

    /// Return this latitude as a raw angle value in radians.
    #[inline]
    pub fn rad(self) -> Float {
        self.0.rad()
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
}

impl Convertible for Lat {
    type Unit = AngleUnit;

    #[inline]
    fn val(self, unit: AngleUnit) -> Float {
        self.0.val(unit)
    }
}

impl Add<Angle> for Lat {
    type Output = Self;

    fn add(self, rhs: Angle) -> Self::Output {
        Self::new(self.0 + rhs)
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
        self.0 = (self.0 + rhs).clamped(Angle::M_PI_2, Angle::PI_2)
    }
}

impl Sub for Lat {
    type Output = Angle;

    fn sub(self, rhs: Self) -> Self::Output {
        self.0 - rhs.0
    }
}

impl Sub<Angle> for Lat {
    type Output = Self;

    fn sub(self, rhs: Angle) -> Self::Output {
        Self::new(self.0 - rhs)
    }
}

impl SubAssign<Angle> for Lat {
    fn sub_assign(&mut self, rhs: Angle) {
        self.0 = (self.0 - rhs).clamped(Angle::M_PI_2, Angle::PI_2);
    }
}

impl Div<Float> for Lat {
    type Output = Lat;

    fn div(self, rhs: Float) -> Self::Output {
        Lat(self.0 / rhs)
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
    use crate::math::Float;
    use crate::units::angle::{DEG, RAD};
    use approx::assert_abs_diff_eq;
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
            let normalized = Lon::new(a).normalize();
            assert!(
                (normalized.angle() - e).rad().abs() < 1e-15,
                "case {}: expected {}, got {}",
                t,
                e,
                normalized
            );
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
        let a = Lon::new(10. * DEG);

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
