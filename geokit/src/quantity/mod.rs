//! The quantity module defines types, constants and methods to work with various quantities, eg angle, length...

/// The angle module defines constants to convert various angle units to radians:
/// ```
/// use geokit::quantity::angle::units::DEG;
/// let one_deg_in_rad = 1. * DEG;
/// ```
///
/// All raw angle values in functions/methods parameters or return values, variables are ***in radians***.
/// When working with angles which have a special meaning and/or whose valid representation is constrained
/// are represented by a wrapper types, eg [Lon], [Lat] and [Azimuth]:
/// ```
/// use std::f64::consts::FRAC_PI_2;
/// use geokit::cs::geodetic::Lon;
/// let east = Lon::new(FRAC_PI_2);
/// ```
#[macro_use]
pub mod angle;

/// The length modue defines constants to convert various length units to meters:
/// ```
/// use geokit::quantity::length::units::FT;
/// let one_us_foot_in_meters = 1. * FT;
/// ```
///
/// All raw length values in functions/methods parameters or return values, variables are **in meters**.
pub mod length;
pub mod scale;
