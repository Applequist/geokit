/// The length modue defines constants to convert various length units to meters:
/// ```
/// use geokit::quantity::length::units::FT;
/// let one_us_foot_in_meters = 1. * FT;
/// ```
///
/// All raw length values in functions/methods parameters or return values, variables are **in meters**.
pub mod length;

/// The scale module defines
pub mod scale;
