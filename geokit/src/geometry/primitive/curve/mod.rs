use super::Primitive;
use crate::{crs::Crs, math::fp::Float, quantities::length::Length};
use line_string::LineString;

pub trait GenericCurve {
    /// Returns the position's coordinates of the first point of the curve.
    fn start(&self) -> &[Float];

    /// Returns the position's coordinates of the last point of the curve.
    fn end(&self) -> &[Float];

    /// Returns the length of this curve in the given [Crs](crate::crs::Crs).
    fn length(&self, crs: &dyn Crs) -> Length;

    /// Returns the parameter for the `start()` point. For a [Curve] it must always be 0.
    /// For a [CurveSegment] within a [Curve], it should be the [param] of the segment's
    /// starting point in the parent [Curve].
    fn start_param(&self) -> Length;

    /// Returns the parameter for the `end()` point. For a [Curve] it must always be the length
    /// of the [Curve].
    /// For a [CurveSegment] within a [Curve], it should be the `param` of the segment's
    /// ending point in the parent [Curve].
    fn end_param(&self) -> Length;

    /// Returns the position on the [GenericCurve] at the given distance from the start.
    ///
    /// This represents the curve parameterization by arc length.
    fn param(&self, s: Length) -> &[Float];

    /// Constructs a [LineString] approximation of the curve where the control points are
    /// no more that `max_distance` apart and/or points on the line string are less than
    /// `max_offset` off from the original curve.
    fn as_line_string(
        &self,
        max_distance: Option<Length>,
        max_offset: Option<Length>,
    ) -> LineString;
}

/// Base trait for geometric curves.
///
/// A [Curve] is the image of an open interval by a continuous mapping.
pub trait Curve: Primitive + GenericCurve {}

pub mod line_string;
