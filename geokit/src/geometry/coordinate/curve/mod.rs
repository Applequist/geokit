use crate::{
    crs::Crs, geometry::primitive::curve::CurveBoundary, math::fp::Float,
    quantities::length::Length,
};
use dyn_clone::DynClone;
use line_string::LineString;
use std::fmt::Debug;

use super::Pos;

/// A [ParameterizedCurve] defines a curve parameterization.
pub trait ParameterizedCurve {
    fn coord_dim(&self) -> usize;

    /// Returns the start [position](Pos) of the curve.
    fn start(&self) -> &Pos;

    /// Returns the end [position](Pos)of the curve.
    fn end(&self) -> &Pos;

    /// Returns the length of this curve.
    fn length(&self) -> Length;

    /// Returns the [position](Pos) on this curve at the given distance from the [Self::start()].
    ///
    /// This represents the curve parameterization by arc length.
    /// The following must hold:
    /// - `c.param(0) == c.start()`
    /// - `c.param(c.length()) == c.end()`
    fn param(&self, s: Length) -> Vec<Float>;

    /// Constructs a [LineString] approximation of the curve where the control points are
    /// no more that `max_distance` apart and/or points on the line string are less than
    /// `max_offset` off from the original curve.
    fn as_line_string(
        &self,
        crs: &dyn Crs,
        max_distance: Option<Length>,
        max_offset: Option<Length>,
    ) -> LineString;
}

pub enum CurveInterpolation {
    /// Linear interpolation returns positions on a straight line between
    /// consecutive pairs of control points.
    Linear,
    /// Geodesic interpolation returns positions on a geodesic line between
    /// consecutive pairs of control points.
    Geodesic,
}

/// A [CurveSegment] defines an homogeneous segment, i.e. section, of a [Curve].
///
/// [Curve]: crate::geometry::primitive::curve::Curve
pub trait CurveSegment: ParameterizedCurve + Debug + DynClone {
    /// Returns the curve interpolation method used for this segment.
    fn interpolation(&self) -> CurveInterpolation;

    /// Returns the boundary of this segment.
    fn curve_boundary(&self) -> CurveBoundary;
}

dyn_clone::clone_trait_object!(CurveSegment);

pub mod line_string;
