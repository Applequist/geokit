use super::{Curve, GenericCurve};
use crate::{
    crs::Crs,
    geometry::{
        Geometry,
        complex::CurveBoundary,
        primitive::{Primitive, point::Point},
    },
    math::fp::Float,
    quantities::length::Length,
};

#[derive(Clone, Debug, PartialEq)]
pub struct LineString {
    coords_dim: usize,
    pos: Vec<Float>,
}

impl LineString {
    /// Returns the number of control points defining this [LineString]
    pub(crate) fn len(&self) -> usize {
        self.pos.len() / self.coords_dim
    }

    /// Returns the coordinates of the control point's position
    /// at index `index`.
    pub(crate) fn pos(&self, index: usize) -> &[Float] {
        let start = self.coords_dim * index;
        let end = start + self.coords_dim;
        &self.pos[start..end]
    }
}

impl Geometry for LineString {
    fn dim(&self) -> usize {
        1
    }

    fn coord_dim(&self) -> usize {
        self.coords_dim
    }

    fn is_cycle(&self) -> bool {
        let start = self.pos(0);
        let end = self.pos(self.len() - 1);
        start == end
    }

    fn boundary(&self) -> Option<Box<dyn Geometry>> {
        if self.is_cycle() {
            None
        } else {
            Some(Box::new(CurveBoundary {
                start: Point::new(self.pos(0), self.coords_dim),
                end: Point::new(self.pos(self.len() - 1), self.coords_dim),
            }))
        }
    }
}

impl Primitive for LineString {}

impl GenericCurve for LineString {
    fn start(&self) -> &[Float] {
        todo!()
    }

    fn end(&self) -> &[Float] {
        todo!()
    }

    fn length(&self, crs: &dyn Crs) -> Length {
        todo!()
    }

    fn start_param(&self) -> Length {
        todo!()
    }

    fn end_param(&self) -> Length {
        todo!()
    }

    fn param(&self, s: Length) -> &[Float] {
        todo!()
    }

    fn as_line_string(
        &self,
        max_distance: Option<Length>,
        max_offset: Option<Length>,
    ) -> LineString {
        todo!()
    }
}

impl Curve for LineString {}

#[derive(Clone, Debug)]
pub struct LineStringBuilder<const D: usize> {
    coords: Vec<Float>,
}

impl<const D: usize> LineStringBuilder<D> {
    /// Creates a new [LineStringBuilder] with the initial line
    /// segment from `start` to `to`.
    pub fn new(start: [Float; D], to: [Float; D]) -> Self {
        let mut coords = vec![];
        coords.extend(start);
        coords.extend(to);
        Self { coords }
    }

    /// Adds a new line segment for the end point of the last segment
    /// to `p`.
    pub fn line_to<P: AsRef<[Float]>>(mut self, p: P) -> Self {
        let s = &p.as_ref()[..D];
        self.coords.extend(s);
        self
    }

    /// Adds a copy of the start position to close the [LineString]
    /// and make it.
    pub fn close(mut self) -> LineString {
        let mut s = [0.; D];
        s.copy_from_slice(&self.coords[..D]);
        self.line_to(s).make()
    }

    /// Returns the [LineString] made from this builder.
    pub fn make(self) -> LineString {
        LineString {
            coords_dim: D,
            pos: self.coords,
        }
    }
}
