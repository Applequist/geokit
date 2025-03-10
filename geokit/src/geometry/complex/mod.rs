//! Defines [Complex] geometric objects.

use super::{Geometry, primitive::point::Point};

/// A [Complex] is a set of [Primitive](crate::geometry::primitive) whose interiors are disjoints,
/// and such that for each primitive in the complex, there exists a set of primitive in the complex
/// whose union is the boundary of that primitive.
pub trait Complex: Geometry {}

pub trait PrimitiveBoundary: Complex {}

/// A [Complex] used to describe a [Curve](crate::geometry::primitive::curve)'s boundary if any.
#[derive(Copy, Clone, Debug)]
pub struct CurveBoundary {
    pub start: Point,
    pub end: Point,
}

impl Geometry for CurveBoundary {
    fn dim(&self) -> usize {
        0
    }

    fn coord_dim(&self) -> usize {
        self.start.dim()
    }

    fn is_cycle(&self) -> bool {
        self.start == self.end
    }

    fn boundary(&self) -> Option<Box<dyn Geometry>> {
        None
    }
}

impl Complex for CurveBoundary {}
