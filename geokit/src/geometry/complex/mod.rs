//! Defines [Complex] geometric objects.

use super::Geometry;

/// A [Complex] is a set of [Primitive](crate::geometry::primitive) whose interiors are disjoints,
/// and such that for each primitive in the complex, there exists a set of primitive in the complex
/// whose union is the boundary of that primitive.
pub trait Complex: Geometry {}
