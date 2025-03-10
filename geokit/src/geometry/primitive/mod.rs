//! [Primitive] is the base trait for geometric primitives, i.e. basic building
//! blocks for this model of geometric objects.

use super::Geometry;

/// Base trait for geometric *primitives*, i.e. geometries that serve as building
/// blocks for this model of geometric objects as set of positions.
///
/// Primitives are *open*, i.e. they do not include their boundary (if they have one).
///
/// Any geometry that is used to describe a feature is a collection of [Primitive].
/// A collection of geometry may or may not be a [Complex](crate::geometry::complex)
pub trait Primitive: Geometry {}

pub mod curve;
pub mod point;
