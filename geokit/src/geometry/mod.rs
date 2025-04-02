//! This module provides geometries.

use primitive::Boundary;
use std::any::Any;
use std::fmt::Debug;

/// Defines the various *geometry type*, a combination of
/// geometric dimension and coordinate dimension.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GeometryType {
    Empty,
    Point(usize),
    Curve(usize),
    Surface(usize),
    Solid,
}

impl GeometryType {
    /// Returns the geometric dimension of a geometry type.
    fn dim(&self) -> Option<usize> {
        match self {
            GeometryType::Empty => None,
            GeometryType::Point(_) => Some(0),
            GeometryType::Curve(_) => Some(1),
            GeometryType::Surface(_) => Some(2),
            GeometryType::Solid => Some(3),
        }
    }

    /// Returns the coordinate dimension of a geometry type.
    fn coord_dim(&self) -> Option<usize> {
        match self {
            GeometryType::Empty => None,
            GeometryType::Point(d) | GeometryType::Curve(d) | GeometryType::Surface(d) => Some(*d),
            GeometryType::Solid => Some(3),
        }
    }
}

/// [Geometry] is the base trait for all **geometric objects** that are finite representations
/// of possibly infinite set of [positions](coordinate::Pos) in a given [Crs](crate::crs::Crs).
pub trait Geometry: Any + Debug {
    /// Returns the [GeometryType] of this geometry.
    fn geometry_type(&self) -> GeometryType;

    /// Returns whether this geometry is *simple*.
    ///
    /// A geomtry is simple each point in its interior has a neighborhood that is
    /// isomorphic to a n-sphere where n is the [dimension](Geometry::dim) of this
    /// geometry.
    /// In simple terms, it is simple if it has no self-intersection or self-tangency.
    ///
    /// WARN: This library assumes that all geometries it deals with are simple. Hence the
    /// default implementation returns `true`.
    fn is_simple(&self) -> bool {
        true
    }

    /// Returns `true` if the boundary of this geometry is empty after topological
    /// simplification.
    ///
    /// A geometry is a cycle if it is isomorphic to a geometry that is the boundary
    /// of a region of some euclidian space:
    /// - a point is a cycle,
    /// - a curve is a cycle if it is closed, ie. isomorphic to a circle,
    /// - a surface is a cycle if it is isomorphic to a sphere or a torus.
    /// - a solid in dimension 3 is never a cycle.
    fn is_cycle(&self) -> bool;

    /// Returns a geometry that contains all the points on the boundary
    /// of this geometry if it has one.
    fn boundary(&self) -> Option<Box<dyn Boundary>>;
}

pub mod complex;
pub mod coordinate;
pub mod empty;
pub mod primitive;
