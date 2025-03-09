//! This module provides geometries.

use std::any::Any;

/// The base trait for all **geometric objects**, finite representation of possibly
/// infinite set of positions in space.
pub trait Geometry: Any {
    /// Returns the *inherent* dimension of this geometry, which should be
    /// less than or equal to its coordinates dimension.
    ///
    /// The dimension of a collection of geometry should be the highest dimension
    /// of any of its elements.
    fn dim(&self) -> usize;

    /// Returns the dimension of the coordinates that defines this geometry,
    /// which should be the same of the dimension of the CRS it is
    /// defined in.
    fn coord_dim(&self) -> usize;

    /// Returns whether this geometry is *simple*.
    ///
    /// A geomtry is simple each point in its interior has a neighborhood that is
    /// isomorphic to a n-sphere where n is the [dimension](Geometry::dim) of this
    /// geometry.
    /// In simple terms, it is simple if it has no self-intersection or self-tangency.
    ///
    /// This library assumes that all geometries it deals with are simple. Hence the
    /// default implementation returns `true`.
    fn is_simple(&self) -> bool {
        true
    }

    /// Returns `true` if the boundary of this geometry is empty after topological
    /// simplification.
    ///
    /// A geometry is a cycle if it is isomorphic to a geometry that is the boundary
    /// of a region of some euclidian space: A point is a cycle, a curve is a cycle
    /// if it is closed, ie. isomorphic to a circle, a surface is a cycle if it is
    /// isomorphic to a sphere or a torus. A solid in dimension 3 is never a cycle.
    fn is_cycle(&self) -> bool;

    /// Returns a geometry that contains all the points on the boundary
    /// of this geometry if it has one.
    fn boundary(&self) -> Option<Box<dyn Geometry>>;
}

pub mod primitive;
