use std::ops::{Deref, DerefMut};

use super::Primitive;
use crate::{
    geometry::{Geometry, GeometryType, coordinate::Pos, primitive::Boundary},
    math::fp::Float,
};

/// A [Point] is a 0-dimensional [Primitive] geometric object
/// containing a single [position](crate::geometry::Pos)
/// and whose boundary is empty.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Point {
    coord_dim: usize,
    coord: [Float; 3],
}

impl Point {
    /// Creates a new [Point] at the given [position](Pos) by copying
    /// the position's coordinates.
    ///
    /// # Panics
    ///
    /// This method panics if the [Pos] length is not 2 or 3.
    pub fn new(pos: &Pos) -> Self {
        let coord_dim = pos.len();
        assert!(coord_dim == 2 || coord_dim == 3);
        let mut coord = [0.0; 3];
        coord[..coord_dim].clone_from_slice(&pos);
        Self { coord_dim, coord }
    }
}

impl Geometry for Point {
    fn geometry_type(&self) -> GeometryType {
        GeometryType::Point(self.coord_dim)
    }

    fn is_cycle(&self) -> bool {
        true
    }

    fn boundary(&self) -> Option<Box<dyn Boundary>> {
        None
    }
}

impl Primitive for Point {}

impl Deref for Point {
    type Target = Pos;

    fn deref(&self) -> &Self::Target {
        &self.coord[0..self.coord_dim]
    }
}

impl DerefMut for Point {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.coord[0..self.coord_dim]
    }
}

impl<T> AsRef<T> for Point
where
    T: ?Sized,
    <Point as Deref>::Target: AsRef<T>,
{
    fn as_ref(&self) -> &T {
        self.deref().as_ref()
    }
}

impl<T> AsMut<T> for Point
where
    <Point as Deref>::Target: AsMut<T>,
{
    fn as_mut(&mut self) -> &mut T {
        self.deref_mut().as_mut()
    }
}

#[cfg(test)]
mod tests {
    use super::Point;
    use crate::geometry::{Geometry, GeometryType};
    use crate::math::fp::Float;
    use std::any::Any;

    #[test]
    fn size() {
        let p2d: Point = Point::new(&[1.0, 2.0]);
        assert_eq!(
            size_of_val(&p2d),
            std::mem::size_of::<usize>() + 3 * std::mem::size_of::<Float>()
        );

        let p3d: Point = Point::new(&[1., 2., 3.]);
        assert_eq!(
            size_of_val(&p3d),
            std::mem::size_of::<usize>() + 3 * std::mem::size_of::<Float>()
        );
    }

    #[test]
    fn type_id() {
        let p2d: &Point = &Point::new(&[1.0, 2.0]);
        assert_eq!(p2d.type_id(), std::any::TypeId::of::<Point>());

        let p3d: &Point = &Point::new(&[1., 2., 3.]);
        assert_eq!(p3d.type_id(), std::any::TypeId::of::<Point>());
    }

    #[test]
    fn geometry() {
        let p2d: Point = Point::new(&[1., 2.]);
        assert_eq!(p2d.geometry_type(), GeometryType::Point(2));
        assert_eq!(p2d.is_cycle(), true);
        assert!(p2d.boundary().is_none());

        let p3d: Point = Point::new(&[1., 2., 3.]);
        assert_eq!(p3d.geometry_type(), GeometryType::Point(3));
        assert_eq!(p3d.is_cycle(), true);
        assert!(p3d.boundary().is_none());
    }
}
