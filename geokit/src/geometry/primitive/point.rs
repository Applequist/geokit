use crate::{geometry::Geometry, math::fp::Float};

pub struct Point {
    dim: usize,
    coords: [Float; 3],
}

impl Point {
    pub fn new(buf: &[Float], dim: usize) -> Self {
        let mut coords = [0.0; 3];
        coords[..dim].clone_from_slice(&buf[..dim]);
        Self { dim, coords }
    }
}

impl From<[Float; 2]> for Point {
    fn from(value: [Float; 2]) -> Self {
        Self::new(&value, 2)
    }
}

impl From<[Float; 3]> for Point {
    fn from(value: [Float; 3]) -> Self {
        Self::new(&value, 3)
    }
}

impl Geometry for Point {
    fn dim(&self) -> usize {
        0
    }

    fn coord_dim(&self) -> usize {
        self.dim
    }

    fn is_cycle(&self) -> bool {
        true
    }

    fn boundary(&self) -> Option<Box<dyn Geometry>> {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::Point;
    use crate::geometry::Geometry;
    use crate::math::fp::Float;
    use std::any::Any;

    #[test]
    fn size() {
        let p2d: Point = [1.0, 2.0].into();
        assert_eq!(
            size_of_val(&p2d),
            std::mem::size_of::<usize>() + 3 * std::mem::size_of::<Float>()
        );

        let p3d: Point = [1., 2., 3.].into();
        assert_eq!(
            size_of_val(&p3d),
            std::mem::size_of::<usize>() + 3 * std::mem::size_of::<Float>()
        );
    }

    #[test]
    fn type_id() {
        let p2d: &Point = &[1.0, 2.0].into();
        assert_eq!(p2d.type_id(), std::any::TypeId::of::<Point>());

        let p3d: &Point = &[1., 2., 3.].into();
        assert_eq!(p3d.type_id(), std::any::TypeId::of::<Point>());
    }

    #[test]
    fn geometry() {
        let p2d: Point = [1., 2.].into();
        assert_eq!(p2d.dim(), 0);
        assert_eq!(p2d.coord_dim(), 2);
        assert_eq!(p2d.is_cycle(), true);
        assert!(p2d.boundary().is_none());

        let p3d: Point = [1., 2., 3.].into();
        assert_eq!(p3d.dim(), 0);
        assert_eq!(p3d.coord_dim(), 3);
        assert_eq!(p3d.is_cycle(), true);
        assert!(p3d.boundary().is_none());
    }
}
