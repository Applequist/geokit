use std::any::Any;

use super::complex::Complex;
use crate::geometry::{Geometry, GeometryType, primitive::Boundary};

/// The special empty geometry.
#[derive(Clone, Copy, Debug)]
pub struct Empty;

impl Geometry for Empty {
    fn geometry_type(&self) -> GeometryType {
        GeometryType::Empty
    }

    fn is_cycle(&self) -> bool {
        true
    }

    fn boundary(&self) -> Box<dyn Boundary> {
        Box::new(Empty)
    }
}

impl Complex for Empty {}
impl Boundary for Empty {}
