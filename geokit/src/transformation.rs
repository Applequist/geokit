use crate::coord::Coord3D;
use std::fmt::Debug;

pub enum TransformationError {}

pub trait Transformation: Debug {
    fn apply(&self, input: &Coord3D) -> Result<Coord3D, TransformationError>;
    fn boxed_clone(&self) -> Box<dyn Transformation>;
}

pub trait InversibleTransformation: Transformation {
    fn inverse(&self) -> Box<dyn Transformation>;
    fn boxed_clone(&self) -> Box<dyn InversibleTransformation>;
}

impl Clone for Box<dyn Transformation> {
    fn clone(&self) -> Box<dyn Transformation> {
        Transformation::boxed_clone(&**self)
    }
}

impl Clone for Box<dyn InversibleTransformation> {
    fn clone(&self) -> Box<dyn InversibleTransformation> {
        InversibleTransformation::boxed_clone(&**self)
    }
}
