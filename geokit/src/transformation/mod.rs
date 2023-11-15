use std::fmt::Debug;

#[derive(Debug)]
pub enum TransformationError {}

pub trait Transformation: Debug {
    fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<(), TransformationError>;
    fn boxed_clone(&self) -> Box<dyn Transformation>;
}

pub trait InvertibleTransformation: Transformation {
    fn inverse(&self) -> Box<dyn InvertibleTransformation>;
    fn boxed_clone(&self) -> Box<dyn InvertibleTransformation>;
}

impl Clone for Box<dyn Transformation> {
    fn clone(&self) -> Box<dyn Transformation> {
        Transformation::boxed_clone(&**self)
    }
}

impl Clone for Box<dyn InvertibleTransformation> {
    fn clone(&self) -> Box<dyn InvertibleTransformation> {
        InvertibleTransformation::boxed_clone(&**self)
    }
}

#[derive(Debug)]
pub struct Identity;

impl Transformation for Identity {
    fn apply(&self, input: &[f64], output: &mut [f64]) -> Result<(), TransformationError> {
        output.copy_from_slice(input);
        Ok(())
    }

    fn boxed_clone(&self) -> Box<dyn Transformation> {
        Box::new(Identity)
    }
}

impl InvertibleTransformation for Identity {
    fn inverse(&self) -> Box<dyn InvertibleTransformation> {
        Box::new(Identity)
    }

    fn boxed_clone(&self) -> Box<dyn InvertibleTransformation> {
        Box::new(Identity)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn box_clone() {
        let mut out = [0.; 3];

        let boxed_transfo: Box<dyn Transformation> = Box::new(Identity);
        let _cloned_as_transfo: Box<dyn Transformation> = boxed_transfo.clone();
        _cloned_as_transfo.apply(&[1., 2., 3.], &mut out).unwrap();
        assert_eq!(
            out,
            [1., 2., 3.],
            "Expected {:?}. Got {:?}.",
            [1., 2., 3.],
            out
        );

        let boxed_inv_transfo: Box<dyn InvertibleTransformation> = Box::new(Identity);
        let _cloned_as_inv_transfo: Box<dyn InvertibleTransformation> = boxed_inv_transfo.clone();
        let inv_t = _cloned_as_inv_transfo.inverse();
        inv_t.apply(&[0., 0., 0.], &mut out).unwrap();
        assert_eq!(out, [0.; 3], "Expected {:?}. Got {:?}", [0.; 3], out);
    }
}
