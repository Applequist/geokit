use crate::math::Float;

pub enum OperationError {
    InputOutOfBounds,
    OutputOutOfBounds,
}

pub trait Operation {
    fn transform(&self, input: &[Float], output: &mut [Float]) -> Result<(), OperationError>;
}

