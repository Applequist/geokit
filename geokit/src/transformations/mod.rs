use crate::math::fp::Float;
use derive_more::derive::Display;
use thiserror::Error;

#[derive(Debug, Error, Display)]
pub enum TransformationError {
    OutOfBounds,
}

/// Trait for coordinates transformation between a *source* and a *destination* [coordinates reference systems](crate::crs::Crs)
pub trait CrsTransformation {
    /// Transforms *source* coordinates from `src` into *destination* coordinates in `dst`.
    fn fwd(&self, src: &[Float], dst: &mut [Float]) -> Result<(), TransformationError>;

    /// Transforms *destination* coordinates from `dst` into *source* coordinates in `src`.
    fn bwd(&self, dst: &[Float], src: &mut [Float]) -> Result<(), TransformationError>;
}

pub mod xyz;
