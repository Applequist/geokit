use crate::{
    math::fp::Float,
    quantities::length::Length,
    units::length::{LengthUnit, M},
};

/// Errors used to compare approximate equality between cartesian coordinates.
///
/// See [XYZ::approx_eq]
/// See [ENH::approx_eq]
#[derive(Copy, Clone, Debug)]
pub struct CartesianErrors((Float, LengthUnit));

impl CartesianErrors {
    pub fn tiny() -> Self {
        CartesianErrors((1e-4, M))
    }

    pub fn small() -> Self {
        CartesianErrors((1e-2, M))
    }

    #[inline]
    pub fn length(&self) -> Length {
        self.0 .0 * self.0 .1
    }
}

impl Default for CartesianErrors {
    fn default() -> Self {
        CartesianErrors::tiny()
    }
}

pub mod geocentric;
pub mod projected;
