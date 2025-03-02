use crate::{
    math::fp::Float,
    quantities::length::Length,
    units::length::{LengthUnit, M},
};

/// Errors used to compare approximate equality between cartesian coordinates.
///
/// See [approx_eq](crate::cs::cartesian::geocentric::XYZ::approx_eq)
#[derive(Copy, Clone, Debug)]
pub struct CartesianTolerance {
    pub all: (Float, LengthUnit),
}

impl CartesianTolerance {
    /// [CartesianTolerance::tiny()] allows a `1e-4 m` error on each axis,
    /// which translates to less that a millimeters positional error.
    pub fn tiny() -> Self {
        CartesianTolerance { all: (1e-4, M) }
    }

    /// [CartesianTolerance::small()] allows a `1e-3 m` error on each axis,
    /// which translates to less than a centimeter positional error.
    pub fn small() -> Self {
        CartesianTolerance { all: (1e-3, M) }
    }

    /// Returns the maximum error allowed on each axis.
    #[inline]
    pub fn length(&self) -> Length {
        self.all.0 * self.all.1
    }

    /// Converts this [CartesianTolerance] into the given unit.
    pub(crate) fn convert_to(&self, unit: LengthUnit) -> Self {
        let all = if self.all.1 == unit {
            self.all
        } else {
            (self.length().val(unit), unit)
        };

        Self { all }
    }
}

impl Default for CartesianTolerance {
    fn default() -> Self {
        CartesianTolerance::tiny()
    }
}

pub mod geocentric;
pub mod projected;
