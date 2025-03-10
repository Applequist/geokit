use crate::quantities::length::Length;

/// Errors used to compare approximate equality between normalizeld cartesian coordinates.
///
/// See [approx_eq](crate::cs::cartesian::geocentric::XYZ::approx_eq)
#[derive(Copy, Clone, Debug)]
pub struct CartesianTolerance {
    pub all: Length,
}

impl CartesianTolerance {
    /// [CartesianTolerance::tiny()] allows a `1e-4 m` error on each axis,
    /// which translates to less that a millimeters positional error.
    pub const fn tiny() -> Self {
        CartesianTolerance { all: Length::TINY }
    }

    /// [CartesianTolerance::small()] allows a `1e-3 m` error on each axis,
    /// which translates to less than a centimeter positional error.
    pub const fn small() -> Self {
        CartesianTolerance { all: Length::SMALL }
    }
}

impl Default for CartesianTolerance {
    fn default() -> Self {
        CartesianTolerance::tiny()
    }
}

pub mod geocentric;
pub mod projected;
