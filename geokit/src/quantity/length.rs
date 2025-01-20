//! This module defines some types used to express length in various quantity.

use derive_more::derive::Display;

/// A 'Length' value.
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd, Display)]
#[display("{} m", _0)]
pub struct Length(f64);

impl Length {
    pub fn m(&self) -> f64 {
        self.0
    }
}

pub mod units {

    use std::ops::Mul;

    use super::Length;

    /// A length unit is basically a conversion to meters.
    #[derive(Debug, Copy, Clone, PartialEq)]
    pub struct LengthUnit(pub f64, pub f64);

    impl LengthUnit {
        pub fn m_per_unit(&self) -> f64 {
            self.0 / self.1
        }
    }

    impl Mul<LengthUnit> for f64 {
        type Output = Length;

        fn mul(self, rhs: LengthUnit) -> Self::Output {
            Length(self * rhs.0 / rhs.1)
        }
    }

    impl Mul<f64> for LengthUnit {
        type Output = Length;

        fn mul(self, rhs: f64) -> Self::Output {
            rhs * self
        }
    }

    /// 1 meter. The reference length unit.
    pub const M: LengthUnit = LengthUnit(1.0, 1.0);

    /// 1 international foot: `1 ft = 0.3048 m`.
    pub const FT: LengthUnit = LengthUnit(0.304_8, 1.0);

    /// 1 US survey foott: `1 ft = 0.3048006... m`.
    /// Beginning of Jan 1st 2023, the US survey foot should be avoided,
    /// except for historic and legacy applications.
    /// It has been superseded by the international foot, see [FT]
    pub const US_FT: LengthUnit = LengthUnit(1_200.0, 3_937.);
}
