/// A length unit is basically a conversion to meters.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct LengthUnit(pub f64, pub f64);

impl LengthUnit {
    #[inline]
    pub const fn m_per_unit(&self) -> f64 {
        self.0 / self.1
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
