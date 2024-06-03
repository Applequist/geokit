//! This module defines some types used to express length in various quantity.

pub mod units {
    /// 1 meter. The reference length unit.
    pub const M: f64 = 1.0;

    /// 1 international foot: `1 ft = 0.3048 m`.
    pub const FT: f64 = 0.3048;

    /// 1 US survey foott: `1 ft = 0.3048006... m`.
    /// Beginning of Jan 1st 2023, the US survey foot should be avoided,
    /// except for historic and legacy applications.
    /// It has been superseded by the international foot, see [FT]
    pub const US_FT: f64 = 1200. / 3937.;
}
