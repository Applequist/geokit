use crate::cs::azimuth::Azimuth;
use crate::cs::geodetic::{Lat, Lon};
use crate::math::Float;
use crate::units::length::LengthUnit;
use approx::AbsDiffEq;
use derive_more::derive::{Add, AddAssign, Display, Neg, Sub, SubAssign};
use std::fmt::{Display, Formatter};
use std::ops::{Div, Mul};

/// [CurveLength] is a length measure along a curve in meters.
#[derive(
    Debug, Copy, Clone, PartialEq, PartialOrd, Display, Add, AddAssign, Sub, SubAssign, Neg,
)]
#[display("{} m", _0)]
pub struct CurveLength(Float);

impl CurveLength {
    pub const ZERO: CurveLength = CurveLength(0.0);

    pub fn new(qty: Float, unit: LengthUnit) -> Self {
        Self(qty * unit.m_per_unit())
    }

    pub fn val(self, unit: LengthUnit) -> Float {
        self.0 / unit.m_per_unit()
    }

    pub fn m(self) -> Float {
        self.0
    }
}

impl Mul<Float> for CurveLength {
    type Output = CurveLength;

    fn mul(self, rhs: Float) -> Self::Output {
        Self(self.0 * rhs)
    }
}

impl Mul<CurveLength> for Float {
    type Output = CurveLength;

    fn mul(self, rhs: CurveLength) -> Self::Output {
        rhs * self
    }
}

impl Div for CurveLength {
    type Output = Float;

    fn div(self, rhs: Self) -> Self::Output {
        self.0 / rhs.0
    }
}

impl Div<Float> for CurveLength {
    type Output = CurveLength;

    fn div(self, rhs: Float) -> Self::Output {
        Self(self.0 / rhs)
    }
}

impl AbsDiffEq for CurveLength {
    type Epsilon = Float;

    fn default_epsilon() -> Self::Epsilon {
        Float::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Geodesic {
    pub p1: (Lon, Lat),
    pub alpha1: Azimuth,
    pub p2: (Lon, Lat),
    pub alpha2: Azimuth,
    pub s: CurveLength,
}

impl Default for Geodesic {
    fn default() -> Self {
        Geodesic {
            p1: (Lon::ZERO, Lat::ZERO),
            alpha1: Azimuth::NORTH,
            p2: (Lon::ZERO, Lat::ZERO),
            alpha2: Azimuth::NORTH,
            s: CurveLength::ZERO,
        }
    }
}

impl Display for Geodesic {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "From: lon = {} lat = {} az = {}\n",
            self.p1.0, self.p1.1, self.alpha1
        )?;
        write!(
            f,
            "To  : lon = {} lat = {} az = {}\n",
            self.p2.0, self.p2.1, self.alpha2
        )?;
        // WARN: longitude difference meaning is still not clear
        write!(f, "L   : {}\n", (self.p2.0 - self.p1.0).to_dms())?;
        write!(f, "s   : {:0.3}\n", self.s)
    }
}

pub trait GeodesicSolver {
    /// Compute the coordinates and forward azimuth of the point at `s12` meters away
    /// from `p1` following the geodesic in the azimuth `alpha1` at `p1`.
    ///
    /// Some implementations may not converge for some inputs where the resulting geodesic starts
    /// and finishes at *near* antipodal points.
    ///
    /// # Parameters
    ///
    /// - ellipsoid: the ellipsoid on which the geodesic is solved
    /// - p1: the **normalized geodetic** coordinates of the starting point
    /// - alpha1: the azimuth **in radians** of the geodesic at `p1`,
    /// - s12: the distance **in meters** along the geodesic from `p1` of the returned point coordinates.
    /// TODO: s12 is not really a r1::Length! More like a s2::ArcLength
    fn solve_direct(
        &self,
        p1: (Lon, Lat),
        alpha1: Azimuth,
        s12: CurveLength,
    ) -> Result<Geodesic, &'static str>;

    /// Compute the azimuths and length in meters of the geodesic from `p1` to `p2`.
    ///
    /// Some implementations may not converge when inputs are *near* antipodal points.
    fn solve_inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Result<Geodesic, &'static str>;
}

pub mod karney;
pub mod rapp;
pub mod vincenty;

#[cfg(test)]
mod tests {

    use crate::cs::azimuth::Azimuth;
    use crate::cs::geodetic::{Lat, Lon};
    use crate::geodesy::ellipsoid::consts;
    use crate::geodesy::geodesics::{CurveLength, Geodesic};
    use crate::geodesy::Ellipsoid;
    use crate::units::angle::DEG;
    use crate::units::length::M;

    /// errors in the 5th decimal of a second
    pub struct DirectDeltas {
        pub delta_lat2: f64,
        pub delta_delta_lon: f64,
        pub delta_alpha2: f64,
    }

    pub struct InverseDeltas {
        /// errors in the 5th decimal of a second
        pub delta_alpha1: f64,
        /// errors in the 5th decimal of a second
        pub delta_alpha2: f64,
        /// errors in millimeters
        pub delta_s: f64,
    }

    pub struct LineData {
        pub ellipsoid: Ellipsoid,
        pub geodesic: Geodesic,
    }

    /// From Vincenty - Direct and inverse solutions of geodesics on the ellipsoid with
    ///                 application of nested equations
    /// Table II - Results of solutions
    pub fn vincenty_lines() -> Vec<LineData> {
        vec![
            // Line (a)
            LineData {
                ellipsoid: consts::BESSEL,
                geodesic: Geodesic {
                    p1: (Lon::ZERO, Lat::dms(55., 45., 0.)),
                    alpha1: Azimuth::dms(96., 36., 8.79960),
                    p2: (Lon::dms(108., 13., 0.), Lat::dms(-33., 26., 0.)),
                    alpha2: Azimuth::dms(137., 52., 22.01454),
                    s: CurveLength::new(14_110_526.170, M),
                },
            },
            // Line (b)
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::ZERO, Lat::dms(37., 19., 54.95367)),
                    alpha1: Azimuth::dms(95., 27., 59.63089),
                    p2: (Lon::dms(41., 28., 35.50729), Lat::dms(26., 7., 42.83946)),
                    alpha2: Azimuth::dms(118., 5., 58.96161),
                    s: CurveLength::new(4_085_966.703, M),
                },
            },
            // Line (c)
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::ZERO, Lat::dms(35., 16., 11.24862)),
                    alpha1: Azimuth::dms(15., 44., 23.74850),
                    p2: (Lon::dms(137., 47., 28.31435), Lat::dms(67., 22., 14.77638)),
                    alpha2: Azimuth::dms(144., 55., 39.92147),
                    s: CurveLength::new(8_084_823.839, M),
                },
            },
            // Line (d)
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::ZERO, Lat::dms(1., 0., 0.)),
                    alpha1: Azimuth::dms(89., 0., 0.),
                    p2: (Lon::dms(179., 17., 48.02997), Lat::dms(-0., 59., 53.83076)),
                    alpha2: Azimuth::dms(91., 0., 6.11733),
                    s: CurveLength::new(19960000.0, M),
                },
            },
            // Line (e)
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::ZERO, Lat::dms(1., 0., 0.)),
                    alpha1: Azimuth::dms(4., 59., 59.99995),
                    p2: (Lon::dms(179., 46., 17.84244), Lat::dms(1., 1., 15.18952)),
                    alpha2: Azimuth::dms(174., 59., 59.88481),
                    s: CurveLength::new(19_780_006.558, M),
                },
            },
        ]
    }

    pub fn vincenty_direct_deltas() -> Vec<DirectDeltas> {
        vec![
            // Line (a)
            DirectDeltas {
                delta_lat2: -1.2,
                delta_delta_lon: 0.7,
                delta_alpha2: -1.2,
            },
            // Line (b)
            DirectDeltas {
                delta_lat2: -0.7,
                delta_delta_lon: 1.2,
                delta_alpha2: 0.5,
            },
            // Line (c)
            DirectDeltas {
                delta_lat2: -2.0,
                delta_delta_lon: 2.9,
                delta_alpha2: 3.0,
            },
            // Line (d)
            DirectDeltas {
                delta_lat2: -0.2,
                delta_delta_lon: 0.6,
                delta_alpha2: -0.3,
            },
            // Line (e)
            DirectDeltas {
                delta_lat2: 2.5,
                delta_delta_lon: -0.2,
                delta_alpha2: -0.3,
            },
        ]
    }

    pub fn vincenty_inverse_deltas() -> Vec<InverseDeltas> {
        vec![
            // Line (a)
            InverseDeltas {
                delta_alpha1: -0.4,
                delta_alpha2: -0.5,
                delta_s: -0.4,
            },
            // Line (b)
            InverseDeltas {
                delta_alpha1: -0.2,
                delta_alpha2: -0.2,
                delta_s: -0.4,
            },
            // Line (c)
            InverseDeltas {
                delta_alpha1: -0.2,
                delta_alpha2: 0.3,
                delta_s: -0.7,
            },
            // Line (d)
            InverseDeltas {
                delta_alpha1: -102.9,
                delta_alpha2: 102.9,
                delta_s: -0.2,
            },
            // Line (e)
            InverseDeltas {
                delta_alpha1: 0.4,
                delta_alpha2: -0.8,
                delta_s: 0.8,
            },
        ]
    }
    // From Rapp - Geometric Geodesy 1.71 Standard Test Lines

    pub fn standard_lines() -> Vec<LineData> {
        vec![
            // Line 1
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::new(0. * DEG), Lat::dms(37., 19., 54.95367)),
                    alpha1: Azimuth::dms(95., 27., 59.630888),
                    p2: (Lon::dms(41., 28., 35.50729), Lat::dms(26., 7., 42.83946)),
                    alpha2: Azimuth::dms(118., 5., 58.961608),
                    s: CurveLength::new(4_085_966.7026, M),
                },
            },
            // Line 2
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::new(0. * DEG), Lat::dms(35., 16., 11.24862)),
                    alpha1: Azimuth::dms(15., 44., 23.748498),
                    p2: (Lon::dms(137., 47., 28.31435), Lat::dms(67., 22., 14.77638)),
                    alpha2: Azimuth::dms(144., 55., 39.921473),
                    s: CurveLength::new(8_084_823.8383, M),
                },
            },
            // Line 3
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::new(0. * DEG), Lat::dms(1., 0., 0.)),
                    alpha1: Azimuth::dms(88., 59., 59.998970),
                    p2: (Lon::dms(179., 17., 48.02997), Lat::dms(-0., 59., 53.83076)),
                    alpha2: Azimuth::dms(91., 0., 6.118357),
                    s: CurveLength::new(19_959_999.9998, M),
                },
            },
            // Line 4
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::new(0. * DEG), Lat::dms(1., 0., 0.)),
                    alpha1: Azimuth::dms(4., 59., 59.999953),
                    p2: (Lon::dms(179., 46., 17.84244), Lat::dms(1., 1., 15.18952)),
                    alpha2: Azimuth::dms(174., 59., 59.884804),
                    s: CurveLength::new(19_780_006.5588, M),
                },
            },
            // Line 5
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::new(0. * DEG), Lat::dms(41., 41., 45.88000)),
                    alpha1: Azimuth::dms(52., 40., 39.390667),
                    p2: (Lon::dms(0., 0., 0.56000), Lat::dms(41., 41., 46.20000)),
                    alpha2: Azimuth::dms(52., 40., 39.763168),
                    s: CurveLength::new(16.2839751, M),
                },
            },
            // Line 6
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::new(0. * DEG), Lat::dms(30., 0., 0.)),
                    alpha1: Azimuth::dms(45., 0., 0.000004),
                    p2: (Lon::dms(116., 19., 16.68843), Lat::dms(37., 53., 32.46584)),
                    alpha2: Azimuth::dms(129., 8., 12.326010),
                    s: CurveLength::new(10_002_499.9999, M),
                },
            },
            // Line 7
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::new(0. * DEG), Lat::dms(37., 0., 0.)),
                    alpha1: Azimuth::dms(195., 0., 0.),
                    p2: (Lon::dms(-2., 37., 39.52918), Lat::dms(28., 15., 36.69535)),
                    alpha2: Azimuth::dms(193., 34., 43.74060),
                    s: CurveLength::new(1_000_000.0, M),
                },
            },
        ]
    }

    // From Rapp - Geometric Geodesy Table 1.3
    pub fn antipodal_lines() -> Vec<LineData> {
        vec![
            // Line A
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::new(0.0 * DEG), Lat::dms(41., 41., 45.88)),
                    alpha1: Azimuth::dms(179., 58., 49.16255),
                    p2: (Lon::dms(179., 59., 59.99985), Lat::dms(-41., 41., 46.20)),
                    alpha2: Azimuth::dms(0., 1., 10.8376),
                    s: CurveLength::new(20_004_566.7228, M),
                },
            },
            // Line B
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::new(0. * DEG), Lat::new(0. * DEG)),
                    alpha1: Azimuth::dms(29., 59., 59.9999),
                    p2: (Lon::new(180. * DEG), Lat::new(0. * DEG)),
                    alpha2: Azimuth::new(150. * DEG),
                    s: CurveLength::new(19_996_147.4168, M),
                },
            },
            // Line C
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::new(0. * DEG), Lat::new(30. * DEG)),
                    alpha1: Azimuth::dms(39., 24., 51.8058),
                    p2: (Lon::new(180. * DEG), Lat::new(-30. * DEG)),
                    alpha2: Azimuth::dms(140., 35., 8.1942),
                    s: CurveLength::new(19_994_364.6069, M),
                },
            },
            // Line D
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::new(0. * DEG), Lat::new(60. * DEG)),
                    alpha1: Azimuth::dms(2., 11., 51.0700),
                    p2: (Lon::dms(17., 58., 53.03674), Lat::dms(-59., 59., 0.)),
                    alpha2: Azimuth::dms(150., 49., 6.8680),
                    s: CurveLength::new(20_000_433.9629, M),
                },
            },
            // Line E
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::new(0. * DEG), Lat::new(30. * DEG)),
                    alpha1: Azimuth::dms(16., 2., 28.3389),
                    p2: (Lon::dms(179., 56., 41.64754), Lat::dms(-29., 50., 0.)),
                    alpha2: Azimuth::dms(163., 59., 10.3369),
                    s: CurveLength::new(19_983_420.1536, M),
                },
            },
            // Line F
            LineData {
                ellipsoid: consts::INTL,
                geodesic: Geodesic {
                    p1: (Lon::new(0. * DEG), Lat::new(30. * DEG)),
                    alpha1: Azimuth::dms(18., 38., 12.5568),
                    p2: (Lon::dms(179., 58., 3.57082), Lat::dms(-29., 55., 0.)),
                    alpha2: Azimuth::dms(161., 22., 45.4373),
                    s: CurveLength::new(19_992_241.7634, M),
                },
            },
        ]
    }
}
