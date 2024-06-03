use num::complex::ComplexFloat;

use crate::cs::geodetic::{Lat, Lon};
use crate::cs::Azimuth;

#[derive(Copy, Clone, Debug, PartialEq, Default)]
pub struct Geodesic {
    pub p1: (Lon, Lat),
    pub alpha1: Azimuth,
    pub p2: (Lon, Lat),
    pub alpha2: Azimuth,
    pub s: f64,
}

pub trait GeodesicSolver {
    /// Compute the coordinates and forward azimuth of the point at `s12` meters away
    /// from `p1` following the geodesic in the azimuth `alpha1` at `p1`.
    ///
    /// The algorithm is taken from Karney - Algorithms for geodesics.
    ///
    /// # Parameters
    ///
    /// - ellipsoid: the ellipsoid on which the geodesic is solved
    /// - p1: the **normalized geodetic** coordinates of the starting point
    /// - alpha1: the azimuth **in radians** of the geodesic at `p1`,
    /// - s12: the distance **in meters** along the geodesic from `p1` of the returned point coordinates.
    fn solve_direct(&self, p1: (Lon, Lat), alpha1: Azimuth, s12: f64) -> Geodesic;

    fn solve_inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Geodesic;
}

pub mod vincenty;
pub mod rapp;
pub mod karney;

#[cfg(test)]
mod tests {

    use crate::cs::geodetic::{Lat, Lon};
    use crate::cs::Azimuth;
    use crate::geodesy::geodesics::Geodesic;
    use crate::quantity::angle::dms;
    use crate::quantity::angle::units::DEG;

    // From Rapp - Geometric Geodesy 1.71 Standard Test Lines
    pub fn standard_lines() -> Vec<Geodesic> {
        vec![
            // Line 1
            Geodesic {
                p1: ( Lon::new(0.), Lat::new(dms(37., 19., 54.95367))),
                alpha1: Azimuth::new(dms(95., 27., 59.630888)),
                p2: ( Lon::new(dms(41., 28., 35.50729)), Lat::new(dms(26., 7., 42.83946))),
                alpha2: Azimuth::new(dms(118., 5., 58.961608)),
                s: 4_085_966.7026,
            },
            // Line 2
            Geodesic {
                p1: (
                    Lon::new(0.),
                    Lat::new(dms(35., 16., 11.24862))),
                    alpha1: Azimuth::new(dms(15., 44., 23.748498)),
                p2: (
                    Lon::new(dms(137., 47., 28.31435)),
                    Lat::new(dms(67., 22., 14.77638))),
                    alpha2: Azimuth::new(dms(144., 55., 39.921473)),
                s: 8_084_823.8383,
            },
            // Line 3
            Geodesic {
                p1: (
                    Lon::new(0.),
                    Lat::new(dms(1., 0., 0.))),
                    alpha1: Azimuth::new(dms(88., 59., 59.998970)),
                p2: (
                    Lon::new(dms(179., 17., 48.02997)),
                    Lat::new(dms(-0., 59., 53.83076))),
                    alpha2: Azimuth::new(dms(91., 0., 6.118357)),
                s: 19_959_999.9998,
            },
            // Line 4
            Geodesic {
                p1: (
                    Lon::new(0.),
                    Lat::new(dms(1., 0., 0.))),
                    alpha1: Azimuth::new(dms(4., 59., 59.999953)),
                p2: (
                    Lon::new(dms(179., 46., 17.84244)),
                    Lat::new(dms(1., 1., 15.18952))),
                    alpha2: Azimuth::new(dms(174., 59., 59.884804)),
                s: 19_780_006.5588,
            },
            // Line 5
            Geodesic {
                p1: (
                    Lon::new(0.),
                    Lat::new(dms(41., 41., 45.88000))),
                    alpha1: Azimuth::new(dms(52., 40., 39.390667)),
                p2: (
                    Lon::new(dms(0., 0., 0.56000)),
                    Lat::new(dms(41., 41., 46.20000))),
                    alpha2: Azimuth::new(dms(52., 40., 39.763168)),
                s: 16.2839751,
            },
            // Line 6
            Geodesic {
                p1: (
                    Lon::new(0.),
                    Lat::new(dms(30., 0., 0.))),
                    alpha1: Azimuth::new(dms(45., 0., 0.000004)),
                p2: (
                    Lon::new(dms(116., 19., 16.68843)),
                    Lat::new(dms(37., 53., 32.46584))),
                    alpha2: Azimuth::new(dms(129., 8., 12.326010)),
                s: 10_002_499.9999,
            },
            // Line 7
            Geodesic {
                p1: (
                    Lon::new(0.),
                    Lat::new(dms(37., 0., 0.))),
                    alpha1: Azimuth::new(dms(195., 0., 0.)),
                p2: (
                    Lon::new(dms(-2., 37., 39.52918)),
                    Lat::new(dms(28., 15., 36.69535))),
                    alpha2: Azimuth::new(dms(193., 34., 43.74060)),
                s: 1_000_000.0,
            },
        ]
    }

    // From Rapp - Geometric Geodesy Table 1.3
    pub fn antipodal_lines() -> Vec<Geodesic> {
        vec![
            // Line A
            Geodesic {
                p1: (
                    Lon::new(0.0),
                    Lat::new(dms(41., 41., 45.88))),
                    alpha1: Azimuth::new(dms(179., 58., 49.16255)),
                p2: (
                    Lon::new(dms(179., 59., 59.99985)),
                    Lat::new(dms(-41., 41., 46.20))),
                    alpha2: Azimuth::new(dms(0., 1., 10.8376)),
                s: 20_004_566.7228,
            },
            // Line B
            Geodesic {
                p1: (
                    Lon::new(0.),
                    Lat::new(0. * DEG)),
                    alpha1: Azimuth::new(dms(29., 59., 59.9999)),
                p2: (
                    Lon::new(180. * DEG),
                    Lat::new(0. * DEG)),
                    alpha2: Azimuth::new(150. * DEG),
                s: 19_996_147.4168,
            },
            // Line C
            Geodesic {
                p1: (
                    Lon::new(0.),
                    Lat::new(30. * DEG)),
                    alpha1: Azimuth::new(dms(39., 24., 51.8058)),
                p2: (
                    Lon::new(180. * DEG),
                    Lat::new(-30. * DEG)),
                    alpha2: Azimuth::new(dms(140., 35., 8.1942)),
                s: 19_994_364.6069,
            },
            // Line D
            Geodesic {
                p1: (
                    Lon::new(0.),
                    Lat::new(60. * DEG)),
                    alpha1: Azimuth::new(dms(29., 11., 51.0700)),
                p2: (
                    Lon::new(dms(179., 58., 53.03674)),
                    Lat::new(dms(-59., 59., 0.))),
                    alpha2: Azimuth::new(dms(150., 49., 6.8680)),
                s: 20_000_433.9629,
            },
            // Line E
            Geodesic {
                p1: (
                    Lon::new(0.),
                    Lat::new(30. * DEG)),
                    alpha1: Azimuth::new(dms(16., 2., 28.3389)),
                p2: (
                    Lon::new(dms(179., 56., 41.64754)),
                    Lat::new(dms(-29., 50., 0.))),
                    alpha2: Azimuth::new(dms(163., 59., 10.3369)),
                s: 19_983_420.1536,
            },
            // Line F
            Geodesic {
                p1: (
                    Lon::new(0.),
                    Lat::new(30. * DEG)),
                    alpha1: Azimuth::new(dms(18., 38., 12.5568)),
                p2: (
                    Lon::new(dms(179., 58., 3.57082)),
                    Lat::new(dms(-29., 55., 0.))),
                    alpha2: Azimuth::new(dms(161., 22., 45.4373)),
                s: 19_992_241.7634,
            },
        ]
    }
}
