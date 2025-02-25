use crate::cs::azimuth::Azimuth;
use crate::cs::geodetic::{Lat, Lon};
use crate::quantities::angle::Angle;
use crate::quantities::length::Length;
use approx::AbsDiffEq;
use derive_more::derive::Display;

/// Geodesic segment data.
#[derive(Copy, Clone, Debug, PartialEq, Display)]
#[display("From: (lon: {}, lat: {}, az: {})\nTo:   (lon: {}, lat: {}, az: {})\nL:   {}\ns:   {}", p1.0, p1.0, alpha1, p2.0, p2.1, alpha2, p2.0 - p1.0, s)]
pub struct Geodesic {
    pub p1: (Lon, Lat),
    pub alpha1: Azimuth,
    pub p2: (Lon, Lat),
    pub alpha2: Azimuth,
    pub s: Length,
}

impl Default for Geodesic {
    fn default() -> Self {
        Geodesic {
            p1: (Lon::ZERO, Lat::ZERO),
            alpha1: Azimuth::NORTH,
            p2: (Lon::ZERO, Lat::ZERO),
            alpha2: Azimuth::NORTH,
            s: Length::ZERO,
        }
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
    fn solve_direct(
        &self,
        p1: (Lon, Lat),
        alpha1: Azimuth,
        s12: Length,
    ) -> Result<Geodesic, &'static str>;

    /// Compute the azimuths and length in meters of the geodesic from `p1` to `p2`.
    ///
    /// Some implementations may not converge when inputs are *near* antipodal points.
    fn solve_inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Result<Geodesic, &'static str>;
}

/// Maximum absolute difference in direct geodesic computation results.
pub struct DirectErrors {
    /// Maximum absolute difference in longitude.
    lon: Angle,
    /// Maximum absolute difference in latitude.
    lat: Angle,
    /// Maximum absolute difference in azimuth
    alpha: Angle,
}

impl Default for DirectErrors {
    fn default() -> Self {
        DirectErrors {
            lon: Angle::default_epsilon(),
            lat: Angle::default_epsilon(),
            alpha: Angle::default_epsilon(),
        }
    }
}

/// Maximum absolute difference in inverse geodesic computation results.
pub struct InverseErrors {
    /// Maximum absolute difference in starting azimuth.
    alpha1: Angle,
    /// Maximum absolute difference in ending azimuth.
    alpha2: Angle,
    /// Maximum absolute difference in geodesic distance.
    s: Length,
}

impl Default for InverseErrors {
    fn default() -> Self {
        InverseErrors {
            alpha1: Angle::default_epsilon(),
            alpha2: Angle::default_epsilon(),
            s: Length::default_epsilon(),
        }
    }
}

/// Checks whether a direct geodesic computation is within the error bounds of the expected
/// result.
pub fn check_direct(computed: &Geodesic, expected: &Geodesic, err: &DirectErrors) {
    let has_lon_error = !computed.p2.0.abs_diff_eq(&expected.p2.0, err.lon);
    let has_lat_error = !computed.p2.1.abs_diff_eq(&expected.p2.1, err.lat);
    let has_az_error = !computed.alpha2.abs_diff_eq(&expected.alpha2, err.alpha);

    if has_lon_error || has_lat_error || has_az_error {
        print_test_case(computed, expected);
        if has_lon_error {
            println!(
                "Longitude error: | {} - {} | = {:e} > {:e}",
                computed.p2.0,
                expected.p2.0,
                (computed.p2.0 - expected.p2.0).abs().deg(),
                err.lon.deg()
            );
        } else {
            println!(
                "Longitude ok:      {} = {} +/- {:e}",
                computed.p2.0,
                expected.p2.0,
                err.lon.deg()
            );
        }
        if has_lat_error {
            println!(
                "Latitude error: | {} - {} | = {:e} > {:e}",
                computed.p2.1,
                expected.p2.1,
                (computed.p2.1 - expected.p2.1).abs().deg(),
                err.lat.deg()
            );
        } else {
            println!(
                "Latitude ok:      {} = {} +/- {:e}",
                computed.p2.1,
                expected.p2.1,
                err.lat.deg()
            );
        }
        if has_az_error {
            println!(
                "Azimuth error: | {} - {} | = {:e} > {:e}",
                computed.alpha2,
                expected.alpha2,
                (computed.alpha2 - expected.alpha2).abs().deg(),
                err.alpha.deg()
            );
        } else {
            println!(
                "Azimuth ok:      {} = {} +/- {:e}",
                computed.alpha2,
                expected.alpha2,
                err.alpha.deg()
            );
        }
        assert!(false);
    }
}

/// Check whether an inverse geodesic computation is within the error bounds of the expected
/// result.
pub fn check_inverse(computed: &Geodesic, expected: &Geodesic, err: &InverseErrors) {
    let has_alpha1_err = !computed.alpha1.abs_diff_eq(&expected.alpha1, err.alpha1);
    let has_alpha2_err = !computed.alpha2.abs_diff_eq(&expected.alpha2, err.alpha2);
    let has_s_err = !computed.s.abs_diff_eq(&expected.s, err.s);

    if has_alpha1_err || has_alpha2_err || has_s_err {
        print_test_case(computed, expected);
        if has_alpha1_err {
            println!(
                "Alpha1 error: | {} - {} | = {:e} > {:e}",
                computed.alpha1,
                expected.alpha1,
                (computed.alpha1 - expected.alpha1).abs().deg(),
                err.alpha1.deg()
            );
        } else {
            println!(
                "Alpha1 ok:     {} = {} +/- {:e}",
                computed.alpha1,
                expected.alpha1,
                err.alpha1.deg(),
            );
        }
        if has_alpha2_err {
            println!(
                "Alpha2 error: | {} - {} | = {:e} > {:e}",
                computed.alpha2,
                expected.alpha2,
                (computed.alpha2 - expected.alpha2).abs().deg(),
                err.alpha2.deg()
            );
        } else {
            println!(
                "Alpha2 ok:     {} = {} +/- {:e}",
                computed.alpha2,
                expected.alpha2,
                err.alpha2.deg(),
            );
        }
        if has_s_err {
            println!(
                "S error: | {} - {} | = {:e} m > {:e} m",
                computed.s,
                expected.s,
                (computed.s - expected.s).abs().m(),
                err.s.m(),
            );
        } else {
            println!(
                "S ok      {} = {} +/- {:e} m",
                computed.s,
                expected.s,
                err.s.m()
            );
        }
        assert!(false)
    }
}

fn print_test_case(computed: &Geodesic, expected: &Geodesic) {
    println!("------------------------ FAILURE --------------------------------");
    println!("Computed: ");
    println!("{}", computed);
    println!("Expected: ");
    println!("{}", expected);
    println!("");
}

pub mod karney;
pub mod rapp;
pub mod vincenty;

#[cfg(test)]
mod tests {
    use crate::cs::azimuth::Azimuth;
    use crate::cs::geodetic::{Lat, Lon};
    use crate::geodesy::ellipsoid::consts::WGS84;
    use crate::geodesy::ellipsoid::{self, consts};
    use crate::geodesy::geodesics::Geodesic;
    use crate::geodesy::Ellipsoid;
    use crate::math::fp::Float;
    use crate::quantities::length::Length;
    use crate::units::angle::DEG;
    use crate::units::length::M;

    pub struct LineData {
        pub ellipsoid: Ellipsoid,
        pub testcases: Vec<Geodesic>,
    }

    /// From Rapp - Geometric Geodesy 1.71 Standard Test Lines Tables 1.1 and 1.2
    pub fn standard_lines() -> LineData {
        LineData {
            ellipsoid: consts::INTL,
            testcases: vec![
                // Line 1
                Geodesic {
                    p1: (Lon::ZERO, Lat::dms(37., 19., 54.95367)),
                    alpha1: Azimuth::dms(95., 27., 59.630888),
                    p2: (Lon::dms(41., 28., 35.50729), Lat::dms(26., 7., 42.83946)),
                    alpha2: Azimuth::dms(118., 5., 58.961608),
                    s: Length::new(4_085_966.7026, M),
                },
                // Line 2
                Geodesic {
                    p1: (Lon::ZERO, Lat::dms(35., 16., 11.24862)),
                    alpha1: Azimuth::dms(15., 44., 23.748498),
                    p2: (Lon::dms(137., 47., 28.31435), Lat::dms(67., 22., 14.77638)),
                    alpha2: Azimuth::dms(144., 55., 39.921473),
                    s: Length::new(8_084_823.8383, M),
                },
                // Line 3
                Geodesic {
                    p1: (Lon::ZERO, Lat::dms(1., 0., 0.)),
                    alpha1: Azimuth::dms(88., 59., 59.998970),
                    p2: (Lon::dms(179., 17., 48.02997), Lat::dms(-0., 59., 53.83076)),
                    alpha2: Azimuth::dms(91., 0., 6.118357),
                    s: Length::new(19_959_999.9998, M),
                },
                // Line 4
                Geodesic {
                    p1: (Lon::ZERO, Lat::dms(1., 0., 0.)),
                    alpha1: Azimuth::dms(4., 59., 59.999953),
                    p2: (Lon::dms(179., 46., 17.84244), Lat::dms(1., 1., 15.18952)),
                    alpha2: Azimuth::dms(174., 59., 59.884804),
                    s: Length::new(19_780_006.5588, M),
                },
                Geodesic {
                    p1: (Lon::ZERO, Lat::dms(41., 41., 45.88000)),
                    alpha1: Azimuth::dms(52., 40., 39.390667),
                    p2: (Lon::dms(0., 0., 0.56000), Lat::dms(41., 41., 46.20000)),
                    alpha2: Azimuth::dms(52., 40., 39.763168),
                    s: Length::new(16.2839751, M),
                },
                Geodesic {
                    p1: (Lon::ZERO, Lat::dms(30., 0., 0.)),
                    alpha1: Azimuth::dms(45., 0., 0.000004),
                    p2: (Lon::dms(116., 19., 16.68843), Lat::dms(37., 53., 32.46584)),
                    alpha2: Azimuth::dms(129., 8., 12.326010),
                    s: Length::new(10_002_499.9999, M),
                },
                Geodesic {
                    p1: (Lon::ZERO, Lat::dms(37., 0., 0.)),
                    alpha1: Azimuth::dms(195., 0., 0.),
                    p2: (Lon::dms(-2., 37., 39.52918), Lat::dms(28., 15., 36.69535)),
                    alpha2: Azimuth::dms(193., 34., 43.74060),
                    s: Length::new(1_000_000.0, M),
                },
            ],
        }
    }

    pub fn equatorial_lines() -> LineData {
        LineData {
            ellipsoid: WGS84,
            testcases: vec![
                Geodesic {
                    p1: (Lon::ZERO, Lat::ZERO),
                    alpha1: Azimuth::EAST,
                    p2: (Lon::new(0.17966306 * DEG), Lat::ZERO),
                    alpha2: Azimuth::EAST,
                    s: Length::new(20_000.0, M),
                },
                Geodesic {
                    p1: (Lon::new(170.0 * DEG), Lat::ZERO),
                    alpha1: Azimuth::new(90.0 * DEG),
                    p2: (Lon::new(-172.03369432 * DEG), Lat::ZERO),
                    alpha2: Azimuth::EAST,
                    s: Length::new(2_000_000.0, M),
                },
                Geodesic {
                    p1: (Lon::new(-10. * DEG), Lat::ZERO),
                    alpha1: Azimuth::EAST,
                    p2: (Lon::new(10. * DEG), Lat::ZERO),
                    alpha2: Azimuth::EAST,
                    s: Length::new(2_226_389.816, M),
                },
                Geodesic {
                    p1: (Lon::new(10. * DEG), Lat::ZERO),
                    alpha1: Azimuth::WEST,
                    p2: (Lon::new(-10. * DEG), Lat::ZERO),
                    alpha2: Azimuth::WEST,
                    s: Length::new(2_226_389.816, M),
                },
                Geodesic {
                    p1: (Lon::new(170. * DEG), Lat::ZERO),
                    alpha1: Azimuth::EAST,
                    p2: (Lon::new(-170. * DEG), Lat::ZERO),
                    alpha2: Azimuth::EAST,
                    s: Length::new(2_226_389.816, M),
                },
            ],
        }
    }

    pub fn meridional_lines() -> LineData {
        LineData {
            ellipsoid: WGS84,
            testcases: vec![
                Geodesic {
                    p1: (Lon::ZERO, Lat::new(-10. * DEG)),
                    alpha1: Azimuth::NORTH,
                    p2: (Lon::ZERO, Lat::new(8.08583903 * DEG)),
                    alpha2: Azimuth::NORTH,
                    s: 2_000_000.0 * M,
                },
                Geodesic {
                    p1: (Lon::ZERO, Lat::new(80. * DEG)),
                    alpha1: Azimuth::NORTH,
                    p2: (Lon::MAX, Lat::new(82.09240627 * DEG)),
                    alpha2: Azimuth::SOUTH,
                    s: Length::new(2_000_000.0, M),
                },
                Geodesic {
                    p1: (Lon::ZERO, Lat::new(-10. * DEG)),
                    alpha1: Azimuth::NORTH,
                    p2: (Lon::ZERO, Lat::new(10. * DEG)),
                    alpha2: Azimuth::NORTH,
                    s: 2_211_709.666 * M,
                },
                Geodesic {
                    p1: (Lon::ZERO, Lat::new(10. * DEG)),
                    alpha1: Azimuth::SOUTH,
                    p2: (Lon::ZERO, Lat::new(-10. * DEG)),
                    alpha2: Azimuth::SOUTH,
                    s: 2_211_709.666 * M,
                },
                Geodesic {
                    p1: (Lon::ZERO, Lat::new(80. * DEG)),
                    alpha1: Azimuth::NORTH,
                    p2: (Lon::new(180. * DEG), Lat::new(80. * DEG)),
                    alpha2: Azimuth::SOUTH,
                    s: 2_233_651.715 * M,
                },
            ],
        }
    }

    /// From Rapp - Geometric Geodesy Table 1.3
    pub fn antipodal_lines() -> LineData {
        LineData {
            ellipsoid: consts::INTL,
            testcases: vec![
                // Line A
                Geodesic {
                    p1: (Lon::ZERO, Lat::dms(41., 41., 45.88)),
                    alpha1: Azimuth::dms(179., 58., 49.16255),
                    p2: (Lon::dms(179., 59., 59.44), Lat::dms(-41., 41., 46.20)),
                    alpha2: Azimuth::dms(0., 1., 10.8376),
                    s: Length::new(20_004_566.7228, M),
                },
                Geodesic {
                    p1: (Lon::ZERO, Lat::ZERO),
                    alpha1: Azimuth::dms(29., 59., 59.9999),
                    p2: (Lon::dms(179., 41., 49.78063), Lat::ZERO),
                    alpha2: Azimuth::new(150. * DEG),
                    s: Length::new(19_996_147.4168, M),
                },
                Geodesic {
                    p1: (Lon::ZERO, Lat::new(30. * DEG)),
                    alpha1: Azimuth::dms(39., 24., 51.8058),
                    p2: (Lon::dms(179., 40., 0.), Lat::new(-30. * DEG)),
                    alpha2: Azimuth::dms(140., 35., 8.1942),
                    s: Length::new(19_994_364.6069, M),
                },
                Geodesic {
                    p1: (Lon::ZERO, Lat::new(60. * DEG)),
                    alpha1: Azimuth::dms(29., 11., 51.0700),
                    p2: (Lon::dms(179., 50., 0.), Lat::dms(-59., 59., 0.)),
                    alpha2: Azimuth::dms(150., 49., 6.8680),
                    s: Length::new(20_000_433.9629, M),
                },
                Geodesic {
                    p1: (Lon::ZERO, Lat::new(30. * DEG)),
                    alpha1: Azimuth::dms(16., 2., 28.3389),
                    p2: (Lon::dms(179., 48., 0.), Lat::dms(-29., 50., 0.)),
                    alpha2: Azimuth::dms(163., 59., 10.3369),
                    s: Length::new(19_983_420.1536, M),
                },
                Geodesic {
                    p1: (Lon::ZERO, Lat::new(30. * DEG)),
                    alpha1: Azimuth::dms(18., 38., 12.5568),
                    p2: (Lon::dms(179., 48., 0.), Lat::dms(-29., 55., 0.)),
                    alpha2: Azimuth::dms(161., 22., 45.4373),
                    s: Length::new(19_992_241.7634, M),
                },
            ],
        }
    }

    pub fn geographiclib_lines() -> LineData {
        LineData {
            ellipsoid: ellipsoid::consts::WGS84,
            testcases: GEOGRAPHICLIB_TESTCASES
                .iter()
                .map(|t| Geodesic {
                    p1: (Lon::new(t[1] * DEG), Lat::new(t[0] * DEG)),
                    alpha1: Azimuth::new(t[2] * DEG),
                    p2: (Lon::new(t[4] * DEG), Lat::new(t[3] * DEG)),
                    alpha2: Azimuth::new(t[5] * DEG),
                    s: Length::new(t[6], M),
                })
                .collect::<Vec<_>>(),
        }
    }

    /// [0]: lat1 in degrees
    /// [1]: lon1 in degrees
    /// [2]: az1 in degrees
    /// [3]: lat2 in degrees
    /// [4]: lon2 in degrees
    /// [5]: az2 in degrees
    /// [6]: s12 in meters
    /// [7]: a12
    /// [8]: m12
    /// [9]: M12
    /// [10]: M21
    /// [11]: S12
    type GeographicLibTestCase = [Float; 12];

    /// From <https://github.com/geographiclib/geographiclib/blob/main/tests/geodtest.cpp#L27>
    /// use WGS84 ellipsoid
    const GEOGRAPHICLIB_TESTCASES: [GeographicLibTestCase; 20] = [
        [
            35.60777,
            -139.44815,
            111.098748429560326,
            -11.17491,
            -69.95921,
            129.289270889708762,
            8935244.5604818305,
            80.50729714281974,
            6273170.2055303837,
            0.16606318447386067,
            0.16479116945612937,
            12841384694976.432,
        ],
        [
            55.52454,
            106.05087,
            22.020059880982801,
            77.03196,
            197.18234,
            109.112041110671519,
            4105086.1713924406,
            36.892740690445894,
            3828869.3344387607,
            0.80076349608092607,
            0.80101006984201008,
            61674961290615.615,
        ],
        [
            -21.97856,
            142.59065,
            -32.44456876433189,
            41.84138,
            98.56635,
            -41.84359951440466,
            8394328.894657671,
            75.62930491011522,
            6161154.5773110616,
            0.24816339233950381,
            0.24930251203627892,
            -6637997720646.717,
        ],
        [
            -66.99028,
            112.2363,
            173.73491240878403,
            -12.70631,
            285.90344,
            2.512956620913668,
            11150344.2312080241,
            100.278634181155759,
            6289939.5670446687,
            -0.17199490274700385,
            -0.17722569526345708,
            -121287239862139.744,
        ],
        [
            -17.42761,
            173.34268,
            -159.033557661192928,
            -15.84784,
            5.93557,
            -20.787484651536988,
            16076603.1631180673,
            144.640108810286253,
            3732902.1583877189,
            -0.81273638700070476,
            -0.81299800519154474,
            97825992354058.708,
        ],
        [
            32.84994,
            48.28919,
            150.492927788121982,
            -56.28556,
            202.29132,
            48.113449399816759,
            16727068.9438164461,
            150.565799985466607,
            3147838.1910180939,
            -0.87334918086923126,
            -0.86505036767110637,
            -72445258525585.010,
        ],
        [
            6.96833,
            52.74123,
            92.581585386317712,
            -7.39675,
            206.17291,
            90.721692165923907,
            17102477.2496958388,
            154.147366239113561,
            2772035.6169917581,
            -0.89991282520302447,
            -0.89986892177110739,
            -1311796973197.995,
        ],
        [
            -50.56724,
            -16.30485,
            -105.439679907590164,
            -33.56571,
            -94.97412,
            -47.348547835650331,
            6455670.5118668696,
            58.083719495371259,
            5409150.7979815838,
            0.53053508035997263,
            0.52988722644436602,
            41071447902810.047,
        ],
        [
            -58.93002,
            -8.90775,
            140.965397902500679,
            -8.91104,
            133.13503,
            19.255429433416599,
            11756066.0219864627,
            105.755691241406877,
            6151101.2270708536,
            -0.26548622269867183,
            -0.27068483874510741,
            -86143460552774.735,
        ],
        [
            -68.82867,
            -74.28391,
            93.774347763114881,
            -50.63005,
            -8.36685,
            34.65564085411343,
            3956936.926063544,
            35.572254987389284,
            3708890.9544062657,
            0.81443963736383502,
            0.81420859815358342,
            -41845309450093.787,
        ],
        [
            -10.62672,
            -32.0898,
            -86.426713286747751,
            5.883,
            -134.31681,
            -80.473780971034875,
            11470869.3864563009,
            103.387395634504061,
            6184411.6622659713,
            -0.23138683500430237,
            -0.23155097622286792,
            4198803992123.548,
        ],
        [
            -21.76221,
            166.90563,
            29.319421206936428,
            48.72884,
            213.97627,
            43.508671946410168,
            9098627.3986554915,
            81.963476716121964,
            6299240.9166992283,
            0.13965943368590333,
            0.14152969707656796,
            10024709850277.476,
        ],
        [
            -19.79938,
            -174.47484,
            71.167275780171533,
            -11.99349,
            -154.35109,
            65.589099775199228,
            2319004.8601169389,
            20.896611684802389,
            2267960.8703918325,
            0.93427001867125849,
            0.93424887135032789,
            -3935477535005.785,
        ],
        [
            -11.95887,
            -116.94513,
            92.712619830452549,
            4.57352,
            7.16501,
            78.64960934409585,
            13834722.5801401374,
            124.688684161089762,
            5228093.177931598,
            -0.56879356755666463,
            -0.56918731952397221,
            -9919582785894.853,
        ],
        [
            -87.85331,
            85.66836,
            -65.120313040242748,
            66.48646,
            16.09921,
            -4.888658719272296,
            17286615.3147144645,
            155.58592449699137,
            2635887.4729110181,
            -0.90697975771398578,
            -0.91095608883042767,
            42667211366919.534,
        ],
        [
            1.74708,
            128.32011,
            -101.584843631173858,
            -11.16617,
            11.87109,
            -86.325793296437476,
            12942901.1241347408,
            116.650512484301857,
            5682744.8413270572,
            -0.44857868222697644,
            -0.44824490340007729,
            10763055294345.653,
        ],
        [
            -25.72959,
            -144.90758,
            -153.647468693117198,
            -57.70581,
            -269.17879,
            -48.343983158876487,
            9413446.7452453107,
            84.664533838404295,
            6356176.6898881281,
            0.09492245755254703,
            0.09737058264766572,
            74515122850712.444,
        ],
        [
            -41.22777,
            122.32875,
            14.285113402275739,
            -7.57291,
            130.37946,
            10.805303085187369,
            3812686.035106021,
            34.34330804743883,
            3588703.8812128856,
            0.82605222593217889,
            0.82572158200920196,
            -2456961531057.857,
        ],
        [
            11.01307,
            138.25278,
            79.43682622782374,
            6.62726,
            247.05981,
            103.708090215522657,
            11911190.819018408,
            107.341669954114577,
            6070904.722786735,
            -0.29767608923657404,
            -0.29785143390252321,
            17121631423099.696,
        ],
        [
            -29.47124,
            95.14681,
            -163.779130441688382,
            -27.46601,
            -69.15955,
            -15.909335945554969,
            13487015.8381145492,
            121.294026715742277,
            5481428.9945736388,
            -0.51527225545373252,
            -0.51556587964721788,
            104679964020340.318,
        ],
    ];
}
