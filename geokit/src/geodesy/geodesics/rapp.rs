#![allow(non_snake_case)]

use crate::cs::azimuth::Azimuth;
use crate::cs::geodetic::{Lat, Lon};
use crate::geodesy::geodesics::{Geodesic, GeodesicSolver};
use crate::geodesy::Ellipsoid;
use crate::math::fp::PI;
use crate::math::polynomial::Polynomial;
use crate::quantities::length::Length;
use crate::units::angle::RAD;
use crate::units::length::M;

pub struct RappIterativeGeodisicSolver<'e> {
    ellipsoid: &'e Ellipsoid,
    axs: [Polynomial<4>; 4],
}

impl<'e> RappIterativeGeodisicSolver<'e> {
    const MAX_ITERATIONS: usize = 30;

    pub fn new(ellipsoid: &'e Ellipsoid) -> Self {
        // From eq 1.56 to eval at cos_alpha_sq
        let f = ellipsoid.f();
        let axs = [
            // A0
            Polynomial::new([
                1.,                                    // cos^0 alpha
                -(f / 4. * (1. + f * (1. + f))),       // cos^2 alpha
                3. * f * f / 16. * (1. + 9. * f / 4.), // cos^4 alpha
                -25. * f * f * f / 128.,               // cos^6 alpha
            ]),
            // A2
            Polynomial::new([
                0.,                               // cos^0 alpha
                f / 4. * (1. + f * (1. + f)),     // cos^2 alpha
                -f * f / 4. * (1. + 9. * f / 4.), // cos^4 alpha
                75. * f * f * f / 256.,           // cos^6 alpha
            ]),
            // A4
            Polynomial::new([
                0.,
                0.,
                f * f / 32. * (1. + 9. * f / 4.),
                -15. * f * f * f / 256.,
            ]),
            // A6
            Polynomial::new([0., 0., 0., 5. * f * f * f / 768.]),
        ];
        Self { ellipsoid, axs }
    }

    pub fn direct(
        &self,
        p1: (Lon, Lat),
        alpha1: Azimuth,
        s12: Length,
    ) -> Result<Geodesic, &'static str> {
        let (lon1, lat1) = p1;
        let beta1 = self.ellipsoid.reduced_latitude(lat1);
        let (sin_beta1, cos_beta1) = beta1.sin_cos();
        let (sin_alpha1, cos_alpha1) = alpha1.sin_cos();
        let sin_alpha = sin_alpha1 * cos_beta1;
        let cos_alpha_sq = 1. - sin_alpha.powi(2);

        let u_sq = self.ellipsoid.e_prime_sq() * cos_alpha_sq;
        let b0 =
            1. + u_sq * (1. / 4. + u_sq * (-3. / 64. + u_sq * (5. / 256. + u_sq * -175. / 16384.)));
        let b2 = u_sq * (-1. / 4. + u_sq * (1. / 16. + u_sq * (-15. / 512. + u_sq * 35. / 2048.)));
        let b4 = u_sq * u_sq * (-1. / 128. + u_sq * (3. / 512. + u_sq * -35. / 8192.));
        let b6 = u_sq * u_sq * u_sq * (-1. / 1536. + u_sq * 5. / 6144.);

        let sigma_0 = s12.m() / (self.ellipsoid.b() * b0).m();

        let mut sigma = sigma_0;
        let (mut sin_sigma, mut cos_sigma);
        let mut sin_2_sigma;
        let mut sin_3_sigma;
        let mut sigma_1;
        let mut cos_2_sigma_m;
        let mut cos_4_sigma_m;
        let mut cos_6_sigma_m;
        let mut iter_count = 0;
        loop {
            iter_count += 1;
            (sin_sigma, cos_sigma) = sigma.sin_cos();
            sin_2_sigma = 2. * sin_sigma * cos_sigma;
            sin_3_sigma = 3. * sin_sigma * cos_sigma.powi(2) - sin_sigma.powi(3);

            // Eq (1.21)
            sigma_1 = beta1.tan().atan2(cos_alpha1);

            let two_sigma_m = 2. * sigma_1 + sigma;
            cos_2_sigma_m = two_sigma_m.cos();
            cos_4_sigma_m = 2. * cos_2_sigma_m.powi(2) - 1.;
            cos_6_sigma_m =
                cos_2_sigma_m.powi(3) - 3. * cos_2_sigma_m * (1. - cos_2_sigma_m.powi(2));

            let next_sigma = sigma_0
                - b2 / b0 * sin_sigma * cos_2_sigma_m
                - b4 / b0 * sin_2_sigma * cos_4_sigma_m
                - b6 / b0 * sin_3_sigma * cos_6_sigma_m;
            let sigma_old = sigma;
            sigma = next_sigma;
            if (sigma - sigma_old).abs() < 0.5e-14 {
                break;
            }
            if iter_count == Self::MAX_ITERATIONS {
                return Err("Direct geodesic computation is not converging!");
            }
        }

        let alpha = sin_alpha.asin();
        // From Eq (1.22)
        let beta2 = ((sigma_1 + sigma).sin() * alpha.cos()).asin() * RAD;
        let (sin_beta2, cos_beta2) = beta2.sin_cos();
        // geodetic latitude from reduced latitude
        let lat2 = beta2.tan().atan2(1. - self.ellipsoid.f());

        // Eq (1.76)
        let sin_lambda = sin_sigma * sin_alpha1 / cos_beta2;
        // From Eq (1.57)
        let cos_lambda = (cos_sigma - sin_beta1 * sin_beta2) / (cos_beta1 * cos_beta2);
        let lambda = sin_lambda.atan2(cos_lambda);

        // Eq (1.56)
        let f = self.ellipsoid.f();
        let big_as = self.axs.map(|p| p.fast_eval_at(cos_alpha_sq));
        let lambda_minus_L = f
            * sin_alpha
            * (big_as[0] * sigma
                + big_as[1] * sin_sigma * cos_2_sigma_m
                + big_as[2] * sin_2_sigma * cos_4_sigma_m
                + big_as[3] * sin_3_sigma * cos_6_sigma_m);

        let L = lambda - lambda_minus_L;

        let alpha21 = if s12.m() < 1000.0 {
            let sin_lambda_2_sq = (1. - cos_lambda) / 2.;
            // Eq (1.73)
            (sin_lambda * cos_beta1)
                .atan2((beta2 - beta1).sin() - 2. * cos_beta1 * sin_beta2 * sin_lambda_2_sq)
        } else {
            // Eq (1.71)
            (sin_lambda * cos_beta1)
                .atan2(sin_beta2 * cos_beta1 * cos_lambda - sin_beta1 * cos_beta2)
        };

        Ok(Geodesic {
            p1,
            alpha1,
            p2: (lon1 + L * RAD, Lat::new(lat2 * RAD)),
            alpha2: Azimuth::new(alpha21 * RAD),
            s: s12,
        })
    }

    pub fn inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Result<Geodesic, &'static str> {
        // Init
        let f = self.ellipsoid.f();
        let (lon1, lat1) = p1;
        let (lon2, lat2) = p2;
        let beta1 = self.ellipsoid.reduced_latitude(lat1);
        let beta2 = self.ellipsoid.reduced_latitude(lat2);
        let (sin_beta1, cos_beta1) = beta1.sin_cos();
        let (sin_beta2, cos_beta2) = beta2.sin_cos();

        let L = (lon2 - lon1).rad();
        let mut lambda = L;
        let (mut sin_lambda, mut cos_lambda);
        let mut sigma;
        let mut sin_sigma;
        let mut sin_2_sigma;
        let mut sin_3_sigma;
        let mut sin_4_sigma;
        let mut sin_alpha;
        let mut cos_alpha_sq;
        let mut cos_2_sigma_m = 0.0;
        let mut cos_4_sigma_m = 0.0;
        let mut cos_6_sigma_m = 0.0;
        let mut cos_8_sigma_m = 0.0;
        let mut iter_count = 0;
        loop {
            iter_count += 1;
            (sin_lambda, cos_lambda) = lambda.sin_cos();
            // Eq (1.57)
            let cos_sigma = sin_beta1 * sin_beta2 + cos_beta1 * cos_beta2 * cos_lambda;
            // Eq (1.58)
            sin_sigma = ((sin_lambda * cos_beta2).powi(2)
                + (sin_beta2 * cos_beta1 - sin_beta1 * cos_beta2 * cos_lambda).powi(2))
            .sqrt();
            sigma = sin_sigma.atan2(cos_sigma);
            sin_2_sigma = 2. * sin_sigma * cos_sigma;
            sin_3_sigma = 3. * sin_sigma * cos_sigma.powi(2) - sin_sigma.powi(3);
            sin_4_sigma = 4. * sin_sigma * cos_sigma.powi(3) - 4. * sin_sigma.powi(3) * cos_sigma;

            // Eq (1.61)
            sin_alpha = cos_beta1 * cos_beta2 * sin_lambda / sin_sigma;
            cos_alpha_sq = 1. - sin_alpha.powi(2);

            // If the geodesic (p1, p2) is an equatorial line Eq (18) becomes indeterminate but
            // big_as becomes [1., 0., 0., ...]
            let lambda_minus_L = if cos_alpha_sq == 0.0 {
                f * sin_alpha * sigma
            } else {
                // Eq (1.66)
                cos_2_sigma_m = cos_sigma - 2. * sin_beta1 * sin_beta2 / cos_alpha_sq;
                cos_4_sigma_m = 2. * cos_2_sigma_m.powi(2) - 1.;
                cos_6_sigma_m =
                    cos_2_sigma_m.powi(3) - 3. * cos_2_sigma_m * (1. - cos_2_sigma_m.powi(2));
                cos_8_sigma_m = cos_2_sigma_m.powi(4)
                    - 6. * cos_2_sigma_m.powi(2) * sin_2_sigma.powi(2)
                    + sin_2_sigma.powi(4);

                // Eq (1.56)
                let big_as = self.axs.map(|p| p.eval_at(cos_alpha_sq));
                f * sin_alpha
                    * (big_as[0] * sigma
                        + big_as[1] * sin_sigma * cos_2_sigma_m
                        + big_as[2] * sin_2_sigma * cos_4_sigma_m
                        + big_as[3] * sin_3_sigma * cos_6_sigma_m)
            };
            let lambda_prev = lambda;
            lambda = L + lambda_minus_L;
            if (lambda - lambda_prev).abs() < 0.5e-14 {
                break;
            }
            if iter_count == Self::MAX_ITERATIONS {
                return Err("Inverse geodesic computation is not converging!");
            }
        }

        // If the geodesic (p1, p2) is an equatorial line b0 = 1 and bj = 0 j > 1, and s expression
        // is simplified.
        let s = if cos_alpha_sq == 0.0 {
            Length::new((self.ellipsoid.b() * sigma).m(), M)
        } else {
            let u_sq = self.ellipsoid.e_prime_sq() * cos_alpha_sq;
            let b0 = 1.
                + u_sq
                    * (1. / 4. + u_sq * (-3. / 64. + u_sq * (5. / 256. + u_sq * -175. / 16384.)));
            let b2 =
                u_sq * (-1. / 4. + u_sq * (1. / 16. + u_sq * (-15. / 512. + u_sq * 35. / 2048.)));
            let b4 = u_sq * u_sq * (-1. / 128. + u_sq * (3. / 512. + u_sq * -35. / 8192.));
            let b6 = u_sq * u_sq * u_sq * (-1. / 1536. + u_sq * 5. / 6144.);
            let b8 = u_sq * u_sq * u_sq * u_sq * (-5. / 65536.);

            let s_m = (self.ellipsoid.b()
                * (b0 * sigma
                    + b2 * sin_sigma * cos_2_sigma_m
                    + b4 * sin_2_sigma * cos_4_sigma_m
                    + b6 * sin_3_sigma * cos_6_sigma_m
                    + b8 * sin_4_sigma * cos_8_sigma_m))
                .m();
            Length::new(s_m, M)
        };
        let sin_alpha1 = sin_alpha / cos_beta1;
        let sin_alpha2 = sin_alpha / cos_beta2;
        let sin_lambda_2_sq = (1. - cos_lambda) / 2.;

        // If the geodesic (p1, p2) is a meridional line, 2 cases:
        let (alpha1, alpha2) = if lambda == 0. {
            // lambda = 0 -> no pole crossing
            if lat1 < lat2 {
                (0., 0.)
            } else {
                (PI, PI)
            }
        } else if lambda == PI {
            // lambda = PI -> pole crossing
            if lat1 >= Lat::ZERO {
                // crossing north pole by going north then south
                (0., PI)
            } else {
                // crossing south pole by going south then north
                (PI, 0.)
            }
        } else {
            let (tan_alpha12, tan_alpha21) = if s.m() < 10_000.0 {
                (
                    // Eq (1.72)
                    lambda.sin() * cos_beta2
                        / ((beta2 - beta1).sin() + 2. * sin_beta1 * cos_beta2 * sin_lambda_2_sq),
                    // Eq (1.73)
                    lambda.sin() * cos_beta1
                        / ((beta2 - beta1).sin() - 2. * cos_beta1 * sin_beta2 * sin_lambda_2_sq),
                )
            } else {
                (
                    // Eq (1.70)
                    lambda.sin() * cos_beta2
                        / (sin_beta2 * cos_beta1 - lambda.cos() * sin_beta1 * cos_beta2),
                    // Eq (1.71)
                    lambda.sin() * cos_beta1
                        / (sin_beta2 * cos_beta1 * lambda.cos() - sin_beta1 * cos_beta2),
                )
            };
            (
                sin_alpha1.atan2(sin_alpha1 / tan_alpha12),
                sin_alpha2.atan2(sin_alpha2 / tan_alpha21),
            )
        };

        Ok(Geodesic {
            p1: (lon1, lat1),
            alpha1: Azimuth::new(alpha1 * RAD),
            p2: (lon2, lat2),
            alpha2: Azimuth::new(alpha2 * RAD),
            s,
        })
    }
}

impl<'e> GeodesicSolver for RappIterativeGeodisicSolver<'e> {
    fn solve_direct(
        &self,
        p1: (Lon, Lat),
        alpha1: Azimuth,
        s12: Length,
    ) -> Result<Geodesic, &'static str> {
        self.direct(p1, alpha1, s12)
    }

    fn solve_inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Result<Geodesic, &'static str> {
        self.inverse(p1, p2)
    }
}

#[cfg(test)]
mod tests {
    use crate::geodesy::geodesics::rapp::RappIterativeGeodisicSolver;
    use crate::geodesy::geodesics::tests::{
        antipodal_lines, equatorial_lines, geographiclib_lines, meridional_lines, standard_lines,
        LineData,
    };
    use crate::geodesy::geodesics::{
        check_direct, check_inverse, DirectErrors, GeodesicSolver, InverseErrors,
    };

    fn test_on(tset: LineData, err_direct: &DirectErrors, err_inverse: &InverseErrors) {
        let solver = RappIterativeGeodisicSolver::new(&tset.ellipsoid);
        for tcase in tset.testcases.into_iter() {
            let direct = solver
                .solve_direct(tcase.p1, tcase.alpha1, tcase.s)
                .unwrap();
            check_direct(&direct, &tcase, err_direct);

            let inverse = solver.solve_inverse(tcase.p1, tcase.p2).unwrap();
            check_inverse(&inverse, &tcase, err_inverse);
        }
    }

    #[test]
    fn on_geographiclib_lines() {
        let tset = geographiclib_lines();
        test_on(tset, &DirectErrors::default(), &InverseErrors::default());
    }

    #[test]
    fn on_standard_lines() {
        let tset = standard_lines();
        test_on(tset, &DirectErrors::default(), &InverseErrors::default());
    }

    #[test]
    fn on_equatorial_lines() {
        let tset = equatorial_lines();
        test_on(tset, &DirectErrors::default(), &InverseErrors::default());
    }

    #[test]
    fn on_meridional_lines() {
        let tset = meridional_lines();
        test_on(tset, &DirectErrors::default(), &InverseErrors::default());
    }

    #[ignore = "Known to fail"]
    #[test]
    fn on_antipodal_lines() {
        let tset = antipodal_lines();
        test_on(tset, &DirectErrors::default(), &InverseErrors::default());
    }
}
