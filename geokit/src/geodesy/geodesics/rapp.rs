use crate::cs::geodetic::{Lat, Lon};
use crate::cs::Azimuth;
use crate::geodesy::geodesics::{Geodesic, GeodesicSolver};
use crate::geodesy::Ellipsoid;
use crate::math::polynomial::Polynomial;

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
                1.,                                     // cos^0 alpha
                -(f / 4. * (1. + f * (1. + f))),        // cos^2 alpha
                3. * f * f / 16. * (1. + 9. * f / 4.),  // cos^4 alpha
                -25. * f * f * f / 128.,                // cos^6 alpha
            ]),
            // A2
            Polynomial::new([
                0.,                              // cos^0 alpha
                f / 4. * (1. + f * (1. + f)),    // cos^2 alpha
                -f * f / 4. * (1. + 9. * f / 4.),// cos^4 alpha
                75. * f * f * f / 256.,          // cos^6 alpha
            ]),
            // A4
            Polynomial::new([
                0.,
                0.,
                f * f / 32. * (1. + 9. * f / 4.),
                -15. * f * f * f / 256.,
            ]),
            // A6
            Polynomial::new([
                0.,
                0.,
                0.,
                5. * f * f * f / 768.,
            ]),
        ];
        Self { ellipsoid, axs }
    }

    pub fn direct(&self, p1: (Lon, Lat), alpha1: Azimuth, s12: f64) -> Geodesic {
        let (lon1, lat1) = p1;
        let beta1 = self.ellipsoid.reduced_latitude(lat1.rad());
        let (sin_beta1, cos_beta1) = beta1.sin_cos();
        let (sin_alpha1, cos_alpha1) = alpha1.rad().sin_cos();
        let sin_alpha = sin_alpha1 * cos_beta1;
        let cos_alpha_sq = 1. - sin_alpha.powi(2);

        let u_sq = self.ellipsoid.e_prime_sq() * cos_alpha_sq;
        let b0 = 1. + u_sq * (1. / 4. + u_sq * (-3. / 64. + u_sq * (5. / 256. + u_sq * -175. / 16384.)));
        let b2 = u_sq * (-1. / 4. + u_sq * (1. / 16. + u_sq * (-15. / 512. + u_sq * 35. / 2048.)));
        let b4 = u_sq * u_sq * (-1. / 128. + u_sq * (3. / 512. + u_sq * -35. / 8192.));
        let b6 = u_sq * u_sq * u_sq * (-1. / 1536. + u_sq * 5. / 6144.);

        let sigma_0 = s12 / (self.ellipsoid.b() * b0);

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
            cos_6_sigma_m = cos_2_sigma_m.powi(3) - 3. * cos_2_sigma_m * (1. - cos_2_sigma_m.powi(2));

            let next_sigma = sigma_0 - b2 / b0 * sin_sigma * cos_2_sigma_m - b4/ b0 * sin_2_sigma * cos_4_sigma_m - b6 / b0 * sin_3_sigma * cos_6_sigma_m;
            let sigma_old = sigma;
            sigma = next_sigma;
            if (sigma - sigma_old).abs() < 0.5e-14 {
                break;
            }
            if iter_count == Self::MAX_ITERATIONS {
                panic!("Direct geodesic computation is not converging after {} iterations!", Self::MAX_ITERATIONS);
            }
        }

        let alpha = sin_alpha.asin();
        // From Eq (1.22)
        let beta2 = ((sigma_1 + sigma).sin() * alpha.cos()).asin();
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
        let big_as = self.axs.map(|p| p.eval_at(cos_alpha_sq));
        let lambda_minus_L = f * sin_alpha * (
            big_as[0] * sigma +
                big_as[1] * sin_sigma * cos_2_sigma_m +
                big_as[2] * sin_2_sigma * cos_4_sigma_m +
                big_as[3] * sin_3_sigma * cos_6_sigma_m
        );

        let L = lambda - lambda_minus_L;

        let alpha21 = if s12 < 1000.0 {
            let sin_lambda_2_sq = (1. - cos_lambda) / 2.;
            // Eq (1.73)
            (sin_lambda * cos_beta1).atan2(
                 (beta2 - beta1).sin()
                - 2. * cos_beta1 * sin_beta2 * sin_lambda_2_sq)
        } else {
            // Eq (1.71)
            (sin_lambda * cos_beta1).atan2(
                sin_beta2 * cos_beta1 * cos_lambda - sin_beta1 * cos_beta2)
        };

        Geodesic {
            p1,
            alpha1,
            p2: (lon1 + L, Lat::new(lat2)),
            alpha2: Azimuth::new(alpha21),
            s: s12
        }
    }


    pub fn inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Geodesic {
        // Init
        let (lon1, lat1) = p1;
        let (lon2, lat2) = p2;
        let beta1 = self.ellipsoid.reduced_latitude(lat1.rad());
        let beta2 = self.ellipsoid.reduced_latitude(lat2.rad());
        let (sin_beta1, cos_beta1) = beta1.sin_cos();
        let (sin_beta2, cos_beta2) = beta2.sin_cos();

        let L = (lon2 - lon1).length().unwrap_or(0.);
        let mut lambda = L;
        let (mut sin_lambda, mut cos_lambda);
        let mut sigma;
        let mut sin_sigma;
        let mut sin_2_sigma;
        let mut sin_3_sigma;
        let mut sin_4_sigma;
        let mut sin_alpha;
        let mut cos_alpha_sq;
        let mut cos_2_sigma_m;
        let mut cos_4_sigma_m;
        let mut cos_6_sigma_m;
        let mut cos_8_sigma_m;
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

            // Eq (1.66)
            cos_2_sigma_m = cos_sigma - 2. * sin_beta1 * sin_beta2 / cos_alpha_sq;
            cos_4_sigma_m = 2. * cos_2_sigma_m.powi(2) - 1.;
            cos_6_sigma_m = cos_2_sigma_m.powi(3) - 3. * cos_2_sigma_m * (1. - cos_2_sigma_m.powi(2));
            cos_8_sigma_m = cos_2_sigma_m.powi(4) - 6. * cos_2_sigma_m.powi(2) * sin_2_sigma.powi(2) + sin_2_sigma.powi(4);

            // Eq (1.56)
            let f = self.ellipsoid.f();
            let big_as = self.axs.map(|p| p.eval_at(cos_alpha_sq));
            let lambda_minus_L = f * sin_alpha * (
                big_as[0] * sigma +
                    big_as[1] * sin_sigma * cos_2_sigma_m +
                    big_as[2] * sin_2_sigma * cos_4_sigma_m +
                    big_as[3] * sin_3_sigma * cos_6_sigma_m
            );

            let lambda_old = lambda;
            lambda = L + lambda_minus_L;
            if (lambda - lambda_old).abs() < 0.5e-14 {
                break;
            }
            if iter_count == Self::MAX_ITERATIONS {
                panic!("Inverse geodesic computation is not converging after {} iterations!", Self::MAX_ITERATIONS);
            }
        }

        let u_sq = self.ellipsoid.e_prime_sq() * cos_alpha_sq;
        let b0 = 1. + u_sq * (1. / 4. + u_sq * (-3. / 64. + u_sq * (5. / 256. + u_sq * -175. / 16384.)));
        let b2 = u_sq * (-1. / 4. + u_sq * (1. / 16. + u_sq * (-15. / 512. + u_sq * 35. / 2048.)));
        let b4 = u_sq * u_sq * (-1. / 128. + u_sq * (3. / 512. + u_sq * -35. / 8192.));
        let b6 = u_sq * u_sq * u_sq * (-1. / 1536. + u_sq * 5. / 6144.);
        let b8 = u_sq * u_sq * u_sq * u_sq * (-5. / 65536.);

        let s = self.ellipsoid.b() * (
            b0 * sigma +
                b2 * sin_sigma * cos_2_sigma_m +
                b4 * sin_2_sigma * cos_4_sigma_m +
                b6 * sin_3_sigma * cos_6_sigma_m +
                b8 * sin_4_sigma * cos_8_sigma_m
        );
        let sin_alpha1 = sin_alpha / cos_beta1;
        let sin_alpha2 = sin_alpha / cos_beta2;
        let sin_lambda_2_sq = (1. - cos_lambda) / 2.;

        let tan_alpha12 = if s < 10_000.0 {
            // Eq (1.72)
            lambda.sin() * cos_beta2
                / ((beta2 - beta1).sin()
                    + 2. * sin_beta1 * cos_beta2 * sin_lambda_2_sq)
        } else {
            // Eq (1.70)
            lambda.sin() * cos_beta2
                / (sin_beta2 * cos_beta1 - lambda.cos() * sin_beta1 * cos_beta2)
        };
        let tan_alpha21 = if s < 10_000.0 {
            // Eq (1.73)
            lambda.sin() * cos_beta1
                / ((beta2 - beta1).sin()
                    - 2. * cos_beta1 * sin_beta2 * sin_lambda_2_sq)
        } else {
            // Eq (1.71)
            lambda.sin() * cos_beta1
                / (sin_beta2 * cos_beta1 * lambda.cos() - sin_beta1 * cos_beta2)
        };

        let az1 = sin_alpha1.atan2(sin_alpha1 / tan_alpha12);
        let az2 = sin_alpha2.atan2(sin_alpha2 / tan_alpha21);

        Geodesic {
            p1: (lon1, lat1),
            alpha1: Azimuth::new(az1),
            p2: (lon2, lat2),
            alpha2: Azimuth::new(az2),
            s,
        }
    }

}

impl<'e> GeodesicSolver for RappIterativeGeodisicSolver<'e> {
    fn solve_direct(&self, p1: (Lon, Lat), alpha1: Azimuth, s12: f64) -> Geodesic {
        self.direct(p1, alpha1, s12)
    }

    fn solve_inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Geodesic {
        self.inverse(p1, p2)
    }
}

#[cfg(test)]
mod tests {
    use crate::geodesy::geodesics::rapp::RappIterativeGeodisicSolver;
    use crate::geodesy::geodesics::tests::{antipodal_lines, standard_lines, vincenty_direct_deltas, vincenty_inverse_deltas, vincenty_lines};
    use crate::geodesy::geodesics::GeodesicSolver;
    use approx::assert_abs_diff_eq;
    use crate::quantity::angle::DMS;

    #[test]
    fn direct_vincenty_lines() {
        let chars = ['a', 'b', 'c', 'd', 'e', 'f'];
        for (ix, (input, delta)) in vincenty_lines()
            .into_iter()
            .zip(vincenty_direct_deltas())
            .enumerate()
        {
            let solver = RappIterativeGeodisicSolver::new(&input.ellipsoid);
            let computed =
                solver.solve_direct(input.geodesic.p1, input.geodesic.alpha1, input.geodesic.s);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Vincenty Line ({})", chars[ix]);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_lon_dms = DMS::from_rad(computed.p2.0.rad() - input.geodesic.p2.0.rad());
            println!("error on lon (res - exp) = {}", diff_lon_dms);

            let diff_lat_dms = DMS::from_rad(computed.p2.1.rad() - input.geodesic.p2.1.rad());
            println!("error on lat (res - exp) = {}", diff_lat_dms);

            let diff_az_dms = DMS::from_rad(computed.alpha2.rad() - input.geodesic.alpha2.rad());
            println!("error on az (res - exp) = {}", diff_az_dms);

            assert_eq!(diff_lon_dms.deg(), 0.0);
            assert_eq!(diff_lon_dms.min(), 0.0);
            assert!(diff_lon_dms.sec() < delta.delta_delta_lon.abs());

            assert_eq!(diff_lat_dms.deg(), 0.0);
            assert_eq!(diff_lat_dms.min(), 0.0);
            assert!(diff_lat_dms.sec() < delta.delta_lat2.abs());

            assert_eq!(diff_az_dms.deg(), 0.0);
            assert_eq!(diff_az_dms.min(), 0.0);
            assert!(diff_az_dms.sec() < delta.delta_alpha2.abs());
        }
    }

    #[test]
    fn direct_standard_lines() {
        // Direct problem
        for (ix, input) in standard_lines().iter().enumerate() {
            let solver = RappIterativeGeodisicSolver::new(&input.ellipsoid);
            let computed = solver.solve_direct(input.geodesic.p1, input.geodesic.alpha1, input.geodesic.s);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Standard Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_lon_dms = DMS::from_rad(computed.p2.0.rad() - input.geodesic.p2.0.rad());
            println!("error on lon (res - exp) = {}", diff_lon_dms);

            let diff_lat_dms = DMS::from_rad(computed.p2.1.rad() - input.geodesic.p2.1.rad());
            println!("error on lat (res - exp) = {}", diff_lat_dms);

            let diff_az_dms = DMS::from_rad(computed.alpha2.rad() - input.geodesic.alpha2.rad());
            println!("error on az (res - exp) = {}", diff_az_dms);

            assert_abs_diff_eq!(computed.p2.0, input.geodesic.p2.0, epsilon = 2e-9);
            assert_abs_diff_eq!(computed.p2.1, input.geodesic.p2.1, epsilon = 2e-9);
            assert_abs_diff_eq!(computed.alpha2, input.geodesic.alpha2, epsilon = 1e-7);
        }
    }

    #[test]
    fn direct_antipodal_lines() {
        for (ix, input) in antipodal_lines().iter().enumerate() {
            let solver = RappIterativeGeodisicSolver::new(&input.ellipsoid);
            let computed = solver.solve_direct(input.geodesic.p1, input.geodesic.alpha1, input.geodesic.s);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Antipodal Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_lon_dms = DMS::from_rad(computed.p2.0.rad() - input.geodesic.p2.0.rad());
            println!("error on lon (res - exp) = {}", diff_lon_dms);

            let diff_lat_dms = DMS::from_rad(computed.p2.1.rad() - input.geodesic.p2.1.rad());
            println!("error on lat (res - exp) = {}", diff_lat_dms);

            let diff_az_dms = DMS::from_rad(computed.alpha2.rad() - input.geodesic.alpha2.rad());
            println!("error on az (res - exp) = {}", diff_az_dms);

            assert_abs_diff_eq!(computed.p2.0, input.geodesic.p2.0, epsilon = 2e-9);
            assert_abs_diff_eq!(computed.p2.1, input.geodesic.p2.1, epsilon = 2e-9);
            assert_abs_diff_eq!(computed.alpha2, input.geodesic.alpha2, epsilon = 1e-7);
        }
    }

    #[test]
    fn inverse_vincenty_lines() {
        // Vincenty lines
        let chars = ['a', 'b', 'c', 'd', 'e', 'f'];
        for (ix, (input, delta)) in vincenty_lines()
            .into_iter()
            .zip(vincenty_inverse_deltas())
            .enumerate()
        {
            let solver = RappIterativeGeodisicSolver::new(&input.ellipsoid);
            let computed = solver.solve_inverse(input.geodesic.p1, input.geodesic.p2);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Vincenty Line ({})", chars[ix]);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_az1_dms = DMS::from_rad(computed.alpha1.rad() - input.geodesic.alpha1.rad());
            println!("error on az1 (computed - input) = {}", diff_az1_dms);

            let diff_az2_dms = DMS::from_rad(computed.alpha2.rad() - input.geodesic.alpha2.rad());
            println!("error on az2 (computed - input) = {}", diff_az2_dms);

            let diff_s_m = computed.s - input.geodesic.s;
            println!("error on s (computed - input) = {}", diff_s_m);

            assert_eq!(diff_az1_dms.deg(), 0.0);
            assert_eq!(diff_az1_dms.min(), 0.0);
            assert!(diff_az1_dms.sec() < delta.delta_alpha1.abs());

            assert_eq!(diff_az2_dms.deg(), 0.0);
            assert_eq!(diff_az2_dms.min(), 0.0);
            assert!(diff_az2_dms.sec() < delta.delta_alpha2.abs());

            assert!(diff_s_m < delta.delta_s.abs());
        }
    }

    #[test]
    fn inverse_standard_lines() {
        // Standard lines
        for (ix, input) in standard_lines().iter().enumerate() {
            let solver = RappIterativeGeodisicSolver::new(&input.ellipsoid);
            let computed = solver.solve_inverse(input.geodesic.p1, input.geodesic.p2);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Standard Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_az1_dms = DMS::from_rad(computed.alpha1.rad() - input.geodesic.alpha1.rad());
            println!("error on az1 (computed - input) = {}", diff_az1_dms);

            let diff_az2_dms = DMS::from_rad(computed.alpha2.rad() - input.geodesic.alpha2.rad());
            println!("error on az2 (computed - input) = {}", diff_az2_dms);

            let diff_s_m = computed.s - input.geodesic.s;
            println!("error on s (computed - input) = {}", diff_s_m);

            // assert_abs_diff_eq!(computed.alpha1, input.geodesic.alpha1, epsilon = 2e-10);
            // assert_abs_diff_eq!(computed.alpha2, input.geodesic.alpha2, epsilon = 1e-10);
            // assert_abs_diff_eq!(computed.s, input.geodesic.s, epsilon = 1e-3);
        }
    }

    #[test]
    fn inverse_antipodal_lines() {
        // Anti-podal lines
        for (ix, input) in antipodal_lines().iter().enumerate() {
            let solver = RappIterativeGeodisicSolver::new(&input.ellipsoid);
            let computed = solver.solve_inverse(input.geodesic.p1, input.geodesic.p2);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Antipodal Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!();

            let diff_az1_dms = DMS::from_rad(computed.alpha1.rad() - input.geodesic.alpha1.rad());
            println!("error on az1 (computed - input) = {}", diff_az1_dms);

            let diff_az2_dms = DMS::from_rad(computed.alpha2.rad() - input.geodesic.alpha2.rad());
            println!("error on az2 (computed - input) = {}", diff_az2_dms);

            let diff_s_m = computed.s - input.geodesic.s;
            println!("error on s (computed - input) = {}", diff_s_m);

            // assert_abs_diff_eq!(computed.alpha1, input.geodesic.alpha1, epsilon = 2e-10);
            // assert_abs_diff_eq!(computed.alpha2, input.geodesic.alpha2, epsilon = 1e-10);
            // assert_abs_diff_eq!(computed.s, input.geodesic.s, epsilon = 1e-3);
        }
    }
}
