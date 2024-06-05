use crate::cs::geodetic::{Lat, Lon};
use crate::cs::Azimuth;
use crate::geodesy::geodesics::{Geodesic, GeodesicSolver};
use crate::geodesy::Ellipsoid;

pub struct VincentyGeodesicSolver<'e> {
    ellipsoid: &'e Ellipsoid,
}

impl<'e> VincentyGeodesicSolver<'e> {
    const MAX_ITERATIONS: usize = 30;

    pub fn new(ellipsoid: &'e Ellipsoid) -> Self {
        Self { ellipsoid }
    }

    fn direct(&self, p1: (Lon, Lat), alpha1: Azimuth, s12: f64) -> Geodesic {
        let f = self.ellipsoid.f();
        let (lon1, lat1) = p1;
        let tan_beta1 = (1. - f) * lat1.rad().tan();
        let cos_beta1 = 1. / (1. + tan_beta1 * tan_beta1).sqrt(); // cos_beta1 >= 0.
        let sin_beta1 = tan_beta1 * cos_beta1;

        // Eq (1): indeterminate for equatorial lines !!!
        let (sin_alpha1, cos_alpha1) = alpha1.rad().sin_cos();
        let sigma1 = tan_beta1.atan2(cos_alpha1);

        // Eq (2)
        let sin_alpha = cos_beta1 * sin_alpha1;
        let cos_alpha_sq = 1. - sin_alpha.powi(2);

        // Eq (3) & (4)
        let u_sq = cos_alpha_sq * self.ellipsoid.e_prime_sq();
        let A = 1. + (u_sq / 16384.) * (4096. + u_sq * (-768. + u_sq * (320. - 175. * u_sq)));
        let B = (u_sq / 1024.) * (256. + u_sq * (-128. + u_sq * (74. - 47. * u_sq)));

        let sigma_0 = s12 / (self.ellipsoid.b() * A);
        let mut sigma = sigma_0;
        let (mut sin_sigma, mut cos_sigma) = sigma.sin_cos();
        let mut cos_2_sigma_m = (2. * sigma1 + sigma).cos();
        let mut iter_count = 0;
        loop {
            iter_count += 1;

            let sigma_prev = sigma;

            // Eq (6)
            let delta_sigma = B * sin_sigma * (
                cos_2_sigma_m + (B / 4.) * (
                    cos_sigma * (-1. + 2. * cos_2_sigma_m.powi(2)) - (B / 6.) * cos_2_sigma_m * (-3. + 4. * sin_sigma.powi(2)) * (-3. + 4. * cos_2_sigma_m.powi(2))
                )
            );
            sigma = sigma_0 + delta_sigma;
            (sin_sigma, cos_sigma) = sigma.sin_cos();
            // Eq (5)
            cos_2_sigma_m = (2. * sigma1 + sigma).cos();

            if (sigma - sigma_prev).abs() < 1e-12 {
                break;
            }
            if iter_count == Self::MAX_ITERATIONS {
                panic!("Direct geodesic computation is not converging after {} iterations!", Self::MAX_ITERATIONS);
            }
        }

        // Eq (8)
        let lat2 = (sin_beta1 * cos_sigma + cos_beta1 * sin_sigma * cos_alpha1).atan2(
            (1. - f) * (sin_alpha.powi(2) + (sin_beta1 * sin_sigma - cos_beta1 * cos_sigma * cos_alpha1).powi(2)).sqrt()
        );

        // Eq (9)
        let lambda = (sin_sigma * sin_alpha1)
            .atan2((cos_beta1 * cos_sigma - sin_beta1 * sin_sigma * cos_alpha1));

        // Eq (10)
        let C = (f / 16.) * cos_alpha_sq * (4. + f * (4. - 3. * cos_alpha_sq));

        // Eq (11) decomposed
        let lambda_minus_L = (1. - C)
            * f
            * sin_alpha
            * (sigma
                + C * sin_sigma * (cos_2_sigma_m + C * cos_sigma * (-1. + 2. * cos_2_sigma_m.powi(2))));
        let L = lambda - lambda_minus_L;

        let alpha2 = sin_alpha.atan2(-sin_beta1 * sin_sigma + cos_beta1 * cos_sigma * cos_alpha1);

        Geodesic {
            p1: (lon1, lat1),
            alpha1,
            p2: (lon1 + L, Lat::new(lat2)),
            alpha2: Azimuth::new(alpha2),
            s: s12,
        }
    }

    fn inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Geodesic {
        let (lon1, lat1) = p1;
        let (lon2, lat2) = p2;
        let beta1 = self.ellipsoid.reduced_latitude(lat1.rad());
        let beta2 = self.ellipsoid.reduced_latitude(lat2.rad());
        let (sin_beta1, cos_beta1) = beta1.sin_cos();
        let (sin_beta2, cos_beta2) = beta2.sin_cos();

        // Eq (13) 1st approximation
        let L = (lon2 - lon1).length().unwrap_or(0.0);
        let mut lambda = L;
        let (mut sin_lambda, mut cos_lambda) = lambda.sin_cos();
        let mut cos_alpha_sq;
        let mut sigma;
        let (mut sin_sigma, mut cos_sigma);
        let mut cos_2_sigma_m;
        let mut iter_count = 0;
        loop {
            iter_count += 1;
            // Eq (14)
            let sin_sigma_sq = (cos_beta2 * sin_lambda).powi(2)
                + (cos_beta1 * sin_beta2 - sin_beta1 * cos_beta2 * cos_lambda).powi(2);

            // Eq (15)
            let cos_sigma_i = sin_beta1 * sin_beta2 + cos_beta1 * cos_beta2 * cos_lambda;

            // Eq (16): tan(sigma) = sin(sigma) / cos(sigma)
            sigma = sin_sigma_sq.sqrt().atan2(cos_sigma_i);
            (sin_sigma, cos_sigma) = sigma.sin_cos();

            // Eq (17)
            let sin_alpha = cos_beta1 * cos_beta2 * sin_lambda / sin_sigma;
            cos_alpha_sq = 1. - sin_alpha.powi(2);

            // Eq (18)
            cos_2_sigma_m = cos_sigma - 2. * sin_beta1 * sin_beta2 / cos_alpha_sq;

            // Update lambda
            let lambda_prev = lambda;
            let f = self.ellipsoid.f();
            // Eq (10)
            let c = f / 16. * cos_alpha_sq * (4. + f * (4. - 3. * cos_alpha_sq));
            // Eq (11)
            let lambda_minus_L = (1. - c)
                * f
                * sin_alpha
                * (sigma
                    + c * sin_sigma * (cos_2_sigma_m + c * cos_sigma * (-1. + 2. * cos_2_sigma_m.powi(2))));
            lambda = lambda_minus_L + L;
            (sin_lambda, cos_lambda) = lambda.sin_cos();
            if (lambda - lambda_prev).abs() < 1e-12 {
                break;
            }
            if iter_count == Self::MAX_ITERATIONS {
                panic!("Inverse geodesic computation is not converging after {} iterations!", Self::MAX_ITERATIONS);
            }
        }

        let u_sq = cos_alpha_sq * self.ellipsoid.e_prime_sq();
        // Eq (3)
        let A = 1. + u_sq / 16384. * (4096. + u_sq * (-768. + u_sq * (320. - 175. * u_sq)));
        // Eq (4)
        let B = u_sq / 1024. * (256. + u_sq * (-128. + u_sq * (74. - 47. * u_sq)));

        // Eq (6)
        let delta_sigma = B * sin_sigma
            * (cos_2_sigma_m
                + B / 4.
                    * (cos_sigma * (-1. + 2. * cos_2_sigma_m.powi(2))
                        - B / 6.
                            * cos_2_sigma_m
                            * (-3. + 4. * sin_sigma.powi(2))
                            * (-3. + 4. * cos_2_sigma_m.powi(2))));
        // Eq (19)
        let s = self.ellipsoid.b() * A * (sigma - delta_sigma);

        // Eq (20)
        let alpha1 = (cos_beta2 * sin_lambda)
            .atan2(cos_beta1 * sin_beta2 - sin_beta1 * cos_beta2 * cos_lambda);
        // Eq (21)
        let alpha2 = (cos_beta1 * sin_lambda)
            .atan2(-sin_beta1 * cos_beta2 + cos_beta1 * sin_beta2 * cos_lambda);

        Geodesic {
            p1,
            alpha1: Azimuth::new(alpha1),
            p2,
            alpha2: Azimuth::new(alpha2),
            s,
        }
    }
}

impl<'e> GeodesicSolver for VincentyGeodesicSolver<'e> {
    fn solve_direct(&self, p1: (Lon, Lat), alpha1: Azimuth, s12: f64) -> Geodesic {
        self.direct(p1, alpha1, s12)
    }

    fn solve_inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Geodesic {
        self.inverse(p1, p2)
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use crate::geodesy::geodesics::tests::{vincenty_lines, DirectDeltas, InverseDeltas, LineData, vincenty_direct_deltas, vincenty_inverse_deltas, standard_lines, antipodal_lines};
    use crate::geodesy::geodesics::vincenty::VincentyGeodesicSolver;
    use crate::geodesy::geodesics::{GeodesicSolver};
    use crate::geodesy::geodesics::rapp::RappIterativeGeodisicSolver;
    use crate::quantity::angle::DMS;

    #[test]
    fn direct_vincenty_lines() {
        let chars = ['a', 'b', 'c', 'd', 'e', 'f'];
        for (ix, (input, delta)) in vincenty_lines()
            .into_iter()
            .zip(vincenty_direct_deltas())
            .enumerate()
        {
            let solver = VincentyGeodesicSolver::new(&input.ellipsoid);
            let computed =
                solver.solve_direct(input.geodesic.p1, input.geodesic.alpha1, input.geodesic.s);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Vincenty Line ({})", chars[ix]);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!("");

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
        for (ix, input) in standard_lines().iter().enumerate() {
            let solver = VincentyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver.solve_direct(input.geodesic.p1, input.geodesic.alpha1, input.geodesic.s);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Standard Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!("");

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
            let solver = VincentyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver.solve_direct(input.geodesic.p1, input.geodesic.alpha1, input.geodesic.s);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Antipodal Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!("");

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
            let solver = VincentyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver.solve_inverse(input.geodesic.p1, input.geodesic.p2);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Vincenty Line ({})", chars[ix]);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!("");

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
            let solver = VincentyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver.solve_inverse(input.geodesic.p1, input.geodesic.p2);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Standard Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!("");

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
            let solver = VincentyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver.solve_inverse(input.geodesic.p1, input.geodesic.p2);
            println!("-----------------------------------------------------------------------------------");
            println!("Input: Antipodal Line ({})", ix);
            println!("{}", input.geodesic);
            println!("Computed: ");
            println!("{}", computed);
            println!("");

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
