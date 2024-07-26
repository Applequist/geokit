#![allow(non_snake_case)]

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

    fn direct(&self, p1: (Lon, Lat), alpha1: Azimuth, s12: f64) -> Result<Geodesic, &'static str> {
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
            let delta_sigma = B
                * sin_sigma
                * (cos_2_sigma_m
                    + (B / 4.)
                        * (cos_sigma * (-1. + 2. * cos_2_sigma_m.powi(2))
                            - (B / 6.)
                                * cos_2_sigma_m
                                * (-3. + 4. * sin_sigma.powi(2))
                                * (-3. + 4. * cos_2_sigma_m.powi(2))));
            sigma = sigma_0 + delta_sigma;
            (sin_sigma, cos_sigma) = sigma.sin_cos();
            // Eq (5)
            cos_2_sigma_m = (2. * sigma1 + sigma).cos();

            if (sigma - sigma_prev).abs() < 1e-12 {
                break;
            }
            if iter_count == Self::MAX_ITERATIONS {
                return Err("Direct geodesic computation is not converging!");
            }
        }

        // Eq (8)
        let lat2 = (sin_beta1 * cos_sigma + cos_beta1 * sin_sigma * cos_alpha1).atan2(
            (1. - f)
                * (sin_alpha.powi(2)
                    + (sin_beta1 * sin_sigma - cos_beta1 * cos_sigma * cos_alpha1).powi(2))
                .sqrt(),
        );

        // Eq (9)
        let lambda = (sin_sigma * sin_alpha1)
            .atan2(cos_beta1 * cos_sigma - sin_beta1 * sin_sigma * cos_alpha1);

        // Eq (10)
        let C = (f / 16.) * cos_alpha_sq * (4. + f * (4. - 3. * cos_alpha_sq));

        // Eq (11) decomposed
        let lambda_minus_L = (1. - C)
            * f
            * sin_alpha
            * (sigma
                + C * sin_sigma
                    * (cos_2_sigma_m + C * cos_sigma * (-1. + 2. * cos_2_sigma_m.powi(2))));
        let L = lambda - lambda_minus_L;

        let alpha2 = sin_alpha.atan2(-sin_beta1 * sin_sigma + cos_beta1 * cos_sigma * cos_alpha1);

        Ok(Geodesic {
            p1: (lon1, lat1),
            alpha1,
            p2: (lon1 + L, Lat::new(lat2)),
            alpha2: Azimuth::new(alpha2),
            s: s12,
        })
    }

    fn inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Result<Geodesic, &'static str> {
        let f = self.ellipsoid.f();
        let (lon1, lat1) = p1;
        let (lon2, lat2) = p2;
        let beta1 = self.ellipsoid.reduced_latitude(lat1.rad());
        let beta2 = self.ellipsoid.reduced_latitude(lat2.rad());
        let (sin_beta1, cos_beta1) = beta1.sin_cos();
        let (sin_beta2, cos_beta2) = beta2.sin_cos();

        // Eq (13) 1st approximation
        let L = (lon2 - lon1).length();
        let mut lambda = L;
        let (mut sin_lambda, mut cos_lambda) = lambda.sin_cos();
        let mut cos_alpha_sq;
        let mut sigma;
        let (mut sin_sigma, mut cos_sigma);
        let mut cos_2_sigma_m = 0.0;
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

            // If the geodesic (p1, p2) is an equatorial line Eq (18) becomes indeterminate but
            // B and C become 0.0
            let lambda_minus_L = if cos_alpha_sq == 0.0 {
                f * sin_alpha * sigma
            } else {
                // Eq (18)
                cos_2_sigma_m = cos_sigma - 2. * sin_beta1 * sin_beta2 / cos_alpha_sq;

                // Eq (10)
                let upper_c = f / 16. * cos_alpha_sq * (4. + f * (4. - 3. * cos_alpha_sq));
                // Eq (11)
                (1. - upper_c)
                    * f
                    * sin_alpha
                    * (sigma
                        + upper_c
                            * sin_sigma
                            * (cos_2_sigma_m
                                + upper_c * cos_sigma * (-1. + 2. * cos_2_sigma_m.powi(2))))
            };

            // Update lambda
            let lambda_prev = lambda;
            lambda = lambda_minus_L + L;
            (sin_lambda, cos_lambda) = lambda.sin_cos();
            if (lambda - lambda_prev).abs() < 1e-12 {
                break;
            }
            if iter_count == Self::MAX_ITERATIONS {
                return Err("Inverse geodesic computation is not converging!");
            }
        }

        let s = if cos_alpha_sq == 0.0 {
            // A = 1. B = 0. -> delta_sigma = 0 s = b * A(=1) * sigma
            self.ellipsoid.b() * sigma
        } else {
            let u_sq = cos_alpha_sq * self.ellipsoid.e_prime_sq();
            // Eq (3)
            let A = 1. + u_sq / 16384. * (4096. + u_sq * (-768. + u_sq * (320. - 175. * u_sq)));
            // Eq (4)
            let B = u_sq / 1024. * (256. + u_sq * (-128. + u_sq * (74. - 47. * u_sq)));

            // Eq (6)
            let delta_sigma = B
                * sin_sigma
                * (cos_2_sigma_m
                    + B / 4.
                        * (cos_sigma * (-1. + 2. * cos_2_sigma_m.powi(2))
                            - B / 6.
                                * cos_2_sigma_m
                                * (-3. + 4. * sin_sigma.powi(2))
                                * (-3. + 4. * cos_2_sigma_m.powi(2))));

            // Eq (19)
            self.ellipsoid.b() * A * (sigma - delta_sigma)
        };

        // Eq (20)
        let alpha1 = (cos_beta2 * sin_lambda)
            .atan2(cos_beta1 * sin_beta2 - sin_beta1 * cos_beta2 * cos_lambda);
        // Eq (21)
        let alpha2 = (cos_beta1 * sin_lambda)
            .atan2(-sin_beta1 * cos_beta2 + cos_beta1 * sin_beta2 * cos_lambda);

        Ok(Geodesic {
            p1,
            alpha1: Azimuth::new(alpha1),
            p2,
            alpha2: Azimuth::new(alpha2),
            s,
        })
    }
}

impl<'e> GeodesicSolver for VincentyGeodesicSolver<'e> {
    fn solve_direct(
        &self,
        p1: (Lon, Lat),
        alpha1: Azimuth,
        s12: f64,
    ) -> Result<Geodesic, &'static str> {
        self.direct(p1, alpha1, s12)
    }

    fn solve_inverse(&self, p1: (Lon, Lat), p2: (Lon, Lat)) -> Result<Geodesic, &'static str> {
        self.inverse(p1, p2)
    }
}

#[cfg(test)]
mod tests {
    use crate::cs::geodetic::{Lat, Lon};
    use crate::cs::Azimuth;
    use crate::geodesy::ellipsoid::consts;
    use crate::geodesy::geodesics::tests::{
        antipodal_lines, standard_lines, vincenty_direct_deltas, vincenty_inverse_deltas,
        vincenty_lines,
    };
    use crate::geodesy::geodesics::vincenty::VincentyGeodesicSolver;
    use crate::geodesy::geodesics::GeodesicSolver;
    use crate::quantity::angle::units::DEG;
    use crate::quantity::angle::DMS;
    use approx::assert_abs_diff_eq;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn direct_vincenty_lines() {
        let chars = ['a', 'b', 'c', 'd', 'e', 'f'];
        for (ix, (input, delta)) in vincenty_lines()
            .into_iter()
            .zip(vincenty_direct_deltas())
            .enumerate()
        {
            let solver = VincentyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver
                .solve_direct(input.geodesic.p1, input.geodesic.alpha1, input.geodesic.s)
                .unwrap();
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
        for (ix, input) in standard_lines().iter().enumerate() {
            let solver = VincentyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver
                .solve_direct(input.geodesic.p1, input.geodesic.alpha1, input.geodesic.s)
                .unwrap();
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

            assert_abs_diff_eq!(computed.p2.0, input.geodesic.p2.0, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.p2.1, input.geodesic.p2.1, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.alpha2, input.geodesic.alpha2, epsilon = 1e-10);
        }
    }

    #[test]
    fn direct_antipodal_lines() {
        for (ix, input) in antipodal_lines().iter().enumerate() {
            let solver = VincentyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver
                .solve_direct(input.geodesic.p1, input.geodesic.alpha1, input.geodesic.s)
                .unwrap();
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

            assert_abs_diff_eq!(computed.p2.0, input.geodesic.p2.0, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.p2.1, input.geodesic.p2.1, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.alpha2, input.geodesic.alpha2, epsilon = 1e-10);
        }
    }

    #[test]
    fn direct_equator_lines() {
        let wgs84 = consts::WGS84;
        let solver = VincentyGeodesicSolver::new(&wgs84);
        let computed = solver
            .solve_direct(
                (Lon::new(0.0), Lat::new(0.0)),
                Azimuth::new(90.0 * DEG),
                20_000.0,
            )
            .unwrap();
        assert_abs_diff_eq!(computed.p2.0.rad(), 0.17966306 * DEG, epsilon = 1e-8);
        assert_abs_diff_eq!(computed.p2.1.rad(), 0.0, epsilon = 1e-8);
        assert_abs_diff_eq!(computed.alpha2.rad(), FRAC_PI_2, epsilon = 1e-8);

        let computed = solver
            .solve_direct(
                (Lon::new(170.0 * DEG), Lat::new(0.0)),
                Azimuth::new(90.0 * DEG),
                2_000_000.0,
            )
            .unwrap();
        assert_abs_diff_eq!(
            computed.p2.0,
            Lon::new(-172.03369432 * DEG),
            epsilon = 1e-10
        );
        assert_abs_diff_eq!(computed.p2.1, Lat::new(0.0), epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2.rad(), FRAC_PI_2, epsilon = 1e-8);
    }

    #[test]
    fn direct_meridian_lines() {
        let wgs84 = consts::WGS84;
        let solver = VincentyGeodesicSolver::new(&wgs84);

        let computed = solver
            .solve_direct(
                (Lon::new(0.), Lat::new(-10. * DEG)),
                Azimuth::new(0.),
                2_000_000.0,
            )
            .unwrap();
        assert_abs_diff_eq!(computed.p2.0.rad(), 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.p2.1.rad(), 8.08583903 * DEG, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2.rad(), 0.0, epsilon = 1e-8);

        let computed = solver
            .solve_direct(
                (Lon::new(0.), Lat::new(80. * DEG)),
                Azimuth::new(0.),
                2_000_000.0,
            )
            .unwrap();
        assert_abs_diff_eq!(computed.p2.0.rad(), PI, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.p2.1.rad(), 82.09240627 * DEG, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2.rad(), PI, epsilon = 1e-8);
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
            let computed = solver
                .solve_inverse(input.geodesic.p1, input.geodesic.p2)
                .unwrap();
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
            let solver = VincentyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver
                .solve_inverse(input.geodesic.p1, input.geodesic.p2)
                .unwrap();
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

            assert_abs_diff_eq!(computed.alpha1, input.geodesic.alpha1, epsilon = 1e-9);
            assert_abs_diff_eq!(computed.alpha2, input.geodesic.alpha2, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.s, input.geodesic.s, epsilon = 1e-3);
        }
    }

    #[test]
    fn inverse_antipodal_lines() {
        // Anti-podal lines
        for (ix, input) in antipodal_lines().iter().enumerate() {
            let solver = VincentyGeodesicSolver::new(&input.ellipsoid);
            let computed = solver
                .solve_inverse(input.geodesic.p1, input.geodesic.p2)
                .unwrap();
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

            assert_abs_diff_eq!(computed.alpha1, input.geodesic.alpha1, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.alpha2, input.geodesic.alpha2, epsilon = 1e-10);
            assert_abs_diff_eq!(computed.s, input.geodesic.s, epsilon = 1e-3);
        }
    }

    #[test]
    fn inverse_equator_lines() {
        let wgs84 = consts::WGS84;
        let solver = VincentyGeodesicSolver::new(&wgs84);
        let computed = solver
            .solve_inverse(
                (Lon::new(-10. * DEG), Lat::new(0.)),
                (Lon::new(10. * DEG), Lat::new(0.)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1, Azimuth::EAST, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2, Azimuth::EAST, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, 2_226_389.816, epsilon = 1e-3);

        let computed = solver
            .solve_inverse(
                (Lon::new(10. * DEG), Lat::new(0.)),
                (Lon::new(-10. * DEG), Lat::new(0.)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1, Azimuth::WEST, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2, Azimuth::WEST, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, 2_226_389.816, epsilon = 1e-3);

        let computed = solver
            .solve_inverse(
                (Lon::new(170. * DEG), Lat::new(0.)),
                (Lon::new(-170. * DEG), Lat::new(0.)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1, Azimuth::EAST, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2, Azimuth::EAST, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, 2_226_389.816, epsilon = 1e-3);
    }

    #[test]
    fn inverse_meridian_lines() {
        let wgs84 = consts::WGS84;
        let solver = VincentyGeodesicSolver::new(&wgs84);
        let computed = solver
            .solve_inverse(
                (Lon::new(0. * DEG), Lat::new(-10. * DEG)),
                (Lon::new(0.), Lat::new(10. * DEG)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1, Azimuth::NORTH, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2, Azimuth::NORTH, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, 2_211_709.666, epsilon = 1e-3);

        let computed = solver
            .solve_inverse(
                (Lon::new(0. * DEG), Lat::new(10. * DEG)),
                (Lon::new(0. * DEG), Lat::new(-10. * DEG)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1, Azimuth::SOUTH, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2, Azimuth::SOUTH, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, 2_211_709.666, epsilon = 1e-3);

        let computed = solver
            .solve_inverse(
                (Lon::new(0.), Lat::new(80. * DEG)),
                (Lon::new(180. * DEG), Lat::new(80. * DEG)),
            )
            .unwrap();
        assert_abs_diff_eq!(computed.alpha1, Azimuth::NORTH, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.alpha2, Azimuth::SOUTH, epsilon = 1e-10);
        assert_abs_diff_eq!(computed.s, 2_233_651.715, epsilon = 1e-3);
    }
}
