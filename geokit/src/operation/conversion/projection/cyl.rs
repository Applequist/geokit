use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

use crate::{
    geodesy::Ellipsoid,
    operation::{self, Operation},
};
use crate::math::polynomial::Polynomial;

// Polynomials for Mercator backward projection.
lazy_static! {
    static ref P2: Polynomial<5> = Polynomial::new([0.0, 0.5, 5.0 / 24.0, 1.0 / 12.0, 13.0 / 360.0]);
    static ref P4: Polynomial<5> = Polynomial::new([0.0, 0.0, 7.0 / 48.0, 29.0 / 240.0, 811.0 / 11520.0]);
    static ref P6: Polynomial<5> = Polynomial::new([0.0, 0.0, 0.0, 7.0 / 120.0, 81.0 / 1120.0]);
    static ref P8: Polynomial<5> = Polynomial::new([0.0, 0.0, 0.0, 0.0, 4279.0 / 161280.0]);
}

/// The [Mercator] map projection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Mercator {
    a: f64,
    e: f64,
    e2: f64,
    lon0: f64,
    lat0: f64,
    k0: f64,
    false_easting: f64,
    false_northing: f64,
}

impl Mercator {
    /// Creates a new [Mercator] projection instance known as *Mercator (1SP)*. The projection is
    /// defined with the equator as the single standard parallel and the scale at the equator also
    /// defined. False grid coordinates are applied at the *natural origin* of the projection, the
    /// intersection of the equator and the *longitude of origin* `lon0`.
    /// Know as **EPSG:9804**.
    pub(crate) fn new_1_sp(
        ellipsoid: &Ellipsoid,
        lon0: f64,
        k0: f64,
        false_easting: f64,
        false_northing: f64,
    ) -> Self {
        Self {
            a: ellipsoid.a(),
            e: ellipsoid.e(),
            e2: ellipsoid.e_sq(),
            lon0,
            lat0: 0.0,
            k0,
            false_easting,
            false_northing,
        }
    }

    /// Creates new [Mercator] projection instance defined through the latitude of 2 parallel
    /// equidistant either side of the equator upon which the grid scale is true. False grid
    /// coordinates are applied at the *natural origin* of the projection, the intersection of the
    /// equator and the *longitude of origin* `lon0`.
    /// Known as EPSG:9805
    pub(crate) fn new_2_sp(
        ellipsoid: &Ellipsoid,
        lon0: f64,
        lat0: f64,
        false_easting: f64,
        false_northing: f64,
    ) -> Self {
        let (sin_lat1, cos_lat1) = lat0.abs().sin_cos();
        Self {
            a: ellipsoid.a(),
            e: ellipsoid.e(),
            e2: ellipsoid.e_sq(),
            lon0,
            lat0: lat0.abs(),
            k0: cos_lat1 / (1.0 - ellipsoid.e_sq() * sin_lat1 * sin_lat1).sqrt(),
            false_easting,
            false_northing,
        }
    }
}

/// Eval a polynom
// pub fn eval_polynom(coef: &[f64], x: f64) -> f64 {
//     coef.iter().rev().fold(0.0, |acc, coef_i| x * acc + coef_i)
// }

impl Operation for Mercator {
    fn in_dim(&self) -> usize {
        3
    }

    fn out_dim(&self) -> usize {
        3
    }

    fn apply_fwd(&self, input: &[f64], output: &mut [f64]) -> operation::Result<()> {
        let sin_lat = input[1].sin();
        let e_sin_lat = self.e * sin_lat;
        let lat1 = FRAC_PI_4 + input[1] / 2.0;
        let r = ((1.0 - e_sin_lat) / (1.0 + e_sin_lat)).powf(self.e / 2.0);

        output[0] = self.false_easting + self.a * self.k0 * (input[0] - self.lon0);
        output[1] = self.false_northing + self.a * self.k0 * (lat1.tan() * r).ln();
        output[2] = input[2];

        Ok(())
    }

    fn apply_bwd(&self, input: &[f64], output: &mut [f64]) -> operation::Result<()> {
        let t = ((self.false_northing - input[1]) / (self.a * self.k0)).exp();
        let xi = FRAC_PI_2 - 2.0 * t.atan();
        let sin_2xi = (2.0 * xi).sin();
        let sin_4xi = (4.0 * xi).sin();
        let sin_6xi = (6.0 * xi).sin();
        let sin_8xi = (8.0 * xi).sin();

        let e2 = self.e * self.e;

        output[0] = self.lon0 + (input[0] - self.false_easting) / (self.a * self.k0);
        output[1] = xi
            + P2.eval_at(e2) * sin_2xi
            + P4.eval_at(e2) * sin_4xi
            + P6.eval_at(e2) * sin_6xi
            + P8.eval_at(e2) * sin_8xi;
        output[2] = input[2];

        Ok(())
    }
}

/// The [TransverseMercator] projection.
/// Known as EPSG:9807.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TransverseMercator {
    e: f64,
    lon0: f64,
    k0: f64,
    false_easting: f64,
    false_northing: f64,
    upper_b: f64,
    h: [f64; 4],
    h_p: [f64; 4],
    m0: f64,
}

impl TransverseMercator {
    const UTM_ZONE_WIDTH_DEG: f64 = 6.0;
    const UTM_ZONE_COUNT: f64 = 360.0 / Self::UTM_ZONE_WIDTH_DEG;
    const UTM_ZONE_1_WEST_LIMIT_DEG: f64 = -177.0;
    const UTM_LAT0: f64 = 0.0;
    const UTM_K0: f64 = 0.9996;
    const UTM_FE: f64 = 500_000.0;
    const UTM_FN_NORTH: f64 = 0.0;
    const UTM_FN_SOUTH: f64 = 10_000_000.0;

    fn utm_zone(lon: f64) -> u8 {
        let lon_deg = lon.to_degrees();
        let delta_lon =
            ((lon_deg - Self::UTM_ZONE_1_WEST_LIMIT_DEG) / Self::UTM_ZONE_WIDTH_DEG).floor();
        let z = delta_lon - (delta_lon / Self::UTM_ZONE_COUNT).floor() * Self::UTM_ZONE_COUNT;
        z as u8
    }

    fn utm_lon0(zone: u8) -> f64 {
        let lon0_deg = Self::UTM_ZONE_1_WEST_LIMIT_DEG + (zone as f64) * Self::UTM_ZONE_WIDTH_DEG
            - Self::UTM_ZONE_WIDTH_DEG / 2.0;
        lon0_deg.to_radians()
    }

    pub(crate) fn new_utm_north(ellipsoid: &Ellipsoid, zone: u8) -> Self {
        Self::new(
            ellipsoid,
            Self::utm_lon0(zone),
            Self::UTM_LAT0,
            Self::UTM_K0,
            Self::UTM_FE,
            Self::UTM_FN_NORTH,
        )
    }

    pub(crate) fn new_utm_south(ellipsoid: &Ellipsoid, zone: u8) -> Self {
        Self::new(
            ellipsoid,
            Self::utm_lon0(zone),
            Self::UTM_LAT0,
            Self::UTM_K0,
            Self::UTM_FE,
            Self::UTM_FN_SOUTH,
        )
    }

    pub(crate) fn new(
        ellipsoid: &Ellipsoid,
        lon0: f64,
        lat0: f64,
        k0: f64,
        false_easting: f64,
        false_northing: f64,
    ) -> Self {
        let e = ellipsoid.e();
        let f = 1.0 / ellipsoid.invf();
        let n = f / (2.0 - f);
        let upper_b =
            ellipsoid.a() / (1.0 + n) * Polynomial::new([1.0, 0.0, 0.25, 0.0, 1.0 / 64.0]).eval_at(n);
        let h = [
            Polynomial::new([0.0, 0.5, -2.0 / 3.0, 5.0 / 16.0, 41.0 / 180.0]).eval_at(n),
            Polynomial::new([0.0, 0.0, 13.0 / 48.0, -3.0 / 5.0, 557.0 / 1440.0]).eval_at(n),
            Polynomial::new([0.0, 0.0, 0.0, 61.0 / 240.0, -103.0 / 140.0]).eval_at(n),
            Polynomial::new([0.0, 0.0, 0.0, 0.0, 49561.0 / 161280.0]).eval_at(n),
        ];

        let h_p = [
            Polynomial::new([0.0, 0.5, -2.0 / 3.0, 37.0 / 96.0, -1.0 / 360.0]).eval_at(n),
            Polynomial::new([0.0, 0.0, 1.0 / 48.0, 1.0 / 15.0, -437.0 / 1440.0]).eval_at(n),
            Polynomial::new([0.0, 0.0, 0.0, 17.0 / 480.0, -37.0 / 840.0]).eval_at(n),
            Polynomial::new([0.0, 0.0, 0.0, 0.0, 4397.0 / 161280.0]).eval_at(n),
        ];

        let m0 = Self::m0(lat0, upper_b, e, h);

        Self {
            e,
            lon0,
            k0,
            false_easting,
            false_northing,
            upper_b,
            h,
            h_p,
            m0,
        }
    }

    fn q(e: f64, lat: f64) -> f64 {
        lat.tan().asinh() - e * (e * lat.sin()).atanh()
    }

    fn beta(q: f64) -> f64 {
        q.sinh().atan()
    }

    fn eta_0(beta: f64, lon: f64, lon0: f64) -> f64 {
        (beta.cos() * (lon - lon0).sin()).atanh()
    }

    fn m0(lat0: f64, upper_b: f64, e: f64, h: [f64; 4]) -> f64 {
        if lat0 == 0.0 {
            0.0
        } else if lat0.abs() == FRAC_PI_2 {
            upper_b * lat0
        } else {
            let q0 = Self::q(e, lat0);
            let beta0 = Self::beta(q0);
            let xi_o0 = beta0; // From note: simplified beta0.sin().asin() to beta0;
            let xi_os = [
                xi_o0,
                h[0] * (2.0 * xi_o0).sin(),
                h[1] * (4.0 * xi_o0).sin(),
                h[2] * (6.0 * xi_o0).sin(),
                h[3] * (8.0 * xi_o0).sin(),
            ];
            let xi_o = xi_os.into_iter().sum::<f64>();
            upper_b * xi_o
        }
    }
}

impl Operation for TransverseMercator {
    fn in_dim(&self) -> usize {
        3
    }

    fn out_dim(&self) -> usize {
        3
    }

    fn apply_fwd(&self, input: &[f64], output: &mut [f64]) -> operation::Result<()> {
        let q = Self::q(self.e, input[1]);
        let beta = Self::beta(q);
        let eta_0 = Self::eta_0(beta, input[0], self.lon0);
        let xi_0 = (beta.sin() * eta_0.cosh()).asin();
        let xis = [
            xi_0,
            self.h[0] * (2.0 * xi_0).sin() * (2.0 * eta_0).cosh(),
            self.h[1] * (4.0 * xi_0).sin() * (4.0 * eta_0).cosh(),
            self.h[2] * (6.0 * xi_0).sin() * (6.0 * eta_0).cosh(),
            self.h[3] * (8.0 * xi_0).sin() * (8.0 * eta_0).cosh(),
        ];
        let xi = xis.into_iter().sum::<f64>();
        let etas = [
            self.h[0] * (2.0 * xi_0).cos() * (2.0 * eta_0).sinh(),
            self.h[1] * (4.0 * xi_0).cos() * (4.0 * eta_0).sinh(),
            self.h[2] * (6.0 * xi_0).cos() * (6.0 * eta_0).sinh(),
            self.h[3] * (8.0 * xi_0).cos() * (8.0 * eta_0).sinh(),
        ];
        let eta = eta_0 + etas.into_iter().sum::<f64>();

        output[0] = self.false_easting + self.k0 * self.upper_b * eta;
        output[1] = self.false_northing + self.k0 * (self.upper_b * xi - self.m0);
        output[2] = input[2];

        Ok(())
    }

    fn apply_bwd(&self, input: &[f64], output: &mut [f64]) -> operation::Result<()> {
        let eta_p = (input[0] - self.false_easting) / (self.upper_b * self.k0);
        let xi_p = (input[1] - self.false_northing + self.k0 * self.m0) / (self.upper_b * self.k0);

        let xi_ps = [
            self.h_p[0] * (2.0 * xi_p).sin() * (2.0 * eta_p).cosh(),
            self.h_p[1] * (4.0 * xi_p).sin() * (4.0 * eta_p).cosh(),
            self.h_p[2] * (6.0 * xi_p).sin() * (6.0 * eta_p).cosh(),
            self.h_p[3] * (8.0 * xi_p).sin() * (8.0 * eta_p).cosh(),
        ];
        let xi_p0 = xi_p - xi_ps.into_iter().sum::<f64>();

        let eta_ps = [
            self.h_p[0] * (2.0 * xi_p).cos() * (2.0 * eta_p).sinh(),
            self.h_p[1] * (4.0 * xi_p).cos() * (4.0 * eta_p).sinh(),
            self.h_p[2] * (6.0 * xi_p).cos() * (6.0 * eta_p).sinh(),
            self.h_p[3] * (8.0 * xi_p).cos() * (8.0 * eta_p).sinh(),
        ];
        let eta_p0 = eta_p - eta_ps.into_iter().sum::<f64>();

        let beta_p = (xi_p0.sin() / eta_p0.cosh()).asin();
        let q_p = beta_p.tan().asinh();

        fn iter_fn(e: f64, q_p: f64, q_pp: f64) -> f64 {
            q_p + (e * (e * q_pp.tanh()).atanh())
        }

        let mut q_pp = iter_fn(self.e, q_p, q_p);
        loop {
            let old_q_pp = q_pp;
            q_pp = iter_fn(self.e, q_p, q_pp);
            if (q_pp - old_q_pp).abs() < 1e-10 {
                break;
            }
        }

        output[0] = self.lon0 + (eta_p0.tanh() / beta_p.cos()).asin();
        output[1] = q_pp.sinh().atan();
        output[2] = input[2];

        Ok(())
    }
}

/// [WebMercator](epsg:1024) also known as 'Pseudo-Mercator' is a projection method
/// used by some popular web mapping and visualisation applications.
/// Strictly speaking the name is misleading as it is **NOT** a Mercator projection.
#[derive(Debug, Clone, PartialEq)]
pub struct WebMercator {
    a: f64,
    /// longitude of *natural origin* in radians.
    lon0: f64,
    /// latitude of *natural origin* in radians.
    lat0: f64,
    /// False easting in meters.
    false_easting: f64,
    /// False norhting in meters.
    false_northing: f64,
}

impl WebMercator {
    /// Creates a new [WebMercator] projection instance.
    pub(crate) fn new(
        ellipsoid: &Ellipsoid,
        lon0: f64,
        lat0: f64,
        false_easting: f64,
        false_northing: f64,
    ) -> Self {
        Self {
            a: ellipsoid.a(),
            lon0,
            lat0,
            false_easting,
            false_northing,
        }
    }
}

impl Operation for WebMercator {
    fn in_dim(&self) -> usize {
        3
    }

    fn out_dim(&self) -> usize {
        3
    }

    fn apply_fwd(&self, input: &[f64], output: &mut [f64]) -> operation::Result<()> {
        output[0] = self.false_easting + self.a * (input[0] - self.lon0);
        output[1] = self.false_northing + self.a * (input[1] / 2.0 + FRAC_PI_4).tan().ln();
        output[2] = input[2];
        Ok(())
    }

    fn apply_bwd(&self, input: &[f64], output: &mut [f64]) -> operation::Result<()> {
        let d = (self.false_northing - input[1]) / self.a;
        output[0] = self.lon0 + (input[0] - self.false_easting) / self.a;
        output[1] = FRAC_PI_2 - 2.0 * d.exp().atan();
        output[2] = input[2];
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use crate::{geodesy::ellipsoid, operation::Operation};

    use super::Mercator;
    use super::TransverseMercator;
    use super::WebMercator;

    #[test]
    fn web_mercator() {
        let proj = WebMercator::new(&ellipsoid::consts::WGS84, 0.0, 0.0, 0.0, 0.0);
        let mut output = [0.0; 3];
        proj.apply_fwd(&[-1.751147016, 0.42554246, 0.0], &mut output)
            .unwrap();
        assert_abs_diff_eq!(output[0], -11_169_055.58, epsilon = 1e-2);
        assert_abs_diff_eq!(output[1], 2_800_000.0, epsilon = 1e-2);

        proj.apply_bwd(&[-11_169_055.58, 2_810_000.0, 0.0], &mut output)
            .unwrap();
        assert_abs_diff_eq!(output[0], -1.751147016, epsilon = 1e-9);
        assert_abs_diff_eq!(output[1], 0.426970023, epsilon = 1e-9);
    }

    #[test]
    fn mercator_1_sp() {
        // From EPSG Guidance Note 7 Part 2 3.5.1 Mercator
        let proj = Mercator::new_1_sp(
            &ellipsoid::consts::BESSEL,
            110.0_f64.to_radians(),
            0.997,
            3_900_000.0,
            900_000.0,
        );
        let mut output = [0.0; 3];
        proj.apply_fwd(
            &[120.0_f64.to_radians(), -3.0_f64.to_radians(), 0.0],
            &mut output,
        )
        .unwrap();
        assert_abs_diff_eq!(output[0], 5009726.58, epsilon = 2e-2);
        assert_abs_diff_eq!(output[1], 569150.82, epsilon = 2e-2);

        proj.apply_bwd(&[5009726.58, 569150.82, 0.0], &mut output)
            .unwrap();
        assert_abs_diff_eq!(output[0].to_degrees(), 120.0, epsilon = 3e-8);
        assert_abs_diff_eq!(output[1].to_degrees(), -3.0, epsilon = 3e-8);
    }

    #[test]
    fn mercator_2_sp() {
        // From EPSG Guidance Note 7 Part 2 3.5.1 Mercator
        let proj = Mercator::new_2_sp(
            &ellipsoid::consts::KRASS,
            51.0_f64.to_radians(),
            42.0_f64.to_radians(),
            0.0,
            0.0,
        );
        let mut output = [0.0; 3];
        proj.apply_fwd(
            &[53.0_f64.to_radians(), 53.0_f64.to_radians(), 0.0],
            &mut output,
        )
        .unwrap();
        assert_abs_diff_eq!(165704.29, output[0], epsilon = 4e-3);
        assert_abs_diff_eq!(5171848.07, output[1], epsilon = 3e-3);

        proj.apply_bwd(&[165704.29, 5171848.07, 0.0], &mut output)
            .unwrap();
        assert_abs_diff_eq!(output[0].to_degrees(), 53.0, epsilon = 5e-8);
        assert_abs_diff_eq!(output[1].to_degrees(), 53.0, epsilon = 4e-8);
    }

    #[test]
    fn transverse_mercator() {
        // From EPSG Guidance Note 7 part 2
        let proj = TransverseMercator::new(
            &ellipsoid::consts::AIRY,
            -2.0_f64.to_radians(), // 2 W
            49.0_f64.to_radians(), // 49 N
            0.9996012717,
            400_000.0,
            -100_000.0,
        );

        let enh = proj
            .fwd_new(&[0.5_f64.to_radians(), 50.5_f64.to_radians(), 0.0])
            .unwrap();
        assert_abs_diff_eq!(enh[0], 577_274.99, epsilon = 1e-2);
        assert_abs_diff_eq!(enh[1], 69_740.5, epsilon = 1e-2);

        let llh = proj.bwd_new(&[577_274.99, 69_740.50, 0.0]).unwrap();
        assert_abs_diff_eq!(llh[0], 0.5_f64.to_radians(), epsilon = 2e-9);
        assert_abs_diff_eq!(llh[1], 50.5_f64.to_radians(), epsilon = 2e-9);
    }

    #[test]
    fn transverse_mercator_bound() {
        let proj = TransverseMercator::new(&ellipsoid::consts::WGS84, 0.0, 0.0, 0.9996, 0.0, 0.0);

        let north_pole = proj.fwd_new(&[0., 90.0_f64.to_radians(), 0.0]).unwrap();
        println!("north_pole = {north_pole:?}");
        let east_north = proj
            .fwd_new(&[90f64.to_radians(), 10.0f64.to_radians(), 0.0])
            .unwrap();
        println!("90E 10N = {east_north:?}");
        let west_north = proj
            .fwd_new(&[-90f64.to_radians(), 10f64.to_radians(), 0.0])
            .unwrap();
        println!("90W 10N = {west_north:?}");

        let south_pole = proj.fwd_new(&[0., -90.0_f64.to_radians(), 0.0]).unwrap();
        println!("south_pole = {south_pole:?}");

        let med_90: Vec<f64> = (-9..=9)
            .flat_map(|i| vec![90.0_f64.to_radians(), (i as f64 * 10.0).to_radians(), 0.0])
            .collect();
        let _proj_med_90 = proj.fwd_new(&med_90).unwrap();
    }

    #[test]
    fn transverse_mercator_roundtrip() {
        let proj = TransverseMercator::new(&ellipsoid::consts::WGS84, 0.0, 0.0, 0.9996, 0.0, 0.0);

        let north_pole = vec![10f64.to_radians(), 90.0_f64.to_radians(), 0.0];
        let north_pole_enh = proj.fwd_new(&north_pole).unwrap();
        println!("north_pole (enh) = {north_pole_enh:?}");

        let north_pole_llh: Vec<f64> = proj.bwd_new(&north_pole_enh).unwrap();
        println!("north_pole (llh)= {north_pole_llh:?}");

        assert_abs_diff_eq!(&north_pole_llh[..], &north_pole[..], epsilon = 1e-9);

        let south_pole = vec![-10f64.to_radians(), -90.0_f64.to_radians(), 0.0];
        let south_pole_enh = proj.fwd_new(&south_pole).unwrap();
        println!("south_pole (enh) = {south_pole_enh:?}");

        let south_pole_llh: Vec<f64> = proj.bwd_new(&south_pole_enh).unwrap();
        println!("south_pole (llh)= {south_pole_llh:?}");

        assert_abs_diff_eq!(&south_pole_llh[..], &south_pole[..], epsilon = 1e-9);
    }
}
