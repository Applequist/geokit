use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

use crate::{
    geodesy::Ellipsoid,
    operation::{self, Operation},
};

// Polynoms for Mercator backward projection.
lazy_static! {
    static ref P2: Vec<f64> = vec![0.0, 0.5, 5.0 / 24.0, 1.0 / 12.0, 13.0 / 360.0];
    static ref P4: Vec<f64> = vec![0.0, 0.0, 7.0 / 48.0, 29.0 / 240.0, 811.0 / 11520.0];
    static ref P6: Vec<f64> = vec![0.0, 0.0, 0.0, 7.0 / 120.0, 81.0 / 1120.0];
    static ref P8: Vec<f64> = vec![0.0, 0.0, 0.0, 0.0, 4279.0 / 161280.0];
}

/// The [Mercator] map projection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mercator {
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
    pub fn new_1_sp(
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
    pub fn new_2_sp(
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

pub fn eval_polynom(coef: &[f64], x: f64) -> f64 {
    coef.iter().rev().fold(0.0, |v, a| x * v + a)
}

impl Operation for Mercator {
    fn in_dim(&self) -> usize {
        3
    }

    fn out_dim(&self) -> usize {
        3
    }

    fn fwd(&self, input: &[f64], output: &mut [f64]) -> operation::Result<()> {
        let sin_lat = input[1].sin();
        let e_sin_lat = self.e * sin_lat;
        let lat1 = FRAC_PI_4 + input[1] / 2.0;
        let r = ((1.0 - e_sin_lat) / (1.0 + e_sin_lat)).powf(self.e / 2.0);

        output[0] = self.false_easting + self.a * self.k0 * (input[0] - self.lon0);
        output[1] = self.false_northing + self.a * self.k0 * (lat1.tan() * r).ln();
        output[2] = input[2];

        Ok(())
    }

    fn bwd(&self, input: &[f64], output: &mut [f64]) -> operation::Result<()> {
        let t = ((self.false_northing - input[1]) / (self.a * self.k0)).exp();
        let xi = FRAC_PI_2 - 2.0 * t.atan();
        let sin_2xi = (2.0 * xi).sin();
        let sin_4xi = (4.0 * xi).sin();
        let sin_6xi = (6.0 * xi).sin();
        let sin_8xi = (8.0 * xi).sin();

        let e2 = self.e * self.e;

        output[0] = self.lon0 + (input[0] - self.false_easting) / (self.a * self.k0);
        output[1] = xi
            + eval_polynom(&P2, e2) * sin_2xi
            + eval_polynom(&P4, e2) * sin_4xi
            + eval_polynom(&P6, e2) * sin_6xi
            + eval_polynom(&P8, e2) * sin_8xi;
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
    pub fn new(
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

    fn fwd(&self, input: &[f64], output: &mut [f64]) -> crate::operation::Result<()> {
        output[0] = self.false_easting + self.a * (input[0] - self.lon0);
        output[1] = self.false_northing + self.a * (input[1] / 2.0 + FRAC_PI_4).tan().ln();
        output[2] = input[2];
        Ok(())
    }

    fn bwd(&self, input: &[f64], output: &mut [f64]) -> crate::operation::Result<()> {
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

    use super::eval_polynom;
    use super::Mercator;
    use super::WebMercator;

    #[test]
    fn test_eval_polynom() {
        assert_eq!(eval_polynom(&[1.0, 0.0, 2.0], 1.0), 3.0);
        assert_eq!(eval_polynom(&[1.0, 0.0, 2.0], 2.0), 9.0);
    }

    #[test]
    fn web_mercator() {
        let proj = WebMercator::new(&ellipsoid::consts::WGS84, 0.0, 0.0, 0.0, 0.0);
        let mut output = [0.0; 3];
        proj.fwd(&[-1.751147016, 0.42554246, 0.0], &mut output)
            .unwrap();
        assert_abs_diff_eq!(output[0], -11_169_055.58, epsilon = 1e-2);
        assert_abs_diff_eq!(output[1], 2_800_000.0, epsilon = 1e-2);

        proj.bwd(&[-11_169_055.58, 2_810_000.0, 0.0], &mut output)
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
        proj.fwd(
            &[120.0_f64.to_radians(), -3.0_f64.to_radians(), 0.0],
            &mut output,
        )
        .unwrap();
        assert_abs_diff_eq!(output[0], 5009726.58, epsilon = 2e-2);
        assert_abs_diff_eq!(output[1], 569150.82, epsilon = 2e-2);

        proj.bwd(&[5009726.58, 569150.82, 0.0], &mut output)
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
        proj.fwd(
            &[53.0_f64.to_radians(), 53.0_f64.to_radians(), 0.0],
            &mut output,
        )
        .unwrap();
        assert_abs_diff_eq!(165704.29, output[0], epsilon = 4e-3);
        assert_abs_diff_eq!(5171848.07, output[1], epsilon = 3e-3);

        proj.bwd(&[165704.29, 5171848.07, 0.0], &mut output)
            .unwrap();
        assert_abs_diff_eq!(output[0].to_degrees(), 53.0, epsilon = 5e-8);
        assert_abs_diff_eq!(output[1].to_degrees(), 53.0, epsilon = 4e-8);
    }
}
