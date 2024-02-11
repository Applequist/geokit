use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

use crate::{geodesy::Ellipsoid, operation::DynOperation};

/// [WebMercator](epsg:1024) also known as 'Pseudo-Mercator' is a projection method
/// used by some popular web mapping and visualisation applications.
/// Strictly speaking the name is misleading as it is **NOT** a Mercator projection.
#[derive(Debug, Clone, PartialEq)]
pub struct WebMercator {
    ellipsoid: Ellipsoid,
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
        ellipsoid: Ellipsoid,
        lon0: f64,
        lat0: f64,
        false_easting: f64,
        false_northing: f64,
    ) -> Self {
        Self {
            ellipsoid,
            lon0,
            lat0,
            false_easting,
            false_northing,
        }
    }
}

impl DynOperation for WebMercator {
    fn fwd_in_dim(&self) -> usize {
        3
    }

    fn fwd_out_dim(&self) -> usize {
        3
    }

    fn fwd(&self, input: &[f64], output: &mut [f64]) -> crate::operation::Result<()> {
        output[0] = self.false_easting + self.ellipsoid.a() * (input[0] - self.lon0);
        output[1] =
            self.false_northing + self.ellipsoid.a() * (input[1] / 2.0 + FRAC_PI_4).tan().ln();
        output[2] = input[2];
        Ok(())
    }

    fn is_invertible(&self) -> bool {
        true
    }

    fn bwd(&self, input: &[f64], output: &mut [f64]) -> crate::operation::Result<()> {
        let d = (self.false_northing - input[1]) / self.ellipsoid.a();
        output[0] = self.lon0 + (input[0] - self.false_easting) / self.ellipsoid.a();
        output[1] = FRAC_PI_2 - 2.0 * d.exp().atan();
        output[2] = input[2];
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use crate::{geodesy::ellipsoid, operation::DynOperation};

    use super::WebMercator;

    #[test]
    fn web_mercator_fwd() {
        let ellipsoid = ellipsoid::consts::WGS84;
        let proj = WebMercator::new(ellipsoid, 0.0, 0.0, 0.0, 0.0);
        let mut output = [0.0; 3];
        proj.fwd(&[-1.751147016, 0.42554246, 0.0], &mut output)
            .unwrap();
        assert_abs_diff_eq!(output[0], -11_169_055.58, epsilon = 1e-2);
        assert_abs_diff_eq!(output[1], 2_800_000.0, epsilon = 1e-2);
    }

    #[test]
    fn web_mercator_bwd() {
        let ellipsoid = ellipsoid::consts::WGS84;
        let proj = WebMercator::new(ellipsoid, 0.0, 0.0, 0.0, 0.0);
        let mut output = [0.0; 3];
        proj.bwd(&[-11_169_055.58, 2_810_000.0, 0.0], &mut output)
            .unwrap();
        assert_abs_diff_eq!(output[0], -1.751147016, epsilon = 1e-9);
        assert_abs_diff_eq!(output[1], 0.426970023, epsilon = 1e-9);
    }
}
