use crate::{
    crs::{GeocentricAxes, GeodeticAxes, ProjectedAxes},
    geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian},
};

use super::{Operation, Result};

type ToOrd = (usize, f64);

/// [Normalization] is used to normalize(fwd)/denormalize(bwd) coordinates by:
/// - reordering coordinates
/// - converting coordinates quantity
/// - zeroing input coordinates in excess or output coordinates in excess.
///
/// Normalized geocentric coordinates are (x, y, z) measured in metres.
/// Normalized geographic coordinates are using the `GeodeticAxes::EastNorthUp(1.0, 1.0)` axes..
#[derive(Debug, Clone)]
pub struct Normalization(Vec<ToOrd>);

impl Normalization {
    pub fn identity(dim: usize) -> Self {
        let to_ord = (0..dim).map(|ord| (ord, 1.0)).collect::<Vec<_>>();
        Normalization(to_ord)
    }
}

impl From<GeocentricAxes> for Normalization {
    fn from(_value: GeocentricAxes) -> Self {
        Normalization::identity(3)
    }
}

impl From<GeodeticAxes> for Normalization {
    fn from(value: GeodeticAxes) -> Self {
        let to_ord = match value {
            GeodeticAxes::EastNorthUp {
                angle_unit,
                height_unit,
            } => vec![
                (0, angle_unit.rad_per_unit()),
                (1, angle_unit.rad_per_unit()),
                (2, height_unit),
            ],
            GeodeticAxes::EastNorth { angle_unit } => {
                vec![
                    (0, angle_unit.rad_per_unit()),
                    (1, angle_unit.rad_per_unit()),
                ]
            }
            GeodeticAxes::NorthEastUp {
                angle_unit,
                height_unit,
            } => vec![
                (1, angle_unit.rad_per_unit()),
                (0, angle_unit.rad_per_unit()),
                (2, height_unit),
            ],
            GeodeticAxes::NorthEast { angle_unit } => {
                vec![
                    (1, angle_unit.rad_per_unit()),
                    (0, angle_unit.rad_per_unit()),
                ]
            }
            GeodeticAxes::NorthWest { angle_unit } => {
                vec![
                    (1, -angle_unit.rad_per_unit()),
                    (0, angle_unit.rad_per_unit()),
                ]
            }
        };
        Normalization(to_ord)
    }
}

impl From<ProjectedAxes> for Normalization {
    fn from(value: ProjectedAxes) -> Self {
        let to_ord = match value {
            ProjectedAxes::EastNorthUp {
                horiz_unit,
                height_unit,
            } => vec![(0, horiz_unit), (1, horiz_unit), (2, height_unit)],
            ProjectedAxes::EastNorth { horiz_unit } => {
                vec![(0, horiz_unit), (1, horiz_unit)]
            }
        };
        Normalization(to_ord)
    }
}

impl Operation for Normalization {
    fn in_dim(&self) -> usize {
        self.0.len()
    }

    fn out_dim(&self) -> usize {
        3
    }

    fn apply_fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        output.copy_from_slice(&[0.; 3]);
        for (nix, (ix, f)) in self.0.iter().enumerate() {
            output[nix] = input[*ix] * f;
        }
        Ok(())
    }

    fn apply_bwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        for (nix, (ix, f)) in self.0.iter().enumerate() {
            output[*ix] = input[nix] / f;
        }
        Ok(())
    }
}

/// [GeogToGeoc] converts **normalized geographic coordinates** to **normalized geocentric coordinates**
/// between a [GeographicCrs] and a [GeocentricCrs] that share the same [GeodeticDatum].
///
/// As mentioned in the EPSG guidance note 7 part 2, paragraph 4.1.1, this transformation
/// first transform to/from non-greenwich base geographic coordinates into greenwich-base ones
/// before/after transformation to/from geocentric coordinates.
///
/// This is epsg:9602.
#[derive(Debug, Clone)]
pub struct GeogToGeoc {
    ellipsoid: Ellipsoid,
    prime_meridian: PrimeMeridian,
}

impl GeogToGeoc {
    pub fn new(datum: &GeodeticDatum) -> Self {
        Self {
            ellipsoid: datum.ellipsoid().clone(),
            prime_meridian: datum.prime_meridian().clone(),
        }
    }

    /// Convert **normalized geodetic coordinates** (lon in rad, lat in rad, height in meters)
    /// into **normalized geocentric coordinates** (x, y, z) all in meters.
    fn llh_to_xyz(&self, llh: &[f64], xyz: &mut [f64]) {
        let lon = llh[0];
        let lat = llh[1];
        let h = llh[2];

        let v = self.ellipsoid.prime_vertical_radius(lat);
        let (sin_lon, cos_lon) = lon.sin_cos();
        let (sin_lat, cos_lat) = lat.sin_cos();

        xyz[0] = (v + h) * cos_lat * cos_lon;
        xyz[1] = (v + h) * cos_lat * sin_lon;
        xyz[2] = (v * (1.0 - self.ellipsoid.e_sq()) + h) * sin_lat;
    }

    /// Convert **normalized geocentric coordinates** (x, y, z) in meters
    /// into **normalized geodetic coordinates** (lon in rad, lat in rad, height in meters)
    /// using Heiskanen and Moritz iterative method.
    fn xyz_to_llh(&self, xyz: &[f64], llh: &mut [f64]) {
        let x = xyz[0];
        let y = xyz[1];
        let z = xyz[2];

        let a2 = self.ellipsoid.a_sq();
        let b2 = self.ellipsoid.b_sq();
        let e2 = self.ellipsoid.e_sq();

        let lon = y.atan2(x);

        let p = x.hypot(y);
        let mut lat = z.atan2(p * (1.0 - e2));
        let (sin_lat, cos_lat) = lat.sin_cos();
        let n = a2 / (a2 * cos_lat * cos_lat + b2 * sin_lat * sin_lat).sqrt();
        let mut h = p / cos_lat - n;
        loop {
            let next_lat = z.atan2(p * (1.0 - e2 * n / (n + h)));
            let (sin_nlat, cos_nlat) = next_lat.sin_cos();
            let next_n = a2 / ((a2 * cos_nlat * cos_nlat) + b2 * sin_nlat * sin_nlat).sqrt();
            let next_h = p / cos_nlat - next_n;
            let delta_lat = (lat - next_lat).abs();
            let delta_h = (h - next_h).abs();
            lat = next_lat;
            h = next_h;
            if delta_lat < 0.5e-5 && delta_h < 0.5e-3 {
                break;
            }
        }

        llh[0] = lon;
        llh[1] = lat;
        llh[2] = h;
    }

    /// Convert a **normalized longitude** with the prime meridian as origin into
    /// a **normalized longitude** with the Greenwich prime meridian as origin.
    #[inline]
    fn convert_lon_to_gw(&self, lon: f64) -> f64 {
        // FIX: What if we cross the antimeridian?
        lon + self.prime_meridian.lon()
    }

    /// Convert a **normalized longitude** with the Greenwich prime meridian as origin into
    /// a **normalized longitude** with this prime meridian as origin.
    #[inline]
    fn convert_lon_from_gw(&self, lon: f64) -> f64 {
        // FIX: What if we cross the antimeridian?
        lon - self.prime_meridian.lon()
    }
}

impl Operation for GeogToGeoc {
    fn in_dim(&self) -> usize {
        3
    }

    fn out_dim(&self) -> usize {
        3
    }

    /// Convert geographic coordinates to geocentric coordinates.
    fn apply_fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        let mut llh_cpy = [0.0; 3];
        llh_cpy.copy_from_slice(input);
        llh_cpy[0] = self.convert_lon_to_gw(input[0]);
        self.llh_to_xyz(&llh_cpy, output);
        Ok(())
    }

    /// Convert geocentric coordinates to geographic coordinates.
    fn apply_bwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.xyz_to_llh(input, output);
        output[0] = self.convert_lon_from_gw(output[0]);
        Ok(())
    }
}

pub mod projection;

#[cfg(test)]
mod tests {
    use crate::crs::GeodeticAxes;
    use crate::operation::conversion::Normalization;
    use crate::operation::Operation;
    use crate::quantity::angle::units::Deg;

    #[test]
    fn normalization() {
        let latlondeg = GeodeticAxes::NorthWest {
            angle_unit: Deg::UNIT,
        };
        let t = Normalization::from(latlondeg);
        assert_eq!(t.in_dim(), 2);
        assert_eq!(t.out_dim(), 3);

        let i = [10.0, 110.0];
        let mut o = [0.0; 3];
        t.apply_fwd(&i, &mut o).unwrap();
        assert_eq!(o, [-110.0f64.to_radians(), 10.0f64.to_radians(), 0.0]);
    }

    #[test]
    #[ignore = "Not Implemented"]
    fn geog_to_geoc_fwd() {
        unimplemented!()
    }

    #[test]
    #[ignore = "Not Implemented"]
    fn geog_to_geoc_bwd() {
        unimplemented!()
    }
}
