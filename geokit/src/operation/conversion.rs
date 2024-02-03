use crate::{
    crs::geographic::GeodeticAxes,
    geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian},
};

use super::DynOperation;
use super::Result;

pub type ToOrd = (usize, f64);

/// [Normalization] is used to normalize(fwd)/denormalize(bwd) coordinates by:
/// - reordering coordinates
/// - converting coordinates units
/// - zeroing input coordinates in excess or output coordinates in excess.
///
/// Normalized geocentric coordinates are (x, y, z) measured in metres.
/// Normalized geographic coordinates are (lon, lat(, height)) measured in radians (and metres).
#[derive(Debug, Clone)]
pub struct Normalization(Vec<ToOrd>);

impl From<&GeodeticAxes> for Normalization {
    fn from(value: &GeodeticAxes) -> Self {
        let to_ord = match *value {
            GeodeticAxes::EastNorthUp {
                angle_unit,
                height_unit,
            } => vec![(0, angle_unit), (1, angle_unit), (2, height_unit)],
            GeodeticAxes::EastNorth { angle_unit } => vec![(0, angle_unit), (1, angle_unit)],
            GeodeticAxes::NorthEastUp {
                angle_unit,
                height_unit,
            } => vec![(1, angle_unit), (0, angle_unit), (2, height_unit)],
            GeodeticAxes::NorthEast { angle_unit } => vec![(1, angle_unit), (0, angle_unit)],
            GeodeticAxes::NorthWest { angle_unit } => vec![(1, -angle_unit), (0, angle_unit)],
        };
        Normalization(to_ord)
    }
}

impl DynOperation for Normalization {
    fn fwd_in_dim(&self) -> usize {
        self.0.len()
    }

    fn fwd_out_dim(&self) -> usize {
        3
    }

    fn fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        output.copy_from_slice(&[0.; 3]);
        for (nix, (ix, f)) in self.0.iter().enumerate() {
            output[nix] = input[*ix] * f;
        }
        Ok(())
    }

    fn is_invertible(&self) -> bool {
        true
    }

    fn bwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        for (nix, (ix, f)) in self.0.iter().enumerate() {
            output[*ix] = input[nix] / f;
        }
        Ok(())
    }
}

/// `GeogToGeoc` converts **normalized geographic coordinates** to **normalized geocentric coordinates**
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
}

impl DynOperation for GeogToGeoc {
    fn fwd_in_dim(&self) -> usize {
        3
    }

    fn fwd_out_dim(&self) -> usize {
        3
    }

    /// Convert geographic coordinates to geocentric coordinates.
    fn fwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        let mut llh_cpy = [0.0; 3];
        llh_cpy.copy_from_slice(input);
        llh_cpy[0] = self.prime_meridian.convert_lon_to_gw(input[0]);
        self.ellipsoid.llh_to_xyz(&llh_cpy, output);
        Ok(())
    }

    /// Convert geocentric coordinates to geographic coordinates.
    fn bwd(&self, input: &[f64], output: &mut [f64]) -> Result<()> {
        self.ellipsoid.xyz_to_llh(input, output);
        output[0] = self.prime_meridian.convert_lon_from_gw(output[0]);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::crs::{GeodeticAxes, GeographicCrs};
    use crate::operation::conversion::Normalization;
    use crate::operation::DynOperation;

    use crate::geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian};

    #[test]
    fn normalization() {
        let latlondeg = GeographicCrs::new(
            "WGS84",
            GeodeticDatum::new(
                "WGS84",
                Ellipsoid::from_ab("WGS84", 1., 0.99),
                PrimeMeridian::new("Greenwich", 0.0),
                None,
            ),
            GeodeticAxes::NorthWest {
                angle_unit: 1.0_f64.to_radians(),
            },
        );
        let t: Normalization = latlondeg.axes().into();
        assert_eq!(t.fwd_in_dim(), 2);
        assert_eq!(t.fwd_out_dim(), 3);

        let i = [10.0, 110.0];
        let mut o = [0.0; 3];
        t.fwd(&i, &mut o).unwrap();
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
