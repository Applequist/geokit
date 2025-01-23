use super::{Operation, Result};
use crate::{
    cs::{
        cartesian::{GeocentricAxes, ProjectedAxes},
        geodetic::GeodeticAxes,
    },
    geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian},
    math::Float,
};

type ToOrd = (usize, Float);

/// [Normalization] is used to normalize(fwd)/denormalize(bwd) coordinates by:
/// - reordering coordinates
/// - converting coordinates quantity
/// - zeroing input coordinates in excess or output coordinates in excess.
///
/// Normalized geocentric coordinates are (x, y, z) measured in metres.
/// Normalized geographic coordinates are using the `GeodeticAxes::EastNorthUp(RAD, M)` axes..
#[deprecated(note = "Use the axes")]
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
                (2, height_unit.m_per_unit()),
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
                (2, height_unit.m_per_unit()),
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
            } => vec![
                (0, horiz_unit.m_per_unit()),
                (1, horiz_unit.m_per_unit()),
                (2, height_unit.m_per_unit()),
            ],
            ProjectedAxes::EastNorth { horiz_unit } => {
                vec![(0, horiz_unit.m_per_unit()), (1, horiz_unit.m_per_unit())]
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

    fn apply_fwd(&self, input: &[Float], output: &mut [Float]) -> Result<()> {
        output.copy_from_slice(&[0.; 3]);
        for (nix, (ix, f)) in self.0.iter().enumerate() {
            output[nix] = input[*ix] * f;
        }
        Ok(())
    }

    fn apply_bwd(&self, input: &[Float], output: &mut [Float]) -> Result<()> {
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

#[deprecated(note = "Use the GeodeticDatum")]
impl GeogToGeoc {
    pub fn new(datum: &GeodeticDatum) -> Self {
        Self {
            ellipsoid: datum.ellipsoid().clone(),
            prime_meridian: datum.prime_meridian().clone(),
        }
    }

    ///// Convert **normalized geodetic coordinates** (lon in rad, lat in rad, height in meters)
    ///// into **normalized geocentric coordinates** (x, y, z) all in meters.
    #[deprecated(note = "llh_to_xyz is now done in GeodeticDatum")]
    fn llh_to_xyz(&self, llh: &[Float], xyz: &mut [Float]) {
        panic!("llh_to_xyz is now done in GeodeticDatum");
    }

    /// Convert **normalized geocentric coordinates** (x, y, z) in meters
    /// into **normalized geodetic coordinates** (lon in rad, lat in rad, height in meters)
    /// using Heiskanen and Moritz iterative method.
    #[deprecated(note = "xyz_to_llh is now done in GeodeticDatum")]
    fn xyz_to_llh(&self, xyz: &[Float], llh: &mut [Float]) {
        panic!("xyz_to_llh is now done in GeodeticDatum");
    }

    /// Convert a **normalized longitude** with the prime meridian as origin into
    /// a **normalized longitude** with the Greenwich prime meridian as origin.
    #[deprecated(note = "done in llh_to_xyz")]
    fn convert_lon_to_gw(&self, lon: Float) -> Float {
        panic!("conversion to Greenwich-based longitude is done in llh_to_xyz");
    }

    /// Convert a **normalized longitude** with the Greenwich prime meridian as origin into
    /// a **normalized longitude** with this prime meridian as origin.
    #[deprecated(note = "done in xyz_to_llh")]
    fn convert_lon_from_gw(&self, lon: Float) -> Float {
        panic!("conversion from Greenwich-based longitude is done in llh_to_xyz");
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
    fn apply_fwd(&self, input: &[Float], output: &mut [Float]) -> Result<()> {
        let mut llh_cpy = [0.0; 3];
        llh_cpy.copy_from_slice(input);
        llh_cpy[0] = self.convert_lon_to_gw(input[0]);
        self.llh_to_xyz(&llh_cpy, output);
        Ok(())
    }

    /// Convert geocentric coordinates to geographic coordinates.
    fn apply_bwd(&self, input: &[Float], output: &mut [Float]) -> Result<()> {
        self.xyz_to_llh(input, output);
        output[0] = self.convert_lon_from_gw(output[0]);
        Ok(())
    }
}

pub mod projection;

#[cfg(test)]
mod tests {
    use crate::cs::geodetic::GeodeticAxes;
    use crate::operation::conversion::Normalization;
    use crate::operation::Operation;
    use crate::units::angle::DEG;

    #[test]
    fn normalization() {
        let latlondeg = GeodeticAxes::NorthWest { angle_unit: DEG };
        let t = Normalization::from(latlondeg);
        assert_eq!(t.in_dim(), 2);
        assert_eq!(t.out_dim(), 3);

        let i = [10.0, 110.0];
        let mut o = [0.0; 3];
        t.apply_fwd(&i, &mut o).unwrap();
        // FIX: remove the need to use f64 suffixes
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
