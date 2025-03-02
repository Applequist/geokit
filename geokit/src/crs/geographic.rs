use super::Crs;
use crate::{
    cs::{
        cartesian::XYZ,
        geodetic::{GeodeticAxes, GeodeticErrors},
    },
    geodesy::GeodeticDatum,
    math::fp::Float,
    transformations::{ToXYZTransformation, ToXYZTransformationProvider, TransformationError},
};
use smol_str::SmolStr;

/// A [GeographicCrs] is a **2D or 3D geodetic coordinates reference system** in which
/// coordinates are made up of longitude, latitude and optionally ellipsoidal height in various order, direction
/// and quantity.
#[derive(Debug, Clone, PartialEq)]
pub struct GeographicCrs {
    pub id: SmolStr,
    pub datum: GeodeticDatum,
    pub axes: GeodeticAxes,
}

impl Crs for GeographicCrs {
    fn id(&self) -> &str {
        &self.id
    }
}

impl ToXYZTransformationProvider for GeographicCrs {
    fn to_xyz_transformation<'a>(&self) -> Box<dyn ToXYZTransformation + 'a> {
        Box::new(self.clone())
    }
}

impl ToXYZTransformation for GeographicCrs {
    fn to_xyz(&self, coords: &[Float]) -> Result<XYZ, TransformationError> {
        let llh = self.axes.normalize(coords);
        Ok(self.datum.llh_to_xyz(llh))
    }

    fn from_xyz(&self, xyz: XYZ, coords: &mut [Float]) -> Result<(), TransformationError> {
        let llh = self.datum.xyz_to_llh(xyz);
        Ok(self.axes.denormalize(llh, coords))
    }
}

#[cfg(test)]
mod tests {
    use super::GeographicCrs;
    use crate::{
        geodesy::{ellipsoid::consts::WGS84, prime_meridian::consts::GREENWICH, GeodeticDatum},
        units::{angle::DEG, length::M},
    };
    use std::any::{Any, TypeId};

    #[test]
    fn test_any() {
        let crs: &dyn Any = &GeographicCrs {
            id: "WGS84".into(),
            datum: GeodeticDatum::new("WGS84", WGS84, GREENWICH),
            axes: crate::cs::geodetic::GeodeticAxes::EastNorthUp {
                angle_unit: DEG,
                height_unit: M,
            },
        };

        assert_eq!(crs.type_id(), TypeId::of::<GeographicCrs>());
        assert!(crs.is::<GeographicCrs>());
    }
}
