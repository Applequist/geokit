use super::Crs;
use crate::cs::cartesian::geocentric::XYZ;
use crate::cs::cartesian::projected::{ProjectedAxes, ProjectedTolerance};
use crate::quantities::length::Length;
use crate::transformations::{
    ToXYZTransformation, ToXYZTransformationProvider, TransformationError,
};
use crate::units::length::M;
use crate::{
    geodesy::GeodeticDatum,
    math::fp::Float,
    projections::{Projection, ProjectionError, ProjectionSpec},
};
use approx::AbsDiffEq;
use smol_str::SmolStr;

/// A [ProjectedCrs] is a **2D/3D cartesian coordinates reference system** derived from a
/// 2D or 3D [Geographic] Crs using a projection.
/// Coordinates are made up of easting, northing and optionally ellipsoidal height in various order,
/// direction and quantity.
#[derive(Debug, Clone, PartialEq)]
pub struct ProjectedCrs {
    pub id: SmolStr,
    pub datum: GeodeticDatum,
    pub axes: ProjectedAxes,
    pub projection: ProjectionSpec,
}

impl Crs for ProjectedCrs {
    type Tolerance = ProjectedTolerance;

    fn id(&self) -> &str {
        &self.id
    }

    fn dim(&self) -> usize {
        self.axes.dim()
    }

    fn approx_eq(&self, a: &[Float], b: &[Float], err: Self::Tolerance) -> bool {
        let (horiz_unit, height_unit) = match self.axes {
            ProjectedAxes::EastNorthUp {
                horiz_unit,
                height_unit,
            } => (horiz_unit, height_unit),
            ProjectedAxes::EastNorth { horiz_unit } => (horiz_unit, M),
        };
        let converted_tol = err.convert_to(horiz_unit, height_unit);

        match self.axes {
            ProjectedAxes::EastNorthUp {
                horiz_unit: _,
                height_unit: _,
            } => {
                a[0].abs_diff_eq(&b[0], converted_tol.horiz.0)
                    && a[1].abs_diff_eq(&b[1], converted_tol.horiz.0)
                    && a[2].abs_diff_eq(&b[2], converted_tol.height.0)
            }
            ProjectedAxes::EastNorth { horiz_unit: _ } => {
                a[0].abs_diff_eq(&b[0], converted_tol.horiz.0)
                    && a[1].abs_diff_eq(&b[1], converted_tol.horiz.0)
            }
        }
    }

    fn dist(&self, a: &[Float], b: &[Float]) -> Result<Length, &'static str> {
        let na = self.axes.normalize(a);
        let nb = self.axes.normalize(b);
        Ok(na.dist_to(nb))
    }
}

impl ToXYZTransformationProvider for ProjectedCrs {
    fn to_xyz_transformation<'a>(&self) -> Box<dyn ToXYZTransformation + 'a> {
        Box::new(ProjectedToXYZ {
            datum: self.datum.clone(),
            axes: self.axes,
            projection: self.projection.applied_to(self.datum.ellipsoid()),
        })
    }
}

struct ProjectedToXYZ<'a> {
    pub datum: GeodeticDatum,
    pub axes: ProjectedAxes,
    pub projection: Box<dyn Projection + 'a>,
}

impl From<ProjectionError> for TransformationError {
    fn from(_err: ProjectionError) -> Self {
        // TODO: do proper error conversion
        Self::OutOfBounds
    }
}

impl<'a> ToXYZTransformation for ProjectedToXYZ<'a> {
    fn to_xyz(&self, coords: &[Float]) -> Result<XYZ, TransformationError> {
        let enh = self.axes.normalize(coords);
        let llh = self.projection.unproj(enh)?;
        Ok(self.datum.llh_to_xyz(llh))
    }

    fn from_xyz(&self, xyz: XYZ, coords: &mut [Float]) -> Result<(), TransformationError> {
        let llh = self.datum.xyz_to_llh(xyz);
        let enh = self.projection.proj(llh)?;
        Ok(self.axes.denormalize(enh, coords))
    }
}

#[cfg(test)]
mod tests {
    use super::ProjectedCrs;
    use crate::{
        crs::Crs,
        cs::cartesian::projected::{ProjectedAxes, ProjectedTolerance},
        geodesy::geodetic_datum::consts::WGS84,
        units::length::M,
    };
    use std::any::{Any, TypeId};

    #[test]
    fn type_id() {
        let crs: &dyn Any = &ProjectedCrs {
            id: "UTM zone1".into(),
            datum: WGS84,
            axes: ProjectedAxes::EastNorthUp {
                horiz_unit: M,
                height_unit: M,
            },
            projection: crate::projections::ProjectionSpec::UTMNorth { zone: 1 },
        };

        assert_eq!(crs.type_id(), TypeId::of::<ProjectedCrs>());
        assert!(crs.is::<ProjectedCrs>());
    }

    #[test]
    fn approx_eq() {
        let crs = ProjectedCrs {
            id: "UTM Zon1".into(),
            datum: WGS84,
            axes: ProjectedAxes::EastNorth { horiz_unit: M },
            projection: crate::projections::ProjectionSpec::UTMNorth { zone: 1 },
        };

        assert!(crs.approx_eq(&[1., 1. + 1e-3], &[1., 1.], ProjectedTolerance::small()));
        assert!(!crs.approx_eq(&[1., 1. + 1e-1], &[1., 1.], ProjectedTolerance::small()));
    }

    #[test]
    fn dist() {
        todo!()
    }
}
