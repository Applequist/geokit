use super::Crs;
use crate::transformations::{
    ToXYZTransformation, ToXYZTransformationProvider, TransformationError,
};
use crate::{
    cs::cartesian::{ProjectedAxes, XYZ},
    geodesy::GeodeticDatum,
    math::fp::Float,
    projections::{Projection, ProjectionError, ProjectionSpec},
};
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
    fn id(&self) -> &str {
        &self.id
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
