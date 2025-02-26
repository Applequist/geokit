use crate::crs::{GeocentricCrs, GeographicCrs, ProjectedCrs};
use crate::cs::cartesian::{ProjectedAxes, XYZ};
use crate::geodesy::GeodeticDatum;
use crate::math::fp::Float;
use crate::projections::{Projection, ProjectionError};
use derive_more::derive::Display;
use thiserror::Error;

pub mod xyz;

#[derive(Debug, Error, Display)]
pub enum TransformationError {
    OutOfBounds,
}

/// Trait for Coordinates transformation between two coordinates reference systems.
pub trait CrsTransformation {
    fn src_to_dst(&self, src: &[Float], dst: &mut [Float]) -> Result<(), TransformationError>;
    fn dst_to_src(&self, dst: &[Float], src: &mut [Float]) -> Result<(), TransformationError>;
}

pub trait ToXYZTransformation {
    fn to_xyz(&self, coords: &[Float]) -> Result<XYZ, TransformationError>;
    fn from_xyz(&self, xyz: XYZ, coords: &mut [Float]) -> Result<(), TransformationError>;
}

/// A trait for ***normalized geocentric coordinates*** transformations.
pub trait XYZTransformation {
    fn src_to_dst(&self, src: XYZ) -> XYZ;
    fn dst_to_src(&self, dst: XYZ) -> XYZ;
}

/// Coordinates transformation between two CRS via [normalized geocentric coordinates][XYZ].
pub struct CrsXYZTransformation<'a> {
    src_to_xyz: Box<dyn ToXYZTransformation + 'a>,
    src_xyz_to_ref_xyz: Box<dyn XYZTransformation + 'a>,
    dst_xyz_to_ref_xyz: Box<dyn XYZTransformation + 'a>,
    dst_to_xyz: Box<dyn ToXYZTransformation + 'a>,
}

impl<'a> CrsXYZTransformation<'a> {
    pub fn new<S: ProvideToXYZTransformation, D: ProvideToXYZTransformation>(
        src: S,
        dst: D,
        src_to_ref: impl XYZTransformation + 'a,
        dst_to_ref: impl XYZTransformation + 'a,
    ) -> Self {
        CrsXYZTransformation {
            src_to_xyz: src.to_xyz_transformation(),
            src_xyz_to_ref_xyz: Box::new(src_to_ref),
            dst_xyz_to_ref_xyz: Box::new(dst_to_ref),
            dst_to_xyz: dst.to_xyz_transformation(),
        }
    }
}

impl<'a> CrsTransformation for CrsXYZTransformation<'a> {
    fn src_to_dst(&self, src: &[Float], dst: &mut [Float]) -> Result<(), TransformationError> {
        let src_xyz = self.src_to_xyz.to_xyz(src)?;
        let ref_xyz = self.src_xyz_to_ref_xyz.src_to_dst(src_xyz);
        let dst_xyz = self.dst_xyz_to_ref_xyz.dst_to_src(ref_xyz);
        self.dst_to_xyz.from_xyz(dst_xyz, dst)
    }

    fn dst_to_src(&self, dst: &[Float], src: &mut [Float]) -> Result<(), TransformationError> {
        let dst_xyz = self.dst_to_xyz.to_xyz(dst)?;
        let ref_xyz = self.dst_xyz_to_ref_xyz.src_to_dst(dst_xyz);
        let src_xyz = self.src_xyz_to_ref_xyz.dst_to_src(ref_xyz);
        self.src_to_xyz.from_xyz(src_xyz, src)
    }
}

pub trait ProvideToXYZTransformation {
    // NOTE: as of v1.84, this doesn't work for our use case:
    // fn to_xyz_transformation(&self) -> impl ToXYZTransformation;
    // `use<...>` precise capturing syntax is currently not allowed in return-position
    // `impl Trait` in traits
    // See issue #130044 <https://github.com/rust-lang/rust/issues/130044>
    fn to_xyz_transformation<'a>(&self) -> Box<dyn ToXYZTransformation + 'a>;
}

impl ProvideToXYZTransformation for GeocentricCrs {
    fn to_xyz_transformation<'a>(&self) -> Box<dyn ToXYZTransformation + 'a> {
        Box::new(self.clone())
    }
}

impl ToXYZTransformation for GeocentricCrs {
    fn to_xyz(&self, coords: &[Float]) -> Result<XYZ, TransformationError> {
        Ok(self.axes.normalize(coords))
    }

    fn from_xyz(&self, xyz: XYZ, coords: &mut [Float]) -> Result<(), TransformationError> {
        Ok(self.axes.denormalize(xyz, coords))
    }
}

impl ProvideToXYZTransformation for GeographicCrs {
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

impl ProvideToXYZTransformation for ProjectedCrs {
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
    fn from(alue: ProjectionError) -> Self {
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
