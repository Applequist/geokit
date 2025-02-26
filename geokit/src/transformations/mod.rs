//! This module defines [Crs](Coordinates Reference Systems).
//! A [Crs] ties a coordinates system to the Earth using a [GeodeticDatum](datum) and allows to
//! unambiguously assign coordinates to location on Earth.

use crate::cs::cartesian::XYZ;
use crate::math::fp::Float;
use derive_more::derive::Display;
use thiserror::Error;

#[derive(Debug, Error, Display)]
pub enum TransformationError {
    OutOfBounds,
}

/// Trait for coordinates transformation between two [coordinates reference systems](Crs)
pub trait CrsTransformation {
    fn src_to_dst(&self, src: &[Float], dst: &mut [Float]) -> Result<(), TransformationError>;
    fn dst_to_src(&self, dst: &[Float], src: &mut [Float]) -> Result<(), TransformationError>;
}

pub trait ToXYZTransformationProvider {
    // NOTE: as of v1.84, this doesn't work for our use case:
    // fn to_xyz_transformation(&self) -> impl ToXYZTransformation;
    // `use<...>` precise capturing syntax is currently not allowed in return-position
    // `impl Trait` in traits
    // See issue #130044 <https://github.com/rust-lang/rust/issues/130044>
    fn to_xyz_transformation<'a>(&self) -> Box<dyn ToXYZTransformation + 'a>;
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

/// Coordinates transformation between two CRS via [normalized geocentric coordinates](XYZ).
pub struct CrsXYZTransformation<'a> {
    src_to_xyz: Box<dyn ToXYZTransformation + 'a>,
    src_xyz_to_ref_xyz: Box<dyn XYZTransformation + 'a>,
    dst_xyz_to_ref_xyz: Box<dyn XYZTransformation + 'a>,
    dst_to_xyz: Box<dyn ToXYZTransformation + 'a>,
}

impl<'a> CrsXYZTransformation<'a> {
    pub fn new<S: ToXYZTransformationProvider, D: ToXYZTransformationProvider>(
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

pub mod xyz;
