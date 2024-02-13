use std::any::Any;

use crate::{
    crs::{projected::ProjectedCrs, Crs, GeocentricCrs, GeographicCrs},
    geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian},
    operation::{identity, Operation},
};

pub trait GeodesyProvider {
    fn ellipsoid_ids(&self) -> impl Iterator<Item = &str>;
    fn ellipsoid(&self, id: &str) -> Option<Ellipsoid>;

    fn prime_meridian_ids(&self) -> impl Iterator<Item = &str>;
    fn prime_meridian(&self, id: &str) -> Option<PrimeMeridian>;

    fn datum_ids(&self) -> impl Iterator<Item = &str>;
    fn datum(&self, id: &str) -> Option<GeodeticDatum>;
}

pub trait CrsProvider {
    fn crs_ids(&self) -> impl Iterator<Item = &str>;
    fn crs(&self, id: &str) -> Option<Box<dyn Crs>>;
}

pub type TransformationPath = (Box<dyn Operation>, Box<dyn Operation>);
pub enum TransformationProviderError {
    NoTransformationPath,
}

pub trait TransformationProvider<S, D> {
    fn transformation(from: S, to: D) -> Result<TransformationPath, TransformationProviderError>;
}

pub struct DefaultTransformationProvider;

impl TransformationProvider<&dyn Any, &dyn Any> for DefaultTransformationProvider {
    fn transformation(
        from: &dyn Any,
        to: &dyn Any,
    ) -> Result<TransformationPath, TransformationProviderError> {
        // TODO: Add the missing datum transformation
        if let Some(proj_src) = from.downcast_ref::<ProjectedCrs>() {
            let (src_to_geoc, geoc_to_src) = proj_src.to_geoc();
            if let Some(proj_dst) = to.downcast_ref::<ProjectedCrs>() {
                let (dst_to_geoc, geoc_to_dst) = proj_dst.to_geoc();
                Ok((
                    Box::new(src_to_geoc.and_then(geoc_to_dst)),
                    Box::new(dst_to_geoc.and_then(geoc_to_src)),
                ))
            } else if let Some(geog_dst) = to.downcast_ref::<GeographicCrs>() {
                let (dst_to_geoc, geoc_to_dst) = geog_dst.to_geoc();
                Ok((
                    Box::new(src_to_geoc.and_then(geoc_to_dst)),
                    Box::new(dst_to_geoc.and_then(geoc_to_src)),
                ))
            } else if let Some(_geoc_crs) = to.downcast_ref::<GeocentricCrs>() {
                Ok((Box::new(src_to_geoc), Box::new(geoc_to_src)))
            } else {
                Err(TransformationProviderError::NoTransformationPath)
            }
        } else if let Some(geog_src) = from.downcast_ref::<GeographicCrs>() {
            let (src_to_geoc, geoc_to_src) = geog_src.to_geoc();
            if let Some(proj_dst) = to.downcast_ref::<ProjectedCrs>() {
                let (dst_to_geoc, geoc_to_dst) = proj_dst.to_geoc();
                Ok((
                    Box::new(src_to_geoc.and_then(geoc_to_dst)),
                    Box::new(dst_to_geoc.and_then(geoc_to_src)),
                ))
            } else if let Some(geog_dst) = to.downcast_ref::<GeographicCrs>() {
                let (dst_to_geoc, geoc_to_dst) = geog_dst.to_geoc();
                Ok((
                    Box::new(src_to_geoc.and_then(geoc_to_dst)),
                    Box::new(dst_to_geoc.and_then(geoc_to_src)),
                ))
            } else if let Some(_geoc_crs) = to.downcast_ref::<GeocentricCrs>() {
                Ok((Box::new(src_to_geoc), Box::new(geoc_to_src)))
            } else {
                Err(TransformationProviderError::NoTransformationPath)
            }
        } else if let Some(_geoc_src) = from.downcast_ref::<GeocentricCrs>() {
            if let Some(proj_dst) = to.downcast_ref::<ProjectedCrs>() {
                let (dst_to_geoc, geoc_to_dst) = proj_dst.to_geoc();
                Ok((Box::new(geoc_to_dst), Box::new(dst_to_geoc)))
            } else if let Some(geog_dst) = to.downcast_ref::<GeographicCrs>() {
                let (dst_to_geoc, geoc_to_dst) = geog_dst.to_geoc();
                Ok((Box::new(geoc_to_dst), Box::new(dst_to_geoc)))
            } else if let Some(_geoc_crs) = to.downcast_ref::<GeocentricCrs>() {
                Ok((Box::new(identity::<3>()), Box::new(identity::<3>())))
            } else {
                Err(TransformationProviderError::NoTransformationPath)
            }
        } else {
            Err(TransformationProviderError::NoTransformationPath)
        }
    }
}
