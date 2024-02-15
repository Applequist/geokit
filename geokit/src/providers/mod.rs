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
        src: &dyn Any,
        dst: &dyn Any,
    ) -> Result<TransformationPath, TransformationProviderError> {
        if let Some(proj_src) = src.downcast_ref::<ProjectedCrs>() {
            let src_ref_datum_id = proj_src.datum().ref_datum_id();
            let (src_to_geoc, geoc_to_src) = proj_src.to_geoc();
            if let Some(proj_dst) = dst.downcast_ref::<ProjectedCrs>() {
                let dst_ref_datum_id = proj_dst.datum().ref_datum_id();
                if src_ref_datum_id != dst_ref_datum_id {
                    return Err(TransformationProviderError::NoTransformationPath);
                }
                let (dst_to_geoc, geoc_to_dst) = proj_dst.to_geoc();
                return Ok((
                    Box::new(src_to_geoc.and_then(geoc_to_dst)),
                    Box::new(dst_to_geoc.and_then(geoc_to_src)),
                ));
            }
            if let Some(geog_dst) = dst.downcast_ref::<GeographicCrs>() {
                let dst_ref_datum_id = geog_dst.datum().ref_datum_id();
                if src_ref_datum_id != dst_ref_datum_id {
                    return Err(TransformationProviderError::NoTransformationPath);
                }
                let (dst_to_geoc, geoc_to_dst) = geog_dst.to_geoc();
                return Ok((
                    Box::new(src_to_geoc.and_then(geoc_to_dst)),
                    Box::new(dst_to_geoc.and_then(geoc_to_src)),
                ));
            }
            if let Some(geoc_dst) = dst.downcast_ref::<GeocentricCrs>() {
                let dst_ref_datum_id = geoc_dst.datum().ref_datum_id();
                if src_ref_datum_id != dst_ref_datum_id {
                    return Err(TransformationProviderError::NoTransformationPath);
                }
                return Ok((Box::new(src_to_geoc), Box::new(geoc_to_src)));
            }
            Err(TransformationProviderError::NoTransformationPath)
        } else if let Some(geog_src) = src.downcast_ref::<GeographicCrs>() {
            let src_ref_datum_id = geog_src.datum().ref_datum_id();
            let (src_to_geoc, geoc_to_src) = geog_src.to_geoc();
            if let Some(proj_dst) = dst.downcast_ref::<ProjectedCrs>() {
                let dst_ref_datum_id = proj_dst.datum().ref_datum_id();
                if src_ref_datum_id != dst_ref_datum_id {
                    return Err(TransformationProviderError::NoTransformationPath);
                }
                let (dst_to_geoc, geoc_to_dst) = proj_dst.to_geoc();
                return Ok((
                    Box::new(src_to_geoc.and_then(geoc_to_dst)),
                    Box::new(dst_to_geoc.and_then(geoc_to_src)),
                ));
            }
            if let Some(geog_dst) = dst.downcast_ref::<GeographicCrs>() {
                let dst_ref_datum_id = geog_dst.datum().ref_datum_id();
                if src_ref_datum_id != dst_ref_datum_id {
                    return Err(TransformationProviderError::NoTransformationPath);
                }
                let (dst_to_geoc, geoc_to_dst) = geog_dst.to_geoc();
                return Ok((
                    Box::new(src_to_geoc.and_then(geoc_to_dst)),
                    Box::new(dst_to_geoc.and_then(geoc_to_src)),
                ));
            }
            if let Some(geoc_dst) = dst.downcast_ref::<GeocentricCrs>() {
                let dst_ref_datum_id = geoc_dst.datum().ref_datum_id();
                if src_ref_datum_id != dst_ref_datum_id {
                    return Err(TransformationProviderError::NoTransformationPath);
                }
                return Ok((Box::new(src_to_geoc), Box::new(geoc_to_src)));
            }
            Err(TransformationProviderError::NoTransformationPath)
        } else if let Some(geoc_src) = src.downcast_ref::<GeocentricCrs>() {
            let src_ref_datum_id = geoc_src.datum().ref_datum_id();
            if let Some(proj_dst) = dst.downcast_ref::<ProjectedCrs>() {
                let dst_ref_datum_id = proj_dst.datum().ref_datum_id();
                if src_ref_datum_id != dst_ref_datum_id {
                    return Err(TransformationProviderError::NoTransformationPath);
                }
                let (dst_to_geoc, geoc_to_dst) = proj_dst.to_geoc();
                return Ok((Box::new(geoc_to_dst), Box::new(dst_to_geoc)));
            }
            if let Some(geog_dst) = dst.downcast_ref::<GeographicCrs>() {
                let dst_ref_datum_id = geog_dst.datum().ref_datum_id();
                if src_ref_datum_id != dst_ref_datum_id {
                    return Err(TransformationProviderError::NoTransformationPath);
                }
                let (dst_to_geoc, geoc_to_dst) = geog_dst.to_geoc();
                return Ok((Box::new(geoc_to_dst), Box::new(dst_to_geoc)));
            }
            if let Some(geoc_dst) = dst.downcast_ref::<GeocentricCrs>() {
                let dst_ref_datum_id = geoc_dst.datum().ref_datum_id();
                if src_ref_datum_id != dst_ref_datum_id {
                    return Err(TransformationProviderError::NoTransformationPath);
                }
                return Ok((Box::new(identity::<3>()), Box::new(identity::<3>())));
            }
            Err(TransformationProviderError::NoTransformationPath)
        } else {
            Err(TransformationProviderError::NoTransformationPath)
        }
    }
}
