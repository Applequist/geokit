use crate::{
    crs::Crs,
    geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian},
    operation::Operation,
    tag::Tag,
};

pub trait GeodesyProvider {
    fn ellipsoid_ids(&self) -> Vec<Tag>;
    fn ellipsoid<I: Into<Tag>>(&self, id: I) -> Option<Ellipsoid>;

    fn prime_meridian_ids(&self) -> Vec<Tag>;
    fn prime_meridian<I: Into<Tag>>(&self, id: I) -> Option<PrimeMeridian>;

    fn datum_ids(&self) -> Vec<Tag>;
    fn datum<I: Into<Tag>>(&self, id: I) -> Option<GeodeticDatum>;
}

pub trait CrsProvider {
    fn crs_ids(&self) -> Vec<Tag>;
    fn crs(&self, id: Tag) -> Option<Box<dyn Crs>>;
}

pub type TransformationPath = (Box<dyn Operation>, Box<dyn Operation>);
pub enum TransformationProviderError {
    NoTransformationPath,
}

pub trait TransformationProvider<S, D> {
    fn transformation(from: &S, to: &D) -> Result<TransformationPath, TransformationProviderError>;
}

pub mod default;
