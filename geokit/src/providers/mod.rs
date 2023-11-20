use crate::{
    crs::Crs,
    geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian},
    id::Id,
    transformation::Transformation,
};

pub trait GeodesyProvider {
    fn ellipsoid_ids(&self) -> Vec<Id>;
    fn ellipsoid<I: Into<Id>>(&self, id: I) -> Option<Ellipsoid>;

    fn prime_meridian_ids(&self) -> Vec<Id>;
    fn prime_meridian<I: Into<Id>>(&self, id: I) -> Option<PrimeMeridian>;

    fn datum_ids(&self) -> Vec<Id>;
    fn datum<I: Into<Id>>(&self, id: I) -> Option<GeodeticDatum>;
}

pub trait CrsProvider {
    fn crs_ids(&self) -> Vec<Id>;
    fn crs(&self, id: Id) -> Option<Box<dyn Crs>>;
}

pub type TransformationPath = (Box<dyn Transformation>, Box<dyn Transformation>);
pub enum TransformationProviderError {
    NoTransformationPath,
}

pub trait TransformationProvider<S, D> {
    fn transformation(from: &S, to: &D) -> Result<TransformationPath, TransformationProviderError>;
}

pub mod default;
