use crate::{
    cs::{
        cartesian::ENH,
        geodetic::{Lat, Lon, LLH},
    },
    geodesy::Ellipsoid,
    math::fp::Float,
    quantities::length::Length,
};
use cyl::{mercator::Mercator, transverse_mercator::TransverseMercator, web_mercator::WebMercator};
use derive_more::derive::Display;
use thiserror::Error;

#[derive(Error, Debug, Display)]
pub enum ProjectionError {
    LLHOutOfBounds,
    ENHOutOfBounds,
}

pub trait Projection {
    fn proj(&self, llh: LLH) -> Result<ENH, ProjectionError>;
    fn unproj(&self, enh: ENH) -> Result<LLH, ProjectionError>;
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProjectionSpec {
    Mercator1SP {
        lon0: Lon,
        k0: Float,
        false_easting: Length,
        false_northing: Length,
    },
    Mercator2SP {
        lon0: Lon,
        lat0: Lat,
        false_easting: Length,
        false_northing: Length,
    },
    UTMNorth {
        zone: u8,
    },
    UTMSouth {
        zone: u8,
    },
    TransverseMercator {
        lon0: Lon,
        lat0: Lat,
        k0: Float,
        false_easting: Length,
        false_northing: Length,
    },
    WebMercator {
        lon0: Lon,
        lat0: Lat,
        false_easting: Length,
        false_northing: Length,
    },
}

mod cyl;

impl ProjectionSpec {
    pub fn projection(&self, ellipsoid: &Ellipsoid) -> Box<dyn Projection> {
        match *self {
            Self::Mercator1SP {
                lon0,
                k0,
                false_easting,
                false_northing,
            } => Box::new(Mercator::new_1_sp(
                ellipsoid,
                lon0,
                k0,
                false_easting,
                false_northing,
            )),
            Self::Mercator2SP {
                lon0,
                lat0,
                false_easting,
                false_northing,
            } => Box::new(Mercator::new_2_sp(
                ellipsoid,
                lon0,
                lat0,
                false_easting,
                false_northing,
            )),
            Self::UTMNorth { zone } => Box::new(TransverseMercator::new_utm_north(ellipsoid, zone)),
            Self::UTMSouth { zone } => Box::new(TransverseMercator::new_utm_south(ellipsoid, zone)),
            Self::TransverseMercator {
                lon0,
                lat0,
                k0,
                false_easting,
                false_northing,
            } => Box::new(TransverseMercator::new(
                ellipsoid,
                lon0,
                lat0,
                k0,
                false_easting,
                false_northing,
            )),
            Self::WebMercator {
                lon0,
                lat0,
                false_easting,
                false_northing,
            } => Box::new(WebMercator::new(
                ellipsoid,
                lon0,
                lat0,
                false_easting,
                false_northing,
            )),
        }
    }
}
