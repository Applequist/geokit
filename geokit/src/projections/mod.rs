use derive_more::derive::Display;
use thiserror::Error;

use crate::{
    cs::{
        cartesian::ENH,
        geodetic::{Lat, Lon, LLH},
    },
    math::Float,
    quantities::length::Length,
};
//use cyl::{Mercator, TransverseMercator, WebMercator};

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

#[derive(Error, Debug, Display)]
pub enum ProjectionError {
    InputOutOfBounds,
    OutputOutOfBounds,
}

pub trait Projection {
    fn proj(&self, llh: LLH) -> Result<ENH, ProjectionError>;
    fn unproj(&self, enh: ENH) -> Result<LLH, ProjectionError>;
}

/// TODO: Move that into transformation somewhere!!!
//impl ProjectionSpec {
//    pub fn projection(&self, ellipsoid: &Ellipsoid) -> Box<dyn Operation> {
//        match *self {
//            Self::Mercator1SP {
//                lon0,
//                k0,
//                false_easting,
//                false_northing,
//            } => Mercator::new_1_sp(
//                ellipsoid,
//                lon0.rad(),
//                k0,
//                false_easting.m(),
//                false_northing.m(),
//            )
//            .boxed(),
//            Self::Mercator2SP {
//                lon0,
//                lat0,
//                false_easting,
//                false_northing,
//            } => Mercator::new_2_sp(
//                ellipsoid,
//                lon0.rad(),
//                lat0.rad(),
//                false_easting.m(),
//                false_northing.m(),
//            )
//            .boxed(),
//            Self::UTMNorth { zone } => TransverseMercator::new_utm_north(ellipsoid, zone).boxed(),
//            Self::UTMSouth { zone } => TransverseMercator::new_utm_south(ellipsoid, zone).boxed(),
//            Self::TransverseMercator {
//                lon0,
//                lat0,
//                k0,
//                false_easting,
//                false_northing,
//            } => TransverseMercator::new(
//                ellipsoid,
//                lon0.rad(),
//                lat0.rad(),
//                k0,
//                false_easting.m(),
//                false_northing.m(),
//            )
//            .boxed(),
//            Self::WebMercator {
//                lon0,
//                lat0,
//                false_easting,
//                false_northing,
//            } => WebMercator::new(
//                ellipsoid,
//                lon0.rad(),
//                lat0.rad(),
//                false_easting.m(),
//                false_northing.m(),
//            )
//            .boxed(),
//        }
//    }
// }
pub mod cyl;
