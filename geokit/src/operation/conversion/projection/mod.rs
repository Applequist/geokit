use cyl::{Mercator, TransverseMercator, WebMercator};

use crate::{
    cs::geodetic::{Lat, Lon},
    geodesy::Ellipsoid,
    operation::Operation,
    quantities::length::Length,
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProjectionSpec {
    Mercator1SP {
        lon0: Lon,
        k0: f64,
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
        k0: f64,
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

/// TODO: Move that into transformation somewhere!!!
impl ProjectionSpec {
    pub fn projection(&self, ellipsoid: &Ellipsoid) -> Box<dyn Operation> {
        match *self {
            Self::Mercator1SP {
                lon0,
                k0,
                false_easting,
                false_northing,
            } => Mercator::new_1_sp(
                ellipsoid,
                lon0.rad(),
                k0,
                false_easting.m(),
                false_northing.m(),
            )
            .boxed(),
            Self::Mercator2SP {
                lon0,
                lat0,
                false_easting,
                false_northing,
            } => Mercator::new_2_sp(
                ellipsoid,
                lon0.rad(),
                lat0.rad(),
                false_easting.m(),
                false_northing.m(),
            )
            .boxed(),
            Self::UTMNorth { zone } => TransverseMercator::new_utm_north(ellipsoid, zone).boxed(),
            Self::UTMSouth { zone } => TransverseMercator::new_utm_south(ellipsoid, zone).boxed(),
            Self::TransverseMercator {
                lon0,
                lat0,
                k0,
                false_easting,
                false_northing,
            } => TransverseMercator::new(
                ellipsoid,
                lon0.rad(),
                lat0.rad(),
                k0,
                false_easting.m(),
                false_northing.m(),
            )
            .boxed(),
            Self::WebMercator {
                lon0,
                lat0,
                false_easting,
                false_northing,
            } => WebMercator::new(
                ellipsoid,
                lon0.rad(),
                lat0.rad(),
                false_easting.m(),
                false_northing.m(),
            )
            .boxed(),
        }
    }
}
pub mod cyl;
