use std::fmt::Debug;

use smol_str::SmolStr;

use crate::{
    geodesy::{Ellipsoid, GeodeticDatum},
    operation::{
        conversion::{projection::WebMercator, GeogToGeoc, Normalization},
        Bwd, DynOperation, Fwd, Operation,
    },
};

use super::Crs;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProjectedAxes {
    EastNorthUp { horiz_unit: f64, height_unit: f64 },
    EastNorth { horiz_unit: f64 },
}

impl ProjectedAxes {
    pub fn dim(&self) -> usize {
        match self {
            ProjectedAxes::EastNorthUp {
                horiz_unit: _,
                height_unit: _,
            } => 3,
            ProjectedAxes::EastNorth { horiz_unit: _ } => 2,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProjectionSpec {
    WebMercator {
        lon0: f64,
        lat0: f64,
        false_easting: f64,
        false_northing: f64,
    },
}

impl ProjectionSpec {
    pub fn projection(&self, ellipsoid: &Ellipsoid) -> Box<dyn DynOperation> {
        match *self {
            Self::WebMercator {
                lon0,
                lat0,
                false_easting,
                false_northing,
            } => Box::new(WebMercator::new(
                ellipsoid.clone(),
                lon0,
                lat0,
                false_easting,
                false_northing,
            )),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ProjectedCrs {
    id: SmolStr,
    datum: GeodeticDatum,
    axes: ProjectedAxes,
    projection: ProjectionSpec,
}

impl ProjectedCrs {
    pub fn new<S: AsRef<str>>(
        id: S,
        datum: GeodeticDatum,
        axes: ProjectedAxes,
        projection: ProjectionSpec,
    ) -> Self {
        Self {
            id: SmolStr::new(id),
            datum,
            axes,
            projection,
        }
    }

    pub fn id(&self) -> &str {
        self.id.as_str()
    }

    #[inline]
    pub fn dim(&self) -> usize {
        self.axes.dim()
    }

    #[inline]
    pub fn datum(&self) -> &GeodeticDatum {
        &self.datum
    }

    #[inline]
    pub fn axes(&self) -> ProjectedAxes {
        self.axes
    }

    pub fn projection(&self) -> Box<dyn DynOperation> {
        self.projection.projection(self.datum().ellipsoid())
    }

    pub fn to_geoc(&self) -> (impl Operation, impl Operation) {
        let fwd = Fwd(Normalization::from(self.axes()))
            .and_then(Bwd(self.projection()))
            .and_then(Fwd(GeogToGeoc::new(self.datum())));
        let bwd = Bwd(GeogToGeoc::new(self.datum()))
            .and_then(Fwd(self.projection()))
            .and_then(Bwd(Normalization::from(self.axes())));
        (fwd, bwd)
    }
}

impl Crs for ProjectedCrs {
    fn is_normalized(&self) -> bool {
        if let ProjectedAxes::EastNorthUp {
            horiz_unit,
            height_unit,
        } = self.axes()
        {
            horiz_unit == 1.0 && height_unit == 1.0
        } else {
            false
        }
    }
}
