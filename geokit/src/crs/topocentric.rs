use crate::coord::{Coord3D, CoordSpace};
use crate::crs::{geocentric::GeocentricCrs, geodetic::Geodetic3DCrs, Crs};
use crate::geodesy::GeodeticDatum;
use crate::id::Id;
use crate::transformation::{Identity, InversibleTransformation};

/// A `TopocentricCrs` is a **3D cartesian coordinates reference system** whose origin is specified
/// as a geodetic location in a base 3D geodetic CRS and axes are derived from the base CRS.
pub struct TopocentricCrs {
    id: Id,
    /// The base 3D geodetic CRS.
    base_crs: Geodetic3DCrs,
    /// The origin of the 3D cartesian frame given in the base CRS. The axes are derive from the base CRS at the given
    /// location.
    origin: Coord3D,
    /// The length unit used for the coordinates in this CRS.
    length_unit: f64,
}

impl TopocentricCrs {
    /// Create a new 3D topocentric CRS.
    pub fn new(id: Id, base_crs: Geodetic3DCrs, origin: Coord3D, length_unit: f64) -> Self {
        Self {
            id,
            base_crs,
            origin,
            length_unit,
        }
    }
}

impl Crs for TopocentricCrs {
    fn id(&self) -> &Id {
        &self.id
    }

    fn dim(&self) -> usize {
        3
    }

    fn kind(&self) -> CoordSpace {
        CoordSpace::Topocentric
    }

    fn datum(&self) -> &GeodeticDatum {
        self.base_crs.datum()
    }

    fn lower(&self) -> Option<(Box<dyn Crs>, Box<dyn InversibleTransformation>)> {
        Some((
            Box::new(GeocentricCrs::new(
                Id::name("n/a"),
                self.base_crs.datum().clone(),
            )),
            // FIX: Replace with proper transformation
            Box::new(Identity),
        ))
    }
}
