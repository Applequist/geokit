use crate::crs::{geodetic::GeodeticCrs, CoordSpace, Crs};
use crate::geodesy::GeodeticDatum;
use crate::id::Id;
use crate::transformation::{CoordScaling, Identity, Transformation};

use super::geocentric::GeocentricCrs;

/// A `TopocentricCrs` is a **3D cartesian coordinates reference system** whose origin is specified
/// as a geodetic location in a base 3D geodetic CRS and axes are:
/// - X pointing East from the origin,
/// - Y pointing North from the origin
/// - Z pointing up at the origin
/// All use the given length unit.
#[derive(Debug, Clone, PartialEq)]
pub struct TopocentricCrs {
    id: Id,
    /// The base geodetic CRS.
    base_crs: GeodeticCrs,
    /// The origin of the 3D cartesian frame given in the base CRS. The axes are derive from the base CRS at the given
    /// location.
    origin: [f64; 3],
    /// The length unit used for the coordinates in this CRS.
    length_unit: f64,
}

impl TopocentricCrs {
    /// Create a new 3D topocentric CRS.
    pub fn new(id: Id, base_crs: GeodeticCrs, origin: [f64; 3], length_unit: f64) -> Self {
        Self {
            id,
            base_crs,
            origin,
            length_unit,
        }
    }

    pub fn geoc_crs(&self) -> (GeocentricCrs, Identity<3>, Identity<3>) {
        // TODO: Transformation part.
        let (crs, base_to_geoc, _geoc_to_base) = self.base_crs.geoc_crs();
        // Transform origin to geocentric coordinates
        let mut origin_geoc = [0.; 3];
        base_to_geoc.apply(&self.origin, &mut origin_geoc).unwrap();

        // TODO: Compute the mat to go from the local topocentric frame to the geocentric frame.
        // then apply normalization then affine transform
        (crs, Identity::<3>, Identity::<3>)
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

    fn normalized(
        &self,
    ) -> (
        Box<dyn Crs>,
        Box<dyn Transformation>,
        Box<dyn Transformation>,
    ) {
        (
            Box::new(TopocentricCrs::new(
                Id::name(format!("{} normalized", self.id())),
                self.base_crs.clone(),
                self.origin,
                1.0,
            )),
            CoordScaling::from((3, self.length_unit)).boxed(),
            CoordScaling::from((3, 1.0 / self.length_unit)).boxed(),
        )
    }

    fn lowered(
        &self,
    ) -> Option<(
        Box<dyn Crs>,
        Box<dyn Transformation>,
        Box<dyn Transformation>,
    )> {
        let (crs, to, from) = self.geoc_crs();
        // TODO: build TopocentricToGeocentric and GeocentricToTopocentric
        Some((Box::new(crs), to.boxed(), from.boxed()))
    }
}
