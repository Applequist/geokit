use crate::crs::{geocentric::GeocentricCrs, geodetic::GeodeticCrs, CoordSpace, Crs};
use crate::geodesy::{GeodeticDatum, PrimeMeridian};
use crate::id::Id;
use crate::transformation::{Identity, Transformation};

/// A `TopocentricCrs` is a **3D cartesian coordinates reference system** whose origin is specified
/// as a geodetic location in a base 3D geodetic CRS and axes are derived from the base CRS.
#[derive(Debug, Clone, PartialEq)]
pub struct TopocentricCrs {
    id: Id,
    /// The base geodetic CRS.
    base_crs: GeodeticCrs,
    /// The origin of the 3D cartesian frame given in the base CRS. The axes are derive from the base CRS at the given
    /// location.
    origin: (f64, f64, f64),
    /// The length unit used for the coordinates in this CRS.
    length_unit: f64,
}

impl TopocentricCrs {
    /// Create a new 3D topocentric CRS.
    pub fn new(id: Id, base_crs: GeodeticCrs, origin: (f64, f64, f64), length_unit: f64) -> Self {
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

    fn normalized(
        &self,
    ) -> (
        Box<dyn Crs>,
        Box<dyn Transformation>,
        Box<dyn Transformation>,
    ) {
        // FIX: replace with proper values
        (
            Box::<GeocentricCrs>::default(),
            Identity::<3>.boxed(),
            Identity::<3>.boxed(),
        )
    }

    fn lowered(
        &self,
    ) -> Option<(
        Box<dyn Crs>,
        Box<dyn Transformation>,
        Box<dyn Transformation>,
    )> {
        Some((
            Box::new(GeocentricCrs::new(
                Id::name(format!("Lowered from {}", self.id())),
                GeodeticDatum::new(
                    Id::name(format!("Derived from {}", self.datum().id())),
                    self.datum().ellipsoid(),
                    PrimeMeridian::default(),
                    None,
                ),
            )),
            // FIX: Replace with proper transformation
            Identity::<3>.boxed(),
            Identity::<3>.boxed(),
        ))
    }
}
