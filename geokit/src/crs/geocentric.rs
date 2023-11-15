use crate::crs::{CoordSpace, Crs};
use crate::geodesy::GeodeticDatum;
use crate::id::Id;
use crate::transformation::InvertibleTransformation;

/// A `GeocentricCrs` is a **3D cartesian coordinates reference system** in which
/// coordinates are given by distance **in meters** along the following axes:
/// - X: axis from the center of the datum's ellipsoid in the equatorial and prime meridian plane,
/// - Y: axis from the center of the datum's ellipsoid in the equatorial plane and 90 degrees
/// meridian plane (east)
/// - Z: axis from the center of the datum's ellipsoid through the north pole.
#[derive(Debug, Clone, PartialEq)]
pub struct GeocentricCrs {
    id: Id,
    datum: GeodeticDatum,
}

impl GeocentricCrs {
    /// Creates a new [`GeocentricCrs`].
    pub fn new(id: Id, datum: GeodeticDatum) -> Self {
        Self { id, datum }
    }
}

impl Crs for GeocentricCrs {
    fn id(&self) -> &Id {
        &self.id
    }

    fn dim(&self) -> usize {
        3
    }

    fn kind(&self) -> CoordSpace {
        CoordSpace::Geocentric
    }

    fn datum(&self) -> &GeodeticDatum {
        &self.datum
    }

    fn lower(&self) -> Option<(Box<dyn Crs>, Box<dyn InvertibleTransformation>)> {
        None
    }
}

impl Default for GeocentricCrs {
    fn default() -> Self {
        GeocentricCrs::new(Id::name("WGS 84 (geocentric)"), GeodeticDatum::default())
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::geodesy::*;
    use crate::id::Id;
    #[test]
    fn clone() {
        let geoc = GeocentricCrs::default();
        let cpy = geoc.clone();
        assert_eq!(geoc, cpy);
    }

    #[test]
    fn partial_eq() {
        let geoc = GeocentricCrs::default();
        let cpy = geoc.clone();
        assert!(geoc.eq(&cpy));
        assert!(!geoc.ne(&cpy));

        let mut different_id = geoc.clone();
        different_id.id = Id::name("WGS 84.1");
        assert_ne!(geoc, different_id);

        let mut different_datum = geoc.clone();
        different_datum.datum = GeodeticDatum::new(
            Id::name("WGS 84.1"),
            Ellipsoid::default(),
            PrimeMeridian::default(),
            None,
        );
        assert_ne!(geoc, different_datum);
    }
}
