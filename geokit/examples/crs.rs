use geokit::crs::{Crs, GeocentricAxes, GeodeticAxes, ProjectedAxes, ProjectionSpec};
use geokit::cs::geodetic::{Lat, Lon};
use geokit::geodesy::geodetic_datum::DatumTransformation;
use geokit::geodesy::{ellipsoid, prime_meridian, Ellipsoid, GeodeticDatum, PrimeMeridian};
use geokit::quantity::angle::units::Deg;
use geokit::quantity::length::units::M;

fn main() {
    let geoc = Crs::Geocentric {
        id: "Geocentric CRS Example".into(),
        datum: GeodeticDatum::new(
            "WGS84",
            Ellipsoid::from_ainvf("WGS84", 6_378_137.0 * M, 298.257223563),
            PrimeMeridian::new("Greenwich", Lon::new(0.0 * Deg)),
            None,
        ),
        axes: GeocentricAxes::XYZ,
    };

    println!("CRS: {:#?}", geoc);

    let geog = Crs::Geographic {
        id: "Geographic 3D Crs Example".into(),
        datum: GeodeticDatum::new(
            "GGRS87",
            ellipsoid::consts::GRS80,
            prime_meridian::consts::GREENWICH,
            Some((
                "WGS84",
                DatumTransformation::GeocentricTranslation {
                    tx: (-199.87 * M).m(),
                    ty: (74.79 * M).m(),
                    tz: (246.64 * M).m(),
                },
            )),
        ),
        axes: GeodeticAxes::EastNorthUp {
            angle_unit: Deg::UNIT,
            height_unit: M,
        },
    };
    println!("CRS: {:#?}", geog);

    let proj = Crs::Projected {
        id: "GRS80 + UTM 0".into(),
        datum: GeodeticDatum::new(
            "n/a",
            ellipsoid::consts::GRS80,
            prime_meridian::consts::GREENWICH,
            None,
        ),
        axes: ProjectedAxes::EastNorth { horiz_unit: M },
        projection: ProjectionSpec::TransverseMercator {
            lon0: Lon::zero(),
            lat0: Lat::zero(),
            k0: 0.9996,
            false_easting: 500_000.0 * M,
            false_northing: 0.0 * M,
        },
    };
    println!("Destination CRS: {:#?}", proj);
}
