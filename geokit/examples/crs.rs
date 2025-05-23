use geokit::crs::{geocentric::GeocentricCrs, geographic::GeographicCrs, projected::ProjectedCrs};
use geokit::cs::cartesian::{geocentric::GeocentricAxes, projected::ProjectedAxes};
use geokit::cs::geodetic::{GeodeticAxes, Lat, Lon};
use geokit::geodesy::{ellipsoid, prime_meridian, Ellipsoid, GeodeticDatum, PrimeMeridian};
use geokit::projections::ProjectionSpec;
use geokit::units::angle::DEG;
use geokit::units::length::M;

fn main() {
    let geoc = GeocentricCrs {
        id: "Geocentric CRS Example".into(),
        datum: GeodeticDatum::new(
            "WGS84",
            Ellipsoid::from_ainvf("WGS84", 6_378_137.0 * M, 298.257223563),
            PrimeMeridian::new("Greenwich", Lon::new(0.0 * DEG)),
        ),
        axes: GeocentricAxes::XYZ,
    };

    println!("{:#?}", geoc);

    let geog = GeographicCrs {
        id: "Geographic 3D Crs Example".into(),
        datum: GeodeticDatum::new(
            "GGRS87",
            ellipsoid::consts::GRS80,
            prime_meridian::consts::GREENWICH,
        ),
        axes: GeodeticAxes::EastNorthUp {
            angle_unit: DEG,
            height_unit: M,
        },
    };
    println!("{:#?}", geog);

    let proj = ProjectedCrs {
        id: "GRS80 + UTM 0".into(),
        datum: GeodeticDatum::new(
            "n/a",
            ellipsoid::consts::GRS80,
            prime_meridian::consts::GREENWICH,
        ),
        axes: ProjectedAxes::EastNorth { horiz_unit: M },
        projection: ProjectionSpec::TransverseMercator {
            lon0: Lon::ZERO,
            lat0: Lat::ZERO,
            k0: 0.9996,
            false_easting: 500_000.0 * M,
            false_northing: 0.0 * M,
        },
    };
    println!("{:#?}", proj);
}
