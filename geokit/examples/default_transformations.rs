use geokit::crs::{geographic::GeographicCrs, projected::ProjectedCrs};
use geokit::cs::cartesian::ProjectedAxes;
use geokit::cs::geodetic::{GeodeticAxes, Lat, Lon};
use geokit::geodesy::{ellipsoid, prime_meridian, GeodeticDatum};
use geokit::projections::ProjectionSpec;
use geokit::transformations::xyz::XYZIdentity;
use geokit::transformations::{CrsTransformation, CrsXYZTransformation};
use geokit::units::angle::DEG;
use geokit::units::length::M;

fn main() {
    let src = GeographicCrs {
        id: "GRS80 (Geographic 3D)".into(),
        datum: GeodeticDatum::new(
            "n/a",
            ellipsoid::consts::GRS80,
            prime_meridian::consts::GREENWICH,
        ),
        axes: GeodeticAxes::EastNorthUp {
            angle_unit: DEG,
            height_unit: M,
        },
    };
    println!("Source CRS: {:#?}", src);

    let dst = ProjectedCrs {
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
    println!("Destination CRS: {:#?}", dst);

    let tx = CrsXYZTransformation::new(src, dst, XYZIdentity, XYZIdentity);

    let mut src_pt = [-10.0, -90.0, 0.0];
    let mut dst_pt = [0.0; 2];
    tx.src_to_dst(&src_pt, &mut dst_pt).unwrap();
    println!("{src_pt:?} --- src_to_dst ---> {dst_pt:?}");

    tx.dst_to_src(&dst_pt, &mut src_pt).unwrap();
    println!("{dst_pt:?} --- dst_to_src ---> {src_pt:?}");
}
