use geokit::cs::geodetic::{Lat, Lon};
use geokit::quantity::angle::units::DEG;
use geokit::quantity::length::units::M;
use geokit::{
    crs::{Crs::Geographic, Crs::Projected, GeodeticAxes, ProjectedAxes, ProjectionSpec},
    geodesy::{ellipsoid, prime_meridian, GeodeticDatum},
    operation::Operation,
    providers::{DefaultTransformationProvider, TransformationProvider},
};

fn main() {
    let src = Geographic {
        id: "GRS80 (Geographic 3D)".into(),
        datum: GeodeticDatum::new(
            "n/a",
            ellipsoid::consts::GRS80,
            prime_meridian::consts::GREENWICH,
            None,
        ),
        axes: GeodeticAxes::EastNorthUp {
            angle_unit: DEG,
            height_unit: M,
        },
    };
    println!("Source CRS: {:#?}", src);

    let dst = Projected {
        id: "GRS80 + UTM 0".into(),
        datum: GeodeticDatum::new(
            "n/a",
            ellipsoid::consts::GRS80,
            prime_meridian::consts::GREENWICH,
            None,
        ),
        axes: ProjectedAxes::EastNorth { horiz_unit: M },
        projection: ProjectionSpec::TransverseMercator {
            lon0: Lon::new(0.0 * DEG),
            lat0: Lat::new(0.0),
            k0: 0.9996,
            false_easting: 500_000.0 * M,
            false_northing: 0.0 * M,
        },
    };
    println!("Destination CRS: {:#?}", dst);

    let provider = DefaultTransformationProvider;
    let (src_to_dst, dst_to_src) = provider.transformation(&src, &dst).unwrap();

    let src_pt = vec![-10.0, -90.0, 0.0];
    let dst_pt = src_to_dst.fwd_new(&src_pt).unwrap();
    println!("{src_pt:?} --- src_to_dst ---> {dst_pt:?}");

    let dst_src_pt = dst_to_src.fwd_new(&dst_pt).unwrap();
    println!("{dst_pt:?} --- dst_to_src ---> {dst_src_pt:?}");
    // let src_dst_bwd_pt = src_to_dst.bwd_new(&dst_pt).unwrap();
    // println!("{dst_pt:?} --- src_to_dst.bwd ---> {src_dst_bwd_pt:?}");
}
