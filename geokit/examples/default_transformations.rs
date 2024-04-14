use geokit::{
    crs::{Crs::Geographic, Crs::Projected, GeodeticAxes, ProjectedAxes, ProjectionSpec},
    geodesy::{ellipsoid, prime_meridian, GeodeticDatum},
    operation::Operation,
    providers::{DefaultTransformationProvider, TransformationProvider},
};
use geokit::units::angle::{Angle, DEG};
use geokit::units::length::{Length, M};

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
            angle_unit: DEG.to_radians(), // degrees
            height_unit: M.to_meters(),                 // metres
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
        axes: ProjectedAxes::EastNorth { horiz_unit: 1.0 },
        projection: ProjectionSpec::TransverseMercator {
            lon0: 0.0,
            lat0: 0.0,
            k0: 0.9996,
            false_easting: 500_000.0,
            false_northing: 0.0,
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
    let src_dst_bwd_pt = src_to_dst.bwd_new(&dst_pt).unwrap();
    println!("{dst_pt:?} --- src_to_dst.bwd ---> {src_dst_bwd_pt:?}");
}
