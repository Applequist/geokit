use approx::assert_abs_diff_eq;
use geokit::units::angle::{Angle, DEG};
use geokit::units::length::{Length, M};
use geokit::{
    crs::{Crs::Geographic, Crs::Projected, GeodeticAxes, ProjectedAxes, ProjectionSpec},
    geodesy::{ellipsoid, prime_meridian, GeodeticDatum},
    operation::{self, Operation},
    providers::{DefaultTransformationProvider, TransformationProvider},
};
use std::default::Default;

fn dist(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), 3);
    a.iter()
        .zip(b.iter())
        .fold(0.0, |mut acc, e| {
            let diff = e.0 - e.1;
            acc += diff * diff;
            acc
        })
        .sqrt()
}

mod input;

#[test]
fn llh_to_utm0() -> operation::Result<()> {
    let llh_orig = input::read_coords("tests/data/llh_grid_restricted.txt", 3).collect::<Vec<_>>();
    let count_llh = llh_orig.len();
    println!("Read {count_llh} 2D coordinates");

    let enh_orig =
        input::read_coords("tests/data/grs80_utm0_restricted.txt", 3).collect::<Vec<_>>();
    let count_enh = enh_orig.len();
    println!("Read {count_enh} 3D coordinates");

    assert_eq!(
        count_llh, count_enh,
        "Expected #llh == #enh. Got {} != {}",
        count_llh, count_enh
    );

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
            height_unit: M.to_meters(),   // metres
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
        axes: ProjectedAxes::EastNorth {
            horiz_unit: M.to_meters(),
        },
        projection: ProjectionSpec::TransverseMercator {
            lon0: 0.0,
            lat0: 0.0,
            k0: 0.9996,
            false_easting: 500_000.0,
            false_northing: 0.0,
        },
    };

    let provider = DefaultTransformationProvider;
    let (src_to_dst, dst_to_src) = provider.transformation(&src, &dst).unwrap();

    let src_pt = vec![-10.0, -90.0, 0.0];
    let dst_pt = src_to_dst.fwd_new(&src_pt).unwrap();
    println!("{src_pt:?} --- src_to_dst ---> {dst_pt:?}");

    let dst_src_pt = dst_to_src.fwd_new(&dst_pt).unwrap();
    println!("{dst_pt:?} --- dst_to_src ---> {dst_src_pt:?}");
    let src_dst_bwd_pt = src_to_dst.bwd_new(&dst_pt).unwrap();
    println!("{dst_pt:?} --- src_to_dst.bwd ---> {src_dst_bwd_pt:?}");
    assert_abs_diff_eq!(&dst_src_pt[..], &src_dst_bwd_pt[..]);

    // Allocating storage for transformed coordinates.
    println!("Converting llh coordinates to enh coordinates...");
    let mut enh = vec![vec![0.; 3]; count_llh];
    let mut trans_count = 0;
    for (i, o) in llh_orig.iter().zip(enh.iter_mut()) {
        src_to_dst.apply_fwd(i, o)?;
        trans_count += 1;
    }
    assert_eq!(
        trans_count, count_llh,
        "Expected #ops = {count_llh}. Got {trans_count}"
    );

    println!("Checking for equality...");
    for (p_enh, p_enh_orig) in enh.iter().zip(enh_orig.iter()) {
        let err = dist(p_enh, p_enh_orig);
        assert!(
            err < 1e-3,
            "Expected dist(actual, expected) < 1e-3. Got {err}"
        );
    }

    println!("Converting enh coordinates to llh coordinates...");
    let mut llh = vec![vec![0.; 3]; count_enh];
    let mut trans_count = 0;
    for (i, o) in enh.iter().zip(llh.iter_mut()) {
        dst_to_src.apply_fwd(i, o)?;
        trans_count += 1;
    }
    assert_eq!(
        trans_count, count_enh,
        "Expected #ops = {count_enh}. Got {trans_count}"
    );

    println!("Checking for equality...");
    for (p_llh, p_llh_orig) in llh.iter().zip(llh_orig.iter()) {
        assert_abs_diff_eq!(p_llh.as_slice(), p_llh_orig.as_slice(), epsilon = 1e-4);
    }

    Ok(())
}
