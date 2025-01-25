use approx::assert_abs_diff_eq;
use geokit::crs::Crs;
use geokit::cs::cartesian::GeocentricAxes;
use geokit::cs::geodetic::GeodeticAxes;
use geokit::geodesy::geodetic_datum;
use geokit::units::angle::DEG;
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
fn llh_to_xyz() -> Result<(), String> {
    let ll_orig = input::read_coords("tests/data/llh_grid.txt", 2).collect::<Vec<_>>();
    let count_ll = ll_orig.len(); // -1 to account for the header line
    println!("Read {count_ll} 2D coordinates");

    let xyz_orig = input::read_coords("tests/data/xyz_grid.txt", 3).collect::<Vec<_>>();
    let count_xyz = xyz_orig.len(); // -1 to account for the header line
    println!("Read {count_xyz} 3D coordinates");

    assert_eq!(
        count_ll, count_xyz,
        "Expected #ll == #xyz. Got {} != {}",
        count_ll, count_xyz
    );

    let src = Crs::Geographic {
        id: "WGS84".into(),
        datum: geodetic_datum::consts::WGS84,
        axes: GeodeticAxes::EastNorth { angle_unit: DEG },
    };
    println!("Source CRS: {:#?}", src);

    let dst = Crs::Geocentric {
        id: "WGS84 Geocentric".into(),
        datum: geodetic_datum::consts::WGS84,
        axes: GeocentricAxes::XYZ,
    };

    //let provider = DefaultTransformationProvider;
    //let (src_to_dst, dst_to_src) = provider.transformation(&src, &dst).unwrap();
    //
    //// Allocating storage for transformed coordinates.
    //println!("Converting ll coordinates to xyz coordinates...");
    //let mut xyz = vec![vec![0.0; 3]; count_ll];
    //let mut trans_count = 0;
    //for (i, o) in ll_orig.iter().zip(xyz.iter_mut()) {
    //    src_to_dst.apply_fwd(i, o)?;
    //    trans_count += 1;
    //}
    //assert_eq!(
    //    trans_count, count_ll,
    //    "Expected #ops = {count_ll}. Got {trans_count}"
    //);
    //
    //println!("Checking for equality...");
    //for (p_xyz, p_xyz_orig) in xyz.iter().zip(xyz_orig.iter()) {
    //    let err = dist(p_xyz, p_xyz_orig);
    //    assert!(
    //        err < 1e-4,
    //        "Expected dist(actual, expected) < 1e-4. Got {err}"
    //    );
    //}
    //
    //println!("Converting xyz coordinates to ll coordinates...");
    //let mut ll = vec![vec![0.0; 2]; count_xyz];
    //let mut trans_count = 0;
    //for (xyz, ll) in xyz_orig.iter().zip(ll.iter_mut()) {
    //    dst_to_src.apply_fwd(xyz, ll)?;
    //    trans_count += 1;
    //}
    //
    //assert_eq!(
    //    trans_count, count_xyz,
    //    "Expected #ops = {count_xyz}. Got {trans_count}"
    //);
    //
    //println!("Checking for equality...");
    //for (p_ll, p_ll_orig) in ll.iter().zip(ll_orig.iter()) {
    //    assert_abs_diff_eq!(p_ll.as_slice(), p_ll_orig.as_slice(), epsilon = 1e-4);
    //}

    Ok(())
}
