use approx::assert_abs_diff_eq;
use geokit::crs::{GeographicCrs, ProjectedCrs};
use geokit::cs::cartesian::ProjectedAxes;
use geokit::cs::geodetic::{GeodeticAxes, Lat, Lon};
use geokit::geodesy::{ellipsoid, prime_meridian, GeodeticDatum};
use geokit::projections::ProjectionSpec;
use geokit::transformations::xyz::XYZIdentity;
use geokit::transformations::{CrsTransformation, CrsXYZTransformation};
use geokit::units::angle::DEG;
use geokit::units::length::M;

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
fn llh_to_utm0() -> Result<(), String> {
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

    let tx = CrsXYZTransformation::new(src, dst, XYZIdentity, XYZIdentity);

    let mut src_pt = vec![-10.0, -90.0, 0.0];
    let mut dst_pt = [0.0; 3];
    tx.src_to_dst(&src_pt, &mut dst_pt).unwrap();
    println!("{src_pt:?} --- src_to_dst ---> {dst_pt:?}");

    tx.dst_to_src(&dst_pt, &mut src_pt).unwrap();
    println!("{dst_pt:?} --- dst_to_src ---> {src_pt:?}");

    // Allocating storage for transformed coordinates.
    println!("Converting llh coordinates to enh coordinates...");
    let mut enh = vec![vec![0.; 3]; count_llh];
    let mut trans_count = 0;
    for (i, o) in llh_orig.iter().zip(enh.iter_mut()) {
        tx.src_to_dst(i, o).unwrap();
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
        tx.dst_to_src(i, o).unwrap();
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
