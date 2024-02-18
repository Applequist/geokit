use approx::assert_abs_diff_eq;
use geokit::{
    crs::{Crs::Geographic, GeodeticAxes},
    geodesy::geodetic_datum,
    operation::apply_seq,
};
use regex::Regex;
use std::default::Default;

fn dist(a: &[f64], b: &[f64]) -> f64 {
    assert!(a.len() == 3);
    a.iter()
        .zip(b.iter())
        .fold(0.0, |mut acc, e| {
            let diff = e.0 - e.1;
            acc += diff * diff;
            acc
        })
        .sqrt()
}

#[test]
fn llh_to_xyz() {
    let re = Regex::new(r"\s+").unwrap();

    let ll_grid = include_str!("data/ll_grid.txt");

    let ll_orig: Vec<f64> = ll_grid
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .flat_map(|l| re.split(l.trim()).map(|s| s.parse::<f64>().unwrap()))
        .collect();

    let count_ll = ll_grid.lines().count() - 1; // -1 to account for the header line
    println!("Read {count_ll} 2D coordinates");

    let xyz_grid = include_str!("data/xyz_grid.txt");
    let xyz_orig: Vec<f64> = xyz_grid
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .flat_map(|l| re.split(l.trim()).map(|s| s.parse::<f64>().unwrap()))
        .collect();
    let count_xyz = xyz_grid.lines().count() - 1; // -1 to account for the header line
    println!("Read {count_xyz} 3D coordinates");

    assert_eq!(
        count_ll, count_xyz,
        "Expected #ll == #xyz. Got {} != {}",
        count_ll, count_xyz
    );

    let src = Geographic {
        id: "WGS84".into(),
        datum: geodetic_datum::consts::WGS84,
        axes: GeodeticAxes::EastNorth {
            angle_unit: 1.0_f64.to_radians(),
        },
    };
    println!("Source CRS: {:#?}", src);

    let (to_geoc, from_geoc) = src.to_wgs84_geoc();

    // Allocating storage for transformed coordinates.
    println!("Converting ll coordinates to xyz coordinates...");
    let mut xyz = vec![0.; count_ll * 3];
    let trans_count = apply_seq(to_geoc, &ll_orig, &mut xyz).unwrap();
    assert_eq!(
        trans_count, count_ll,
        "Expected #ops = {count_ll}. Got {trans_count}"
    );

    println!("Checking for equality...");
    for (p_xyz, p_xyz_orig) in xyz.chunks_exact(3).zip(xyz_orig.chunks_exact(3)) {
        let err = dist(p_xyz, p_xyz_orig);
        assert!(
            err < 1e-4,
            "Expected dist(actual, expected) < 1e-4. Got {err}"
        );
    }

    println!("Converting xyz coordinates to ll coordinates...");
    let mut ll = vec![0.; count_xyz * 2];
    let trans_count = apply_seq(from_geoc, &xyz_orig, &mut ll).unwrap();
    assert_eq!(
        trans_count, count_xyz,
        "Expected #ops = {count_xyz}. Got {trans_count}"
    );

    println!("Checking for equality...");
    for (p_ll, p_ll_orig) in ll.chunks_exact(2).zip(ll_orig.chunks_exact(2)) {
        assert_abs_diff_eq!(p_ll, p_ll_orig, epsilon = 1e-4);
    }
}
