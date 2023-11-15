use geokit::geodesy::Ellipsoid;
use itertools::Itertools;
use regex::Regex;
use std::default::Default;

#[test]
fn llh_to_xyz() {
    let ellipsoid = Ellipsoid::default();

    let re = Regex::new(r"\s+").unwrap();
    let ll_grid = include_str!("data/ll_grid.txt");
    let llh_it = ll_grid
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .flat_map(|l| {
            re.split(l.trim())
                .map(|s| s.parse::<f64>().unwrap().to_radians())
                .chain(std::iter::once(0.))
        })
        .tuples::<(_, _, _)>();

    let xyz_grid = include_str!("data/xyz_grid.txt");
    let xyz_it = xyz_grid
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .flat_map(|l| re.split(l.trim()).map(|s| s.parse::<f64>().unwrap()))
        .tuples::<(_, _, _)>();

    for (llh, xyz) in llh_it.zip(xyz_it) {
        let computed_xyz = ellipsoid.llh_to_xyz(&llh);
        assert!(
            (xyz.0 - computed_xyz.0).abs() < 1e-3
                && (xyz.1 - computed_xyz.1).abs() < 1e-3
                && (xyz.2 - computed_xyz.2).abs() < 1e-3,
            "{:?} -> {:?} != {:?}",
            llh,
            computed_xyz,
            xyz
        );

        let computed_llh = ellipsoid.xyz_to_llh(&xyz);
        assert!(
            (llh.0 - computed_llh.0).abs() < 1e-3
                && (llh.1 - computed_llh.1).abs() < 1e-3
                && (llh.2 - computed_llh.2).abs() < 1e-3,
            "{:?} -> {:?} != {:?}",
            xyz,
            computed_llh,
            llh
        );
    }
}
