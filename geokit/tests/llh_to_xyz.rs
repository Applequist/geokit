use geokit::geodesy::Ellipsoid;
use regex::Regex;
use std::default::Default;

#[test]
fn llh_to_xyz() {
    let ellipsoid = Ellipsoid::default();

    let re = Regex::new(r"\s+").unwrap();
    let ll_grid = include_str!("data/ll_grid.txt");
    let llh_flat: Vec<f64> = ll_grid
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .flat_map(|l| {
            re.split(l.trim())
                .map(|s| s.parse::<f64>().unwrap().to_radians())
                .chain(std::iter::once(0.))
        })
        .collect();

    let xyz_grid = include_str!("data/xyz_grid.txt");
    let xyz_flat: Vec<f64> = xyz_grid
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .flat_map(|l| re.split(l.trim()).map(|s| s.parse::<f64>().unwrap()))
        .collect();

    for (llh, xyz) in llh_flat.chunks(3).zip(xyz_flat.chunks(3)) {
        let mut computed_xyz = [0.; 3];
        ellipsoid.llh_to_xyz(llh, &mut computed_xyz);
        assert!(
            (xyz[0] - computed_xyz[0]).abs() < 1e-3
                && (xyz[1] - computed_xyz[1]).abs() < 1e-3
                && (xyz[2] - computed_xyz[2]).abs() < 1e-3,
            "{:?} -> {:?} != {:?}",
            llh,
            computed_xyz,
            xyz
        );

        let mut computed_llh = [0.; 3];
        ellipsoid.xyz_to_llh(xyz, &mut computed_llh);
        assert!(
            (llh[0] - computed_llh[0]).abs() < 1e-3
                && (llh[1] - computed_llh[1]).abs() < 1e-3
                && (llh[2] - computed_llh[2]).abs() < 1e-3,
            "{:?} -> {:?} != {:?}",
            xyz,
            computed_llh,
            llh
        );
    }
}
