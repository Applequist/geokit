use approx::assert_abs_diff_eq;
use geokit::crs::{GeocentricCrs, GeographicCrs};
use geokit::cs::cartesian::GeocentricAxes;
use geokit::cs::geodetic::GeodeticAxes;
use geokit::geodesy::geodetic_datum;
use geokit::transformations::xyz::XYZIdentity;
use geokit::transformations::{CrsTransformation, CrsXYZTransformation};
use geokit::units::angle::DEG;

const XYZLLH: &'static str = include_str!("data/xyz_to_llh_10_000.txt");

fn read_conversions() -> impl Iterator<Item = (XYZ, LLH)> {
    XYZLLH
        .lines()
        .filter(|&l| !l.is_empty() && !l.starts_with("#"))
        .map(move |l| {
            let xyzllh = l
                .split(',')
                .map(|s| s.trim().parse::<Float>().unwrap())
                .collect::<Vec<_>>();
            let xyz = XYZ {
                x: xyzllh[0] * M,
                y: xyzllh[1] * M,
                z: xyzllh[2] * M,
            };
            let llh = LLH {
                lon: Lon::new(xyzllh[3] * DEG),
                lat: Lat::new(xyzllh[4] * DEG),
                height: xyzllh[5] * M,
            };
            (xyz, llh)
        })
}

#[test]
fn xyz_to_llh_10_000() {
    let datum = WGS84;

    for (xyz, expected_llh) in read_conversions() {
        let llh = datum.xyz_to_llh(xyz);
        llh.approx_eq(&expected_llh, GeodeticErrors::default());
    }
}
