use geokit::{
    cs::{
        cartesian::{
            geocentric::{approx_eq_xyz, XYZ},
            CartesianTolerance,
        },
        geodetic::{approx_eq_llh, GeodeticTolerance, Lat, Lon, LLH},
    },
    geodesy::geodetic_datum::consts::WGS84,
    math::fp::Float,
    units::{angle::DEG, length::M},
};

const XYZLLH: &'static str = include_str!("data/xyz_to_llh_10_000.txt");

fn read_conversions() -> impl Iterator<Item = (XYZ, LLH)> {
    XYZLLH
        .lines()
        .filter(|&l| !l.is_empty() && !l.starts_with("#"))
        .map(|l| {
            let xyzllh = l
                .split(" ")
                .filter(|s| !s.is_empty())
                .map(|s| s.parse::<Float>().unwrap())
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

    for (xyz, llh) in read_conversions() {
        let computed_xyz = datum.llh_to_xyz(llh);
        assert!(approx_eq_xyz(
            computed_xyz,
            xyz,
            CartesianTolerance::default()
        ));
        let computed_llh = datum.xyz_to_llh(xyz);
        assert!(approx_eq_llh(
            computed_llh,
            llh,
            GeodeticTolerance::default()
        ));
    }
}
