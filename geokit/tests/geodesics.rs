use flate2::read::GzDecoder;
use geokit::{
    cs::{
        azimuth::Azimuth,
        geodetic::{Lat, Lon},
    },
    geodesy::{
        ellipsoid::consts::WGS84,
        geodesics::{
            check_direct, karney::KarneyGeodesicSolver, DirectErrors, Geodesic, GeodesicSolver,
        },
    },
    math::fp::Float,
    units::{angle::DEG, length::M},
};
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

fn read_geographiclib_testcases(path: &str) -> impl Iterator<Item = Geodesic> {
    let file = File::open(path).unwrap();
    let reader = BufReader::new(GzDecoder::new(file));
    reader.lines().map(|l| {
        let numbers = l
            .unwrap()
            .split(" ")
            .map(|n| n.parse::<Float>().unwrap())
            .collect::<Vec<_>>();
        println!("Parsed numbers = {:?}", numbers);

        Geodesic {
            p1: (Lon::new(numbers[1] * DEG), Lat::new(numbers[0] * DEG)),
            alpha1: Azimuth::new(numbers[2] * DEG),
            p2: (Lon::new(numbers[4] * DEG), Lat::new(numbers[3] * DEG)),
            alpha2: Azimuth::new(numbers[5] * DEG),
            s: numbers[6] * M,
        }
    })
}

#[test]
fn geographiclib_geodtest_short() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/GeodTest-short.dat.gz"
    );
    let testcases = read_geographiclib_testcases(path);
    let ellipsoid = WGS84;
    let solver = KarneyGeodesicSolver::new(&ellipsoid);
    for tc in testcases {
        let computed_direct = solver.solve_direct(tc.p1, tc.alpha1, tc.s).unwrap();
        check_direct(&computed_direct, &tc, &DirectErrors::default());

        //let computed_inverse = solver.solve_inverse(tc.p1, tc.p2).unwrap();
        //check_inverse(&computed_inverse, &tc, &InverseErrors::default());
    }
}
