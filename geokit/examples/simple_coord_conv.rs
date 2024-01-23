use geokit::{
    crs::{GeodeticAxes, GeographicCrs},
    geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian},
    operation::{
        conversion::{GeogToGeoc, Normalization},
        Bwd, Fwd, Operation,
    },
};
use regex::Regex;

fn main() {
    let re = Regex::new(r"\s+").unwrap();

    let ll_grid = std::fs::read_to_string("../tests/data/ll_grid.txt")
        .expect("Unable to read '../tests/data/ll_grid.txt'");
    let count = ll_grid.lines().count();
    println!("Read {} 2D coordinates", count);

    let ll: Vec<f64> = ll_grid
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .flat_map(|l| re.split(l.trim()).map(|s| s.parse::<f64>().unwrap()))
        .collect();

    let src = GeographicCrs::new(
        "WGS84",
        GeodeticDatum::new(
            "WGS84",
            Ellipsoid::default(),
            PrimeMeridian::default(),
            None,
        ),
        GeodeticAxes::EastNorth {
            angle_unit: 1.0_f64.to_radians(),
        },
    );
    println!("Source CRS: {:#?}", src);

    let norm = Normalization::from(src.axes());
    let geog_to_geoc = GeogToGeoc::new(src.datum());
    let to_geoc = Fwd(norm.clone()).and_then(Fwd(geog_to_geoc));
    let from_geoc = Bwd(geog_to_geoc).and_then(Bwd(norm));

    assert_eq!(
        to_geoc.in_dim(),
        2,
        "Expected to_geoc input dim == 2. Got {}",
        to_geoc.in_dim()
    );
    assert_eq!(
        to_geoc.out_dim(),
        3,
        "Expected to_geoc output dim == 3. Got {}",
        to_geoc.out_dim()
    );
    println!("Transformation from source to destination: {:#?}", to_geoc);

    // Allocating storage for transformed coordinates.
    let mut xyz = vec![0.; count * 3];
    let trans_count = to_geoc.apply_seq(&ll, &mut xyz).unwrap();
    assert_eq!(trans_count, count);

    for xyz in xyz.chunks_exact(3).take(5) {
        println!("{:?}", xyz);
    }

    let mut ll_back = vec![0.; count * 2];
    let trans_count = from_geoc.apply_seq(&xyz, &mut ll_back).unwrap();
    assert_eq!(trans_count, count);

    for ll in ll_back.chunks_exact(2).take(5) {
        println!("{:?}", ll);
    }
}
