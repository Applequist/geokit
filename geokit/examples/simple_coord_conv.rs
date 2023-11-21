use geokit::{
    crs::geodetic::{GeodeticAxes, GeodeticCrs},
    geodesy::{Ellipsoid, GeodeticDatum, PrimeMeridian},
    id::Id,
    transformation::Transformation,
};
use regex::Regex;

fn main() {
    let re = Regex::new(r"\s+").unwrap();

    let ll_grid = std::fs::read_to_string("../tests/data/ll_grid.txt").unwrap();
    let count = ll_grid.lines().count();
    println!("Read {} 2D coordinates", count);

    let ll: Vec<f64> = ll_grid
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .flat_map(|l| re.split(l.trim()).map(|s| s.parse::<f64>().unwrap()))
        .collect();

    let src = GeodeticCrs::new(
        Id::name("WGS84"),
        GeodeticDatum::new(
            Id::name("WGS84"),
            Ellipsoid::default(),
            PrimeMeridian::default(),
            None,
        ),
        GeodeticAxes::EastNorth {
            angle_unit: 1.0_f64.to_radians(),
        },
    );
    println!("Source CRS: {:#?}", src);

    let (dst, to_geoc, _from_geoc) = src.geoc_crs();
    println!("Destination CRS: {:#?}", dst);

    assert_eq!(to_geoc.in_dim(), 2);
    assert_eq!(to_geoc.out_dim(), 3);
    println!("Transformation from source to destination: {:#?}", to_geoc);

    // Allocating storage for transformed coordinates.
    let mut xyz = vec![0.; count * 3];
    let trans_count = to_geoc.apply_seq(&ll, &mut xyz).unwrap();

    assert_eq!(trans_count, count);

    for xyz in xyz.chunks_exact(3).take(5) {
        println!("{:?}", xyz);
    }
}
