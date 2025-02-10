use geokit::{geodesy::ellipsoid::consts::WGS84, math::Float};
use rand_distr::{Distribution, UnitSphere};

/// `cargo run sample_s2 | cct -I +proj=cart +ellps=WGS84`
///
fn main() {
    let distr = UnitSphere;
    let height_distr = rand::distr::Uniform::new_inclusive(-150, 10_000).unwrap();
    let mut rng = rand::rng();
    for _ in 0..100 {
        let [x, y, z]: [Float; 3] = distr.sample(&mut rng);
        //let lon = y.atan2(x).to_degrees();
        //let lat = z.asin().to_degrees();
        //let h = height_distr.sample(&mut rng);
        println!(
            "{:17.8} {:17.8} {:17.8}",
            x * WGS84.a().m(),
            y * WGS84.a().m(),
            z * WGS84.b().m()
        );
        //let lon = y.atan2(x);
        //let lat = z.asin();
        //let height = height_distr.sample(&mut rand::rng());
        //println!("{} {} {}", lon, lat, height);
    }
}
