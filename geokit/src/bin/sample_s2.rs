use clap::{Parser, ValueEnum};
use geokit::math::fp::PI;
use rand::Rng;
use rand_distr::{Distribution, Uniform, UnitSphere};

#[derive(Parser)]
#[command(name = "sample_s2")]
#[command(version = "1.0")]
#[command(about = "Sample points on the sphere", long_about = None)]
struct Cli {
    /// The number of samples.
    count: u32,

    /// The type of sampled coordinates.
    #[arg(value_enum)]
    #[arg(default_value_t = CoordType::XYZ)]
    coord_type: CoordType,

    /// The ellipsoid semi-major axis in meters.
    #[arg(short = 'a', long = "semi-major-axis")]
    #[arg(default_value_t = 6_378_137.)]
    a: f64,

    /// The ellipsoid's inverse flattening.
    #[arg(long)]
    #[arg(default_value_t = 298.257_223_563)]
    inverse_flattening: f64,

    /// Length unit of the output coordinates.
    #[arg(long)]
    #[arg(default_value_t = 1.0)]
    meter_per_unit: f64,

    /// Angle unit of the output coordinates.
    #[arg(long)]
    #[arg(default_value_t = PI / 180.)]
    rad_per_unit: f64,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum CoordType {
    XYZ,
    LonLat,
    LonLatHeight,
}
/// `cargo run sample_s2 | cct -I +proj=cart +ellps=WGS84`
fn main() {
    let cmd = Cli::parse();

    let rng = rand::rng();
    match cmd.coord_type {
        CoordType::XYZ => sample_s2(rng)
            .map(move |xyz| s2_to_xyz(cmd.a, cmd.inverse_flattening, cmd.meter_per_unit, xyz))
            .take(cmd.count as usize)
            .for_each(|xyz| {
                let [x, y, z] = xyz;
                println!("{x:17.8} {y:17.8} {z:17.8}");
            }),
        CoordType::LonLat => sample_s2(rng)
            .map(move |xyz| s2_to_ll(cmd.inverse_flattening, cmd.rad_per_unit, xyz))
            .take(cmd.count as usize)
            .for_each(|ll| {
                let [lon, lat] = ll;
                println!("{lon:17.12} {lat:16.12}");
            }),
        CoordType::LonLatHeight => {
            let h_rng = rand::rng();
            let height_distr = Uniform::new_inclusive(-1500., 30_000.)
                .unwrap()
                .sample_iter(h_rng);
            sample_s2(rng)
                .map(move |xyz| s2_to_ll(cmd.inverse_flattening, cmd.rad_per_unit, xyz))
                .zip(height_distr)
                .take(cmd.count as usize)
                .map(|(ll, h)| [ll[0], ll[1], h])
                .for_each(|llh| {
                    let [lon, lat, h] = llh;
                    println!("{lon:17.12} {lat:16.12} {h:9.3}");
                })
        }
    }

    //let rng = rand::rng();
    //let max_lat = sample_ll_on_s2(WGS84.f(), rng)
    //    .take(5000)
    //    .fold(0.0_f64, |acc, ll| acc.max(ll[1].abs()));
    //println!("{max_lat}");
}

fn sample_s2(rng: impl Rng) -> impl Iterator<Item = [f64; 3]> {
    UnitSphere.sample_iter(rng)
}

fn s2_to_xyz(a: f64, inverse_flattening: f64, meter_per_unit: f64, s2: [f64; 3]) -> [f64; 3] {
    let [x, y, z] = s2;
    [
        a * x / meter_per_unit,
        a * y / meter_per_unit,
        (a - a / inverse_flattening) * z / meter_per_unit,
    ]
}

fn s2_to_ll(inverse_flattening: f64, rad_per_unit: f64, s2: [f64; 3]) -> [f64; 2] {
    let [x, y, z] = s2;
    let lon = y.atan2(x) / rad_per_unit;
    let lat = z.atan2((1. - 1. / inverse_flattening) * ((1. - z) * (1. + z)).sqrt()) / rad_per_unit;
    [lon, lat]
}
