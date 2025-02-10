use regex::Regex;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub fn read_coords(path: &str) -> impl Iterator<Item = (XYZ, LLH)> {
    println!("{}", std::env::current_dir().unwrap().display());
    let file = File::open(path).unwrap_or_else(|err| panic!("Unable to read {}: {}", path, err));
    let buf_reader = BufReader::new(file);
    let re = Regex::new(r"\s+").unwrap();
    buf_reader
        .lines()
        .filter(|l| {
            let s = l.as_deref().unwrap();
            !s.is_empty() && !s.starts_with("#")
        })
        .map(move |l| {
            l.map(|s| {
                let xyzllh = re
                    .split(s.trim())
                    .map(|s| s.parse::<Float>().unwrap())
                    .take(coord_dim)
                    .collect::<Vec<_>>();
                let xyz = XYZ {
                    x: xyzllh[0] * M,
                    y: xyzllh[1] * M,
                    z: xyzllh[2] * M,
                };
                let llh = LLH {
                    lon: Lon::new(xyzllh[3] * DEG),
                    lat: Lat::new(xyzllh[4] * DEG),
                    heigght: xyzllh[5] * M,
                };
                (xyz, llh)
            })
            .unwrap()
        })
}
