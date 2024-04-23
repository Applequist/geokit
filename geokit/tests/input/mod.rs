use core::f64;
use regex::Regex;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub fn read_coords(path: &str, coord_dim: usize) -> impl Iterator<Item = Vec<f64>> {
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
                re.split(s.trim())
                    .map(|s| s.parse::<f64>().unwrap())
                    .take(coord_dim)
                    .collect::<Vec<_>>()
            })
            .unwrap()
        })
}
