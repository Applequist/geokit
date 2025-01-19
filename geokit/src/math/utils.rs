/// Wrap 'val' into [-r.abs(), r.abs()]
pub fn wrap(val: f64, r: f64) -> f64 {
    debug_assert!(r != 0.0, "Expected r != 0. Got {}", r);
    let mut x = val;
    let r_abs = r.abs();
    while x > r_abs {
        x -= 2. * r_abs;
    }
    while x < -r_abs {
        x += 2. * r_abs;
    }
    x
}

pub fn iter_fn<const D: usize>(
    init: [f64; D],
    deltas: &dyn Fn([f64; D]) -> [f64; D],
    epsilon: [f64; D],
) -> [f64; D] {
    let mut res = init;
    loop {
        let ds = deltas(res);
        for i in 0..D {
            res[i] += ds[i];
        }
        if ds.iter().zip(epsilon.iter()).all(|(d, e)| d.abs() < *e) {
            break;
        }
    }
    res
}

#[cfg(test)]
mod test {
    #[test]
    fn wrap() {
        use super::wrap;
        assert_eq!(wrap(-180.0, 180.0), -180.0);
        assert_eq!(wrap(0.0, 180.0), 0.0);
        assert_eq!(wrap(180.0, 180.0), 180.0);
        assert_eq!(wrap(360.0, 180.0), 0.0);
    }
}
