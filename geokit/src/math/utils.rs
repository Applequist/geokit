use super::Float;

pub fn remainder(val: Float, r: Float) -> Float {
    val - (val / r).round_ties_even() * r
}

pub fn iter_fn<const D: usize>(
    init: [Float; D],
    deltas: &dyn Fn([Float; D]) -> [Float; D],
    epsilon: [Float; D],
) -> [Float; D] {
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
        use super::remainder;
        assert_eq!(remainder(-180.0, 360.0), -180.0);
        assert_eq!(remainder(0.0, 360.0), 0.0);
        assert_eq!(remainder(180.0, 360.0), 180.0);
        assert_eq!(remainder(360.0, 360.0), 0.0);
    }
}
