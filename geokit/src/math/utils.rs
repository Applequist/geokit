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
