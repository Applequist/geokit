use super::fp::Float;

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
