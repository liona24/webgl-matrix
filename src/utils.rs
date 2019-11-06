// we redfine it, because of no_std
// pub const EPSILON : f32 = 1.19209290e-07_f32;
pub const EPSILON: f32 = 1e-5_f32;

/// Checks if two sequences of numbers are equal up to EPSILON precision.
pub fn almost_eq(a: &[f32], b: &[f32]) -> bool {
    if a.len() == b.len() {
        a.iter()
            .zip(b.iter())
            .all(|(&ai, &bi)| ai - EPSILON <= bi && ai + EPSILON >= bi)
    } else {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn almost_eq_if_equal() {
        let a = [1., 2., 3., 4.];
        let b = [1., 2., 3., 4.];

        assert!(almost_eq(&a, &b));
    }

    #[test]
    fn almost_eq_if_not_equal() {
        let a = [1., 2., 3., 4.];
        let b = [1., 2.001, 3., 4.];

        assert!(!almost_eq(&a, &b));

        let a = [1., 2., 3., 4.];
        let b = [1., 2., 2.9999, 4.];

        assert!(!almost_eq(&a, &b));
    }
}
