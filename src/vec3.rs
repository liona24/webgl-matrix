use crate::mat3::Mat3;
use crate::vector::Vector;
use std::f32;

pub type Vec3 = [f32; 3];

impl Vector for Vec3 {
    type VectorType = Vec3;
    type MatrixType = Mat3;

    fn zeros() -> Self::VectorType {
        [0., 0., 0.]
    }

    fn ones() -> Self::VectorType {
        [1., 1., 1.]
    }

    fn mul(mut self, rhs: &[f32]) -> Self::VectorType {
        debug_assert!(rhs.len() >= 3);

        self[0] *= rhs[0];
        self[1] *= rhs[1];
        self[2] *= rhs[2];

        self
    }

    fn mul_matrix_left(mut self, lhs: &Self::MatrixType) -> Self::VectorType {
        let x = self[0];
        let y = self[1];
        let z = self[2];

        self[0] = lhs[0] * x + lhs[1] * y + lhs[2] * z;
        self[1] = lhs[3] * x + lhs[4] * y + lhs[5] * z;
        self[2] = lhs[6] * x + lhs[7] * y + lhs[8] * z;

        self
    }

    fn mul_matrix(mut self, rhs: &Self::MatrixType) -> Self::VectorType {
        let x = self[0];
        let y = self[1];
        let z = self[2];

        self[0] = rhs[0] * x + rhs[3] * y + rhs[6] * z;
        self[1] = rhs[1] * x + rhs[4] * y + rhs[7] * z;
        self[2] = rhs[2] * x + rhs[5] * y + rhs[8] * z;

        self
    }

    fn add(mut self, rhs: &[f32]) -> Self::VectorType {
        debug_assert!(rhs.len() >= 3);

        self[0] += rhs[0];
        self[1] += rhs[1];
        self[2] += rhs[2];

        self
    }

    fn sub(mut self, rhs: &[f32]) -> Self::VectorType {
        debug_assert!(rhs.len() >= 3);

        self[0] -= rhs[0];
        self[1] -= rhs[1];
        self[2] -= rhs[2];

        self
    }

    fn scale(mut self, factor: f32) -> Self::VectorType {
        self[0] *= factor;
        self[1] *= factor;
        self[2] *= factor;

        self
    }

    fn mag(&self) -> f32 {
        self.mag2().sqrt()
    }

    fn mag2(&self) -> f32 {
        let x = self[0];
        let y = self[1];
        let z = self[2];

        x * x + y * y + z * z
    }

    fn dot(&self, rhs: &[f32]) -> f32 {
        let x1 = self[0];
        let y1 = self[1];
        let z1 = self[2];

        let x2 = rhs[0];
        let y2 = rhs[1];
        let z2 = rhs[2];

        x1 * x2 + y1 * y2 + z1 * z2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::almost_eq;

    #[test]
    fn vec3_mul_matrix_left() {
        let a = [1., 2., 3., 4., 5., 6., 7., 8., 9.];
        let b = [11., 12., 13.];

        let c = b.mul_matrix_left(&a);
        assert_eq!(c, [74., 182., 290.]);
    }

    #[test]
    fn vec3_mul_matrix() {
        let a = [1., 2., 3., 4., 5., 6., 7., 8., 9.];
        let b = [11., 12., 13.];

        let c = b.mul_matrix(&a);
        assert_eq!(c, [150., 186., 222.]);
    }

    #[test]
    fn vec3_add() {
        let a = [1., 2., 3.];
        let b = [-1., -2., -3.];

        assert_eq!(a.add(&b), [0., 0., 0.]);
    }

    #[test]
    fn vec3_sub() {
        let a = [1., 2., 3.];
        let b = [1., 2., 3.];

        assert_eq!(a.sub(&b), [0., 0., 0.]);
    }

    #[test]
    fn vec3_mul() {
        let a = [1., 2., 3.];
        let b = [2., 3., 4.];

        assert_eq!(a.mul(&b), [2., 6., 12.]);
    }

    #[test]
    fn vec3_scale() {
        let a = [1., 2., 3.];

        assert_eq!(a.scale(3.), [3., 6., 9.]);
    }

    #[test]
    fn vec3_dot() {
        let a = [1., 2., 3.];
        let b = [2., 3., 4.];

        assert_eq!(a.dot(&b), 2. + 6. + 12.);
    }
    #[test]
    fn vec3_mag() {
        let b = [2., 3., 4.];
        assert!(almost_eq(&[b.mag()], &[5.385164807]));
    }
    #[test]
    fn vec3_mag2() {
        let b = [2., 3., 4.];
        assert!(almost_eq(&[b.mag2()], &[29.]));
    }
}
