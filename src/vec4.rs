#[cfg(feature = "Matrix4")]
use crate::mat4::Mat4;
#[cfg(feature = "Matrix4")]
use crate::vector::MulVectorMatrix;
use crate::vector::Vector;
use std::f32;

pub type Vec4 = [f32; 4];

impl Vector for Vec4 {
    type VectorType = Vec4;

    fn zeros() -> Self::VectorType {
        [0., 0., 0., 0.]
    }

    fn ones() -> Self::VectorType {
        [1., 1., 1., 1.]
    }

    fn mul(mut self, rhs: &[f32]) -> Self::VectorType {
        debug_assert!(rhs.len() >= 4);

        self[0] *= rhs[0];
        self[1] *= rhs[1];
        self[2] *= rhs[2];
        self[3] *= rhs[3];

        self
    }

    fn add(mut self, rhs: &[f32]) -> Self::VectorType {
        debug_assert!(rhs.len() >= 4);

        self[0] += rhs[0];
        self[1] += rhs[1];
        self[2] += rhs[2];
        self[3] += rhs[3];

        self
    }

    fn sub(mut self, rhs: &[f32]) -> Self::VectorType {
        debug_assert!(rhs.len() >= 4);

        self[0] -= rhs[0];
        self[1] -= rhs[1];
        self[2] -= rhs[2];
        self[3] -= rhs[3];

        self
    }

    fn scale(mut self, factor: f32) -> Self::VectorType {
        self[0] *= factor;
        self[1] *= factor;
        self[2] *= factor;
        self[3] *= factor;

        self
    }

    fn mag(&self) -> f32 {
        self.mag2().sqrt()
    }

    fn mag2(&self) -> f32 {
        let x = self[0];
        let y = self[1];
        let z = self[2];
        let w = self[3];

        x * x + y * y + z * z + w * w
    }

    fn dot(&self, rhs: &[f32]) -> f32 {
        let x1 = self[0];
        let y1 = self[1];
        let z1 = self[2];
        let w1 = self[3];

        let x2 = rhs[0];
        let y2 = rhs[1];
        let z2 = rhs[2];
        let w2 = rhs[3];

        x1 * x2 + y1 * y2 + z1 * z2 + w1 * w2
    }
}

#[cfg(feature = "Matrix4")]
impl MulVectorMatrix for Vec4 {
    type VectorType = Vec4;
    type MatrixType = Mat4;

    fn mul_matrix_left(mut self, lhs: &Self::MatrixType) -> Self::VectorType {
        let x = self[0];
        let y = self[1];
        let z = self[2];
        let w = self[3];

        self[0] = lhs[0] * x + lhs[1] * y + lhs[2] * z + lhs[3] * w;
        self[1] = lhs[4] * x + lhs[5] * y + lhs[6] * z + lhs[7] * w;
        self[2] = lhs[8] * x + lhs[9] * y + lhs[10] * z + lhs[11] * w;
        self[3] = lhs[12] * x + lhs[13] * y + lhs[14] * z + lhs[15] * w;

        self
    }

    fn mul_matrix(mut self, rhs: &Self::MatrixType) -> Self::VectorType {
        let x = self[0];
        let y = self[1];
        let z = self[2];
        let w = self[3];

        self[0] = rhs[0] * x + rhs[4] * y + rhs[8] * z + rhs[12] * w;
        self[1] = rhs[1] * x + rhs[5] * y + rhs[9] * z + rhs[13] * w;
        self[2] = rhs[2] * x + rhs[6] * y + rhs[10] * z + rhs[14] * w;
        self[3] = rhs[3] * x + rhs[7] * y + rhs[11] * z + rhs[15] * w;

        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::almost_eq;

    #[test]
    #[cfg(feature = "Matrix4")]
    fn vec4_mul_matrix_left() {
        let a = [
            1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.,
        ];
        let b = [17., 18., 19., 20.];

        let c = b.mul_matrix_left(&a);
        assert_eq!(c, [190., 486., 782., 1078.]);
    }

    #[test]
    #[cfg(feature = "Matrix4")]
    fn vec4_mul_matrix() {
        let a = [
            1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.,
        ];
        let b = [17., 18., 19., 20.];

        let c = b.mul_matrix(&a);
        assert_eq!(c, [538., 612., 686., 760.]);
    }

    #[test]
    fn vec4_add() {
        let a = [1., 2., 3., 4.];
        let b = [-1., -2., -3., -4.];

        assert_eq!(a.add(&b), [0., 0., 0., 0.]);
    }

    #[test]
    fn vec4_sub() {
        let a = [1., 2., 3., 4.];
        let b = [1., 2., 3., 4.];

        assert_eq!(a.sub(&b), [0., 0., 0., 0.]);
    }

    #[test]
    fn vec4_mul() {
        let a = [1., 2., 3., 4.];
        let b = [2., 3., 4., 5.];

        assert_eq!(a.mul(&b), [2., 6., 12., 20.]);
    }

    #[test]
    fn vec4_scale() {
        let a = [1., 2., 3., 4.];

        assert_eq!(a.scale(3.), [3., 6., 9., 12.]);
    }

    #[test]
    fn vec4_dot() {
        let a = [1., 2., 3., 4.];
        let b = [2., 3., 4., 5.];

        assert_eq!(a.dot(&b), 2. + 6. + 12. + 20.);
    }
    #[test]
    fn vec4_mag() {
        let b = [2., 3., 4., 5.];
        assert!(almost_eq(&[b.mag()], &[7.3484693]));
    }
    #[test]
    fn vec4_mag2() {
        let b = [2., 3., 4., 5.];
        assert!(almost_eq(&[b.mag2()], &[54.]));
    }
}
