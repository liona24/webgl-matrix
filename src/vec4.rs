#[cfg(feature = "Matrix4")]
use crate::mat4::Mat4;
use crate::slice_ops::*;
#[cfg(feature = "Matrix4")]
use crate::vector::MulVectorMatrix;
use crate::vector::Vector;
use std::f32;

pub type Vec4 = [f32; 4];

impl_vector!(Vec4, 4);

#[cfg(feature = "Matrix4")]
impl MulVectorMatrix<Mat4> for Vec4 {
    type VectorType = Vec4;

    fn mul_matrix_left(&self, lhs: &Mat4) -> Self::VectorType {
        let x = self[0];
        let y = self[1];
        let z = self[2];
        let w = self[3];

        [
            lhs[0] * x + lhs[1] * y + lhs[2] * z + lhs[3] * w,
            lhs[4] * x + lhs[5] * y + lhs[6] * z + lhs[7] * w,
            lhs[8] * x + lhs[9] * y + lhs[10] * z + lhs[11] * w,
            lhs[12] * x + lhs[13] * y + lhs[14] * z + lhs[15] * w,
        ]
    }

    fn mul_matrix(&self, rhs: &Mat4) -> Self::VectorType {
        let x = self[0];
        let y = self[1];
        let z = self[2];
        let w = self[3];

        [
            rhs[0] * x + rhs[4] * y + rhs[8] * z + rhs[12] * w,
            rhs[1] * x + rhs[5] * y + rhs[9] * z + rhs[13] * w,
            rhs[2] * x + rhs[6] * y + rhs[10] * z + rhs[14] * w,
            rhs[3] * x + rhs[7] * y + rhs[11] * z + rhs[15] * w,
        ]
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
