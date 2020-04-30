#[cfg(feature = "Matrix3")]
use crate::mat3::Mat3;
#[cfg(feature = "Matrix4")]
use crate::mat4::Mat4;
use crate::slice_ops::*;
#[cfg(feature = "Matrix4")]
use crate::vec4::Vec4;
#[cfg(any(feature = "Matrix3", feature = "Matrix4"))]
use crate::vector::MulVectorMatrix;
use crate::vector::Vector;
use std::f32;

pub type Vec3 = [f32; 3];

impl_vector!(Vec3, 3);

#[cfg(feature = "Matrix3")]
impl MulVectorMatrix<Mat3> for Vec3 {
    type VectorType = Vec3;

    fn mul_matrix_left(&self, lhs: &Mat3) -> Self::VectorType {
        let x = self[0];
        let y = self[1];
        let z = self[2];

        [
            lhs[0] * x + lhs[1] * y + lhs[2] * z,
            lhs[3] * x + lhs[4] * y + lhs[5] * z,
            lhs[6] * x + lhs[7] * y + lhs[8] * z,
        ]
    }

    fn mul_matrix(&self, rhs: &Mat3) -> Self::VectorType {
        let x = self[0];
        let y = self[1];
        let z = self[2];

        [
            rhs[0] * x + rhs[3] * y + rhs[6] * z,
            rhs[1] * x + rhs[4] * y + rhs[7] * z,
            rhs[2] * x + rhs[5] * y + rhs[8] * z,
        ]
    }
}

#[cfg(feature = "Matrix4")]
impl MulVectorMatrix<Mat4> for Vec3 {
    type VectorType = Vec4;

    /// Interprets `self` as a column vector with the 4th component equal to 1 and multiplies the given matrix
    /// from the left-hand-side, i.e. `lhs * [...self, 1.0]`
    fn mul_matrix_left(&self, lhs: &Mat4) -> Self::VectorType {
        let x = self[0];
        let y = self[1];
        let z = self[2];
        // let w = 1.0

        [
            lhs[0] * x + lhs[1] * y + lhs[2] * z + lhs[3],
            lhs[4] * x + lhs[5] * y + lhs[6] * z + lhs[7],
            lhs[8] * x + lhs[9] * y + lhs[10] * z + lhs[11],
            lhs[12] * x + lhs[13] * y + lhs[14] * z + lhs[15],
        ]
    }

    /// Interprets `self` as a row vector with the 4th component equal to 1 and multiplies the given matrix
    /// from the right-hand-side, i.e. `[...self, 1.0] * rhs`
    fn mul_matrix(&self, rhs: &Mat4) -> Self::VectorType {
        let x = self[0];
        let y = self[1];
        let z = self[2];
        // let w = 1.0;

        [
            rhs[0] * x + rhs[4] * y + rhs[8] * z + rhs[12],
            rhs[1] * x + rhs[5] * y + rhs[9] * z + rhs[13],
            rhs[2] * x + rhs[6] * y + rhs[10] * z + rhs[14],
            rhs[3] * x + rhs[7] * y + rhs[11] * z + rhs[15],
        ]
    }
}

pub trait CrossProduct {
    fn cross(&self, v: &Vec3) -> Vec3;
}

impl CrossProduct for Vec3 {
    fn cross(&self, v: &Self) -> Self {
        let u = self;

        /*
         *  i  -j   k
         * u[0] u[1] u[2]
         * v[0] v[1] v[2]
         */
        [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::almost_eq;

    #[test]
    #[cfg(feature = "Matrix3")]
    fn vec3_mul_matrix_left() {
        let a = [1., 2., 3., 4., 5., 6., 7., 8., 9.];
        let b = [11., 12., 13.];

        let c = b.mul_matrix_left(&a);
        assert_eq!(c, [74., 182., 290.]);
    }

    #[test]
    #[cfg(feature = "Matrix3")]
    fn vec3_mul_matrix() {
        let a = [1., 2., 3., 4., 5., 6., 7., 8., 9.];
        let b = [11., 12., 13.];

        let c = b.mul_matrix(&a);
        assert_eq!(c, [150., 186., 222.]);
    }
    #[test]
    #[cfg(feature = "Matrix4")]
    fn vec3_mul_matrix4() {
        let a = [
            1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.,
        ];
        let b = [17., 18., 19.];

        let c = b.mul_matrix(&a);
        assert_eq!(c, [291., 346., 401., 456.]);
    }

    #[test]
    #[cfg(feature = "Matrix4")]
    fn vec3_mul_matrix4_left() {
        let a = [
            1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.,
        ];
        let b = [17., 18., 19.];

        let c = b.mul_matrix_left(&a);
        assert_eq!(c, [114., 334., 554., 774.]);
    }

    #[test]
    fn vec3_is_immutable() {
        let b = [2., 3., 4.];
        let c = [3., 2., 3.];

        let _d = b.add(&c);
        assert_eq!(b, [2., 3., 4.]);
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

    #[test]
    fn vec3_cross1() {
        let a = [1., 0., 0.];
        let b = [0., 1., 0.];

        // Right-handed
        assert_eq!(a.cross(&b), [0., 0., 1.]);

        // Left-handed
        assert_eq!(b.cross(&a), [0., 0., -1.]);
    }

    #[test]
    fn vec3_cross2() {
        let a = [2., 3., 4.];
        let b = [5., 6., 7.];

        // Right-handed
        assert_eq!(a.cross(&b), [-3., 6., -3.]);

        // Left-handed
        assert_eq!(b.cross(&a), [3., -6., 3.]);
    }
}
