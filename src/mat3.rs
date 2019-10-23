use crate::matrix::Matrix;
use crate::utils::EPSILON;

pub type Mat3 = [f32; 9];
pub type Vec3 = [f32; 3];

impl Matrix for Mat3 {
    type MatrixType = Mat3;
    type VectorType = Vec3;

    fn zeros() -> Self {
        [0., 0., 0., 0., 0., 0., 0., 0., 0.]
    }
    fn ones() -> Self {
        [1., 1., 1., 1., 1., 1., 1., 1., 1.]
    }
    fn identity() -> Self {
        [1., 0., 0., 0., 1., 0., 0., 0., 1.]
    }

    fn copy_to(&self, dst: &mut Self) {
        dst[0] = self[0];
        dst[1] = self[1];
        dst[2] = self[2];
        dst[3] = self[3];
        dst[4] = self[4];
        dst[5] = self[5];
        dst[6] = self[6];
        dst[7] = self[7];
        dst[8] = self[8];
    }

    fn transpose(mut self) -> Self {
        let v01 = self[1];
        let v02 = self[2];
        let v12 = self[5];
        self[1] = self[3];
        self[2] = self[6];
        self[3] = v01;
        self[5] = self[7];
        self[6] = v02;
        self[7] = v12;

        self
    }

    fn mul(mut self, rhs: &Self) -> Self {
        let lhs00 = self[0];
        let lhs01 = self[1];
        let lhs02 = self[2];
        let lhs10 = self[3];
        let lhs11 = self[4];
        let lhs12 = self[5];
        let lhs20 = self[6];
        let lhs21 = self[7];
        let lhs22 = self[8];

        let rhs00 = rhs[0];
        let rhs01 = rhs[1];
        let rhs02 = rhs[2];
        let rhs10 = rhs[3];
        let rhs11 = rhs[4];
        let rhs12 = rhs[5];
        let rhs20 = rhs[6];
        let rhs21 = rhs[7];
        let rhs22 = rhs[8];

        self[0] = lhs00 * rhs00 + lhs01 * rhs10 + lhs02 * rhs20;
        self[1] = lhs00 * rhs01 + lhs01 * rhs11 + lhs02 * rhs21;
        self[2] = lhs00 * rhs02 + lhs01 * rhs12 + lhs02 * rhs22;
        self[3] = lhs10 * rhs00 + lhs11 * rhs10 + lhs12 * rhs20;
        self[4] = lhs10 * rhs01 + lhs11 * rhs11 + lhs12 * rhs21;
        self[5] = lhs10 * rhs02 + lhs11 * rhs12 + lhs12 * rhs22;
        self[6] = lhs20 * rhs00 + lhs21 * rhs10 + lhs22 * rhs20;
        self[7] = lhs20 * rhs01 + lhs21 * rhs11 + lhs22 * rhs21;
        self[8] = lhs20 * rhs02 + lhs21 * rhs12 + lhs22 * rhs22;

        self
    }

    fn mul_vector(&self, rhs: &Vec3) -> Vec3 {
        let x = rhs[0];
        let y = rhs[1];
        let w = rhs[2];
        [
            self[0] * x + self[1] * y + self[2] * w,
            self[3] * x + self[4] * y + self[5] * w,
            self[6] * x + self[7] * y + self[8] * w,
        ]
    }
    fn mul_vector_left(&self, lhs: &Vec3) -> Vec3 {
        let x = lhs[0];
        let y = lhs[1];
        let w = lhs[2];
        [
            self[0] * x + self[3] * y + self[6] * w,
            self[1] * x + self[4] * y + self[7] * w,
            self[2] * x + self[5] * y + self[8] * w,
        ]
    }

    fn add(mut self, rhs: &Self) -> Self {
        self[0] += rhs[0];
        self[1] += rhs[1];
        self[2] += rhs[2];
        self[3] += rhs[3];
        self[4] += rhs[4];
        self[5] += rhs[5];
        self[6] += rhs[6];
        self[7] += rhs[7];
        self[8] += rhs[8];

        self
    }

    fn sub(mut self, rhs: &Self) -> Self {
        self[0] -= rhs[0];
        self[1] -= rhs[1];
        self[2] -= rhs[2];
        self[3] -= rhs[3];
        self[4] -= rhs[4];
        self[5] -= rhs[5];
        self[6] -= rhs[6];
        self[7] -= rhs[7];
        self[8] -= rhs[8];

        self
    }

    fn scale(mut self, factor: f32) -> Self {
        self[0] *= factor;
        self[1] *= factor;
        self[2] *= factor;
        self[3] *= factor;
        self[4] *= factor;
        self[5] *= factor;
        self[6] *= factor;
        self[7] *= factor;
        self[8] *= factor;

        self
    }

    fn inverse(mut self) -> Option<Self> {
        let v00 = self[0];
        let v01 = self[1];
        let v02 = self[2];
        let v10 = self[3];
        let v11 = self[4];
        let v12 = self[5];
        let v20 = self[6];
        let v21 = self[7];
        let v22 = self[8];

        let tmp01 = v22 * v11 - v12 * v21;
        let tmp11 = -v22 * v10 + v12 * v20;
        let tmp21 = v21 * v10 - v11 * v20;

        let det = v00 * tmp01 + v01 * tmp11 + v02 * tmp21;

        if det.abs() <= EPSILON {
            return None;
        }

        let det_inv = 1.0 / det;

        self[0] = tmp01 * det_inv;
        self[1] = (-v22 * v01 + v02 * v21) * det_inv;
        self[2] = (v12 * v01 - v02 * v11) * det_inv;
        self[3] = tmp11 * det_inv;
        self[4] = (v22 * v00 - v02 * v20) * det_inv;
        self[5] = (-v12 * v00 + v02 * v10) * det_inv;
        self[6] = tmp21 * det_inv;
        self[7] = (-v21 * v00 + v01 * v20) * det_inv;
        self[8] = (v11 * v00 - v01 * v10) * det_inv;

        Some(self)
    }

    fn det(&self) -> f32 {
        let v00 = self[0];
        let v01 = self[1];
        let v02 = self[2];
        let v10 = self[3];
        let v11 = self[4];
        let v12 = self[5];
        let v20 = self[6];
        let v21 = self[7];
        let v22 = self[8];

        v00 * (v22 * v11 - v12 * v21)
            + v01 * (-v22 * v10 + v12 * v20)
            + v02 * (v21 * v10 - v11 * v20)
    }

    fn adjugate(mut self) -> Self {
        let v00 = self[0];
        let v01 = self[1];
        let v02 = self[2];
        let v10 = self[3];
        let v11 = self[4];
        let v12 = self[5];
        let v20 = self[6];
        let v21 = self[7];
        let v22 = self[8];

        self[0] = v11 * v22 - v12 * v21;
        self[1] = v02 * v21 - v01 * v22;
        self[2] = v01 * v12 - v02 * v11;
        self[3] = v12 * v20 - v10 * v22;
        self[4] = v00 * v22 - v02 * v20;
        self[5] = v02 * v10 - v00 * v12;
        self[6] = v10 * v21 - v11 * v20;
        self[7] = v01 * v20 - v00 * v21;
        self[8] = v00 * v11 - v01 * v10;

        self
    }

    fn translate(mut self, direction: &Vec3) -> Self {
        let x = direction[0] / direction[2];
        let y = direction[1] / direction[2];

        self[6] += x * self[0] + y * self[3];
        self[7] += x * self[1] + y * self[4];
        self[8] += x * self[2] + y * self[5];

        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::almost_eq;

    #[test]
    fn mat3_zeros() {
        let zeros = Mat3::zeros();
        assert!(zeros.iter().all(|&x| x == 0.0));
    }

    #[test]
    fn mat3_ones() {
        let ones = Mat3::ones();
        assert!(ones.iter().all(|&x| x == 1.0));
    }

    #[test]
    fn mat3_identity() {
        let i = Mat3::identity();
        assert_eq!(i[0], 1.0);
        assert_eq!(i[1], 0.0);
        assert_eq!(i[2], 0.0);

        assert_eq!(i[3], 0.0);
        assert_eq!(i[4], 1.0);
        assert_eq!(i[5], 0.0);

        assert_eq!(i[6], 0.0);
        assert_eq!(i[7], 0.0);
        assert_eq!(i[8], 1.0);
    }

    #[test]
    fn mat3_copy_to() {
        let mut a = Mat3::zeros();
        let b = Mat3::ones();

        b.copy_to(&mut a);
        assert!(a.iter().all(|&x| x == 1.0));
    }

    #[test]
    fn mat3_transpose() {
        // 1 2 3
        // 4 5 6
        // 7 8 9
        //
        // --->
        //
        // 1 4 7
        // 2 5 8
        // 3 6 9
        let a = ([1., 2., 3., 4., 5., 6., 7., 8., 9.] as Mat3).transpose();
        let b: Mat3 = [1., 4., 7., 2., 5., 8., 3., 6., 9.];

        assert_eq!(a, b);
    }

    #[test]
    fn mat3_mul() {
        let a: Mat3 = [1., 2., 3., 4., 5., 6., 7., 8., 9.];
        let b: Mat3 = [11., 12., 13., 14., 15., 16., 17., 18., 19.];

        let c: Mat3 = [90., 96., 102., 216., 231., 246., 342., 366., 390.];

        assert_eq!(a.mul(&b), c);
    }

    #[test]
    fn mat3_mul_identity() {
        let a = [1., 2., 3., 4., 5., 6., 7., 8., 9.];
        let b = Mat3::identity();

        assert_eq!(b.mul(&a), a);
    }

    #[test]
    fn mat3_add() {
        let a = [1., 2., 3., 4., 5., 6., 7., 8., 9.];
        let b = [11., 12., 13., 14., 15., 16., 17., 18., 19.];

        let c = [12., 14., 16., 18., 20., 22., 24., 26., 28.];

        assert_eq!(a.add(&b), c);
    }

    #[test]
    fn mat3_sub() {
        let a = [9., 8., 7., 6., 5., 4., 3., 2., 1.];
        let b = [11., 12., 13., 14., 15., 16., 17., 18., 19.];

        let c = [-2., -4., -6., -8., -10., -12., -14., -16., -18.];

        assert_eq!(a.sub(&b), c);
    }

    #[test]
    fn mat3_scale() {
        let a = [9., 8., 7., 6., 5., 4., 3., 2., 1.];
        let b = [18., 16., 14., 12., 10., 8., 6., 4., 2.];

        assert_eq!(a.scale(2.0), b);
    }

    #[test]
    fn mat3_inverse_valid() {
        let a = [1., 3., 2., 4., 2., 8., 9., 2., 7.];
        let b = a.clone();

        let a = a.inverse().expect("Inverse should exist");

        let inv = [
            -0.01818, -0.15455, 0.18182, 0.4, -0.1, 0.0, -0.09091, 0.22727, -0.09091,
        ];
        assert!(almost_eq(&inv, &a));

        assert!(almost_eq(&a.mul(&b), &Mat3::identity()));
    }

    #[test]
    fn mat3_inverse_invalid() {
        let a = [9., 8., 7., 6., 5., 4., 3., 2., 1.];
        assert_eq!(a.inverse(), None);
    }

    #[test]
    fn mat3_det() {
        let a = [1., 3., 2., 4., 2., 8., 9., 2., 7.];
        assert_eq!(a.det(), 110.0);
    }

    #[test]
    fn mat3_adjugate() {
        let a = [1., 3., 2., 4., 2., 8., 9., 2., 7.];
        let b = [-2., -17., 20., 44., -11., 0., -10., 25., -10.];
        assert_eq!(a.adjugate(), b);
    }

    #[test]
    fn mat3_mul_vector() {
        let a = [1., 2., 3., 4., 5., 6., 7., 8., 9.];
        let b = [11., 12., 13.];

        let c = a.mul_vector(&b);
        assert_eq!(c, [74., 182., 290.]);
    }

    #[test]
    fn mat3_mul_vector_left() {
        let a = [1., 2., 3., 4., 5., 6., 7., 8., 9.];
        let b = [11., 12., 13.];

        let c = a.mul_vector_left(&b);
        assert_eq!(c, [150., 186., 222.]);
    }

    #[test]
    fn mat3_translate() {
        let d = [3., -5., 1.];
        let m = Mat3::identity().translate(&d);

        let a = [-3., 5., 1.];
        assert_eq!(m.mul_vector_left(&a), [0., 0., 1.]);
    }
}
