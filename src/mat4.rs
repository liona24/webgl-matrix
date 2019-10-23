use crate::matrix::Matrix;
use crate::utils::EPSILON;
use std::f32;

pub type Mat4 = [f32; 16];

impl Matrix for Mat4 {
    type Type = Self;

    fn zeros() -> Self {
        [0.; 16]
    }
    fn ones() -> Self {
        [1.; 16]
    }
    fn identity() -> Self {
        [
            1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1.,
        ]
    }

    fn copy_to(&self, dst: &mut Self) {
        for i in 0..16 {
            dst[i] = self[i];
        }
    }

    fn transpose(mut self) -> Self {
        let v01 = self[1];
        let v02 = self[2];
        let v03 = self[3];
        let v12 = self[6];
        let v13 = self[7];
        let v23 = self[11];

        self[1] = self[4];
        self[2] = self[8];
        self[3] = self[12];
        self[4] = v01;
        self[6] = self[9];
        self[7] = self[13];
        self[8] = v02;
        self[9] = v12;
        self[11] = self[14];
        self[12] = v03;
        self[13] = v13;
        self[14] = v23;

        self
    }

    fn mul(mut self, rhs: &Self) -> Self {
        let r00 = rhs[0];
        let r01 = rhs[1];
        let r02 = rhs[2];
        let r03 = rhs[3];
        let r10 = rhs[4];
        let r11 = rhs[5];
        let r12 = rhs[6];
        let r13 = rhs[7];
        let r20 = rhs[8];
        let r21 = rhs[9];
        let r22 = rhs[10];
        let r23 = rhs[11];
        let r30 = rhs[12];
        let r31 = rhs[13];
        let r32 = rhs[14];
        let r33 = rhs[15];

        let mut v0 = self[0];
        let mut v1 = self[1];
        let mut v2 = self[2];
        let mut v3 = self[3];
        self[0] = v0 * r00 + v1 * r10 + v2 * r20 + v3 * r30;
        self[1] = v0 * r01 + v1 * r11 + v2 * r21 + v3 * r31;
        self[2] = v0 * r02 + v1 * r12 + v2 * r22 + v3 * r32;
        self[3] = v0 * r03 + v1 * r13 + v2 * r23 + v3 * r33;

        v0 = self[4];
        v1 = self[5];
        v2 = self[6];
        v3 = self[7];
        self[4] = v0 * r00 + v1 * r10 + v2 * r20 + v3 * r30;
        self[5] = v0 * r01 + v1 * r11 + v2 * r21 + v3 * r31;
        self[6] = v0 * r02 + v1 * r12 + v2 * r22 + v3 * r32;
        self[7] = v0 * r03 + v1 * r13 + v2 * r23 + v3 * r33;

        v0 = self[8];
        v1 = self[9];
        v2 = self[10];
        v3 = self[11];
        self[8] = v0 * r00 + v1 * r10 + v2 * r20 + v3 * r30;
        self[9] = v0 * r01 + v1 * r11 + v2 * r21 + v3 * r31;
        self[10] = v0 * r02 + v1 * r12 + v2 * r22 + v3 * r32;
        self[11] = v0 * r03 + v1 * r13 + v2 * r23 + v3 * r33;

        v0 = self[12];
        v1 = self[13];
        v2 = self[14];
        v3 = self[15];
        self[12] = v0 * r00 + v1 * r10 + v2 * r20 + v3 * r30;
        self[13] = v0 * r01 + v1 * r11 + v2 * r21 + v3 * r31;
        self[14] = v0 * r02 + v1 * r12 + v2 * r22 + v3 * r32;
        self[15] = v0 * r03 + v1 * r13 + v2 * r23 + v3 * r33;

        self
    }
    fn add(mut self, rhs: &Self) -> Self {
        for i in 0..16 {
            self[i] += rhs[i];
        }

        self
    }
    fn sub(mut self, rhs: &Self) -> Self {
        for i in 0..16 {
            self[i] -= rhs[i];
        }

        self
    }

    fn scale(mut self, factor: f32) -> Self {
        for i in 0..16 {
            self[i] *= factor;
        }

        self
    }

    fn inverse(mut self) -> Option<Self> {
        let v00 = self[0];
        let v01 = self[1];
        let v02 = self[2];
        let v03 = self[3];
        let v10 = self[4];
        let v11 = self[5];
        let v12 = self[6];
        let v13 = self[7];
        let v20 = self[8];
        let v21 = self[9];
        let v22 = self[10];
        let v23 = self[11];
        let v30 = self[12];
        let v31 = self[13];
        let v32 = self[14];
        let v33 = self[15];

        let tmp00 = v00 * v11 - v01 * v10;
        let tmp01 = v00 * v12 - v02 * v10;
        let tmp02 = v00 * v13 - v03 * v10;
        let tmp03 = v01 * v12 - v02 * v11;
        let tmp04 = v01 * v13 - v03 * v11;
        let tmp05 = v02 * v13 - v03 * v12;
        let tmp06 = v20 * v31 - v21 * v30;
        let tmp07 = v20 * v32 - v22 * v30;
        let tmp08 = v20 * v33 - v23 * v30;
        let tmp09 = v21 * v32 - v22 * v31;
        let tmp10 = v21 * v33 - v23 * v31;
        let tmp11 = v22 * v33 - v23 * v32;

        let det = tmp00 * tmp11 - tmp01 * tmp10 + tmp02 * tmp09 + tmp03 * tmp08 - tmp04 * tmp07
            + tmp05 * tmp06;

        if det.abs() <= EPSILON {
            return None;
        }
        let det_inv = 1.0 / det;

        self[0] = (v11 * tmp11 - v12 * tmp10 + v13 * tmp09) * det_inv;
        self[1] = (v02 * tmp10 - v01 * tmp11 - v03 * tmp09) * det_inv;
        self[2] = (v31 * tmp05 - v32 * tmp04 + v33 * tmp03) * det_inv;
        self[3] = (v22 * tmp04 - v21 * tmp05 - v23 * tmp03) * det_inv;
        self[4] = (v12 * tmp08 - v10 * tmp11 - v13 * tmp07) * det_inv;
        self[5] = (v00 * tmp11 - v02 * tmp08 + v03 * tmp07) * det_inv;
        self[6] = (v32 * tmp02 - v30 * tmp05 - v33 * tmp01) * det_inv;
        self[7] = (v20 * tmp05 - v22 * tmp02 + v23 * tmp01) * det_inv;
        self[8] = (v10 * tmp10 - v11 * tmp08 + v13 * tmp06) * det_inv;
        self[9] = (v01 * tmp08 - v00 * tmp10 - v03 * tmp06) * det_inv;
        self[10] = (v30 * tmp04 - v31 * tmp02 + v33 * tmp00) * det_inv;
        self[11] = (v21 * tmp02 - v20 * tmp04 - v23 * tmp00) * det_inv;
        self[12] = (v11 * tmp07 - v10 * tmp09 - v12 * tmp06) * det_inv;
        self[13] = (v00 * tmp09 - v01 * tmp07 + v02 * tmp06) * det_inv;
        self[14] = (v31 * tmp01 - v30 * tmp03 - v32 * tmp00) * det_inv;
        self[15] = (v20 * tmp03 - v21 * tmp01 + v22 * tmp00) * det_inv;

        Some(self)
    }

    fn det(&self) -> f32 {
        let v00 = self[0];
        let v01 = self[1];
        let v02 = self[2];
        let v03 = self[3];
        let v10 = self[4];
        let v11 = self[5];
        let v12 = self[6];
        let v13 = self[7];
        let v20 = self[8];
        let v21 = self[9];
        let v22 = self[10];
        let v23 = self[11];
        let v30 = self[12];
        let v31 = self[13];
        let v32 = self[14];
        let v33 = self[15];

        let tmp00 = v00 * v11 - v01 * v10;
        let tmp01 = v00 * v12 - v02 * v10;
        let tmp02 = v00 * v13 - v03 * v10;
        let tmp03 = v01 * v12 - v02 * v11;
        let tmp04 = v01 * v13 - v03 * v11;
        let tmp05 = v02 * v13 - v03 * v12;
        let tmp06 = v20 * v31 - v21 * v30;
        let tmp07 = v20 * v32 - v22 * v30;
        let tmp08 = v20 * v33 - v23 * v30;
        let tmp09 = v21 * v32 - v22 * v31;
        let tmp10 = v21 * v33 - v23 * v31;
        let tmp11 = v22 * v33 - v23 * v32;

        tmp00 * tmp11 - tmp01 * tmp10 + tmp02 * tmp09 + tmp03 * tmp08 - tmp04 * tmp07
            + tmp05 * tmp06
    }

    fn adjugate(mut self) -> Self {
        let v00 = self[0];
        let v01 = self[1];
        let v02 = self[2];
        let v03 = self[3];
        let v10 = self[4];
        let v11 = self[5];
        let v12 = self[6];
        let v13 = self[7];
        let v20 = self[8];
        let v21 = self[9];
        let v22 = self[10];
        let v23 = self[11];
        let v30 = self[12];
        let v31 = self[13];
        let v32 = self[14];
        let v33 = self[15];

        self[0] = v11 * (v22 * v33 - v23 * v32) - v21 * (v12 * v33 - v13 * v32)
            + v31 * (v12 * v23 - v13 * v22);
        self[1] = -(v01 * (v22 * v33 - v23 * v32) - v21 * (v02 * v33 - v03 * v32)
            + v31 * (v02 * v23 - v03 * v22));
        self[2] = v01 * (v12 * v33 - v13 * v32) - v11 * (v02 * v33 - v03 * v32)
            + v31 * (v02 * v13 - v03 * v12);
        self[3] = -(v01 * (v12 * v23 - v13 * v22) - v11 * (v02 * v23 - v03 * v22)
            + v21 * (v02 * v13 - v03 * v12));
        self[4] = -(v10 * (v22 * v33 - v23 * v32) - v20 * (v12 * v33 - v13 * v32)
            + v30 * (v12 * v23 - v13 * v22));
        self[5] = v00 * (v22 * v33 - v23 * v32) - v20 * (v02 * v33 - v03 * v32)
            + v30 * (v02 * v23 - v03 * v22);
        self[6] = -(v00 * (v12 * v33 - v13 * v32) - v10 * (v02 * v33 - v03 * v32)
            + v30 * (v02 * v13 - v03 * v12));
        self[7] = v00 * (v12 * v23 - v13 * v22) - v10 * (v02 * v23 - v03 * v22)
            + v20 * (v02 * v13 - v03 * v12);
        self[8] = v10 * (v21 * v33 - v23 * v31) - v20 * (v11 * v33 - v13 * v31)
            + v30 * (v11 * v23 - v13 * v21);
        self[9] = -(v00 * (v21 * v33 - v23 * v31) - v20 * (v01 * v33 - v03 * v31)
            + v30 * (v01 * v23 - v03 * v21));
        self[10] = v00 * (v11 * v33 - v13 * v31) - v10 * (v01 * v33 - v03 * v31)
            + v30 * (v01 * v13 - v03 * v11);
        self[11] = -(v00 * (v11 * v23 - v13 * v21) - v10 * (v01 * v23 - v03 * v21)
            + v20 * (v01 * v13 - v03 * v11));
        self[12] = -(v10 * (v21 * v32 - v22 * v31) - v20 * (v11 * v32 - v12 * v31)
            + v30 * (v11 * v22 - v12 * v21));
        self[13] = v00 * (v21 * v32 - v22 * v31) - v20 * (v01 * v32 - v02 * v31)
            + v30 * (v01 * v22 - v02 * v21);
        self[14] = -(v00 * (v11 * v32 - v12 * v31) - v10 * (v01 * v32 - v02 * v31)
            + v30 * (v01 * v12 - v02 * v11));
        self[15] = v00 * (v11 * v22 - v12 * v21) - v10 * (v01 * v22 - v02 * v21)
            + v20 * (v01 * v12 - v02 * v11);

        self
    }
}

pub trait ProjectionMatrix {
    fn create_perspective(fov_y: f32, aspect_ratio: f32, near: f32, far: f32) -> Mat4;
    fn create_perspective_from_viewport(
        vp_left: f32,
        vp_right: f32,
        vp_top: f32,
        vp_bot: f32,
        near: f32,
        far: f32,
    ) -> Mat4;

    fn create_orthogonal_from_viewport(
        vp_left: f32,
        vp_right: f32,
        vp_top: f32,
        vp_bot: f32,
        near: f32,
        far: f32,
    ) -> Mat4;
}

impl ProjectionMatrix for Mat4 {
    fn create_perspective(fov_y: f32, aspect_ratio: f32, near: f32, far: f32) -> Self {
        let f = 1. / (fov_y / 2.).tan();
        let nf = 1. / (near - far);
        [
            f / aspect_ratio,
            0.,
            0.,
            0.,
            0.,
            f,
            0.,
            0.,
            0.,
            0.,
            (far + near) * nf,
            -1.,
            0.,
            0.,
            2. * far * near * nf,
            0.,
        ]
    }
    fn create_perspective_from_viewport(
        vp_left: f32,
        vp_right: f32,
        vp_top: f32,
        vp_bot: f32,
        near: f32,
        far: f32,
    ) -> Self {
        let wi = 1. / (vp_right - vp_left);
        let hi = 1. / (vp_top - vp_bot);
        let nf = 1. / (near - far);

        [
            near * 2. * wi,
            0.,
            0.,
            0.,
            0.,
            near * 2. * hi,
            0.,
            0.,
            (vp_right + vp_left) * wi,
            (vp_top + vp_bot) * hi,
            (far + near) * nf,
            -1.,
            0.,
            0.,
            far * near * 2. * nf,
            0.,
        ]
    }

    fn create_orthogonal_from_viewport(
        vp_left: f32,
        vp_right: f32,
        vp_top: f32,
        vp_bot: f32,
        near: f32,
        far: f32,
    ) -> Self {
        let wi = 1. / (vp_right - vp_left);
        let hi = 1. / (vp_top - vp_bot);
        let nf = 1. / (near - far);

        [
            2. * wi,
            0.,
            0.,
            0.,
            0.,
            2. * hi,
            0.,
            0.,
            0.,
            0.,
            2. * nf,
            0.,
            -(vp_left + vp_right) * wi,
            -(vp_top + vp_bot) * hi,
            (far + near) * nf,
            1.,
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::almost_eq;

    #[test]
    fn mat4_zeros() {
        let zeros = Mat4::zeros();
        assert!(zeros.iter().all(|&x| x == 0.0));
    }

    #[test]
    fn mat4_ones() {
        let zeros = Mat4::ones();
        assert!(zeros.iter().all(|&x| x == 1.0));
    }

    #[test]
    fn mat4_identity() {
        let i = Mat4::identity();
        assert_eq!(i[0], 1.0);
        assert_eq!(i[1], 0.0);
        assert_eq!(i[2], 0.0);
        assert_eq!(i[3], 0.0);

        assert_eq!(i[4], 0.0);
        assert_eq!(i[5], 1.0);
        assert_eq!(i[6], 0.0);
        assert_eq!(i[7], 0.0);

        assert_eq!(i[8], 0.0);
        assert_eq!(i[9], 0.0);
        assert_eq!(i[10], 1.0);
        assert_eq!(i[11], 0.0);

        assert_eq!(i[12], 0.0);
        assert_eq!(i[13], 0.0);
        assert_eq!(i[14], 0.0);
        assert_eq!(i[15], 1.0);
    }

    #[test]
    fn mat4_copy_to() {
        let mut a = Mat4::zeros();
        let b = Mat4::ones();

        b.copy_to(&mut a);
        assert!(a.iter().all(|&x| x == 1.0));
    }

    #[test]
    fn mat4_transpose() {
        //  1  2  3  4
        //  5  6  7  8
        //  9 10 11 12
        // 13 14 15 16
        //
        // --->
        //
        //  1  5  9 13
        //  2  6 10 14
        //  3  7 11 15
        //  4  8 12 16
        let a = ([
            1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.,
        ] as Mat4)
            .transpose();
        let b = [
            1., 5., 9., 13., 2., 6., 10., 14., 3., 7., 11., 15., 4., 8., 12., 16.,
        ];

        assert_eq!(a, b);
    }

    #[test]
    fn mat4_mul() {
        let a = [
            1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.,
        ];
        let b = [
            11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26.,
        ];

        let c = [
            190., 200., 210., 220., 462., 488., 514., 540., 734., 776., 818., 860., 1006., 1064.,
            1122., 1180.,
        ];

        assert_eq!(a.mul(&b), c);
    }

    #[test]
    fn mat4_mul_identity() {
        let a = [
            1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.,
        ];
        let b = Mat4::identity();

        assert_eq!(b.mul(&a), a);
    }

    #[test]
    fn mat4_add() {
        let a = [
            1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.,
        ];
        let b = [
            11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26.,
        ];

        let c = [
            12., 14., 16., 18., 20., 22., 24., 26., 28., 30., 32., 34., 36., 38., 40., 42.,
        ];

        assert_eq!(a.add(&b), c);
    }

    #[test]
    fn mat4_sub() {
        let a = [
            16., 15., 14., 13., 12., 11., 10., 9., 8., 7., 6., 5., 4., 3., 2., 1.,
        ];
        let b = [
            11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26.,
        ];

        let c = [
            5., 3., 1., -1., -3., -5., -7., -9., -11., -13., -15., -17., -19., -21., -23., -25.,
        ];

        assert_eq!(a.sub(&b), c);
    }

    #[test]
    fn mat4_scale() {
        let a = [
            16., 15., 14., 13., 12., 11., 10., 9., 8., 7., 6., 5., 4., 3., 2., 1.,
        ];
        let b = [
            32., 30., 28., 26., 24., 22., 20., 18., 16., 14., 12., 10., 8., 6., 4., 2.,
        ];

        assert_eq!(a.scale(2.0), b);
    }

    #[test]
    fn mat4_inverse_valid() {
        let a = [
            3., 4., 1., 2., 3., 6., 10., 12., 2., 7., 3., 14., 16., 4., 8., 18.,
        ];
        let b = a.clone();

        let a = a.inverse().expect("Inverse should exist");

        let inv = [
            0.0976342, -0.0401802, -0.0548254, 0.0585805, 0.2482163, 0.0093879, 0.0221555,
            -0.0510702, -0.0168982, 0.1415697, -0.1058956, -0.0101389, -0.1344348, -0.0292903,
            0.0908750, 0.0193391,
        ];

        assert!(almost_eq(&inv, &a));

        assert!(almost_eq(&a.mul(&b), &Mat4::identity()));
    }

    #[test]
    fn mat4_inverse_invalid() {
        let a = [
            1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.,
        ];
        assert_eq!(a.inverse(), None);
    }

    #[test]
    fn mat4_det() {
        let a = [
            3., 4., 1., 2., 3., 6., 10., 12., 2., 7., 3., 14., 16., 4., 8., 18.,
        ];
        assert_eq!(a.det(), -5326.0);
    }

    #[test]
    fn mat4_adjugate() {
        let a = [
            3., 4., 1., 2., 3., 6., 10., 12., 2., 7., 3., 14., 16., 4., 8., 18.,
        ];
        let b = [
            -520., 214., 292., -312., -1322., -50., -118., 272., 90., -754., 564., 54., 716., 156.,
            -484., -103.,
        ];
        assert_eq!(a.adjugate(), b);
    }
}
