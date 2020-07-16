webgl-matrix
============

[![Build Status](https://dev.azure.com/ackermlion/webgl-matrix/_apis/build/status/liona24.webgl-matrix?branchName=master)](https://dev.azure.com/ackermlion/webgl-matrix/_build/latest?definitionId=1&branchName=master)
[![codecov](https://codecov.io/gh/liona24/webgl-matrix/branch/master/graph/badge.svg)](https://codecov.io/gh/liona24/webgl-matrix)

## Getting Started

This is a lightweight matrix / vector library meant for usage with WebGL.

At the core this library only consists of the traits `Matrix` and `Vector`.

Available features:
* `Matrix4`: 4x4 matrix operations (includes *Vector4*)
* `Matrix3`: 3x3 matrix operations (includes *Vector3*)
* `Vector4`: 4-dimensional vector operations
* `Vector3`: 3-dimensional vector operations
* `SliceOps`: Low level slice operations such as addition, subtraction, scaling etc.

## Examples

All the types are simple arrays. You may also just use slices as operands.

```rust
use webgl_matrix::{Matrix, Vector, ProjectionMatrix, Mat4, Vec4, Mat3, Vec3};

fn main() {
    // all the default operations available
    let mut B = [1., 2., 3.,
                 4., 5., 6.,
                 7., 8., 9.];
    let b = Vec3::ones();
    // Matrix operations are in-place
    B.inverse();
    B.transpose();
    // ..

    // Some basic vector operations
    let c = B.mul_vector_left(&b);
    let mag = c.mag(); // magnitude
    let d = c.scale(5.);
    let e = c.add(&b);

    // Or fancier transformations
    B.translate(&[1., 2., 3.]);

    let A = Mat4::identity();
    // operate on slices
    let b = [1., 2., 3., 4., 5., 6., 7.];

    // with automatic homogenous coordinate expansion
    let c = A.mul_vector(&b[0..=2]);
    // or using all four coordinates
    let d = A.mul_vector(&b[3..]);

    // create projection matrices (left, right, bot, top, near, far)
    let P = Mat4::create_perspective_from_viewport(0., 1., 0., 1., 0.1, 10.);
}
```
