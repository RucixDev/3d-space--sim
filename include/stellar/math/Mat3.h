#pragma once

#include "stellar/math/Vec3.h"

#include <algorithm>
#include <cmath>

namespace stellar::math {

// Column-major 3x3 matrix (OpenGL-friendly).
//
// Indexing matches Mat4d:
//   m[col * 3 + row]
//
// This is intentionally minimal and header-only to keep it usable across the
// project (sim/math/tools) without introducing heavy dependencies.
struct Mat3d {
  double m[9]{};

  constexpr Mat3d() = default;

  bool isFinite() const {
    for (double v : m) {
      if (!std::isfinite(v)) return false;
    }
    return true;
  }

  static Mat3d zero() { return Mat3d{}; }

  static Mat3d identity() {
    Mat3d r{};
    r.m[0] = r.m[4] = r.m[8] = 1.0;
    return r;
  }

  static Mat3d fromColumns(const Vec3d& c0, const Vec3d& c1, const Vec3d& c2) {
    Mat3d r{};
    r.m[0] = c0.x;
    r.m[1] = c0.y;
    r.m[2] = c0.z;

    r.m[3] = c1.x;
    r.m[4] = c1.y;
    r.m[5] = c1.z;

    r.m[6] = c2.x;
    r.m[7] = c2.y;
    r.m[8] = c2.z;
    return r;
  }

  // Outer product: a * b^T.
  //
  // Element form:
  //   M_ij = a_i * b_j
  //
  // In column-major storage: column j is a * b_j.
  static Mat3d outerProduct(const Vec3d& a, const Vec3d& b) {
    return fromColumns(a * b.x, a * b.y, a * b.z);
  }

  double trace() const { return m[0] + m[4] + m[8]; }

  // Matrix-vector product (treating v as a column vector).
  Vec3d operator*(const Vec3d& v) const {
    return {
      m[0] * v.x + m[3] * v.y + m[6] * v.z,
      m[1] * v.x + m[4] * v.y + m[7] * v.z,
      m[2] * v.x + m[5] * v.y + m[8] * v.z,
    };
  }

  Mat3d operator+(const Mat3d& o) const {
    Mat3d r{};
    for (int i = 0; i < 9; ++i) r.m[i] = m[i] + o.m[i];
    return r;
  }

  Mat3d operator-(const Mat3d& o) const {
    Mat3d r{};
    for (int i = 0; i < 9; ++i) r.m[i] = m[i] - o.m[i];
    return r;
  }

  Mat3d operator*(double s) const {
    Mat3d r{};
    for (int i = 0; i < 9; ++i) r.m[i] = m[i] * s;
    return r;
  }

  Mat3d& operator+=(const Mat3d& o) {
    for (int i = 0; i < 9; ++i) m[i] += o.m[i];
    return *this;
  }

  Mat3d& operator-=(const Mat3d& o) {
    for (int i = 0; i < 9; ++i) m[i] -= o.m[i];
    return *this;
  }

  Mat3d& operator*=(double s) {
    for (int i = 0; i < 9; ++i) m[i] *= s;
    return *this;
  }

  void addToDiagonal(double s) {
    m[0] += s;
    m[4] += s;
    m[8] += s;
  }

  double determinant() const {
    const double a00 = m[0], a01 = m[3], a02 = m[6];
    const double a10 = m[1], a11 = m[4], a12 = m[7];
    const double a20 = m[2], a21 = m[5], a22 = m[8];

    return a00 * (a11 * a22 - a12 * a21) -
           a01 * (a10 * a22 - a12 * a20) +
           a02 * (a10 * a21 - a11 * a20);
  }

  // Compute the inverse matrix.
  //
  // Returns false if:
  //  - the matrix is non-finite
  //  - the determinant is too close to zero (|det| <= detEps)
  bool inverse(Mat3d& out, double detEps = 1e-18) const {
    if (!isFinite()) return false;

    const double a00 = m[0], a01 = m[3], a02 = m[6];
    const double a10 = m[1], a11 = m[4], a12 = m[7];
    const double a20 = m[2], a21 = m[5], a22 = m[8];

    const double c00 = (a11 * a22 - a12 * a21);
    const double c01 = -(a10 * a22 - a12 * a20);
    const double c02 = (a10 * a21 - a11 * a20);

    const double c10 = -(a01 * a22 - a02 * a21);
    const double c11 = (a00 * a22 - a02 * a20);
    const double c12 = -(a00 * a21 - a01 * a20);

    const double c20 = (a01 * a12 - a02 * a11);
    const double c21 = -(a00 * a12 - a02 * a10);
    const double c22 = (a00 * a11 - a01 * a10);

    const double det = a00 * c00 + a01 * c01 + a02 * c02;

    if (!std::isfinite(det) || std::abs(det) <= detEps) return false;

    const double invDet = 1.0 / det;

    // Inverse = adjugate / det. Adjugate = cofactor^T.
    const double inv00 = c00 * invDet;
    const double inv01 = c10 * invDet;
    const double inv02 = c20 * invDet;

    const double inv10 = c01 * invDet;
    const double inv11 = c11 * invDet;
    const double inv12 = c21 * invDet;

    const double inv20 = c02 * invDet;
    const double inv21 = c12 * invDet;
    const double inv22 = c22 * invDet;

    out = Mat3d{};
    out.m[0] = inv00;
    out.m[3] = inv01;
    out.m[6] = inv02;

    out.m[1] = inv10;
    out.m[4] = inv11;
    out.m[7] = inv12;

    out.m[2] = inv20;
    out.m[5] = inv21;
    out.m[8] = inv22;
    return true;
  }

  // Solve A*x = b for x using matrix inversion.
  bool solve(Vec3d& outX, const Vec3d& b, double detEps = 1e-18) const {
    Mat3d inv{};
    if (!inverse(inv, detEps)) return false;
    outX = inv * b;
    return true;
  }
};

inline Mat3d operator*(double s, const Mat3d& M) { return M * s; }

} // namespace stellar::math
