#pragma once

#include "stellar/math/Vec3.h"

#include <cmath>

namespace stellar::math {

// Minimal 4x4 matrix (column-major) for rendering.
// Compatible with OpenGL's expected memory layout.
struct Mat4f {
  float m[16]{}; // column-major: m[col*4 + row]

  static Mat4f identity() {
    Mat4f out{};
    out.m[0] = 1.0f;
    out.m[5] = 1.0f;
    out.m[10] = 1.0f;
    out.m[15] = 1.0f;
    return out;
  }

  static Mat4f perspective(float fovYRad, float aspect, float zNear, float zFar) {
    const float f = 1.0f / std::tan(fovYRad * 0.5f);

    Mat4f out{};
    out.m[0] = f / aspect;
    out.m[5] = f;
    out.m[10] = (zFar + zNear) / (zNear - zFar);
    out.m[11] = -1.0f;
    out.m[14] = (2.0f * zFar * zNear) / (zNear - zFar);
    return out;
  }

  static Mat4f lookAt(const Vec3d& eye, const Vec3d& center, const Vec3d& up) {
    const Vec3d f = (center - eye).normalized();
    const Vec3d s = cross(f, up).normalized();
    const Vec3d u = cross(s, f);

    Mat4f out = Mat4f::identity();

    out.m[0] = static_cast<float>(s.x);
    out.m[1] = static_cast<float>(u.x);
    out.m[2] = static_cast<float>(-f.x);

    out.m[4] = static_cast<float>(s.y);
    out.m[5] = static_cast<float>(u.y);
    out.m[6] = static_cast<float>(-f.y);

    out.m[8] = static_cast<float>(s.z);
    out.m[9] = static_cast<float>(u.z);
    out.m[10] = static_cast<float>(-f.z);

    out.m[12] = static_cast<float>(-dot(s, eye));
    out.m[13] = static_cast<float>(-dot(u, eye));
    out.m[14] = static_cast<float>(dot(f, eye));

    return out;
  }

  Mat4f operator*(const Mat4f& rhs) const {
    Mat4f out{};
    for (int col = 0; col < 4; ++col) {
      for (int row = 0; row < 4; ++row) {
        float sum = 0.0f;
        for (int k = 0; k < 4; ++k) {
          const float a = m[k * 4 + row];     // (row, k)
          const float b = rhs.m[col * 4 + k]; // (k, col)
          sum += a * b;
        }
        out.m[col * 4 + row] = sum;
      }
    }
    return out;
  }
};

} // namespace stellar::math
