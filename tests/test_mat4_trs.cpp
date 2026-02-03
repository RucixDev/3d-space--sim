#include "stellar/math/Mat4.h"

#include <cmath>
#include <iostream>

static bool approx(double a, double b, double eps) {
  return std::abs(a - b) <= eps;
}

static bool approxVec3(const stellar::math::Vec3d& a, const stellar::math::Vec3d& b, double eps) {
  return approx(a.x, b.x, eps) && approx(a.y, b.y, eps) && approx(a.z, b.z, eps);
}

static bool approxMat(const stellar::math::Mat4d& a, const stellar::math::Mat4d& b, double eps) {
  for (int i = 0; i < 16; ++i) {
    if (!approx(a.m[i], b.m[i], eps)) return false;
  }
  return true;
}

static bool approxSameRotation(const stellar::math::Quatd& a, const stellar::math::Quatd& b, double minAbsDot) {
  const auto an = a.normalized();
  const auto bn = b.normalized();
  const double d = std::abs(stellar::math::Quatd::dot(an, bn));
  return d >= minAbsDot;
}

int test_mat4_trs() {
  using namespace stellar::math;
  int fails = 0;

  const Vec3d t{1.25, -2.0, 3.5};
  const Quatd q = Quatd::fromEulerYXZ(degToRad(30.0), degToRad(-10.0), degToRad(5.0));

  // rigidTransform should match translation(t) * rotation(q).
  {
    const Mat4d a = Mat4d::translation(t) * Mat4d::rotation(q);
    const Mat4d b = Mat4d::rigidTransform(t, q);

    if (!a.isFinite() || !b.isFinite() || !approxMat(a, b, 1e-12)) {
      std::cerr << "[test_mat4_trs] rigidTransform: matrix mismatch\n";
      ++fails;
    }

    const Vec3d p{3.0, -1.0, 7.0};
    const Vec3d ap = a.transformPointAffine(p);
    const Vec3d bp = b.transformPointAffine(p);
    if (!approxVec3(ap, bp, 1e-12)) {
      std::cerr << "[test_mat4_trs] rigidTransform: point transform mismatch\n";
      ++fails;
    }

    const Vec3d v{0.25, 2.0, -0.5};
    const Vec3d av = a.transformVector(v);
    const Vec3d bv = b.transformVector(v);
    if (!approxVec3(av, bv, 1e-12)) {
      std::cerr << "[test_mat4_trs] rigidTransform: vector transform mismatch\n";
      ++fails;
    }
  }

  // trs should match translation(t) * rotation(q) * scale(s).
  {
    const Vec3d s{2.0, 0.5, 3.0};
    const Mat4d a = Mat4d::translation(t) * Mat4d::rotation(q) * Mat4d::scale(s);
    const Mat4d b = Mat4d::trs(t, q, s);

    if (!a.isFinite() || !b.isFinite()) {
      std::cerr << "[test_mat4_trs] trs: non-finite matrix\n";
      ++fails;
    }

    const Vec3d p{3.0, -1.0, 7.0};
    const Vec3d ap = a.transformPointAffine(p);
    const Vec3d bp = b.transformPointAffine(p);
    if (!approxVec3(ap, bp, 1e-12)) {
      std::cerr << "[test_mat4_trs] trs: point transform mismatch\n";
      ++fails;
    }

    const Vec3d v{-2.0, 1.25, 0.75};
    const Vec3d av = a.transformVector(v);
    const Vec3d bv = b.transformVector(v);
    if (!approxVec3(av, bv, 1e-12)) {
      std::cerr << "[test_mat4_trs] trs: vector transform mismatch\n";
      ++fails;
    }

    if (!approxVec3(b.translationPart(), t, 1e-12)) {
      std::cerr << "[test_mat4_trs] trs: translationPart mismatch\n";
      ++fails;
    }

    const Quatd qr = b.rotationPart();
    if (!approxSameRotation(q, qr, 0.999999)) {
      std::cerr << "[test_mat4_trs] trs: rotationPart mismatch\n";
      ++fails;
    }
  }

  return fails;
}
