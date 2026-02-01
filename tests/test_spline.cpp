#include "test_harness.h"

#include "stellar/math/Spline.h"

#include <cmath>
#include <vector>

static bool approxEq(double a, double b, double eps = 1e-8) {
  return std::abs(a - b) <= eps;
}

static bool nearVec(const stellar::math::Vec3d& a, const stellar::math::Vec3d& b, double eps = 1e-8) {
  return approxEq(a.x, b.x, eps) && approxEq(a.y, b.y, eps) && approxEq(a.z, b.z, eps);
}

static double coeffVar(const std::vector<double>& v) {
  if (v.empty()) return 0.0;
  double mean = 0.0;
  for (double x : v) mean += x;
  mean /= (double)v.size();
  if (!(mean > 1e-12)) return 0.0;

  double var = 0.0;
  for (double x : v) {
    const double d = x - mean;
    var += d * d;
  }
  var /= (double)v.size();
  const double stddev = std::sqrt(std::max(0.0, var));
  return stddev / mean;
}

int test_spline() {
  int failures = 0;

  using stellar::math::Vec3d;
  using stellar::math::ArcLengthTable;
  using stellar::math::buildArcLengthTable;
  using stellar::math::catmullRomCentripetal;
  using stellar::math::catmullRomUniform;
  using stellar::math::reparamByArcLength;

  // A sharp-ish corner / S-curve, chosen to produce noticeable speed variation
  // under the raw parameterization.
  const Vec3d p0{-10.0, 0.0, 0.0};
  const Vec3d p1{0.0, 0.0, 0.0};
  const Vec3d p2{0.0, 10.0, 0.0};
  const Vec3d p3{10.0, 0.0, 0.0};

  // Endpoint interpolation invariants.
  CHECK(nearVec(catmullRomUniform(p0, p1, p2, p3, 0.0), p1));
  CHECK(nearVec(catmullRomUniform(p0, p1, p2, p3, 1.0), p2));
  CHECK(nearVec(catmullRomCentripetal(p0, p1, p2, p3, 0.0), p1));
  CHECK(nearVec(catmullRomCentripetal(p0, p1, p2, p3, 1.0), p2));

  // Arc-length table sanity.
  const ArcLengthTable tbl = buildArcLengthTable(
      [&](double u) { return catmullRomCentripetal(p0, p1, p2, p3, u); },
      96);
  CHECK(tbl.samples >= 2);
  CHECK(tbl.total > 0.0);

  // Constant-speed check: arc-length reparameterization should reduce speed
  // variation when sampling uniformly in time.
  const int steps = 60;

  std::vector<double> dRaw;
  std::vector<double> dArc;
  dRaw.reserve((std::size_t)steps);
  dArc.reserve((std::size_t)steps);

  Vec3d prevRaw = catmullRomCentripetal(p0, p1, p2, p3, 0.0);
  Vec3d prevArc = prevRaw;

  for (int i = 1; i <= steps; ++i) {
    const double uTime = (double)i / (double)steps;

    const Vec3d pRaw = catmullRomCentripetal(p0, p1, p2, p3, uTime);
    dRaw.push_back((pRaw - prevRaw).length());
    prevRaw = pRaw;

    const double uArc = reparamByArcLength(tbl, uTime);
    const Vec3d pA = catmullRomCentripetal(p0, p1, p2, p3, uArc);
    dArc.push_back((pA - prevArc).length());
    prevArc = pA;
  }

  const double cvRaw = coeffVar(dRaw);
  const double cvArc = coeffVar(dArc);

  // Expect arc-length mapping to reduce speed variability.
  CHECK(cvArc < cvRaw);

  // Also expect it to be reasonably uniform (not perfect, it's still a sampled table).
  CHECK(cvArc < 0.35);

  return failures;
}
