#include "test_harness.h"

// Starfield generation is CPU-only and can be tested even in headless builds.
#include "stellar/render/Starfield.h"

#include <algorithm>
#include <cmath>
#include <vector>

static double vlen(const stellar::render::PointVertex& v) {
  const double x = (double)v.px;
  const double y = (double)v.py;
  const double z = (double)v.pz;
  return std::sqrt(x * x + y * y + z * z);
}

static double dist(const stellar::render::PointVertex& a, const stellar::render::PointVertex& b) {
  const double dx = (double)a.px - (double)b.px;
  const double dy = (double)a.py - (double)b.py;
  const double dz = (double)a.pz - (double)b.pz;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
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
  return std::sqrt(std::max(0.0, var)) / mean;
}

int test_starfield() {
  int failures = 0;

  stellar::render::Starfield sf;
  sf.setRadius(1.0);

  stellar::render::Starfield::Settings s;
  s.distribution = stellar::render::Starfield::Distribution::RelaxedFibonacci;
  s.randomRotation = true;
  s.jitter01 = 0.12;
  s.relaxIterations = 6;
  s.relaxStrength = 0.25;
  sf.setSettings(s);

  sf.regenerate(1234567ULL, 512);
  sf.update({0.0, 0.0, 0.0}, 0.0);

  const auto& pts = sf.points();
  CHECK(pts.size() == 512);

  // Unit-length sanity (update uses radius=1).
  double meanLen = 0.0;
  for (const auto& p : pts) meanLen += vlen(p);
  meanLen /= (double)pts.size();
  CHECK(std::abs(meanLen - 1.0) < 1e-3);

  // Nearest-neighbor statistics (should be reasonably uniform; not clumped).
  std::vector<double> nn;
  nn.reserve(pts.size());

  for (std::size_t i = 0; i < pts.size(); ++i) {
    double best = 1.0e9;
    for (std::size_t j = 0; j < pts.size(); ++j) {
      if (i == j) continue;
      best = std::min(best, dist(pts[i], pts[j]));
    }
    nn.push_back(best);
  }

  const double minNN = *std::min_element(nn.begin(), nn.end());
  double meanNN = 0.0;
  for (double x : nn) meanNN += x;
  meanNN /= (double)nn.size();

  const double cvNN = coeffVar(nn);

  // With 512 points, the expected spacing chord length is ~0.156. We just
  // check broad sanity bounds so the test is robust to small algorithm tweaks.
  CHECK(minNN > 0.05);
  CHECK(meanNN > 0.10);
  CHECK(cvNN < 0.40);

  return failures;
}
