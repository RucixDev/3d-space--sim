#include "stellar/sim/BearingTrack.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

static bool approxVec3(const math::Vec3d& a, const math::Vec3d& b, double eps = 1e-6) {
  return approx(a.x, b.x, eps) && approx(a.y, b.y, eps) && approx(a.z, b.z, eps);
}

int test_geometry_bearing_triangulation() {
  int fails = 0;

  // --- Two bearings with baseline should triangulate a unique point ---
  {
    sim::BearingTrack3d tr{};
    sim::BearingTrackParams p{};
    p.observationHalfLifeSec = 0.0; // no forgetting in this deterministic test
    p.solveRegularization = 0.0;
    p.determinantEps = 1.0e-14;
    p.minEffectiveWeight = 1.0; // allow solve as soon as geometry is solvable
    p.sigmaMinKm = 0.0;
    p.sigmaMaxKm = 1.0e9;

    const math::Vec3d target{100.0, -50.0, 20.0};

    const math::Vec3d o0{0.0, 0.0, 0.0};
    const math::Vec3d d0 = (target - o0).normalized();

    const auto r0 = sim::updateBearingTrack(tr, 0.0, true, o0, d0, 1.0, p);
    if (r0.valid) {
      std::cerr << "[test_geometry_bearing_triangulation] single bearing unexpectedly produced a valid solution.\n";
      ++fails;
    }

    const math::Vec3d o1{50.0, 20.0, -10.0};
    const math::Vec3d d1 = (target - o1).normalized();

    const auto r1 = sim::updateBearingTrack(tr, 1.0, true, o1, d1, 1.0, p);

    if (!r1.valid) {
      std::cerr << "[test_geometry_bearing_triangulation] triangulation did not become valid after two bearings.\n";
      ++fails;
    } else if (!approxVec3(r1.posKm, target, 1e-6)) {
      std::cerr << "[test_geometry_bearing_triangulation] bad triangulated position. got=(" << r1.posKm.x << "," << r1.posKm.y
                << "," << r1.posKm.z << ") expected=(" << target.x << "," << target.y << "," << target.z << ")\n";
      ++fails;
    }

    if (r1.sigmaKm > 1e-6) {
      std::cerr << "[test_geometry_bearing_triangulation] expected near-zero sigma for perfect bearings. sigmaKm=" << r1.sigmaKm
                << "\n";
      ++fails;
    }

    // Coasting without a measurement should keep the estimate valid.
    const auto r2 = sim::updateBearingTrack(tr, 1.0, false, o1, d1, 0.0, p);
    if (!r2.valid) {
      std::cerr << "[test_geometry_bearing_triangulation] track became invalid when coasting without measurements.\n";
      ++fails;
    } else if (!approxVec3(r2.posKm, target, 1e-6)) {
      std::cerr << "[test_geometry_bearing_triangulation] coasting changed position unexpectedly. got=(" << r2.posKm.x << "," << r2.posKm.y
                << "," << r2.posKm.z << ")\n";
      ++fails;
    }
  }

  // --- Repeating a single bearing (no baseline) should never become solvable ---
  {
    sim::BearingTrack3d tr{};
    sim::BearingTrackParams p{};
    p.observationHalfLifeSec = 0.0;
    p.solveRegularization = 1.0e-6; // even with regularization, determinant gate should block
    p.determinantEps = 1.0e-12;
    p.minEffectiveWeight = 0.0;
    p.sigmaMinKm = 0.0;

    const math::Vec3d target{100.0, 0.0, 0.0};
    const math::Vec3d o{0.0, 0.0, 0.0};
    const math::Vec3d d = (target - o).normalized();

    bool everValid = false;
    for (int i = 0; i < 6; ++i) {
      const auto r = sim::updateBearingTrack(tr, 1.0, true, o, d, 1.0, p);
      everValid = everValid || r.valid;
    }

    if (everValid) {
      std::cerr << "[test_geometry_bearing_triangulation] repeated single bearing unexpectedly produced a valid solution.\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_geometry_bearing_triangulation] PASS\n";
  }

  return fails;
}
