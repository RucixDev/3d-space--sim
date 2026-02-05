#include "stellar/sim/KinematicTrack.h"

#include "test_harness.h"

#include <cmath>

using namespace stellar;

namespace {

static bool approx(double a, double b, double eps) {
  return std::abs(a - b) <= eps;
}

static bool approxVec(const math::Vec3d& a, const math::Vec3d& b, double eps) {
  return approx(a.x, b.x, eps) && approx(a.y, b.y, eps) && approx(a.z, b.z, eps);
}

} // namespace

int test_kinematic_track() {
  int failures = 0;

  sim::KinematicTrackParams p;
  // Fix gains so the test behaves deterministically across platforms.
  p.minAlpha = p.maxAlpha = 0.75;
  p.minBeta = p.maxBeta = 0.12;
  p.maxSpeedKmS = 20.0;

  p.sigmaInitKm = 10.0;
  p.sigmaMinKm = 0.1;
  p.sigmaMaxKm = 1000.0;
  p.sigmaGrowthKmPerSec = 1.0;
  p.sigmaShrinkFactor = 0.50;

  sim::KinematicTrack3d tr{};
  CHECK(!tr.initialized);

  // ---- Initialization ----
  {
    const math::Vec3d z0{0.0, 0.0, 0.0};
    const auto r = sim::updateKinematicTrack(tr, /*dtSec=*/0.0, /*hasMeasurement=*/true, z0, /*strength01=*/1.0, p);
    CHECK(tr.initialized);
    CHECK(approxVec(r.posKm, z0, 1e-12));
    CHECK(approx(r.ageSinceMeasSec, 0.0, 1e-12));
    CHECK(r.sigmaKm >= p.sigmaMinKm);
  }

  // ---- Converge on a constant-velocity target ----
  // True motion: x(t) = t km, v = 1 km/s.
  {
    const math::Vec3d vTrue{1.0, 0.0, 0.0};
    for (int t = 1; t <= 12; ++t) {
      const math::Vec3d z = vTrue * (double)t;
      sim::updateKinematicTrack(tr, /*dtSec=*/1.0, /*hasMeasurement=*/true, z, /*strength01=*/1.0, p);
    }

    // After enough measurements the velocity estimate should be close.
    CHECK(std::abs(tr.velKmS.x - 1.0) < 0.35);
    // And the position should be near the last measurement.
    CHECK(std::abs(tr.posKm.x - 12.0) < 0.75);
  }

  // ---- Dropout: coast using the estimated velocity and grow uncertainty ----
  {
    const double sigmaBefore = tr.sigmaKm;
    const double xBefore = tr.posKm.x;
    const double vx = tr.velKmS.x;

    for (int i = 0; i < 5; ++i) {
      sim::updateKinematicTrack(tr, /*dtSec=*/1.0, /*hasMeasurement=*/false, /*measPosKm=*/math::Vec3d{}, /*strength01=*/0.0, p);
    }

    CHECK(tr.ageSinceMeasSec > 0.0);
    CHECK(tr.sigmaKm > sigmaBefore);

    // Prediction should move forward roughly along +X.
    const double dx = tr.posKm.x - xBefore;
    CHECK(dx > 0.25 * std::abs(vx) * 5.0);
  }

  // ---- Reacquire: shrink uncertainty and reset age counter ----
  {
    const double sigmaBefore = tr.sigmaKm;
    const math::Vec3d z{17.0, 0.0, 0.0};
    sim::updateKinematicTrack(tr, /*dtSec=*/0.0, /*hasMeasurement=*/true, z, /*strength01=*/1.0, p);

    CHECK(approx(tr.ageSinceMeasSec, 0.0, 1e-12));
    CHECK(tr.sigmaKm < sigmaBefore);
  }

  return failures;
}
