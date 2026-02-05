#include "test_harness.h"

#include "stellar/sim/CollisionWarning.h"

#include <cmath>
#include <vector>

using stellar::math::Vec3d;
using stellar::sim::CollisionWarningLevel;
using stellar::sim::CollisionWarningParams;
using stellar::sim::computeCollisionWarning;
using stellar::sim::ProximityFieldKm;
using stellar::sim::SphereObstacleKm;

static bool feq(double a, double b, double eps = 1e-6) {
  return std::fabs(a - b) <= eps;
}

int test_collision_warning() {
  int failures = 0;

  ProximityFieldKm field;
  {
    std::vector<SphereObstacleKm> obs;
    SphereObstacleKm s;
    s.centerKm = Vec3d{10.0, 0.0, 0.0};
    s.radiusKm = 1.0;
    s.hardness01 = 1.0;
    obs.push_back(s);
    field.build(std::move(obs));
  }

  CollisionWarningParams p;
  p.horizonSec = 30.0;
  p.padKm = 0.0;
  p.cautionTtiSec = 14.0;
  p.dangerTtiSec = 5.0;
  p.minSpeedKmS = 0.01;
  p.useStopDistance = true;
  p.stopMarginFactor = 0.10;

  // Case 1: moving toward the sphere => impact.
  {
    const Vec3d pos{0.0, 0.0, 0.0};
    const Vec3d vel{1.0, 0.0, 0.0}; // km/s
    const double maxDecel = 0.2;    // km/s^2

    const auto r = computeCollisionWarning(field, pos, vel, maxDecel, p);
    CHECK(r.hasImpact);
    CHECK(r.obstacleId == 0);
    CHECK(r.impactDistKm > 8.9 && r.impactDistKm < 9.1);  // surface at x=9
    CHECK(r.ttiSec > 8.9 && r.ttiSec < 9.1);
    CHECK(r.level == CollisionWarningLevel::Caution);
    CHECK(r.hazard01 > 0.1);
    CHECK(r.canStopBeforeImpact);
  }

  // Case 2: too little decel => can't stop before impact => danger.
  {
    const Vec3d pos{0.0, 0.0, 0.0};
    const Vec3d vel{1.0, 0.0, 0.0};
    const double maxDecel = 0.05; // stop dist = 10 km

    const auto r = computeCollisionWarning(field, pos, vel, maxDecel, p);
    CHECK(r.hasImpact);
    CHECK(!r.canStopBeforeImpact);
    CHECK(r.level == CollisionWarningLevel::Danger);
    CHECK(r.hazard01 > 0.5);
  }

  // Case 3: moving away => no impact.
  {
    const Vec3d pos{0.0, 0.0, 0.0};
    const Vec3d vel{-1.0, 0.0, 0.0};
    const double maxDecel = 0.2;

    const auto r = computeCollisionWarning(field, pos, vel, maxDecel, p);
    CHECK(!r.hasImpact);
    CHECK(r.level == CollisionWarningLevel::None);
    CHECK(feq(r.hazard01, 0.0));
  }

  // Case 4: ignoreId should suppress the hit.
  {
    const Vec3d pos{0.0, 0.0, 0.0};
    const Vec3d vel{1.0, 0.0, 0.0};

    const auto r = computeCollisionWarning(field, pos, vel, /*maxDecel=*/0.2, p, /*ignoreId=*/0);
    CHECK(!r.hasImpact);
    CHECK(r.level == CollisionWarningLevel::None);
  }

  return failures;
}
