#include "stellar/sim/AgentAvoidance.h"

#include "test_harness.h"

#include "stellar/core/Hash.h"

#include <algorithm>
#include <cmath>
#include <vector>

using namespace stellar;

namespace {

static math::Vec3d stablePerp(const math::Vec3d& dirUnit) {
  const math::Vec3d axis = (std::abs(dirUnit.y) < 0.90) ? math::Vec3d{0, 1, 0} : math::Vec3d{1, 0, 0};
  auto v = math::cross(axis, dirUnit);
  auto n = v.normalized();
  if (n.lengthSq() > 1e-12) return n;
  v = math::cross(math::Vec3d{0, 0, 1}, dirUnit);
  return v.normalized();
}

} // namespace

int test_agent_avoidance() {
  using stellar::math::Vec3d;

  int failures = 0;

  sim::AgentAvoidanceParams p{};
  p.enabled = true;
  p.selfRadiusKm = 0.05;
  p.padKm = 0.0;
  p.horizonSec = 20.0;
  p.nearMissExtraKm = 5.0;
  p.minRelSpeedKmS = 0.0;
  p.selfVelBlend01 = 0.5;
  p.strength = 3.0;
  p.maxSteerDeg = 80.0;

  const Vec3d selfPosKm{0, 0, 0};
  const Vec3d selfVelKmS{1, 0, 0};
  const Vec3d desiredDir{1, 0, 0};
  const double desiredSpeedKmS = 1.0;

  std::vector<sim::AgentSphere> neigh;
  {
    sim::AgentSphere n{};
    n.id = 42;
    n.posKm = {10, 0, 0};
    n.velKmS = {-1, 0, 0};
    n.radiusKm = 0.05;
    n.hardness01 = 1.0;
    neigh.push_back(n);
  }

  const auto ar = sim::steerAvoidAgents(selfPosKm, selfVelKmS, desiredDir, desiredSpeedKmS, neigh, 123, p);
  CHECK(ar.steering);
  CHECK(ar.safeDirUnit.lengthSq() > 0.9);
  CHECK(std::abs(ar.safeDirUnit.y) + std::abs(ar.safeDirUnit.z) > 1e-3);

  // Compare closest-approach distance using desired direction vs safe direction.
  auto minSepKm = [&](const Vec3d& vSelfKmS) -> double {
    const auto& n = neigh[0];
    const Vec3d dp = n.posKm - selfPosKm;
    const Vec3d dv = n.velKmS - vSelfKmS;

    double t = 0.0;
    const double dv2 = stellar::math::dot(dv, dv);
    if (dv2 > 1e-12) {
      t = -stellar::math::dot(dp, dv) / dv2;
      t = std::clamp(t, 0.0, p.horizonSec);
    }

    const Vec3d sep = dp + dv * t;
    return sep.length();
  };

  const double sepDesired = minSepKm(desiredDir.normalized() * desiredSpeedKmS);
  const double sepSafe = minSepKm(ar.safeDirUnit.normalized() * desiredSpeedKmS);
  CHECK(sepSafe > sepDesired + 1e-6);

  const double rInfl = p.selfRadiusKm + p.padKm + neigh[0].radiusKm;
  CHECK(sepSafe > rInfl);

  // Determinism: same seed => same output.
  const auto ar2 = sim::steerAvoidAgents(selfPosKm, selfVelKmS, desiredDir, desiredSpeedKmS, neigh, 123, p);
  CHECK(stellar::math::dot(ar.safeDirUnit.normalized(), ar2.safeDirUnit.normalized()) > 0.999999);

  // Tie-break sanity: in a perfectly symmetric head-on case (sep ~ 0 at closest approach),
  // the direction of the perpendicular escape is chosen deterministically from (seed, neighborId).
  {
    const auto perp = stablePerp(desiredDir.normalized());
    const auto h = core::hashCombine(123, neigh[0].id);
    const double s = math::dot(ar.safeDirUnit.normalized(), perp.normalized());
    CHECK(std::abs(s) > 1e-6);
    if (h & 1ull) CHECK(s < 0.0);
    else CHECK(s > 0.0);
  }

  // No neighbors => no steering.
  std::vector<sim::AgentSphere> empty;
  const auto ar4 = sim::steerAvoidAgents(selfPosKm, selfVelKmS, desiredDir, desiredSpeedKmS, empty, 0, p);
  CHECK(!ar4.steering);
  CHECK(stellar::math::dot(ar4.safeDirUnit.normalized(), desiredDir.normalized()) > 0.999999);

  return failures;
}
