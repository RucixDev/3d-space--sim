#include "stellar/sim/ThreatAvoidance.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-6) {
  return std::fabs(a - b) <= eps;
}

int test_threat_avoidance() {
  int fails = 0;

  // ---------------------------------------------------------------------------
  // Collision-only: moving toward an inflated sphere should produce a lateral
  // jink direction and recommend braking when hazard is high.
  // ---------------------------------------------------------------------------
  {
    sim::ProximityFieldKm field;
    std::vector<sim::SphereObstacleKm> obs;
    obs.push_back(sim::SphereObstacleKm{{0, 0, 0}, 1.0, 1.0});
    field.build(obs);

    sim::Ship ship;
    ship.setPositionKm({1.0, 0.0, 5.0});
    ship.setVelocityKmS({0.0, 0.0, -1.0});
    ship.setOrientation({1, 0, 0, 0});

    sim::ThreatAvoidanceParams p{};
    p.missileEnable = false;
    p.collisionEngageHazard01 = 0.2;      // engage aggressively (test)
    p.collisionBrakeEngageHazard01 = 0.5; // brake for danger
    p.maxThrust01 = 1.0;

    const auto r = sim::computeThreatAvoidance(
      ship,
      &field,
      /*maxDecelKmS2=*/0.5,
      /*missiles=*/nullptr,
      /*missileCount=*/0,
      sim::CombatTargetKind::Ship,
      /*targetId=*/0,
      /*seed=*/0x1234ull,
      p);

    if (!r.active || !r.collisionActive) {
      std::cerr << "[test_threat_avoidance] expected collision avoid to be active.\n";
      ++fails;
    } else {
      if (!r.input.brake) {
        std::cerr << "[test_threat_avoidance] expected brake recommendation under high hazard.\n";
        ++fails;
      }

      // Jink should be lateral to the current velocity.
      const math::Vec3d vDir = ship.velocityKmS().normalized();
      const double dv = std::fabs(math::dot(r.dirWorld, vDir));
      if (dv > 1e-4) {
        std::cerr << "[test_threat_avoidance] expected lateral jink (dot=" << dv << ")\n";
        ++fails;
      }

      // In this geometry, the lateral component should have a positive X.
      if (!(r.dirWorld.x > 0.05)) {
        std::cerr << "[test_threat_avoidance] expected jink to bias +X. got x=" << r.dirWorld.x << "\n";
        ++fails;
      }

      const double len = r.dirWorld.length();
      if (!approx(len, 1.0, 1e-6)) {
        std::cerr << "[test_threat_avoidance] expected normalized dirWorld. len=" << len << "\n";
        ++fails;
      }

      // With identity orientation, thrustLocal should match dirWorld scaled by thrust01.
      if (std::fabs(r.input.thrustLocal.x - r.dirWorld.x * r.thrust01) > 1e-6) {
        std::cerr << "[test_threat_avoidance] thrustLocal not aligned with dirWorld (x mismatch).\n";
        ++fails;
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Missile-only: nearest inbound missile should trigger a lateral evasion dir.
  // ---------------------------------------------------------------------------
  {
    sim::Ship ship;
    ship.setPositionKm({0.0, 0.0, 0.0});
    ship.setVelocityKmS({0.0, 0.0, 0.0});
    ship.setOrientation({1, 0, 0, 0});

    sim::Missile m{};
    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Player;
    m.targetId = 42;
    m.ttlSimSec = 60.0;
    m.posKm = {0.0, -120.0, 100.0};
    m.velKmS = {0.0, 0.0, -30.0};
    m.seeker = sim::MissileSeekerType::Radar;
    m.radarDopplerNotchKmS = 0.75;

    sim::ThreatAvoidanceParams p{};
    p.collisionEnable = false;
    p.missileEngageTtiSec = 25.0;
    p.missileEvasion.enforceLateralToLos = true;

    const auto r = sim::computeThreatAvoidance(
      ship,
      /*obstacles=*/nullptr,
      /*maxDecelKmS2=*/0.5,
      &m,
      /*missileCount=*/1,
      sim::CombatTargetKind::Player,
      /*targetId=*/42,
      /*seed=*/0x9876ull,
      p);

    if (!r.active || !r.missileActive || !r.missilePlan.valid) {
      std::cerr << "[test_threat_avoidance] expected missile evasion to be active/valid.\n";
      ++fails;
    } else {
      // Evasion dir should be lateral to LOS (LOS from target -> missile).
      const math::Vec3d los = math::safeNormalized(m.posKm - ship.positionKm(), math::Vec3d{0, 0, 1}, 1e-12);
      const double dl = std::fabs(math::dot(r.dirWorld, los));
      if (dl > 1e-4) {
        std::cerr << "[test_threat_avoidance] expected lateral-to-LOS evasion (dot=" << dl << ")\n";
        ++fails;
      }

      const double len = r.dirWorld.length();
      if (!approx(len, 1.0, 1e-6)) {
        std::cerr << "[test_threat_avoidance] expected normalized dirWorld under missile. len=" << len << "\n";
        ++fails;
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Blend: when both threats are present, the result should still be normalized
  // and active.
  // ---------------------------------------------------------------------------
  {
    sim::ProximityFieldKm field;
    std::vector<sim::SphereObstacleKm> obs;
    obs.push_back(sim::SphereObstacleKm{{0, 0, 0}, 1.0, 1.0});
    field.build(obs);

    sim::Ship ship;
    ship.setPositionKm({1.0, 0.0, 5.0});
    ship.setVelocityKmS({0.0, 0.0, -1.2});
    ship.setOrientation({1, 0, 0, 0});

    sim::Missile m{};
    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 7;
    m.ttlSimSec = 60.0;
    m.posKm = {-80.0, 0.0, 120.0};
    m.velKmS = {0.0, 0.0, -28.0};
    m.seeker = sim::MissileSeekerType::Heat;

    sim::ThreatAvoidanceParams p{};
    p.collisionEngageHazard01 = 0.2;
    p.missileEngageTtiSec = 30.0;
    p.maxThrust01 = 1.0;
    p.missileWeight = 0.75;
    p.collisionWeight = 1.0;

    const auto r = sim::computeThreatAvoidance(
      ship,
      &field,
      /*maxDecelKmS2=*/0.5,
      &m,
      /*missileCount=*/1,
      sim::CombatTargetKind::Ship,
      /*targetId=*/7,
      /*seed=*/0x5555ull,
      p);

    if (!r.active) {
      std::cerr << "[test_threat_avoidance] expected blended avoidance to be active.\n";
      ++fails;
    } else {
      const double len = r.dirWorld.length();
      if (!approx(len, 1.0, 1e-6)) {
        std::cerr << "[test_threat_avoidance] expected normalized blended dirWorld. len=" << len << "\n";
        ++fails;
      }
    }
  }

  return fails;
}
