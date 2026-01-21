#include "stellar/sim/MissileDefense.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_missile_defense() {
  int fails = 0;

  const sim::CombatTargetKind kind = sim::CombatTargetKind::Ship;
  const core::u64 targetId = 123;
  const math::Vec3d targetPos{0, 0, 10};
  const math::Vec3d targetVel{0, 0, 0};

  {
    std::vector<sim::Missile> missiles;
    missiles.push_back(sim::Missile{});

    auto& m = missiles.back();
    m.hasTarget = true;
    m.targetKind = kind;
    m.targetId = targetId;
    m.seeker = sim::MissileSeekerType::Heat;
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 5};
    m.ttlSimSec = 10.0;

    const auto t = sim::nearestInboundMissile(
      missiles.data(), missiles.size(), kind, targetId, targetPos, targetVel, sim::MissileThreatParams{});

    if (!t.inbound) {
      std::cerr << "FAIL: expected inbound missile\n";
      fails++;
    }
    if (!approx(t.distKm, 10.0, 1e-6)) {
      std::cerr << "FAIL: distKm expected 10, got " << t.distKm << "\n";
      fails++;
    }
    if (!approx(t.closingKmS, 5.0, 1e-6)) {
      std::cerr << "FAIL: closingKmS expected 5, got " << t.closingKmS << "\n";
      fails++;
    }
    if (!approx(t.ttiSec, 2.0, 1e-6)) {
      std::cerr << "FAIL: ttiSec expected 2, got " << t.ttiSec << "\n";
      fails++;
    }
  }

  {
    std::vector<sim::Missile> missiles;
    missiles.push_back(sim::Missile{});

    auto& m = missiles.back();
    m.hasTarget = true;
    m.targetKind = kind;
    m.targetId = targetId;
    m.seeker = sim::MissileSeekerType::Heat;
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, -5};
    m.ttlSimSec = 10.0;

    const auto t = sim::nearestInboundMissile(
      missiles.data(), missiles.size(), kind, targetId, targetPos, targetVel, sim::MissileThreatParams{});

    if (t.inbound) {
      std::cerr << "FAIL: expected NOT inbound missile (moving away)\n";
      fails++;
    }
  }

  {
    std::vector<sim::Missile> missiles;

    // Slow/long threat.
    {
      sim::Missile m{};
      m.hasTarget = true;
      m.targetKind = kind;
      m.targetId = targetId;
      m.seeker = sim::MissileSeekerType::Heat;
      m.posKm = {0, 0, 0};
      m.velKmS = {0, 0, 2};
      m.ttlSimSec = 10.0;
      missiles.push_back(m);
    }

    // Fast/near threat.
    {
      sim::Missile m{};
      m.hasTarget = true;
      m.targetKind = kind;
      m.targetId = targetId;
      m.seeker = sim::MissileSeekerType::Radar;
      m.posKm = {0, 0, 7};
      m.velKmS = {0, 0, 3};
      m.ttlSimSec = 10.0;
      missiles.push_back(m);
    }

    const auto t = sim::nearestInboundMissile(
      missiles.data(), missiles.size(), kind, targetId, targetPos, targetVel, sim::MissileThreatParams{});
    if (!t.inbound || t.seeker != sim::MissileSeekerType::Radar || t.missileIndex != 1) {
      std::cerr << "FAIL: expected nearest inbound radar missile (index 1)\n";
      fails++;
    }
  }


  // Evasion plan: offset pass -> should push away from the predicted closest-approach point.
  {
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 1};

    const math::Vec3d tgtPos{1, 0, 10};
    const math::Vec3d tgtVel{0, 0, 0};
    const auto plan = sim::planMissileEvasion(m, tgtPos, tgtVel, /*seed=*/42);
    if (!plan.valid) {
      std::cerr << "FAIL: expected valid evasion plan (offset case)\n";
      fails++;
    } else {
      if (!approx(plan.tClosestSec, 10.0, 1e-6)) {
        std::cerr << "FAIL: tClosestSec mismatch. got=" << plan.tClosestSec << " expected=10\n";
        fails++;
      }
      if (!approx(plan.missDistanceKm, 1.0, 1e-6)) {
        std::cerr << "FAIL: missDistanceKm mismatch. got=" << plan.missDistanceKm << " expected=1\n";
        fails++;
      }

      const math::Vec3d los = (m.posKm - tgtPos).normalized();
      const double lateralDot = std::fabs(math::dot(plan.dirWorld, los));
      if (lateralDot > 1e-6) {
        std::cerr << "FAIL: expected dirWorld lateral to LOS. |dot|=" << lateralDot << "\n";
        fails++;
      }

      // In this setup the missile passes on the target's -X side, so +X is a good escape direction.
      if (plan.dirWorld.x < 0.5) {
        std::cerr << "FAIL: expected dirWorld.x positive (move away). got=" << plan.dirWorld.x << "\n";
        fails++;
      }
    }
  }

  // Evasion plan: head-on -> seeded perpendicular direction.
  {
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.velKmS = {0, 0, 1};

    const math::Vec3d tgtPos{0, 0, 10};
    const math::Vec3d tgtVel{0, 0, 0};
    const auto p1 = sim::planMissileEvasion(m, tgtPos, tgtVel, /*seed=*/123);
    const auto p2 = sim::planMissileEvasion(m, tgtPos, tgtVel, /*seed=*/123);
    if (!p1.valid) {
      std::cerr << "FAIL: expected valid evasion plan (head-on)\n";
      fails++;
    } else {
      const double axisDot = std::fabs(math::dot(p1.dirWorld, math::Vec3d{0, 0, 1}));
      if (axisDot > 1e-6) {
        std::cerr << "FAIL: expected dirWorld ⟂ approach axis. |dot|=" << axisDot << "\n";
        fails++;
      }
    }

    // Deterministic for identical seeds.
    if (!approx(p1.dirWorld.x, p2.dirWorld.x, 1e-12) ||
        !approx(p1.dirWorld.y, p2.dirWorld.y, 1e-12) ||
        !approx(p1.dirWorld.z, p2.dirWorld.z, 1e-12)) {
      std::cerr << "FAIL: expected deterministic dirWorld for identical seed\n";
      fails++;
    }
  }
  return fails;
}
