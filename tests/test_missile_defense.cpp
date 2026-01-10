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

  return fails;
}
