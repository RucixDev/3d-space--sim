#include "stellar/sim/Signals.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/Universe.h"

#include <cmath>
#include <iostream>

static bool nearly(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

static bool vecNearly(const stellar::math::Vec3d& a, const stellar::math::Vec3d& b, double eps = 1e-6) {
  return nearly(a.x, b.x, eps) && nearly(a.y, b.y, eps) && nearly(a.z, b.z, eps);
}

int test_signals() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;

  const core::u64 seed = 424242ull;
  Universe u(seed);

  const auto stubs = u.queryNearby({0,0,0}, 35.0, 8);
  if (stubs.empty()) {
    std::cerr << "[test_signals] universe returned no stubs\n";
    return 1;
  }

  const auto& stub = stubs.front();
  const auto& sys = u.getSystem(stub.id, &stub);

  const double t = 1000.25;

  SignalGenParams p{};
  p.resourceFieldCount = 1;
  p.includeDailyDerelict = true;
  p.includeDistress = true;
  p.distressPerDay = 2;
  p.distressTtlDays = 1.0;

  std::vector<Mission> missions;
  std::vector<core::u64> resolved;

  const auto a = generateSystemSignals(seed, sys, t, missions, resolved, p);
  const auto b = generateSystemSignals(seed, sys, t, missions, resolved, p);

  // Basic determinism check.
  if (a.sites.size() != b.sites.size()) {
    std::cerr << "[test_signals] sites size mismatch\n";
    ++fails;
  }

  const std::size_t n = std::min(a.sites.size(), b.sites.size());
  for (std::size_t i = 0; i < n; ++i) {
    const auto& sa = a.sites[i];
    const auto& sb = b.sites[i];

    if (sa.id != sb.id || sa.kind != sb.kind || sa.resolved != sb.resolved) {
      std::cerr << "[test_signals] site mismatch (id/kind/resolved)\n";
      ++fails;
      break;
    }
    if (!nearly(sa.expireDay, sb.expireDay, 1e-9)) {
      std::cerr << "[test_signals] site mismatch (expireDay)\n";
      ++fails;
      break;
    }
    if (!vecNearly(sa.posKm, sb.posKm, 1e-6)) {
      std::cerr << "[test_signals] site mismatch (posKm)\n";
      ++fails;
      break;
    }
    if (sa.kind == SignalKind::Distress) {
      if (!sa.hasDistressPlan || !sb.hasDistressPlan) {
        std::cerr << "[test_signals] distress missing plan\n";
        ++fails;
        break;
      }
      if (sa.distress.scenario != sb.distress.scenario || sa.distress.needCommodity != sb.distress.needCommodity) {
        std::cerr << "[test_signals] distress plan mismatch (scenario/commodity)\n";
        ++fails;
        break;
      }
      if (!nearly(sa.distress.needUnits, sb.distress.needUnits) || !nearly(sa.distress.rewardCr, sb.distress.rewardCr)) {
        std::cerr << "[test_signals] distress plan mismatch (numeric)\n";
        ++fails;
        break;
      }
    }
  }

  // Resource field plan should be deterministic too.
  if (a.resourceFields.fields.size() != b.resourceFields.fields.size() || a.resourceFields.asteroids.size() != b.resourceFields.asteroids.size()) {
    std::cerr << "[test_signals] resource field plan size mismatch\n";
    ++fails;
  }
  if (!a.resourceFields.fields.empty() && !b.resourceFields.fields.empty()) {
    if (a.resourceFields.fields.front().id != b.resourceFields.fields.front().id) {
      std::cerr << "[test_signals] resource field id mismatch\n";
      ++fails;
    }
  }

  // Resolved signals should flip the resolved flag (daily derelict).
  {
    core::u64 derelictId = 0;
    for (const auto& s : a.sites) {
      if (s.kind == SignalKind::Derelict) { derelictId = s.id; break; }
    }
    if (derelictId == 0) {
      std::cerr << "[test_signals] expected a derelict signal\n";
      ++fails;
    } else {
      resolved.push_back(derelictId);
      const auto c = generateSystemSignals(seed, sys, t, missions, resolved, p);
      bool found = false;
      for (const auto& s : c.sites) {
        if (s.id == derelictId && s.kind == SignalKind::Derelict) {
          found = true;
          if (!s.resolved) {
            std::cerr << "[test_signals] resolved derelict not marked resolved\n";
            ++fails;
          }
          break;
        }
      }
      if (!found) {
        std::cerr << "[test_signals] resolved derelict disappeared (should remain, but marked)\n";
        ++fails;
      }
    }
  }

  // Mission salvage sites should appear when we have an active salvage mission.
  {
    Mission m{};
    m.id = 999ull;
    m.type = MissionType::Salvage;
    m.toSystem = sys.stub.id;
    m.toStation = sys.stations.empty() ? 0 : sys.stations.front().id;
    m.deadlineDay = t + 3.0;
    m.targetNpcId = 0x8000000000000000ull | 12345ull;
    m.scanned = false;

    missions.push_back(m);

    const auto d = generateSystemSignals(seed, sys, t, missions, resolved, p);
    bool found = false;
    for (const auto& s : d.sites) {
      if (s.kind == SignalKind::MissionSalvage && s.id == m.targetNpcId) {
        found = true;
        if (s.missionId != m.id) {
          std::cerr << "[test_signals] mission salvage site did not preserve missionId\n";
          ++fails;
        }
        break;
      }
    }
    if (!found) {
      std::cerr << "[test_signals] missing mission salvage site\n";
      ++fails;
    }
  }

  return fails;
}
