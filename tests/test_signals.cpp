#include "stellar/sim/Signals.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/Traffic.h"
#include "stellar/sim/TrafficLedger.h"
#include "stellar/sim/Universe.h"

#include <cmath>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

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

  // Pick a system that actually has stations (required by signal generation).
  const StarSystem* sysPtr = nullptr;
  const SystemStub* stubPtr = nullptr;
  for (const auto& s : stubs) {
    const auto& candidate = u.getSystem(s.id, &s);
    if (!candidate.stations.empty()) {
      sysPtr = &candidate;
      stubPtr = &s;
      break;
    }
  }
  if (!sysPtr || !stubPtr) {
    std::cerr << "[test_signals] no nearby system with stations\n";
    return 1;
  }

  const auto& stub = *stubPtr;
  const auto& sys = *sysPtr;

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

    // Optional traffic convoy payloads should be deterministic too.
    if (sa.hasTrafficConvoy != sb.hasTrafficConvoy) {
      std::cerr << "[test_signals] site mismatch (hasTrafficConvoy)\n";
      ++fails;
      break;
    }
    if (sa.hasTrafficConvoy) {
      const auto& ca = sa.trafficConvoy;
      const auto& cb = sb.trafficConvoy;
      if (ca.id != cb.id || ca.systemId != cb.systemId || ca.fromStation != cb.fromStation || ca.toStation != cb.toStation || ca.commodity != cb.commodity) {
        std::cerr << "[test_signals] traffic convoy mismatch (ids/route/commodity)\n";
        ++fails;
        break;
      }
      if (!nearly(ca.units, cb.units) || !nearly(ca.departDay, cb.departDay) || !nearly(ca.arriveDay, cb.arriveDay)) {
        std::cerr << "[test_signals] traffic convoy mismatch (numeric)\n";
        ++fails;
        break;
      }

      const auto& saState = sa.trafficState;
      const auto& sbState = sb.trafficState;
      if (saState.active != sbState.active) {
        std::cerr << "[test_signals] traffic state mismatch (active)\n";
        ++fails;
        break;
      }
      if (!nearly(saState.progress01, sbState.progress01) || !nearly(saState.distKm, sbState.distKm) || !nearly(saState.speedKmS, sbState.speedKmS)) {
        std::cerr << "[test_signals] traffic state mismatch (numeric)\n";
        ++fails;
        break;
      }
      if (!vecNearly(saState.velKmS, sbState.velKmS, 1e-6)) {
        std::cerr << "[test_signals] traffic state mismatch (velKmS)\n";
        ++fails;
        break;
      }
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

  // Traffic convoy sites should appear when enabled (requires >=2 stations).
  {
    const StarSystem* trafficSys = nullptr;
    for (const auto& s : stubs) {
      const auto& candidate = u.getSystem(s.id, &s);
      if (candidate.stations.size() >= 2) {
        trafficSys = &candidate;
        break;
      }
    }

    if (trafficSys) {
      SignalGenParams pt{};
      pt.resourceFieldCount = 0;
      pt.includeDailyDerelict = false;
      pt.includeDistress = false;
      pt.includeTrafficConvoys = true;

      // Increase the likelihood of at least one convoy by bumping traffic density.
      pt.trafficLaneParams.convoysPerDayBase = 6;
      pt.trafficLaneParams.convoysPerStation = 3;
      pt.trafficLaneParams.maxConvoysPerDay = 24;
      pt.trafficLaneParams.includeInactive = true;
      pt.trafficLaneParams.genWindowDays = 0;
      pt.trafficLaneParams.minDurationDays = 0.25;
      pt.trafficLaneParams.maxDurationDays = 0.75;

      std::vector<Mission> emptyMissions;
      std::vector<core::u64> emptyResolved;

      const auto all = generateSystemSignals(seed, *trafficSys, t, emptyMissions, emptyResolved, pt);
      const SignalSite* first = nullptr;
      for (const auto& s : all.sites) {
        if (s.kind == SignalKind::TrafficConvoy && s.hasTrafficConvoy) {
          first = &s;
          break;
        }
      }

      if (!first) {
        std::cerr << "[test_signals] expected traffic convoy sites when enabled\n";
        ++fails;
      } else {
        // Pick a time inside this convoy so it must be active.
        const double midT = 0.5 * (first->trafficConvoy.departDay + first->trafficConvoy.arriveDay);

        pt.trafficLaneParams.includeInactive = false;
        pt.trafficLaneParams.genWindowDays = 1;

        const auto mid = generateSystemSignals(seed, *trafficSys, midT, emptyMissions, emptyResolved, pt);
        bool found = false;
        for (const auto& s : mid.sites) {
          if (s.kind != SignalKind::TrafficConvoy) continue;
          if (s.id != first->id) continue;
          if (!s.hasTrafficConvoy) continue;
          found = true;
          if (!s.trafficState.active) {
            std::cerr << "[test_signals] expected chosen convoy to be active at midT\n";
            ++fails;
          }
          if (s.trafficState.progress01 < -1e-6 || s.trafficState.progress01 > 1.0 + 1e-6) {
            std::cerr << "[test_signals] unexpected convoy progress01\n";
            ++fails;
          }
          if (!nearly(s.expireDay, s.trafficConvoy.arriveDay, 1e-9)) {
            std::cerr << "[test_signals] convoy signal expireDay should match arriveDay\n";
            ++fails;
          }
          break;
        }
        if (!found) {
          std::cerr << "[test_signals] expected to find chosen convoy in midT plan\n";
          ++fails;
        }
      }
    }
  }

  // If the caller provides a TrafficLedger recorded from simulateNpcTradeTraffic(...),
  // signal generation should prefer *those* convoys (bridging the market nudger to
  // visible traffic sites).
  {
    const StarSystem* ledgerSys = nullptr;
    TrafficLedger ledger;
    std::unordered_map<SystemId, int> stamps;

    // Find a nearby multi-station system that produces at least one recorded shipment.
    for (const auto& s : stubs) {
      const auto& candidate = u.getSystem(s.id, &s);
      if (candidate.stations.size() < 2) continue;

      ledger.clear();
      stamps.clear();

      simulateNpcTradeTraffic(u, candidate, t, stamps, /*kMaxBackfillDays=*/14, &ledger);
      if (!ledger.shipments.empty()) {
        ledgerSys = &candidate;
        break;
      }
    }

    if (ledgerSys) {
      SignalGenParams pt{};
      pt.resourceFieldCount = 0;
      pt.includeDailyDerelict = false;
      pt.includeDistress = false;
      pt.includeTrafficConvoys = true;

      // We want to see the full day's manifest so we can pick a specific convoy.
      pt.trafficLaneParams.includeInactive = true;
      pt.trafficLaneParams.genWindowDays = 0;

      std::vector<Mission> emptyMissions;
      std::vector<core::u64> emptyResolved;

      const auto all = generateSystemSignals(seed, *ledgerSys, t, emptyMissions, emptyResolved, pt, &ledger);

      // The plan should include at least one convoy derived from the ledger.
      std::unordered_set<core::u64> shipmentIds;
      shipmentIds.reserve(ledger.shipments.size() * 2 + 1);
      for (const auto& sh : ledger.shipments) shipmentIds.insert(sh.id);

      const SignalSite* first = nullptr;
      for (const auto& s : all.sites) {
        if (s.kind != SignalKind::TrafficConvoy) continue;
        if (!s.hasTrafficConvoy) continue;
        if (shipmentIds.find(s.id) == shipmentIds.end()) {
          std::cerr << "[test_signals] traffic convoy id not found in ledger shipments\n";
          ++fails;
          break;
        }
        first = &s;
        break;
      }

      if (!first) {
        std::cerr << "[test_signals] expected ledger-backed traffic convoys\n";
        ++fails;
      } else {
        // Verify basic payload mapping back to the matching shipment.
        const TrafficShipment* match = nullptr;
        for (const auto& sh : ledger.shipments) {
          if (sh.id == first->id) { match = &sh; break; }
        }
        if (!match) {
          std::cerr << "[test_signals] could not find matching ledger shipment for convoy\n";
          ++fails;
        } else {
          if (first->trafficConvoy.commodity != match->commodity) {
            std::cerr << "[test_signals] ledger convoy commodity mismatch\n";
            ++fails;
          }
          if (!nearly(first->trafficConvoy.units, match->units, 1e-9)) {
            std::cerr << "[test_signals] ledger convoy units mismatch\n";
            ++fails;
          }
        }

        // Choose a time inside the convoy window so it must be active.
        const double midT = 0.5 * (first->trafficConvoy.departDay + first->trafficConvoy.arriveDay);

        pt.trafficLaneParams.includeInactive = false;
        pt.trafficLaneParams.genWindowDays = 1;

        const auto mid = generateSystemSignals(seed, *ledgerSys, midT, emptyMissions, emptyResolved, pt, &ledger);
        bool found = false;
        for (const auto& s : mid.sites) {
          if (s.kind != SignalKind::TrafficConvoy) continue;
          if (s.id != first->id) continue;
          if (!s.hasTrafficConvoy) continue;
          found = true;
          if (!s.trafficState.active) {
            std::cerr << "[test_signals] expected ledger convoy to be active at midT\n";
            ++fails;
          }
          break;
        }
        if (!found) {
          std::cerr << "[test_signals] expected to find ledger convoy in midT plan\n";
          ++fails;
        }
      }
    }
  }

  return fails;
}
