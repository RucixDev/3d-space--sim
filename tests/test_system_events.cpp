#include "stellar/sim/SystemEvents.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) { return std::fabs(a - b) <= eps; }

int test_system_events() {
  int fails = 0;

  const core::u64 seed = 123456789ULL;
  const sim::SystemId sysId = 42;

  sim::SystemSecurityProfile base{};
  base.security01 = 0.35;
  base.piracy01 = 0.70;
  base.traffic01 = 0.45;
  base.contest01 = 0.20;
  base.traits.authority = 0.60;
  base.traits.stability = 0.40;
  base.traits.wealth = 0.35;
  base.traits.tech = 0.55;
  base.traits.corruption = 0.65;

  sim::SystemEventParams p{};
  p.cycleDays = 5.0;
  p.eventChance = 1.0; // always roll an event (kind may still be None by weight).

  const auto evA = sim::generateSystemEvent(seed, sysId, 10.2, base, p);
  const auto evB = sim::generateSystemEvent(seed, sysId, 10.8, base, p);
  const auto evA2 = sim::generateSystemEvent(seed, sysId, 10.2, base, p);
  const auto evNext = sim::generateSystemEvent(seed, sysId, 15.1, base, p);

  // Deterministic for identical inputs.
  if (evA.kind != evA2.kind || !approx(evA.severity01, evA2.severity01) ||
      !approx(evA.securityDelta, evA2.securityDelta) || !approx(evA.piracyDelta, evA2.piracyDelta) ||
      !approx(evA.trafficDelta, evA2.trafficDelta)) {
    std::cerr << "FAIL: SystemEvents not deterministic for identical inputs.\n";
    ++fails;
  }

  // Stable within a cycle.
  if (!approx(evA.startDay, evB.startDay) || !approx(evA.endDay, evB.endDay) || evA.kind != evB.kind ||
      !approx(evA.severity01, evB.severity01)) {
    std::cerr << "FAIL: SystemEvents not stable within a cycle.\n";
    ++fails;
  }

  // Boundary changes cycle interval bookkeeping.
  if (approx(evA.startDay, evNext.startDay)) {
    std::cerr << "FAIL: SystemEvents cycle start did not advance across boundary.\n";
    ++fails;
  }

  // Deltas should always be clamped to configured maxAbsDelta when active.
  if (evA.active) {
    if (std::fabs(evA.securityDelta) > p.maxAbsDelta + 1e-12 || std::fabs(evA.piracyDelta) > p.maxAbsDelta + 1e-12 ||
        std::fabs(evA.trafficDelta) > p.maxAbsDelta + 1e-12) {
      std::cerr << "FAIL: SystemEvents deltas exceed maxAbsDelta.\n";
      ++fails;
    }
  } else {
    // If inactive, deltas should be zeros (or extremely close).
    if (!approx(evA.securityDelta, 0.0) || !approx(evA.piracyDelta, 0.0) || !approx(evA.trafficDelta, 0.0)) {
      std::cerr << "FAIL: SystemEvents inactive event should have zero deltas.\n";
      ++fails;
    }
  }

  // Applying an event should clamp profile channels to [0,1].
  sim::SystemEvent manual{};
  manual.active = true;
  manual.kind = sim::SystemEventKind::PirateRaid;
  manual.securityDelta = -9.0;
  manual.piracyDelta = +9.0;
  manual.trafficDelta = -9.0;

  const auto applied = sim::applySystemEventToProfile(base, manual);
  if (applied.security01 < 0.0 || applied.security01 > 1.0 || applied.piracy01 < 0.0 || applied.piracy01 > 1.0 ||
      applied.traffic01 < 0.0 || applied.traffic01 > 1.0) {
    std::cerr << "FAIL: applySystemEventToProfile did not clamp to [0,1].\n";
    ++fails;
  }

  return fails;
}
