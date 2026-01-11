#include "stellar/sim/SystemConditions.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) { return std::fabs(a - b) <= eps; }

int test_system_conditions() {
  int fails = 0;

  const core::u64 seed = 123456789ULL;

  sim::StarSystem sys{};
  sys.stub.id = 42;
  sys.stub.name = "TestSystem";
  sys.stub.posLy = {0.0, 0.0, 0.0};

  // One station => stable controlling faction id.
  sim::Station st{};
  st.id = 1;
  st.name = "Test Station";
  st.factionId = 7;
  sys.stations.push_back(st);

  const double t = 100.0;

  sim::SystemSecurityDynamicsParams dynParams{};
  sim::SystemEventParams evParams{};
  evParams.eventChance = 0.0; // deterministic off

  // --- Baseline (no dynamics, no events) ---
  {
    const auto snap = sim::snapshotSystemConditions(seed, sys, t, nullptr, dynParams, evParams);
    if (snap.systemId != sys.stub.id) {
      std::cerr << "FAIL: SystemConditionsSnapshot systemId mismatch.\n";
      ++fails;
    }

    if (snap.hasDynamics) {
      std::cerr << "FAIL: Baseline snapshot should not have dynamics.\n";
      ++fails;
    }

    if (snap.event.active) {
      std::cerr << "FAIL: Baseline snapshot should not have an active event.\n";
      ++fails;
    }

    if (!approx(snap.base.security01, snap.afterDynamics.security01) ||
        !approx(snap.base.piracy01, snap.afterDynamics.piracy01) ||
        !approx(snap.base.traffic01, snap.afterDynamics.traffic01)) {
      std::cerr << "FAIL: afterDynamics should equal base when no deltaState is supplied.\n";
      ++fails;
    }

    if (!approx(snap.base.security01, snap.effective.security01) ||
        !approx(snap.base.piracy01, snap.effective.piracy01) ||
        !approx(snap.base.traffic01, snap.effective.traffic01)) {
      std::cerr << "FAIL: effective should equal base when no events and no dynamics.\n";
      ++fails;
    }
  }

  // --- With dynamics (but no events) ---
  sim::SystemSecurityDeltaState delta{};
  delta.systemId = sys.stub.id;
  delta.lastUpdateDay = 0.0;
  delta.securityDelta = +0.12;
  delta.piracyDelta = -0.08;
  delta.trafficDelta = +0.05;

  {
    const auto snap = sim::snapshotSystemConditions(seed, sys, t, &delta, dynParams, evParams);

    if (!snap.hasDynamics) {
      std::cerr << "FAIL: Snapshot should report hasDynamics when deltaState is provided.\n";
      ++fails;
    }

    const auto base = sim::systemSecurityProfile(seed, sys);
    const auto expectedAfter = sim::applySystemSecurityDelta(base, delta, t, dynParams);
    const auto expectedDyn = sim::decayedSystemSecurityDelta(delta, t, dynParams);

    if (!approx(snap.afterDynamics.security01, expectedAfter.security01) ||
        !approx(snap.afterDynamics.piracy01, expectedAfter.piracy01) ||
        !approx(snap.afterDynamics.traffic01, expectedAfter.traffic01)) {
      std::cerr << "FAIL: afterDynamics doesn't match applySystemSecurityDelta.\n";
      ++fails;
    }

    if (!approx(snap.dynamicsNow.securityDelta, expectedDyn.securityDelta) ||
        !approx(snap.dynamicsNow.piracyDelta, expectedDyn.piracyDelta) ||
        !approx(snap.dynamicsNow.trafficDelta, expectedDyn.trafficDelta) ||
        !approx(snap.dynamicsNow.lastUpdateDay, expectedDyn.lastUpdateDay)) {
      std::cerr << "FAIL: dynamicsNow doesn't match decayedSystemSecurityDelta.\n";
      ++fails;
    }

    // No events => effective == afterDynamics.
    if (!approx(snap.effective.security01, snap.afterDynamics.security01) ||
        !approx(snap.effective.piracy01, snap.afterDynamics.piracy01) ||
        !approx(snap.effective.traffic01, snap.afterDynamics.traffic01)) {
      std::cerr << "FAIL: effective should equal afterDynamics when events are disabled.\n";
      ++fails;
    }
  }

  // --- Event wiring check (compare against direct SystemEvents calls) ---
  {
    sim::SystemEventParams p = evParams;
    p.eventChance = 1.0; // always roll an event slot (may still pick "None")
    p.cycleDays = 5.0;

    const auto base = sim::systemSecurityProfile(seed, sys);
    const auto ev = sim::generateSystemEvent(seed, sys.stub.id, t, base, p);
    const auto expectedEff = sim::applySystemEventToProfile(base, ev);

    const auto snap = sim::snapshotSystemConditions(seed, sys, t, nullptr, dynParams, p);

    if (snap.event.active != ev.active || snap.event.kind != ev.kind ||
        !approx(snap.event.startDay, ev.startDay) || !approx(snap.event.endDay, ev.endDay) ||
        !approx(snap.event.severity01, ev.severity01) ||
        !approx(snap.event.securityDelta, ev.securityDelta) ||
        !approx(snap.event.piracyDelta, ev.piracyDelta) ||
        !approx(snap.event.trafficDelta, ev.trafficDelta)) {
      std::cerr << "FAIL: snapshotSystemConditions event doesn't match generateSystemEvent.\n";
      ++fails;
    }

    if (!approx(snap.effective.security01, expectedEff.security01) ||
        !approx(snap.effective.piracy01, expectedEff.piracy01) ||
        !approx(snap.effective.traffic01, expectedEff.traffic01)) {
      std::cerr << "FAIL: snapshotSystemConditions effective doesn't match applySystemEventToProfile.\n";
      ++fails;
    }

    sim::SystemEvent outEv{};
    const auto eff = sim::effectiveSystemSecurityProfile(seed, sys, t, nullptr, dynParams, p, &outEv);
    if (!approx(eff.security01, snap.effective.security01) ||
        !approx(eff.piracy01, snap.effective.piracy01) ||
        !approx(eff.traffic01, snap.effective.traffic01) ||
        outEv.kind != snap.event.kind ||
        outEv.active != snap.event.active) {
      std::cerr << "FAIL: effectiveSystemSecurityProfile doesn't match snapshotSystemConditions.\n";
      ++fails;
    }
  }

  return fails;
}
