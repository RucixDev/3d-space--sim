#include "test_harness.h"

#include "stellar/core/Assert.h"
#include "stellar/sim/Signals.h"
#include "stellar/sim/Universe.h"
#include "stellar/sim/SystemConditions.h"

#include <algorithm>
#include <optional>

#include <vector>

using namespace stellar;

// Validates that the optional SystemConditionsSnapshot hook in generateSystemSignals()
// can make signals "event-reactive" without impacting the baseline deterministic output
// when conditions are not provided.
int test_signal_event_reactivity() {
  const core::u64 seed = 0xC0FFEE1234ULL;

  sim::Universe u(seed);
  u.initGalaxy(sim::GalaxyParams{});

  // Pick a real generated system that has at least one station (signals need an anchor).
  std::optional<sim::SystemStub> stubOpt;
  {
    auto stubs = u.queryNearby(math::Vec3d{0.0, 0.0, 0.0}, /*radiusLy=*/50.0, /*maxResults=*/512);
    for (const auto& s : stubs) {
      const auto& sys = u.getSystem(s.id, &s);
      if (!sys.stations.empty()) {
        stubOpt = s;
        break;
      }
    }
  }
  CHECK(stubOpt.has_value());
  const auto stub = *stubOpt;
  const auto& sys = u.getSystem(stub.id, &stub);

  const double timeDays = 1234.25;
  const std::vector<sim::Mission> missions{};
  const std::vector<core::u64> resolved{};

  sim::SignalGenParams gen{};
  gen.resourceFieldCount = 0;
  gen.includeDailyDerelict = true;
  gen.includeDistress = true;
  gen.distressPerDay = 1;
  gen.distressTtlDays = 1.0;
  gen.includeTrafficConvoys = false;

  // Baseline (no conditions): should exactly match the requested distress density.
  const auto base = sim::generateSystemSignals(seed, sys, timeDays, missions, resolved, gen);

  int baseDistress = 0;
  int baseDerelicts = 0;
  for (const auto& s : base.sites) {
    if (s.kind == sim::SignalKind::Distress) ++baseDistress;
    if (s.kind == sim::SignalKind::Derelict) ++baseDerelicts;
  }
  CHECK_EQ(baseDistress, 1);
  CHECK_EQ(baseDerelicts, 1); // the daily derelict

  // Now fabricate a conditions snapshot with a severe pirate raid.
  sim::SystemConditionsSnapshot cond{};
  cond.systemId = sys.stub.id;
  cond.base = sim::systemSecurityProfile(seed, sys);
  cond.effective = cond.base;

  // Nudge the effective knobs into a plausible raid regime (more piracy, less security),
  // but keep them clamped for sanity.
  cond.effective.piracy01 = std::clamp(cond.effective.piracy01 + 0.35, 0.0, 1.0);
  cond.effective.security01 = std::clamp(cond.effective.security01 - 0.25, 0.0, 1.0);

  cond.event.active = true;
  cond.event.kind = sim::SystemEventKind::PirateRaid;
  cond.event.severity01 = 1.0;

  const auto boosted = sim::generateSystemSignals(seed, sys, timeDays, missions, resolved, gen, /*trafficLedger=*/nullptr, &cond);

  int boostedDistress = 0;
  int boostedDerelicts = 0;
  for (const auto& s : boosted.sites) {
    if (s.kind == sim::SignalKind::Distress) ++boostedDistress;
    if (s.kind == sim::SignalKind::Derelict) ++boostedDerelicts;
  }

  // PirateRaid: distressPerDay += 1 + round(sev*3) => +4 at sev=1 => 5 total.
  CHECK_EQ(boostedDistress, 5);

  // PirateRaid: extraDerelictsPerDay = ceil(sev*2) => 2 extra, plus the usual daily derelict.
  CHECK_EQ(boostedDerelicts, 3);

  return 0;
}
