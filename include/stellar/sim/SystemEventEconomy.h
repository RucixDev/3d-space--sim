#pragma once

#include "stellar/econ/Economy.h"
#include "stellar/sim/SystemEvents.h"

#include <array>
#include <string>

namespace stellar::sim {

// Tuning for how system events influence station economies.
//
// This layer is intentionally implemented as an additive per-day inventory flow
// (units/day) so that existing market quotes, trading, scanners, and tests can
// keep using the station's base StationEconomyModel without needing to thread an
// "effective model" through the whole codebase.
struct SystemEventEconomyParams {
  // Baseline scale: how much inventory can be pushed per day at severity=1 and weight=1,
  // as a fraction of desired stock.
  double flowScaleFracOfDesiredPerDay{0.10};

  // Clamp for absolute extra net flow: prevents extreme values if desired stock is huge.
  double maxAbsFlowFracOfDesiredPerDay{0.25};

  // One-shot inventory shocks applied at event start (cycle boundary).
  // Pirate raids + unrest remove a fraction of *current* inventory.
  double pirateRaidTheftFrac{0.22};
  double civilUnrestLossFrac{0.18};

  // Crackdowns / breakthroughs can inject inventory (fraction of desired stock).
  double securityCrackdownDumpFrac{0.22};
  double researchBreakthroughSurplusFrac{0.20};

  // Safety clamp.
  double maxSeverity01{1.0};
};

// Returns an additive net flow per commodity (units/day) that should be added to
// (productionPerDay - consumptionPerDay) when advancing the economy.
//
// Positive values increase inventory (tend to lower prices).
// Negative values decrease inventory (tend to raise prices).
std::array<double, econ::kCommodityCount> systemEventExtraNetPerDay(
    const econ::StationEconomyModel& model,
    const SystemEvent& ev,
    const SystemEventEconomyParams& params = {});

// Applies a one-shot inventory shock at the *start* of an event.
// Intended to be called exactly once when a new event begins.
void applySystemEventInventoryShock(econ::StationEconomyState& state,
                                   const econ::StationEconomyModel& model,
                                   const SystemEvent& ev,
                                   const SystemEventEconomyParams& params = {});

// Human-friendly, compact summary of expected market impact for bulletins / UI.
// Example:
//   "Market impact: Food↑ Water↑ Medicine↑ | Luxury↓ Stimulants↓"
std::string systemEventEconomySummary(const SystemEvent& ev,
                                     int maxUp = 3,
                                     int maxDown = 3);

} // namespace stellar::sim
