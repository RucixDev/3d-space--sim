#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/System.h"
#include "stellar/sim/Universe.h"

#include <unordered_map>
#include <vector>

namespace stellar::sim {

// Simulate low-cost "NPC trader" traffic for the given system.
//
// This nudges station inventories (and therefore prices) by moving commodity units
// between stations that naturally produce/consume complementary goods.
//
// Design goals:
//  - Deterministic across runs (seeded by universe.seed(), system id, and day stamp).
//  - Independent of update cadence (only steps integer day stamps).
//  - Bounded CPU cost (backfills at most kMaxBackfillDays on large time jumps).
//
// lastTrafficDayBySystem holds the last simulated day stamp per system.
void simulateNpcTradeTraffic(Universe& universe,
                             const StarSystem& system,
                             double timeDays,
                             std::unordered_map<SystemId, int>& lastTrafficDayBySystem,
                             int kMaxBackfillDays = 14);

// SaveGame-friendly overload.
//
// Convenience for callers that already store traffic stamps in the SaveGame format.
// The vector will be updated in-place (the matching system entry is added/updated).
void simulateNpcTradeTraffic(Universe& universe,
                             const StarSystem& system,
                             double timeDays,
                             std::vector<SystemTrafficStamp>& trafficStamps,
                             int kMaxBackfillDays = 14);

} // namespace stellar::sim
