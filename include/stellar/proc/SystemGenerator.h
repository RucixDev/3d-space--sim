#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/Faction.h"
#include "stellar/sim/System.h"

#include <vector>

namespace stellar::proc {

// New (Round 23): explicit universe seed so the generated system can consume
// galaxy-coherent TradeProfiles.
sim::StarSystem generateSystem(core::u64 universeSeed,
                               const sim::SystemStub& stub,
                               const std::vector<sim::Faction>& factions);

// Legacy overload kept for external callers; uses `stub.seed` as the universe
// seed for trade profiling.
sim::StarSystem generateSystem(const sim::SystemStub& stub, const std::vector<sim::Faction>& factions);

} // namespace stellar::proc
