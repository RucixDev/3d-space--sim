#pragma once

#include "stellar/sim/Faction.h"
#include "stellar/sim/System.h"

#include <vector>

namespace stellar::proc {

sim::StarSystem generateSystem(const sim::SystemStub& stub, const std::vector<sim::Faction>& factions);

} // namespace stellar::proc
