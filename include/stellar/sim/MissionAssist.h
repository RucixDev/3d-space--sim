#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/econ/Economy.h"
#include "stellar/sim/SaveGame.h"

#include <cstddef>
#include <vector>

namespace stellar::sim {

class Universe;
struct StarSystem;
struct Station;

// -----------------------------------------------------------------------------
// Mission logistics helpers
// -----------------------------------------------------------------------------
//
// The prototype has a rich mission board, but once a contract is accepted the
// player can still get stuck in "logistics limbo":
//   - a buy-to-deliver mission where the commodity isn't available at the
//     current station
//   - losing mission cargo to piracy / jettison / mistakes
//
// This module provides a deterministic, headless "cargo sourcing" scan that
// searches nearby systems for stations that sell the required commodity.
// It is intentionally lightweight (no pathfinding); the game layer can decide
// whether to plot a route or filter candidates by jump range.

struct MissionCargoSourceParams {
  // Nearby-system query settings.
  double searchRadiusLy{220.0};
  std::size_t maxSystems{160};

  // Candidate trimming.
  std::size_t maxResults{8};

  // Market quote tuning.
  double bidAskSpread{0.10};

  // If false, the current system is excluded from the scan.
  bool includeCurrentSystem{true};

  // If true, only keep stations with inventory >= missingUnits.
  // If false, stations with partial inventory can still be suggested.
  bool requireEnoughInventory{false};
};

struct MissionCargoSourceCandidate {
  SystemId systemId{0};
  StationId stationId{0};
  econ::StationType stationType{econ::StationType::Outpost};

  // Straight-line distance from the origin system to the candidate system.
  double distanceLy{0.0};

  // Market estimate at `timeDays`.
  double inventoryUnits{0.0};
  double askCr{0.0};
  double feeRate{0.0};
  double askEffCr{0.0}; // askCr * (1 + feeRate)
};

struct MissionCargoSourcePlan {
  // True if the scan was executed (it can still return zero candidates).
  bool ok{false};

  econ::CommodityId commodity{econ::CommodityId::Food};
  double missingUnits{0.0};
  double missingMassKg{0.0};

  // Sorted "best first".
  std::vector<MissionCargoSourceCandidate> candidates{};
};

// Find stations that can supply (buy) the given commodity.
//
// Deterministic given:
//   - universe seed
//   - timeDays
//   - origin system
//   - params
MissionCargoSourcePlan planMissionCargoSourcing(Universe& universe,
                                                const StarSystem& originSystem,
                                                double timeDays,
                                                econ::CommodityId commodity,
                                                double missingUnits,
                                                const MissionCargoSourceParams& params = {});

// Convenience helper: compute missingUnits from a mission + current cargo, and
// run planMissionCargoSourcing() when the mission implies a commodity objective.
//
// Returns ok=false when the mission doesn't require cargo sourcing.
MissionCargoSourcePlan planMissionCargoSourcingForMission(Universe& universe,
                                                          const StarSystem& originSystem,
                                                          double timeDays,
                                                          const Mission& mission,
                                                          const std::array<double, econ::kCommodityCount>& cargo,
                                                          const MissionCargoSourceParams& params = {});

} // namespace stellar::sim
