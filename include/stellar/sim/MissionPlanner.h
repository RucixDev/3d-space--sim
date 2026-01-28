#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/SystemSecurityDynamics.h"

#include <limits>
#include <span>
#include <vector>

namespace stellar::sim {

class Universe;
struct StarSystem;

// Heuristic mission itinerary planner.
//
// This module produces a deterministic "best next stops" list for a set of
// missions given the player's current location and jump range.
//
// Design goals:
//  - lightweight: intended for real-time UI (Mission Control window)
//  - deterministic: stable tie-breaking for identical inputs
//  - headless: no ImGui dependency
//
// NOTE: This is *not* an optimal TSP solver. It uses greedy selection with
// route-cost estimates from NavRouteBatch.

struct MissionObjective {
  SystemId systemId{0};
  StationId stationId{0};

  // True when the objective is a "site" inside the system rather than a dock
  // action (e.g. mission salvage site / bounty target area).
  bool isSite{false};
};

// Return the next objective for a mission.
//
// - Multi-leg missions (MultiDelivery/Passenger) return the via stop when leg==0.
// - Salvage missions return the salvage site (system-level) until m.scanned is true.
// - Bounty missions return a system-level objective.
MissionObjective missionNextObjective(const Mission& m);

struct MissionItineraryParams {
  // Navigation constraints.
  double maxJumpLy{28.0};
  double costPerJump{1.0};
  double costPerLy{0.0};

  // Optional risk-aware navigation.
  //
  // If > 0, route planning will prefer safer paths (lower piracy/contest/
  // low-security) even if they are longer.
  //
  // cost model: legCost = costPerJump + costPerLy*d + navRiskWeightPerLy*avgRisk01*d
  // where avgRisk01 is the average endpoint travel-risk for the leg.
  double navRiskWeightPerLy{0.0};

  // Route batching search radius.
  // If <= 0, an automatic radius is chosen based on remaining objectives.
  double queryRadiusLy{0.0};
  std::size_t maxSystems{1800};

  // Stop limit (UI sanity).
  int maxStops{12};

  // Greedy scoring weights.
  double rewardWeight{1.0};
  double riskWeight{0.45};
  double urgencyWeight{0.35};

  // ETA model (rough travel-time estimate, in *simulation* seconds).
  //
  // This enables deadline-aware scoring to reason about whether an objective is
  // still feasible given a (very approximate) travel budget.
  //
  // Notes:
  //  - This is NOT a physics-accurate estimator.
  //  - The defaults intentionally include a per-stop overhead to account for
  //    in-system travel/menu/docking friction.
  bool etaAwareUrgency{true};
  double etaSecondsPerJump{45.0};
  double etaSecondsPerLy{0.0};
  double etaSecondsPerStop{600.0};
  double etaSecondsPerSite{420.0};

  // When true, objectives are grouped by system (stationId differences are
  // ignored). This can reduce plan length, but may hide that multiple docks are
  // required inside the same system.
  bool groupBySystem{false};
};

struct MissionItineraryStop {
  MissionObjective objective{};

  // Mission ids covered by this stop.
  std::vector<core::u64> missionIds;

  // Aggregates for UI.
  double rewardCr{0.0};
  double avgRisk01{0.0};
  double earliestDeadlineDay{0.0}; // 0 == none

  // Travel from previous stop.
  bool reachable{true};
  int hopsFromPrev{0};
  double distanceLyFromPrev{0.0};
  double costFromPrev{0.0};

  // ETA analytics (derived from MissionItineraryParams::eta* fields).
  //
  // etaDay is an *absolute* simulation timestamp.
  // etaSlackHours is relative to earliestDeadlineDay (positive = on time).
  double etaTravelDaysFromPrev{0.0};
  double etaDay{0.0};
  double etaSlackHours{std::numeric_limits<double>::infinity()};
};

struct MissionItineraryResult {
  bool ok{false};
  SystemId startSystemId{0};
  std::vector<MissionItineraryStop> stops;

  int unreachableStops{0};
  double totalCost{0.0};
  double totalRewardCr{0.0};
  int totalHops{0};
  double totalDistanceLy{0.0};
};

// Compute a deterministic itinerary for the given missions.
//
// `playerRepWithFaction` is an optional map (factionId -> reputation) used when
// scoring mission risk. If not provided, 0 rep is assumed.
MissionItineraryResult planMissionItinerary(Universe& universe,
                                           const StarSystem& currentSystem,
                                           double timeDays,
                                           std::span<const Mission> missions,
                                           std::span<const FactionReputation> playerRepWithFaction = {},
                                           std::span<const SystemSecurityDeltaState> securityDeltas = {},
                                           const MissionItineraryParams& params = {});

} // namespace stellar::sim