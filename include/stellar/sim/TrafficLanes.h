#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/System.h"

#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Traffic Lanes / Convoys (prototype)
// -----------------------------------------------------------------------------
//
// The core library already has a headless station-economy + "ambient" NPC traffic
// model (sim/Traffic.*) that nudges inventories over time.
//
// This module adds a *visible* counterpart: deterministic, on-rails "convoys"
// that travel between stations inside a star system.
//
// Design goals:
//  - Deterministic schedule (seed + system + dayStamp) so tools/tests/UI agree.
//  - Headless: no renderer or app dependencies.
//  - Cheap to query: systems have few stations, so O(stations^2 * commodities) is fine.
//
// NOTE: This is intentionally a lightweight prototype. The game can use this
// output to spawn actual trader/escort contacts, provide interdiction targets,
// or render "traffic lane" UI hints.

struct TrafficLaneParams {
  // Approximate convoy count per day.
  int convoysPerDayBase{2};
  int convoysPerStation{2};
  int maxConvoysPerDay{16};

  // Randomized trip duration targets (days). Duration is clamped by the speed
  // limits below.
  double minDurationDays{0.04}; // ~1.0h
  double maxDurationDays{0.35}; // ~8.4h

  // Speed clamps (km/s). For reference, the prototype supercruise max is ~18,000.
  double speedMinKmS{250.0};
  double speedMaxKmS{18000.0};

  // Lane arc styling: add a smooth sideways offset that is zero at endpoints.
  // This makes the route visually distinct from a straight line.
  double arcMinKm{8000.0};
  double arcMaxKm{120000.0};
  double arcMaxFracOfDistance{0.25};

  // When generating convoys around a query time, include schedules for the
  // surrounding +/- window days (helps capture trips that straddle midnight).
  int genWindowDays{1};

  // If true, generateTrafficConvoys() returns all scheduled convoys in the
  // window (with state evaluated at timeDays). If false, only active convoys
  // are returned.
  bool includeInactive{false};
};

// A deterministic shipment traveling between two stations.
struct TrafficConvoy {
  core::u64 id{0};
  SystemId systemId{0};

  StationId fromStation{0};
  StationId toStation{0};

  // Faction for law/security flavor (defaults to origin station faction).
  core::u32 factionId{0};

  econ::CommodityId commodity{econ::CommodityId::Food};
  double units{0.0};

  // Absolute schedule (timeDays).
  double departDay{0.0};
  double arriveDay{0.0};
};

// Evaluated state of a convoy at a given time.
struct TrafficConvoyState {
  bool active{false};
  double progress01{0.0}; // clamped 0..1

  double distKm{0.0};
  double speedKmS{0.0};

  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};
  math::Vec3d dir{0, 0, 1};
};

struct TrafficConvoyView {
  TrafficConvoy convoy{};
  TrafficConvoyState state{};
};

// Generate a deterministic convoy schedule for a specific day stamp.
//
// NOTE: The returned convoys may arrive after the dayStamp depending on station
// distance and speed clamps.
std::vector<TrafficConvoy> generateTrafficConvoysForDay(core::u64 universeSeed,
                                                        const StarSystem& system,
                                                        int dayStamp,
                                                        const TrafficLaneParams& params = {});

// Evaluate a convoy's position/velocity at a time. The returned state is safe
// to compute even when timeDays is outside [departDay, arriveDay] (progress is
// clamped and active=false).
TrafficConvoyState evaluateTrafficConvoy(const TrafficConvoy& convoy,
                                         const StarSystem& system,
                                         double timeDays,
                                         const TrafficLaneParams& params = {});

// Convenience: generate convoys around timeDays (dayStamp +/- genWindowDays)
// and return each convoy with evaluated state.
std::vector<TrafficConvoyView> generateTrafficConvoys(core::u64 universeSeed,
                                                      const StarSystem& system,
                                                      double timeDays,
                                                      const TrafficLaneParams& params = {});

} // namespace stellar::sim
