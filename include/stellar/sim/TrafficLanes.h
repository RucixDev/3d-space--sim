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

  // ---------------------------------------------------------------------------
  // Lane corridor modeling
  // ---------------------------------------------------------------------------
  // If true, lane geometry is derived primarily from the unordered station pair
  // (systemId, min(from,to), max(from,to)). This makes multiple convoys between
  // the same endpoints share a stable "corridor".
  //
  // If false, lane geometry is derived from the convoy id (legacy behavior), so
  // each convoy can have a distinct arc plane/offset.
  bool bundleByStationPair{true};

  // If true, opposite-direction travel between the same station pair uses a
  // mirrored arc ("dual carriageway") by flipping the arc side sign.
  bool dualCarriageway{true};

  // Optional per-convoy arc magnitude jitter (fraction of base arc magnitude).
  // This keeps multiple same-route convoys from perfectly overlapping while
  // still forming a coherent corridor. 0 disables.
  double arcJitterFrac{0.0};
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

// Precomputed lane/corridor geometry for a specific convoy schedule.
//
// This is a purely geometric helper: it does not own simulation state.
// It exists to:
//  - make lane evaluation/sample cheaper (compute once, evaluate many)
//  - expose stable "lane keys" so UI/tools can group convoys into corridors
struct TrafficLaneGeometry {
  // Stable key used to group lanes into corridors.
  //
  // When TrafficLaneParams::bundleByStationPair==true, this is derived from the
  // unordered station pair (systemId, min(from,to), max(from,to)).
  // When false, this may be derived from the convoy id (legacy behavior).
  core::u64 laneKey{0};

  SystemId systemId{0};
  StationId fromStation{0};
  StationId toStation{0};

  // +1 or -1 for dual-carriageway lanes (sign applied to `side`).
  int directionSign{1};

  // Schedule snapshot.
  double departDay{0.0};
  double arriveDay{0.0};

  // Endpoints (km) evaluated at (departDay/arriveDay).
  math::Vec3d p0Km{0, 0, 0};
  math::Vec3d p1Km{0, 0, 0};

  // Orthonormal frame:
  //  - dir: along the chord from p0->p1
  //  - side: perpendicular direction used for lane arc offsets
  //  - up: completes the basis (useful for lane "ribbon" jitter)
  math::Vec3d dir{0, 0, 1};
  math::Vec3d side{1, 0, 0};
  math::Vec3d up{0, 1, 0};

  double distKm{0.0};
  double durSec{0.0};

  // Base arc amplitude (km). Applied as sin^2(pi*t) so endpoints are zero.
  double arcMagKm{0.0};

  // Jitter applied to arcMagKm for this specific convoy (km). Typically 0.
  double arcJitterKm{0.0};
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

// Compute the lane geometry used for a convoy.
//
// This is deterministic from (systemId, station ids, convoy schedule) and the
// lane parameters.
TrafficLaneGeometry computeTrafficLaneGeometry(const TrafficConvoy& convoy,
                                               const StarSystem& system,
                                               const TrafficLaneParams& params = {});

// Evaluate a convoy's position/velocity at a time. The returned state is safe
// to compute even when timeDays is outside [departDay, arriveDay] (progress is
// clamped and active=false).
TrafficConvoyState evaluateTrafficConvoy(const TrafficConvoy& convoy,
                                         const StarSystem& system,
                                         double timeDays,
                                         const TrafficLaneParams& params = {});

// Evaluate using a precomputed lane geometry (faster when sampling many times).
TrafficConvoyState evaluateTrafficConvoy(const TrafficLaneGeometry& lane,
                                         double timeDays,
                                         const TrafficLaneParams& params = {});

// Sample a lane path as positions in km.
//
// Returns `segments + 1` points (including both endpoints) when segments >= 1.
std::vector<math::Vec3d> sampleTrafficLanePathKm(const TrafficLaneGeometry& lane,
                                                 int segments,
                                                 const TrafficLaneParams& params = {});

// Convenience: generate convoys around timeDays (dayStamp +/- genWindowDays)
// and return each convoy with evaluated state.
std::vector<TrafficConvoyView> generateTrafficConvoys(core::u64 universeSeed,
                                                      const StarSystem& system,
                                                      double timeDays,
                                                      const TrafficLaneParams& params = {});

} // namespace stellar::sim
