#pragma once

#include "stellar/sim/TrafficLedger.h"
#include "stellar/sim/TrafficLanes.h"

#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Traffic Convoy Layer (Ledger → Convoys)
// -----------------------------------------------------------------------------
//
// The sim has two related systems:
//  - `sim::simulateNpcTradeTraffic(...)` which mutates station inventories
//    deterministically using small station-to-station shipments.
//  - `sim::TrafficLanes` which can generate deterministic "on-rails" convoys.
//
// `TrafficLedger` records the *actual* shipments chosen by the ambient
// simulation.
//
// This module bridges those recorded shipments to the convoy primitives so
// tooling (and eventually gameplay) can render/track them as moving objects.

// Lossless conversion (copies ids + schedule + cargo metadata).
TrafficConvoy convoyFromShipment(const TrafficShipment& s);

// Convert recorded shipments to convoy views (with evaluated state).
//
// windowDays controls which day stamps are returned (timeDays +/- windowDays).
// If includeInactive is false, only shipments whose schedules are active at
// timeDays are returned.
std::vector<TrafficConvoyView> generateTrafficConvoysFromLedger(const TrafficLedger& ledger,
                                                                const StarSystem& system,
                                                                double timeDays,
                                                                int windowDays,
                                                                bool includeInactive = false,
                                                                const TrafficLaneParams& laneParams = {});

// Sample a convoy path as positions in km.
//
// Returns `segments + 1` points (including both endpoints) when segments >= 1.
std::vector<math::Vec3d> sampleTrafficConvoyPathKm(const TrafficConvoy& convoy,
                                                   const StarSystem& system,
                                                   int segments,
                                                   const TrafficLaneParams& laneParams = {});

} // namespace stellar::sim
