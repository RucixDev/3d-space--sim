#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/sim/System.h"

#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Traffic Ledger (instrumentation)
// -----------------------------------------------------------------------------
//
// The core sim already supports a low-cost background market nudger:
//   sim::simulateNpcTradeTraffic(...)
//
// This header adds an *optional* instrumentation layer that can record the
// station-to-station shipments chosen by that model.
//
// The intent is to provide a bridge from "invisible inventory deltas" to
// gameplay/UI primitives:
//  - tooling can print daily shipment manifests
//  - the game can spawn visible convoy contacts/signals
//  - tests can verify determinism and continuity across seeds
//
// NOTE: Recording shipments is designed to be side-effect free.
// It must not consume RNG from the traffic simulation itself.

struct TrafficLedgerParams {
  // Keep only shipments whose dayStamp is within this many days of the most
  // recent prune() call.
  int keepDays{8};

  // Scheduling model (used only for metadata / visualization).
  // Duration targets are clamped by the speed limits below.
  double minDurationDays{0.04}; // ~1.0h
  double maxDurationDays{0.35}; // ~8.4h

  // Speed clamps (km/s).
  double speedMinKmS{250.0};
  double speedMaxKmS{18000.0};
};

// A single station-to-station shipment chosen by the ambient traffic model.
//
// Units are "commodity units" (same units used by station inventories).
// Times are in absolute "timeDays".
struct TrafficShipment {
  core::u64 id{0};
  SystemId systemId{0};
  int dayStamp{0};

  StationId fromStation{0};
  StationId toStation{0};
  core::u32 factionId{0};

  econ::CommodityId commodity{econ::CommodityId::Food};
  double units{0.0};

  // Optional schedule metadata for visualization.
  double departDay{0.0};
  double arriveDay{0.0};
  double distKm{0.0};
  double speedKmS{0.0};
};

// Create a deterministic shipment record with schedule metadata.
//
// IMPORTANT: This function is deterministic and does not depend on any mutable
// simulation state beyond the provided system/stations.
TrafficShipment makeNpcTradeShipment(core::u64 universeSeed,
                                    const StarSystem& system,
                                    int dayStamp,
                                    int runIndex,
                                    const Station& from,
                                    const Station& to,
                                    econ::CommodityId commodity,
                                    double units,
                                    const TrafficLedgerParams& params = {});

// Returns true when the shipment's schedule metadata appears to be valid.
//
// This is primarily used to support backward compatibility: older/corrupted save files
// may contain TrafficShipment records without depart/arrive metadata, which would
// otherwise prevent those shipments from being replayed as moving convoys.
bool shipmentScheduleLooksValid(const TrafficShipment& s);

// If the shipment's schedule metadata is missing/corrupt, rebuild it deterministically
// from (shipment id + system id + dayStamp) using the same model as makeNpcTradeShipment(...).
//
// Returns false if the referenced stations cannot be resolved inside `system`.
bool hydrateShipmentScheduleFromId(TrafficShipment& s,
                                  const StarSystem& system,
                                  const TrafficLedgerParams& params = {});


// Simple in-memory log of recent shipments.
struct TrafficLedger {
  TrafficLedgerParams params{};
  std::vector<TrafficShipment> shipments{};

  void clear();
  void prune(double timeDays);
  void record(const TrafficShipment& s);
  void record(TrafficShipment&& s);

  // Query shipments for a system in a +/- day window around timeDays.
  //
  // If includeInactive==false, only shipments whose schedule is "active"
  // at timeDays are returned.
  std::vector<TrafficShipment> query(SystemId systemId,
                                     double timeDays,
                                     int windowDays,
                                     bool includeInactive = false) const;
};

} // namespace stellar::sim
