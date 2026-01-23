#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/sim/ShipLoadout.h"

namespace stellar::sim {

// -----------------------------------------------------------------------------
// ShipScan — deterministic "active scan" report generation
// -----------------------------------------------------------------------------
//
// The renderer/game layer already implements a time-based "scan" interaction.
// Historically, finishing a scan was only meaningful for a small handful of
// missions (distress / bounty scan).
//
// This module turns that scan action into a reusable, headless mechanic by
// generating a bounded, deterministic "scan report" (threat estimate, cargo
// estimate, and EW hints) from a snapshot of target state.
//
// Notes:
//  - The scan report is deliberately *imperfect*: it introduces bounded noise
//    that shrinks as scan quality rises.
//  - The noise is deterministic per (universeSeed, targetId), so tests and
//    replays stay stable.

struct ShipScanInput {
  core::u64 targetId{0};

  // Loadout
  ShipHullClass hullClass{ShipHullClass::Scout};
  int thrusterMk{1};
  int shieldMk{1};
  int distributorMk{1};
  WeaponType weapon{WeaponType::BeamLaser};

  // Current condition
  double hullFrac{1.0};   // 0..1
  double shieldFrac{1.0}; // 0..1

  // AI skill in [0,1] (used as a mild threat multiplier).
  double aiSkill01{0.5};

  // Electronic warfare — arbitrary "jammer power" used by ElectronicWarfare.
  double jammerPower{0.0};

  // Cargo snapshot.
  // If cargoUnits>0, cargoValue is derived from commodity base price.
  // Otherwise cargoValueCr is used as a fallback abstraction.
  econ::CommodityId cargoCommodity{econ::CommodityId::Food};
  double cargoUnits{0.0};
  double cargoValueCr{0.0};

  // Sensor confidence observed during the scan.
  // In stellar_game this maps nicely to SensorTrackResult::strength01.
  double scanStrength01{1.0}; // 0..1
};

struct ShipScanReport {
  core::u64 targetId{0};

  // Overall scan quality / confidence [0,1].
  double quality01{0.0};

  // Threat rating [0,1] derived from loadout, skill, and current condition.
  double threat01{0.0};

  // Identification.
  bool hullKnown{false};
  ShipHullClass hullClass{ShipHullClass::Scout};

  bool weaponKnown{false};
  WeaponType weapon{WeaponType::BeamLaser};

  bool healthKnown{false};
  double hullFrac{1.0};
  double shieldFrac{1.0};

  // Cargo estimate.
  bool cargoDetected{false};         // any meaningful cargo signature
  bool cargoKnown{false};            // have a numeric value estimate
  bool cargoCommodityKnown{false};   // commodity classification is reliable
  bool cargoUnitsKnown{false};       // unit estimate is reliable

  econ::CommodityId cargoCommodity{econ::CommodityId::Food};
  double cargoUnitsEst{0.0};

  // Credits estimate with a bounded error band.
  double cargoValueEstCr{0.0};
  double cargoValueMinCr{0.0};
  double cargoValueMaxCr{0.0};

  // Convenience: 0..1 confidence for cargo estimate (0 => unknown).
  double cargoConfidence01{0.0};

  // Electronic warfare hints.
  bool jammerSuspected{false};
  bool jammerDetected{false};
  double jammerStrength01{0.0}; // 0..1 (only meaningful if detected)
};

// Convert an instantaneous scan strength and jammer power into a scan quality.
//
// This is intentionally conservative: at the edge of the radar you can still
// complete a scan, but the report will be noisy and may hide details.
//
// Returns value in [0,1].
double shipScanQuality01(double scanStrength01, double jammerPower);

// Approximate threat rating in [0,1] from a target snapshot.
double shipThreatRating01(const ShipScanInput& in);

// Build a scan report.
//
// universeSeed is used only to make bounded noise deterministic.
ShipScanReport computeShipScanReport(core::u64 universeSeed, const ShipScanInput& in);

} // namespace stellar::sim
