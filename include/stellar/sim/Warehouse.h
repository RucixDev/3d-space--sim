#pragma once

#include "stellar/econ/Commodity.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/System.h"

#include <array>
#include <vector>

namespace stellar::sim {

// Simple per-station cargo storage ("warehouse") system.
//
// The game can use this to let players deposit cargo into station storage and
// withdraw it later. Storage accrues fees over time based on:
//  - stored mass (kg)
//  - station type
//  - station tariff snapshot (Station::feeRate)
//  - player reputation (discounts / penalties)
//
// This module is intentionally deterministic and lightweight; it does not
// interact with the station economy model (inventories/prices).

struct WarehouseResult {
  bool ok{false};
  // Optional short machine-friendly reason on failure.
  // Common values:
  //  - "invalid"
  //  - "no_storage"
  //  - "empty"
  //  - "capacity"
  //  - "cargo_full"
  //  - "fees_due"
  const char* reason{nullptr};

  // Operation side-effects.
  double unitsMoved{0.0};
  double creditsPaid{0.0};
};

StationStorage* findStorage(std::vector<StationStorage>& storage, StationId stationId);
const StationStorage* findStorage(const std::vector<StationStorage>& storage, StationId stationId);

// Returns an existing station storage entry or creates a new one.
// When creating, captures a snapshot of the station fields required for fee computation.
StationStorage& getOrCreateStorage(std::vector<StationStorage>& storage,
                                   const Station& station,
                                   double timeDays);

// Total stored cargo mass (kg).
double storageMassKg(const StationStorage& entry);

// Maximum allowed storage mass (kg) for this station type.
// This is a gameplay-facing soft cap; it does not affect the station economy.
double storageCapacityKg(const StationStorage& entry);

// Fee model: credits per kg per day.
double storageFeeRateCrPerKgPerDay(const StationStorage& entry, double rep);

// Estimated daily storage fee for the current stored cargo.
double estimateStorageDailyFeeCr(const StationStorage& entry, double rep);

// Accrues storage fees for this entry up to `timeDays`.
// Updates entry.lastFeeDay.
void accrueStorageFees(StationStorage& entry, double timeDays, double rep);

// Pay some (or all) storage fees due at the station.
WarehouseResult payStorageFees(std::vector<StationStorage>& storage,
                               const Station& station,
                               double timeDays,
                               double rep,
                               double& credits,
                               double amountCr);

// Move cargo from ship to station storage.
WarehouseResult depositToStorage(std::vector<StationStorage>& storage,
                                 const Station& station,
                                 double timeDays,
                                 double rep,
                                 std::array<double, econ::kCommodityCount>& shipCargo,
                                 econ::CommodityId commodity,
                                 double units);

// Move cargo from station storage to ship.
// If fees are due and the player can afford them, the withdrawal will automatically
// pay the outstanding amount first (so withdrawals are "one click").
WarehouseResult withdrawFromStorage(std::vector<StationStorage>& storage,
                                   const Station& station,
                                   double timeDays,
                                   double rep,
                                   std::array<double, econ::kCommodityCount>& shipCargo,
                                   double& credits,
                                   double cargoCapacityKg,
                                   econ::CommodityId commodity,
                                   double units);

// Remove empty (no cargo, no fees) entries.
void pruneEmptyStorage(std::vector<StationStorage>& storage,
                       double epsUnits = 1e-6,
                       double epsFees = 1e-2);

} // namespace stellar::sim
