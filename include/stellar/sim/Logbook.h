#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Economy.h"
#include "stellar/sim/Celestial.h"

#include <cstddef>
#include <vector>

namespace stellar::sim {

// Persistent player logbook entry.
//
// This is intentionally compact and save-friendly:
//  - only numeric fields (no variable-length strings)
//  - stable IDs so entries can be referenced across save/load
//
// The UI can enrich display text by looking up system/station names from Universe.
enum class LogbookEntryKind : core::u8 {
  StarScan = 0,
  PlanetScan = 1,
  StationScan = 2,
  SignalScan = 3,
  AsteroidProspect = 4,
  SystemSurveyBonus = 5,
  Count
};

const char* logbookEntryKindName(LogbookEntryKind k);

// A single discovery / record in the player's logbook.
struct LogbookEntry {
  // Stable unique key for this entry (typically the scan key).
  core::u64 key{0};

  LogbookEntryKind kind{LogbookEntryKind::SignalScan};

  // Optional sub-kind (meaning depends on kind):
  //  - PlanetScan: (core::u8)PlanetType
  //  - StationScan: (core::u8)StationType
  //  - SignalScan: (core::u8)SignalKind
  core::u8 subKind{0};

  // Location context.
  SystemId systemId{0};
  StationId stationId{0};

  // Optional object id:
  //  - SignalScan: signal site id
  //  - AsteroidProspect: asteroid id
  //  - PlanetScan: planet index (stored as u64)
  core::u64 objectId{0};

  // Optional commodity context (used by asteroid prospecting).
  econ::CommodityId commodity{econ::CommodityId::Food};
  double units{0.0};

  // When discovered.
  double discoveredDay{0.0};

  // Base exploration value (credits). For non-exploration entries this can be 0.
  double valueCr{0.0};

  // True once this entry has been sold via an exploration data broker.
  bool sold{false};
};

// Sum of unsold base value across entries.
double logbookUnsoldValueCr(const std::vector<LogbookEntry>& entries);

// Mark all unsold entries as sold; returns the base value sold.
double logbookMarkAllSold(std::vector<LogbookEntry>& entries);

// Mark only entries whose key is in `keys` as sold; returns the base value sold.
double logbookMarkSoldByKeys(std::vector<LogbookEntry>& entries,
                             const std::vector<core::u64>& keys);

// Keep the logbook bounded. By default, preferentially evicts the oldest SOLD
// entries first, then oldest overall.
void pruneLogbook(std::vector<LogbookEntry>& entries,
                  std::size_t maxEntries = 2000);

} // namespace stellar::sim
