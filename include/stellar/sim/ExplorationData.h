#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Economy.h"
#include "stellar/sim/Celestial.h"
#include "stellar/sim/Logbook.h"
#include "stellar/sim/Signals.h"
#include "stellar/sim/System.h"

#include <cstddef>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Exploration Data
// -----------------------------------------------------------------------------
//
// The exploration loop produces:
//   - a stable 64-bit scan/discovery key (stored in SaveGame::scannedKeys)
//   - a LogbookEntry (stored in SaveGame::logbook)
//   - a base credit value (stored in LogbookEntry::valueCr)
//
// The rendering/gameplay layer is free to decide *when* to record a scan.
// This module only defines stable identifiers + value formulas.
//
// NOTE: These keys are player-local. They are not meant to be globally unique
// across all possible universes; they are stable within a save because they are
// derived from deterministic ids (SystemId/StationId/signalId/etc).

// Scan/discovery keys ----------------------------------------------------------
core::u64 scanKeyStar(SystemId systemId);
core::u64 scanKeyPlanet(SystemId systemId, std::size_t planetIndex);
core::u64 scanKeyStation(StationId stationId);
core::u64 scanKeySignal(core::u64 signalId);
core::u64 scanKeyAsteroid(core::u64 asteroidId);
core::u64 scanKeySystemComplete(SystemId systemId);

// Base value formulas (credits) ------------------------------------------------
// These are deliberately lightweight and "arcade plausible": they create a
// consistent reward ladder without requiring a full astrophysics model.

double scanValueStar(const Star& star);
double scanValuePlanet(const Planet& planet);
double scanValueStation(const Station& station);
double scanValueSignal(SignalKind kind);
double scanValueAsteroidProspect();
double scanValueSystemSurveyBonus(int planetCount, int stationCount);

// Optional broker-side multiplier helpers -------------------------------------
//
// The in-game broker UI may layer premiums/penalties on top of the base values.
// We expose a small helper to keep the multiplier logic deterministic and
// testable, while still letting the UI decide how it presents groupings.

struct ExplorationDataBrokerParams {
  // Distance-based premium.
  bool enableDistancePremium{true};
  double distanceScaleLy{250.0};     // distance that yields near-max premium
  double maxDistancePremium{0.30};   // +30% max

  // Faction alignment bonus/penalty.
  double sameFactionBonus{0.10};     // +10%
  double otherFactionPenalty{0.05};  // -5%

  // Station-type demand shaping.
  // If enabled, the station's economy type biases what discoveries it values.
  bool enableStationDemand{true};
  double demandStrength{0.35};       // 0 = off, 1 = full demand table

  // Final clamp.
  double minMultiplier{0.70};
  double maxMultiplier{1.65};
};

// A station-demand multiplier for a specific logbook entry kind.
//
// This is intentionally heuristic (not simulation-heavy): it makes selling data
// at a Research station feel different from selling at a Refinery, without
// requiring new economy state.
double explorationDataStationDemandMultiplier(econ::StationType stationType,
                                             LogbookEntryKind kind,
                                             core::u8 subKind);

// Full broker payout multiplier for a discovery originating in `scanSystem` and
// sold at `saleStation` in `saleSystem`.
double explorationDataBrokerMultiplier(const ExplorationDataBrokerParams& params,
                                       const SystemStub& saleSystem,
                                       const Station& saleStation,
                                       const SystemStub& scanSystem,
                                       LogbookEntryKind kind,
                                       core::u8 subKind);

} // namespace stellar::sim
