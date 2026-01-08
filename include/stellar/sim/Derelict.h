#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/sim/Celestial.h"

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Derelict encounters
// -----------------------------------------------------------------------------
// A deterministic generator for lightweight "signal" derelict salvage encounters.
//
// This mirrors the design of sim::Distress:
//  - headless + deterministic
//  - small data model suitable for game/runtime + tooling/tests
//  - does NOT spawn ships or manage runtime state

enum class DerelictScenario : core::u8 {
  CivilianWreck = 0,
  SmugglerCache = 1,
  MilitaryHulk = 2,
  PirateTrap = 3,
};

inline const char* derelictScenarioName(DerelictScenario s) {
  switch (s) {
    case DerelictScenario::CivilianWreck: return "Civilian";
    case DerelictScenario::SmugglerCache: return "Smuggler";
    case DerelictScenario::MilitaryHulk: return "Military";
    case DerelictScenario::PirateTrap: return "Trap";
    default: return "Derelict";
  }
}

// Small set of size buckets used only for shaping loot and threat.
// 0 = small, 1 = medium, 2 = large.
inline const char* derelictWreckClassName(core::u8 wreckClass) {
  switch (wreckClass) {
    case 0: return "Small";
    case 1: return "Medium";
    case 2: return "Large";
    default: return "Unknown";
  }
}

struct DerelictPlan {
  DerelictScenario scenario{DerelictScenario::CivilianWreck};
  core::u8 wreckClass{0};

  // Salvage (cargo pods drifting near the wreck).
  bool hasSalvage{true};
  econ::CommodityId salvageCommodity{econ::CommodityId::Machinery};
  double salvageUnits{0.0};
  int salvagePods{0};

  // Some derelicts still contain an intact data core. The prototype represents
  // this as a small Electronics pod you can scoop.
  bool hasDataCore{false};
  econ::CommodityId dataCommodity{econ::CommodityId::Electronics};
  double dataUnits{0.0};

  // Optional hostile response ("salvage ambush").
  bool ambush{false};
  int pirateCount{0};

  // 0..1-ish risk score used only for tuning/telemetry/UI.
  double risk{0.0};
};

// Generate a deterministic derelict plan.
//
// Inputs:
//  - universeSeed: global seed
//  - systemId: used for locality
//  - signalId: unique per signal instance
//  - timeDays: only used when includeDayStamp==true (integer day mixed in)
//  - piracy/security/contest: 0..1 system knobs (see SecurityModel)
//  - missionSite: true when this derelict is spawned for a mission objective
//  - includeDayStamp: when true, mixes the integer day into the RNG seed. This
//    should be enabled for "daily" derelicts but disabled for mission sites so
//    their contents don't change mid-mission.
DerelictPlan planDerelictEncounter(core::u64 universeSeed,
                                   SystemId systemId,
                                   core::u64 signalId,
                                   double timeDays,
                                   double piracy01,
                                   double security01,
                                   double contest01,
                                   bool missionSite,
                                   bool includeDayStamp = true);

} // namespace stellar::sim
