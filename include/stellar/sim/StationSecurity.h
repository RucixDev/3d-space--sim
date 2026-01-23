#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Law.h"
#include "stellar/sim/System.h"

#include <optional>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Station security / traffic control (headless)
// -----------------------------------------------------------------------------
//
// A lightweight enforcement model for the station "no-fire" zone and basic
// approach rules.
//
// Design goals:
//  - Deterministic and UI-agnostic (usable by tests and tools).
//  - Avoid per-frame spam (accumulators + cooldowns).
//  - Graduated response: warning -> fine -> bounty.
//
// The app is responsible for actually applying fines/bounties and for displaying
// toasts/comms. This module only produces recommended enforcement events.

enum class StationOffenseKind : core::u8 {
  None = 0,
  WeaponDischarge,
  Speeding,
  Trespass,
};

enum class StationOffenseAction : core::u8 {
  None = 0,
  Warning,
  Fine,
  Bounty,
};

const char* stationOffenseKindName(StationOffenseKind k);
const char* stationOffenseActionName(StationOffenseAction a);

struct StationOffenseEvent {
  StationOffenseKind kind{StationOffenseKind::None};
  StationOffenseAction action{StationOffenseAction::None};

  // 0..1 severity hint (used for tuning and UI copy).
  double severity01{0.0};

  // Recommended numbers (0 if not applicable).
  double fineCr{0.0};
  double bountyCr{0.0};
  double repPenalty{0.0}; // negative values penalize reputation

  // Strike counter AFTER this event is applied.
  int strikes{0};

  // Optional telemetry (mainly for speeding).
  double speedLimitKmS{0.0};
  double measuredSpeedKmS{0.0};
};

struct StationOffenseTracker {
  int strikes{0};
  double lastStrikeDays{0.0};
  double cooldownUntilDays{0.0};
};

struct StationSecurityState {
  StationOffenseTracker weapons{};
  StationOffenseTracker speeding{};
  StationOffenseTracker trespass{};

  // Accumulators for continuous offenses.
  double speedingAccumSec{0.0};
  double trespassAccumSec{0.0};
};

struct StationSecurityParams {
  // Zone geometry.
  double noFireZoneRadiusFactor{25.0};
  double speedZoneRadiusFactor{12.0};

  // Speed envelope: speedLimit is interpolated between near/far multipliers.
  double speedZoneInnerFactor{2.0};
  double speedLimitNearMult{1.10};
  double speedLimitFarMult{7.0};

  // Allow a little overspeed before we start accumulating a speeding event.
  double speedToleranceFrac{0.18};

  // Strike timing.
  double strikeResetSec{12.0 * 60.0};

  // Cooldowns / trigger durations.
  double weaponEventCooldownSec{1.25};
  double genericCooldownSec{28.0};
  double speedingTriggerSec{2.0};
  double trespassTriggerSec{1.0};

  // Escalation thresholds.
  int strikesToFine{2};
  int strikesToBounty{3};

  // Baseline fine schedule (law profile further modulates).
  double baseFineWeaponCr{550.0};
  double baseFineSpeedingCr{180.0};
  double baseFineTrespassCr{420.0};

  // Baseline bounty schedule (law profile further modulates).
  double baseBountyWeaponCr{950.0};
  double baseBountyTrespassCr{700.0};
};

// Convenience helpers (pure, deterministic).
double stationNoFireZoneRadiusKm(const Station& st, const StationSecurityParams& params = {});
double stationSpeedZoneRadiusKm(const Station& st, const StationSecurityParams& params = {});
double stationSpeedLimitKmS(const Station& st, const StationSecurityParams& params, double distanceKm);

// True if the ship is inside the mail-slot tunnel volume (station-local space).
// This is intentionally a loose geometric test (no clearance / alignment checks).
bool insideStationSlotTunnel(const Station& st, const math::Vec3d& relLocalKm);

// Evaluate station security rules.
//
// Returns at most one event per call (priority: weapons > trespass > speeding)
// to keep UI/ledger spam manageable.
std::optional<StationOffenseEvent> updateStationSecurity(StationSecurityState& ioState,
                                                        const StationSecurityParams& params,
                                                        const Station& st,
                                                        const LawProfile& law,
                                                        const math::Vec3d& stationPosKm,
                                                        const math::Vec3d& stationVelKmS,
                                                        const math::Quatd& stationOrient,
                                                        const math::Vec3d& shipPosKm,
                                                        const math::Vec3d& shipVelKmS,
                                                        bool hasClearance,
                                                        bool shipDocked,
                                                        bool shipFiredWeaponThisTick,
                                                        double timeDays,
                                                        double dtSimSec);

} // namespace stellar::sim
