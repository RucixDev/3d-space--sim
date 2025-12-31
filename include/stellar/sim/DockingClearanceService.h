#pragma once

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/core/Types.h"
#include "stellar/sim/System.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Docking clearance (headless)
// -----------------------------------------------------------------------------
//
// The SDL prototype historically handled docking clearance requests directly in
// apps/stellar_game/main.cpp using the game's global RNG.
//
// This caused clearance outcomes to be subtly influenced by unrelated random
// events (cosmetics, loot rolls, etc.). This module extracts the clearance rules
// into the core library and makes decisions deterministic per:
//   (universeSeed, stationId, requestCount)
//
// The app still owns the UI/toasts; this module just updates a small state blob.
//
// Behavior model:
//  - You can request clearance when within station.commsRangeKm.
//  - If granted, it remains valid for grantDurationSec.
//  - If denied, a denial cooldown applies for denyCooldownSec.
//  - Auto-request mode adds an additional retryIntervalSec throttle.
//
// traffic01 is a congestion hint in [0,1] (0 = empty, 1 = very busy). It reduces
// the probability of being granted.

struct DockingClearanceState {
  bool granted{false};
  double expiresDays{0.0};
  double cooldownUntilDays{0.0};
  double lastRequestDays{0.0};

  // Counts actual request attempts (in range, not throttled/cooldown).
  // Used to derive a stable per-request RNG seed.
  core::u32 requestCount{0};
};

enum class DockingClearanceStatus : core::u8 {
  OutOfRange = 0,
  Throttled  = 1, // cooldown or retry interval not elapsed
  Granted    = 2,
  Denied     = 3,
};

struct DockingClearanceDecision {
  DockingClearanceStatus status{DockingClearanceStatus::OutOfRange};

  // True if the player has a currently valid clearance (either pre-existing or
  // granted by this call).
  bool hasClearance{false};

  // The probability used for the grant roll (0..1). 0 for non-attempt outcomes.
  double pGrant{0.0};
};

struct DockingClearanceParams {
  // Baseline probability of clearance grant before traffic adjustments.
  double baseGrantProb{0.82};

  // Clamp bounds for the final probability.
  double minGrantProb{0.05};
  double maxGrantProb{0.98};

  // Validity window for a granted clearance.
  double grantDurationSec{12.0 * 60.0};

  // Denial cooldown window (simulates traffic / channel busy).
  double denyCooldownSec{90.0};

  // Auto-request throttle window.
  double retryIntervalSec{6.0};

  // Final pGrant = clamp(base + typeBias - traffic01 * trafficPenalty, min, max)
  double trafficPenalty{0.35};

  // If false, calling requestDockingClearance() while already granted does
  // nothing (and cannot revoke/replace clearance).
  bool allowRefresh{false};
};

inline bool dockingClearanceValid(const DockingClearanceState& s, double timeDays) {
  return s.granted && (timeDays <= s.expiresDays);
}

// Estimate a deterministic station congestion factor in [0,1].
//
// The app can pass nearbyShips (e.g. number of local contacts near the station)
// to nudge the congestion up/down.
//
// This function is intentionally stable and independent of the app's RNG.
// It's designed to be used for gameplay feel only (not a physical simulation).
double estimateDockingTraffic01(core::u64 universeSeed,
                                const Station& st,
                                double timeDays,
                                int nearbyShips = 0);

// Request docking clearance.
//
// Returns:
//  - OutOfRange: ship is outside comms range (state unchanged)
//  - Throttled: request is on cooldown (state unchanged)
//  - Granted: clearance granted and state updated
//  - Denied: clearance denied and cooldown applied
DockingClearanceDecision requestDockingClearance(core::u64 universeSeed,
                                                const Station& st,
                                                double timeDays,
                                                double distanceKm,
                                                DockingClearanceState& ioState,
                                                double traffic01 = 0.0,
                                                const DockingClearanceParams& params = DockingClearanceParams{});

// Auto-request docking clearance.
//
// Like requestDockingClearance() but also enforces the retryIntervalSec throttle
// so it can be called each frame by an autopilot without spamming.
DockingClearanceDecision autoRequestDockingClearance(core::u64 universeSeed,
                                                     const Station& st,
                                                     double timeDays,
                                                     double distanceKm,
                                                     DockingClearanceState& ioState,
                                                     double traffic01 = 0.0,
                                                     const DockingClearanceParams& params = DockingClearanceParams{});

} // namespace stellar::sim
