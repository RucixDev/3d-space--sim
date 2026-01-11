#pragma once

#include "stellar/core/Types.h"

#include "stellar/sim/SecurityModel.h"
#include "stellar/sim/System.h"
#include "stellar/sim/SystemEvents.h"
#include "stellar/sim/SystemSecurityDynamics.h"

namespace stellar::sim {

// A structured view of the system's "conditions" layers.
//
// Layers:
//  1) Base (deterministic): systemSecurityProfile(universeSeed, sys)
//  2) Dynamics (persistent): decaying deltas applied to the base profile
//  3) Event (deterministic): time-based "weather" applied on top of base+dynamics
struct SystemConditionsSnapshot {
  SystemId systemId{0};

  // Deterministic baseline values.
  SystemSecurityProfile base{};

  // Base profile with the (optional) persistent/decaying delta applied.
  SystemSecurityProfile afterDynamics{};

  // The delta values after being decayed to `timeDays`.
  // (lastUpdateDay is set to `timeDays`.)
  SystemSecurityDeltaState dynamicsNow{};
  bool hasDynamics{false};

  // Deterministic system event (may be inactive).
  SystemEvent event{};

  // Final effective values used by gameplay systems.
  SystemSecurityProfile effective{};
};

// Compute a snapshot of system conditions.
//
// `deltaState` may be null if the system has no persistent dynamics.
SystemConditionsSnapshot snapshotSystemConditions(core::u64 universeSeed,
                                                  const StarSystem& sys,
                                                  double timeDays,
                                                  const SystemSecurityDeltaState* deltaState,
                                                  const SystemSecurityDynamicsParams& dynParams = {},
                                                  const SystemEventParams& evParams = {});

// Convenience: compute just the effective security profile.
// If `outEvent` is non-null, it receives the generated event.
SystemSecurityProfile effectiveSystemSecurityProfile(core::u64 universeSeed,
                                                     const StarSystem& sys,
                                                     double timeDays,
                                                     const SystemSecurityDeltaState* deltaState,
                                                     const SystemSecurityDynamicsParams& dynParams = {},
                                                     const SystemEventParams& evParams = {},
                                                     SystemEvent* outEvent = nullptr);

} // namespace stellar::sim
