#include "stellar/sim/SystemConditions.h"

#include "stellar/sim/SystemEvents.h"
#include "stellar/sim/SystemSecurityDynamics.h"

namespace stellar::sim {

SystemConditionsSnapshot snapshotSystemConditions(core::u64 universeSeed,
                                                  const StarSystem& sys,
                                                  double timeDays,
                                                  const SystemSecurityDeltaState* deltaState,
                                                  const SystemSecurityDynamicsParams& dynParams,
                                                  const SystemEventParams& evParams) {
  SystemConditionsSnapshot out{};
  out.systemId = sys.stub.id;

  out.base = systemSecurityProfile(universeSeed, sys);
  out.afterDynamics = out.base;
  out.hasDynamics = (deltaState != nullptr);

  if (deltaState) {
    out.dynamicsNow = decayedSystemSecurityDelta(*deltaState, timeDays, dynParams);
    out.afterDynamics = applySystemSecurityDelta(out.base, *deltaState, timeDays, dynParams);
  } else {
    out.dynamicsNow = SystemSecurityDeltaState{};
  }

  // Generate the event from the baseline-with-dynamics profile (so events can
  // react to persistent player-driven deltas).
  out.event = generateSystemEvent(universeSeed, sys.stub.id, timeDays, out.afterDynamics, evParams);
  out.effective = applySystemEventToProfile(out.afterDynamics, out.event);
  return out;
}

SystemSecurityProfile effectiveSystemSecurityProfile(core::u64 universeSeed,
                                                     const StarSystem& sys,
                                                     double timeDays,
                                                     const SystemSecurityDeltaState* deltaState,
                                                     const SystemSecurityDynamicsParams& dynParams,
                                                     const SystemEventParams& evParams,
                                                     SystemEvent* outEvent) {
  const auto snap = snapshotSystemConditions(universeSeed, sys, timeDays, deltaState, dynParams, evParams);
  if (outEvent) *outEvent = snap.event;
  return snap.effective;
}

} // namespace stellar::sim