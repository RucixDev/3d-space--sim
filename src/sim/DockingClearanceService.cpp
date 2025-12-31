#include "stellar/sim/DockingClearanceService.h"

#include "stellar/core/Hash.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static double clamp01(double x) {
  return std::clamp(x, 0.0, 1.0);
}

static double stationCapacityFactor(econ::StationType type) {
  using econ::StationType;
  switch (type) {
    case StationType::Outpost:        return 0.75;
    case StationType::Agricultural:   return 0.95;
    case StationType::Mining:         return 0.90;
    case StationType::Refinery:       return 1.00;
    case StationType::Industrial:     return 1.05;
    case StationType::Research:       return 0.95;
    case StationType::TradeHub:       return 1.25;
    case StationType::Shipyard:       return 1.20;
    default:                          return 1.00;
  }
}

static double stationTypeGrantBias(econ::StationType type) {
  using econ::StationType;
  switch (type) {
    case StationType::Outpost:   return -0.06;
    case StationType::Research:  return -0.02;
    case StationType::TradeHub:  return +0.05;
    case StationType::Shipyard:  return +0.04;
    default:                     return 0.0;
  }
}

double estimateDockingTraffic01(core::u64 universeSeed,
                                const Station& st,
                                double timeDays,
                                int nearbyShips) {
  constexpr double kTwoPi = 6.283185307179586476925286766559;

  const core::u64 salt = core::fnv1a64("DockingTraffic");
  core::SplitMix64 rng(core::hashCombine(core::hashCombine(universeSeed, salt), (core::u64)st.id));

  // Baseline congestion curve (smooth daily cycle).
  const double base = rng.range(0.18, 0.48); // average
  const double amp = rng.range(0.08, 0.28);  // daily swing
  const double phase = rng.range(0.0, kTwoPi);

  const double cyc = 0.5 + 0.5 * std::sin(kTwoPi * timeDays + phase); // 0..1
  const double baseline = clamp01(base + (cyc - 0.5) * amp);

  // Nearby traffic nudges congestion upward relative to station capacity.
  const double cap = stationCapacityFactor(st.type);
  const double crowd = clamp01(static_cast<double>(std::max(0, nearbyShips)) / (8.0 * cap));

  return clamp01(baseline * 0.70 + crowd * 0.55);
}

DockingClearanceDecision requestDockingClearance(core::u64 universeSeed,
                                                const Station& st,
                                                double timeDays,
                                                double distanceKm,
                                                DockingClearanceState& ioState,
                                                double traffic01,
                                                const DockingClearanceParams& params) {
  DockingClearanceDecision out{};

  const bool valid = dockingClearanceValid(ioState, timeDays);
  if (valid && !params.allowRefresh) {
    out.status = DockingClearanceStatus::Granted;
    out.hasClearance = true;
    return out;
  }

  if (distanceKm > st.commsRangeKm) {
    out.status = DockingClearanceStatus::OutOfRange;
    out.hasClearance = valid;
    return out;
  }

  if (timeDays < ioState.cooldownUntilDays) {
    out.status = DockingClearanceStatus::Throttled;
    out.hasClearance = valid;
    return out;
  }

  const double traffic = clamp01(traffic01);
  double p = params.baseGrantProb + stationTypeGrantBias(st.type) - traffic * params.trafficPenalty;
  p = std::clamp(p, params.minGrantProb, params.maxGrantProb);
  out.pGrant = p;

  // Stable seed per request attempt.
  const core::u64 salt = core::fnv1a64("DockingClearance");
  core::u64 s = core::hashCombine(universeSeed, salt);
  s = core::hashCombine(s, (core::u64)st.id);
  s = core::hashCombine(s, (core::u64)ioState.requestCount);

  core::SplitMix64 rng(s);

  ioState.lastRequestDays = timeDays;
  ioState.requestCount += 1;

  const bool grant = rng.chance(p);
  if (grant) {
    ioState.granted = true;
    ioState.expiresDays = timeDays + params.grantDurationSec / 86400.0;
    ioState.cooldownUntilDays = timeDays;

    out.status = DockingClearanceStatus::Granted;
    out.hasClearance = true;
    return out;
  }

  // Denied: apply cooldown and revoke.
  ioState.granted = false;
  ioState.expiresDays = 0.0;
  ioState.cooldownUntilDays = timeDays + params.denyCooldownSec / 86400.0;

  out.status = DockingClearanceStatus::Denied;
  out.hasClearance = false;
  return out;
}

DockingClearanceDecision autoRequestDockingClearance(core::u64 universeSeed,
                                                     const Station& st,
                                                     double timeDays,
                                                     double distanceKm,
                                                     DockingClearanceState& ioState,
                                                     double traffic01,
                                                     const DockingClearanceParams& params) {
  DockingClearanceDecision out{};

  const bool valid = dockingClearanceValid(ioState, timeDays);
  if (valid && !params.allowRefresh) {
    out.status = DockingClearanceStatus::Granted;
    out.hasClearance = true;
    return out;
  }

  if (distanceKm > st.commsRangeKm) {
    out.status = DockingClearanceStatus::OutOfRange;
    out.hasClearance = valid;
    return out;
  }

  if (timeDays < ioState.cooldownUntilDays) {
    out.status = DockingClearanceStatus::Throttled;
    out.hasClearance = valid;
    return out;
  }

  const double retryDays = params.retryIntervalSec / 86400.0;
  if (retryDays > 0.0 && (timeDays - ioState.lastRequestDays) < retryDays) {
    out.status = DockingClearanceStatus::Throttled;
    out.hasClearance = valid;
    return out;
  }

  return requestDockingClearance(universeSeed, st, timeDays, distanceKm, ioState, traffic01, params);
}

} // namespace stellar::sim
