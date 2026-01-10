#pragma once

#include "stellar/sim/Celestial.h"

namespace stellar::sim {

// -----------------------------------------------------------------------------
// System Security Dynamics (lightweight, persistent, headless)
// -----------------------------------------------------------------------------
//
// The deterministic SystemSecurityProfile (see sim/SecurityModel) gives a stable
// baseline per system. This module adds a small *persistent delta* per system so
// player actions can nudge local conditions (security / piracy / traffic) over
// time.
//
// Design goals:
//   - save-friendly (plain POD)
//   - deterministic (no RNG)
//   - robust (clamps + exponential decay back toward baseline)
//
// The deltas are intended to be applied additively:
//   effective = clamp01(baseline + delta)
// where delta decays exponentially toward 0.

struct SystemSecurityDeltaState {
  SystemId systemId{0};

  // Additive deltas applied to SystemSecurityProfile fields.
  // Typical range is small (e.g. +/- 0.05). Values are clamped when impulses are applied.
  double securityDelta{0.0};
  double piracyDelta{0.0};
  double trafficDelta{0.0};

  // Timestamp in simulation days when this delta was last updated/decayed.
  double lastUpdateDay{0.0};
};

struct SystemSecurityDynamicsParams {
  // Half-life (in days) for exponential decay toward baseline.
  // Example: halfLife=10 means the delta halves every ~10 in-game days.
  double securityHalfLifeDays{10.0};
  double piracyHalfLifeDays{10.0};
  double trafficHalfLifeDays{10.0};

  // Clamp deltas to avoid extreme values.
  double maxAbsDelta{0.35};

  // Deltas smaller than this (after decay-to-now) can be pruned.
  double negligibleAbs{1e-4};
};

// Exponential decay factor in [0,1].
//
// halfLifeDays:
//  - dt == halfLife => factor ~ 0.5
//  - dt == 0        => factor == 1
//  - halfLife <= 0  => factor == 0 for dt > 0
double decayFactorDays(double dtDays, double halfLifeDays);

// Return a decayed copy of the delta state at nowDays.
SystemSecurityDeltaState decayedSystemSecurityDelta(const SystemSecurityDeltaState& st,
                                                    double nowDays,
                                                    const SystemSecurityDynamicsParams& params);

// Decay the delta state in-place to nowDays.
void decaySystemSecurityDeltaInPlace(SystemSecurityDeltaState& st,
                                     double nowDays,
                                     const SystemSecurityDynamicsParams& params);

// Apply an additive impulse at nowDays (decays-to-now first, then adds, then clamps).
void applySystemSecurityImpulse(SystemSecurityDeltaState& st,
                                double nowDays,
                                double dSecurity,
                                double dPiracy,
                                double dTraffic,
                                const SystemSecurityDynamicsParams& params);

// True if the delta is small enough (after decaying to nowDays) that it can be discarded.
bool isSystemSecurityDeltaNegligible(const SystemSecurityDeltaState& st,
                                     double nowDays,
                                     const SystemSecurityDynamicsParams& params);

// Forward declaration (defined in sim/SecurityModel.h).
struct SystemSecurityProfile;

// Apply a (decayed) delta to a baseline security profile.
SystemSecurityProfile applySystemSecurityDelta(const SystemSecurityProfile& base,
                                               const SystemSecurityDeltaState& st,
                                               double nowDays,
                                               const SystemSecurityDynamicsParams& params);

} // namespace stellar::sim
