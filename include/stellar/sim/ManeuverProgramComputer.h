#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/Gravity.h"
#include "stellar/sim/ManeuverComputer.h"
#include "stellar/sim/System.h"

#include <string>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// ManeuverProgramComputer (headless)
// -----------------------------------------------------------------------------
//
// The ManeuverComputer executes an *instantaneous* delta-v as a centered,
// continuous burn.
//
// The ManeuverProgramComputer computes those maneuver nodes from high-level
// goals.
//
// This is intentionally conservative: it assumes Keplerian two-body dynamics
// around a single reference body. That makes it fast and predictable for tools,
// unit tests, and early gameplay.

enum class ManeuverProgramKind : core::u8 {
  CircularizeAtPeriapsis = 0,
  CircularizeAtApoapsis = 1,
};

struct ManeuverProgramParams {
  // If true, multiply the reference body's mu by gravityScale.
  // (GravityParams::scale scales accelerations; for two-body math this is
  // equivalent to scaling mu.)
  bool applyGravityScale{true};

  // If |dv| is below this threshold, treat the plan as "already circular".
  double minDvKmS{1e-6}; // 1 mm/s
};

struct ManeuverProgramResult {
  bool valid{false};
  std::string reason{};

  GravityBody refBody{};
  ManeuverPlan plan{};

  // Diagnostics
  double timeToNodeSec{0.0};
  double dvKmS{0.0};
  double nodeRadiusKm{0.0};
  double nodeAltitudeKm{0.0};
  double nodeSpeedKmS{0.0};
  double targetCircularSpeedKmS{0.0};
};

// Build a single maneuver plan for a circularization burn relative to `refBody`.
//
// - For elliptic orbits: targets apoapsis/periapsis based on `kind`.
// - For already-circular orbits (within minDvKmS): returns a zero-dv plan at now.
ManeuverProgramResult planCircularize(const Ship& ship,
                                      double nowTimeDays,
                                      const GravityBody& refBody,
                                      ManeuverProgramKind kind,
                                      double gravityScale = 1.0,
                                      const ManeuverProgramParams& params = {});

// Convenience overload: chooses the dominant gravity body as the reference.
// Returns valid=false if no bodies were enabled.
ManeuverProgramResult planCircularizeDominantBody(const StarSystem& sys,
                                                  double nowTimeDays,
                                                  const Ship& ship,
                                                  const GravityParams& gravityParams,
                                                  ManeuverProgramKind kind,
                                                  const ManeuverProgramParams& params = {});

} // namespace stellar::sim
