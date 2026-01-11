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

  // Additional programs used by UI/tools.
  //
  // Notes:
  //  - These are *planning* helpers (impulsive/Keplerian). Execution is still
  //    handled by sim::ManeuverComputer.
  //  - Some programs require extra inputs (target radius, etc.), so the enum is
  //    mainly used for reporting/debug.
  SetApoapsisAtPeriapsis = 2,
  SetPeriapsisAtApoapsis = 3,
  EscapeNow = 4,
  AlignPlaneAtAscendingNode = 5,
  AlignPlaneAtDescendingNode = 6,
};

struct ManeuverProgramParams {
  // If true, multiply the reference body's mu by gravityScale.
  // (GravityParams::scale scales accelerations; for two-body math this is
  // equivalent to scaling mu.)
  bool applyGravityScale{true};

  // If |dv| is below this threshold, treat the plan as a no-op.
  //
  // For Circularize plans, this also avoids the circular-orbit degeneracy
  // (periapsis/apoapsis are not unique).
  double minDvKmS{1e-6}; // 1 mm/s
};

struct ManeuverProgramResult {
  bool valid{false};
  std::string reason{};

  ManeuverProgramKind program{ManeuverProgramKind::CircularizeAtPeriapsis};

  GravityBody refBody{};
  ManeuverPlan plan{};

  // Diagnostics
  double timeToNodeSec{0.0};
  double dvKmS{0.0};
  double nodeRadiusKm{0.0};
  double nodeAltitudeKm{0.0};
  double nodeSpeedKmS{0.0};

  // Program-specific target at the node (if applicable).
  // - Circularize: targetSpeedKmS == circular speed at nodeRadiusKm
  // - SetApo/SetPeri: targetRadiusKm == requested opposite apsis radius
  // - Escape: targetSpeedKmS == escape speed at nodeRadiusKm
  double targetRadiusKm{0.0};
  double targetSpeedKmS{0.0};

  // Back-compat name for the circularize program (filled when applicable).
  double targetCircularSpeedKmS{0.0};

  // Body-centric state at the burn node (relative to refBody). Useful for
  // computing RTN bases in UI.
  math::Vec3d nodeRelPosKm{0, 0, 0};
  math::Vec3d nodeRelVelKmS{0, 0, 0};
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

// Plan a tangential burn at periapsis to set the orbit's apoapsis radius.
//
// `targetApoapsisKm` is clamped to be >= periapsis radius and >= refBody.radiusKm.
ManeuverProgramResult planSetApoapsisAtPeriapsis(const Ship& ship,
                                                 double nowTimeDays,
                                                 const GravityBody& refBody,
                                                 double targetApoapsisKm,
                                                 double gravityScale = 1.0,
                                                 const ManeuverProgramParams& params = {});

// Plan a tangential burn at apoapsis to set the orbit's periapsis radius.
//
// `targetPeriapsisKm` is clamped to be <= apoapsis radius and >= refBody.radiusKm.
ManeuverProgramResult planSetPeriapsisAtApoapsis(const Ship& ship,
                                                 double nowTimeDays,
                                                 const GravityBody& refBody,
                                                 double targetPeriapsisKm,
                                                 double gravityScale = 1.0,
                                                 const ManeuverProgramParams& params = {});

// Plan a prograde burn at the current position to reach escape energy.
// If the current speed is already >= escape speed (within minDvKmS), returns a
// zero-dv plan at now.
ManeuverProgramResult planEscapeNow(const Ship& ship,
                                    double nowTimeDays,
                                    const GravityBody& refBody,
                                    double gravityScale = 1.0,
                                    const ManeuverProgramParams& params = {});

// Plan an inclination/plane change burn to align the orbit plane to the
// reference XY plane (normal ±Z).
//
// - ascendingNode: true => ascending node, false => descending node.
// - forcePrograde: if false, will preserve retrograde orientation (target normal
//   will match the current orbit normal sign).
ManeuverProgramResult planAlignPlaneAtNode(const Ship& ship,
                                           double nowTimeDays,
                                           const GravityBody& refBody,
                                           bool ascendingNode,
                                           bool forcePrograde = false,
                                           double gravityScale = 1.0,
                                           const ManeuverProgramParams& params = {});

inline ManeuverProgramResult planAlignPlaneAtAscendingNode(const Ship& ship,
                                                           double nowTimeDays,
                                                           const GravityBody& refBody,
                                                           bool forcePrograde = false,
                                                           double gravityScale = 1.0,
                                                           const ManeuverProgramParams& params = {}) {
  return planAlignPlaneAtNode(ship, nowTimeDays, refBody, /*ascendingNode=*/true, forcePrograde, gravityScale, params);
}

inline ManeuverProgramResult planAlignPlaneAtDescendingNode(const Ship& ship,
                                                            double nowTimeDays,
                                                            const GravityBody& refBody,
                                                            bool forcePrograde = false,
                                                            double gravityScale = 1.0,
                                                            const ManeuverProgramParams& params = {}) {
  return planAlignPlaneAtNode(ship, nowTimeDays, refBody, /*ascendingNode=*/false, forcePrograde, gravityScale, params);
}

// Convenience overload: chooses the dominant gravity body as the reference.
// Returns valid=false if no bodies were enabled.
ManeuverProgramResult planCircularizeDominantBody(const StarSystem& sys,
                                                  double nowTimeDays,
                                                  const Ship& ship,
                                                  const GravityParams& gravityParams,
                                                  ManeuverProgramKind kind,
                                                  const ManeuverProgramParams& params = {});

} // namespace stellar::sim
