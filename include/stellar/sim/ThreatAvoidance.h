#pragma once

#include "stellar/core/Types.h"

#include "stellar/math/Vec3.h"

#include "stellar/sim/CollisionWarning.h"
#include "stellar/sim/MissileDefense.h"
#include "stellar/sim/ProximityField.h"
#include "stellar/sim/Ship.h"

#include <cstddef>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// ThreatAvoidance (headless)
// -----------------------------------------------------------------------------
//
// A tiny deterministic helper that blends two cockpit-level "oh no" signals:
//   1) Collision prediction along the current velocity vector (CollisionWarning)
//   2) Nearest inbound missile detection + lateral jink direction (MissileDefense)
//
// It outputs a *suggested* ShipInput overlay (thrustLocal + optional brake/boost)
// that game/UI layers can apply as:
//   - a soft assist (add and clamp)
//   - an autopilot override
//   - an NPC defense behavior
//
// The module is intentionally conservative and avoids global state.

struct ThreatAvoidanceParams {
  // --- Collision avoidance ---
  bool collisionEnable{true};
  CollisionWarningParams collision{};

  // Engage lateral jink when collision.hazard01 >= this.
  double collisionEngageHazard01{0.55};

  // Engage brake when collision.hazard01 >= this.
  double collisionBrakeEngageHazard01{0.85};

  // If the ship is moving slower than this, we don't bother producing a jink.
  double collisionMinSpeedForJinkKmS{0.08};

  // --- Missile evasion ---
  bool missileEnable{true};
  MissileThreatParams missileThreat{};
  MissileEvasionParams missileEvasion{};

  // Engage missile evasion when threat.ttiSec <= this.
  double missileEngageTtiSec{12.0};

  // --- Blending / weighting ---
  //
  // The collision weight is scaled by collision.hazard01.
  // The missile weight is scaled by a ramp based on threat.ttiSec.
  double collisionWeight{1.0};
  double missileWeight{1.0};

  // --- Output shaping ---
  // Clamp final thrust magnitude.
  double maxThrust01{1.0};

  // If true, project the final thrust direction into the plane perpendicular
  // to the current velocity. This tends to produce classic "side jinks" that
  // increase miss distance without accelerating deeper into a collision.
  bool preferLateralJink{true};

  // Optional boost recommendation when the situation is extremely urgent.
  bool allowBoost{true};
  double boostEngageHazard01{0.95};
  double boostEngageMissileTtiSec{3.0};
};

struct ThreatAvoidanceResult {
  bool active{false};
  bool collisionActive{false};
  bool missileActive{false};

  // World-space unit direction for the suggested thrust.
  math::Vec3d dirWorld{0, 0, 0};

  // Suggested thrust magnitude in [0,1].
  double thrust01{0.0};

  // Suggested control overlay.
  ShipInput input{};

  // Diagnostics for UI/telemetry.
  CollisionWarningResult collision{};
  MissileThreatSummary missileThreat{};
  MissileEvasionPlan missilePlan{};
};

// Compute a deterministic threat-avoidance suggestion for `ship`.
//
// - obstacles: optional; when provided, collision is predicted against its spheres.
// - missiles: optional; when provided, detects the nearest inbound missile that
//             targets (targetKind,targetId).
// - seed: used for deterministic tie-breaking / degenerate perpendicular picks.
//
// The returned ShipInput is intended as an overlay (caller may add to existing
// input and clamp).
ThreatAvoidanceResult computeThreatAvoidance(const Ship& ship,
                                            const ProximityFieldKm* obstacles,
                                            double maxDecelKmS2,
                                            const Missile* missiles,
                                            std::size_t missileCount,
                                            CombatTargetKind targetKind,
                                            core::u64 targetId,
                                            core::u64 seed,
                                            const ThreatAvoidanceParams& params = {},
                                            int ignoreObstacleId = -1);

} // namespace stellar::sim
