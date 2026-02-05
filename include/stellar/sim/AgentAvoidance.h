#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"

#include <span>

namespace stellar::sim {

// Deterministic local collision-avoidance steering against *moving* spherical agents.
//
// This helper is intentionally lightweight: it does not do global path planning.
// It biases a desired travel direction away from neighbors by predicting the
// closest approach between two agents under constant-velocity motion over a
// finite horizon. The output is a unit direction suitable for turning a point
// goal into an "avoidance-biased" chase point.
//
// The implementation is inspired by the velocity-obstacle family of techniques
// (RVO/ORCA), but keeps the math simple for a game/simulation context.

struct AgentSphere {
  core::u64 id{0};

  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};

  // Radius of the agent (km).
  double radiusKm{0.0};

  // Hardness (0..1) acts as a steering weight: 0 = ignored, 1 = full repulsion.
  double hardness01{1.0};
};

struct AgentAvoidanceParams {
  bool enabled{true};

  // Inflates every neighbor sphere by (selfRadiusKm + padKm) before computing clearance.
  double selfRadiusKm{0.03};
  double padKm{0.05};

  // Horizon for closest-approach prediction.
  double horizonSec{12.0};

  // Begin steering when predicted clearance falls below this (km).
  double nearMissExtraKm{0.35};

  // Treat very small relative speeds as static.
  double minRelSpeedKmS{0.01};

  // Blend factor [0,1] between current velocity and the commanded velocity
  // implied by (desiredDir, desiredSpeed) for prediction.
  double selfVelBlend01{0.55};

  // Unitless steering strength.
  double strength{1.25};

  // Clamp steering cone relative to desiredDirUnitWorld.
  double maxSteerDeg{35.0};
};

struct AgentAvoidanceResult {
  math::Vec3d desiredDirUnit{0, 0, 1};
  math::Vec3d safeDirUnit{0, 0, 1};
  bool steering{false};

  core::u64 threatId{0};
  double threatTtiSec{0.0};
  double threatClearanceKm{0.0};     // negative => intersects inflated spheres
  double threatDistAtClosestKm{0.0}; // center-to-center distance at closest approach
};

// Compute an avoidance-biased travel direction given dynamic neighbors.
//
// desiredDirUnitWorld:
//   Desired travel direction (world space). It is normalized internally.
//
// desiredSpeedKmS:
//   Used only for a rough prediction of our near-future motion. It does not clamp
//   the output.
//
// seed:
//   Used solely to break symmetry in perfectly head-on or overlapping cases.
AgentAvoidanceResult steerAvoidAgents(const math::Vec3d& selfPosKm,
                                     const math::Vec3d& selfVelKmS,
                                     const math::Vec3d& desiredDirUnitWorld,
                                     double desiredSpeedKmS,
                                     std::span<const AgentSphere> neighbors,
                                     core::u64 seed,
                                     const AgentAvoidanceParams& params);

} // namespace stellar::sim
