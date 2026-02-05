#pragma once

#include "stellar/math/Bvh.h"
#include "stellar/math/Vec3.h"

#include <cstddef>
#include <vector>
#include <utility>

namespace stellar::sim {

// A simple spherical obstacle used for proximity queries and steering.
// Units: center in km, radius in km.
struct SphereObstacleKm {
  math::Vec3d centerKm{0, 0, 0};
  double radiusKm{0.0};

  // Hardness (0..1) acts as a steering weight: 0 = ignored by avoidance,
  // 1 = full repulsion.
  double hardness01{1.0};
};

struct ProximityHitKm {
  bool hit{false};
  int id{-1};

  // Distance along the ray to the first intersection point.
  double tKm{0.0};

  // If the hit came from a velocity prediction helper, this is tKm / |vel|.
  // Otherwise it is left at 0.
  double tSec{0.0};

  math::Vec3d pointKm{0, 0, 0};

  // Snapshot of the obstacle at the time of the query.
  double obstacleRadiusKm{0.0};
  double obstacleHardness01{1.0};
};

// BVH-backed proximity queries over spherical obstacles.
//
// This is meant for “cockpit-level” helpers:
//  - collision warning (time-to-impact along current velocity)
//  - local obstacle avoidance steering (nav assist)
//
// The field is deterministic and safe for headless / server-side use.
class ProximityFieldKm {
public:
  void clear();
  bool empty() const { return obstacles_.empty(); }
  std::size_t obstacleCount() const { return obstacles_.size(); }

  // Build a BVH over the provided obstacles.
  //
  // Important: the field preserves obstacle ordering. The BVH stores obstacle
  // indices as IDs, so callers can safely use their own vector indices as
  // ignore IDs (e.g., ignore the target asteroid by index).
  void build(std::vector<SphereObstacleKm> obstacles, std::size_t leafSize = 6);

  const SphereObstacleKm* obstacle(int id) const;

  // Broadphase query helper.
  template <class Fn>
  void queryAabb(const math::Aabb3d& query, Fn&& onHit) const {
    bvh_.queryAabb(query, std::forward<Fn>(onHit));
  }

  // Broadphase frustum query helper (returns obstacle IDs whose AABBs intersect the frustum).
  template <class Fn>
  void queryFrustum(const math::Frustumd& frustum, Fn&& onHit) const {
    bvh_.queryFrustum(frustum, std::forward<Fn>(onHit));
  }

  // Raycast against obstacles (inflated by padKm). Returns the closest hit along
  // [0, maxDistKm].
  ProximityHitKm raycastFirstHit(const math::Vec3d& originKm,
                                const math::Vec3d& dirUnitWorld,
                                double maxDistKm,
                                double padKm = 0.0,
                                int ignoreId = -1,
                                int maxSphereTests = 256) const;

  // Predict the earliest impact along linear motion with constant velocity.
  //
  // - velKmS: world-space velocity
  // - maxTimeSec: horizon
  // - padKm: inflates obstacles (use ship radius + safety margin)
  ProximityHitKm predictLinearImpact(const math::Vec3d& originKm,
                                     const math::Vec3d& velKmS,
                                     double maxTimeSec,
                                     double padKm = 0.0,
                                     int ignoreId = -1,
                                     int maxSphereTests = 256) const;

private:
  static bool raySphereFirstHitKm_(const math::Vec3d& originKm,
                                  const math::Vec3d& dirUnitWorld,
                                  const math::Vec3d& centerKm,
                                  double radiusKm,
                                  double& outTKm);

  std::vector<SphereObstacleKm> obstacles_;
  math::Bvh3d bvh_;
};

struct AvoidanceParamsKm {
  bool enabled{true};

  // Inflates all obstacles by (shipRadiusKm + padKm) before computing clearance.
  double shipRadiusKm{0.03};
  double padKm{0.05};

  // Lookahead distance = base + speed * lookaheadTimeSec.
  double lookaheadTimeSec{8.0};
  double lookaheadBaseKm{0.0};
  double minSpeedForLookaheadKmS{0.05};

  // Begin steering when clearance to the inflated obstacle falls below this.
  // (0 => steer only when intersecting the inflated sphere.)
  double nearMissExtraKm{0.35};

  // Unitless steering strength. Higher values bias harder away from obstacles.
  double strength{1.25};

  // Clamp the steering cone relative to the desired direction.
  double maxSteerDeg{35.0};
};

struct AvoidanceResultKm {
  math::Vec3d desiredDirUnit{0, 0, 1};
  math::Vec3d safeDirUnit{0, 0, 1};
  bool steering{false};

  // Most threatening obstacle (smallest clearance along the lookahead segment).
  int threatId{-1};
  double threatAlongKm{0.0};
  double threatClearanceKm{0.0}; // negative => intersects inflated sphere
  double lookaheadKm{0.0};
};



struct TangentBypassParamsKm {
  bool enabled{true};

  // Inflates obstacles by (shipRadiusKm + padKm) and additionally by
  // nearMissExtraKm to keep a comfortable clearance while bypassing.
  double shipRadiusKm{0.03};
  double padKm{0.05};
  double nearMissExtraKm{0.35};

  // Place the waypoint beyond the obstacle exit by:
  //   tBeyond = tExit + clearanceKm * aheadClearanceMult + aheadExtraKm,
  // clamped to the remaining distance to the goal.
  double aheadClearanceMult{1.0};
  double aheadExtraKm{0.0};

  // Optional clamp on distance from start -> waypoint (0 = no clamp).
  double maxWaypointDistKm{0.0};
};

struct TangentBypassResultKm {
  bool used{false};
  math::Vec3d waypointKm{0, 0, 0};

  int obstacleId{-1};

  // -1 or +1 (relative to the planner's internal stable frame), 0 if unused.
  int sideSign{0};
};

// Plan a local “tangent bypass” waypoint around the nearest blocking obstacle
// on the straight-line segment [pos -> goal].
//
// This is NOT global path planning: it returns a single waypoint that helps the
// caller steer around a spherical obstacle deterministically without jitter.
//
// preferredSideSign can be set to -1/+1 to keep a stable side choice across
// frames. Pass 0 to let the planner pick the best side.
TangentBypassResultKm planTangentBypassWaypoint(const ProximityFieldKm& field,
                                               const math::Vec3d& posKm,
                                               const math::Vec3d& goalKm,
                                               const TangentBypassParamsKm& params,
                                               int ignoreId = -1,
                                               int preferredSideSign = 0,
                                               int maxSphereTests = 256);

// Compute a deterministic steering direction that biases away from obstacles.
//
// This is a lightweight "potential field" style helper: it improves local
// robustness without doing global path planning.
AvoidanceResultKm steerAvoidObstacles(const ProximityFieldKm& field,
                                     const math::Vec3d& posKm,
                                     const math::Vec3d& velKmS,
                                     const math::Vec3d& desiredDirUnitWorld,
                                     double desiredSpeedKmS,
                                     const AvoidanceParamsKm& params,
                                     int ignoreId = -1);

} // namespace stellar::sim
