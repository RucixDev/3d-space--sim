#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Ship.h"
#include "stellar/sim/ShipLoadout.h"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Combat (headless helpers)
// -----------------------------------------------------------------------------
//
// The SDL prototype historically kept simple combat utilities (damage application,
// ray hits, projectile stepping) inside apps/stellar_game/main.cpp.
//
// This module extracts those mechanics into the core library so:
//  - tests can validate combat math deterministically
//  - future apps (stellar_sandbox, server sims) can reuse the same rules
//  - main.cpp stays focused on orchestration/rendering
//
// Design goals:
//  - no renderer/UI dependencies
//  - deterministic + easy to reason about
//  - stable API that can grow (status effects, resistances, etc.)

enum class CombatTargetKind : core::u8 {
  Ship = 0,
  Asteroid = 1,
  Player = 2,
};

struct SphereTarget {
  CombatTargetKind kind{CombatTargetKind::Ship};
  std::size_t index{0};
  core::u64 id{0};

  math::Vec3d centerKm{0, 0, 0};
  double radiusKm{1.0};

  // Optional aim-cone filter for soft aim assist.
  // If >= -0.5, a target is only considered when
  //   dot(rayDir, (center-origin).normalized()) >= minAimCos.
  // Use -1.0 to disable the filter.
  double minAimCos{-1.0};
};

struct RaycastHit {
  bool hit{false};
  double tKm{0.0}; // distance along the ray to first intersection
  math::Vec3d pointKm{0, 0, 0};

  CombatTargetKind kind{CombatTargetKind::Ship};
  std::size_t index{0};
  core::u64 id{0};
};

// Returns true if the ray intersects the sphere. If so, outTEnterKm is set to the
// distance along the ray to the entry point (clamped to >= 0).
bool raySphereIntersectKm(const math::Vec3d& originKm,
                          const math::Vec3d& dirNormalized,
                          const math::Vec3d& centerKm,
                          double radiusKm,
                          double& outTEnterKm);

// Raycast against a list of spherical targets and return the nearest intersection
// within maxRangeKm. Targets may specify an optional minAimCos filter.
RaycastHit raycastNearestSphereKm(const math::Vec3d& originKm,
                                  const math::Vec3d& dirNormalized,
                                  double maxRangeKm,
                                  const SphereTarget* targets,
                                  std::size_t targetCount);

// Segment-sphere hit test (useful for fast-moving projectiles).
bool segmentHitsSphereKm(const math::Vec3d& aKm,
                         const math::Vec3d& bKm,
                         const math::Vec3d& centerKm,
                         double radiusKm);

// Apply damage to (shield, hull) with shields absorbing first.
inline void applyDamage(double dmg, double& shield, double& hull) {
  dmg = std::max(0.0, dmg);
  if (shield > 0.0) {
    const double s = std::min(shield, dmg);
    shield -= s;
    dmg -= s;
  }
  if (dmg > 0.0) {
    hull -= dmg;
  }
}

// Visual event for beam-style weapons. Units are kilometers (sim space).
struct BeamEvent {
  math::Vec3d aKm{0, 0, 0};
  math::Vec3d bKm{0, 0, 0};
  float r{1.0f}, g{1.0f}, b{1.0f};
};

// Ballistic projectile (kinetic cannons / slugs). Units are kilometers.
struct Projectile {
  math::Vec3d prevKm{0, 0, 0};
  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};

  float r{1.0f}, g{1.0f}, b{1.0f};

  double ttlSimSec{0.0};
  double radiusKm{450.0};
  double dmg{0.0};

  bool fromPlayer{false};
  core::u64 shooterId{0};
};

struct ProjectileHit {
  CombatTargetKind kind{CombatTargetKind::Ship};
  std::size_t targetIndex{0};
  core::u64 targetId{0};

  math::Vec3d pointKm{0, 0, 0};

  double dmg{0.0};
  bool fromPlayer{false};
  core::u64 shooterId{0};
};

// Move projectiles forward by dtSim seconds, perform segment collision checks
// against spherical targets, and emit hit events.
//
// The caller decides what each SphereTarget represents (contacts, asteroids,
// player) and what to do with hits.
void stepProjectiles(std::vector<Projectile>& projectiles,
                     double dtSim,
                     const SphereTarget* targets,
                     std::size_t targetCount,
                     std::vector<ProjectileHit>& outHits);

// Weapon firing output. This does not apply damage; it only emits the geometry
// and hit metadata required for gameplay code to resolve the effect.
struct FireResult {
  bool fired{false};

  // Cooldown handling (sim seconds). newCooldownSimSec is either the unchanged
  // cooldown (if not fired) or weaponDef.cooldownSimSec (if fired).
  double newCooldownSimSec{0.0};

  // Heat added for this shot (caller applies clamping).
  double heatDelta{0.0};

  // Beam/projectile payload.
  bool hasBeam{false};
  BeamEvent beam{};

  bool hasProjectile{false};
  Projectile projectile{};

  // For beams: whether something was hit and, if so, which target.
  bool hit{false};
  CombatTargetKind hitKind{CombatTargetKind::Ship};
  std::size_t hitIndex{0};
  core::u64 hitId{0};
  math::Vec3d hitPointKm{0, 0, 0};
};

// Try to fire the selected weapon.
//
// - cooldownSimSec is the current remaining cooldown (sim seconds).
// - distributorMk is used to scale heat generation (better distributors run cooler).
// - shooterId is copied into spawned projectiles so collision code can avoid
//   immediate self-hits.
// - For beam weapons, provide a target array so the beam can stop on the nearest
//   hit (or max range if none).
FireResult tryFireWeapon(Ship& shooter,
                         WeaponType weapon,
                         double cooldownSimSec,
                         int distributorMk,
                         core::u64 shooterId,
                         bool fromPlayer,
                         const SphereTarget* beamTargets,
                         std::size_t beamTargetCount);

// Compute the heat delta for firing a weapon with a given distributor Mk.
double weaponHeatDelta(const WeaponDef& w, int distributorMk);

} // namespace stellar::sim
