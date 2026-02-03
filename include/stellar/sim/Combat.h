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
  Decoy = 3,
};

// Missile seeker channel for simple countermeasure logic.
// Heat seekers are attracted to flares; Radar seekers are attracted to chaff/jammers.
enum class MissileSeekerType : core::u8 {
  Heat = 0,
  Radar = 1,
};

// Missile guidance mode.
//
// LeadPursuit:
//  - Simple "lead solve" each frame + turn-rate-limited steering.
//  - Works well for arcade feel and integrates nicely with decoy overrides.
//
// ProNav:
//  - Proportional Navigation (PN / Pro-Nav) guidance using a navigation constant.
//  - Tends to produce cleaner intercepts against maneuvering targets.
//
// NOTE: This is intentionally a small enum that is *not* persisted.
// Keep explicit values stable anyway for safety.
enum class MissileGuidance : core::u8 {
  LeadPursuit = 0,
  ProNav = 1,
};

struct SphereTarget {
  CombatTargetKind kind{CombatTargetKind::Ship};
  std::size_t index{0};
  core::u64 id{0};

  math::Vec3d centerKm{0, 0, 0};
  // Optional kinematics for guidance / lead solves.
  // Leave at zero if unknown.
  math::Vec3d velKmS{0, 0, 0};
  double radiusKm{1.0};

  // Optional aim-cone filter for soft aim assist.
  // If >= -0.5, a target is only considered when
  //   dot(rayDir, (center-origin).normalized()) >= minAimCos.
  // Use -1.0 to disable the filter.
  double minAimCos{-1.0};

  // Optional decoy signatures for missile seekers.
  // Only consulted when kind == CombatTargetKind::Decoy.
  double decoyHeat{0.0};
  double decoyRadar{0.0};
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

// Guided projectile (simple homing missile). Units are kilometers.
//
// The game treats missiles as a constant-speed projectile that can steer toward
// a target with a limited turn rate.
struct Missile {
  math::Vec3d prevKm{0, 0, 0};
  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};

  float r{1.0f}, g{1.0f}, b{1.0f};

  double ttlSimSec{0.0};
  double radiusKm{650.0};
  double dmg{0.0};

  // Explosion radius for splash damage.
  double blastRadiusKm{0.0};

  // Optional directional blast model.
  //
  // Some warheads focus fragments forward along the missile velocity vector at
  // detonation. When blastDirectionalStrength > 0, splash damage is multiplied
  // by a factor derived from the cosine between the missile's travel direction
  // at detonation and the direction from the detonation point to the target.
  //
  // - blastDirectionalStrength == 0: isotropic (legacy behavior)
  // - blastDirectionalStrength == 1: full directional weighting
  //
  // blastDirectionalMinFactor clamps the minimum damage multiplier for targets
  // directly behind the missile (cos = -1).
  double blastDirectionalStrength{0.0};   // 0..1
  double blastDirectionalMinFactor{0.25}; // 0..1

  // Optional blast occlusion by asteroids (line-of-sight).
  //
  // When enabled, asteroid spheres between the detonation point and a candidate
  // target attenuate splash damage. This allows ships to "hide" behind
  // rocks from nearby explosions.
  bool blastRequireLineOfSight{false};
  // Optional occlusion padding added to asteroid radii when testing blast LOS (km).
  double blastLineOfSightOcclusionPadKm{0.0};
  // Damage multiplier applied when occluded by an asteroid (0 => fully blocked).
  double blastOccludedFactor{0.0};

  // Optional proximity fuse standoff distance (km).
  //
  // When > 0, a missile may detonate without a direct collision when it comes within
  //   (target.radiusKm + missile.radiusKm + proximityFuseKm)
  // of a target. Detonation time is chosen near closest approach during the step.
  double proximityFuseKm{0.0};

  // Max steering rate (rad/s) toward the desired aim direction.
  double turnRateRadS{0.0};

  // Optional missile motor model (boost/coast).
  //
  // By default missiles are constant-speed projectiles. If you populate these
  // fields, the missile can accelerate along its current velocity direction
  // while its motor is burning.
  //
  // Acceleration only occurs when:
  //   - thrustAccelKmS2 > 0
  //   - motorBurnRemainingSimSec > 0
  //   - maxSpeedKmS > currentSpeedKmS
  //
  // These parameters are *not* persisted; they are a lightweight sim-side
  // tuning tool for gameplay feel.
  double thrustAccelKmS2{0.0};
  double maxSpeedKmS{0.0};
  double motorBurnRemainingSimSec{0.0};

  // Optional lateral acceleration limit (km/s^2).
  //
  // When > 0, this caps turning on top of turnRateRadS by ensuring:
  //   speedKmS * effectiveTurnRateRadS <= maxLateralAccelKmS2
  //
  // This is useful when missiles can accelerate: it preserves a consistent
  // "G-limit" instead of allowing higher speeds to implicitly create more
  // turning authority.
  double maxLateralAccelKmS2{0.0};

  // Optional turn acceleration limit (rad/s^2).
  //
  // When > 0, missiles do not instantly achieve their commanded turn rate. Instead,
  // they maintain an internal angular-velocity state that can change by at most this
  // rate per second (a lightweight "autopilot lag" model).
  //
  // This makes last-second maneuvers more meaningful at close range and helps avoid
  // fully instantaneous heading snaps when using high turn rates.
  //
  // 0.0 disables turn acceleration limiting (default).
  double maxTurnAccelRadS2{0.0};

  // Internal: current turn angular velocity in world space (rad/s).
  // Not serialized/persisted. Used only when maxTurnAccelRadS2 > 0.
  bool hasTurnOmega{false};
  math::Vec3d turnOmegaWorld{0, 0, 0};

  // Optional seeker activation delay (sim seconds).
  //
  // When > 0, the missile behaves like it has a midcourse "datalink" phase:
  //   - it continues to guide toward its locked target even if the target is
  //     outside the seeker's field-of-view
  //   - it ignores countermeasure / decoy attraction until the seeker activates
  //
  // This produces a simple two-phase behavior for radar-style missiles:
  // midcourse guidance -> terminal active seeker.
  //
  // 0.0 means the seeker is active immediately (default).
  double seekerActivationSimSec{0.0};

  // Optional midcourse datalink maximum range (km).
  //
  // When > 0 and the seeker is not yet active (seekerActivationSimSec > 0), the missile
  // only receives target updates if the launching platform (shooterId) is within this
  // range of the locked target (and passes optional LOS gating).
  //
  // When 0, midcourse guidance behaves like a perfect datalink (legacy behavior).
  double datalinkRangeKm{0.0};

  // Optional: line-of-sight requirement for midcourse datalink updates.
  //
  // When true and datalinkRangeKm > 0, target updates before seeker activation require
  // that the shooter has a clear line-of-sight to the target (not occluded by asteroids).
  bool datalinkRequireLineOfSight{false};
  // Optional occlusion padding added to asteroid radii when testing datalink LOS (km).
  double datalinkOcclusionPadKm{0.0};

  // Optional autonomous target acquisition range (km).
  //
  // When > 0 and the seeker is active, if the missile does not currently have a
  // trackable target (and target memory has expired), it may acquire a new
  // Ship/Player target within this range and within the seeker FOV.
  //
  // 0.0 disables reacquisition (default).
  double autoAcquireRangeKm{0.0};

  // Seeker tuning (optional): used for target tracking + countermeasure / decoy attraction.
  MissileSeekerType seeker{MissileSeekerType::Heat};
  // Seeker field-of-view cosine (dot(fwdDir, toDir)).
  // 1.0 = only straight ahead, 0.0 = 90 degrees, -1.0 = any direction.
  double seekerFovCos{0.0};

  // Optional: seeker gimbal slew rate (rad/s).
  //
  // When > 0, the missile maintains an internal seeker line-of-sight direction
  // that slews toward its current tracked guide at this maximum rate.
  //
  // This primarily exists to make seeker lock quality more dynamic:
  // high line-of-sight rates (close-range merges, sharp target maneuvers, or
  // decoy pulls) can temporarily reduce track quality and make countermeasures
  // and evasive flying more effective.
  //
  // 0.0 disables slew limiting (instant pointing).
  double seekerSlewRateRadS{0.0};

  // Internal: current seeker pointing direction in world space.
  // Not serialized/persisted.
  bool hasSeekerDir{false};
  math::Vec3d seekerDirWorld{0, 0, 1};
  // How resistant the missile is to decoys. A decoy must exceed
  //   targetScore * decoyResistance
  // to override guidance.
  double decoyResistance{1.0};

  // Optional: line-of-sight requirement for seeker tracking.
  //
  // When true, once the seeker is active the missile can only track/compare
  // targets and decoys that are not occluded by asteroids (sphere targets of
  // kind Asteroid) between the missile and the candidate.
  //
  // This creates a simple, deterministic "terrain masking" mechanic:
  // breaking LOS behind an asteroid can temporarily break lock (subject to
  // targetMemorySimSec).
  bool requireLineOfSight{false};
  // Optional occlusion padding added to asteroid radii when testing LOS (km).
  double lineOfSightOcclusionPadKm{0.0};

  // Optional: radar doppler notch (km/s).
  //
  // When seeker == Radar and this value is > 0, the missile considers a target
  // untrackable while the absolute radial velocity (along the line-of-sight)
  // is below this threshold. This is a deterministic approximation of doppler
  // notch / "beaming" behavior.
  double radarDopplerNotchKmS{0.0};

  // Optional decoy discrimination gates (seeker-active only).
  //
  // These parameters help seekers (especially radar) reject decoys that are
  // clearly not co-located with the true target in angle / range / radial
  // velocity. They are intentionally simple and deterministic:
  //
  //  - decoyAngleGateCos: decoy must satisfy
  //        dot(toTargetDir, toDecoyDir) >= decoyAngleGateCos
  //    when this value is >= -0.5 (disabled by default).
  //
  //  - decoyRangeGateKm: decoy score is weighted by
  //        clamp01(1 - |rangeDecoy - rangeTarget| / decoyRangeGateKm)
  //    when > 0 (disabled by default).
  //
  //  - decoyDopplerGateKmS: decoy score is weighted by
  //        clamp01(1 - |vrDecoy - vrTarget| / decoyDopplerGateKmS)
  //    when > 0 (disabled by default).
  //
  // Together these approximate range/angle/velocity track gates and make
  // countermeasures feel more believable without adding non-deterministic
  // randomness.
  double decoyAngleGateCos{-1.0};
  double decoyRangeGateKm{0.0};
  double decoyDopplerGateKmS{0.0};

  // Optional: close-range decoy "burn-through" (radar seekers).
  //
  // In practice, many radar countermeasures (noise jamming, some deceptive
  // techniques) become less effective at short range as the seeker has higher
  // SNR / better discrimination. This is often discussed as "burn-through".
  //
  // When enabled (decoyBurnThroughRangeKm > 0) and seeker == Radar, decoy
  // attraction is reduced as the missile approaches the current target track.
  // The decoy score is multiplied by:
  //
  //   lerp(decoyBurnThroughMinFactor, 1,
  //        clamp01(targetRangeKm / decoyBurnThroughRangeKm))
  //
  // 0 disables the feature.
  double decoyBurnThroughRangeKm{0.0};

  // Lower bound for the burn-through multiplier at zero range.
  // Typical values: 0.02..0.2 (clamped to [0,1]).
  double decoyBurnThroughMinFactor{0.05};

  // -------------------------------------------------------------------------
  // Optional: asteroid avoidance (deterministic)
  // -------------------------------------------------------------------------
  //
  // When enabled (asteroidAvoidanceStrength > 0), the missile performs a
  // forward-looking collision check against asteroid spheres and biases its
  // commanded heading away from imminent intersections. This is a small, fully
  // deterministic "potential-field" style steering nudge — not a full path
  // planner — intended to keep missiles from suiciding into dense belts.
  //
  // Strength is a unitless blend factor (0 disables).
  double asteroidAvoidanceStrength{0.0};

  // Lookahead horizon used for the collision predictor (sim seconds).
  // If <= 0, a conservative default (~1s) is used.
  double asteroidAvoidanceLookaheadSimSec{0.0};

  // Extra padding (km) added to asteroid radius for the predictor.
  double asteroidAvoidancePadKm{0.0};


  // -------------------------------------------------------------------------
  // Optional: missile swarm / salvo coordination (deterministic)
  // -------------------------------------------------------------------------
  //
  // When enabled, missiles can apply a small steering bias based on nearby
  // friendly missiles that are pursuing the same target. This produces a
  // light-weight "swarm" effect that:
  //   - reduces missile stacking (multiple missiles occupying the same path)
  //   - enables simple multi-axis attacks without an explicit planner
  //
  // The behavior is a deterministic variant of classic flocking rules:
  // separation / cohesion / alignment, applied as a bounded heading bias.
  //
  // All strengths are unitless blend weights (0 disables each term).
  double swarmSeparationStrength{0.0};
  double swarmCohesionStrength{0.0};
  double swarmAlignmentStrength{0.0};

  // Desired minimum spacing between missiles (km). When <= 0, separation is
  // treated as disabled even if swarmSeparationStrength > 0.
  double swarmSeparationKm{0.0};

  // Neighbor consideration radius (km). When <= 0, a conservative default is
  // derived from swarmSeparationKm.
  double swarmNeighborRangeKm{0.0};

  // Maximum steering angle applied by swarm bias (radians). When <= 0, the
  // bias is limited to a fraction of the missile's per-step turn capability.
  double swarmMaxSteerRad{0.0};



  // Optional: decoy commitment duration (sim seconds).
  //
  // When > 0, if a decoy overrides the true target, the missile will remain
  // committed to that decoy for at least this duration before it can switch
  // back. This prevents rapid "thrashing" when scores are close and makes
  // countermeasures feel more meaningful.
  double decoyCommitSimSec{0.0};

  // Internal: current committed decoy target (not serialized/persisted).
  core::u64 committedDecoyId{0};
  double decoyCommitRemainingSimSec{0.0};

  // Optional: target kinematic memory (sim seconds).
  //
  // When > 0, if the missile temporarily cannot track its locked target (e.g.
  // due to seeker field-of-view limits), it will continue guiding toward the
  // last known target position/velocity for up to this duration.
  //
  // This is a lightweight, deterministic approximation of inertial midcourse
  // guidance and makes "lock breaks" less binary while still allowing evasive
  // maneuvers to matter.
  double targetMemorySimSec{0.0};

  // Internal: last known target state (not serialized/persisted).
  bool hasLastKnownTarget{false};
  math::Vec3d lastKnownTargetPosKm{0, 0, 0};
  math::Vec3d lastKnownTargetVelKmS{0, 0, 0};
  double lastKnownTargetAgeSimSec{0.0};

  // Internal: missile age (sim seconds). Used for seeker activation timing.
  double ageSimSec{0.0};


  // Optional: seeker track quality (deterministic).
  //
  // When enabled and the seeker is active, the missile maintains a [0,1] scalar:
  //   - rises toward 1 while the locked target is directly trackable
  //   - decays toward 0 while guiding on memory / without fresh measurements
  //
  // Track quality modulates effective decoy resistance to create believable counterplay:
  // losing lock temporarily makes countermeasures more likely to work.
  //
  // effectiveResist = decoyResistance * lerp(trackQualityResistFloor, 1.0, trackQuality)
  bool enableTrackQuality{false};
  double trackQuality{1.0};
  // Half-life to rise from current value toward 1 when directly tracking.
  double trackQualityRiseHalfLifeSimSec{0.15};
  // Half-life to decay toward 0 when not directly tracking.
  double trackQualityFallHalfLifeSimSec{0.25};
  // Minimum fraction of decoyResistance when trackQuality==0. Clamped to [0,1].
  double trackQualityResistFloor{0.25};


  // Optional: augmented proportional navigation (APN) target acceleration feed-forward.
  //
  // When guidance == ProNav and apnTargetAccelGain > 0, the missile adds a feed-forward
  // term proportional to the estimated target acceleration component that is
  // perpendicular to the current line-of-sight (and projected to be normal to the
  // missile velocity for energy-conserving steering).
  //
  // A common choice is apnTargetAccelGain = navConstant * 0.5 (classical APN scaling).
  double apnTargetAccelGain{0.0};
  // Optional magnitude clamp for the estimated target acceleration used for APN (km/s^2).
  // 0 disables clamping.
  double apnMaxTargetAccelKmS2{0.0};
  // Optional half-life (sim seconds) for filtering the APN target acceleration estimate.
  // 0 disables filtering (raw finite-difference estimate).
  double apnAccelHalfLifeSimSec{0.0};

  // Internal: APN target acceleration estimation state (not serialized/persisted).
  bool hasTargetVelSample{false};
  math::Vec3d targetVelSampleKmS{0, 0, 0};
  double targetVelSampleAgeSimSec{0.0};
  math::Vec3d targetAccelEstKmS2{0, 0, 0};


  // Guidance mode.
  MissileGuidance guidance{MissileGuidance::LeadPursuit};
  // Navigation constant for ProNav guidance (typical real-world values: ~3-5).
  // Only used when guidance == ProNav.
  double navConstant{3.5};

  bool fromPlayer{false};
  core::u64 shooterId{0};

  // Target lock.
  bool hasTarget{false};
  CombatTargetKind targetKind{CombatTargetKind::Ship};
  core::u64 targetId{0};
};

struct MissileDetonation {
  math::Vec3d pointKm{0, 0, 0};
  double blastRadiusKm{0.0};
  double baseDmg{0.0};
  bool fromPlayer{false};
  core::u64 shooterId{0};
};

struct MissileHit {
  CombatTargetKind kind{CombatTargetKind::Ship};
  std::size_t targetIndex{0};
  core::u64 targetId{0};
  math::Vec3d pointKm{0, 0, 0};
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

// Move missiles forward by dtSim seconds, steer toward locked targets, and emit
// detonation events + splash hit events.
//
// Missiles are consumed when they detonate.
void stepMissiles(std::vector<Missile>& missiles,
                  double dtSim,
                  const SphereTarget* targets,
                  std::size_t targetCount,
                  std::vector<MissileDetonation>& outDetonations,
                  std::vector<MissileHit>& outHits);

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

  bool hasMissile{false};
  Missile missile{};

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
