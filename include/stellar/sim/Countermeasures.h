#pragma once

#include "stellar/core/Random.h"
#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Combat.h"

#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Countermeasures (flares / chaff) — headless
// -----------------------------------------------------------------------------
//
// A small, deterministic countermeasure model intended to plug into the existing
// missile guidance system in Combat.
//
// The general idea is:
//  - countermeasures are spawned as short-lived, drifting "decoy" spheres
//  - they expose decoyHeat / decoyRadar signatures via SphereTarget fields
//  - stepMissiles() can (optionally) bias toward these decoys based on seeker type
//
// This keeps the core sim deterministic and rendering/UI-independent.

enum class CountermeasureType : core::u8 {
  Flare = 0, // attracts heat seekers
  Chaff = 1, // attracts radar seekers
};

struct Countermeasure {
  core::u64 id{0};
  CountermeasureType type{CountermeasureType::Flare};

  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};

  double radiusKm{0.15};

  // Remaining time-to-live (sim seconds).
  double ttlSimSec{0.0};
  double ttlMaxSimSec{0.0};

  // Seeker attraction strengths (arbitrary units).
  // These are mapped into SphereTarget::decoyHeat/decoyRadar and can decay with TTL.
  double heatStrength{0.0};
  double radarStrength{0.0};
};

struct CountermeasureBurstParams {
  int count{8};

  // How long each decoy persists in sim seconds.
  double ttlSimSec{6.0};

  // Size of each decoy collision sphere (km).
  double radiusKm{0.15};

  // Ejection speed relative to the ship (km/s).
  double ejectSpeedKmS{0.10};

  // Scatter amount in local right/up directions (unitless).
  // 0 -> all decoys eject exactly opposite forward
  // ~0.35 -> moderately wide spread cone
  double spread{0.35};

  // Base strengths (arbitrary units).
  double heatStrength{8.0};
  double radarStrength{8.0};
};

// Spawn a burst of countermeasures from a ship-like frame.
//
// - ioNextId must be maintained by the caller for stable unique ids.
// - shipRight/up/forward should be unit-ish vectors (they will be normalized).
void spawnCountermeasureBurst(std::vector<Countermeasure>& ioCountermeasures,
                              core::u64& ioNextId,
                              CountermeasureType type,
                              const math::Vec3d& originKm,
                              const math::Vec3d& baseVelKmS,
                              const math::Vec3d& shipRight,
                              const math::Vec3d& shipUp,
                              const math::Vec3d& shipForward,
                              const CountermeasureBurstParams& params = {});

// Integrate countermeasures forward in time and remove expired items.
void stepCountermeasures(std::vector<Countermeasure>& ioCountermeasures, double dtSim);

// Append countermeasures as Combat SphereTargets so missiles can steer to them and/or collide.
// Caller controls the lifetime and ordering of ioTargets.
void appendCountermeasureTargets(const std::vector<Countermeasure>& countermeasures,
                                 std::vector<SphereTarget>& ioTargets);

} // namespace stellar::sim
