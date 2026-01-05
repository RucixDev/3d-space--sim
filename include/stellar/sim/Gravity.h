#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/System.h"
#include "stellar/sim/Units.h"

#include <algorithm>
#include <string>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Gravity (deterministic, headless)
// -----------------------------------------------------------------------------
//
// The core sim uses:
//   - km / km/s / km/s^2 for moment-to-moment ship dynamics
//   - AU / days for orbital elements (see Orbit.h + Units.h)
//
// This module adds a *small* amount of Newtonian gravity so local-space flight
// can optionally feel more "space-y" near planets without forcing a full
// high-fidelity n-body integrator.
//
// Design goals:
//   - deterministic + headless (usable in tests/tools)
//   - cheap enough to call per-frame (star + a few planets)
//   - safe: includes softening/clamping to avoid singular accelerations

// Physical constants (converted to km/kg/s units).
constexpr double kGravG_Km3KgS2 = 6.67430e-20; // 6.67430e-11 m^3/(kg*s^2) * 1e-9

constexpr double kMassSunKg = 1.98847e30;
constexpr double kMassEarthKg = 5.9722e24;

constexpr double kRadiusSunKm = 695700.0;
constexpr double kRadiusEarthKm = 6371.0;

constexpr double kStandardGravityMS2 = 9.80665;

inline double muFromMassKg(double massKg) {
  return kGravG_Km3KgS2 * massKg;
}

inline double muFromSolarMass(double massSol) {
  return muFromMassKg(massSol * kMassSunKg);
}

inline double muFromEarthMass(double massEarth) {
  return muFromMassKg(massEarth * kMassEarthKg);
}

inline double radiusKmFromSolarRadius(double radiusSol) {
  return radiusSol * kRadiusSunKm;
}

inline double radiusKmFromEarthRadius(double radiusEarth) {
  return radiusEarth * kRadiusEarthKm;
}

struct GravityParams {
  bool includeStar{true};
  bool includePlanets{true};

  // Multiply the resulting acceleration vector by this scalar.
  double scale{1.0};

  // Softening radius is bodyRadiusKm * softeningRadiusScale.
  // This prevents singular accelerations *inside* body volumes.
  double softeningRadiusScale{1.05};

  // If > 0, clamp magnitude of resulting acceleration (km/s^2).
  double maxAccelKmS2{0.0};
};

struct GravityBody {
  enum class Kind : core::u8 { Star = 0, Planet = 1 };

  Kind kind{Kind::Planet};

  // Stable-ish identifier:
  //  - star: systemId
  //  - planet: planet index within StarSystem::planets
  core::u64 id{0};

  // Debug label for UI/tools.
  std::string name{};

  // Inertial frame (local system) at timeDays.
  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};

  // Gravitational parameter (mu = G*M), in km^3/s^2.
  double muKm3S2{0.0};

  // Physical radius (km) for softening/diagnostics.
  double radiusKm{0.0};
};

double muStarKm3S2(const Star& star);
double muPlanetKm3S2(const Planet& planet);

double radiusStarKm(const Star& star);
double radiusPlanetKm(const Planet& planet);

// Fill `out` with gravitational bodies for the system at `timeDays`.
// `out` is cleared.
void gatherGravityBodies(const StarSystem& sys,
                         double timeDays,
                         std::vector<GravityBody>& out,
                         const GravityParams& params = {});

// Acceleration from a single body at a query position.
//
// minRadiusKm is used as a softening distance to avoid singularities.
math::Vec3d gravityAccelFromBodyKmS2(const math::Vec3d& bodyPosKm,
                                    double muKm3S2,
                                    const math::Vec3d& queryPosKm,
                                    double minRadiusKm);

// Total gravitational acceleration at `posKm` due to system star and planets.
math::Vec3d systemGravityAccelKmS2(const StarSystem& sys,
                                  double timeDays,
                                  const math::Vec3d& posKm,
                                  const GravityParams& params = {});

struct DominantGravityBody {
  bool valid{false};
  GravityBody body{};

  // Accel contribution from this body (after params.scale, before global clamp).
  math::Vec3d accelKmS2{0, 0, 0};
  double accelMagKmS2{0.0};
};

// Which body dominates gravity at `posKm` (largest |a| contribution).
// Returns {valid=false} if no bodies were included.
DominantGravityBody dominantGravityBody(const StarSystem& sys,
                                       double timeDays,
                                       const math::Vec3d& posKm,
                                       const GravityParams& params = {});

} // namespace stellar::sim
