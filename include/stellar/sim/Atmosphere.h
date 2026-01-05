#pragma once

#include "stellar/math/Vec3.h"
#include "stellar/sim/System.h"

#include <cstddef>
#include <string>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Atmosphere (deterministic, headless)
// -----------------------------------------------------------------------------
//
// This module provides a lightweight, gameplay-friendly atmosphere model:
//  - Planet-specific parameters derived deterministically from PlanetType and
//    rough mass/radius.
//  - Exponential density falloff: rho(h) = rho0 * exp(-h / H)
//  - Simple drag based on dynamic pressure (q = 0.5 * rho * v^2).
//
// The goal is *feel* and determinism rather than high-fidelity aerodynamics.

struct AtmosphereModel {
  bool hasAtmosphere{false};
  double seaLevelDensityKgM3{0.0};
  double scaleHeightKm{0.0};
  double topAltitudeKm{0.0}; // hard cutoff (km above surface)
};

// Deterministic heuristic atmosphere model from PlanetType + mass/radius.
AtmosphereModel atmosphereModelForPlanet(const Planet& p);

// Density at altitude above surface (km). Returns 0 if outside the atmosphere.
double atmosphereDensityKgM3(const AtmosphereModel& model, double altitudeKm);

struct AtmosphereParams {
  bool includePlanets{true};

  // Gameplay tuning: multiply computed density by this factor.
  double densityScale{1.0};

  // Drag: a = q * Cd * A / m
  double dragCd{0.9};
  double referenceAreaM2{80.0}; // effective cross-sectional area

  // Clamp magnitude of drag decel (km/s^2). 0 disables clamping.
  double maxDecelKmS2{0.75};

  // Heating model: heatRate = (q[kPa] * heatPerKPaSec) * heatingScale
  bool applyHeating{true};
  double heatPerKPaSec{0.06};
  double heatingScale{1.0};

  // Optional guard: if > 0, clamp relative speed used for q calculation.
  // (Useful as a safety net if the ship somehow enters atmosphere at extreme
  // speeds due to timewarp/supercruise edge cases.)
  double maxRelSpeedKmS{0.0};
};

struct AtmosphereSample {
  bool inAtmosphere{false};

  // Dominant planet (index in StarSystem::planets).
  std::size_t planetIndex{0};
  std::string planetName{};

  double planetRadiusKm{0.0};
  double altitudeKm{0.0};

  // Local properties.
  double densityKgM3{0.0};
  double dynamicPressurePa{0.0};

  // Relative velocity to the atmosphere frame (planet translational velocity).
  math::Vec3d relVelKmS{0, 0, 0};
  double relSpeedKmS{0.0};

  // Outputs.
  math::Vec3d dragAccelKmS2{0, 0, 0};
  double heatingHeatPerSec{0.0};
};

// Compute atmosphere drag (and optional heating rate) at the ship state.
//
// The returned sample is for the single most relevant planet atmosphere
// (highest dynamic pressure). If no atmosphere applies, inAtmosphere=false.
AtmosphereSample sampleSystemAtmosphere(const StarSystem& sys,
                                        double timeDays,
                                        const math::Vec3d& shipPosKm,
                                        const math::Vec3d& shipVelKmS,
                                        double shipMassKg,
                                        const AtmosphereParams& params = {});

} // namespace stellar::sim
