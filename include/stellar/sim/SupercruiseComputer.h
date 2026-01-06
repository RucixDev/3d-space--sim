#pragma once

#include "stellar/math/Vec3.h"
#include "stellar/sim/Ship.h"

namespace stellar::sim {

// Core guidance logic for "supercruise" (fast in-system travel).
//
// The game/app can wrap this with its own state machine (charging/cooldown,
// interdictions, encounters, etc.). This module focuses on:
//  - Computing a safe-drop window (time-to-arrival heuristic).
//  - Producing translational + rotational ShipInput to approach a destination.
//  - Providing a braking-distance-aware speed limit so nav-assist can actually
//    reach the safe drop window under finite acceleration.

struct SupercruiseParams {
  // "6/7 second rule" target time-to-arrival used by the safe-drop heuristic.
  double safeTtaSec{7.0};
  // Acceptable +/- slack around safeTtaSec.
  double safeWindowSlackSec{2.0};
  // Minimum positive closing speed required to consider the window valid.
  double minClosingKmS{0.05};

  // Additional safety gates for the drop window.
  //
  // The classic "6/7 second rule" is mostly about along-track closing speed, but
  // in practice you also want to avoid dropping while sliding past the target.
  //
  // Lateral velocity is measured orthogonal to the ship→destination line.
  // A value of 0.6 corresponds to ~31 degrees off-axis (tan(angle)).
  double maxLateralFrac{0.6};

  // Supercruise handling caps (used by callers to set ship caps while active).
  double accelCapKmS2{6.0};
  double angularCapRadS2{1.2};

  // Speed model.
  double maxSpeedKmS{18000.0};
  double assistSpeedMinKmS{60.0};
  double manualSpeedMinKmS{90.0};
  double manualSpeedDistGain{0.0008};

  // Guidance tuning.
  double velTimeConstantSec{6.0}; // smaller = more aggressive
  double faceGain{1.6};           // yaw/pitch gain

  // If true, nav-assist uses a braking-distance speed limit so it can
  // realistically decelerate to the safe-drop speed at the drop radius.
  bool useBrakingDistanceLimit{true};
};

struct SupercruiseHud {
  bool safeDropReady{false};
  double distKm{0.0};
  double closingKmS{0.0};
  double ttaSec{0.0};

  // Lateral relative speed (orthogonal to the approach line) and predicted miss
  // distance if we kept current relative velocity.
  double lateralKmS{0.0};
  double missKm{0.0};

  // Diagnostics for tuning / UI.
  double desiredSpeedKmS{0.0};
  double speedLimitKmS{0.0};
};

struct SupercruiseGuidanceResult {
  ShipInput input{};
  SupercruiseHud hud{};

  double recommendedMaxLinearAccelKmS2{6.0};
  double recommendedMaxAngularAccelRadS2{1.2};

  // If true, the caller should drop to normal space this frame.
  bool dropNow{false};
  // If dropNow is true, indicates whether it's an emergency drop.
  bool emergencyDrop{false};

  // Direction from the ship to destination (world-space, normalized when valid).
  math::Vec3d dirToDest{0, 0, 1};
};

// Computes supercruise guidance toward destPosKm/destVelKmS.
//
// - navAssistEnabled: if true, the controller enforces a safe approach profile
//   and will request an auto-drop when safeDropReady.
// - dropRequested: user requested a manual drop *this frame*.
// - interdicted: if true, rotational input will not be overridden (caller can
//   let the player steer), but translational guidance is still provided.
SupercruiseGuidanceResult guideSupercruise(const Ship& ship,
                                          const math::Vec3d& destPosKm,
                                          const math::Vec3d& destVelKmS,
                                          double dropRadiusKm,
                                          bool navAssistEnabled,
                                          bool dropRequested,
                                          bool interdicted,
                                          const SupercruiseParams& params = {});

} // namespace stellar::sim
