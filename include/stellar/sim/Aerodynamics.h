#pragma once

#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Atmosphere.h"

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Aerodynamics (deterministic, gameplay-friendly)
// -----------------------------------------------------------------------------
//
// Atmosphere.cpp provides:
//   - exponential density model (rho)
//   - dynamic pressure (q)
//   - a simple body-drag acceleration
//
// This module layers a *lift + stability* model on top of AtmosphereSample.
//
// Intent:
//   - Make close-to-planet flight feel like "flight" rather than "draggy space".
//   - Keep it deterministic + cheap enough for per-frame use.
//   - Be tunable from UI and safe under extreme time scales.
//
// Notes:
//   - This is not a CFD solver; it's a compact, controllable approximation.
//   - Output is expressed as accelerations (km/s^2, rad/s^2) so it can be fed
//     directly into Ship integration.

struct AerodynamicsParams {
  bool enabled{false};

  // Ignore ultra-low q to prevent jitter near vacuum.
  double minDynamicPressurePa{25.0};

  // Effective lifting area (m^2). (Think: wings + body lift.)
  double wingAreaM2{120.0};

  // Thin-airfoil-ish slope (dCL/dalpha) in 1/rad.
  double liftSlopePerRad{5.0};

  // Maximum absolute CL (pre-stall clamp).
  double clMax{1.25};

  // Stall onset (degrees). Past this, lift rolls off and drag rises.
  double stallAngleDeg{18.0};

  // Induced drag factor: CDi = k * CL^2.
  double inducedDragK{0.06};

  // Extra drag added in deep stall (dimensionless Cd).
  double stallDragCd{1.10};

  // Lift axis in ship-local coordinates (defaults to +Y).
  math::Vec3d liftAxisLocal{0, 1, 0};

  // Safety clamps.
  double maxLiftAccelKmS2{0.75};      // 0 = no clamp
  double maxExtraDragAccelKmS2{0.75}; // 0 = no clamp
  double maxAngAccelRadS2{6.0};       // 0 = no clamp

  // Control/stability is scaled by q / qRef (clamped).
  double qRefPa{20000.0};             // 20 kPa (~200 mbar)
  double qScaleMax{6.0};              // max q/qRef

  // Optional: treat pilot torque input as "control surfaces" in atmosphere.
  bool controlSurfaces{true};
  // Angular acceleration authority (rad/s^2) at qRef.
  math::Vec3d controlAngAccelRadS2{2.0, 1.2, 2.6}; // pitch,yaw,roll

  // Aerodynamic rotational damping (per-second) at qRef.
  // Produces angAccel = -omega * dampingPerS * (q/qRef).
  math::Vec3d dampingPerS{1.6, 1.2, 2.0}; // pitch,yaw,roll

  // "Weathervane" alignment toward the velocity vector (rad/s^2) at qRef.
  // Helps aircraft-like feel: the nose naturally wants to follow the airflow.
  bool alignToVelocity{true};
  double alignGain{3.0}; // multiplies alpha/beta (rad)
};

struct AerodynamicsSample {
  bool active{false};

  double dynamicPressurePa{0.0};

  // Angle of attack (alpha) and sideslip (beta), radians.
  double alphaRad{0.0};
  double betaRad{0.0};

  // Lift coefficient (after stall shaping).
  double cl{0.0};

  // Extra drag coefficient contributed by induced + stall drag.
  double cdExtra{0.0};

  // 0..1 stall blend factor (0 = no stall, 1 = deep stall).
  double stall01{0.0};

  // Linear accelerations (world space).
  math::Vec3d liftAccelKmS2{0, 0, 0};
  math::Vec3d extraDragAccelKmS2{0, 0, 0};

  // Angular acceleration (ship-local).
  math::Vec3d angAccelLocalRadS2{0, 0, 0};
};

// Compute aerodynamic lift/stability for a ship in a sampled atmosphere.
//
// Inputs:
//   - atmo: from sampleSystemAtmosphere (contains q, relVel)
//   - shipOrient: ship attitude (local -> world)
//   - shipAngVelLocalRadS: ship angular velocity in local axes (rad/s)
//   - shipMassKg: mass used for forces (include cargo if desired)
//   - controlInputLocal: optional pilot control vector in [-1,1] per axis
//                        (pitch,yaw,roll), used when controlSurfaces=true
//
// Output:
//   - AerodynamicsSample with accelerations in km/s^2 and rad/s^2.
AerodynamicsSample computeAerodynamics(const AtmosphereSample& atmo,
                                      const math::Quatd& shipOrient,
                                      const math::Vec3d& shipAngVelLocalRadS,
                                      double shipMassKg,
                                      const AerodynamicsParams& params,
                                      math::Vec3d controlInputLocal = {0, 0, 0});

} // namespace stellar::sim
