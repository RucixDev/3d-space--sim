#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Ship.h"

namespace stellar::sim {

// -----------------------------------------------------------------------------
// ManeuverComputer — execute an "instantaneous" maneuver node as a continuous burn
// -----------------------------------------------------------------------------
//
// The trajectory preview system supports an instantaneous delta-v (ManeuverNode) for
// planning. This computer turns the planned delta-v into a continuous burn by:
//
//  - Pointing the ship's +Z (forward) axis along the desired delta-v.
//  - Starting the burn ~half the estimated burn duration early (classic centered burn).
//  - Throttling down on the final step to avoid overshooting the requested delta-v.
//
// Notes / limitations (by design for the prototype):
//  - Assumes thrust is applied along the ship's forward axis (6DOF translation is
//    available, but the guidance aligns the ship instead of using lateral thrust).
//  - Tracks achieved delta-v by projecting the ship's actual velocity change onto
//    the desired burn direction. Environmental accelerations (gravity) can therefore
//    perturb the estimate slightly.
//  - Designed to be deterministic + headless; safe to use in tests and tools.

enum class ManeuverComputerPhase : core::u8 {
  Off = 0,
  Orient = 1,   // rotating toward burn vector, waiting for start time
  Burn = 2,     // actively thrusting
  Complete = 3,
  Aborted = 4,
};

struct ManeuverPlan {
  // Absolute time in simulation days when the burn is centered.
  double nodeTimeDays{0.0};

  // Desired delta-v in world space (km/s).
  math::Vec3d deltaVWorldKmS{0, 0, 0};
};

struct ManeuverComputerParams {
  // Attitude control gain (maps orientation error -> normalized torque input).
  double faceGain{2.4};

  // Require alignment before burning (degrees).
  double alignToleranceDeg{4.0};

  // Extra lead time added to the centered-burn start time (seconds).
  // Positive values start earlier, negative values start later.
  double extraLeadTimeSec{0.0};

  // Delta-v completion threshold (km/s). Default: 2 m/s.
  double dvToleranceKmS{0.002};

  // If true, disables velocity dampers while burning (recommended).
  bool disableDampersDuringBurn{true};

  // If true, requests boost during the burn. (Gameplay may still clamp boost
  // availability based on capacitor state.)
  bool allowBoost{false};

  // If the burn is missed by more than this many seconds, abort (0 disables).
  double abortAfterMissedSec{45.0};
};

struct ManeuverComputerOutput {
  ShipInput input{};

  ManeuverComputerPhase phase{ManeuverComputerPhase::Off};

  // Diagnostics (seconds, km/s, degrees)
  double timeToNodeSec{0.0};
  double burnDurationSec{0.0};
  double burnStartLeadSec{0.0};
  double dvTotalKmS{0.0};
  double dvRemainingKmS{0.0};
  double alignmentErrorDeg{0.0};

  bool finished{false};
  bool aborted{false};
};

class ManeuverComputer {
public:
  void disengage();

  // Engage the computer with a fully specified plan.
  // The plan is copied, and internal progress state is reset.
  void engage(const Ship& ship, const ManeuverPlan& plan);

  bool active() const {
    return phase_ != ManeuverComputerPhase::Off &&
           phase_ != ManeuverComputerPhase::Complete &&
           phase_ != ManeuverComputerPhase::Aborted;
  }

  ManeuverComputerPhase phase() const { return phase_; }
  const ManeuverPlan& plan() const { return plan_; }

  // Step the maneuver computer and return recommended ship inputs.
  ManeuverComputerOutput update(const Ship& ship,
                                double nowTimeDays,
                                double dtSimSec,
                                const ManeuverComputerParams& params = {});

private:
  ManeuverPlan plan_{};
  ManeuverComputerPhase phase_{ManeuverComputerPhase::Off};

  math::Vec3d burnDirWorld_{0, 0, 1};
  double dvTotalKmS_{0.0};

  // Captured at burn start so we can estimate achieved delta-v from ship state.
  math::Vec3d velAtBurnStartKmS_{0, 0, 0};
};

} // namespace stellar::sim
