#pragma once

#include "stellar/core/Types.h"

#include <algorithm>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// ThermalSystem
// -----------------------------------------------------------------------------
// A small, deterministic ship thermal model.
//
// Notes:
//  - This intentionally keeps the original prototype "heat" feel: linear heating/cooling
//    in real time, with a simple overheat consequence.
//  - It is extracted into the core library so both the game and headless tools/tests
//    can share the logic.
//
// Units:
//  - heat is an abstract scalar in [heatMin..heatMax]
//  - rates are "heat units per real second"

enum class ThermalFsdState : core::u8 {
  Idle = 0,
  Charging = 1,
  Jumping = 2,
};

struct ThermalParams {
  // Heat clamp range.
  double heatMin{0.0};
  double heatMax{120.0};

  // Base cooling rates (for a ship with heatCoolRate == referenceCoolRate).
  double baseCoolUndocked{10.0};
  double baseCoolDocked{20.0};
  double referenceCoolRate{10.0};

  // Continuous heating sources (when undocked).
  double heatPerBoostSec{18.0};
  double heatPerSupercruiseSec{6.0};
  double heatPerFsdChargeSec{12.0};
  double heatPerFsdJumpSec{16.0};

  // Overheat behavior.
  double overheatStart{100.0};

  // Hull damage model (matches the original prototype):
  //   hullDamage = (heat - overheatStart) * overheatHullDamageRate * hullMax * dtReal
  double overheatHullDamageRate{0.0008};
};

struct ThermalInputs {
  // Real-time timestep in seconds.
  double dtReal{0.0};

  bool docked{false};
  bool supercruiseActive{false};
  ThermalFsdState fsd{ThermalFsdState::Idle};

  // Fraction of this frame where boost was applied (0..1).
  double boostAppliedFrac{0.0};

  // Instantaneous heat added this frame (weapon shots, emergency drops, etc).
  // Applied before cooling/heating integration.
  double heatImpulse{0.0};

  // Continuous external heat input (heat units per real second).
  // Used for star proximity heating, fuel scooping, etc.
  double externalHeatPerSec{0.0};

  // Ship cooling stat (see ShipLoadout::ShipDerivedStats::heatCoolRate).
  // Values around 10..15 for early-game ships.
  double heatCoolRate{10.0};

  // Hull max, used to scale overheat damage. If <= 0, damage is suppressed.
  double hullMax{0.0};
};

struct ThermalStepResult {
  double heat{0.0};
  double hullDamage{0.0};
  bool overheated{false};

  // Optional debug values.
  double heatInRate{0.0};
  double coolRate{0.0};
};

// Advance the thermal model by one real-time step.
//
// The model is intentionally linear:
//   heat += (heatInRate - coolRate) * dtReal
//
// where:
//   heatInRate is determined by inputs (boost/supercruise/fsd)
//   coolRate scales from baseCool*heatCoolRate/referenceCoolRate
//
// Returns new heat (clamped) and any hull damage due to overheating.
ThermalStepResult stepThermal(double currentHeat,
                              const ThermalInputs& in,
                              const ThermalParams& params = {});

} // namespace stellar::sim
