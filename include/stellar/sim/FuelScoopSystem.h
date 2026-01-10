#pragma once

#include "stellar/sim/Celestial.h"

namespace stellar::sim {

// -----------------------------------------------------------------------------
// FuelScoopSystem (headless)
// -----------------------------------------------------------------------------
//
// Deterministic helper used by the game layer for:
//  - computing star proximity heating (simple inverse-square falloff)
//  - computing fuel scooping rates when a fuel scoop module is deployed
//
// Notes:
//  - The caller is responsible for integrating rates over dtSim/dtReal.
//  - Fuel and heat are intentionally separate knobs so the game can keep
//    fuel progression (sim-time) and heat feedback (real-time) consistent with
//    the rest of the prototype.

struct FuelScoopParams {
  // Range window expressed in multiples of star radius:
  //  distR = distanceFromStarCenterKm / starRadiusKm
  double minRangeR{1.05};
  double maxRangeR{8.0};

  // Star heat is applied out to this range (beyond, it becomes 0).
  double starHeatRangeR{12.0};

  // Reference heat (heat points per real second) at distR == 2.
  double starHeatPerSecAt2R{12.0};

  // Hard clamp to keep the thermal system stable extremely close to the star.
  double starHeatMaxPerSec{120.0};

  // Maximum scoop fuel rate (fuel units per simulated second) for Mk1
  // at minRangeR around a 1.0 luminosity star.
  double scoopFuelPerSimSecMk1{0.012};

  // Additional heat from active scooping (heat points per real second) for Mk1
  // at minRangeR.
  double scoopHeatPerSecAtMaxMk1{6.0};

  // Mk multipliers.
  double scoopMk2Mult{1.6};
  double scoopMk3Mult{2.3};

  // Shapes the distance falloff (higher means "needs to be closer").
  double scoopFalloffPower{1.35};

  // Hard clamp on fuel rate to avoid explosive tuning issues.
  double scoopFuelMaxPerSimSec{0.06};
};

struct FuelScoopRates {
  // Whether the ship is within the nominal scooping window [minRangeR,maxRangeR].
  bool inRange{false};

  // True when an installed scoop is deployed and the ship is in-range.
  bool active{false};

  double distKm{0.0};
  double distR{0.0};

  // 0..1 factor used for UI (1 = very close / strong, 0 = far / none).
  double scoopFactor01{0.0};

  // Fuel gain rate (per simulated second). Only >0 when active.
  double fuelPerSimSec{0.0};

  // Star proximity heating (per real second). Can be >0 even without a scoop.
  double starHeatPerSec{0.0};

  // Additional heat from actively scooping (per real second). Only >0 when active.
  double scoopHeatPerSec{0.0};
};

FuelScoopRates computeFuelScoopRates(const Star& star,
                                     double distKmFromStarCenter,
                                     int fuelScoopMk,
                                     bool deployed,
                                     const FuelScoopParams& params = FuelScoopParams{});

} // namespace stellar::sim
