#include "stellar/sim/ThermalSystem.h"

#include <cmath>

namespace stellar::sim {

ThermalStepResult stepThermal(double currentHeat,
                              const ThermalInputs& in,
                              const ThermalParams& params) {
  ThermalStepResult out{};

  const double dt = std::max(0.0, in.dtReal);

  double heat = currentHeat;
  if (!std::isfinite(heat)) heat = 0.0;

  // Apply instantaneous heat sources first.
  if (std::isfinite(in.heatImpulse)) {
    heat += in.heatImpulse;
  }

  // Compute continuous heat-in rate (per real second).
  double heatIn = 0.0;
  if (!in.docked) {
    const double boostFrac = std::clamp(in.boostAppliedFrac, 0.0, 1.0);
    heatIn += params.heatPerBoostSec * boostFrac;
    if (in.supercruiseActive) heatIn += params.heatPerSupercruiseSec;
    if (in.fsd == ThermalFsdState::Charging) heatIn += params.heatPerFsdChargeSec;
    if (in.fsd == ThermalFsdState::Jumping) heatIn += params.heatPerFsdJumpSec;
  }

  // Cooling rate scales with ship cooling stat.
  const double baseCool = in.docked ? params.baseCoolDocked : params.baseCoolUndocked;
  const double ref = std::max(1e-6, params.referenceCoolRate);
  const double coolStat = std::max(0.0, in.heatCoolRate);
  const double coolRate = baseCool * (coolStat / ref);

  if (dt > 0.0) {
    heat += (heatIn - coolRate) * dt;
  }

  heat = std::clamp(heat, params.heatMin, params.heatMax);

  // Overheat damage.
  double hullDamage = 0.0;
  const double hullMax = std::max(0.0, in.hullMax);
  if (hullMax > 0.0 && heat > params.overheatStart && dt > 0.0) {
    const double over = heat - params.overheatStart;
    hullDamage = std::max(0.0, over) * params.overheatHullDamageRate * hullMax * dt;
  }

  out.heat = heat;
  out.hullDamage = hullDamage;
  out.overheated = (heat > params.overheatStart + 1e-9);
  out.heatInRate = heatIn;
  out.coolRate = coolRate;
  return out;
}

} // namespace stellar::sim
