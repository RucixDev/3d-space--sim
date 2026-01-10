#include "stellar/sim/FuelScoopSystem.h"

#include "stellar/sim/Gravity.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

FuelScoopRates computeFuelScoopRates(const Star& star,
                                     double distKmFromStarCenter,
                                     int fuelScoopMk,
                                     bool deployed,
                                     const FuelScoopParams& params) {
  FuelScoopRates out{};
  out.distKm = std::max(0.0, distKmFromStarCenter);

  const double rStarKm = std::max(1.0, radiusStarKm(star));
  out.distR = out.distKm / rStarKm;

  // Star heat (simple inverse-square falloff).
  if (out.distR > 1e-9 && out.distR <= std::max(0.0, params.starHeatRangeR)) {
    const double lum = std::clamp(star.luminositySol, 0.05, 50.0);
    const double lumMul = std::sqrt(lum);

    const double inv = 2.0 / out.distR;
    const double heat = params.starHeatPerSecAt2R * lumMul * inv * inv;
    out.starHeatPerSec = std::clamp(heat, 0.0, std::max(0.0, params.starHeatMaxPerSec));
  }

  const int mk = std::clamp(fuelScoopMk, 0, 3);
  out.inRange = (out.distR >= params.minRangeR && out.distR <= params.maxRangeR);

  if (!deployed || mk <= 0 || !out.inRange) {
    return out;
  }

  out.active = true;

  // 0..1, where 1 = minRangeR and 0 = maxRangeR.
  const double denom = std::max(1e-6, params.maxRangeR - params.minRangeR);
  const double t = std::clamp((out.distR - params.minRangeR) / denom, 0.0, 1.0);
  out.scoopFactor01 = 1.0 - t;

  const double power = std::max(0.05, params.scoopFalloffPower);
  const double shaped = std::pow(out.scoopFactor01, power);

  double mkMult = 1.0;
  if (mk == 2) mkMult = params.scoopMk2Mult;
  else if (mk >= 3) mkMult = params.scoopMk3Mult;

  // Scooping is a "hardware" feature, but we loosely scale it with star luminosity so
  // dim stars feel weaker and bright stars feel more rewarding (and hotter).
  const double lum = std::clamp(star.luminositySol, 0.05, 50.0);
  const double lumMul = std::sqrt(lum);

  out.fuelPerSimSec = params.scoopFuelPerSimSecMk1 * mkMult * shaped * lumMul;
  out.fuelPerSimSec = std::clamp(out.fuelPerSimSec, 0.0, params.scoopFuelMaxPerSimSec);

  out.scoopHeatPerSec = params.scoopHeatPerSecAtMaxMk1 * mkMult * shaped;

  return out;
}

} // namespace stellar::sim
