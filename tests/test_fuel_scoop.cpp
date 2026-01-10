#include "stellar/sim/FuelScoopSystem.h"

#include "stellar/sim/Gravity.h"

#include <cassert>

int test_fuel_scoop() {
  using namespace stellar;

  sim::Star star{};
  star.radiusSol = 1.0;
  star.luminositySol = 1.0;

  const double rKm = sim::radiusStarKm(star);

  sim::FuelScoopParams p{};
  p.minRangeR = 1.05;
  p.maxRangeR = 8.0;
  p.starHeatRangeR = 12.0;

  // Star heat: inverse-square style falloff (2R should be hotter than 4R).
  const auto h2 = sim::computeFuelScoopRates(star, rKm * 2.0, 0, false, p);
  const auto h4 = sim::computeFuelScoopRates(star, rKm * 4.0, 0, false, p);
  assert(h2.starHeatPerSec > h4.starHeatPerSec);
  assert(h4.starHeatPerSec > 0.0);

  // Out of scoop range: no fuel even if deployed.
  const auto far = sim::computeFuelScoopRates(star, rKm * 20.0, 3, true, p);
  assert(!far.inRange);
  assert(!far.active);
  assert(far.fuelPerSimSec == 0.0);

  // In range but not deployed: no fuel.
  const auto notDeployed = sim::computeFuelScoopRates(star, rKm * 2.0, 3, false, p);
  assert(notDeployed.inRange);
  assert(!notDeployed.active);
  assert(notDeployed.fuelPerSimSec == 0.0);

  // Mk3 should scoop faster (and generate more scoop heat) than Mk1 at the same distance.
  const auto mk1 = sim::computeFuelScoopRates(star, rKm * 2.0, 1, true, p);
  const auto mk3 = sim::computeFuelScoopRates(star, rKm * 2.0, 3, true, p);
  assert(mk1.active && mk3.active);
  assert(mk3.fuelPerSimSec > mk1.fuelPerSimSec);
  assert(mk3.scoopHeatPerSec > mk1.scoopHeatPerSec);

  // Closer should scoop faster than farther (within the window).
  const auto close = sim::computeFuelScoopRates(star, rKm * 1.2, 2, true, p);
  const auto mid = sim::computeFuelScoopRates(star, rKm * 6.0, 2, true, p);
  assert(close.fuelPerSimSec > mid.fuelPerSimSec);
  assert(close.scoopFactor01 > mid.scoopFactor01);

  return 0;
}
