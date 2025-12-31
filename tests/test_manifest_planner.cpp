#include "stellar/econ/RoutePlanner.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps) { return std::abs(a - b) <= eps; }

int test_manifest_planner() {
  int fails = 0;
  auto expect = [&](bool ok, const char* msg) {
    if (!ok) {
      std::cerr << "[test_manifest_planner] FAIL: " << msg << "\n";
      ++fails;
    }
  };

  // Build two station economies:
  //  - Only Luxury and Food are profitable (dest is undersupplied, origin is oversupplied).
  //  - Luxury has limited supply so the manifest must mix in Food to fill remaining capacity.
  econ::StationEconomyModel fromM{};
  econ::StationEconomyModel toM{};

  fromM.priceVolatility = 1.0;
  toM.priceVolatility = 1.0;
  fromM.shockVolatility = 0.0;
  toM.shockVolatility = 0.0;

  econ::StationEconomyState fromS{};
  econ::StationEconomyState toS{};

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    fromM.capacity[i] = 300.0;
    toM.capacity[i] = 300.0;
    fromM.desiredStock[i] = 100.0;
    toM.desiredStock[i] = 100.0;

    fromS.inventory[i] = 100.0;
    toS.inventory[i] = 100.0;
  }

  // Luxury: cheap at origin, expensive at destination, but origin supply limited to 10 units.
  const auto lux = econ::CommodityId::Luxury;
  const std::size_t luxI = (std::size_t)lux;
  fromM.capacity[luxI] = 10.0;
  fromM.desiredStock[luxI] = 1.0;
  fromS.inventory[luxI] = 10.0;

  toM.capacity[luxI] = 300.0;
  toM.desiredStock[luxI] = 100.0;
  toS.inventory[luxI] = 0.0;

  // Food: cheap at origin (oversupplied), expensive at destination (undersupplied).
  const auto food = econ::CommodityId::Food;
  const std::size_t foodI = (std::size_t)food;
  fromM.capacity[foodI] = 300.0;
  fromM.desiredStock[foodI] = 50.0;
  fromS.inventory[foodI] = 200.0;

  toM.capacity[foodI] = 300.0;
  toM.desiredStock[foodI] = 100.0;
  toS.inventory[foodI] = 0.0;

  econ::CargoManifestParams p{};
  p.cargoCapacityKg = 10.0;
  p.bidAskSpread = 0.10;
  p.fromFeeRate = 0.0;
  p.toFeeRate = 0.0;
  p.stepKg = 1.0;
  p.maxBuyCreditsCr = 0.0;
  p.simulatePriceImpact = true;

  const auto plan = econ::bestManifestForCargo(fromS, fromM, toS, toM, p);

  expect(plan.netProfitCr > 0.0, "plan should be profitable");
  expect(approx(plan.cargoFilledKg, 10.0, 1e-3), "plan should fill the hold (10 kg)");

  double luxUnits = 0.0;
  double foodUnits = 0.0;
  for (const auto& ln : plan.lines) {
    if (ln.commodity == lux) luxUnits = ln.units;
    if (ln.commodity == food) foodUnits = ln.units;
  }

  // With stepKg=1, Luxury mass=0.2 => 5 units/step, limited to 10 units (2 kg).
  // The remaining 8 kg should be filled by Food (mass=1 => 8 units).
  expect(approx(luxUnits, 10.0, 1e-2), "luxury units should hit the supply cap (10)");
  expect(approx(foodUnits, 8.0, 1e-2), "food units should fill remaining mass (8)");

  // Credit-limited planning: ensure we don't exceed maxBuyCreditsCr.
  econ::CargoManifestParams p2 = p;
  p2.maxBuyCreditsCr = 50.0; // small budget
  const auto plan2 = econ::bestManifestForCargo(fromS, fromM, toS, toM, p2);

  expect(plan2.netBuyCr <= 50.0 + 1e-6, "credit-limited plan must respect maxBuyCreditsCr");
  expect(plan2.cargoFilledKg <= p2.cargoCapacityKg + 1e-6, "credit-limited plan must respect cargo capacity");
  expect(plan2.netProfitCr >= 0.0, "credit-limited plan should not produce negative profit");

  if (fails == 0) {
    std::cout << "[test_manifest_planner] PASS\n";
  }
  return fails;
}
