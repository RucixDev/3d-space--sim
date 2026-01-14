#include "stellar/econ/RoutePlanner.h"

#include "test_harness.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static constexpr std::size_t cidx(econ::CommodityId id) {
  return static_cast<std::size_t>(id);
}

static const econ::CargoManifestLine* findLine(const econ::CargoManifestPlan& plan, econ::CommodityId id) {
  for (const auto& l : plan.lines) {
    if (l.commodity == id) return &l;
  }
  return nullptr;
}

int test_manifest_planner_beam() {
  // Crafted setup:
  //  - Luxury has enormous profit/kg but is so expensive that a tight credit cap
  //    prevents filling the hold.
  //  - Food is cheaper; filling the hold yields higher total profit.
  // Greedy (profit/kg) gets baited into buying only Luxury.
  // Beam search optimizes total profit under the credit constraint.

  econ::StationEconomyModel fromM{};
  econ::StationEconomyModel toM{};
  fromM.capacity.fill(1000.0);
  toM.capacity.fill(1000.0);
  fromM.desiredStock.fill(100.0);
  toM.desiredStock.fill(100.0);
  fromM.priceVolatility = 1.0;
  toM.priceVolatility = 1.0;
  fromM.shockVolatility = 0.0;
  toM.shockVolatility = 0.0;

  econ::StationEconomyState fromS{};
  econ::StationEconomyState toS{};
  fromS.inventory.fill(100.0);
  toS.inventory.fill(100.0);

  // Make Food cheaper at origin and more expensive at destination.
  fromS.inventory[cidx(econ::CommodityId::Food)] = 130.0; // oversupplied => cheaper
  toS.inventory[cidx(econ::CommodityId::Food)] = 50.0;    // scarce => more expensive

  // Make Luxury profitable but expensive.
  fromS.inventory[cidx(econ::CommodityId::Luxury)] = 120.0; // slightly oversupplied
  toS.inventory[cidx(econ::CommodityId::Luxury)] = 80.0;    // somewhat scarce

  econ::CargoManifestParams p{};
  p.cargoCapacityKg = 10.0;
  p.bidAskSpread = 0.10;
  p.fromFeeRate = 0.0;
  p.toFeeRate = 0.0;
  p.stepKg = 1.0;
  p.maxBuyCreditsCr = 200.0; // tight credit cap
  p.simulatePriceImpact = false;
  p.planner = econ::CargoManifestPlanner::Greedy;

  const auto greedy = econ::bestManifestForCargo(fromS, fromM, toS, toM, p);

  p.planner = econ::CargoManifestPlanner::BeamSearch;
  p.beamWidth = 8;
  const auto beam = econ::bestManifestForCargo(fromS, fromM, toS, toM, p);

  int failures = 0;

  CHECK(greedy.netProfitCr > 1.0);
  CHECK(beam.netProfitCr > greedy.netProfitCr + 1.0);

  const auto* greedyLux = findLine(greedy, econ::CommodityId::Luxury);
  CHECK(greedyLux != nullptr);

  const auto* beamFood = findLine(beam, econ::CommodityId::Food);
  CHECK(beamFood != nullptr);

  // Beam should usually fill substantially more mass than greedy in this setup.
  CHECK(beam.cargoFilledKg > greedy.cargoFilledKg + 1.0);

  if (failures == 0) {
    std::cout << "[test_manifest_planner_beam] PASS\n";
  }
  return failures;
}
