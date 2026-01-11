#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"

#include <array>
#include <vector>

namespace stellar::sim {

// A single line item in a jettison plan.
struct CargoJettisonLine {
  econ::CommodityId commodity{econ::CommodityId::Food};
  double units{0.0};
  double valueCr{0.0};
  double massKg{0.0};
};

// A plan describing which commodities to jettison to reach a required cargo value.
//
// Intended use cases:
//  - pirate extortion ("drop cargo worth X")
//  - emergency mass shedding for jump range
//
// Notes:
//  - This planner treats cargo as discrete whole units (it floors available units).
//  - It minimizes overpay first (plannedValue - requiredValue), then minimizes total mass.
//  - The returned plan may include mission-reserved cargo only when allowUsingReserved=true.
struct CargoJettisonPlan {
  bool success{false};        // plannedValueCr >= requiredValueCr
  bool usedReserved{false};   // plan dips into reservedUnits

  double requiredValueCr{0.0};
  double plannedValueCr{0.0};
  double plannedMassKg{0.0};

  // If success==false, plannedValueCr is the maximum achievable with the allowed cargo.
  std::vector<CargoJettisonLine> lines{};
};

// Plan which cargo to jettison to satisfy a required value.
//
// cargoUnits:     current ship hold units per commodity
// reservedUnits:  mission-reserved units per commodity (can be empty)
// requiredValueCr: value target to reach (credits)
// allowUsingReserved: if true, reserved cargo may be used if needed
CargoJettisonPlan planCargoJettisonForValue(
  const std::array<double, econ::kCommodityCount>& cargoUnits,
  const std::array<double, econ::kCommodityCount>& reservedUnits,
  double requiredValueCr,
  bool allowUsingReserved);

} // namespace stellar::sim
