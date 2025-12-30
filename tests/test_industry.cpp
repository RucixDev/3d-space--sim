#include "stellar/sim/Industry.h"

#include "stellar/core/Types.h"

#include <cmath>
#include <iostream>

namespace {

static bool nearly(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

} // namespace

int test_industry() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;
  using namespace stellar::econ;

  // Basic availability checks.
  {
    const auto refinery = availableIndustryRecipes(StationType::Refinery);
    if (refinery.empty()) {
      std::cerr << "[industry] expected refinery to have recipes\n";
      ++fails;
    }
    const auto outpost = availableIndustryRecipes(StationType::Outpost);
    if (!outpost.empty()) {
      std::cerr << "[industry] expected outpost to have no recipes\n";
      ++fails;
    }
  }

  // Determinism: station modifiers and quotes must be stable across calls.
  {
    const StationId stId = 123456789ull;
    const auto m1 = stationIndustryModifiers(stId);
    const auto m2 = stationIndustryModifiers(stId);
    if (!nearly(m1.yieldMul, m2.yieldMul) || !nearly(m1.speedMul, m2.speedMul) || !nearly(m1.feeMul, m2.feeMul)) {
      std::cerr << "[industry] stationIndustryModifiers not deterministic\n";
      ++fails;
    }

    const auto* r = findIndustryRecipe(IndustryRecipeId::SmeltOre);
    if (!r) {
      std::cerr << "[industry] missing SmeltOre recipe\n";
      ++fails;
    } else {
      const auto q1 = quoteIndustryOrder(*r, stId, StationType::Refinery, 10.0, 0.05, 25.0);
      const auto q2 = quoteIndustryOrder(*r, stId, StationType::Refinery, 10.0, 0.05, 25.0);
      if (!nearly(q1.outputUnits, q2.outputUnits) || !nearly(q1.timeDays, q2.timeDays) || !nearly(q1.serviceFeeCr, q2.serviceFeeCr)) {
        std::cerr << "[industry] quoteIndustryOrder not deterministic\n";
        ++fails;
      }
      if (q1.inputA != CommodityId::Ore || q1.output != CommodityId::Metals) {
        std::cerr << "[industry] SmeltOre I/O mismatch\n";
        ++fails;
      }
      if (!(q1.inputAUnits > 0.0) || !(q1.outputUnits > 0.0) || !(q1.timeDays > 0.0)) {
        std::cerr << "[industry] expected positive quote values\n";
        ++fails;
      }
    }
  }

  // Clamping: negative batches should result in a zero-quote.
  {
    const StationId stId = 42;
    const auto* r = findIndustryRecipe(IndustryRecipeId::SmeltOre);
    if (r) {
      const auto q = quoteIndustryOrder(*r, stId, StationType::Refinery, -1.0, 0.05, 0.0);
      if (!nearly(q.batches, 0.0) || !nearly(q.outputUnits, 0.0) || !nearly(q.timeDays, 0.0)) {
        std::cerr << "[industry] negative batches did not clamp to zero\n";
        ++fails;
      }
    }
  }

  return fails;
}
