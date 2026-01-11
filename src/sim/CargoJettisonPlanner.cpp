#include "stellar/sim/CargoJettisonPlanner.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace stellar::sim {

namespace {

struct Item {
  int valueCr{0};
  double massKg{0.0};
  int commodityIdx{0};
  int units{0};
};

static int toIntUnits(double units) {
  // Cargo is generally integral in gameplay, but we store it as double.
  // We floor to keep plans stable and to avoid planning fractional pods.
  return (int)std::floor(std::max(0.0, units) + 1e-6);
}

static int toIntCredits(double cr) {
  // Credits are effectively integral for this use case.
  return (int)std::ceil(std::max(0.0, cr) - 1e-9);
}

} // namespace


CargoJettisonPlan planCargoJettisonForValue(
  const std::array<double, econ::kCommodityCount>& cargoUnits,
  const std::array<double, econ::kCommodityCount>& reservedUnits,
  double requiredValueCr,
  bool allowUsingReserved) {

  CargoJettisonPlan out{};
  out.requiredValueCr = std::max(0.0, requiredValueCr);

  const int req = toIntCredits(requiredValueCr);
  if (req <= 0) {
    out.success = true;
    out.plannedValueCr = 0.0;
    out.plannedMassKg = 0.0;
    return out;
  }

  std::array<int, econ::kCommodityCount> have{};
  std::array<int, econ::kCommodityCount> free{};
  have.fill(0);
  free.fill(0);

  int totalValue = 0;
  int maxUnitValue = 0;

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const int h = toIntUnits(cargoUnits[i]);
    const int r = toIntUnits(reservedUnits[i]);
    const int f = std::max(0, h - r);
    have[i] = h;
    free[i] = f;

    const int price = (int)std::lround(econ::commodityDef((econ::CommodityId)i).basePrice);
    const int avail = allowUsingReserved ? h : f;
    if (avail <= 0 || price <= 0) continue;
    totalValue += avail * price;
    maxUnitValue = std::max(maxUnitValue, price);
  }

  if (totalValue <= 0) {
    out.success = false;
    out.plannedValueCr = 0.0;
    out.plannedMassKg = 0.0;
    return out;
  }

  // If we can't satisfy the demand, return a "dump everything" plan.
  if (totalValue < req) {
    out.success = false;
    out.plannedValueCr = (double)totalValue;
    out.plannedMassKg = 0.0;
    for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
      const int avail = allowUsingReserved ? have[i] : free[i];
      if (avail <= 0) continue;
      const auto cid = (econ::CommodityId)i;
      const auto& def = econ::commodityDef(cid);
      CargoJettisonLine ln{};
      ln.commodity = cid;
      ln.units = (double)avail;
      ln.valueCr = ln.units * def.basePrice;
      ln.massKg = ln.units * def.massKg;
      out.plannedMassKg += ln.massKg;
      out.lines.push_back(ln);

      if (allowUsingReserved && avail > free[i]) out.usedReserved = true;
    }

    std::sort(out.lines.begin(), out.lines.end(), [](const CargoJettisonLine& a, const CargoJettisonLine& b) {
      return a.valueCr > b.valueCr;
    });

    return out;
  }

  // Bounded knapsack via binary splitting.
  std::vector<Item> items;
  items.reserve(64);
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const int avail = allowUsingReserved ? have[i] : free[i];
    if (avail <= 0) continue;

    const auto cid = (econ::CommodityId)i;
    const auto& def = econ::commodityDef(cid);
    const int price = (int)std::lround(def.basePrice);
    if (price <= 0) continue;

    int remaining = avail;
    int chunk = 1;
    while (remaining > 0) {
      const int take = std::min(chunk, remaining);
      Item it{};
      it.valueCr = take * price;
      it.massKg = (double)take * def.massKg;
      it.commodityIdx = (int)i;
      it.units = take;
      items.push_back(it);

      remaining -= take;
      chunk <<= 1;
    }
  }

  // Search window: overpay is typically bounded by the max unit value, but use a
  // generous margin in case the available values are sparse.
  int cap = std::min(totalValue, req + std::max(1, maxUnitValue) * 50);
  cap = std::max(cap, req);

  auto solveWithCap = [&](int capValue,
                          std::vector<double>& dp,
                          std::vector<int>& prevV,
                          std::vector<int>& prevItem) -> int {
    const double INF = std::numeric_limits<double>::infinity();
    dp.assign((std::size_t)capValue + 1, INF);
    prevV.assign((std::size_t)capValue + 1, -1);
    prevItem.assign((std::size_t)capValue + 1, -1);
    dp[0] = 0.0;

    for (int idx = 0; idx < (int)items.size(); ++idx) {
      const Item& it = items[(std::size_t)idx];
      if (it.valueCr <= 0) continue;

      for (int v = capValue; v >= it.valueCr; --v) {
        const double prev = dp[(std::size_t)(v - it.valueCr)];
        if (!std::isfinite(prev)) continue;
        const double cand = prev + it.massKg;
        const double cur = dp[(std::size_t)v];
        if (cand + 1e-12 < cur) {
          dp[(std::size_t)v] = cand;
          prevV[(std::size_t)v] = v - it.valueCr;
          prevItem[(std::size_t)v] = idx;
        }
      }
    }

    for (int v = req; v <= capValue; ++v) {
      if (std::isfinite(dp[(std::size_t)v])) return v;
    }
    return -1;
  };

  std::vector<double> dp;
  std::vector<int> prevV;
  std::vector<int> prevItem;
  int bestV = solveWithCap(cap, dp, prevV, prevItem);

  // Fallback: if the window was too small, expand to totalValue.
  if (bestV < 0 && cap < totalValue) {
    cap = totalValue;
    bestV = solveWithCap(cap, dp, prevV, prevItem);
  }

  if (bestV < 0) {
    // Should not happen if totalValue >= req, but keep a safe fallback.
    out.success = false;
    out.plannedValueCr = 0.0;
    out.plannedMassKg = 0.0;
    return out;
  }

  std::array<int, econ::kCommodityCount> chosen{};
  chosen.fill(0);

  int v = bestV;
  while (v > 0) {
    const int itemIdx = prevItem[(std::size_t)v];
    if (itemIdx < 0) break;
    const Item& it = items[(std::size_t)itemIdx];
    chosen[(std::size_t)it.commodityIdx] += it.units;
    v = prevV[(std::size_t)v];
    if (v < 0) break;
  }

  out.success = (bestV >= req);
  out.plannedValueCr = (double)bestV;
  out.plannedMassKg = 0.0;

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const int u = chosen[i];
    if (u <= 0) continue;
    const auto cid = (econ::CommodityId)i;
    const auto& def = econ::commodityDef(cid);

    CargoJettisonLine ln{};
    ln.commodity = cid;
    ln.units = (double)u;
    ln.valueCr = ln.units * def.basePrice;
    ln.massKg = ln.units * def.massKg;
    out.plannedMassKg += ln.massKg;
    out.lines.push_back(ln);

    if (allowUsingReserved && u > free[i]) out.usedReserved = true;
  }

  std::sort(out.lines.begin(), out.lines.end(), [](const CargoJettisonLine& a, const CargoJettisonLine& b) {
    // Sort by value first, then by mass (descending value is easier to read in UI).
    if (std::abs(a.valueCr - b.valueCr) > 1e-6) return a.valueCr > b.valueCr;
    return a.massKg > b.massKg;
  });

  return out;
}

} // namespace stellar::sim
