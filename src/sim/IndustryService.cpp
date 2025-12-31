#include "stellar/sim/IndustryService.h"

#include "stellar/econ/Cargo.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace stellar::sim {
namespace {

static double finiteOr(double v, double fallback = 0.0) {
  return std::isfinite(v) ? v : fallback;
}

static double stationQueueEndDay(const std::vector<IndustryOrder>& orders, StationId stationId) {
  double end = -std::numeric_limits<double>::infinity();
  for (const auto& o : orders) {
    if (o.stationId != stationId) continue;
    const double rd = finiteOr(o.readyDay, -std::numeric_limits<double>::infinity());
    end = std::max(end, rd);
  }
  return end;
}

static constexpr double kEpsUnits = 0.01;

} // namespace

IndustrySubmitResult submitIndustryOrder(core::u64& nextOrderId,
                                        std::vector<IndustryOrder>& orders,
                                        std::array<double, econ::kCommodityCount>& shipCargo,
                                        double& credits,
                                        const Station& station,
                                        double timeDays,
                                        double rep,
                                        const IndustryRecipeDef& recipe,
                                        int batches,
                                        double effectiveStationFeeRate) {
  IndustrySubmitResult out{};

  timeDays = finiteOr(timeDays, 0.0);
  rep = finiteOr(rep, 0.0);
  credits = finiteOr(credits, 0.0);

  if (batches <= 0) {
    out.ok = false;
    out.reason = "no_batches";
    return out;
  }

  const double batchesD = (double)batches;
  out.quote = quoteIndustryOrder(recipe, station.id, station.type, batchesD, effectiveStationFeeRate, rep);
  if (out.quote.batches <= 0.0 || out.quote.outputUnits <= 0.0 || out.quote.timeDays <= 0.0) {
    out.ok = false;
    out.reason = "no_service";
    return out;
  }

  // Input checks.
  const std::size_t ia = static_cast<std::size_t>(out.quote.inputA);
  const std::size_t ib = static_cast<std::size_t>(out.quote.inputB);

  const double haveA = finiteOr(shipCargo[ia], 0.0);
  const double needA = finiteOr(out.quote.inputAUnits, 0.0);
  if (haveA + 1e-6 < needA) {
    out.ok = false;
    out.reason = "need_inputs";
    return out;
  }

  if (out.quote.inputBUnits > 0.0) {
    const double haveB = finiteOr(shipCargo[ib], 0.0);
    const double needB = finiteOr(out.quote.inputBUnits, 0.0);
    if (haveB + 1e-6 < needB) {
      out.ok = false;
      out.reason = "need_inputs";
      return out;
    }
  }

  const double fee = std::max(0.0, finiteOr(out.quote.serviceFeeCr, 0.0));
  if (credits + 1e-6 < fee) {
    out.ok = false;
    out.reason = "need_credits";
    return out;
  }

  // Deduct inputs + fee.
  shipCargo[ia] = std::max(0.0, haveA - needA);
  if (out.quote.inputBUnits > 0.0) {
    shipCargo[ib] = std::max(0.0, finiteOr(shipCargo[ib], 0.0) - finiteOr(out.quote.inputBUnits, 0.0));
  }

  credits -= fee;
  out.creditsPaid = fee;

  // Queueing: station processes one order at a time.
  const double queueEnd = stationQueueEndDay(orders, station.id);
  const double start = std::max(timeDays, queueEnd);

  IndustryOrder o{};
  o.id = nextOrderId;
  nextOrderId = std::max<core::u64>(1, nextOrderId + 1);
  o.recipe = recipe.id;
  o.stationId = station.id;
  o.inputA = out.quote.inputA;
  o.inputAUnits = out.quote.inputAUnits;
  o.inputB = out.quote.inputB;
  o.inputBUnits = out.quote.inputBUnits;
  o.output = out.quote.output;
  o.outputUnits = out.quote.outputUnits;
  o.submittedDay = timeDays;
  o.readyDay = start + out.quote.timeDays;
  o.claimed = false;

  orders.push_back(o);

  out.order = o;
  out.ok = true;
  out.reason = nullptr;
  return out;
}

IndustryClaimResult claimIndustryOrderToCargo(IndustryOrder& order,
                                             const Station& station,
                                             double timeDays,
                                             std::array<double, econ::kCommodityCount>& shipCargo,
                                             double cargoCapacityKg) {
  IndustryClaimResult out{};
  timeDays = finiteOr(timeDays, 0.0);
  cargoCapacityKg = std::max(0.0, finiteOr(cargoCapacityKg, 0.0));

  if (order.stationId != station.id) {
    out.ok = false;
    out.reason = "wrong_station";
    return out;
  }
  if (order.claimed) {
    out.ok = false;
    out.reason = "claimed";
    out.completed = true;
    return out;
  }
  if (timeDays < finiteOr(order.readyDay, 0.0) - 1e-6) {
    out.ok = false;
    out.reason = "not_ready";
    return out;
  }

  double remaining = std::max(0.0, finiteOr(order.outputUnits, 0.0));
  if (remaining <= kEpsUnits) {
    order.outputUnits = 0.0;
    order.claimed = true;
    out.ok = false;
    out.reason = "empty";
    out.completed = true;
    return out;
  }

  const double outMassKg = std::max(0.0, econ::commodityDef(order.output).massKg);
  const double usedKg = std::max(0.0, econ::cargoMassKg(shipCargo));
  const double freeKg = std::max(0.0, cargoCapacityKg - usedKg);

  const double maxUnitsByMass = (outMassKg > 1e-9) ? (freeKg / outMassKg) : remaining;
  const double moved = std::min(remaining, std::max(0.0, maxUnitsByMass));
  if (moved <= kEpsUnits) {
    out.ok = false;
    out.reason = "cargo_full";
    return out;
  }

  shipCargo[static_cast<std::size_t>(order.output)] += moved;
  remaining -= moved;
  order.outputUnits = remaining;

  out.ok = true;
  out.reason = nullptr;
  out.unitsMoved = moved;

  if (order.outputUnits <= kEpsUnits) {
    order.outputUnits = 0.0;
    order.claimed = true;
    out.completed = true;
  }

  return out;
}

IndustryClaimResult moveIndustryOrderOutputToWarehouse(IndustryOrder& order,
                                                       const Station& station,
                                                       double timeDays,
                                                       double rep,
                                                       std::vector<StationStorage>& stationStorage) {
  IndustryClaimResult out{};
  timeDays = finiteOr(timeDays, 0.0);
  rep = finiteOr(rep, 0.0);

  if (order.stationId != station.id) {
    out.ok = false;
    out.reason = "wrong_station";
    return out;
  }
  if (order.claimed) {
    out.ok = false;
    out.reason = "claimed";
    out.completed = true;
    return out;
  }
  if (timeDays < finiteOr(order.readyDay, 0.0) - 1e-6) {
    out.ok = false;
    out.reason = "not_ready";
    return out;
  }

  const double remaining = std::max(0.0, finiteOr(order.outputUnits, 0.0));
  if (remaining <= kEpsUnits) {
    order.outputUnits = 0.0;
    order.claimed = true;
    out.ok = false;
    out.reason = "empty";
    out.completed = true;
    return out;
  }

  auto& e = getOrCreateStorage(stationStorage, station, timeDays);
  accrueStorageFees(e, timeDays, rep);
  e.cargo[static_cast<std::size_t>(order.output)] += remaining;

  order.outputUnits = 0.0;
  order.claimed = true;

  out.ok = true;
  out.reason = nullptr;
  out.unitsMoved = remaining;
  out.completed = true;
  return out;
}

void pruneClaimedIndustryOrders(std::vector<IndustryOrder>& orders) {
  orders.erase(std::remove_if(orders.begin(), orders.end(),
                              [](const IndustryOrder& o) { return o.claimed; }),
               orders.end());
}

} // namespace stellar::sim
