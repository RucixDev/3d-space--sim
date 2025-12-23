#include "stellar/econ/RoutePlanner.h"

#include <algorithm>

namespace stellar::econ {

static constexpr std::size_t idx(CommodityId id) { return static_cast<std::size_t>(id); }

std::vector<RouteOpportunity> bestRoutes(const StationEconomyState& fromState,
                                        const StationEconomyModel& fromModel,
                                        const StationEconomyState& toState,
                                        const StationEconomyModel& toModel,
                                        double bidAskSpread,
                                        std::size_t maxResults) {
  std::vector<RouteOpportunity> out;
  out.reserve(kCommodityCount);

  for (std::size_t i = 0; i < kCommodityCount; ++i) {
    const CommodityId cid = static_cast<CommodityId>(i);
    const auto qFrom = quote(fromState, fromModel, cid, bidAskSpread);
    const auto qTo   = quote(toState, toModel, cid, bidAskSpread);

    const double profit = qTo.bid - qFrom.ask;
    if (profit > 0.0) {
      RouteOpportunity r{};
      r.commodity = cid;
      r.profitPerUnit = profit;
      r.buyPrice = qFrom.ask;
      r.sellPrice = qTo.bid;
      out.push_back(r);
    }
  }

  std::sort(out.begin(), out.end(), [](const RouteOpportunity& a, const RouteOpportunity& b) {
    return a.profitPerUnit > b.profitPerUnit;
  });

  if (out.size() > maxResults) out.resize(maxResults);
  return out;
}

} // namespace stellar::econ
