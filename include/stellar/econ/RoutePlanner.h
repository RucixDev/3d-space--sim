#pragma once

#include "stellar/econ/Market.h"

#include <vector>

namespace stellar::econ {

struct RouteOpportunity {
  CommodityId commodity{};
  double profitPerUnit{0.0};
  double buyPrice{0.0};
  double sellPrice{0.0};
};

std::vector<RouteOpportunity> bestRoutes(const StationEconomyState& fromState,
                                        const StationEconomyModel& fromModel,
                                        const StationEconomyState& toState,
                                        const StationEconomyModel& toModel,
                                        double bidAskSpread = 0.10,
                                        std::size_t maxResults = 5);

} // namespace stellar::econ
