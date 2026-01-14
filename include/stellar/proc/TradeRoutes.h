#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/proc/HyperlaneRouter.h"
#include "stellar/proc/TradeProfile.h"
#include "stellar/sim/Celestial.h"

#include <vector>

namespace stellar::proc {

// Procedural, profile-driven trade route suggestions.
//
// This is intentionally a *macro* signal — it does not inspect live station inventories
// or prices. Instead it uses the deterministic TradeProfile export/import scores to
// estimate which nearby systems are likely trade partners.
//
// Primary use-cases:
//  - A "galaxy econ" lab window / map overlay.
//  - Seeding NPC logistics / trade lane simulation in later rounds.
//  - Biasing mission generation toward plausible origin->destination pairs.

struct TradeRoute {
  // For an export route: origin -> otherSystem.
  // For an import route: otherSystem -> origin.
  sim::SystemId otherSystem{0};

  // Commodity that best matches the route direction.
  econ::CommodityId commodity{econ::CommodityId::Food};

  // Dimensionless score (higher = stronger macro signal).
  double score{0.0};

  // Route distance metric in light-years.
  //
  // - In the classic overload, this is the straight-line distance.
  // - In the hyperlane-aware overload, this is an *effective travel cost* (ly-like)
  //   computed by routing over the hyperlane network (distance modulated by
  //   bandwidth and risk).
  double distanceLy{0.0};

  // Straight-line system-to-system distance for reference.
  double directDistanceLy{0.0};

  // Hyperlane path metrics (all zero if the route was scored geometrically).
  double laneDistanceLy{0.0};          // physical lane path length (sum of edge distances)
  int laneHops{0};                     // hop count along the lane path
  double laneBottleneckBandwidth01{0.0}; // min bandwidth along the path

  // Rough risk proxy in [0,1] (currently based on lawlessness).
  double risk01{0.0};
};

struct TradeRouteParams {
  // Maximum number of routes to return for exports and for imports.
  std::size_t maxRoutes{8};

  // Optional filter for candidates farther than this. <= 0 means "no filter".
  double maxDistanceLy{0.0};

  // How strongly distance suppresses trade affinity.
  // Typical values: 1.0 .. 2.0
  double distanceExponent{1.35};

  // Softening factor to avoid blowing up extremely short distances.
  // Effective distance = max(distanceLy, distanceSofteningLy)
  double distanceSofteningLy{4.0};

  // If true, same-faction pairs receive a small score bonus.
  bool sameFactionBonus{true};

  // If true, routes that have extremely low commodity affinity are dropped.
  bool dropWeakRoutes{true};

  // Minimum raw commodity affinity (export*import) required when dropWeakRoutes is enabled.
  double minAffinity{0.02};
};

struct TradeRouteSuggestions {
  // Macro export routes: origin -> other
  std::vector<TradeRoute> exports;

  // Macro import routes: other -> origin
  std::vector<TradeRoute> imports;
};

// Compute trade route suggestions using precomputed TradeProfiles.
//
// Requirements:
//  - candidates.size() == candidateProfiles.size()
//
// Notes:
//  - The origin stub may appear in candidates; it is automatically ignored.
//  - Output order is deterministic (score desc, distance asc, id asc).
TradeRouteSuggestions suggestTradeRoutes(const sim::SystemStub& origin,
                                        const TradeProfile& originProfile,
                                        const std::vector<sim::SystemStub>& candidates,
                                        const std::vector<TradeProfile>& candidateProfiles,
                                        const TradeRouteParams& params = {});

// Hyperlane-aware overload: distance/risk are computed along the supplied
// lane network using HyperlaneRouter.
//
// Notes:
//  - Unreachable candidates are skipped.
//  - TradeRoute::distanceLy is set to the lane travel cost (see struct docs).
TradeRouteSuggestions suggestTradeRoutes(const sim::SystemStub& origin,
                                        const TradeProfile& originProfile,
                                        const std::vector<sim::SystemStub>& candidates,
                                        const std::vector<TradeProfile>& candidateProfiles,
                                        const HyperlaneNetwork& hyperlanes,
                                        const TradeRouteParams& params = {},
                                        const HyperlaneTravelParams& travelParams = {});

// Convenience overload: computes TradeProfiles internally.
TradeRouteSuggestions suggestTradeRoutes(core::u64 universeSeed,
                                        const sim::SystemStub& origin,
                                        const std::vector<sim::SystemStub>& candidates,
                                        const TradeRouteParams& params = {});

// Hyperlane-aware convenience overload (computes TradeProfiles internally).
TradeRouteSuggestions suggestTradeRoutes(core::u64 universeSeed,
                                        const sim::SystemStub& origin,
                                        const std::vector<sim::SystemStub>& candidates,
                                        const HyperlaneNetwork& hyperlanes,
                                        const TradeRouteParams& params = {},
                                        const HyperlaneTravelParams& travelParams = {});

} // namespace stellar::proc
