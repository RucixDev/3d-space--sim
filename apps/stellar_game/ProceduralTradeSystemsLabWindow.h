#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/proc/HyperlaneRouter.h"
#include "stellar/proc/TradeFlow.h"
#include "stellar/proc/TradeProfile.h"
#include "stellar/sim/Celestial.h"

#include <functional>
#include <string>
#include <vector>

namespace stellar::sim {
class Universe;
struct StarSystem;
} // namespace stellar::sim

namespace stellar::game {

struct ProceduralTradeSystemsLabWindowState {
  bool open{false};

  // Query controls.
  bool centerOnCurrentSystem{true};
  math::Vec3d centerLy{0.0, 0.0, 0.0};
  double radiusLy{500.0};
  int maxSystems{256};

  // Display controls.
  int topNCommodities{3};
  int sortMode{0}; // 0=distance, 1=hub, 2=wealth, 3=lawlessness, 4=population
  bool showExportsImports{true};
  bool showMacroFactors{true};

  // Macro trade-route suggestion settings (profile-driven; read-only signal).
  bool showTradeRoutes{true};
  bool useHyperlaneRouting{false};
  proc::HyperlaneParams hyperlaneParams{};
  proc::HyperlaneTravelParams hyperlaneTravel{};
  double lastHyperlaneMs{0.0};
  proc::HyperlaneNetwork hyperlanes{};
  int tradeRouteMax{8};
  double tradeRouteMaxDistanceLy{0.0};
  double tradeRouteDistanceExponent{1.35};

  // Interstellar traffic / corridor estimation (macro; routed on hyperlanes).
  bool showInterstellarTraffic{true};
  proc::TradeFlowParams flowParams{};
  bool flowDirty{true};
  double lastFlowMs{0.0};
  proc::TradeFlowNetwork flowNet{};

  // Macro route economics (profile-driven price/profit estimates).
  bool showRouteEconomics{true};
  bool showTradeLoops{true};
  double econBidAskSpread{0.10};
  double econBuyFeeRate{0.00};
  double econSellFeeRate{0.00};
  double econCargoKg{420.0};
  int econMaxLoops{8};

  // Selected system row (index into stubs/profiles). -1 = none.
  int selectedIndex{-1};

  // Cached results.
  bool dirty{true};
  core::u64 lastUniverseSeed{0};
  math::Vec3d lastCenterLy{0.0, 0.0, 0.0};
  double lastRadiusLy{0.0};
  int lastMaxSystems{0};
  double lastBuildMs{0.0};

  std::vector<sim::SystemStub> stubs{};
  std::vector<proc::TradeProfile> profiles{};
  std::vector<double> distsLy{};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

void drawProceduralTradeSystemsLabWindow(ProceduralTradeSystemsLabWindowState& state,
                                        sim::Universe& universe,
                                        const sim::StarSystem* currentSystem,
                                        float timeSec,
                                        const ToastFn& toast);

} // namespace stellar::game
