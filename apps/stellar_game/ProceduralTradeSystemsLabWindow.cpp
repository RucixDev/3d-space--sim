#include "ProceduralTradeSystemsLabWindow.h"

#include "stellar/econ/Commodity.h"
#include "stellar/proc/GalaxyHazards.h"
#include "stellar/proc/TradeIntel.h"
#include "stellar/proc/TradeRoutes.h"
#include "stellar/sim/Faction.h"
#include "stellar/sim/Universe.h"

#include <imgui.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace stellar::game {

using Clock = std::chrono::high_resolution_clock;

static const char* sortModeName(int mode) {
  switch (mode) {
    case 0: return "Distance";
    case 1: return "Hub";
    case 2: return "Wealth";
    case 3: return "Lawlessness";
    case 4: return "Population";
    default: return "Distance";
  }
}

static std::string topCommodities(const std::array<double, econ::kCommodityCount>& scores,
                                 int topN) {
  topN = std::clamp(topN, 1, (int)econ::kCommodityCount);
  std::vector<std::pair<double, econ::CommodityId>> v;
  v.reserve(econ::kCommodityCount);
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    v.emplace_back(scores[i], static_cast<econ::CommodityId>(i));
  }

  std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) {
    if (a.first != b.first) return a.first > b.first;
    return (std::size_t)a.second < (std::size_t)b.second;
  });

  std::ostringstream oss;
  for (int i = 0; i < topN; ++i) {
    if (i) oss << " ";
    oss << econ::commodityCode(v[i].second);
  }
  return oss.str();
}

static std::string fmtCr(double cr) {
  const double a = std::abs(cr);
  std::ostringstream oss;
  oss.setf(std::ios::fixed);
  if (a >= 1e9) {
    oss << std::setprecision(2) << (cr / 1e9) << "B";
  } else if (a >= 1e6) {
    oss << std::setprecision(2) << (cr / 1e6) << "M";
  } else if (a >= 1e3) {
    oss << std::setprecision(1) << (cr / 1e3) << "k";
  } else {
    oss << std::setprecision(0) << cr;
  }
  return oss.str();
}

static const sim::Faction* findFaction(core::u32 id, const std::vector<sim::Faction>& factions) {
  if (id == 0) return nullptr;
  for (const auto& f : factions) {
    if (f.id == id) return &f;
  }
  return nullptr;
}

static const sim::SystemStub* findStub(sim::SystemId id, const std::vector<sim::SystemStub>& stubs) {
  for (const auto& s : stubs) {
    if (s.id == id) return &s;
  }
  return nullptr;
}

static const proc::HyperlaneEdge* findLaneEdge(sim::SystemId a, sim::SystemId b, const proc::HyperlaneNetwork& net) {
  if (a > b) std::swap(a, b);
  for (const auto& e : net.edges) {
    if (e.a == a && e.b == b) return &e;
  }
  return nullptr;
}

static double distLy(const math::Vec3d& a, const math::Vec3d& b) {
  const double dx = a.x - b.x;
  const double dy = a.y - b.y;
  const double dz = a.z - b.z;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

static void rebuild(ProceduralTradeSystemsLabWindowState& st, sim::Universe& u, const math::Vec3d& centerLy) {
  const auto t0 = Clock::now();
  sim::SystemId selectedId = 0;
  if (st.selectedIndex >= 0 && (std::size_t)st.selectedIndex < st.stubs.size()) {
    selectedId = st.stubs[(std::size_t)st.selectedIndex].id;
  }


  st.stubs.clear();
  st.profiles.clear();
  st.distsLy.clear();

  const int maxN = std::clamp(st.maxSystems, 1, 10'000);
  const double radius = std::max(10.0, st.radiusLy);

  auto stubs = u.queryNearby(centerLy, radius, (std::size_t)maxN);

  st.stubs = std::move(stubs);
  st.profiles.reserve(st.stubs.size());
  st.distsLy.reserve(st.stubs.size());

  for (const auto& s : st.stubs) {
    st.profiles.push_back(proc::generateTradeProfile(u.seed(), s));
    st.distsLy.push_back(distLy(centerLy, s.posLy));
  }

  // Optional sorting.
  std::vector<std::size_t> idx(st.stubs.size());
  for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = i;

  const auto keyFor = [&](std::size_t i) -> double {
    const auto& p = st.profiles[i];
    switch (st.sortMode) {
      case 1: return p.hub;
      case 2: return p.wealth;
      case 3: return p.lawlessness;
      case 4: return p.population;
      default: return -st.distsLy[i]; // distance ascending
    }
  };

  if (st.sortMode == 0) {
    // Already ordered by deterministic distance from Universe::queryNearby.
  } else {
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
      const double ka = keyFor(a);
      const double kb = keyFor(b);
      if (ka != kb) return ka > kb;
      // Tie-breaker: nearer first, then id.
      if (st.distsLy[a] != st.distsLy[b]) return st.distsLy[a] < st.distsLy[b];
      return st.stubs[a].id < st.stubs[b].id;
    });

    std::vector<sim::SystemStub> stubsSorted;
    std::vector<proc::TradeProfile> profilesSorted;
    std::vector<double> distsSorted;
    stubsSorted.reserve(idx.size());
    profilesSorted.reserve(idx.size());
    distsSorted.reserve(idx.size());

    for (std::size_t j : idx) {
      stubsSorted.push_back(st.stubs[j]);
      profilesSorted.push_back(st.profiles[j]);
      distsSorted.push_back(st.distsLy[j]);
    }

    st.stubs = std::move(stubsSorted);
    st.profiles = std::move(profilesSorted);
    st.distsLy = std::move(distsSorted);
  }

  // Preserve selection across rebuilds by id (sorting/filtering may re-order rows).
  st.selectedIndex = -1;
  if (selectedId != 0) {
    for (std::size_t i = 0; i < st.stubs.size(); ++i) {
      if (st.stubs[i].id == selectedId) {
        st.selectedIndex = (int)i;
        break;
      }
    }
  }

  // Build a hyperlane network over the queried systems when enabled.
  st.hyperlanes.edges.clear();
  st.lastHyperlaneMs = 0.0;

  // Interstellar traffic cache depends on the hyperlane graph and travel params.
  st.flowNet = proc::TradeFlowNetwork{};
  st.lastFlowMs = 0.0;
  st.flowDirty = true;
  if (st.useHyperlaneRouting && st.stubs.size() >= 2) {
    const auto tLane0 = Clock::now();
    std::vector<sim::SystemStub> laneStubs = st.stubs;
    std::sort(laneStubs.begin(), laneStubs.end(), [](const auto& a, const auto& b) {
      return a.id < b.id;
    });
    st.hyperlanes = proc::generateHyperlaneNetwork(u.seed(), laneStubs, st.hyperlaneParams);
    st.lastHyperlaneMs = std::chrono::duration<double, std::milli>(Clock::now() - tLane0).count();
  }

  st.lastBuildMs = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
  st.dirty = false;
}

void drawProceduralTradeSystemsLabWindow(ProceduralTradeSystemsLabWindowState& st,
                                        sim::Universe& universe,
                                        const sim::StarSystem* currentSystem,
                                        float /*timeSec*/,
                                        const ToastFn& toast) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(1180, 760), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Procedural Trade Systems", &st.open)) {
    ImGui::End();
    return;
  }

  bool dirty = false;

  if (ImGui::BeginTable("##trade_sys_lab_layout", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV)) {
    ImGui::TableNextColumn();

    ImGui::TextUnformatted("Inspect galaxy-coherent procedural trade profiles.");
    ImGui::TextDisabled(
        "This window is a read-only probe: it doesn't mutate the live economy. (Round 23+: profiles do shape generation.)");
    ImGui::Separator();

    // ----- Query controls -----
    dirty |= ImGui::Checkbox("Center on current system", &st.centerOnCurrentSystem);

    math::Vec3d derivedCenter = st.centerLy;
    if (st.centerOnCurrentSystem && currentSystem) {
      derivedCenter = currentSystem->stub.posLy;
      ImGui::Text("Center: %s", currentSystem->stub.name.c_str());
      ImGui::TextDisabled("(%.0f, %.0f, %.0f) ly", derivedCenter.x, derivedCenter.y, derivedCenter.z);
    } else {
      dirty |= ImGui::DragScalar("Center X (ly)", ImGuiDataType_Double, &st.centerLy.x, 10.0, nullptr, nullptr, "%.0f");
      dirty |= ImGui::DragScalar("Center Y (ly)", ImGuiDataType_Double, &st.centerLy.y, 10.0, nullptr, nullptr, "%.0f");
      dirty |= ImGui::DragScalar("Center Z (ly)", ImGuiDataType_Double, &st.centerLy.z, 10.0, nullptr, nullptr, "%.0f");
      derivedCenter = st.centerLy;
    }

    dirty |= ImGui::DragScalar("Radius (ly)", ImGuiDataType_Double, &st.radiusLy, 5.0, nullptr, nullptr, "%.0f");
    dirty |= ImGui::SliderInt("Max systems", &st.maxSystems, 16, 2048);

    ImGui::Separator();
    ImGui::TextUnformatted("Display");

    dirty |= ImGui::SliderInt("Top commodities", &st.topNCommodities, 1, (int)econ::kCommodityCount);

    const char* kSortModes[] = {"Distance", "Hub", "Wealth", "Lawlessness", "Population"};
    dirty |= ImGui::Combo("Sort", &st.sortMode, kSortModes, (int)(sizeof(kSortModes) / sizeof(kSortModes[0])));

    dirty |= ImGui::Checkbox("Show exports/imports", &st.showExportsImports);
    dirty |= ImGui::Checkbox("Show macro factors", &st.showMacroFactors);

    ImGui::Separator();
    ImGui::TextUnformatted("Macro trade routes");
    ImGui::Checkbox("Show trade routes", &st.showTradeRoutes);
    if (st.showTradeRoutes) {
      st.tradeRouteMax = std::clamp(st.tradeRouteMax, 1, 24);
      ImGui::SliderInt("Max routes", &st.tradeRouteMax, 1, 24);

      if (st.tradeRouteMaxDistanceLy < 0.0) st.tradeRouteMaxDistanceLy = 0.0;
      ImGui::DragScalar("Route max dist (ly)", ImGuiDataType_Double, &st.tradeRouteMaxDistanceLy, 10.0, nullptr, nullptr, "%.0f");

      const double minExp = 0.50;
      const double maxExp = 2.50;
      if (st.tradeRouteDistanceExponent < minExp) st.tradeRouteDistanceExponent = minExp;
      if (st.tradeRouteDistanceExponent > maxExp) st.tradeRouteDistanceExponent = maxExp;
      ImGui::SliderScalar("Distance exponent", ImGuiDataType_Double, &st.tradeRouteDistanceExponent, &minExp, &maxExp, "%.2f");

      ImGui::Separator();
      dirty |= ImGui::Checkbox("Use hyperlane routing", &st.useHyperlaneRouting);
      if (st.useHyperlaneRouting) {
        ImGui::TextDisabled("Routes are constrained by the local hyperlane graph.");

        // Hyperlane graph generation params (requires rebuild).
        if (st.hyperlaneParams.maxNeighborDistanceLy < 1.0) st.hyperlaneParams.maxNeighborDistanceLy = 1.0;
        dirty |= ImGui::DragScalar("Lane neighbor dist (ly)", ImGuiDataType_Double, &st.hyperlaneParams.maxNeighborDistanceLy,
                                  1.0, nullptr, nullptr, "%.1f");

        st.hyperlaneParams.neighborK = std::clamp(st.hyperlaneParams.neighborK, 1, 12);
        dirty |= ImGui::SliderInt("Lane neighbor K", &st.hyperlaneParams.neighborK, 1, 12);

        dirty |= ImGui::Checkbox("Lane force connected", &st.hyperlaneParams.forceConnected);

        st.hyperlaneParams.minDegree = std::clamp(st.hyperlaneParams.minDegree, 0, 6);
        dirty |= ImGui::SliderInt("Lane min degree", &st.hyperlaneParams.minDegree, 0, 6);

        const double chanceMin = 0.0;
        const double chanceMax = 0.75;
        st.hyperlaneParams.extraEdgeChance = std::clamp(st.hyperlaneParams.extraEdgeChance, chanceMin, chanceMax);
        dirty |= ImGui::SliderScalar("Lane extra edge chance", ImGuiDataType_Double, &st.hyperlaneParams.extraEdgeChance,
                                    &chanceMin, &chanceMax, "%.2f");

        if (st.hyperlaneParams.regionCellSizeLy < 0.0) st.hyperlaneParams.regionCellSizeLy = 0.0;
        dirty |= ImGui::DragScalar("Lane region cell size (ly)", ImGuiDataType_Double, &st.hyperlaneParams.regionCellSizeLy,
                                  25.0, nullptr, nullptr, "%.0f");

        // Optional hazard modulation (requires rebuild).
        ImGui::Separator();
        dirty |= ImGui::Checkbox("Lane hazards / space weather", &st.hyperlaneParams.hazards.enabled);
        if (st.hyperlaneParams.hazards.enabled) {
          ImGui::TextDisabled("Hazards bias topology + risk + bandwidth using a smooth galaxy-scale noise field.");

          // Field motion.
          dirty |= ImGui::DragScalar("Hazard time (days)", ImGuiDataType_Double, &st.hyperlaneParams.hazards.timeDays,
                                    1.0, nullptr, nullptr, "%.0f");

          if (st.hyperlaneParams.hazards.driftLyPerDay < 0.0) st.hyperlaneParams.hazards.driftLyPerDay = 0.0;
          dirty |= ImGui::DragScalar("Hazard drift (ly/day)", ImGuiDataType_Double, &st.hyperlaneParams.hazards.driftLyPerDay,
                                    0.05, nullptr, nullptr, "%.2f");

          // Field structure.
          if (st.hyperlaneParams.hazards.nebulaScaleLy < 1.0) st.hyperlaneParams.hazards.nebulaScaleLy = 1.0;
          if (st.hyperlaneParams.hazards.stormScaleLy < 1.0) st.hyperlaneParams.hazards.stormScaleLy = 1.0;
          dirty |= ImGui::DragScalar("Nebula scale (ly)", ImGuiDataType_Double, &st.hyperlaneParams.hazards.nebulaScaleLy,
                                    25.0, nullptr, nullptr, "%.0f");
          dirty |= ImGui::DragScalar("Storm scale (ly)", ImGuiDataType_Double, &st.hyperlaneParams.hazards.stormScaleLy,
                                    10.0, nullptr, nullptr, "%.0f");

          st.hyperlaneParams.hazards.edgeSamples = std::clamp(st.hyperlaneParams.hazards.edgeSamples, 1, 11);
          dirty |= ImGui::SliderInt("Edge hazard samples", &st.hyperlaneParams.hazards.edgeSamples, 1, 11);

          // Strength knobs.
          const double costMin = 0.0;
          const double costMax = 2.0;
          st.hyperlaneParams.hazards.hazardToCost = std::clamp(st.hyperlaneParams.hazards.hazardToCost, costMin, costMax);
          dirty |= ImGui::SliderScalar("Hazard->topology", ImGuiDataType_Double, &st.hyperlaneParams.hazards.hazardToCost,
                                      &costMin, &costMax, "%.2f");

          const double riskMin2 = 0.0;
          const double riskMax2 = 1.0;
          st.hyperlaneParams.hazards.hazardToRisk = std::clamp(st.hyperlaneParams.hazards.hazardToRisk, riskMin2, riskMax2);
          dirty |= ImGui::SliderScalar("Hazard->risk", ImGuiDataType_Double, &st.hyperlaneParams.hazards.hazardToRisk,
                                      &riskMin2, &riskMax2, "%.2f");

          const double bwMin2 = 0.0;
          const double bwMax2 = 0.95;
          st.hyperlaneParams.hazards.hazardToBandwidth =
              std::clamp(st.hyperlaneParams.hazards.hazardToBandwidth, bwMin2, bwMax2);
          dirty |= ImGui::SliderScalar("Hazard->bw penalty", ImGuiDataType_Double, &st.hyperlaneParams.hazards.hazardToBandwidth,
                                      &bwMin2, &bwMax2, "%.2f");
        }

        ImGui::TextDisabled("Hyperlanes: %zu edges (%.2f ms)", st.hyperlanes.edges.size(), st.lastHyperlaneMs);

        // Travel-cost shaping is immediate (no rebuild needed).
        ImGui::Separator();
        ImGui::TextUnformatted("Lane travel cost model");

        bool travelChanged = false;

        const double riskMin = 0.0;
        const double riskMax = 2.0;
        st.hyperlaneTravel.riskWeight = std::clamp(st.hyperlaneTravel.riskWeight, riskMin, riskMax);
        travelChanged |= ImGui::SliderScalar("Risk weight", ImGuiDataType_Double, &st.hyperlaneTravel.riskWeight, &riskMin, &riskMax, "%.2f");

        const double biasMin = 0.0;
        const double biasMax = 1.0;
        st.hyperlaneTravel.bandwidthBias = std::clamp(st.hyperlaneTravel.bandwidthBias, biasMin, biasMax);
        travelChanged |= ImGui::SliderScalar("Bandwidth bias", ImGuiDataType_Double, &st.hyperlaneTravel.bandwidthBias, &biasMin, &biasMax, "%.2f");

        const double minBwMin = 0.05;
        const double minBwMax = 1.0;
        st.hyperlaneTravel.minBandwidthFactor = std::clamp(st.hyperlaneTravel.minBandwidthFactor, minBwMin, minBwMax);
        travelChanged |= ImGui::SliderScalar("Min bandwidth factor", ImGuiDataType_Double, &st.hyperlaneTravel.minBandwidthFactor, &minBwMin, &minBwMax, "%.2f");

        if (travelChanged) {
          st.flowDirty = true;
        }

        // ----- Macro corridor / traffic field -----
        ImGui::Separator();
        ImGui::TextUnformatted("Interstellar traffic field");

        bool flowParamsChanged = false;
        flowParamsChanged |= ImGui::Checkbox("Compute traffic corridors", &st.showInterstellarTraffic);

        if (st.showInterstellarTraffic) {
          // Sampling controls.
          int fullLimit = (int)st.flowParams.fullPairLimit;
          fullLimit = std::clamp(fullLimit, 32, 4096);
          flowParamsChanged |= ImGui::SliderInt("Full-pairs limit", &fullLimit, 32, 2048);
          st.flowParams.fullPairLimit = (std::size_t)fullLimit;

          st.flowParams.sampleSourceCount = std::clamp(st.flowParams.sampleSourceCount, 4, 512);
          flowParamsChanged |= ImGui::SliderInt("Sample sources", &st.flowParams.sampleSourceCount, 4, 256);

          st.flowParams.samplePartnersPerSource = std::clamp(st.flowParams.samplePartnersPerSource, 0, 4096);
          flowParamsChanged |= ImGui::SliderInt("Partners per source", &st.flowParams.samplePartnersPerSource, 0, 256);

          // Model parameters.
          const double geMin = 0.50;
          const double geMax = 2.50;
          st.flowParams.gravityExponent = std::clamp(st.flowParams.gravityExponent, geMin, geMax);
          flowParamsChanged |= ImGui::SliderScalar("Gravity exponent", ImGuiDataType_Double, &st.flowParams.gravityExponent, &geMin, &geMax, "%.2f");

          const double cwMin = 0.0;
          const double cwMax = 1.0;
          st.flowParams.commodityCompatibilityWeight =
              std::clamp(st.flowParams.commodityCompatibilityWeight, cwMin, cwMax);
          flowParamsChanged |= ImGui::SliderScalar("Commodity compatibility", ImGuiDataType_Double,
                                                &st.flowParams.commodityCompatibilityWeight, &cwMin, &cwMax, "%.2f");

          ImGui::TextDisabled("Traffic: %zu nodes, %zu lanes (%.2f ms)", st.flowNet.nodes.size(), st.flowNet.edges.size(), st.lastFlowMs);
        } else {
          ImGui::TextDisabled("Traffic computation disabled.");
        }

        if (flowParamsChanged) {
          st.flowDirty = true;
        }
      }

      ImGui::Separator();
      ImGui::TextUnformatted("Macro route economics");
      ImGui::Checkbox("Show economics", &st.showRouteEconomics);
      if (st.showRouteEconomics) {
        const double spreadMin = 0.0;
        const double spreadMax = 0.35;
        st.econBidAskSpread = std::clamp(st.econBidAskSpread, spreadMin, spreadMax);
        ImGui::SliderScalar("Bid/ask spread", ImGuiDataType_Double, &st.econBidAskSpread, &spreadMin, &spreadMax, "%.2f");

        const double feeMin = 0.0;
        const double feeMax = 0.25;
        st.econBuyFeeRate = std::clamp(st.econBuyFeeRate, feeMin, feeMax);
        st.econSellFeeRate = std::clamp(st.econSellFeeRate, feeMin, feeMax);
        ImGui::SliderScalar("Buy fee", ImGuiDataType_Double, &st.econBuyFeeRate, &feeMin, &feeMax, "%.2f");
        ImGui::SliderScalar("Sell fee", ImGuiDataType_Double, &st.econSellFeeRate, &feeMin, &feeMax, "%.2f");

        if (st.econCargoKg < 0.0) st.econCargoKg = 0.0;
        ImGui::DragScalar("Cargo capacity (kg)", ImGuiDataType_Double, &st.econCargoKg, 10.0, nullptr, nullptr, "%.0f");

        ImGui::Checkbox("Show round-trip loops", &st.showTradeLoops);
        if (st.showTradeLoops) {
          st.econMaxLoops = std::clamp(st.econMaxLoops, 1, 24);
          ImGui::SliderInt("Max loops", &st.econMaxLoops, 1, 24);
        }
      }
    }


    if (ImGui::Button("Refresh")) {
      st.dirty = true;
    }

    ImGui::Separator();
    ImGui::Text("Universe seed: %llu", (unsigned long long)universe.seed());
    ImGui::Text("Systems: %zu", st.stubs.size());
    ImGui::Text("Build: %.2f ms", st.lastBuildMs);

    if (toast && ImGui::Button("Copy center coords")) {
      std::ostringstream oss;
      oss << derivedCenter.x << ", " << derivedCenter.y << ", " << derivedCenter.z;
      ImGui::SetClipboardText(oss.str().c_str());
      toast("Copied center coords to clipboard", 1.6);
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Selection");
    if (st.selectedIndex >= 0 && (std::size_t)st.selectedIndex < st.stubs.size()) {
      const auto& sel = st.stubs[(std::size_t)st.selectedIndex];
      const auto& sp = st.profiles[(std::size_t)st.selectedIndex];

      ImGui::Text("%s", sel.name.c_str());
      ImGui::SameLine();
      ImGui::TextDisabled("(%.0f ly from center)", st.distsLy[(std::size_t)st.selectedIndex]);

      if (const sim::Faction* f = findFaction(sel.factionId, universe.factions())) {
        ImGui::Text("Faction: %s", f->name.c_str());
      } else {
        ImGui::TextDisabled("Faction: Independent");
      }

      if (st.useHyperlaneRouting && st.hyperlaneParams.hazards.enabled) {
        proc::GalaxyHazardsParams haz{};
        haz.timeDays = st.hyperlaneParams.hazards.timeDays;
        haz.driftLyPerDay = st.hyperlaneParams.hazards.driftLyPerDay;
        haz.nebulaScaleLy = st.hyperlaneParams.hazards.nebulaScaleLy;
        haz.stormScaleLy = st.hyperlaneParams.hazards.stormScaleLy;
        haz.regionCellSizeLy = (st.hyperlaneParams.regionCellSizeLy > 0.0) ? st.hyperlaneParams.regionCellSizeLy : 0.0;

        const proc::GalaxyHazardSample hs = proc::sampleGalaxyHazards(universe.seed(), sel.posLy, haz);
        ImGui::Text("Hazards: %s", proc::galaxyHazardKindName(hs.kind));
        ImGui::TextDisabled("haz=%.0f%% | nebula=%.0f%% | storm=%.0f%%", hs.hazard01 * 100.0, hs.nebula01 * 100.0, hs.storm01 * 100.0);
      }

      if (st.showExportsImports) {
        const std::string ex = topCommodities(sp.exportScore, st.topNCommodities);
        const std::string im = topCommodities(sp.importScore, st.topNCommodities);
        ImGui::Text("Top exports: %s", ex.c_str());
        ImGui::Text("Top imports: %s", im.c_str());
      }

      if (st.showInterstellarTraffic && st.useHyperlaneRouting && !st.flowNet.nodes.empty()) {
        // Selected-system traffic summary.
        double traffic = 0.0;
        double traffic01 = 0.0;
        bool has = false;
        for (const auto& n : st.flowNet.nodes) {
          if (n.id == sel.id) {
            traffic = n.traffic;
            traffic01 = n.traffic01;
            has = true;
            break;
          }
        }

        ImGui::Separator();
        ImGui::TextUnformatted("Traffic / corridors");
        if (has) {
          ImGui::Text("System traffic: %.0f (%.0f%% of local max)", traffic, traffic01 * 100.0);
        } else {
          ImGui::TextDisabled("System traffic: (n/a)");
        }

        if (ImGui::CollapsingHeader("Top traffic hubs", ImGuiTreeNodeFlags_DefaultOpen)) {
          const int top = (int)std::min<std::size_t>(st.flowNet.nodes.size(), 10);
          for (int i = 0; i < top; ++i) {
            const auto& n = st.flowNet.nodes[(std::size_t)i];
            const auto* s = findStub(n.id, st.stubs);
            const char* name = s ? s->name.c_str() : "?";
            ImGui::BulletText("#%d %s (%.0f%%)", i + 1, name, n.traffic01 * 100.0);
          }
        }

        if (ImGui::CollapsingHeader("Top trade corridors", ImGuiTreeNodeFlags_DefaultOpen)) {
          const int top = (int)std::min<std::size_t>(st.flowNet.edges.size(), 10);
          for (int i = 0; i < top; ++i) {
            const auto& e = st.flowNet.edges[(std::size_t)i];
            const auto* sa = findStub(e.a, st.stubs);
            const auto* sb = findStub(e.b, st.stubs);
            const char* na = sa ? sa->name.c_str() : "?";
            const char* nb = sb ? sb->name.c_str() : "?";
            const proc::HyperlaneEdge* lane = findLaneEdge(e.a, e.b, st.hyperlanes);
            if (lane) {
              double haz01 = 0.0;
              if (st.hyperlaneParams.hazards.enabled && sa && sb) {
                proc::GalaxyHazardsParams haz{};
                haz.timeDays = st.hyperlaneParams.hazards.timeDays;
                haz.driftLyPerDay = st.hyperlaneParams.hazards.driftLyPerDay;
                haz.nebulaScaleLy = st.hyperlaneParams.hazards.nebulaScaleLy;
                haz.stormScaleLy = st.hyperlaneParams.hazards.stormScaleLy;
                haz.regionCellSizeLy = (st.hyperlaneParams.regionCellSizeLy > 0.0) ? st.hyperlaneParams.regionCellSizeLy : 0.0;
                haz01 = proc::sampleGalaxyHazardAvgOnSegment(universe.seed(), sa->posLy, sb->posLy, haz,
                                                           st.hyperlaneParams.hazards.edgeSamples);
              }

              if (st.hyperlaneParams.hazards.enabled) {
                ImGui::BulletText("#%d %s <-> %s (%.0f%% flow | d=%.0f ly, bw=%.0f%%, risk=%.0f%%, haz=%.0f%%)",
                                  i + 1, na, nb,
                                  e.flow01 * 100.0,
                                  lane->distanceLy,
                                  lane->bandwidth01 * 100.0,
                                  lane->risk01 * 100.0,
                                  haz01 * 100.0);
              } else {
                ImGui::BulletText("#%d %s <-> %s (%.0f%% flow | d=%.0f ly, bw=%.0f%%, risk=%.0f%%)",
                                  i + 1, na, nb,
                                  e.flow01 * 100.0,
                                  lane->distanceLy,
                                  lane->bandwidth01 * 100.0,
                                  lane->risk01 * 100.0);
              }
            } else {
              ImGui::BulletText("#%d %s <-> %s (%.0f%% flow)", i + 1, na, nb, e.flow01 * 100.0);
            }
          }
        }
      }

      if (st.showTradeRoutes) {
        proc::TradeRouteParams rp{};
        rp.maxRoutes = (std::size_t)std::clamp(st.tradeRouteMax, 1, 24);
        rp.maxDistanceLy = st.tradeRouteMaxDistanceLy;
        rp.distanceExponent = st.tradeRouteDistanceExponent;

        auto nameFor = [&](sim::SystemId id) -> const char* {
          for (const auto& s : st.stubs) {
            if (s.id == id) return s.name.c_str();
          }
          return "?";
        };

        if (st.showRouteEconomics) {
          proc::TradeIntelParams ip{};
          ip.bidAskSpread = st.econBidAskSpread;
          ip.buyFeeRate = st.econBuyFeeRate;
          ip.sellFeeRate = st.econSellFeeRate;
          ip.cargoCapacityKg = st.econCargoKg;
          ip.maxLoops = (std::size_t)std::clamp(st.econMaxLoops, 1, 24);

          const bool useLanes = st.useHyperlaneRouting && !st.hyperlanes.edges.empty();
          const auto rep = useLanes ? proc::buildTradeIntel(sel, sp, st.stubs, st.profiles, st.hyperlanes, rp, st.hyperlaneTravel, ip)
                                    : proc::buildTradeIntel(sel, sp, st.stubs, st.profiles, rp, ip);

          ImGui::TextDisabled("Macro route intel (profile-driven prices; not live inventories):");

          if (!rep.exports.empty()) {
            ImGui::TextUnformatted("Exports");
            for (const auto& e : rep.exports) {
              const auto code = econ::commodityCode(e.commodity);
              const std::string buy = fmtCr(e.buyAskCr);
              const std::string sell = fmtCr(e.sellBidCr);
              const std::string dU = fmtCr(e.profitPerUnitCr);
              const std::string dTrip = fmtCr(e.profitForCargoCr);
              if (useLanes) {
                ImGui::BulletText("%.*s -> %s (cost %.0f, lane %.0f, direct %.0f, hops %d, bw %.0f%%, s=%.3g, risk=%.0f%% | buy %s, sell %s, Δ/u %s, Δ/trip %s)",
                                  (int)code.size(), code.data(), nameFor(e.otherSystem),
                                  e.distanceLy, e.laneDistanceLy, e.directDistanceLy,
                                  e.laneHops, e.laneBottleneckBandwidth01 * 100.0,
                                  e.score, e.risk01 * 100.0,
                                  buy.c_str(), sell.c_str(), dU.c_str(), dTrip.c_str());
              } else {
                ImGui::BulletText("%.*s -> %s (%.0f ly, s=%.3g, risk=%.0f%% | buy %s, sell %s, Δ/u %s, Δ/trip %s)",
                                  (int)code.size(), code.data(), nameFor(e.otherSystem),
                                  e.distanceLy, e.score, e.risk01 * 100.0,
                                  buy.c_str(), sell.c_str(), dU.c_str(), dTrip.c_str());
              }
            }
          } else {
            ImGui::TextDisabled("Exports: (none)");
          }

          if (!rep.imports.empty()) {
            ImGui::TextUnformatted("Imports");
            for (const auto& e : rep.imports) {
              const auto code = econ::commodityCode(e.commodity);
              const std::string buy = fmtCr(e.buyAskCr);
              const std::string sell = fmtCr(e.sellBidCr);
              const std::string dU = fmtCr(e.profitPerUnitCr);
              const std::string dTrip = fmtCr(e.profitForCargoCr);
              if (useLanes) {
                ImGui::BulletText("%.*s <- %s (cost %.0f, lane %.0f, direct %.0f, hops %d, bw %.0f%%, s=%.3g, risk=%.0f%% | buy %s, sell %s, Δ/u %s, Δ/trip %s)",
                                  (int)code.size(), code.data(), nameFor(e.otherSystem),
                                  e.distanceLy, e.laneDistanceLy, e.directDistanceLy,
                                  e.laneHops, e.laneBottleneckBandwidth01 * 100.0,
                                  e.score, e.risk01 * 100.0,
                                  buy.c_str(), sell.c_str(), dU.c_str(), dTrip.c_str());
              } else {
                ImGui::BulletText("%.*s <- %s (%.0f ly, s=%.3g, risk=%.0f%% | buy %s, sell %s, Δ/u %s, Δ/trip %s)",
                                  (int)code.size(), code.data(), nameFor(e.otherSystem),
                                  e.distanceLy, e.score, e.risk01 * 100.0,
                                  buy.c_str(), sell.c_str(), dU.c_str(), dTrip.c_str());
              }
            }
          } else {
            ImGui::TextDisabled("Imports: (none)");
          }

          if (st.showTradeLoops) {
            ImGui::Separator();
            ImGui::TextUnformatted("Round-trip loops");
            if (!rep.loops.empty()) {
              for (const auto& l : rep.loops) {
                const auto outCode = econ::commodityCode(l.outLeg.commodity);
                const auto backCode = econ::commodityCode(l.backLeg.commodity);
                const std::string total = fmtCr(l.roundTripProfitCr);
                const std::string outP = fmtCr(l.outLeg.profitForCargoCr);
                const std::string backP = fmtCr(l.backLeg.profitForCargoCr);
                if (useLanes) {
                  const double laneRt = l.outLeg.laneDistanceLy + l.backLeg.laneDistanceLy;
                  const double directRt = l.outLeg.directDistanceLy + l.backLeg.directDistanceLy;
                  const int hopsRt = l.outLeg.laneHops + l.backLeg.laneHops;
                  const double bwRt = 100.0 * std::min(l.outLeg.laneBottleneckBandwidth01, l.backLeg.laneBottleneckBandwidth01);
                  ImGui::BulletText("%s: out %.*s (%s), back %.*s (%s) => %s (cost %.0f, lane %.0f, direct %.0f, hops %d, bw %.0f%%, risk %.0f%%)",
                                    nameFor(l.otherSystem),
                                    (int)outCode.size(), outCode.data(), outP.c_str(),
                                    (int)backCode.size(), backCode.data(), backP.c_str(),
                                    total.c_str(), l.roundTripDistanceLy, laneRt, directRt, hopsRt, bwRt, l.avgRisk01 * 100.0);
                } else {
                  ImGui::BulletText("%s: out %.*s (%s), back %.*s (%s) => %s (%.0f ly, risk=%.0f%%)",
                                    nameFor(l.otherSystem),
                                    (int)outCode.size(), outCode.data(), outP.c_str(),
                                    (int)backCode.size(), backCode.data(), backP.c_str(),
                                    total.c_str(), l.roundTripDistanceLy, l.avgRisk01 * 100.0);
                }
              }
            } else {
              ImGui::TextDisabled("(no loops found in current route set)");
            }
          }
        } else {
          const bool useLanes = st.useHyperlaneRouting && !st.hyperlanes.edges.empty();
          const auto routes = useLanes ? proc::suggestTradeRoutes(sel, sp, st.stubs, st.profiles, st.hyperlanes, rp, st.hyperlaneTravel)
                                       : proc::suggestTradeRoutes(sel, sp, st.stubs, st.profiles, rp);
          ImGui::TextDisabled("Macro route signal (profile-driven, not live prices):");

          if (!routes.exports.empty()) {
            ImGui::TextUnformatted("Exports");
            for (const auto& r : routes.exports) {
              const auto code = econ::commodityCode(r.commodity);
              if (useLanes) {
                ImGui::BulletText("%.*s -> %s (cost %.0f, lane %.0f, direct %.0f, hops %d, bw %.0f%%, s=%.3g, risk=%.0f%%)",
                                  (int)code.size(), code.data(), nameFor(r.otherSystem),
                                  r.distanceLy, r.laneDistanceLy, r.directDistanceLy,
                                  r.laneHops, r.laneBottleneckBandwidth01 * 100.0,
                                  r.score, r.risk01 * 100.0);
              } else {
                ImGui::BulletText("%.*s -> %s (%.0f ly, s=%.3g, risk=%.0f%%)",
                                  (int)code.size(), code.data(), nameFor(r.otherSystem),
                                  r.distanceLy, r.score, r.risk01 * 100.0);
              }
            }
          } else {
            ImGui::TextDisabled("Exports: (none)");
          }

          if (!routes.imports.empty()) {
            ImGui::TextUnformatted("Imports");
            for (const auto& r : routes.imports) {
              const auto code = econ::commodityCode(r.commodity);
              if (useLanes) {
                ImGui::BulletText("%.*s <- %s (cost %.0f, lane %.0f, direct %.0f, hops %d, bw %.0f%%, s=%.3g, risk=%.0f%%)",
                                  (int)code.size(), code.data(), nameFor(r.otherSystem),
                                  r.distanceLy, r.laneDistanceLy, r.directDistanceLy,
                                  r.laneHops, r.laneBottleneckBandwidth01 * 100.0,
                                  r.score, r.risk01 * 100.0);
              } else {
                ImGui::BulletText("%.*s <- %s (%.0f ly, s=%.3g, risk=%.0f%%)",
                                  (int)code.size(), code.data(), nameFor(r.otherSystem),
                                  r.distanceLy, r.score, r.risk01 * 100.0);
              }
            }
          } else {
            ImGui::TextDisabled("Imports: (none)");
          }
        }
      } else {
        ImGui::TextDisabled("Trade routes disabled (enable above).");
      }
    } else {
      ImGui::TextDisabled("Click a system row to inspect its macro trade routes.");
    }

    ImGui::TableNextColumn();

    // ----- Results table -----
    if (dirty) st.dirty = true;

    const bool seedChanged = (st.lastUniverseSeed != universe.seed());
    if (seedChanged) {
      st.dirty = true;
    }

    if (!st.centerOnCurrentSystem) {
      // st.centerLy is authoritative.
    } else if (currentSystem) {
      // Center depends on current system.
    } else {
      // No system loaded; keep previous center.
    }

    const bool centerChanged = (st.lastCenterLy.x != derivedCenter.x) || (st.lastCenterLy.y != derivedCenter.y) ||
                               (st.lastCenterLy.z != derivedCenter.z);
    const bool radiusChanged = (st.lastRadiusLy != st.radiusLy);
    const bool maxChanged = (st.lastMaxSystems != st.maxSystems);

    if (st.dirty || centerChanged || radiusChanged || maxChanged) {
      st.lastUniverseSeed = universe.seed();
      st.lastCenterLy = derivedCenter;
      st.lastRadiusLy = st.radiusLy;
      st.lastMaxSystems = st.maxSystems;
      rebuild(st, universe, derivedCenter);
    }

    // Recompute the macro traffic field when requested.
    if (st.showInterstellarTraffic && st.useHyperlaneRouting && !st.hyperlanes.edges.empty() && !st.stubs.empty()) {
      st.flowParams.travelParams = st.hyperlaneTravel;
      if (st.flowDirty) {
        const auto tFlow0 = Clock::now();
        st.flowNet = proc::computeTradeFlow(universe.seed(), st.stubs, st.profiles, st.hyperlanes, st.flowParams);
        st.lastFlowMs = std::chrono::duration<double, std::milli>(Clock::now() - tFlow0).count();
        st.flowDirty = false;
      }
    }

    ImGui::TextUnformatted("Nearby systems");
    ImGui::SameLine();
    ImGui::TextDisabled("(sorted by %s)", sortModeName(st.sortMode));

    ImGuiTableFlags flags = ImGuiTableFlags_ScrollY | ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_Resizable;
    if (ImGui::BeginTable("##trade_sys_table", st.showExportsImports ? 11 : 9, flags, ImVec2(0.0f, ImGui::GetContentRegionAvail().y))) {
      ImGui::TableSetupScrollFreeze(0, 1);
      ImGui::TableSetupColumn("Name", ImGuiTableColumnFlags_WidthFixed, 180.0f);
      ImGui::TableSetupColumn("Dist", ImGuiTableColumnFlags_WidthFixed, 60.0f);
      ImGui::TableSetupColumn("Star", ImGuiTableColumnFlags_WidthFixed, 45.0f);
      ImGui::TableSetupColumn("Faction", ImGuiTableColumnFlags_WidthFixed, 120.0f);
      ImGui::TableSetupColumn("Stations", ImGuiTableColumnFlags_WidthFixed, 65.0f);

      if (st.showMacroFactors) {
        ImGui::TableSetupColumn("Hub", ImGuiTableColumnFlags_WidthFixed, 55.0f);
        ImGui::TableSetupColumn("Pop", ImGuiTableColumnFlags_WidthFixed, 55.0f);
        ImGui::TableSetupColumn("Wealth", ImGuiTableColumnFlags_WidthFixed, 60.0f);
        ImGui::TableSetupColumn("Law", ImGuiTableColumnFlags_WidthFixed, 55.0f);
      } else {
        ImGui::TableSetupColumn("Hub", ImGuiTableColumnFlags_DefaultHide, 55.0f);
        ImGui::TableSetupColumn("Pop", ImGuiTableColumnFlags_DefaultHide, 55.0f);
        ImGui::TableSetupColumn("Wealth", ImGuiTableColumnFlags_DefaultHide, 60.0f);
        ImGui::TableSetupColumn("Law", ImGuiTableColumnFlags_DefaultHide, 55.0f);
      }

      if (st.showExportsImports) {
        ImGui::TableSetupColumn("Exports", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableSetupColumn("Imports", ImGuiTableColumnFlags_WidthFixed, 120.0f);
      }

      ImGui::TableHeadersRow();

      const auto& factions = universe.factions();

      for (std::size_t i = 0; i < st.stubs.size(); ++i) {
        const auto& s = st.stubs[i];
        const auto& p = st.profiles[i];

        ImGui::TableNextRow();

        ImGui::TableSetColumnIndex(0);
        ImGui::PushID((int)i);
        const bool selected = (st.selectedIndex == (int)i);
        if (ImGui::Selectable(s.name.c_str(), selected, ImGuiSelectableFlags_SpanAllColumns)) {
          st.selectedIndex = (int)i;
        }
        ImGui::PopID();

        ImGui::TableSetColumnIndex(1);
        ImGui::Text("%.0f", st.distsLy[i]);

        ImGui::TableSetColumnIndex(2);
        ImGui::TextUnformatted([&]() {
          switch (s.primaryClass) {
            case sim::StarClass::O: return "O";
            case sim::StarClass::B: return "B";
            case sim::StarClass::A: return "A";
            case sim::StarClass::F: return "F";
            case sim::StarClass::G: return "G";
            case sim::StarClass::K: return "K";
            case sim::StarClass::M: return "M";
            default: return "?";
          }
        }());

        ImGui::TableSetColumnIndex(3);
        if (const sim::Faction* f = findFaction(s.factionId, factions)) {
          ImGui::TextUnformatted(f->name.c_str());
        } else {
          ImGui::TextDisabled("Independent");
        }

        ImGui::TableSetColumnIndex(4);
        ImGui::Text("%d", s.stationCount);

        ImGui::TableSetColumnIndex(5);
        ImGui::Text("%.2f", p.hub);

        ImGui::TableSetColumnIndex(6);
        ImGui::Text("%.2f", p.population);

        ImGui::TableSetColumnIndex(7);
        ImGui::Text("%.2f", p.wealth);

        ImGui::TableSetColumnIndex(8);
        ImGui::Text("%.2f", p.lawlessness);

        int col = 9;
        if (st.showExportsImports) {
          ImGui::TableSetColumnIndex(col++);
          const std::string ex = topCommodities(p.exportScore, st.topNCommodities);
          ImGui::TextUnformatted(ex.c_str());

          ImGui::TableSetColumnIndex(col++);
          const std::string im = topCommodities(p.importScore, st.topNCommodities);
          ImGui::TextUnformatted(im.c_str());
        }
      }

      ImGui::EndTable();
    }

    ImGui::EndTable();
  }

  ImGui::End();
}

} // namespace stellar::game
