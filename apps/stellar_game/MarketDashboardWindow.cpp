#include "MarketDashboardWindow.h"

#include "stellar/econ/Market.h"
#include "stellar/math/Math.h"
#include "stellar/ui/FuzzySearch.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>

namespace stellar::game {

static constexpr std::size_t cidx(stellar::econ::CommodityId id) {
  return static_cast<std::size_t>(id);
}

static const char* stationTypeName(stellar::econ::StationType t) {
  using stellar::econ::StationType;
  switch (t) {
    case StationType::Outpost:       return "Outpost";
    case StationType::Agricultural:  return "Agricultural";
    case StationType::Mining:        return "Mining";
    case StationType::Refinery:      return "Refinery";
    case StationType::Industrial:    return "Industrial";
    case StationType::Research:      return "Research";
    case StationType::TradeHub:      return "Trade Hub";
    case StationType::Shipyard:      return "Shipyard";
    case StationType::Count:         break;
  }
  return "?";
}

static stellar::econ::CommodityId clampCommodity(int idx) {
  const int maxIdx = (int)stellar::econ::kCommodityCount - 1;
  idx = std::clamp(idx, 0, std::max(0, maxIdx));
  return (stellar::econ::CommodityId)idx;
}

static int cacheStampFor(double timeDays) {
  // Align with station economy sampling cadence (default 0.25 days).
  if (!std::isfinite(timeDays)) timeDays = 0.0;
  return (int)std::floor(timeDays * 4.0);
}

static double clamp01(double x) {
  if (!std::isfinite(x)) return 0.0;
  return std::clamp(x, 0.0, 1.0);
}

static void rebuildRows(MarketDashboardWindowState& st, const MarketDashboardContext& ctx) {
  st.rows.clear();

  if (!ctx.currentSystem) return;

  const stellar::econ::CommodityId commodity = clampCommodity(st.commodityIndex);
  const double bidAskSpread = 0.10;
  const stellar::math::Vec3d origin = ctx.currentSystem->stub.posLy;

  auto addStation = [&](const stellar::sim::StarSystem& sys,
                        const stellar::sim::Station& station,
                        double distLy) {
    auto& econState = ctx.universe.stationEconomy(station, ctx.timeDays);
    const auto q = stellar::econ::quote(econState, station.economyModel, commodity, bidAskSpread);
    const double fee = ctx.effectiveFeeRate ? clamp01(ctx.effectiveFeeRate(station)) : clamp01(station.feeRate);

    MarketDashboardRow row;
    row.systemId = sys.stub.id;
    row.stationId = station.id;
    row.systemName = sys.stub.name;
    row.stationName = station.name;
    row.stationType = station.type;
    row.distanceLy = distLy;
    row.mid = q.mid;
    row.askEff = q.ask * (1.0 + fee);
    row.bidEff = q.bid * (1.0 - fee);
    row.inventory = q.inventory;
    row.capacity = std::max(0.0, station.economyModel.capacity[cidx(commodity)]);
    row.feeRate = fee;

    // Trend (best-effort) from last 2 samples.
    const auto& hist = econState.history[cidx(commodity)];
    if (hist.size() >= 2) {
      const auto& a = hist[hist.size() - 2];
      const auto& b = hist[hist.size() - 1];
      const double dt = std::max(1e-9, b.day - a.day);
      row.trendPerDay = (b.price - a.price) / dt;
      row.lastSampleDay = b.day;
    } else if (!hist.empty()) {
      row.trendPerDay = 0.0;
      row.lastSampleDay = hist.back().day;
    } else {
      row.trendPerDay = 0.0;
      row.lastSampleDay = econState.lastUpdateDay;
    }

    st.rows.push_back(std::move(row));
  };

  if (!st.scopeNearby) {
    // Current system only.
    const auto& sys = *ctx.currentSystem;
    for (const auto& station : sys.stations) {
      addStation(sys, station, 0.0);
    }
  } else {
    const auto stubs = ctx.universe.queryNearby(origin, std::max(0.0, st.radiusLy), (std::size_t)std::max(1, st.maxSystems));
    for (const auto& stub : stubs) {
      if (!st.includeCurrentSystem && stub.id == ctx.currentSystem->stub.id) continue;

      const auto& sys = ctx.universe.getSystem(stub.id, &stub);
      if (sys.stations.empty()) continue;

      const double distLy = stellar::math::length(stub.posLy - origin);
      for (const auto& station : sys.stations) {
        addStation(sys, station, distLy);
      }
    }
  }

  // Sorting.
  const int mode = st.sortMode;
  std::stable_sort(st.rows.begin(), st.rows.end(), [&](const MarketDashboardRow& a, const MarketDashboardRow& b) {
    switch (mode) {
      case 0: // Best Buy (ask asc)
        if (a.askEff != b.askEff) return a.askEff < b.askEff;
        break;
      case 1: // Best Sell (bid desc)
        if (a.bidEff != b.bidEff) return a.bidEff > b.bidEff;
        break;
      case 2: // Distance
        if (a.distanceLy != b.distanceLy) return a.distanceLy < b.distanceLy;
        break;
      case 3: // Name
        if (a.systemName != b.systemName) return a.systemName < b.systemName;
        if (a.stationName != b.stationName) return a.stationName < b.stationName;
        break;
      default:
        break;
    }
    // Deterministic tie-break.
    if (a.systemId != b.systemId) return a.systemId < b.systemId;
    return a.stationId < b.stationId;
  });
}

static bool matchesFilter(const MarketDashboardRow& r, const char* filter, int* outScore = nullptr) {
  if (!filter || filter[0] == 0) {
    if (outScore) *outScore = 0;
    return true;
  }

  const std::string hay = r.systemName + " " + r.stationName;
  const auto m = stellar::ui::fuzzyMatch(filter, hay);
  const int score = m.score;
  if (outScore) *outScore = score;
  return score >= 0;
}

static void drawTrend(double trendPerDay) {
  if (!std::isfinite(trendPerDay)) {
    ImGui::TextDisabled("?");
    return;
  }
  if (std::abs(trendPerDay) < 1e-6) {
    ImGui::TextDisabled("~");
    return;
  }
  const char* arrow = (trendPerDay > 0.0) ? "^" : "v";
  ImGui::Text("%s %.2f", arrow, std::abs(trendPerDay));
}

static const MarketDashboardRow* findRow(const MarketDashboardWindowState& st,
                                        stellar::sim::SystemId sysId,
                                        stellar::sim::StationId stId) {
  for (const auto& r : st.rows) {
    if (r.systemId == sysId && r.stationId == stId) return &r;
  }
  return nullptr;
}

static const stellar::sim::Station* findStation(const stellar::sim::StarSystem& sys, stellar::sim::StationId id) {
  for (const auto& st : sys.stations) {
    if (st.id == id) return &st;
  }
  return nullptr;
}

static void drawHistoryPanel(MarketDashboardWindowState& st, const MarketDashboardContext& ctx) {
  if (!st.showHistoryPanel) return;
  if (!ctx.currentSystem) return;
  if (st.selectedStation == 0 || st.selectedSystem == 0) return;

  const auto* row = findRow(st, st.selectedSystem, st.selectedStation);
  if (!row) return;

  const stellar::sim::StarSystem* sysPtr = nullptr;
  if (ctx.currentSystem && ctx.currentSystem->stub.id == st.selectedSystem) {
    sysPtr = ctx.currentSystem;
  } else {
    sysPtr = &ctx.universe.getSystem(st.selectedSystem);
  }
  if (!sysPtr) return;

  const auto* station = findStation(*sysPtr, st.selectedStation);
  if (!station) return;

  const stellar::econ::CommodityId commodity = clampCommodity(st.commodityIndex);
  auto& econState = ctx.universe.stationEconomy(*station, ctx.timeDays);
  const auto& hist = econState.history[cidx(commodity)];

  ImGui::SeparatorText("Price history");
  ImGui::Text("%s / %s", row->systemName.c_str(), row->stationName.c_str());
  ImGui::TextDisabled("%s, fee %.1f%%", stationTypeName(row->stationType), row->feeRate * 100.0);

  // Build plot array.
  std::vector<float> ys;
  ys.reserve(hist.size());
  float ymin = std::numeric_limits<float>::infinity();
  float ymax = -std::numeric_limits<float>::infinity();
  for (const auto& p : hist) {
    const float v = (float)p.price;
    ys.push_back(v);
    ymin = std::min(ymin, v);
    ymax = std::max(ymax, v);
  }
  if (!ys.empty()) {
    // Pad ranges a bit.
    const float pad = (ymax - ymin) * 0.05f;
    const float a = std::isfinite(pad) ? pad : 0.0f;
    ymin = std::max(0.0f, ymin - a);
    ymax = ymax + a;
  } else {
    ymin = 0.0f;
    ymax = 1.0f;
  }

  const auto q = stellar::econ::quote(econState, station->economyModel, commodity, 0.10);
  const std::string overlay = std::string("mid ") + std::to_string((int)std::round(q.mid));

  ImGui::PlotLines("##mid_history",
                   ys.empty() ? nullptr : ys.data(),
                   (int)ys.size(),
                   0,
                   overlay.c_str(),
                   ymin,
                   ymax,
                   ImVec2(0, 120));

  const double askEff = q.ask * (1.0 + row->feeRate);
  const double bidEff = q.bid * (1.0 - row->feeRate);
  ImGui::Text("Ask (buy): %.1f", askEff);
  ImGui::Text("Bid (sell): %.1f", bidEff);
  ImGui::Text("Inventory: %.0f / %.0f", row->inventory, row->capacity);
  if (!hist.empty()) {
    const double age = std::max(0.0, ctx.timeDays - hist.back().day);
    ImGui::TextDisabled("Last sample: %.2f days ago", age);
  }
}

void drawMarketDashboardWindow(MarketDashboardWindowState& st, const MarketDashboardContext& ctx) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(900, 650), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Market Dashboard", &st.open)) {
    ImGui::End();
    return;
  }

  if (!ctx.currentSystem) {
    ImGui::TextDisabled("No current system.");
    ImGui::End();
    return;
  }

  // Commodity selector.
  {
    const stellar::econ::CommodityId commodity = clampCommodity(st.commodityIndex);
    const auto& def = stellar::econ::commodityDef(commodity);
    std::string label = std::string(def.code) + " - " + def.name;

    if (ImGui::BeginCombo("Commodity", label.c_str())) {
      for (std::size_t i = 0; i < stellar::econ::kCommodityCount; ++i) {
        const auto id = (stellar::econ::CommodityId)i;
        const auto& d = stellar::econ::commodityDef(id);
        const bool sel = (int)i == st.commodityIndex;
        std::string item = std::string(d.code) + " - " + d.name;
        if (ImGui::Selectable(item.c_str(), sel)) {
          st.commodityIndex = (int)i;
        }
        if (sel) ImGui::SetItemDefaultFocus();
      }
      ImGui::EndCombo();
    }
  }

  // Scope controls.
  {
    bool cur = !st.scopeNearby;
    if (ImGui::RadioButton("Current system", cur)) st.scopeNearby = false;
    ImGui::SameLine();
    bool near = st.scopeNearby;
    if (ImGui::RadioButton("Nearby systems", near)) st.scopeNearby = true;

    if (st.scopeNearby) {
      const double minLy = 25.0;
      const double maxLy = 1500.0;
      ImGui::SliderScalar("Radius (ly)", ImGuiDataType_Double, &st.radiusLy, &minLy, &maxLy, "%.0f", ImGuiSliderFlags_Logarithmic);
      ImGui::SliderInt("Max systems", &st.maxSystems, 8, 1024);
      ImGui::Checkbox("Include current system", &st.includeCurrentSystem);
    }
  }

  // Sort/filter.
  {
    const char* sortItems[] = {"Best Buy", "Best Sell", "Distance", "Name"};
    ImGui::Combo("Sort", &st.sortMode, sortItems, IM_ARRAYSIZE(sortItems));
    ImGui::SameLine();
    ImGui::Checkbox("History", &st.showHistoryPanel);
    ImGui::InputTextWithHint("Filter", "station or system", st.filter, IM_ARRAYSIZE(st.filter));
  }

  // Refresh cache.
  const int stamp = cacheStampFor(ctx.timeDays);
  const bool needsRebuild =
    st.cacheStamp != stamp ||
    st.cacheCommodity != st.commodityIndex ||
    st.cacheScopeNearby != st.scopeNearby ||
    st.cacheRadiusLy != st.radiusLy ||
    st.cacheMaxSystems != st.maxSystems ||
    st.cacheIncludeCurrentSystem != st.includeCurrentSystem;

  if (needsRebuild) {
    st.cacheStamp = stamp;
    st.cacheCommodity = st.commodityIndex;
    st.cacheScopeNearby = st.scopeNearby;
    st.cacheRadiusLy = st.radiusLy;
    st.cacheMaxSystems = st.maxSystems;
    st.cacheIncludeCurrentSystem = st.includeCurrentSystem;
    rebuildRows(st, ctx);

    // Keep selection stable if possible.
    if (st.selectedStation != 0 && st.selectedSystem != 0) {
      if (!findRow(st, st.selectedSystem, st.selectedStation)) {
        st.selectedStation = 0;
        st.selectedSystem = 0;
      }
    }
  }

  // Compute best buy/sell across currently matching filter.
  const MarketDashboardRow* bestBuy = nullptr;
  const MarketDashboardRow* bestSell = nullptr;
  int bestBuyScore = std::numeric_limits<int>::min();
  int bestSellScore = std::numeric_limits<int>::min();

  for (const auto& r : st.rows) {
    int score = 0;
    if (!matchesFilter(r, st.filter, &score)) continue;

    if (!bestBuy || r.askEff < bestBuy->askEff || (r.askEff == bestBuy->askEff && score > bestBuyScore)) {
      bestBuy = &r;
      bestBuyScore = score;
    }
    if (!bestSell || r.bidEff > bestSell->bidEff || (r.bidEff == bestSell->bidEff && score > bestSellScore)) {
      bestSell = &r;
      bestSellScore = score;
    }
  }

  if (bestBuy && bestSell) {
    const double profitPerUnit = bestSell->bidEff - bestBuy->askEff;
    ImGui::SeparatorText("Quick arbitrage");
    ImGui::Text("Best buy:  %s / %s  (%.1f)", bestBuy->systemName.c_str(), bestBuy->stationName.c_str(), bestBuy->askEff);
    ImGui::Text("Best sell: %s / %s  (%.1f)", bestSell->systemName.c_str(), bestSell->stationName.c_str(), bestSell->bidEff);
    if (profitPerUnit > 0.0) {
      ImGui::Text("Profit/unit: %.1f", profitPerUnit);

      const auto commodity = clampCommodity(st.commodityIndex);
      const auto& def = stellar::econ::commodityDef(commodity);
      const double emptyHoldKg = std::max(0.0, ctx.cargoCapacityKg);
      const double freeHoldKg = std::max(0.0, ctx.cargoCapacityKg - ctx.cargoUsedKg);
      const double unitsEmptyHold = (def.massKg > 1e-9) ? (emptyHoldKg / def.massKg) : 0.0;
      const double unitsFreeHold = (def.massKg > 1e-9) ? (freeHoldKg / def.massKg) : 0.0;

      const double unitsByCredits = (bestBuy->askEff > 1e-9) ? (ctx.playerCreditsCr / bestBuy->askEff) : 0.0;
      const double unitsCap = std::max(0.0, bestBuy->inventory);

      const double unitsEmpty = std::max(0.0, std::min({unitsEmptyHold, unitsByCredits, unitsCap}));
      const double unitsFree = std::max(0.0, std::min({unitsFreeHold, unitsByCredits, unitsCap}));

      ImGui::TextDisabled("Est. units (empty hold): %.0f  | profit: %.0f",
                          unitsEmpty,
                          unitsEmpty * profitPerUnit);
      ImGui::TextDisabled("Est. units (free space): %.0f  | profit: %.0f",
                          unitsFree,
                          unitsFree * profitPerUnit);
    } else {
      ImGui::TextDisabled("No profitable arbitrage in current scan.");
    }
  }

  ImGui::SeparatorText("Stations");

  // Results table.
  const bool showSystemCol = st.scopeNearby;
  const bool showDistCol = st.scopeNearby;
  const int cols = showSystemCol ? 9 : 7;
  const float tableH = st.showHistoryPanel ? 260.0f : 420.0f;

  ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_Resizable |
                          ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_ScrollY;

  if (ImGui::BeginTable("##market_table", cols, flags, ImVec2(0, tableH))) {
    int col = 0;
    if (showSystemCol) ImGui::TableSetupColumn("System", ImGuiTableColumnFlags_WidthStretch);
    ImGui::TableSetupColumn("Station", ImGuiTableColumnFlags_WidthStretch);
    ImGui::TableSetupColumn("Type", ImGuiTableColumnFlags_WidthFixed);
    if (showDistCol) ImGui::TableSetupColumn("Dist", ImGuiTableColumnFlags_WidthFixed);
    ImGui::TableSetupColumn("Ask", ImGuiTableColumnFlags_WidthFixed);
    ImGui::TableSetupColumn("Bid", ImGuiTableColumnFlags_WidthFixed);
    ImGui::TableSetupColumn("Inv", ImGuiTableColumnFlags_WidthFixed);
    ImGui::TableSetupColumn("Trend", ImGuiTableColumnFlags_WidthFixed);
    ImGui::TableSetupColumn("Actions", ImGuiTableColumnFlags_WidthFixed);
    ImGui::TableHeadersRow();

    for (const auto& r : st.rows) {
      if (!matchesFilter(r, st.filter)) continue;

      ImGui::TableNextRow();
      col = 0;

      const bool selected = (r.systemId == st.selectedSystem && r.stationId == st.selectedStation);

      if (showSystemCol) {
        ImGui::TableSetColumnIndex(col++);
        ImGui::TextUnformatted(r.systemName.c_str());
      }

      ImGui::TableSetColumnIndex(col++);
      if (ImGui::Selectable(r.stationName.c_str(), selected,
                            ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap)) {
        st.selectedSystem = r.systemId;
        st.selectedStation = r.stationId;
      }

      ImGui::TableSetColumnIndex(col++);
      ImGui::TextUnformatted(stationTypeName(r.stationType));

      if (showDistCol) {
        ImGui::TableSetColumnIndex(col++);
        ImGui::Text("%.0f", r.distanceLy);
      }

      ImGui::TableSetColumnIndex(col++);
      ImGui::Text("%.1f", r.askEff);

      ImGui::TableSetColumnIndex(col++);
      ImGui::Text("%.1f", r.bidEff);

      ImGui::TableSetColumnIndex(col++);
      ImGui::Text("%.0f", r.inventory);

      ImGui::TableSetColumnIndex(col++);
      drawTrend(r.trendPerDay);

      ImGui::TableSetColumnIndex(col++);
      ImGui::PushID((int)r.stationId);
      if (ctx.routeToStation && ImGui::SmallButton("Route")) {
        ctx.routeToStation(r.systemId, r.stationId);
      }
      if (ctx.targetStation && ctx.currentSystem && r.systemId == ctx.currentSystem->stub.id) {
        ImGui::SameLine();
        if (ImGui::SmallButton("Target")) {
          ctx.targetStation(r.stationId);
        }
      }
      ImGui::PopID();
    }

    ImGui::EndTable();
  }

  drawHistoryPanel(st, ctx);

  ImGui::End();
}

} // namespace stellar::game
