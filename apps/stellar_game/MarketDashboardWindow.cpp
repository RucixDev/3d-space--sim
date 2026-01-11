#include "MarketDashboardWindow.h"

#include "stellar/econ/Market.h"
#include "stellar/sim/MarketAnalysis.h"
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

    // Trend/volatility over a configurable window.
    const auto& hist = econState.history[cidx(commodity)];
    const auto stats = stellar::sim::analyzePriceHistory(hist, ctx.timeDays, st.trendWindowDays);
    if (stats.valid) {
      row.trendValid = true;
      row.trendPerDay = stats.slopePerDay;
      row.trendR2 = stats.r2;
      row.changePct = stats.pctChange;
      row.volatilityPerDay = stats.volatilityPerDay;
      row.lastSampleDay = stats.lastDay;
    } else {
      row.trendValid = false;
      row.trendPerDay = 0.0;
      row.trendR2 = 0.0;
      row.changePct = 0.0;
      row.volatilityPerDay = 0.0;
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

static void drawTrendCell(const MarketDashboardRow& r, double windowDays, double nowDay) {
  if (!r.trendValid || !std::isfinite(r.changePct)) {
    ImGui::TextDisabled("?");
    return;
  }

  const double pct = r.changePct;
  const double absPct = std::abs(pct);

  const char* arrow = "~";
  if (absPct >= 0.05) arrow = (pct > 0.0) ? "^" : "v";

  ImGui::Text("%s %.1f%%", arrow, absPct);

  if (ImGui::IsItemHovered()) {
    ImGui::BeginTooltip();
    ImGui::Text("Trend window: %.0f d", std::max(0.0, windowDays));
    if (std::isfinite(r.trendPerDay)) {
      ImGui::Text("Slope: %.2f cr/u/day", r.trendPerDay);
    }
    if (std::isfinite(r.trendR2)) {
      ImGui::Text("Fit R^2: %.2f", r.trendR2);
    }
    if (std::isfinite(r.volatilityPerDay)) {
      ImGui::Text("Volatility: %.2f%%/day", r.volatilityPerDay * 100.0);
    }
    if (std::isfinite(nowDay) && std::isfinite(r.lastSampleDay)) {
      const double age = std::max(0.0, nowDay - r.lastSampleDay);
      ImGui::TextDisabled("Last sample: %.2f days ago", age);
    }
    ImGui::EndTooltip();
  }
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


struct PlotPt {
  double day{0.0};
  double price{0.0};
};

static void drawPriceHistoryPlot(const std::vector<stellar::econ::PricePoint>& hist,
                                double nowDay,
                                double windowDays,
                                const char* id,
                                double referencePrice) {
  const double minDay = (windowDays > 0.0) ? (nowDay - windowDays) : -std::numeric_limits<double>::infinity();

  std::vector<PlotPt> pts;
  pts.reserve(hist.size());
  for (const auto& p : hist) {
    if (!std::isfinite(p.day) || !std::isfinite(p.price)) continue;
    if (p.day < minDay) continue;
    if (p.day > nowDay) continue;
    pts.push_back(PlotPt{p.day, p.price});
  }

  if (pts.size() < 2) {
    ImGui::TextDisabled("No history yet (time needs to advance).");
    return;
  }

  double minP = std::numeric_limits<double>::infinity();
  double maxP = -std::numeric_limits<double>::infinity();
  for (const auto& p : pts) {
    minP = std::min(minP, p.price);
    maxP = std::max(maxP, p.price);
  }

  if (!std::isfinite(minP) || !std::isfinite(maxP)) {
    ImGui::TextDisabled("History contains non-finite values.");
    return;
  }

  // Pad range a bit so flat-ish curves are visible.
  double range = maxP - minP;
  if (range < 1e-9) range = std::max(1.0, std::abs(maxP) * 0.05);
  const double pad = range * 0.05;
  minP = std::max(0.0, minP - pad);
  maxP = maxP + pad;
  range = std::max(1e-9, maxP - minP);

  const double t0 = pts.front().day;
  double t1 = pts.back().day;
  double trange = t1 - t0;
  if (trange < 1e-9) {
    t1 = t0 + 1.0;
    trange = 1.0;
  }

  const float h = 140.0f;
  ImVec2 size(ImGui::GetContentRegionAvail().x, h);
  if (size.x < 64.0f) size.x = 64.0f;

  const ImVec2 p0 = ImGui::GetCursorScreenPos();
  const ImVec2 p1 = ImVec2(p0.x + size.x, p0.y + size.y);

  ImGui::InvisibleButton(id, size);
  const bool hovered = ImGui::IsItemHovered();

  ImDrawList* draw = ImGui::GetWindowDrawList();
  const ImU32 colBg = ImGui::GetColorU32(ImGuiCol_FrameBg);
  const ImU32 colBorder = ImGui::GetColorU32(ImGuiCol_Border);
  const ImU32 colLine = ImGui::GetColorU32(ImGuiCol_PlotLines);
  const ImU32 colLineHover = ImGui::GetColorU32(ImGuiCol_PlotLinesHovered);

  draw->AddRectFilled(p0, p1, colBg, 4.0f);
  draw->AddRect(p0, p1, colBorder, 4.0f);

  // Draw a reference line at the current (mid) price.
  if (std::isfinite(referencePrice)) {
    const double yn = (referencePrice - minP) / range;
    const float y = p0.y + (float)((1.0 - yn) * (double)size.y);
    if (y >= p0.y && y <= p1.y) {
      draw->AddLine(ImVec2(p0.x, y), ImVec2(p1.x, y), colBorder, 1.0f);
    }
  }

  std::vector<ImVec2> poly;
  poly.reserve(pts.size());
  for (const auto& p : pts) {
    const double xn = (p.day - t0) / trange;
    const double yn = (p.price - minP) / range;
    const float x = p0.x + (float)(xn * (double)size.x);
    const float y = p0.y + (float)((1.0 - yn) * (double)size.y);
    poly.push_back(ImVec2(x, y));
  }

  if (poly.size() >= 2) {
    draw->AddPolyline(poly.data(), (int)poly.size(), colLine, false, 2.0f);
  }

  if (hovered && !poly.empty()) {
    const float mx = ImGui::GetIO().MousePos.x;
    int best = 0;
    float bestDx = std::abs(poly[0].x - mx);
    for (int i = 1; i < (int)poly.size(); ++i) {
      const float dx = std::abs(poly[i].x - mx);
      if (dx < bestDx) {
        bestDx = dx;
        best = i;
      }
    }

    const auto& pt = pts[(std::size_t)best];
    const ImVec2 mp = poly[(std::size_t)best];

    draw->AddLine(ImVec2(mp.x, p0.y), ImVec2(mp.x, p1.y), colLineHover, 1.0f);
    draw->AddCircleFilled(mp, 3.0f, colLineHover);

    ImGui::BeginTooltip();
    ImGui::Text("Day: %.2f (%.2f d ago)", pt.day, std::max(0.0, nowDay - pt.day));
    ImGui::Text("Price: %.1f", pt.price);
    if (std::isfinite(referencePrice) && referencePrice > 1e-9) {
      const double pct = (pt.price - referencePrice) / referencePrice * 100.0;
      ImGui::Text("vs now: %+0.1f%%", pct);
    }
    ImGui::EndTooltip();
  }

  // Range hint below plot.
  ImGui::TextDisabled("Samples: %zu   days [%.1f .. %.1f]", pts.size(), t0, t1);
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

  const double bidAskSpread = 0.10;
  const auto q = stellar::econ::quote(econState, station->economyModel, commodity, bidAskSpread);

  ImGui::SeparatorText("Price analytics");
  ImGui::Text("%s / %s", row->systemName.c_str(), row->stationName.c_str());
  ImGui::TextDisabled("%s, fee %.1f%%", stationTypeName(row->stationType), row->feeRate * 100.0);

  // Window controls.
  {
    ImGui::TextDisabled("History window:");
    ImGui::SameLine();
    if (ImGui::SmallButton("7d")) st.historyWindowDays = 7.0;
    ImGui::SameLine();
    if (ImGui::SmallButton("30d")) st.historyWindowDays = 30.0;
    ImGui::SameLine();
    if (ImGui::SmallButton("90d")) st.historyWindowDays = 90.0;
    ImGui::SameLine();
    if (ImGui::SmallButton("All")) st.historyWindowDays = 0.0;

    double w = (st.historyWindowDays <= 0.0) ? 30.0 : st.historyWindowDays;
    const double minW = 2.0;
    const double maxW = 180.0;
    if (ImGui::SliderScalar("##hist_window", ImGuiDataType_Double, &w, &minW, &maxW, "%.0f d", ImGuiSliderFlags_Logarithmic)) {
      st.historyWindowDays = w;
    }
  }

  // Plot.
  drawPriceHistoryPlot(hist, ctx.timeDays, st.historyWindowDays, "##mid_history_plot", q.mid);

  // Trade-relevant numbers (incl. fees).
  const double askEff = q.ask * (1.0 + row->feeRate);
  const double bidEff = q.bid * (1.0 - row->feeRate);

  ImGui::Text("Mid: %.1f  |  Ask (buy): %.1f  |  Bid (sell): %.1f", q.mid, askEff, bidEff);
  ImGui::Text("Inventory: %.0f / %.0f", row->inventory, row->capacity);

  // Window stats.
  const auto stt = stellar::sim::analyzePriceHistory(hist, ctx.timeDays, st.historyWindowDays);
  if (stt.valid) {
    const double span = std::max(0.0, stt.lastDay - stt.firstDay);
    const double z = (stt.stdDevPrice > 1e-9) ? ((q.mid - stt.meanPrice) / stt.stdDevPrice) : 0.0;

    ImGui::Separator();
    ImGui::Text("Window (%.0f d, %zu samples): min %.1f  mean %.1f  max %.1f",
                span, stt.samples, stt.minPrice, stt.meanPrice, stt.maxPrice);
    ImGui::Text("Δ: %.2f%%  |  slope: %.2f cr/u/day  |  vol: %.2f%%/day  |  z: %.2f",
                stt.pctChange, stt.slopePerDay, stt.volatilityPerDay * 100.0, z);
    ImGui::TextDisabled("Fit R^2: %.2f", stt.r2);
  } else {
    ImGui::TextDisabled("Not enough history in the selected window.");
  }

  // Naive forecast (expected drift only).
  {
    const double mid1d = stellar::sim::forecastMidPrice(econState, station->economyModel, commodity, 1.0);
    const double pct = (q.mid > 1e-9) ? ((mid1d - q.mid) / q.mid * 100.0) : 0.0;

    const std::size_t i = cidx(commodity);
    const double netPerDay = station->economyModel.productionPerDay[i] - station->economyModel.consumptionPerDay[i];

    ImGui::Separator();
    ImGui::Text("Expected mid (1d, mean shock=0): %.1f (%+.2f%%)", mid1d, pct);
    ImGui::TextDisabled("Net flow: %+0.1f units/day (prod %.1f - cons %.1f)",
                        netPerDay,
                        station->economyModel.productionPerDay[i],
                        station->economyModel.consumptionPerDay[i]);
  }

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

    const double minTw = 2.0;
    const double maxTw = 30.0;
    ImGui::SliderScalar("Trend window (days)", ImGuiDataType_Double, &st.trendWindowDays, &minTw, &maxTw, "%.0f d", ImGuiSliderFlags_Logarithmic);
  }

  // Refresh cache.
  const int stamp = cacheStampFor(ctx.timeDays);
  const bool needsRebuild =
    st.cacheStamp != stamp ||
    st.cacheCommodity != st.commodityIndex ||
    st.cacheScopeNearby != st.scopeNearby ||
    st.cacheRadiusLy != st.radiusLy ||
    st.cacheMaxSystems != st.maxSystems ||
    st.cacheIncludeCurrentSystem != st.includeCurrentSystem ||
    st.cacheTrendWindowDays != st.trendWindowDays;

  if (needsRebuild) {
    st.cacheStamp = stamp;
    st.cacheCommodity = st.commodityIndex;
    st.cacheScopeNearby = st.scopeNearby;
    st.cacheRadiusLy = st.radiusLy;
    st.cacheMaxSystems = st.maxSystems;
    st.cacheIncludeCurrentSystem = st.includeCurrentSystem;
    st.cacheTrendWindowDays = st.trendWindowDays;
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
      drawTrendCell(r, st.trendWindowDays, ctx.timeDays);

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
