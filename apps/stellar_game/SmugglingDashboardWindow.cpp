#include "SmugglingDashboardWindow.h"

#include "stellar/econ/Commodity.h"
#include "stellar/sim/SystemEvents.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstring>
#include <string>

#include <imgui.h>

namespace stellar::game {

static bool icontains(std::string_view hay, std::string_view needle) {
  if (needle.empty()) return true;
  if (hay.empty()) return false;

  auto lower = [](unsigned char c) -> unsigned char { return (unsigned char)std::tolower(c); };

  const std::size_t n = needle.size();
  for (std::size_t i = 0; i + n <= hay.size(); ++i) {
    bool ok = true;
    for (std::size_t j = 0; j < n; ++j) {
      if (lower((unsigned char)hay[i + j]) != lower((unsigned char)needle[j])) {
        ok = false;
        break;
      }
    }
    if (ok) return true;
  }
  return false;
}

static const char* availabilityName(stellar::sim::SmugglingAvailabilityMode m) {
  using stellar::sim::SmugglingAvailabilityMode;
  switch (m) {
    case SmugglingAvailabilityMode::TodayOnly: return "Today only";
    case SmugglingAvailabilityMode::Expected:  return "Expected access";
    case SmugglingAvailabilityMode::Ignore:    return "Ignore availability";
    default:                                   return "Today only";
  }
}

static const char* scoreModeName(stellar::sim::SmugglingScoreMode m) {
  using stellar::sim::SmugglingScoreMode;
  switch (m) {
    case SmugglingScoreMode::ExpectedProfit:      return "Expected profit";
    case SmugglingScoreMode::RiskAdjusted:        return "Risk-adjusted";
    case SmugglingScoreMode::CleanProfit:         return "Clean profit";
    case SmugglingScoreMode::ExpectedProfitPerLy: return "Expected / ly";
    default:                                      return "Expected profit";
  }
}

static bool beginComboEnum(const char* label,
                           stellar::sim::SmugglingAvailabilityMode& v) {
  const char* cur = availabilityName(v);
  bool changed = false;
  if (ImGui::BeginCombo(label, cur)) {
    using stellar::sim::SmugglingAvailabilityMode;
    for (const auto m : {SmugglingAvailabilityMode::TodayOnly,
                         SmugglingAvailabilityMode::Expected,
                         SmugglingAvailabilityMode::Ignore}) {
      const bool isSel = (m == v);
      if (ImGui::Selectable(availabilityName(m), isSel)) {
        v = m;
        changed = true;
      }
      if (isSel) ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }
  return changed;
}

static bool beginComboEnum(const char* label,
                           stellar::sim::SmugglingScoreMode& v) {
  const char* cur = scoreModeName(v);
  bool changed = false;
  if (ImGui::BeginCombo(label, cur)) {
    using stellar::sim::SmugglingScoreMode;
    for (const auto m : {SmugglingScoreMode::ExpectedProfit,
                         SmugglingScoreMode::RiskAdjusted,
                         SmugglingScoreMode::CleanProfit,
                         SmugglingScoreMode::ExpectedProfitPerLy}) {
      const bool isSel = (m == v);
      if (ImGui::Selectable(scoreModeName(m), isSel)) {
        v = m;
        changed = true;
      }
      if (isSel) ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }
  return changed;
}

static void maybeToast(const SmugglingDashboardContext& ctx, std::string_view msg, double ttlSec) {
  if (ctx.toast) ctx.toast(msg, ttlSec);
}

static stellar::sim::StationId pickDefaultOrigin(const SmugglingDashboardWindowState& st,
                                                 const SmugglingDashboardContext& ctx) {
  if (st.originStationId != 0) return st.originStationId;
  if (ctx.preferredOriginStationId != 0) return ctx.preferredOriginStationId;
  if (ctx.currentSystem && !ctx.currentSystem->stations.empty()) return ctx.currentSystem->stations.front().id;
  return 0;
}

static bool cacheMatches(const SmugglingDashboardWindowState& st,
                         stellar::sim::StationId originStationId,
                         int dayStamp) {
  if (st.cacheDay != dayStamp) return false;
  if (st.cacheOriginStationId != originStationId) return false;
  if (std::fabs(st.cacheRadiusLy - st.radiusLy) > 1e-6) return false;
  if (st.cacheMaxSystems != st.maxSystems) return false;
  if (st.cacheIncludeSameSystem != st.includeSameSystem) return false;
  if (st.cacheUseFreeHold != st.useFreeHold) return false;
  if (std::fabs(st.cacheBidAskSpread - st.bidAskSpread) > 1e-6) return false;
  if (st.cacheRequireOriginLegal != st.requireOriginLegal) return false;
  if (st.cacheUseLiveSystemConditions != st.useLiveSystemConditions) return false;
  if (st.cacheAvailability != st.availability) return false;
  if (st.cacheScoreMode != st.scoreMode) return false;
  if (std::fabs(st.cacheRiskLambda - st.riskLambda) > 1e-6) return false;
  if (std::fabs(st.cacheMinScoreCr - st.minScoreCr) > 1e-6) return false;
  if (st.cacheMaxResults != st.maxResults) return false;
  if (st.cachePerStationLimit != st.perStationLimit) return false;
  return true;
}

static void updateCache(SmugglingDashboardWindowState& st,
                        stellar::sim::StationId originStationId,
                        int dayStamp) {
  st.cacheDay = dayStamp;
  st.cacheOriginStationId = originStationId;
  st.cacheRadiusLy = st.radiusLy;
  st.cacheMaxSystems = st.maxSystems;
  st.cacheIncludeSameSystem = st.includeSameSystem;
  st.cacheUseFreeHold = st.useFreeHold;
  st.cacheBidAskSpread = st.bidAskSpread;
  st.cacheRequireOriginLegal = st.requireOriginLegal;
  st.cacheUseLiveSystemConditions = st.useLiveSystemConditions;
  st.cacheAvailability = st.availability;
  st.cacheScoreMode = st.scoreMode;
  st.cacheRiskLambda = st.riskLambda;
  st.cacheMinScoreCr = st.minScoreCr;
  st.cacheMaxResults = st.maxResults;
  st.cachePerStationLimit = st.perStationLimit;
}

static void doScan(SmugglingDashboardWindowState& st,
                   const SmugglingDashboardContext& ctx,
                   const stellar::sim::StarSystem& originSys,
                   const stellar::sim::Station& originStation) {
  using Clock = std::chrono::high_resolution_clock;
  const auto t0 = Clock::now();

  stellar::sim::SmugglingScanParams p{};
  p.radiusLy = st.radiusLy;
  p.maxSystems = (std::size_t)std::max(1, st.maxSystems);
  p.includeSameSystem = st.includeSameSystem;
  p.useFreeHold = st.useFreeHold;
  p.cargoCapacityKg = ctx.cargoCapacityKg;
  p.cargoUsedKg = ctx.cargoUsedKg;

  p.playerHeat = ctx.playerHeat;
  p.smuggleHoldMk = ctx.smuggleHoldMk;

  // Use origin faction rep as a sane fallback, but prefer per-faction rep when available.
  if (ctx.reputationForFaction) {
    p.playerRep = ctx.reputationForFaction(originStation.factionId);
    p.repForFaction = ctx.reputationForFaction;
  }

  p.bidAskSpread = st.bidAskSpread;
  p.requireOriginLegal = st.requireOriginLegal;

  p.availability = st.availability;
  p.scoreMode = st.scoreMode;
  p.riskLambda = st.riskLambda;
  p.minScoreCr = st.minScoreCr;
  p.maxResults = std::max(1, st.maxResults);
  p.perStationLimit = std::max(1, st.perStationLimit);

  p.useLiveSystemConditions = st.useLiveSystemConditions;

  // Candidates: query once so we can show scale & avoid repeated allocations.
  const auto candidates = ctx.universe.queryNearby(originSys.stub.posLy,
                                                   p.radiusLy,
                                                   std::max<std::size_t>(1, p.maxSystems));

  // Fee model (origin buy only).
  stellar::sim::SmugglingFeeRateFn feeFn = ctx.effectiveFeeRate;
  st.results = stellar::sim::scanSmugglingOpportunities(ctx.universe,
                                                        originSys.stub,
                                                        originStation,
                                                        ctx.timeDays,
                                                        candidates,
                                                        p,
                                                        std::move(feeFn));

  const auto t1 = Clock::now();
  st.lastComputeMs = std::chrono::duration<double, std::milli>(t1 - t0).count();

  maybeToast(ctx,
             "Smuggling scan: " + std::to_string(st.results.size()) + " opportunities",
             1.6);
}

static void drawResultTable(SmugglingDashboardWindowState& st,
                            const SmugglingDashboardContext& ctx) {
  const std::string_view filter = std::string_view(st.filter, std::strlen(st.filter));

  const ImGuiTableFlags flags =
      ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersV | ImGuiTableFlags_RowBg |
      ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable |
      ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingFixedFit;

  if (!ImGui::BeginTable("##SmugglingOpps", 11, flags, ImVec2(0.0f, 360.0f))) return;

  ImGui::TableSetupScrollFreeze(0, 1);
  ImGui::TableSetupColumn("Score");
  ImGui::TableSetupColumn("Commodity");
  ImGui::TableSetupColumn("Units");
  ImGui::TableSetupColumn("Buy net");
  ImGui::TableSetupColumn("Sell BM");
  ImGui::TableSetupColumn("Expected");
  ImGui::TableSetupColumn("Clean");
  ImGui::TableSetupColumn("Sting");
  ImGui::TableSetupColumn("Dist");
  ImGui::TableSetupColumn("To");
  ImGui::TableSetupColumn("Action");
  ImGui::TableHeadersRow();

  for (std::size_t i = 0; i < st.results.size(); ++i) {
    const auto& r = st.results[i];

    const std::string_view sysName = r.toSystemName;
    const std::string_view stName = r.toStationName;
    const std::string_view comName = stellar::econ::commodityDef(r.commodity).name;

    if (!filter.empty()) {
      if (!icontains(sysName, filter) && !icontains(stName, filter) && !icontains(comName, filter)) {
        continue;
      }
    }

    ImGui::PushID((int)i);
    ImGui::TableNextRow();

    ImGui::TableSetColumnIndex(0);
    ImGui::Text("%.0f", r.scoreCr);

    ImGui::TableSetColumnIndex(1);
    ImGui::TextUnformatted(comName.data());

    ImGui::TableSetColumnIndex(2);
    ImGui::Text("%.0f", r.unitsPossible);

    ImGui::TableSetColumnIndex(3);
    ImGui::Text("%.2f", r.buyAskNetCr);

    ImGui::TableSetColumnIndex(4);
    ImGui::Text("%.2f", r.sellBidCr);

    ImGui::TableSetColumnIndex(5);
    ImGui::Text("%.0f", r.expectedProfitCr);

    ImGui::TableSetColumnIndex(6);
    ImGui::Text("%.0f", r.cleanProfitCr);

    ImGui::TableSetColumnIndex(7);
    ImGui::Text("%.1f%%", r.stingChance * 100.0);

    ImGui::TableSetColumnIndex(8);
    ImGui::Text("%.1f", r.distanceLy);

    ImGui::TableSetColumnIndex(9);
    ImGui::Text("%s:%s", r.toSystemName.c_str(), r.toStationName.c_str());

    ImGui::TableSetColumnIndex(10);
    const bool hasGo = (bool)ctx.goToStation;
    const bool hasPlot = (bool)ctx.routeToStation || hasGo;

    if (hasPlot) {
      if (ImGui::SmallButton("Plot")) {
        if (ctx.goToStation) ctx.goToStation(r.toSystem, r.toStation, /*armAutoRun=*/false);
        else if (ctx.routeToStation) ctx.routeToStation(r.toSystem, r.toStation);
      }
      ImGui::SameLine();
      ImGui::BeginDisabled(!hasGo);
      if (ImGui::SmallButton("Go")) {
        if (ctx.goToStation) ctx.goToStation(r.toSystem, r.toStation, /*armAutoRun=*/true);
      }
      ImGui::EndDisabled();
    } else {
      ImGui::TextDisabled("-");
    }

    // Tooltip with deeper context (BM access, conditions, per-unit edges).
    if (ImGui::IsItemHovered() || ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
      ImGui::BeginTooltip();
      ImGui::TextUnformatted("Smuggling opportunity");
      ImGui::Separator();
      ImGui::Text("To: %s:%s", r.toSystemName.c_str(), r.toStationName.c_str());
      ImGui::Text("Commodity: %s", comName.data());
      ImGui::Text("Units: %.0f  (mass %.1f kg/u)", r.unitsPossible, r.unitMassKg);
      ImGui::Text("Buy ask: %.2f  (net %.2f)", r.buyAskCr, r.buyAskNetCr);
      ImGui::Text("Sell bid: %.2f  (official %.2f)", r.sellBidCr, r.officialSellBidCr);
      ImGui::Text("Expected profit: %.0f  | Clean: %.0f", r.expectedProfitCr, r.cleanProfitCr);
      ImGui::Text("Sting: %.1f%%  Fine: %.0f", r.stingChance * 100.0, r.fineCr);
      ImGui::Text("BM access: %.0f%%  (%s)", r.blackMarketAccess01 * 100.0, r.blackMarketAvailable ? "available" : "unavailable");

      if (r.systemEventKind != stellar::sim::SystemEventKind::None && r.systemEventSeverity01 > 1e-6) {
        ImGui::Text("Event: %s (%.0f%%)", stellar::sim::systemEventKindName(r.systemEventKind), r.systemEventSeverity01 * 100.0);
      }
      ImGui::Text("Security: %.0f%%  Piracy: %.0f%%  Traffic: %.0f%%",
                  r.systemSecurity01 * 100.0,
                  r.systemPiracy01 * 100.0,
                  r.systemTraffic01 * 100.0);
      ImGui::EndTooltip();
    }

    ImGui::PopID();
  }

  ImGui::EndTable();
}

void drawSmugglingDashboardWindow(SmugglingDashboardWindowState& st, const SmugglingDashboardContext& ctx) {
  if (!st.open) return;

  bool open = st.open;
  if (!ImGui::Begin("Smuggling Dashboard", &open)) {
    ImGui::End();
    st.open = open;
    return;
  }
  st.open = open;

  if (!ctx.currentSystem) {
    ImGui::TextDisabled("No current system.");
    ImGui::End();
    return;
  }

  const auto& sys = *ctx.currentSystem;

  // Origin station picker.
  stellar::sim::StationId originId = pickDefaultOrigin(st, ctx);
  if (originId == 0) {
    ImGui::TextDisabled("No stations in the current system.");
    ImGui::End();
    return;
  }

  int originIdx = -1;
  for (std::size_t i = 0; i < sys.stations.size(); ++i) {
    if (sys.stations[i].id == originId) {
      originIdx = (int)i;
      break;
    }
  }
  if (originIdx < 0) originIdx = 0;
  originId = sys.stations[(std::size_t)originIdx].id;

  // Keep state in sync.
  st.originStationId = originId;

  const auto& originStation = sys.stations[(std::size_t)originIdx];

  if (ImGui::BeginCombo("Origin station", originStation.name.c_str())) {
    for (std::size_t i = 0; i < sys.stations.size(); ++i) {
      const bool isSel = ((int)i == originIdx);
      if (ImGui::Selectable(sys.stations[i].name.c_str(), isSel)) {
        st.originStationId = sys.stations[i].id;
        originIdx = (int)i;
      }
      if (isSel) ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }

  // Scan controls.
  int dayStamp = (int)std::floor(ctx.timeDays);

  bool changed = false;
  changed |= ImGui::Checkbox("Auto-refresh daily", &st.autoRefreshDaily);
  ImGui::SameLine();
  changed |= ImGui::Checkbox("Advanced", &st.showAdvanced);

  changed |= ImGui::InputDouble("Radius (ly)", &st.radiusLy, 25.0, 100.0, "%.0f");
  changed |= ImGui::InputInt("Max systems", &st.maxSystems);
  changed |= ImGui::InputInt("Max results", &st.maxResults);
  changed |= ImGui::InputInt("Per-station cap", &st.perStationLimit);

  changed |= ImGui::Checkbox("Include same system", &st.includeSameSystem);
  changed |= ImGui::Checkbox("Require origin legal", &st.requireOriginLegal);

  changed |= ImGui::Checkbox("Use free hold", &st.useFreeHold);
  ImGui::SameLine();
  changed |= ImGui::Checkbox("Live conditions", &st.useLiveSystemConditions);

  changed |= beginComboEnum("Availability", st.availability);
  changed |= beginComboEnum("Score mode", st.scoreMode);

  if (st.showAdvanced) {
    changed |= ImGui::InputDouble("Bid/ask spread", &st.bidAskSpread, 0.01, 0.05, "%.2f");
    changed |= ImGui::InputDouble("Risk lambda", &st.riskLambda, 0.05, 0.25, "%.2f");
    changed |= ImGui::InputDouble("Min score (cr)", &st.minScoreCr, 50.0, 250.0, "%.0f");
  }

  // Clamp some inputs gently.
  st.radiusLy = std::max(1.0, st.radiusLy);
  st.maxSystems = std::clamp(st.maxSystems, 1, 4096);
  st.maxResults = std::clamp(st.maxResults, 1, 256);
  st.perStationLimit = std::clamp(st.perStationLimit, 1, 16);
  st.bidAskSpread = std::clamp(st.bidAskSpread, 0.0, 0.50);
  st.riskLambda = std::max(0.0, st.riskLambda);
  st.minScoreCr = std::max(0.0, st.minScoreCr);

  const bool cachedOk = cacheMatches(st, st.originStationId, dayStamp);

  bool wantsScan = false;
  if (ImGui::Button("Scan now") || (st.autoRefreshDaily && !cachedOk) || (changed && st.autoRefreshDaily)) {
    wantsScan = true;
  }

  ImGui::SameLine();
  if (ImGui::Button("Clear results")) {
    st.results.clear();
  }

  ImGui::SameLine();
  ImGui::TextDisabled("Last scan: %.1f ms | %zu opps", st.lastComputeMs, st.results.size());

  if (wantsScan) {
    doScan(st, ctx, sys, originStation);
    updateCache(st, st.originStationId, dayStamp);
  }

  ImGui::Separator();

  // Filter box.
  ImGui::InputTextWithHint("Filter", "system / station / commodity", st.filter, sizeof(st.filter));

  drawResultTable(st, ctx);

  ImGui::End();
}

} // namespace stellar::game
