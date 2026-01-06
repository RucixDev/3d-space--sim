#include "TradePlannerWindow.h"

#include "stellar/econ/Commodity.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include <imgui.h>

namespace stellar::game {

namespace {

int cacheStampFor(double timeDays) {
  if (!std::isfinite(timeDays)) return -999999;
  // The economy updates at sub-day resolution. A quarter-day stamp keeps the
  // UI feeling responsive without refreshing every frame.
  return (int)std::floor(timeDays * 4.0);
}

const char* modeLabel(TradePlannerWindowState::Mode m) {
  return (m == TradePlannerWindowState::Mode::Runs) ? "Runs" : "Loops";
}

const char* scoreModeLabel(stellar::sim::TradeRunScoreMode m) {
  using M = stellar::sim::TradeRunScoreMode;
  switch (m) {
    case M::TotalProfit: return "Total profit";
    case M::ProfitPerLy: return "Profit / ly";
    case M::ProfitPerHop: return "Profit / hop";
    case M::ProfitPerCost: return "Profit / cost";
    default: return "Score";
  }
}

const char* stationTypeName(stellar::econ::StationType t) {
  using stellar::econ::StationType;
  switch (t) {
    case StationType::Outpost: return "Outpost";
    case StationType::Agricultural: return "Agricultural";
    case StationType::Mining: return "Mining";
    case StationType::Refinery: return "Refinery";
    case StationType::Industrial: return "Industrial";
    case StationType::Research: return "Research";
    case StationType::TradeHub: return "Trade Hub";
    case StationType::Shipyard: return "Shipyard";
    default: return "Station";
  }
}

std::string stopChainForRun(const stellar::sim::TradeRun& r) {
  std::ostringstream oss;
  if (r.legs.empty()) return "";

  // Start at the first leg's origin.
  oss << r.legs.front().fromSystemName << ":" << r.legs.front().fromStationName;
  for (const auto& leg : r.legs) {
    oss << " -> " << leg.toSystemName << ":" << leg.toStationName;
  }
  return oss.str();
}

std::string stopChainForLoop(const stellar::sim::TradeLoop& l) {
  std::ostringstream oss;
  if (l.legs.empty()) return "";

  oss << l.legs.front().fromSystemName << ":" << l.legs.front().fromStationName;
  for (const auto& leg : l.legs) {
    oss << " -> " << leg.toSystemName << ":" << leg.toStationName;
  }
  return oss.str();
}

bool exportRunsCsv(const char* path, const std::vector<stellar::sim::TradeRun>& runs) {
  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f) return false;

  f << "index,total_profit_cr,profit_per_ly,profit_per_hop,profit_per_cost,total_distance_ly,total_hops,legs,stop_chain\n";
  for (std::size_t i = 0; i < runs.size(); ++i) {
    const auto& r = runs[i];
    std::string chain = stopChainForRun(r);
    // Escape quotes for a single CSV string field.
    for (char& c : chain) {
      if (c == '"') c = '\'';
    }
    f << i
      << "," << r.totalProfitCr
      << "," << r.profitPerLy
      << "," << r.profitPerHop
      << "," << r.profitPerCost
      << "," << r.totalRouteDistanceLy
      << "," << r.totalHops
      << "," << r.legs.size()
      << ",\"" << chain << "\"\n";
  }
  return true;
}

bool exportLoopsCsv(const char* path, const std::vector<stellar::sim::TradeLoop>& loops) {
  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f) return false;

  f << "index,total_profit_cr,profit_per_ly,total_distance_ly,legs,stop_chain\n";
  for (std::size_t i = 0; i < loops.size(); ++i) {
    const auto& l = loops[i];
    std::string chain = stopChainForLoop(l);
    for (char& c : chain) {
      if (c == '"') c = '\'';
    }
    f << i
      << "," << l.totalProfitCr
      << "," << l.profitPerLy
      << "," << l.totalDistanceLy
      << "," << l.legs.size()
      << ",\"" << chain << "\"\n";
  }
  return true;
}

void toastMaybe(const TradePlannerContext& ctx, std::string_view msg, double ttl) {
  if (ctx.toast) ctx.toast(msg, ttl);
}

// Lazily (re)create a job system if the requested thread count changes.
void ensureJobSystem(TradePlannerWindowState& st) {
  if (!st.useParallel) {
    st.jobs.reset();
    st.jobsThreadCount = 0;
    return;
  }

  std::size_t desired = 0;
  if (st.threads > 0) desired = (std::size_t)std::max(1, st.threads);

  // If desired==0 we let JobSystem pick a default thread count.
  const bool needNew = !st.jobs || (desired != 0 && st.jobsThreadCount != desired);
  if (needNew) {
    if (desired != 0) {
      st.jobs = std::make_unique<stellar::core::JobSystem>(desired);
    } else {
      st.jobs = std::make_unique<stellar::core::JobSystem>();
    }
    st.jobsThreadCount = st.jobs ? st.jobs->threadCount() : 0;
  }
}

const stellar::sim::Station* findStationById(const stellar::sim::StarSystem& sys, stellar::sim::StationId id) {
  for (const auto& st : sys.stations) {
    if (st.id == id) return &st;
  }
  return nullptr;
}

// Draw a compact manifest summary (top N lines).
void drawManifestSummary(const stellar::econ::CargoManifestPlan& plan, int maxLines) {
  if (plan.lines.empty()) {
    ImGui::TextDisabled("No profitable cargo mix");
    return;
  }

  const int n = std::min(maxLines, (int)plan.lines.size());
  for (int i = 0; i < n; ++i) {
    const auto& ln = plan.lines[(std::size_t)i];
    const auto def = stellar::econ::commodityDef(ln.commodity);
    ImGui::BulletText("%s: %.0f u (%.0f kg)  profit %.0f cr",
                      def.name,
                      ln.units,
                      ln.massKg,
                      ln.netProfitCr);
  }
  if ((int)plan.lines.size() > n) {
    ImGui::TextDisabled("... +%d more", (int)plan.lines.size() - n);
  }
}

} // namespace

void drawTradePlannerWindow(TradePlannerWindowState& st, const TradePlannerContext& ctx) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(820, 640), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Trade Planner", &st.open)) {
    ImGui::End();
    return;
  }

  if (!ctx.currentSystem) {
    ImGui::TextDisabled("No current system.");
    ImGui::End();
    return;
  }

  const auto& sys = *ctx.currentSystem;
  const auto originStub = sys.stub;

  // Choose / validate origin station selection.
  if (st.fromStationId == 0 && ctx.preferredFromStationId != 0) {
    st.fromStationId = ctx.preferredFromStationId;
  }
  if (st.fromStationId == 0 && !sys.stations.empty()) {
    st.fromStationId = sys.stations.front().id;
  }

  const stellar::sim::Station* originStation = findStationById(sys, st.fromStationId);
  if (!originStation && !sys.stations.empty()) {
    st.fromStationId = sys.stations.front().id;
    originStation = &sys.stations.front();
  }

  // Mode selector.
  {
    int modeInt = (int)st.mode;
    const char* items[] = {"Runs (multi-leg)", "Loops (cycle)"};
    if (ImGui::Combo("Mode", &modeInt, items, IM_ARRAYSIZE(items))) {
      st.mode = (TradePlannerWindowState::Mode)modeInt;
    }
  }

  // Origin station.
  if (!originStation) {
    ImGui::TextDisabled("No stations in this system.");
    ImGui::End();
    return;
  }

  {
    std::string preview = originStation->name;
    preview += " (";
    preview += stationTypeName(originStation->type);
    preview += ")";

    if (ImGui::BeginCombo("From station", preview.c_str())) {
      for (const auto& s : sys.stations) {
        const bool selected = (s.id == st.fromStationId);
        std::string label = s.name + " (" + std::string(stationTypeName(s.type)) + ")";
        if (ImGui::Selectable(label.c_str(), selected)) {
          st.fromStationId = s.id;
          originStation = &s;
        }
        if (selected) ImGui::SetItemDefaultFocus();
      }
      ImGui::EndCombo();
    }
  }

  const double feeEff = ctx.effectiveFeeRate ? ctx.effectiveFeeRate(*originStation) : 0.0;
  ImGui::TextDisabled("Fees (rep-adjusted): %.1f%%", feeEff * 100.0);

  ImGui::Separator();

  // ---- Scope ----
  {
    float radius = (float)st.radiusLy;
    if (ImGui::SliderFloat("Search radius (ly)", &radius, 50.0f, 800.0f, "%.0f")) {
      st.radiusLy = (double)radius;
    }

    int maxSys = std::max(8, st.maxSystems);
    if (ImGui::SliderInt("Max systems", &maxSys, 8, 512)) {
      st.maxSystems = maxSys;
    }

    int maxStations = std::max(16, st.maxStations);
    if (ImGui::SliderInt("Max stations", &maxStations, 16, 512)) {
      st.maxStations = maxStations;
    }

    ImGui::Checkbox("Include same system", &st.includeSameSystem);
  }

  // ---- Cargo assumptions ----
  {
    const double freeKg = std::max(0.0, ctx.cargoCapacityKg - ctx.cargoUsedKg);
    ImGui::Text("Hold: %.0f / %.0f kg", ctx.cargoUsedKg, ctx.cargoCapacityKg);
    ImGui::SameLine();
    ImGui::TextDisabled("(free %.0f kg)", freeKg);

    ImGui::Checkbox("Use free hold", &st.useFreeHold);
    ImGui::SameLine();
    ImGui::Checkbox("Limit by credits", &st.limitByCredits);

    float spread = (float)st.bidAskSpread;
    if (ImGui::SliderFloat("Bid/ask spread", &spread, 0.0f, 0.20f, "%.2f")) {
      st.bidAskSpread = std::clamp((double)spread, 0.0, 0.5);
    }

    float stepKg = (float)st.stepKg;
    if (ImGui::SliderFloat("Manifest step (kg)", &stepKg, 0.1f, 5.0f, "%.1f")) {
      st.stepKg = std::max(0.05, (double)stepKg);
    }

    ImGui::Checkbox("Simulate price impact", &st.simulatePriceImpact);
  }

  // ---- Profit filters ----
  {
    float minLeg = (float)st.minLegProfitCr;
    if (ImGui::SliderFloat("Min leg profit (cr)", &minLeg, 0.0f, 20000.0f, "%.0f")) {
      st.minLegProfitCr = std::max(0.0, (double)minLeg);
    }

    float minTotal = (float)st.minTotalProfitCr;
    if (ImGui::SliderFloat("Min total profit (cr)", &minTotal, 0.0f, 80000.0f, "%.0f")) {
      st.minTotalProfitCr = std::max(0.0, (double)minTotal);
    }

    int maxRes = std::max(1, st.maxResults);
    if (ImGui::SliderInt("Max results", &maxRes, 1, 40)) {
      st.maxResults = maxRes;
    }
  }

  ImGui::Separator();

  // ---- Mode-specific tuning ----
  if (st.mode == TradePlannerWindowState::Mode::Runs) {
    int legs = std::clamp(st.runLegs, 1, 4);
    if (ImGui::SliderInt("Legs", &legs, 1, 4)) st.runLegs = legs;

    int beam = std::clamp(st.beamWidth, 4, 128);
    if (ImGui::SliderInt("Beam width", &beam, 4, 128)) st.beamWidth = beam;

    int legCand = std::clamp(st.maxLegCandidates, 4, 48);
    if (ImGui::SliderInt("Candidates / leg", &legCand, 4, 48)) st.maxLegCandidates = legCand;

    ImGui::Checkbox("Loopless", &st.loopless);

    ImGui::Separator();
    ImGui::Checkbox("Enforce jump range", &st.enforceJumpRange);
    if (st.enforceJumpRange) {
      float jr = (float)st.jumpRangeLy;
      if (ImGui::SliderFloat("Jump range (ly)", &jr, 2.0f, 40.0f, "%.1f")) {
        st.jumpRangeLy = std::max(0.1, (double)jr);
      }
      float costHop = (float)st.routeCostPerJump;
      if (ImGui::SliderFloat("Route cost / hop", &costHop, 0.0f, 10.0f, "%.2f")) {
        st.routeCostPerJump = std::max(0.0, (double)costHop);
      }
      float costLy = (float)st.routeCostPerLy;
      if (ImGui::SliderFloat("Route cost / ly", &costLy, 0.0f, 2.0f, "%.2f")) {
        st.routeCostPerLy = std::max(0.0, (double)costLy);
      }
    }

    // Score mode.
    {
      using M = stellar::sim::TradeRunScoreMode;
      const M modes[] = {M::TotalProfit, M::ProfitPerLy, M::ProfitPerHop, M::ProfitPerCost};
      int idx = std::clamp(st.scoreMode, 0, (int)(IM_ARRAYSIZE(modes) - 1));
      const char* labels[] = {"Total profit", "Profit / ly", "Profit / hop", "Profit / cost"};
      if (ImGui::Combo("Score", &idx, labels, IM_ARRAYSIZE(labels))) {
        st.scoreMode = idx;
      }
      ImGui::SameLine();
      ImGui::TextDisabled("(%s)", scoreModeLabel(modes[idx]));
    }
  } else {
    int legs = std::clamp(st.loopLegs, 2, 3);
    if (ImGui::SliderInt("Loop legs", &legs, 2, 3)) st.loopLegs = legs;
  }

  ImGui::Separator();

  // ---- Performance ----
  {
    ImGui::Checkbox("Use parallel scan", &st.useParallel);
    ImGui::SameLine();
    ImGui::Checkbox("Auto refresh", &st.autoRefresh);

    if (st.useParallel) {
      int th = st.threads;
      if (ImGui::SliderInt("Threads (0=auto)", &th, 0, 16)) {
        st.threads = std::max(0, th);
      }
      if (st.jobs) {
        ImGui::SameLine();
        ImGui::TextDisabled("active: %zu", st.jobsThreadCount);
      }
    }
  }

  // ---- Actions ----
  const int stampNow = cacheStampFor(ctx.timeDays);
  const double capNow = ctx.cargoCapacityKg;
  const double usedNow = ctx.cargoUsedKg;
  const double creditsNow = ctx.playerCreditsCr;
  const bool systemChanged = (st.cacheSystemId != originStub.id);

  const bool cacheStale = systemChanged || st.cacheMode != st.mode || st.cacheStamp != stampNow ||
                          st.cacheCargoCapKg != capNow || st.cacheCargoUsedKg != usedNow ||
                          st.cacheCreditsCr != creditsNow;

  if (cacheStale) {
    ImGui::TextDisabled("Results are stale (%s%s).",
                        systemChanged ? "system changed" : "inputs changed",
                        (st.autoRefresh ? ", auto" : ""));
  }

  const bool canScan = originStation && (ctx.cargoCapacityKg > 1e-6);

  if (!canScan) {
    ImGui::TextDisabled("Cannot scan: no origin station or cargo capacity.");
  }

  bool doScan = false;
  if (ImGui::Button("Scan")) doScan = true;
  ImGui::SameLine();
  if (ImGui::Button("Clear")) {
    st.runs.clear();
    st.loops.clear();
    st.cacheStamp = -999999;
    st.cacheSystemId = 0;
  }

  if (st.autoRefresh && cacheStale && canScan) {
    doScan = true;
  }

  ImGui::SameLine();
  ImGui::TextDisabled("Last scan: %.1f ms", st.lastComputeMs);

  // Export.
  {
    ImGui::InputText("Export path", st.exportPath, IM_ARRAYSIZE(st.exportPath));
    ImGui::SameLine();
    if (ImGui::Button("Export CSV")) {
      bool ok = false;
      if (st.mode == TradePlannerWindowState::Mode::Runs) {
        ok = exportRunsCsv(st.exportPath, st.runs);
      } else {
        ok = exportLoopsCsv(st.exportPath, st.loops);
      }
      if (ok) {
        toastMaybe(ctx, std::string("Exported ") + modeLabel(st.mode) + " to '" + st.exportPath + "'.", 3.0);
      } else {
        toastMaybe(ctx, std::string("Failed to export to '") + st.exportPath + "'.", 4.0);
      }
    }
  }

  // ---- Execute scan (synchronously; uses internal parallelism if enabled) ----
  if (doScan && canScan) {
    // Build manifest params.
    stellar::econ::CargoManifestParams mp;
    const double freeKg = std::max(0.0, ctx.cargoCapacityKg - ctx.cargoUsedKg);
    mp.cargoCapacityKg = st.useFreeHold ? freeKg : ctx.cargoCapacityKg;
    mp.bidAskSpread = st.bidAskSpread;
    mp.stepKg = st.stepKg;
    mp.simulatePriceImpact = st.simulatePriceImpact;
    mp.maxBuyCreditsCr = (st.limitByCredits ? ctx.playerCreditsCr : 0.0);

    // If there is effectively no free capacity, bail quickly.
    if (mp.cargoCapacityKg <= 1e-6) {
      toastMaybe(ctx, "No free cargo capacity for trade planning.", 3.0);
    } else {
      auto start = std::chrono::steady_clock::now();

      auto feeFn = [&](const stellar::sim::Station& s) {
        return ctx.effectiveFeeRate ? ctx.effectiveFeeRate(s) : 0.0;
      };

      if (st.useParallel) ensureJobSystem(st);

      if (st.mode == TradePlannerWindowState::Mode::Runs) {
        stellar::sim::TradeRunScanParams p;
        p.legs = (std::size_t)std::clamp(st.runLegs, 1, 6);
        p.beamWidth = (std::size_t)std::clamp(st.beamWidth, 1, 256);
        p.maxLegCandidates = (std::size_t)std::clamp(st.maxLegCandidates, 1, 256);
        p.loopless = st.loopless;
        p.includeSameSystem = st.includeSameSystem;
        p.maxResults = (std::size_t)std::clamp(st.maxResults, 1, 128);
        p.maxStations = (std::size_t)std::clamp(st.maxStations, 2, 4096);

        p.minLegProfitCr = st.minLegProfitCr;
        p.minRunProfitCr = st.minTotalProfitCr;

        p.jumpRangeLy = st.enforceJumpRange ? std::max(0.1, st.jumpRangeLy) : 0.0;
        p.routeCostPerJump = std::max(0.0, st.routeCostPerJump);
        p.routeCostPerLy = std::max(0.0, st.routeCostPerLy);

        using SM = stellar::sim::TradeRunScoreMode;
        const SM modes[] = {SM::TotalProfit, SM::ProfitPerLy, SM::ProfitPerHop, SM::ProfitPerCost};
        const int idx = std::clamp(st.scoreMode, 0, (int)(IM_ARRAYSIZE(modes) - 1));
        p.scoreMode = modes[idx];

        p.manifest = mp;

        if (st.useParallel && st.jobs) {
          st.runs = stellar::sim::planTradeRunsParallel(*st.jobs,
                                                       ctx.universe,
                                                       originStub,
                                                       *originStation,
                                                       ctx.timeDays,
                                                       st.radiusLy,
                                                       (std::size_t)std::max(8, st.maxSystems),
                                                       p,
                                                       feeFn);
        } else {
          st.runs = stellar::sim::planTradeRuns(ctx.universe,
                                               originStub,
                                               *originStation,
                                               ctx.timeDays,
                                               st.radiusLy,
                                               (std::size_t)std::max(8, st.maxSystems),
                                               p,
                                               feeFn);
        }
        st.loops.clear();
      } else {
        stellar::sim::TradeLoopScanParams p;
        p.legs = (std::size_t)std::clamp(st.loopLegs, 2, 3);
        p.includeSameSystem = st.includeSameSystem;
        p.maxResults = (std::size_t)std::clamp(st.maxResults, 1, 128);
        p.maxStations = (std::size_t)std::clamp(st.maxStations, 2, 4096);
        p.minLegProfitCr = st.minLegProfitCr;
        p.minLoopProfitCr = st.minTotalProfitCr;
        p.manifest = mp;

        if (st.useParallel && st.jobs) {
          st.loops = stellar::sim::scanTradeLoopsParallel(*st.jobs,
                                                         ctx.universe,
                                                         originStub,
                                                         *originStation,
                                                         ctx.timeDays,
                                                         st.radiusLy,
                                                         (std::size_t)std::max(8, st.maxSystems),
                                                         p,
                                                         feeFn);
        } else {
          st.loops = stellar::sim::scanTradeLoops(ctx.universe,
                                                 originStub,
                                                 *originStation,
                                                 ctx.timeDays,
                                                 st.radiusLy,
                                                 (std::size_t)std::max(8, st.maxSystems),
                                                 p,
                                                 feeFn);
        }
        st.runs.clear();
      }

      const auto end = std::chrono::steady_clock::now();
      st.lastComputeMs = std::chrono::duration<double, std::milli>(end - start).count();

      st.cacheSystemId = originStub.id;
      st.cacheStamp = stampNow;
      st.cacheCargoCapKg = capNow;
      st.cacheCargoUsedKg = usedNow;
      st.cacheCreditsCr = creditsNow;
      st.cacheMode = st.mode;

      toastMaybe(ctx, std::string("Trade ") + modeLabel(st.mode) + " scan complete.", 2.2);
    }
  }

  ImGui::Separator();

  // ---- Results ----
  if (st.mode == TradePlannerWindowState::Mode::Runs) {
    ImGui::Text("Runs: %d", (int)st.runs.size());
    if (st.runs.empty()) {
      ImGui::TextDisabled("No runs found.");
    } else {
      if (ImGui::BeginTable("trade_runs", 7, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, ImVec2(0, 260))) {
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 26.0f);
        ImGui::TableSetupColumn("Profit", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("/ly", ImGuiTableColumnFlags_WidthFixed, 60.0f);
        ImGui::TableSetupColumn("/hop", ImGuiTableColumnFlags_WidthFixed, 60.0f);
        ImGui::TableSetupColumn("Hops", ImGuiTableColumnFlags_WidthFixed, 40.0f);
        ImGui::TableSetupColumn("Dist (ly)", ImGuiTableColumnFlags_WidthFixed, 70.0f);
        ImGui::TableSetupColumn("Route", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableHeadersRow();

        for (std::size_t i = 0; i < st.runs.size(); ++i) {
          const auto& r = st.runs[i];
          ImGui::TableNextRow();

          ImGui::TableNextColumn();
          const std::string nodeId = "##run_" + std::to_string(i);
          const bool open = ImGui::TreeNodeEx(nodeId.c_str(), ImGuiTreeNodeFlags_SpanFullWidth | ImGuiTreeNodeFlags_NoTreePushOnOpen, "%zu", i + 1);

          ImGui::TableNextColumn();
          ImGui::Text("%.0f", r.totalProfitCr);
          ImGui::TableNextColumn();
          ImGui::Text("%.0f", r.profitPerLy);
          ImGui::TableNextColumn();
          ImGui::Text("%.0f", r.profitPerHop);
          ImGui::TableNextColumn();
          ImGui::Text("%d", (int)r.totalHops);
          ImGui::TableNextColumn();
        ImGui::Text("%.1f", r.totalRouteDistanceLy);
          ImGui::TableNextColumn();
          std::string chain = stopChainForRun(r);
          ImGui::TextUnformatted(chain.c_str());

          if (open) {
            ImGui::TreePush(nodeId.c_str());
            ImGui::Separator();
            if (ImGui::SmallButton(("Copy##run" + std::to_string(i)).c_str())) {
              ImGui::SetClipboardText(chain.c_str());
              toastMaybe(ctx, "Run copied to clipboard.", 1.8);
            }
            ImGui::SameLine();
            if (ImGui::SmallButton(("Plot final##run" + std::to_string(i)).c_str())) {
              if (ctx.routeToStation && !r.legs.empty()) {
                const auto& last = r.legs.back();
                ctx.routeToStation(last.toSystem, last.toStation);
              }
            }

            for (std::size_t li = 0; li < r.legs.size(); ++li) {
              const auto& leg = r.legs[li];
              ImGui::PushID((int)li);
              ImGui::Text("Leg %zu: %s:%s  ->  %s:%s", li + 1,
                          leg.fromSystemName.c_str(), leg.fromStationName.c_str(),
                          leg.toSystemName.c_str(), leg.toStationName.c_str());
              ImGui::TextDisabled("profit %.0f cr | dist %.2f ly | hops %d",
                                  leg.manifest.netProfitCr, leg.routeDistanceLy, (int)leg.routeHops);

              if (ctx.routeToStation && ImGui::SmallButton("Plot leg")) {
                ctx.routeToStation(leg.toSystem, leg.toStation);
              }
              ImGui::SameLine();
              if (ImGui::SmallButton("Copy manifest")) {
                std::ostringstream oss;
                oss << "Leg " << (li + 1) << ": "
                    << leg.fromSystemName << ":" << leg.fromStationName
                    << " -> " << leg.toSystemName << ":" << leg.toStationName
                    << "\n";
                oss << "Profit: " << leg.manifest.netProfitCr << " cr\n";
                for (const auto& ln : leg.manifest.lines) {
                  const auto def = stellar::econ::commodityDef(ln.commodity);
                  oss << "- " << def.name << ": " << ln.units << " u (" << ln.massKg << " kg) profit " << ln.netProfitCr << " cr\n";
                }
                ImGui::SetClipboardText(oss.str().c_str());
                toastMaybe(ctx, "Manifest copied to clipboard.", 1.8);
              }

              drawManifestSummary(leg.manifest, 3);
              ImGui::Separator();
              ImGui::PopID();
            }

            ImGui::TreePop();
          }
        }

        ImGui::EndTable();
      }
    }
  } else {
    ImGui::Text("Loops: %d", (int)st.loops.size());
    if (st.loops.empty()) {
      ImGui::TextDisabled("No loops found.");
    } else {
      if (ImGui::BeginTable("trade_loops", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, ImVec2(0, 260))) {
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 26.0f);
        ImGui::TableSetupColumn("Profit", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("/ly", ImGuiTableColumnFlags_WidthFixed, 60.0f);
        ImGui::TableSetupColumn("Dist (ly)", ImGuiTableColumnFlags_WidthFixed, 70.0f);
        ImGui::TableSetupColumn("Route", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableHeadersRow();

        for (std::size_t i = 0; i < st.loops.size(); ++i) {
          const auto& l = st.loops[i];
          ImGui::TableNextRow();

          ImGui::TableNextColumn();
          const std::string nodeId = "##loop_" + std::to_string(i);
          const bool open = ImGui::TreeNodeEx(nodeId.c_str(), ImGuiTreeNodeFlags_SpanFullWidth | ImGuiTreeNodeFlags_NoTreePushOnOpen, "%zu", i + 1);

          ImGui::TableNextColumn();
          ImGui::Text("%.0f", l.totalProfitCr);
          ImGui::TableNextColumn();
          ImGui::Text("%.0f", l.profitPerLy);
          ImGui::TableNextColumn();
          ImGui::Text("%.1f", l.totalDistanceLy);
          ImGui::TableNextColumn();
          std::string chain = stopChainForLoop(l);
          ImGui::TextUnformatted(chain.c_str());

          if (open) {
            ImGui::TreePush(nodeId.c_str());
            ImGui::Separator();
            if (ImGui::SmallButton(("Copy##loop" + std::to_string(i)).c_str())) {
              ImGui::SetClipboardText(chain.c_str());
              toastMaybe(ctx, "Loop copied to clipboard.", 1.8);
            }
            ImGui::SameLine();
            if (ImGui::SmallButton(("Plot first##loop" + std::to_string(i)).c_str())) {
              if (ctx.routeToStation && !l.legs.empty()) {
                ctx.routeToStation(l.legs.front().toSystem, l.legs.front().toStation);
              }
            }

            for (std::size_t li = 0; li < l.legs.size(); ++li) {
              const auto& leg = l.legs[li];
              ImGui::PushID((int)li);
              ImGui::Text("Leg %zu: %s:%s  ->  %s:%s", li + 1,
                          leg.fromSystemName.c_str(), leg.fromStationName.c_str(),
                          leg.toSystemName.c_str(), leg.toStationName.c_str());
              ImGui::TextDisabled("profit %.0f cr | dist %.2f ly",
                                  leg.manifest.netProfitCr, leg.distanceLy);

              if (ctx.routeToStation && ImGui::SmallButton("Plot leg")) {
                ctx.routeToStation(leg.toSystem, leg.toStation);
              }
              drawManifestSummary(leg.manifest, 3);
              ImGui::Separator();
              ImGui::PopID();
            }

            ImGui::TreePop();
          }
        }

        ImGui::EndTable();
      }
    }
  }

  ImGui::End();
}

} // namespace stellar::game
