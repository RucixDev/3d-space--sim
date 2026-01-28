#include "MissionControlWindow.h"

#include "stellar/core/Hash.h"
#include "stellar/sim/MissionBriefing.h"
#include "stellar/ui/TextFx.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

#include <imgui.h>

namespace stellar::game {

namespace {

static ui::textfx::ProgramCache g_textFxCache(1024);

int cacheStampFor(double timeDays) {
  if (!std::isfinite(timeDays)) return -999999;
  return (int)std::floor(timeDays * 4.0);
}

const char* missionTypeLabel(stellar::sim::MissionType t) {
  using stellar::sim::MissionType;
  switch (t) {
    case MissionType::Delivery: return "Delivery";
    case MissionType::MultiDelivery: return "Multi Delivery";
    case MissionType::Courier: return "Courier";
    case MissionType::Smuggle: return "Smuggle";
    case MissionType::Passenger: return "Passenger";
    case MissionType::BountyKill: return "Bounty (Kill)";
    case MissionType::BountyScan: return "Bounty (Scan)";
    case MissionType::Salvage: return "Salvage";
    case MissionType::Escort: return "Escort";
    default: return "Mission";
  }
}

const stellar::sim::Mission* findMissionById(std::span<const stellar::sim::Mission> missions,
                                             stellar::core::u64 id) {
  if (id == 0) return nullptr;
  for (const auto& m : missions) {
    if (m.id == id) return &m;
  }
  return nullptr;
}

const stellar::sim::Station* findStationById(const stellar::sim::StarSystem& sys, stellar::sim::StationId id) {
  if (id == 0) return nullptr;
  for (const auto& st : sys.stations) {
    if (st.id == id) return &st;
  }
  return nullptr;
}

std::string formatDeadline(double deadlineDay, double nowDay) {
  if (deadlineDay <= 0.0) return {};
  const double hrs = (deadlineDay - nowDay) * 24.0;
  if (!std::isfinite(hrs)) return {};
  char buf[64];
  if (hrs < 0.0) {
    std::snprintf(buf, sizeof(buf), "OVERDUE (%.1fh)", hrs);
  } else if (hrs < 48.0) {
    std::snprintf(buf, sizeof(buf), "%.1fh", hrs);
  } else {
    std::snprintf(buf, sizeof(buf), "%.1fd", hrs / 24.0);
  }
  return buf;
}

stellar::core::u64 computePlanCacheKey(const MissionControlWindowState& st,
                                       const MissionControlContext& ctx,
                                       std::span<const stellar::sim::Mission> missions) {
  stellar::core::u64 h = stellar::core::fnv1a64("MissionControlPlanV2");
  h = stellar::core::hashCombine(h, (stellar::core::u64)(ctx.currentSystem ? ctx.currentSystem->stub.id : 0));
  h = stellar::core::hashCombine(h, (stellar::core::u64)cacheStampFor(ctx.timeDays));

  h = stellar::core::hashCombine(h, (stellar::core::u64)st.routeMode);
  h = stellar::core::hashCombine(h, (stellar::core::u64)(st.groupBySystem ? 1 : 0));
  h = stellar::core::hashCombine(h, (stellar::core::u64)st.maxSystems);
  h = stellar::core::hashCombine(h, (stellar::core::u64)st.maxStops);

  // Quantize doubles for cache stability.
  auto q = [](double v, double s) -> stellar::core::u64 {
    if (!std::isfinite(v)) v = 0.0;
    return (stellar::core::u64)std::llround(v * s);
  };
  h = stellar::core::hashCombine(h, q(st.maxJumpLyOverride, 1000.0));
  h = stellar::core::hashCombine(h, q(st.queryRadiusLy, 100.0));
  h = stellar::core::hashCombine(h, q(st.rewardWeight, 1000.0));
  h = stellar::core::hashCombine(h, q(st.riskWeight, 1000.0));
  h = stellar::core::hashCombine(h, q(st.urgencyWeight, 1000.0));

  // ETA model.
  h = stellar::core::hashCombine(h, (stellar::core::u64)(st.etaAwareUrgency ? 1 : 0));
  h = stellar::core::hashCombine(h, q(st.etaSecondsPerJump, 10.0));
  h = stellar::core::hashCombine(h, q(st.etaSecondsPerLy, 100.0));
  h = stellar::core::hashCombine(h, q(st.etaSecondsPerStop, 10.0));
  h = stellar::core::hashCombine(h, q(st.etaSecondsPerSite, 10.0));

  // Selection.
  std::vector<stellar::core::u64> ids;
  ids.reserve(st.selectedMissionIds.size());
  for (auto id : st.selectedMissionIds) ids.push_back(id);
  std::sort(ids.begin(), ids.end());
  for (auto id : ids) {
    const auto* m = findMissionById(missions, id);
    if (!m) continue;
    h = stellar::core::hashCombine(h, id);
    h = stellar::core::hashCombine(h, (stellar::core::u64)(m->completed ? 1 : 0));
    h = stellar::core::hashCombine(h, (stellar::core::u64)(m->failed ? 1 : 0));
    h = stellar::core::hashCombine(h, (stellar::core::u64)m->leg);
    h = stellar::core::hashCombine(h, (stellar::core::u64)(m->scanned ? 1 : 0));
  }
  return h;
}

void drawFxWrapped(std::string_view markup, stellar::core::u64 seed, float fontScale, float timeRealSec) {
  ui::textfx::DrawParams dp;
  dp.baseColor = ImGui::GetColorU32(ImGuiCol_Text);
  dp.seed = seed;
  dp.wrapWidthPx = ImGui::GetContentRegionAvail().x;
  dp.fontSizePx = ImGui::GetFontSize() * std::max(0.75f, fontScale);

  ImDrawList* dl = ImGui::GetWindowDrawList();
  const ImVec2 pos = ImGui::GetCursorScreenPos();
  const ui::textfx::Program& prog = g_textFxCache.get(markup);
  const ImVec2 sz = ui::textfx::CalcSize(prog, dp);
  ui::textfx::Draw(dl, pos, prog, timeRealSec, dp);
  ImGui::Dummy(ImVec2(std::max(1.0f, sz.x), std::max(1.0f, sz.y)));
}

} // namespace

void tickMissionControl(MissionControlWindowState& st, const MissionControlContext& ctx) {
  if (!st.runner.active) return;
  if (!ctx.currentSystem || !ctx.missions) {
    st.runner.active = false;
    st.runner.stops.clear();
    return;
  }

  const auto missionsSpan = std::span<const stellar::sim::Mission>(*ctx.missions);

  // Emit a one-shot "armed" toast.
  if (!st.runner.startToastEmitted) {
    st.runner.startToastEmitted = true;
    if (ctx.toast && st.runner.stopIndex >= 0 && (std::size_t)st.runner.stopIndex < st.runner.stops.size()) {
      const auto& stop = st.runner.stops[(std::size_t)st.runner.stopIndex];
      const auto& sys = ctx.universe.getSystem(stop.objective.systemId);
      std::string msg = "Mission itinerary armed. Next: " + sys.stub.name;
      if (stop.objective.stationId != 0) {
        if (const auto* stn = findStationById(sys, stop.objective.stationId)) {
          msg += ":" + stn->name;
        }
      } else if (stop.objective.isSite) {
        msg += " (site)";
      }
      ctx.toast(msg, 2.2);
    }
  }

  if (st.runner.stopIndex < 0) st.runner.stopIndex = 0;
  if (st.runner.stops.empty() || (std::size_t)st.runner.stopIndex >= st.runner.stops.size()) {
    st.runner.active = false;
    st.runner.stops.clear();
    if (ctx.toast) ctx.toast("Mission itinerary complete.", 2.0);
    return;
  }

  // Replot on demand (only once per transition).
  if (st.runner.replotPending) {
    st.runner.replotPending = false;
    const auto& stop = st.runner.stops[(std::size_t)st.runner.stopIndex];

    if (stop.objective.stationId != 0) {
      if (ctx.goToStation) {
        ctx.goToStation(stop.objective.systemId, stop.objective.stationId, st.runner.armAutoRun);
      } else if (ctx.routeToStation) {
        ctx.routeToStation(stop.objective.systemId, stop.objective.stationId);
      } else if (ctx.plotRouteToSystem) {
        ctx.plotRouteToSystem(stop.objective.systemId);
      }
    } else {
      if (ctx.plotRouteToSystem) {
        ctx.plotRouteToSystem(stop.objective.systemId);
      }
    }
  }

  // Determine whether the current stop is still "active" for any of its missions.
  const auto& curStop = st.runner.stops[(std::size_t)st.runner.stopIndex];
  bool stillRelevant = false;
  for (auto id : curStop.missionIds) {
    const auto* m = findMissionById(missionsSpan, id);
    if (!m) continue;
    if (m->completed || m->failed) continue;
    const auto objNow = stellar::sim::missionNextObjective(*m);
    const bool same = (objNow.systemId == curStop.objective.systemId) &&
                      (objNow.stationId == curStop.objective.stationId) &&
                      (objNow.isSite == curStop.objective.isSite);
    if (same) {
      stillRelevant = true;
      break;
    }
  }

  if (stillRelevant) return;

  // Stop completed -> advance.
  st.runner.stopIndex += 1;
  st.runner.lastAdvanceAtDays = ctx.timeDays;

  if ((std::size_t)st.runner.stopIndex >= st.runner.stops.size()) {
    st.runner.active = false;
    st.runner.stops.clear();
    if (ctx.toast) ctx.toast("Mission itinerary complete.", 2.0);
    return;
  }

  st.runner.replotPending = st.runner.autoPlotNextStop;

  // Optionally auto-track a mission that belongs to the next stop.
  if (st.runner.autoTrackOnAdvance && ctx.trackMission) {
    const auto& nextStop = st.runner.stops[(std::size_t)st.runner.stopIndex];
    stellar::core::u64 best = 0;
    double bestDeadline = std::numeric_limits<double>::infinity();
    for (auto id : nextStop.missionIds) {
      const auto* m = findMissionById(missionsSpan, id);
      if (!m) continue;
      if (m->completed || m->failed) continue;
      const double hrs = (m->deadlineDay > 0.0) ? (m->deadlineDay - ctx.timeDays) * 24.0 : std::numeric_limits<double>::infinity();
      if (best == 0 || hrs < bestDeadline) {
        best = id;
        bestDeadline = hrs;
      }
    }
    if (best != 0) ctx.trackMission(best);
  }
}

void drawMissionControlWindow(MissionControlWindowState& st, const MissionControlContext& ctx) {
  // Keep runner alive even if the UI is open.
  tickMissionControl(st, ctx);

  if (!st.open) return;

  if (ImGui::Begin("Mission Control", &st.open)) {
    if (!ctx.currentSystem || !ctx.missions) {
      ImGui::TextDisabled("No active universe context.");
      ImGui::End();
      return;
    }

    const auto missionsSpan = std::span<const stellar::sim::Mission>(*ctx.missions);

    // Initialize selection to all active missions.
    if (!st.selectionInitialized) {
      st.selectedMissionIds.clear();
      for (const auto& m : missionsSpan) {
        if (m.completed || m.failed) continue;
        st.selectedMissionIds.insert(m.id);
      }
      st.selectionInitialized = true;
      if (st.focusedMissionId == 0) st.focusedMissionId = ctx.trackedMissionId;
    }

    // Prune stale selection ids.
    {
      std::unordered_set<stellar::core::u64> present;
      present.reserve(missionsSpan.size());
      for (const auto& m : missionsSpan) present.insert(m.id);
      for (auto it = st.selectedMissionIds.begin(); it != st.selectedMissionIds.end();) {
        if (present.find(*it) == present.end()) it = st.selectedMissionIds.erase(it);
        else ++it;
      }
    }

    // --- Controls / Planner settings ---
    if (ImGui::CollapsingHeader("Planner", ImGuiTreeNodeFlags_DefaultOpen)) {
      auto sliderDouble = [](const char* label, double* v, double vMin, double vMax, const char* fmt) -> bool {
        double mn = vMin;
        double mx = vMax;
        return ImGui::SliderScalar(label, ImGuiDataType_Double, v, &mn, &mx, fmt);
      };

      ImGui::Checkbox("Group by system", &st.groupBySystem);
      ImGui::SameLine();
      ImGui::Checkbox("Auto rebuild", &st.autoRebuild);

      ImGui::SetNextItemWidth(180);
      ImGui::Combo("Route", &st.routeMode,
                   "Min hops\0"
                   "Min distance\0\0");

      ImGui::SetNextItemWidth(180);
      ImGui::SliderInt("Max stops", &st.maxStops, 1, 32);
      ImGui::SetNextItemWidth(180);
      ImGui::SliderInt("Max systems", &st.maxSystems, 64, 4096);

      const double defaultJump = std::max(0.0, ctx.maxJumpLy);
      ImGui::SetNextItemWidth(180);
      sliderDouble("Jump range (ly)", &st.maxJumpLyOverride, 0.0, std::max(5.0, defaultJump * 1.5), "%.1f");
      if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("0 uses current ship jump range (%.1f ly)", defaultJump);
      }

      ImGui::SetNextItemWidth(180);
      sliderDouble("Batch radius (ly)", &st.queryRadiusLy, 0.0, 2400.0, "%.0f");
      if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("0 picks an automatic radius based on objectives.");
      }

      ImGui::Separator();
      ImGui::TextDisabled("Score weights");
      ImGui::SetNextItemWidth(220);
      sliderDouble("Reward", &st.rewardWeight, 0.0, 3.0, "%.2f");
      ImGui::SetNextItemWidth(220);
      sliderDouble("Risk penalty", &st.riskWeight, 0.0, 2.0, "%.2f");
      ImGui::SetNextItemWidth(220);
      sliderDouble("Urgency", &st.urgencyWeight, 0.0, 2.0, "%.2f");

      ImGui::Separator();
      ImGui::TextDisabled("ETA model (deadline-aware urgency + analytics)");
      ImGui::Checkbox("ETA-aware urgency", &st.etaAwareUrgency);
      if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("When enabled, the planner discounts mission deadlines by an estimated travel time.\n"
                          "This reduces plans that look good on paper but are impossible to finish before expiry.");
      }
      ImGui::SetNextItemWidth(220);
      sliderDouble("ETA sec / jump", &st.etaSecondsPerJump, 0.0, 900.0, "%.0f");
      ImGui::SetNextItemWidth(220);
      sliderDouble("ETA sec / ly", &st.etaSecondsPerLy, 0.0, 120.0, "%.1f");
      ImGui::SetNextItemWidth(220);
      sliderDouble("ETA sec / stop", &st.etaSecondsPerStop, 0.0, 7200.0, "%.0f");
      ImGui::SetNextItemWidth(220);
      sliderDouble("ETA sec / site", &st.etaSecondsPerSite, 0.0, 7200.0, "%.0f");

      if (ImGui::Button("Rebuild plan")) {
        st.cacheKey = 0;
      }
    }

    // --- Mission selection + details ---
    if (ImGui::CollapsingHeader("Missions", ImGuiTreeNodeFlags_DefaultOpen)) {
      if (ImGui::Button("Select all active")) {
        st.selectedMissionIds.clear();
        for (const auto& m : missionsSpan) {
          if (m.completed || m.failed) continue;
          st.selectedMissionIds.insert(m.id);
        }
      }
      ImGui::SameLine();
      if (ImGui::Button("Select none")) {
        st.selectedMissionIds.clear();
      }
      ImGui::SameLine();
      if (ImGui::Button("Select tracked") && ctx.trackedMissionId != 0) {
        st.selectedMissionIds.clear();
        st.selectedMissionIds.insert(ctx.trackedMissionId);
      }

      ImGui::SameLine();
      ImGui::Checkbox("Show completed", &st.includeCompleted);
      ImGui::SameLine();
      ImGui::Checkbox("Show failed", &st.includeFailed);

      if (ImGui::BeginTable("mission_list", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable)) {
        ImGui::TableSetupColumn("Use", ImGuiTableColumnFlags_WidthFixed, 36);
        ImGui::TableSetupColumn("Type", ImGuiTableColumnFlags_WidthFixed, 110);
        ImGui::TableSetupColumn("Reward", ImGuiTableColumnFlags_WidthFixed, 70);
        ImGui::TableSetupColumn("Deadline", ImGuiTableColumnFlags_WidthFixed, 70);
        ImGui::TableSetupColumn("Objective", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Track", ImGuiTableColumnFlags_WidthFixed, 46);
        ImGui::TableHeadersRow();

        for (const auto& m : missionsSpan) {
          if (!st.includeCompleted && m.completed) continue;
          if (!st.includeFailed && m.failed) continue;

          ImGui::TableNextRow();

          const bool selected = st.selectedMissionIds.find(m.id) != st.selectedMissionIds.end();

          ImGui::TableNextColumn();
          bool use = selected;
          ImGui::PushID((int)m.id);
          if (ImGui::Checkbox("##use", &use)) {
            if (use) st.selectedMissionIds.insert(m.id);
            else st.selectedMissionIds.erase(m.id);
          }

          ImGui::TableNextColumn();
          ImGui::TextUnformatted(missionTypeLabel(m.type));

          ImGui::TableNextColumn();
          ImGui::Text("%.0f", std::max(0.0, m.reward));

          ImGui::TableNextColumn();
          const auto deadline = formatDeadline(m.deadlineDay, ctx.timeDays);
          if (!deadline.empty()) ImGui::TextUnformatted(deadline.c_str());
          else ImGui::TextDisabled("-");

          ImGui::TableNextColumn();
          const auto obj = stellar::sim::missionNextObjective(m);
          const auto& sys = ctx.universe.getSystem(obj.systemId);
          std::string label = sys.stub.name;
          if (obj.stationId != 0) {
            if (const auto* stn = findStationById(sys, obj.stationId)) {
              label += ":" + stn->name;
            }
          } else if (obj.isSite) {
            label += " (site)";
          }
          const bool focused = (st.focusedMissionId == m.id);
          if (ImGui::Selectable(label.c_str(), focused, ImGuiSelectableFlags_SpanAllColumns)) {
            st.focusedMissionId = m.id;
          }

          ImGui::TableNextColumn();
          if (ctx.trackMission) {
            const bool isTracked = (ctx.trackedMissionId == m.id);
            ImGui::BeginDisabled(isTracked);
            if (ImGui::SmallButton("Track")) {
              ctx.trackMission(m.id);
            }
            ImGui::EndDisabled();
          } else {
            ImGui::TextDisabled("-");
          }

          ImGui::PopID();
        }

        ImGui::EndTable();
      }

      // Focused mission briefing.
      const auto* fm = findMissionById(missionsSpan, st.focusedMissionId);
      if (fm) {
        ImGui::Separator();
        ImGui::TextDisabled("Focused mission");
        ImGui::SameLine();
        ImGui::Text("%s", missionTypeLabel(fm->type));
        ImGui::SameLine();
        ImGui::TextDisabled("(id %llu)", (unsigned long long)fm->id);

        // Origin station context.
        const stellar::sim::StarSystem* originSys = nullptr;
        const stellar::sim::Station* originSt = nullptr;
        if (fm->fromSystem != 0) {
          originSys = &ctx.universe.getSystem(fm->fromSystem);
          originSt = findStationById(*originSys, fm->fromStation);
          if (!originSt && originSys && !originSys->stations.empty()) originSt = &originSys->stations.front();
        }
        if (!originSys) originSys = ctx.currentSystem;
        stellar::sim::Station dummy{};
        if (!originSt) {
          if (originSys && !originSys->stations.empty()) originSt = &originSys->stations.front();
          else originSt = &dummy;
        }

        // Faction rep lookup.
        double rep = 0.0;
        for (const auto& r : ctx.playerRepWithFaction) {
            if (r.factionId == fm->factionId) { rep = r.rep; break; }
        }

        stellar::sim::MissionBriefingParams bp;
        bp.useMarkup = true;
        bp.includeRiskHints = true;
        const auto brief = stellar::sim::generateMissionBriefing(ctx.universe,
                                                                 *originSys,
                                                                 *originSt,
                                                                 ctx.timeDays,
                                                                 rep,
                                                                 *fm,
                                                                 ctx.securityDeltas,
                                                                 bp);

        const stellar::core::u64 seed = stellar::core::hashCombine(stellar::core::fnv1a64("mission_control_brief"), (stellar::core::u64)fm->id);
        drawFxWrapped(brief.titleMarkup, seed, 1.25f, (float)ctx.timeRealSec);
        ImGui::Spacing();
        drawFxWrapped(brief.synopsisMarkup, seed + 1, 1.0f, (float)ctx.timeRealSec);
        ImGui::Spacing();
        for (std::size_t i = 0; i < brief.bulletsMarkup.size(); ++i) {
          drawFxWrapped(brief.bulletsMarkup[i], seed + 2 + (stellar::core::u64)i, 0.95f, (float)ctx.timeRealSec);
        }

        ImGui::Separator();
        ImGui::TextDisabled("Risk");
        ImGui::ProgressBar((float)brief.risk.overall01, ImVec2(-1, 0), stellar::sim::riskTierName(brief.risk.overall01));
        if (brief.risk.lawRisk01 > 1e-3) {
          ImGui::ProgressBar((float)brief.risk.lawRisk01, ImVec2(-1, 0), "Law");
        }
        if (brief.risk.combatRisk01 > 1e-3) {
          ImGui::ProgressBar((float)brief.risk.combatRisk01, ImVec2(-1, 0), "Combat");
        }
        if (brief.risk.danger01 > 1e-3) {
          ImGui::ProgressBar((float)brief.risk.danger01, ImVec2(-1, 0), "Danger");
        }
      }
    }

    // --- Itinerary ---
    if (ImGui::CollapsingHeader("Itinerary", ImGuiTreeNodeFlags_DefaultOpen)) {
      // Build a filtered mission set.
      std::vector<stellar::sim::Mission> selected;
      selected.reserve(st.selectedMissionIds.size());
      for (const auto& m : missionsSpan) {
        if (st.selectedMissionIds.find(m.id) == st.selectedMissionIds.end()) continue;
        if (!st.includeCompleted && m.completed) continue;
        if (!st.includeFailed && m.failed) continue;
        selected.push_back(m);
      }

      // Compute plan if needed.
      const stellar::core::u64 key = computePlanCacheKey(st, ctx, missionsSpan);
      if (st.autoRebuild && key != st.cacheKey) {
        st.cacheKey = key;

        stellar::sim::MissionItineraryParams pp;
        pp.groupBySystem = st.groupBySystem;
        pp.maxSystems = (std::size_t)std::max(1, st.maxSystems);
        pp.maxStops = std::clamp(st.maxStops, 1, 64);
        pp.maxJumpLy = (st.maxJumpLyOverride > 1e-6) ? st.maxJumpLyOverride : std::max(0.0, ctx.maxJumpLy);
        pp.queryRadiusLy = st.queryRadiusLy;
        pp.rewardWeight = st.rewardWeight;
        pp.riskWeight = st.riskWeight;
        pp.urgencyWeight = st.urgencyWeight;

        pp.etaAwareUrgency = st.etaAwareUrgency;
        pp.etaSecondsPerJump = st.etaSecondsPerJump;
        pp.etaSecondsPerLy = st.etaSecondsPerLy;
        pp.etaSecondsPerStop = st.etaSecondsPerStop;
        pp.etaSecondsPerSite = st.etaSecondsPerSite;

        if (st.routeMode == 1) {
          pp.costPerJump = 0.0;
          pp.costPerLy = 1.0;
        } else {
          pp.costPerJump = 1.0;
          pp.costPerLy = 0.0;
        }

        st.plan = stellar::sim::planMissionItinerary(ctx.universe,
                                                    *ctx.currentSystem,
                                                    ctx.timeDays,
                                                    selected,
                                                    ctx.playerRepWithFaction,
                                                    ctx.securityDeltas,
                                                    pp);
      }

      if (selected.empty()) {
        ImGui::TextDisabled("No missions selected.");
      } else if (st.plan.stops.empty()) {
        ImGui::TextDisabled("No itinerary stops (all selected missions may already be complete).\nTry selecting active missions.");
      } else {
        std::string etaEnd;
        if (!st.plan.stops.empty() && std::isfinite(st.plan.stops.back().etaDay)) {
          const double etaHrs = (st.plan.stops.back().etaDay - ctx.timeDays) * 24.0;
          char buf[48];
          if (etaHrs < 48.0) std::snprintf(buf, sizeof(buf), "+%.1fh", std::max(0.0, etaHrs));
          else std::snprintf(buf, sizeof(buf), "+%.1fd", std::max(0.0, etaHrs) / 24.0);
          etaEnd = buf;
        } else {
          etaEnd = "-";
        }

        ImGui::TextDisabled("Stops: %d | Cost: %.1f | Hops: %d | Dist: %.1fly | ETA: %s",
                            (int)st.plan.stops.size(),
                            st.plan.totalCost,
                            st.plan.totalHops,
                            st.plan.totalDistanceLy,
                            etaEnd.c_str());

        if (st.plan.unreachableStops > 0) {
          ImGui::SameLine();
          ImGui::TextColored(ImVec4(1, 0.6f, 0.3f, 1), "(unreachable: %d)", st.plan.unreachableStops);
        }

        // Runner controls.
        {
          const bool canArm = !st.groupBySystem;
          ImGui::BeginDisabled(!canArm);
          if (!st.runner.active) {
            if (ImGui::Button("Arm runner")) {
              st.runner.active = true;
              st.runner.armAutoRun = true;
              st.runner.autoPlotNextStop = true;
              st.runner.autoTrackOnAdvance = true;
              st.runner.stopIndex = 0;
              st.runner.startedAtDays = ctx.timeDays;
              st.runner.lastAdvanceAtDays = ctx.timeDays;
              st.runner.startToastEmitted = false;
              st.runner.replotPending = true;
              st.runner.stops = st.plan.stops;
            }
          } else {
            if (ImGui::Button("Cancel runner")) {
              st.runner.active = false;
              st.runner.stops.clear();
              if (ctx.toast) ctx.toast("Mission itinerary cancelled.", 2.0);
            }
          }
          ImGui::EndDisabled();
          if (!canArm && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
            ImGui::SetTooltip("Runner requires 'Group by system' to be OFF (needs concrete station/site stops).");
          }

          if (st.runner.active) {
            ImGui::SameLine();
            ImGui::TextDisabled("Runner: stop %d/%d", st.runner.stopIndex + 1, (int)st.runner.stops.size());
          }
        }

        if (ImGui::BeginTable("itinerary", 8, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable)) {
          ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 26);
          ImGui::TableSetupColumn("Destination", ImGuiTableColumnFlags_WidthStretch);
          ImGui::TableSetupColumn("Missions", ImGuiTableColumnFlags_WidthFixed, 60);
          ImGui::TableSetupColumn("Reward", ImGuiTableColumnFlags_WidthFixed, 70);
          ImGui::TableSetupColumn("Risk", ImGuiTableColumnFlags_WidthFixed, 70);
          ImGui::TableSetupColumn("Travel", ImGuiTableColumnFlags_WidthFixed, 70);
          ImGui::TableSetupColumn("ETA", ImGuiTableColumnFlags_WidthFixed, 86);
          ImGui::TableSetupColumn("Go", ImGuiTableColumnFlags_WidthFixed, 44);
          ImGui::TableHeadersRow();

          for (std::size_t i = 0; i < st.plan.stops.size(); ++i) {
            const auto& stop = st.plan.stops[i];
            const auto& sys = ctx.universe.getSystem(stop.objective.systemId);

            std::string dst = sys.stub.name;
            if (stop.objective.stationId != 0) {
              if (const auto* stn = findStationById(sys, stop.objective.stationId)) {
                dst += ":" + stn->name;
              }
            } else if (stop.objective.isSite) {
              dst += " (site)";
            }

            ImGui::TableNextRow();

            ImGui::TableNextColumn();
            ImGui::Text("%d", (int)i + 1);

            ImGui::TableNextColumn();
            ImGui::TextUnformatted(dst.c_str());
            if (!stop.reachable) {
              ImGui::SameLine();
              ImGui::TextColored(ImVec4(1, 0.6f, 0.3f, 1), "(est)");
            }
            if (stop.earliestDeadlineDay > 0.0) {
              ImGui::SameLine();
              const auto d = formatDeadline(stop.earliestDeadlineDay, ctx.timeDays);
              if (!d.empty()) {
                ImGui::TextDisabled("[%s]", d.c_str());
              }
            }

            ImGui::TableNextColumn();
            ImGui::Text("%d", (int)stop.missionIds.size());

            ImGui::TableNextColumn();
            ImGui::Text("%.0f", std::max(0.0, stop.rewardCr));

            ImGui::TableNextColumn();
            ImGui::Text("%s", stellar::sim::riskTierName(stop.avgRisk01));

            ImGui::TableNextColumn();
            if (st.routeMode == 0) {
              ImGui::Text("%d hops", stop.hopsFromPrev);
            } else {
              ImGui::Text("%.0fly", stop.distanceLyFromPrev);
            }

            ImGui::TableNextColumn();
            {
              const double etaHrs = (stop.etaDay - ctx.timeDays) * 24.0;
              if (!std::isfinite(etaHrs)) {
                ImGui::TextDisabled("-");
              } else {
                char buf[48];
                if (etaHrs < 48.0) std::snprintf(buf, sizeof(buf), "+%.1fh", std::max(0.0, etaHrs));
                else std::snprintf(buf, sizeof(buf), "+%.1fd", std::max(0.0, etaHrs) / 24.0);
                ImGui::TextUnformatted(buf);

                if (stop.earliestDeadlineDay > 0.0 && std::isfinite(stop.etaSlackHours)) {
                  char b2[64];
                  if (stop.etaSlackHours < 0.0) {
                    std::snprintf(b2, sizeof(b2), "late %.1fh", -stop.etaSlackHours);
                    ImGui::TextColored(ImVec4(1, 0.25f, 0.25f, 1), "%s", b2);
                  } else if (stop.etaSlackHours < 6.0) {
                    std::snprintf(b2, sizeof(b2), "slack %.1fh", stop.etaSlackHours);
                    ImGui::TextColored(ImVec4(1, 0.85f, 0.25f, 1), "%s", b2);
                  } else {
                    std::snprintf(b2, sizeof(b2), "slack %.1fh", stop.etaSlackHours);
                    ImGui::TextDisabled("%s", b2);
                  }
                }
              }
            }

            ImGui::TableNextColumn();
            ImGui::PushID((int)i);
            if (ImGui::SmallButton("Go")) {
              if (stop.objective.stationId != 0) {
                if (ctx.goToStation) ctx.goToStation(stop.objective.systemId, stop.objective.stationId, /*armAutoRun=*/false);
                else if (ctx.routeToStation) ctx.routeToStation(stop.objective.systemId, stop.objective.stationId);
                else if (ctx.plotRouteToSystem) ctx.plotRouteToSystem(stop.objective.systemId);
              } else {
                if (ctx.plotRouteToSystem) ctx.plotRouteToSystem(stop.objective.systemId);
              }
            }
            ImGui::PopID();

            // Optional per-stop tooltip.
            if (ImGui::IsItemHovered()) {
              ImGui::BeginTooltip();
              ImGui::TextDisabled("Stop %d", (int)i + 1);
              ImGui::Text("Missions: %d", (int)stop.missionIds.size());
              ImGui::Text("Reward: %.0f", stop.rewardCr);
              ImGui::Text("Risk: %.0f%% (%s)", stop.avgRisk01 * 100.0, stellar::sim::riskTierName(stop.avgRisk01));
              if (std::isfinite(stop.etaDay)) {
                const double etaHrs = (stop.etaDay - ctx.timeDays) * 24.0;
                if (std::isfinite(etaHrs)) {
                  ImGui::Text("ETA: +%.1fh", std::max(0.0, etaHrs));
                }

                if (stop.earliestDeadlineDay > 0.0 && std::isfinite(stop.etaSlackHours)) {
                  ImGui::Text("Deadline slack: %.1fh", stop.etaSlackHours);
                }
              }
              if (stop.earliestDeadlineDay > 0.0) {
                ImGui::Text("Deadline: %s", formatDeadline(stop.earliestDeadlineDay, ctx.timeDays).c_str());
              }
              ImGui::Separator();
              ImGui::TextDisabled("ETA model: %.0fs/jump, %.1fs/ly, %.0fs/stop, %.0fs/site", st.etaSecondsPerJump, st.etaSecondsPerLy, st.etaSecondsPerStop, st.etaSecondsPerSite);
              ImGui::EndTooltip();
            }
          }

          ImGui::EndTable();
        }
      }
    }
  }
  ImGui::End();
}

} // namespace stellar::game
