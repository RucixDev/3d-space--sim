#include "CopilotWindow.h"

#include "stellar/econ/Commodity.h"
#include "stellar/sim/Mission.h"
#include "stellar/sim/MissionAssist.h"
#include "stellar/sim/System.h"
#include "stellar/sim/Units.h"
#include "stellar/sim/Universe.h"

#include <imgui.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace stellar::game {

namespace {

enum class Severity : int {
  Info = 0,
  Warn = 1,
  Critical = 2,
};

struct SuggestionAction {
  std::string label;
  std::string tooltip;
  std::function<void()> fn;
  bool disabled{false};
};

struct Suggestion {
  std::string id;
  Severity severity{Severity::Info};
  std::string title;
  std::string body;
  std::vector<SuggestionAction> actions;
};

static const char* severityGlyph(Severity s) {
  switch (s) {
    case Severity::Info: return "i";
    case Severity::Warn: return "!";
    case Severity::Critical: return "!!";
  }
  return "?";
}

static bool containsCaseInsensitive(std::string_view hay, std::string_view needle) {
  if (needle.empty()) return true;
  if (hay.size() < needle.size()) return false;

  auto lower = [](unsigned char c) -> unsigned char {
    return (unsigned char)std::tolower(c);
  };

  for (std::size_t i = 0; i + needle.size() <= hay.size(); ++i) {
    bool ok = true;
    for (std::size_t j = 0; j < needle.size(); ++j) {
      if (lower((unsigned char)hay[i + j]) != lower((unsigned char)needle[j])) {
        ok = false;
        break;
      }
    }
    if (ok) return true;
  }
  return false;
}

static std::pair<sim::SystemId, sim::StationId> missionNextStop(const sim::Mission& m) {
  // Multi-leg deliveries: go to the via stop first.
  const bool hasViaStop = (m.viaSystem != 0 && m.viaStation != 0 && m.leg == 0
                           && (m.type == sim::MissionType::MultiDelivery || m.type == sim::MissionType::Passenger));
  if (hasViaStop) {
    return {m.viaSystem, m.viaStation};
  }

  // Salvage jobs: before the site is visited, the objective is the mission site (not a station).
  if (m.type == sim::MissionType::Salvage && !m.scanned) {
    return {m.toSystem, 0};
  }

  // Default: the destination station (if any).
  return {m.toSystem, m.toStation};
}

static std::string systemName(sim::Universe& universe, sim::SystemId systemId) {
  if (systemId == 0) return "?";
  return universe.getSystem(systemId).stub.name;
}

static std::string stationName(sim::Universe& universe, sim::SystemId systemId, sim::StationId stationId) {
  if (systemId == 0 || stationId == 0) return "?";
  const auto& sys = universe.getSystem(systemId);
  for (const auto& st : sys.stations) {
    if (st.id == stationId) return st.name;
  }
  return "?";
}

struct NearestStation {
  sim::StationId id{0};
  double distKm{0.0};
};

static NearestStation nearestStationInSystem(const sim::StarSystem& sys, const math::Vec3d& shipPosKm, double timeDays) {
  NearestStation out{};
  double best = std::numeric_limits<double>::infinity();

  for (const auto& st : sys.stations) {
    const math::Vec3d stPos = sim::stationPosKm(st, timeDays);
    const double d = (shipPosKm - stPos).length();
    if (d < best) {
      best = d;
      out.id = st.id;
      out.distKm = d;
    }
  }

  return out;
}

static std::string formatDeadline(const sim::Mission& m, double timeDays) {
  if (m.deadlineDay <= 0.0) return "—";

  const double hrsLeft = (m.deadlineDay - timeDays) * 24.0;
  if (!std::isfinite(hrsLeft)) return "?";

  if (hrsLeft < 0.0) {
    return "OVERDUE";
  }
  if (hrsLeft < 48.0) {
    const int h = (int)std::ceil(hrsLeft);
    return std::to_string(std::max(0, h)) + " h";
  }

  const int d = (int)std::ceil(hrsLeft / 24.0);
  return std::to_string(std::max(0, d)) + " d";
}

static Severity missionDeadlineSeverity(const sim::Mission& m, double timeDays) {
  if (m.deadlineDay <= 0.0) return Severity::Info;
  const double hrsLeft = (m.deadlineDay - timeDays) * 24.0;
  if (!std::isfinite(hrsLeft)) return Severity::Warn;
  if (hrsLeft <= 0.0) return Severity::Critical;
  if (hrsLeft <= 6.0) return Severity::Critical;
  if (hrsLeft <= 24.0) return Severity::Warn;
  return Severity::Info;
}

static void addAction(Suggestion& s, std::string label, std::string tooltip, std::function<void()> fn, bool disabled = false) {
  if (!fn) return;
  SuggestionAction a;
  a.label = std::move(label);
  a.tooltip = std::move(tooltip);
  a.fn = std::move(fn);
  a.disabled = disabled;
  s.actions.push_back(std::move(a));
}

static void addShipRecommendations(std::vector<Suggestion>& out, const CopilotContext& ctx) {
  const double hullMax = std::max(0.0, ctx.hullMax);
  const double shieldMax = std::max(0.0, ctx.shieldMax);
  const double fuelMax = std::max(0.0, ctx.fuelMax);

  // Hull
  if (hullMax > 1e-6) {
    const double hullFrac = std::clamp(ctx.hull / hullMax, 0.0, 1.0);
    if (hullFrac < 0.40 && (hullMax - ctx.hull) > 0.5) {
      Suggestion s;
      s.id = "ship.hull.low";
      s.severity = (hullFrac < 0.20) ? Severity::Critical : Severity::Warn;
      s.title = "Hull low (" + std::to_string((int)std::round(hullFrac * 100.0)) + "%)";
      s.body = ctx.docked ? "Docked: open Services to repair." : "Consider docking to repair before taking risks.";

      if (ctx.docked) {
        addAction(s,
                  "Services",
                  "Open station Services (repairs/refuel).",
                  ctx.openStationServices);
      } else if (ctx.currentSystem && ctx.goToStation) {
        const auto ns = nearestStationInSystem(*ctx.currentSystem, ctx.shipPosKm, ctx.timeDays);
        if (ns.id != 0) {
          const sim::SystemId sysId = ctx.currentSystem->stub.id;
          addAction(s,
                    "Go nearest",
                    "Arm auto-run to the nearest station.",
                    [cb = ctx.goToStation, sysId, stId = ns.id]() { cb(sysId, stId, true); });
        }
      }

      out.push_back(std::move(s));
    }
  }

  // Shield
  if (shieldMax > 1e-6) {
    const double shFrac = std::clamp(ctx.shield / shieldMax, 0.0, 1.0);
    if (!ctx.docked && shFrac < 0.15 && shieldMax > 1.0) {
      Suggestion s;
      s.id = "ship.shield.low";
      s.severity = (shFrac < 0.06) ? Severity::Critical : Severity::Warn;
      s.title = "Shields low (" + std::to_string((int)std::round(shFrac * 100.0)) + "%)";
      s.body = "Break line-of-sight and let shields regenerate.";
      out.push_back(std::move(s));
    }
  }

  // Fuel
  if (fuelMax > 1e-6) {
    const double fFrac = std::clamp(ctx.fuel / fuelMax, 0.0, 1.0);
    if (fFrac < 0.25) {
      Suggestion s;
      s.id = "ship.fuel.low";
      s.severity = (fFrac < 0.10) ? Severity::Critical : Severity::Warn;
      s.title = "Fuel low (" + std::to_string((int)std::round(fFrac * 100.0)) + "%)";
      s.body = ctx.docked ? "Docked: open Services to refuel." : "Consider docking to refuel before long trips.";

      if (ctx.docked) {
        addAction(s,
                  "Services",
                  "Open station Services (repairs/refuel).",
                  ctx.openStationServices);
      } else if (ctx.currentSystem && ctx.goToStation) {
        const auto ns = nearestStationInSystem(*ctx.currentSystem, ctx.shipPosKm, ctx.timeDays);
        if (ns.id != 0) {
          const sim::SystemId sysId = ctx.currentSystem->stub.id;
          addAction(s,
                    "Go nearest",
                    "Arm auto-run to the nearest station.",
                    [cb = ctx.goToStation, sysId, stId = ns.id]() { cb(sysId, stId, true); });
        }
      }

      if (ctx.fuelScoopMk > 0) {
        addAction(s,
                  "Fuel scoop",
                  "You have a fuel scoop installed. In some ships you can refuel near stars.",
                  [t = ctx.toast]() {
                    if (t) t("Fuel scoop installed: skim a star (carefully) to refuel.", 3.0);
                  });
      }

      out.push_back(std::move(s));
    }
  }

  // Heat
  if (std::isfinite(ctx.heat) && ctx.heat > 92.0) {
    Suggestion s;
    s.id = "ship.heat.high";
    s.severity = (ctx.heat > 98.0) ? Severity::Critical : Severity::Warn;
    s.title = "Heat high (" + std::to_string((int)std::round(ctx.heat)) + "%)";
    s.body = "Throttle down and avoid boosting; heat increases signature and can cause damage.";
    out.push_back(std::move(s));
  }
}

static const sim::Mission* findTrackedMission(const CopilotContext& ctx) {
  if (!ctx.missions || ctx.trackedMissionId == 0) return nullptr;
  for (const auto& m : *ctx.missions) {
    if (m.id == ctx.trackedMissionId && !m.completed && !m.failed) return &m;
  }
  return nullptr;
}

static const sim::Mission* selectMostUrgentMission(const std::vector<sim::Mission>& missions, double timeDays) {
  const sim::Mission* best = nullptr;
  double bestScore = std::numeric_limits<double>::infinity();

  for (const auto& m : missions) {
    if (m.completed || m.failed) continue;

    double score = 1e9;
    if (m.deadlineDay > 0.0) {
      score = m.deadlineDay - timeDays;
    }

    if (score < bestScore) {
      bestScore = score;
      best = &m;
    }
  }

  // Fallback: first active mission.
  if (!best) {
    for (const auto& m : missions) {
      if (!m.completed && !m.failed) return &m;
    }
  }

  return best;
}

static void addMissionRecommendations(std::vector<Suggestion>& out, const CopilotContext& ctx) {
  if (!ctx.missions) return;

  // Track something if nothing is tracked.
  const sim::Mission* tracked = findTrackedMission(ctx);
  if (!tracked) {
    const sim::Mission* urgent = selectMostUrgentMission(*ctx.missions, ctx.timeDays);
    if (urgent) {
      Suggestion s;
      s.id = "mission.track";
      s.severity = missionDeadlineSeverity(*urgent, ctx.timeDays);
      s.title = "Track a mission";
      s.body = std::string("Most urgent: ") + sim::missionTypeName(urgent->type) + " (" + formatDeadline(*urgent, ctx.timeDays) + ")";

      addAction(s,
                "Track",
                "Set this mission as tracked (Objective HUD + recommendations).",
                [cb = ctx.trackMission, id = urgent->id]() {
                  if (cb) cb(id);
                },
                !ctx.trackMission);

      addAction(s,
                "Missions",
                "Open the Missions window.",
                ctx.openMissionsWindow);

      out.push_back(std::move(s));
    }

    return;
  }

  // Tracked mission summary.
  {
    const auto [nextSys, nextSt] = missionNextStop(*tracked);

    Suggestion s;
    s.id = "mission.tracked";
    s.severity = missionDeadlineSeverity(*tracked, ctx.timeDays);

    std::string dest = systemName(ctx.universe, nextSys);
    if (nextSt != 0) {
      dest += " / " + stationName(ctx.universe, nextSys, nextSt);
    } else {
      dest += " (site)";
    }

    s.title = std::string("Tracked: ") + sim::missionTypeName(tracked->type) + " → " + dest;
    s.body = "Deadline " + formatDeadline(*tracked, ctx.timeDays) + " | Reward " + std::to_string((int)std::round(tracked->reward)) + " cr";

    // Travel actions.
    if (nextSys != 0) {
      if (!ctx.currentSystem || nextSys != ctx.currentSystem->stub.id) {
        // Cross-system.
        addAction(s,
                  "Plot",
                  "Plot an FSD route to the destination system.",
                  [cb = ctx.plotRouteToSystem, sysId = nextSys]() {
                    if (cb) cb(sysId);
                  },
                  !ctx.plotRouteToSystem);

        addAction(s,
                  "Go (auto)",
                  "Arm auto-run to travel and dock.",
                  [cb = ctx.goToStation, sysId = nextSys, stId = nextSt]() {
                    if (cb) cb(sysId, stId, true);
                  },
                  !ctx.goToStation);
      } else {
        // In-system.
        if (nextSt != 0) {
          addAction(s,
                    "Target",
                    "Target the destination station.",
                    [cb = ctx.targetStation, stId = nextSt]() {
                      if (cb) cb(stId);
                    },
                    !ctx.targetStation);

          addAction(s,
                    "Dock",
                    "Request docking clearance at the destination station.",
                    [cb = ctx.requestDocking, stId = nextSt]() {
                      if (cb) cb(stId);
                    },
                    !ctx.requestDocking);

          addAction(s,
                    "Auto-dock",
                    "Engage the docking computer/autopilot (requires a station target).",
                    ctx.engageDockingComputer);
        }
      }
    }

    addAction(s,
              "Missions",
              "Open the Missions window.",
              ctx.openMissionsWindow);

    out.push_back(std::move(s));
  }

  // Cargo sourcing for delivery/smuggle missions.
  if (ctx.currentSystem && ctx.cargo) {
    const bool cargoMission = (tracked->type == sim::MissionType::Delivery
                               || tracked->type == sim::MissionType::MultiDelivery
                               || tracked->type == sim::MissionType::Smuggle);

    if (cargoMission && tracked->units > 0.0) {
      const std::size_t ci = (std::size_t)tracked->commodity;
      const double have = (*ctx.cargo)[ci];
      const double need = tracked->units;
      const double missing = std::max(0.0, need - have);

      if (missing > 0.5) {
        sim::MissionCargoSourceParams p{};
        p.searchRadiusLy = 160.0;
        p.maxResults = 6;
        p.includeCurrentSystem = true;

        const auto plan = sim::planMissionCargoSourcingForMission(ctx.universe,
                                                                  *ctx.currentSystem,
                                                                  ctx.timeDays,
                                                                  *tracked,
                                                                  *ctx.cargo,
                                                                  p);

        Suggestion s;
        s.id = "mission.cargo." + std::to_string((unsigned long long)tracked->id);
        s.severity = Severity::Warn;
        const auto& def = econ::commodityDef(tracked->commodity);
        s.title = "Missing " + std::to_string((int)std::round(missing)) + " units of " + std::string(def.name);

        if (plan.ok && !plan.candidates.empty()) {
          const auto& c = plan.candidates.front();
          const std::string sysName = systemName(ctx.universe, c.systemId);
          const std::string stName = stationName(ctx.universe, c.systemId, c.stationId);

          s.body = "Best supplier: " + sysName + " / " + stName + " (" + std::to_string((int)std::round(c.distanceLy)) + " ly)";

          addAction(s,
                    "Go (auto)",
                    "Arm auto-run to the suggested supplier station.",
                    [cb = ctx.goToStation, sysId = c.systemId, stId = c.stationId]() {
                      if (cb) cb(sysId, stId, true);
                    },
                    !ctx.goToStation);

          addAction(s,
                    "Plot",
                    "Plot an FSD route to the supplier system.",
                    [cb = ctx.routeToStation, sysId = c.systemId, stId = c.stationId]() {
                      if (cb) cb(sysId, stId);
                    },
                    !ctx.routeToStation);
        } else {
          s.body = "No nearby supplier found. Try a Trade Hub / Industrial station or widen the search.";
        }

        addAction(s,
                  "Trade Planner",
                  "Open the Trade Planner (commodity scouting).",
                  ctx.openTradePlanner);

        out.push_back(std::move(s));
      }
    }
  }
}

static void drawSuggestionTable(CopilotWindowState& st,
                                const CopilotContext& ctx,
                                const std::vector<Suggestion>& suggestions) {
  const double now = ctx.timeRealSec;

  // Prune expired dismissals.
  for (auto it = st.dismissedUntilRealSec.begin(); it != st.dismissedUntilRealSec.end(); ) {
    if (now >= it->second) it = st.dismissedUntilRealSec.erase(it);
    else ++it;
  }

  std::string_view filter{};
  if (st.filter[0] != '\0') {
    filter = std::string_view(st.filter, std::strlen(st.filter));
  }

  if (ImGui::BeginTable("##copilot_recs", 3,
                       ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_BordersInnerH
                       | ImGuiTableFlags_SizingStretchProp)) {
    ImGui::TableSetupColumn("!", ImGuiTableColumnFlags_WidthFixed, 26.0f);
    ImGui::TableSetupColumn("Recommendation", ImGuiTableColumnFlags_WidthStretch);
    ImGui::TableSetupColumn("Actions", ImGuiTableColumnFlags_WidthFixed, 220.0f);
    ImGui::TableHeadersRow();

    for (const auto& s : suggestions) {
      const auto it = st.dismissedUntilRealSec.find(s.id);
      const bool dismissed = (it != st.dismissedUntilRealSec.end() && now < it->second);
      if (dismissed && !st.showDismissed) continue;

      if (!filter.empty()) {
        if (!containsCaseInsensitive(s.title, filter) && !containsCaseInsensitive(s.body, filter)) {
          continue;
        }
      }

      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      ImGui::TextUnformatted(severityGlyph(s.severity));

      ImGui::TableSetColumnIndex(1);
      ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + ImGui::GetColumnWidth());
      ImGui::TextUnformatted(s.title.c_str());
      if (!s.body.empty()) {
        ImGui::TextDisabled("%s", s.body.c_str());
      }
      ImGui::PopTextWrapPos();

      ImGui::TableSetColumnIndex(2);
      ImGui::PushID(s.id.c_str());

      // Actions (first 3 inline, remainder via popup).
      const int inlineMax = 3;
      int shown = 0;

      for (int i = 0; i < (int)s.actions.size(); ++i) {
        if (shown >= inlineMax) break;
        const auto& a = s.actions[(std::size_t)i];

        ImGui::BeginDisabled(a.disabled);
        ImGui::PushID(i);
        if (ImGui::SmallButton(a.label.c_str())) {
          if (!a.disabled && a.fn) a.fn();
        }
        ImGui::PopID();
        ImGui::EndDisabled();

        if (!a.tooltip.empty() && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
          ImGui::SetTooltip("%s", a.tooltip.c_str());
        }

        ImGui::SameLine();
        shown++;
      }

      const bool hasMore = ((int)s.actions.size() > inlineMax);
      if (hasMore) {
        if (ImGui::SmallButton("...")) {
          ImGui::OpenPopup("more_actions");
        }
        if (ImGui::BeginPopup("more_actions")) {
          for (std::size_t i = (std::size_t)inlineMax; i < s.actions.size(); ++i) {
            const auto& a = s.actions[i];
            ImGui::BeginDisabled(a.disabled);
            if (ImGui::MenuItem(a.label.c_str())) {
              if (!a.disabled && a.fn) a.fn();
            }
            ImGui::EndDisabled();
            if (!a.tooltip.empty() && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
              ImGui::SetTooltip("%s", a.tooltip.c_str());
            }
          }
          ImGui::EndPopup();
        }
        ImGui::SameLine();
      }

      // Dismiss button.
      if (ImGui::SmallButton("Dismiss")) {
        st.dismissedUntilRealSec[s.id] = now + (double)std::max(0.0f, st.dismissTtlSec);
      }
      if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("Hide this recommendation for %.0f seconds", st.dismissTtlSec);
      }

      ImGui::PopID();
    }

    ImGui::EndTable();
  }
}

static void drawShipStatusBars(const CopilotContext& ctx, bool compact) {
  auto bar = [&](const char* label, double value, double maxValue) {
    const double denom = std::max(1e-6, maxValue);
    const float frac = (float)std::clamp(value / denom, 0.0, 1.0);
    std::string overlay = std::to_string((int)std::round(value)) + "/" + std::to_string((int)std::round(maxValue));
    ImGui::ProgressBar(frac, ImVec2(-FLT_MIN, 0.0f), overlay.c_str());
    if (!compact) {
      ImGui::SameLine();
      ImGui::TextUnformatted(label);
    }
  };

  if (ctx.hullMax > 0.0) {
    bar("Hull", ctx.hull, ctx.hullMax);
  }
  if (ctx.shieldMax > 0.0) {
    bar("Shield", ctx.shield, ctx.shieldMax);
  }
  if (ctx.fuelMax > 0.0) {
    bar("Fuel", ctx.fuel, ctx.fuelMax);
  }

  if (std::isfinite(ctx.heat)) {
    const float frac = (float)std::clamp(ctx.heat / 100.0, 0.0, 1.0);
    std::string overlay = std::to_string((int)std::round(ctx.heat)) + "%";
    ImGui::ProgressBar(frac, ImVec2(-FLT_MIN, 0.0f), overlay.c_str());
    if (!compact) {
      ImGui::SameLine();
      ImGui::TextUnformatted("Heat");
    }
  }
}

} // namespace

void drawCopilotWindow(CopilotWindowState& st, const CopilotContext& ctx) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(520.0f, 620.0f), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Copilot", &st.open)) {
    ImGui::End();
    return;
  }

  ImGui::TextDisabled("Context-aware recommendations and quick actions.");

  // Header controls.
  ImGui::Checkbox("Ship status", &st.showShipStatus);
  ImGui::SameLine();
  ImGui::Checkbox("Recommendations", &st.showRecommendations);
  ImGui::SameLine();
  ImGui::Checkbox("Missions", &st.showMissions);
  ImGui::SameLine();
  ImGui::Checkbox("Progression", &st.showProgression);

  ImGui::Checkbox("Compact", &st.compact);
  ImGui::SameLine();
  ImGui::Checkbox("Show dismissed", &st.showDismissed);

  ImGui::SetNextItemWidth(200.0f);
  ImGui::InputTextWithHint("##copilot_filter", "Filter (mission, fuel, dock...)", st.filter, sizeof(st.filter));

  ImGui::SameLine();
  if (ImGui::SmallButton("Clear")) {
    st.filter[0] = '\0';
  }

  ImGui::SameLine();
  ImGui::SetNextItemWidth(140.0f);
  ImGui::SliderFloat("Dismiss (sec)", &st.dismissTtlSec, 10.0f, 900.0f, "%.0f");

  ImGui::SameLine();
  if (ImGui::SmallButton("Reset dismissals")) {
    st.dismissedUntilRealSec.clear();
  }

  // System summary.
  if (ctx.currentSystem) {
    ImGui::Separator();
    ImGui::Text("System: %s", ctx.currentSystem->stub.name.c_str());
    if (ctx.docked && ctx.dockedStationId != 0) {
      ImGui::SameLine();
      ImGui::TextDisabled("(Docked)");
    } else if (ctx.supercruiseActive) {
      ImGui::SameLine();
      ImGui::TextDisabled("(Supercruise)");
    } else if (ctx.fsdBusy) {
      ImGui::SameLine();
      ImGui::TextDisabled("(FSD)");
    }
  }

  if (st.showShipStatus) {
    ImGui::SeparatorText("Ship");
    drawShipStatusBars(ctx, st.compact);

    if (ctx.docked && ctx.openStationServices) {
      if (ImGui::SmallButton("Open Services")) {
        ctx.openStationServices();
      }
      if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("Repairs, refuel, legal services, shipyard, etc.");
      }
    }
  }

  if (st.showProgression && !ctx.progression.title.empty()) {
    ImGui::SeparatorText("Progression");
    ImGui::TextUnformatted(ctx.progression.title.c_str());
    if (!ctx.progression.detail.empty()) {
      ImGui::TextDisabled("%s", ctx.progression.detail.c_str());
    }

    const float p = std::clamp(ctx.progression.progress01, 0.0f, 1.0f);
    ImGui::ProgressBar(p, ImVec2(-FLT_MIN, 0.0f), "");

    if (ctx.openProgressionWindow) {
      if (ImGui::SmallButton("Open Progression")) {
        ctx.openProgressionWindow();
      }
    }
  }

  // Build recommendations.
  std::vector<Suggestion> recs;
  recs.reserve(12);
  addShipRecommendations(recs, ctx);
  addMissionRecommendations(recs, ctx);

  // Sort critical first.
  std::stable_sort(recs.begin(), recs.end(), [](const Suggestion& a, const Suggestion& b) {
    return (int)a.severity > (int)b.severity;
  });

  if (st.showRecommendations) {
    ImGui::SeparatorText("Recommendations");

    if (recs.empty()) {
      ImGui::TextDisabled("No immediate recommendations.");
    } else {
      drawSuggestionTable(st, ctx, recs);
    }
  }

  if (st.showMissions && ctx.missions) {
    ImGui::SeparatorText("Active Missions");

    // Build a sorted list of active missions.
    std::vector<const sim::Mission*> active;
    active.reserve(ctx.missions->size());
    for (const auto& m : *ctx.missions) {
      if (!m.completed && !m.failed) active.push_back(&m);
    }

    if (active.empty()) {
      ImGui::TextDisabled("No active missions.");
    } else {
      std::stable_sort(active.begin(), active.end(), [&](const sim::Mission* a, const sim::Mission* b) {
        const double da = (a->deadlineDay > 0.0) ? (a->deadlineDay - ctx.timeDays) : 1e9;
        const double db = (b->deadlineDay > 0.0) ? (b->deadlineDay - ctx.timeDays) : 1e9;
        return da < db;
      });

      if (ImGui::BeginTable("##copilot_missions", 4,
                           ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_BordersInnerH
                           | ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("Type", ImGuiTableColumnFlags_WidthFixed, 110.0f);
        ImGui::TableSetupColumn("Destination", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Due", ImGuiTableColumnFlags_WidthFixed, 70.0f);
        ImGui::TableSetupColumn("Track", ImGuiTableColumnFlags_WidthFixed, 72.0f);
        ImGui::TableHeadersRow();

        int row = 0;
        for (const auto* m : active) {
          ImGui::TableNextRow();
          ImGui::TableSetColumnIndex(0);
          ImGui::TextUnformatted(sim::missionTypeName(m->type));

          ImGui::TableSetColumnIndex(1);
          const auto [nextSys, nextSt] = missionNextStop(*m);
          std::string dest = systemName(ctx.universe, nextSys);
          if (nextSt != 0) dest += " / " + stationName(ctx.universe, nextSys, nextSt);
          else if (nextSys != 0) dest += " (site)";
          ImGui::TextUnformatted(dest.c_str());

          ImGui::TableSetColumnIndex(2);
          ImGui::TextUnformatted(formatDeadline(*m, ctx.timeDays).c_str());

          ImGui::TableSetColumnIndex(3);
          ImGui::PushID(row++);
          const bool isTracked = (ctx.trackedMissionId != 0 && ctx.trackedMissionId == m->id);
          if (isTracked) {
            ImGui::TextDisabled("Tracked");
          } else {
            ImGui::BeginDisabled(!ctx.trackMission);
            if (ImGui::SmallButton("Track")) {
              if (ctx.trackMission) ctx.trackMission(m->id);
            }
            ImGui::EndDisabled();
          }
          ImGui::PopID();
        }

        ImGui::EndTable();
      }

      if (ctx.openMissionsWindow) {
        if (ImGui::SmallButton("Open Missions")) {
          ctx.openMissionsWindow();
        }
      }
    }
  }

  ImGui::End();
}

} // namespace stellar::game
