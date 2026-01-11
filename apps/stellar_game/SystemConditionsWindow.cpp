#include "SystemConditionsWindow.h"

#include "stellar/sim/SystemConditions.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <utility>
#include <string>

#include <imgui.h>

namespace stellar::game {
namespace {

static std::string_view factionName(const stellar::sim::Universe& u, stellar::core::u32 factionId) {
  if (factionId == 0) return "Independent";
  const auto& f = u.factions();
  for (const auto& fac : f) {
    if (fac.id == factionId) return fac.name;
  }
  return "Unknown";
}

static double systemDistanceLy(const stellar::sim::StarSystem& a, const stellar::sim::StarSystem& b) {
  return (a.stub.posLy - b.stub.posLy).length();
}

static const stellar::sim::SystemSecurityDeltaState* findDeltaState(
    const std::unordered_map<stellar::sim::SystemId, stellar::sim::SystemSecurityDeltaState>* map,
    stellar::sim::SystemId id) {
  if (!map) return nullptr;
  const auto it = map->find(id);
  if (it == map->end()) return nullptr;
  return &it->second;
}

static const char* sortKeyName(SystemConditionsWindowState::SortKey k) {
  switch (k) {
    case SystemConditionsWindowState::SortKey::Distance: return "Distance";
    case SystemConditionsWindowState::SortKey::Security: return "Security";
    case SystemConditionsWindowState::SortKey::Piracy: return "Piracy";
    case SystemConditionsWindowState::SortKey::Traffic: return "Traffic";
    case SystemConditionsWindowState::SortKey::EventSeverity: return "Event severity";
  }
  return "Distance";
}

} // namespace

void drawSystemConditionsWindow(SystemConditionsWindowState& st, const SystemConditionsWindowContext& ctx) {
  if (!st.open) return;

  if (!ImGui::Begin("System Conditions", &st.open)) {
    ImGui::End();
    return;
  }

  if (!ctx.currentSystem) {
    ImGui::TextDisabled("No current system.");
    ImGui::End();
    return;
  }

  const int dayStamp = (int)std::floor(ctx.timeDays);

  // --- Current system header ---
  const auto* deltaNow = findDeltaState(ctx.systemSecurityDeltaBySystem, ctx.currentSystem->stub.id);
  const stellar::sim::SystemConditionsSnapshot snap = stellar::sim::snapshotSystemConditions(
      ctx.universe.seed(), *ctx.currentSystem, ctx.timeDays, deltaNow, ctx.dynamicsParams, ctx.eventParams);

  ImGui::Text("System: %s (id=%llu)", ctx.currentSystem->stub.name.c_str(),
              (unsigned long long)ctx.currentSystem->stub.id);
  ImGui::Text("Controlling faction: %s",
              std::string(factionName(ctx.universe, snap.base.controllingFactionId)).c_str());

  ImGui::Separator();

  // Event summary.
  if (snap.event.active) {
    const double endsIn = snap.event.endDay - ctx.timeDays;
    ImGui::Text("Event: %s (severity %.2f, ends in %.1f days)",
                stellar::sim::systemEventKindName(snap.event.kind),
                snap.event.severity01,
                endsIn);
    ImGui::TextDisabled("Deltas: security %+0.2f  piracy %+0.2f  traffic %+0.2f",
                        snap.event.securityDelta, snap.event.piracyDelta, snap.event.trafficDelta);
  } else {
    ImGui::Text("Event: none");
  }

  // Profile summary.
  ImGui::Text("Base:      security %.2f  piracy %.2f  traffic %.2f",
              snap.base.security01, snap.base.piracy01, snap.base.traffic01);
  if (snap.hasDynamics) {
    ImGui::Text("Dynamics:  security %+0.2f  piracy %+0.2f  traffic %+0.2f",
                snap.dynamicsNow.securityDelta, snap.dynamicsNow.piracyDelta, snap.dynamicsNow.trafficDelta);
  } else {
    ImGui::TextDisabled("Dynamics: none");
  }
  ImGui::Text("Effective: security %.2f  piracy %.2f  traffic %.2f",
              snap.effective.security01, snap.effective.piracy01, snap.effective.traffic01);

  ImGui::Separator();

  // --- Controls ---
  ImGui::TextDisabled("Nearby systems");
  ImGui::SetNextItemWidth(220.0f);
  ImGui::SliderFloat("Radius (ly)", &st.radiusLy, 20.0f, 1600.0f, "%.0f");
  ImGui::SetNextItemWidth(220.0f);
  ImGui::SliderInt("Max results", &st.maxResults, 16, 256);

  ImGui::Checkbox("Only active events", &st.onlyActiveEvents);
  ImGui::SameLine();
  ImGui::Checkbox("Include current system", &st.includeCurrentSystem);

  ImGui::SetNextItemWidth(220.0f);
  int sortIdx = (int)st.sortKey;
  const char* sortItems[] = {"Distance", "Security", "Piracy", "Traffic", "Event severity"};
  if (ImGui::Combo("Sort", &sortIdx, sortItems, (int)(sizeof(sortItems) / sizeof(sortItems[0])))) {
    st.sortKey = (SystemConditionsWindowState::SortKey)sortIdx;
  }
  ImGui::SameLine();
  ImGui::Checkbox("Desc", &st.sortDescending);

  // Cache invalidation.
  const bool cacheDirty =
      (st.cacheFromSystem != ctx.currentSystem->stub.id) ||
      (st.cacheDayStamp != dayStamp) ||
      (std::abs(st.cacheRadiusLy - st.radiusLy) > 1e-6f) ||
      (st.cacheMaxResults != st.maxResults) ||
      (st.cacheOnlyActiveEvents != st.onlyActiveEvents) ||
      (st.cacheIncludeCurrentSystem != st.includeCurrentSystem) ||
      (st.cacheSortKey != st.sortKey) ||
      (st.cacheSortDescending != st.sortDescending);

  if (cacheDirty) {
    const auto t0 = std::chrono::high_resolution_clock::now();

    st.cacheRows.clear();
    st.cacheRows.reserve((std::size_t)std::max(0, st.maxResults));

    const double radius = std::max(0.0, (double)st.radiusLy);
    const std::size_t maxRes = (std::size_t)std::max(1, st.maxResults);

    auto nearby = ctx.universe.queryNearby(ctx.currentSystem->stub.posLy, radius, maxRes);

    for (const auto& stub : nearby) {
      if (!st.includeCurrentSystem && stub.id == ctx.currentSystem->stub.id) continue;

      const auto& sys = ctx.universe.getSystem(stub.id);
      const auto* d = findDeltaState(ctx.systemSecurityDeltaBySystem, stub.id);
      const auto snapRow = stellar::sim::snapshotSystemConditions(
          ctx.universe.seed(), sys, ctx.timeDays, d, ctx.dynamicsParams, ctx.eventParams);

      if (st.onlyActiveEvents && !snapRow.event.active) continue;

      SystemConditionsWindowState::Row row{};
      row.id = stub.id;
      row.distanceLy = systemDistanceLy(*ctx.currentSystem, sys);
      row.controllingFactionId = snapRow.base.controllingFactionId;
      row.event = snapRow.event;
      row.base = snapRow.base;
      row.effective = snapRow.effective;
      st.cacheRows.push_back(std::move(row));
    }

    auto keyValue = [&](const SystemConditionsWindowState::Row& r) -> double {
      switch (st.sortKey) {
        case SystemConditionsWindowState::SortKey::Distance: return r.distanceLy;
        case SystemConditionsWindowState::SortKey::Security: return r.effective.security01;
        case SystemConditionsWindowState::SortKey::Piracy: return r.effective.piracy01;
        case SystemConditionsWindowState::SortKey::Traffic: return r.effective.traffic01;
        case SystemConditionsWindowState::SortKey::EventSeverity: return r.event.active ? r.event.severity01 : 0.0;
      }
      return r.distanceLy;
    };

    std::sort(st.cacheRows.begin(), st.cacheRows.end(),
              [&](const auto& a, const auto& b) {
                const double ka = keyValue(a);
                const double kb = keyValue(b);
                if (std::abs(ka - kb) > 1e-12) {
                  return st.sortDescending ? (ka > kb) : (ka < kb);
                }
                if (std::abs(a.distanceLy - b.distanceLy) > 1e-12) {
                  return a.distanceLy < b.distanceLy;
                }
                return a.id < b.id;
              });

    const auto t1 = std::chrono::high_resolution_clock::now();
    st.lastComputeMs = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // Update cache key.
    st.cacheFromSystem = ctx.currentSystem->stub.id;
    st.cacheDayStamp = dayStamp;
    st.cacheRadiusLy = st.radiusLy;
    st.cacheMaxResults = st.maxResults;
    st.cacheOnlyActiveEvents = st.onlyActiveEvents;
    st.cacheIncludeCurrentSystem = st.includeCurrentSystem;
    st.cacheSortKey = st.sortKey;
    st.cacheSortDescending = st.sortDescending;
  }

  ImGui::TextDisabled("Computed in %.2f ms  |  Sort: %s%s",
                      st.lastComputeMs,
                      sortKeyName(st.sortKey),
                      st.sortDescending ? " (desc)" : "");

  // --- Table ---
  const ImGuiTableFlags flags =
      ImGuiTableFlags_BordersInnerV |
      ImGuiTableFlags_BordersOuter |
      ImGuiTableFlags_RowBg |
      ImGuiTableFlags_ScrollY |
      ImGuiTableFlags_Resizable;

  if (ImGui::BeginTable("##sys_conditions", 10, flags, ImVec2(0.0f, 360.0f))) {
    ImGui::TableSetupColumn("System", ImGuiTableColumnFlags_WidthStretch);
    ImGui::TableSetupColumn("Dist", ImGuiTableColumnFlags_WidthFixed, 55.0f);
    ImGui::TableSetupColumn("Faction", ImGuiTableColumnFlags_WidthStretch);
    ImGui::TableSetupColumn("Event", ImGuiTableColumnFlags_WidthStretch);
    ImGui::TableSetupColumn("Ends", ImGuiTableColumnFlags_WidthFixed, 55.0f);
    ImGui::TableSetupColumn("Sec", ImGuiTableColumnFlags_WidthFixed, 45.0f);
    ImGui::TableSetupColumn("Pir", ImGuiTableColumnFlags_WidthFixed, 45.0f);
    ImGui::TableSetupColumn("Trf", ImGuiTableColumnFlags_WidthFixed, 45.0f);
    ImGui::TableSetupColumn("Plot", ImGuiTableColumnFlags_WidthFixed, 44.0f);
    ImGui::TableSetupColumn("Target", ImGuiTableColumnFlags_WidthFixed, 55.0f);
    ImGui::TableHeadersRow();

    for (const auto& r : st.cacheRows) {
      const auto& sys = ctx.universe.getSystem(r.id);

      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      ImGui::TextUnformatted(sys.stub.name.c_str());

      ImGui::TableSetColumnIndex(1);
      ImGui::Text("%.0f", r.distanceLy);

      ImGui::TableSetColumnIndex(2);
      ImGui::TextUnformatted(std::string(factionName(ctx.universe, r.controllingFactionId)).c_str());

      ImGui::TableSetColumnIndex(3);
      if (r.event.active) {
        ImGui::Text("%s (%.2f)", stellar::sim::systemEventKindName(r.event.kind), r.event.severity01);
        if (ImGui::IsItemHovered()) {
          ImGui::BeginTooltip();
          ImGui::Text("%s", stellar::sim::systemEventKindName(r.event.kind));
          ImGui::Separator();
          ImGui::Text("Severity: %.2f", r.event.severity01);
          ImGui::Text("Deltas: sec %+0.2f  pir %+0.2f  trf %+0.2f",
                      r.event.securityDelta, r.event.piracyDelta, r.event.trafficDelta);
          ImGui::EndTooltip();
        }
      } else {
        ImGui::TextDisabled("—");
      }

      ImGui::TableSetColumnIndex(4);
      if (r.event.active) {
        const double endsIn = r.event.endDay - ctx.timeDays;
        ImGui::Text("%.1f", endsIn);
      } else {
        ImGui::TextDisabled("—");
      }

      ImGui::TableSetColumnIndex(5);
      ImGui::Text("%.2f", r.effective.security01);

      ImGui::TableSetColumnIndex(6);
      ImGui::Text("%.2f", r.effective.piracy01);

      ImGui::TableSetColumnIndex(7);
      ImGui::Text("%.2f", r.effective.traffic01);

      ImGui::TableSetColumnIndex(8);
      {
        const bool enabled = (bool)ctx.plotRouteToSystem;
        if (!enabled) ImGui::BeginDisabled();
        const std::string label = "Plot##" + std::to_string((unsigned long long)r.id);
        if (ImGui::SmallButton(label.c_str())) {
          const bool ok = ctx.plotRouteToSystem(r.id);
          if (!ok && ctx.toast) ctx.toast("Couldn't plot route.", 2.0);
        }
        if (!enabled) ImGui::EndDisabled();
      }

      ImGui::TableSetColumnIndex(9);
      {
        const bool enabled = (bool)ctx.targetSystem;
        if (!enabled) ImGui::BeginDisabled();
        const std::string label = "Select##" + std::to_string((unsigned long long)r.id);
        if (ImGui::SmallButton(label.c_str())) {
          ctx.targetSystem(r.id);
        }
        if (!enabled) ImGui::EndDisabled();
      }
    }

    ImGui::EndTable();
  }

  ImGui::End();
}

} // namespace stellar::game
