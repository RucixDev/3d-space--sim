#include "TrafficLanesWindow.h"

#include "stellar/econ/Commodity.h"
#include "stellar/math/Math.h"
#include "stellar/sim/TrafficConvoyLayer.h"
#include "stellar/sim/Universe.h"
#include "stellar/sim/Units.h"

#include "imgui.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace stellar::game {
namespace {

using stellar::core::u64;
using stellar::math::Vec3d;
using stellar::sim::Station;
using stellar::sim::StationId;
using stellar::sim::StarSystem;
using stellar::sim::TrafficConvoy;
using stellar::sim::TrafficConvoyView;
using stellar::sim::TrafficLaneParams;

static const Station* findStationById(const StarSystem& sys, StationId id) {
  for (const auto& st : sys.stations) {
    if (st.id == id) return &st;
  }
  return nullptr;
}

static std::string fmtSeconds(double sec) {
  if (!std::isfinite(sec)) return "--";
  sec = std::max(0.0, sec);

  const int s = (int)std::llround(sec);
  const int hh = s / 3600;
  const int mm = (s % 3600) / 60;
  const int ss = s % 60;

  char buf[64];
  if (hh > 0) std::snprintf(buf, sizeof(buf), "%dh %02dm %02ds", hh, mm, ss);
  else if (mm > 0) std::snprintf(buf, sizeof(buf), "%dm %02ds", mm, ss);
  else std::snprintf(buf, sizeof(buf), "%ds", ss);
  return buf;
}

struct RelMetrics {
  double distKm{0.0};
  double closingKmS{0.0};
  double tcaSec{0.0};
  double missKm{0.0};
};

// Relative-motion helper for "how hard is it to intercept this convoy?"
static RelMetrics relMetrics(const Vec3d& playerPosKm,
                             const Vec3d& playerVelKmS,
                             const Vec3d& targetPosKm,
                             const Vec3d& targetVelKmS) {
  RelMetrics out{};

  const Vec3d r = targetPosKm - playerPosKm;
  const Vec3d v = targetVelKmS - playerVelKmS;

  out.distKm = r.length();

  const double rLen = out.distKm;
  if (rLen > 1e-9) {
    out.closingKmS = -math::dot(r / rLen, v);
  } else {
    out.closingKmS = 0.0;
  }

  const double v2 = v.lengthSq();
  if (v2 > 1e-12) {
    // t* = argmin |r + v t|  =>  t* = -dot(r,v)/|v|^2
    const double t = -math::dot(r, v) / v2;
    out.tcaSec = std::max(0.0, t);
    const Vec3d rAt = r + v * out.tcaSec;
    out.missKm = rAt.length();
  } else {
    out.tcaSec = 0.0;
    out.missKm = out.distKm;
  }

  return out;
}

static ImVec2 toScreen2D(const Vec3d& pKm,
                         const Vec3d& centerKm,
                         double scalePxPerKm,
                         const ImVec2& originPx,
                         const ImVec2& sizePx) {
  const double x = (pKm.x - centerKm.x) * scalePxPerKm;
  const double y = (pKm.y - centerKm.y) * scalePxPerKm;
  // Put +Y up.
  return ImVec2(originPx.x + sizePx.x * 0.5f + (float)x,
                originPx.y + sizePx.y * 0.5f - (float)y);
}

static void expandBounds(Vec3d& mn, Vec3d& mx, const Vec3d& p) {
  mn.x = std::min(mn.x, p.x); mn.y = std::min(mn.y, p.y); mn.z = std::min(mn.z, p.z);
  mx.x = std::max(mx.x, p.x); mx.y = std::max(mx.y, p.y); mx.z = std::max(mx.z, p.z);
}

static void drawLaneMap(const TrafficConvoyView& v,
                        const StarSystem& sys,
                        const TrafficLaneParams& laneParams,
                        int segments,
                        const Vec3d& playerPosKm,
                        const ImVec2& sizePx) {
  ImDrawList* dl = ImGui::GetWindowDrawList();
  const ImVec2 origin = ImGui::GetCursorScreenPos();

  const ImU32 colBg = ImGui::GetColorU32(ImGuiCol_FrameBg);
  const ImU32 colGrid = ImGui::GetColorU32(ImVec4(1, 1, 1, 0.06f));
  const ImU32 colPath = ImGui::GetColorU32(ImVec4(0.25f, 0.70f, 1.0f, 0.65f));
  const ImU32 colStation = ImGui::GetColorU32(ImVec4(1, 1, 1, 0.85f));
  const ImU32 colStation2 = ImGui::GetColorU32(ImVec4(1, 1, 1, 0.25f));
  const ImU32 colConvoy = ImGui::GetColorU32(ImVec4(1.0f, 0.70f, 0.25f, 0.95f));
  const ImU32 colPlayer = ImGui::GetColorU32(ImVec4(0.60f, 1.0f, 0.60f, 0.95f));

  dl->AddRectFilled(origin, ImVec2(origin.x + sizePx.x, origin.y + sizePx.y), colBg, 6.0f);
  dl->AddRect(origin, ImVec2(origin.x + sizePx.x, origin.y + sizePx.y),
              ImGui::GetColorU32(ImGuiCol_Border), 6.0f);

  // Sample the path polyline.
  std::vector<Vec3d> pts = stellar::sim::sampleTrafficConvoyPathKm(v.convoy, sys, std::max(4, segments), laneParams);

  // Bounds: stations + path + convoy + player (XY projection).
  Vec3d mn{+1e30, +1e30, +1e30};
  Vec3d mx{-1e30, -1e30, -1e30};

  for (const auto& st : sys.stations) {
    const Vec3d p = stellar::sim::stationPosKm(st, v.convoy.departDay);
    expandBounds(mn, mx, p);
  }
  for (const auto& p : pts) expandBounds(mn, mx, p);
  expandBounds(mn, mx, v.state.posKm);
  expandBounds(mn, mx, playerPosKm);

  const Vec3d center{(mn.x + mx.x) * 0.5, (mn.y + mx.y) * 0.5, 0.0};
  const double spanX = std::max(1.0, mx.x - mn.x);
  const double spanY = std::max(1.0, mx.y - mn.y);
  const double span = std::max(spanX, spanY);

  const double pad = 0.18;
  const double scale = (std::min((double)sizePx.x, (double)sizePx.y) * (0.5 - pad)) / span;

  // Grid (simple crosshair + box thirds).
  {
    const ImVec2 c = ImVec2(origin.x + sizePx.x * 0.5f, origin.y + sizePx.y * 0.5f);
    dl->AddLine(ImVec2(origin.x, c.y), ImVec2(origin.x + sizePx.x, c.y), colGrid, 1.0f);
    dl->AddLine(ImVec2(c.x, origin.y), ImVec2(c.x, origin.y + sizePx.y), colGrid, 1.0f);

    for (int i = 1; i <= 2; ++i) {
      const float fx = origin.x + sizePx.x * (float)i / 3.0f;
      const float fy = origin.y + sizePx.y * (float)i / 3.0f;
      dl->AddLine(ImVec2(fx, origin.y), ImVec2(fx, origin.y + sizePx.y), colGrid, 1.0f);
      dl->AddLine(ImVec2(origin.x, fy), ImVec2(origin.x + sizePx.x, fy), colGrid, 1.0f);
    }
  }

  // Path
  for (std::size_t i = 1; i < pts.size(); ++i) {
    const ImVec2 a = toScreen2D(pts[i - 1], center, scale, origin, sizePx);
    const ImVec2 b = toScreen2D(pts[i], center, scale, origin, sizePx);
    dl->AddLine(a, b, colPath, 2.0f);
  }

  // Stations
  for (const auto& st : sys.stations) {
    const Vec3d p = stellar::sim::stationPosKm(st, v.convoy.departDay);
    const ImVec2 sp = toScreen2D(p, center, scale, origin, sizePx);
    dl->AddCircleFilled(sp, 4.0f, colStation);
    dl->AddCircle(sp, 11.0f, colStation2, 0, 1.0f);
    dl->AddText(ImVec2(sp.x + 6.0f, sp.y - 12.0f), colStation, st.name.c_str());
  }

  // Player
  {
    const ImVec2 sp = toScreen2D(playerPosKm, center, scale, origin, sizePx);
    dl->AddCircleFilled(sp, 4.0f, colPlayer);
    dl->AddText(ImVec2(sp.x + 6.0f, sp.y - 12.0f), colPlayer, "You");
  }

  // Convoy
  {
    const ImVec2 sp = toScreen2D(v.state.posKm, center, scale, origin, sizePx);
    dl->AddCircleFilled(sp, 5.0f, colConvoy);
    dl->AddText(ImVec2(sp.x + 6.0f, sp.y - 12.0f), colConvoy, "Convoy");
  }
}

static std::vector<TrafficConvoyView> buildTrafficViews(const TrafficLanesContext& ctx, bool includeInactive) {
  std::vector<TrafficConvoyView> out;
  if (!ctx.currentSystem) return out;

  const StarSystem& sys = *ctx.currentSystem;
  const TrafficLaneParams laneParams = ctx.laneParams;

  if (ctx.trafficLedger) {
    out = stellar::sim::generateTrafficConvoysFromLedger(*ctx.trafficLedger,
                                                         sys,
                                                         ctx.timeDays,
                                                         laneParams.genWindowDays,
                                                         includeInactive,
                                                         laneParams);
  } else {
    // Deterministic fallback: convoys are scheduled from the universe seed.
    TrafficLaneParams p = laneParams;
    p.includeInactive = includeInactive;
    out = stellar::sim::generateTrafficConvoys(ctx.universe.seed(), sys, ctx.timeDays, p);
  }

  // Optional filter (e.g., convoy already interdicted in save).
  if (ctx.isConvoySuppressed) {
    out.erase(std::remove_if(out.begin(), out.end(),
                             [&](const TrafficConvoyView& v) { return ctx.isConvoySuppressed(v.convoy.id); }),
              out.end());
  }

  return out;
}

} // namespace

void drawTrafficLanesWindow(TrafficLanesWindowState& st, const TrafficLanesContext& ctx) {
  if (!st.open) return;

  if (!ImGui::Begin("Traffic Lanes", &st.open)) {
    ImGui::End();
    return;
  }

  if (!ctx.currentSystem) {
    ImGui::TextDisabled("No system loaded.");
    ImGui::End();
    return;
  }

  const StarSystem& sys = *ctx.currentSystem;

  ImGui::TextDisabled("System: %s | stations %d", sys.stub.name.c_str(), (int)sys.stations.size());
  ImGui::Separator();

  ImGui::Checkbox("Show inactive (preview)", &st.showInactive);
  ImGui::SameLine();
  ImGui::Checkbox("Show lane map", &st.showMap);

  ImGui::SameLine();
  ImGui::SetNextItemWidth(110.0f);
  ImGui::DragInt("Path segments", &st.pathSegments, 1.0f, 8, 256);

  ImGui::SameLine();
  ImGui::SetNextItemWidth(110.0f);
  ImGui::DragFloat("Map height", &st.mapHeightPx, 1.0f, 160.0f, 900.0f, "%.0f px");

  // Build the convoy view list.
  std::vector<TrafficConvoyView> views = buildTrafficViews(ctx, st.showInactive);

  // Sort by: active first, then distance-to-player.
  std::sort(views.begin(), views.end(),
            [&](const TrafficConvoyView& a, const TrafficConvoyView& b) {
              if (a.state.active != b.state.active) return a.state.active > b.state.active;
              const double da = (a.state.posKm - ctx.playerPosKm).lengthSq();
              const double db = (b.state.posKm - ctx.playerPosKm).lengthSq();
              return da < db;
            });

  ImGui::TextDisabled("Convoys in window: %d", (int)views.size());

  // Table
  const ImGuiTableFlags flags = ImGuiTableFlags_BordersV
                             | ImGuiTableFlags_BordersOuterH
                             | ImGuiTableFlags_RowBg
                             | ImGuiTableFlags_ScrollY
                             | ImGuiTableFlags_Resizable;

  const float tableHeight = st.showMap ? std::max(180.0f, st.mapHeightPx * 0.55f) : 320.0f;
  if (ImGui::BeginTable("##traffic_convoys", 9, flags, ImVec2(0, tableHeight))) {
    ImGui::TableSetupScrollFreeze(0, 1);
    ImGui::TableSetupColumn("Status", ImGuiTableColumnFlags_WidthFixed, 62.0f);
    ImGui::TableSetupColumn("Commodity");
    ImGui::TableSetupColumn("Units", ImGuiTableColumnFlags_WidthFixed, 70.0f);
    ImGui::TableSetupColumn("From");
    ImGui::TableSetupColumn("To");
    ImGui::TableSetupColumn("ETA", ImGuiTableColumnFlags_WidthFixed, 90.0f);
    ImGui::TableSetupColumn("Dist", ImGuiTableColumnFlags_WidthFixed, 80.0f);
    ImGui::TableSetupColumn("TCA / miss", ImGuiTableColumnFlags_WidthFixed, 120.0f);
    ImGui::TableSetupColumn("Action", ImGuiTableColumnFlags_WidthFixed, 70.0f);
    ImGui::TableHeadersRow();

    for (const auto& v : views) {
      const bool selected = (st.selectedConvoyId != 0 && st.selectedConvoyId == v.convoy.id);

      const auto def = stellar::econ::commodityDef(v.convoy.commodity);
      const Station* fromSt = findStationById(sys, v.convoy.fromStation);
      const Station* toSt = findStationById(sys, v.convoy.toStation);

      const RelMetrics rm = relMetrics(ctx.playerPosKm, ctx.playerVelKmS, v.state.posKm, v.state.velKmS);

      const double etaSec = std::max(0.0, (v.convoy.arriveDay - ctx.timeDays) * 86400.0);

      ImGui::TableNextRow();
      ImGui::TableNextColumn();

      // Row selection click target (first column).
      {
        char idbuf[64];
        std::snprintf(idbuf, sizeof(idbuf), "%s##sel_%llx",
                      v.state.active ? "Active" : "Idle",
                      (unsigned long long)v.convoy.id);
        if (ImGui::Selectable(idbuf, selected, ImGuiSelectableFlags_SpanAllColumns)) {
          st.selectedConvoyId = v.convoy.id;
        }
      }

      ImGui::TableNextColumn();
      ImGui::TextUnformatted(def.name);

      ImGui::TableNextColumn();
      ImGui::Text("%.0f", std::round(std::max(0.0, v.convoy.units)));

      ImGui::TableNextColumn();
      ImGui::TextUnformatted(fromSt ? fromSt->name.c_str() : "--");

      ImGui::TableNextColumn();
      ImGui::TextUnformatted(toSt ? toSt->name.c_str() : "--");

      ImGui::TableNextColumn();
      ImGui::TextUnformatted(fmtSeconds(etaSec).c_str());

      ImGui::TableNextColumn();
      ImGui::Text("%.0f km", std::round(rm.distKm));

      ImGui::TableNextColumn();
      {
        if (rm.tcaSec <= 0.0) {
          ImGui::TextDisabled("--");
        } else {
          char buf[96];
          std::snprintf(buf, sizeof(buf), "%s / %.0f",
                        fmtSeconds(rm.tcaSec).c_str(),
                        std::round(rm.missKm));
          ImGui::TextUnformatted(buf);
        }
      }

      ImGui::TableNextColumn();
      ImGui::BeginDisabled(!ctx.targetTrafficConvoy);
      {
        char b[64];
        std::snprintf(b, sizeof(b), "Target##t_%llx", (unsigned long long)v.convoy.id);
        if (ImGui::SmallButton(b)) {
          ctx.targetTrafficConvoy(v.convoy.id);
          st.selectedConvoyId = v.convoy.id;
          if (ctx.toast) ctx.toast("Targeting convoy signal.", 1.8);
        }
      }
      ImGui::EndDisabled();
    }

    ImGui::EndTable();
  }

  // Selected detail
  const TrafficConvoyView* sel = nullptr;
  for (const auto& v : views) {
    if (v.convoy.id == st.selectedConvoyId) { sel = &v; break; }
  }

  ImGui::Separator();
  if (!sel) {
    ImGui::TextDisabled("Select a convoy row to see details.");
    ImGui::End();
    return;
  }

  const auto def = stellar::econ::commodityDef(sel->convoy.commodity);
  const Station* fromSt = findStationById(sys, sel->convoy.fromStation);
  const Station* toSt = findStationById(sys, sel->convoy.toStation);

  ImGui::Text("Convoy: %s", def.name);
  ImGui::TextDisabled("ID: %llx", (unsigned long long)sel->convoy.id);
  ImGui::TextDisabled("Route: %s -> %s",
                      fromSt ? fromSt->name.c_str() : "--",
                      toSt ? toSt->name.c_str() : "--");

  const double departSec = (sel->convoy.departDay - ctx.timeDays) * 86400.0;
  const double arriveSec = (sel->convoy.arriveDay - ctx.timeDays) * 86400.0;
  ImGui::TextDisabled("Schedule: dep %s | arr %s",
                      fmtSeconds(departSec).c_str(),
                      fmtSeconds(arriveSec).c_str());

  ImGui::TextDisabled("State: %s | progress %.0f%% | speed %.0f km/s",
                      sel->state.active ? "active" : "inactive",
                      std::round(sel->state.progress01 * 100.0),
                      std::round(sel->state.speedKmS));

  const RelMetrics rm = relMetrics(ctx.playerPosKm, ctx.playerVelKmS, sel->state.posKm, sel->state.velKmS);
  ImGui::TextDisabled("Relative: dist %.0f km | closing %.0f km/s | miss %.0f km @ %s",
                      std::round(rm.distKm),
                      std::round(rm.closingKmS),
                      std::round(rm.missKm),
                      fmtSeconds(rm.tcaSec).c_str());

  // Actions
  ImGui::BeginDisabled(!ctx.targetTrafficConvoy);
  if (ImGui::Button("Target this convoy")) {
    ctx.targetTrafficConvoy(sel->convoy.id);
    if (ctx.toast) ctx.toast("Target set: convoy signal.", 1.8);
  }
  ImGui::EndDisabled();

  ImGui::SameLine();
  ImGui::BeginDisabled(!ctx.requestSupercruise);
  if (ImGui::Button("Supercruise to target")) {
    if (ctx.targetTrafficConvoy) ctx.targetTrafficConvoy(sel->convoy.id);
    ctx.requestSupercruise();
  }
  ImGui::EndDisabled();

  if (st.showMap) {
    ImGui::Separator();
    ImGui::TextDisabled("Lane map (XY projection)");
    const ImVec2 size = ImVec2(ImGui::GetContentRegionAvail().x,
                               std::max(120.0f, st.mapHeightPx));
    drawLaneMap(*sel, sys, ctx.laneParams, st.pathSegments, ctx.playerPosKm, size);
    ImGui::InvisibleButton("##tlmap_canvas", size);
  }

  ImGui::End();
}

} // namespace stellar::game
