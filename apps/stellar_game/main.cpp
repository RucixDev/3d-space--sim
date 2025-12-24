#include "stellar/core/Log.h"
#include "stellar/core/Random.h"
#include "stellar/econ/Market.h"
#include "stellar/econ/RoutePlanner.h"
#include "stellar/render/Camera.h"
#include "stellar/render/Gl.h"
#include "stellar/render/LineRenderer.h"
#include "stellar/render/Mesh.h"
#include "stellar/render/MeshRenderer.h"
#include "stellar/render/PointRenderer.h"
#include "stellar/render/Texture.h"
#include "stellar/sim/Orbit.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/Ship.h"
#include "stellar/sim/Universe.h"

#include <SDL.h>
#include <SDL_opengl.h>

#include <backends/imgui_impl_opengl3.h>
#include <backends/imgui_impl_sdl2.h>
#include <imgui.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

using namespace stellar;

static constexpr double kAU_KM = 149597870.7;
static constexpr double kSOLAR_RADIUS_KM = 695700.0;
static constexpr double kEARTH_RADIUS_KM = 6371.0;
static constexpr double kRENDER_UNIT_KM = 1.0e6; // 1 unit = 1 million km

static void matToFloat(const math::Mat4d& m, float out[16]) {
  for (int i = 0; i < 16; ++i) out[i] = static_cast<float>(m.m[i]);
}

static const char* starClassName(sim::StarClass c) {
  switch (c) {
    case sim::StarClass::O: return "O";
    case sim::StarClass::B: return "B";
    case sim::StarClass::A: return "A";
    case sim::StarClass::F: return "F";
    case sim::StarClass::G: return "G";
    case sim::StarClass::K: return "K";
    case sim::StarClass::M: return "M";
    default: return "?";
  }
}

static const char* stationTypeName(econ::StationType t) {
  switch (t) {
    case econ::StationType::Outpost: return "Outpost";
    case econ::StationType::Agricultural: return "Agricultural";
    case econ::StationType::Mining: return "Mining";
    case econ::StationType::Refinery: return "Refinery";
    case econ::StationType::Industrial: return "Industrial";
    case econ::StationType::Research: return "Research";
    case econ::StationType::TradeHub: return "Trade Hub";
    case econ::StationType::Shipyard: return "Shipyard";
    default: return "?";
  }
}

static double clampd(double v, double a, double b) {
  return std::max(a, std::min(v, b));
}

static math::Vec3d safeNorm(const math::Vec3d& v, const math::Vec3d& fallback = {0, 0, 1}) {
  const double len = v.length();
  if (len <= 1e-12) return fallback;
  return v * (1.0 / len);
}

static math::Vec3d orbitPosKm(const sim::OrbitElements& el, double timeDays) {
  return sim::orbitPosition3DAU(el, timeDays) * kAU_KM;
}

static math::Vec3d orbitVelKmS(const sim::OrbitElements& el, double timeDays) {
  // Finite-difference (central) derivative.
  const double dtDays = 1.0e-4; // 8.64 seconds
  const double dtSec = dtDays * 86400.0;
  const math::Vec3d p0 = orbitPosKm(el, timeDays - dtDays);
  const math::Vec3d p1 = orbitPosKm(el, timeDays + dtDays);
  return (p1 - p0) * (1.0 / (2.0 * dtSec));
}

static bool projectToScreen(const math::Vec3d& worldU,
                            const math::Mat4d& view,
                            const math::Mat4d& proj,
                            int screenW,
                            int screenH,
                            ImVec2& outPx,
                            bool& outBehind) {
  // clip = proj * view * vec4(world,1)
  const math::Mat4d vp = proj * view;

  const double x = vp.m[0] * worldU.x + vp.m[4] * worldU.y + vp.m[8] * worldU.z + vp.m[12];
  const double y = vp.m[1] * worldU.x + vp.m[5] * worldU.y + vp.m[9] * worldU.z + vp.m[13];
  const double z = vp.m[2] * worldU.x + vp.m[6] * worldU.y + vp.m[10] * worldU.z + vp.m[14];
  const double w = vp.m[3] * worldU.x + vp.m[7] * worldU.y + vp.m[11] * worldU.z + vp.m[15];

  if (std::abs(w) < 1e-9) {
    outPx = ImVec2(0, 0);
    outBehind = true;
    return false;
  }

  const double ndcX = x / w;
  const double ndcY = y / w;
  const double ndcZ = z / w;

  outBehind = (ndcZ < 0.0);

  // Convert to pixels
  const float px = (float)((ndcX * 0.5 + 0.5) * (double)screenW);
  const float py = (float)((-ndcY * 0.5 + 0.5) * (double)screenH);
  outPx = ImVec2(px, py);

  return (ndcX >= -1.2 && ndcX <= 1.2 && ndcY >= -1.2 && ndcY <= 1.2);
}

struct Toast {
  std::string text;
  double untilSec{0.0};
};

static void pushToast(std::vector<Toast>& toasts, const std::string& msg, double nowSec, double durSec = 2.5) {
  toasts.push_back({msg, nowSec + durSec});
}

static void drawToasts(std::vector<Toast>& toasts, double nowSec) {
  // Drop expired
  toasts.erase(std::remove_if(toasts.begin(), toasts.end(), [&](const Toast& t) { return nowSec > t.untilSec; }),
               toasts.end());
  if (toasts.empty()) return;

  ImGuiViewport* vp = ImGui::GetMainViewport();
  ImDrawList* dl = ImGui::GetForegroundDrawList(vp);

  const float pad = 10.0f;
  ImVec2 p = ImVec2(vp->Pos.x + pad, vp->Pos.y + pad);

  for (const auto& t : toasts) {
    const float alpha = (float)clampd((t.untilSec - nowSec) / 0.35, 0.0, 1.0);
    const ImU32 bg = IM_COL32(20, 20, 24, (int)(200 * alpha));
    const ImU32 fg = IM_COL32(220, 220, 235, (int)(255 * alpha));

    const ImVec2 sz = ImGui::CalcTextSize(t.text.c_str());
    const ImVec2 boxMin = p;
    const ImVec2 boxMax = ImVec2(p.x + sz.x + 18.0f, p.y + sz.y + 12.0f);

    dl->AddRectFilled(boxMin, boxMax, bg, 6.0f);
    dl->AddRect(boxMin, boxMax, IM_COL32(80, 80, 95, (int)(200 * alpha)), 6.0f);
    dl->AddText(ImVec2(p.x + 9.0f, p.y + 6.0f), fg, t.text.c_str());

    p.y += (sz.y + 16.0f);
  }
}

enum class FlightMode : int {
  Normal = 0,
  Supercruise = 1,
  Docked = 2,
};

enum class TargetKind : int {
  None = 0,
  Station = 1,
  Planet = 2,
};

struct NavTarget {
  TargetKind kind{TargetKind::None};
  int index{-1};
};

enum class ClearanceState : int {
  None = 0,
  Granted = 1,
  Denied = 2,
};

struct DockingClearance {
  sim::StationId stationId{0};
  ClearanceState state{ClearanceState::None};
  double untilTimeDays{0.0};
};

static std::optional<int> findStationIndexById(const sim::StarSystem& sys, sim::StationId id) {
  for (std::size_t i = 0; i < sys.stations.size(); ++i) {
    if (sys.stations[i].id == id) return (int)i;
  }
  return std::nullopt;
}

static std::optional<int> findNearestStationIndex(const sim::StarSystem& sys, double timeDays, const math::Vec3d& shipPosKm) {
  if (sys.stations.empty()) return std::nullopt;
  double bestD2 = 1e300;
  int best = 0;
  for (std::size_t i = 0; i < sys.stations.size(); ++i) {
    const auto& st = sys.stations[i];
    const math::Vec3d p = orbitPosKm(st.orbit, timeDays);
    const double d2 = (p - shipPosKm).lengthSq();
    if (d2 < bestD2) {
      bestD2 = d2;
      best = (int)i;
    }
  }
  return best;
}

static void spawnNearStation(sim::Ship& ship, const sim::StarSystem& sys, double timeDays, int stationIndex) {
  if (sys.stations.empty()) {
    ship.setPositionKm({0, 0, -8000.0});
    ship.setVelocityKmS({0, 0, 0});
    ship.setAngularVelocityRadS({0, 0, 0});
    ship.setOrientation(math::Quatd::identity());
    return;
  }

  stationIndex = std::clamp(stationIndex, 0, (int)sys.stations.size() - 1);
  const auto& st = sys.stations[(std::size_t)stationIndex];

  const math::Vec3d stPos = orbitPosKm(st.orbit, timeDays);
  const math::Vec3d stVel = orbitVelKmS(st.orbit, timeDays);
  const math::Vec3d approachDir = safeNorm(stPos);

  // Spawn outside the corridor, facing the station.
  const double spawnAxial = st.docking.corridorLengthKm + std::max(10.0, st.docking.commsRangeKm * 0.5);
  ship.setPositionKm(stPos + approachDir * spawnAxial);
  ship.setVelocityKmS(stVel);
  ship.setAngularVelocityRadS({0, 0, 0});

  // Point nose along -approachDir (towards station).
  // Build a quaternion from a "look" direction by aligning local +Z to desired.
  const math::Vec3d fwd = safeNorm(-approachDir);
  const math::Vec3d up = {0, 1, 0};
  const math::Vec3d right = safeNorm(math::cross(up, fwd), {1, 0, 0});
  const math::Vec3d up2 = math::cross(fwd, right);

  // Convert basis to quaternion (column-major rotation matrix).
  const double m00 = right.x, m01 = up2.x, m02 = fwd.x;
  const double m10 = right.y, m11 = up2.y, m12 = fwd.y;
  const double m20 = right.z, m21 = up2.z, m22 = fwd.z;

  const double tr = m00 + m11 + m22;
  math::Quatd q;
  if (tr > 0.0) {
    const double S = std::sqrt(tr + 1.0) * 2.0;
    q.w = 0.25 * S;
    q.x = (m21 - m12) / S;
    q.y = (m02 - m20) / S;
    q.z = (m10 - m01) / S;
  } else if ((m00 > m11) && (m00 > m22)) {
    const double S = std::sqrt(1.0 + m00 - m11 - m22) * 2.0;
    q.w = (m21 - m12) / S;
    q.x = 0.25 * S;
    q.y = (m01 + m10) / S;
    q.z = (m02 + m20) / S;
  } else if (m11 > m22) {
    const double S = std::sqrt(1.0 + m11 - m00 - m22) * 2.0;
    q.w = (m02 - m20) / S;
    q.x = (m01 + m10) / S;
    q.y = 0.25 * S;
    q.z = (m12 + m21) / S;
  } else {
    const double S = std::sqrt(1.0 + m22 - m00 - m11) * 2.0;
    q.w = (m10 - m01) / S;
    q.x = (m02 + m20) / S;
    q.y = (m12 + m21) / S;
    q.z = 0.25 * S;
  }

  ship.setOrientation(q.normalized());
}

static bool stationDockingGeometry(const sim::Station& st,
                                  double timeDays,
                                  math::Vec3d& outPosKm,
                                  math::Vec3d& outVelKmS,
                                  math::Vec3d& outApproachDir) {
  outPosKm = orbitPosKm(st.orbit, timeDays);
  outVelKmS = orbitVelKmS(st.orbit, timeDays);
  outApproachDir = safeNorm(outPosKm);
  return true;
}

static bool isInDockingCorridor(const sim::Station& st,
                               const math::Vec3d& stationPosKm,
                               const math::Vec3d& approachDir,
                               const math::Vec3d& shipPosKm,
                               double& outAxial,
                               double& outRadial) {
  const math::Vec3d rel = shipPosKm - stationPosKm;
  const double axial = math::dot(rel, approachDir);
  const math::Vec3d radialV = rel - approachDir * axial;
  const double radial = radialV.length();

  outAxial = axial;
  outRadial = radial;

  return axial >= 0.0 && axial <= st.docking.corridorLengthKm && radial <= st.docking.corridorRadiusKm;
}

static bool canDockNow(const sim::Station& st,
                       const math::Vec3d& stationPosKm,
                       const math::Vec3d& stationVelKmS,
                       const math::Vec3d& approachDir,
                       const sim::Ship& ship,
                       const DockingClearance& clearance,
                       double timeDays,
                       double& outDistKm,
                       double& outRelSpeedKmS,
                       double& outAlignCos,
                       bool& outInCorridor) {
  const math::Vec3d rel = ship.positionKm() - stationPosKm;
  const math::Vec3d relV = ship.velocityKmS() - stationVelKmS;
  const double dist = rel.length();
  const double relSpeed = relV.length();

  double axial = 0.0;
  double radial = 0.0;
  const bool inCorridor = isInDockingCorridor(st, stationPosKm, approachDir, ship.positionKm(), axial, radial);

  // Ship should point "in" along -approachDir.
  const math::Vec3d requiredFwd = -approachDir;
  const double alignCos = math::dot(safeNorm(ship.forward()), safeNorm(requiredFwd));

  const bool clearanceOk = (clearance.state == ClearanceState::Granted &&
                            clearance.stationId == st.id &&
                            timeDays <= clearance.untilTimeDays);

  outDistKm = dist;
  outRelSpeedKmS = relSpeed;
  outAlignCos = alignCos;
  outInCorridor = inCorridor;

  return clearanceOk && inCorridor && dist <= st.docking.dockRangeKm && relSpeed <= st.docking.speedLimitKmS &&
         alignCos >= st.docking.alignCosMin;
}

static void stepShipSubdiv(sim::Ship& ship, double dtSeconds, const sim::ShipInput& input) {
  const double maxStep = 0.05; // seconds
  if (dtSeconds <= 0.0) return;
  const int steps = std::clamp((int)std::ceil(dtSeconds / maxStep), 1, 256);
  const double h = dtSeconds / (double)steps;
  for (int i = 0; i < steps; ++i) ship.step(h, input);
}

static double maxAllowedTimeScale(const sim::StarSystem& sys,
                                 double timeDays,
                                 const sim::Ship& ship,
                                 FlightMode mode) {
  if (mode == FlightMode::Supercruise) return 1.0;
  if (mode == FlightMode::Docked) return 12000.0;

  const math::Vec3d pos = ship.positionKm();
  const double speed = ship.velocityKmS().length();

  // Very high speeds: keep time scale low so controls remain usable.
  if (speed > 20.0) return 1.0;
  if (speed > 8.0) return 5.0;

  // Distance to star surface.
  const double starRadKm = sys.star.radiusSol * kSOLAR_RADIUS_KM;
  double minDistKm = std::max(0.0, pos.length() - starRadKm);

  // Planets
  for (const auto& p : sys.planets) {
    const math::Vec3d ppos = orbitPosKm(p.orbit, timeDays);
    const double rad = p.radiusEarth * kEARTH_RADIUS_KM;
    const double d = std::max(0.0, (ppos - pos).length() - rad);
    minDistKm = std::min(minDistKm, d);
  }

  // Stations
  for (const auto& st : sys.stations) {
    const math::Vec3d spos = orbitPosKm(st.orbit, timeDays);
    const double d = std::max(0.0, (spos - pos).length() - st.radiusKm);
    minDistKm = std::min(minDistKm, d);
  }

  // Safety tiers (km)
  if (minDistKm < 20.0) return 1.0;
  if (minDistKm < 80.0) return 10.0;
  if (minDistKm < 300.0) return 50.0;
  if (minDistKm < 1500.0) return 250.0;
  return 12000.0;
}

static bool beginStationSelectionUI(const sim::StarSystem& sys, int& stationIndex, bool docked, std::string_view dockedName) {
  bool changed = false;

  ImGui::Begin("Dock / Market");

  ImGui::Text("System: %s  (Star %s, planets %d, stations %d)",
              sys.stub.name.c_str(),
              starClassName(sys.stub.primaryClass),
              sys.stub.planetCount,
              sys.stub.stationCount);

  if (docked) {
    ImGui::TextColored(ImVec4(0.3f, 0.9f, 0.4f, 1.0f), "Docked at: %s", std::string(dockedName).c_str());
  } else {
    ImGui::TextDisabled("Not docked");
  }

  if (!sys.stations.empty()) {
    std::vector<const char*> names;
    names.reserve(sys.stations.size());
    for (const auto& st : sys.stations) names.push_back(st.name.c_str());

    int old = stationIndex;
    ImGui::Combo("Station", &stationIndex, names.data(), (int)names.size());
    changed = (old != stationIndex);

    const auto& st = sys.stations[(std::size_t)stationIndex];
    ImGui::SameLine();
    ImGui::TextDisabled("(%s, fee %.1f%%)", stationTypeName(st.type), st.feeRate * 100.0);
  } else {
    ImGui::Text("No stations in system.");
  }

  ImGui::End();
  return changed;
}

static double supercruiseMaxSpeedKmS(double minDistToMassKm) {
  // A simple "gravity well" curve:
  // - close to masses, allow a minimum high-speed cruise
  // - far away, asymptotically approach a big cap
  const double base = 200.0;       // km/s (near masses)
  const double cap = 100000.0;     // km/s (deep space)
  const double scale = 20.0e6;     // km
  const double t = 1.0 - std::exp(-std::max(0.0, minDistToMassKm) / scale);
  return base + (cap - base) * t;
}

static double minDistanceToMassKm(const sim::StarSystem& sys, double timeDays, const math::Vec3d& shipPosKm) {
  // Star at origin.
  const double starRadKm = sys.star.radiusSol * kSOLAR_RADIUS_KM;
  double minD = std::max(0.0, shipPosKm.length() - starRadKm);

  // Planets.
  for (const auto& p : sys.planets) {
    const math::Vec3d ppos = orbitPosKm(p.orbit, timeDays);
    const double rad = p.radiusEarth * kEARTH_RADIUS_KM;
    const double d = std::max(0.0, (ppos - shipPosKm).length() - rad);
    minD = std::min(minD, d);
  }

  // Stations.
  for (const auto& st : sys.stations) {
    const math::Vec3d spos = orbitPosKm(st.orbit, timeDays);
    const double d = std::max(0.0, (spos - shipPosKm).length() - st.radiusKm);
    minD = std::min(minD, d);
  }

  return minD;
}

static double stationDropDistanceKm(const sim::Station& st) {
  // Drop outside comms range so you're not instantly inside the "traffic" bubble.
  return std::max(80.0, st.docking.commsRangeKm * 4.0);
}

static double planetDropDistanceKm(const sim::Planet& p) {
  const double rad = p.radiusEarth * kEARTH_RADIUS_KM;
  return std::max(1500.0, rad * 6.0);
}

static bool getTargetState(const sim::StarSystem& sys,
                           const NavTarget& tgt,
                           double timeDays,
                           math::Vec3d& outPosKm,
                           math::Vec3d& outVelKmS,
                           std::string& outName) {
  outPosKm = {0, 0, 0};
  outVelKmS = {0, 0, 0};
  outName.clear();

  if (tgt.kind == TargetKind::Station) {
    if (tgt.index < 0 || tgt.index >= (int)sys.stations.size()) return false;
    const auto& st = sys.stations[(std::size_t)tgt.index];
    outPosKm = orbitPosKm(st.orbit, timeDays);
    outVelKmS = orbitVelKmS(st.orbit, timeDays);
    outName = st.name;
    return true;
  }
  if (tgt.kind == TargetKind::Planet) {
    if (tgt.index < 0 || tgt.index >= (int)sys.planets.size()) return false;
    const auto& p = sys.planets[(std::size_t)tgt.index];
    outPosKm = orbitPosKm(p.orbit, timeDays);
    outVelKmS = orbitVelKmS(p.orbit, timeDays);
    outName = p.name;
    return true;
  }
  return false;
}

static NavTarget cycleTarget(const sim::StarSystem& sys, const NavTarget& cur) {
  // Order: stations then planets.
  const int nStations = (int)sys.stations.size();
  const int nPlanets = (int)sys.planets.size();
  const int total = nStations + nPlanets;
  if (total <= 0) return {};

  int curFlat = -1;
  if (cur.kind == TargetKind::Station) curFlat = cur.index;
  if (cur.kind == TargetKind::Planet) curFlat = nStations + cur.index;

  int nextFlat = (curFlat + 1) % total;
  NavTarget out;
  if (nextFlat < nStations) {
    out.kind = TargetKind::Station;
    out.index = nextFlat;
  } else {
    out.kind = TargetKind::Planet;
    out.index = nextFlat - nStations;
  }
  return out;
}

static void requestDockingClearance(const sim::StarSystem& sys,
                                   int stationIndex,
                                   double timeDays,
                                   const math::Vec3d& shipPosKm,
                                   const math::Vec3d& shipVelKmS,
                                   DockingClearance& clearance,
                                   std::vector<Toast>& toasts) {
  if (stationIndex < 0 || stationIndex >= (int)sys.stations.size()) {
    pushToast(toasts, "No station selected.", ImGui::GetTime());
    return;
  }

  const auto& st = sys.stations[(std::size_t)stationIndex];
  const math::Vec3d spos = orbitPosKm(st.orbit, timeDays);
  const math::Vec3d svel = orbitVelKmS(st.orbit, timeDays);
  const double dist = (spos - shipPosKm).length();

  // Expire old
  if (clearance.state != ClearanceState::None && timeDays > clearance.untilTimeDays) {
    clearance = {};
  }

  if (dist > st.docking.commsRangeKm) {
    pushToast(toasts, "Docking request denied: out of range.", ImGui::GetTime());
    clearance.stationId = st.id;
    clearance.state = ClearanceState::Denied;
    clearance.untilTimeDays = timeDays + (st.docking.deniedCooldownSec / 86400.0);
    return;
  }

  if (clearance.state == ClearanceState::Denied && clearance.stationId == st.id && timeDays < clearance.untilTimeDays) {
    pushToast(toasts, "Docking request denied: try again soon.", ImGui::GetTime());
    return;
  }

  // Simple traffic model:
  // - prefer granting when you are slow relative to the station
  // - some chance of denial to keep it from being automatic
  const double relSpeed = (shipVelKmS - svel).length();

  double pGrant = 0.85;
  if (st.type == econ::StationType::Outpost) pGrant = 0.90;
  if (st.type == econ::StationType::Shipyard) pGrant = 0.80;

  if (relSpeed > st.docking.speedLimitKmS * 3.0) {
    pGrant = 0.15; // you're screaming in too fast
  }

  core::SplitMix64 rng(core::hashCombine((core::u64)st.id, (core::u64)std::floor(timeDays * 1440.0)));
  const double roll = rng.nextDouble();

  if (roll < pGrant) {
    clearance.stationId = st.id;
    clearance.state = ClearanceState::Granted;
    clearance.untilTimeDays = timeDays + (st.docking.clearanceDurationMin / 1440.0);
    pushToast(toasts, "Docking clearance granted.", ImGui::GetTime());
  } else {
    clearance.stationId = st.id;
    clearance.state = ClearanceState::Denied;
    clearance.untilTimeDays = timeDays + (st.docking.deniedCooldownSec / 86400.0);
    pushToast(toasts, "Docking clearance denied: hold position.", ImGui::GetTime());
  }
}

static bool undock(sim::Ship& ship,
                   const sim::StarSystem& sys,
                   double timeDays,
                   int dockedStationIndex,
                   std::vector<Toast>& toasts) {
  if (sys.stations.empty() || dockedStationIndex < 0 || dockedStationIndex >= (int)sys.stations.size()) return false;

  const auto& st = sys.stations[(std::size_t)dockedStationIndex];
  math::Vec3d spos, svel, approach;
  stationDockingGeometry(st, timeDays, spos, svel, approach);

  // Pop out slightly along the approach axis.
  ship.setPositionKm(spos + approach * (st.docking.corridorLengthKm + 5.0));
  ship.setVelocityKmS(svel + (-approach) * 0.2);
  ship.setAngularVelocityRadS({0, 0, 0});

  pushToast(toasts, "Undocked.", ImGui::GetTime());
  return true;
}

static bool attemptDock(sim::Ship& ship,
                        const sim::StarSystem& sys,
                        double timeDays,
                        DockingClearance& clearance,
                        sim::StationId& outDockedStation,
                        int& outDockedStationIndex,
                        std::vector<Toast>& toasts) {
  if (sys.stations.empty()) {
    pushToast(toasts, "No stations in this system.", ImGui::GetTime());
    return false;
  }

  // Find a station you're close to and that passes corridor checks.
  int best = -1;
  double bestDist = 1e300;

  for (int i = 0; i < (int)sys.stations.size(); ++i) {
    const auto& st = sys.stations[(std::size_t)i];
    math::Vec3d spos, svel, approach;
    stationDockingGeometry(st, timeDays, spos, svel, approach);

    double distKm = 0.0, relSpeed = 0.0, alignCos = 0.0;
    bool inCorridor = false;
    const bool ok = canDockNow(st, spos, svel, approach, ship, clearance, timeDays, distKm, relSpeed, alignCos, inCorridor);
    if (ok && distKm < bestDist) {
      bestDist = distKm;
      best = i;
    }
  }

  if (best < 0) {
    // Provide a more informative message: check the nearest station and report why.
    auto nearest = findNearestStationIndex(sys, timeDays, ship.positionKm());
    if (!nearest) {
      pushToast(toasts, "No station nearby.", ImGui::GetTime());
      return false;
    }

    const auto& st = sys.stations[(std::size_t)*nearest];
    math::Vec3d spos, svel, approach;
    stationDockingGeometry(st, timeDays, spos, svel, approach);

    double distKm = 0.0, relSpeed = 0.0, alignCos = 0.0;
    bool inCorridor = false;
    (void)canDockNow(st, spos, svel, approach, ship, clearance, timeDays, distKm, relSpeed, alignCos, inCorridor);

    if (!(clearance.state == ClearanceState::Granted && clearance.stationId == st.id && timeDays <= clearance.untilTimeDays)) {
      pushToast(toasts, "Docking failed: no clearance (press L to request).", ImGui::GetTime());
    } else if (!inCorridor) {
      pushToast(toasts, "Docking failed: not in approach corridor.", ImGui::GetTime());
    } else if (distKm > st.docking.dockRangeKm) {
      pushToast(toasts, "Docking failed: too far.", ImGui::GetTime());
    } else if (relSpeed > st.docking.speedLimitKmS) {
      pushToast(toasts, "Docking failed: too fast.", ImGui::GetTime());
    } else if (alignCos < st.docking.alignCosMin) {
      pushToast(toasts, "Docking failed: not aligned.", ImGui::GetTime());
    } else {
      pushToast(toasts, "Docking failed.", ImGui::GetTime());
    }
    return false;
  }

  const auto& st = sys.stations[(std::size_t)best];
  math::Vec3d spos, svel, approach;
  stationDockingGeometry(st, timeDays, spos, svel, approach);

  // Snap to docking position and match station velocity.
  ship.setPositionKm(spos + approach * st.docking.dockRangeKm);
  ship.setVelocityKmS(svel);
  ship.setAngularVelocityRadS({0, 0, 0});

  outDockedStation = st.id;
  outDockedStationIndex = best;

  pushToast(toasts, "Docked.", ImGui::GetTime());
  return true;
}

static sim::ShipInput autopilotStationInput(const sim::Station& st,
                                           double timeDays,
                                           const sim::Ship& ship,
                                           const DockingClearance& clearance,
                                           bool& outWantsDock) {
  outWantsDock = false;

  math::Vec3d spos, svel, approach;
  stationDockingGeometry(st, timeDays, spos, svel, approach);

  const math::Vec3d shipPos = ship.positionKm();
  const math::Vec3d shipVel = ship.velocityKmS();

  const math::Vec3d rel = shipPos - spos;
  const math::Vec3d relV = shipVel - svel;

  const double dist = rel.length();
  const double relSpeed = relV.length();

  // Corridor terms
  double axial = 0.0, radial = 0.0;
  const bool inCorridor = isInDockingCorridor(st, spos, approach, shipPos, axial, radial);

  // If we don't have clearance yet, hold outside the corridor.
  const bool clearanceOk =
      (clearance.state == ClearanceState::Granted && clearance.stationId == st.id && timeDays <= clearance.untilTimeDays);

  double targetAxial = st.docking.corridorLengthKm;
  if (clearanceOk) {
    targetAxial = st.docking.dockRangeKm;
  } else {
    // waiting point: just outside corridor
    targetAxial = st.docking.corridorLengthKm + std::max(5.0, st.docking.corridorRadiusKm * 1.2);
  }

  const math::Vec3d desiredPos = spos + approach * targetAxial;

  // Desired relative velocity: reduce position error with a time constant.
  const double tau = clearanceOk ? 2.0 : 4.0;
  math::Vec3d desiredRelVel = (desiredPos - shipPos) * (1.0 / std::max(0.5, tau));

  // Speed limits
  double maxRelSpeed = 6.0;
  if (clearanceOk && inCorridor) maxRelSpeed = st.docking.speedLimitKmS;

  // Smoothly reduce max rel speed as we get close.
  if (clearanceOk) {
    const double stopDist = std::max(1e-6, dist);
    const double a = std::max(0.01, ship.maxLinearAccelKmS2());
    const double vStop = std::sqrt(2.0 * a * stopDist);
    maxRelSpeed = std::min(maxRelSpeed, vStop);
  }

  if (desiredRelVel.length() > maxRelSpeed) {
    desiredRelVel = desiredRelVel.normalized() * maxRelSpeed;
  }

  const math::Vec3d desiredVel = svel + desiredRelVel;

  // PD on velocity
  const double tauV = 0.75;
  math::Vec3d desiredAcc = (desiredVel - shipVel) * (1.0 / std::max(0.1, tauV));

  // Clamp to ship capabilities (convert to thruster inputs).
  const double aCap = ship.maxLinearAccelKmS2();
  if (desiredAcc.length() > aCap) desiredAcc = desiredAcc.normalized() * aCap;

  // Convert world accel -> local thruster demand.
  const math::Vec3d accLocal = ship.orientation().conjugate().rotate(desiredAcc);
  sim::ShipInput out{};
  out.thrustLocal = {accLocal.x / std::max(1e-6, aCap), accLocal.y / std::max(1e-6, aCap), accLocal.z / std::max(1e-6, aCap)};

  // Attitude: align forward with -approachDir when near corridor.
  const math::Vec3d desiredFwdWorld = safeNorm(spos - shipPos); // point at station
  const math::Vec3d desiredFwdLocal = safeNorm(ship.orientation().conjugate().rotate(desiredFwdWorld));

  // Error axis between local forward and desired direction.
  const math::Vec3d fLocal = {0, 0, 1};
  const math::Vec3d axis = math::cross(fLocal, desiredFwdLocal);
  const double gain = 3.0;
  out.torqueLocal = axis * gain;

  // Clamp inputs
  out.thrustLocal.x = clampd(out.thrustLocal.x, -1.0, 1.0);
  out.thrustLocal.y = clampd(out.thrustLocal.y, -1.0, 1.0);
  out.thrustLocal.z = clampd(out.thrustLocal.z, -1.0, 1.0);

  out.torqueLocal.x = clampd(out.torqueLocal.x, -1.0, 1.0);
  out.torqueLocal.y = clampd(out.torqueLocal.y, -1.0, 1.0);
  out.torqueLocal.z = clampd(out.torqueLocal.z, -1.0, 1.0);

  out.dampers = true;

  // If we're basically ready, request docking.
  if (clearanceOk && dist <= st.docking.dockRangeKm * 1.25 && relSpeed <= st.docking.speedLimitKmS * 1.1) {
    outWantsDock = true;
  }

  return out;
}

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;

  core::setLogLevel(core::LogLevel::Info);

  if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) != 0) {
    core::log(core::LogLevel::Error, std::string("SDL_Init failed: ") + SDL_GetError());
    return 1;
  }

  // GL 3.3 core
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, 0);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
  SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

  SDL_Window* window = SDL_CreateWindow("Stellar Forge (prototype)",
                                        SDL_WINDOWPOS_CENTERED,
                                        SDL_WINDOWPOS_CENTERED,
                                        1280,
                                        720,
                                        SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);

  if (!window) {
    core::log(core::LogLevel::Error, std::string("SDL_CreateWindow failed: ") + SDL_GetError());
    return 1;
  }

  SDL_GLContext glContext = SDL_GL_CreateContext(window);
  SDL_GL_MakeCurrent(window, glContext);
  SDL_GL_SetSwapInterval(1);

  if (!render::gl::load()) {
    core::log(core::LogLevel::Error, "Failed to load OpenGL functions.");
    return 1;
  }

  core::log(core::LogLevel::Info, std::string("OpenGL: ") + render::gl::glVersionString());

  glEnable(GL_DEPTH_TEST);

  // --- ImGui setup ---
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
  ImGui::StyleColorsDark();

  ImGui_ImplSDL2_InitForOpenGL(window, glContext);
  ImGui_ImplOpenGL3_Init("#version 330 core");

  // --- Render assets ---
  render::Mesh sphere = render::Mesh::makeUvSphere(48, 24);
  render::Mesh cube = render::Mesh::makeCube();

  render::Texture2D checker;
  checker.createChecker(256, 256, 16);

  render::MeshRenderer meshRenderer;
  std::string err;
  if (!meshRenderer.init(&err)) {
    core::log(core::LogLevel::Error, err);
    return 1;
  }
  meshRenderer.setTexture(&checker);

  render::LineRenderer lineRenderer;
  if (!lineRenderer.init(&err)) {
    core::log(core::LogLevel::Error, err);
    return 1;
  }

  render::PointRenderer pointRenderer;
  if (!pointRenderer.init(&err)) {
    core::log(core::LogLevel::Error, err);
    return 1;
  }

  // --- Universe / sim state ---
  core::u64 seed = 1337;
  sim::Universe universe(seed);

  sim::SystemStub currentStub{};
  if (auto s = universe.findClosestSystem({0, 0, 0}, 80.0)) {
    currentStub = *s;
  } else {
    // Should be rare; fall back
    currentStub = sim::SystemStub{};
    currentStub.id = 0;
    currentStub.seed = seed;
    currentStub.name = "Origin";
    currentStub.posLy = {0, 0, 0};
    currentStub.planetCount = 6;
    currentStub.stationCount = 1;
  }

  const sim::StarSystem* currentSystem = &universe.getSystem(currentStub.id, &currentStub);

  sim::Ship ship;
  ship.setMaxLinearAccelKmS2(0.08);
  ship.setMaxAngularAccelRadS2(1.2);

  double timeDays = 0.0;

  // Time compression (sim seconds per real second)
  static const std::array<double, 11> kTimeLevels = {1.0, 2.0, 5.0, 10.0, 30.0, 60.0, 120.0, 300.0, 600.0, 2000.0, 12000.0};
  int timeLevelIndex = 5; // 60x
  double requestedTimeScale = kTimeLevels[(std::size_t)timeLevelIndex];
  double appliedTimeScale = requestedTimeScale;
  bool paused = false;

  // Player economy
  double credits = 2500.0;
  std::array<double, econ::kCommodityCount> cargo{};

  int selectedStationIndex = 0;

  // Nav / flight state
  FlightMode flightMode = FlightMode::Normal;
  sim::StationId dockedStationId = 0;
  int dockedStationIndex = -1;

  NavTarget target{};

  // Approach autopilot (normal flight)
  bool autopilotEnabled = false;

  // Supercruise state
  bool scAssistEnabled = false; // 7-second rule throttle hold + auto-drop
  double scThrottle = 1.0;
  double scSpeedKmS = 0.0;

  // Docking clearance
  DockingClearance clearance{};

  // UI state
  bool showGalaxy = true;
  bool showShip = true;
  bool showEconomy = true;
  bool showNav = true;
  bool showHud = true;

  bool cockpitCamera = false;

  // Toasts
  std::vector<Toast> toasts;

  // Save/load
  const std::string savePath = "savegame.txt";

  // Pre-generate starfield directions
  std::vector<math::Vec3d> starDirs;
  std::vector<math::Vec3d> starCols;
  std::vector<float> starSizes;
  {
    core::SplitMix64 rng(0xBADC0FFEEULL);
    const int n = 2200;
    starDirs.reserve(n);
    starCols.reserve(n);
    starSizes.reserve(n);
    for (int i = 0; i < n; ++i) {
      // random direction
      const double z = rng.range(-1.0, 1.0);
      const double t = rng.range(0.0, 2.0 * math::kPi);
      const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
      const double x = r * std::cos(t);
      const double y = r * std::sin(t);
      starDirs.push_back(safeNorm({x, y, z}));

      // subtle color variance
      const double c = rng.range(0.75, 1.0);
      const double warm = rng.range(0.0, 1.0);
      const math::Vec3d col = {c * (0.75 + 0.25 * warm), c * (0.78 + 0.18 * warm), c};
      starCols.push_back(col);

      const float sz = (float)rng.range(1.0, 2.5);
      starSizes.push_back(sz);
    }
  }

  // Initial spawn near first station (if any).
  if (!currentSystem->stations.empty()) {
    spawnNearStation(ship, *currentSystem, timeDays, 0);
  } else {
    ship.setPositionKm({0, 0, -8000.0});
  }

  bool running = true;
  auto last = std::chrono::high_resolution_clock::now();

  SDL_SetRelativeMouseMode(SDL_FALSE);

  while (running) {
    // Timing
    auto now = std::chrono::high_resolution_clock::now();
    const double dtReal = std::chrono::duration<double>(now - last).count();
    last = now;

    // Events
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      ImGui_ImplSDL2_ProcessEvent(&event);

      if (event.type == SDL_QUIT) running = false;
      if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_CLOSE) running = false;

      if (event.type == SDL_KEYDOWN && !event.key.repeat) {
        const bool captureKeys = io.WantCaptureKeyboard;

        if (event.key.keysym.sym == SDLK_ESCAPE) running = false;

        if (event.key.keysym.sym == SDLK_F5) {
          sim::SaveGame s{};
          s.version = 1;
          s.seed = universe.seed();
          s.timeDays = timeDays;
          s.currentSystem = currentSystem->stub.id;
          s.dockedStation = dockedStationId;

          s.shipPosKm = ship.positionKm();
          s.shipVelKmS = ship.velocityKmS();
          s.shipOrient = ship.orientation();
          s.shipAngVelRadS = ship.angularVelocityRadS();

          s.credits = credits;
          s.cargo = cargo;
          s.stationOverrides = universe.exportStationOverrides();

          if (sim::saveToFile(s, savePath)) {
            core::log(core::LogLevel::Info, "Saved to " + savePath);
            pushToast(toasts, "Saved.", ImGui::GetTime());
          }
        }

        if (event.key.keysym.sym == SDLK_F9) {
          sim::SaveGame s{};
          if (sim::loadFromFile(savePath, s)) {
            universe = sim::Universe(s.seed);
            universe.importStationOverrides(s.stationOverrides);

            timeDays = s.timeDays;

            const sim::StarSystem& sys = universe.getSystem(s.currentSystem);
            currentStub = sys.stub;
            currentSystem = &sys;

            ship.setPositionKm(s.shipPosKm);
            ship.setVelocityKmS(s.shipVelKmS);
            ship.setOrientation(s.shipOrient);
            ship.setAngularVelocityRadS(s.shipAngVelRadS);

            credits = s.credits;
            cargo = s.cargo;

            dockedStationId = s.dockedStation;
            flightMode = (dockedStationId != 0) ? FlightMode::Docked : FlightMode::Normal;

            dockedStationIndex = -1;
            if (dockedStationId != 0) {
              if (auto idx = findStationIndexById(*currentSystem, dockedStationId)) {
                dockedStationIndex = *idx;
                selectedStationIndex = *idx;
              }
            }

            // Reset other transient state
            target = {};
            autopilotEnabled = false;
            scAssistEnabled = false;
            scThrottle = 1.0;
            scSpeedKmS = 0.0;
            clearance = {};

            core::log(core::LogLevel::Info, "Loaded " + savePath);
            pushToast(toasts, "Loaded.", ImGui::GetTime());
          }
        }

        if (event.key.keysym.sym == SDLK_TAB && !captureKeys) showGalaxy = !showGalaxy;
        if (event.key.keysym.sym == SDLK_F1) showShip = !showShip;
        if (event.key.keysym.sym == SDLK_F2) showEconomy = !showEconomy;
        if (event.key.keysym.sym == SDLK_F3) showNav = !showNav;
        if (event.key.keysym.sym == SDLK_F4) showHud = !showHud;
        if (event.key.keysym.sym == SDLK_SPACE && !captureKeys) paused = !paused;

        if (!captureKeys) {
          if (event.key.keysym.sym == SDLK_v) cockpitCamera = !cockpitCamera;

          if (event.key.keysym.sym == SDLK_t) {
            target = cycleTarget(*currentSystem, target);
          }

          if (event.key.keysym.sym == SDLK_y) {
            target = {};
          }

          if (event.key.keysym.sym == SDLK_p) {
            if (flightMode == FlightMode::Supercruise) {
              scAssistEnabled = !scAssistEnabled;
              pushToast(toasts, scAssistEnabled ? "Supercruise assist ON" : "Supercruise assist OFF", ImGui::GetTime());
            } else {
              autopilotEnabled = !autopilotEnabled;
              pushToast(toasts, autopilotEnabled ? "Autopilot ON" : "Autopilot OFF", ImGui::GetTime());
            }
          }

          if (event.key.keysym.sym == SDLK_j) {
            if (flightMode == FlightMode::Docked) {
              pushToast(toasts, "Cannot enter supercruise while docked.", ImGui::GetTime());
            } else if (flightMode == FlightMode::Supercruise) {
              // drop out
              flightMode = FlightMode::Normal;
              // set a sane normal-space speed
              ship.setVelocityKmS(ship.forward() * 4.0);
              scSpeedKmS = 0.0;
              pushToast(toasts, "Exited supercruise.", ImGui::GetTime());
            } else {
              const double minD = minDistanceToMassKm(*currentSystem, timeDays, ship.positionKm());
              if (minD < 150.0) {
                pushToast(toasts, "Mass lock: too close to a body.", ImGui::GetTime());
              } else {
                flightMode = FlightMode::Supercruise;
                scThrottle = 1.0;
                scSpeedKmS = std::max(400.0, ship.velocityKmS().length());
                pushToast(toasts, "Entered supercruise.", ImGui::GetTime());
              }
            }
          }

          if (event.key.keysym.sym == SDLK_l) {
            int stIdx = -1;
            if (target.kind == TargetKind::Station) {
              stIdx = target.index;
            } else {
              auto n = findNearestStationIndex(*currentSystem, timeDays, ship.positionKm());
              if (n) stIdx = *n;
            }

            if (stIdx >= 0) {
              requestDockingClearance(*currentSystem,
                                     stIdx,
                                     timeDays,
                                     ship.positionKm(),
                                     ship.velocityKmS(),
                                     clearance,
                                     toasts);
            } else {
              pushToast(toasts, "No station to request docking with.", ImGui::GetTime());
            }
          }

          if (event.key.keysym.sym == SDLK_g) {
            if (flightMode == FlightMode::Docked) {
              if (undock(ship, *currentSystem, timeDays, dockedStationIndex, toasts)) {
                flightMode = FlightMode::Normal;
                dockedStationId = 0;
                dockedStationIndex = -1;
                // Clear clearance; require re-request next time.
                clearance = {};
              }
            } else {
              sim::StationId did = 0;
              int didIdx = -1;
              if (attemptDock(ship, *currentSystem, timeDays, clearance, did, didIdx, toasts)) {
                flightMode = FlightMode::Docked;
                dockedStationId = did;
                dockedStationIndex = didIdx;
                selectedStationIndex = didIdx;
                autopilotEnabled = false;
                scAssistEnabled = false;
              }
            }
          }

          if (event.key.keysym.sym == SDLK_PAGEUP) {
            timeLevelIndex = std::min<int>(timeLevelIndex + 1, (int)kTimeLevels.size() - 1);
          }
          if (event.key.keysym.sym == SDLK_PAGEDOWN) {
            timeLevelIndex = std::max<int>(timeLevelIndex - 1, 0);
          }
          if (event.key.keysym.sym == SDLK_HOME) {
            timeLevelIndex = 0;
          }
          if (event.key.keysym.sym == SDLK_END) {
            timeLevelIndex = (int)kTimeLevels.size() - 1;
          }
        }
      }
    }

    // Compute requested/applied time scale
    requestedTimeScale = kTimeLevels[(std::size_t)timeLevelIndex];
    const bool forceTime = (SDL_GetKeyboardState(nullptr)[SDL_SCANCODE_LCTRL] != 0) ||
                           (SDL_GetKeyboardState(nullptr)[SDL_SCANCODE_RCTRL] != 0);

    const double maxAllowed = maxAllowedTimeScale(*currentSystem, timeDays, ship, flightMode);
    appliedTimeScale = forceTime ? requestedTimeScale : std::min(requestedTimeScale, maxAllowed);

    // Expire clearance
    if (clearance.state != ClearanceState::None && timeDays > clearance.untilTimeDays) {
      clearance = {};
    }

    // Input
    sim::ShipInput input{};
    const Uint8* keys = SDL_GetKeyboardState(nullptr);
    const bool captureKeys = io.WantCaptureKeyboard;

    if (!captureKeys) {
      if (flightMode == FlightMode::Normal) {
        input.thrustLocal.z += (keys[SDL_SCANCODE_W] ? 1.0 : 0.0);
        input.thrustLocal.z -= (keys[SDL_SCANCODE_S] ? 1.0 : 0.0);

        input.thrustLocal.x += (keys[SDL_SCANCODE_D] ? 1.0 : 0.0);
        input.thrustLocal.x -= (keys[SDL_SCANCODE_A] ? 1.0 : 0.0);

        input.thrustLocal.y += (keys[SDL_SCANCODE_R] ? 1.0 : 0.0);
        input.thrustLocal.y -= (keys[SDL_SCANCODE_F] ? 1.0 : 0.0);

        input.torqueLocal.x += (keys[SDL_SCANCODE_UP] ? 1.0 : 0.0);
        input.torqueLocal.x -= (keys[SDL_SCANCODE_DOWN] ? 1.0 : 0.0);

        input.torqueLocal.y += (keys[SDL_SCANCODE_RIGHT] ? 1.0 : 0.0);
        input.torqueLocal.y -= (keys[SDL_SCANCODE_LEFT] ? 1.0 : 0.0);

        input.torqueLocal.z += (keys[SDL_SCANCODE_E] ? 1.0 : 0.0);
        input.torqueLocal.z -= (keys[SDL_SCANCODE_Q] ? 1.0 : 0.0);

        input.boost = keys[SDL_SCANCODE_LSHIFT] != 0;
        input.brake = keys[SDL_SCANCODE_X] != 0;

        static bool dampers = true;
        if (keys[SDL_SCANCODE_Z]) dampers = true;
        if (keys[SDL_SCANCODE_C]) dampers = false;
        input.dampers = dampers;
      } else if (flightMode == FlightMode::Supercruise) {
        // Rotation still works.
        input.torqueLocal.x += (keys[SDL_SCANCODE_UP] ? 1.0 : 0.0);
        input.torqueLocal.x -= (keys[SDL_SCANCODE_DOWN] ? 1.0 : 0.0);

        input.torqueLocal.y += (keys[SDL_SCANCODE_RIGHT] ? 1.0 : 0.0);
        input.torqueLocal.y -= (keys[SDL_SCANCODE_LEFT] ? 1.0 : 0.0);

        input.torqueLocal.z += (keys[SDL_SCANCODE_E] ? 1.0 : 0.0);
        input.torqueLocal.z -= (keys[SDL_SCANCODE_Q] ? 1.0 : 0.0);

        input.dampers = true;

        // Manual throttle when not using SC assist
        if (!scAssistEnabled) {
          const double d = 0.6 * dtReal;
          if (keys[SDL_SCANCODE_W]) scThrottle = clampd(scThrottle + d, 0.0, 1.0);
          if (keys[SDL_SCANCODE_S]) scThrottle = clampd(scThrottle - d, 0.0, 1.0);
          if (keys[SDL_SCANCODE_X]) scThrottle = 0.0;
        }
      }
    }

    // Autopilot override (normal flight)
    bool wantsAutoDock = false;
    if (!paused && flightMode == FlightMode::Normal && autopilotEnabled && target.kind == TargetKind::Station &&
        target.index >= 0 && target.index < (int)currentSystem->stations.size()) {
      const auto& st = currentSystem->stations[(std::size_t)target.index];

      // Auto-request docking when in comms range
      const math::Vec3d spos = orbitPosKm(st.orbit, timeDays);
      const double dist = (spos - ship.positionKm()).length();
      if (dist <= st.docking.commsRangeKm &&
          !(clearance.state == ClearanceState::Granted && clearance.stationId == st.id && timeDays <= clearance.untilTimeDays)) {
        requestDockingClearance(*currentSystem,
                               target.index,
                               timeDays,
                               ship.positionKm(),
                               ship.velocityKmS(),
                               clearance,
                               toasts);
      }

      input = autopilotStationInput(st, timeDays, ship, clearance, wantsAutoDock);
    }

    // Simulation step
    if (!paused) {
      const double dtSim = dtReal * appliedTimeScale;

      if (flightMode == FlightMode::Normal) {
        stepShipSubdiv(ship, dtSim, input);
        timeDays += dtSim / 86400.0;

        if (wantsAutoDock) {
          sim::StationId did = 0;
          int didIdx = -1;
          if (attemptDock(ship, *currentSystem, timeDays, clearance, did, didIdx, toasts)) {
            flightMode = FlightMode::Docked;
            dockedStationId = did;
            dockedStationIndex = didIdx;
            selectedStationIndex = didIdx;
            autopilotEnabled = false;
            scAssistEnabled = false;
          }
        }
      } else if (flightMode == FlightMode::Docked) {
        // Advance time, keep ship attached.
        timeDays += (dtSim / 86400.0);

        if (dockedStationIndex >= 0 && dockedStationIndex < (int)currentSystem->stations.size()) {
          const auto& st = currentSystem->stations[(std::size_t)dockedStationIndex];
          math::Vec3d spos, svel, approach;
          stationDockingGeometry(st, timeDays, spos, svel, approach);
          ship.setPositionKm(spos + approach * st.docking.dockRangeKm);
          ship.setVelocityKmS(svel);
        }
      } else if (flightMode == FlightMode::Supercruise) {
        // In supercruise, we do NOT apply time acceleration (to avoid insane stacking).
        timeDays += dtReal / 86400.0;

        // Rotate ship (no linear motion through ship.step)
        ship.setVelocityKmS({0, 0, 0});
        ship.setPositionKm(ship.positionKm());
        stepShipSubdiv(ship, dtReal, input);

        // Steering assist: keep nose pointed at target when SC assist is enabled.
        math::Vec3d targetPos{}, targetVel{};
        std::string targetName;
        const bool hasTarget = getTargetState(*currentSystem, target, timeDays, targetPos, targetVel, targetName);

        if (scAssistEnabled && hasTarget) {
          const math::Vec3d dir = safeNorm(targetPos - ship.positionKm());
          const math::Vec3d dirLocal = safeNorm(ship.orientation().conjugate().rotate(dir));
          const math::Vec3d fLocal = {0, 0, 1};
          const math::Vec3d axis = math::cross(fLocal, dirLocal);
          sim::ShipInput steer{};
          steer.dampers = true;
          steer.torqueLocal = axis * 2.5;
          steer.torqueLocal.x = clampd(steer.torqueLocal.x, -1.0, 1.0);
          steer.torqueLocal.y = clampd(steer.torqueLocal.y, -1.0, 1.0);
          steer.torqueLocal.z = clampd(steer.torqueLocal.z, -1.0, 1.0);
          stepShipSubdiv(ship, dtReal, steer);

          // 7-second-rule throttle hold (aim to keep time-to-target ~ 7s)
          const double distKm = (targetPos - ship.positionKm()).length();
          const double etaTarget = 7.0;

          const double minD = minDistanceToMassKm(*currentSystem, timeDays, ship.positionKm());
          const double maxV = supercruiseMaxSpeedKmS(minD);

          const double desiredSpeed = clampd(distKm / std::max(etaTarget, 0.5), 0.0, maxV);
          scThrottle = clampd(desiredSpeed / std::max(1.0, maxV), 0.0, 1.0);
        }

        // Speed update
        const double minD = minDistanceToMassKm(*currentSystem, timeDays, ship.positionKm());
        const double maxV = supercruiseMaxSpeedKmS(minD);
        const double desired = scThrottle * maxV;
        const double k = 1.5; // response
        const double a = 1.0 - std::exp(-dtReal * k);
        scSpeedKmS = scSpeedKmS + (desired - scSpeedKmS) * a;

        // Move forward at scSpeed
        ship.setVelocityKmS(ship.forward() * scSpeedKmS);
        ship.setPositionKm(ship.positionKm() + ship.velocityKmS() * dtReal);

        // Auto-drop near target
        if (scAssistEnabled && hasTarget) {
          const math::Vec3d toT = targetPos - ship.positionKm();
          const double distKm = toT.length();
          const double alignCos = math::dot(safeNorm(ship.forward()), safeNorm(toT));

          double dropDist = 1500.0;
          if (target.kind == TargetKind::Station && target.index >= 0 && target.index < (int)currentSystem->stations.size()) {
            dropDist = stationDropDistanceKm(currentSystem->stations[(std::size_t)target.index]);
          }
          if (target.kind == TargetKind::Planet && target.index >= 0 && target.index < (int)currentSystem->planets.size()) {
            dropDist = planetDropDistanceKm(currentSystem->planets[(std::size_t)target.index]);
          }

          if (distKm <= dropDist && alignCos >= 0.985) {
            flightMode = FlightMode::Normal;
            // Give a manageable starting velocity in normal space.
            ship.setVelocityKmS(ship.forward() * 6.0);
            scSpeedKmS = 0.0;
            pushToast(toasts, "Auto-drop complete.", ImGui::GetTime());

            // Optional: continue into normal-space autopilot for station approaches.
            if (target.kind == TargetKind::Station) autopilotEnabled = true;
          }
        }
      }
    }

    // ---- Camera follow ----
    render::Camera cam;
    int w = 1280, h = 720;
    SDL_GetWindowSize(window, &w, &h);
    const double aspect = (h > 0) ? (double)w / (double)h : 16.0 / 9.0;

    cam.setPerspective(math::degToRad(60.0), aspect, 0.01, 20000.0);

    const math::Vec3d shipPosU = ship.positionKm() * (1.0 / kRENDER_UNIT_KM);

    if (!cockpitCamera) {
      const math::Vec3d back = ship.forward() * (-8.0);
      const math::Vec3d up = ship.up() * (3.0);
      cam.setPosition(shipPosU + back + up);
      cam.setOrientation(ship.orientation());
    } else {
      cam.setPosition(shipPosU + ship.up() * 0.2);
      cam.setOrientation(ship.orientation());
    }

    const math::Mat4d view = cam.viewMatrix();
    const math::Mat4d proj = cam.projectionMatrix();

    float viewF[16], projF[16];
    matToFloat(view, viewF);
    matToFloat(proj, projF);

    meshRenderer.setViewProj(viewF, projF);
    lineRenderer.setViewProj(viewF, projF);
    pointRenderer.setViewProj(viewF, projF);

    // ---- Starfield (camera-anchored) ----
    std::vector<render::PointVertex> stars;
    stars.reserve(starDirs.size());

    const double starRadiusU = 16000.0; // within far plane
    for (std::size_t i = 0; i < starDirs.size(); ++i) {
      const math::Vec3d pU = cam.position() + starDirs[i] * starRadiusU;
      const auto& c = starCols[i];
      stars.push_back({(float)pU.x, (float)pU.y, (float)pU.z, (float)c.x, (float)c.y, (float)c.z, starSizes[i]});
    }

    // ---- Build instances (star + planets + stations + ship) ----
    std::vector<render::InstanceData> spheres;
    spheres.reserve(1 + currentSystem->planets.size());

    // Star at origin
    {
      const double starRadiusKm = currentSystem->star.radiusSol * kSOLAR_RADIUS_KM;
      const float starScale = (float)std::max(0.8, (starRadiusKm / kRENDER_UNIT_KM) * 3.0);
      spheres.push_back({0, 0, 0, starScale, 1, 0, 0, 0, 1.0f, 0.95f, 0.75f});
    }

    // Planets
    for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
      const auto& p = currentSystem->planets[i];
      const math::Vec3d posKm = orbitPosKm(p.orbit, timeDays);
      const math::Vec3d posU = posKm * (1.0 / kRENDER_UNIT_KM);

      const double radiusKm = p.radiusEarth * kEARTH_RADIUS_KM;
      const float scale = (float)std::max(0.25, (radiusKm / kRENDER_UNIT_KM) * 200.0);

      // Simple color palette by type
      float cr = 0.6f, cg = 0.6f, cb = 0.6f;
      switch (p.type) {
        case sim::PlanetType::Rocky: cr = 0.6f; cg = 0.55f; cb = 0.5f; break;
        case sim::PlanetType::Desert: cr = 0.8f; cg = 0.7f; cb = 0.35f; break;
        case sim::PlanetType::Ocean: cr = 0.25f; cg = 0.45f; cb = 0.85f; break;
        case sim::PlanetType::Ice: cr = 0.7f; cg = 0.85f; cb = 0.95f; break;
        case sim::PlanetType::GasGiant: cr = 0.7f; cg = 0.55f; cb = 0.35f; break;
        default: break;
      }

      spheres.push_back({(float)posU.x,
                         (float)posU.y,
                         (float)posU.z,
                         scale,
                         1,
                         0,
                         0,
                         0,
                         cr,
                         cg,
                         cb});
    }

    // Station cubes
    std::vector<render::InstanceData> cubes;
    cubes.reserve(currentSystem->stations.size() + 1);

    for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
      const auto& st = currentSystem->stations[i];
      const math::Vec3d posKm = orbitPosKm(st.orbit, timeDays);
      const math::Vec3d posU = posKm * (1.0 / kRENDER_UNIT_KM);

      float cr = 0.7f, cg = 0.8f, cb = 0.9f;
      if ((int)i == selectedStationIndex) {
        cr = 1.0f;
        cg = 0.75f;
        cb = 0.35f;
      }
      if (target.kind == TargetKind::Station && target.index == (int)i) {
        cr = 1.0f;
        cg = 0.35f;
        cb = 0.35f;
      }
      if (dockedStationId != 0 && st.id == dockedStationId) {
        cr = 0.3f;
        cg = 0.95f;
        cb = 0.4f;
      }

      // Use a visual scale factor so stations are visible even though radiusKm is small.
      const float scale = (float)std::max(0.008, (st.radiusKm / kRENDER_UNIT_KM) * 6000.0);

      cubes.push_back({(float)posU.x,
                       (float)posU.y,
                       (float)posU.z,
                       scale,
                       1,
                       0,
                       0,
                       0,
                       cr,
                       cg,
                       cb});
    }

    // Ship cube (oriented)
    {
      render::InstanceData shipInst{};
      shipInst.px = (float)shipPosU.x;
      shipInst.py = (float)shipPosU.y;
      shipInst.pz = (float)shipPosU.z;
      shipInst.scale = 0.0065f;
      shipInst.qw = (float)ship.orientation().w;
      shipInst.qx = (float)ship.orientation().x;
      shipInst.qy = (float)ship.orientation().y;
      shipInst.qz = (float)ship.orientation().z;

      shipInst.cr = 0.9f;
      shipInst.cg = 0.9f;
      shipInst.cb = 1.0f;
      cubes.push_back(shipInst);
    }

    // Orbit lines
    std::vector<render::LineVertex> orbitLines;
    orbitLines.reserve(currentSystem->planets.size() * 128);

    for (const auto& p : currentSystem->planets) {
      const int seg = 96;
      math::Vec3d prev{};
      for (int s = 0; s <= seg; ++s) {
        const double t = (double)s / (double)seg * p.orbit.periodDays;
        const math::Vec3d posKm = orbitPosKm(p.orbit, t);
        const math::Vec3d posU = posKm * (1.0 / kRENDER_UNIT_KM);

        if (s > 0) {
          orbitLines.push_back({(float)prev.x, (float)prev.y, (float)prev.z, 0.25f, 0.25f, 0.25f});
          orbitLines.push_back({(float)posU.x, (float)posU.y, (float)posU.z, 0.25f, 0.25f, 0.25f});
        }
        prev = posU;
      }
    }

    // Station corridor line for selected station
    std::vector<render::LineVertex> corridorLines;
    if (!currentSystem->stations.empty() && selectedStationIndex >= 0 && selectedStationIndex < (int)currentSystem->stations.size()) {
      const auto& st = currentSystem->stations[(std::size_t)selectedStationIndex];
      math::Vec3d spos, svel, approach;
      stationDockingGeometry(st, timeDays, spos, svel, approach);
      const math::Vec3d a = (spos)* (1.0 / kRENDER_UNIT_KM);
      const math::Vec3d b = (spos + approach * st.docking.corridorLengthKm) * (1.0 / kRENDER_UNIT_KM);

      corridorLines.push_back({(float)a.x, (float)a.y, (float)a.z, 0.8f, 0.7f, 0.2f});
      corridorLines.push_back({(float)b.x, (float)b.y, (float)b.z, 0.8f, 0.7f, 0.2f});
    }

    // ---- Render ----
    glViewport(0, 0, w, h);
    glClearColor(0.01f, 0.01f, 0.02f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    pointRenderer.drawPoints(stars);

    lineRenderer.drawLines(orbitLines);
    lineRenderer.drawLines(corridorLines);

    // Star + planets
    meshRenderer.setMesh(&sphere);
    meshRenderer.drawInstances(spheres);

    // Stations + ship
    meshRenderer.setMesh(&cube);
    meshRenderer.drawInstances(cubes);

    // ---- UI ----
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplSDL2_NewFrame(window);
    ImGui::NewFrame();

    ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());

    // HUD overlay
    if (showHud) {
      ImGuiViewport* vp = ImGui::GetMainViewport();
      ImDrawList* dl = ImGui::GetForegroundDrawList(vp);

      const ImVec2 center = ImVec2(vp->Pos.x + vp->Size.x * 0.5f, vp->Pos.y + vp->Size.y * 0.5f);
      dl->AddLine(ImVec2(center.x - 10, center.y), ImVec2(center.x + 10, center.y), IM_COL32(220, 220, 235, 160), 1.0f);
      dl->AddLine(ImVec2(center.x, center.y - 10), ImVec2(center.x, center.y + 10), IM_COL32(220, 220, 235, 160), 1.0f);

      // Target marker
      math::Vec3d tposKm{}, tvelKmS{};
      std::string tname;
      if (getTargetState(*currentSystem, target, timeDays, tposKm, tvelKmS, tname)) {
        const math::Vec3d tposU = tposKm * (1.0 / kRENDER_UNIT_KM);
        ImVec2 p2;
        bool behind = false;
        projectToScreen(tposU, view, proj, (int)vp->Size.x, (int)vp->Size.y, p2, behind);

        // clamp to viewport
        ImVec2 marker = ImVec2(vp->Pos.x + p2.x, vp->Pos.y + p2.y);
        if (behind) {
          // flip to edge
          marker.x = vp->Pos.x + vp->Size.x - 30.0f;
          marker.y = vp->Pos.y + 30.0f;
        }

        dl->AddCircle(marker, 10.0f, IM_COL32(255, 120, 120, 220), 0, 2.0f);
        const double distKm = (tposKm - ship.positionKm()).length();

        char buf[256];
        if (flightMode == FlightMode::Supercruise) {
          const double eta = distKm / std::max(1.0, std::max(scSpeedKmS, 1.0));
          std::snprintf(buf, sizeof(buf), "%s  %.0f km  ETA %.1fs", tname.c_str(), distKm, eta);
        } else {
          std::snprintf(buf, sizeof(buf), "%s  %.0f km", tname.c_str(), distKm);
        }

        dl->AddText(ImVec2(marker.x + 14.0f, marker.y - 8.0f), IM_COL32(235, 235, 245, 220), buf);
      }

      // Mode line
      const char* modeStr = (flightMode == FlightMode::Normal) ? "NORMAL" : (flightMode == FlightMode::Supercruise ? "SUPERCRUISE" : "DOCKED");
      char mline[256];
      std::snprintf(mline,
                    sizeof(mline),
                    "%s  |v| %.2f km/s  Time x%.0f%s",
                    modeStr,
                    ship.velocityKmS().length(),
                    appliedTimeScale,
                    paused ? " [PAUSED]" : "");
      dl->AddText(ImVec2(vp->Pos.x + 12.0f, vp->Pos.y + vp->Size.y - 26.0f), IM_COL32(235, 235, 245, 220), mline);

      // Clearance line
      if (clearance.state == ClearanceState::Granted) {
        dl->AddText(ImVec2(vp->Pos.x + 12.0f, vp->Pos.y + vp->Size.y - 46.0f), IM_COL32(120, 255, 150, 220), "Docking clearance: GRANTED");
      } else if (clearance.state == ClearanceState::Denied) {
        dl->AddText(ImVec2(vp->Pos.x + 12.0f, vp->Pos.y + vp->Size.y - 46.0f), IM_COL32(255, 120, 120, 220), "Docking clearance: DENIED");
      }
    }

    // Flight window
    if (showShip) {
      ImGui::Begin("Ship / Flight");

      const auto pos = ship.positionKm();
      const auto vel = ship.velocityKmS();
      const auto wv = ship.angularVelocityRadS();

      ImGui::Text("Time: %.4f days", timeDays);
      ImGui::Text("Time compression: requested x%.0f, applied x%.0f%s", requestedTimeScale, appliedTimeScale, (appliedTimeScale < requestedTimeScale ? " (clamped)" : ""));
      if ((SDL_GetKeyboardState(nullptr)[SDL_SCANCODE_LCTRL] != 0) || (SDL_GetKeyboardState(nullptr)[SDL_SCANCODE_RCTRL] != 0)) {
        ImGui::TextDisabled("CTRL held: force time compression (unsafe)");
      }

      ImGui::Separator();
      ImGui::Text("Pos (km):   [%.1f %.1f %.1f]", pos.x, pos.y, pos.z);
      ImGui::Text("Vel (km/s): [%.3f %.3f %.3f] |v|=%.3f", vel.x, vel.y, vel.z, vel.length());
      ImGui::Text("AngVel (rad/s): [%.3f %.3f %.3f]", wv.x, wv.y, wv.z);

      if (flightMode == FlightMode::Supercruise) {
        const double minD = minDistanceToMassKm(*currentSystem, timeDays, ship.positionKm());
        const double maxV = supercruiseMaxSpeedKmS(minD);
        ImGui::Separator();
        ImGui::Text("Supercruise speed: %.0f km/s  (max %.0f km/s)", scSpeedKmS, maxV);
        ImGui::Text("Throttle: %.0f%%", scThrottle * 100.0);
        ImGui::Text("SC assist: %s (P toggles)", scAssistEnabled ? "ON" : "OFF");
      }

      ImGui::Separator();
      ImGui::TextDisabled("Controls:");
      ImGui::BulletText("Translate: WASD + R/F (normal)");
      ImGui::BulletText("Rotate: Arrow keys + Q/E roll");
      ImGui::BulletText("Boost: LShift   Brake: X");
      ImGui::BulletText("Dampers: Z (on) / C (off)");
      ImGui::BulletText("Target: T cycle, Y clear");
      ImGui::BulletText("Autopilot/SC assist: P");
      ImGui::BulletText("Request docking: L");
      ImGui::BulletText("Dock/Undock: G");
      ImGui::BulletText("Supercruise: J");
      ImGui::BulletText("Time compression: PgUp/PgDn, Home=1x, End=max (hold CTRL to force)");
      ImGui::BulletText("Camera: V");
      ImGui::BulletText("Pause: Space   Save: F5   Load: F9");

      ImGui::End();
    }

    // Nav window
    if (showNav) {
      ImGui::Begin("Nav / Targets");

      if (ImGui::CollapsingHeader("Stations", ImGuiTreeNodeFlags_DefaultOpen)) {
        for (int i = 0; i < (int)currentSystem->stations.size(); ++i) {
          const auto& st = currentSystem->stations[(std::size_t)i];
          const bool isSel = (target.kind == TargetKind::Station && target.index == i);
          if (ImGui::Selectable(st.name.c_str(), isSel)) {
            target.kind = TargetKind::Station;
            target.index = i;
            selectedStationIndex = i;
          }
        }
      }

      if (ImGui::CollapsingHeader("Planets", ImGuiTreeNodeFlags_DefaultOpen)) {
        for (int i = 0; i < (int)currentSystem->planets.size(); ++i) {
          const auto& p = currentSystem->planets[(std::size_t)i];
          const bool isSel = (target.kind == TargetKind::Planet && target.index == i);
          if (ImGui::Selectable(p.name.c_str(), isSel)) {
            target.kind = TargetKind::Planet;
            target.index = i;
          }
        }
      }

      ImGui::Separator();
      ImGui::TextDisabled("Tip: target a station and hit P for approach autopilot.");

      ImGui::End();
    }

    // Economy window
    if (showEconomy) {
      std::string_view dockedName = "";
      if (flightMode == FlightMode::Docked && dockedStationIndex >= 0 && dockedStationIndex < (int)currentSystem->stations.size()) {
        dockedName = currentSystem->stations[(std::size_t)dockedStationIndex].name;
      }

      beginStationSelectionUI(*currentSystem, selectedStationIndex, flightMode == FlightMode::Docked, dockedName);

      if (!currentSystem->stations.empty()) {
        selectedStationIndex = std::clamp(selectedStationIndex, 0, (int)currentSystem->stations.size() - 1);
        const auto& station = currentSystem->stations[(std::size_t)selectedStationIndex];
        auto& stEcon = universe.stationEconomy(station, timeDays);

        const bool dockedHere = (flightMode == FlightMode::Docked && dockedStationId == station.id);

        ImGui::Begin("Market Details");
        ImGui::Text("Credits: %.2f", credits);
        ImGui::Text("Trading: %s", dockedHere ? "ENABLED (docked)" : "DISABLED (dock to trade)");

        static int selectedCommodity = 0;
        ImGui::SliderInt("Plot commodity", &selectedCommodity, 0, (int)econ::kCommodityCount - 1);

        if (ImGui::BeginTable("market", 6, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
          ImGui::TableSetupColumn("Commodity");
          ImGui::TableSetupColumn("Inv");
          ImGui::TableSetupColumn("Ask");
          ImGui::TableSetupColumn("Bid");
          ImGui::TableSetupColumn("Cargo");
          ImGui::TableSetupColumn("Trade");
          ImGui::TableHeadersRow();

          for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
            const auto cid = (econ::CommodityId)i;
            const auto q = econ::quote(stEcon, station.economyModel, cid, 0.10);

            ImGui::TableNextRow();

            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%s", std::string(econ::commodityName(cid)).c_str());

            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%.0f", q.inventory);

            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%.2f", q.ask);

            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%.2f", q.bid);

            ImGui::TableSetColumnIndex(4);
            ImGui::Text("%.0f", cargo[i]);

            ImGui::TableSetColumnIndex(5);
            ImGui::PushID((int)i);

            static float qty[(int)econ::kCommodityCount] = {};
            if (qty[i] <= 0.0f) qty[i] = 10.0f;
            ImGui::SetNextItemWidth(70);
            ImGui::InputFloat("##qty", &qty[i], 1.0f, 10.0f, "%.0f");

            ImGui::SameLine();
            if (!dockedHere) ImGui::BeginDisabled();
            if (ImGui::SmallButton("Buy")) {
              auto tr = econ::buy(stEcon, station.economyModel, cid, qty[i], credits, 0.10, station.feeRate);
              if (tr.ok) cargo[i] += qty[i];
            }
            if (!dockedHere) ImGui::EndDisabled();

            ImGui::SameLine();
            if (!dockedHere) ImGui::BeginDisabled();
            if (ImGui::SmallButton("Sell")) {
              const double sellUnits = std::min<double>(qty[i], cargo[i]);
              if (sellUnits > 0.0) {
                auto tr = econ::sell(stEcon, station.economyModel, cid, sellUnits, credits, 0.10, station.feeRate);
                if (tr.ok) cargo[i] -= sellUnits;
              }
            }
            if (!dockedHere) ImGui::EndDisabled();

            ImGui::PopID();
          }

          ImGui::EndTable();
        }

        // Price history plot for selected commodity
        const std::size_t cidx = (std::size_t)selectedCommodity;
        const auto& hist = stEcon.history[cidx];
        if (!hist.empty()) {
          std::vector<float> vals;
          vals.reserve(hist.size());
          for (const auto& p : hist) vals.push_back((float)p.price);

          ImGui::PlotLines("Price history", vals.data(), (int)vals.size(), 0, nullptr, 0.0f, 0.0f, ImVec2(0, 120));
        } else {
          ImGui::TextDisabled("No history yet (time needs to advance).");
        }

        ImGui::End();
      }
    }

    // Galaxy window
    if (showGalaxy) {
      ImGui::Begin("Galaxy / Streaming");

      const auto center = currentSystem->stub.posLy;
      static float radius = 200.0f;
      ImGui::SliderFloat("Query radius (ly)", &radius, 20.0f, 1200.0f);

      auto nearby = universe.queryNearby(center, radius, 128);

      ImGui::Text("Nearby systems: %d", (int)nearby.size());

      const ImVec2 canvasSize = ImVec2(420, 420);
      ImGui::BeginChild("map", canvasSize, true, ImGuiWindowFlags_NoScrollbar);

      ImDrawList* draw = ImGui::GetWindowDrawList();
      const ImVec2 p0 = ImGui::GetCursorScreenPos();
      const ImVec2 p1 = ImVec2(p0.x + canvasSize.x, p0.y + canvasSize.y);
      const ImVec2 centerPx = ImVec2((p0.x + p1.x) * 0.5f, (p0.y + p1.y) * 0.5f);

      draw->AddRectFilled(p0, p1, IM_COL32(10, 10, 14, 255));
      draw->AddRect(p0, p1, IM_COL32(80, 80, 95, 255));

      auto toPx = [&](const math::Vec3d& posLy) -> ImVec2 {
        const math::Vec3d d = posLy - center;
        const float sx = (float)(d.x / (double)radius) * (canvasSize.x * 0.5f);
        const float sy = (float)(d.y / (double)radius) * (canvasSize.y * 0.5f);
        return ImVec2(centerPx.x + sx, centerPx.y + sy);
      };

      // Star lanes: connect each system to 3 nearest neighbors
      const int k = 3;
      for (std::size_t i = 0; i < nearby.size(); ++i) {
        struct N {
          std::size_t j;
          double d2;
        };
        std::vector<N> ns;
        ns.reserve(nearby.size());
        for (std::size_t j = 0; j < nearby.size(); ++j)
          if (j != i) {
            const auto di = nearby[j].posLy - nearby[i].posLy;
            const double d2 = di.x * di.x + di.y * di.y + di.z * di.z;
            ns.push_back({j, d2});
          }
        std::sort(ns.begin(), ns.end(), [](const N& a, const N& b) { return a.d2 < b.d2; });
        const int count = std::min<int>(k, (int)ns.size());

        const ImVec2 a = toPx(nearby[i].posLy);
        for (int n = 0; n < count; ++n) {
          const ImVec2 b = toPx(nearby[ns[n].j].posLy);
          draw->AddLine(a, b, IM_COL32(50, 80, 120, 100), 1.0f);
        }
      }

      static sim::SystemId selected = 0;
      for (const auto& s : nearby) {
        const ImVec2 p = toPx(s.posLy);
        const bool isCurrent = (s.id == currentSystem->stub.id);
        const bool isSel = (s.id == selected);

        ImU32 col = isCurrent ? IM_COL32(255, 240, 160, 255) : IM_COL32(170, 170, 190, 255);
        if (s.factionId != 0) col = IM_COL32(160, 220, 170, 255);
        if (isSel) col = IM_COL32(255, 120, 120, 255);

        draw->AddCircleFilled(p, isCurrent ? 5.5f : 4.0f, col);

        const float rClick = 6.0f;
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
          const ImVec2 mp = ImGui::GetIO().MousePos;
          const float dx = mp.x - p.x;
          const float dy = mp.y - p.y;
          if (dx * dx + dy * dy <= rClick * rClick) selected = s.id;
        }
      }

      ImGui::EndChild();

      if (selected != 0 && selected != currentSystem->stub.id) {
        if (ImGui::Button("Jump to selected system")) {
          auto it = std::find_if(nearby.begin(), nearby.end(), [&](const sim::SystemStub& s) { return s.id == selected; });
          if (it != nearby.end()) {
            currentStub = *it;
            const auto& sys = universe.getSystem(currentStub.id, &currentStub);
            currentSystem = &sys;

            // spawn near a station for immediate gameplay
            if (!currentSystem->stations.empty()) {
              spawnNearStation(ship, *currentSystem, timeDays, 0);
              selectedStationIndex = 0;
            } else {
              ship.setPositionKm({0, 0, -8000.0});
              ship.setVelocityKmS({0, 0, 0});
              ship.setAngularVelocityRadS({0, 0, 0});
            }

            flightMode = FlightMode::Normal;
            dockedStationId = 0;
            dockedStationIndex = -1;

            target = {};
            autopilotEnabled = false;
            scAssistEnabled = false;
            clearance = {};

            pushToast(toasts, "Arrived in new system.", ImGui::GetTime());
          }
        }
      }

      ImGui::TextDisabled("Tip: TAB toggles this window, F1 Flight, F2 Economy, F3 Nav, F4 HUD");

      ImGui::End();
    }

    // Toasts
    drawToasts(toasts, ImGui::GetTime());

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    SDL_GL_SwapWindow(window);
  }

  // Cleanup
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplSDL2_Shutdown();
  ImGui::DestroyContext();

  SDL_GL_DeleteContext(glContext);
  SDL_DestroyWindow(window);
  SDL_Quit();

  return 0;
}
