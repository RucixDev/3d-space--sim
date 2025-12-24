#include "stellar/core/Hash.h"
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
#include "stellar/sim/Mission.h"
#include "stellar/sim/Orbit.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/Ship.h"
#include "stellar/sim/Universe.h"

#include <SDL.h>
#include <SDL_opengl.h>

#include <backends/imgui_impl_opengl3.h>
#include <backends/imgui_impl_sdl2.h>
#include <imgui.h>

#include <array>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <optional>
#include <string>
#include <vector>

using namespace stellar;

// ---
// Units
// ---
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

static const char* missionTypeName(sim::MissionType t) {
  switch (t) {
    case sim::MissionType::Courier: return "Courier";
    case sim::MissionType::Delivery: return "Delivery";
    case sim::MissionType::BountyScan: return "Bounty Scan";
    default: return "?";
  }
}

// ---
// Math helpers
// ---

static math::Quatd quatFromBasis(const math::Vec3d& right,
                                 const math::Vec3d& up,
                                 const math::Vec3d& forward) {
  // Columns are basis vectors.
  const double m00 = right.x,  m01 = up.x,  m02 = forward.x;
  const double m10 = right.y,  m11 = up.y,  m12 = forward.y;
  const double m20 = right.z,  m21 = up.z,  m22 = forward.z;

  const double tr = m00 + m11 + m22;
  double qw, qx, qy, qz;
  if (tr > 0.0) {
    const double s = std::sqrt(tr + 1.0) * 2.0;
    qw = 0.25 * s;
    qx = (m21 - m12) / s;
    qy = (m02 - m20) / s;
    qz = (m10 - m01) / s;
  } else if (m00 > m11 && m00 > m22) {
    const double s = std::sqrt(1.0 + m00 - m11 - m22) * 2.0;
    qw = (m21 - m12) / s;
    qx = 0.25 * s;
    qy = (m01 + m10) / s;
    qz = (m02 + m20) / s;
  } else if (m11 > m22) {
    const double s = std::sqrt(1.0 + m11 - m00 - m22) * 2.0;
    qw = (m02 - m20) / s;
    qx = (m01 + m10) / s;
    qy = 0.25 * s;
    qz = (m12 + m21) / s;
  } else {
    const double s = std::sqrt(1.0 + m22 - m00 - m11) * 2.0;
    qw = (m10 - m01) / s;
    qx = (m02 + m20) / s;
    qy = (m12 + m21) / s;
    qz = 0.25 * s;
  }
  return math::Quatd{qw, qx, qy, qz}.normalized();
}

static math::Quatd quatLookRotation(const math::Vec3d& forward, const math::Vec3d& up) {
  math::Vec3d f = forward.normalized();
  if (f.length() < 1e-9) return math::Quatd::identity();

  math::Vec3d u = up.normalized();
  if (u.length() < 1e-9) u = {0,1,0};

  math::Vec3d r = math::cross(u, f).normalized();
  if (r.length() < 1e-9) {
    // forward parallel to up; pick an arbitrary right.
    r = math::cross({1,0,0}, f).normalized();
    if (r.length() < 1e-9) r = math::cross({0,0,1}, f).normalized();
  }
  u = math::cross(f, r);
  return quatFromBasis(r, u, f);
}

static math::Vec3d mulMat4Vec4(const math::Mat4d& m, const math::Vec3d& v, double w, double& outW) {
  // Column-major m; vector is (x,y,z,w).
  const double x = v.x, y = v.y, z = v.z;
  double rx = m.m[0]*x + m.m[4]*y + m.m[8]*z + m.m[12]*w;
  double ry = m.m[1]*x + m.m[5]*y + m.m[9]*z + m.m[13]*w;
  double rz = m.m[2]*x + m.m[6]*y + m.m[10]*z + m.m[14]*w;
  outW     = m.m[3]*x + m.m[7]*y + m.m[11]*z + m.m[15]*w;
  return {rx, ry, rz};
}

static bool projectToScreen(const math::Vec3d& worldPos,
                            const math::Mat4d& view,
                            const math::Mat4d& proj,
                            int screenW,
                            int screenH,
                            ImVec2& out) {
  double w1 = 1.0;
  const math::Vec3d vpos = mulMat4Vec4(view, worldPos, 1.0, w1);
  (void)w1;
  double w2 = 1.0;
  const math::Vec3d cpos = mulMat4Vec4(proj, vpos, 1.0, w2);
  if (w2 <= 1e-6) return false;
  const double ndcX = cpos.x / w2;
  const double ndcY = cpos.y / w2;
  const double ndcZ = cpos.z / w2;
  if (ndcZ < -1.0 || ndcZ > 1.0) return false;
  out.x = static_cast<float>((ndcX * 0.5 + 0.5) * screenW);
  out.y = static_cast<float>((1.0 - (ndcY * 0.5 + 0.5)) * screenH);
  return true;
}

// ---
// Nav targets
// ---

enum class NavTargetType {
  None = 0,
  Planet,
  Station,
};

struct NavTarget {
  NavTargetType type{NavTargetType::None};
  int index{-1};
};

static bool navTargetValid(const NavTarget& t, const sim::StarSystem& sys) {
  if (t.type == NavTargetType::Planet) return t.index >= 0 && t.index < (int)sys.planets.size();
  if (t.type == NavTargetType::Station) return t.index >= 0 && t.index < (int)sys.stations.size();
  return false;
}

static math::Vec3d orbitPosKm(const sim::OrbitElements& o, double timeDays) {
  const math::Vec3d au = sim::orbitPosition3DAU(o, timeDays);
  return au * kAU_KM;
}

static math::Vec3d orbitVelKmS(const sim::OrbitElements& o, double timeDays) {
  // Numerical derivative (1 second step).
  const double dtDays = 1.0 / 86400.0;
  const math::Vec3d p0 = orbitPosKm(o, timeDays);
  const math::Vec3d p1 = orbitPosKm(o, timeDays + dtDays);
  return (p1 - p0) * (1.0 / 1.0); // km/s
}

struct TargetKinematics {
  math::Vec3d posKm{0,0,0};
  math::Vec3d velKmS{0,0,0};
};

static std::optional<TargetKinematics> getTargetKinematics(const NavTarget& t,
                                                           const sim::StarSystem& sys,
                                                           double timeDays) {
  if (!navTargetValid(t, sys)) return std::nullopt;
  TargetKinematics k{};
  if (t.type == NavTargetType::Planet) {
    const auto& p = sys.planets[(std::size_t)t.index];
    k.posKm = orbitPosKm(p.orbit, timeDays);
    k.velKmS = orbitVelKmS(p.orbit, timeDays);
    return k;
  }
  if (t.type == NavTargetType::Station) {
    const auto& s = sys.stations[(std::size_t)t.index];
    k.posKm = orbitPosKm(s.orbit, timeDays);
    k.velKmS = orbitVelKmS(s.orbit, timeDays);
    return k;
  }
  return std::nullopt;
}

static const char* navTargetName(const NavTarget& t, const sim::StarSystem& sys) {
  if (!navTargetValid(t, sys)) return "None";
  if (t.type == NavTargetType::Planet) return sys.planets[(std::size_t)t.index].name.c_str();
  if (t.type == NavTargetType::Station) return sys.stations[(std::size_t)t.index].name.c_str();
  return "None";
}

// ---
// Toast / notifications
// ---

struct Toast {
  std::string text;
  float ttl{0.0f};
};

static void pushToast(std::vector<Toast>& toasts, std::string msg, float seconds = 3.0f) {
  toasts.push_back(Toast{std::move(msg), seconds});
}

// ---
// Mission helpers
// ---

static std::string formatMissionSummary(const sim::Mission& m, const sim::Universe& uni) {
  const auto& destSys = uni.getSystem(m.destSystem);
  std::string s = std::string(missionTypeName(m.type)) + ": " + destSys.stub.name;
  if (m.type == sim::MissionType::Delivery) {
    s += " (" + std::to_string((int)std::round(m.units)) + " " + std::string(econ::commodityName(m.commodity)) + ")";
  }
  return s;
}

static double missionDistanceLy(const sim::Universe& uni, const sim::Mission& m) {
  const auto& os = uni.getSystem(m.originSystem);
  const auto& ds = uni.getSystem(m.destSystem);
  return (ds.stub.posLy - os.stub.posLy).length();
}

static std::vector<sim::Mission> generateMissionOffers(const sim::Universe& uni,
                                                       const sim::StarSystem& originSys,
                                                       const sim::Station& originStation,
                                                       double timeDays) {
  // Deterministic mission offers refreshed daily.
  const core::u64 day = (core::u64)std::floor(timeDays);
  const core::u64 seed = core::hashCombine((core::u64)originStation.id, day);
  core::SplitMix64 rng(seed);

  std::vector<sim::Mission> out;
  out.reserve(8);

  // Candidate nearby systems.
  const auto near = uni.queryNearby(originSys.stub.posLy, 120.0, 48);
  if (near.empty()) return out;

  auto pickDest = [&]() -> sim::SystemId {
    for (int tries = 0; tries < 12; ++tries) {
      const auto& s = near[(std::size_t)rng.range<int>(0, (int)near.size() - 1)];
      if (s.id != originSys.stub.id) return s.id;
    }
    return originSys.stub.id;
  };

  auto pickStationId = [&](sim::SystemId sid) -> sim::StationId {
    const auto& sys = uni.getSystem(sid);
    if (sys.stations.empty()) return 0;
    const int si = rng.range<int>(0, (int)sys.stations.size() - 1);
    return sys.stations[(std::size_t)si].id;
  };

  const int offerCount = 6;
  for (int i = 0; i < offerCount; ++i) {
    sim::Mission m{};
    m.id = rng.nextU64();
    m.originSystem = originSys.stub.id;
    m.originStation = originStation.id;

    // Mix mission types.
    const double r = rng.nextDouble();
    if (r < 0.40) m.type = sim::MissionType::Courier;
    else if (r < 0.80) m.type = sim::MissionType::Delivery;
    else m.type = sim::MissionType::BountyScan;

    m.destSystem = pickDest();
    m.destStation = pickStationId(m.destSystem);
    if (m.destStation == 0) {
      // If destination has no stations, fallback to origin system.
      m.destSystem = originSys.stub.id;
      m.destStation = originStation.id;
    }

    const double distLy = missionDistanceLy(uni, m);
    const double base = 250.0 + distLy * 20.0;
    const double urgencyDays = rng.range(1.0, 6.0);
    m.expiryDay = (rng.chance(0.65) ? (std::floor(timeDays) + urgencyDays) : 0.0);

    if (m.type == sim::MissionType::Courier) {
      m.rewardCredits = base * rng.range(0.9, 1.3);
    } else if (m.type == sim::MissionType::Delivery) {
      const int cid = rng.range<int>(0, (int)econ::CommodityId::Count - 1);
      m.commodity = static_cast<econ::CommodityId>(cid);
      m.units = (double)rng.range<int>(4, 24);
      const double cargoValue = econ::commodityDef(m.commodity).basePrice * m.units;
      m.rewardCredits = (base + cargoValue * 0.25) * rng.range(0.9, 1.35);
    } else {
      // Bounty scan
      m.rewardCredits = (base * 0.75 + 120.0) * rng.range(0.9, 1.4);
      m.scanned = false;
    }

    out.push_back(std::move(m));
  }

  return out;
}

// ---
// Flight / docking state
// ---

enum class FlightMode {
  Normal = 0,
  Supercruise
};

struct DockingClearance {
  sim::StationId stationId{0};
  bool granted{false};
  float ttlSec{0.0f};
  float requestCooldownSec{0.0f};
  std::string lastMessage;
};

static void tickClearance(DockingClearance& c, float dtSec) {
  if (c.ttlSec > 0.0f) {
    c.ttlSec -= dtSec;
    if (c.ttlSec <= 0.0f) {
      c.granted = false;
      c.stationId = 0;
      c.ttlSec = 0.0f;
    }
  }
  if (c.requestCooldownSec > 0.0f) c.requestCooldownSec -= dtSec;
}

static bool requestDocking(const sim::Station& st,
                           const math::Vec3d& shipPosKm,
                           const math::Vec3d& stPosKm,
                           DockingClearance& c,
                           core::SplitMix64& rng,
                           std::vector<Toast>& toasts) {
  if (c.requestCooldownSec > 0.0f) return false;
  const double distKm = (shipPosKm - stPosKm).length();
  if (distKm > st.commsRangeKm) {
    pushToast(toasts, "Too far for docking request (out of comms range).", 2.5f);
    return false;
  }

  // Simple traffic heuristic: busier stations more likely to deny.
  double base = 0.85;
  if (st.type == econ::StationType::TradeHub) base = 0.75;
  if (st.type == econ::StationType::Shipyard) base = 0.78;
  if (st.type == econ::StationType::Outpost) base = 0.90;
  base -= st.feeRate * 0.15; // high-tax faction slightly less friendly
  base = std::clamp(base, 0.55, 0.95);

  const bool ok = rng.chance(base);
  c.stationId = st.id;
  c.granted = ok;
  c.ttlSec = ok ? 120.0f : 0.0f;
  c.requestCooldownSec = ok ? 2.5f : 7.5f;
  c.lastMessage = ok ? "Docking granted" : "Docking denied";
  pushToast(toasts, ok ? "Docking request accepted." : "Docking request denied.", 2.5f);
  return ok;
}

static double clamp01(double v) {
  return std::max(0.0, std::min(1.0, v));
}

static void drawHelpOverlay(bool* pOpen) {
  if (!ImGui::Begin("Help / Controls", pOpen, ImGuiWindowFlags_AlwaysAutoResize)) {
    ImGui::End();
    return;
  }

  ImGui::TextUnformatted("Flight:");
  ImGui::BulletText("W/S/A/D: forward/back/left/right thrust");
  ImGui::BulletText("Q/E: down/up thrust");
  ImGui::BulletText("Mouse: pitch/yaw (hold RMB)");
  ImGui::BulletText("Z/X: roll");
  ImGui::BulletText("Space: boost   Left Shift: brake");
  ImGui::BulletText("P: toggle approach autopilot");
  ImGui::BulletText("J: toggle supercruise");
  ImGui::BulletText("V: toggle supercruise approach assist (7s rule)");
  ImGui::Separator();
  ImGui::TextUnformatted("Navigation / docking:");
  ImGui::BulletText("T: cycle targets (stations -> planets)");
  ImGui::BulletText("L: request docking clearance (target station)");
  ImGui::BulletText("G: dock/undock");
  ImGui::BulletText("K: scan (for bounty missions / target scan)");
  ImGui::Separator();
  ImGui::TextUnformatted("Time warp:");
  ImGui::BulletText("[ / ]: decrease/increase time warp");
  ImGui::Separator();
  ImGui::TextUnformatted("UI:");
  ImGui::BulletText("F1: Flight HUD   F2: Dock/Market   F3: Galaxy/FSD   F4: Missions");
  ImGui::BulletText("F5: Save   F9: Load savegame.txt");

  ImGui::End();
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

  SDL_Window* window = SDL_CreateWindow(
      "Stellar Forge (prototype)",
      SDL_WINDOWPOS_CENTERED,
      SDL_WINDOWPOS_CENTERED,
      1280,
      720,
      SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
  if (!window) {
    core::log(core::LogLevel::Error, std::string("SDL_CreateWindow failed: ") + SDL_GetError());
    SDL_Quit();
    return 1;
  }

  SDL_GLContext glctx = SDL_GL_CreateContext(window);
  SDL_GL_MakeCurrent(window, glctx);
  SDL_GL_SetSwapInterval(1);

  if (!render::gl::loadGLFunctions(SDL_GL_GetProcAddress)) {
    core::log(core::LogLevel::Error, "Failed to load GL functions");
    SDL_GL_DeleteContext(glctx);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 1;
  }

  // ImGui
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
  ImGui::StyleColorsDark();
  ImGui_ImplSDL2_InitForOpenGL(window, glctx);
  ImGui_ImplOpenGL3_Init("#version 330");

  // Render resources
  render::Mesh sphere = render::Mesh::makeSphere(1.0f, 24, 16);
  render::Mesh cube = render::Mesh::makeCube();
  render::Texture2D whiteTex;
  whiteTex.createChecker(64, 64);

  render::MeshRenderer meshRenderer;
  std::string err;
  if (!meshRenderer.init(&err)) {
    core::log(core::LogLevel::Error, "MeshRenderer init failed: " + err);
  }
  meshRenderer.setTexture(&whiteTex);

  render::PointRenderer pointRenderer;
  if (!pointRenderer.init(&err)) {
    core::log(core::LogLevel::Error, "PointRenderer init failed: " + err);
  }
  render::LineRenderer lineRenderer;
  if (!lineRenderer.init(&err)) {
    core::log(core::LogLevel::Error, "LineRenderer init failed: " + err);
  }

  // Camera
  render::Camera cam;
  cam.setPosition({0, 0.02, -0.08});
  cam.setTarget({0, 0, 0});

  // ---
  // Load save
  // ---
  sim::SaveGame save;
  const std::string savePath = "savegame.txt";
  const bool hadSave = sim::loadFromFile(savePath, save);

  core::u64 seed = hadSave ? save.seed : 12345ull;
  double timeDays = hadSave ? save.timeDays : 0.0;
  double credits = hadSave ? save.credits : 1000.0;
  std::array<double, econ::kCommodityCount> cargo = hadSave ? save.cargo : std::array<double, econ::kCommodityCount>{};

  double hull = hadSave ? save.hull : 100.0;
  double hullMax = hadSave ? save.hullMax : 100.0;
  double fuel = hadSave ? save.fuel : 20.0;
  double fuelMax = hadSave ? save.fuelMax : 20.0;
  double cargoCapacity = hadSave ? save.cargoCapacity : 120.0;
  double fsdCooldownSec = hadSave ? save.fsdCooldownSec : 0.0;
  std::vector<sim::Mission> missions = hadSave ? save.missions : std::vector<sim::Mission>{};
  std::vector<sim::StationEconomyOverride> stationOverrides = hadSave ? save.stationOverrides : std::vector<sim::StationEconomyOverride>{};

  sim::Universe universe(seed);
  // apply station overrides loaded from save
  universe.importStationOverrides(stationOverrides);

  // Load / generate current system
  sim::SystemId currentSystemId = hadSave ? save.currentSystem : 0;
  const sim::StarSystem* currentSystem = nullptr;
  if (currentSystemId == 0) {
    // pick a deterministic nearby system near origin
    const auto list = universe.queryNearby({0,0,0}, 30.0, 16);
    currentSystemId = list.empty() ? universe.getSystem(0).stub.id : list.front().id;
  }
  currentSystem = &universe.getSystem(currentSystemId);

  // Ship
  sim::Ship ship;
  if (hadSave) {
    ship.setPositionKm(save.shipPosKm);
    ship.setVelocityKmS(save.shipVelKmS);
    ship.setOrientation(save.shipOrient);
    ship.setAngularVelocityRadS(save.shipAngVelRadS);
  } else {
    // Start near the first station if possible, else near star.
    if (!currentSystem->stations.empty()) {
      const auto& st = currentSystem->stations.front();
      const math::Vec3d stPos = orbitPosKm(st.orbit, timeDays);
      const math::Vec3d axis = (stPos.length() > 1e-6) ? stPos.normalized() : math::Vec3d{0,0,1};
      ship.setPositionKm(stPos + axis * (st.corridorLengthKm + 30.0));
      ship.setVelocityKmS(orbitVelKmS(st.orbit, timeDays));
      ship.setOrientation(quatLookRotation(-axis, {0,1,0}));
    } else {
      ship.setPositionKm({0, 0, -8000});
      ship.setVelocityKmS({0, 0, 0});
      ship.setOrientation(math::Quatd::identity());
    }
  }

  // Docked state
  sim::StationId dockedStationId = hadSave ? save.dockedStation : 0;
  bool isDocked = (dockedStationId != 0);

  // Selected station index (for UI)
  int selectedStationIndex = 0;
  if (!currentSystem->stations.empty()) {
    // If docked, select that station.
    if (isDocked) {
      for (int i = 0; i < (int)currentSystem->stations.size(); ++i) {
        if (currentSystem->stations[(std::size_t)i].id == dockedStationId) {
          selectedStationIndex = i;
          break;
        }
      }
    }
  }

  // Nav target
  NavTarget target;
  if (!currentSystem->stations.empty()) {
    target.type = NavTargetType::Station;
    target.index = selectedStationIndex;
  } else if (!currentSystem->planets.empty()) {
    target.type = NavTargetType::Planet;
    target.index = 0;
  }

  // Flight mode
  FlightMode flightMode = FlightMode::Normal;
  double supercruiseSpeedKmS = 0.0;
  double supercruiseThrottle = 0.0;
  bool supercruiseAssist = true;
  bool supercruiseAutoDrop = true;
  const double scMaxSpeedKmS = 6000.0;

  // Autopilot
  bool approachAutopilot = false;

  // Time warp levels (sim seconds per real second)
  const std::vector<double> timeWarpLevels = {1.0, 5.0, 10.0, 50.0, 200.0, 1000.0, 5000.0};
  int timeWarpIndex = 0;
  bool paused = false;

  // UI toggles
  bool showHUD = true;
  bool showDock = true;
  bool showGalaxy = false;
  bool showMissions = true;
  bool showHelp = false;

  // Galaxy selection
  math::Vec3d galQueryPos = currentSystem->stub.posLy;
  float galRadius = 50.0f;
  int galMax = 32;
  std::vector<sim::SystemStub> galList;
  int selectedSystemIndex = 0;
  int selectedDestStationIndex = 0;

  // Docking clearance
  DockingClearance clearance;

  // Toasts
  std::vector<Toast> toasts;

  // RNG for gameplay (requests, etc.)
  core::SplitMix64 gameRng(core::hashCombine(seed, 0xA5A5A5A5ull));

  // ---
  // Main loop
  // ---
  bool running = true;
  auto last = std::chrono::high_resolution_clock::now();

  // Mouse flight
  bool mouseLook = false;
  float mouseSensitivity = 0.0025f;

  while (running) {
    // timing
    auto now = std::chrono::high_resolution_clock::now();
    double dtReal = std::chrono::duration<double>(now - last).count();
    last = now;
    if (dtReal > 0.25) dtReal = 0.25; // avoid huge steps on debug breaks

    // events
    SDL_Event e;
    while (SDL_PollEvent(&e)) {
      ImGui_ImplSDL2_ProcessEvent(&e);
      if (e.type == SDL_QUIT) running = false;
      if (e.type == SDL_WINDOWEVENT && e.window.event == SDL_WINDOWEVENT_CLOSE && e.window.windowID == SDL_GetWindowID(window))
        running = false;

      if (e.type == SDL_KEYDOWN && !e.key.repeat) {
        const SDL_Keycode kc = e.key.keysym.sym;

        if (kc == SDLK_ESCAPE) {
          running = false;
        } else if (kc == SDLK_F1) {
          showHUD = !showHUD;
        } else if (kc == SDLK_F2) {
          showDock = !showDock;
        } else if (kc == SDLK_F3) {
          showGalaxy = !showGalaxy;
        } else if (kc == SDLK_F4) {
          showMissions = !showMissions;
        } else if (kc == SDLK_F5) {
          // Save
          sim::SaveGame s{};
          s.seed = seed;
          s.timeDays = timeDays;
          s.currentSystem = currentSystemId;
          s.dockedStation = isDocked ? dockedStationId : 0;
          s.shipPosKm = ship.positionKm();
          s.shipVelKmS = ship.velocityKmS();
          s.shipOrient = ship.orientation();
          s.shipAngVelRadS = ship.angularVelocityRadS();
          s.credits = credits;
          s.cargo = cargo;
          s.hull = hull;
          s.hullMax = hullMax;
          s.fuel = fuel;
          s.fuelMax = fuelMax;
          s.cargoCapacity = cargoCapacity;
          s.fsdCooldownSec = fsdCooldownSec;
          s.missions = missions;

          // persist station overrides
          s.stationOverrides = universe.exportStationOverrides();

          if (sim::saveToFile(s, savePath)) {
            pushToast(toasts, "Saved savegame.txt", 2.0f);
          } else {
            pushToast(toasts, "Save failed", 2.0f);
          }
        } else if (kc == SDLK_F9) {
          sim::SaveGame s2;
          if (sim::loadFromFile(savePath, s2)) {
            seed = s2.seed;
            timeDays = s2.timeDays;
            credits = s2.credits;
            cargo = s2.cargo;
            hull = s2.hull;
            hullMax = s2.hullMax;
            fuel = s2.fuel;
            fuelMax = s2.fuelMax;
            cargoCapacity = s2.cargoCapacity;
            fsdCooldownSec = s2.fsdCooldownSec;
            missions = s2.missions;
            stationOverrides = s2.stationOverrides;

            universe = sim::Universe(seed);
            universe.importStationOverrides(stationOverrides);

            currentSystemId = s2.currentSystem;
            currentSystem = &universe.getSystem(currentSystemId);

            ship.setPositionKm(s2.shipPosKm);
            ship.setVelocityKmS(s2.shipVelKmS);
            ship.setOrientation(s2.shipOrient);
            ship.setAngularVelocityRadS(s2.shipAngVelRadS);

            dockedStationId = s2.dockedStation;
            isDocked = (dockedStationId != 0);

            pushToast(toasts, "Loaded savegame.txt", 2.0f);
          } else {
            pushToast(toasts, "Load failed", 2.0f);
          }
        } else if (kc == SDLK_h) {
          showHelp = !showHelp;
        } else if (kc == SDLK_p) {
          approachAutopilot = !approachAutopilot;
          pushToast(toasts, approachAutopilot ? "Approach autopilot: ON" : "Approach autopilot: OFF", 1.5f);
        } else if (kc == SDLK_v) {
          supercruiseAssist = !supercruiseAssist;
          pushToast(toasts, supercruiseAssist ? "Supercruise assist: ON" : "Supercruise assist: OFF", 1.5f);
        } else if (kc == SDLK_j) {
          if (isDocked) {
            pushToast(toasts, "Can't enter supercruise while docked.", 2.0f);
          } else {
            if (flightMode == FlightMode::Normal) {
              // Engage supercruise
              flightMode = FlightMode::Supercruise;
              supercruiseSpeedKmS = ship.velocityKmS().length();
              supercruiseThrottle = clamp01(supercruiseSpeedKmS / 50.0);
              pushToast(toasts, "Supercruise: ENGAGED", 2.0f);
            } else {
              // Drop to normal space
              flightMode = FlightMode::Normal;
              // Keep some forward velocity but clamp
              const double v = std::min(1.5, supercruiseSpeedKmS * 0.002);
              ship.setVelocityKmS(ship.forward() * v);
              supercruiseSpeedKmS = 0.0;
              supercruiseThrottle = 0.0;
              pushToast(toasts, "Supercruise: DISENGAGED", 2.0f);
            }
          }
        } else if (kc == SDLK_g) {
          // Dock / undock
          if (isDocked) {
            isDocked = false;
            dockedStationId = 0;
            pushToast(toasts, "Undocked", 2.0f);
          } else {
            // attempt docking at target station
            if (target.type != NavTargetType::Station || !navTargetValid(target, *currentSystem)) {
              pushToast(toasts, "No station targeted.", 2.0f);
            } else {
              const auto& st = currentSystem->stations[(std::size_t)target.index];
              const math::Vec3d stPos = orbitPosKm(st.orbit, timeDays);
              const math::Vec3d stVel = orbitVelKmS(st.orbit, timeDays);
              const math::Vec3d relP = ship.positionKm() - stPos;
              const math::Vec3d relV = ship.velocityKmS() - stVel;
              const double dist = relP.length();
              const double speed = relV.length();
              const math::Vec3d axis = (stPos.length() > 1e-6) ? stPos.normalized() : math::Vec3d{0,0,1};
              const double align = math::dot(ship.forward().normalized(), (-axis));
              const bool inRange = dist < (st.radiusKm + 6.0);
              const bool slowEnough = speed < 0.08;
              const bool alignedEnough = align > st.corridorAlignCos;
              const bool hasClearance = clearance.granted && clearance.stationId == st.id;
              if (!hasClearance) {
                pushToast(toasts, "Docking not permitted: request clearance (L)", 2.5f);
              } else if (!inRange) {
                pushToast(toasts, "Too far to dock. Approach the station.", 2.0f);
              } else if (!slowEnough) {
                pushToast(toasts, "Reduce relative speed to dock.", 2.0f);
              } else if (!alignedEnough) {
                pushToast(toasts, "Align with approach corridor to dock.", 2.0f);
              } else {
                isDocked = true;
                dockedStationId = st.id;
                selectedStationIndex = target.index;
                // snap ship to docking position
                const math::Vec3d dockPos = stPos + axis * (st.radiusKm + 2.0);
                ship.setPositionKm(dockPos);
                ship.setVelocityKmS(stVel);
                ship.setOrientation(quatLookRotation(-axis, {0,1,0}));
                pushToast(toasts, "Docked", 2.0f);
              }
            }
          }
        } else if (kc == SDLK_l) {
          // Request docking
          if (target.type == NavTargetType::Station && navTargetValid(target, *currentSystem)) {
            const auto& st = currentSystem->stations[(std::size_t)target.index];
            const math::Vec3d stPos = orbitPosKm(st.orbit, timeDays);
            requestDocking(st, ship.positionKm(), stPos, clearance, gameRng, toasts);
          } else {
            pushToast(toasts, "Target a station first.", 2.0f);
          }
        } else if (kc == SDLK_k) {
          // Scan
          bool scannedSomething = false;
          if (target.type == NavTargetType::Station && navTargetValid(target, *currentSystem)) {
            const auto& st = currentSystem->stations[(std::size_t)target.index];
            const math::Vec3d stPos = orbitPosKm(st.orbit, timeDays);
            const double dist = (ship.positionKm() - stPos).length();
            if (dist < std::max(40.0, st.radiusKm + 20.0)) {
              // complete any matching bounty mission
              for (auto& m : missions) {
                if (m.type == sim::MissionType::BountyScan && !m.scanned && m.destSystem == currentSystemId && m.destStation == st.id) {
                  m.scanned = true;
                  credits += m.rewardCredits;
                  pushToast(toasts, "Bounty scan complete! Paid " + std::to_string((int)std::round(m.rewardCredits)) + " cr", 3.0f);
                  scannedSomething = true;
                }
              }
              if (!scannedSomething) {
                pushToast(toasts, "Scan complete (no active bounty here).", 2.0f);
              }
            } else {
              pushToast(toasts, "Too far to scan target.", 2.0f);
            }
          }
          if (!scannedSomething && (target.type != NavTargetType::Station)) {
            pushToast(toasts, "Nothing scannable targeted.", 2.0f);
          }
        } else if (kc == SDLK_t) {
          // Cycle nav targets: stations then planets
          if (!currentSystem->stations.empty()) {
            if (target.type != NavTargetType::Station) {
              target.type = NavTargetType::Station;
              target.index = 0;
            } else {
              target.index = (target.index + 1) % (int)currentSystem->stations.size();
            }
          } else if (!currentSystem->planets.empty()) {
            target.type = NavTargetType::Planet;
            target.index = (target.index + 1) % (int)currentSystem->planets.size();
          } else {
            target.type = NavTargetType::None;
            target.index = -1;
          }
        } else if (kc == SDLK_LEFTBRACKET) {
          timeWarpIndex = std::max(0, timeWarpIndex - 1);
        } else if (kc == SDLK_RIGHTBRACKET) {
          timeWarpIndex = std::min((int)timeWarpLevels.size() - 1, timeWarpIndex + 1);
        } else if (kc == SDLK_SPACE) {
          // handled as continuous in key state
        } else if (kc == SDLK_BACKQUOTE) {
          paused = !paused;
        }
      }

      if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_RIGHT) {
        mouseLook = true;
        SDL_SetRelativeMouseMode(SDL_TRUE);
      }
      if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_RIGHT) {
        mouseLook = false;
        SDL_SetRelativeMouseMode(SDL_FALSE);
      }
    }

    // Determine allowed time warp.
    double warp = paused ? 0.0 : timeWarpLevels[(std::size_t)timeWarpIndex];
    if (flightMode == FlightMode::Supercruise) {
      // supercruise already provides the travel speed; keep time warp at 1 for stability
      warp = std::min(warp, 1.0);
    }
    if (isDocked) {
      // allow time warp while docked, but cap to prevent crazy economy jumps by accident
      warp = std::min(warp, 200.0);
    }

    // Safety clamps: if close to a station or moving fast in normal space, clamp.
    if (flightMode == FlightMode::Normal && !isDocked) {
      // No high warp if actively thrusting/turning unless autopilot.
      const Uint8* ks = SDL_GetKeyboardState(nullptr);
      const bool manualInput = ks[SDL_SCANCODE_W] || ks[SDL_SCANCODE_S] || ks[SDL_SCANCODE_A] || ks[SDL_SCANCODE_D] ||
                               ks[SDL_SCANCODE_Q] || ks[SDL_SCANCODE_E] || ks[SDL_SCANCODE_Z] || ks[SDL_SCANCODE_X] || mouseLook;
      if (manualInput && !approachAutopilot) warp = std::min(warp, 10.0);

      const double speed = ship.velocityKmS().length();
      if (speed > 2.0) warp = std::min(warp, 10.0);

      // Distance to nearest station
      double nearestStationKm = 1e30;
      for (const auto& st : currentSystem->stations) {
        const math::Vec3d p = orbitPosKm(st.orbit, timeDays);
        const double d = (ship.positionKm() - p).length();
        nearestStationKm = std::min(nearestStationKm, d);
      }
      if (nearestStationKm < 600.0) warp = std::min(warp, 10.0);
      if (nearestStationKm < 200.0) warp = std::min(warp, 1.0);
    }

    const double dtSim = dtReal * warp;
    timeDays += dtSim / 86400.0;
    if (fsdCooldownSec > 0.0) fsdCooldownSec = std::max(0.0, fsdCooldownSec - dtReal);

    // Tick docking clearance timers (real-time, not warped)
    tickClearance(clearance, (float)dtReal);

    // Update toasts
    for (auto it = toasts.begin(); it != toasts.end();) {
      it->ttl -= (float)dtReal;
      if (it->ttl <= 0.0f) it = toasts.erase(it);
      else ++it;
    }

    // Update system pointer if needed
    currentSystem = &universe.getSystem(currentSystemId);

    // Kinematics for docked station attachment
    std::optional<TargetKinematics> dockKin;
    if (isDocked) {
      // Find dock station in current system
      for (int si = 0; si < (int)currentSystem->stations.size(); ++si) {
        if (currentSystem->stations[(std::size_t)si].id == dockedStationId) {
          const auto& st = currentSystem->stations[(std::size_t)si];
          dockKin = TargetKinematics{orbitPosKm(st.orbit, timeDays), orbitVelKmS(st.orbit, timeDays)};
          break;
        }
      }
      if (dockKin) {
        const math::Vec3d axis = (dockKin->posKm.length() > 1e-6) ? dockKin->posKm.normalized() : math::Vec3d{0,0,1};
        ship.setPositionKm(dockKin->posKm + axis * 10.0);
        ship.setVelocityKmS(dockKin->velKmS);
        ship.setAngularVelocityRadS({0,0,0});
      }
    }

    // ---
    // Input -> ship controls
    // ---

    sim::ShipInput input;
    input.dampers = true;

    // Keyboard state
    const Uint8* kstate = SDL_GetKeyboardState(nullptr);
    const bool boost = kstate[SDL_SCANCODE_SPACE];
    const bool brake = kstate[SDL_SCANCODE_LSHIFT] || kstate[SDL_SCANCODE_RSHIFT];
    input.boost = boost;
    input.brake = brake;

    // Translation
    if (!isDocked && flightMode == FlightMode::Normal) {
      if (kstate[SDL_SCANCODE_W]) input.thrustLocal.z += 1.0;
      if (kstate[SDL_SCANCODE_S]) input.thrustLocal.z -= 1.0;
      if (kstate[SDL_SCANCODE_D]) input.thrustLocal.x += 1.0;
      if (kstate[SDL_SCANCODE_A]) input.thrustLocal.x -= 1.0;
      if (kstate[SDL_SCANCODE_E]) input.thrustLocal.y += 1.0;
      if (kstate[SDL_SCANCODE_Q]) input.thrustLocal.y -= 1.0;
    }

    // Rotation
    if (!isDocked) {
      if (kstate[SDL_SCANCODE_Z]) input.torqueLocal.z += 1.0;
      if (kstate[SDL_SCANCODE_X]) input.torqueLocal.z -= 1.0;
    }

    // Mouse look (pitch/yaw)
    if (mouseLook && !isDocked) {
      int mx = 0, my = 0;
      SDL_GetRelativeMouseState(&mx, &my);
      input.torqueLocal.x += (double)(-my) * mouseSensitivity;
      input.torqueLocal.y += (double)(mx) * mouseSensitivity;
    }

    // Normalize inputs
    if (input.thrustLocal.length() > 1.0) input.thrustLocal = input.thrustLocal.normalized();
    if (input.torqueLocal.length() > 1.0) input.torqueLocal = input.torqueLocal.normalized();

    // ---
    // Autopilot + supercruise flight model
    // ---

    if (!isDocked) {
      if (flightMode == FlightMode::Normal) {
        // Approach autopilot
        if (approachAutopilot && target.type == NavTargetType::Station && navTargetValid(target, *currentSystem)) {
          const auto& st = currentSystem->stations[(std::size_t)target.index];
          const math::Vec3d stPos = orbitPosKm(st.orbit, timeDays);
          const math::Vec3d stVel = orbitVelKmS(st.orbit, timeDays);

          // Auto-request clearance once in comms range.
          const double distKm = (ship.positionKm() - stPos).length();
          if (!(clearance.granted && clearance.stationId == st.id) && distKm < st.commsRangeKm) {
            requestDocking(st, ship.positionKm(), stPos, clearance, gameRng, toasts);
          }

          const math::Vec3d axis = (stPos.length() > 1e-6) ? stPos.normalized() : math::Vec3d{0,0,1};
          const math::Vec3d entryPos = stPos + axis * st.corridorLengthKm;

          // Corridor metrics
          const math::Vec3d rel = ship.positionKm() - stPos;
          const double t = math::dot(rel, axis); // km from station along axis
          const math::Vec3d lateral = rel - axis * t;
          const double latDist = lateral.length();

          // Choose a moving target point (moves with station).
          math::Vec3d desiredPos = stPos;
          if (t < 0.0) {
            desiredPos = entryPos; // overshot behind station
          } else if (t > st.corridorLengthKm + 40.0) {
            desiredPos = entryPos;
          } else {
            // Within corridor length band
            if (latDist > st.corridorRadiusKm * 0.85) {
              desiredPos = stPos + axis * std::clamp(t, 0.0, st.corridorLengthKm);
            } else {
              desiredPos = stPos;
            }
          }

          // Desired velocity: match station + approach along axis when in corridor.
          math::Vec3d desiredVel = stVel;
          const bool inCorridor = (t >= 0.0 && t <= st.corridorLengthKm && latDist <= st.corridorRadiusKm);
          if (inCorridor) {
            const double distToDock = std::max(0.0, t);
            const double approachSpd = std::min(st.corridorSpeedLimitKmS, 0.02 + distToDock * 0.001);
            desiredVel = stVel + (-axis) * approachSpd;
          }

          // PD controller in world space
          const math::Vec3d posErr = ship.positionKm() - desiredPos;
          const math::Vec3d velErr = ship.velocityKmS() - desiredVel;
          const math::Vec3d desiredAccelWorld = (-posErr * 0.08) + (-velErr * 0.65);

          // Convert to local thrust command
          const double maxA = ship.maxLinearAccelKmS2();
          const math::Vec3d accelClamped = desiredAccelWorld.length() > maxA ? desiredAccelWorld.normalized() * maxA : desiredAccelWorld;
          const math::Quatd q = ship.orientation();
          const math::Vec3d accelLocal = q.conjugate().rotate(accelClamped) * (1.0 / std::max(1e-6, maxA));
          input.thrustLocal = accelLocal;
          if (input.thrustLocal.length() > 1.0) input.thrustLocal = input.thrustLocal.normalized();

          // Orientation: look toward station / entry
          math::Vec3d lookDir = (desiredPos - ship.positionKm());
          if (lookDir.length() < 1e-3) lookDir = -axis;
          lookDir = lookDir.normalized();
          const math::Vec3d localDir = q.conjugate().rotate(lookDir);
          const double yawErr = std::atan2(localDir.x, localDir.z);
          const double pitchErr = -std::atan2(localDir.y, localDir.z);
          input.torqueLocal.x = pitchErr * 2.2;
          input.torqueLocal.y = yawErr * 2.2;
          input.torqueLocal.z *= 0.2; // discourage roll
          if (input.torqueLocal.length() > 1.0) input.torqueLocal = input.torqueLocal.normalized();

          // Auto-dock when conditions met
          const math::Vec3d stVelNow = stVel;
          const math::Vec3d relP = ship.positionKm() - stPos;
          const math::Vec3d relV = ship.velocityKmS() - stVelNow;
          const double dist = relP.length();
          const double speed = relV.length();
          const double align = math::dot(ship.forward().normalized(), (-axis));
          const bool hasClearance = clearance.granted && clearance.stationId == st.id;
          if (hasClearance && dist < (st.radiusKm + 4.0) && speed < 0.06 && align > st.corridorAlignCos) {
            isDocked = true;
            dockedStationId = st.id;
            selectedStationIndex = target.index;
            ship.setPositionKm(stPos + axis * (st.radiusKm + 2.0));
            ship.setVelocityKmS(stVelNow);
            ship.setOrientation(quatLookRotation(-axis, {0,1,0}));
            pushToast(toasts, "Autopilot: Docked", 2.0f);
          }
        }

        // Step full physics
        ship.step(dtSim, input);
      } else {
        // Supercruise: simplified linear model + keep angular dynamics.
        ship.stepAngularOnly(dtReal, input);

        // Determine target distances for assist / max speed.
        double nearestBodyKm = (ship.positionKm()).length();
        for (const auto& p : currentSystem->planets) {
          const double d = (ship.positionKm() - orbitPosKm(p.orbit, timeDays)).length();
          nearestBodyKm = std::min(nearestBodyKm, d);
        }
        for (const auto& st : currentSystem->stations) {
          const double d = (ship.positionKm() - orbitPosKm(st.orbit, timeDays)).length();
          nearestBodyKm = std::min(nearestBodyKm, d);
        }

        const double maxSpeed = std::min(scMaxSpeedKmS,
                                         std::max(40.0, 30.0 + nearestBodyKm * 0.002));

        // Throttle controls (W/S) in supercruise
        if (kstate[SDL_SCANCODE_W]) supercruiseThrottle = clamp01(supercruiseThrottle + dtReal * 0.35);
        if (kstate[SDL_SCANCODE_S]) supercruiseThrottle = clamp01(supercruiseThrottle - dtReal * 0.35);

        // 7-second rule assist when approaching target
        if (supercruiseAssist && navTargetValid(target, *currentSystem)) {
          const auto tk = getTargetKinematics(target, *currentSystem, timeDays);
          if (tk) {
            const math::Vec3d toT = tk->posKm - ship.positionKm();
            const double distKm = toT.length();
            const double desiredT = 7.0;
            const double desiredSpeed = std::clamp(distKm / desiredT, 0.0, maxSpeed);

            // Smoothly drive speed toward desiredSpeed.
            const double err = desiredSpeed - supercruiseSpeedKmS;
            supercruiseSpeedKmS += err * dtReal * 1.5;
            supercruiseSpeedKmS = std::clamp(supercruiseSpeedKmS, 0.0, maxSpeed);
            supercruiseThrottle = clamp01(supercruiseSpeedKmS / maxSpeed);

            // Auto-drop when close
            double dropDistKm = 0.0;
            if (target.type == NavTargetType::Station) {
              const auto& st = currentSystem->stations[(std::size_t)target.index];
              dropDistKm = std::max(250.0, st.corridorLengthKm * 1.4);
            } else {
              dropDistKm = 5000.0;
            }
            if (supercruiseAutoDrop && distKm < dropDistKm && supercruiseSpeedKmS < 150.0) {
              flightMode = FlightMode::Normal;
              // Match target velocity and place near approach entry for station
              ship.setVelocityKmS(tk->velKmS + ship.forward() * 0.5);
              supercruiseSpeedKmS = 0.0;
              supercruiseThrottle = 0.0;
              pushToast(toasts, "Supercruise: auto-drop", 2.0f);
            }
          }
        } else {
          // Manual throttle -> speed
          const double targetSpeed = supercruiseThrottle * maxSpeed;
          supercruiseSpeedKmS += (targetSpeed - supercruiseSpeedKmS) * dtReal * 1.2;
          supercruiseSpeedKmS = std::clamp(supercruiseSpeedKmS, 0.0, maxSpeed);
        }

        ship.setVelocityKmS(ship.forward() * supercruiseSpeedKmS);
        ship.setPositionKm(ship.positionKm() + ship.velocityKmS() * dtReal);
      }
    }

    // ---
    // Simple station collision + speeding penalty in corridor
    // ---
    if (!isDocked && flightMode == FlightMode::Normal) {
      for (const auto& st : currentSystem->stations) {
        const math::Vec3d stPos = orbitPosKm(st.orbit, timeDays);
        const math::Vec3d stVel = orbitVelKmS(st.orbit, timeDays);
        math::Vec3d rel = ship.positionKm() - stPos;
        const double dist = rel.length();
        if (dist < st.radiusKm) {
          const math::Vec3d n = (dist > 1e-6) ? (rel * (1.0 / dist)) : math::Vec3d{0,0,1};
          // bounce in station frame
          math::Vec3d vRel = ship.velocityKmS() - stVel;
          const double impact = std::max(0.0, math::dot(vRel, -n));
          vRel = vRel + n * (impact * 1.8);
          ship.setVelocityKmS(stVel + vRel);
          ship.setPositionKm(stPos + n * (st.radiusKm + 0.25));

          const double dmg = std::clamp(impact * 18.0, 1.0, 35.0);
          hull = std::max(0.0, hull - dmg);
          pushToast(toasts, "Impact! Hull -" + std::to_string((int)std::round(dmg)) + "%", 2.0f);
        }

        // Speeding warning inside corridor
        const math::Vec3d axis = (stPos.length() > 1e-6) ? stPos.normalized() : math::Vec3d{0,0,1};
        const double t = math::dot(rel, axis);
        const math::Vec3d lat = rel - axis * t;
        if (t >= 0.0 && t <= st.corridorLengthKm && lat.length() <= st.corridorRadiusKm) {
          const double relSpd = (ship.velocityKmS() - stVel).length();
          if (relSpd > st.corridorSpeedLimitKmS * 1.15) {
            // gentle hull wear + small fine
            hull = std::max(0.0, hull - dtSim * 0.002);
            credits = std::max(0.0, credits - dtSim * 0.03 * (1.0 + st.feeRate));
          }
        }
      }
    }

    // Mission expiry + completion for delivery/courier
    for (auto it = missions.begin(); it != missions.end();) {
      bool remove = false;
      if (it->expiryDay > 0.0 && timeDays > it->expiryDay) {
        pushToast(toasts, "Mission failed (expired): " + formatMissionSummary(*it, universe), 3.0f);
        remove = true;
      }

      if (!remove) {
        if (it->type == sim::MissionType::Courier) {
          if (isDocked && dockedStationId == it->destStation && currentSystemId == it->destSystem) {
            credits += it->rewardCredits;
            pushToast(toasts, "Courier delivered! +" + std::to_string((int)std::round(it->rewardCredits)) + " cr", 3.0f);
            remove = true;
          }
        } else if (it->type == sim::MissionType::Delivery) {
          if (isDocked && dockedStationId == it->destStation && currentSystemId == it->destSystem) {
            const std::size_t cid = (std::size_t)it->commodity;
            if (cid < cargo.size() && cargo[cid] >= it->units) {
              cargo[cid] -= it->units;
              credits += it->rewardCredits;
              pushToast(toasts, "Delivery complete! +" + std::to_string((int)std::round(it->rewardCredits)) + " cr", 3.0f);
              remove = true;
            }
          }
        } else if (it->type == sim::MissionType::BountyScan) {
          if (it->scanned) {
            // paid on scan, just remove
            remove = true;
          }
        }
      }

      if (remove) it = missions.erase(it);
      else ++it;
    }

    // Clamp stats
    hullMax = std::max(1.0, hullMax);
    hull = std::clamp(hull, 0.0, hullMax);
    fuelMax = std::max(0.0, fuelMax);
    fuel = std::clamp(fuel, 0.0, fuelMax);
    cargoCapacity = std::max(0.0, cargoCapacity);

    // ---
    // Rendering
    // ---

    int displayW = 1280, displayH = 720;
    SDL_GL_GetDrawableSize(window, &displayW, &displayH);
    glViewport(0, 0, displayW, displayH);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.02f, 0.02f, 0.03f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Camera follow ship
    const math::Vec3d shipRender = ship.positionKm() * (1.0 / kRENDER_UNIT_KM);
    const math::Vec3d camPos = shipRender + ship.forward() * (-0.03) + ship.up() * 0.015;
    cam.setPosition(camPos);
    cam.setTarget(shipRender + ship.forward() * 0.02);

    const math::Mat4d view = cam.view();
    const math::Mat4d proj = math::Mat4d::perspective(60.0 * math::degToRad(1.0),
                                                      (double)displayW / (double)displayH,
                                                      0.00001,
                                                      200.0);
    float vf[16], pf[16];
    matToFloat(view, vf);
    matToFloat(proj, pf);
    meshRenderer.setViewProj(vf, pf);
    pointRenderer.setViewProj(vf, pf);
    lineRenderer.setViewProj(vf, pf);

    // Instances: star + planets + stations + ship
    std::vector<render::InstanceData> instances;
    instances.reserve(1 + currentSystem->planets.size() + currentSystem->stations.size() + 1);

    auto pushSphere = [&](const math::Vec3d& posKm, double radiusKm, const math::Vec3d& color, const math::Quatd& q) {
      const math::Vec3d p = posKm * (1.0 / kRENDER_UNIT_KM);
      render::InstanceData inst{};
      inst.px = (float)p.x;
      inst.py = (float)p.y;
      inst.pz = (float)p.z;
      inst.scale = (float)(radiusKm / kRENDER_UNIT_KM);
      inst.qx = (float)q.x;
      inst.qy = (float)q.y;
      inst.qz = (float)q.z;
      inst.qw = (float)q.w;
      inst.cr = (float)color.x;
      inst.cg = (float)color.y;
      inst.cb = (float)color.z;
      instances.push_back(inst);
    };

    // Star at origin
    const double starR = currentSystem->star.radiusSol * kSOLAR_RADIUS_KM;
    pushSphere({0,0,0}, starR, {1.0, 0.95, 0.8}, math::Quatd::identity());

    // Planets
    for (const auto& p : currentSystem->planets) {
      const math::Vec3d pPos = orbitPosKm(p.orbit, timeDays);
      const double r = p.radiusEarth * kEARTH_RADIUS_KM;
      math::Vec3d col{0.55,0.55,0.65};
      switch (p.type) {
        case sim::PlanetType::Rocky: col = {0.55,0.52,0.48}; break;
        case sim::PlanetType::Ocean: col = {0.35,0.45,0.75}; break;
        case sim::PlanetType::Ice: col = {0.75,0.80,0.92}; break;
        case sim::PlanetType::Desert: col = {0.75,0.65,0.35}; break;
        case sim::PlanetType::GasGiant: col = {0.65,0.55,0.75}; break;
        default: break;
      }
      pushSphere(pPos, r, col, math::Quatd::identity());
    }

    // Stations as cubes
    std::vector<render::InstanceData> stationInstances;
    stationInstances.reserve(currentSystem->stations.size());
    for (const auto& st : currentSystem->stations) {
      const math::Vec3d stPos = orbitPosKm(st.orbit, timeDays);
      const math::Vec3d p = stPos * (1.0 / kRENDER_UNIT_KM);
      render::InstanceData inst{};
      inst.px = (float)p.x;
      inst.py = (float)p.y;
      inst.pz = (float)p.z;
      inst.scale = (float)(st.radiusKm / kRENDER_UNIT_KM);
      inst.qx = 0.0f; inst.qy = 0.0f; inst.qz = 0.0f; inst.qw = 1.0f;
      // color by station type
      math::Vec3d col{0.75,0.75,0.75};
      if (st.type == econ::StationType::Mining) col = {0.70,0.68,0.55};
      if (st.type == econ::StationType::Refinery) col = {0.70,0.55,0.55};
      if (st.type == econ::StationType::TradeHub) col = {0.55,0.70,0.70};
      if (st.type == econ::StationType::Shipyard) col = {0.70,0.60,0.80};
      inst.cr = (float)col.x;
      inst.cg = (float)col.y;
      inst.cb = (float)col.z;
      stationInstances.push_back(inst);
    }

    // Ship as cube (tiny)
    render::InstanceData shipInst{};
    shipInst.px = (float)shipRender.x;
    shipInst.py = (float)shipRender.y;
    shipInst.pz = (float)shipRender.z;
    shipInst.scale = 3e-6f; // 3 km-ish visual scale
    shipInst.qx = (float)ship.orientation().x;
    shipInst.qy = (float)ship.orientation().y;
    shipInst.qz = (float)ship.orientation().z;
    shipInst.qw = (float)ship.orientation().w;
    shipInst.cr = 0.9f;
    shipInst.cg = 0.9f;
    shipInst.cb = 1.0f;

    // Draw spheres
    meshRenderer.setMesh(&sphere);
    meshRenderer.drawInstances(instances);

    // Draw stations + ship
    meshRenderer.setMesh(&cube);
    // append ship
    stationInstances.push_back(shipInst);
    meshRenderer.drawInstances(stationInstances);

    // Draw approach corridors (lines)
    std::vector<render::LineVertex> corridorLines;
    corridorLines.reserve(currentSystem->stations.size() * 2);
    for (const auto& st : currentSystem->stations) {
      const math::Vec3d stPos = orbitPosKm(st.orbit, timeDays);
      const math::Vec3d axis = (stPos.length() > 1e-6) ? stPos.normalized() : math::Vec3d{0,0,1};
      const math::Vec3d entry = stPos + axis * st.corridorLengthKm;
      const math::Vec3d a = stPos * (1.0 / kRENDER_UNIT_KM);
      const math::Vec3d b = entry * (1.0 / kRENDER_UNIT_KM);

      render::LineVertex v0{(float)a.x,(float)a.y,(float)a.z, 0.25f,0.9f,0.25f};
      render::LineVertex v1{(float)b.x,(float)b.y,(float)b.z, 0.25f,0.9f,0.25f};
      corridorLines.push_back(v0);
      corridorLines.push_back(v1);
    }
    lineRenderer.drawLines(corridorLines);

    // ---
    // ImGui frame
    // ---
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplSDL2_NewFrame();
    ImGui::NewFrame();

    if (showHelp) drawHelpOverlay(&showHelp);

    // HUD window
    if (showHUD) {
      ImGui::Begin("Flight HUD", &showHUD);
      ImGui::Text("System: %s (Star %s)", currentSystem->stub.name.c_str(), starClassName(currentSystem->stub.primaryClass));
      ImGui::Text("Mode: %s", (flightMode == FlightMode::Normal) ? "Normal" : "Supercruise");

      const double speed = ship.velocityKmS().length();
      ImGui::Text("Speed: %.3f km/s", speed);
      ImGui::Text("Hull: %.0f / %.0f", hull, hullMax);
      ImGui::Text("Fuel: %.1f / %.1f", fuel, fuelMax);
      ImGui::Text("Credits: %.0f", credits);

      const double warpDisplay = warp;
      ImGui::Text("Time: %.2f days   Warp: x%.0f %s", timeDays, warpDisplay, paused ? "[PAUSED]" : "");

      ImGui::Separator();
      ImGui::Text("Target: %s", navTargetName(target, *currentSystem));
      if (navTargetValid(target, *currentSystem)) {
        const auto tk = getTargetKinematics(target, *currentSystem, timeDays);
        if (tk) {
          const double dist = (tk->posKm - ship.positionKm()).length();
          ImGui::Text("Range: %.1f km", dist);
        }
      }

      ImGui::Checkbox("Approach autopilot (P)", &approachAutopilot);
      ImGui::Checkbox("Supercruise assist (V)", &supercruiseAssist);
      ImGui::Checkbox("Supercruise auto-drop", &supercruiseAutoDrop);

      ImGui::Separator();
      if (target.type == NavTargetType::Station && navTargetValid(target, *currentSystem)) {
        const auto& st = currentSystem->stations[(std::size_t)target.index];
        const math::Vec3d stPos = orbitPosKm(st.orbit, timeDays);
        const double d = (ship.positionKm() - stPos).length();
        ImGui::Text("Docking: %s", isDocked ? "Docked" : "In space");
        ImGui::Text("Comms range: %.0f km   Dist: %.0f km", st.commsRangeKm, d);
        ImGui::Text("Clearance: %s", (clearance.granted && clearance.stationId == st.id) ? "GRANTED" : "None");
        if (clearance.granted && clearance.stationId == st.id) {
          ImGui::SameLine();
          ImGui::TextDisabled("(%.0fs)", clearance.ttlSec);
        }
      }

      ImGui::End();
    }

    // Dock / Market / Services
    if (showDock) {
      ImGui::Begin("Dock / Market / Services", &showDock);
      ImGui::Text("System: %s", currentSystem->stub.name.c_str());

      if (!currentSystem->stations.empty()) {
        std::vector<const char*> names;
        names.reserve(currentSystem->stations.size());
        for (const auto& st : currentSystem->stations) names.push_back(st.name.c_str());
        ImGui::Combo("Station", &selectedStationIndex, names.data(), (int)names.size());

        const auto& st = currentSystem->stations[(std::size_t)selectedStationIndex];
        ImGui::SameLine();
        ImGui::TextDisabled("(%s, fee %.1f%%)", stationTypeName(st.type), st.feeRate * 100.0);

        const bool dockedHere = isDocked && dockedStationId == st.id;
        ImGui::Text("Status: %s", dockedHere ? "Docked" : "Not docked");

        // Services
        ImGui::SeparatorText("Services");
        if (!dockedHere) {
          ImGui::TextDisabled("Dock at this station to access services.");
        } else {
          // Repairs
          const double dmg = hullMax - hull;
          const double repairCost = dmg * 12.0 * (1.0 + st.feeRate);
          ImGui::Text("Repair: %.0f%% damage", dmg);
          ImGui::SameLine();
          ImGui::TextDisabled("(cost %.0f cr)", repairCost);
          if (ImGui::Button("Repair to 100%")) {
            if (dmg <= 0.0) {
              pushToast(toasts, "Hull already at 100%", 1.5f);
            } else if (credits < repairCost) {
              pushToast(toasts, "Insufficient credits for repairs.", 2.0f);
            } else {
              credits -= repairCost;
              hull = hullMax;
              pushToast(toasts, "Repaired hull.", 2.0f);
            }
          }

          // Refuel
          const double needFuel = fuelMax - fuel;
          if (needFuel > 0.0) {
            auto& stState = universe.stationEconomy(st, timeDays);
            const auto q = econ::quote(stState, st.economyModel, econ::CommodityId::Fuel);
            const double unitCost = q.ask * (1.0 + st.feeRate);
            const double fillCost = needFuel * unitCost;
            ImGui::Text("Refuel: need %.1f units", needFuel);
            ImGui::SameLine();
            ImGui::TextDisabled("(est cost %.0f cr)", fillCost);
            if (ImGui::Button("Refuel to full")) {
              double buyUnits = needFuel;
              // Buy from station economy, but put into fuel tank instead of cargo.
              auto res = econ::buy(stState, st.economyModel, econ::CommodityId::Fuel, buyUnits, credits, 0.10, st.feeRate);
              if (!res.ok) {
                pushToast(toasts, std::string("Refuel failed: ") + (res.reason ? res.reason : ""), 2.5f);
              } else {
                fuel = std::min(fuelMax, fuel + res.unitsDelta);
                pushToast(toasts, "Refueled.", 2.0f);
              }
            }
          } else {
            ImGui::Text("Refuel: tank full");
          }

          // Shipyard upgrades
          if (st.type == econ::StationType::Shipyard) {
            ImGui::SeparatorText("Upgrades (Shipyard)");
            const double cargoUpCost = 2500.0 * (1.0 + st.feeRate);
            if (ImGui::Button("Upgrade cargo bay (+20 units)")) {
              if (credits < cargoUpCost) pushToast(toasts, "Not enough credits.", 2.0f);
              else {
                credits -= cargoUpCost;
                cargoCapacity += 20.0;
                pushToast(toasts, "Cargo bay upgraded.", 2.0f);
              }
            }
            ImGui::SameLine();
            ImGui::TextDisabled("(%.0f cr)", cargoUpCost);

            const double fuelUpCost = 2200.0 * (1.0 + st.feeRate);
            if (ImGui::Button("Upgrade fuel tank (+5 units)")) {
              if (credits < fuelUpCost) pushToast(toasts, "Not enough credits.", 2.0f);
              else {
                credits -= fuelUpCost;
                fuelMax += 5.0;
                pushToast(toasts, "Fuel tank upgraded.", 2.0f);
              }
            }
            ImGui::SameLine();
            ImGui::TextDisabled("(%.0f cr)", fuelUpCost);
          }
        }

        // Market (gated by docking)
        ImGui::SeparatorText("Market");
        auto& stState = universe.stationEconomy(st, timeDays);

        // Cargo capacity
        double cargoUsed = 0.0;
        for (double u : cargo) cargoUsed += u;
        ImGui::Text("Cargo: %.0f / %.0f units", cargoUsed, cargoCapacity);

        static int commodityIdx = 0;
        std::vector<const char*> cnames;
        cnames.reserve(econ::kCommodityCount);
        for (std::size_t i = 0; i < econ::kCommodityCount; ++i) cnames.push_back(econ::commodityTable()[i].name);
        ImGui::Combo("Commodity", &commodityIdx, cnames.data(), (int)cnames.size());

        const auto cid = static_cast<econ::CommodityId>(commodityIdx);
        const auto q = econ::quote(stState, st.economyModel, cid);
        ImGui::Text("Ask: %.1f  Bid: %.1f  Inv: %.0f", q.ask, q.bid, q.inventory);
        ImGui::Text("You: %.0f units", cargo[(std::size_t)cid]);

        static float units = 1.0f;
        ImGui::SliderFloat("Units", &units, 1.0f, 25.0f, "%.0f");

        const bool dockedHere = isDocked && dockedStationId == st.id;
        if (!dockedHere) {
          ImGui::TextDisabled("Dock to trade.");
        } else {
          if (ImGui::Button("Buy")) {
            if (cargoUsed + units > cargoCapacity) {
              pushToast(toasts, "Not enough cargo space.", 2.0f);
            } else {
              auto res = econ::buy(stState, st.economyModel, cid, units, credits, 0.10, st.feeRate);
              if (res.ok) {
                cargo[(std::size_t)cid] += res.unitsDelta;
                pushToast(toasts, "Bought cargo.", 2.0f);
              } else {
                pushToast(toasts, std::string("Buy failed: ") + (res.reason ? res.reason : ""), 2.5f);
              }
            }
          }
          ImGui::SameLine();
          if (ImGui::Button("Sell")) {
            const double have = cargo[(std::size_t)cid];
            if (have < units) {
              pushToast(toasts, "Not enough units to sell.", 2.0f);
            } else {
              auto res = econ::sell(stState, st.economyModel, cid, units, credits, 0.10, st.feeRate);
              if (res.ok) {
                cargo[(std::size_t)cid] += res.unitsDelta; // negative
                pushToast(toasts, "Sold cargo.", 2.0f);
              } else {
                pushToast(toasts, std::string("Sell failed: ") + (res.reason ? res.reason : ""), 2.5f);
              }
            }
          }
        }
      }
      ImGui::End();
    }

    // Galaxy / FSD jump
    if (showGalaxy) {
      ImGui::Begin("Galaxy / FSD", &showGalaxy);
      ImGui::Text("Current system: %s", currentSystem->stub.name.c_str());
      ImGui::Text("FSD Fuel: %.1f/%.1f   Cooldown: %.0fs", fuel, fuelMax, fsdCooldownSec);
      ImGui::Text("Note: hyperspace jumps require docking.");

      ImGui::InputDouble("Query X (ly)", &galQueryPos.x);
      ImGui::InputDouble("Query Y (ly)", &galQueryPos.y);
      ImGui::InputDouble("Query Z (ly)", &galQueryPos.z);
      ImGui::SliderFloat("Radius (ly)", &galRadius, 10.0f, 200.0f);
      ImGui::SliderInt("Max systems", &galMax, 8, 80);

      if (ImGui::Button("Query")) {
        galList = universe.queryNearby(galQueryPos, galRadius, galMax);
        selectedSystemIndex = 0;
        selectedDestStationIndex = 0;
      }

      if (!galList.empty()) {
        std::vector<const char*> names;
        names.reserve(galList.size());
        for (const auto& s : galList) names.push_back(s.name.c_str());
        ImGui::Combo("Destination system", &selectedSystemIndex, names.data(), (int)names.size());

        const auto& destStub = galList[(std::size_t)selectedSystemIndex];
        const double distLy = (destStub.posLy - currentSystem->stub.posLy).length();
        ImGui::Text("Distance: %.2f ly", distLy);

        const auto& destSys = universe.getSystem(destStub.id, &destStub);
        if (!destSys.stations.empty()) {
          std::vector<const char*> stNames;
          stNames.reserve(destSys.stations.size());
          for (const auto& st : destSys.stations) stNames.push_back(st.name.c_str());
          ImGui::Combo("Arrival station", &selectedDestStationIndex, stNames.data(), (int)stNames.size());
        } else {
          ImGui::TextDisabled("Destination has no stations.");
          selectedDestStationIndex = 0;
        }

        // Simple FSD model
        const double fuelPerLy = 0.35;
        const double minFuel = 1.0;
        const double fuelCost = std::max(minFuel, distLy * fuelPerLy);
        const double maxJumpLy = 70.0;
        const bool inRange = distLy <= maxJumpLy;

        ImGui::Text("Fuel cost: %.1f (%.2f/unit ly), Max range: %.0f ly", fuelCost, fuelPerLy, maxJumpLy);

        bool dockedSomewhere = isDocked;
        if (!dockedSomewhere) {
          ImGui::TextDisabled("Dock to engage hyperspace.");
        }
        if (fsdCooldownSec > 0.0) {
          ImGui::TextDisabled("FSD cooling down.");
        }
        if (!inRange) {
          ImGui::TextDisabled("Out of range.");
        }
        if (fuel < fuelCost) {
          ImGui::TextDisabled("Not enough fuel.");
        }

        if (ImGui::Button("Engage hyperspace jump")) {
          if (!isDocked) {
            pushToast(toasts, "You must be docked to engage hyperspace.", 2.5f);
          } else if (fsdCooldownSec > 0.0) {
            pushToast(toasts, "FSD cooling down.", 2.0f);
          } else if (!inRange) {
            pushToast(toasts, "Destination out of range.", 2.0f);
          } else if (fuel < fuelCost) {
            pushToast(toasts, "Insufficient fuel.", 2.0f);
          } else {
            fuel -= fuelCost;
            fsdCooldownSec = 18.0;

            // Travel time
            const double travelDays = std::max(0.002, distLy * 0.0015);
            timeDays += travelDays;

            // Switch system
            currentSystemId = destStub.id;
            currentSystem = &universe.getSystem(currentSystemId, &destStub);
            isDocked = false;
            dockedStationId = 0;
            approachAutopilot = false;
            flightMode = FlightMode::Normal;

            // Spawn near selected station
            math::Vec3d spawnPos{0,0,0};
            math::Vec3d spawnVel{0,0,0};
            math::Vec3d axis{0,0,1};
            if (!currentSystem->stations.empty()) {
              const int si = std::clamp(selectedDestStationIndex, 0, (int)currentSystem->stations.size() - 1);
              const auto& st = currentSystem->stations[(std::size_t)si];
              const math::Vec3d stPos = orbitPosKm(st.orbit, timeDays);
              const math::Vec3d stVel = orbitVelKmS(st.orbit, timeDays);
              axis = (stPos.length() > 1e-6) ? stPos.normalized() : math::Vec3d{0,0,1};
              spawnPos = stPos + axis * (st.corridorLengthKm + 80.0);
              spawnVel = stVel;

              // set nav target to arrival station
              target.type = NavTargetType::Station;
              target.index = si;
            }
            ship.setPositionKm(spawnPos);
            ship.setVelocityKmS(spawnVel);
            ship.setOrientation(quatLookRotation(-axis, {0,1,0}));
            ship.setAngularVelocityRadS({0,0,0});

            pushToast(toasts, "Hyperspace jump complete.", 3.0f);
          }
        }
      }
      ImGui::End();
    }

    // Missions
    if (showMissions) {
      ImGui::Begin("Missions", &showMissions);
      ImGui::Text("Active missions: %d", (int)missions.size());
      if (missions.empty()) {
        ImGui::TextDisabled("No active missions.");
      }
      for (int i = 0; i < (int)missions.size(); ++i) {
        auto& m = missions[(std::size_t)i];
        const auto& ds = universe.getSystem(m.destSystem);
        ImGui::Separator();
        ImGui::Text("%s", missionTypeName(m.type));
        ImGui::Text("Destination: %s", ds.stub.name.c_str());
        if (m.type == sim::MissionType::Delivery) {
          ImGui::Text("Cargo: %.0f %s", m.units, econ::commodityName(m.commodity).data());
        }
        ImGui::Text("Reward: %.0f cr", m.rewardCredits);
        if (m.expiryDay > 0.0) ImGui::Text("Expires day: %.0f", m.expiryDay);

        if (ImGui::Button(("Set nav target##" + std::to_string(i)).c_str())) {
          if (m.destSystem == currentSystemId) {
            // set target to station if in this system
            for (int si = 0; si < (int)currentSystem->stations.size(); ++si) {
              if (currentSystem->stations[(std::size_t)si].id == m.destStation) {
                target.type = NavTargetType::Station;
                target.index = si;
                break;
              }
            }
          } else {
            showGalaxy = true;
            pushToast(toasts, "Open Galaxy/FSD to jump closer.", 2.0f);
          }
        }
      }

      // Mission board when docked
      if (isDocked && !currentSystem->stations.empty()) {
        const sim::Station* dockSt = nullptr;
        for (const auto& st : currentSystem->stations) if (st.id == dockedStationId) { dockSt = &st; break; }
        if (dockSt) {
          ImGui::SeparatorText("Mission board");
          const auto offers = generateMissionOffers(universe, *currentSystem, *dockSt, timeDays);
          int shown = 0;
          for (const auto& m : offers) {
            // Don't show offers that we already accepted (same id)
            bool already = false;
            for (const auto& a : missions) if (a.id == m.id) { already = true; break; }
            if (already) continue;
            ImGui::Separator();
            ImGui::Text("%s", missionTypeName(m.type));
            const auto& ds = universe.getSystem(m.destSystem);
            ImGui::Text("To: %s", ds.stub.name.c_str());
            if (m.type == sim::MissionType::Delivery) {
              ImGui::Text("Deliver: %.0f %s", m.units, econ::commodityName(m.commodity).data());
            }
            ImGui::Text("Reward: %.0f cr", m.rewardCredits);
            if (m.expiryDay > 0.0) ImGui::Text("Expires day: %.0f", m.expiryDay);

            if (ImGui::Button(("Accept##" + std::to_string((long long)m.id)).c_str())) {
              if ((int)missions.size() >= 8) {
                pushToast(toasts, "Mission log full.", 2.0f);
              } else {
                missions.push_back(m);
                pushToast(toasts, "Mission accepted.", 2.0f);
              }
            }
            shown++;
            if (shown >= 6) break;
          }
        }
      } else {
        ImGui::SeparatorText("Mission board");
        ImGui::TextDisabled("Dock at a station to accept missions.");
      }

      ImGui::End();
    }

    // ---
    // Overlay: target marker + toasts
    // ---
    {
      ImDrawList* dl = ImGui::GetForegroundDrawList();
      if (navTargetValid(target, *currentSystem)) {
        const auto tk = getTargetKinematics(target, *currentSystem, timeDays);
        if (tk) {
          ImVec2 sp;
          const math::Vec3d tpos = tk->posKm * (1.0 / kRENDER_UNIT_KM);
          if (projectToScreen(tpos, view, proj, displayW, displayH, sp)) {
            dl->AddCircle(sp, 10.0f, IM_COL32(0, 255, 0, 220), 24, 2.0f);
            dl->AddText(ImVec2(sp.x + 12, sp.y - 10), IM_COL32(220, 255, 220, 230), navTargetName(target, *currentSystem));
          }
        }
      }

      // Toast stack (top-left)
      float y = 10.0f;
      for (const auto& t : toasts) {
        dl->AddText(ImVec2(10.0f, y), IM_COL32(255, 255, 255, 220), t.text.c_str());
        y += 18.0f;
      }
    }

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    SDL_GL_SwapWindow(window);
  }

  // Cleanup
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplSDL2_Shutdown();
  ImGui::DestroyContext();

  SDL_GL_DeleteContext(glctx);
  SDL_DestroyWindow(window);
  SDL_Quit();

  return 0;
}
