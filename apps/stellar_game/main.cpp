#include "stellar/core/Log.h"
#include "stellar/core/SplitMix64.h"
#include "stellar/econ/Market.h"
#include "stellar/econ/RoutePlanner.h"
#include "stellar/math/Math.h"
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

#include <imgui.h>
#include <backends/imgui_impl_opengl3.h>
#include <backends/imgui_impl_sdl2.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

using namespace stellar;

static constexpr double kAU_KM = 149597870.7;
static constexpr double kSOLAR_RADIUS_KM = 695700.0;
static constexpr double kEARTH_RADIUS_KM = 6371.0;
static constexpr double kRENDER_UNIT_KM = 1.0e6; // 1 unit = 1 million km

static void matToFloat(const math::Mat4d& m, float out[16]) {
  for (int i = 0; i < 16; ++i) out[i] = static_cast<float>(m.m[i]);
}

static double clampd(double v, double lo, double hi) {
  return std::max(lo, std::min(hi, v));
}

static math::Vec3d clampLen(const math::Vec3d& v, double maxLen) {
  const double len = v.length();
  if (len <= maxLen || len <= 1e-12) return v;
  return v * (maxLen / len);
}

static math::Vec3d safeNormalize(const math::Vec3d& v, const math::Vec3d& fallback = {0, 0, 1}) {
  const double len = v.length();
  if (len <= 1e-12) return fallback;
  return v * (1.0 / len);
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

// Orbit state (pos + vel) in kilometers / km/s.
static void orbitStateKm(const sim::OrbitElements& el,
                         double timeDays,
                         math::Vec3d* outPosKm,
                         math::Vec3d* outVelKmS) {
  const math::Vec3d p0Km = sim::orbitPosition3DAU(el, timeDays) * kAU_KM;
  if (outPosKm) *outPosKm = p0Km;

  if (outVelKmS) {
    // Finite difference over a small interval (seconds) to estimate orbital velocity.
    const double dtSec = 60.0;
    const double dtDays = dtSec / 86400.0;
    const math::Vec3d p1Km = sim::orbitPosition3DAU(el, timeDays + dtDays) * kAU_KM;
    *outVelKmS = (p1Km - p0Km) * (1.0 / dtSec);
  }
}

enum class TargetKind : int {
  None = 0,
  Station,
  Planet,
};

struct NavTarget {
  TargetKind kind{TargetKind::None};
  std::size_t index{0};

  bool isValid(const sim::StarSystem& sys) const {
    switch (kind) {
      case TargetKind::Station: return index < sys.stations.size();
      case TargetKind::Planet: return index < sys.planets.size();
      default: return false;
    }
  }

  void clear() {
    kind = TargetKind::None;
    index = 0;
  }
};

struct BodyState {
  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};
  double radiusKm{0.0};
  const char* kindLabel{""};
  std::string name;
};

static bool getTargetState(const sim::StarSystem& sys,
                           const NavTarget& t,
                           double timeDays,
                           BodyState* out) {
  if (!out) return false;

  if (!t.isValid(sys)) return false;

  if (t.kind == TargetKind::Station) {
    const auto& st = sys.stations[t.index];
    orbitStateKm(st.orbit, timeDays, &out->posKm, &out->velKmS);
    out->radiusKm = st.radiusKm;
    out->kindLabel = "Station";
    out->name = st.name;
    return true;
  }
  if (t.kind == TargetKind::Planet) {
    const auto& p = sys.planets[t.index];
    orbitStateKm(p.orbit, timeDays, &out->posKm, &out->velKmS);
    out->radiusKm = p.radiusEarth * kEARTH_RADIUS_KM;
    out->kindLabel = "Planet";
    out->name = p.name;
    return true;
  }

  return false;
}

struct CorridorInfo {
  bool inside{false};
  double axialKm{0.0};
  double lateralKm{0.0};
  double cosToAxis{0.0};
  math::Vec3d axisOut{0, 0, 1};
};

static CorridorInfo stationCorridorInfo(const sim::Station& st,
                                        const math::Vec3d& shipPosKm,
                                        const math::Vec3d& stPosKm) {
  CorridorInfo info{};
  info.axisOut = safeNormalize(stPosKm, {0, 0, 1});

  const math::Vec3d rel = shipPosKm - stPosKm;
  const double dist = rel.length();

  if (dist <= 1e-9) {
    info.inside = true;
    return info;
  }

  const math::Vec3d relN = rel * (1.0 / dist);
  info.cosToAxis = math::dot(relN, info.axisOut);
  info.axialKm = math::dot(rel, info.axisOut);
  const math::Vec3d relPerp = rel - info.axisOut * info.axialKm;
  info.lateralKm = relPerp.length();

  const double cosHalf = std::cos(st.dockingCorridorHalfAngleRad);
  const bool correctSide = (info.axialKm > 0.0);
  const bool withinCone = (info.cosToAxis >= cosHalf);
  const bool withinLen = (info.axialKm <= st.dockingCorridorLengthKm);
  info.inside = correctSide && withinCone && withinLen;
  return info;
}

struct DockCheck {
  bool ok{false};
  bool inCorridor{false};
  double distKm{0.0};
  double relSpeedKmS{0.0};
  double cosAlign{0.0};
};

static DockCheck canDock(const sim::Ship& ship,
                         const sim::Station& st,
                         const math::Vec3d& stPosKm,
                         const math::Vec3d& stVelKmS) {
  DockCheck c{};
  const math::Vec3d rel = ship.positionKm() - stPosKm;
  c.distKm = rel.length();

  const CorridorInfo corridor = stationCorridorInfo(st, ship.positionKm(), stPosKm);
  c.inCorridor = corridor.inside;

  const math::Vec3d vRel = ship.velocityKmS() - stVelKmS;
  c.relSpeedKmS = vRel.length();

  const math::Vec3d dirToStation = safeNormalize(stPosKm - ship.positionKm(), ship.forward());
  c.cosAlign = math::dot(ship.forward(), dirToStation);

  const double cosAlignReq = std::cos(stellar::math::degToRad(12.0));

  c.ok = (c.distKm <= st.dockingRangeKm) &&
         c.inCorridor &&
         (c.relSpeedKmS <= st.dockingSpeedLimitKmS) &&
         (c.cosAlign >= cosAlignReq);
  return c;
}

// Project a world-space position (in render units) to screen coordinates.
// Returns true if the point is in front of the camera (w>0). The marker may still be offscreen.
static bool projectToScreen(const math::Vec3d& posU,
                            const math::Mat4d& view,
                            const math::Mat4d& proj,
                            const ImVec2& vpPos,
                            const ImVec2& vpSize,
                            ImVec2* outPx,
                            bool* outOnScreen = nullptr,
                            ImVec2* outNdc = nullptr) {
  const math::Mat4d vp = proj * view;

  // Column-major matrix * column vec4
  const double x = vp.m[0] * posU.x + vp.m[4] * posU.y + vp.m[8] * posU.z + vp.m[12] * 1.0;
  const double y = vp.m[1] * posU.x + vp.m[5] * posU.y + vp.m[9] * posU.z + vp.m[13] * 1.0;
  const double z = vp.m[2] * posU.x + vp.m[6] * posU.y + vp.m[10] * posU.z + vp.m[14] * 1.0;
  const double w = vp.m[3] * posU.x + vp.m[7] * posU.y + vp.m[11] * posU.z + vp.m[15] * 1.0;

  (void)z;
  if (w <= 1e-9) {
    if (outOnScreen) *outOnScreen = false;
    return false;
  }

  const double ndcX = x / w;
  const double ndcY = y / w;

  if (outNdc) *outNdc = ImVec2((float)ndcX, (float)ndcY);

  const bool on = (ndcX >= -1.0 && ndcX <= 1.0 && ndcY >= -1.0 && ndcY <= 1.0);
  if (outOnScreen) *outOnScreen = on;

  // Map NDC -> pixel
  const float px = vpPos.x + (float)((ndcX * 0.5 + 0.5) * vpSize.x);
  const float py = vpPos.y + (float)((1.0 - (ndcY * 0.5 + 0.5)) * vpSize.y);
  if (outPx) *outPx = ImVec2(px, py);
  return true;
}

struct Toast {
  std::string msg;
  double timer{0.0};

  void show(std::string m, double seconds = 2.5) {
    msg = std::move(m);
    timer = seconds;
  }
};

// Pioneer-style time compression levels (sim seconds per real second).
static constexpr std::array<double, 8> kTimeLevels = {
    1.0,   // 1x
    5.0,   // 5x
    10.0,  // 10x
    50.0,  // 50x
    100.0, // 100x
    500.0, // 500x
    1000.0, // 1000x
    5000.0, // 5000x
};

static int maxTimeLevelIndex() {
  return (int)kTimeLevels.size() - 1;
}

static double timeScaleForIndex(int idx) {
  idx = std::max(0, std::min(idx, maxTimeLevelIndex()));
  return kTimeLevels[(std::size_t)idx];
}

static int allowedTimeLevelIndex(const sim::StarSystem& sys,
                                 const sim::Ship& ship,
                                 double timeDays,
                                 bool docked,
                                 bool supercruise,
                                 bool autopilot,
                                 bool manualInputActive) {
  if (supercruise) return 0; // keep supercruise stable / deterministic
  if (docked) return maxTimeLevelIndex();

  const math::Vec3d shipPos = ship.positionKm();
  const double shipSpeed = ship.velocityKmS().length();

  // Find nearest distance-to-surface among star/planets/stations.
  double nearestSurfaceKm = 1.0e30;

  // Star
  {
    const double rStarKm = sys.star.radiusSol * kSOLAR_RADIUS_KM;
    const double d = shipPos.length() - rStarKm;
    nearestSurfaceKm = std::min(nearestSurfaceKm, d);
  }

  // Planets
  for (const auto& p : sys.planets) {
    math::Vec3d pPosKm{};
    orbitStateKm(p.orbit, timeDays, &pPosKm, nullptr);
    const double rKm = p.radiusEarth * kEARTH_RADIUS_KM;
    const double d = (shipPos - pPosKm).length() - rKm;
    nearestSurfaceKm = std::min(nearestSurfaceKm, d);
  }

  // Stations
  for (const auto& st : sys.stations) {
    math::Vec3d stPosKm{};
    orbitStateKm(st.orbit, timeDays, &stPosKm, nullptr);
    const double d = (shipPos - stPosKm).length() - st.radiusKm;
    nearestSurfaceKm = std::min(nearestSurfaceKm, d);
  }

  nearestSurfaceKm = std::max(0.0, nearestSurfaceKm);

  int maxIdx = maxTimeLevelIndex();

  // Safety clamp by proximity.
  if (nearestSurfaceKm < 200.0) {
    maxIdx = 0;
  } else if (nearestSurfaceKm < 1500.0) {
    maxIdx = std::min(maxIdx, 1); // 5x
  } else if (nearestSurfaceKm < 8000.0) {
    maxIdx = std::min(maxIdx, 2); // 10x
  } else if (nearestSurfaceKm < 50000.0) {
    maxIdx = std::min(maxIdx, 3); // 50x
  } else if (nearestSurfaceKm < 250000.0) {
    maxIdx = std::min(maxIdx, 4); // 100x
  }

  // Safety clamp by speed.
  if (shipSpeed > 5.0) {
    maxIdx = std::min(maxIdx, 1);
  } else if (shipSpeed > 2.0) {
    maxIdx = std::min(maxIdx, 2);
  } else if (shipSpeed > 0.8) {
    maxIdx = std::min(maxIdx, 3);
  }

  // Manual control should keep time compression modest.
  if (manualInputActive && !autopilot) {
    maxIdx = std::min(maxIdx, 2); // 10x max while actively flying
  }

  return std::max(0, maxIdx);
}

static double supercruiseMaxSpeedKmS(const sim::StarSystem& sys,
                                     const math::Vec3d& shipPosKm,
                                     double timeDays) {
  // Very simple "gravity well" limiter: max speed is proportional to the
  // distance to the nearest body's surface, clamped.
  double nearestSurfaceKm = 1.0e30;

  // Star
  {
    const double rStarKm = sys.star.radiusSol * kSOLAR_RADIUS_KM;
    nearestSurfaceKm = std::min(nearestSurfaceKm, shipPosKm.length() - rStarKm);
  }

  // Planets
  for (const auto& p : sys.planets) {
    math::Vec3d pPosKm{};
    orbitStateKm(p.orbit, timeDays, &pPosKm, nullptr);
    const double rKm = p.radiusEarth * kEARTH_RADIUS_KM;
    nearestSurfaceKm = std::min(nearestSurfaceKm, (shipPosKm - pPosKm).length() - rKm);
  }

  // Stations
  for (const auto& st : sys.stations) {
    math::Vec3d stPosKm{};
    orbitStateKm(st.orbit, timeDays, &stPosKm, nullptr);
    nearestSurfaceKm = std::min(nearestSurfaceKm, (shipPosKm - stPosKm).length() - st.radiusKm);
  }

  nearestSurfaceKm = std::max(0.0, nearestSurfaceKm);

  const double minKmS = 5.0;
  const double maxKmS = 20000.0;
  const double scaled = nearestSurfaceKm * 0.01; // 100,000 km -> 1000 km/s
  return clampd(minKmS + scaled, minKmS, maxKmS);
}

static bool canEnterSupercruise(const sim::StarSystem& sys,
                                const math::Vec3d& shipPosKm,
                                double timeDays,
                                std::string* outReason) {
  // "Mass lock": disallow if too close to a body.
  double nearestSurfaceKm = 1.0e30;

  // Star
  {
    const double rStarKm = sys.star.radiusSol * kSOLAR_RADIUS_KM;
    nearestSurfaceKm = std::min(nearestSurfaceKm, shipPosKm.length() - rStarKm);
  }
  // Planets
  for (const auto& p : sys.planets) {
    math::Vec3d pPosKm{};
    orbitStateKm(p.orbit, timeDays, &pPosKm, nullptr);
    const double rKm = p.radiusEarth * kEARTH_RADIUS_KM;
    nearestSurfaceKm = std::min(nearestSurfaceKm, (shipPosKm - pPosKm).length() - rKm);
  }
  // Stations
  for (const auto& st : sys.stations) {
    math::Vec3d stPosKm{};
    orbitStateKm(st.orbit, timeDays, &stPosKm, nullptr);
    nearestSurfaceKm = std::min(nearestSurfaceKm, (shipPosKm - stPosKm).length() - st.radiusKm);
  }

  nearestSurfaceKm = std::max(0.0, nearestSurfaceKm);

  if (nearestSurfaceKm < 2500.0) {
    if (outReason) *outReason = "Too close to a body (mass lock).";
    return false;
  }
  return true;
}

static bool beginStationSelectHUD(const sim::StarSystem& sys, int& stationIndex, const char* title) {
  bool changed = false;

  ImGui::Begin(title);

  ImGui::Text("System: %s  (Star %s, planets %d, stations %d)",
              sys.stub.name.c_str(),
              starClassName(sys.stub.primaryClass),
              sys.stub.planetCount,
              sys.stub.stationCount);

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
  bool paused = false;

  int requestedTimeIdx = 2; // 10x default (feels less static)
  int appliedTimeIdx = requestedTimeIdx;

  // Docking state
  bool docked = false;
  sim::StationId dockedStationId = 0;
  int dockedStationIndex = 0; // also used by market selection

  // Navigation / autopilot
  NavTarget navTarget{};
  bool autopilot = false;

  // Supercruise
  bool supercruise = false;
  double scThrottle = 0.0;        // 0..1
  double scSpeedKmS = 0.0;        // actual
  bool scAutoDrop = true;

  // Player economy
  double credits = 2500.0;
  std::array<double, econ::kCommodityCount> cargo{};

  // Save/load
  const std::string savePath = "savegame.txt";

  // UI state
  bool showGalaxy = true;
  bool showShip = true;
  bool showEconomy = true;
  bool showNav = true;
  bool showHUD = true;

  // Toasts
  Toast toast;

  // Starfield (directions)
  struct StarDir {
    math::Vec3d dir;
    float size;
    float c;
  };
  std::vector<StarDir> starDirs;
  {
    core::SplitMix64 rng(0xC0FFEEULL);
    const int n = 2600;
    starDirs.reserve(n);
    for (int i = 0; i < n; ++i) {
      // Random direction on sphere
      const double z = rng.range(-1.0, 1.0);
      const double a = rng.range(0.0, 2.0 * math::kPi);
      const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
      const math::Vec3d dir{r * std::cos(a), r * std::sin(a), z};
      const float size = (float)rng.range(1.0, 2.5);
      const float c = (float)rng.range(0.75, 1.0);
      starDirs.push_back({dir, size, c});
    }
  }

  // Spawn near first station (if present) so it's immediately playable.
  {
    if (!currentSystem->stations.empty()) {
      dockedStationIndex = 0;
      const auto& st = currentSystem->stations[0];
      math::Vec3d stPosKm{}, stVelKmS{};
      orbitStateKm(st.orbit, timeDays, &stPosKm, &stVelKmS);
      const math::Vec3d axisOut = safeNormalize(stPosKm, {0, 0, 1});
      const double spawnDist = std::max(1200.0, st.dockingCorridorLengthKm * 0.85);
      ship.setPositionKm(stPosKm + axisOut * spawnDist);
      ship.setVelocityKmS(stVelKmS);
      navTarget.kind = TargetKind::Station;
      navTarget.index = 0;
      toast.show("Spawned near " + st.name + " (use T to cycle targets, P for autopilot)", 3.0);
    } else {
      ship.setPositionKm({0, 0, -8000.0});
    }
  }

  bool running = true;
  auto last = std::chrono::high_resolution_clock::now();

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
        const SDL_Keycode key = event.key.keysym.sym;

        if (key == SDLK_ESCAPE) running = false;

        if (key == SDLK_SPACE) paused = !paused;

        if (key == SDLK_F5) {
          sim::SaveGame s{};
          s.seed = universe.seed();
          s.timeDays = timeDays;
          s.currentSystem = currentSystem->stub.id;
          s.dockedStation = docked ? dockedStationId : 0;

          s.shipPosKm = ship.positionKm();
          s.shipVelKmS = ship.velocityKmS();
          s.shipOrient = ship.orientation();
          s.shipAngVelRadS = ship.angularVelocityRadS();

          s.credits = credits;
          s.cargo = cargo;
          s.stationOverrides = universe.exportStationOverrides();

          if (sim::saveToFile(s, savePath)) {
            toast.show("Saved to " + savePath);
          }
        }

        if (key == SDLK_F9) {
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

            docked = (s.dockedStation != 0);
            dockedStationId = s.dockedStation;
            dockedStationIndex = 0;
            for (std::size_t i = 0; i < sys.stations.size(); ++i) {
              if (sys.stations[i].id == dockedStationId) dockedStationIndex = (int)i;
            }

            // Reset transient flight modes
            autopilot = false;
            supercruise = false;
            scThrottle = 0.0;
            scSpeedKmS = 0.0;

            navTarget.clear();
            if (!sys.stations.empty()) {
              navTarget.kind = TargetKind::Station;
              navTarget.index = (std::size_t)dockedStationIndex;
            }

            toast.show("Loaded " + savePath);
          }
        }

        if (key == SDLK_TAB) showGalaxy = !showGalaxy;
        if (key == SDLK_F1) showShip = !showShip;
        if (key == SDLK_F2) showEconomy = !showEconomy;
        if (key == SDLK_F3) showNav = !showNav;
        if (key == SDLK_F4) showHUD = !showHUD;

        // Time compression hotkeys
        if (key == SDLK_PAGEUP) requestedTimeIdx = std::min(requestedTimeIdx + 1, maxTimeLevelIndex());
        if (key == SDLK_PAGEDOWN) requestedTimeIdx = std::max(requestedTimeIdx - 1, 0);
        if (key == SDLK_HOME) requestedTimeIdx = 0;
        if (key == SDLK_END) requestedTimeIdx = maxTimeLevelIndex();

        // Targeting
        if (key == SDLK_t) {
          // Cycle: stations then planets.
          std::vector<NavTarget> list;
          list.reserve(currentSystem->stations.size() + currentSystem->planets.size());
          for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) list.push_back({TargetKind::Station, i});
          for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) list.push_back({TargetKind::Planet, i});

          if (list.empty()) {
            navTarget.clear();
          } else {
            int idx = 0;
            for (int i = 0; i < (int)list.size(); ++i) {
              if (navTarget.kind == list[(std::size_t)i].kind && navTarget.index == list[(std::size_t)i].index) {
                idx = i;
                break;
              }
            }
            idx = (idx + 1) % (int)list.size();
            navTarget = list[(std::size_t)idx];

            BodyState bs{};
            if (getTargetState(*currentSystem, navTarget, timeDays, &bs)) {
              toast.show(std::string("Target: ") + bs.kindLabel + " - " + bs.name);
            }
          }
        }

        if (key == SDLK_y) {
          navTarget.clear();
          toast.show("Target cleared");
        }

        // Autopilot
        if (key == SDLK_p) {
          autopilot = !autopilot;
          if (autopilot) toast.show("Autopilot ON");
          else toast.show("Autopilot OFF");
        }

        // Dock / undock
        if (key == SDLK_g) {
          if (docked) {
            // Undock: keep station velocity, add a small push outward.
            if (dockedStationIndex >= 0 && (std::size_t)dockedStationIndex < currentSystem->stations.size()) {
              const auto& st = currentSystem->stations[(std::size_t)dockedStationIndex];
              math::Vec3d stPosKm{}, stVelKmS{};
              orbitStateKm(st.orbit, timeDays, &stPosKm, &stVelKmS);
              const math::Vec3d axisOut = safeNormalize(stPosKm, {0, 0, 1});
              ship.setPositionKm(stPosKm + axisOut * (st.radiusKm + 5.0));
              ship.setVelocityKmS(stVelKmS + axisOut * 0.02);
            }
            docked = false;
            dockedStationId = 0;
            toast.show("Undocked");
          } else {
            // Attempt docking at closest station.
            int best = -1;
            double bestD = 1.0e30;
            for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
              math::Vec3d stPosKm{};
              orbitStateKm(currentSystem->stations[i].orbit, timeDays, &stPosKm, nullptr);
              const double d = (ship.positionKm() - stPosKm).length();
              if (d < bestD) {
                bestD = d;
                best = (int)i;
              }
            }

            if (best >= 0) {
              const auto& st = currentSystem->stations[(std::size_t)best];
              math::Vec3d stPosKm{}, stVelKmS{};
              orbitStateKm(st.orbit, timeDays, &stPosKm, &stVelKmS);
              const DockCheck c = canDock(ship, st, stPosKm, stVelKmS);
              if (c.ok) {
                docked = true;
                dockedStationId = st.id;
                dockedStationIndex = best;
                ship.setVelocityKmS(stVelKmS);
                toast.show("Docked at " + st.name);
              } else {
                toast.show("Docking conditions not met (align, slow down, be in corridor)");
              }
            }
          }
        }

        // Supercruise
        if (key == SDLK_j) {
          if (supercruise) {
            supercruise = false;
            scSpeedKmS = 0.0;
            scThrottle = 0.0;
            toast.show("Supercruise OFF");
          } else {
            if (docked) {
              toast.show("Cannot engage supercruise while docked");
            } else {
              std::string reason;
              if (canEnterSupercruise(*currentSystem, ship.positionKm(), timeDays, &reason)) {
                supercruise = true;
                autopilot = false;
                paused = false;
                // Start with a modest speed.
                scSpeedKmS = std::max(20.0, ship.velocityKmS().length());
                scThrottle = 0.35;
                toast.show("Supercruise ON");
              } else {
                toast.show("Supercruise blocked: " + reason);
              }
            }
          }
        }
      }
    }

    // Input
    const Uint8* keys = SDL_GetKeyboardState(nullptr);
    const bool captureKeys = io.WantCaptureKeyboard;

    // Manual input intent (used both for direct flight and for canceling autopilot).
    sim::ShipInput manual{};
    bool manualActive = false;

    if (!captureKeys) {
      // Rotation always available (except while docked)
      if (!docked) {
        manual.torqueLocal.x += (keys[SDL_SCANCODE_UP] ? 1.0 : 0.0);
        manual.torqueLocal.x -= (keys[SDL_SCANCODE_DOWN] ? 1.0 : 0.0);
        manual.torqueLocal.y += (keys[SDL_SCANCODE_RIGHT] ? 1.0 : 0.0);
        manual.torqueLocal.y -= (keys[SDL_SCANCODE_LEFT] ? 1.0 : 0.0);
        manual.torqueLocal.z += (keys[SDL_SCANCODE_E] ? 1.0 : 0.0);
        manual.torqueLocal.z -= (keys[SDL_SCANCODE_Q] ? 1.0 : 0.0);
      }

      manual.boost = keys[SDL_SCANCODE_LSHIFT] != 0;
      manual.brake = keys[SDL_SCANCODE_X] != 0;

      static bool dampers = true;
      if (keys[SDL_SCANCODE_Z]) dampers = true;
      if (keys[SDL_SCANCODE_C]) dampers = false;
      manual.dampers = dampers;

      if (!supercruise && !docked) {
        manual.thrustLocal.z += (keys[SDL_SCANCODE_W] ? 1.0 : 0.0);
        manual.thrustLocal.z -= (keys[SDL_SCANCODE_S] ? 1.0 : 0.0);
        manual.thrustLocal.x += (keys[SDL_SCANCODE_D] ? 1.0 : 0.0);
        manual.thrustLocal.x -= (keys[SDL_SCANCODE_A] ? 1.0 : 0.0);
        manual.thrustLocal.y += (keys[SDL_SCANCODE_R] ? 1.0 : 0.0);
        manual.thrustLocal.y -= (keys[SDL_SCANCODE_F] ? 1.0 : 0.0);
      }

      if (supercruise && !docked) {
        const double rate = 0.5; // per second
        if (keys[SDL_SCANCODE_W]) scThrottle = clampd(scThrottle + rate * dtReal, 0.0, 1.0);
        if (keys[SDL_SCANCODE_S]) scThrottle = clampd(scThrottle - rate * dtReal, 0.0, 1.0);
        if (keys[SDL_SCANCODE_X]) scThrottle = 0.0;
      }

      manualActive = (manual.thrustLocal.lengthSq() > 1e-6) || (manual.torqueLocal.lengthSq() > 1e-6);
    }

    // Cancel autopilot on manual input.
    if (autopilot && manualActive && !captureKeys) {
      autopilot = false;
      toast.show("Autopilot canceled (manual input)");
    }

    // Time compression clamp / application
    const bool forceTime = (!captureKeys) && (keys[SDL_SCANCODE_LCTRL] || keys[SDL_SCANCODE_RCTRL]);
    const int allowedIdx = allowedTimeLevelIndex(*currentSystem, ship, timeDays, docked, supercruise, autopilot, manualActive);
    appliedTimeIdx = forceTime ? requestedTimeIdx : std::min(requestedTimeIdx, allowedIdx);
    const double timeScale = paused ? 0.0 : timeScaleForIndex(appliedTimeIdx);
    const double dtSim = dtReal * timeScale;

    // --- Simulation update ---
    if (toast.timer > 0.0) toast.timer = std::max(0.0, toast.timer - dtReal);

    if (!paused && dtSim > 0.0) {
      // Advance "world time"
      timeDays += dtSim / 86400.0;

      if (docked) {
        // Snap ship to docked station's orbit state
        if (dockedStationIndex >= 0 && (std::size_t)dockedStationIndex < currentSystem->stations.size()) {
          const auto& st = currentSystem->stations[(std::size_t)dockedStationIndex];
          math::Vec3d stPosKm{}, stVelKmS{};
          orbitStateKm(st.orbit, timeDays, &stPosKm, &stVelKmS);
          const math::Vec3d axisOut = safeNormalize(stPosKm, {0, 0, 1});
          ship.setPositionKm(stPosKm + axisOut * (st.radiusKm + 3.0));
          ship.setVelocityKmS(stVelKmS);
        }
      } else if (supercruise) {
        // Supercruise: simplified translation, keep rotation thrusters.
        const double maxKmS = supercruiseMaxSpeedKmS(*currentSystem, ship.positionKm(), timeDays);
        const double desired = scThrottle * maxKmS;
        const double accel = 0.8; // 1/s smoothing
        scSpeedKmS += (desired - scSpeedKmS) * (1.0 - std::exp(-accel * dtSim));

        // Angular integration using ship physics (but freeze linear).
        sim::ShipInput ang = manual;
        ang.thrustLocal = {0, 0, 0};
        // Don't let dampers kill our (synthetic) supercruise speed.
        ship.setVelocityKmS({0, 0, 0});
        ship.step(dtSim, ang);

        const math::Vec3d newPos = ship.positionKm() + ship.forward() * (scSpeedKmS * dtSim);
        ship.setPositionKm(newPos);
        ship.setVelocityKmS(ship.forward() * scSpeedKmS);

        // Auto-drop near target.
        if (scAutoDrop && navTarget.isValid(*currentSystem)) {
          BodyState bs{};
          if (getTargetState(*currentSystem, navTarget, timeDays, &bs)) {
            const double dist = (bs.posKm - ship.positionKm()).length();
            double dropDist = 0.0;
            if (navTarget.kind == TargetKind::Station) {
              const auto& st = currentSystem->stations[navTarget.index];
              dropDist = std::max(5000.0, st.dockingCorridorLengthKm * 0.8);
            } else {
              dropDist = std::max(20000.0, bs.radiusKm * 3.0);
            }

            if (dist <= dropDist) {
              supercruise = false;
              scSpeedKmS = 0.0;
              scThrottle = 0.0;
              // Drop with the target's orbital velocity to avoid ridiculous relative speeds.
              ship.setVelocityKmS(bs.velKmS);
              toast.show("Auto-drop near " + bs.name);
            }
          }
        }
      } else {
        // Normal flight: manual or autopilot input.
        sim::ShipInput in = manual;

        bool autoDockNow = false;

        if (autopilot && navTarget.isValid(*currentSystem)) {
          BodyState bs{};
          if (getTargetState(*currentSystem, navTarget, timeDays, &bs)) {
            // Autopilot behavior depends on target.
            if (navTarget.kind == TargetKind::Station) {
              const auto& st = currentSystem->stations[navTarget.index];
              const math::Vec3d stPos = bs.posKm;
              const math::Vec3d stVel = bs.velKmS;

              const math::Vec3d shipPos = ship.positionKm();
              const math::Vec3d shipVel = ship.velocityKmS();

              const CorridorInfo corridor = stationCorridorInfo(st, shipPos, stPos);
              const math::Vec3d axisOut = corridor.axisOut;

              const math::Vec3d rel = shipPos - stPos;
              const double dist = rel.length();
              const double axial = corridor.axialKm;
              const math::Vec3d relPerp = rel - axisOut * axial;
              const double lateral = relPerp.length();

              // Waypoint: corridor entry point (outward end)
              const math::Vec3d entryPoint = stPos + axisOut * st.dockingCorridorLengthKm;
              const math::Vec3d toEntry = entryPoint - shipPos;

              // Desired relative speed profile (stop at dockingRange)
              const double maxA = std::max(0.01, ship.maxLinearAccelKmS2());
              const double stopDist = std::max(0.0, dist - st.dockingRangeKm);
              const double vStop = std::sqrt(2.0 * maxA * stopDist);

              const double farCap = std::max(2.0, st.dockingSpeedLimitKmS * 12.0);
              const double speedCap = corridor.inside ? st.dockingSpeedLimitKmS : farCap;
              const double desiredSpeed = std::min(speedCap, vStop);

              math::Vec3d desiredRelVel{0, 0, 0};
              if (!corridor.inside || axial < 0.0) {
                // Get to corridor entry point.
                const math::Vec3d dir = safeNormalize(toEntry, -axisOut);
                const double s = std::min(farCap, toEntry.length() * 0.02); // ~2% of distance per second
                desiredRelVel = dir * s;
              } else {
                // In corridor: go inward and correct lateral offset.
                desiredRelVel = (-axisOut) * desiredSpeed;
                if (lateral > 1e-3) {
                  const math::Vec3d latDir = (-relPerp) * (1.0 / lateral);
                  const double latSpeed = std::min(desiredSpeed * 0.6 + 0.02, lateral * 0.02);
                  desiredRelVel += latDir * latSpeed;
                }
              }

              const math::Vec3d desiredVelWorld = stVel + desiredRelVel;

              // Velocity controller -> acceleration -> thrusters.
              const math::Vec3d velErr = desiredVelWorld - shipVel;
              const math::Vec3d desiredAccelWorld = clampLen(velErr * 0.7, maxA);
              const math::Vec3d accelLocal = ship.orientation().conjugate().rotate(desiredAccelWorld);
              in.thrustLocal = clampLen(accelLocal * (1.0 / maxA), 1.0);
              in.dampers = true;

              // Orientation controller: point toward station.
              const math::Vec3d desiredFwd = safeNormalize(stPos - shipPos, ship.forward());
              const math::Vec3d errAxisWorld = math::cross(ship.forward(), desiredFwd);
              const math::Vec3d errAxisLocal = ship.orientation().conjugate().rotate(errAxisWorld);
              in.torqueLocal = clampLen(errAxisLocal * 4.0, 1.0);

              // Auto-dock when conditions met.
              const DockCheck dc = canDock(ship, st, stPos, stVel);
              autoDockNow = dc.ok;
            } else {
              // Planet target: just match its orbital velocity and point at it.
              const math::Vec3d shipPos = ship.positionKm();
              const math::Vec3d shipVel = ship.velocityKmS();
              const math::Vec3d to = bs.posKm - shipPos;
              const double dist = to.length();
              const double safeDist = std::max(5000.0, bs.radiusKm * 2.0);

              const double maxA = std::max(0.01, ship.maxLinearAccelKmS2());
              const double stopDist = std::max(0.0, dist - safeDist);
              const double vStop = std::sqrt(2.0 * maxA * stopDist);
              const double desiredSpeed = std::min(5.0, vStop);
              const math::Vec3d desiredVelWorld = bs.velKmS + safeNormalize(to, ship.forward()) * desiredSpeed;

              const math::Vec3d velErr = desiredVelWorld - shipVel;
              const math::Vec3d desiredAccelWorld = clampLen(velErr * 0.5, maxA);
              const math::Vec3d accelLocal = ship.orientation().conjugate().rotate(desiredAccelWorld);
              in.thrustLocal = clampLen(accelLocal * (1.0 / maxA), 1.0);
              in.dampers = true;

              const math::Vec3d desiredFwd = safeNormalize(to, ship.forward());
              const math::Vec3d errAxisWorld = math::cross(ship.forward(), desiredFwd);
              const math::Vec3d errAxisLocal = ship.orientation().conjugate().rotate(errAxisWorld);
              in.torqueLocal = clampLen(errAxisLocal * 3.5, 1.0);
            }
          }
        }

        // Integrate ship simulation. Cap per-step to avoid huge dt when time compression is high.
        const double maxStep = 1.0; // seconds
        double remain = dtSim;
        while (remain > 0.0) {
          const double step = std::min(remain, maxStep);
          ship.step(step, in);
          remain -= step;
        }

        if (autoDockNow && navTarget.kind == TargetKind::Station && navTarget.index < currentSystem->stations.size()) {
          const auto& st = currentSystem->stations[navTarget.index];
          docked = true;
          dockedStationId = st.id;
          dockedStationIndex = (int)navTarget.index;
          toast.show("Auto-docked at " + st.name);
        }
      }
    }

    // ---- Camera (third-person follow) ----
    render::Camera cam;
    int w = 1280, h = 720;
    SDL_GetWindowSize(window, &w, &h);
    const double aspect = (h > 0) ? (double)w / (double)h : 16.0 / 9.0;

    cam.setPerspective(math::degToRad(60.0), aspect, 0.01, 20000.0);

    const math::Vec3d shipPosU = ship.positionKm() * (1.0 / kRENDER_UNIT_KM);
    const math::Vec3d back = ship.forward() * (-6.0);
    const math::Vec3d up = ship.up() * (2.0);

    cam.setPosition(shipPosU + back + up);
    cam.setOrientation(ship.orientation());

    const math::Mat4d view = cam.viewMatrix();
    const math::Mat4d proj = cam.projectionMatrix();

    float viewF[16], projF[16];
    matToFloat(view, viewF);
    matToFloat(proj, projF);

    meshRenderer.setViewProj(viewF, projF);
    lineRenderer.setViewProj(viewF, projF);
    pointRenderer.setViewProj(viewF, projF);

    // ---- Build instances (star + planets) ----
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
      math::Vec3d posKm{};
      orbitStateKm(p.orbit, timeDays, &posKm, nullptr);
      const math::Vec3d posU = posKm * (1.0 / kRENDER_UNIT_KM);

      const double radiusKm = p.radiusEarth * kEARTH_RADIUS_KM;
      const float scale = (float)std::max(0.25, (radiusKm / kRENDER_UNIT_KM) * 200.0);

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

    // ---- Orbit lines (planets + stations) ----
    std::vector<render::LineVertex> orbitLines;
    orbitLines.reserve((currentSystem->planets.size() + currentSystem->stations.size()) * 128);

    auto appendOrbit = [&](const sim::OrbitElements& el, float r, float g, float b) {
      const int seg = 96;
      math::Vec3d prev{};
      for (int s = 0; s <= seg; ++s) {
        const double t = (double)s / (double)seg * el.periodDays;
        const math::Vec3d posU = (sim::orbitPosition3DAU(el, t) * kAU_KM) * (1.0 / kRENDER_UNIT_KM);
        if (s > 0) {
          orbitLines.push_back({(float)prev.x, (float)prev.y, (float)prev.z, r, g, b});
          orbitLines.push_back({(float)posU.x, (float)posU.y, (float)posU.z, r, g, b});
        }
        prev = posU;
      }
    };

    for (const auto& p : currentSystem->planets) appendOrbit(p.orbit, 0.25f, 0.25f, 0.25f);
    for (const auto& st : currentSystem->stations) appendOrbit(st.orbit, 0.22f, 0.32f, 0.42f);

    // ---- Cubes: stations + ship ----
    std::vector<render::InstanceData> cubes;
    cubes.reserve(1 + currentSystem->stations.size());

    // Stations
    for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
      const auto& st = currentSystem->stations[i];
      math::Vec3d stPosKm{};
      orbitStateKm(st.orbit, timeDays, &stPosKm, nullptr);
      const math::Vec3d stPosU = stPosKm * (1.0 / kRENDER_UNIT_KM);
      const float scale = (float)std::max(0.12, (st.radiusKm / kRENDER_UNIT_KM) * 20000.0);

      float cr = 0.75f, cg = 0.75f, cb = 0.85f;
      if (navTarget.kind == TargetKind::Station && navTarget.index == i) {
        cr = 1.0f;
        cg = 0.55f;
        cb = 0.55f;
      }
      if (docked && (int)i == dockedStationIndex) {
        cr = 0.55f;
        cg = 1.0f;
        cb = 0.55f;
      }

      // Identity rotation (cube is symmetric anyway)
      cubes.push_back({(float)stPosU.x, (float)stPosU.y, (float)stPosU.z, scale, 1, 0, 0, 0, cr, cg, cb});

      // Corridor visualization for targeted station
      if (navTarget.kind == TargetKind::Station && navTarget.index == i) {
        const math::Vec3d axisOut = safeNormalize(stPosKm, {0, 0, 1});
        const math::Vec3d aU = axisOut * (1.0 / kRENDER_UNIT_KM);
        const math::Vec3d startU = stPosU;
        const math::Vec3d endU = stPosU + aU * st.dockingCorridorLengthKm;

        orbitLines.push_back({(float)startU.x, (float)startU.y, (float)startU.z, 0.9f, 0.35f, 0.2f});
        orbitLines.push_back({(float)endU.x, (float)endU.y, (float)endU.z, 0.9f, 0.35f, 0.2f});
      }
    }

    // Ship
    {
      const auto q = ship.orientation();
      cubes.push_back({(float)shipPosU.x,
                       (float)shipPosU.y,
                       (float)shipPosU.z,
                       0.35f,
                       (float)q.w,
                       (float)q.x,
                       (float)q.y,
                       (float)q.z,
                       0.9f,
                       0.9f,
                       1.0f});
    }

    // ---- Starfield points (attached to camera) ----
    std::vector<render::PointVertex> stars;
    stars.reserve(starDirs.size());
    {
      const math::Vec3d camPosU = cam.position();
      const double r = 12000.0;
      for (const auto& sdir : starDirs) {
        const math::Vec3d p = camPosU + sdir.dir * r;
        stars.push_back({(float)p.x, (float)p.y, (float)p.z, sdir.c, sdir.c, sdir.c, sdir.size});
      }
    }

    // ---- Render ---
    glViewport(0, 0, w, h);
    glClearColor(0.01f, 0.01f, 0.02f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Stars (points)
    pointRenderer.drawPoints(stars);

    // Orbits + corridor line (lines)
    lineRenderer.drawLines(orbitLines);

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

    if (showShip) {
      ImGui::Begin("Ship / Flight");

      const auto pos = ship.positionKm();
      const auto vel = ship.velocityKmS();
      const auto wv = ship.angularVelocityRadS();

      ImGui::Text("Time: %.2f days  %s", timeDays, paused ? "[PAUSED]" : "");
      ImGui::Text("Time compression: requested x%.0f | applied x%.0f | allowed max x%.0f %s",
                  timeScaleForIndex(requestedTimeIdx),
                  timeScaleForIndex(appliedTimeIdx),
                  timeScaleForIndex(allowedIdx),
                  forceTime ? "[FORCED]" : "");

      ImGui::Text("Flight mode: %s%s", docked ? "DOCKED" : (supercruise ? "SUPERCRUISE" : "NORMAL"),
                  autopilot ? " + AUTOPILOT" : "");

      ImGui::Separator();
      ImGui::Text("Pos (km):   [%.1f %.1f %.1f]", pos.x, pos.y, pos.z);
      ImGui::Text("Vel (km/s): [%.3f %.3f %.3f] |v|=%.3f", vel.x, vel.y, vel.z, vel.length());
      ImGui::Text("AngVel (rad/s): [%.3f %.3f %.3f]", wv.x, wv.y, wv.z);

      if (supercruise) {
        ImGui::Separator();
        ImGui::Text("Supercruise speed: %.0f km/s  throttle: %.2f", scSpeedKmS, scThrottle);
        ImGui::Checkbox("Auto-drop", &scAutoDrop);
      }

      ImGui::Separator();
      ImGui::TextDisabled("Controls:");
      ImGui::BulletText("Translate: WASD + R/F (up/down)");
      ImGui::BulletText("Rotate: Arrow keys + Q/E roll");
      ImGui::BulletText("Boost: LShift   Brake: X");
      ImGui::BulletText("Dampers: Z (on) / C (off)");
      ImGui::BulletText("Target cycle: T  Clear: Y");
      ImGui::BulletText("Autopilot: P   Dock/Undock: G");
      ImGui::BulletText("Supercruise toggle: J (W/S throttle while active)");
      ImGui::BulletText("Time compression: PgUp/PgDn, Home=1x, End=max (hold Ctrl to force)");
      ImGui::BulletText("Pause: Space   Save: F5   Load: F9");
      ImGui::BulletText("Windows: TAB Galaxy, F1 Flight, F2 Economy, F3 Nav, F4 HUD");

      ImGui::End();
    }

    if (showNav) {
      ImGui::Begin("System / Nav");

      if (navTarget.kind == TargetKind::None) {
        ImGui::Text("Target: (none)");
      } else {
        BodyState bs{};
        if (getTargetState(*currentSystem, navTarget, timeDays, &bs)) {
          const double dist = (bs.posKm - ship.positionKm()).length();
          const double relSpeed = (ship.velocityKmS() - bs.velKmS).length();
          ImGui::Text("Target: %s - %s", bs.kindLabel, bs.name.c_str());
          ImGui::Text("Distance: %.0f km", dist);
          ImGui::Text("Relative speed: %.3f km/s", relSpeed);

          if (navTarget.kind == TargetKind::Station) {
            const auto& st = currentSystem->stations[navTarget.index];
            const DockCheck dc = canDock(ship, st, bs.posKm, bs.velKmS);
            ImGui::Text("Docking: %s | corridor=%s | vRel=%.3f (limit %.3f)",
                        dc.ok ? "OK" : "NO",
                        dc.inCorridor ? "IN" : "OUT",
                        dc.relSpeedKmS,
                        st.dockingSpeedLimitKmS);
          }
        }
      }

      ImGui::Separator();
      ImGui::Checkbox("Autopilot", &autopilot);
      ImGui::SameLine();
      ImGui::Checkbox("Supercruise", &supercruise);
      if (ImGui::IsItemHovered()) ImGui::SetTooltip("Use J hotkey for safety checks");

      ImGui::Separator();
      ImGui::Text("Stations:");
      for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
        const auto& st = currentSystem->stations[i];
        bool selected = (navTarget.kind == TargetKind::Station && navTarget.index == i);
        if (ImGui::Selectable(st.name.c_str(), selected)) {
          navTarget.kind = TargetKind::Station;
          navTarget.index = i;
        }
      }

      ImGui::Separator();
      ImGui::Text("Planets:");
      for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
        const auto& p = currentSystem->planets[i];
        bool selected = (navTarget.kind == TargetKind::Planet && navTarget.index == i);
        if (ImGui::Selectable(p.name.c_str(), selected)) {
          navTarget.kind = TargetKind::Planet;
          navTarget.index = i;
        }
      }

      ImGui::End();
    }

    if (showEconomy) {
      beginStationSelectHUD(*currentSystem, dockedStationIndex, "Dock / Market");

      ImGui::Begin("Dock / Status");
      if (docked) {
        const auto& st = currentSystem->stations[(std::size_t)dockedStationIndex];
        ImGui::Text("Docked at: %s", st.name.c_str());
      } else {
        ImGui::TextDisabled("Not docked. Fly to a station, align in corridor, slow down, press G.");
      }
      ImGui::End();

      if (!currentSystem->stations.empty()) {
        const auto& station = currentSystem->stations[(std::size_t)dockedStationIndex];
        auto& stEcon = universe.stationEconomy(station, timeDays);

        ImGui::Begin("Market Details");
        ImGui::Text("Credits: %.2f", credits);

        const bool atThisStation = docked && (station.id == dockedStationId);
        if (!atThisStation) {
          ImGui::TextDisabled("Trade is disabled unless you are docked at the selected station.");
        }

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

            if (!atThisStation) ImGui::BeginDisabled(true);

            ImGui::SameLine();
            if (ImGui::SmallButton("Buy")) {
              auto tr = econ::buy(stEcon, station.economyModel, cid, qty[i], credits, 0.10, station.feeRate);
              if (tr.ok) cargo[i] += qty[i];
            }

            ImGui::SameLine();
            if (ImGui::SmallButton("Sell")) {
              const double sellUnits = std::min<double>(qty[i], cargo[i]);
              if (sellUnits > 0.0) {
                auto tr = econ::sell(stEcon, station.economyModel, cid, sellUnits, credits, 0.10, station.feeRate);
                if (tr.ok) cargo[i] -= sellUnits;
              }
            }

            if (!atThisStation) ImGui::EndDisabled();

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
          ImGui::TextDisabled("No history yet (time needs to advance). Use PgUp/PgDn time compression.");
        }

        ImGui::End();
      }
    }

    if (showGalaxy) {
      ImGui::Begin("Galaxy / Streaming");

      const auto center = currentSystem->stub.posLy;
      static float radius = 200.0f;
      ImGui::SliderFloat("Query radius (ly)", &radius, 20.0f, 1200.0f);

      auto nearby = universe.queryNearby(center, radius, 128);
      ImGui::Text("Nearby systems: %d", (int)nearby.size());

      // Mini-map canvas
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

            // Reset transient state
            docked = false;
            dockedStationId = 0;
            dockedStationIndex = 0;
            autopilot = false;
            supercruise = false;
            scThrottle = 0.0;
            scSpeedKmS = 0.0;

            navTarget.clear();

            // Spawn near first station if available.
            if (!sys.stations.empty()) {
              const auto& st = sys.stations[0];
              math::Vec3d stPosKm{}, stVelKmS{};
              orbitStateKm(st.orbit, timeDays, &stPosKm, &stVelKmS);
              const math::Vec3d axisOut = safeNormalize(stPosKm, {0, 0, 1});
              ship.setPositionKm(stPosKm + axisOut * std::max(1200.0, st.dockingCorridorLengthKm * 0.85));
              ship.setVelocityKmS(stVelKmS);
              navTarget.kind = TargetKind::Station;
              navTarget.index = 0;
            } else {
              ship.setPositionKm({0, 0, -8000.0});
              ship.setVelocityKmS({0, 0, 0});
            }

            toast.show("Jumped to " + sys.stub.name);
          }
        }
      }

      ImGui::TextDisabled("Tip: TAB toggles this window, F1 Flight, F2 Economy, F3 Nav, F4 HUD");

      ImGui::End();
    }

    // HUD overlay (target marker + crosshair)
    if (showHUD) {
      ImDrawList* fg = ImGui::GetForegroundDrawList();
      const ImGuiViewport* vp = ImGui::GetMainViewport();
      const ImVec2 vpPos = vp->Pos;
      const ImVec2 vpSize = vp->Size;
      const ImVec2 center = ImVec2(vpPos.x + vpSize.x * 0.5f, vpPos.y + vpSize.y * 0.5f);

      // Crosshair
      const float s = 6.0f;
      fg->AddLine(ImVec2(center.x - s, center.y), ImVec2(center.x + s, center.y), IM_COL32(220, 220, 240, 170), 1.0f);
      fg->AddLine(ImVec2(center.x, center.y - s), ImVec2(center.x, center.y + s), IM_COL32(220, 220, 240, 170), 1.0f);

      // Target marker
      if (navTarget.isValid(*currentSystem)) {
        BodyState bs{};
        if (getTargetState(*currentSystem, navTarget, timeDays, &bs)) {
          const math::Vec3d targetU = bs.posKm * (1.0 / kRENDER_UNIT_KM);
          ImVec2 px{}, ndc{};
          bool onScreen = false;
          if (projectToScreen(targetU, view, proj, vpPos, vpSize, &px, &onScreen, &ndc)) {
            ImVec2 drawPos = px;
            if (!onScreen) {
              // Clamp to edge.
              const float ax = ndc.x;
              const float ay = ndc.y;
              const float m = std::max(std::abs(ax), std::abs(ay));
              const float kEdge = 0.92f;
              const float sx = (m > 1e-5f) ? (ax / m) * kEdge : 0.0f;
              const float sy = (m > 1e-5f) ? (ay / m) * kEdge : 0.0f;
              const float px2 = vpPos.x + (sx * 0.5f + 0.5f) * vpSize.x;
              const float py2 = vpPos.y + (1.0f - (sy * 0.5f + 0.5f)) * vpSize.y;
              drawPos = ImVec2(px2, py2);
            }

            fg->AddCircle(drawPos, 10.0f, IM_COL32(255, 140, 80, 220), 0, 1.5f);
            const double dist = (bs.posKm - ship.positionKm()).length();
            const std::string label = bs.name + "  " + std::to_string((long long)dist) + " km";
            fg->AddText(ImVec2(drawPos.x + 12.0f, drawPos.y - 6.0f), IM_COL32(255, 200, 160, 220), label.c_str());
          }
        }
      }

      // Mini status top-left
      ImGui::SetNextWindowPos(ImVec2(vpPos.x + 10.0f, vpPos.y + 10.0f));
      ImGui::SetNextWindowBgAlpha(0.35f);
      ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize |
                               ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing |
                               ImGuiWindowFlags_NoNav;
      ImGui::Begin("##hud", nullptr, flags);
      ImGui::Text("%s", currentSystem->stub.name.c_str());
      ImGui::Text("x%.0f %s", timeScale, paused ? "PAUSED" : "");
      if (autopilot) ImGui::Text("AUTOPILOT");
      if (supercruise) ImGui::Text("SUPERCRUISE %.0f km/s", scSpeedKmS);
      if (docked) ImGui::Text("DOCKED");
      ImGui::End();
    }

    // Toast overlay
    if (toast.timer > 0.0 && !toast.msg.empty()) {
      const ImGuiViewport* vp = ImGui::GetMainViewport();
      ImGui::SetNextWindowPos(ImVec2(vp->Pos.x + 10.0f, vp->Pos.y + vp->Size.y - 60.0f));
      ImGui::SetNextWindowBgAlpha(0.55f);
      ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize |
                               ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing |
                               ImGuiWindowFlags_NoNav;
      ImGui::Begin("##toast", nullptr, flags);
      ImGui::TextUnformatted(toast.msg.c_str());
      ImGui::End();
    }

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
