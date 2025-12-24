#include "stellar/core/Log.h"
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

static constexpr double kC_KM_S = 299792.458; // speed of light (km/s)

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

static bool beginDockedHUD(const sim::StarSystem& sys, int& stationIndex) {
  bool changed = false;

  ImGui::Begin("Dock / Market");

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

static math::Vec3d planetPosKm(const sim::Planet& p, double timeDays) {
  const math::Vec3d posAU = sim::orbitPosition3DAU(p.orbit, timeDays);
  return posAU * kAU_KM;
}

static double planetRadiusKm(const sim::Planet& p) {
  return p.radiusEarth * kEARTH_RADIUS_KM;
}

static math::Vec3d stationPosKm(const sim::Station& st, double timeDays) {
  const math::Vec3d posAU = sim::orbitPosition3DAU(st.orbit, timeDays);
  return posAU * kAU_KM;
}

static math::Vec3d stationOutDir(const math::Vec3d& stationPosKm) {
  // Station "out" direction points away from the primary star at origin.
  const math::Vec3d d = stationPosKm.normalized();
  return (d.lengthSq() > 1e-12) ? d : math::Vec3d{0,0,1};
}

struct ProximityInfo {
  double distToStarKm{0.0};

  double nearestPlanetDistKm{1.0e30};
  double nearestPlanetRadiusKm{0.0};
  int nearestPlanetIndex{-1};

  double nearestStationDistKm{1.0e30};
  int nearestStationIndex{-1};
};

static ProximityInfo computeProximity(const sim::StarSystem& sys,
                                     const math::Vec3d& shipPosKm,
                                     double timeDays) {
  ProximityInfo p{};
  p.distToStarKm = shipPosKm.length();

  for (std::size_t i = 0; i < sys.planets.size(); ++i) {
    const double rKm = planetRadiusKm(sys.planets[i]);
    const double dKm = (shipPosKm - planetPosKm(sys.planets[i], timeDays)).length();
    if (dKm < p.nearestPlanetDistKm) {
      p.nearestPlanetDistKm = dKm;
      p.nearestPlanetRadiusKm = rKm;
      p.nearestPlanetIndex = (int)i;
    }
  }

  for (std::size_t i = 0; i < sys.stations.size(); ++i) {
    const double dKm = (shipPosKm - stationPosKm(sys.stations[i], timeDays)).length();
    if (dKm < p.nearestStationDistKm) {
      p.nearestStationDistKm = dKm;
      p.nearestStationIndex = (int)i;
    }
  }

  return p;
}

static double moveToward(double current, double target, double maxDelta) {
  if (current < target) return std::min(target, current + maxDelta);
  if (current > target) return std::max(target, current - maxDelta);
  return current;
}

static double allowedTimeScale(const ProximityInfo& prox,
                              double shipSpeedKmS,
                              bool docked,
                              bool inSupercruise) {
  // Conservative defaults. This is a *prototype safety system*:
  // - close to planets/stations: keep time compression low
  // - in supercruise: avoid extremely high compression (overshoot / instability)
  if (docked) return 200000.0;

  double limit = 10000.0;
  if (inSupercruise) limit = std::min(limit, 200.0);

  // Close to stations: keep it low.
  if (prox.nearestStationIndex >= 0) {
    if (prox.nearestStationDistKm < 250.0) limit = std::min(limit, 5.0);
    else if (prox.nearestStationDistKm < 2500.0) limit = std::min(limit, 20.0);
    else if (prox.nearestStationDistKm < 20000.0) limit = std::min(limit, 200.0);
  }

  // Close to planets: keep it low (use radius-scaled bands).
  if (prox.nearestPlanetIndex >= 0 && prox.nearestPlanetRadiusKm > 0.0) {
    const double r = prox.nearestPlanetRadiusKm;
    const double d = prox.nearestPlanetDistKm;
    if (d < r * 5.0) limit = std::min(limit, 5.0);
    else if (d < r * 20.0) limit = std::min(limit, 20.0);
    else if (d < r * 80.0) limit = std::min(limit, 200.0);
  }

  // High speed: reduce max time compression.
  if (shipSpeedKmS > 2.0) limit = std::min(limit, 200.0);
  if (shipSpeedKmS > 20.0) limit = std::min(limit, 20.0);

  // Star proximity also matters.
  if (prox.distToStarKm < 250000.0) limit = std::min(limit, 20.0);
  if (prox.distToStarKm < 50000.0) limit = std::min(limit, 5.0);

  return std::max(1.0, limit);
}

static double supercruiseMaxSpeedKmS(const ProximityInfo& prox) {
  // Simple "gravity well" scaling inspired by Elite's notion that supercruise
  // max speed depends heavily on distance from large masses.
  // We compute the nearest "mass" distance and scale a capped power-law.
  const double d = std::max(1.0, std::min({prox.distToStarKm, prox.nearestPlanetDistKm, prox.nearestStationDistKm}));

  // Tunable curve: ~50 km/s at 1,000 km, rising toward a cap as distance grows.
  const double base = 50.0;
  const double scale = 1000.0;
  const double exp = 1.2;
  const double cap = 300000.0; // 1c in km/s

  const double v = base * std::pow(d / scale, exp);
  return std::min(cap, std::max(50.0, v));
}

int main(int argc, char** argv) {
  (void)argc; (void)argv;

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
      SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
      1280, 720,
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
  render::Mesh cube   = render::Mesh::makeCube();

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

  // --- Universe / sim state ---
  core::u64 seed = 1337;
  sim::Universe universe(seed);

  sim::SystemStub currentStub{};
  if (auto s = universe.findClosestSystem({0,0,0}, 80.0)) {
    currentStub = *s;
  } else {
    // Should be rare; fall back
    currentStub = sim::SystemStub{};
    currentStub.id = 0;
    currentStub.seed = seed;
    currentStub.name = "Origin";
    currentStub.posLy = {0,0,0};
    currentStub.planetCount = 6;
    currentStub.stationCount = 1;
  }

  const sim::StarSystem* currentSystem = &universe.getSystem(currentStub.id, &currentStub);

  sim::Ship ship;
  // Spawn near the first station (when present) to reach an early-playable docking loop quickly.
  double timeDays = 0.0;
  if (!currentSystem->stations.empty()) {
    const auto& st0 = currentSystem->stations[0];
    const math::Vec3d stPosKm = stationPosKm(st0, timeDays);
    const math::Vec3d outDir = stationOutDir(stPosKm);
    ship.setPositionKm(stPosKm + outDir * (st0.radiusKm + 25.0));
  } else {
    ship.setPositionKm({0, 0, -8000.0}); // fallback
  }
  ship.setMaxLinearAccelKmS2(0.08);
  ship.setMaxAngularAccelRadS2(1.2);

  // --- Time acceleration ---
  // Time compression is expressed as "sim seconds per real second".
  // (1 = real-time, 10 = 10x, etc.)
  static constexpr std::array<double, 12> kTimeLevels = {
    1.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0
  };
  int timeLevelIndex = 2; // default ~10x
  double timeScale = kTimeLevels[(std::size_t)timeLevelIndex];
  bool paused = false;

  // --- Flight modes ---
  bool docked = false;
  sim::StationId dockedStationId = 0;

  bool supercruise = false;
  double supercruiseThrottle = 0.0; // 0..1
  double supercruiseSpeedKmS = 0.0;

  // Player economy
  double credits = 2500.0;
  std::array<double, econ::kCommodityCount> cargo{};
  int selectedStationIndex = 0;
  int dockedStationIndex = 0; // valid when `docked == true`

  // Save/load
  const std::string savePath = "savegame.txt";

  // UI state
  bool showGalaxy = true;
  bool showShip = true;
  bool showEconomy = true;

  // Destination for route planner
  sim::SystemId routeToSystem = 0;
  sim::StationId routeToStation = 0;

  bool running = true;
  auto last = std::chrono::high_resolution_clock::now();

  // One-line HUD/toast message for player feedback.
  std::string toast;
  double toastTimer = 0.0;
  auto setToast = [&](const std::string& msg) {
    toast = msg;
    toastTimer = 4.0;
    core::log(core::LogLevel::Info, msg);
  };

  bool requestDockToggle = false;
  bool requestSupercruiseToggle = false;

  SDL_SetRelativeMouseMode(SDL_FALSE);

  while (running) {
    // Timing
    auto now = std::chrono::high_resolution_clock::now();
    const double dt = std::chrono::duration<double>(now - last).count();
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
          s.seed = universe.seed();
          s.timeDays = timeDays;
          s.currentSystem = currentSystem->stub.id;
          // Only persist the *actual* docked station.
          s.dockedStation = docked ? dockedStationId : 0;

          s.shipPosKm = ship.positionKm();
          s.shipVelKmS = ship.velocityKmS();
          s.shipOrient = ship.orientation();
          s.shipAngVelRadS = ship.angularVelocityRadS();

          s.credits = credits;
          s.cargo = cargo;
          s.stationOverrides = universe.exportStationOverrides();

          if (sim::saveToFile(s, savePath)) {
            core::log(core::LogLevel::Info, "Saved to " + savePath);
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

            // Restore docking state.
            docked = (s.dockedStation != 0);
            dockedStationId = s.dockedStation;
            dockedStationIndex = 0;
            for (std::size_t i = 0; i < sys.stations.size(); ++i) {
              if (sys.stations[i].id == s.dockedStation) dockedStationIndex = (int)i;
            }
            selectedStationIndex = dockedStationIndex;

            // Supercruise is transient.
            supercruise = false;
            supercruiseThrottle = 0.0;
            supercruiseSpeedKmS = 0.0;

            if (docked && !sys.stations.empty()) {
              const auto& st = sys.stations[(std::size_t)dockedStationIndex];
              ship.setPositionKm(stationPosKm(st, timeDays));
              ship.setVelocityKmS({0,0,0});
              ship.setAngularVelocityRadS({0,0,0});
            }

            core::log(core::LogLevel::Info, "Loaded " + savePath);
          }
        }

        if (!captureKeys) {
          if (event.key.keysym.sym == SDLK_TAB) showGalaxy = !showGalaxy;
          if (event.key.keysym.sym == SDLK_F1) showShip = !showShip;
          if (event.key.keysym.sym == SDLK_F2) showEconomy = !showEconomy;
          if (event.key.keysym.sym == SDLK_SPACE) paused = !paused;

          // Time compression (PageUp/PageDown like classic space sims).
          if (event.key.keysym.sym == SDLK_PAGEUP) {
            timeLevelIndex = std::min(timeLevelIndex + 1, (int)kTimeLevels.size() - 1);
          }
          if (event.key.keysym.sym == SDLK_PAGEDOWN) {
            timeLevelIndex = std::max(timeLevelIndex - 1, 0);
          }
          if (event.key.keysym.sym == SDLK_HOME) timeLevelIndex = 0;
          if (event.key.keysym.sym == SDLK_END) timeLevelIndex = (int)kTimeLevels.size() - 1;

          // Dock / undock.
          if (event.key.keysym.sym == SDLK_g) requestDockToggle = true;

          // Supercruise.
          if (event.key.keysym.sym == SDLK_j) requestSupercruiseToggle = true;
        }
      }
    }

    // --- Keyboard state ---
    sim::ShipInput input{};
    const Uint8* keys = SDL_GetKeyboardState(nullptr);
    const bool captureKeys = io.WantCaptureKeyboard;
    const bool ctrlHeld = (keys[SDL_SCANCODE_LCTRL] || keys[SDL_SCANCODE_RCTRL]);

    // --- Time accel (with safety clamp) ---
    const double requestedTimeScale = kTimeLevels[(std::size_t)timeLevelIndex];
    const ProximityInfo prox = computeProximity(*currentSystem, ship.positionKm(), timeDays);
    const double speedKmS = ship.velocityKmS().length();
    const double allowed = allowedTimeScale(prox, speedKmS, docked, supercruise);
    timeScale = ctrlHeld ? requestedTimeScale : std::min(requestedTimeScale, allowed);

    // --- Docking corridor status for the *selected* station ---
    const sim::Station* selStation = nullptr;
    math::Vec3d selStationPosKm{0,0,0};
    math::Vec3d selStationOut{0,0,1};
    double selStationDistKm = 1.0e30;
    double corridorAlongKm = 0.0;
    double corridorRadialKm = 0.0;
    double corridorAlignCos = -1.0;
    bool inCorridor = false;
    bool alignedForDock = false;
    bool speeding = false;
    bool canDock = false;

    if (!currentSystem->stations.empty()) {
      // When docked, your target/market station is the station you're docked at.
      if (docked) selectedStationIndex = dockedStationIndex;

      const int rawIndex = docked ? dockedStationIndex : selectedStationIndex;
      const int clampedIndex = std::max(0, std::min(rawIndex, (int)currentSystem->stations.size() - 1));
      if (docked) dockedStationIndex = clampedIndex;
      else selectedStationIndex = clampedIndex;

      selStation = &currentSystem->stations[(std::size_t)clampedIndex];
      selStationPosKm = stationPosKm(*selStation, timeDays);
      selStationOut = stationOutDir(selStationPosKm);

      const math::Vec3d rel = ship.positionKm() - selStationPosKm;
      selStationDistKm = rel.length();

      corridorAlongKm = math::dot(rel, selStationOut);
      const math::Vec3d radialV = rel - selStationOut * corridorAlongKm;
      corridorRadialKm = radialV.length();

      // Ship should be facing *toward* the station when approaching (i.e. -out).
      corridorAlignCos = math::dot(ship.forward(), selStationOut * (-1.0));
      const double cosReq = std::cos(math::degToRad(selStation->corridorAlignHalfAngleDeg));
      alignedForDock = corridorAlignCos >= cosReq;

      inCorridor = (corridorAlongKm > 0.0 && corridorAlongKm < selStation->corridorLengthKm &&
                    corridorRadialKm < selStation->corridorRadiusKm);
      speeding = inCorridor && (speedKmS > selStation->corridorSpeedLimitKmS);

      // Dock when close, inside corridor, aligned, and under speed limit.
      const double dockDist = std::max(selStation->radiusKm * 1.2, selStation->radiusKm + 5.0);
      canDock = (!supercruise && !docked && inCorridor && alignedForDock &&
                 speedKmS <= selStation->corridorSpeedLimitKmS && selStationDistKm <= dockDist);
    }

    // --- Handle one-shot requests ---
    if (requestDockToggle) {
      requestDockToggle = false;

      if (!selStation) {
        setToast("No station selected.");
      } else if (supercruise) {
        setToast("Cannot dock while in supercruise.");
      } else if (docked) {
        // Undock: place ship just outside the corridor, gently pushing away from the star.
        docked = false;
        dockedStationId = 0;
        ship.setVelocityKmS(selStationOut * (selStation->corridorSpeedLimitKmS * 0.5));
        ship.setAngularVelocityRadS({0,0,0});
        ship.setPositionKm(selStationPosKm + selStationOut * (selStation->radiusKm + 25.0));
        setToast("Undocked from " + selStation->name);
      } else {
        if (canDock) {
          docked = true;
          dockedStationId = selStation->id;
          dockedStationIndex = selectedStationIndex;
          selectedStationIndex = dockedStationIndex;
          ship.setVelocityKmS({0,0,0});
          ship.setAngularVelocityRadS({0,0,0});
          ship.setPositionKm(selStationPosKm);
          setToast("Docked at " + selStation->name);
        } else {
          if (!inCorridor) setToast("Docking failed: not in approach corridor.");
          else if (!alignedForDock) setToast("Docking failed: align with corridor." );
          else if (speedKmS > selStation->corridorSpeedLimitKmS) setToast("Docking failed: SPEEDING (slow down)." );
          else setToast("Docking failed: get closer to the station." );
        }
      }
    }

    if (requestSupercruiseToggle) {
      requestSupercruiseToggle = false;

      if (docked) {
        setToast("Cannot enter supercruise while docked.");
      } else if (supercruise) {
        supercruise = false;
        // Drop-out: clamp to a sane speed.
        const double exitSpeed = std::min(2.0, supercruiseSpeedKmS);
        ship.setVelocityKmS(ship.forward() * exitSpeed);
        supercruiseSpeedKmS = exitSpeed;
        setToast("Exited supercruise");
      } else {
        // Simple mass-lock style constraint.
        const bool tooCloseStation = (prox.nearestStationDistKm < 8000.0);
        const bool tooClosePlanet = (prox.nearestPlanetIndex >= 0 && prox.nearestPlanetRadiusKm > 0.0 &&
                                     prox.nearestPlanetDistKm < prox.nearestPlanetRadiusKm * 80.0);
        const bool tooCloseStar = (prox.distToStarKm < 400000.0);

        if (tooCloseStation || tooClosePlanet || tooCloseStar) {
          setToast("Mass-locked: move farther from star/planet/station to enter supercruise.");
        } else {
          supercruise = true;
          supercruiseThrottle = 0.50;
          supercruiseSpeedKmS = std::max(1.0, ship.velocityKmS().length());
          ship.setVelocityKmS(ship.forward() * supercruiseSpeedKmS);
          setToast("Entered supercruise");
        }
      }
    }

    // --- Build flight input by mode ---
    if (!captureKeys) {
      // Rotation is always available in flight.
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

      if (!docked) {
        if (!supercruise) {
          // Normal 6DOF translation.
          input.thrustLocal.z += (keys[SDL_SCANCODE_W] ? 1.0 : 0.0);
          input.thrustLocal.z -= (keys[SDL_SCANCODE_S] ? 1.0 : 0.0);

          input.thrustLocal.x += (keys[SDL_SCANCODE_D] ? 1.0 : 0.0);
          input.thrustLocal.x -= (keys[SDL_SCANCODE_A] ? 1.0 : 0.0);

          input.thrustLocal.y += (keys[SDL_SCANCODE_R] ? 1.0 : 0.0);
          input.thrustLocal.y -= (keys[SDL_SCANCODE_F] ? 1.0 : 0.0);
        } else {
          // Supercruise uses a throttle (W/S) instead of direct translation.
          const double throttleRate = 0.55;
          if (keys[SDL_SCANCODE_W]) supercruiseThrottle += dt * throttleRate;
          if (keys[SDL_SCANCODE_S]) supercruiseThrottle -= dt * throttleRate;
          if (keys[SDL_SCANCODE_X]) supercruiseThrottle -= dt * (throttleRate * 2.5);
          supercruiseThrottle = math::clamp(supercruiseThrottle, 0.0, 1.0);

          // Prevent ship dampers from fighting the supercruise velocity.
          input.dampers = true;
        }
      }
    }

    // Toast lifetime
    if (toastTimer > 0.0) {
      toastTimer -= dt;
      if (toastTimer < 0.0) toastTimer = 0.0;
    }

    // --- Simulation step (time compression applies to ship + orbits + economy) ---
    if (!paused) {
      double remaining = dt * timeScale;
      const double maxStep = (timeScale > 1000.0) ? 2.0 : 0.5; // keep loops bounded

      while (remaining > 0.0) {
        const double step = std::min(remaining, maxStep);

        if (docked && selStation) {
          // While docked, the ship is "attached" to the station.
          timeDays += step / 86400.0;
          const math::Vec3d pKm = stationPosKm(*selStation, timeDays);
          ship.setPositionKm(pKm);
          ship.setVelocityKmS({0,0,0});
          ship.setAngularVelocityRadS({0,0,0});
        } else if (supercruise) {
          // Supercruise: simplified model.
          // 1) Integrate attitude only (use Ship's angular physics).
          ship.setVelocityKmS({0,0,0});
          ship.step(step, input);

          // 2) Update speed toward throttle*max.
          const double vMax = supercruiseMaxSpeedKmS(prox);
          const double targetSpeed = supercruiseThrottle * vMax;
          const double accel = 200.0; // km/s^2, arcade-ish but responsive
          supercruiseSpeedKmS = moveToward(supercruiseSpeedKmS, targetSpeed, accel * step);

          ship.setVelocityKmS(ship.forward() * supercruiseSpeedKmS);
          ship.setPositionKm(ship.positionKm() + ship.velocityKmS() * step);

          timeDays += step / 86400.0;
        } else {
          ship.step(step, input);
          timeDays += step / 86400.0;
        }

        remaining -= step;
      }
    }

    // ---- Camera follow (third-person) ----
    render::Camera cam;
    int w = 1280, h = 720;
    SDL_GetWindowSize(window, &w, &h);
    const double aspect = (h > 0) ? (double)w / (double)h : 16.0/9.0;

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

    // ---- Build instances (star + planets + ship) ----
    std::vector<render::InstanceData> spheres;
    spheres.reserve(1 + currentSystem->planets.size());

    // Star at origin
    {
      const double starRadiusKm = currentSystem->star.radiusSol * kSOLAR_RADIUS_KM;
      const float starScale = (float)std::max(0.8, (starRadiusKm / kRENDER_UNIT_KM) * 3.0);
      spheres.push_back({0,0,0, starScale, 1.0f, 0.95f, 0.75f});
    }

    // Planets
    for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
      const auto& p = currentSystem->planets[i];
      const math::Vec3d posAU = sim::orbitPosition3DAU(p.orbit, timeDays);
      const math::Vec3d posKm = posAU * kAU_KM;
      const math::Vec3d posU = posKm * (1.0 / kRENDER_UNIT_KM);

      const double radiusKm = p.radiusEarth * kEARTH_RADIUS_KM;
      const float scale = (float)std::max(0.25, (radiusKm / kRENDER_UNIT_KM) * 200.0);

      // Simple color palette by type
      float cr=0.6f, cg=0.6f, cb=0.6f;
      switch (p.type) {
        case sim::PlanetType::Rocky: cr=0.6f; cg=0.55f; cb=0.5f; break;
        case sim::PlanetType::Desert: cr=0.8f; cg=0.7f; cb=0.35f; break;
        case sim::PlanetType::Ocean: cr=0.25f; cg=0.45f; cb=0.85f; break;
        case sim::PlanetType::Ice: cr=0.7f; cg=0.85f; cb=0.95f; break;
        case sim::PlanetType::GasGiant: cr=0.7f; cg=0.55f; cb=0.35f; break;
        default: break;
      }

      spheres.push_back({(float)posU.x, (float)posU.y, (float)posU.z, scale, cr,cg,cb});
    }

    // Orbit lines
    std::vector<render::LineVertex> orbitLines;
    orbitLines.reserve(currentSystem->planets.size() * 128);

    for (const auto& p : currentSystem->planets) {
      const int seg = 96;
      math::Vec3d prev{};
      for (int s = 0; s <= seg; ++s) {
        const double t = (double)s / (double)seg * p.orbit.periodDays;
        const math::Vec3d posAU = sim::orbitPosition3DAU(p.orbit, t);
        const math::Vec3d posU = (posAU * kAU_KM) * (1.0 / kRENDER_UNIT_KM);

        if (s > 0) {
          orbitLines.push_back({(float)prev.x,(float)prev.y,(float)prev.z, 0.25f,0.25f,0.25f});
          orbitLines.push_back({(float)posU.x,(float)posU.y,(float)posU.z, 0.25f,0.25f,0.25f});
        }
        prev = posU;
      }
    }

    // Cubes: ship + stations (prototype visuals)
    std::vector<render::InstanceData> cubes;
    cubes.reserve(1 + currentSystem->stations.size());

    // Ship
    {
      const float s = supercruise ? 0.40f : (docked ? 0.30f : 0.35f);
      const float r = supercruise ? 0.7f : 0.9f;
      const float g = docked ? 0.95f : 0.9f;
      const float b = 1.0f;
      cubes.push_back(render::InstanceData{
        (float)shipPosU.x, (float)shipPosU.y, (float)shipPosU.z,
        s,
        r, g, b
      });
    }

    // Stations
    for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
      const auto& st = currentSystem->stations[i];
      const math::Vec3d pKm = stationPosKm(st, timeDays);
      const math::Vec3d pU = pKm * (1.0 / kRENDER_UNIT_KM);

      const bool isDockedStation = docked && ((int)i == dockedStationIndex);
      const bool isSelectedStation = ((int)i == selectedStationIndex);

      float cr = 0.6f, cg = 0.6f, cb = 0.65f;
      if (isSelectedStation) { cr = 0.95f; cg = 0.85f; cb = 0.25f; }
      if (isDockedStation) { cr = 0.25f; cg = 0.95f; cb = 0.35f; }

      cubes.push_back(render::InstanceData{
        (float)pU.x, (float)pU.y, (float)pU.z,
        0.55f,
        cr, cg, cb
      });
    }

    // ---- Render ---
    glViewport(0, 0, w, h);
    glClearColor(0.01f, 0.01f, 0.02f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Orbits (lines)
    lineRenderer.drawLines(orbitLines);

    // Star + planets
    meshRenderer.setMesh(&sphere);
    meshRenderer.drawInstances(spheres);

    // Ship + stations (as cubes)
    meshRenderer.setMesh(&cube);
    meshRenderer.drawInstances(cubes);

    // ---- UI ----
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplSDL2_NewFrame(window);
    ImGui::NewFrame();

    // Main dockspace
    ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());

    if (showShip) {
      ImGui::Begin("Ship / Flight");

      const auto pos = ship.positionKm();
      const auto vel = ship.velocityKmS();
      const auto wv  = ship.angularVelocityRadS();

      // Flight mode banner
      if (docked && !currentSystem->stations.empty()) {
        const auto& st = currentSystem->stations[(std::size_t)dockedStationIndex];
        ImGui::TextColored(ImVec4(0.25f,0.95f,0.35f,1.0f), "DOCKED: %s", st.name.c_str());
      } else if (supercruise) {
        ImGui::TextColored(ImVec4(0.7f,0.9f,1.0f,1.0f), "SUPERCRUISE");
      } else {
        ImGui::Text("NORMAL FLIGHT");
      }

      // Time accel readout
      ImGui::Text("Time: %.2f days", timeDays);
      if (paused) {
        ImGui::SameLine();
        ImGui::TextColored(ImVec4(1.0f,0.6f,0.2f,1.0f), "[PAUSED]");
      }
      ImGui::Text("Time compression: x%.0f (requested x%.0f, allowed x%.0f)%s",
                  timeScale,
                  requestedTimeScale,
                  allowed,
                  ctrlHeld ? "  [CTRL OVERRIDE]" : "");
      if (!ctrlHeld && requestedTimeScale > allowed) {
        ImGui::TextColored(ImVec4(1.0f,0.75f,0.3f,1.0f),
                           "Limited for safety (near traffic/masses). Hold CTRL to override.");
      }

      ImGui::Separator();
      ImGui::Text("Pos (km):   [%.1f %.1f %.1f]", pos.x, pos.y, pos.z);
      ImGui::Text("Vel (km/s): [%.3f %.3f %.3f] |v|=%.3f", vel.x, vel.y, vel.z, vel.length());
      ImGui::Text("AngVel (rad/s): [%.3f %.3f %.3f]", wv.x, wv.y, wv.z);

      ImGui::Separator();

      if (selStation) {
        ImGui::Text("Target station: %s", selStation->name.c_str());
        ImGui::Text("Dist: %.1f km | Corridor: %s", selStationDistKm, inCorridor ? "IN" : "OUT");
        ImGui::Text("Corridor along: %.1f km  radial: %.1f km", corridorAlongKm, corridorRadialKm);
        const double cosReq = std::cos(math::degToRad(selStation->corridorAlignHalfAngleDeg));
        ImGui::Text("Alignment cos: %.3f (req >= %.3f)  %s", corridorAlignCos, cosReq,
                    alignedForDock ? "[OK]" : "[NOT ALIGNED]");
        ImGui::Text("Speed: %.3f km/s  limit: %.3f km/s  %s", speedKmS, selStation->corridorSpeedLimitKmS,
                    speeding ? "[SPEEDING]" : "");
        if (!docked && !supercruise) {
          ImGui::Text("Docking: %s  (press G)", canDock ? "READY" : "NOT READY");
        }
      }

      if (supercruise) {
        ImGui::Separator();
        const double vMax = supercruiseMaxSpeedKmS(prox);
        ImGui::Text("Supercruise throttle: %.0f%%", supercruiseThrottle * 100.0);
        ImGui::Text("Supercruise speed: %.1f km/s (max %.1f km/s)", supercruiseSpeedKmS, vMax);
      }

      ImGui::Separator();
      ImGui::TextDisabled("Controls:");
      ImGui::BulletText("Rotate: Arrow keys + Q/E roll");
      ImGui::BulletText("Translate (normal): WASD + R/F");
      ImGui::BulletText("Dock/Undock: G  | Supercruise toggle: J");
      ImGui::BulletText("Supercruise throttle: W/S (X = fast decel)");
      ImGui::BulletText("Boost: LShift   Brake: X");
      ImGui::BulletText("Dampers: Z (on) / C (off)");
      ImGui::BulletText("Time accel: PageUp/PageDown  (Home=1x, End=max, CTRL=override)");
      ImGui::BulletText("Pause: Space   Save: F5   Load: F9");

      ImGui::End();
    }

    if (showEconomy) {
      beginDockedHUD(*currentSystem, selectedStationIndex);

      if (!currentSystem->stations.empty()) {
        selectedStationIndex = std::max(0, std::min(selectedStationIndex, (int)currentSystem->stations.size() - 1));
        const auto& station = currentSystem->stations[(std::size_t)selectedStationIndex];
        auto& stEcon = universe.stationEconomy(station, timeDays);
        const bool atThisStation = docked && (dockedStationId == station.id);

        ImGui::Begin("Market Details");
        ImGui::Text("Credits: %.2f", credits);
        ImGui::SameLine();
        ImGui::TextDisabled("| %s", atThisStation ? "DOCKED" : "IN FLIGHT (dock to trade)");

        static int selectedCommodity = 0;
        ImGui::SliderInt("Plot commodity", &selectedCommodity, 0, (int)econ::kCommodityCount - 1);

        // Table
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

            static float qty[ (int)econ::kCommodityCount ] = {};
            if (qty[i] <= 0.0f) qty[i] = 10.0f;
            ImGui::SetNextItemWidth(70);
            ImGui::InputFloat("##qty", &qty[i], 1.0f, 10.0f, "%.0f");

            ImGui::SameLine();
            if (ImGui::SmallButton("Buy")) {
              if (!atThisStation) {
                setToast("You must dock at " + station.name + " to trade.");
              } else {
                auto tr = econ::buy(stEcon, station.economyModel, cid, qty[i], credits, 0.10, station.feeRate);
                if (tr.ok) cargo[i] += qty[i];
              }
            }

            ImGui::SameLine();
            if (ImGui::SmallButton("Sell")) {
              if (!atThisStation) {
                setToast("You must dock at " + station.name + " to trade.");
              } else {
                const double sellUnits = std::min<double>(qty[i], cargo[i]);
                if (sellUnits > 0.0) {
                  auto tr = econ::sell(stEcon, station.economyModel, cid, sellUnits, credits, 0.10, station.feeRate);
                  if (tr.ok) cargo[i] -= sellUnits;
                }
              }
            }

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

        ImGui::Separator();
        ImGui::Text("Route planner (profit/unit):");

        // Choose destination among nearby systems
        auto nearby = universe.queryNearby(currentSystem->stub.posLy, 200.0, 64);

        static int destIndex = 0;
        if (!nearby.empty()) {
          destIndex = std::min(destIndex, (int)nearby.size() - 1);

          std::vector<std::string> labels;
          labels.reserve(nearby.size());
          for (const auto& s : nearby) {
            const math::Vec3d d = s.posLy - currentSystem->stub.posLy;
            const double dist = std::sqrt(d.lengthSq());
            labels.push_back(s.name + " (" + std::to_string(dist) + " ly)");
          }
          std::vector<const char*> cstr;
          cstr.reserve(labels.size());
          for (auto& l : labels) cstr.push_back(l.c_str());

          ImGui::Combo("Destination system", &destIndex, cstr.data(), (int)cstr.size());

          const auto& dstStub = nearby[(std::size_t)destIndex];
          const auto& dstSys = universe.getSystem(dstStub.id, &dstStub);
          if (!dstSys.stations.empty()) {
            static int dstStationIndex = 0;
            dstStationIndex = std::min(dstStationIndex, (int)dstSys.stations.size() - 1);

            std::vector<const char*> stNames;
            stNames.reserve(dstSys.stations.size());
            for (const auto& st : dstSys.stations) stNames.push_back(st.name.c_str());
            ImGui::Combo("Destination station", &dstStationIndex, stNames.data(), (int)stNames.size());

            const auto& dstStation = dstSys.stations[(std::size_t)dstStationIndex];

            auto& dstEcon = universe.stationEconomy(dstStation, timeDays);

            const auto routes = econ::bestRoutes(stEcon, station.economyModel, dstEcon, dstStation.economyModel, 0.10, 5);
            if (routes.empty()) {
              ImGui::TextDisabled("No profitable routes right now.");
            } else {
              for (const auto& r : routes) {
                ImGui::BulletText("%s: +%.2f (buy %.2f → sell %.2f)",
                                  std::string(econ::commodityName(r.commodity)).c_str(),
                                  r.profitPerUnit, r.buyPrice, r.sellPrice);
              }
            }
          }
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

      // Background
      draw->AddRectFilled(p0, p1, IM_COL32(10, 10, 14, 255));
      draw->AddRect(p0, p1, IM_COL32(80, 80, 95, 255));

      auto toPx = [&](const math::Vec3d& posLy) -> ImVec2 {
        const math::Vec3d d = posLy - center;
        const float sx = (float)(d.x / (double)radius) * (canvasSize.x * 0.5f);
        const float sy = (float)(d.y / (double)radius) * (canvasSize.y * 0.5f);
        return ImVec2(centerPx.x + sx, centerPx.y + sy);
      };

      // Star lanes: connect each system to 3 nearest neighbors in XY
      const int k = 3;
      for (std::size_t i = 0; i < nearby.size(); ++i) {
        // find neighbors
        struct N { std::size_t j; double d2; };
        std::vector<N> ns;
        ns.reserve(nearby.size());
        for (std::size_t j = 0; j < nearby.size(); ++j) if (j != i) {
          const auto di = nearby[j].posLy - nearby[i].posLy;
          const double d2 = di.x*di.x + di.y*di.y + di.z*di.z;
          ns.push_back({j, d2});
        }
        std::sort(ns.begin(), ns.end(), [](const N& a, const N& b){ return a.d2 < b.d2; });
        const int count = std::min<int>(k, (int)ns.size());

        const ImVec2 a = toPx(nearby[i].posLy);
        for (int n = 0; n < count; ++n) {
          const ImVec2 b = toPx(nearby[ns[n].j].posLy);
          draw->AddLine(a, b, IM_COL32(50, 80, 120, 100), 1.0f);
        }
      }

      // Systems
      static sim::SystemId selected = 0;
      for (const auto& s : nearby) {
        const ImVec2 p = toPx(s.posLy);
        const bool isCurrent = (s.id == currentSystem->stub.id);
        const bool isSel = (s.id == selected);

        ImU32 col = isCurrent ? IM_COL32(255, 240, 160, 255) : IM_COL32(170, 170, 190, 255);
        if (s.factionId != 0) col = IM_COL32(160, 220, 170, 255);
        if (isSel) col = IM_COL32(255, 120, 120, 255);

        draw->AddCircleFilled(p, isCurrent ? 5.5f : 4.0f, col);

        // Click detection
        const float rClick = 6.0f;
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
          const ImVec2 mp = ImGui::GetIO().MousePos;
          const float dx = mp.x - p.x;
          const float dy = mp.y - p.y;
          if (dx*dx + dy*dy <= rClick*rClick) selected = s.id;
        }
      }

      ImGui::EndChild();

      // List + jump
      if (selected != 0 && selected != currentSystem->stub.id) {
        if (ImGui::Button("Jump to selected system")) {
          // Find stub in list
          auto it = std::find_if(nearby.begin(), nearby.end(), [&](const sim::SystemStub& s){ return s.id == selected; });
          if (it != nearby.end()) {
            currentStub = *it;
            const auto& sys = universe.getSystem(currentStub.id, &currentStub);
            currentSystem = &sys;

            // Reset flight state.
            docked = false;
            dockedStationId = 0;
            dockedStationIndex = 0;
            selectedStationIndex = 0;
            supercruise = false;
            supercruiseThrottle = 0.0;
            supercruiseSpeedKmS = 0.0;

            ship.setVelocityKmS({0,0,0});
            ship.setAngularVelocityRadS({0,0,0});

            // Spawn near the first station if available, otherwise near origin.
            if (!currentSystem->stations.empty()) {
              const auto& st0 = currentSystem->stations[0];
              const math::Vec3d pKm = stationPosKm(st0, timeDays);
              const math::Vec3d outDir = stationOutDir(pKm);
              ship.setPositionKm(pKm + outDir * (st0.radiusKm + 25.0));
            } else {
              ship.setPositionKm({0, 0, -8000.0});
            }

            setToast("Jumped to " + currentSystem->stub.name);
          }
        }
      }

      ImGui::TextDisabled("Tip: TAB toggles this window, F1 Flight, F2 Economy");

      ImGui::End();
    }

    // Toast / HUD message (top-left)
    if (toastTimer > 0.0 && !toast.empty()) {
      ImGui::SetNextWindowPos(ImVec2(20, 20), ImGuiCond_Always);
      ImGui::SetNextWindowBgAlpha(0.35f);
      ImGui::Begin("##toast", nullptr,
                   ImGuiWindowFlags_NoDecoration |
                   ImGuiWindowFlags_AlwaysAutoResize |
                   ImGuiWindowFlags_NoMove |
                   ImGuiWindowFlags_NoSavedSettings |
                   ImGuiWindowFlags_NoFocusOnAppearing |
                   ImGuiWindowFlags_NoNav);
      ImGui::TextUnformatted(toast.c_str());
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
