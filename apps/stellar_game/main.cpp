#include "stellar/core/Log.h"
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

static const char* planetTypeName(sim::PlanetType t) {
  switch (t) {
    case sim::PlanetType::Rocky: return "Rocky";
    case sim::PlanetType::Desert: return "Desert";
    case sim::PlanetType::Ocean: return "Ocean";
    case sim::PlanetType::Ice: return "Ice";
    case sim::PlanetType::GasGiant: return "Gas Giant";
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

static void stationColorRgb(econ::StationType t, float& r, float& g, float& b) {
  // Visually distinct station colors.
  switch (t) {
    case econ::StationType::Outpost:       r=0.70f; g=0.70f; b=0.75f; break;
    case econ::StationType::Agricultural: r=0.35f; g=0.85f; b=0.45f; break;
    case econ::StationType::Mining:       r=0.85f; g=0.55f; b=0.25f; break;
    case econ::StationType::Refinery:     r=0.90f; g=0.85f; b=0.30f; break;
    case econ::StationType::Industrial:   r=0.35f; g=0.55f; b=0.95f; break;
    case econ::StationType::Research:     r=0.80f; g=0.40f; b=0.95f; break;
    case econ::StationType::TradeHub:     r=0.25f; g=0.85f; b=0.95f; break;
    case econ::StationType::Shipyard:     r=0.95f; g=0.35f; b=0.35f; break;
    default:                              r=0.80f; g=0.80f; b=0.85f; break;
  }
}

static math::Vec3d planetPositionKm(const sim::Planet& p, double timeDays) {
  const math::Vec3d posAU = sim::orbitPosition3DAU(p.orbit, timeDays);
  return posAU * kAU_KM;
}

static math::Vec3d stationPositionKm(const sim::Station& st, double timeDays) {
  const math::Vec3d posAU = sim::orbitPosition3DAU(st.orbit, timeDays);
  return posAU * kAU_KM;
}

static math::Vec3d clampMagnitude(const math::Vec3d& v, double maxLen) {
  const double len = v.length();
  if (len <= maxLen || len <= 1e-12) return v;
  return v * (maxLen / len);
}

enum class TargetKind {
  None,
  Station,
  Planet,
};

struct TargetRef {
  TargetKind kind{TargetKind::None};
  int index{-1};
};

static bool isValidTarget(const TargetRef& t, const sim::StarSystem& sys) {
  if (t.kind == TargetKind::Station) return t.index >= 0 && t.index < (int)sys.stations.size();
  if (t.kind == TargetKind::Planet)  return t.index >= 0 && t.index < (int)sys.planets.size();
  return false;
}

static const char* targetKindName(TargetKind k) {
  switch (k) {
    case TargetKind::Station: return "Station";
    case TargetKind::Planet:  return "Planet";
    default:                  return "None";
  }
}

enum class CameraMode {
  Chase,
  Cockpit,
};

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

  render::PointRenderer pointRenderer;
  if (!pointRenderer.init(&err)) {
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

  // --- Time ---
  double timeDays = 0.0;
  double timeScale = 60.0; // simulated seconds per real second
  bool paused = false;

  // --- Player ship ---
  sim::Ship ship;
  ship.setMaxLinearAccelKmS2(0.08);
  ship.setMaxAngularAccelRadS2(1.2);

  // Docking state (gameplay)
  bool docked = false;
  int dockedStationIndex = 0;
  int viewStationIndex = 0; // which station market is shown for

  // Targeting / autopilot
  TargetRef target{};
  bool autopilotEnabled = false;

  // Camera
  CameraMode camMode = CameraMode::Chase;
  bool mouseFlight = false;
  float mouseSens = 0.0035f;
  double chaseDistanceU = 7.5;
  double chaseHeightU = 2.8;

  // Spawn: start docked at first station if possible (makes the build immediately playable).
  if (!currentSystem->stations.empty()) {
    docked = true;
    dockedStationIndex = 0;
    viewStationIndex = 0;
    const math::Vec3d stPosKm = stationPositionKm(currentSystem->stations[0], timeDays);
    ship.setPositionKm(stPosKm + math::Vec3d{0, 0, -(currentSystem->stations[0].radiusKm + 40.0)});
    ship.setVelocityKmS({0,0,0});
    ship.setAngularVelocityRadS({0,0,0});
  } else {
    ship.setPositionKm({0, 0, -8000.0}); // 8000 km behind origin
  }

  // Player economy
  double credits = 2500.0;
  std::array<double, econ::kCommodityCount> cargo{};

  // Save/load
  const std::string savePath = "savegame.txt";

  // UI state
  bool showGalaxy = true;
  bool showShip = true;
  bool showEconomy = true;
  bool showSystem = true;
  bool showHud = true;

  // Galaxy nav selection
  sim::SystemId galaxySelectedSystem = currentSystem->stub.id;
  int galaxySelectedStationIndex = 0;

  // Starfield (generated once; rendered with view translation removed for an "infinite" feel).
  std::vector<render::PointVertex> starfield;
  {
    starfield.reserve(2000);
    core::SplitMix64 rng(seed ^ 0xA1B2C3D4u);
    for (int i = 0; i < 2000; ++i) {
      // random point on sphere
      const double z = rng.range(-1.0, 1.0);
      const double a = rng.range(0.0, 2.0 * math::kPi);
      const double r = std::sqrt(std::max(0.0, 1.0 - z*z));
      const double radius = rng.range(3500.0, 9000.0); // in render units

      const float px = (float)(std::cos(a) * r * radius);
      const float py = (float)(z * radius);
      const float pz = (float)(std::sin(a) * r * radius);

      const float tint = (float)rng.range(0.85, 1.0);
      const float warm = (float)rng.range(0.0, 1.0);
      const float cr = tint * (0.80f + 0.20f*warm);
      const float cg = tint * (0.80f + 0.10f*(1.0f-warm));
      const float cb = tint * (0.90f + 0.10f*(1.0f-warm));
      const float size = (float)rng.range(1.0, 2.2);

      starfield.push_back({px,py,pz, cr,cg,cb, size});
    }
  }

  bool running = true;
  auto last = std::chrono::high_resolution_clock::now();

  SDL_SetRelativeMouseMode(SDL_FALSE);

  bool requestDockToggle = false;

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
          s.dockedStation = (docked && !currentSystem->stations.empty())
            ? currentSystem->stations[(std::size_t)dockedStationIndex].id
            : 0;

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

            // Restore docking
            docked = false;
            dockedStationIndex = 0;
            if (s.dockedStation != 0 && !sys.stations.empty()) {
              for (std::size_t i = 0; i < sys.stations.size(); ++i) {
                if (sys.stations[i].id == s.dockedStation) {
                  docked = true;
                  dockedStationIndex = (int)i;
                  break;
                }
              }
            }

            if (!sys.stations.empty()) {
              viewStationIndex = std::clamp(viewStationIndex, 0, (int)sys.stations.size() - 1);
              if (docked) viewStationIndex = dockedStationIndex;
            } else {
              viewStationIndex = 0;
            }

            // If docked, snap to station position.
            if (docked && !sys.stations.empty()) {
              const auto& st = sys.stations[(std::size_t)dockedStationIndex];
              const math::Vec3d stPosKm = stationPositionKm(st, timeDays);
              ship.setPositionKm(stPosKm + math::Vec3d{0, 0, -(st.radiusKm + 40.0)});
              ship.setVelocityKmS({0,0,0});
              ship.setAngularVelocityRadS({0,0,0});
            }

            // Reset navigation state
            autopilotEnabled = false;
            target = TargetRef{};
            galaxySelectedSystem = currentSystem->stub.id;
            galaxySelectedStationIndex = 0;

            core::log(core::LogLevel::Info, "Loaded " + savePath);
          }
        }

        if (event.key.keysym.sym == SDLK_TAB) showGalaxy = !showGalaxy;
        if (event.key.keysym.sym == SDLK_F1) showShip = !showShip;
        if (event.key.keysym.sym == SDLK_F2) showEconomy = !showEconomy;
        if (event.key.keysym.sym == SDLK_F3) showSystem = !showSystem;
        if (event.key.keysym.sym == SDLK_F4) showHud = !showHud;
        if (event.key.keysym.sym == SDLK_SPACE) paused = !paused;

        if (!captureKeys) {
          if (event.key.keysym.sym == SDLK_g) {
            requestDockToggle = true;
          }

          if (event.key.keysym.sym == SDLK_t) {
            // Cycle target through stations, then planets.
            const int nStations = (int)currentSystem->stations.size();
            const int nPlanets  = (int)currentSystem->planets.size();
            const int total = nStations + nPlanets;
            if (total == 0) {
              target = TargetRef{};
            } else {
              int cursor = -1;
              if (target.kind == TargetKind::Station) cursor = target.index;
              if (target.kind == TargetKind::Planet)  cursor = nStations + target.index;

              cursor = (cursor + 1) % total;
              if (cursor < nStations) {
                target.kind = TargetKind::Station;
                target.index = cursor;
              } else {
                target.kind = TargetKind::Planet;
                target.index = cursor - nStations;
              }
            }
          }

          if (event.key.keysym.sym == SDLK_p) {
            autopilotEnabled = !autopilotEnabled;
          }

          if (event.key.keysym.sym == SDLK_v) {
            camMode = (camMode == CameraMode::Chase) ? CameraMode::Cockpit : CameraMode::Chase;
          }

          if (event.key.keysym.sym == SDLK_m) {
            mouseFlight = !mouseFlight;
            SDL_SetRelativeMouseMode(mouseFlight ? SDL_TRUE : SDL_FALSE);
          }

          if (event.key.keysym.sym == SDLK_LEFTBRACKET) {
            timeScale = std::max(0.0, timeScale * 0.5);
          }
          if (event.key.keysym.sym == SDLK_RIGHTBRACKET) {
            timeScale = std::min(50000.0, timeScale * 2.0);
          }
        }
      }
    }

    // Build per-frame planet/station positions for gameplay + UI
    std::vector<math::Vec3d> planetPosKm(currentSystem->planets.size());
    for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
      planetPosKm[i] = planetPositionKm(currentSystem->planets[i], timeDays);
    }

    std::vector<math::Vec3d> stationPosKm(currentSystem->stations.size());
    for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
      stationPosKm[i] = stationPositionKm(currentSystem->stations[i], timeDays);
    }

    // Ship speed for control / autopilot logic
    double shipSpeed = ship.velocityKmS().length();

    auto dockingRangeKmFor = [&](int stIdx) -> double {
      if (stIdx < 0 || stIdx >= (int)currentSystem->stations.size()) return 0.0;
      const double r = currentSystem->stations[(std::size_t)stIdx].radiusKm;
      return std::max(120.0, r * 35.0);
    };

    // Input (6DOF)
    sim::ShipInput input{};
    const Uint8* keys = SDL_GetKeyboardState(nullptr);

    const bool captureKeys = io.WantCaptureKeyboard;

    if (!captureKeys && !docked) {
      // Dampers toggle: hold Z or C
      static bool dampers = true;
      if (keys[SDL_SCANCODE_Z]) dampers = true;
      if (keys[SDL_SCANCODE_C]) dampers = false;
      input.dampers = dampers;

      // Manual control only if autopilot is off
      if (!autopilotEnabled || !isValidTarget(target, *currentSystem)) {
        // Translate
        input.thrustLocal.z += (keys[SDL_SCANCODE_W] ? 1.0 : 0.0);
        input.thrustLocal.z -= (keys[SDL_SCANCODE_S] ? 1.0 : 0.0);

        input.thrustLocal.x += (keys[SDL_SCANCODE_D] ? 1.0 : 0.0);
        input.thrustLocal.x -= (keys[SDL_SCANCODE_A] ? 1.0 : 0.0);

        input.thrustLocal.y += (keys[SDL_SCANCODE_R] ? 1.0 : 0.0);
        input.thrustLocal.y -= (keys[SDL_SCANCODE_F] ? 1.0 : 0.0);

        // Rotate (keyboard)
        input.torqueLocal.x += (keys[SDL_SCANCODE_UP] ? 1.0 : 0.0);
        input.torqueLocal.x -= (keys[SDL_SCANCODE_DOWN] ? 1.0 : 0.0);

        input.torqueLocal.y += (keys[SDL_SCANCODE_RIGHT] ? 1.0 : 0.0);
        input.torqueLocal.y -= (keys[SDL_SCANCODE_LEFT] ? 1.0 : 0.0);

        input.torqueLocal.z += (keys[SDL_SCANCODE_E] ? 1.0 : 0.0);
        input.torqueLocal.z -= (keys[SDL_SCANCODE_Q] ? 1.0 : 0.0);

        // Mouse flight (optional)
        if (mouseFlight && !io.WantCaptureMouse) {
          int mx = 0, my = 0;
          SDL_GetRelativeMouseState(&mx, &my);
          input.torqueLocal.y += (double)mx * (double)mouseSens; // yaw
          input.torqueLocal.x += (double)(-my) * (double)mouseSens; // pitch
        }

        input.boost = keys[SDL_SCANCODE_LSHIFT] != 0;
        input.brake = keys[SDL_SCANCODE_X] != 0;
      } else {
        // --- Autopilot ---
        const math::Vec3d shipPos = ship.positionKm();
        const math::Vec3d shipVel = ship.velocityKmS();

        math::Vec3d tgtPos{0,0,0};
        if (target.kind == TargetKind::Station) {
          tgtPos = stationPosKm[(std::size_t)target.index];
        } else {
          tgtPos = planetPosKm[(std::size_t)target.index];
        }

        const math::Vec3d to = tgtPos - shipPos;
        const double dist = to.length();
        const math::Vec3d dir = (dist > 1e-6) ? (to * (1.0 / dist)) : math::Vec3d{0,0,1};

        // Desired approach speed decreases as we get closer.
        const double vMax = (target.kind == TargetKind::Planet) ? 40.0 : 14.0; // km/s
        double vDesired = std::clamp(dist * 0.02, 0.0, vMax);
        if (target.kind == TargetKind::Station) {
          vDesired = std::min(vDesired, std::max(1.0, dist * 0.01));
        }

        const math::Vec3d velDesired = dir * vDesired;
        math::Vec3d accelCmd = (velDesired - shipVel) * 0.75; // ~1/s

        const double linCap = ship.maxLinearAccelKmS2();
        const bool boost = dist > 25000.0;
        const double cap = linCap * (boost ? 1.8 : 1.0);
        accelCmd = clampMagnitude(accelCmd, cap);

        // Convert world accel to local thrust command.
        const math::Vec3d accelLocal = ship.orientation().conjugate().rotate(accelCmd);
        math::Vec3d thr = accelLocal * (1.0 / std::max(1e-9, linCap));
        thr.x = std::clamp(thr.x, -1.0, 1.0);
        thr.y = std::clamp(thr.y, -1.0, 1.0);
        thr.z = std::clamp(thr.z, -1.0, 1.0);
        input.thrustLocal = thr;

        // Gentle orientation assist: point the ship forward towards the target.
        {
          const math::Vec3d fwd = ship.forward().normalized();
          const math::Vec3d axisW = cross(fwd, dir);
          const double s = axisW.length();
          if (s > 1e-6) {
            const double ang = std::asin(std::clamp(s, 0.0, 1.0));
            const math::Vec3d axisL = ship.orientation().conjugate().rotate(axisW.normalized());
            math::Vec3d torque = axisL * std::clamp(ang * 1.2, 0.0, 1.0);
            torque.x = std::clamp(torque.x, -1.0, 1.0);
            torque.y = std::clamp(torque.y, -1.0, 1.0);
            torque.z = std::clamp(torque.z, -1.0, 1.0);
            input.torqueLocal = torque;
          }
        }

        input.dampers = true;
        input.boost = boost;
        input.brake = (dist < 500.0) && (shipSpeed > 2.5);

        // Auto-dock when we're close enough.
        const double autoDockRangeKm = dockingRangeKmFor(target.index);
        if (target.kind == TargetKind::Station && dist <= autoDockRangeKm && shipSpeed <= 1.0) {
          requestDockToggle = true;
        }
      }
    }

    // Step sim
    if (!paused) {
      if (!docked) {
        ship.step(dt, input);
      }
      timeDays += (dt * timeScale) / 86400.0;
    }

    // Recompute station/planet positions after time update (time affects orbit positions).
    for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
      planetPosKm[i] = planetPositionKm(currentSystem->planets[i], timeDays);
    }
    for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
      stationPosKm[i] = stationPositionKm(currentSystem->stations[i], timeDays);
    }

    // Nearest station (for docking prompt / manual docking)
    int nearestStation = -1;
    double nearestStationDistKm = 1e300;
    for (std::size_t i = 0; i < stationPosKm.size(); ++i) {
      const double d = (stationPosKm[i] - ship.positionKm()).length();
      if (d < nearestStationDistKm) {
        nearestStationDistKm = d;
        nearestStation = (int)i;
      }
    }

    // Update ship speed post-sim-step for UI / docking checks.
    shipSpeed = ship.velocityKmS().length();

    const double dockRangeKm = dockingRangeKmFor(nearestStation);
    const bool canDockNow = (!docked)
      && (nearestStation >= 0)
      && (nearestStationDistKm <= dockRangeKm)
      && (shipSpeed <= 2.0)
      && (!paused);

    // Handle docking request (G or autopilot).
    if (requestDockToggle && !captureKeys) {
      if (docked) {
        // Undock
        if (!currentSystem->stations.empty()) {
          const auto& st = currentSystem->stations[(std::size_t)dockedStationIndex];
          const math::Vec3d stPos = stationPosKm[(std::size_t)dockedStationIndex];
          math::Vec3d dir = stPos.normalized();
          if (dir.length() < 1e-6) dir = {0,0,1};

          docked = false;
          ship.setPositionKm(stPos + dir * (st.radiusKm + 120.0));
          ship.setVelocityKmS({0,0,0});
          ship.setAngularVelocityRadS({0,0,0});
        } else {
          docked = false;
        }
      } else {
        // Dock
        if (canDockNow && nearestStation >= 0) {
          docked = true;
          dockedStationIndex = nearestStation;
          viewStationIndex = dockedStationIndex;

          const auto& st = currentSystem->stations[(std::size_t)dockedStationIndex];
          const math::Vec3d stPos = stationPosKm[(std::size_t)dockedStationIndex];
          ship.setPositionKm(stPos + math::Vec3d{0, 0, -(st.radiusKm + 40.0)});
          ship.setVelocityKmS({0,0,0});
          ship.setAngularVelocityRadS({0,0,0});
        }
      }

      requestDockToggle = false;
    }

    // If docked, keep ship "attached" to the station (station moves on its orbit).
    if (docked && !currentSystem->stations.empty()) {
      dockedStationIndex = std::clamp(dockedStationIndex, 0, (int)currentSystem->stations.size() - 1);
      const auto& st = currentSystem->stations[(std::size_t)dockedStationIndex];
      const math::Vec3d stPos = stationPosKm[(std::size_t)dockedStationIndex];
      ship.setPositionKm(stPos + math::Vec3d{0, 0, -(st.radiusKm + 40.0)});
      ship.setVelocityKmS({0,0,0});
      ship.setAngularVelocityRadS({0,0,0});
      autopilotEnabled = false;
    }

    // Soft collision with star/planets (prevents drifting through bodies).
    {
      const math::Vec3d shipPos = ship.positionKm();
      const double starRadiusKm = currentSystem->star.radiusSol * kSOLAR_RADIUS_KM;
      const double dStar = shipPos.length();
      if (dStar < starRadiusKm * 1.02) {
        const math::Vec3d n = (dStar > 1e-9) ? (shipPos * (1.0 / dStar)) : math::Vec3d{0,0,1};
        ship.setPositionKm(n * (starRadiusKm * 1.05));
        ship.setVelocityKmS({0,0,0});
      }

      for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
        const auto& p = currentSystem->planets[i];
        const double rKm = p.radiusEarth * kEARTH_RADIUS_KM;
        const math::Vec3d dp = shipPos - planetPosKm[i];
        const double d = dp.length();
        if (d < rKm * 1.02) {
          const math::Vec3d n = (d > 1e-9) ? (dp * (1.0 / d)) : math::Vec3d{0,0,1};
          ship.setPositionKm(planetPosKm[i] + n * (rKm * 1.06));
          ship.setVelocityKmS({0,0,0});
        }
      }
    }

    // ---- Camera ----
    int w = 1280, h = 720;
    SDL_GetWindowSize(window, &w, &h);
    const double aspect = (h > 0) ? (double)w / (double)h : 16.0/9.0;

    render::Camera cam;
    cam.setPerspective(math::degToRad(60.0), aspect, 0.01, 20000.0);

    const math::Vec3d shipPosU = ship.positionKm() * (1.0 / kRENDER_UNIT_KM);

    math::Mat4d view{};
    if (camMode == CameraMode::Chase) {
      const math::Vec3d eye = shipPosU + ship.forward() * (-chaseDistanceU) + math::Vec3d{0,1,0} * chaseHeightU;
      view = math::Mat4d::lookAt(eye, shipPosU, {0,1,0});
    } else {
      const math::Vec3d eye = shipPosU + ship.forward() * 0.10 + ship.up() * 0.02;
      view = math::Mat4d::lookAt(eye, eye + ship.forward(), ship.up());
    }

    const math::Mat4d proj = cam.projectionMatrix();

    float viewF[16], projF[16];
    matToFloat(view, viewF);
    matToFloat(proj, projF);

    meshRenderer.setViewProj(viewF, projF);
    lineRenderer.setViewProj(viewF, projF);

    // Starfield uses a view matrix without translation.
    float viewNoTrans[16];
    std::memcpy(viewNoTrans, viewF, sizeof(viewF));
    viewNoTrans[12] = viewNoTrans[13] = viewNoTrans[14] = 0.0f;
    pointRenderer.setViewProj(viewNoTrans, projF);

    // ---- Build instances (star + planets) ----
    std::vector<render::InstanceData> spheres;
    spheres.reserve(1 + currentSystem->planets.size());

    // Star at origin
    {
      const double starRadiusKm = currentSystem->star.radiusSol * kSOLAR_RADIUS_KM;
      const float starScale = (float)std::max(0.8, (starRadiusKm / kRENDER_UNIT_KM) * 3.0);
      spheres.push_back({0,0,0, starScale, 1,0,0,0, 1.0f, 0.95f, 0.75f});
    }

    // Planets
    for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
      const auto& p = currentSystem->planets[i];
      const math::Vec3d posU = planetPosKm[i] * (1.0 / kRENDER_UNIT_KM);

      const double radiusKm = p.radiusEarth * kEARTH_RADIUS_KM;
      const float scale = (float)std::max(0.25, (radiusKm / kRENDER_UNIT_KM) * 200.0);

      // Simple color palette by type
      float cr=0.6f, cg=0.6f, cb=0.6f;
      switch (p.type) {
        case sim::PlanetType::Rocky:   cr=0.6f; cg=0.55f; cb=0.5f; break;
        case sim::PlanetType::Desert:  cr=0.8f; cg=0.7f;  cb=0.35f; break;
        case sim::PlanetType::Ocean:   cr=0.25f;cg=0.45f; cb=0.85f; break;
        case sim::PlanetType::Ice:     cr=0.7f; cg=0.85f; cb=0.95f; break;
        case sim::PlanetType::GasGiant:cr=0.7f; cg=0.55f; cb=0.35f; break;
        default: break;
      }

      spheres.push_back({(float)posU.x, (float)posU.y, (float)posU.z,
                         scale,
                         1,0,0,0,
                         cr,cg,cb});
    }

    // Orbit lines
    std::vector<render::LineVertex> lines;
    lines.reserve(currentSystem->planets.size() * 128 + 8);

    for (const auto& p : currentSystem->planets) {
      const int seg = 96;
      math::Vec3d prev{};
      for (int s = 0; s <= seg; ++s) {
        const double t = (double)s / (double)seg * p.orbit.periodDays;
        const math::Vec3d posAU = sim::orbitPosition3DAU(p.orbit, t);
        const math::Vec3d posU = (posAU * kAU_KM) * (1.0 / kRENDER_UNIT_KM);

        if (s > 0) {
          lines.push_back({(float)prev.x,(float)prev.y,(float)prev.z, 0.25f,0.25f,0.25f});
          lines.push_back({(float)posU.x,(float)posU.y,(float)posU.z, 0.25f,0.25f,0.25f});
        }
        prev = posU;
      }
    }

    // Target line (ship -> target)
    if (isValidTarget(target, *currentSystem)) {
      math::Vec3d tgtU{0,0,0};
      if (target.kind == TargetKind::Station) tgtU = stationPosKm[(std::size_t)target.index] * (1.0 / kRENDER_UNIT_KM);
      if (target.kind == TargetKind::Planet)  tgtU = planetPosKm[(std::size_t)target.index] * (1.0 / kRENDER_UNIT_KM);

      lines.push_back({(float)shipPosU.x, (float)shipPosU.y, (float)shipPosU.z, 0.2f, 0.95f, 0.2f});
      lines.push_back({(float)tgtU.x, (float)tgtU.y, (float)tgtU.z, 0.2f, 0.95f, 0.2f});
    }

    // Stations + ship
    std::vector<render::InstanceData> cubes;
    cubes.reserve(currentSystem->stations.size() + 1);

    for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
      const auto& st = currentSystem->stations[i];
      const math::Vec3d posU = stationPosKm[i] * (1.0 / kRENDER_UNIT_KM);

      float cr, cg, cb;
      stationColorRgb(st.type, cr, cg, cb);

      if ((int)i == dockedStationIndex && docked) {
        // highlight docked station
        cr = std::min(1.0f, cr * 1.2f);
        cg = std::min(1.0f, cg * 1.2f);
        cb = std::min(1.0f, cb * 1.2f);
      }

      const float base = 0.25f;
      const float size = base + (float)std::clamp(st.radiusKm / 30.0, 0.0, 1.0) * 0.35f;

      cubes.push_back({(float)posU.x, (float)posU.y, (float)posU.z,
                       size,
                       1,0,0,0,
                       cr,cg,cb});
    }

    // Ship cube instance
    {
      const auto q = ship.orientation().normalized();
      cubes.push_back({
        (float)shipPosU.x, (float)shipPosU.y, (float)shipPosU.z,
        0.35f,
        (float)q.w, (float)q.x, (float)q.y, (float)q.z,
        0.90f, 0.90f, 1.0f
      });
    }

    // ---- Render ---
    glViewport(0, 0, w, h);
    glClearColor(0.01f, 0.01f, 0.02f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Starfield
    glDepthMask(GL_FALSE);
    pointRenderer.drawPoints(starfield);
    glDepthMask(GL_TRUE);

    // Orbits + target line
    lineRenderer.drawLines(lines);

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

    // Main dockspace
    ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());

    // HUD overlay
    if (showHud) {
      ImGuiWindowFlags hudFlags = ImGuiWindowFlags_NoDecoration
        | ImGuiWindowFlags_AlwaysAutoResize
        | ImGuiWindowFlags_NoSavedSettings
        | ImGuiWindowFlags_NoFocusOnAppearing
        | ImGuiWindowFlags_NoNav
        | ImGuiWindowFlags_NoMove;

      ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_Always);
      ImGui::Begin("HUD", nullptr, hudFlags);

      ImGui::Text("%s  (%s)", currentSystem->stub.name.c_str(), starClassName(currentSystem->stub.primaryClass));
      ImGui::Text("Time: %.2f d  x%.1f %s", timeDays, timeScale, paused ? "[PAUSED]" : "");

      const auto pos = ship.positionKm();
      const auto vel = ship.velocityKmS();
      ImGui::Text("Speed: %.2f km/s", vel.length());

      ImGui::Separator();

      if (docked && !currentSystem->stations.empty()) {
        const auto& st = currentSystem->stations[(std::size_t)dockedStationIndex];
        ImGui::TextColored(ImVec4(0.7f,1.0f,0.7f,1.0f), "Docked: %s", st.name.c_str());
      } else {
        ImGui::Text("Docked: no");
        if (nearestStation >= 0) {
          const auto& st = currentSystem->stations[(std::size_t)nearestStation];
          ImGui::Text("Nearest: %s (%.0f km)", st.name.c_str(), nearestStationDistKm);
          if (canDockNow) {
            ImGui::TextColored(ImVec4(1.0f,0.9f,0.4f,1.0f), "Press G to dock");
          }
        }
      }

      if (isValidTarget(target, *currentSystem)) {
        math::Vec3d tgtPos{};
        std::string name;
        if (target.kind == TargetKind::Station) {
          const auto& st = currentSystem->stations[(std::size_t)target.index];
          name = st.name;
          tgtPos = stationPosKm[(std::size_t)target.index];
        } else {
          const auto& p = currentSystem->planets[(std::size_t)target.index];
          name = p.name;
          tgtPos = planetPosKm[(std::size_t)target.index];
        }
        const double dist = (tgtPos - ship.positionKm()).length();

        ImGui::Separator();
        ImGui::Text("Target: %s", targetKindName(target.kind));
        ImGui::Text("  %s", name.c_str());
        ImGui::Text("  Dist: %.0f km", dist);
      } else {
        ImGui::Separator();
        ImGui::TextDisabled("Target: none (press T)");
      }

      ImGui::Text("Autopilot: %s (P)", (autopilotEnabled && isValidTarget(target, *currentSystem) && !docked) ? "ON" : "OFF");

      const double fuelUnits = cargo[(std::size_t)econ::CommodityId::Fuel];
      ImGui::Text("Credits: %.0f   Fuel: %.0f", credits, fuelUnits);

      ImGui::End();

      // Crosshair
      ImDrawList* fg = ImGui::GetForegroundDrawList();
      const ImVec2 c = ImGui::GetMainViewport()->GetCenter();
      const float s = 8.0f;
      fg->AddLine(ImVec2(c.x - s, c.y), ImVec2(c.x + s, c.y), IM_COL32(200, 200, 220, 160), 1.0f);
      fg->AddLine(ImVec2(c.x, c.y - s), ImVec2(c.x, c.y + s), IM_COL32(200, 200, 220, 160), 1.0f);
    }

    if (showShip) {
      ImGui::Begin("Ship / Flight");

      const auto pos = ship.positionKm();
      const auto vel = ship.velocityKmS();
      const auto wv  = ship.angularVelocityRadS();

      ImGui::Text("Pos (km):   [%.1f %.1f %.1f]", pos.x, pos.y, pos.z);
      ImGui::Text("Vel (km/s): [%.3f %.3f %.3f] |v|=%.3f", vel.x, vel.y, vel.z, vel.length());
      ImGui::Text("AngVel (rad/s): [%.3f %.3f %.3f]", wv.x, wv.y, wv.z);

      ImGui::Separator();
      ImGui::SliderFloat("Time scale (sim sec / real sec)", (float*)&timeScale, 0.0f, 50000.0f);

      ImGui::Separator();
      ImGui::TextDisabled("Controls:");
      ImGui::BulletText("Translate: WASD + R/F");
      ImGui::BulletText("Rotate: Arrow keys + Q/E roll");
      ImGui::BulletText("Boost: LShift   Brake: X");
      ImGui::BulletText("Dampers: Z (on) / C (off)");
      ImGui::BulletText("Dock/Undock: G (near station)");
      ImGui::BulletText("Target cycle: T   Autopilot: P");
      ImGui::BulletText("Camera: V (chase/cockpit)");
      ImGui::BulletText("Mouse flight: M (toggle)");
      ImGui::BulletText("Save: F5   Load: F9");
      ImGui::BulletText("Toggle windows: TAB=Galaxy, F1=Flight, F2=Economy, F3=System, F4=HUD");

      ImGui::End();
    }

    if (showSystem) {
      ImGui::Begin("System / Nav");

      ImGui::Text("System: %s", currentSystem->stub.name.c_str());
      ImGui::Text("Star class: %s  Planets: %d  Stations: %d", starClassName(currentSystem->stub.primaryClass),
                  (int)currentSystem->planets.size(), (int)currentSystem->stations.size());

      if (!currentSystem->stations.empty()) {
        ImGui::Separator();
        ImGui::Text("Stations:");

        for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
          const auto& st = currentSystem->stations[i];
          const double dist = (stationPosKm[i] - ship.positionKm()).length();

          ImGui::PushID((int)i);
          const bool isTgt = (target.kind == TargetKind::Station && target.index == (int)i);
          if (ImGui::Selectable((st.name + std::string(" (") + stationTypeName(st.type) + ")").c_str(), isTgt)) {
            target.kind = TargetKind::Station;
            target.index = (int)i;
          }
          ImGui::SameLine();
          ImGui::TextDisabled("%.0f km", dist);
          ImGui::PopID();
        }
      }

      if (!currentSystem->planets.empty()) {
        ImGui::Separator();
        ImGui::Text("Planets:");

        for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
          const auto& p = currentSystem->planets[i];
          const double dist = (planetPosKm[i] - ship.positionKm()).length();

          ImGui::PushID((int)(1000 + i));
          const bool isTgt = (target.kind == TargetKind::Planet && target.index == (int)i);
          if (ImGui::Selectable((p.name + std::string(" (") + planetTypeName(p.type) + ")").c_str(), isTgt)) {
            target.kind = TargetKind::Planet;
            target.index = (int)i;
          }
          ImGui::SameLine();
          ImGui::TextDisabled("%.0f km", dist);
          ImGui::PopID();
        }
      }

      ImGui::End();
    }

    if (showEconomy) {
      ImGui::Begin("Docking / Market");

      ImGui::Text("System: %s", currentSystem->stub.name.c_str());

      if (docked && !currentSystem->stations.empty()) {
        const auto& st = currentSystem->stations[(std::size_t)dockedStationIndex];
        ImGui::TextColored(ImVec4(0.7f,1.0f,0.7f,1.0f), "Docked at: %s", st.name.c_str());
        ImGui::SameLine();
        if (ImGui::Button("Undock")) requestDockToggle = true;
      } else {
        ImGui::TextDisabled("Not docked");
        if (nearestStation >= 0) {
          const auto& st = currentSystem->stations[(std::size_t)nearestStation];
          ImGui::Text("Nearest: %s (%.0f km)", st.name.c_str(), nearestStationDistKm);
          ImGui::Text("Dock window: %.0f km  Speed: %.2f km/s", dockRangeKm, shipSpeed);
          if (canDockNow) {
            ImGui::TextColored(ImVec4(1.0f,0.9f,0.4f,1.0f), "Press G to dock");
          } else {
            ImGui::TextDisabled("To dock: get closer and reduce speed below 2 km/s");
          }
        } else {
          ImGui::TextDisabled("No stations in system.");
        }
      }

      ImGui::Separator();

      if (!currentSystem->stations.empty()) {
        viewStationIndex = std::clamp(viewStationIndex, 0, (int)currentSystem->stations.size() - 1);

        std::vector<const char*> names;
        names.reserve(currentSystem->stations.size());
        for (const auto& st : currentSystem->stations) names.push_back(st.name.c_str());

        ImGui::Combo("View station", &viewStationIndex, names.data(), (int)names.size());
        const auto& st = currentSystem->stations[(std::size_t)viewStationIndex];
        ImGui::SameLine();
        ImGui::TextDisabled("(%s, fee %.1f%%)", stationTypeName(st.type), st.feeRate * 100.0);
      }

      ImGui::End();

      if (!currentSystem->stations.empty()) {
        const bool canTrade = docked && (viewStationIndex == dockedStationIndex);
        const auto& station = currentSystem->stations[(std::size_t)viewStationIndex];
        auto& stEcon = universe.stationEconomy(station, timeDays);

        ImGui::Begin("Market Details");
        ImGui::Text("Credits: %.2f", credits);

        if (!canTrade) {
          ImGui::TextDisabled("(To trade here: dock at this station)");
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

            ImGui::SameLine();
            ImGui::BeginDisabled(!canTrade);
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
            ImGui::EndDisabled();

            ImGui::PopID();
          }

          ImGui::EndTable();
        }

        const std::size_t cidx = (std::size_t)selectedCommodity;
        const auto& hist = stEcon.history[cidx];
        if (!hist.empty()) {
          std::vector<float> vals;
          vals.reserve(hist.size());
          for (const auto& p : hist) vals.push_back((float)p.price);

          ImGui::PlotLines("Price history", vals.data(), (int)vals.size(), 0, nullptr, 0.0f, 0.0f, ImVec2(0, 120));
        } else {
          ImGui::TextDisabled("No history yet (time needs to advance)."
                              " Tip: use [ and ] to change time scale.");
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
      for (const auto& s : nearby) {
        const ImVec2 p = toPx(s.posLy);
        const bool isCurrent = (s.id == currentSystem->stub.id);
        const bool isSel = (s.id == galaxySelectedSystem);

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
          if (dx*dx + dy*dy <= rClick*rClick) {
            galaxySelectedSystem = s.id;
            galaxySelectedStationIndex = 0;
          }
        }
      }

      ImGui::EndChild();

      // Jump UI (requires docking)
      if (galaxySelectedSystem != currentSystem->stub.id) {
        auto it = std::find_if(nearby.begin(), nearby.end(), [&](const sim::SystemStub& s){ return s.id == galaxySelectedSystem; });
        if (it != nearby.end()) {
          const sim::SystemStub& dstStub = *it;
          const math::Vec3d d = dstStub.posLy - currentSystem->stub.posLy;
          const double distLy = std::sqrt(d.lengthSq());

          const int fuelNeeded = std::max(1, (int)std::ceil(distLy / 20.0));
          const double fuelHave = cargo[(std::size_t)econ::CommodityId::Fuel];

          const auto& dstSys = universe.getSystem(dstStub.id, &dstStub);
          if (!dstSys.stations.empty()) {
            galaxySelectedStationIndex = std::clamp(galaxySelectedStationIndex, 0, (int)dstSys.stations.size() - 1);
            std::vector<const char*> stNames;
            stNames.reserve(dstSys.stations.size());
            for (const auto& st : dstSys.stations) stNames.push_back(st.name.c_str());
            ImGui::Combo("Arrive at station", &galaxySelectedStationIndex, stNames.data(), (int)stNames.size());
          } else {
            galaxySelectedStationIndex = 0;
            ImGui::TextDisabled("Destination has no stations");
          }

          ImGui::Text("Selected: %s  (%.1f ly)", dstStub.name.c_str(), distLy);
          ImGui::Text("Jump cost: Fuel %d (you have %.0f)", fuelNeeded, fuelHave);

          const bool canJump = docked && (fuelHave >= (double)fuelNeeded);
          if (!docked) {
            ImGui::TextDisabled("To hyperjump: dock at a station");
          }

          ImGui::BeginDisabled(!canJump);
          if (ImGui::Button("Hyperjump")) {
            // Consume fuel
            cargo[(std::size_t)econ::CommodityId::Fuel] = std::max(0.0, fuelHave - (double)fuelNeeded);

            // Travel time (very rough)
            timeDays += distLy * 0.15;

            // Switch system
            currentStub = dstStub;
            const auto& sys = universe.getSystem(currentStub.id, &currentStub);
            currentSystem = &sys;

            // Reset ship & dock at selected destination station if possible
            ship.setVelocityKmS({0,0,0});
            ship.setAngularVelocityRadS({0,0,0});

            docked = false;
            dockedStationIndex = 0;
            viewStationIndex = 0;

            if (!sys.stations.empty()) {
              docked = true;
              dockedStationIndex = std::clamp(galaxySelectedStationIndex, 0, (int)sys.stations.size() - 1);
              viewStationIndex = dockedStationIndex;

              const auto& st = sys.stations[(std::size_t)dockedStationIndex];
              const math::Vec3d stPosKm = stationPositionKm(st, timeDays);
              ship.setPositionKm(stPosKm + math::Vec3d{0,0,-(st.radiusKm + 40.0)});
            } else {
              ship.setPositionKm({0,0,-8000.0});
            }

            // Reset nav
            autopilotEnabled = false;
            target = TargetRef{};
            galaxySelectedSystem = currentSystem->stub.id;
            galaxySelectedStationIndex = 0;
          }
          ImGui::EndDisabled();
        }
      }

      ImGui::TextDisabled("Tip: TAB toggles this window");

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
