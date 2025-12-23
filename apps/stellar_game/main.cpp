#include "stellar/core/Log.h"
#include "stellar/math/Math.h"
#include "stellar/render/Camera.h"
#include "stellar/render/Gl.h"
#include "stellar/render/PointRenderer.h"
#include "stellar/sim/Faction.h"
#include "stellar/sim/Market.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/Ship.h"
#include "stellar/sim/Universe.h"

#include <SDL.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace {

constexpr double kSecondsPerDay = 86400.0;

void printControls() {
  std::cout
    << "\nStellarForge (prototype) controls:\n"
    << "  Mouse              Look (yaw/pitch)\n"
    << "  W/S                Thrust forward/back\n"
    << "  A/D                Strafe left/right\n"
    << "  Q/E                Down/up\n"
    << "  Space              Brake\n"
    << "  Shift              Boost\n"
    << "  N / P              Next / Previous system index (streamed)\n"
    << "  J                  Random jump (new system index)\n"
    << "  M                  Print market for current system\n"
    << "  [ / ]              Select commodity\n"
    << "  , / .              Buy 1 / Sell 1 of selected commodity\n"
    << "  F5 / F9            Save / Load\n"
    << "  ESC                Toggle mouse capture\n"
    << "  Ctrl+C / Window X  Quit\n\n";
}

void printSystemBanner(
  std::size_t systemIndex,
  const stellar::sim::StarSystem& sys,
  const stellar::sim::Faction& fac)
{
  std::cout
    << "\n=== System #" << systemIndex << " ===\n"
    << "Star: " << sys.primary.name
    << " | Planets: " << sys.planets.size()
    << " | Faction: " << fac.name
    << "\n";
}

void printMarket(
  const stellar::sim::Market& market,
  const stellar::sim::CargoHold& cargo,
  double credits,
  std::size_t selected)
{
  std::cout << "\n--- Market ---\n";
  std::cout << "Credits: " << std::fixed << std::setprecision(2) << credits << "\n";
  for (std::size_t i = 0; i < stellar::sim::kCommodityCount; ++i) {
    const auto c = static_cast<stellar::sim::Commodity>(i);
    const auto* offer = market.find(c);
    if (!offer) continue;
    const bool isSel = (i == selected);
    std::cout
      << (isSel ? "> " : "  ")
      << std::setw(12) << stellar::sim::commodityName(c)
      << "  price=" << std::setw(8) << std::fixed << std::setprecision(2) << offer->price
      << "  supply=" << std::setw(6) << offer->supply
      << "  you=" << cargo[c]
      << "\n";
  }
  std::cout << "---------------\n";
}

stellar::render::PointVertex bodyVertex(const stellar::math::Vec3d& posAU, stellar::sim::PlanetType type) {
  stellar::render::PointVertex v;
  v.px = static_cast<float>(posAU.x);
  v.py = static_cast<float>(posAU.y);
  v.pz = static_cast<float>(posAU.z);

  switch (type) {
    case stellar::sim::PlanetType::Rocky:
      v.r = 0.75f; v.g = 0.75f; v.b = 0.78f; v.size = 6.0f;
      break;
    case stellar::sim::PlanetType::Ice:
      v.r = 0.55f; v.g = 0.75f; v.b = 1.0f; v.size = 6.0f;
      break;
    case stellar::sim::PlanetType::GasGiant:
      v.r = 1.0f; v.g = 0.72f; v.b = 0.35f; v.size = 10.0f;
      break;
    case stellar::sim::PlanetType::IceGiant:
      v.r = 0.55f; v.g = 0.85f; v.b = 0.95f; v.size = 8.0f;
      break;
    case stellar::sim::PlanetType::AsteroidBelt:
      v.r = 0.75f; v.g = 0.62f; v.b = 0.45f; v.size = 4.0f;
      break;
  }
  return v;
}

} // namespace

int main(int argc, char** argv) {
  stellar::core::setLogLevel(stellar::core::LogLevel::Info);

  std::filesystem::path savePath = "savegame.txt";
  bool tryLoad = true;
  stellar::core::u64 seed = 1;
  std::size_t factionCount = 32;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto requireValue = [&]() -> std::optional<std::string> {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << arg << "\n";
        return std::nullopt;
      }
      return std::string(argv[++i]);
    };

    if (arg == "--help" || arg == "-h") {
      std::cout
        << "StellarForge real-time prototype\n\n"
        << "Options:\n"
        << "  --seed <u64>        Universe seed (only used if no save is loaded)\n"
        << "  --save <path>       Save file path (default: savegame.txt)\n"
        << "  --no-load           Don't auto-load an existing save\n"
        << "  --factions <n>      Faction count (default: 32)\n";
      return 0;
    } else if (arg == "--seed") {
      const auto v = requireValue();
      if (!v) return 2;
      seed = static_cast<stellar::core::u64>(std::stoull(*v));
    } else if (arg == "--save") {
      const auto v = requireValue();
      if (!v) return 2;
      savePath = *v;
    } else if (arg == "--no-load") {
      tryLoad = false;
    } else if (arg == "--factions") {
      const auto v = requireValue();
      if (!v) return 2;
      factionCount = static_cast<std::size_t>(std::stoull(*v));
      if (factionCount == 0) factionCount = 1;
    } else {
      std::cerr << "Unknown option: " << arg << "\n";
      return 2;
    }
  }

  // --- Load or initialize game state ---
  stellar::sim::SaveGame save;
  if (tryLoad && std::filesystem::exists(savePath)) {
    std::string err;
    if (auto loaded = stellar::sim::readSave(savePath, &err)) {
      save = *loaded;
      seed = save.universeSeed;
      std::cout << "Loaded save from " << savePath << " (seed=" << seed << ")\n";
    } else {
      std::cerr << "Failed to load save: " << err << "\n";
    }
  }
  save.universeSeed = seed;
  if (save.selectedCommodity >= stellar::sim::kCommodityCount) {
    save.selectedCommodity = 0;
  }

  // --- SDL2 + OpenGL init ---
  SDL_SetMainReady();
  if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS) != 0) {
    std::cerr << "SDL_Init failed: " << SDL_GetError() << "\n";
    return 1;
  }

  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

#if defined(__APPLE__)
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG);
#endif

  SDL_Window* window = SDL_CreateWindow(
    "StellarForge",
    SDL_WINDOWPOS_CENTERED,
    SDL_WINDOWPOS_CENTERED,
    1280,
    720,
    SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);

  if (!window) {
    std::cerr << "SDL_CreateWindow failed: " << SDL_GetError() << "\n";
    SDL_Quit();
    return 1;
  }

  SDL_GLContext glctx = SDL_GL_CreateContext(window);
  if (!glctx) {
    std::cerr << "SDL_GL_CreateContext failed: " << SDL_GetError() << "\n";
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 1;
  }

  SDL_GL_MakeCurrent(window, glctx);
  SDL_GL_SetSwapInterval(1);

  std::string glErr;
  if (!stellar::render::gl::load(SDL_GL_GetProcAddress, &glErr)) {
    std::cerr << glErr << "\n";
    SDL_GL_DeleteContext(glctx);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 1;
  }

  stellar::render::PointRenderer renderer;
  std::string rErr;
  if (!renderer.init(&rErr)) {
    std::cerr << "Renderer init failed: " << rErr << "\n";
    SDL_GL_DeleteContext(glctx);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 1;
  }

  int winW = 1280;
  int winH = 720;
  SDL_GetWindowSize(window, &winW, &winH);
  renderer.resize(winW, winH);

  // --- Sim state ---
  stellar::sim::UniverseConfig ucfg;
  ucfg.seed = seed;
  ucfg.systemCacheSize = 512;
  ucfg.systemCountHint = 0;
  stellar::sim::Universe universe(ucfg);

  stellar::sim::FactionGenerator factions(seed, factionCount);
  stellar::sim::MarketGenerator markets(seed);

  stellar::sim::Ship ship;
  ship.positionAU = save.ship.positionAU;
  ship.velocityAUPerDay = save.ship.velocityAUPerDay;
  ship.yawDeg = save.ship.yawDeg;
  ship.pitchDeg = save.ship.pitchDeg;

  double simTimeDays = save.simTimeDays;
  std::size_t systemIndex = save.currentSystemIndex;

  auto currentSystem = universe.system(systemIndex);
  auto currentFaction = factions.controllingFaction(currentSystem->id);
  printSystemBanner(systemIndex, *currentSystem, currentFaction);

  printControls();

  bool running = true;
  bool mouseCaptured = true;
  SDL_SetRelativeMouseMode(SDL_TRUE);
  SDL_ShowCursor(SDL_DISABLE);

  stellar::core::SplitMix64 jumpRng(stellar::core::deriveSeed(seed, "jumps"));

  const std::uint64_t freq = SDL_GetPerformanceFrequency();
  std::uint64_t last = SDL_GetPerformanceCounter();

  while (running) {
    // dt
    const std::uint64_t now = SDL_GetPerformanceCounter();
    double dt = static_cast<double>(now - last) / static_cast<double>(freq);
    last = now;
    dt = std::clamp(dt, 0.0, 0.1);

    stellar::sim::ShipControls controls;

    SDL_Event e;
    while (SDL_PollEvent(&e)) {
      if (e.type == SDL_QUIT) {
        running = false;
      } else if (e.type == SDL_WINDOWEVENT) {
        if (e.window.event == SDL_WINDOWEVENT_SIZE_CHANGED) {
          renderer.resize(e.window.data1, e.window.data2);
          winW = e.window.data1;
          winH = e.window.data2;
        }
      } else if (e.type == SDL_KEYDOWN && !e.key.repeat) {
        const SDL_Keycode k = e.key.keysym.sym;
        if (k == SDLK_ESCAPE) {
          mouseCaptured = !mouseCaptured;
          SDL_SetRelativeMouseMode(mouseCaptured ? SDL_TRUE : SDL_FALSE);
          SDL_ShowCursor(mouseCaptured ? SDL_DISABLE : SDL_ENABLE);
        } else if (k == SDLK_F5) {
          save.version = 1;
          save.universeSeed = seed;
          save.currentSystemIndex = systemIndex;
          save.simTimeDays = simTimeDays;
          save.ship.positionAU = ship.positionAU;
          save.ship.velocityAUPerDay = ship.velocityAUPerDay;
          save.ship.yawDeg = ship.yawDeg;
          save.ship.pitchDeg = ship.pitchDeg;
          std::string err;
          if (stellar::sim::writeSave(save, savePath, &err)) {
            std::cout << "Saved to " << savePath << "\n";
          } else {
            std::cerr << "Save failed: " << err << "\n";
          }
        } else if (k == SDLK_F9) {
          std::string err;
          if (auto loaded = stellar::sim::readSave(savePath, &err)) {
            save = *loaded;
            seed = save.universeSeed;
            // Rebuild generators from seed.
            ucfg.seed = seed;
            universe = stellar::sim::Universe(ucfg);
            factions = stellar::sim::FactionGenerator(seed, factionCount);
            markets = stellar::sim::MarketGenerator(seed);

            simTimeDays = save.simTimeDays;
            systemIndex = save.currentSystemIndex;
            ship.positionAU = save.ship.positionAU;
            ship.velocityAUPerDay = save.ship.velocityAUPerDay;
            ship.yawDeg = save.ship.yawDeg;
            ship.pitchDeg = save.ship.pitchDeg;
            if (save.selectedCommodity >= stellar::sim::kCommodityCount) save.selectedCommodity = 0;

            currentSystem = universe.system(systemIndex);
            currentFaction = factions.controllingFaction(currentSystem->id);
            printSystemBanner(systemIndex, *currentSystem, currentFaction);
          } else {
            std::cerr << "Load failed: " << err << "\n";
          }
        } else if (k == SDLK_n) {
          ++systemIndex;
          currentSystem = universe.system(systemIndex);
          currentFaction = factions.controllingFaction(currentSystem->id);
          ship.positionAU = {0.0, 0.0, 5.0};
          ship.velocityAUPerDay = {0.0, 0.0, 0.0};
          printSystemBanner(systemIndex, *currentSystem, currentFaction);
        } else if (k == SDLK_p) {
          if (systemIndex > 0) {
            --systemIndex;
            currentSystem = universe.system(systemIndex);
            currentFaction = factions.controllingFaction(currentSystem->id);
            ship.positionAU = {0.0, 0.0, 5.0};
            ship.velocityAUPerDay = {0.0, 0.0, 0.0};
            printSystemBanner(systemIndex, *currentSystem, currentFaction);
          }
        } else if (k == SDLK_j) {
          systemIndex = static_cast<std::size_t>(jumpRng.nextU64());
          currentSystem = universe.system(systemIndex);
          currentFaction = factions.controllingFaction(currentSystem->id);
          ship.positionAU = {0.0, 0.0, 5.0};
          ship.velocityAUPerDay = {0.0, 0.0, 0.0};
          printSystemBanner(systemIndex, *currentSystem, currentFaction);
        } else if (k == SDLK_m) {
          const auto day = static_cast<std::int64_t>(std::floor(simTimeDays));
          const auto market = markets.generate(*currentSystem, currentFaction, day);
          printMarket(market, save.cargo, save.credits, save.selectedCommodity);
        } else if (k == SDLK_LEFTBRACKET) {
          if (save.selectedCommodity > 0) --save.selectedCommodity;
        } else if (k == SDLK_RIGHTBRACKET) {
          if (save.selectedCommodity + 1 < stellar::sim::kCommodityCount) ++save.selectedCommodity;
        } else if (k == SDLK_COMMA || k == SDLK_PERIOD) {
          const auto day = static_cast<std::int64_t>(std::floor(simTimeDays));
          const auto market = markets.generate(*currentSystem, currentFaction, day);
          const auto c = static_cast<stellar::sim::Commodity>(save.selectedCommodity);
          if (const auto* offer = market.find(c)) {
            if (k == SDLK_COMMA) {
              // buy 1
              if (save.credits >= offer->price) {
                save.credits -= offer->price;
                save.cargo[c] += 1;
              }
            } else {
              // sell 1
              if (save.cargo[c] > 0) {
                save.cargo[c] -= 1;
                save.credits += offer->price;
              }
            }
          }
        }
      } else if (e.type == SDL_MOUSEMOTION) {
        if (mouseCaptured) {
          const double sens = 0.08;
          controls.yawDeltaDeg += static_cast<double>(e.motion.xrel) * sens;
          controls.pitchDeltaDeg += static_cast<double>(-e.motion.yrel) * sens;
        }
      }
    }

    // Continuous controls
    const std::uint8_t* keys = SDL_GetKeyboardState(nullptr);
    controls.thrustForward = (keys[SDL_SCANCODE_W] ? 1.0 : 0.0) - (keys[SDL_SCANCODE_S] ? 1.0 : 0.0);
    controls.thrustRight = (keys[SDL_SCANCODE_D] ? 1.0 : 0.0) - (keys[SDL_SCANCODE_A] ? 1.0 : 0.0);
    controls.thrustUp = (keys[SDL_SCANCODE_E] ? 1.0 : 0.0) - (keys[SDL_SCANCODE_Q] ? 1.0 : 0.0);
    controls.brake = keys[SDL_SCANCODE_SPACE] != 0;
    controls.boost = (keys[SDL_SCANCODE_LSHIFT] || keys[SDL_SCANCODE_RSHIFT]);

    // Step sim
    ship.step(controls, dt);
    simTimeDays += dt / kSecondsPerDay;

    // Camera: third-person behind the ship.
    const auto shipFwd = ship.forward();
    const auto shipUp = ship.up();

    const stellar::math::Vec3d camPos = ship.positionAU - shipFwd * 15.0 + shipUp * 6.0;
    const stellar::math::Vec3d camTarget = ship.positionAU + shipFwd * 2.0;
    const auto dir = (camTarget - camPos).normalized();
    const double camYaw = stellar::math::radToDeg(std::atan2(dir.x, dir.z));
    const double camPitch = stellar::math::radToDeg(std::asin(std::clamp(dir.y, -1.0, 1.0)));

    stellar::render::Camera camera;
    camera.position = camPos;
    camera.yawDeg = camYaw;
    camera.pitchDeg = camPitch;
    camera.nearPlane = 0.01;
    camera.farPlane = 2000.0;

    // Build points: star + planets + ship.
    std::vector<stellar::render::PointVertex> pts;
    pts.reserve(2 + currentSystem->planets.size());

    // Star at origin
    pts.push_back(stellar::render::PointVertex{
      .px = 0.0f,
      .py = 0.0f,
      .pz = 0.0f,
      .r = 1.0f,
      .g = 0.95f,
      .b = 0.65f,
      .size = 16.0f,
    });

    for (const auto& p : currentSystem->planets) {
      const auto pos = p.orbit.positionAU(simTimeDays);
      pts.push_back(bodyVertex(pos, p.type));
    }

    // Ship
    pts.push_back(stellar::render::PointVertex{
      .px = static_cast<float>(ship.positionAU.x),
      .py = static_cast<float>(ship.positionAU.y),
      .pz = static_cast<float>(ship.positionAU.z),
      .r = 1.0f,
      .g = 1.0f,
      .b = 1.0f,
      .size = 7.0f,
    });

    const float aspect = (winH > 0) ? static_cast<float>(winW) / static_cast<float>(winH) : 1.0f;
    renderer.beginFrame(0.02f, 0.02f, 0.04f, 1.0f);
    renderer.drawPoints(pts, camera.viewProj(aspect));
    SDL_GL_SwapWindow(window);
  }

  renderer.shutdown();
  SDL_GL_DeleteContext(glctx);
  SDL_DestroyWindow(window);
  SDL_Quit();
  return 0;
}
