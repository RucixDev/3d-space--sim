#include "stellar/core/Hash.h"
#include "stellar/core/Log.h"
#include "stellar/core/Random.h"
#include "stellar/econ/Commodity.h"
#include "stellar/econ/Market.h"
#include "stellar/math/Math.h"
#include "stellar/math/Mat4.h"
#include "stellar/math/Vec3.h"
#include "stellar/proc/NameGenerator.h"
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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Dear ImGui
#include "imgui.h"
#include "imgui_impl_opengl3.h"
#include "imgui_impl_sdl2.h"

// ------------ Constants / helpers ------------

static constexpr double kAuKm = 149597870.7;
static constexpr double kLyKm = 9.4607304725808e12;
static constexpr double kSolarRadiusKm = 695700.0;

static double auToKm(double au) { return au * kAuKm; }

static stellar::math::Vec3d auToKm(const stellar::math::Vec3d& p) {
  return {p.x * kAuKm, p.y * kAuKm, p.z * kAuKm};
}

struct Vec4d {
  double x{0}, y{0}, z{0}, w{1};
};

static Vec4d mul(const stellar::math::Mat4d& m, const Vec4d& v) {
  // Column-major
  return {
    m.m[0] * v.x + m.m[4] * v.y + m.m[8] * v.z + m.m[12] * v.w,
    m.m[1] * v.x + m.m[5] * v.y + m.m[9] * v.z + m.m[13] * v.w,
    m.m[2] * v.x + m.m[6] * v.y + m.m[10] * v.z + m.m[14] * v.w,
    m.m[3] * v.x + m.m[7] * v.y + m.m[11] * v.z + m.m[15] * v.w,
  };
}

static stellar::math::Vec3d planetPosKm(const stellar::sim::Planet& p, double timeDays) {
  return auToKm(stellar::sim::orbitPosition3DAU(p.orbit, timeDays));
}

static stellar::math::Vec3d stationPosKm(const stellar::sim::Station& st, double timeDays) {
  return auToKm(stellar::sim::orbitPosition3DAU(st.orbit, timeDays));
}

static stellar::math::Vec3d orbitVelKmS(const stellar::sim::OrbitElements& orbit, double timeDays) {
  // Finite-difference velocity.
  const double dtDays = 10.0 / 86400.0; // 10 seconds
  const stellar::math::Vec3d p0 = auToKm(stellar::sim::orbitPosition3DAU(orbit, timeDays));
  const stellar::math::Vec3d p1 = auToKm(stellar::sim::orbitPosition3DAU(orbit, timeDays + dtDays));
  const stellar::math::Vec3d vKmPerDay = (p1 - p0) / dtDays;
  return vKmPerDay / 86400.0;
}

static ImVec2 worldToScreen(const stellar::math::Mat4d& view,
                            const stellar::math::Mat4d& proj,
                            const stellar::math::Vec3d& world,
                            const ImVec2& viewportSize,
                            bool* outBehind = nullptr) {
  const Vec4d clip = mul(proj, mul(view, Vec4d{world.x, world.y, world.z, 1.0}));
  if (outBehind) *outBehind = (clip.w <= 0.0);

  if (std::abs(clip.w) < 1e-9) return {-10000, -10000};

  const double ndcX = clip.x / clip.w;
  const double ndcY = clip.y / clip.w;

  // NDC (-1..1) to screen
  const float sx = static_cast<float>((ndcX * 0.5 + 0.5) * viewportSize.x);
  const float sy = static_cast<float>((-ndcY * 0.5 + 0.5) * viewportSize.y);
  return {sx, sy};
}

static stellar::render::Texture2D makeCheckerTexture(int w, int h) {
  std::vector<std::uint8_t> data(static_cast<std::size_t>(w * h * 4));
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      const bool c = ((x / 8) % 2) ^ ((y / 8) % 2);
      const std::uint8_t v = c ? 200 : 60;
      std::size_t idx = static_cast<std::size_t>((y * w + x) * 4);
      data[idx + 0] = v;
      data[idx + 1] = v;
      data[idx + 2] = v;
      data[idx + 3] = 255;
    }
  }
  stellar::render::Texture2D tex;
  tex.uploadRGBA(w, h, data.data());
  return tex;
}

static const char* stationTypeName(stellar::econ::StationType t) {
  using stellar::econ::StationType;
  switch (t) {
    case StationType::Outpost:      return "Outpost";
    case StationType::Agricultural: return "Agricultural";
    case StationType::Mining:       return "Mining";
    case StationType::Refinery:     return "Refinery";
    case StationType::Industrial:   return "Industrial";
    case StationType::Research:     return "Research";
    case StationType::TradeHub:     return "Trade Hub";
    case StationType::Shipyard:     return "Shipyard";
    default:                        return "Unknown";
  }
}

static ImVec4 stationTypeColor(stellar::econ::StationType t) {
  using stellar::econ::StationType;
  switch (t) {
    case StationType::Outpost:      return {0.7f, 0.7f, 0.7f, 1.0f};
    case StationType::Agricultural: return {0.4f, 0.85f, 0.4f, 1.0f};
    case StationType::Mining:       return {0.75f, 0.55f, 0.35f, 1.0f};
    case StationType::Refinery:     return {0.8f, 0.7f, 0.2f, 1.0f};
    case StationType::Industrial:   return {0.35f, 0.55f, 0.9f, 1.0f};
    case StationType::Research:     return {0.7f, 0.35f, 0.9f, 1.0f};
    case StationType::TradeHub:     return {0.9f, 0.65f, 0.2f, 1.0f};
    case StationType::Shipyard:     return {0.95f, 0.35f, 0.35f, 1.0f};
    default:                        return {1,1,1,1};
  }
}

// ------------ Gameplay types ------------

enum class NavTargetType : std::uint8_t {
  None = 0,
  Planet,
  Station,
  Npc,
};

struct NavTarget {
  NavTargetType type{NavTargetType::None};
  int planetIndex{-1};
  stellar::sim::StationId stationId{0};
  std::uint64_t npcId{0};
};

enum class NpcRole : std::uint8_t {
  Trader = 0,
  Pirate,
  Police,
};

enum class NpcState : std::uint8_t {
  Docked = 0,
  Travelling,
  Loiter,
};

struct NpcShip {
  std::uint64_t id{0};
  std::string name;
  NpcRole role{NpcRole::Trader};

  // Position is always in "normal space" km coordinates.
  stellar::math::Vec3d posKm{0,0,0};

  // Simple travel plan between stations (for traffic).
  stellar::sim::StationId origin{0};
  stellar::sim::StationId dest{0};
  double travelT{0.0};
  double travelDurationS{1.0};
  stellar::math::Vec3d travelStart{0,0,0};
  stellar::math::Vec3d travelEnd{0,0,0};

  NpcState state{NpcState::Docked};
  double stateTimerS{0.0};

  // Trading payload
  stellar::econ::CommodityId cargoCommodity{stellar::econ::CommodityId::Food};
  double cargoUnits{0.0};
};

struct DockingClearance {
  bool granted{false};
  double expiresDay{0.0};
};

struct StationTraffic {
  // A very small traffic gate: only one ship can be on final approach at a time.
  double occupiedUntilDay{0.0};
};

enum class FlightMode : std::uint8_t {
  Normal = 0,
  Supercruise,
  Hyperspace,
  Docked,
};

struct InterdictionEvent {
  bool active{false};
  double timeLeftS{0.0};
  double escapeMeter{0.0}; // 0..1
  stellar::math::Vec3d escapeDir{0,0,1};
  std::uint64_t pirateId{0};
};

struct Toast {
  std::string text;
  double timeLeftS{0.0};
};

static void pushToast(std::deque<Toast>& toasts, std::string text, double ttlS = 4.0) {
  toasts.push_front(Toast{std::move(text), ttlS});
  while (toasts.size() > 6) toasts.pop_back();
}

static double cargoMassKg(const std::array<double, stellar::econ::kCommodityCount>& cargo) {
  double m = 0.0;
  for (std::size_t i = 0; i < stellar::econ::kCommodityCount; ++i) {
    const auto id = static_cast<stellar::econ::CommodityId>(i);
    m += cargo[i] * stellar::econ::commodityDef(id).massKg;
  }
  return m;
}

static double cargoValueEstimate(const stellar::sim::StarSystem& sys,
                                 stellar::sim::Universe& uni,
                                 double timeDays,
                                 const std::array<double, stellar::econ::kCommodityCount>& cargo) {
  // Very rough estimate: use mid prices at the first station (if any) as a reference.
  if (sys.stations.empty()) return 0.0;
  const auto& st = sys.stations.front();
  auto& econ = uni.stationEconomy(st, timeDays);

  double v = 0.0;
  for (std::size_t i = 0; i < stellar::econ::kCommodityCount; ++i) {
    const auto id = static_cast<stellar::econ::CommodityId>(i);
    if (cargo[i] <= 0.0) continue;
    const double p = stellar::econ::midPrice(econ, st.economyModel, id);
    v += cargo[i] * p;
  }
  return v;
}

static std::optional<const stellar::sim::Station*> findStationById(const stellar::sim::StarSystem& sys,
                                                                   stellar::sim::StationId id) {
  for (const auto& st : sys.stations) {
    if (st.id == id) return &st;
  }
  return std::nullopt;
}

static std::optional<int> findPlanetIndexByName(const stellar::sim::StarSystem& sys, const std::string& name) {
  for (int i = 0; i < static_cast<int>(sys.planets.size()); ++i) {
    if (sys.planets[static_cast<std::size_t>(i)].name == name) return i;
  }
  return std::nullopt;
}

static const char* npcRoleName(NpcRole r) {
  switch (r) {
    case NpcRole::Trader: return "Trader";
    case NpcRole::Pirate: return "Pirate";
    case NpcRole::Police: return "Security";
    default: return "Contact";
  }
}

// Basic deterministic NPC name.
static std::string makeNpcName(std::uint64_t seed, NpcRole role) {
  stellar::proc::NameGenerator ng(seed);
  std::string base = ng.systemName();
  if (base.size() > 12) base.resize(12);
  switch (role) {
    case NpcRole::Trader: return "Hauler " + base;
    case NpcRole::Pirate: return "Pirate " + base;
    case NpcRole::Police: return "Sec " + base;
    default: return base;
  }
}

// ------------ Mission generation ------------

static const char* missionTypeName(stellar::sim::MissionType t) {
  using stellar::sim::MissionType;
  switch (t) {
    case MissionType::Courier: return "Courier";
    case MissionType::Delivery: return "Delivery";
    case MissionType::BountyScan: return "Bounty Scan";
    default: return "Mission";
  }
}

static std::vector<stellar::sim::Mission> generateMissionBoard(
    const stellar::sim::SystemStub& hereStub,
    const stellar::sim::StarSystem& hereSys,
    const stellar::sim::Station& fromStation,
    stellar::sim::Universe& uni,
    double timeDays,
    int count,
    double maxRadiusLy) {

  std::vector<stellar::sim::Mission> out;
  out.reserve(static_cast<std::size_t>(count));

  // Deterministic-ish seed: station id + day
  const std::uint64_t dayKey = static_cast<std::uint64_t>(std::floor(timeDays));
  stellar::core::SplitMix64 rng(stellar::core::hashCombine(static_cast<std::uint64_t>(fromStation.id), dayKey));

  auto nearby = uni.queryNearby(hereStub.posLy, maxRadiusLy, 96);
  if (nearby.empty()) return out;

  // Helper: pick a random destination station (possibly in another system).
  auto pickDestStation = [&]() -> std::optional<std::pair<stellar::sim::SystemStub, const stellar::sim::Station*>> {
    for (int tries = 0; tries < 12; ++tries) {
      const auto& dstStub = nearby[static_cast<std::size_t>(rng.rangeInt(0, static_cast<int>(nearby.size()) - 1))];
      const auto& dstSys = uni.getSystem(dstStub.id, &dstStub);
      if (dstSys.stations.empty()) continue;

      const auto& dstStation = dstSys.stations[static_cast<std::size_t>(rng.rangeInt(0, static_cast<int>(dstSys.stations.size()) - 1))];
      // Avoid "from -> same" loops.
      if (dstStub.id == hereStub.id && dstStation.id == fromStation.id) continue;
      return std::make_pair(dstStub, &dstStation);
    }
    return std::nullopt;
  };

  // Use current economy snapshot for availability weighting.
  auto& fromEcon = uni.stationEconomy(fromStation, timeDays);

  for (int i = 0; i < count; ++i) {
    auto dst = pickDestStation();
    if (!dst) break;

    const auto& dstStub = dst->first;
    const auto* dstStation = dst->second;
    const double distLy = (dstStub.posLy - hereStub.posLy).length();

    const double r = rng.nextDouble();
    stellar::sim::Mission m{};
    m.type = (r < 0.45) ? stellar::sim::MissionType::Courier
         : (r < 0.85) ? stellar::sim::MissionType::Delivery
                      : stellar::sim::MissionType::BountyScan;

    m.fromSystem = hereStub.id;
    m.fromStation = fromStation.id;
    m.toSystem = dstStub.id;
    m.toStation = dstStation->id;

    // Deadline scales with distance.
    m.deadlineDay = timeDays + std::max(1.0, distLy * 0.6) + rng.range(1.0, 6.0);

    // Base reward scales with distance.
    const double base = 120.0 + distLy * 22.0;

    if (m.type == stellar::sim::MissionType::Courier) {
      m.units = 0.0;
      m.reward = base * rng.range(0.9, 1.25);
      m.cargoProvided = false;
    } else if (m.type == stellar::sim::MissionType::Delivery) {
      // Choose a commodity the station actually has in stock.
      // Prefer "surplus" goods (inventory above desiredStock).
      int bestCid = -1;
      double bestScore = -1e9;
      for (std::size_t cid = 0; cid < stellar::econ::kCommodityCount; ++cid) {
        const double inv = fromEcon.inventory[cid];
        const double desired = fromStation.economyModel.desiredStock[cid];
        if (inv < 8.0) continue;
        const double surplus = inv - desired;
        const double score = surplus + rng.range(-5.0, 5.0);
        if (score > bestScore) {
          bestScore = score;
          bestCid = static_cast<int>(cid);
        }
      }

      if (bestCid < 0) {
        // Fallback: Food.
        bestCid = static_cast<int>(stellar::econ::CommodityId::Food);
      }
      m.commodity = static_cast<stellar::econ::CommodityId>(bestCid);

      // Units are limited by station stock, so scarcity matters.
      const double inv = fromEcon.inventory[static_cast<std::size_t>(m.commodity)];
      m.units = std::clamp(rng.range(6.0, 26.0), 4.0, std::max(4.0, inv * 0.25));
      m.cargoProvided = true;

      // Delivery is worth more than courier.
      m.reward = (base + m.units * 6.0) * rng.range(1.0, 1.35);
    } else {
      // Bounty scan: spawn a "wanted" ship in the target system.
      m.units = 0.0;
      m.cargoProvided = false;
      m.reward = (base + 180.0) * rng.range(1.0, 1.35);

      // Target id will be filled when accepted.
      m.targetNpcId = 0;
    }

    out.push_back(std::move(m));
  }

  return out;
}

// ------------ Trading helper ------------

static std::optional<std::pair<stellar::econ::CommodityId, double>> pickBestTrade(
    stellar::sim::Universe& uni,
    const stellar::sim::Station& origin,
    const stellar::sim::Station& dest,
    double timeDays,
    double maxUnits) {

  auto& e0 = uni.stationEconomy(origin, timeDays);
  auto& e1 = uni.stationEconomy(dest, timeDays);

  double bestProfitPerUnit = 0.0;
  stellar::econ::CommodityId bestId = stellar::econ::CommodityId::Food;

  for (std::size_t i = 0; i < stellar::econ::kCommodityCount; ++i) {
    auto id = static_cast<stellar::econ::CommodityId>(i);
    const auto q0 = stellar::econ::quote(e0, origin.economyModel, id);
    const auto q1 = stellar::econ::quote(e1, dest.economyModel, id);

    if (q0.ask <= 0.0 || q1.bid <= 0.0) continue;
    if (q0.inventory < 2.0) continue;

    const double profit = q1.bid - q0.ask;
    if (profit > bestProfitPerUnit) {
      bestProfitPerUnit = profit;
      bestId = id;
    }
  }

  if (bestProfitPerUnit <= 0.1) return std::nullopt;

  // Small-ish NPC loads to keep stations from being instantly drained.
  auto& e0b = uni.stationEconomy(origin, timeDays);
  const auto q0b = stellar::econ::quote(e0b, origin.economyModel, bestId);
  const double units = std::clamp(maxUnits, 2.0, std::max(2.0, q0b.inventory * 0.18));

  return std::make_pair(bestId, units);
}

// ------------ Main ------------

int main(int /*argc*/, char** /*argv*/) {
  // SDL
  if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_EVENTS) != 0) {
    std::fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
    return 1;
  }

  SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, 0);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
  SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

  SDL_Window* window = SDL_CreateWindow(
    "Stellar Forge (Early Playable)",
    SDL_WINDOWPOS_CENTERED,
    SDL_WINDOWPOS_CENTERED,
    1280,
    720,
    SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE
  );
  if (!window) {
    std::fprintf(stderr, "SDL_CreateWindow failed: %s\n", SDL_GetError());
    SDL_Quit();
    return 1;
  }

  SDL_GLContext gl_ctx = SDL_GL_CreateContext(window);
  SDL_GL_MakeCurrent(window, gl_ctx);
  SDL_GL_SetSwapInterval(1);

  if (!stellar::render::gl::init()) {
    std::fprintf(stderr, "Failed to load OpenGL functions.\n");
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 1;
  }

  // ImGui
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  (void)io;
  ImGui::StyleColorsDark();

  ImGui_ImplSDL2_InitForOpenGL(window, gl_ctx);
  ImGui_ImplOpenGL3_Init("#version 330");

  // Renderers
  std::string err;
  stellar::render::MeshRenderer meshRenderer;
  stellar::render::LineRenderer lineRenderer;
  stellar::render::PointRenderer pointRenderer;

  if (!meshRenderer.init(&err) || !lineRenderer.init(&err) || !pointRenderer.init(&err)) {
    std::fprintf(stderr, "Renderer init failed: %s\n", err.c_str());
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 1;
  }

  stellar::render::Mesh cube = stellar::render::Mesh::makeCube();
  stellar::render::Mesh sphere = stellar::render::Mesh::makeUvSphere(32, 16);
  stellar::render::Texture2D checker = makeCheckerTexture(64, 64);

  meshRenderer.setTexture(&checker);

  // Game state
  constexpr const char* kSavePath = "savegame.txt";
  stellar::sim::SaveGame save{};
  (void)stellar::sim::loadFromFile(kSavePath, save);

  // Universe
  if (save.seed == 0) save.seed = 1337;
  stellar::sim::Universe universe(save.seed);
  if (!save.stationOverrides.empty()) universe.importStationOverrides(save.stationOverrides);

  // Choose a system if missing.
  std::vector<stellar::sim::SystemStub> nearby;
  if (save.currentSystem == 0) {
    nearby = universe.queryNearby({0, 0, 0}, 45.0, 64);
    if (nearby.empty()) {
      std::fprintf(stderr, "No systems generated.\n");
      SDL_DestroyWindow(window);
      SDL_Quit();
      return 1;
    }
    save.currentSystem = nearby.front().id;
  }

  // Keep a cached stub for current system (for mission board + galaxy map distance).
  // If it isn't in the current nearby list, try to find it.
  if (nearby.empty()) {
    nearby = universe.queryNearby({0, 0, 0}, 80.0, 256);
  }
  auto findStub = [&](stellar::sim::SystemId id) -> std::optional<stellar::sim::SystemStub> {
    for (const auto& s : nearby) if (s.id == id) return s;
    // Fallback: search around the origin; this is cheap enough for now.
    auto more = universe.queryNearby({0, 0, 0}, 200.0, 512);
    for (const auto& s : more) if (s.id == id) return s;
    return std::nullopt;
  };

  std::optional<stellar::sim::SystemStub> currentStubOpt = findStub(save.currentSystem);
  stellar::sim::SystemStub currentStub = currentStubOpt.value_or(stellar::sim::SystemStub{});

  const stellar::sim::StarSystem* currentSys = &universe.getSystem(save.currentSystem, currentStubOpt ? &currentStub : nullptr);

  // Player ship
  stellar::sim::Ship ship;
  ship.setPositionKm(save.shipPosKm);
  ship.setVelocityKmS(save.shipVelKmS);
  ship.setOrientation(save.shipOrient);
  ship.setAngularVelocityRadS(save.shipAngVelRadS);

  // If no meaningful saved position, spawn near first station.
  if (ship.positionKm().length() < 1e-6 && !currentSys->stations.empty()) {
    const auto& st = currentSys->stations.front();
    const auto p = stationPosKm(st, save.timeDays);
    const auto axis = p.normalized();
    ship.setPositionKm(p + axis * (st.docking.corridorLengthKm * 0.85));
    ship.setVelocityKmS(orbitVelKmS(st.orbit, save.timeDays));
    ship.setOrientation(stellar::math::Quatd::fromAxisAngle({0, 1, 0}, 0));
  }

  // Core gameplay
  FlightMode mode = (save.dockedStation != 0) ? FlightMode::Docked : FlightMode::Normal;
  NavTarget navTarget{};
  bool autopilot = false;
  bool dampers = true;

  // Time acceleration: Pioneer-style levels.
  std::array<double, 5> timeWarpLevels{1.0, 10.0, 100.0, 1000.0, 10000.0};
  int timeWarpIndex = 0;

  // Supercruise
  double scThrottle = 0.0;         // 0..1
  bool scAssist = true;            // 7-second rule assist
  double scSpeedKmS = 0.0;
  double scMaxSpeedKmS = 220000.0; // ~0.73c (game scale)

  // FSD hyperspace
  double fuel = save.fuel;
  double fuelMax = save.fuelMax;
  double hull = save.hull;
  double cargoCapacityKg = save.cargoCapacityKg;
  double fsdReadyDay = save.fsdReadyDay;

  // Thermal
  bool scoopEnabled = true;
  double heat = 0.0; // 0..1+

  // Docking
  stellar::sim::StationId dockedStationId = save.dockedStation;
  std::unordered_map<stellar::sim::StationId, DockingClearance> clearance{};
  std::unordered_map<stellar::sim::StationId, StationTraffic> stationTraffic{};

  // Missions
  std::vector<stellar::sim::Mission> activeMissions = save.missions;
  std::uint64_t nextMissionId = save.nextMissionId;

  // Mission board cache
  stellar::sim::StationId boardStationId = 0;
  std::uint64_t boardDayKey = 0;
  std::vector<stellar::sim::Mission> boardMissions;

  // NPCs
  std::vector<NpcShip> npcs;

  auto regenNpcsForSystem = [&](const stellar::sim::StarSystem& sys, const stellar::sim::SystemStub& stub) {
    npcs.clear();
    stationTraffic.clear();
    clearance.clear();

    // Pre-create traffic state for stations.
    for (const auto& st : sys.stations) {
      stationTraffic[st.id] = StationTraffic{};
      clearance[st.id] = DockingClearance{};
    }

    stellar::core::SplitMix64 rng(stellar::core::hashCombine(static_cast<std::uint64_t>(stub.seed), 0xC0FFEEULL));

    const int trafficCount = std::clamp(static_cast<int>(sys.stations.size()) * 3, 5, 18);
    const int pirateCount = std::clamp(static_cast<int>(sys.stations.size()) / 2, 1, 4);
    const int policeCount = std::clamp(static_cast<int>(sys.stations.size()) / 3, 0, 2);

    auto pickStation = [&]() -> std::optional<const stellar::sim::Station*> {
      if (sys.stations.empty()) return std::nullopt;
      return &sys.stations[static_cast<std::size_t>(rng.rangeInt(0, static_cast<int>(sys.stations.size()) - 1))];
    };

    auto makeNpc = [&](NpcRole role) {
      NpcShip n{};
      n.id = stellar::core::hashCombine(static_cast<std::uint64_t>(stub.id), rng.nextU64());
      n.role = role;
      n.name = makeNpcName(stellar::core::hashCombine(static_cast<std::uint64_t>(n.id), 0xABCDEFULL), role);

      auto st0 = pickStation();
      auto st1 = pickStation();
      if (st0 && st1) {
        n.origin = (*st0)->id;
        n.dest = (*st1)->id;
        if (n.dest == n.origin && sys.stations.size() > 1) {
          n.dest = sys.stations[(static_cast<std::size_t>(rng.rangeInt(0, static_cast<int>(sys.stations.size()) - 1)) + 1) % sys.stations.size()].id;
        }
      }

      n.state = NpcState::Docked;
      n.stateTimerS = rng.range(2.0, 20.0);
      n.posKm = {0, 0, 0};

      npcs.push_back(std::move(n));
    };

    for (int i = 0; i < trafficCount; ++i) makeNpc(NpcRole::Trader);
    for (int i = 0; i < pirateCount; ++i) makeNpc(NpcRole::Pirate);
    for (int i = 0; i < policeCount; ++i) makeNpc(NpcRole::Police);

    // Place non-traders in loose orbits around the star.
    for (auto& n : npcs) {
      if (n.role == NpcRole::Trader) continue;
      const double rAU = rng.range(0.15, 1.8);
      const double ang = rng.range(0.0, 2.0 * stellar::math::kPi);
      n.posKm = {std::cos(ang) * auToKm(rAU), rng.range(-0.02, 0.02) * auToKm(rAU), std::sin(ang) * auToKm(rAU)};
      n.state = NpcState::Loiter;
      n.stateTimerS = rng.range(10.0, 45.0);
    }
  };

  regenNpcsForSystem(*currentSys, currentStub);

  // Interdiction
  InterdictionEvent interdiction{};

  // Toasts
  std::deque<Toast> toasts;

  // UI toggles
  bool showHelp = true;
  bool showNav = true;
  bool showMarket = true;
  bool showGalaxy = true;
  bool showMissions = true;

  // Remember last selected system in galaxy map.
  stellar::sim::SystemId galaxySelectedSystem = 0;

  // Random for moment-to-moment events
  stellar::core::SplitMix64 sessionRng(stellar::core::hashCombine(save.seed, 0x12345678ULL));

  // Render camera
  stellar::render::Camera cam;
  cam.setPerspective(stellar::math::degToRad(70.0), 1280.0 / 720.0, 0.01, 25000.0);

  // Main loop
  bool running = true;
  std::uint64_t lastTicks = SDL_GetPerformanceCounter();

  auto isDocked = [&]() { return mode == FlightMode::Docked && dockedStationId != 0; };

  auto setMode = [&](FlightMode m) {
    // Don't stack time-warp with other flight modes.
    if (m != FlightMode::Normal) {
      timeWarpIndex = 0;
    }

    mode = m;
    if (mode == FlightMode::Docked) {
      scSpeedKmS = 0.0;
      scThrottle = 0.0;
      autopilot = false;
      interdiction.active = false;
    }
  };

  auto currentSystemName = [&]() -> std::string {
    return currentSys ? currentSys->stub.name : std::string("?");
  };

  auto canTimeWarp = [&]() {
    if (mode != FlightMode::Normal) return false;

    if (autopilot) return false;

    // If player is thrusting hard, don't allow time accel.
    // (We apply time accel to physics; this avoids accidental slingshots.)
    const auto v = ship.velocityKmS();
    const double speed = v.length();
    (void)speed;

    // Must be reasonably far from stations & planets.
    for (const auto& st : currentSys->stations) {
      const auto p = stationPosKm(st, save.timeDays);
      if ((p - ship.positionKm()).length() < 2000.0) return false;
    }

    for (const auto& pl : currentSys->planets) {
      const auto p = planetPosKm(pl, save.timeDays);
      if ((p - ship.positionKm()).length() < 8000.0) return false;
    }
    return true;
  };

  auto clearTargetIfInvalid = [&]() {
    if (navTarget.type == NavTargetType::Station) {
      auto st = findStationById(*currentSys, navTarget.stationId);
      if (!st) navTarget = NavTarget{};
    } else if (navTarget.type == NavTargetType::Planet) {
      if (navTarget.planetIndex < 0 || navTarget.planetIndex >= static_cast<int>(currentSys->planets.size())) navTarget = NavTarget{};
    } else if (navTarget.type == NavTargetType::Npc) {
      const bool exists = std::any_of(npcs.begin(), npcs.end(), [&](const NpcShip& n) { return n.id == navTarget.npcId; });
      if (!exists) navTarget = NavTarget{};
    }
  };

  auto targetWorldPosKm = [&]() -> std::optional<stellar::math::Vec3d> {
    clearTargetIfInvalid();
    if (navTarget.type == NavTargetType::Station) {
      auto st = findStationById(*currentSys, navTarget.stationId);
      if (!st) return std::nullopt;
      return stationPosKm(**st, save.timeDays);
    }
    if (navTarget.type == NavTargetType::Planet) {
      if (navTarget.planetIndex < 0 || navTarget.planetIndex >= static_cast<int>(currentSys->planets.size())) return std::nullopt;
      return planetPosKm(currentSys->planets[static_cast<std::size_t>(navTarget.planetIndex)], save.timeDays);
    }
    if (navTarget.type == NavTargetType::Npc) {
      for (const auto& n : npcs) {
        if (n.id == navTarget.npcId) return n.posKm;
      }
      return std::nullopt;
    }
    return std::nullopt;
  };

  auto targetName = [&]() -> std::string {
    clearTargetIfInvalid();
    if (navTarget.type == NavTargetType::Station) {
      auto st = findStationById(*currentSys, navTarget.stationId);
      if (!st) return "(none)";
      return std::string("[Station] ") + (*st)->name;
    }
    if (navTarget.type == NavTargetType::Planet) {
      if (navTarget.planetIndex < 0 || navTarget.planetIndex >= static_cast<int>(currentSys->planets.size())) return "(none)";
      return std::string("[Planet] ") + currentSys->planets[static_cast<std::size_t>(navTarget.planetIndex)].name;
    }
    if (navTarget.type == NavTargetType::Npc) {
      for (const auto& n : npcs) if (n.id == navTarget.npcId) return std::string("[Contact] ") + n.name;
      return "(none)";
    }
    return "(none)";
  };

  auto tryDock = [&](const stellar::sim::Station& st) {
    const auto stPos = stationPosKm(st, save.timeDays);
    const auto stVel = orbitVelKmS(st.orbit, save.timeDays);

    const auto axis = stPos.normalized();
    const auto inbound = axis * -1.0;

    const auto r = ship.positionKm() - stPos;
    const double dist = r.length();

    // Within physical docking radius.
    if (dist > st.docking.radiusKm + 2.0) return false;

    // Must have clearance.
    auto& cl = clearance[st.id];
    if (!cl.granted || save.timeDays > cl.expiresDay) {
      return false;
    }

    // Must be slow.
    const double relSpeed = (ship.velocityKmS() - stVel).length();
    if (relSpeed > 0.03) return false;

    // Must be aligned.
    const double align = stellar::math::dot(ship.forward(), inbound);
    if (align < st.docking.alignCos) return false;

    dockedStationId = st.id;
    save.dockedStation = st.id;
    setMode(FlightMode::Docked);

    // Snap to a "docking port" point.
    ship.setPositionKm(stPos + axis * (st.docking.radiusKm + 1.0));
    ship.setVelocityKmS(stVel);

    pushToast(toasts, "Docked at " + st.name);
    return true;
  };

  auto requestDocking = [&](const stellar::sim::Station& st) {
    auto& traf = stationTraffic[st.id];
    auto& cl = clearance[st.id];

    const auto stPos = stationPosKm(st, save.timeDays);
    const double dist = (stPos - ship.positionKm()).length();
    if (dist > st.docking.commsRangeKm) {
      pushToast(toasts, "Out of comms range.");
      return;
    }

    // If currently occupied by another ship, deny.
    if (save.timeDays < traf.occupiedUntilDay) {
      pushToast(toasts, "Docking denied: traffic." );
      return;
    }

    traf.occupiedUntilDay = save.timeDays + (90.0 / 86400.0); // reserve ~90 seconds
    cl.granted = true;
    cl.expiresDay = save.timeDays + (6.0 / 1440.0); // 6 minutes

    pushToast(toasts, "Docking clearance granted.");
  };

  auto undock = [&]() {
    if (!isDocked()) return;

    auto stOpt = findStationById(*currentSys, dockedStationId);
    if (!stOpt) {
      dockedStationId = 0;
      save.dockedStation = 0;
      setMode(FlightMode::Normal);
      return;
    }

    const auto& st = **stOpt;
    const auto stPos = stationPosKm(st, save.timeDays);
    const auto stVel = orbitVelKmS(st.orbit, save.timeDays);
    const auto axis = stPos.normalized();

    ship.setPositionKm(stPos + axis * (st.docking.corridorLengthKm * 0.9));
    ship.setVelocityKmS(stVel);

    dockedStationId = 0;
    save.dockedStation = 0;
    setMode(FlightMode::Normal);

    pushToast(toasts, "Undocked.");
  };

  auto dropFromSupercruise = [&](bool hardDrop) {
    if (mode != FlightMode::Supercruise) return;

    // Convert supercruise speed to a smaller "normal" velocity.
    const double dropSpeed = std::min(6.0, std::max(0.0, scSpeedKmS * 0.00005));
    ship.setVelocityKmS(ship.forward() * dropSpeed);

    scThrottle = 0.0;
    scSpeedKmS = 0.0;
    setMode(FlightMode::Normal);

    if (hardDrop) {
      // Minor hull damage.
      hull = std::max(0.0, hull - 0.04);
      pushToast(toasts, "Emergency drop! Hull stressed.");
    } else {
      pushToast(toasts, "Dropped to normal space.");
    }
  };

  auto engageSupercruise = [&]() {
    if (mode != FlightMode::Normal) return;
    if (isDocked()) return;

    // Don't allow supercruise too close to a station.
    for (const auto& st : currentSys->stations) {
      const auto p = stationPosKm(st, save.timeDays);
      if ((p - ship.positionKm()).length() < 2000.0) {
        pushToast(toasts, "Too close to station to engage supercruise.");
        return;
      }
    }

    setMode(FlightMode::Supercruise);
    scThrottle = std::max(scThrottle, 0.1);
    pushToast(toasts, "Supercruise engaged.");
  };

  auto jumpFuelPerLy = [&]() {
    // Base cost per ly, scaled a bit by cargo mass.
    const double base = 1.0;
    const double massFactor = 1.0 + cargoMassKg(save.cargo) / std::max(1.0, cargoCapacityKg) * 0.35;
    return base * massFactor;
  };

  auto maxJumpRangeLy = [&]() {
    return (fuel <= 0.0) ? 0.0 : fuel / std::max(1e-6, jumpFuelPerLy());
  };

  auto arriveInSystemNearStation = [&](stellar::sim::SystemId sysId, const stellar::sim::SystemStub& stub) {
    currentStub = stub;
    currentSys = &universe.getSystem(sysId, &currentStub);
    save.currentSystem = sysId;

    regenNpcsForSystem(*currentSys, currentStub);

    // Put the player near the first station (or at origin if none).
    if (!currentSys->stations.empty()) {
      const auto& st = currentSys->stations.front();
      const auto p = stationPosKm(st, save.timeDays);
      const auto axis = p.normalized();
      ship.setPositionKm(p + axis * (st.docking.corridorLengthKm * 0.95));
      ship.setVelocityKmS(orbitVelKmS(st.orbit, save.timeDays));
      ship.setAngularVelocityRadS({0, 0, 0});
    } else {
      ship.setPositionKm({0, 0, auToKm(0.2)});
      ship.setVelocityKmS({0, 0, 0});
    }

    dockedStationId = 0;
    save.dockedStation = 0;
    setMode(FlightMode::Normal);

    pushToast(toasts, "Arrived in " + currentSystemName());
  };

  while (running) {
    // dt
    std::uint64_t nowTicks = SDL_GetPerformanceCounter();
    const double dt = static_cast<double>(nowTicks - lastTicks) / static_cast<double>(SDL_GetPerformanceFrequency());
    lastTicks = nowTicks;

    // Events
    SDL_Event e;
    while (SDL_PollEvent(&e)) {
      ImGui_ImplSDL2_ProcessEvent(&e);
      if (e.type == SDL_QUIT) running = false;
      if (e.type == SDL_WINDOWEVENT && e.window.event == SDL_WINDOWEVENT_CLOSE) running = false;

      if (e.type == SDL_KEYDOWN && !e.key.repeat) {
        const SDL_Keycode key = e.key.keysym.sym;
        if (key == SDLK_ESCAPE) running = false;

        // UI toggles
        if (key == SDLK_F1) showHelp = !showHelp;
        if (key == SDLK_F2) showNav = !showNav;
        if (key == SDLK_F3) showMarket = !showMarket;
        if (key == SDLK_F4) showGalaxy = !showGalaxy;
        if (key == SDLK_F5) showMissions = !showMissions;

        // Flight assist
        if (key == SDLK_f) {
          dampers = !dampers;
          pushToast(toasts, dampers ? "Flight assist: ON" : "Flight assist: OFF");
        }

        // Time warp levels (1..5)
        if (key >= SDLK_1 && key <= SDLK_5) {
          const int idx = static_cast<int>(key - SDLK_1);
          if (idx >= 0 && idx < static_cast<int>(timeWarpLevels.size())) {
            if (idx == 0 || canTimeWarp()) {
              timeWarpIndex = idx;
              pushToast(toasts, "Time accel: x" + std::to_string(static_cast<int>(timeWarpLevels[static_cast<std::size_t>(timeWarpIndex)])));
            } else {
              timeWarpIndex = 0;
              pushToast(toasts, "Time accel unavailable (not safe).");
            }
          }
        }

        // Docking
        if (key == SDLK_l) {
          // Request docking at target station, or nearest station.
          const stellar::sim::Station* st = nullptr;
          if (navTarget.type == NavTargetType::Station) {
            auto s = findStationById(*currentSys, navTarget.stationId);
            if (s) st = *s;
          }
          if (!st) {
            // nearest
            double best = std::numeric_limits<double>::max();
            for (const auto& cand : currentSys->stations) {
              const double d = (stationPosKm(cand, save.timeDays) - ship.positionKm()).length();
              if (d < best) {
                best = d;
                st = &cand;
              }
            }
          }
          if (st) requestDocking(*st);
        }

        if (key == SDLK_u) {
          undock();
        }

        // Autopilot
        if (key == SDLK_p) {
          autopilot = !autopilot;
          if (autopilot) {
            timeWarpIndex = 0;
          }
          pushToast(toasts, autopilot ? "Autopilot: ON" : "Autopilot: OFF");
        }

        // Supercruise
        if (key == SDLK_b) {
          if (mode == FlightMode::Supercruise) {
            dropFromSupercruise(false);
          } else if (mode == FlightMode::Normal) {
            engageSupercruise();
          }
        }

        if (key == SDLK_z) {
          scAssist = !scAssist;
          pushToast(toasts, scAssist ? "Supercruise assist: ON" : "Supercruise assist: OFF");
        }

        // Fuel scoop
        if (key == SDLK_o) {
          scoopEnabled = !scoopEnabled;
          pushToast(toasts, scoopEnabled ? "Fuel scoop: ON" : "Fuel scoop: OFF");
        }

        // Target cycling
        if (key == SDLK_t) {
          if (!currentSys->stations.empty()) {
            // cycle
            std::size_t idx = 0;
            if (navTarget.type == NavTargetType::Station) {
              for (std::size_t i = 0; i < currentSys->stations.size(); ++i) {
                if (currentSys->stations[i].id == navTarget.stationId) {
                  idx = (i + 1) % currentSys->stations.size();
                  break;
                }
              }
            }
            navTarget.type = NavTargetType::Station;
            navTarget.stationId = currentSys->stations[idx].id;
            pushToast(toasts, "Target: " + currentSys->stations[idx].name);
          }
        }

        if (key == SDLK_y) {
          if (!currentSys->planets.empty()) {
            int idx = 0;
            if (navTarget.type == NavTargetType::Planet) idx = (navTarget.planetIndex + 1) % static_cast<int>(currentSys->planets.size());
            navTarget.type = NavTargetType::Planet;
            navTarget.planetIndex = idx;
            pushToast(toasts, "Target: " + currentSys->planets[static_cast<std::size_t>(idx)].name);
          }
        }

        if (key == SDLK_g) {
          navTarget = NavTarget{};
          pushToast(toasts, "Target cleared.");
        }
      }
    }

    // Advance simulation time.
    const double simDt = dt * timeWarpLevels[static_cast<std::size_t>(timeWarpIndex)];
    save.timeDays += simDt / 86400.0;

    // Update toast timers.
    for (auto& t : toasts) t.timeLeftS -= dt;
    while (!toasts.empty() && toasts.back().timeLeftS <= 0.0) toasts.pop_back();

    // If in docked mode, keep ship glued.
    if (isDocked()) {
      auto stOpt = findStationById(*currentSys, dockedStationId);
      if (stOpt) {
        const auto& st = **stOpt;
        const auto stPos = stationPosKm(st, save.timeDays);
        const auto axis = stPos.normalized();
        const auto stVel = orbitVelKmS(st.orbit, save.timeDays);
        ship.setPositionKm(stPos + axis * (st.docking.radiusKm + 1.0));
        ship.setVelocityKmS(stVel);
      }
    }

    // Interdiction logic (supercruise only)
    if (mode == FlightMode::Supercruise && !interdiction.active) {
      const double cargoV = cargoValueEstimate(*currentSys, universe, save.timeDays, save.cargo);
      if (cargoV > 250.0 && sessionRng.nextDouble() < (dt * 0.015)) {
        // Chance per second.
        interdiction.active = true;
        interdiction.timeLeftS = 10.0;
        interdiction.escapeMeter = 0.0;
        // Random escape direction in world space.
        const double a = sessionRng.range(0.0, 2.0 * stellar::math::kPi);
        const double b = sessionRng.range(-0.35, 0.35);
        interdiction.escapeDir = stellar::math::Vec3d{std::cos(a) * std::cos(b), std::sin(b), std::sin(a) * std::cos(b)}.normalized();

        interdiction.pirateId = stellar::core::hashCombine(static_cast<std::uint64_t>(save.timeDays * 1000.0), sessionRng.nextU64());
        pushToast(toasts, "INTERDICTION WARNING!");

        // Spawn a fake pirate contact near the player.
        NpcShip pirate{};
        pirate.id = interdiction.pirateId;
        pirate.role = NpcRole::Pirate;
        pirate.name = makeNpcName(pirate.id, NpcRole::Pirate);
        pirate.state = NpcState::Loiter;
        pirate.posKm = ship.positionKm() - ship.forward() * 5000.0;
        npcs.push_back(std::move(pirate));
      }
    }

    if (interdiction.active) {
      interdiction.timeLeftS -= dt;

      const bool resist = (SDL_GetKeyboardState(nullptr)[SDL_SCANCODE_R] != 0);
      const double align = stellar::math::dot(ship.forward(), interdiction.escapeDir);

      if (resist && align > 0.92) {
        interdiction.escapeMeter += dt * (0.20 + 0.55 * (align - 0.92) / 0.08);
      } else {
        interdiction.escapeMeter -= dt * 0.08;
      }
      interdiction.escapeMeter = std::clamp(interdiction.escapeMeter, 0.0, 1.0);

      if (interdiction.escapeMeter >= 1.0) {
        interdiction.active = false;
        pushToast(toasts, "Interdiction evaded.");
      } else if (interdiction.timeLeftS <= 0.0) {
        interdiction.active = false;
        // Forced drop.
        dropFromSupercruise(true);
      }
    }

    // NPC update (simple)
    for (auto& n : npcs) {
      if (n.role != NpcRole::Trader) {
        // Simple loiter timer.
        if (n.state == NpcState::Loiter) {
          n.stateTimerS -= simDt;
          if (n.stateTimerS <= 0.0) {
            n.stateTimerS = sessionRng.range(8.0, 40.0);
            // Slight drift
            n.posKm += stellar::math::Vec3d{sessionRng.range(-1.0, 1.0), sessionRng.range(-0.3, 0.3), sessionRng.range(-1.0, 1.0)} * 40.0;
          }
        }
        continue;
      }

      if (n.state == NpcState::Docked) {
        n.stateTimerS -= simDt;
        if (n.stateTimerS <= 0.0 && n.origin != 0 && n.dest != 0 && n.origin != n.dest) {
          auto s0 = findStationById(*currentSys, n.origin);
          auto s1 = findStationById(*currentSys, n.dest);
          if (!s0 || !s1) {
            n.stateTimerS = sessionRng.range(5.0, 15.0);
            continue;
          }

          const auto p0 = stationPosKm(**s0, save.timeDays);
          const auto p1 = stationPosKm(**s1, save.timeDays);

          // Do a small economy trade to add real scarcity.
          if (auto pick = pickBestTrade(universe, **s0, **s1, save.timeDays, sessionRng.range(4.0, 20.0))) {
            n.cargoCommodity = pick->first;
            n.cargoUnits = pick->second;

            auto& e0 = universe.stationEconomy(**s0, save.timeDays);
            double dummyCredits = 1e12;
            (void)stellar::econ::buy(e0, (*s0)->economyModel, n.cargoCommodity, n.cargoUnits, dummyCredits, 0.0, (*s0)->feeRate);
          } else {
            n.cargoUnits = 0.0;
          }

          n.travelStart = p0;
          n.travelEnd = p1;
          const double dist = (p1 - p0).length();
          const double speed = sessionRng.range(12000.0, 38000.0); // km/s (supercruise-ish)
          n.travelDurationS = std::max(2.0, dist / speed);
          n.travelT = 0.0;
          n.state = NpcState::Travelling;

          // Spawn just outside origin station.
          n.posKm = p0 + (p0.normalized()) * 200.0;
        }
      } else if (n.state == NpcState::Travelling) {
        n.travelT += simDt;
        const double u = std::clamp(n.travelT / std::max(1e-6, n.travelDurationS), 0.0, 1.0);
        n.posKm = n.travelStart * (1.0 - u) + n.travelEnd * u;

        if (u >= 1.0 - 1e-6) {
          // Arrive; sell cargo.
          auto s1 = findStationById(*currentSys, n.dest);
          if (s1 && n.cargoUnits > 0.0) {
            auto& e1 = universe.stationEconomy(**s1, save.timeDays);
            double dummyCredits = 0.0;
            (void)stellar::econ::sell(e1, (*s1)->economyModel, n.cargoCommodity, n.cargoUnits, dummyCredits, 0.0, (*s1)->feeRate);
          }

          // Swap origin/dest and dock for a while.
          std::swap(n.origin, n.dest);
          n.state = NpcState::Docked;
          n.stateTimerS = sessionRng.range(3.0, 18.0);
        }
      }
    }

    // Player update
    const Uint8* keys = SDL_GetKeyboardState(nullptr);

    if (!isDocked()) {
      // Player input (manual)
      stellar::sim::ShipInput input{};
      input.dampers = dampers;

      // Translation
      //  W/S forward/back, A/D strafe left/right, Space/Ctrl up/down.
      const double fb = (keys[SDL_SCANCODE_W] ? 1.0 : 0.0) - (keys[SDL_SCANCODE_S] ? 1.0 : 0.0);
      const double lr = (keys[SDL_SCANCODE_D] ? 1.0 : 0.0) - (keys[SDL_SCANCODE_A] ? 1.0 : 0.0);
      const double ud = (keys[SDL_SCANCODE_SPACE] ? 1.0 : 0.0) - (keys[SDL_SCANCODE_LCTRL] ? 1.0 : 0.0);

      // Rotation
      const double yaw = (keys[SDL_SCANCODE_RIGHT] ? 1.0 : 0.0) - (keys[SDL_SCANCODE_LEFT] ? 1.0 : 0.0);
      const double pitch = (keys[SDL_SCANCODE_UP] ? 1.0 : 0.0) - (keys[SDL_SCANCODE_DOWN] ? 1.0 : 0.0);
      const double roll = (keys[SDL_SCANCODE_E] ? 1.0 : 0.0) - (keys[SDL_SCANCODE_Q] ? 1.0 : 0.0);

      // If the player starts maneuvering, drop time acceleration back to 1x.
      const bool manualInput =
        (std::abs(fb) > 0.01) || (std::abs(lr) > 0.01) || (std::abs(ud) > 0.01) ||
        (std::abs(yaw) > 0.01) || (std::abs(pitch) > 0.01) || (std::abs(roll) > 0.01) ||
        (keys[SDL_SCANCODE_LSHIFT] != 0) || (keys[SDL_SCANCODE_X] != 0);

      if (manualInput && timeWarpIndex > 0) {
        timeWarpIndex = 0;
        pushToast(toasts, "Time accel disengaged (manual input).");
      }

      input.thrustLocal = { lr, ud, fb };
      input.torqueLocal = { pitch, yaw, roll };

      input.boost = keys[SDL_SCANCODE_LSHIFT] != 0;
      input.brake = keys[SDL_SCANCODE_X] != 0;

      // Autopilot overrides input.
      if (autopilot && navTarget.type != NavTargetType::None) {
        if (auto tPos = targetWorldPosKm()) {
          // Simple station approach autopilot.
          if (navTarget.type == NavTargetType::Station) {
            auto stOpt = findStationById(*currentSys, navTarget.stationId);
            if (stOpt) {
              const auto& st = **stOpt;
              const auto stPos = stationPosKm(st, save.timeDays);
              const auto stVel = orbitVelKmS(st.orbit, save.timeDays);
              const auto axis = stPos.normalized();
              const auto inbound = axis * -1.0;

              const auto r = ship.positionKm() - stPos;
              const double t = stellar::math::dot(r, axis);
              const double lateralSq = std::max(0.0, r.lengthSq() - t * t);
              const double lateral = std::sqrt(lateralSq);

              const bool inCorridor = (t >= 0.0 && t <= st.docking.corridorLengthKm && lateral <= st.docking.corridorRadiusKm);

              stellar::math::Vec3d desiredPos = stPos + axis * (st.docking.corridorLengthKm * 0.85);
              if (inCorridor) {
                // Move down the corridor.
                desiredPos = stPos + axis * std::max(0.0, t - 180.0);
              }

              // Desired velocity: relative to station.
              const auto relPos = desiredPos - ship.positionKm();
              const auto relVel = (ship.velocityKmS() - stVel);

              // PD control
              const double kp = inCorridor ? 0.015 : 0.010;
              const double kd = 0.28;
              stellar::math::Vec3d desiredAcc = relPos * kp - relVel * kd;

              // Clamp
              const double maxA = ship.maxLinearAccelKmS2() * (input.boost ? 1.8 : 1.0);
              const double aLen = desiredAcc.length();
              if (aLen > maxA && aLen > 1e-9) desiredAcc *= (maxA / aLen);

              // Convert to local thrust
              const auto fwd = ship.forward();
              const auto right = ship.right();
              const auto up = ship.up();

              const double ax = stellar::math::dot(desiredAcc, right) / std::max(1e-6, ship.maxLinearAccelKmS2());
              const double ay = stellar::math::dot(desiredAcc, up) / std::max(1e-6, ship.maxLinearAccelKmS2());
              const double az = stellar::math::dot(desiredAcc, fwd) / std::max(1e-6, ship.maxLinearAccelKmS2());

              input.thrustLocal = { std::clamp(ax, -1.0, 1.0), std::clamp(ay, -1.0, 1.0), std::clamp(az, -1.0, 1.0) };

              // Orientation: align to inbound axis once close.
              const double wantAlign = inCorridor ? 1.0 : 0.0;
              const double curAlign = stellar::math::dot(ship.forward(), inbound);
              if (wantAlign > 0.5 && curAlign < 0.999) {
                const auto axisRot = stellar::math::cross(ship.forward(), inbound);
                const double axisLen = axisRot.length();
                if (axisLen > 1e-6) {
                  const auto nAxis = axisRot / axisLen;
                  // Use axis components as torque hints (very simple).
                  const auto local = stellar::math::Vec3d{
                    stellar::math::dot(nAxis, right),
                    stellar::math::dot(nAxis, up),
                    stellar::math::dot(nAxis, fwd)
                  };
                  input.torqueLocal = { std::clamp(local.x * 2.0, -1.0, 1.0), std::clamp(local.y * 2.0, -1.0, 1.0), std::clamp(local.z * 1.0, -1.0, 1.0) };
                }
              }

              // If we have clearance and we are basically at the station, dock.
              if (tryDock(st)) {
                input.thrustLocal = {0, 0, 0};
                input.torqueLocal = {0, 0, 0};
              }

              // Respect speed limit in corridor.
              if (inCorridor) {
                const double relSp = relVel.length();
                if (relSp > st.docking.speedLimitKmS) {
                  input.brake = true;
                }
              }
            }
          }
        }
      }

      // Supercruise uses simplified physics for translation.
      if (mode == FlightMode::Supercruise) {
        // Throttle control.
        const double dThrottle = ((keys[SDL_SCANCODE_PAGEUP] ? 1.0 : 0.0) - (keys[SDL_SCANCODE_PAGEDOWN] ? 1.0 : 0.0)) * dt * 0.35;
        scThrottle = std::clamp(scThrottle + dThrottle, 0.0, 1.0);

        // 7-second rule style assist
        if (scAssist) {
          if (auto tPos = targetWorldPosKm()) {
            const double dist = (*tPos - ship.positionKm()).length();
            const double desiredTimeS = 7.0;
            double desiredSpeed = dist / std::max(0.5, desiredTimeS);
            desiredSpeed = std::clamp(desiredSpeed, 2000.0, scMaxSpeedKmS);

            // Convert desired speed to a throttle fraction.
            const double targetThrottle = std::clamp(desiredSpeed / scMaxSpeedKmS, 0.0, 1.0);
            const double rate = 0.35;
            scThrottle += (targetThrottle - scThrottle) * std::clamp(dt * rate, 0.0, 1.0);
          }
        }

        // Speed response
        const double desiredSpeed = scThrottle * scMaxSpeedKmS;
        scSpeedKmS += (desiredSpeed - scSpeedKmS) * std::clamp(dt * 0.5, 0.0, 1.0);

        // Drop automatically when close to target.
        if (auto tPos = targetWorldPosKm()) {
          const double dist = (*tPos - ship.positionKm()).length();
          if (dist < 1200.0 && scSpeedKmS < 9000.0) {
            dropFromSupercruise(false);
          }
        }

        // Very basic "mass lock": if near a station, clamp speed.
        for (const auto& st : currentSys->stations) {
          const double d = (stationPosKm(st, save.timeDays) - ship.positionKm()).length();
          if (d < 5000.0) {
            scSpeedKmS = std::min(scSpeedKmS, 12000.0);
          }
        }

        // Use ship.step only for rotation.
        input.thrustLocal = {0, 0, 0};
        ship.step(simDt, input);

        // Translate along forward.
        ship.setPositionKm(ship.positionKm() + ship.forward() * (scSpeedKmS * simDt));

      } else {
        // Normal flight.
        ship.step(simDt, input);
      }
    }

    // Fuel scoop + heat model
    {
      const double starRadiusKm = currentSys->star.radiusSol * kSolarRadiusKm;
      const double distToStar = ship.positionKm().length();

      const double scoopMin = starRadiusKm * 2.0;
      const double scoopMax = starRadiusKm * 55.0;

      // Heat rises when close.
      const double closeness = std::clamp(1.0 - (distToStar - scoopMin) / std::max(1.0, scoopMax - scoopMin), 0.0, 1.0);
      const double heatGain = closeness * closeness * (scoopEnabled ? 0.35 : 0.22);

      // Heat also rises a bit with very high supercruise speeds.
      const double scHeat = (mode == FlightMode::Supercruise) ? std::clamp(scSpeedKmS / scMaxSpeedKmS, 0.0, 1.0) * 0.08 : 0.0;

      heat += (heatGain + scHeat) * dt;
      heat -= 0.18 * dt;
      heat = std::max(0.0, heat);

      // Scoop fuel if enabled and in range.
      if (scoopEnabled && distToStar > scoopMin && distToStar < scoopMax) {
        const double rate = 1.2 + 1.8 * closeness; // units / second
        fuel = std::min(fuelMax, fuel + rate * dt);
      }

      if (heat > 1.0) {
        hull = std::max(0.0, hull - (heat - 1.0) * 0.035 * dt);
        if (mode == FlightMode::Supercruise) {
          // Overheat forces drop.
          dropFromSupercruise(true);
        }
      }
    }

    // Auto-dock if physically right on top of a station.
    if (mode == FlightMode::Normal && !isDocked()) {
      for (const auto& st : currentSys->stations) {
        if (tryDock(st)) break;
      }
    }

    // If hull is dead, reset.
    if (hull <= 0.0) {
      hull = 1.0;
      fuel = std::max(fuelMax * 0.35, fuel);
      save.credits = std::max(0.0, save.credits - 200.0);
      pushToast(toasts, "Ship destroyed! Emergency recovery fee paid.");

      // Respawn at first station.
      if (!currentSys->stations.empty()) {
        const auto& st = currentSys->stations.front();
        const auto p = stationPosKm(st, save.timeDays);
        ship.setPositionKm(p + p.normalized() * (st.docking.corridorLengthKm * 0.9));
        ship.setVelocityKmS(orbitVelKmS(st.orbit, save.timeDays));
      }
      setMode(FlightMode::Normal);
    }

    // Mission completion / failure checks.
    for (auto& m : activeMissions) {
      if (m.completed || m.failed) continue;
      if (save.timeDays > m.deadlineDay) {
        m.failed = true;
        pushToast(toasts, "Mission failed (deadline).", 6.0);
        continue;
      }

      if (isDocked() && dockedStationId == m.toStation && save.currentSystem == m.toSystem) {
        if (m.type == stellar::sim::MissionType::Courier) {
          m.completed = true;
          save.credits += m.reward;
          pushToast(toasts, "Courier completed. +" + std::to_string(static_cast<int>(m.reward)) + " cr");
        } else if (m.type == stellar::sim::MissionType::Delivery) {
          const std::size_t cid = static_cast<std::size_t>(m.commodity);
          if (save.cargo[cid] >= m.units) {
            save.cargo[cid] -= m.units;
            m.completed = true;
            save.credits += m.reward;
            pushToast(toasts, "Delivery completed. +" + std::to_string(static_cast<int>(m.reward)) + " cr");
          }
        }
      }

      if (m.type == stellar::sim::MissionType::BountyScan && !m.completed && !m.failed) {
        // Scan completes in space: be near target NPC and hold V.
        if (save.currentSystem == m.toSystem && m.targetNpcId != 0) {
          for (const auto& n : npcs) {
            if (n.id != m.targetNpcId) continue;
            const double d = (n.posKm - ship.positionKm()).length();
            if (d < 12.0) {
              if (keys[SDL_SCANCODE_V]) {
                // immediate for now
                m.completed = true;
                save.credits += m.reward;
                pushToast(toasts, "Bounty scan complete. +" + std::to_string(static_cast<int>(m.reward)) + " cr");
              }
            }
          }
        }
      }
    }

    // Save persistent values back into save struct each frame.
    save.shipPosKm = ship.positionKm();
    save.shipVelKmS = ship.velocityKmS();
    save.shipOrient = ship.orientation();
    save.shipAngVelRadS = ship.angularVelocityRadS();

    save.fuel = fuel;
    save.fuelMax = fuelMax;
    save.hull = hull;
    save.cargoCapacityKg = cargoCapacityKg;
    save.fsdReadyDay = fsdReadyDay;
    save.missions = activeMissions;
    save.nextMissionId = nextMissionId;

    // Update station overrides snapshot occasionally.
    save.stationOverrides = universe.exportStationOverrides();

    // ---------------- Rendering ----------------

    int w = 1280, h = 720;
    SDL_GetWindowSize(window, &w, &h);
    cam.setPerspective(stellar::math::degToRad(70.0), static_cast<double>(w) / static_cast<double>(h), 0.01, 25000.0);

    // Camera: third-person behind ship
    const double worldScale = 1.0 / 1e6; // 1 unit = 1e6 km
    const auto shipPosU = ship.positionKm() * worldScale;

    const auto forward = ship.forward();
    const auto camPosKm = ship.positionKm() - forward * 220.0 + stellar::math::Vec3d{0, 60.0, 0};
    cam.setPosition(camPosKm * worldScale);
    cam.setOrientation(ship.orientation());

    const auto view = cam.viewMatrix();
    const auto proj = cam.projectionMatrix();

    float viewF[16];
    float projF[16];
    for (int i = 0; i < 16; ++i) {
      viewF[i] = static_cast<float>(view.m[i]);
      projF[i] = static_cast<float>(proj.m[i]);
    }

    meshRenderer.setViewProj(viewF, projF);
    lineRenderer.setViewProj(viewF, projF);
    pointRenderer.setViewProj(viewF, projF);

    // Instances
    std::vector<stellar::render::InstanceData> cubes;
    std::vector<stellar::render::InstanceData> spheres;
    cubes.reserve(256);
    spheres.reserve(512);

    // Star
    {
      stellar::render::InstanceData s{};
      s.px = 0.0f;
      s.py = 0.0f;
      s.pz = 0.0f;
      const double r = currentSys->star.radiusSol * kSolarRadiusKm * worldScale;
      s.scale = static_cast<float>(std::clamp(r, 0.002, 1.5));
      s.cr = 1.0f;
      s.cg = 0.95f;
      s.cb = 0.75f;
      spheres.push_back(s);
    }

    // Planets
    for (const auto& p : currentSys->planets) {
      const auto pk = planetPosKm(p, save.timeDays) * worldScale;

      stellar::render::InstanceData inst{};
      inst.px = static_cast<float>(pk.x);
      inst.py = static_cast<float>(pk.y);
      inst.pz = static_cast<float>(pk.z);
      inst.scale = static_cast<float>(0.0025 + 0.0012 * p.radiusEarth);

      // color by type
      switch (p.type) {
        case stellar::sim::PlanetType::Rocky:   inst.cr = 0.7f; inst.cg = 0.5f; inst.cb = 0.4f; break;
        case stellar::sim::PlanetType::Desert:  inst.cr = 0.9f; inst.cg = 0.7f; inst.cb = 0.3f; break;
        case stellar::sim::PlanetType::Ocean:   inst.cr = 0.2f; inst.cg = 0.5f; inst.cb = 0.9f; break;
        case stellar::sim::PlanetType::Ice:     inst.cr = 0.8f; inst.cg = 0.9f; inst.cb = 1.0f; break;
        case stellar::sim::PlanetType::GasGiant:inst.cr = 0.6f; inst.cg = 0.8f; inst.cb = 0.5f; break;
        default: inst.cr = inst.cg = inst.cb = 1.0f; break;
      }

      spheres.push_back(inst);
    }

    // Stations
    for (const auto& st : currentSys->stations) {
      const auto sk = stationPosKm(st, save.timeDays) * worldScale;
      const auto c = stationTypeColor(st.type);

      stellar::render::InstanceData inst{};
      inst.px = static_cast<float>(sk.x);
      inst.py = static_cast<float>(sk.y);
      inst.pz = static_cast<float>(sk.z);
      inst.scale = static_cast<float>(std::clamp(st.docking.radiusKm * worldScale * 0.7, 0.0012, 0.03));
      inst.cr = c.x;
      inst.cg = c.y;
      inst.cb = c.z;
      cubes.push_back(inst);
    }

    // NPC contacts
    for (const auto& n : npcs) {
      const auto nk = n.posKm * worldScale;
      stellar::render::InstanceData inst{};
      inst.px = static_cast<float>(nk.x);
      inst.py = static_cast<float>(nk.y);
      inst.pz = static_cast<float>(nk.z);
      inst.scale = 0.0012f;

      if (n.role == NpcRole::Trader) { inst.cr = 0.2f; inst.cg = 0.9f; inst.cb = 0.8f; }
      else if (n.role == NpcRole::Pirate) { inst.cr = 0.95f; inst.cg = 0.2f; inst.cb = 0.2f; }
      else { inst.cr = 0.25f; inst.cg = 0.55f; inst.cb = 0.95f; }

      cubes.push_back(inst);
    }

    // Player ship
    {
      stellar::render::InstanceData inst{};
      inst.px = static_cast<float>(shipPosU.x);
      inst.py = static_cast<float>(shipPosU.y);
      inst.pz = static_cast<float>(shipPosU.z);
      inst.scale = 0.0016f;
      inst.cr = 1.0f;
      inst.cg = 1.0f;
      inst.cb = 1.0f;
      cubes.push_back(inst);
    }

    // Docking corridor lines (target station)
    std::vector<stellar::render::LineVertex> lines;
    if (navTarget.type == NavTargetType::Station) {
      auto stOpt = findStationById(*currentSys, navTarget.stationId);
      if (stOpt) {
        const auto& st = **stOpt;
        const auto p = stationPosKm(st, save.timeDays);
        const auto axis = p.normalized();

        // Draw a simple center line.
        const auto p0 = p;
        const auto p1 = p + axis * st.docking.corridorLengthKm;

        stellar::render::LineVertex a{};
        a.px = static_cast<float>((p0.x * worldScale));
        a.py = static_cast<float>((p0.y * worldScale));
        a.pz = static_cast<float>((p0.z * worldScale));
        a.cr = 1.0f; a.cg = 0.9f; a.cb = 0.2f;

        stellar::render::LineVertex b{};
        b.px = static_cast<float>((p1.x * worldScale));
        b.py = static_cast<float>((p1.y * worldScale));
        b.pz = static_cast<float>((p1.z * worldScale));
        b.cr = 1.0f; b.cg = 0.9f; b.cb = 0.2f;

        lines.push_back(a);
        lines.push_back(b);
      }
    }

    // Starfield points (cheap, deterministic)
    std::vector<stellar::render::PointVertex> stars;
    {
      constexpr int kStars = 600;
      stars.reserve(kStars);
      stellar::core::SplitMix64 rng(stellar::core::hashCombine(save.seed, 0xDEADBEEF));
      for (int i = 0; i < kStars; ++i) {
        const double r = rng.range(4000.0, 20000.0);
        const double a = rng.range(0.0, 2.0 * stellar::math::kPi);
        const double b = rng.range(-0.85, 0.85);

        const double x = std::cos(a) * std::cos(b) * r;
        const double y = std::sin(b) * r;
        const double z = std::sin(a) * std::cos(b) * r;

        stellar::render::PointVertex pv{};
        pv.px = static_cast<float>(x);
        pv.py = static_cast<float>(y);
        pv.pz = static_cast<float>(z);
        const float c = static_cast<float>(rng.range(0.7, 1.0));
        pv.cr = c;
        pv.cg = c;
        pv.cb = c;
        pv.size = static_cast<float>(rng.range(1.0, 2.2));
        stars.push_back(pv);
      }
    }

    // Clear
    stellar::render::gl::Viewport(0, 0, w, h);
    stellar::render::gl::Enable(GL_DEPTH_TEST);
    stellar::render::gl::ClearColor(0.02f, 0.03f, 0.05f, 1.0f);
    stellar::render::gl::Clear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    pointRenderer.drawPoints(stars);

    // Draw planets + star
    meshRenderer.setMesh(&sphere);
    meshRenderer.drawInstances(spheres);

    // Draw cubes
    meshRenderer.setMesh(&cube);
    meshRenderer.drawInstances(cubes);

    // Draw corridor lines
    if (!lines.empty()) lineRenderer.drawLines(lines);

    // ---------------- ImGui ----------------
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplSDL2_NewFrame(window);
    ImGui::NewFrame();

    const ImVec2 vpSize = ImGui::GetIO().DisplaySize;

    // HUD overlay
    {
      ImGui::SetNextWindowPos({10, 10});
      ImGui::SetNextWindowBgAlpha(0.35f);
      ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_AlwaysAutoResize;
      ImGui::Begin("HUD", nullptr, flags);

      const double speed = ship.velocityKmS().length();
      const char* modeName = (mode == FlightMode::Docked) ? "Docked" : (mode == FlightMode::Supercruise) ? "Supercruise" : "Normal";

      ImGui::Text("System: %s", currentSystemName().c_str());
      ImGui::Text("Mode: %s", modeName);
      ImGui::Text("Time: day %.2f (x%.0f)", save.timeDays, timeWarpLevels[static_cast<std::size_t>(timeWarpIndex)]);
      ImGui::Separator();

      ImGui::Text("Speed: %.2f km/s", speed);
      if (mode == FlightMode::Supercruise) {
        ImGui::Text("SC Speed: %.0f km/s  Throttle: %.0f%%", scSpeedKmS, scThrottle * 100.0);
      }

      ImGui::Text("Fuel: %.1f / %.1f", fuel, fuelMax);
      ImGui::Text("Heat: %.0f%%", heat * 100.0);
      ImGui::Text("Hull: %.0f%%", hull * 100.0);

      const double cm = cargoMassKg(save.cargo);
      ImGui::Text("Cargo: %.1f / %.1f kg", cm, cargoCapacityKg);
      ImGui::Text("Credits: %.0f", save.credits);

      ImGui::Separator();
      ImGui::Text("Target: %s", targetName().c_str());
      if (auto tPos = targetWorldPosKm()) {
        const double d = (*tPos - ship.positionKm()).length();
        ImGui::Text("Dist: %.0f km", d);
        if (mode == FlightMode::Supercruise && scSpeedKmS > 1.0) {
          ImGui::Text("TtT: %.1f s", d / scSpeedKmS);
        }
      }

      if (interdiction.active) {
        ImGui::Separator();
        ImGui::TextColored({1.0f, 0.25f, 0.25f, 1.0f}, "INTERDICTION!");
        ImGui::Text("Hold R and align with escape dir.");
        ImGui::ProgressBar(static_cast<float>(interdiction.escapeMeter), ImVec2(180, 0));
        ImGui::Text("%.1f s", interdiction.timeLeftS);
      }

      ImGui::End();
    }

    // Toasts (bottom-left)
    {
      ImGui::SetNextWindowPos({10, vpSize.y - 10}, 0, {0, 1});
      ImGui::SetNextWindowBgAlpha(0.0f);
      ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_AlwaysAutoResize;
      ImGui::Begin("Toasts", nullptr, flags);
      for (const auto& t : toasts) {
        ImGui::TextUnformatted(t.text.c_str());
      }
      ImGui::End();
    }

    if (showHelp) {
      ImGui::Begin("Help / Controls", &showHelp);
      ImGui::Text("Flight:");
      ImGui::BulletText("W/S = forward/back thrusters");
      ImGui::BulletText("A/D = strafe left/right");
      ImGui::BulletText("Space/Ctrl = up/down");
      ImGui::BulletText("Arrow keys = pitch/yaw");
      ImGui::BulletText("Q/E = roll");
      ImGui::BulletText("Shift = boost, X = brake, F = flight assist");
      ImGui::Separator();
      ImGui::Text("Nav / progression:");
      ImGui::BulletText("T = cycle station targets, Y = cycle planets, G = clear target");
      ImGui::BulletText("B = toggle Supercruise, Z = toggle SC Assist (7s rule)");
      ImGui::BulletText("O = toggle Fuel Scoop");
      ImGui::BulletText("L = request docking clearance, P = toggle autopilot, U = undock");
      ImGui::BulletText("1..5 = time accel levels (restricted)");
      ImGui::BulletText("Hold V near bounty target to scan");
      ImGui::Separator();
      ImGui::Text("UI: F1..F5 toggles panels.");
      ImGui::End();
    }

    // Nav panel
    if (showNav) {
      ImGui::Begin("Navigation", &showNav);
      ImGui::Text("Stations:");
      for (const auto& st : currentSys->stations) {
        const bool selected = (navTarget.type == NavTargetType::Station && navTarget.stationId == st.id);
        if (ImGui::Selectable((st.name + "##st").c_str(), selected)) {
          navTarget.type = NavTargetType::Station;
          navTarget.stationId = st.id;
        }
      }
      ImGui::Separator();
      ImGui::Text("Planets:");
      for (int i = 0; i < static_cast<int>(currentSys->planets.size()); ++i) {
        const auto& p = currentSys->planets[static_cast<std::size_t>(i)];
        const bool selected = (navTarget.type == NavTargetType::Planet && navTarget.planetIndex == i);
        if (ImGui::Selectable((p.name + "##pl").c_str(), selected)) {
          navTarget.type = NavTargetType::Planet;
          navTarget.planetIndex = i;
        }
      }
      ImGui::Separator();
      ImGui::Text("Contacts:");
      for (const auto& n : npcs) {
        const bool selected = (navTarget.type == NavTargetType::Npc && navTarget.npcId == n.id);
        std::string label = n.name + " (" + npcRoleName(n.role) + ")";
        if (ImGui::Selectable((label + "##npc").c_str(), selected)) {
          navTarget.type = NavTargetType::Npc;
          navTarget.npcId = n.id;
        }
      }
      ImGui::Separator();
      ImGui::Checkbox("Autopilot (P)", &autopilot);
      ImGui::Checkbox("Supercruise Assist (Z)", &scAssist);
      ImGui::Checkbox("Fuel Scoop (O)", &scoopEnabled);
      ImGui::End();
    }

    // Station services: market / repairs / refuel / upgrades
    if (showMarket) {
      ImGui::Begin("Station Services", &showMarket);

      if (!isDocked()) {
        ImGui::TextColored({1, 0.6f, 0.2f, 1}, "Dock at a station to access services.");
        ImGui::End();
      } else {
        auto stOpt = findStationById(*currentSys, dockedStationId);
        if (!stOpt) {
          ImGui::Text("Docked station not found.");
          ImGui::End();
        } else {
          const auto& st = **stOpt;
          ImGui::Text("%s (%s)", st.name.c_str(), stationTypeName(st.type));
          ImGui::Text("Fee: %.1f%%", st.feeRate * 100.0);

          auto& econState = universe.stationEconomy(st, save.timeDays);

          if (ImGui::CollapsingHeader("Market", ImGuiTreeNodeFlags_DefaultOpen)) {
            if (ImGui::BeginTable("market", 7, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
              ImGui::TableSetupColumn("Commodity");
              ImGui::TableSetupColumn("Inv");
              ImGui::TableSetupColumn("Bid");
              ImGui::TableSetupColumn("Ask");
              ImGui::TableSetupColumn("You");
              ImGui::TableSetupColumn("Buy");
              ImGui::TableSetupColumn("Sell");
              ImGui::TableHeadersRow();

              for (std::size_t i = 0; i < stellar::econ::kCommodityCount; ++i) {
                const auto id = static_cast<stellar::econ::CommodityId>(i);
                const auto q = stellar::econ::quote(econState, st.economyModel, id);

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::TextUnformatted(stellar::econ::commodityName(id).data());
                ImGui::TableSetColumnIndex(1);
                ImGui::Text("%.1f", q.inventory);
                ImGui::TableSetColumnIndex(2);
                ImGui::Text("%.1f", q.bid);
                ImGui::TableSetColumnIndex(3);
                ImGui::Text("%.1f", q.ask);
                ImGui::TableSetColumnIndex(4);
                ImGui::Text("%.1f", save.cargo[i]);

                // Buy 1
                ImGui::TableSetColumnIndex(5);
                std::string buyId = std::string("Buy##") + std::to_string(i);
                if (ImGui::SmallButton(buyId.c_str())) {
                  const double mass = stellar::econ::commodityDef(id).massKg;
                  if (cargoMassKg(save.cargo) + mass > cargoCapacityKg) {
                    pushToast(toasts, "Cargo hold full.");
                  } else {
                    auto tr = stellar::econ::buy(econState, st.economyModel, id, 1.0, save.credits, 0.10, st.feeRate);
                    if (!tr.ok) pushToast(toasts, tr.reason ? tr.reason : "Buy failed");
                    else save.cargo[i] += 1.0;
                  }
                }

                // Sell 1
                ImGui::TableSetColumnIndex(6);
                std::string sellId = std::string("Sell##") + std::to_string(i);
                if (ImGui::SmallButton(sellId.c_str())) {
                  if (save.cargo[i] < 1.0) {
                    pushToast(toasts, "You don't have that cargo.");
                  } else {
                    auto tr = stellar::econ::sell(econState, st.economyModel, id, 1.0, save.credits, 0.10, st.feeRate);
                    if (!tr.ok) pushToast(toasts, tr.reason ? tr.reason : "Sell failed");
                    else save.cargo[i] -= 1.0;
                  }
                }
              }

              ImGui::EndTable();
            }
          }

          if (ImGui::CollapsingHeader("Refuel & Repairs", ImGuiTreeNodeFlags_DefaultOpen)) {
            // Refuel uses station Fuel inventory.
            const auto fuelQuote = stellar::econ::quote(econState, st.economyModel, stellar::econ::CommodityId::Fuel);
            const double need = std::max(0.0, fuelMax - fuel);
            const double maxBuy = std::min(need, fuelQuote.inventory);

            ImGui::Text("Fuel price: ask %.1f  (station inv %.1f)", fuelQuote.ask, fuelQuote.inventory);
            if (ImGui::Button("Refuel 5 units")) {
              const double amount = std::min(5.0, maxBuy);
              if (amount <= 0.0) {
                pushToast(toasts, "No fuel available.");
              } else {
                double credits = save.credits;
                auto tr = stellar::econ::buy(econState, st.economyModel, stellar::econ::CommodityId::Fuel, amount, credits, 0.10, st.feeRate);
                if (!tr.ok) pushToast(toasts, tr.reason ? tr.reason : "Refuel failed");
                else {
                  save.credits = credits;
                  fuel = std::min(fuelMax, fuel + amount);
                  pushToast(toasts, "Refueled.");
                }
              }
            }

            // Repairs consume Metals from station to gate repair capability.
            const auto metQuote = stellar::econ::quote(econState, st.economyModel, stellar::econ::CommodityId::Metals);
            const double missing = std::max(0.0, 1.0 - hull);
            const double metalPerHull = 6.0; // units per 100% hull
            const double needMetal = missing * metalPerHull;
            const double canMetal = std::min(needMetal, metQuote.inventory);

            ImGui::Text("Metals inv: %.1f (repair needs %.1f)", metQuote.inventory, needMetal);
            if (ImGui::Button("Repair 10%")) {
              const double wantHull = 0.10;
              const double needM = wantHull * metalPerHull;
              if (hull >= 0.999) {
                pushToast(toasts, "Hull already full.");
              } else if (metQuote.inventory < needM) {
                pushToast(toasts, "Not enough Metals for repairs.");
              } else {
                // Consume metals.
                double dummyCredits = 1e12;
                auto tr = stellar::econ::buy(econState, st.economyModel, stellar::econ::CommodityId::Metals, needM, dummyCredits, 0.0, st.feeRate);
                if (!tr.ok) {
                  pushToast(toasts, "Repair inventory error.");
                } else {
                  const double cost = 120.0 * wantHull * (1.0 + st.feeRate * 2.0);
                  if (save.credits < cost) {
                    // Refund metals if player can't pay
                    (void)stellar::econ::sell(econState, st.economyModel, stellar::econ::CommodityId::Metals, needM, dummyCredits, 0.0, st.feeRate);
                    pushToast(toasts, "Insufficient credits.");
                  } else {
                    save.credits -= cost;
                    hull = std::min(1.0, hull + wantHull);
                    pushToast(toasts, "Repaired hull.");
                  }
                }
              }
            }
          }

          if (ImGui::CollapsingHeader("Upgrades")) {
            ImGui::Text("(Early placeholder upgrade loop)");

            // Upgrades are now lightly gated by station inventory so scarcity matters.
            {
              const double cost = 650.0 * (1.0 + st.feeRate);
              const double parts = 6.0;
              const auto partId = stellar::econ::CommodityId::Machinery;
              const auto q = stellar::econ::quote(econState, st.economyModel, partId);

              ImGui::Text("Cargo racks require %.0f x Machinery (inv %.1f)", parts, q.inventory);
              if (ImGui::Button("Buy cargo racks (+80 kg)")) {
                if (save.credits < cost) {
                  pushToast(toasts, "Not enough credits.");
                } else if (q.inventory < parts) {
                  pushToast(toasts, "Station lacks Machinery parts.");
                } else {
                  double dummyCredits = 1e12;
                  auto tr = stellar::econ::buy(econState, st.economyModel, partId, parts, dummyCredits, 0.0, st.feeRate);
                  if (!tr.ok) {
                    pushToast(toasts, "Upgrade failed (inventory).");
                  } else {
                    save.credits -= cost;
                    cargoCapacityKg += 80.0;
                    pushToast(toasts, "Cargo capacity upgraded.");
                  }
                }
              }
            }

            {
              const double cost = 550.0 * (1.0 + st.feeRate);
              const double parts = 5.0;
              const auto partId = stellar::econ::CommodityId::Metals;
              const auto q = stellar::econ::quote(econState, st.economyModel, partId);
              ImGui::Text("Aux tank requires %.0f x Metals (inv %.1f)", parts, q.inventory);

              if (ImGui::Button("Buy auxiliary tank (+20 fuel)")) {
                if (save.credits < cost) {
                  pushToast(toasts, "Not enough credits.");
                } else if (q.inventory < parts) {
                  pushToast(toasts, "Station lacks Metals parts.");
                } else {
                  double dummyCredits = 1e12;
                  auto tr = stellar::econ::buy(econState, st.economyModel, partId, parts, dummyCredits, 0.0, st.feeRate);
                  if (!tr.ok) {
                    pushToast(toasts, "Upgrade failed (inventory).");
                  } else {
                    save.credits -= cost;
                    fuelMax += 20.0;
                    fuel = std::min(fuelMax, fuel + 10.0);
                    pushToast(toasts, "Fuel tank upgraded.");
                  }
                }
              }
            }
          }

          ImGui::End();
        }
      }
    }

    // Missions panel
    if (showMissions) {
      ImGui::Begin("Missions", &showMissions);

      if (!isDocked()) {
        ImGui::TextColored({1, 0.6f, 0.2f, 1}, "Dock to pick up missions.");
      } else {
        auto stOpt = findStationById(*currentSys, dockedStationId);
        if (stOpt) {
          const auto& st = **stOpt;
          const std::uint64_t dayKey = static_cast<std::uint64_t>(std::floor(save.timeDays));
          if (boardStationId != st.id || boardDayKey != dayKey) {
            boardStationId = st.id;
            boardDayKey = dayKey;
            boardMissions = generateMissionBoard(currentStub, *currentSys, st, universe, save.timeDays, 6, 35.0);
          }

          ImGui::Text("Available:");
          for (std::size_t i = 0; i < boardMissions.size(); ++i) {
            auto& m = boardMissions[i];

            // Gate delivery missions by stock.
            bool canAccept = true;
            if (m.type == stellar::sim::MissionType::Delivery && m.cargoProvided) {
              auto& econState = universe.stationEconomy(st, save.timeDays);
              const double inv = econState.inventory[static_cast<std::size_t>(m.commodity)];
              if (inv < m.units) canAccept = false;
              if (cargoMassKg(save.cargo) + m.units * stellar::econ::commodityDef(m.commodity).massKg > cargoCapacityKg) canAccept = false;
            }

            ImGui::Separator();
            ImGui::Text("%s -> station %llu", missionTypeName(m.type), static_cast<unsigned long long>(m.toStation));
            if (m.type == stellar::sim::MissionType::Delivery) {
              ImGui::Text("Cargo: %.0f x %s", m.units, stellar::econ::commodityName(m.commodity).data());
            }
            ImGui::Text("Reward: %.0f   Deadline: day %.1f", m.reward, m.deadlineDay);

            ImGui::BeginDisabled(!canAccept);
            std::string btn = std::string("Accept##") + std::to_string(i);
            if (ImGui::Button(btn.c_str())) {
              stellar::sim::Mission a = m;
              a.id = nextMissionId++;

              if (a.type == stellar::sim::MissionType::Delivery && a.cargoProvided) {
                // Load cargo from station inventory right now.
                auto& econState = universe.stationEconomy(st, save.timeDays);
                double dummyCredits = 1e12;
                auto tr = stellar::econ::buy(econState, st.economyModel, a.commodity, a.units, dummyCredits, 0.0, st.feeRate);
                if (!tr.ok) {
                  pushToast(toasts, "Mission cargo unavailable.");
                } else {
                  save.cargo[static_cast<std::size_t>(a.commodity)] += a.units;
                  pushToast(toasts, "Mission accepted: delivery cargo loaded.");
                  activeMissions.push_back(std::move(a));
                }
              } else if (a.type == stellar::sim::MissionType::BountyScan) {
                // Spawn a deterministic target id.
                a.targetNpcId = stellar::core::hashCombine(a.id, 0xB00B1EULL);

                // Place the target in the destination system on arrival.
                // For now we spawn it immediately if target is current system.
                if (a.toSystem == save.currentSystem) {
                  NpcShip tgt{};
                  tgt.id = a.targetNpcId;
                  tgt.role = NpcRole::Pirate;
                  tgt.name = "Wanted " + makeNpcName(tgt.id, NpcRole::Pirate);
                  tgt.state = NpcState::Loiter;
                  tgt.posKm = ship.positionKm() + stellar::math::Vec3d{1200.0, 0.0, 1200.0};
                  npcs.push_back(std::move(tgt));
                }

                pushToast(toasts, "Mission accepted: bounty scan.");
                activeMissions.push_back(std::move(a));
              } else {
                pushToast(toasts, "Mission accepted.");
                activeMissions.push_back(std::move(a));
              }
            }
            ImGui::EndDisabled();
            if (!canAccept) {
              ImGui::TextColored({1, 0.5f, 0.3f, 1}, "Unavailable (scarcity or cargo capacity).");
            }
          }
        }
      }

      ImGui::Separator();
      ImGui::Text("Active:");
      for (auto& m : activeMissions) {
        if (m.failed) {
          ImGui::TextColored({1, 0.3f, 0.3f, 1}, "FAILED: %s (id %llu)", missionTypeName(m.type), static_cast<unsigned long long>(m.id));
        } else if (m.completed) {
          ImGui::TextColored({0.3f, 1, 0.4f, 1}, "DONE: %s (id %llu)", missionTypeName(m.type), static_cast<unsigned long long>(m.id));
        } else {
          ImGui::Text("%s (id %llu) -> station %llu (deadline %.1f)",
                      missionTypeName(m.type),
                      static_cast<unsigned long long>(m.id),
                      static_cast<unsigned long long>(m.toStation),
                      m.deadlineDay);
        }
      }

      ImGui::End();
    }

    // Galaxy map / FSD jump
    if (showGalaxy) {
      ImGui::Begin("Galaxy Map", &showGalaxy);

      ImGui::Text("Jump range: %.1f ly (fuel %.1f, cost/ly %.2f)", maxJumpRangeLy(), fuel, jumpFuelPerLy());
      const bool fsdReady = (save.timeDays >= fsdReadyDay);
      if (!fsdReady) {
        ImGui::TextColored({1, 0.6f, 0.2f, 1}, "FSD cooldown... ready at day %.2f", fsdReadyDay);
      }

      // Nearby systems around current.
      auto sysList = universe.queryNearby(currentStub.posLy, 60.0, 64);
      if (sysList.empty()) {
        ImGui::Text("No nearby systems.");
      } else {
        if (ImGui::BeginTable("gal", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, {0, 300})) {
          ImGui::TableSetupColumn("Name");
          ImGui::TableSetupColumn("Dist (ly)");
          ImGui::TableSetupColumn("Planets");
          ImGui::TableSetupColumn("Stations");
          ImGui::TableHeadersRow();

          for (const auto& s : sysList) {
            if (s.id == currentStub.id) continue;
            const double distLy = (s.posLy - currentStub.posLy).length();
            const bool inRange = distLy <= maxJumpRangeLy();

            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            const bool selected = (galaxySelectedSystem == s.id);
            ImGui::BeginDisabled(!inRange);
            if (ImGui::Selectable((s.name + "##sys").c_str(), selected, ImGuiSelectableFlags_SpanAllColumns)) {
              galaxySelectedSystem = s.id;
            }
            ImGui::EndDisabled();
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%.1f", distLy);
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%d", s.planetCount);
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%d", s.stationCount);
          }

          ImGui::EndTable();
        }

        if (galaxySelectedSystem != 0) {
          auto dstStubOpt = [&]() -> std::optional<stellar::sim::SystemStub> {
            for (const auto& s : sysList) if (s.id == galaxySelectedSystem) return s;
            return std::nullopt;
          }();

          if (dstStubOpt) {
            const double distLy = (dstStubOpt->posLy - currentStub.posLy).length();
            const double fuelNeed = distLy * jumpFuelPerLy();

            ImGui::Separator();
            ImGui::Text("Selected: %s (%.1f ly)", dstStubOpt->name.c_str(), distLy);
            ImGui::Text("Fuel needed: %.1f", fuelNeed);

            ImGui::BeginDisabled(!(isDocked() && fsdReady && distLy <= maxJumpRangeLy() && fuel >= fuelNeed));
            if (ImGui::Button("Hyperspace Jump")) {
              fuel -= fuelNeed;
              fsdReadyDay = save.timeDays + (3.0 / 1440.0); // 3 min cooldown

              // Simple "jump" (no animation yet).
              arriveInSystemNearStation(dstStubOpt->id, *dstStubOpt);

              // If any bounty mission targets the destination, spawn its target contact.
              for (const auto& m : activeMissions) {
                if (m.type == stellar::sim::MissionType::BountyScan && !m.completed && !m.failed && m.toSystem == dstStubOpt->id && m.targetNpcId != 0) {
                  NpcShip tgt{};
                  tgt.id = m.targetNpcId;
                  tgt.role = NpcRole::Pirate;
                  tgt.name = "Wanted " + makeNpcName(tgt.id, NpcRole::Pirate);
                  tgt.state = NpcState::Loiter;
                  tgt.posKm = ship.positionKm() + stellar::math::Vec3d{1500.0, 0.0, 1000.0};
                  npcs.push_back(std::move(tgt));
                }
              }

              pushToast(toasts, "FSD jump complete.");
            }
            ImGui::EndDisabled();

            if (!isDocked()) {
              ImGui::TextColored({1, 0.6f, 0.2f, 1}, "(Dock to jump)" );
            }
          }
        }
      }

      ImGui::End();
    }

    // Simple on-screen marker for current target
    if (navTarget.type != NavTargetType::None) {
      if (auto tPos = targetWorldPosKm()) {
        bool behind = false;
        const ImVec2 pt = worldToScreen(view, proj, (*tPos) * worldScale, vpSize, &behind);
        if (!behind && pt.x >= -100 && pt.x <= vpSize.x + 100 && pt.y >= -100 && pt.y <= vpSize.y + 100) {
          auto* dl = ImGui::GetForegroundDrawList();
          dl->AddCircle(pt, 8.0f, IM_COL32(255, 220, 80, 220), 16, 2.0f);
          dl->AddText({pt.x + 10.0f, pt.y - 10.0f}, IM_COL32(255, 220, 80, 220), targetName().c_str());
        }
      }
    }

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    SDL_GL_SwapWindow(window);

    // Autosave
    static double saveTimer = 0.0;
    saveTimer += dt;
    if (saveTimer > 2.0) {
      saveTimer = 0.0;
      (void)stellar::sim::saveToFile(save, kSavePath);
    }
  }

  // Final save
  (void)stellar::sim::saveToFile(save, kSavePath);

  // Shutdown
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplSDL2_Shutdown();
  ImGui::DestroyContext();

  SDL_GL_DeleteContext(gl_ctx);
  SDL_DestroyWindow(window);
  SDL_Quit();

  return 0;
}
