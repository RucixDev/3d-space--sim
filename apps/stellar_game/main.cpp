// On Windows, SDL may redefine main() to SDL_main unless SDL2main is linked.
// We provide our own entry point, so prevent SDL from overriding it.
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define SDL_MAIN_HANDLED
#endif

#include "stellar/core/Log.h"
#include "stellar/core/Random.h"
#include "stellar/core/Hash.h"
#include "stellar/math/Vec2.h"
#include "stellar/proc/Noise.h"
#include "stellar/econ/Market.h"
#include "stellar/econ/RoutePlanner.h"
#include "stellar/render/Camera.h"
#include "stellar/render/Gl.h"
#include "stellar/render/LineRenderer.h"
#include "stellar/render/PointRenderer.h"
#include "stellar/render/Mesh.h"
#include "stellar/render/MeshRenderer.h"
#include "stellar/render/Texture.h"
#include "stellar/render/ProceduralSprite.h"
#include "stellar/render/ProceduralPlanet.h"
#include "stellar/render/ProceduralLivery.h"
#include "stellar/render/RenderTarget.h"
#include "stellar/render/AtmosphereRenderer.h"
#include "stellar/render/Starfield.h"
#include "stellar/render/Nebula.h"
#include "stellar/render/ParticleSystem.h"
#include "stellar/render/PostFX.h"
#include "stellar/ui/HudLayout.h"
#include "stellar/ui/HudSettings.h"
#include "stellar/ui/Livery.h"
#include "stellar/ui/UiSettings.h"
#include "stellar/ui/Bookmarks.h"
#include "stellar/ui/FuzzySearch.h"
#include "stellar/sim/Orbit.h"
#include "stellar/sim/MissionLogic.h"
#include "stellar/sim/Contraband.h"
#include "stellar/sim/Law.h"
#include "stellar/sim/PoliceScan.h"
#include "stellar/sim/Distress.h"
#include "stellar/sim/ResourceField.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/WorldIds.h"
#include "stellar/sim/Warehouse.h"
#include "stellar/sim/IndustryService.h"
#include "stellar/sim/TradeScanner.h"
#include "stellar/sim/Ship.h"
#include "stellar/sim/ShipLoadout.h"
#include "stellar/sim/Combat.h"
#include "stellar/sim/Ballistics.h"
#include "stellar/sim/PowerDistributor.h"
#include "stellar/sim/FlightController.h"
#include "stellar/sim/SupercruiseComputer.h"
#include "stellar/sim/Interdiction.h"
#include "stellar/sim/EncounterDirector.h"
#include "stellar/sim/Docking.h"
#include "stellar/sim/DockingComputer.h"
#include "stellar/sim/DockingClearanceService.h"
#include "stellar/sim/ThermalSystem.h"
#include "stellar/sim/Traffic.h"
#include "stellar/sim/NavRoute.h"
#include "stellar/sim/Universe.h"

#include "ControlsConfig.h"
#include "CommandPalette.h"
#include "ControlsWindow.h"
#include "ConsoleWindow.h"

#include <SDL.h>
#include <SDL_opengl.h>

// Some SDL configurations still define main -> SDL_main; ensure our main symbol remains intact.
#ifdef main
#undef main
#endif

#include <imgui.h>
#ifdef IMGUI_HAS_DOCK
#include <imgui_internal.h>
#endif
#include "stellar/ui/ImGuiCompat.h"
#include <backends/imgui_impl_opengl3.h>
#include <backends/imgui_impl_sdl2.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <cstdio>
#include <cctype>
#include <cstdint>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <optional>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace stellar;

// Convenience aliases for frequently used gameplay enums.
// (The game prototype uses these heavily; keep the names short in this TU.)
using sim::ShipHullClass;
using sim::WeaponType;

static constexpr double kAU_KM = 149597870.7;
static constexpr double kSOLAR_RADIUS_KM = 695700.0;
static constexpr double kEARTH_RADIUS_KM = 6371.0;

// Rendering scale:
// The sim uses kilometers. For rendering, we scale down by this factor.
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

static void starClassRgb(sim::StarClass cls, float& r, float& g, float& b) {
  auto set255 = [&](int rr, int gg, int bb) {
    r = rr / 255.0f;
    g = gg / 255.0f;
    b = bb / 255.0f;
  };
  switch (cls) {
    case sim::StarClass::O: set255(120, 180, 255); break;
    case sim::StarClass::B: set255(155, 205, 255); break;
    case sim::StarClass::A: set255(240, 240, 255); break;
    case sim::StarClass::F: set255(255, 248, 220); break;
    case sim::StarClass::G: set255(255, 236, 175); break;
    case sim::StarClass::K: set255(255, 200, 140); break;
    case sim::StarClass::M: set255(255, 165, 140); break;
    default: set255(200, 200, 220); break;
  }
}
static const char* planetTypeName(sim::PlanetType t) {
  switch (t) {
    case sim::PlanetType::Rocky:    return "Rocky";
    case sim::PlanetType::Desert:   return "Desert";
    case sim::PlanetType::Ocean:    return "Ocean";
    case sim::PlanetType::Ice:      return "Ice";
    case sim::PlanetType::GasGiant: return "Gas Giant";
    default: return "Unknown";
  }
}

static render::SurfaceKind planetSurfaceKind(sim::PlanetType t) {
  switch (t) {
    case sim::PlanetType::Rocky:    return render::SurfaceKind::Rocky;
    case sim::PlanetType::Desert:   return render::SurfaceKind::Desert;
    case sim::PlanetType::Ocean:    return render::SurfaceKind::Ocean;
    case sim::PlanetType::Ice:      return render::SurfaceKind::Ice;
    case sim::PlanetType::GasGiant: return render::SurfaceKind::GasGiant;
    default: return render::SurfaceKind::Rocky;
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

static double cargoMassKg(const std::array<double, econ::kCommodityCount>& cargo) {
  double kg = 0.0;
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const auto cid = (econ::CommodityId)i;
    const double units = cargo[i];
    if (units <= 0.0) continue;
    kg += units * econ::commodityDef(cid).massKg;
  }
  return kg;
}

static double clampRep(double r) {
  return std::clamp(r, -100.0, 100.0);
}

static double repNorm(double r) {
  return clampRep(r) / 100.0;
}

static double applyRepToFee(double baseFeeRate, double rep) {
  // Positive rep reduces fees, negative rep increases fees.
  const double kMaxEffect = 0.35; // +/- 35%
  const double eff = baseFeeRate * (1.0 - kMaxEffect * repNorm(rep));
  return std::clamp(eff, 0.0, 0.25);
}

static math::Vec3d toRenderU(const math::Vec3d& km) { return km * (1.0 / kRENDER_UNIT_KM); }

static math::Quatd quatFromTo(const math::Vec3d& from, const math::Vec3d& to) {
  math::Vec3d f = from.normalized();
  math::Vec3d t = to.normalized();
  const double c = std::clamp(math::dot(f, t), -1.0, 1.0);

  if (c > 0.999999) return math::Quatd::identity();

  if (c < -0.999999) {
    math::Vec3d axis = math::cross({1,0,0}, f);
    if (axis.lengthSq() < 1e-12) axis = math::cross({0,1,0}, f);
    return math::Quatd::fromAxisAngle(axis, math::kPi);
  }

  const math::Vec3d axis = math::cross(f, t);
  const double s = std::sqrt((1.0 + c) * 2.0);
  const double invs = 1.0 / s;
  return math::Quatd{ s * 0.5, axis.x * invs, axis.y * invs, axis.z * invs }.normalized();
}

static math::Vec3d stationPosKm(const sim::Station& st, double timeDays) {
  const math::Vec3d posAU = sim::orbitPosition3DAU(st.orbit, timeDays);
  return posAU * kAU_KM;
}

static math::Vec3d stationVelKmS(const sim::Station& st, double timeDays) {
  // Analytic Keplerian velocity (avoids numerical jitter in local-frame velocity matching).
  const math::Vec3d velAUPerDay = sim::orbitVelocity3DAU(st.orbit, timeDays);
  return velAUPerDay * (kAU_KM / 86400.0);
}

static math::Vec3d planetPosKm(const sim::Planet& p, double timeDays) {
  const math::Vec3d posAU = sim::orbitPosition3DAU(p.orbit, timeDays);
  return posAU * kAU_KM;
}

static math::Vec3d planetVelKmS(const sim::Planet& p, double timeDays) {
  // Analytic Keplerian velocity.
  const math::Vec3d velAUPerDay = sim::orbitVelocity3DAU(p.orbit, timeDays);
  return velAUPerDay * (kAU_KM / 86400.0);
}

static math::Quatd stationOrient(const sim::Station& st, const math::Vec3d& posKm, double timeDays) {
  // Station local +Z points outward from the slot.
  // We point it away from the star (radial outward) for intuitive docking.
  math::Vec3d outward = posKm.normalized();
  if (outward.lengthSq() < 1e-12) outward = {0,0,1};
  math::Quatd base = quatFromTo({0,0,1}, outward);

  // Add a small spin to make orientation/geometry feel alive.
  // Spin around the station's local +Z axis.
  const double spinRadPerDay = math::degToRad(35.0);
  math::Quatd spin = math::Quatd::fromAxisAngle({0,0,1}, std::fmod(timeDays * spinRadPerDay, 2.0 * math::kPi));
  return (base * spin).normalized();
}

static render::InstanceData makeInst(const math::Vec3d& posU,
                                    const math::Vec3d& scaleU,
                                    const math::Quatd& qWorld,
                                    float cr, float cg, float cb) {
  // shader expects quat as x,y,z,w
  const math::Quatd q = qWorld.normalized();
  return render::InstanceData{
    (float)posU.x, (float)posU.y, (float)posU.z,
    (float)scaleU.x, (float)scaleU.y, (float)scaleU.z,
    (float)q.x, (float)q.y, (float)q.z, (float)q.w,
    cr, cg, cb
  };
}

static render::InstanceData makeInstUniform(const math::Vec3d& posU,
                                           double s,
                                           float cr, float cg, float cb) {
  return makeInst(posU, {s,s,s}, math::Quatd::identity(), cr,cg,cb);
}

static bool projectToScreen(const math::Vec3d& worldU,
                            const math::Mat4d& view,
                            const math::Mat4d& proj,
                            int w, int h,
                            ImVec2& outPx) {
  // clip = proj * view * vec4(world, 1)
  const math::Mat4d vp = proj * view;
  const double x = worldU.x, y = worldU.y, z = worldU.z;

  const double cx = vp.m[0]*x + vp.m[4]*y + vp.m[8]*z + vp.m[12];
  const double cy = vp.m[1]*x + vp.m[5]*y + vp.m[9]*z + vp.m[13];
  const double cz = vp.m[2]*x + vp.m[6]*y + vp.m[10]*z + vp.m[14];
  const double cw = vp.m[3]*x + vp.m[7]*y + vp.m[11]*z + vp.m[15];

  if (cw <= 1e-6) return false;

  const double ndcX = cx / cw;
  const double ndcY = cy / cw;
  // const double ndcZ = cz / cw;

  // off-screen (still allow a bit of slack)
  if (ndcX < -1.2 || ndcX > 1.2 || ndcY < -1.2 || ndcY > 1.2) return false;

  outPx.x = (float)((ndcX * 0.5 + 0.5) * (double)w);
  outPx.y = (float)((-ndcY * 0.5 + 0.5) * (double)h);
  return true;
}

// Like projectToScreen, but always returns a screen position (can be off-screen) as long as
// the point is in front of the camera. Useful for edge-of-screen indicators.
static bool projectToScreenAny(const math::Vec3d& worldU,
                               const math::Mat4d& view,
                               const math::Mat4d& proj,
                               int w, int h,
                               ImVec2& outPx,
                               bool& outOffscreen) {
  const math::Mat4d vp = proj * view;
  const double x = worldU.x, y = worldU.y, z = worldU.z;

  const double cx = vp.m[0]*x + vp.m[4]*y + vp.m[8]*z + vp.m[12];
  const double cy = vp.m[1]*x + vp.m[5]*y + vp.m[9]*z + vp.m[13];
  const double cw = vp.m[3]*x + vp.m[7]*y + vp.m[11]*z + vp.m[15];
  if (cw <= 1e-6) return false;

  const double ndcX = cx / cw;
  const double ndcY = cy / cw;

  outOffscreen = (ndcX < -1.0 || ndcX > 1.0 || ndcY < -1.0 || ndcY > 1.0);
  outPx.x = (float)((ndcX * 0.5 + 0.5) * (double)w);
  outPx.y = (float)((-ndcY * 0.5 + 0.5) * (double)h);
  return true;
}

static bool segmentHitsSphere(const math::Vec3d& aKm,
                             const math::Vec3d& bKm,
                             const math::Vec3d& centerKm,
                             double radiusKm) {
  // Share the core collision helper (extracted from main.cpp in Round 5).
  return sim::segmentHitsSphereKm(aKm, bKm, centerKm, radiusKm);
}

static render::Texture2D makeRadialSpriteTextureRGBA(int size,
                                                     float coreExp,
                                                     float haloExp,
                                                     float haloStrength,
                                                     float spikeStrength,
                                                     int spikes,
                                                     core::u64 seed) {
  // A tiny procedural point-sprite used by the PointRenderer's textured mode.
  // The result is a white-ish sprite; point color is applied via vertex tint.
  std::vector<std::uint8_t> rgba;
  rgba.resize((std::size_t)size * (std::size_t)size * 4);

  auto hash01 = [&](int x, int y) -> float {
    // Cheap per-pixel deterministic noise in [0,1].
    core::u64 h = seed;
    h = core::hashCombine(h, (core::u64)(x * 73856093));
    h = core::hashCombine(h, (core::u64)(y * 19349663));
    h ^= (h >> 33);
    h *= 0xff51afd7ed558ccdULL;
    h ^= (h >> 33);
    // Take high bits.
    const core::u64 v = (h >> 40) & 0xFFFFULL;
    return (float)v / 65535.0f;
  };

  const float inv = (size > 1) ? (1.0f / (float)(size - 1)) : 1.0f;
  for (int y = 0; y < size; ++y) {
    for (int x = 0; x < size; ++x) {
      const float fx = (float)x * inv;
      const float fy = (float)y * inv;
      const float u = fx * 2.0f - 1.0f;
      const float v = fy * 2.0f - 1.0f;
      const float r2 = u*u + v*v;

      // Inside a radius slightly >1 so edges don't get clipped by smoothing.
      float a = 0.0f;
      if (r2 <= 1.25f) {
        const float core = std::exp(-r2 * coreExp);
        const float halo = std::exp(-r2 * haloExp) * haloStrength;
        a = core + halo;

        if (spikeStrength > 1e-4f && spikes > 0) {
          const float ang = std::atan2(v, u);
          const float s = std::abs(std::sin(ang * (float)spikes));
          const float sp = std::pow(s, 6.0f);
          a *= (1.0f - spikeStrength) + spikeStrength * sp;
        }

        // Subtle dithering so the sprite doesn't band when blurred.
        const float n = hash01(x, y) * 2.0f - 1.0f;
        a *= (1.0f + n * 0.06f);
      }

      a = std::clamp(a, 0.0f, 1.0f);
      const std::uint8_t alpha = (std::uint8_t)std::lround(a * 255.0f);

      const std::size_t i = ((std::size_t)y * (std::size_t)size + (std::size_t)x) * 4;
      rgba[i + 0] = 255;
      rgba[i + 1] = 255;
      rgba[i + 2] = 255;
      rgba[i + 3] = alpha;
    }
  }

  render::Texture2D out;
  out.createRGBA(size, size, rgba.data(), true, false, true);
  return out;
}

static render::Texture2D makeCloudSpriteTextureRGBA(int size,
                                                    core::u64 seed,
                                                    float falloffExp = 2.4f,
                                                    float noiseScale = 3.2f,
                                                    int octaves = 5) {
  // A soft nebula puff sprite (white RGB with noisy alpha). Tint is applied via vertex color.
  std::vector<std::uint8_t> rgba;
  rgba.resize((std::size_t)size * (std::size_t)size * 4);

  auto hash01 = [&](int x, int y) -> float {
    // Cheap per-pixel deterministic noise in [0,1].
    core::u64 h = seed;
    h = core::hashCombine(h, (core::u64)(x * 73856093));
    h = core::hashCombine(h, (core::u64)(y * 19349663));
    h ^= (h >> 33);
    h *= 0xff51afd7ed558ccdULL;
    h ^= (h >> 33);
    const core::u64 v = (h >> 40) & 0xFFFFULL;
    return (float)v / 65535.0f;
  };

  const float inv = (size > 1) ? (1.0f / (float)(size - 1)) : 1.0f;
  for (int y = 0; y < size; ++y) {
    for (int x = 0; x < size; ++x) {
      const float fx = (float)x * inv;
      const float fy = (float)y * inv;
      const float u = fx * 2.0f - 1.0f;
      const float v = fy * 2.0f - 1.0f;
      const float r2 = u*u + v*v;

      float a = 0.0f;
      if (r2 <= 1.35f) {
        // Radial falloff, then modulate with fBm noise to get cloudiness.
        const float base = std::exp(-r2 * falloffExp);
        const double n = proc::fbm2D(seed ^ 0xBADC0FFEE0DDF00DULL,
                                     (double)u * (double)noiseScale + 12.3,
                                     (double)v * (double)noiseScale - 7.1,
                                     octaves, 2.05, 0.52);
        const float nf = (float)std::clamp(n, 0.0, 1.0);
        // Bias noise so we get soft holes.
        const float cloud = std::pow(nf, 1.15f);
        a = base * (0.35f + 0.75f * cloud);

        // Subtle dither (prevents banding when blurred/bloomed).
        const float d = (hash01(x, y) * 2.0f - 1.0f) * 0.06f;
        a *= (1.0f + d);
      }

      a = std::clamp(a, 0.0f, 1.0f);
      const std::uint8_t alpha = (std::uint8_t)std::lround(a * 255.0f);

      const std::size_t i = ((std::size_t)y * (std::size_t)size + (std::size_t)x) * 4;
      rgba[i + 0] = 255;
      rgba[i + 1] = 255;
      rgba[i + 2] = 255;
      rgba[i + 3] = alpha;
    }
  }

  render::Texture2D out;
  out.createRGBA(size, size, rgba.data(), true, false, true);
  return out;
}

static bool beginStationSelectorHUD(const sim::StarSystem& sys, int& stationIndex, bool docked, sim::StationId dockedId) {
  bool changed = false;

  ImGui::Begin("Dock / Station");

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

    if (docked) {
      ImGui::Text("Docked at: %s", (st.id == dockedId) ? st.name.c_str() : "(other station)");
    } else {
      ImGui::TextDisabled("Not docked");
    }
  } else {
    ImGui::Text("No stations in system.");
  }

  ImGui::End();
  return changed;
}

struct ToastMsg {
  std::string text;
  double ttl{3.0};      // seconds remaining
  double ttlTotal{3.0}; // initial ttl (for fade)
};

struct ToastHistoryEntry {
  double timeDays{0.0};
  std::string text;
};

static std::deque<ToastHistoryEntry>* gToastHistorySink = nullptr;
static const double* gToastTimeDaysPtr = nullptr;

static void setToastHistorySink(std::deque<ToastHistoryEntry>* sink, const double* timeDaysPtr) {
  gToastHistorySink = sink;
  gToastTimeDaysPtr = timeDaysPtr;
}

static void toast(std::vector<ToastMsg>& toasts, std::string msg, double ttlSec=3.0) {
  if (ttlSec <= 0.0) ttlSec = 0.01;
  toasts.push_back({std::move(msg), ttlSec, ttlSec});

  if (gToastHistorySink && gToastTimeDaysPtr) {
    gToastHistorySink->push_back({*gToastTimeDaysPtr, toasts.back().text});
    static constexpr std::size_t kCap = 1200;
    while (gToastHistorySink->size() > kCap) gToastHistorySink->pop_front();
  }
}

using ClearanceState = sim::DockingClearanceState;

enum class TargetKind : int {
  None = 0,
  Station = 1,
  Planet = 2,
  Contact = 3,
  Star = 4,
  Cargo = 5,
  Asteroid = 6,
  Signal = 7,
};

struct Target {
  TargetKind kind{TargetKind::None};
  std::size_t index{0}; // index into the relevant list for kind (stations/planets/contacts/...) 
};

enum class ContactRole : int { Pirate=0, Trader=1, Police=2 };

static const char* contactRoleName(ContactRole r) {
  switch (r) {
    case ContactRole::Pirate: return "Pirate";
    case ContactRole::Trader: return "Trader";
    case ContactRole::Police: return "Police";
    default: return "?";
  }
}

struct Contact {
  core::u64 id{0};
  std::string name;
  sim::Ship ship{};

  ContactRole role{ContactRole::Pirate};
  core::u32 factionId{0}; // police/trader affiliation; 0 for pirates / independent
  std::size_t homeStationIndex{0}; // best-effort "patrol" anchor (index into system stations)

  bool missionTarget{false};

  // Group behavior (squads)
  core::u64 groupId{0};   // 0 = solo. Non-zero = part of a squad.
  core::u64 leaderId{0};  // 0 = squad leader (or solo). Otherwise, id of leader.

  // Optional behavior relationships
  // - followId: try to stay near this contact (used e.g. for convoy escorts)
  // - attackTargetId: prefer attacking this contact instead of the player (used for mission ambushes)
  core::u64 followId{0};
  core::u64 attackTargetId{0};

  // Short "under fire" window for basic behavior tuning (e.g., morale / evasion).
  double underFireUntilDays{0.0};

  // Behavior flags
  bool hostileToPlayer{false}; // police become hostile when you're wanted / shoot them
  double fleeUntilDays{0.0};   // traders flee after being attacked

  // --- Ship loadout + combat stats ---
  // NPCs use the same core loadout tables as the player.
  ShipHullClass hullClass{ShipHullClass::Scout};
  int thrusterMk{1};
  int shieldMk{1};
  int distributorMk{1};
  WeaponType weapon{WeaponType::BeamLaser};

  double shield{60.0};
  double shieldMax{60.0};
  double hull{70.0};
  double hullMax{70.0};

  // Derived from loadout (points per simulated minute, before pips multipliers).
  double shieldRegenPerSimMin{2.0};

  // Power distributor state (capacitors + pips).
  sim::Pips pips{2,2,2};
  sim::DistributorConfig distributorCfg{};
  sim::DistributorState distributorState{};

  // Basic AI tuning: 0..1 where higher = more accurate and more aggressive firing windows.
  double aiSkill{0.55};

  // Traders can have an approximate loot value (paid on destruction for now).
  double cargoValueCr{0.0};

  // Best-effort value of *current haul* used for ambush scaling.
  // For normal traders this can mirror `cargoValueCr`; for escort convoys it is
  // set from the mission payload at launch.
  double tradeCargoValueCr{0.0};

  // Escort-mission convoy marker (leader ship). Used only for UI/logic; not saved.
  bool escortConvoy{false};

  // --- Distress / rescue (signal encounters) ---
  // Distressed civilian ship that can be helped by delivering cargo.
  // Triggered via a normal contact scan (scanner action).
  bool distressVictim{false};
  core::u64 distressSignalId{0};
  econ::CommodityId distressNeedCommodity{econ::CommodityId::Food};
  double distressNeedUnits{0.0};
  double distressRewardCr{0.0};
  double distressRepReward{0.0};
  core::u32 distressPayerFactionId{0};

  // --- Trader economy simulation (local-system hauling) ---
  // Traders move a small amount of real station inventory between stations.
  // This makes markets feel "alive" and creates/relieves shortages.
  std::size_t tradeDestStationIndex{0};
  econ::CommodityId tradeCommodity{econ::CommodityId::Food};
  double tradeUnits{0.0};        // units currently carried
  double tradeCapacityKg{180.0}; // cargo capacity (kg)
  double tradeCooldownUntilDays{0.0};

  // Traders use a crude "supercruise" when far from the player to keep the in-system economy moving.
  double tradeSupercruiseSpeedKmS{9000.0};
  double tradeSupercruiseDropDistKm{120000.0};

  // Supercruise "shadows": contacts that are kept near the player during supercruise
  // so interdictions can happen even though the full combat sim is paused.
  bool supercruiseShadow{false};
  math::Vec3d shadowOffsetLocalKm{0,0,0};

  double fireCooldown{0.0}; // simulated seconds
  bool alive{true};
};

// --- Space objects (salvage / mining / encounters) ---
struct FloatingCargo {
  core::u64 id{0};
  econ::CommodityId commodity{econ::CommodityId::Food};
  double units{0.0};
  math::Vec3d posKm{0,0,0};
  math::Vec3d velKmS{0,0,0};
  double expireDay{0.0};
  bool fromPlayer{false};
  core::u64 missionId{0}; // non-zero if spawned for a mission objective
};

struct AsteroidNode {
  core::u64 id{0};
  math::Vec3d posKm{0,0,0};
  double radiusKm{2500.0};

  // What this asteroid yields when mined / prospected.
  econ::CommodityId yield{econ::CommodityId::Ore};

  // Units remaining inside the rock.
  double remainingUnits{120.0};

  // Mining can generate fractional chunks; accumulate and emit whole units over time.
  double chunkAccumulator{0.0};
};

enum class SignalType : int {
  Distress = 0,
  Derelict,
  Resource,
};

static const char* signalTypeName(SignalType t) {
  switch (t) {
    case SignalType::Distress: return "Distress Call";
    case SignalType::Derelict: return "Derelict";
    case SignalType::Resource: return "Resource Field";
    default: return "Signal";
  }
}

struct SignalSource {
  core::u64 id{0};
  SignalType type{SignalType::Distress};
  math::Vec3d posKm{0,0,0};
  double expireDay{0.0};

  // Deterministic plan for distress calls (victim request + optional ambush parameters).
  // Filled at spawn time so scan / resolution UI can display stable info.
  bool hasDistressPlan{false};
  sim::DistressPlan distress{};
  core::u64 distressVictimId{0};
  bool distressCompleted{false};

  // Whether the player has "resolved" this site (i.e., we have fired its one-shot content).
  bool resolved{false};

  // Some sites (resource fields) need a second one-shot flag even after 'resolved' for clarity / future expansion.
  bool fieldSpawned{false};

  // Deterministic plan for resource fields (composition/richness).
  // Filled at seed time so scanner / UI can show stable hints.
  bool hasResourcePlan{false};
  sim::ResourceFieldPlan resource{};
};

static void applyDamage(double dmg, double& shield, double& hull) {
  // Keep the prototype's call sites stable, but share the core damage rule.
  sim::applyDamage(dmg, shield, hull);
}

static void emitStationGeometry(const sim::Station& st,
                                const math::Vec3d& stPosKm,
                                const math::Quatd& stQ,
                                std::vector<render::InstanceData>& outCubeInstances) {
  // Render in render-units
  const math::Vec3d posU = toRenderU(stPosKm);

  // Station scale factor: exaggerate a bit so it's readable at prototype camera distances.
  const double s = std::max(0.8, (st.radiusKm / kRENDER_UNIT_KM) * 1800.0);
  const math::Vec3d baseScaleU = {s, s, s};

  auto addPart = [&](const math::Vec3d& localPosU, const math::Vec3d& localScaleU, float r, float g, float b) {
    const math::Vec3d worldPosU = posU + stQ.rotate(localPosU);
    const math::Quatd worldQ = stQ;
    outCubeInstances.push_back(makeInst(worldPosU, localScaleU, worldQ, r,g,b));
  };

  // Central body (elongated)
  addPart({0,0,0}, {baseScaleU.x*0.9, baseScaleU.y*0.9, baseScaleU.z*1.3}, 0.65f,0.68f,0.72f);

  // Ring (8 segments)
  const double ringR = baseScaleU.x * 1.25;
  const double segLen = baseScaleU.z * 0.55;
  for (int i = 0; i < 8; ++i) {
    const double ang = (double)i / 8.0 * 2.0 * math::kPi;
    const double cx = std::cos(ang), sx = std::sin(ang);
    const math::Vec3d lp{ ringR * cx, ringR * sx, 0.0 };
    addPart(lp, {baseScaleU.x*0.22, baseScaleU.y*0.22, segLen}, 0.45f,0.55f,0.70f);
  }

  // Docking "mail slot" frame on +Z face
  const double frameZ = baseScaleU.z * 1.3;
  const double frameW = baseScaleU.x * 1.15;
  const double frameH = baseScaleU.y * 0.60;
  const double thickness = baseScaleU.x * 0.12;

  // Left/right pylons
  addPart({-frameW*0.5, 0, frameZ}, {thickness, frameH, thickness}, 0.85f,0.80f,0.55f);
  addPart({+frameW*0.5, 0, frameZ}, {thickness, frameH, thickness}, 0.85f,0.80f,0.55f);

  // Top/bottom bars
  addPart({0, +frameH*0.5, frameZ}, {frameW, thickness, thickness}, 0.85f,0.80f,0.55f);
  addPart({0, -frameH*0.5, frameZ}, {frameW, thickness, thickness}, 0.85f,0.80f,0.55f);

  // A small "light strip" above slot
  addPart({0, +frameH*0.65, frameZ}, {frameW*0.8, thickness*0.6, thickness*0.6}, 0.25f,0.85f,0.35f);
}

int main(int argc, char** argv) {
  (void)argc; (void)argv;

  core::setLogLevel(core::LogLevel::Info);

#ifdef _WIN32
  SDL_SetMainReady();
#endif

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

  { // Scope GL objects so they are destroyed before SDL_GL_DeleteContext

  // --- ImGui setup ---
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
#ifdef IMGUI_HAS_DOCK
  io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
#endif
  // ---- UI settings (persistent) ----
  // Stored separately from Dear ImGui's imgui.ini so users can tweak look/scale/layout
  // without nuking window positions.
  const std::string uiSettingsPath = ui::defaultUiSettingsPath();
  ui::UiSettings uiSettings = ui::makeDefaultUiSettings();
  bool uiSettingsDirty = false;
  bool showUiSettingsWindow = false;

  (void)ui::loadUiSettingsFromFile(uiSettingsPath, uiSettings);

  // Allow users to keep multiple ImGui layout profiles by pointing ImGui at a
  // different ini file.
  std::string uiIniFilename = uiSettings.imguiIniFile;
  if (uiIniFilename.empty()) uiIniFilename = "imgui.ini";
  io.IniFilename = uiIniFilename.c_str();

  ui::UiTheme uiTheme = uiSettings.theme;

  auto rebuildUiStyle = [&](ui::UiTheme theme, float scale) {
    // Rebuild the style from scratch so switching themes doesn't leave behind
    // old border/rounding sizes, etc.
    ImGuiStyle st = ImGuiStyle();
    switch (theme) {
      case ui::UiTheme::Dark: ImGui::StyleColorsDark(&st); break;
      case ui::UiTheme::Light: ImGui::StyleColorsLight(&st); break;
      case ui::UiTheme::Classic: ImGui::StyleColorsClassic(&st); break;
      case ui::UiTheme::HighContrast: {
        // Start from Dark and push borders/frames toward higher contrast.
        ImGui::StyleColorsDark(&st);
        st.WindowRounding = 0.0f;
        st.FrameRounding = 0.0f;
        st.GrabRounding = 0.0f;
        st.ScrollbarRounding = 0.0f;
        st.WindowBorderSize = 1.0f;
        st.FrameBorderSize = 1.0f;

        ImVec4* c = st.Colors;
        c[ImGuiCol_Text] = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
        c[ImGuiCol_TextDisabled] = ImVec4(0.70f, 0.70f, 0.70f, 1.00f);
        c[ImGuiCol_WindowBg] = ImVec4(0.02f, 0.02f, 0.02f, 1.00f);
        c[ImGuiCol_ChildBg] = ImVec4(0.02f, 0.02f, 0.02f, 0.00f);
        c[ImGuiCol_PopupBg] = ImVec4(0.02f, 0.02f, 0.02f, 0.98f);
        c[ImGuiCol_Border] = ImVec4(0.80f, 0.80f, 0.80f, 0.30f);
        c[ImGuiCol_BorderShadow] = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);

        c[ImGuiCol_FrameBg] = ImVec4(0.10f, 0.10f, 0.10f, 1.00f);
        c[ImGuiCol_FrameBgHovered] = ImVec4(0.16f, 0.16f, 0.16f, 1.00f);
        c[ImGuiCol_FrameBgActive] = ImVec4(0.20f, 0.20f, 0.20f, 1.00f);

        c[ImGuiCol_TitleBg] = ImVec4(0.00f, 0.00f, 0.00f, 1.00f);
        c[ImGuiCol_TitleBgActive] = ImVec4(0.06f, 0.06f, 0.06f, 1.00f);
        c[ImGuiCol_TitleBgCollapsed] = ImVec4(0.00f, 0.00f, 0.00f, 0.80f);

        c[ImGuiCol_MenuBarBg] = ImVec4(0.05f, 0.05f, 0.05f, 1.00f);
        c[ImGuiCol_ScrollbarBg] = ImVec4(0.02f, 0.02f, 0.02f, 0.80f);
        c[ImGuiCol_ScrollbarGrab] = ImVec4(0.60f, 0.60f, 0.60f, 0.50f);
        c[ImGuiCol_ScrollbarGrabHovered] = ImVec4(0.75f, 0.75f, 0.75f, 0.70f);
        c[ImGuiCol_ScrollbarGrabActive] = ImVec4(0.90f, 0.90f, 0.90f, 0.80f);

        c[ImGuiCol_CheckMark] = ImVec4(0.95f, 0.95f, 0.95f, 1.00f);
        c[ImGuiCol_SliderGrab] = ImVec4(0.85f, 0.85f, 0.85f, 0.70f);
        c[ImGuiCol_SliderGrabActive] = ImVec4(0.95f, 0.95f, 0.95f, 0.90f);

        c[ImGuiCol_Button] = ImVec4(0.12f, 0.12f, 0.12f, 1.00f);
        c[ImGuiCol_ButtonHovered] = ImVec4(0.20f, 0.20f, 0.20f, 1.00f);
        c[ImGuiCol_ButtonActive] = ImVec4(0.28f, 0.28f, 0.28f, 1.00f);

        c[ImGuiCol_Header] = ImVec4(0.18f, 0.18f, 0.18f, 1.00f);
        c[ImGuiCol_HeaderHovered] = ImVec4(0.25f, 0.25f, 0.25f, 1.00f);
        c[ImGuiCol_HeaderActive] = ImVec4(0.30f, 0.30f, 0.30f, 1.00f);

        c[ImGuiCol_Separator] = ImVec4(0.85f, 0.85f, 0.85f, 0.25f);
        c[ImGuiCol_SeparatorHovered] = ImVec4(0.95f, 0.95f, 0.95f, 0.45f);
        c[ImGuiCol_SeparatorActive] = ImVec4(1.00f, 1.00f, 1.00f, 0.60f);

        c[ImGuiCol_Tab] = ImVec4(0.08f, 0.08f, 0.08f, 1.00f);
        c[ImGuiCol_TabHovered] = ImVec4(0.20f, 0.20f, 0.20f, 1.00f);
        c[ImGuiCol_TabActive] = ImVec4(0.14f, 0.14f, 0.14f, 1.00f);
      } break;
      default: ImGui::StyleColorsDark(&st); break;
    }

    if (std::abs(scale - 1.0f) > 0.001f) {
      st.ScaleAllSizes(scale);
    }
    ImGui::GetStyle() = st;
  };

  // ---- UI appearance / layout ----
  // NOTE: These are referenced during initial ImGui setup (HiDPI scaling + docking layout
  // seed) and then later during runtime UI menus, so they must live for the full app.
  bool uiAutoScaleFromDpi = uiSettings.autoScaleFromDpi;
  float uiScaleUser = uiSettings.scaleUser; // multiplier (not DPI-derived)
  float uiScaleDpi = 1.0f;                 // derived from SDL_GetDisplayDPI
  float uiScale = 1.0f;                    // effective scale applied to style + fonts
  float uiScaleApplied = 1.0f;             // internal: style already scaled to this value
  bool uiSettingsAutoSaveOnExit = uiSettings.autoSaveOnExit;
  bool showImGuiDemo = false;
  bool showImGuiMetrics = false;

#ifdef IMGUI_HAS_DOCK
  // Docking (requires Dear ImGui docking branch/tag).
  bool uiDockingEnabled = uiSettings.dock.dockingEnabled;
  bool uiDockPassthruCentral = uiSettings.dock.passthruCentral;  // keeps 3D view visible when central node is empty
  bool uiDockLockCentralView = uiSettings.dock.lockCentralView;  // prevents docking into the central node
  bool uiDockResetLayout = false;     // rebuild default layout next frame
  float uiDockLeftRatio = uiSettings.dock.leftRatio;
  float uiDockRightRatio = uiSettings.dock.rightRatio;
  float uiDockBottomRatio = uiSettings.dock.bottomRatio;
#endif

  auto recomputeUiDpiScale = [&]() {
    uiScaleDpi = 1.0f;
    if (!uiAutoScaleFromDpi) return;
    float ddpi = 0.0f;
    if (SDL_GetDisplayDPI(0, &ddpi, nullptr, nullptr) == 0 && ddpi > 0.0f) {
      uiScaleDpi = std::clamp(ddpi / 96.0f, 0.75f, 1.75f);
    }
  };

  auto recomputeUiScale = [&]() {
    const float dpiFactor = uiAutoScaleFromDpi ? uiScaleDpi : 1.0f;
    uiScale = std::clamp(uiScaleUser * dpiFactor, 0.75f, 2.50f);
  };

  // Initial UI scaling + theme.
  recomputeUiDpiScale();
  recomputeUiScale();
  rebuildUiStyle(uiTheme, uiScale);
  io.FontGlobalScale = uiScale;
  uiScaleApplied = uiScale;

  auto applyUiScaleNow = [&]() {
    recomputeUiScale();
    const float denom = std::max(0.001f, uiScaleApplied);
    ImGui::GetStyle().ScaleAllSizes(uiScale / denom);
    io.FontGlobalScale = uiScale;
    uiScaleApplied = uiScale;
  };

  auto syncUiSettingsFromRuntime = [&]() {
    uiSettings.imguiIniFile = uiIniFilename;
    uiSettings.autoScaleFromDpi = uiAutoScaleFromDpi;
    uiSettings.scaleUser = uiScaleUser;
    uiSettings.theme = uiTheme;
    uiSettings.autoSaveOnExit = uiSettingsAutoSaveOnExit;
#ifdef IMGUI_HAS_DOCK
    uiSettings.dock.dockingEnabled = uiDockingEnabled;
    uiSettings.dock.passthruCentral = uiDockPassthruCentral;
    uiSettings.dock.lockCentralView = uiDockLockCentralView;
    uiSettings.dock.leftRatio = uiDockLeftRatio;
    uiSettings.dock.rightRatio = uiDockRightRatio;
    uiSettings.dock.bottomRatio = uiDockBottomRatio;
#endif
  };

  auto applyUiSettingsToRuntime = [&](const ui::UiSettings& s, bool loadImGuiIni) {
    uiTheme = s.theme;
    uiAutoScaleFromDpi = s.autoScaleFromDpi;
    uiScaleUser = s.scaleUser;
    uiSettingsAutoSaveOnExit = s.autoSaveOnExit;

    uiIniFilename = s.imguiIniFile.empty() ? std::string("imgui.ini") : s.imguiIniFile;
    io.IniFilename = uiIniFilename.c_str();

#ifdef IMGUI_HAS_DOCK
    uiDockingEnabled = s.dock.dockingEnabled;
    uiDockPassthruCentral = s.dock.passthruCentral;
    uiDockLockCentralView = s.dock.lockCentralView;
    uiDockLeftRatio = s.dock.leftRatio;
    uiDockRightRatio = s.dock.rightRatio;
    uiDockBottomRatio = s.dock.bottomRatio;
#endif

    recomputeUiDpiScale();
    recomputeUiScale();
    rebuildUiStyle(uiTheme, uiScale);
    io.FontGlobalScale = uiScale;
    uiScaleApplied = uiScale;

    if (loadImGuiIni && io.IniFilename && io.IniFilename[0] != '\0') {
      ImGui::LoadIniSettingsFromDisk(io.IniFilename);
    }
  };

#ifdef IMGUI_HAS_DOCK
  // If this is the first run (no imgui.ini yet), seed a sensible default dock layout.
  if (uiDockingEnabled && io.IniFilename && io.IniFilename[0] != '\0') {
    std::ifstream ini(io.IniFilename);
    if (!ini.good())
      uiDockResetLayout = true;
  }
#endif

  ImGui_ImplSDL2_InitForOpenGL(window, glContext);
  ImGui_ImplOpenGL3_Init("#version 330 core");

  // --- Universe / sim seed ---
  // Keep this deterministic by default so procedural generation is stable.
  core::u64 seed = 1337;

  // --- Render assets ---
  render::Mesh sphere = render::Mesh::makeUvSphere(48, 24);
  render::Mesh cube   = render::Mesh::makeCube();

  render::Texture2D checker;
  checker.createChecker(256, 256, 16);

  std::string err;

  // --- Ship livery (procedural texture) ---
  // Stored in a tiny text file next to the savegame for quick iteration.
  const std::string liveryPath = ui::defaultLiveryPath();
  ui::LiveryConfig liveryCfg = ui::makeDefaultLivery();
  (void)ui::loadLiveryFromFile(liveryPath, liveryCfg);

  // Dedicated RNG so randomizing cosmetics doesn't perturb gameplay/system RNG.
  core::SplitMix64 liveryRng(core::hashCombine(seed, core::fnv1a64("livery")));

  render::Texture2D shipLiveryTex;
  bool shipLiveryDirty = true;
  double shipLiveryLastRegenTimeSec = -1e9;
  float shipLiveryRegenMinIntervalSec = 0.08f;

  auto liveryToDesc = [&](const ui::LiveryConfig& c) -> render::LiveryDesc {
    render::LiveryDesc d{};
    d.pattern = c.pattern;
    d.seed = c.seed;
    d.base[0] = c.base[0]; d.base[1] = c.base[1]; d.base[2] = c.base[2];
    d.accent1[0] = c.accent1[0]; d.accent1[1] = c.accent1[1]; d.accent1[2] = c.accent1[2];
    d.accent2[0] = c.accent2[0]; d.accent2[1] = c.accent2[1]; d.accent2[2] = c.accent2[2];
    d.scale = c.scale;
    d.angleDeg = c.angleDeg;
    d.detail = c.detail;
    d.wear = c.wear;
    d.contrast = c.contrast;
    d.decal = c.decal;
    return d;
  };

  auto rebuildShipLivery = [&]() {
    const render::LiveryDesc desc = liveryToDesc(liveryCfg);
    const int sizePx = std::clamp(liveryCfg.textureSizePx, 64, 2048);
    const auto img = render::generateLiveryTexture(desc, sizePx);
    if (!img.rgba.empty()) {
      // Repeat wrap: UVs on simple meshes (cube) look nicer with repeat than clamp.
      shipLiveryTex.createRGBA(img.w, img.h, img.rgba.data(),
                              /*generateMips=*/true,
                              /*nearestFilter=*/false,
                              /*clampToEdge=*/false);
    }
    shipLiveryDirty = false;
  };

  // Initial livery build.
  rebuildShipLivery();

  // Hangar preview render target (render-to-texture -> ImGui::Image).
  render::RenderTarget2D hangarTarget;
  if (!hangarTarget.init(640, 640, &err)) {
    core::log(core::LogLevel::Warn, err);
    err.clear();
  }
  int hangarPreviewSizePx = 640;
  bool hangarAnimate = true;
  float hangarSpinRadPerSec = 0.35f;
  float hangarYawDeg = 35.0f;
  float hangarPitchDeg = 18.0f;
  float hangarZoom = 3.25f;
  bool liveryAutoSaveOnExit = true;

  // UI icon sprites (procedural, cached as GL textures)
  render::SpriteCache spriteCache;
  spriteCache.setMaxEntries(2048);

  // HUD icon atlas: packs many procedural icons into a single texture so HUD overlays
  // (radar, tactical markers, target reticle) can draw lots of icons without binding
  // a different texture per marker.
  render::SpriteAtlas hudAtlas;
  hudAtlas.init(/*atlasSizePx=*/1024, /*cellSizePx=*/64, /*paddingPx=*/2, /*nearestFilter=*/true);
  hudAtlas.setMaxEntries(256);

  // Procedural surface textures for stars/planets (UV-sphere equirectangular albedo).
  // This is separate from the HUD icon sprites: these textures are intended for 3D meshes.
  render::SurfaceTextureCache surfaceTexCache;
  surfaceTexCache.setMaxEntries(128);

  render::MeshRenderer meshRenderer;
  if (!meshRenderer.init(&err)) {
    core::log(core::LogLevel::Error, err);
    return 1;
  }
  meshRenderer.setTexture(&checker);

  render::AtmosphereRenderer atmosphereRenderer;
  if (!atmosphereRenderer.init(&err)) {
    core::log(core::LogLevel::Error, err);
    return 1;
  }
  atmosphereRenderer.setMesh(&sphere);

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

  // Point-sprite textures (procedural): used for starfield/particles when "textured" mode is enabled.
  // These are small RGBA textures sampled by the PointRenderer shader via gl_PointCoord.
  render::Texture2D starPointSpriteTex = makeRadialSpriteTextureRGBA(
      64,
      /*coreExp*/ 18.0f,
      /*haloExp*/ 4.0f,
      /*haloStrength*/ 0.18f,
      /*spikeStrength*/ 0.35f,
      /*spikes*/ 4,
      seed ^ 0xA55A5AA5u);

  render::Texture2D particlePointSpriteTex = makeRadialSpriteTextureRGBA(
      64,
      /*coreExp*/ 10.0f,
      /*haloExp*/ 2.5f,
      /*haloStrength*/ 0.32f,
      /*spikeStrength*/ 0.0f,
      /*spikes*/ 0,
      seed ^ 0xC0DEC0DEu);

  // Nebula puff sprite texture (procedural alpha clouds; tinted per-vertex).
  render::Texture2D nebulaPointSpriteTex = makeCloudSpriteTextureRGBA(
      128,
      seed ^ 0xF00DF00DF00DF00DULL);

  // --- Post-processing (HDR + bloom + tonemap) ---
  render::PostFX postFx;
  if (!postFx.init(&err)) {
    core::log(core::LogLevel::Error, err);
    return 1;
  }
  render::PostFXSettings postFxSettings{};
  bool postFxAutoWarpFromSpeed = true;
  bool postFxAutoHyperspaceFromFsd = true;
  float postFxFsdWarpBoost = 0.065f;
  float postFxFsdHyperspaceBoost = 1.00f;
  bool hudJumpOverlay = true;

  // --- Universe / sim state ---

  // --- Visual effects (background stars + particles) ---
  bool vfxStarfieldEnabled = true;
  bool vfxStarfieldTextured = true;
  int vfxStarCount = 5200;
  double vfxStarRadiusU = 18000.0;

  // Nebula (large background cloud puffs)
  bool vfxNebulaEnabled = true;
  int vfxNebulaPuffCount = 1400;
  int vfxNebulaVariant = 0;
  int vfxNebulaVariantLast = vfxNebulaVariant;
  double vfxNebulaInnerRadiusU = 9000.0;
  double vfxNebulaOuterRadiusU = 22000.0;
  double vfxNebulaParallax = 0.25;
  float vfxNebulaIntensity = 1.4f;
  float vfxNebulaOpacity = 0.18f;
  float vfxNebulaSizeMinPx = 90.0f;
  float vfxNebulaSizeMaxPx = 320.0f;
  float vfxNebulaBandPower = 1.8f;
  float vfxNebulaBandPowerLast = vfxNebulaBandPower;
  float vfxNebulaTurbulence = 0.35f;
  float vfxNebulaTurbulenceSpeed = 0.10f;

  bool vfxParticlesEnabled = true;
  bool vfxParticlesTextured = true;
  bool vfxThrustersEnabled = true;
  bool vfxImpactsEnabled = true;
  bool vfxExplosionsEnabled = true;
  float vfxParticleIntensity = 1.0f; // global scaler

  // --- World visuals (stars/planets) ---
  bool worldUseProceduralSurfaces = true;
  bool worldUseSurfaceInUi = true; // show surface previews in tooltips
  int worldSurfaceTexWidth = 512;  // equirectangular width (height is width/2)
  bool worldStarUnlit = true;
  float worldStarIntensity = 3.5f; // HDR multiplier; affects bloom

  // Planet secondary layers (visual only).
  bool worldCloudsEnabled = true;
  float worldCloudOpacity = 0.55f;      // alpha multiplier (0..1)
  float worldCloudShellScale = 1.012f;  // sphere scale multiplier
  float worldCloudSpinDegPerSec = 2.0f; // visual rotation speed

  bool worldAtmospheresEnabled = true;
  float worldAtmoIntensity = 1.35f;     // additive color multiplier (HDR friendly)
  float worldAtmoPower = 5.25f;         // rim width
  float worldAtmoShellScale = 1.035f;   // sphere scale multiplier
  float worldAtmoSunLitBoost = 0.85f;   // how much day-side brightens the rim
  float worldAtmoForwardScatter = 0.35f; // when looking toward the star
  bool worldAtmoTintWithStar = true;

  render::Starfield starfield;
  starfield.setRadius(vfxStarRadiusU);
  starfield.regenerate(seed ^ 0xC5A0F1A9u, vfxStarCount);

  render::NebulaField nebula;
  {
    const core::u64 nSeed = core::hashCombine(seed ^ 0xBADC0FFEE0DDF00DULL, (core::u64)vfxNebulaVariant);
    nebula.regenerate(nSeed, vfxNebulaPuffCount, vfxNebulaBandPower);
  }

  render::ParticleSystem particles;
  particles.reseed(seed ^ 0xBADC0FFEu);
  particles.setMaxParticles(14000);

  std::vector<render::PointVertex> particleVerts;

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

  // --- Ship loadout / progression (shared: stellar/sim/ShipLoadout.h) ---
  using ShipHullClass = sim::ShipHullClass;
  using WeaponType = sim::WeaponType;
  using HullDef = sim::HullDef;
  using MkDef = sim::MkDef;
  using WeaponDef = sim::WeaponDef;

  // Local aliases to keep the rest of main.cpp diffs small.
  const auto& kHullDefs = sim::kHullDefs;
  const auto& kThrusters = sim::kThrusters;
  const auto& kShields = sim::kShields;
  const auto& kDistributors = sim::kDistributors;
  const auto& kWeaponDefs = sim::kWeaponDefs;

  auto weaponDef = [&](WeaponType t) -> const WeaponDef& {
    return sim::weaponDef(t);
  };

  ShipHullClass shipHullClass = ShipHullClass::Scout;
  int thrusterMk = 1;
  int shieldMk = 1;
  int distributorMk = 1;
  int smuggleHoldMk = 0; // 0..3 hidden compartments for contraband
  WeaponType weaponPrimary = WeaponType::BeamLaser;
  WeaponType weaponSecondary = WeaponType::Cannon;

  // Derived stats (recomputed from hull/modules)
  double playerHullMax = 100.0;
  double playerShieldMax = 100.0;
  double playerShieldRegenPerSimMin = 2.5; // points per simulated minute
  double playerHeatCoolRate = 10.0;        // heat points per real second
  double playerBaseLinAccelKmS2 = 0.08;
  double playerBaseAngAccelRadS2 = 1.2;

  // Player combat state
  double playerShield = 100.0;
  double playerHull = 100.0;
  double weaponPrimaryCooldown = 0.0;   // simulated seconds
  double weaponSecondaryCooldown = 0.0; // simulated seconds
  // Power distributor (pips + capacitors)
  sim::Pips distributorPips{2, 2, 2};
  sim::DistributorConfig distributorCfg = sim::distributorConfig(shipHullClass, distributorMk);
  sim::DistributorState distributorState = sim::makeFull(distributorCfg);

  // Used to scale heat input when boost is only partially available in a frame.
  double boostAppliedFrac = 0.0;


  // Baseline NPC combat tuning (used for random encounters / ambushes)
  const double npcShieldMax = 80.0;
  const double npcHullMax = 90.0;
  const double npcShieldRegenPerSec = 0.035; // points per simulated second (~2.1 per sim minute)

  // Helper: initialize a contact's loadout-derived stats + distributor.
  // NPCs are not persisted, so it's fine to reconstruct from these knobs at spawn time.
  auto configureContactLoadout = [&](Contact& c,
                                     ShipHullClass hull,
                                     int tMk,
                                     int sMk,
                                     int dMk,
                                     WeaponType weapon,
                                     double hullMul,
                                     double shieldMul,
                                     double regenMul,
                                     double accelMul,
                                     double aiSkill) {
    c.hullClass = hull;
    c.thrusterMk = std::clamp(tMk, 1, 3);
    c.shieldMk = std::clamp(sMk, 1, 3);
    c.distributorMk = std::clamp(dMk, 1, 3);
    c.weapon = weapon;
    c.aiSkill = std::clamp(aiSkill, 0.05, 0.98);

    const sim::ShipDerivedStats ds = sim::computeShipDerivedStats(
      c.hullClass, c.thrusterMk, c.shieldMk, c.distributorMk
    );

    // Use derived accelerations, but NPCs are intentionally a bit slower than the player by default.
    const double am = std::clamp(accelMul, 0.20, 2.0);
    c.ship.setMaxLinearAccelKmS2(ds.baseLinAccelKmS2 * am);
    c.ship.setMaxAngularAccelRadS2(ds.baseAngAccelRadS2 * am);

    c.hullMax = std::max(1.0, ds.hullMax * hullMul);
    c.hull = c.hullMax;
    c.shieldMax = std::max(0.0, ds.shieldMax * shieldMul);
    c.shield = c.shieldMax;
    c.shieldRegenPerSimMin = std::max(0.0, ds.shieldRegenPerSimMin * regenMul);

    c.pips = {2, 2, 2};
    sim::normalizePips(c.pips);
    c.distributorCfg = sim::distributorConfig(c.hullClass, c.distributorMk);
    c.distributorState = sim::makeFull(c.distributorCfg);
  };

  auto recalcPlayerStats = [&]() {
    const sim::ShipDerivedStats ds = sim::computeShipDerivedStats(
      shipHullClass, thrusterMk, shieldMk, distributorMk
    );

    playerHullMax = ds.hullMax;
    playerShieldMax = ds.shieldMax;
    playerShieldRegenPerSimMin = ds.shieldRegenPerSimMin;
    playerHeatCoolRate = ds.heatCoolRate;
    playerBaseLinAccelKmS2 = ds.baseLinAccelKmS2;
    playerBaseAngAccelRadS2 = ds.baseAngAccelRadS2;

    ship.setMaxLinearAccelKmS2(playerBaseLinAccelKmS2);
    ship.setMaxAngularAccelRadS2(playerBaseAngAccelRadS2);

    // Recompute distributor config and preserve capacitor fractions across loadout changes.
    {
      const double engFrac = (distributorCfg.capEng > 1e-9) ? (distributorState.eng / distributorCfg.capEng) : 1.0;
      const double wepFrac = (distributorCfg.capWep > 1e-9) ? (distributorState.wep / distributorCfg.capWep) : 1.0;
      const double sysFrac = (distributorCfg.capSys > 1e-9) ? (distributorState.sys / distributorCfg.capSys) : 1.0;

      distributorCfg = sim::distributorConfig(shipHullClass, distributorMk);
      sim::normalizePips(distributorPips);

      distributorState.eng = std::clamp(engFrac, 0.0, 1.0) * distributorCfg.capEng;
      distributorState.wep = std::clamp(wepFrac, 0.0, 1.0) * distributorCfg.capWep;
      distributorState.sys = std::clamp(sysFrac, 0.0, 1.0) * distributorCfg.capSys;
    }


    playerHull = std::clamp(playerHull, 0.0, playerHullMax);
    playerShield = std::clamp(playerShield, 0.0, playerShieldMax);
  };

  recalcPlayerStats();
  playerShield = playerShieldMax;
  playerHull = playerHullMax;

  // Time
  double timeDays = 0.0;
  double timeScale = 60.0; // simulated seconds per real second
  bool paused = false;

  // Economy
  double credits = 2500.0;
  double explorationDataCr = 0.0; // sellable scan data
  std::array<double, econ::kCommodityCount> cargo{};
  double cargoCapacityKg = 420.0;
  int passengerSeats = 2;
  int selectedStationIndex = 0;

  // Ship meta / progression
  double fuel = 45.0;
  double fuelMax = 45.0;
  double fsdRangeLy = 18.0;
  double fsdReadyDay = 0.0;

  // Heat (0..100). Kept intentionally simple for early gameplay feedback.
  double heat = 0.0;
  // Per-frame accumulator for instantaneous heat spikes (weapons, emergency drops, etc).
  // This is consumed by the ThermalSystem step.
  double heatImpulse = 0.0;

  // Missions + reputation
  core::u64 nextMissionId = 1;
  std::vector<sim::Mission> missions;

  // Industry / fabrication orders (station processing queues)
  core::u64 nextIndustryOrderId = 1;
  std::vector<sim::IndustryOrder> industryOrders;

  // Station warehouse / cargo storage (player-owned, per station).
  std::vector<sim::StationStorage> stationStorage;

  // Mission UI convenience: which mission is currently tracked in the HUD.
  core::u64 trackedMissionId = 0;
  bool missionTrackerAutoPlotNextLeg = true;
  bool objectiveHudEnabled = true; // overlay showing tracked mission objective in-flight
  std::unordered_map<core::u32, double> repByFaction;

  // Smuggling / contraband legality (per faction, deterministic bitmask of illegal commodities)
  std::unordered_map<core::u32, core::u32> illegalMaskByFaction;

  // Per-faction enforcement tuning (scan strictness, fine schedule, corruption).
  // Deterministic from (universe seed, faction id).
  std::unordered_map<core::u32, sim::LawProfile> lawProfileByFaction;

  // Law / crime (per-faction bounties)
  std::unordered_map<core::u32, double> bountyByFaction;
  // Bounty vouchers (earned by destroying criminals; redeem later at stations).
  std::unordered_map<core::u32, double> bountyVoucherByFaction;
  double policeAlertUntilDays = 0.0;
  double policeHeat = 0.0; // soft pursuit intensity (drives police escalation / spawn rate)

  // When police/stations scan you while wanted, you get a short "submit or fight" window.
  struct PoliceDemand {
    bool active{false};
    core::u32 factionId{0};
    double amountCr{0.0};
    double untilDays{0.0};
    std::string sourceName;
  };
  PoliceDemand policeDemand;


  // Pirates sometimes demand cargo instead of immediately opening fire.
  // Drop cargo pods (Ship/Status -> Cargo management -> Jettison) to satisfy them.
  struct PirateDemand {
    bool active{false};
    core::u64 groupId{0};
    core::u64 leaderId{0};
    double requiredValueCr{0.0};
    double deliveredValueCr{0.0};
    double startDays{0.0};
    double untilDays{0.0};
    std::string leaderName;
  };
  PirateDemand pirateDemand;

  // Escort mission runtime bookkeeping (not saved). Escort missions spawn a convoy contact that
  // the player must stay near while traveling between stations. We keep small transient state
  // here to avoid bloating the SaveGame format.
  struct EscortRuntime {
    core::u64 missionId{0};
    core::u64 convoyId{0};

    // If the player strays too far from the convoy for too long, the mission fails.
    double tooFarSec{0.0};

    // One-shot pirate ambush scheduling.
    bool ambushSpawned{false};
    double nextAmbushDays{0.0};
    core::u64 pirateGroupId{0};
  };
  std::vector<EscortRuntime> escortRuntime;

  // Mission-critical bounty targets (bounty scan/kill missions).
  // These are persisted so targets don't reset when leaving/re-entering a system or saving/loading.
  std::vector<sim::BountyTargetState> bountyTargetStates;

  // When contraband is detected by a police scan, you may get a short bribe/compliance window.
  // This is intentionally lightweight for the prototype loop: pay to keep cargo, comply to
  // accept confiscation + fine, or run and take a bounty.
  struct BribeOffer {
    bool active{false};
    core::u32 factionId{0};

    // Only used for police-scan bribes: lets us detect "run out of range" vs. comply.
    core::u64 scannerContactId{0};
    double scannerRangeKm{0.0};

    double amountCr{0.0};
    double fineCr{0.0};
    double illegalValueCr{0.0};

    double startDays{0.0};
    double untilDays{0.0};

    std::string sourceName;
    std::string detail;
    std::array<double, econ::kCommodityCount> scannedIllegal{};
  };
  BribeOffer bribeOffer;

  // Cargo scans (contraband)
  enum class CargoScanSourceKind { Station, Police };
  bool cargoScanActive = false;
  CargoScanSourceKind cargoScanSourceKind = CargoScanSourceKind::Station;
  sim::StationId cargoScanStationId = 0;
  core::u64 cargoScanContactId = 0;
  core::u32 cargoScanFactionId = 0;
  std::string cargoScanSourceName;
  double cargoScanProgressSec = 0.0;
  double cargoScanDurationSec = 0.0;
  double cargoScanRangeKm = 0.0;
  double cargoScanCooldownUntilDays = 0.0;

  // Exploration / discovery
  std::unordered_set<core::u64> scannedKeys; // scanned bodies/stations in the universe (player-local)

  // Cached mission board offers (regenerated when docking / day changes)
  std::vector<sim::Mission> missionOffers;
  sim::StationId missionOffersStationId = 0;
  int missionOffersDayStamp = -1;

  // Background NPC trade traffic stamps (per visited system).
  // Used to deterministically advance low-cost "ambient" station-to-station trade.
  std::unordered_map<sim::SystemId, int> trafficDayStampBySystem;

  // Cached trade helper suggestions (regenerated when docking / day changes)
  struct TradeIdea {
    sim::SystemId toSystem{0};
    sim::StationId toStation{0};
    std::string toSystemName;
    std::string toStationName;
    econ::CommodityId commodity{econ::CommodityId::Food};
    double buyPrice{0.0};
    double sellPrice{0.0};
    double unitsFrom{0.0};
    double unitsToSpace{0.0};
    double unitsPossible{0.0};
    double netProfitPerUnit{0.0};
    double netTripProfit{0.0};
    double profitPerKg{0.0};
    double distanceLy{0.0};
  };

  std::vector<TradeIdea> tradeIdeas;
  sim::StationId tradeFromStationId = 0;
  int tradeIdeasDayStamp = -1;
  double tradeSearchRadiusLy = 200.0;
  double tradeMinNetProfit = 0.0;
  int tradeCommodityFilter = -1; // -1 = any
  int tradeIdeasPerStation = 1;
  bool tradeUseFreeHold = true;
  bool tradeIncludeSameSystem = true;

  // Multi-commodity trade "mix" (cargo manifest) mode.
  bool tradeUseManifest = false;
  std::vector<sim::TradeManifestOpportunity> tradeMixIdeas;
  sim::StationId tradeMixFromStationId = 0;
  int tradeMixDayStamp = -1;
  double tradeMixStepKg = 1.0;
  int tradeMixLinesShown = 3;
  bool tradeMixSimulateImpact = true;


  // Docking state
  bool docked = false;
  sim::StationId dockedStationId = 0;
  std::unordered_map<sim::StationId, ClearanceState> clearances;

  // Contacts (pirates etc.)
  std::vector<Contact> contacts;
  core::SplitMix64 rng(seed ^ 0xC0FFEEu);
  sim::EncounterDirectorState encounterDirector = sim::makeEncounterDirector(seed, timeDays);

  // Salvage / mining / signal sources
  bool cargoScoopDeployed = true;
  double cargoFullToastCooldownUntilDays = 0.0;
  double incomingMissileToastCooldownUntilDays = 0.0;
  std::vector<FloatingCargo> floatingCargo;
  std::vector<AsteroidNode> asteroids;
  std::vector<SignalSource> signals;

  // Deterministic world-state persistence:
  // - resolvedSignalIds: one-shot system-entry signals (e.g., starter derelicts) that have already fired.
  // - asteroidRemainingById: remaining units for deterministically generated asteroids so mining can't be reset
  //   by leaving/re-entering a system.
  std::unordered_set<core::u64> resolvedSignalIds;
  std::unordered_map<core::u64, double> asteroidRemainingById;

  // Tag bit for deterministically-generated world objects.
  // Shared constant lives in sim/WorldIds.h so ids remain stable across apps/tools.
  constexpr core::u64 kDeterministicWorldIdBit = sim::kDeterministicWorldIdBit;

  core::u64 nextWorldObjectId = 1;
  double nextSignalSpawnDays = 0.01;
  double nextSupercruiseShadowSpawnDays = 0.01;

  // Beams (for laser visuals)
  struct Beam { math::Vec3d aU, bU; float r,g,b; double ttl; };
  std::vector<Beam> beams;

  // Projectiles (kinetic cannons / slugs)
  std::vector<sim::Projectile> projectiles;

  // Guided projectiles (homing missiles)
  std::vector<sim::Missile> missiles;

  // Save/load
  const std::string savePath = "savegame.txt";

  // Controls (keybindings)
  std::string controlsPath = game::defaultControlsPath();
  game::ControlsConfig controls = game::makeDefaultControls();
  const game::ControlsConfig controlsDefaults = game::makeDefaultControls();
  bool controlsDirty = false;
  bool controlsAutoSaveOnExit = true;
  game::ControlsWindowState controlsWindow{};
  game::ConsoleWindowState consoleWindow{};

  // UI state
  bool showGalaxy = true;
  bool showShip = true;
  bool showEconomy = true;
  bool showMissions = true;
  bool showContacts = true;
  bool showScanner = true;
  bool showTrade = true;
  bool showGuide = true;
  bool showSprites = false;
  bool showVfx = false;
  bool showPostFx = false;
  bool showWorldVisuals = false;
  bool showHangar = false;
  // controls window state lives in controlsWindow

  // HUD overlays
  bool showRadarHud = true;
  double radarRangeKm = 220000.0;
  int radarMaxBlips = 72;

  // Pirate demand HUD (threat/tribute): this HUD can appear contextually when pirates
  // extort you. Keep a master toggle so players can disable it.
  bool hudThreatOverlayEnabled = true;

  // HUD layout persistence: position + scale for in-flight overlays.
  const std::string hudLayoutPath = ui::defaultHudLayoutPath();
  ui::HudLayout hudLayout = ui::makeDefaultHudLayout();
  bool showHudLayoutWindow = false;
  bool hudLayoutEditMode = false;
  bool hudLayoutAutoSaveOnExit = true;

  // HUD settings persistence: master toggles + tuning (radar/combat/tactical).
  const std::string hudSettingsPath = ui::defaultHudSettingsPath();
  ui::HudSettings hudSettings = ui::makeDefaultHudSettings();
  ui::HudSettings hudSettingsSaved = hudSettings;
  bool showHudSettingsWindow = false;
  bool hudSettingsDirty = false;
  bool hudSettingsAutoSaveOnExit = true;

  // Draw an edge-of-screen arrow for the current target when it is off-screen.
  bool hudOffscreenTargetIndicator = true;

  // Combat HUD symbology (procedural reticle + weapon rings + lead indicator + flight path marker).
  bool hudCombatHud = true;
  bool hudUseProceduralReticle = true;
  bool hudShowWeaponRings = true;
  bool hudShowHeatRing = true;
  bool hudShowDistributorRings = true;
  float hudReticleSizePx = 44.0f;
  float hudReticleAlpha = 0.80f;

  bool hudShowLeadIndicator = true;
  bool hudLeadUseLastFiredWeapon = true;
  float hudLeadSizePx = 22.0f;
  double hudLeadMaxTimeSec = 18.0; // hide very long leads

  bool hudShowFlightPathMarker = true;
  bool hudFlightMarkerUseLocalFrame = true;
  bool hudFlightMarkerClampToEdge = true;
  float hudFlightMarkerSizePx = 22.0f;

  // Tracks last fired weapon for HUD lead indicator (primary LMB vs secondary RMB).
  bool hudLastFiredPrimary = true;

  // World-space icon overlay (screen projection + procedural sprites).
  bool showTacticalOverlay = true;
  bool tacticalShowLabels = true;
  double tacticalRangeKm = 450000.0;
  int tacticalMaxMarkers = 96;
  bool tacticalShowStations = true;
  bool tacticalShowPlanets = true;
  bool tacticalShowContacts = true;
  bool tacticalShowCargo = true;
  bool tacticalShowAsteroids = true;
  bool tacticalShowSignals = true;

  // Load HUD layout file (if present) and apply persisted widget enabled states.
  {
    ui::HudLayout loaded = ui::makeDefaultHudLayout();
    if (ui::loadFromFile(hudLayoutPath, loaded)) {
      hudLayout = loaded;
    }

    showRadarHud = hudLayout.widget(ui::HudWidgetId::Radar).enabled;
    objectiveHudEnabled = hudLayout.widget(ui::HudWidgetId::Objective).enabled;
    hudThreatOverlayEnabled = hudLayout.widget(ui::HudWidgetId::Threat).enabled;
    hudJumpOverlay = hudLayout.widget(ui::HudWidgetId::Jump).enabled;
  }

  auto syncHudSettingsFromRuntime = [&]() {
    hudSettings.autoSaveOnExit = hudSettingsAutoSaveOnExit;

    hudSettings.showRadarHud = showRadarHud;
    hudSettings.objectiveHudEnabled = objectiveHudEnabled;
    hudSettings.threatHudEnabled = hudThreatOverlayEnabled;
    hudSettings.jumpHudEnabled = hudJumpOverlay;

    hudSettings.radarRangeKm = radarRangeKm;
    hudSettings.radarMaxBlips = radarMaxBlips;

    hudSettings.offscreenTargetIndicator = hudOffscreenTargetIndicator;

    hudSettings.combatHudEnabled = hudCombatHud;
    hudSettings.useProceduralReticle = hudUseProceduralReticle;
    hudSettings.showWeaponRings = hudShowWeaponRings;
    hudSettings.showHeatRing = hudShowHeatRing;
    hudSettings.showDistributorRings = hudShowDistributorRings;
    hudSettings.reticleSizePx = hudReticleSizePx;
    hudSettings.reticleAlpha = hudReticleAlpha;

    hudSettings.showLeadIndicator = hudShowLeadIndicator;
    hudSettings.leadUseLastFiredWeapon = hudLeadUseLastFiredWeapon;
    hudSettings.leadSizePx = hudLeadSizePx;
    hudSettings.leadMaxTimeSec = hudLeadMaxTimeSec;

    hudSettings.showFlightPathMarker = hudShowFlightPathMarker;
    hudSettings.flightMarkerUseLocalFrame = hudFlightMarkerUseLocalFrame;
    hudSettings.flightMarkerClampToEdge = hudFlightMarkerClampToEdge;
    hudSettings.flightMarkerSizePx = hudFlightMarkerSizePx;

    hudSettings.tacticalOverlayEnabled = showTacticalOverlay;
    hudSettings.tacticalShowLabels = tacticalShowLabels;
    hudSettings.tacticalRangeKm = tacticalRangeKm;
    hudSettings.tacticalMaxMarkers = tacticalMaxMarkers;
    hudSettings.tacticalShowStations = tacticalShowStations;
    hudSettings.tacticalShowPlanets = tacticalShowPlanets;
    hudSettings.tacticalShowContacts = tacticalShowContacts;
    hudSettings.tacticalShowCargo = tacticalShowCargo;
    hudSettings.tacticalShowAsteroids = tacticalShowAsteroids;
    hudSettings.tacticalShowSignals = tacticalShowSignals;
  };

  auto applyHudSettingsToRuntime = [&](const ui::HudSettings& s) {
    hudSettingsAutoSaveOnExit = s.autoSaveOnExit;

    showRadarHud = s.showRadarHud;
    objectiveHudEnabled = s.objectiveHudEnabled;
    hudThreatOverlayEnabled = s.threatHudEnabled;
    hudJumpOverlay = s.jumpHudEnabled;

    radarRangeKm = s.radarRangeKm;
    radarMaxBlips = s.radarMaxBlips;

    hudOffscreenTargetIndicator = s.offscreenTargetIndicator;

    hudCombatHud = s.combatHudEnabled;
    hudUseProceduralReticle = s.useProceduralReticle;
    hudShowWeaponRings = s.showWeaponRings;
    hudShowHeatRing = s.showHeatRing;
    hudShowDistributorRings = s.showDistributorRings;
    hudReticleSizePx = s.reticleSizePx;
    hudReticleAlpha = s.reticleAlpha;

    hudShowLeadIndicator = s.showLeadIndicator;
    hudLeadUseLastFiredWeapon = s.leadUseLastFiredWeapon;
    hudLeadSizePx = s.leadSizePx;
    hudLeadMaxTimeSec = s.leadMaxTimeSec;

    hudShowFlightPathMarker = s.showFlightPathMarker;
    hudFlightMarkerUseLocalFrame = s.flightMarkerUseLocalFrame;
    hudFlightMarkerClampToEdge = s.flightMarkerClampToEdge;
    hudFlightMarkerSizePx = s.flightMarkerSizePx;

    showTacticalOverlay = s.tacticalOverlayEnabled;
    tacticalShowLabels = s.tacticalShowLabels;
    tacticalRangeKm = s.tacticalRangeKm;
    tacticalMaxMarkers = s.tacticalMaxMarkers;
    tacticalShowStations = s.tacticalShowStations;
    tacticalShowPlanets = s.tacticalShowPlanets;
    tacticalShowContacts = s.tacticalShowContacts;
    tacticalShowCargo = s.tacticalShowCargo;
    tacticalShowAsteroids = s.tacticalShowAsteroids;
    tacticalShowSignals = s.tacticalShowSignals;
  };

  auto hudSettingsEquivalent = [&](const ui::HudSettings& a, const ui::HudSettings& b) -> bool {
    auto deq = [](double x, double y, double eps = 1e-6) { return std::fabs(x - y) <= eps; };
    auto feq = [](float x, float y, float eps = 1e-4f) { return std::fabs(x - y) <= eps; };
    return a.autoSaveOnExit == b.autoSaveOnExit
        && a.showRadarHud == b.showRadarHud
        && a.objectiveHudEnabled == b.objectiveHudEnabled
        && a.threatHudEnabled == b.threatHudEnabled
        && a.jumpHudEnabled == b.jumpHudEnabled
        && deq(a.radarRangeKm, b.radarRangeKm)
        && a.radarMaxBlips == b.radarMaxBlips
        && a.offscreenTargetIndicator == b.offscreenTargetIndicator
        && a.combatHudEnabled == b.combatHudEnabled
        && a.useProceduralReticle == b.useProceduralReticle
        && a.showWeaponRings == b.showWeaponRings
        && a.showHeatRing == b.showHeatRing
        && a.showDistributorRings == b.showDistributorRings
        && feq(a.reticleSizePx, b.reticleSizePx)
        && feq(a.reticleAlpha, b.reticleAlpha)
        && a.showLeadIndicator == b.showLeadIndicator
        && a.leadUseLastFiredWeapon == b.leadUseLastFiredWeapon
        && feq(a.leadSizePx, b.leadSizePx)
        && deq(a.leadMaxTimeSec, b.leadMaxTimeSec)
        && a.showFlightPathMarker == b.showFlightPathMarker
        && a.flightMarkerUseLocalFrame == b.flightMarkerUseLocalFrame
        && a.flightMarkerClampToEdge == b.flightMarkerClampToEdge
        && feq(a.flightMarkerSizePx, b.flightMarkerSizePx)
        && a.tacticalOverlayEnabled == b.tacticalOverlayEnabled
        && a.tacticalShowLabels == b.tacticalShowLabels
        && deq(a.tacticalRangeKm, b.tacticalRangeKm)
        && a.tacticalMaxMarkers == b.tacticalMaxMarkers
        && a.tacticalShowStations == b.tacticalShowStations
        && a.tacticalShowPlanets == b.tacticalShowPlanets
        && a.tacticalShowContacts == b.tacticalShowContacts
        && a.tacticalShowCargo == b.tacticalShowCargo
        && a.tacticalShowAsteroids == b.tacticalShowAsteroids
        && a.tacticalShowSignals == b.tacticalShowSignals;
  };

  // Load HUD settings (if present). Missing file -> keep runtime defaults/layout-driven toggles.
  {
    ui::HudSettings loaded = ui::makeDefaultHudSettings();
    if (ui::loadHudSettingsFromFile(hudSettingsPath, loaded)) {
      hudSettings = loaded;
      applyHudSettingsToRuntime(hudSettings);
    } else {
      syncHudSettingsFromRuntime();
    }
    hudSettingsSaved = hudSettings;
    hudSettingsDirty = false;
  }

  // Load controls (if present). Missing file -> defaults.
  {
    game::loadFromFile(controlsPath, controls);
  }

  // Load bookmarks (if present). Missing file -> empty.
  {
    ui::Bookmarks loaded = ui::makeDefaultBookmarks();
    if (ui::loadBookmarksFromFile(bookmarksPath, loaded)) {
      bookmarks = std::move(loaded);
    }
  }

  // Optional mouse steering (relative mouse mode). Toggle with M.
  bool mouseSteer = false;
  float mouseSensitivity = 0.0025f; // torque intent per pixel
  bool mouseInvertY = false;

  // Flight assistance
  bool autopilot = false;
  sim::DockingComputer dockingComputer;
  bool dockingComputerDisengageOnManualInput = true;
  double dockingComputerManualDeadzone = 0.20;

  // Local reference frame ("space is local" feel near moving bodies)
  bool localFrameEnabled = true;
  math::Vec3d localFrameVelKmS{0,0,0};
  double localFrameBlendTauSec = 1.0; // seconds (real time)

  // "Blue zone" turn assist (max turn rate near ~mid-speed when flight assist is on)
  bool blueZoneTurnAssist = true;

  // Supercruise (Elite-style in-system travel)
  enum class SupercruiseState { Idle, Charging, Active, Cooldown };
  SupercruiseState supercruiseState = SupercruiseState::Idle;

  bool supercruiseAssist = true; // "Nav assist" keeps a safe drop profile
  double supercruiseMaxSpeedKmS = 18000.0;

  double supercruiseChargeRemainingSec = 0.0;   // real seconds
  double supercruiseCooldownRemainingSec = 0.0; // real seconds
  double supercruiseCooldownTotalSec = 0.0;     // real seconds (for UI)
  bool supercruiseDropRequested = false;
  bool supercruiseSafeDropReady = false;
  double supercruiseTtaSec = 0.0;
  double supercruiseDistKm = 0.0;
  double supercruiseClosingKmS = 0.0;

  // Interdiction (supercruise tether minigame)
  sim::InterdictionState interdiction;
  sim::InterdictionParams interdictionParams;
  sim::InterdictionTriggerParams interdictionTriggerParams;
  std::string interdictionPirateName;
  double interdictionPirateStrength = 1.0;
  bool interdictionSubmitRequested = false;
  double interdictionCooldownUntilDays = 0.0;

  // FSD / hyperspace (system-to-system)
  enum class FsdState { Idle, Charging, Jumping };
  FsdState fsdState = FsdState::Idle;
  sim::SystemId fsdTargetSystem{0};
  double fsdChargeRemainingSec = 0.0;
  double fsdTravelRemainingSec = 0.0;
  double fsdTravelTotalSec = 0.0;
  double fsdFuelCost = 0.0;
  double fsdJumpDistanceLy = 0.0;

  // Galaxy navigation / route plotting
  sim::SystemId galaxySelectedSystemId = 0;
  std::vector<sim::SystemId> navRoute;
  std::size_t navRouteHop = 0;
  bool navAutoRun = false;

  // Route planner settings
  enum class NavRouteMode { Hops = 0, Distance = 1, Fuel = 2 };
  NavRouteMode navRouteMode = NavRouteMode::Hops;
  // If true, route edges are constrained by current-fuel jump range (more "what can I do now?").
  // If false, edges are constrained by the ship's max range (useful for planning before refueling).
  bool navConstrainToCurrentFuelRange = true;

  // Planner diagnostics for the currently plotted route.
  sim::RoutePlanStats navRoutePlanStats{};
  bool navRoutePlanStatsValid = false;
  NavRouteMode navRoutePlannedMode = NavRouteMode::Hops;
  bool navRoutePlannedUsedCurrentFuelRange = true;
  double navRoutePlanMaxJumpLy = 0.0;
  double navRoutePlanCostPerJump = 0.0;
  double navRoutePlanCostPerLy = 0.0;

  // Optional helper: when arriving in a system, auto-target a specific station (e.g. from a mission/trade suggestion).
  sim::StationId pendingArrivalTargetStationId = 0;

  // Bookmarks: persistent navigation shortcuts (systems + stations).
  const std::string bookmarksPath = ui::defaultBookmarksPath();
  ui::Bookmarks bookmarks = ui::makeDefaultBookmarks();
  bool showBookmarksWindow = false;
  bool bookmarksDirty = false;
  bool bookmarksAutoSaveOnExit = true;

  // Mission Board: cached route preview for offers (computed when offers/settings change)
  struct MissionOfferRoutePreview {
    bool ok{false};
    int jumps{0};
    double distanceLy{0.0};
    double fuel{0.0};
  };
  bool missionOfferRoutePreviewEnabled = true;
  std::vector<MissionOfferRoutePreview> missionOfferRoutePreview;
  sim::StationId missionOfferRoutePreviewStationId = 0;
  int missionOfferRoutePreviewDayStamp = -1;
  NavRouteMode missionOfferRoutePreviewMode = NavRouteMode::Hops;
  bool missionOfferRoutePreviewConstrainToCurrentFuelRange = true;
  double missionOfferRoutePreviewMaxJumpLy = 0.0;
  sim::SystemId missionOfferRoutePreviewFromSystem = 0;

  // Scanner interaction (bounty scan + exploration scans)
  bool scanning = false;
  Target scanLockedTarget{};
  core::u64 scanLockedId = 0;
  std::string scanLabel;
  double scanProgressSec = 0.0;
  double scanDurationSec = 4.0;
  double scanRangeKm = 80000.0;

  // Target
  Target target{};

  std::vector<ToastMsg> toasts;

  // Persistent notification history (auto-filled via the `toast()` helper).
  std::deque<ToastHistoryEntry> toastHistory;
  bool showNotifications = false;
  int notificationsSelected = -1;
  char notificationsFilter[128]{};

  // UI: command palette (Ctrl+P by default)
  game::CommandPaletteState commandPalette;
  std::vector<game::PaletteItem> commandPaletteItems;

  // Wire the toast history sink now that timeDays + storage exist.
  setToastHistorySink(&toastHistory, &timeDays);

  // Console: hook core logging + built-in commands.
  consoleWindow.simTimeDaysPtr = &timeDays;
  game::consoleAddBuiltins(consoleWindow);
  game::consoleInstallCoreLogSink(consoleWindow);

  auto respawnNearStation = [&](const sim::StarSystem& sys, std::size_t stationIdx) {
    if (sys.stations.empty()) {
      ship.setPositionKm({0,0,-8000.0});
      ship.setVelocityKmS({0,0,0});
      ship.setAngularVelocityRadS({0,0,0});
      ship.setOrientation(math::Quatd::identity());
      return;
    }
    stationIdx = std::min(stationIdx, sys.stations.size() - 1);
    const auto& st = sys.stations[stationIdx];

    const math::Vec3d stPos = stationPosKm(st, timeDays);
    const math::Quatd stQ = stationOrient(st, stPos, timeDays);
    const math::Vec3d axis = stQ.rotate({0,0,1}); // outward
    const double startDistKm = st.radiusKm * 14.0;
    ship.setPositionKm(stPos + axis * startDistKm);
    ship.setVelocityKmS(stationVelKmS(st, timeDays));
    ship.setAngularVelocityRadS({0,0,0});
    // face toward the slot (into station)
    ship.setOrientation(quatFromTo({0,0,1}, -axis));
  };

  auto allocWorldId = [&]() -> core::u64 {
    // Keep IDs stable within a session; collisions with old IDs don't matter since these are transient.
    return nextWorldObjectId++;
  };

  auto randUnit = [&](double maxAbs = 1.0) -> math::Vec3d {
    math::Vec3d v{rng.range(-maxAbs, maxAbs), rng.range(-maxAbs, maxAbs), rng.range(-maxAbs, maxAbs)};
    const double lsq = v.lengthSq();
    if (lsq < 1e-9) return {1, 0, 0};
    return v / std::sqrt(lsq);
  };

  auto seedSystemSpaceObjects = [&](const sim::StarSystem& sys) {
    // Reset transient objects on system entry.
    floatingCargo.clear();
    asteroids.clear();
    signals.clear();
    cargoFullToastCooldownUntilDays = 0.0;
            incomingMissileToastCooldownUntilDays = 0.0;

    // Always seed at least one resource field near a "useful" station so mining is easy to find.
    if (sys.stations.empty()) return;
    std::size_t anchorIdx = 0;
    for (std::size_t i = 0; i < sys.stations.size(); ++i) {
      const auto t = sys.stations[i].type;
      if (t == econ::StationType::Mining || t == econ::StationType::Refinery) { anchorIdx = i; break; }
    }

    const auto& anchor = sys.stations[anchorIdx];
    const math::Vec3d anchorPos = stationPosKm(anchor, timeDays);

    // Deterministic in-system seeding:
    // - Stable IDs mean we can persist depletion/resolution across system transitions.
    // - Positions still evolve with time because anchors (stations) move on their orbits.
    const int dayStamp = (int)std::floor(timeDays);
    const core::u64 sysKey = core::hashCombine(seed, (core::u64)sys.stub.id);

    auto detId = [&](core::u64 typeCode, core::u64 salt) -> core::u64 {
      // Mix: universe seed + system id + type + salt, then tag as deterministic.
      return (core::hashCombine(core::hashCombine(sysKey, typeCode), salt) | kDeterministicWorldIdBit);
    };

    auto randUnitDet = [&](core::SplitMix64& srng) -> math::Vec3d {
      const double x = srng.range(-1.0, 1.0);
      const double y = srng.range(-1.0, 1.0);
      const double z = srng.range(-1.0, 1.0);
      const double len = std::sqrt(std::max(1e-9, x*x + y*y + z*z));
      return {x/len, y/len, z/len};
    };

    // Persistent resource fields (mining). Never expire; depletion is persisted by asteroid id.
    //
    // Generated in the core sim so the math can be unit tested.
    for (int fieldIdx = 0; fieldIdx < sim::kDefaultResourceFieldCount; ++fieldIdx) {
      const auto gen = sim::generateResourceField(seed,
                                                  sys.stub.id,
                                                  anchor.id,
                                                  anchorPos,
                                                  anchor.commsRangeKm,
                                                  timeDays,
                                                  fieldIdx);

      SignalSource s{};
      s.id = gen.field.signalId;
      s.type = SignalType::Resource;
      s.expireDay = 0.0;
      s.resolved = false;
      s.fieldSpawned = true;
      s.posKm = gen.field.posKm;
      s.hasResourcePlan = true;
      s.resource = gen.field;
      signals.push_back(s);

      for (const auto& ap : gen.asteroids) {
        AsteroidNode a{};
        a.id = ap.asteroidId;
        a.posKm = ap.posKm;
        a.radiusKm = ap.radiusKm;
        a.yield = ap.yield;

        const double baseUnits = std::max(0.0, ap.baseUnits);
        auto it = asteroidRemainingById.find(a.id);
        a.remainingUnits = (it != asteroidRemainingById.end()) ? std::clamp(it->second, 0.0, baseUnits)
                                                               : baseUnits;
        a.chunkAccumulator = 0.0;
        asteroids.push_back(a);
      }
    }

    // A daily derelict somewhere nearby for early salvage.
    // Deterministic id + persistence stops "farm by re-entering the system".
    {
      SignalSource s{};
      s.id = detId(/*typeCode=*/2, /*salt=*/(core::u64)dayStamp);
      s.type = SignalType::Derelict;
      s.expireDay = (double)(dayStamp + 1); // end of day
      s.resolved = (resolvedSignalIds.find(s.id) != resolvedSignalIds.end());

      core::SplitMix64 srng(core::hashCombine(s.id, 0xD311E1C7ull));
      const math::Vec3d dir = randUnitDet(srng);
      const double distKm = anchor.commsRangeKm * 1.6 + 190000.0;
      s.posKm = anchorPos + dir * distKm;
      signals.push_back(s);
    }

    // Mission-specific salvage sites: spawn a stable derelict signal for each active Salvage job
    // targeting this system. We only spawn these while the mission site has not been visited yet
    // (m.scanned==false) so players can't farm infinite "mission salvage" by re-entering the system.
    for (const auto& m : missions) {
      if (m.completed || m.failed) continue;
      if (m.type != sim::MissionType::Salvage) continue;
      if (m.toSystem != sys.stub.id) continue;
      if (m.scanned) continue;
      if (m.targetNpcId == 0) continue;

      // Find the mission's station (used as the anchor). Fallback to the system's anchor station.
      const sim::Station* baseSt = nullptr;
      for (const auto& st : sys.stations) {
        if (st.id == m.toStation) { baseSt = &st; break; }
      }
      const math::Vec3d basePos = baseSt ? stationPosKm(*baseSt, timeDays) : anchorPos;
      const double comms = baseSt ? baseSt->commsRangeKm : anchor.commsRangeKm;

      // Deterministic placement: combine the system id + mission "signal id" into a local seed.
      core::SplitMix64 srng(core::hashCombine((core::u64)sys.stub.id, (core::u64)m.targetNpcId));
      auto randUnitDet = [&]() -> math::Vec3d {
        const double x = srng.range(-1.0, 1.0);
        const double y = srng.range(-1.0, 1.0);
        const double z = srng.range(-1.0, 1.0);
        const double len = std::sqrt(std::max(1e-9, x*x + y*y + z*z));
        return {x/len, y/len, z/len};
      };

      const double dist = comms * 1.8 + srng.range(150000.0, 260000.0);

      SignalSource s{};
      s.id = m.targetNpcId; // stable mission signal id
      s.type = SignalType::Derelict;
      s.posKm = basePos + randUnitDet() * dist;

      // Keep it alive at least until the mission deadline (clamped).
      const double ttl = (m.deadlineDay > timeDays) ? std::clamp(m.deadlineDay - timeDays, 0.25, 4.0) : 1.0;
      s.expireDay = timeDays + ttl;
      s.resolved = false;
      signals.push_back(s);
    }
  };

  // Spawn near first station for immediate gameplay.
  respawnNearStation(*currentSystem, 0);
  seedSystemSpaceObjects(*currentSystem);

  galaxySelectedSystemId = currentSystem->stub.id;

  auto findFaction = [&](core::u32 factionId) -> const sim::Faction* {
    for (const auto& f : universe.factions()) {
      if (f.id == factionId) return &f;
    }
    return nullptr;
  };

  auto factionName = [&](core::u32 factionId) -> std::string {
    if (factionId == 0) return "Independent";
    if (auto f = findFaction(factionId)) return f->name;
    return "Faction " + std::to_string(factionId);
  };

  auto getRep = [&](core::u32 factionId) -> double {
    auto it = repByFaction.find(factionId);
    return it == repByFaction.end() ? 0.0 : it->second;
  };

  auto addRep = [&](core::u32 factionId, double delta) {
    if (factionId == 0) return;
    repByFaction[factionId] = clampRep(getRep(factionId) + delta);
  };

auto getBounty = [&](core::u32 factionId) -> double {
  auto it = bountyByFaction.find(factionId);
  return it == bountyByFaction.end() ? 0.0 : std::max(0.0, it->second);
};

auto addBounty = [&](core::u32 factionId, double deltaCr) {
  if (factionId == 0) return;
  bountyByFaction[factionId] = std::max(0.0, getBounty(factionId) + deltaCr);
};

auto clearBounty = [&](core::u32 factionId) {
  if (factionId == 0) return;
  bountyByFaction[factionId] = 0.0;
};

  auto getVoucher = [&](core::u32 factionId) -> double {
    auto it = bountyVoucherByFaction.find(factionId);
    return it == bountyVoucherByFaction.end() ? 0.0 : std::max(0.0, it->second);
  };

  auto addVoucher = [&](core::u32 factionId, double deltaCr) {
    if (factionId == 0) return;
    bountyVoucherByFaction[factionId] = std::max(0.0, getVoucher(factionId) + deltaCr);
  };

  auto clearVoucher = [&](core::u32 factionId) {
    if (factionId == 0) return;
    bountyVoucherByFaction[factionId] = 0.0;
  };

auto commitCrime = [&](core::u32 factionId, double bountyAddCr, double repPenalty, const std::string& reason, bool showToast = true) {
  if (factionId == 0) return;
  addBounty(factionId, bountyAddCr);
  addRep(factionId, repPenalty);
  // Escalation: bigger crimes ramp up response intensity for a while.
  policeHeat = std::clamp(policeHeat + 0.75 + std::min(2.5, bountyAddCr / 1500.0), 0.0, 6.0);

  const double alertSec = 120.0 + 25.0 * policeHeat;
  policeAlertUntilDays = std::max(policeAlertUntilDays, timeDays + (alertSec / 86400.0));

  // Higher heat = faster reinforcements.
  const double respSec = std::clamp(7.0 - 0.65 * policeHeat, 2.0, 7.0);
  encounterDirector.nextPoliceSpawnDays = std::min(encounterDirector.nextPoliceSpawnDays, timeDays + (respSec / 86400.0));
  if (showToast) {
    toast(toasts, "Crime (" + factionName(factionId) + "): " + reason, 3.0);
  }
};


// --- Smuggling / contraband -------------------------------------------------
auto illegalMaskForFaction = [&](core::u32 factionId) -> core::u32 {
  if (factionId == 0) return 0u;

  auto it = illegalMaskByFaction.find(factionId);
  if (it != illegalMaskByFaction.end()) return it->second;

  const core::u32 mask = sim::illegalCommodityMask(universe.seed(), factionId);
  illegalMaskByFaction[factionId] = mask;
  return mask;
};

auto isIllegalCommodity = [&](core::u32 factionId, econ::CommodityId cid) -> bool {
  if (factionId == 0) return false;
  return (illegalMaskForFaction(factionId) & sim::commodityBit(cid)) != 0u;
};

auto illegalListString = [&](core::u32 factionId) -> std::string {
  return sim::illegalCommodityListString(universe.seed(), factionId);
};

auto hasIllegalCargo = [&](core::u32 factionId) -> bool {
  return sim::hasIllegalCargo(universe.seed(), factionId, cargo);
};

auto lawForFaction = [&](core::u32 factionId) -> sim::LawProfile {
  if (factionId == 0) return sim::LawProfile{};

  auto it = lawProfileByFaction.find(factionId);
  if (it != lawProfileByFaction.end()) return it->second;

  const sim::LawProfile prof = sim::lawProfile(universe.seed(), factionId);
  lawProfileByFaction.emplace(factionId, prof);
  return prof;
};

// Apply confiscation + fine + reputation/bounty/alert effects when contraband is enforced.
// `scannedIllegal` is what security saw during the scan (used for messaging + smuggle mission attribution).
auto enforceContraband = [&](core::u32 jurisdiction,
                             const std::string& sourceName,
                             double illegalValueCr,
                             const std::string& detail,
                             const std::array<double, econ::kCommodityCount>& scannedIllegal) {
  if (jurisdiction == 0) return;

  const sim::LawProfile law = lawForFaction(jurisdiction);

  const auto res = sim::enforceContraband(law, credits, cargo, scannedIllegal, illegalValueCr);
  credits = res.creditsAfter;
  cargo = res.cargoAfter;

  // Rep / bounty
  addRep(jurisdiction, res.repPenalty);
  if (res.bountyAddedCr > 1e-6) {
    addBounty(jurisdiction, res.bountyAddedCr);
  }

  // Heightened police attention for a short window.
  policeAlertUntilDays = std::max(policeAlertUntilDays, timeDays + (res.policeAlertSeconds / 86400.0));
  encounterDirector.nextPoliceSpawnDays = std::min(encounterDirector.nextPoliceSpawnDays, timeDays + (res.nextPoliceSpawnDelaySeconds / 86400.0));

  // Raise local security alert (affects patrol spawn rate / response).
  policeHeat = std::clamp(policeHeat + res.policeHeatDelta, 0.0, 6.0);

  std::string msg = "Contraband detected! Confiscated: "
                    + (detail.empty() ? std::string("illegal cargo") : detail)
                    + ". Fine: " + std::to_string((int)std::round(res.fineCr)) + " cr"
                    + (res.unpaidCr > 1e-6 ? " (unpaid -> bounty)" : "")
                    + ". Rep " + std::to_string((int)std::round(res.repPenalty));

  if (res.unpaidCr > 1e-6 && !docked) {
    // Give the player a short "submit or fight" window.
    policeDemand.active = true;
    policeDemand.factionId = jurisdiction;
    policeDemand.untilDays = timeDays + (18.0 / 86400.0);
    policeDemand.amountCr = getBounty(jurisdiction);
    policeDemand.sourceName = sourceName;
    msg += ". Press I to submit/pay.";
  }

  toast(toasts, msg, 4.0);

  // If the player loses contraband that was tied to an active smuggling job,
  // fail those jobs immediately (the contact's package is gone).
  int failedSmuggle = 0;
  for (auto& m : missions) {
    if (m.completed || m.failed) continue;
    if (m.type != sim::MissionType::Smuggle) continue;
    const std::size_t mi = (std::size_t)m.commodity;
    if (mi >= econ::kCommodityCount) continue;
    if (scannedIllegal[mi] <= 1e-6) continue;
    if (cargo[mi] + 1e-6 < m.units) {
      m.failed = true;
      addRep(m.factionId, -4.0);
      ++failedSmuggle;
    }
  }
  if (failedSmuggle > 0) {
    toast(toasts, "Smuggle mission failed: contraband confiscated.", 3.2);
  }
};

// Scan/discovery keys (player-local)
const auto scanKeyStar = [&](sim::SystemId sysId) -> core::u64 {
  return core::hashCombine((core::u64)sysId, 0x53544152ULL); // 'STAR'
};
const auto scanKeyPlanet = [&](sim::SystemId sysId, std::size_t planetIndex) -> core::u64 {
  return core::hashCombine((core::u64)sysId, core::hashCombine(0x504C414EULL, (core::u64)planetIndex)); // 'PLAN'
};
const auto scanKeyStation = [&](sim::StationId stId) -> core::u64 {
  return core::hashCombine((core::u64)stId, 0x53544154ULL); // 'STAT'
};
	const auto scanKeySignal = [&](core::u64 signalId) -> core::u64 {
	  return core::hashCombine((core::u64)signalId, 0x5349474EULL); // 'SIGN'
	};
	const auto scanKeyAsteroid = [&](core::u64 asteroidId) -> core::u64 {
	  return core::hashCombine((core::u64)asteroidId, 0x41535452ULL); // 'ASTR'
	};
const auto scanKeySystemComplete = [&](sim::SystemId sysId) -> core::u64 {
  return core::hashCombine((core::u64)sysId, 0x434F4D50ULL); // 'COMP'
};

  auto effectiveFeeRate = [&](const sim::Station& st) -> double {
    return applyRepToFee(st.feeRate, getRep(st.factionId));
  };

  // FSD / jump parameters
  const double kFsdFuelBase = 2.0;
  const double kFsdFuelPerLy = 0.5;
  const double kFsdChargeSec = 4.0;
  const double kFsdCooldownSec = 25.0;

  // Supercruise parameters (in-system travel mode)
  const double kSupercruiseChargeSec = 3.0;
  const double kSupercruiseCooldownSec = 6.0;
  const double kSupercruiseEmergencyCooldownSec = 14.0;
  const double kSupercruiseSafeTtaSec = 7.0; // the classic "7-second rule"

  auto fsdBaseRangeLy = [&]() -> double {
    const double cap = std::max(1.0, cargoCapacityKg);
    const double load = std::clamp(cargoMassKg(cargo) / cap, 0.0, 1.0);
    // Cargo load reduces effective range a bit (keeps hauling interesting).
    return std::max(0.0, fsdRangeLy * (1.0 - 0.25 * load));
  };

  auto fsdFuelLimitedRangeLy = [&]() -> double {
    if (fuel <= kFsdFuelBase) return 0.0;
    return std::max(0.0, (fuel - kFsdFuelBase) / kFsdFuelPerLy);
  };

  auto fsdCurrentRangeLy = [&]() -> double {
    return std::min(fsdBaseRangeLy(), fsdFuelLimitedRangeLy());
  };

  auto fsdFuelCostFor = [&](double distanceLy) -> double {
    return kFsdFuelBase + distanceLy * kFsdFuelPerLy;
  };

  // ---- Navigation helpers ----
  auto tryTargetStationById = [&](sim::StationId stationId) -> bool {
    if (!currentSystem) return false;
    for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
      if (currentSystem->stations[i].id == stationId) {
        target.kind = TargetKind::Station;
        target.index = i;
        selectedStationIndex = (int)i;
        return true;
      }
    }
    return false;
  };

  auto tryTargetSignalById = [&](core::u64 signalId) -> bool {
    for (std::size_t i = 0; i < signals.size(); ++i) {
      if (signals[i].id == signalId) {
        target.kind = TargetKind::Signal;
        target.index = i;
        return true;
      }
    }
    return false;
  };

  auto plotRouteToSystem = [&](sim::SystemId destSystemId, bool showToast = true) -> bool {
    if (!currentSystem || destSystemId == 0) return false;

    if (destSystemId == currentSystem->stub.id) {
      navRoute.clear();
      navRouteHop = 0;
      navAutoRun = false;
      galaxySelectedSystemId = destSystemId;
      return true;
    }

    const auto& destStub = universe.getSystem(destSystemId).stub;
    const double distLy = (destStub.posLy - currentSystem->stub.posLy).length();

    const double jrMaxLy = navConstrainToCurrentFuelRange ? fsdCurrentRangeLy() : fsdBaseRangeLy();
    if (jrMaxLy <= 0.0) {
      if (showToast) toast(toasts, "No jump range available (refuel or reduce cargo).", 2.6);
      return false;
    }
    double radius = std::clamp(distLy * 1.20 + 60.0, 180.0, 1400.0);

    for (int attempt = 0; attempt < 6; ++attempt) {
      auto nearby = universe.queryNearby(currentSystem->stub.posLy, radius);
      sim::RoutePlanStats stats{};

      double costPerJump = 1.0;
      double costPerLy = 0.0;
      std::vector<sim::SystemId> route;

      if (navRouteMode == NavRouteMode::Hops) {
        costPerJump = 1.0;
        costPerLy = 0.0;
        route = sim::plotRouteAStarHops(nearby, currentSystem->stub.id, destSystemId, jrMaxLy, &stats);
      } else if (navRouteMode == NavRouteMode::Distance) {
        costPerJump = 0.0;
        costPerLy = 1.0;
        route = sim::plotRouteAStarCost(nearby, currentSystem->stub.id, destSystemId, jrMaxLy, costPerJump, costPerLy, &stats);
      } else {
        costPerJump = kFsdFuelBase;
        costPerLy = kFsdFuelPerLy;
        route = sim::plotRouteAStarCost(nearby, currentSystem->stub.id, destSystemId, jrMaxLy, costPerJump, costPerLy, &stats);
      }
      if (!route.empty()) {
        navRoute = std::move(route);
        navRouteHop = 0;
        navAutoRun = false;
        galaxySelectedSystemId = destSystemId;

        navRoutePlanStats = stats;
        navRoutePlanStatsValid = true;
        navRoutePlannedMode = navRouteMode;
        navRoutePlannedUsedCurrentFuelRange = navConstrainToCurrentFuelRange;
        navRoutePlanMaxJumpLy = jrMaxLy;
        navRoutePlanCostPerJump = costPerJump;
        navRoutePlanCostPerLy = costPerLy;

        if (showToast) {
          const double totalDist = sim::routeDistanceLy(nearby, navRoute);
          const double totalFuel = sim::routeCost(nearby, navRoute, kFsdFuelBase, kFsdFuelPerLy);
          toast(toasts,
                "Route plotted (" + std::to_string((int)navRoute.size() - 1) + " jumps, "
                + std::to_string((int)std::round(totalDist)) + " ly, est fuel "
                + std::to_string((int)std::round(totalFuel)) + ").",
                2.4);
        }
        return true;
      }
      radius *= 1.35;
    }

    navRoutePlanStatsValid = false;
    if (showToast) toast(toasts, "No route found (try refueling or upgrading your FSD).", 2.6);
    return false;
  };

  // Console: game commands (stateful helpers).
  // These are registered once at startup and can be invoked from the Console window.
  game::consoleAddCommand(consoleWindow, "pause", "Pause simulation time.",
    [&](game::ConsoleWindowState& c, const std::vector<std::string_view>&) {
      paused = true;
      game::consolePrint(c, core::LogLevel::Info, "Paused.");
    });

  game::consoleAddCommand(consoleWindow, "unpause", "Resume simulation time.",
    [&](game::ConsoleWindowState& c, const std::vector<std::string_view>&) {
      paused = false;
      game::consolePrint(c, core::LogLevel::Info, "Unpaused.");
    });

  game::consoleAddCommand(consoleWindow, "timescale", "Set time scale. Usage: timescale <multiplier>",
    [&](game::ConsoleWindowState& c, const std::vector<std::string_view>& args) {
      if (args.empty()) {
        game::consolePrint(c, core::LogLevel::Info, "Usage: timescale <multiplier>");
        return;
      }
      try {
        const double v = std::stod(std::string(args[0]));
        timeScale = std::clamp(v, 0.0, 200.0);
        paused = (timeScale <= 0.0);
        game::consolePrint(c, core::LogLevel::Info, "Time scale set.");
      } catch (...) {
        game::consolePrint(c, core::LogLevel::Warn, "Invalid number.");
      }
    });

  game::consoleAddCommand(consoleWindow, "credits", "Set credits. Usage: credits <amount>",
    [&](game::ConsoleWindowState& c, const std::vector<std::string_view>& args) {
      if (args.empty()) {
        game::consolePrint(c, core::LogLevel::Info, "Usage: credits <amount>");
        return;
      }
      try {
        const double v = std::stod(std::string(args[0]));
        credits = std::max(0.0, v);
        game::consolePrint(c, core::LogLevel::Info, "Credits updated.");
      } catch (...) {
        game::consolePrint(c, core::LogLevel::Warn, "Invalid number.");
      }
    });

  game::consoleAddCommand(consoleWindow, "window", "Toggle UI windows. Usage: window <name> [on|off|toggle]",
    [&](game::ConsoleWindowState& c, const std::vector<std::string_view>& args) {
      if (args.empty()) {
        game::consolePrint(c, core::LogLevel::Info, "Usage: window <name> [on|off|toggle]");
        game::consolePrint(c, core::LogLevel::Info, "Names: galaxy ship market contacts missions scanner trade guide hangar visuals notifications bookmarks console");
        return;
      }
      auto lower = [&](std::string_view s) {
        std::string out;
        out.reserve(s.size());
        for (char ch : s) out.push_back((char)std::tolower((unsigned char)ch));
        return out;
      };
      const std::string name = lower(args[0]);
      const std::string mode = args.size() >= 2 ? lower(args[1]) : std::string("toggle");
      auto apply = [&](bool& flag) {
        if (mode == "on" || mode == "1" || mode == "true") flag = true;
        else if (mode == "off" || mode == "0" || mode == "false") flag = false;
        else flag = !flag;
      };

      if (name == "galaxy") apply(showGalaxy);
      else if (name == "ship") apply(showShip);
      else if (name == "market") apply(showEconomy);
      else if (name == "contacts") apply(showContacts);
      else if (name == "missions") apply(showMissions);
      else if (name == "scanner") apply(showScanner);
      else if (name == "trade") apply(showTrade);
      else if (name == "guide") apply(showGuide);
      else if (name == "hangar") apply(showHangar);
      else if (name == "visuals" || name == "world" || name == "worldvisuals") apply(showWorldVisuals);
      else if (name == "notifications") apply(showNotifications);
      else if (name == "bookmarks") apply(showBookmarksWindow);
      else if (name == "console") { apply(consoleWindow.open); if (consoleWindow.open) consoleWindow.focusInput = true; }
      else {
        game::consolePrint(c, core::LogLevel::Warn, "Unknown window name.");
        return;
      }
      game::consolePrint(c, core::LogLevel::Info, "OK.");
    });

  game::consoleAddCommand(consoleWindow, "route", "Plot a route to a system id. Usage: route <system_id>",
    [&](game::ConsoleWindowState& c, const std::vector<std::string_view>& args) {
      if (args.empty()) {
        game::consolePrint(c, core::LogLevel::Info, "Usage: route <system_id>");
        return;
      }
      try {
        const auto id = (sim::SystemId)std::stoull(std::string(args[0]));
        if (plotRouteToSystem(id, /*showToast=*/true)) {
          game::consolePrint(c, core::LogLevel::Info, "Route plotted.");
        } else {
          game::consolePrint(c, core::LogLevel::Warn, "Route failed.");
        }
      } catch (...) {
        game::consolePrint(c, core::LogLevel::Warn, "Invalid system id.");
      }
    });


  auto isMassLocked = [&]() -> bool {
    if (!currentSystem) return false;
    const math::Vec3d shipPos = ship.positionKm();

    // Stations (treat as heavy bodies / traffic control).
    for (const auto& st : currentSystem->stations) {
      const math::Vec3d stPos = stationPosKm(st, timeDays);
      const double distKm = (shipPos - stPos).length();
      const double lockKm = st.radiusKm * 12.0;
      if (distKm < lockKm) return true;
    }

    // Planets (simple gravity-well proxy).
    for (const auto& p : currentSystem->planets) {
      const math::Vec3d pPos = sim::orbitPosition3DAU(p.orbit, timeDays) * kAU_KM;
      const double rKm = p.radiusEarth * kEARTH_RADIUS_KM;
      const double distKm = (shipPos - pPos).length();
      const double lockKm = std::max(5000.0, rKm * 20.0);
      if (distKm < lockKm) return true;
    }

    return false;
  };

  auto startFsdJumpTo = [&](sim::SystemId destId) {
    if (destId == 0 || destId == currentSystem->stub.id) {
      toast(toasts, "No destination selected.", 1.8);
      return;
    }
    if (docked) {
      toast(toasts, "Can't jump while docked.", 2.0);
      return;
    }
    if (supercruiseState != SupercruiseState::Idle) {
      toast(toasts, "Disengage supercruise before jumping.", 2.0);
      return;
    }
    if (fsdState != FsdState::Idle) {
      toast(toasts, "FSD already busy.", 2.0);
      return;
    }
    if (timeDays < fsdReadyDay) {
      toast(toasts, "FSD cooling down...", 1.8);
      return;
    }
    if (isMassLocked()) {
      toast(toasts, "Mass-locked: move away from bodies/stations.", 2.2);
      return;
    }

    const sim::StarSystem& destSys = universe.getSystem(destId);
    const double distLy = (destSys.stub.posLy - currentSystem->stub.posLy).length();

    const double rangeLy = fsdBaseRangeLy();
    if (distLy > rangeLy + 1e-9) {
      toast(toasts, "Out of jump range (plot a multi-jump route).", 2.5);
      return;
    }

    const double fuelCost = fsdFuelCostFor(distLy);
    if (fuel < fuelCost) {
      toast(toasts, "Not enough fuel for this jump.", 2.5);
      return;
    }

    // Consume fuel on charge complete (so you can still cancel cleanly).
    fsdTargetSystem = destId;
    fsdJumpDistanceLy = distLy;
    fsdFuelCost = fuelCost;
    fsdChargeRemainingSec = kFsdChargeSec;
    fsdTravelRemainingSec = 0.0;
    fsdTravelTotalSec = 0.0;
    fsdState = FsdState::Charging;

    autopilot = false;
    scanning = false;
    scanProgressSec = 0.0;
    beams.clear();

    toast(toasts, "FSD charging...", 2.0);
  };

  // --- Gameplay helpers ---
  auto cargoValueEstimateCr = [&]() -> double {
    double v = 0.0;
    for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
      const auto cid = (econ::CommodityId)i;
      const double units = cargo[i];
      if (units <= 0.0) continue;
      v += units * econ::commodityDef(cid).basePrice;
    }
    return v;
  };


  auto resolvePirateDemand = [&](bool satisfied, const std::string& msg) {
    if (!pirateDemand.active) return;

    const core::u64 gid = pirateDemand.groupId;

    if (satisfied) {
      // Simulate pirates actually taking the tribute: remove player-jettisoned pods nearby.
      double wantCr = std::max(0.0, pirateDemand.requiredValueCr);
      if (wantCr > 0.0) {
        std::vector<std::pair<double, std::size_t>> candidates;
        candidates.reserve(floatingCargo.size());
        for (std::size_t i = 0; i < floatingCargo.size(); ++i) {
          const auto& pod = floatingCargo[i];
          if (pod.units <= 0.0) continue;
          if (!pod.fromPlayer) continue;
          const double d = (pod.posKm - ship.positionKm()).length();
          if (d > 200000.0) continue;
          candidates.push_back({d, i});
        }

        std::sort(candidates.begin(), candidates.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        for (const auto& it : candidates) {
          if (wantCr <= 0.0) break;
          auto& pod = floatingCargo[it.second];
          const double v = pod.units * econ::commodityDef(pod.commodity).basePrice;
          wantCr -= v;
          pod.units = 0.0;
        }
      }
    }

    if (gid != 0) {
      for (auto& c : contacts) {
        if (!c.alive) continue;
        if (c.role != ContactRole::Pirate) continue;
        if (c.groupId != gid) continue;

        if (satisfied) {
          c.hostileToPlayer = false;
          c.fleeUntilDays = timeDays + (rng.range(70.0, 150.0) / 86400.0);
        } else {
          c.hostileToPlayer = true;
          c.fireCooldown = std::min(c.fireCooldown, 0.10);
        }
      }
    }

    if (!msg.empty()) toast(toasts, msg, 3.6);
    pirateDemand = PirateDemand{};
  };

  auto spawnCargoPod = [&](econ::CommodityId cid, double units, const math::Vec3d& posKm,
                          const math::Vec3d& inheritVelKmS, double scatterMaxKmS,
                          bool fromPlayer = false, core::u64 missionId = 0) {
    if (units <= 0.0) return;
    FloatingCargo pod{};
    pod.id = allocWorldId();
    pod.commodity = cid;
    pod.units = units;
    pod.posKm = posKm;
    pod.velKmS = inheritVelKmS + randUnit() * rng.range(0.0, std::max(0.0, scatterMaxKmS));
    pod.expireDay = timeDays + (rng.range(8.0 * 60.0, 18.0 * 60.0) / 86400.0); // 8–18 min
    pod.fromPlayer = fromPlayer;
    pod.missionId = missionId;
    floatingCargo.push_back(pod);
  };

  auto spawnCargoBurst = [&](econ::CommodityId cid, double totalUnits, const math::Vec3d& posKm,
                            const math::Vec3d& inheritVelKmS, int pods,
                            bool fromPlayer = false, core::u64 missionId = 0) {
    pods = std::max(1, std::min(pods, 6));
    totalUnits = std::max(0.0, totalUnits);
    if (totalUnits <= 0.0) return;
    const double per = totalUnits / (double)pods;
    for (int i = 0; i < pods; ++i) {
      // Small random split so the last pod isn't always tiny.
      const double u = (i == pods - 1) ? (totalUnits - per * (pods - 1)) : (per * rng.range(0.75, 1.25));
      spawnCargoPod(cid, std::max(0.0, u), posKm, inheritVelKmS, 1.2, fromPlayer, missionId);
    }
  };

  bool running = true;
  auto last = std::chrono::high_resolution_clock::now();

  double timeRealSec = 0.0; // for purely visual effects

  SDL_SetRelativeMouseMode(SDL_FALSE);

  while (running) {
    // Timing
    auto now = std::chrono::high_resolution_clock::now();
    const double dtReal = std::chrono::duration<double>(now - last).count();
    last = now;

    timeRealSec += dtReal;

    // Police heat decays in real time (stays stable even if you change timeScale).
    policeHeat = std::max(0.0, policeHeat - dtReal * 0.015);

    // Keep relative mouse mode in sync with our control mode.
    // (Automatically release the mouse when docked so UI interaction is painless.)
    if (docked && mouseSteer) {
      mouseSteer = false;
      SDL_SetRelativeMouseMode(SDL_FALSE);
      SDL_GetRelativeMouseState(nullptr, nullptr); // flush deltas
      toast(toasts, "Mouse steer disabled while docked.", 2.0);
    }

    const SDL_bool wantRelMouse = mouseSteer ? SDL_TRUE : SDL_FALSE;
    if (SDL_GetRelativeMouseMode() != wantRelMouse) {
      SDL_SetRelativeMouseMode(wantRelMouse);
      SDL_GetRelativeMouseState(nullptr, nullptr); // flush deltas
    }


    // Per-frame helper: apply damage from the player to a contact (laser / cannon / projectiles).
    auto playerDamageContact = [&](int idx, double dmg) {
      if (idx < 0 || idx >= (int)contacts.size()) return;
      auto& hit = contacts[(std::size_t)idx];
      if (!hit.alive) return;

      hit.underFireUntilDays = std::max(hit.underFireUntilDays, timeDays + (6.0 / 86400.0));

      // Attacking pirates while they are extorting you voids any "deal".
      if (pirateDemand.active && hit.role == ContactRole::Pirate && hit.groupId != 0 && hit.groupId == pirateDemand.groupId) {
        const std::string src = pirateDemand.leaderName.empty() ? std::string("Pirates") : pirateDemand.leaderName;
        resolvePirateDemand(false, src + ": betrayal! Open fire!");
      }

      // Apply damage
      applyDamage(dmg, hit.shield, hit.hull);

      // VFX: impact sparks at the target.
      if (vfxParticlesEnabled && vfxImpactsEnabled) {
        const math::Vec3d posU = toRenderU(hit.ship.positionKm());
        math::Vec3d n = hit.ship.positionKm() - ship.positionKm();
        if (n.lengthSq() < 1e-12) n = math::Vec3d{0,1,0};
        n = n.normalized();

        const double energy = std::clamp((dmg / 18.0) * (double)vfxParticleIntensity, 0.15, 2.5);
        particles.spawnSparks(posU, n, toRenderU(hit.ship.velocityKmS()), energy);
      }

      // Crimes / reactions (on hit)
      if (hit.alive && hit.role == ContactRole::Trader) {
        hit.fleeUntilDays = timeDays + (180.0 / 86400.0); // flee ~3 minutes
        commitCrime(hit.factionId, 250.0, -5.0, "Assault on trader");
      }
      if (hit.alive && hit.role == ContactRole::Police) {
        hit.hostileToPlayer = true;
        commitCrime(hit.factionId, 600.0, -10.0, "Assault on security");
      }

      // Pirates may disengage when badly hurt (adds a little "morale" to fights).
      if (hit.alive && hit.role == ContactRole::Pirate && !hit.missionTarget) {
        const double hullFrac = (hit.hullMax > 1e-6) ? (hit.hull / hit.hullMax) : 0.0;
        if (hullFrac < 0.25 && timeDays >= hit.fleeUntilDays) {
          hit.fleeUntilDays = timeDays + (110.0 / 86400.0);
        }
      }

      // If destroyed
      if (hit.hull <= 0.0) {
        const core::u64 deadId = hit.id;
        const ContactRole deadRole = hit.role;
        const core::u32 deadFaction = hit.factionId;
        const double lootCr = hit.cargoValueCr;

        hit.alive = false;

	        const math::Vec3d deadPos = hit.ship.positionKm();
	        const math::Vec3d deadVel = hit.ship.velocityKmS();

          // VFX: explosion burst on destruction.
          if (vfxParticlesEnabled && vfxExplosionsEnabled) {
            const double eBase = std::clamp((hit.hullMax / 140.0) * (double)vfxParticleIntensity, 0.45, 2.25);
            particles.spawnExplosion(toRenderU(deadPos), toRenderU(deadVel), eBase);
          }

	        if (deadRole == ContactRole::Pirate) {
	          // Bounty vouchers: redeemed at stations (ties combat -> station loop).
	          const core::u32 authority = currentSystem ? currentSystem->stub.factionId : 0;
	          const double bountyCr = 450.0;
	          if (authority != 0) {
	            addVoucher(authority, bountyCr);
	            toast(toasts, "Pirate destroyed. Bounty voucher +" + std::to_string((int)bountyCr) + " cr", 2.5);
	            addRep(authority, +0.5);
	          } else {
	            credits += bountyCr;
	            toast(toasts, "Pirate destroyed. +" + std::to_string((int)bountyCr) + " cr", 2.5);
	          }

	          // Salvage: pirates drop a few cargo pods.
	          static constexpr econ::CommodityId kPirateLoot[] = {
	            econ::CommodityId::Ore,
	            econ::CommodityId::Metals,
	            econ::CommodityId::Machinery,
	            econ::CommodityId::Electronics,
	            econ::CommodityId::Fuel,
	            econ::CommodityId::Luxury,
	          };
	          const int pods = rng.range<int>(1, 3);
	          for (int i = 0; i < pods; ++i) {
	            const auto cid = kPirateLoot[rng.range<int>(0, (int)std::size(kPirateLoot) - 1)];
	            const double units = (double)rng.range<int>(4, 14);
	            spawnCargoPod(cid, units, deadPos, deadVel, 1.0);
	          }
	        } else if (deadRole == ContactRole::Trader) {
	          // Traders drop their carried commodity as pods (ties piracy -> scooping loop).
	          const double totalUnits = std::max(0.0, hit.tradeUnits * 0.75);
	          if (totalUnits > 0.0) {
	            spawnCargoBurst(hit.tradeCommodity, totalUnits, deadPos, deadVel, rng.range<int>(2, 4));
	          }
	          toast(toasts, "Trader destroyed. Cargo pods ejected. (WANTED!)", 3.0);
	          commitCrime(deadFaction, 1200.0, -18.0, "Murder of trader");
	        } else if (deadRole == ContactRole::Police) {
	          toast(toasts, "Security destroyed. (WANTED!)", 3.0);
	          commitCrime(deadFaction, 2500.0, -35.0, "Murder of security");
	        }

        // Bounty kill missions (pirate targets)
        {
          sim::SaveGame tmp{};
          tmp.credits = credits;
          tmp.missions = missions;

          tmp.reputation.reserve(repByFaction.size());
          for (const auto& kv : repByFaction) {
            tmp.reputation.push_back(sim::FactionReputation{kv.first, kv.second});
          }

          const sim::SystemId sysId = currentSystem ? currentSystem->stub.id : 0;
          const auto res = sim::tryCompleteBountyKill(tmp, sysId, deadId, +2.0);
          if (res.completed > 0) {
            credits = tmp.credits;
            missions = std::move(tmp.missions);

            repByFaction.clear();
            for (const auto& r : tmp.reputation) {
              repByFaction[r.factionId] = r.rep;
            }

            toast(toasts,
                  "Mission complete: bounty target eliminated. +" + std::to_string((int)std::round(res.rewardCr)) + " cr",
                  3.0);
          }
        }
      }
    };

    // Events
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      ImGui_ImplSDL2_ProcessEvent(&event);

      if (event.type == SDL_QUIT) running = false;
      if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_CLOSE) running = false;

      if (event.type == SDL_KEYDOWN && !event.key.repeat) {
        const auto key = [&](const game::KeyChord& chord) { return game::chordMatchesEvent(chord, event.key); };

        // If we're rebinding a control, capture this key press and don't also trigger gameplay actions.
        if (game::handleControlsRebindKeydown(
                controlsWindow, event.key, controls, controlsDirty,
                [&](const std::string& msg, double ttlSec) { toast(toasts, msg, ttlSec); })) {
          continue;
        }

        if (key(controls.actions.commandPalette) && !io.WantCaptureKeyboard) {
          game::openCommandPalette(commandPalette);
        }

        if (key(controls.actions.quit)) running = false;

        if (key(controls.actions.quicksave)) {
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
          s.cargoCapacityKg = cargoCapacityKg;
          s.passengerSeats = passengerSeats;
          s.fuel = fuel;
          s.fuelMax = fuelMax;
          s.fsdRangeLy = fsdRangeLy;
          s.hull = std::clamp(playerHull / std::max(1e-6, playerHullMax), 0.0, 1.0);
          s.shield = std::clamp(playerShield / std::max(1e-6, playerShieldMax), 0.0, 1.0);
          s.heat = std::clamp(heat, 0.0, 200.0);
          // Power distributor
          s.pipsEng = distributorPips.eng;
          s.pipsWep = distributorPips.wep;
          s.pipsSys = distributorPips.sys;
          s.capEngFrac = (distributorCfg.capEng > 1e-9) ? std::clamp(distributorState.eng / distributorCfg.capEng, 0.0, 1.0) : 1.0;
          s.capWepFrac = (distributorCfg.capWep > 1e-9) ? std::clamp(distributorState.wep / distributorCfg.capWep, 0.0, 1.0) : 1.0;
          s.capSysFrac = (distributorCfg.capSys > 1e-9) ? std::clamp(distributorState.sys / distributorCfg.capSys, 0.0, 1.0) : 1.0;

          s.fsdReadyDay = fsdReadyDay;

          // Navigation state (plotted route / auto-run).
          s.navRoute = navRoute;
          s.navRouteHop = (core::u32)std::min<std::size_t>(navRouteHop, (std::size_t)std::numeric_limits<core::u32>::max());
          s.navAutoRun = navAutoRun;
          s.pendingArrivalStation = pendingArrivalTargetStationId;

          // Loadout
          s.shipHull = (core::u8)std::clamp((int)shipHullClass, 0, 2);
          s.thrusterMk = (core::u8)std::clamp(thrusterMk, 1, 3);
          s.shieldMk = (core::u8)std::clamp(shieldMk, 1, 3);
          s.distributorMk = (core::u8)std::clamp(distributorMk, 1, 3);
	          const int maxWeaponIdx = (int)std::size(kWeaponDefs) - 1;
	          s.weaponPrimary = (core::u8)std::clamp((int)weaponPrimary, 0, maxWeaponIdx);
	          s.weaponSecondary = (core::u8)std::clamp((int)weaponSecondary, 0, maxWeaponIdx);
          s.smuggleHoldMk = (core::u8)std::clamp(smuggleHoldMk, 0, 3);

          s.nextMissionId = nextMissionId;
          s.missions = missions;

	          // Escort convoy persistence: store minimal state for the mission-critical convoy
	          // so escort missions survive save/load (including full restart, not just quickload).
	          s.escortConvoys.clear();
	          if (currentSystem) {
	            for (const auto& m : missions) {
	              if (m.type != sim::MissionType::Escort) continue;
	              if (m.completed || m.failed) continue;
	              if (m.leg < 1) continue; // convoy not launched yet
	              if (m.targetNpcId == 0) continue;
	              if (m.fromSystem != currentSystem->stub.id) continue;

	              // Best-effort: locate the actual convoy contact.
	              const Contact* convoy = nullptr;
	              for (const auto& c : contacts) {
	                if (c.id == m.targetNpcId) { convoy = &c; break; }
	              }

	              sim::EscortConvoyState cs{};
	              cs.convoyId = m.targetNpcId;
	              cs.missionId = m.id;
	              cs.systemId = currentSystem->stub.id;
	              cs.fromStation = m.fromStation;
	              cs.toStation = m.toStation;

	              // Defaults if the contact is missing (should be rare, but keeps older saves usable).
	              cs.posKm = ship.positionKm();
	              cs.velKmS = ship.velocityKmS();
	              cs.orient = ship.orientation();
	              cs.angVelRadS = {0,0,0};
	              cs.hullFrac = 1.0;
	              cs.shieldFrac = 1.0;
	              cs.cargoValueCr = 0.0;

	              if (convoy) {
	                cs.posKm = convoy->ship.positionKm();
	                cs.velKmS = convoy->ship.velocityKmS();
	                cs.orient = convoy->ship.orientation();
	                cs.angVelRadS = convoy->ship.angularVelocityRadS();
	                cs.hullFrac = std::clamp(convoy->hull / std::max(1e-6, convoy->hullMax), 0.0, 1.0);
	                cs.shieldFrac = std::clamp(convoy->shield / std::max(1e-6, convoy->shieldMax), 0.0, 1.0);
	                cs.cargoValueCr = convoy->tradeCargoValueCr;
	              }

	              // Restore ambush scheduling / too-far timers where possible.
	              for (const auto& er : escortRuntime) {
	                if (er.missionId == m.id) {
	                  cs.tooFarSec = std::max(0.0, er.tooFarSec);
	                  cs.ambushSpawned = er.ambushSpawned;
	                  cs.nextAmbushDays = er.nextAmbushDays;
	                  break;
	                }
	              }

	              s.escortConvoys.push_back(cs);
	            }
	          }
          // Bounty target persistence (mission-critical NPCs).
          s.bountyTargets = bountyTargetStates;

          s.nextIndustryOrderId = nextIndustryOrderId;
          s.industryOrders = industryOrders;
          s.trackedMissionId = trackedMissionId;

          // Station storage / warehouse
          sim::pruneEmptyStorage(stationStorage);
          s.stationStorage = stationStorage;

          // Mission board cache
          s.missionOffersStationId = missionOffersStationId;
          s.missionOffersDayStamp = missionOffersDayStamp;
          s.missionOffers = missionOffers;

          s.reputation.clear();
          s.reputation.reserve(repByFaction.size());
          for (const auto& kv : repByFaction) {
            sim::FactionReputation r{};
            r.factionId = kv.first;
            r.rep = kv.second;
            s.reputation.push_back(r);
          }
	          s.stationOverrides = universe.exportStationOverrides();

	          // Exploration / law
	          s.explorationDataCr = explorationDataCr;
	          s.scannedKeys.assign(scannedKeys.begin(), scannedKeys.end());

	          // Procedural world persistence (signals / mining depletion)
	          s.resolvedSignalIds.assign(resolvedSignalIds.begin(), resolvedSignalIds.end());
	          std::sort(s.resolvedSignalIds.begin(), s.resolvedSignalIds.end());
	          s.asteroidStates.clear();
	          s.asteroidStates.reserve(asteroidRemainingById.size());
	          for (const auto& kv : asteroidRemainingById) {
	            if ((kv.first & kDeterministicWorldIdBit) == 0) continue;
	            s.asteroidStates.push_back(sim::AsteroidState{kv.first, kv.second});
	          }
	          std::sort(s.asteroidStates.begin(), s.asteroidStates.end(),
	                    [](const sim::AsteroidState& a, const sim::AsteroidState& b) {
	                      return a.asteroidId < b.asteroidId;
	                    });

	          s.bounties.clear();
	          for (const auto& [fid, b] : bountyByFaction) {
	            if (b > 0.0) s.bounties.push_back({fid, b});
	          }
	          s.bountyVouchers.clear();
	          for (const auto& [fid, v] : bountyVoucherByFaction) {
	            if (v > 0.0) s.bountyVouchers.push_back({fid, v});
	          }

          // Background traffic stamps
          s.trafficStamps.clear();
          s.trafficStamps.reserve(trafficDayStampBySystem.size());
          for (const auto& kv : trafficDayStampBySystem) {
            s.trafficStamps.push_back(sim::SystemTrafficStamp{kv.first, kv.second});
          }
          std::sort(s.trafficStamps.begin(), s.trafficStamps.end(),
                    [](const sim::SystemTrafficStamp& a, const sim::SystemTrafficStamp& b) {
                      return a.systemId < b.systemId;
                    });

          if (sim::saveToFile(s, savePath)) {
            toast(toasts, "Saved to " + savePath, 2.5);
          }
        }

        if (key(controls.actions.quickload)) {
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

            cargoCapacityKg = s.cargoCapacityKg;
            passengerSeats = s.passengerSeats;
            fuel = s.fuel;
            fuelMax = s.fuelMax;
            fsdRangeLy = s.fsdRangeLy;
            fsdReadyDay = s.fsdReadyDay;

            // Loadout + derived stats
            shipHullClass = (ShipHullClass)std::clamp((int)s.shipHull, 0, 2);
            thrusterMk = std::clamp((int)s.thrusterMk, 1, 3);
            shieldMk = std::clamp((int)s.shieldMk, 1, 3);
            distributorMk = std::clamp((int)s.distributorMk, 1, 3);
	            const int maxWeaponIdx = (int)std::size(kWeaponDefs) - 1;
	            weaponPrimary = (WeaponType)std::clamp((int)s.weaponPrimary, 0, maxWeaponIdx);
	            weaponSecondary = (WeaponType)std::clamp((int)s.weaponSecondary, 0, maxWeaponIdx);
            smuggleHoldMk = std::clamp((int)s.smuggleHoldMk, 0, 3);
            recalcPlayerStats();

            playerHull = std::clamp(s.hull, 0.0, 1.0) * playerHullMax;
            playerShield = std::clamp(s.shield, 0.0, 1.0) * playerShieldMax;
            heat = std::clamp(s.heat, 0.0, 200.0);
            heatImpulse = 0.0;

            // Power distributor
            distributorPips.eng = std::clamp(s.pipsEng, 0, 4);
            distributorPips.wep = std::clamp(s.pipsWep, 0, 4);
            distributorPips.sys = std::clamp(s.pipsSys, 0, 4);
            sim::normalizePips(distributorPips);

            distributorState.eng = std::clamp(s.capEngFrac, 0.0, 1.0) * distributorCfg.capEng;
            distributorState.wep = std::clamp(s.capWepFrac, 0.0, 1.0) * distributorCfg.capWep;
            distributorState.sys = std::clamp(s.capSysFrac, 0.0, 1.0) * distributorCfg.capSys;


            nextMissionId = s.nextMissionId;
            missions = s.missions;
            nextIndustryOrderId = s.nextIndustryOrderId;
            industryOrders = s.industryOrders;
            trackedMissionId = s.trackedMissionId;

            // Mission board cache
            missionOffersStationId = s.missionOffersStationId;
            missionOffersDayStamp = s.missionOffersDayStamp;
            missionOffers = s.missionOffers;

	            repByFaction.clear();
	            for (const auto& r : s.reputation) repByFaction[r.factionId] = r.rep;

	            // Exploration / law
	            explorationDataCr = s.explorationDataCr;
	            scannedKeys.clear();
	            for (core::u64 k : s.scannedKeys) scannedKeys.insert(k);

	            // Procedural world persistence (signals / mining depletion)
	            resolvedSignalIds.clear();
	            for (core::u64 id : s.resolvedSignalIds) resolvedSignalIds.insert(id);
	            asteroidRemainingById.clear();
	            for (const auto& a : s.asteroidStates) {
	              if ((a.asteroidId & kDeterministicWorldIdBit) == 0) continue;
	              asteroidRemainingById[a.asteroidId] = std::max(0.0, a.remainingUnits);
	            }

	            bountyByFaction.clear();
	            for (const auto& b : s.bounties) bountyByFaction[b.factionId] = b.bountyCr;
	            bountyVoucherByFaction.clear();
	            for (const auto& v : s.bountyVouchers) bountyVoucherByFaction[v.factionId] = v.bountyCr;

                // Background traffic stamps
                trafficDayStampBySystem.clear();
                for (const auto& t : s.trafficStamps) {
                  trafficDayStampBySystem[t.systemId] = t.dayStamp;
                }

                // Station storage / warehouse
                stationStorage = s.stationStorage;
                sim::pruneEmptyStorage(stationStorage);
	            policeAlertUntilDays = 0.0;
	            policeDemand = PoliceDemand{};

            docked = (s.dockedStation != 0);
            dockedStationId = s.dockedStation;
            selectedStationIndex = 0;
            if (docked) {
              for (std::size_t i = 0; i < sys.stations.size(); ++i) {
                if (sys.stations[i].id == s.dockedStation) selectedStationIndex = (int)i;
              }
            }

            // clear transient runtime things
            contacts.clear();
	            escortRuntime.clear();
            bountyTargetStates = s.bountyTargets;
            beams.clear();
            projectiles.clear();
            missiles.clear();

            // VFX
            particles.clear();

            floatingCargo.clear();
            asteroids.clear();
            signals.clear();
            nextWorldObjectId = 1;
            nextSignalSpawnDays = timeDays + 0.01;
            nextSupercruiseShadowSpawnDays = timeDays + 0.01;
            cargoScoopDeployed = true;
            cargoFullToastCooldownUntilDays = 0.0;

            // Reset local-space encounter scheduling after load.
            // (Contacts are transient, so this is purely cadence/QA convenience.)
            encounterDirector = sim::makeEncounterDirector(seed, timeDays);
            encounterDirector.nextPirateSpawnDays = timeDays + (encounterDirector.rng.range(25.0, 55.0) / 86400.0);
            encounterDirector.nextTraderSpawnDays = timeDays + (encounterDirector.rng.range(15.0, 35.0) / 86400.0);
            encounterDirector.nextPoliceSpawnDays = timeDays + (encounterDirector.rng.range(10.0, 25.0) / 86400.0);
            autopilot = false;
            supercruiseState = SupercruiseState::Idle;
            supercruiseChargeRemainingSec = 0.0;
            supercruiseCooldownRemainingSec = 0.0;
      supercruiseCooldownTotalSec = 0.0;
            supercruiseDropRequested = false;
            supercruiseSafeDropReady = false;
            supercruiseTtaSec = 0.0;
            supercruiseDistKm = 0.0;
            supercruiseClosingKmS = 0.0;
            fsdState = FsdState::Idle;
            fsdTargetSystem = 0;
            // Navigation state (plotted route / auto-run / pending arrival target).
            navRoute = s.navRoute;
            navRouteHop = (std::size_t)s.navRouteHop;
            navAutoRun = s.navAutoRun;
            pendingArrivalTargetStationId = s.pendingArrivalStation;

            // Validate/repair route cursor so it matches the loaded current system.
            if (currentSystem) {
              if (navRoute.size() >= 2) {
                std::size_t idx = navRoute.size();
                for (std::size_t i = 0; i < navRoute.size(); ++i) {
                  if (navRoute[i] == currentSystem->stub.id) { idx = i; break; }
                }
                if (idx == navRoute.size()) {
                  navRoute.clear();
                  navRouteHop = 0;
                  navAutoRun = false;
                } else {
                  navRouteHop = idx;
                }
              } else {
                navRoute.clear();
                navRouteHop = 0;
                navAutoRun = false;
              }

              if (!navRoute.empty() && navRouteHop + 1 >= navRoute.size()) {
                navAutoRun = false;
              }
            } else {
              navRoute.clear();
              navRouteHop = 0;
              navAutoRun = false;
            }

            scanning = false;
            scanProgressSec = 0.0;
            scanLockedId = 0;
            scanLabel.clear();
            scanLockedTarget = Target{};
            clearances.clear();
            target = Target{};

            if (!navRoute.empty() && navRouteHop + 1 < navRoute.size()) {
              galaxySelectedSystemId = navRoute[navRouteHop + 1];
            } else {
              galaxySelectedSystemId = currentSystem->stub.id;
            }

            // If we already are in the system for a pending arrival station, target it immediately.
            if (pendingArrivalTargetStationId != 0 && currentSystem) {
              for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
                if (currentSystem->stations[i].id == pendingArrivalTargetStationId) {
                  target.kind = TargetKind::Station;
                  target.index = i;
                  selectedStationIndex = (int)i;
                  pendingArrivalTargetStationId = 0;
                  break;
                }
              }
            }


            // Re-seed deterministic in-system objects so signals/asteroids are present after load.
            seedSystemSpaceObjects(*currentSystem);

	            // Restore escort mission convoys (mission-critical NPCs) after a load.
	            //
	            // Without this, any save/load during an escort would recreate the world without the
	            // convoy contact (contacts are transient), causing the mission to instantly fail.
	            {
	              auto findSavedConvoy = [&](core::u64 convoyId) -> const sim::EscortConvoyState* {
	                for (const auto& c : s.escortConvoys) {
	                  if (c.convoyId == convoyId) return &c;
	                }
	                return nullptr;
	              };

	              auto findStationIndexById = [&](sim::StationId id) -> std::optional<std::size_t> {
	                if (!currentSystem) return std::nullopt;
	                for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
	                  if (currentSystem->stations[i].id == id) return i;
	                }
	                return std::nullopt;
	              };

	              auto findContactById = [&](core::u64 id) -> Contact* {
	                for (auto& c : contacts) {
	                  if (c.id == id) return &c;
	                }
	                return nullptr;
	              };

	              auto spawnEscortConvoy = [&](const sim::Mission& m, const sim::EscortConvoyState* saved) {
	                if (!currentSystem) return;
	
	                const auto fromIdx = findStationIndexById(m.fromStation);
	                const auto toIdx = findStationIndexById(m.toStation);
	                if (!fromIdx || !toIdx) return;
	
		                const auto& originSt = currentSystem->stations[*fromIdx];
	
	                Contact convoy{};
	                convoy.alive = true;
	                convoy.id = m.targetNpcId;
	                convoy.role = ContactRole::Trader;
	                convoy.factionId = m.factionId;
	                convoy.escortConvoy = true;
	                convoy.groupId = convoy.id;
	                convoy.leaderId = 0;
	
	                if (m.units > 0.0) {
	                  convoy.name = std::string("Convoy - ") + econ::commodityDef(m.commodity).name;
	                } else {
	                  convoy.name = "Convoy";
	                }
	
	                configureContactLoadout(convoy,
	                                       ShipHullClass::Hauler,
	                                       /*thrMk=*/1,
	                                       /*shieldMk=*/2,
	                                       /*distMk=*/1,
	                                       /*weapon=*/WeaponType::MiningLaser,
	                                       /*hullMul=*/0.90,
	                                       /*shieldMul=*/0.90,
	                                       /*regenMul=*/1.0,
	                                       /*accelMul=*/0.72,
	                                       /*aiSkill=*/0.50);
	
	                convoy.ship.setMassKg(20000.0);
	
	                if (saved) {
	                  convoy.ship.setPositionKm(saved->posKm);
	                  convoy.ship.setVelocityKmS(saved->velKmS);
	                  convoy.ship.setOrientation(saved->orient);
	                  convoy.ship.setAngularVelocityRadS(saved->angVelRadS);
	                  convoy.hull = convoy.hullMax * std::clamp(saved->hullFrac, 0.0, 1.0);
	                  convoy.shield = convoy.shieldMax * std::clamp(saved->shieldFrac, 0.0, 1.0);
	                  convoy.tradeCargoValueCr = saved->cargoValueCr;
	                  convoy.cargoValueCr = convoy.tradeCargoValueCr;
	                } else {
	                  // Fallback spawn for older saves that predate escort_convoys persistence.
	                  // Put the convoy near the player so the escort can be resumed immediately.
	                  convoy.ship.setPositionKm(ship.positionKm() + randUnit() * rng.range(70000.0, 110000.0));
	                  convoy.ship.setVelocityKmS(ship.velocityKmS());
	                  convoy.hull = convoy.hullMax;
	                  convoy.shield = convoy.shieldMax;
	                }
	
	                // Mission haul: destination is fixed and cargo is already reserved at acceptance.
	                convoy.tradeDestStationIndex = (int)*toIdx;
	                convoy.tradeCommodity = m.commodity;
	                convoy.tradeUnits = m.units;
	                const double cargoMassKg = std::max(0.0, m.units) * econ::commodityDef(m.commodity).massKg;
	                convoy.tradeCapacityKg = std::max(240.0, cargoMassKg * 1.15);
	                convoy.tradeCooldownUntilDays = timeDays + (20.0 / 86400.0);
	                convoy.tradeSupercruiseSpeedKmS = 0.0;
	
	                if (m.units > 0.0 && convoy.tradeCargoValueCr <= 0.0) {
	                  auto& originEcon = universe.stationEconomy(originSt, timeDays);
	                  const econ::MarketQuote q = econ::quote(originEcon, originSt.economyModel, m.commodity);
	                  convoy.tradeCargoValueCr = q.mid * m.units;
	                  convoy.cargoValueCr = convoy.tradeCargoValueCr;
	                }
	
	                contacts.push_back(std::move(convoy));

	                // Spawn a small police escort wing that follows the convoy.
	                const int escortCount = 2;
	                const math::Vec3d basePos = contacts.back().ship.positionKm();
	                const math::Vec3d baseVel = contacts.back().ship.velocityKmS();
	                for (int e = 0; e < escortCount; ++e) {
	                  Contact escort{};
	                  escort.alive = true;
	                  escort.id = allocWorldId();
	                  escort.role = ContactRole::Police;
	                  escort.factionId = m.factionId;
	                  escort.name = "Escort";
	                  escort.followId = m.targetNpcId;
	                  escort.groupId = m.targetNpcId;
	                  escort.leaderId = m.targetNpcId;
	                  configureContactLoadout(escort,
	                                         ShipHullClass::Fighter,
	                                         /*thrMk=*/1,
	                                         /*shieldMk=*/1,
	                                         /*distMk=*/1,
	                                         /*weapon=*/WeaponType::PulseLaser,
	                                         /*hullMul=*/1.10,
	                                         /*shieldMul=*/0.80,
	                                         /*regenMul=*/1.0,
	                                         /*accelMul=*/0.75,
	                                         /*aiSkill=*/0.65);
	                  escort.ship.setMassKg(13000.0);
	                  escort.ship.setPositionKm(basePos + randUnit() * rng.range(2500.0, 4200.0));
	                  escort.ship.setVelocityKmS(baseVel);
	                  contacts.push_back(std::move(escort));
	                }
	              };

	              // Rebuild transient escort mission runtime from the save file.
	              // (The mission list itself is already loaded into `missions` above.)
	              for (const auto& m : missions) {
	                if (m.type != sim::MissionType::Escort) continue;
	                if (m.completed || m.failed) continue;
	                if (m.leg < 1) continue;
	                if (!currentSystem) continue;
	                if (m.fromSystem != currentSystem->stub.id) continue;
	                if (m.targetNpcId == 0) continue;
	
	                const sim::EscortConvoyState* saved = findSavedConvoy(m.targetNpcId);
	
	                // Ensure the convoy contact exists.
	                if (!findContactById(m.targetNpcId)) {
	                  spawnEscortConvoy(m, saved);
	                }
	
	                EscortRuntime er{};
	                er.missionId = m.id;
	                er.convoyId = m.targetNpcId;
	                er.pirateGroupId = (m.targetNpcId != 0)
	                    ? (m.targetNpcId ^ 0x5a5a5a5a5a5a5a5aull)
	                    : std::max<core::u64>(1, rng.nextU64());
	
	                if (saved) {
	                  er.tooFarSec = std::max(0.0, saved->tooFarSec);
	                  er.ambushSpawned = saved->ambushSpawned;
	                  if (saved->nextAmbushDays > 0.0) {
	                    er.nextAmbushDays = saved->nextAmbushDays;
	                  } else {
	                    er.nextAmbushDays = timeDays + (rng.range(70.0, 150.0) / 86400.0);
	                  }
	                } else {
	                  er.tooFarSec = 0.0;
	                  er.ambushSpawned = false;
	                  er.nextAmbushDays = timeDays + (rng.range(70.0, 150.0) / 86400.0);
	                }
	
	                escortRuntime.push_back(er);
	              }
	            }

            toast(toasts, "Loaded " + savePath, 2.5);
          }
        }


        // Window toggles
        if (key(controls.actions.toggleGalaxy) && !io.WantCaptureKeyboard) showGalaxy = !showGalaxy;
        if (key(controls.actions.toggleShip) && !io.WantCaptureKeyboard) showShip = !showShip;
        if (key(controls.actions.toggleMarket) && !io.WantCaptureKeyboard) showEconomy = !showEconomy;
        if (key(controls.actions.toggleMissions) && !io.WantCaptureKeyboard) showMissions = !showMissions;
        if (key(controls.actions.toggleContacts) && !io.WantCaptureKeyboard) showContacts = !showContacts;
        if (key(controls.actions.toggleScanner) && !io.WantCaptureKeyboard) showScanner = !showScanner;
        if (key(controls.actions.toggleTrade) && !io.WantCaptureKeyboard) showTrade = !showTrade;
        if (key(controls.actions.toggleGuide) && !io.WantCaptureKeyboard) showGuide = !showGuide;
        if (key(controls.actions.toggleHangar) && !io.WantCaptureKeyboard) showHangar = !showHangar;
        if (key(controls.actions.toggleWorldVisuals) && !io.WantCaptureKeyboard) showWorldVisuals = !showWorldVisuals;
        if (key(controls.actions.toggleSpriteLab) && !io.WantCaptureKeyboard) showSprites = !showSprites;
        if (key(controls.actions.toggleVfxLab) && !io.WantCaptureKeyboard) showVfx = !showVfx;
        if (key(controls.actions.togglePostFx) && !io.WantCaptureKeyboard) showPostFx = !showPostFx;
        if (key(controls.actions.toggleControlsWindow) && !io.WantCaptureKeyboard) {
          controlsWindow.open = !controlsWindow.open;
          if (controlsWindow.open) controlsWindow.focusFilter = true;
        }

        // HUD layout editor shortcuts
        if (key(controls.actions.hudLayoutToggleEdit) && !io.WantCaptureKeyboard) {
          hudLayoutEditMode = !hudLayoutEditMode;
          showHudLayoutWindow = true;
          toast(toasts,
                std::string("HUD layout edit ") + (hudLayoutEditMode ? "ON" : "OFF") + " (" +
                  game::chordLabel(controls.actions.hudLayoutToggleEdit) + ")",
                1.8);
        }

        if (key(controls.actions.hudLayoutSave) && !io.WantCaptureKeyboard) {
          // Sync enabled toggles into the layout before saving.
          hudLayout.widget(ui::HudWidgetId::Radar).enabled = showRadarHud;
          hudLayout.widget(ui::HudWidgetId::Objective).enabled = objectiveHudEnabled;
          hudLayout.widget(ui::HudWidgetId::Threat).enabled = hudThreatOverlayEnabled;
          hudLayout.widget(ui::HudWidgetId::Jump).enabled = hudJumpOverlay;

          if (ui::saveToFile(hudLayout, hudLayoutPath)) {
            toast(toasts, "Saved HUD layout to " + hudLayoutPath, 2.0);
          } else {
            toast(toasts, "Failed to save HUD layout.", 2.0);
          }
        }

        if (key(controls.actions.hudLayoutLoad) && !io.WantCaptureKeyboard) {
          ui::HudLayout loaded = ui::makeDefaultHudLayout();
          if (ui::loadFromFile(hudLayoutPath, loaded)) {
            hudLayout = loaded;
            showRadarHud = hudLayout.widget(ui::HudWidgetId::Radar).enabled;
            objectiveHudEnabled = hudLayout.widget(ui::HudWidgetId::Objective).enabled;
            hudThreatOverlayEnabled = hudLayout.widget(ui::HudWidgetId::Threat).enabled;
            hudJumpOverlay = hudLayout.widget(ui::HudWidgetId::Jump).enabled;
            toast(toasts, "Loaded HUD layout from " + hudLayoutPath, 2.0);
          } else {
            toast(toasts, "HUD layout file not found (using defaults).", 2.0);
          }
          showHudLayoutWindow = true;
        }

        if (key(controls.actions.hudLayoutReset) && !io.WantCaptureKeyboard) {
          hudLayout = ui::makeDefaultHudLayout();
          showRadarHud = hudLayout.widget(ui::HudWidgetId::Radar).enabled;
          objectiveHudEnabled = hudLayout.widget(ui::HudWidgetId::Objective).enabled;
          hudThreatOverlayEnabled = hudLayout.widget(ui::HudWidgetId::Threat).enabled;
          hudJumpOverlay = hudLayout.widget(ui::HudWidgetId::Jump).enabled;
          toast(toasts, "HUD layout reset to defaults.", 2.0);
          showHudLayoutWindow = true;
        }

        if (key(controls.actions.toggleRadarHud) && !io.WantCaptureKeyboard) {
          showRadarHud = !showRadarHud;
          toast(toasts,
                std::string("Radar HUD ") + (showRadarHud ? "ON" : "OFF") + " (" +
                  game::chordLabel(controls.actions.toggleRadarHud) + ")",
                1.6);
        }

        if (key(controls.actions.toggleTacticalOverlay) && !io.WantCaptureKeyboard) {
          showTacticalOverlay = !showTacticalOverlay;
          toast(toasts,
                std::string("Tactical overlay ") + (showTacticalOverlay ? "ON" : "OFF") + " (" +
                  game::chordLabel(controls.actions.toggleTacticalOverlay) + ")",
                1.6);
        }

        if (key(controls.actions.pause) && !io.WantCaptureKeyboard) paused = !paused;

        if (key(controls.actions.toggleAutopilot) && !io.WantCaptureKeyboard) {
          if (autopilot) {
            autopilot = false;
            dockingComputer.reset();
            toast(toasts, "Docking computer disengaged.", 1.6);
          } else {
            if (docked) {
              toast(toasts, "Already docked.", 1.6);
            } else if (!(currentSystem && target.kind == TargetKind::Station && target.index < currentSystem->stations.size())) {
              toast(toasts, "Target a station to engage docking computer.", 2.0);
            } else if (supercruiseState != SupercruiseState::Idle || fsdState != FsdState::Idle) {
              toast(toasts, "Cannot engage docking computer while in supercruise/FSD.", 2.0);
            } else {
              autopilot = true;
              dockingComputer.reset();
              toast(toasts, "Docking computer engaged.", 1.6);
            }
          }
        }

        if (key(controls.actions.toggleCargoScoop) && !io.WantCaptureKeyboard) {
          cargoScoopDeployed = !cargoScoopDeployed;
          toast(toasts,
                std::string("Cargo scoop ") + (cargoScoopDeployed ? "DEPLOYED" : "RETRACTED") + " (" +
                  game::chordLabel(controls.actions.toggleCargoScoop) + ")",
                1.8);
        }

        if (key(controls.actions.cycleTargets) && !io.WantCaptureKeyboard && !docked && currentSystem) {
	            // Cycle targets across signal sources, cargo pods, and asteroid nodes.
	            struct Cand {
	              TargetKind kind;
	              std::size_t idx;
	              double distKm;
	            };
	            std::vector<Cand> cands;
	            cands.reserve(signals.size() + floatingCargo.size() + asteroids.size());
	            const math::Vec3d p = ship.positionKm();
	            for (std::size_t i = 0; i < signals.size(); ++i) {
	              if (timeDays > signals[i].expireDay) continue;
	              cands.push_back({TargetKind::Signal, i, (signals[i].posKm - p).length()});
	            }
	            for (std::size_t i = 0; i < floatingCargo.size(); ++i) {
	              cands.push_back({TargetKind::Cargo, i, (floatingCargo[i].posKm - p).length()});
	            }
	            for (std::size_t i = 0; i < asteroids.size(); ++i) {
	              if (asteroids[i].remainingUnits <= 0.0) continue;
	              cands.push_back({TargetKind::Asteroid, i, (asteroids[i].posKm - p).length()});
	            }
	            if (cands.empty()) {
	              toast(toasts, "No signal sources / salvage / asteroids detected.", 1.8);
	            } else {
	              std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return a.distKm < b.distKm; });
	              std::size_t pick = 0;
	              for (std::size_t i = 0; i < cands.size(); ++i) {
	                if (target.kind == cands[i].kind && target.index == cands[i].idx) {
	                  pick = (i + 1) % cands.size();
	                  break;
	                }
	              }
	              target.kind = cands[pick].kind;
	              target.index = cands[pick].idx;
	            }
	          }
	        }

	        if (key(controls.actions.complyOrSubmit)) {
          if (!io.WantCaptureKeyboard && !docked && fsdState == FsdState::Idle && supercruiseState == SupercruiseState::Idle) {
            // Contraband scan resolution takes priority over bounty submission.
            if (bribeOffer.active && timeDays < bribeOffer.untilDays) {
              enforceContraband(bribeOffer.factionId,
                                bribeOffer.sourceName.empty() ? "Security" : bribeOffer.sourceName,
                                bribeOffer.illegalValueCr,
                                bribeOffer.detail,
                                bribeOffer.scannedIllegal);
              bribeOffer = BribeOffer{};
            } else if (policeDemand.active && timeDays < policeDemand.untilDays) {
              const core::u32 fid = policeDemand.factionId;
              const double owed = getBounty(fid);
              const double pay = std::min(credits, owed);
              if (pay > 0.0) {
                credits -= pay;
                addBounty(fid, -pay);
              }
              const double remaining = getBounty(fid);
              if (remaining <= 0.0) {
                clearBounty(fid);
                policeAlertUntilDays = timeDays;
                toast(toasts, "Submitted. Bounty paid: CLEAN.", 2.4);
                policeDemand = PoliceDemand{};
              } else {
                toast(toasts,
                      "Submitted. Partial payment " + std::to_string((int)pay) + " cr. Remaining bounty " + std::to_string((int)remaining) + " cr.",
                      3.0);
                // Unable to settle: authorities will escalate immediately.
                policeDemand.active = false;
                policeDemand.untilDays = timeDays;
                policeDemand.amountCr = remaining;
              }
            } else {
              toast(toasts, "No active authority demand to submit to.", 1.8);
            }
          }
        }

        if (key(controls.actions.bribe)) {
          if (!io.WantCaptureKeyboard && !docked && fsdState == FsdState::Idle && supercruiseState == SupercruiseState::Idle) {
            if (bribeOffer.active && timeDays < bribeOffer.untilDays) {
              if (credits + 1e-6 < bribeOffer.amountCr) {
                toast(toasts,
                      "Insufficient credits for bribe (" + std::to_string((int)std::round(bribeOffer.amountCr)) + " cr).",
                      2.6);
              } else {
                credits -= bribeOffer.amountCr;

                // Bribery still generates some local alert, but avoids confiscation + formal bounty.
                policeHeat = std::clamp(policeHeat + std::min(0.9, 0.35 + bribeOffer.illegalValueCr / 12000.0), 0.0, 6.0);
                addRep(bribeOffer.factionId, -1.0);

                // If the player somehow lost the smuggle package while the offer was up, fail it.
                int failedSmuggle = 0;
                for (auto& m : missions) {
                  if (m.completed || m.failed) continue;
                  if (m.type != sim::MissionType::Smuggle) continue;
                  const std::size_t mi = (std::size_t)m.commodity;
                  if (mi >= econ::kCommodityCount) continue;
                  if (bribeOffer.scannedIllegal[mi] <= 1e-6) continue;
                  if (cargo[mi] + 1e-6 < m.units) {
                    m.failed = true;
                    addRep(m.factionId, -4.0);
                    ++failedSmuggle;
                  }
                }
                if (failedSmuggle > 0) {
                  toast(toasts, "Smuggle mission failed: package missing.", 3.2);
                }

                toast(toasts,
                      (bribeOffer.sourceName.empty() ? "Security" : bribeOffer.sourceName) +
                        std::string(": bribe accepted. Cargo scan waived."),
                      3.2);
                bribeOffer = BribeOffer{};
              }
            } else {
              toast(toasts, "No active bribe offer.", 1.8);
            }
          }
        }

        if (key(controls.actions.toggleMouseSteer)) {
          if (!io.WantCaptureKeyboard) {
            mouseSteer = !mouseSteer;
            SDL_SetRelativeMouseMode(mouseSteer ? SDL_TRUE : SDL_FALSE);
            SDL_GetRelativeMouseState(nullptr, nullptr); // flush deltas
            toast(toasts,
                  std::string("Mouse steer ") + (mouseSteer ? "ON" : "OFF") + " (" +
                    game::chordLabel(controls.actions.toggleMouseSteer) + ")",
                  1.8);
          }
        }

        if (key(controls.actions.supercruise)) {
          if (!io.WantCaptureKeyboard && !docked && fsdState == FsdState::Idle) {
            if (supercruiseState == SupercruiseState::Idle) {
              if (target.kind == TargetKind::None) {
                toast(toasts, "No nav target set.", 1.8);
              } else {
                // Soft warning: you *can* engage, but pirates nearby make interdiction more likely.
                double nearestPirateKm = 1e99;
                for (const auto& c : contacts) {
                  if (!c.alive || c.role != ContactRole::Pirate) continue;
                  nearestPirateKm = std::min(nearestPirateKm, (c.ship.positionKm() - ship.positionKm()).length());
                }
                const bool pirateNearby = (nearestPirateKm < 160000.0);

                supercruiseState = SupercruiseState::Charging;
                supercruiseChargeRemainingSec = kSupercruiseChargeSec;
                supercruiseDropRequested = false;
                interdiction = sim::InterdictionState{};
                interdictionSubmitRequested = false;
                interdictionPirateName.clear();
                interdictionPirateStrength = 1.0;
                autopilot = false;
                dockingComputer.reset();
                scanning = false;
                scanProgressSec = 0.0;
                navAutoRun = false;
                toast(toasts,
                      pirateNearby ? "Supercruise charging... (pirate signals nearby)" : "Supercruise charging...",
                      1.8);
              }
            } else if (supercruiseState == SupercruiseState::Charging) {
              supercruiseState = SupercruiseState::Idle;
              supercruiseChargeRemainingSec = 0.0;
              interdiction = sim::InterdictionState{};
              interdictionSubmitRequested = false;
              interdictionPirateName.clear();
              interdictionPirateStrength = 1.0;
              toast(toasts, "Supercruise canceled.", 1.4);
            } else if (supercruiseState == SupercruiseState::Active) {
              if (sim::interdictionInProgress(interdiction)) {
                interdictionSubmitRequested = true;
                toast(toasts, "Submitted to interdiction.", 1.4);
              } else {
                supercruiseDropRequested = true;
                toast(toasts, "Drop requested.", 1.2);
              }
            } else if (supercruiseState == SupercruiseState::Cooldown) {
              toast(toasts, "Supercruise cooling down...", 1.6);
            }
          }
        }

        if (key(controls.actions.scannerAction)) {
            // Scanner action (mission scans + exploration scans)
            if (scanning) {
              scanning = false;
              scanProgressSec = 0.0;
              scanLockedId = 0;
              scanLabel.clear();
              toast(toasts, "Scan cancelled.", 1.5);
            } else {
              if (docked) {
                toast(toasts, "Undock to scan.", 1.8);
              } else if (supercruiseState != SupercruiseState::Idle || fsdState != FsdState::Idle) {
                toast(toasts, "Scanning unavailable in supercruise / hyperspace.", 2.0);
              } else if (target.kind == TargetKind::None) {
                {
                  std::string hint;
                  hint.reserve(128);
                  hint += "No scan target selected. (Use ";
                  hint += game::chordLabel(controls.actions.targetStationCycle);
                  hint += "/";
                  hint += game::chordLabel(controls.actions.targetPlanetCycle);
                  hint += "/";
                  hint += game::chordLabel(controls.actions.targetContactCycle);
                  hint += "/";
                  hint += game::chordLabel(controls.actions.targetStar);
                  hint += " or System Scanner)";
                  toast(toasts, hint, 2.5);
                }
              } else {
                // Lock scan to current target so cycling targets cancels the scan.
                scanLockedTarget = target;
                scanLockedId = 0;
                scanLabel.clear();
                scanProgressSec = 0.0;

                bool ok = false;

                if (target.kind == TargetKind::Contact && target.index < contacts.size()) {
                  const auto& c = contacts[target.index];
                  if (c.alive) {
                    scanLockedId = c.id;
                    scanDurationSec = 4.0;
                    scanRangeKm = 85000.0;
                    scanLabel = "Contact scan: " + c.name;
                    ok = true;
                  }
                } else if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
                  const auto& st = currentSystem->stations[target.index];
                  scanLockedId = st.id;
                  scanDurationSec = 3.0;
                  scanRangeKm = std::max(25000.0, st.commsRangeKm * 0.9);
                  scanLabel = "Station scan: " + st.name;
                  ok = true;
                } else if (target.kind == TargetKind::Planet && target.index < currentSystem->planets.size()) {
                  const auto& p = currentSystem->planets[target.index];
                  scanLockedId = (core::u64)target.index;
                  scanDurationSec = 5.0;
                  const double rKm = p.radiusEarth * kEARTH_RADIUS_KM;
                  scanRangeKm = std::max(200000.0, rKm * 45.0);
                  scanLabel = "Planet scan: " + p.name;
                  ok = true;
                } else if (target.kind == TargetKind::Star) {
                  scanLockedId = currentSystem->stub.id;
                  scanDurationSec = 2.5;
                  scanRangeKm = 1.0e18; // anywhere in-system for now
                  scanLabel = std::string("Star scan: ") + starClassName(currentSystem->star.cls);
                  ok = true;
	                } else if (target.kind == TargetKind::Signal && target.index < signals.size()) {
	                  const auto& s = signals[target.index];
	                  scanLockedId = s.id;
	                  scanDurationSec = (s.type == SignalType::Resource) ? 3.5 : 4.0;
	                  scanRangeKm = 200000.0;
	                  const char* base = signalTypeName(s.type);
	                  if (s.type == SignalType::Resource && s.hasResourcePlan) {
	                    base = sim::resourceFieldKindName(s.resource.kind);
	                  }
	                  scanLabel = std::string("Signal scan: ") + base;
	                  ok = true;
	                } else if (target.kind == TargetKind::Asteroid && target.index < asteroids.size()) {
	                  const auto& a = asteroids[target.index];
	                  scanLockedId = a.id;
	                  scanDurationSec = 4.5;
	                  scanRangeKm = std::max(140000.0, a.radiusKm * 55.0);
	                  scanLabel = std::string("Prospect asteroid: ") + econ::commodityDef(a.yield).name;
	                  ok = true;
                }

                if (ok) {
                  scanning = true;
                  toast(toasts, scanLabel + " (hold steady)...", 2.0);
                } else {
                  toast(toasts, "Invalid scan target.", 2.0);
                }
              }
            }
          }

        if (key(controls.actions.fsdJump)) {
          // If we're mid-charge, allow a clean abort (fuel is consumed when charge completes).
          if (fsdState == FsdState::Charging) {
            fsdState = FsdState::Idle;
            fsdTargetSystem = 0;
            fsdChargeRemainingSec = 0.0;
            fsdTravelRemainingSec = 0.0;
            fsdTravelTotalSec = 0.0;
            fsdFuelCost = 0.0;
            fsdJumpDistanceLy = 0.0;
            navAutoRun = false;
            toast(toasts, "FSD canceled.", 1.6);
          } else {
            // Jump to next hop in the plotted route, otherwise jump to selected system.
            if (!navRoute.empty() && navRouteHop + 1 < navRoute.size()) {
              startFsdJumpTo(navRoute[navRouteHop + 1]);
            } else {
              startFsdJumpTo(galaxySelectedSystemId);
            }
          }
        }
        if (key(controls.actions.targetStationCycle) && !io.WantCaptureKeyboard) {
  // cycle station targets
  if (scanning) { scanning = false; scanProgressSec = 0.0; scanLockedId = 0; scanLabel.clear(); }
  if (!currentSystem->stations.empty()) {
    if (target.kind != TargetKind::Station) {
      target.kind = TargetKind::Station;
      target.index = 0;
    } else {
      target.index = (target.index + 1) % currentSystem->stations.size();
    }
  }
}

        if (key(controls.actions.targetPlanetCycle) && !io.WantCaptureKeyboard) {
  // cycle planet targets
  if (scanning) { scanning = false; scanProgressSec = 0.0; scanLockedId = 0; scanLabel.clear(); }
  if (!currentSystem->planets.empty()) {
    if (target.kind != TargetKind::Planet) {
      target.kind = TargetKind::Planet;
      target.index = 0;
    } else {
      target.index = (target.index + 1) % currentSystem->planets.size();
    }
  }
}

        if (key(controls.actions.targetContactCycle) && !io.WantCaptureKeyboard) {
  // cycle contact targets (alive only)
  if (scanning) { scanning = false; scanProgressSec = 0.0; scanLockedId = 0; scanLabel.clear(); }
  if (!contacts.empty()) {
    std::size_t start = 0;
    if (target.kind == TargetKind::Contact && target.index < contacts.size()) {
      start = (target.index + 1) % contacts.size();
    }
    std::size_t idx = start;
    for (std::size_t step = 0; step < contacts.size(); ++step) {
      if (contacts[idx].alive) {
        target.kind = TargetKind::Contact;
        target.index = idx;
        break;
      }
      idx = (idx + 1) % contacts.size();
    }
  }
}

        if (key(controls.actions.targetStar) && !io.WantCaptureKeyboard) {
  // target system primary star (for scanning / nav reference)
  if (scanning) { scanning = false; scanProgressSec = 0.0; scanLockedId = 0; scanLabel.clear(); }
  target.kind = TargetKind::Star;
  target.index = 0;
}

        if (key(controls.actions.clearTarget) && !io.WantCaptureKeyboard) {
  if (scanning) { scanning = false; scanProgressSec = 0.0; scanLockedId = 0; scanLabel.clear(); }
  target = Target{};
}

        if (key(controls.actions.requestDockingClearance) && !io.WantCaptureKeyboard) {
          // Request docking clearance from targeted station
          if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
            const auto& st = currentSystem->stations[target.index];
            const math::Vec3d stPos = stationPosKm(st, timeDays);
            const double dist = (ship.positionKm() - stPos).length();

            auto& cs = clearances[st.id];

            const bool hadValid = sim::dockingClearanceValid(cs, timeDays);

            // Optional congestion hint: count nearby ships in comms range (excluding the player).
            int nearby = 0;
            for (const auto& c : contacts) {
              if (!c.alive) continue;
              if ((c.ship.positionKm() - stPos).length() <= st.commsRangeKm) ++nearby;
            }

            const double traffic01 = sim::estimateDockingTraffic01(seed, st, timeDays, nearby);
            const auto dec = sim::requestDockingClearance(seed, st, timeDays, dist, cs, traffic01);

            if (hadValid && dec.status == sim::DockingClearanceStatus::Granted) {
              toast(toasts, "Docking clearance already granted.", 2.0);
            } else if (dec.status == sim::DockingClearanceStatus::OutOfRange) {
              toast(toasts, "Out of comms range for clearance.", 2.5);
            } else if (dec.status == sim::DockingClearanceStatus::Throttled) {
              toast(toasts, "Clearance channel busy. Try again soon.", 2.5);
            } else if (dec.status == sim::DockingClearanceStatus::Granted) {
              toast(toasts, "Docking clearance GRANTED.", 3.0);
            } else {
              toast(toasts, "Docking clearance DENIED. (traffic)", 3.0);
            }
          } else {
            toast(toasts, "No station targeted for clearance.", 2.0);
          }
        }

        if (key(controls.actions.dockOrUndock) && !io.WantCaptureKeyboard) {
          if (supercruiseState != SupercruiseState::Idle) {
            toast(toasts, "Cannot dock while in supercruise.", 2.0);
          } else if (docked) {
            // Undock: place ship just outside the slot
            const auto it = std::find_if(currentSystem->stations.begin(), currentSystem->stations.end(),
                                         [&](const sim::Station& s){ return s.id == dockedStationId; });
            if (it != currentSystem->stations.end()) {
              const auto& st = *it;
              const math::Vec3d stPos = stationPosKm(st, timeDays);
              const math::Quatd stQ = stationOrient(st, stPos, timeDays);
              const math::Vec3d axis = stQ.rotate({0,0,1});
              ship.setPositionKm(stPos + axis * (st.radiusKm * 1.8));
              ship.setVelocityKmS(stationVelKmS(st, timeDays));
              ship.setAngularVelocityRadS({0,0,0});
              ship.setOrientation(quatFromTo({0,0,1}, -axis));

              const sim::StationId undockedFromStationId = dockedStationId;
              docked = false;
              dockedStationId = 0;
              toast(toasts, "Undocked.", 2.0);

              // Start any pending escort missions from this station by spawning their convoy
              // contacts (plus a small police escort wing). The mission tracks the convoy via
              // `targetNpcId`.
              for (auto& m : missions) {
                if (m.type != sim::MissionType::Escort) continue;
                if (m.failed || m.completed) continue;
                if (m.leg != 0) continue; // already started
                if (!currentSystem) continue;
                if (m.fromSystem != currentSystem->stub.id || m.fromStation != undockedFromStationId) continue;

                // Ensure we have a unique convoy id even for older saves.
                if (m.targetNpcId == 0) {
                  m.targetNpcId = (allocWorldId() | 0x4000000000000000ull);
                }

                // If the convoy already exists in-world (e.g. after a load), don't duplicate it.
                bool convoyExists = false;
                for (const auto& c : contacts) {
                  if (c.alive && c.id == m.targetNpcId) {
                    convoyExists = true;
                    break;
                  }
                }
                if (convoyExists) {
                  m.leg = 1;
                  continue;
                }

                // Map station ids to indices in the current system.
                auto findStationIndexById = [&](sim::StationId id) -> std::optional<std::size_t> {
                  for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
                    if (currentSystem->stations[i].id == id) return i;
                  }
                  return std::nullopt;
                };

                const auto fromIdx = findStationIndexById(m.fromStation);
                const auto toIdx = findStationIndexById(m.toStation);
                if (!fromIdx || !toIdx) {
                  toast(toasts, "Escort mission failed: station missing in this system.", 3.0);
                  m.failed = true;
                  continue;
                }

                const auto& originSt = currentSystem->stations[*fromIdx];
                const math::Vec3d originPosKm = stationPosKm(originSt, timeDays);
                const math::Vec3d originVelKmS = stationVelKmS(originSt, timeDays);

                // Spawn convoy (trader) contact.
                Contact convoy{};
                convoy.alive = true;
                convoy.id = m.targetNpcId;
                convoy.role = ContactRole::Trader;
                convoy.factionId = m.factionId;
                convoy.escortConvoy = true;
                convoy.groupId = convoy.id;
                convoy.leaderId = 0;

                if (m.units > 0.0) {
                  convoy.name = std::string("Convoy - ") + econ::commodityDef(m.commodity).name;
                } else {
                  convoy.name = "Convoy";
                }

                // Loadout: hauler with stronger shields. Match legacy convoy durability while
                // switching to the unified loadout/distributor model.
                configureContactLoadout(convoy,
                                       ShipHullClass::Hauler,
                                       /*thrMk=*/1,
                                       /*shieldMk=*/2,
                                       /*distMk=*/1,
                                       /*weapon=*/WeaponType::MiningLaser,
                                       /*hullMul=*/0.90,
                                       /*shieldMul=*/0.90,
                                       /*regenMul=*/1.0,
                                       /*accelMul=*/0.72,
                                       /*aiSkill=*/0.50);

                convoy.ship.setMassKg(20000.0);
                convoy.ship.setPositionKm(originPosKm + randUnit() * (originSt.radiusKm + originSt.commsRangeKm + 3500.0));
                convoy.ship.setVelocityKmS(originVelKmS);

                // Mission haul: destination is fixed and cargo is already reserved at acceptance.
                convoy.tradeDestStationIndex = (int)*toIdx;
                convoy.tradeCommodity = m.commodity;
                convoy.tradeUnits = m.units;
                const double cargoMassKg = std::max(0.0, m.units) * econ::commodityDef(m.commodity).massKg;
                convoy.tradeCapacityKg = std::max(240.0, cargoMassKg * 1.15);
                convoy.tradeCooldownUntilDays = timeDays + (20.0 / 86400.0);
                convoy.tradeSupercruiseSpeedKmS = 0.0;

                if (m.units > 0.0) {
                  auto& originEcon = universe.stationEconomy(originSt, timeDays);
                  const econ::MarketQuote q = econ::quote(originEcon, originSt.economyModel, m.commodity);
                  convoy.tradeCargoValueCr = q.mid * m.units;
                  convoy.cargoValueCr = convoy.tradeCargoValueCr;
                }

                contacts.push_back(std::move(convoy));

                // Spawn a small police escort wing that follows the convoy.
                const int escortCount = 2;
                for (int e = 0; e < escortCount; ++e) {
                  Contact escort{};
                  escort.alive = true;
                  escort.id = allocWorldId();
                  escort.role = ContactRole::Police;
                  escort.factionId = m.factionId;
                  escort.name = "Escort";
                  escort.followId = m.targetNpcId;
                  escort.groupId = m.targetNpcId;
                  escort.leaderId = m.targetNpcId;
                  configureContactLoadout(escort,
                                         ShipHullClass::Fighter,
                                         /*thrMk=*/1,
                                         /*shieldMk=*/1,
                                         /*distMk=*/1,
                                         /*weapon=*/WeaponType::PulseLaser,
                                         /*hullMul=*/1.10,
                                         /*shieldMul=*/0.80,
                                         /*regenMul=*/1.0,
                                         /*accelMul=*/0.75,
                                         /*aiSkill=*/0.65);
                  escort.ship.setMassKg(13000.0);
                  escort.ship.setPositionKm(originPosKm + randUnit() * (originSt.radiusKm + originSt.commsRangeKm + 4200.0));
                  escort.ship.setVelocityKmS(originVelKmS);
                  contacts.push_back(std::move(escort));
                }

                m.scanned = false;
                m.leg = 1;

                // Initialize / refresh runtime tracking.
                EscortRuntime* rt = nullptr;
                for (auto& r : escortRuntime) {
                  if (r.missionId == m.id) { rt = &r; break; }
                }
                if (!rt) {
                  escortRuntime.push_back(EscortRuntime{});
                  rt = &escortRuntime.back();
                }
                rt->missionId = m.id;
                rt->convoyId = m.targetNpcId;
                rt->pirateGroupId = (m.targetNpcId ^ 0x5a5a5a5a5a5a5a5aull);
                rt->tooFarSec = 0.0;
                rt->ambushSpawned = false;
                rt->nextAmbushDays = timeDays + (rng.range(70.0, 140.0) / 86400.0);

                toast(toasts, "Convoy launched. Stay close and escort to your destination.", 3.2);
              }
            } else {
              docked = false;
              dockedStationId = 0;
            }
          } else {
            // Dock attempt: must be inside slot tunnel, aligned, and have clearance.
            if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
              const auto& st = currentSystem->stations[target.index];
              auto& cs = clearances[st.id];
              const bool clearanceValid = sim::dockingClearanceValid(cs, timeDays);

              const math::Vec3d stPos = stationPosKm(st, timeDays);
              const math::Quatd stQ = stationOrient(st, stPos, timeDays);
              const math::Vec3d stV = stationVelKmS(st, timeDays);

              const math::Vec3d relWorldKm = ship.positionKm() - stPos;
              const math::Vec3d relLocalKm = stQ.conjugate().rotate(relWorldKm);
              const math::Vec3d relVelWorld = ship.velocityKmS() - stV;
              const math::Vec3d relVelLocal = stQ.conjugate().rotate(relVelWorld);

              const math::Vec3d fwdLocal = stQ.conjugate().rotate(ship.forward());
              const math::Vec3d upLocal = stQ.conjugate().rotate(ship.up());

              if (!clearanceValid) {
                toast(toasts, "No valid clearance. Press L to request.", 2.5);
              } else if (!sim::dockingSlotConditions(st, relLocalKm, relVelLocal, fwdLocal, upLocal, clearanceValid)) {
                toast(toasts, "Docking failed: align and enter the slot under speed limit.", 2.8);
              } else {
                // Docked: lock ship to a point inside hangar.
                docked = true;
                dockedStationId = st.id;
                selectedStationIndex = (int)target.index;

                ship.setVelocityKmS(stV);
                ship.setAngularVelocityRadS({0,0,0});

                toast(toasts, "Docked at " + st.name, 2.5);
              }
            } else {
              toast(toasts, "Target a station (T) before docking.", 2.5);
            }
          }
        }
      }

      if (event.type == SDL_MOUSEBUTTONDOWN &&
          (event.button.button == SDL_BUTTON_LEFT || event.button.button == SDL_BUTTON_RIGHT)) {
        if (!io.WantCaptureMouse) {
          const bool primary = (event.button.button == SDL_BUTTON_LEFT);

          const WeaponType wType = primary ? weaponPrimary : weaponSecondary;
          const WeaponDef& w = weaponDef(wType);

          double& cd = primary ? weaponPrimaryCooldown : weaponSecondaryCooldown;

          if (cd <= 0.0 && !docked && fsdState == FsdState::Idle && supercruiseState == SupercruiseState::Idle) {
            const double capCost = sim::weaponCapacitorCost(w);
            if (distributorState.wep + 1e-9 < capCost) {
              // Not enough weapon capacitor to fire.
              continue;
            }

            // Contact aim cone (soft aim assist).
            // - Pulse lasers are forgiving.
            // - Missiles have a wide lock cone, but still require targets in the forward hemisphere.
            double contactCone = 0.995;
            if (wType == WeaponType::PulseLaser) contactCone = 0.993;
            if (w.guided) contactCone = 0.985;

            // Build target list on demand so combat logic stays reusable.
            // (Beams need asteroids, missiles need lockable ships/player.)
            std::vector<sim::SphereTarget> fireTargets;
            if (w.beam) {
              fireTargets.reserve(contacts.size() + asteroids.size());

              for (std::size_t i = 0; i < contacts.size(); ++i) {
                const auto& c = contacts[i];
                if (!c.alive) continue;
                sim::SphereTarget t{};
                t.kind = sim::CombatTargetKind::Ship;
                t.index = i;
                t.id = c.id;
                t.centerKm = c.ship.positionKm();
                t.velKmS = c.ship.velocityKmS();
                t.radiusKm = 900.0;
                t.minAimCos = contactCone;
                fireTargets.push_back(t);
              }

              for (std::size_t i = 0; i < asteroids.size(); ++i) {
                const auto& a = asteroids[i];
                sim::SphereTarget t{};
                t.kind = sim::CombatTargetKind::Asteroid;
                t.index = i;
                t.id = a.id;
                t.centerKm = a.posKm;
                t.velKmS = {0,0,0};
                t.radiusKm = a.radiusKm;
                t.minAimCos = -1.0;
                fireTargets.push_back(t);
              }
            } else if (w.guided) {
              // Missiles: prefer explicitly targeted contacts; otherwise lock what's in front.
              if (target.kind == TargetKind::Contact && target.index < contacts.size() && contacts[target.index].alive) {
                const auto& c = contacts[target.index];
                sim::SphereTarget t{};
                t.kind = sim::CombatTargetKind::Ship;
                t.index = target.index;
                t.id = c.id;
                t.centerKm = c.ship.positionKm();
                t.velKmS = c.ship.velocityKmS();
                t.radiusKm = 900.0;
                t.minAimCos = 0.0; // wide lock when explicitly targeted
                fireTargets.push_back(t);
              } else {
                fireTargets.reserve(contacts.size());
                for (std::size_t i = 0; i < contacts.size(); ++i) {
                  const auto& c = contacts[i];
                  if (!c.alive) continue;
                  sim::SphereTarget t{};
                  t.kind = sim::CombatTargetKind::Ship;
                  t.index = i;
                  t.id = c.id;
                  t.centerKm = c.ship.positionKm();
                  t.velKmS = c.ship.velocityKmS();
                  t.radiusKm = 900.0;
                  t.minAimCos = contactCone;
                  fireTargets.push_back(t);
                }
              }
            }

            const sim::FireResult fr = sim::tryFireWeapon(
              ship,
              wType,
              cd,
              distributorMk,
              /*shooterId*/0,
              /*fromPlayer*/true,
              fireTargets.empty() ? nullptr : fireTargets.data(),
              fireTargets.size()
            );

            if (fr.fired) {
              distributorState.wep = std::max(0.0, distributorState.wep - capCost);
              // NOTE: cooldown is in *sim* seconds, so it naturally scales with timeScale.
              cd = fr.newCooldownSimSec;

              // Track last fired weapon for HUD lead indicator.
              hudLastFiredPrimary = primary;

              heatImpulse += fr.heatDelta;

              if (fr.hasBeam) {
                beams.push_back({toRenderU(fr.beam.aKm), toRenderU(fr.beam.bKm), fr.beam.r, fr.beam.g, fr.beam.b, 0.10});
              }

              if (fr.hasProjectile) {
                projectiles.push_back(fr.projectile);
              }

              if (fr.hasMissile) {
                missiles.push_back(fr.missile);
              }

              // Resolve beam hit effects (mining/damage). Projectile impacts are handled in the tick.
              if (w.beam && fr.hit) {
                if (fr.hitKind == sim::CombatTargetKind::Asteroid) {
                  if (wType == WeaponType::MiningLaser && fr.hitIndex < asteroids.size()) {
                    auto& a = asteroids[fr.hitIndex];
                    if (a.remainingUnits > 1e-6) {
                      const double chunk = std::min(a.remainingUnits, rng.range(6.0, 14.0));
                      a.remainingUnits -= chunk;

                      // Persist depletion for deterministic asteroid nodes so mining can't be reset
                      // by leaving/re-entering the system.
                      if (a.id & kDeterministicWorldIdBit) {
                        asteroidRemainingById[a.id] = a.remainingUnits;
                      }

                      spawnCargoPod(a.yield, chunk, fr.hitPointKm, math::Vec3d{0,0,0}, 0.45);
                      toast(toasts,
                            std::string("Mined ") + econ::commodityDef(a.yield).name + " x" + std::to_string((int)std::round(chunk)),
                            2.0);
                    } else {
                      toast(toasts, "Asteroid depleted.", 1.5);
                    }
                  }
                } else if (fr.hitKind == sim::CombatTargetKind::Ship) {
                  if (fr.hitIndex < contacts.size()) {
                    playerDamageContact((int)fr.hitIndex, w.dmg);
                  }
                }
              }
            }
          }
        }
      }

    } // while (SDL_PollEvent)

    // Input (6DOF)
    sim::ShipInput input{};
    const Uint8* keys = SDL_GetKeyboardState(nullptr);

    const bool captureKeys = io.WantCaptureKeyboard;
    if (!captureKeys && !docked) {
      input.thrustLocal.z += game::axisValue(controls.axes.thrustForward, keys);
      input.thrustLocal.x += game::axisValue(controls.axes.thrustRight, keys);
      input.thrustLocal.y += game::axisValue(controls.axes.thrustUp, keys);

      input.torqueLocal.x += game::axisValue(controls.axes.pitch, keys);
      input.torqueLocal.y += game::axisValue(controls.axes.yaw, keys);
      input.torqueLocal.z += game::axisValue(controls.axes.roll, keys);

      // Mouse steering (relative mode). Uses local pitch/yaw.
      if (mouseSteer && !io.WantCaptureMouse) {
        int dx = 0, dy = 0;
        SDL_GetRelativeMouseState(&dx, &dy);
        const double sx = (double)mouseSensitivity;
        const double yaw = (double)dx * sx;
        const double pitch = (double)(mouseInvertY ? dy : -dy) * sx;
        input.torqueLocal.y += yaw;
        input.torqueLocal.x += pitch;
      }

      input.boost = game::holdDown(controls.holds.boost, keys);
      input.brake = game::holdDown(controls.holds.brake, keys);

      static bool dampers = true;
      if (game::holdDown(controls.holds.dampersEnable, keys)) dampers = true;
      if (game::holdDown(controls.holds.dampersDisable, keys)) dampers = false;
      input.dampers = dampers;
    }

    // Docking computer: auto approach + auto-clearance + auto-dock through the slot.
    if (autopilot && !docked && fsdState == FsdState::Idle && supercruiseState == SupercruiseState::Idle &&
        currentSystem && target.kind == TargetKind::Station && target.index < currentSystem->stations.size() && !captureKeys) {

      // Optional: disengage if the player provides a clear manual input signal.
      if (dockingComputerDisengageOnManualInput) {
        auto magMax = [](double a, double b, double c) {
          return std::max(std::abs(a), std::max(std::abs(b), std::abs(c)));
        };
        const double thrustMag = magMax(input.thrustLocal.x, input.thrustLocal.y, input.thrustLocal.z);
        const double torqueMag = magMax(input.torqueLocal.x, input.torqueLocal.y, input.torqueLocal.z);
        if (thrustMag > dockingComputerManualDeadzone || torqueMag > dockingComputerManualDeadzone || input.boost || input.brake) {
          autopilot = false;
          dockingComputer.reset();
          toast(toasts, "Docking computer disengaged (manual override).", 1.6);
        }
      }

      if (autopilot) {
        const auto& st = currentSystem->stations[target.index];

        const math::Vec3d stPos = stationPosKm(st, timeDays);
        const math::Quatd stQ = stationOrient(st, stPos, timeDays);
        const math::Vec3d stV = stationVelKmS(st, timeDays);

        auto& cs = clearances[st.id];
        const bool clearanceValid = sim::dockingClearanceValid(cs, timeDays);
        const double distKm = (ship.positionKm() - stPos).length();

        // If we don't have clearance yet, try to request it automatically when in comms range.
        if (!clearanceValid) {
          int nearby = 0;
          for (const auto& c : contacts) {
            if (!c.alive) continue;
            if ((c.ship.positionKm() - stPos).length() <= st.commsRangeKm) ++nearby;
          }

          const double traffic01 = sim::estimateDockingTraffic01(seed, st, timeDays, nearby);
          const auto dec = sim::autoRequestDockingClearance(seed, st, timeDays, distKm, cs, traffic01);

          if (dec.status == sim::DockingClearanceStatus::Granted) {
            toast(toasts, "Docking clearance GRANTED.", 2.2);
          } else if (dec.status == sim::DockingClearanceStatus::Denied) {
            toast(toasts, "Docking clearance DENIED. (traffic)", 2.2);
          }
        }

        const bool clearanceNow = sim::dockingClearanceValid(cs, timeDays);

        const auto dcOut = dockingComputer.update(ship, st, stPos, stQ, stV, clearanceNow);
        input.thrustLocal = dcOut.input.thrustLocal;
        input.torqueLocal = dcOut.input.torqueLocal;
        input.dampers = dcOut.input.dampers;
        input.boost = dcOut.input.boost;
        input.brake = dcOut.input.brake;

        if (dcOut.docked) {
          docked = true;
          dockedStationId = st.id;
          selectedStationIndex = (int)target.index;

          ship.setVelocityKmS(stV);
          ship.setAngularVelocityRadS({0,0,0});

          autopilot = false;
          dockingComputer.reset();

          toast(toasts, "Docked at " + st.name + " (docking computer)", 2.5);
        }
      }
    }

    // Supercruise: fast in-system travel assist to current target.
    // Also: update the "local reference frame" velocity that Flight Assist uses for dampers.
    if (localFrameEnabled && currentSystem) {
      math::Vec3d desiredVel{0,0,0};
      double best = 1e99;

      // Prefer stations when nearby (keeps you "riding along" with their orbital motion).
      for (const auto& st : currentSystem->stations) {
        const math::Vec3d stPos = stationPosKm(st, timeDays);
        const double d = (stPos - ship.positionKm()).length();
        const double influence = std::max(st.commsRangeKm * 1.5, st.approachLengthKm * 2.0 + st.radiusKm * 12.0);
        if (d < influence && d < best) {
          best = d;
          desiredVel = stationVelKmS(st, timeDays);
        }
      }

      // Fall back to planets if no station "wins".
      for (const auto& p : currentSystem->planets) {
        const math::Vec3d pPos = planetPosKm(p, timeDays);
        const double d = (pPos - ship.positionKm()).length();
        const double rKm = p.radiusEarth * 6371.0;
        const double influence = std::max(250000.0, rKm * 80.0);
        if (d < influence && d < best) {
          best = d;
          desiredVel = planetVelKmS(p, timeDays);
        }
      }

      const double tau = std::max(0.05, localFrameBlendTauSec);
      const double alpha = 1.0 - std::exp(-dtReal / tau);
      localFrameVelKmS = localFrameVelKmS + (desiredVel - localFrameVelKmS) * alpha;
    } else {
      localFrameVelKmS = {0,0,0};
    }
    ship.setDampingFrameVelocityKmS(localFrameVelKmS);

    // Supercruise state timers (real time)
    if (supercruiseState == SupercruiseState::Charging) {
      supercruiseChargeRemainingSec = std::max(0.0, supercruiseChargeRemainingSec - dtReal);
      if (supercruiseChargeRemainingSec <= 0.0) {
        supercruiseState = SupercruiseState::Active;
        supercruiseDropRequested = false;
        // Small grace window before interdictions can trigger.
        interdiction = sim::InterdictionState{};
        interdictionSubmitRequested = false;
        interdictionPirateName.clear();
        interdictionPirateStrength = 1.0;
        interdictionCooldownUntilDays = std::max(interdictionCooldownUntilDays, timeDays + (12.0 / 86400.0));
        toast(toasts, "Supercruise engaged.", 1.6);
      }
    } else if (supercruiseState == SupercruiseState::Cooldown) {
      supercruiseCooldownRemainingSec = std::max(0.0, supercruiseCooldownRemainingSec - dtReal);
      if (supercruiseCooldownRemainingSec <= 0.0) {
        supercruiseState = SupercruiseState::Idle;
      }
    }

    // Reset per-frame HUD values
    supercruiseSafeDropReady = false;
    supercruiseTtaSec = 0.0;
    supercruiseDistKm = 0.0;
    supercruiseClosingKmS = 0.0;

    auto endSupercruise = [&](bool emergency, const math::Vec3d& destVelKmS, const math::Vec3d& dirToDest, const char* msg) {
      supercruiseState = SupercruiseState::Cooldown;
      supercruiseCooldownRemainingSec = emergency ? kSupercruiseEmergencyCooldownSec : kSupercruiseCooldownSec;
      supercruiseCooldownTotalSec = supercruiseCooldownRemainingSec;
      supercruiseDropRequested = false;

      // Clear any ongoing interdiction state.
      interdiction = sim::InterdictionState{};
      interdictionSubmitRequested = false;
      interdictionPirateName.clear();
      interdictionPirateStrength = 1.0;

      const double approachSpeed = emergency ? 0.45 : 0.16;
      ship.setVelocityKmS(destVelKmS + dirToDest * approachSpeed);

	      // Detach any supercruise "shadow" contacts so they remain nearby in normal space.
	      for (auto& c : contacts) {
	        if (!c.alive || !c.supercruiseShadow) continue;
	        c.supercruiseShadow = false;
	        c.shadowOffsetLocalKm = {0,0,0};
	        c.ship.setVelocityKmS(ship.velocityKmS() + math::Vec3d{rng.range(-0.08, 0.08), rng.range(-0.08, 0.08), rng.range(-0.08, 0.08)});
	        c.ship.setAngularVelocityRadS({0,0,0});
	      }

      if (emergency) {
        heatImpulse += 25.0;
        playerShield = std::max(0.0, playerShield - playerShieldMax * 0.12);
        playerHull = std::max(0.0, playerHull - playerHullMax * 0.04);
        ship.setAngularVelocityRadS({rng.range(-0.6, 0.6), rng.range(-0.6, 0.6), rng.range(-0.4, 0.4)});
      } else {
        ship.setAngularVelocityRadS({0,0,0});
      }

      toast(toasts, msg, 2.2);
    };

    if (supercruiseState == SupercruiseState::Active && !docked && !captureKeys && fsdState == FsdState::Idle) {
	      // Ensure active escort convoys (and their police escorts) remain nearby in supercruise.
	      if (currentSystem) {
	        std::vector<core::u64> activeConvoys;
	        activeConvoys.reserve(4);
	        for (const auto& m : missions) {
	          if (m.type != sim::MissionType::Escort) continue;
	          if (m.failed || m.completed) continue;
	          if (m.leg < 1) continue;  // Not started yet.
	          if (m.fromSystem != currentSystem->stub.id) continue;
	          if (m.targetNpcId != 0) activeConvoys.push_back(m.targetNpcId);
	        }

	        if (!activeConvoys.empty()) {
	          auto isActiveConvoyId = [&](core::u64 id) -> bool {
	            for (core::u64 cid : activeConvoys) {
	              if (cid == id) return true;
	            }
	            return false;
	          };

	          for (auto& c : contacts) {
	            if (!c.alive) continue;
	            const bool isConvoy = c.escortConvoy && isActiveConvoyId(c.id);
	            const bool isEscort = (c.role == ContactRole::Police) && (c.followId != 0) && isActiveConvoyId(c.followId);
	            if (!isConvoy && !isEscort) continue;

	            c.supercruiseShadow = true;
	            if (c.shadowOffsetLocalKm.lengthSq() < 1.0) {
	              const double lateral = rng.range(-22000.0, 22000.0);
	              const double vertical = rng.range(-9000.0, 9000.0);
	              const double back = isConvoy ? -rng.range(110000.0, 160000.0)
	                                         : -rng.range(140000.0, 200000.0);
	              c.shadowOffsetLocalKm = {lateral, vertical, back};
	            }
	          }
	        }
	      }

	      // Keep supercruise "shadow" contacts tethered to the player.
	      for (auto& c : contacts) {
	        if (!c.alive || !c.supercruiseShadow) continue;
	        const math::Vec3d offWorld = ship.orientation().rotate(c.shadowOffsetLocalKm);
	        c.ship.setPositionKm(ship.positionKm() + offWorld);
	        c.ship.setVelocityKmS(ship.velocityKmS());
	        c.ship.setAngularVelocityRadS({0,0,0});
	      }

	      // Spawn occasional supercruise encounters (pirate shadows on valuable cargo routes;
	      // police shadows if you're wanted in this jurisdiction).
	      if (currentSystem && timeDays >= nextSupercruiseShadowSpawnDays) {
	        nextSupercruiseShadowSpawnDays = timeDays + (rng.range(55.0, 95.0) / 86400.0);
	        const core::u32 jurisdiction = currentSystem->stub.factionId;
	        const double cargoV = cargoValueEstimateCr();
	        const bool wantedHere = (jurisdiction != 0) && (getBounty(jurisdiction) > 0.0);

	        int shadowPirates = 0;
	        int shadowPolice = 0;
	        for (const auto& c : contacts) {
	          if (!c.alive || !c.supercruiseShadow) continue;
	          if (c.role == ContactRole::Pirate) shadowPirates++;
	          if (c.role == ContactRole::Police) shadowPolice++;
	        }

	        // Pirates are more likely to show up the more valuable your cargo is.
	        if (shadowPirates < 1 && cargoV > 2500.0 && timeDays >= interdictionCooldownUntilDays) {
	          const double risk = std::clamp(cargoV / 9000.0, 0.0, 1.5);
	          const double chance = std::clamp(0.20 + 0.35 * risk, 0.0, 0.75);
	          if (rng.nextUnit() < chance) {
	            Contact p{};
	            p.id = allocWorldId();
	            p.role = ContactRole::Pirate;
	            p.name = "Pirate " + std::to_string((int)contacts.size() + 1);
	            configureContactLoadout(p,
	                                   ShipHullClass::Fighter,
	                                   /*thrMk=*/1,
	                                   /*shieldMk=*/1,
	                                   /*distMk=*/1,
	                                   /*weapon=*/WeaponType::Cannon,
	                                   /*hullMul=*/0.95,
	                                   /*shieldMul=*/0.67,
	                                   /*regenMul=*/0.84,
	                                   /*accelMul=*/0.70,
	                                   /*aiSkill=*/0.55);
	            p.hostileToPlayer = true;
	            p.supercruiseShadow = true;
	            // Keep them behind you so they can pressure an interdiction.
	            p.shadowOffsetLocalKm = {rng.range(-25000.0, 25000.0), rng.range(-12000.0, 12000.0), -rng.range(165000.0, 220000.0)};
	            p.ship.setPositionKm(ship.positionKm() + ship.orientation().rotate(p.shadowOffsetLocalKm));
	            p.ship.setVelocityKmS(ship.velocityKmS());
	            p.ship.setOrientation(ship.orientation());
	            p.ship.setAngularVelocityRadS({0,0,0});
	            contacts.push_back(std::move(p));
	            toast(toasts, "Contact: pirate wake signature detected!", 2.8);
	          }
	        }

	        // If you're wanted, law enforcement can follow you through supercruise and be waiting on drop.
	        if (wantedHere && shadowPolice < 1) {
	          if (rng.nextUnit() < 0.55) {
	            Contact pc{};
	            pc.id = allocWorldId();
	            pc.role = ContactRole::Police;
	            pc.factionId = jurisdiction;
	            pc.name = "Police";
	            configureContactLoadout(pc,
	                                   ShipHullClass::Fighter,
	                                   /*thrMk=*/1,
	                                   /*shieldMk=*/1,
	                                   /*distMk=*/1,
	                                   /*weapon=*/WeaponType::PulseLaser,
	                                   /*hullMul=*/0.95,
	                                   /*shieldMul=*/0.75,
	                                   /*regenMul=*/0.90,
	                                   /*accelMul=*/0.75,
	                                   /*aiSkill=*/0.65);
	            pc.supercruiseShadow = true;
	            pc.shadowOffsetLocalKm = {rng.range(-22000.0, 22000.0), rng.range(-10000.0, 10000.0), -rng.range(140000.0, 200000.0)};
	            pc.ship.setPositionKm(ship.positionKm() + ship.orientation().rotate(pc.shadowOffsetLocalKm));
	            pc.ship.setVelocityKmS(ship.velocityKmS());
	            pc.ship.setOrientation(ship.orientation());
	            pc.ship.setAngularVelocityRadS({0,0,0});
	            contacts.push_back(std::move(pc));
	          }
	        }
	      }

	      // Determine destination (station / planet / signal source)
      bool hasDest = false;
      math::Vec3d destPosKm{0,0,0};
      math::Vec3d destVelKmS{0,0,0};
      double dropKm = 0.0;

      if (currentSystem && target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
        const auto& st = currentSystem->stations[target.index];
        destPosKm = stationPosKm(st, timeDays);
        destVelKmS = stationVelKmS(st, timeDays);
        dropKm = std::max(15000.0, st.radiusKm * 3.0);
        hasDest = true;
      } else if (currentSystem && target.kind == TargetKind::Planet && target.index < currentSystem->planets.size()) {
        const auto& p = currentSystem->planets[target.index];
        destPosKm = planetPosKm(p, timeDays);
        destVelKmS = planetVelKmS(p, timeDays);
        const double rKm = p.radiusEarth * 6371.0;
        dropKm = std::max(60000.0, rKm * 12.0);
        hasDest = true;
	      } else if (target.kind == TargetKind::Signal && target.index < signals.size()) {
	        const auto& s = signals[target.index];
	        destPosKm = s.posKm;
	        destVelKmS = {0, 0, 0};
	        dropKm = 70000.0;
	        hasDest = true;
      }

      math::Vec3d dirToDest = ship.forward().normalized();
      if (hasDest) {
        const math::Vec3d rel = destPosKm - ship.positionKm();
        const double dist = rel.length();
        if (dist > 1e-6) dirToDest = rel / dist;
      }

	      // Spawn intermittent signal sources while travelling in supercruise.
	      if (currentSystem && hasDest && timeDays >= nextSignalSpawnDays && signals.size() < 10) {
	        nextSignalSpawnDays = timeDays + (rng.range(55.0, 120.0) / 86400.0);

	        const double r = rng.nextUnit();
	        SignalType t = SignalType::Distress;
	        if (r < 0.55) t = SignalType::Distress;
	        else if (r < 0.82) t = SignalType::Derelict;
	        else t = SignalType::Resource;

	        // Place the signal somewhere ahead along the line to the destination, with some lateral offset.
	        const double aheadKm = rng.range(200000.0, 520000.0);
	        math::Vec3d side = math::cross(dirToDest, randUnit());
	        if (side.length() < 1e-6) side = math::cross(dirToDest, math::Vec3d{0,1,0});
	        side = side.normalized();
	        const double sideKm = rng.range(-140000.0, 140000.0);
	        math::Vec3d pos = ship.positionKm() + dirToDest * aheadKm + side * sideKm;

	        SignalSource s;
	        s.id = allocWorldId();
	        s.type = t;
	        s.posKm = pos;
	        s.expireDay = timeDays + rng.range(0.010, 0.030); // ~14-43 minutes (sim-days)
	        s.resolved = false;
	        s.fieldSpawned = false;

	        // Generate a stable plan for distress calls so scans can reveal useful info
	        // and the on-drop content doesn't change if the player loops back.
	        if (t == SignalType::Distress) {
	          s.hasDistressPlan = true;
	          s.distress = sim::planDistressEncounter(universe.seed(), currentSystem->stub.id, s.id, timeDays, jurisdiction);
	        }
	        signals.push_back(s);

	        toast(toasts, std::string("Signal detected: ") + signalTypeName(t), 3.0);
	      }

	      // Interdiction (supercruise tether minigame)
      double nearestPirateKm = 1e99;
      core::u64 nearestPirateId = 0;
      std::string nearestPirateName;
      double nearestPirateStrength = 1.0;

      auto pirateStrengthFrom = [&](const Contact& c) -> double {
        // Rough difficulty scaling for the interdiction tether.
        // Keep it stable and inexpensive; the detailed combat sim starts after drop.
        double s = 0.85 + 0.70 * std::clamp(c.aiSkill, 0.0, 1.0);
        if (c.hullClass == ShipHullClass::Fighter) s += 0.20;
        s += 0.05 * std::max(0, c.thrusterMk - 1);
        s += 0.03 * std::max(0, c.shieldMk - 1);
        s += 0.02 * std::max(0, c.distributorMk - 1);
        return std::clamp(s, 0.60, 1.80);
      };

      for (const auto& c : contacts) {
        if (!c.alive || c.role != ContactRole::Pirate) continue;
        const double d = (c.ship.positionKm() - ship.positionKm()).length();
        if (d < nearestPirateKm) {
          nearestPirateKm = d;
          nearestPirateId = c.id;
          nearestPirateName = c.name;
          nearestPirateStrength = pirateStrengthFrom(c);
        }
      }

      const bool pirateThreat = (nearestPirateKm < 240000.0);

      // Compute closeness for the active interdictor (if any), otherwise nearest pirate.
      double interdictorKm = nearestPirateKm;
      if (sim::interdictionInProgress(interdiction) && interdiction.pirateId != 0) {
        for (const auto& c : contacts) {
          if (!c.alive || c.id != interdiction.pirateId) continue;
          interdictorKm = (c.ship.positionKm() - ship.positionKm()).length();
          break;
        }
      }
      const double interdictorCloseness = sim::interdictionCloseness01(interdictorKm, 240000.0);

      // "Submit" (H) while the interdiction is ongoing: safe-ish drop with short cooldown.
      if (interdictionSubmitRequested && sim::interdictionInProgress(interdiction)) {
        interdictionSubmitRequested = false;

        // End the minigame (submitted) and drop without emergency effects.
        (void)sim::stepInterdiction(interdiction, 0.0, ship.forward().normalized(), interdictorCloseness,
                                    interdictionPirateStrength, /*submit=*/true, interdictionParams);

        interdictionCooldownUntilDays = timeDays + (90.0 / 86400.0);
        endSupercruise(false, destVelKmS, dirToDest, "Submitted: dropped from supercruise.");
      } else if (!hasDest) {
        endSupercruise(true, localFrameVelKmS, dirToDest, "Supercruise lost target. Emergency drop.");
      } else {
        const double cargoV = cargoValueEstimateCr();

        // Start interdiction warning (random-ish cadence based on proximity + cargo value)
        if (!sim::interdictionInProgress(interdiction) && pirateThreat && timeDays >= interdictionCooldownUntilDays) {
          if (sim::rollInterdictionStart(rng, dtReal, interdictorCloseness, cargoV, nearestPirateStrength, interdictionTriggerParams)) {
            interdiction = sim::beginInterdiction(universe.seed(), nearestPirateId, timeDays,
                                                 ship.forward().normalized(), ship.up().normalized(),
                                                 interdictorCloseness, nearestPirateStrength,
                                                 interdictionParams);
            interdictionPirateName = nearestPirateName;
            interdictionPirateStrength = nearestPirateStrength;

            supercruiseDropRequested = false; // disable manual drop while interdicted
            toast(toasts, "INTERDICTION WARNING! Keep the escape vector centered (H to submit).", 2.6);
          }
        }

        // Step the minigame.
        if (sim::interdictionInProgress(interdiction)) {
          const auto out = sim::stepInterdiction(interdiction, dtReal, ship.forward().normalized(),
                                                 interdictorCloseness, interdictionPirateStrength,
                                                 /*submit=*/false, interdictionParams);

          if (out.beganActive) {
            toast(toasts, "Interdiction engaged! Follow the escape vector.", 2.2);
          }

          if (out.evaded) {
            interdictionCooldownUntilDays = timeDays + (150.0 / 86400.0);
            interdictionPirateName.clear();
            interdictionPirateStrength = 1.0;
            toast(toasts, "Interdiction evaded.", 2.0);
          } else if (out.failed) {
            interdictionCooldownUntilDays = timeDays + (180.0 / 86400.0);
            endSupercruise(true, destVelKmS, dirToDest, "Interdicted! Emergency drop.");
          }
        }

        // Supercruise guidance (safe-drop heuristic + braking-distance-aware speed profile).
        // This extracts the math/control logic from the game loop so it's testable and reusable.
        const bool interdicted = sim::interdictionInProgress(interdiction);

        sim::SupercruiseParams scParams{};
        scParams.safeTtaSec = kSupercruiseSafeTtaSec;
        scParams.safeWindowSlackSec = 2.0;
        scParams.minClosingKmS = 0.05;
        scParams.maxSpeedKmS = supercruiseMaxSpeedKmS;
        scParams.accelCapKmS2 = 6.0;
        scParams.angularCapRadS2 = 1.2;
        scParams.useBrakingDistanceLimit = true;

        const auto sc = sim::guideSupercruise(ship, destPosKm, destVelKmS, dropKm,
                                              /*navAssistEnabled=*/supercruiseAssist,
                                              /*dropRequested=*/supercruiseDropRequested,
                                              /*interdicted=*/interdicted,
                                              scParams);

        supercruiseSafeDropReady = sc.hud.safeDropReady;
        supercruiseTtaSec = sc.hud.ttaSec;
        supercruiseDistKm = sc.hud.distKm;
        supercruiseClosingKmS = sc.hud.closingKmS;

        if (sc.dropNow) {
          if (sc.emergencyDrop) endSupercruise(true, destVelKmS, sc.dirToDest, "EMERGENCY DROP!");
          else endSupercruise(false, destVelKmS, sc.dirToDest, "Supercruise drop.");
        } else {
          input.thrustLocal = sc.input.thrustLocal;
          if (!interdicted) {
            input.torqueLocal = sc.input.torqueLocal;
          }
          input.dampers = sc.input.dampers;
          input.brake = sc.input.brake;
          input.boost = sc.input.boost;

          ship.setMaxLinearAccelKmS2(sc.recommendedMaxLinearAccelKmS2);
          ship.setMaxAngularAccelRadS2(sc.recommendedMaxAngularAccelRadS2);
        }
      }
    }

    if (supercruiseState == SupercruiseState::Active) {
      ship.setMaxLinearAccelKmS2(6.0);
      ship.setMaxAngularAccelRadS2(1.2);
    }

    // Restore handling caps when not in active supercruise (also applies "blue zone" turn assist).
    if (supercruiseState != SupercruiseState::Active) {
      ship.setMaxLinearAccelKmS2(playerBaseLinAccelKmS2);

      double ang = playerBaseAngAccelRadS2;
      if (blueZoneTurnAssist && input.dampers && !docked && fsdState == FsdState::Idle) {
        const double v = (ship.velocityKmS() - localFrameVelKmS).length();
        const double vTerm = std::max(0.001, ship.maxLinearAccelKmS2() / std::max(0.05, ship.dampingLinear()));
        const double x = std::clamp(v / vTerm, 0.0, 1.0);

        // Gaussian-ish bump centered around mid-speed.
        const double center = 0.55;
        const double sigma = 0.25;
        const double bump = std::exp(-std::pow((x - center) / sigma, 2.0));
        const double factor = 0.65 + 0.35 * bump;
        ang *= factor;
      }

      ship.setMaxAngularAccelRadS2(ang);
    }


    // VFX step uses real time (not sim time). When paused, freeze particles.
    const double dtVfx = paused ? 0.0 : dtReal;
    if (vfxParticlesEnabled) {
      particles.update(dtVfx);
    }

    // Sim step
    const double dtSim = dtReal * timeScale;
    if (!paused) {
      // --- FSD (system-to-system travel) ---
      bool fsdJustArrived = false;
      if (fsdState == FsdState::Charging) {
        fsdChargeRemainingSec = std::max(0.0, fsdChargeRemainingSec - dtReal);
        if (fsdChargeRemainingSec <= 0.0) {
          // Consume fuel at the moment we enter hyperspace.
          fuel = std::clamp(fuel - fsdFuelCost, 0.0, fuelMax);
          fsdState = FsdState::Jumping;
          fsdTravelTotalSec = 1.25 + fsdJumpDistanceLy * 0.08;
          fsdTravelRemainingSec = fsdTravelTotalSec;
          toast(toasts, "FSD: entering hyperspace...", 1.8);
        }
      } else if (fsdState == FsdState::Jumping) {
        fsdTravelRemainingSec = std::max(0.0, fsdTravelRemainingSec - dtReal);
        if (fsdTravelRemainingSec <= 0.0) {
          const auto& nextSystem = universe.getSystem(fsdTargetSystem);
          currentStub = nextSystem.stub;
          currentSystem = &nextSystem;

          // Clear transient state.
          docked = false;
          dockedStationId = 0;
          selectedStationIndex = 0;
          clearances.clear();
          contacts.clear();
          beams.clear();
          projectiles.clear();
          missiles.clear();

          // VFX
          particles.clear();
	          floatingCargo.clear();
	          signals.clear();
	          asteroids.clear();
	          policeDemand = PoliceDemand{};
          scanning = false;
          scanProgressSec = 0.0;
          scanLockedId = 0;
          scanLabel.clear();
          scanLockedTarget = Target{};
          target = Target{};
	          seedSystemSpaceObjects(*currentSystem);
	          nextSignalSpawnDays = timeDays + (rng.range(35.0, 70.0) / 86400.0);
	          nextSupercruiseShadowSpawnDays = timeDays + (rng.range(45.0, 90.0) / 86400.0);

          // Spawn near the first station.
          respawnNearStation(*currentSystem, 0);
          galaxySelectedSystemId = currentSystem->stub.id;

          // If we have a pending arrival target station (from missions / trade helper), auto-target it when we reach its system.
          if (pendingArrivalTargetStationId != 0) {
            for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
              if (currentSystem->stations[i].id == pendingArrivalTargetStationId) {
                target.kind = TargetKind::Station;
                target.index = i;
                selectedStationIndex = (int)i;
                toast(toasts, "Target set: " + currentSystem->stations[i].name, 2.0);
                pendingArrivalTargetStationId = 0;
                break;
              }
            }
          }

          // Cooldown (sim time)
          fsdReadyDay = timeDays + (kFsdCooldownSec / 86400.0);

          // Advance route cursor if this jump matches the plotted route.
          if (!navRoute.empty() && navRouteHop + 1 < navRoute.size() && navRoute[navRouteHop + 1] == currentSystem->stub.id) {
            navRouteHop++;
            if (navRouteHop + 1 >= navRoute.size()) {
              navAutoRun = false;
              toast(toasts, "Route complete.", 2.0);
            }
          }
// QoL: auto-select the next hop in the Galaxy UI after each jump.
if (!navRoute.empty() && navRouteHop + 1 < navRoute.size()) {
  galaxySelectedSystemId = navRoute[navRouteHop + 1];
} else {
  galaxySelectedSystemId = currentSystem->stub.id;
}


          toast(toasts, std::string("Arrived in ") + currentSystem->stub.name + ".", 2.0);

          fsdState = FsdState::Idle;
          fsdTargetSystem = 0;
          fsdJumpDistanceLy = 0.0;
          fsdFuelCost = 0.0;
          fsdTravelTotalSec = 0.0;
          fsdJustArrived = true;
        }
      }

      // Auto-run route (hands-free multi-jump). We only auto-trigger when safe.
      // QoL: if cargo/loadout changes make the next hop invalid, attempt a one-shot replot to the
      // final destination instead of immediately disabling auto-run.
      if (navAutoRun && fsdState == FsdState::Idle && timeDays >= fsdReadyDay && !docked && supercruiseState == SupercruiseState::Idle) {
        if (!navRoute.empty() && navRouteHop + 1 < navRoute.size() && !isMassLocked() && currentSystem) {

          auto targetNearestStationForRefuel = [&]() -> bool {
            if (currentSystem->stations.empty()) return false;
            const math::Vec3d p = ship.positionKm();
            double bestD = 1e30;
            std::size_t bestIdx = 0;
            for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
              const auto& st = currentSystem->stations[i];
              const math::Vec3d sp = stationPosKm(st, timeDays);
              const double d = (sp - p).length();
              if (d < bestD) { bestD = d; bestIdx = i; }
            }
            target.kind = TargetKind::Station;
            target.index = bestIdx;
            selectedStationIndex = (int)bestIdx;
            return true;
          };

          sim::SystemId nextHopId = navRoute[navRouteHop + 1];
          double hopDistLy = (universe.getSystem(nextHopId).stub.posLy - currentSystem->stub.posLy).length();

          const double rangeLy = fsdBaseRangeLy();
          if (hopDistLy > rangeLy + 1e-9) {
            const sim::SystemId destSys = navRoute.back();
            const bool wantAuto = navAutoRun;

            // Try to replot the route to the final destination. This can happen if the player
            // changed cargo/loadout after plotting, reducing jump range.
            if (destSys != 0 && destSys != currentSystem->stub.id
                && plotRouteToSystem(destSys, false)
                && !navRoute.empty() && navRouteHop + 1 < navRoute.size()) {

              // plotRouteToSystem resets navAutoRun; restore it.
              navAutoRun = wantAuto;
              nextHopId = navRoute[navRouteHop + 1];
              hopDistLy = (universe.getSystem(nextHopId).stub.posLy - currentSystem->stub.posLy).length();
              toast(toasts, "Auto-run: route replotted (range changed).", 2.4);
            } else {
              navAutoRun = false;
              toast(toasts, "Auto-run stopped: next hop is out of range (replot route).", 2.8);
            }
          }

          if (navAutoRun) {
            const double hopFuel = fsdFuelCostFor(hopDistLy);
            if (fuel + 1e-6 < hopFuel) {
              navAutoRun = false;
              if (targetNearestStationForRefuel()) {
                toast(toasts, "Auto-run stopped: refuel needed (target set).", 2.8);
              } else {
                toast(toasts, "Auto-run stopped: not enough fuel for next hop.", 2.8);
              }
            } else if (hopDistLy > fsdBaseRangeLy() + 1e-9) {
              // Safety net (should be extremely rare).
              navAutoRun = false;
              toast(toasts, "Auto-run stopped: next hop is out of range.", 2.8);
            } else {
              startFsdJumpTo(nextHopId);
            }
          }
        }
      }


      // --- Power distributor ---
      boostAppliedFrac = 0.0;
      sim::stepDistributor(distributorState, distributorCfg, distributorPips, dtSim);

      // --- Ship physics ---
      if (!docked) {
        if (fsdState != FsdState::Idle || fsdJustArrived) {
          sim::ShipInput hold{};
          hold.dampers = true;
          hold.brake = true;
          ship.step(dtSim, hold);
        } else {
          // Apply boost only for the fraction we have ENG capacitor for.
          boostAppliedFrac = input.boost ? sim::consumeBoostFraction(distributorState, distributorCfg, dtSim) : 0.0;
          input.boost = (boostAppliedFrac > 1e-6);

          const double dtBoost = dtSim * boostAppliedFrac;
          const double dtNoBoost = dtSim - dtBoost;

          if (dtBoost > 0.0) {
            sim::ShipInput inBoost = input;
            inBoost.boost = true;
            ship.step(dtBoost, inBoost);
          }
          if (dtNoBoost > 0.0) {
            sim::ShipInput inNoBoost = input;
            inNoBoost.boost = false;
            ship.step(dtNoBoost, inNoBoost);
          }
        }
      } else {
        // While docked, keep ship attached to station.
        auto it = std::find_if(currentSystem->stations.begin(), currentSystem->stations.end(),
                               [&](const sim::Station& s){ return s.id == dockedStationId; });
        if (it != currentSystem->stations.end()) {
          const auto& st = *it;
          const math::Vec3d stPos = stationPosKm(st, timeDays);
          const math::Quatd stQ = stationOrient(st, stPos, timeDays);
          const math::Vec3d stV = stationVelKmS(st, timeDays);

          const double wz = st.radiusKm * 1.10;
          const double dockZ = wz - st.slotDepthKm - st.radiusKm * 0.25;
          const math::Vec3d dockLocal{0,0, dockZ};
          ship.setPositionKm(stPos + stQ.rotate(dockLocal));
          ship.setVelocityKmS(stV);
          ship.setOrientation(stQ * math::Quatd::fromAxisAngle({0,1,0}, math::kPi)); // face outward-ish
        }
      }

      // VFX: thruster plume.
      if (vfxParticlesEnabled && vfxThrustersEnabled && !docked && fsdState == FsdState::Idle) {
        math::Vec3d accelLocal = input.thrustLocal;
        double intensity = std::clamp(accelLocal.length(), 0.0, 1.0);

        // Braking with no explicit thrust still feels like "something" should fire.
        if (input.brake && intensity < 0.12) {
          intensity = 0.35;
        }

        if (input.boost) {
          intensity = std::clamp(intensity + 0.35, 0.0, 1.0);
        }

        intensity *= (double)vfxParticleIntensity;

        if (intensity > 1e-4) {
          math::Vec3d accelWorld{0,0,1};
          if (accelLocal.lengthSq() > 1e-8) {
            accelWorld = ship.orientation().rotate(accelLocal).normalized();
          } else if (input.brake) {
            const math::Vec3d v = ship.velocityKmS() - localFrameVelKmS;
            if (v.lengthSq() > 1e-10) {
              accelWorld = (-v).normalized();
            } else {
              accelWorld = ship.forward().normalized();
            }
          } else {
            accelWorld = ship.forward().normalized();
          }

          const math::Vec3d exhaustDir = (-accelWorld).normalized();
          const math::Vec3d shipPosU = toRenderU(ship.positionKm());

          // Push the emitter slightly toward the exhaust direction so the plume starts behind the ship.
          const math::Vec3d emitterPosU = shipPosU + exhaustDir * 0.75;

          particles.emitThruster(emitterPosU, exhaustDir, std::clamp(intensity, 0.0, 1.0), dtVfx, input.boost);
        }
      }

      // --- Fuel burn (sim-time) ---
      if (!docked && fsdState == FsdState::Idle) {
        if (input.boost) {
          fuel -= dtSim * 0.0025;
        }
        if (supercruiseState == SupercruiseState::Active) {
          const double v = ship.velocityKmS().length();
          fuel -= dtSim * (0.0020 + v * 0.00000015);
        }

        if (fuel < 0.0) fuel = 0.0;
        if (fuel > fuelMax) fuel = fuelMax;

        if (fuel <= 0.0 && supercruiseState == SupercruiseState::Active) {
          supercruiseState = SupercruiseState::Cooldown;
          supercruiseCooldownRemainingSec = kSupercruiseEmergencyCooldownSec;
          supercruiseCooldownTotalSec = supercruiseCooldownRemainingSec;
          supercruiseDropRequested = false;

          heatImpulse += 25.0;
          playerShield = std::max(0.0, playerShield - 10.0);
          playerHull = std::max(0.0, playerHull - 3.0);
          ship.setVelocityKmS(localFrameVelKmS + ship.forward().normalized() * 0.45);
          ship.setAngularVelocityRadS({rng.range(-0.6, 0.6), rng.range(-0.6, 0.6), rng.range(-0.4, 0.4)});

          toast(toasts, "Fuel depleted: EMERGENCY DROP!", 2.4);
        }
      }

const bool combatSimEnabled = (fsdState == FsdState::Idle && supercruiseState == SupercruiseState::Idle);
if (combatSimEnabled) {
  // ---- Spawns / combat / traffic ----

  int alivePirates = 0;
int aliveTraders = 0;
int alivePolice = 0;
int aliveTotal = 0;
for (const auto& c : contacts) {
  if (!c.alive) continue;
  ++aliveTotal;
  if (c.role == ContactRole::Pirate) ++alivePirates;
  if (c.role == ContactRole::Trader) ++aliveTraders;
  if (c.role == ContactRole::Police) ++alivePolice;
}

  const core::u32 localFaction = currentSystem ? currentSystem->stub.factionId : 0;
  const double localRep = getRep(localFaction);
  const double localBounty = getBounty(localFaction);
  const bool playerWantedHere = (localFaction != 0) && (localBounty > 0.0);

  // Simple helper for random directions.
  auto randDir = [&]() -> math::Vec3d {
    return math::Vec3d{rng.range(-1.0,1.0), rng.range(-0.3,0.3), rng.range(-1.0,1.0)}.normalized();
  };

  auto chooseHomeStation = [&]() -> std::size_t {
    if (!currentSystem || currentSystem->stations.empty()) return 0;
    return (std::size_t)(rng.nextU64() % (core::u64)currentSystem->stations.size());
  };

  auto spawnNearStation = [&](std::size_t stIdx, double minDistMul, double maxDistMul, sim::Ship& outShip) {
    if (!currentSystem || currentSystem->stations.empty()) {
      // fallback: near player
      const math::Vec3d d = randDir();
      const double distKm = rng.range(60000.0, 130000.0);
      outShip.setPositionKm(ship.positionKm() + d * distKm);
      outShip.setVelocityKmS(ship.velocityKmS());
      outShip.setOrientation(quatFromTo({0,0,1}, (-d).normalized()));
      return;
    }
    stIdx = std::min(stIdx, currentSystem->stations.size() - 1);
    const auto& st = currentSystem->stations[stIdx];
    const math::Vec3d stPos = stationPosKm(st, timeDays);
    const math::Quatd stQ = stationOrient(st, stPos, timeDays);
    const math::Vec3d axis = stQ.rotate({0,0,1});
    const math::Vec3d tangent = stQ.rotate({1,0,0});

    const double distKm = rng.range(st.radiusKm * minDistMul, st.radiusKm * maxDistMul);
    const math::Vec3d p = stPos + axis * distKm + tangent * rng.range(-st.radiusKm*2.0, st.radiusKm*2.0);
    outShip.setPositionKm(p);
    outShip.setVelocityKmS(stationVelKmS(st, timeDays));
    outShip.setOrientation(quatFromTo({0,0,1}, (stPos - p).normalized()));
  };

  // Pick a different station index in the current system (best-effort).
  auto chooseOtherStation = [&](std::size_t avoidIdx) -> std::size_t {
    if (!currentSystem || currentSystem->stations.empty()) return 0;
    const std::size_t n = currentSystem->stations.size();
    if (n == 1) return 0;

    avoidIdx = std::min(avoidIdx, n - 1);
    // Choose an offset in [1, n-1] and wrap.
    const std::size_t off = 1u + (std::size_t)(rng.nextU32() % (core::u32)(n - 1));
    return (avoidIdx + off) % n;
  };

  // Give a trader a simple local hauling job:
  // - choose a destination station (different from source)
  // - choose a commodity and load units by reducing station inventory
  // This makes NPC traffic actually influence station inventories/prices.
  auto planTraderHaul = [&](Contact& t, std::size_t fromIdx) -> bool {
    if (!currentSystem || currentSystem->stations.size() < 2) return false;

    const std::size_t n = currentSystem->stations.size();
    fromIdx = std::min(fromIdx, n - 1);
    const std::size_t toIdx = chooseOtherStation(fromIdx);

    const auto& fromSt = currentSystem->stations[fromIdx];
    const auto& toSt   = currentSystem->stations[toIdx];

    auto& fromE = universe.stationEconomy(fromSt, timeDays);
    auto& toE   = universe.stationEconomy(toSt, timeDays);

    // Prefer profitable routes (creates believable trade lanes).
    econ::CommodityId cid = econ::CommodityId::Food;
    {
      const auto routes = econ::bestRoutes(fromE, fromSt.economyModel, toE, toSt.economyModel, 0.10, 3);
      if (!routes.empty()) {
        cid = routes.front().commodity;
      } else {
        // If no clear arbitrage, move whichever good is most "surplus" at the source.
        double bestSurplus = -1e99;
        std::size_t bestI = 0;
        for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
          const double surplus = fromE.inventory[i] - fromSt.economyModel.desiredStock[i];
          if (surplus > bestSurplus) {
            bestSurplus = surplus;
            bestI = i;
          }
        }
        cid = (econ::CommodityId)bestI;
      }
    }

    // Compute how many units we can carry.
    const double capKg = std::max(40.0, t.tradeCapacityKg);
    const double massKg = std::max(1e-6, econ::commodityDef(cid).massKg);
    const double wantUnits = std::max(1.0, std::floor(capKg / massKg));

    // Load from source inventory (clamped).
    const double loaded = econ::takeInventory(fromE, fromSt.economyModel, cid, wantUnits);
    if (loaded <= 0.0) {
      // Fallback: try a couple of random commodities to avoid "empty" traders.
      for (int attempt = 0; attempt < 4; ++attempt) {
        const auto rcid = (econ::CommodityId)(rng.nextU32() % (core::u32)econ::kCommodityCount);
        const double m2 = std::max(1e-6, econ::commodityDef(rcid).massKg);
        const double w2 = std::max(1.0, std::floor(capKg / m2));
        const double l2 = econ::takeInventory(fromE, fromSt.economyModel, rcid, w2);
        if (l2 > 0.0) {
          cid = rcid;
          t.tradeUnits = l2;
          t.tradeCommodity = rcid;
          t.tradeDestStationIndex = toIdx;
          const auto q = econ::quote(fromE, fromSt.economyModel, rcid, 0.10);
          t.cargoValueCr = std::max(0.0, l2 * q.mid);
          return true;
        }
      }
      t.tradeUnits = 0.0;
      t.tradeCommodity = cid;
      t.tradeDestStationIndex = toIdx;
      t.cargoValueCr = 0.0;
      return false;
    }

    t.tradeCommodity = cid;
    t.tradeUnits = loaded;
    t.tradeDestStationIndex = toIdx;

    // Approximate cargo value for piracy/loot purposes.
    const auto q = econ::quote(fromE, fromSt.economyModel, cid, 0.10);
    t.cargoValueCr = std::max(0.0, loaded * q.mid);
    return true;
  };

auto spawnPiratePack = [&](int maxCount) -> int {
  maxCount = std::max(1, std::min(maxCount, 3));
  const core::u64 group = std::max<core::u64>(1, rng.nextU64());

  // If the player has cargo, pirates sometimes try extortion first.
  const double playerCargoValueCr = cargoValueEstimateCr();
  const bool canDemand = (playerCargoValueCr > 180.0)
                      && !pirateDemand.active
                      && !docked
                      && (fsdState == FsdState::Idle)
                      && (supercruiseState == SupercruiseState::Idle);
  const bool willDemand = canDemand && (rng.nextUnit() < 0.65);

  // Pick a pack size (bounded by maxCount).
  int count = 1;
  if (maxCount >= 2 && rng.nextUnit() < 0.60) ++count;
  if (maxCount >= 3 && rng.nextUnit() < 0.35) ++count;

  // Spawn somewhere near the player, but not too close.
  const math::Vec3d baseDir = randDir();
  const double baseDistKm = rng.range(55000.0, 120000.0);
  const math::Vec3d basePosKm = ship.positionKm() + baseDir * baseDistKm;

  std::string leaderName;

  core::u64 leaderId = 0;
  for (int k = 0; k < count; ++k) {
    Contact p{};
    p.id = std::max<core::u64>(1, rng.nextU64());
    p.role = ContactRole::Pirate;
    p.hostileToPlayer = !willDemand; // hostile by default; negotiation packs hold fire
    p.groupId = group;
    p.leaderId = (k == 0) ? 0 : leaderId;
    if (k == 0) leaderId = p.id;

    p.name = (k == 0) ? ("Pirate Leader " + std::to_string((int)contacts.size() + 1))
                      : ("Pirate Wing " + std::to_string((int)contacts.size() + 1));
    if (k == 0) leaderName = p.name;

    const bool leader = (k == 0);
    const int tMk = leader ? 2 : ((rng.nextUnit() < 0.20) ? 2 : 1);
    const int dMk = leader ? 2 : 1;

    auto pickPirateWeapon = [&](bool isLeader) -> WeaponType {
      const double r = rng.nextUnit();
      if (isLeader) {
        if (r < 0.40) return WeaponType::Railgun;
        if (r < 0.85) return WeaponType::Cannon;
        return WeaponType::PulseLaser;
      }
      if (r < 0.48) return WeaponType::Cannon;
      if (r < 0.80) return WeaponType::BeamLaser;
      return WeaponType::PulseLaser;
    };

    configureContactLoadout(p,
                            ShipHullClass::Fighter,
                            tMk,
                            /*shieldMk=*/1,
                            dMk,
                            pickPirateWeapon(leader),
                            /*hullMul=*/leader ? 1.05 : 0.95,
                            /*shieldMul=*/leader ? 0.75 : 0.67,
                            /*regenMul=*/leader ? 0.70 : 0.63,
                            /*accelMul=*/leader ? 0.72 : 0.70,
                            /*aiSkill=*/leader ? 0.68 : 0.55);

    // Pirates like to keep weapons hot.
    p.pips = {1, 3, 2};
    sim::normalizePips(p.pips);

    // Slight spread so squads don't overlap perfectly.
    const math::Vec3d spread = randDir() * rng.range(1200.0, 7200.0) + math::Vec3d{0,1,0} * rng.range(-1200.0, 1200.0);
    p.ship.setPositionKm(basePosKm + spread);
    p.ship.setVelocityKmS(ship.velocityKmS());

    const math::Vec3d toP = ship.positionKm() - p.ship.positionKm();
    p.ship.setOrientation(quatFromTo({0,0,1}, (toP.lengthSq() > 1e-9) ? toP.normalized() : math::Vec3d{0,0,1}));

    contacts.push_back(std::move(p));
  }

  toast(toasts, (count > 1) ? "Contacts: pirate pack detected!" : "Contact: pirate detected!", 2.6);

  if (willDemand && leaderId != 0) {
    pirateDemand.active = true;
    pirateDemand.groupId = group;
    pirateDemand.leaderId = leaderId;
    pirateDemand.startDays = timeDays;
    pirateDemand.untilDays = timeDays + (22.0 / 86400.0);
    pirateDemand.requiredValueCr = std::round(std::clamp(160.0 + playerCargoValueCr * rng.range(0.09, 0.18), 200.0, 5000.0) / 10.0) * 10.0;
    pirateDemand.deliveredValueCr = 0.0;
    pirateDemand.leaderName = leaderName.empty() ? "Pirates" : leaderName;

    toast(toasts,
          pirateDemand.leaderName + ": jettison cargo worth ~" + std::to_string((int)std::round(pirateDemand.requiredValueCr))
            + " cr within 22s or be destroyed!",
          4.2);
  }
  return count;
};

auto spawnTrader = [&]() {

    Contact t{};
    t.id = std::max<core::u64>(1, rng.nextU64());
    t.role = ContactRole::Trader;
    t.factionId = localFaction;
    t.homeStationIndex = chooseHomeStation();
    t.name = "Trader " + std::to_string((int)contacts.size() + 1);
    configureContactLoadout(t,
                            ShipHullClass::Hauler,
                            /*thrMk=*/1,
                            /*shieldMk=*/1,
                            /*distMk=*/1,
                            /*weapon=*/WeaponType::MiningLaser,
                            /*hullMul=*/0.56,
                            /*shieldMul=*/0.62,
                            /*regenMul=*/0.0,
                            /*accelMul=*/0.72,
                            /*aiSkill=*/0.35);

    // Traders bias pips into engines for escape.
    t.pips = {4, 1, 1};
    sim::normalizePips(t.pips);

    // Cargo capacity (kg) varies per trader.
    t.tradeCapacityKg = rng.range(120.0, 420.0);

    // Traders have a high-speed shortcut when far from the player (keeps economy moving).
    t.tradeSupercruiseSpeedKmS = rng.range(7000.0, 14000.0);
    t.tradeSupercruiseDropDistKm = rng.range(70000.0, 140000.0);

    // Assign an initial hauling job (if possible) and load from source inventory.
    planTraderHaul(t, t.homeStationIndex);

    spawnNearStation(t.homeStationIndex, 18.0, 30.0, t.ship);

    contacts.push_back(std::move(t));
  };

auto spawnPolicePack = [&](int maxCount) -> int {
  if (localFaction == 0) return 0;
  maxCount = std::max(1, std::min(maxCount, 3));

  const core::u64 group = std::max<core::u64>(1, rng.nextU64());

  // Pack size depends on whether you're wanted / heat is high.
  int count = 1;
  if (maxCount >= 2 && rng.nextUnit() < (playerWantedHere ? 0.75 : 0.35)) ++count;
  if (maxCount >= 3 && rng.nextUnit() < ((playerWantedHere && policeHeat > 2.5) ? 0.55 : 0.20)) ++count;

  // "Tier" bumps stats a bit for high bounties / hot pursuits.
  const int tier = std::clamp(1 + (localBounty > 1200.0) + (localBounty > 4500.0) + (policeHeat > 3.0), 1, 3);
  const double tierMul = 1.0 + 0.15 * (double)(tier - 1);

  core::u64 leaderId = 0;
  for (int k = 0; k < count; ++k) {
    Contact p{};
    p.id = std::max<core::u64>(1, rng.nextU64());
    p.role = ContactRole::Police;
    p.groupId = group;
    p.leaderId = (k == 0) ? 0 : leaderId;
    if (k == 0) leaderId = p.id;

    p.factionId = localFaction;
    p.homeStationIndex = chooseHomeStation();

    if (tier >= 3 && k == 0) {
      p.name = "Interceptor " + std::to_string((int)contacts.size() + 1);
    } else {
      p.name = "Security " + std::to_string((int)contacts.size() + 1);
    }

    const bool leader = (k == 0);
    const WeaponType w = (tier >= 3 && leader) ? WeaponType::Railgun : WeaponType::PulseLaser;

    // Police ships are balanced to roughly match the legacy stats (90 hull/shield baseline),
    // but use the unified loadout + distributor model.
    const int thrMk = 2;
    const int shieldMk = 1;
    const int distMk = tier;

    const sim::ShipDerivedStats dsTmp = sim::computeShipDerivedStats(ShipHullClass::Scout, thrMk, shieldMk, distMk);
    const double targetRegenPerSimMin = (npcShieldRegenPerSec * 60.0) * tierMul;
    const double regenMul = (dsTmp.shieldRegenPerSimMin > 1e-9)
                              ? (targetRegenPerSimMin / dsTmp.shieldRegenPerSimMin)
                              : 0.0;

    configureContactLoadout(p,
                            ShipHullClass::Scout,
                            thrMk,
                            shieldMk,
                            distMk,
                            w,
                            /*hullMul=*/0.90 * tierMul,
                            /*shieldMul=*/0.90 * tierMul,
                            /*regenMul=*/regenMul,
                            /*accelMul=*/0.80 * tierMul,
                            /*aiSkill=*/leader ? (tier >= 3 ? 0.78 : 0.70) : (0.64 + 0.03 * (double)(tier - 1)));

    // Police bias a bit toward SYS for survival.
    p.pips = {2, 2, 2};
    sim::normalizePips(p.pips);

    // Wingmen spawn slightly further out to avoid immediate overlaps.
    spawnNearStation(p.homeStationIndex, 16.0 + 1.5 * (double)k, 22.0 + 2.0 * (double)k, p.ship);

    contacts.push_back(std::move(p));
  }
  return count;
};

  // --- Encounter scheduling (refactored into core) ---
  // The app still owns *how* contacts are instantiated; the core module only
  // decides *when* to request spawns.
  sim::EncounterDirectorCounts ec{};
  ec.alivePirates = alivePirates;
  ec.aliveTraders = aliveTraders;
  ec.alivePolice = alivePolice;
  ec.aliveTotal = aliveTotal;

  sim::EncounterDirectorContext ectx{};
  ectx.timeDays = timeDays;
  ectx.combatSimEnabled = combatSimEnabled;
  ectx.stationCount = (currentSystem && !currentSystem->stations.empty()) ? (int)currentSystem->stations.size() : 0;
  ectx.localFactionId = localFaction;
  ectx.localRep = localRep;
  ectx.localBountyCr = localBounty;
  ectx.policeHeat = policeHeat;
  ectx.policeAlert = (policeAlertUntilDays > timeDays);
  ectx.playerWantedHere = playerWantedHere;

  // Pirates occasionally (baseline threat).
  if (const int cap = sim::planPirateSpawn(encounterDirector, ectx, ec); cap > 0) {
    const int spawned = spawnPiratePack(cap);
    alivePirates += spawned;
    aliveTotal += spawned;
    ec.alivePirates += spawned;
    ec.aliveTotal += spawned;
  }

  // Traders / traffic: gives you something to pirate (but doing so triggers police).
  if (sim::planTraderSpawn(encounterDirector, ectx, ec) > 0) {
    spawnTrader();
    ++aliveTraders;
    ++aliveTotal;
    ++ec.aliveTraders;
    ++ec.aliveTotal;
  }

  // Police / patrols: scale with threat + your legal status.
  if (localFaction != 0) {
    const sim::PoliceSpawnPlan p = sim::planPoliceSpawn(encounterDirector, ectx, ec);
    if (p.spawnMaxCount > 0) {
      const int spawned = spawnPolicePack(p.spawnMaxCount);
      alivePolice += spawned;
      aliveTotal += spawned;
      ec.alivePolice += spawned;
      ec.aliveTotal += spawned;
    }
  }

  // Ensure mission bounty targets exist in their target system.
  if (!docked && !missions.empty()) {
    // Prune cached bounty targets that no longer correspond to an active bounty mission.
    if (!bountyTargetStates.empty()) {
      bountyTargetStates.erase(
        std::remove_if(
          bountyTargetStates.begin(),
          bountyTargetStates.end(),
          [&](const sim::BountyTargetState& bt) {
            for (const auto& m : missions) {
              if (m.completed || m.failed) continue;
              if (m.targetNpcId == 0) continue;
              if (!((m.type == sim::MissionType::BountyScan) || (m.type == sim::MissionType::BountyKill))) continue;
              if (m.targetNpcId == bt.targetId) return false;
            }
            return true;
          }),
        bountyTargetStates.end());
    }

    for (const auto& m : missions) {
      if (m.completed || m.failed) continue;
      if (!((m.type == sim::MissionType::BountyScan) || (m.type == sim::MissionType::BountyKill))) continue;
      if (m.toSystem != currentSystem->stub.id) continue;
      if (m.targetNpcId == 0) continue;

      // Best-effort: anchor the target near a "hideout" station to make the hunt less random.
      std::optional<std::size_t> hideoutIdx{};
      if (currentSystem && !currentSystem->stations.empty()) {
        if (m.toStation != 0) {
          for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
            if (currentSystem->stations[i].id == m.toStation) { hideoutIdx = i; break; }
          }
        }
        if (!hideoutIdx) hideoutIdx = 0;
      }

      bool present = false;
      for (auto& c : contacts) {
        if (c.alive && c.id == m.targetNpcId) {
          c.missionTarget = true;
          if (hideoutIdx) c.homeStationIndex = (int)*hideoutIdx;

          // Update the persisted cache so the target doesn't "reset" if the player leaves and returns.
          sim::BountyTargetState* bt = nullptr;
          for (auto& s : bountyTargetStates) {
            if (s.targetId == m.targetNpcId) { bt = &s; break; }
          }
          if (!bt) {
            bountyTargetStates.push_back(sim::BountyTargetState{});
            bt = &bountyTargetStates.back();
          }
          bt->targetId = m.targetNpcId;
          bt->missionId = m.id;
          bt->systemId = currentSystem ? currentSystem->stub.id : 0;
          bt->hideoutStation = m.toStation;
          bt->posKm = c.ship.positionKm();
          bt->velKmS = c.ship.velocityKmS();
          bt->orient = c.ship.orientation();
          bt->angVelRadS = c.ship.angularVelocityRadS();
          bt->hullFrac = std::clamp(c.hull / std::max(1e-6, c.hullMax), 0.0, 1.0);
          bt->shieldFrac = std::clamp(c.shield / std::max(1e-6, c.shieldMax), 0.0, 1.0);

          present = true;
          break;
        }
      }
      if (present) continue;

      // Spawn a distinct target pirate.
      Contact tgt{};
      tgt.id = m.targetNpcId;
      tgt.role = ContactRole::Pirate;
      tgt.missionTarget = true;
      tgt.name = "Bounty Target";
      if (hideoutIdx) tgt.homeStationIndex = (int)*hideoutIdx;

      // Loadout: kill targets are a bit tougher than scan targets.
      const bool killTarget = (m.type == sim::MissionType::BountyKill);
      configureContactLoadout(tgt,
                              ShipHullClass::Fighter,
                              /*thrMk=*/2,
                              /*shieldMk=*/killTarget ? 2 : 1,
                              /*distMk=*/2,
                              /*weapon=*/killTarget ? WeaponType::Railgun : WeaponType::PulseLaser,
                              /*hullMul=*/killTarget ? 1.12 : 1.02,
                              /*shieldMul=*/killTarget ? 0.95 : 0.75,
                              /*regenMul=*/killTarget ? 0.78 : 0.68,
                              /*accelMul=*/killTarget ? 0.78 : 0.74,
                              /*aiSkill=*/killTarget ? 0.82 : 0.72);

      // Bounty targets stay "cold" until you're in the area (see pirate AI above).
      tgt.hostileToPlayer = false;
      tgt.pips = {1, 3, 2};
      sim::normalizePips(tgt.pips);

      // Find-or-create the cached state record for this target.
      sim::BountyTargetState* bt = nullptr;
      for (auto& s : bountyTargetStates) {
        if (s.targetId == m.targetNpcId) { bt = &s; break; }
      }
      if (!bt) {
        bountyTargetStates.push_back(sim::BountyTargetState{});
        bt = &bountyTargetStates.back();
      }

      const sim::SystemId sysId = currentSystem ? currentSystem->stub.id : 0;
      const bool haveSaved = (bt->targetId == m.targetNpcId) && (bt->systemId == sysId);

      if (haveSaved) {
        tgt.ship.setPositionKm(bt->posKm);
        tgt.ship.setVelocityKmS(bt->velKmS);
        tgt.ship.setOrientation(bt->orient);
        tgt.ship.setAngularVelocityRadS(bt->angVelRadS);
        tgt.hull = tgt.hullMax * std::clamp(bt->hullFrac, 0.0, 1.0);
        tgt.shield = tgt.shieldMax * std::clamp(bt->shieldFrac, 0.0, 1.0);
      } else {
        // Deterministic initial spawn near the hideout station so the mission "search" is repeatable.
        core::u64 seed = core::hashCombine(core::fnv1a64("bounty_spawn"), (core::u64)m.targetNpcId);
        seed = core::hashCombine(seed, (core::u64)m.id);
        seed = core::hashCombine(seed, (core::u64)m.toStation);
        core::SplitMix64 brng(seed);

        if (hideoutIdx && currentSystem && *hideoutIdx < currentSystem->stations.size()) {
          const auto& st = currentSystem->stations[*hideoutIdx];
          const math::Vec3d stPos = stationPosKm(st, timeDays);
          const math::Vec3d stVel = stationVelKmS(st, timeDays);
          const math::Quatd stQ = stationOrient(st, stPos, timeDays);
          const math::Vec3d axis = stQ.rotate({0,0,1});
          const math::Vec3d tangent = stQ.rotate({1,0,0});

          const double distKm = brng.range(st.radiusKm * 22.0, st.radiusKm * 36.0);
          const double lateralKm = brng.range(-st.radiusKm * 3.0, st.radiusKm * 3.0);
          const math::Vec3d p = stPos + axis * distKm + tangent * lateralKm;

          tgt.ship.setPositionKm(p);
          tgt.ship.setVelocityKmS(stVel);
          const math::Vec3d toSt = (stPos - p);
          tgt.ship.setOrientation(quatFromTo({0,0,1}, (toSt.lengthSq() > 1e-9) ? toSt.normalized() : math::Vec3d{0,0,1}));
        } else {
          const math::Vec3d d = randDir();
          const double distKm = brng.range(65000.0, 150000.0);
          tgt.ship.setPositionKm(ship.positionKm() + d * distKm);
          tgt.ship.setVelocityKmS(ship.velocityKmS());
          tgt.ship.setOrientation(quatFromTo({0,0,1}, (-d).normalized()));
        }

        tgt.ship.setAngularVelocityRadS({0,0,0});
        tgt.hull = tgt.hullMax;
        tgt.shield = tgt.shieldMax;

        // Seed the cache so leaving/re-entering doesn't reset the target.
        bt->targetId = m.targetNpcId;
        bt->missionId = m.id;
        bt->systemId = sysId;
        bt->hideoutStation = m.toStation;
        bt->posKm = tgt.ship.positionKm();
        bt->velKmS = tgt.ship.velocityKmS();
        bt->orient = tgt.ship.orientation();
        bt->angVelRadS = tgt.ship.angularVelocityRadS();
        bt->hullFrac = 1.0;
        bt->shieldFrac = 1.0;
      }

      contacts.push_back(std::move(tgt));
      break; // spawn at most one target per frame
    }
  }

  // Find a nearby pirate index (for police).
  auto nearestPirateIndex = [&](const math::Vec3d& fromKm, double maxDistKm) -> std::optional<std::size_t> {
    double bestD = maxDistKm;
    std::optional<std::size_t> best{};
    for (std::size_t i = 0; i < contacts.size(); ++i) {
      const auto& c = contacts[i];
      if (!c.alive || c.role != ContactRole::Pirate) continue;
      const double d = (c.ship.positionKm() - fromKm).length();
      if (d < bestD) { bestD = d; best = i; }
    }
    return best;
  };

  auto chaseTarget = [&](sim::Ship& selfShip,
                         sim::ShipInput& ai,
                         const math::Vec3d& targetPosKm,
                         const math::Vec3d& targetVelKmS,
                         double desiredDistKm,
                         double maxSpeedKmS,
                         double faceGain) {
	    sim::FlightControlParams fp{};
	    fp.maxSpeedKmS = maxSpeedKmS;
	    fp.speedGain = 0.000004;
	    fp.velGain = 1.8;
	    fp.desiredDistKm = desiredDistKm;
	    fp.allowBoost = false;
	    fp.dampers = true;

	    sim::AttitudeControlParams ap{};
	    ap.faceGain = faceGain;
	    ap.alignUp = false;

	    // Intercept-course guidance reduces "tail chasing" when the target has significant
	    // lateral velocity (common in dogfights / pursuit).
	    sim::InterceptCourseParams ic{};
	    ic.enabled = true;
	    ic.maxLeadTimeSec = 180.0;
	    ic.minSpeedKmS = 0.05;
	    ic.useMaxSpeedForSolve = true;

	    const auto out = sim::chaseTargetIntercept(selfShip, targetPosKm, targetVelKmS, fp, ap, ic);
	    ai = out.input;
  };

  // Variant that lets AI keep chasing a point while aiming elsewhere (e.g. projectile lead).
  auto chaseTargetFace = [&](sim::Ship& selfShip,
                             sim::ShipInput& ai,
                             const math::Vec3d& targetPosKm,
                             const math::Vec3d& targetVelKmS,
                             double desiredDistKm,
                             double maxSpeedKmS,
                             double faceGain,
                             const math::Vec3d& desiredForwardWorld,
                             bool allowBoost) {
    sim::FlightControlParams fp{};
    fp.maxSpeedKmS = maxSpeedKmS;
    fp.speedGain = 0.000004;
    fp.velGain = 1.8;
    fp.desiredDistKm = desiredDistKm;
    fp.allowBoost = allowBoost;
    fp.dampers = true;

    sim::AttitudeControlParams ap{};
    ap.faceGain = faceGain;
    ap.alignUp = false;

    sim::InterceptCourseParams ic{};
    ic.enabled = true;
    ic.maxLeadTimeSec = 180.0;
    ic.minSpeedKmS = 0.05;
    ic.useMaxSpeedForSolve = true;

    const auto out = sim::approachTargetIntercept(selfShip, targetPosKm, targetVelKmS, fp, ap, desiredForwardWorld, ic);
    ai = out.input;
  };

  // Contacts AI + combat
  for (auto& c : contacts) {
    if (!c.alive) continue;

    // Keep contact dampers in the same local reference frame as the player.
    c.ship.setDampingFrameVelocityKmS(localFrameVelKmS);

    // Power distributor regen (NPCs share the same capacitor + pip model as the player).
    sim::stepDistributor(c.distributorState, c.distributorCfg, c.pips, dtSim);

    // cooldowns
    c.fireCooldown = std::max(0.0, c.fireCooldown - dtSim);

    // Shields regen (consumes SYS capacitor, scaled by SYS pips).
    if (c.alive && c.hull > 0.0 && c.shieldMax > 0.0 && c.shieldRegenPerSimMin > 0.0 && c.shield < c.shieldMax) {
      const double regenMul = sim::shieldRegenMultiplierFromPips(c.pips.sys);
      const double desired = c.shieldRegenPerSimMin * regenMul * (dtSim / 60.0);
      if (desired > 0.0) {
        const double costPerPt = c.distributorCfg.shieldRegenCostPerPoint;
        const double affordable = (costPerPt > 0.0) ? (c.distributorState.sys / costPerPt) : desired;
        const double actual = std::clamp(desired, 0.0, affordable);
        if (actual > 0.0) {
          c.shield = std::min(c.shieldMax, c.shield + actual);
          if (costPerPt > 0.0) c.distributorState.sys = std::max(0.0, c.distributorState.sys - actual * costPerPt);
        }
      }
    }

    // ---- PIRATES ----
    if (c.role == ContactRole::Pirate) {
      sim::ShipInput ai{};

      // Mission-critical bounty targets shouldn't immediately chase the player across the entire system.
      // If the player is far from the target's "hideout" station, keep the target loitering there.
      if (c.missionTarget && c.attackTargetId == 0 && currentSystem && !currentSystem->stations.empty()) {
        const std::size_t stIdx = std::min<std::size_t>((std::size_t)c.homeStationIndex, currentSystem->stations.size() - 1);
        const auto& st = currentSystem->stations[stIdx];
        const math::Vec3d stPos = stationPosKm(st, timeDays);

        // Engage when the player is plausibly "in the area" (or after being fired upon).
        const double engageDistKm = std::max(90000.0, st.commsRangeKm * 1.25);
        const double playerDistKm = (ship.positionKm() - stPos).length();

        if (playerDistKm > engageDistKm && timeDays >= c.underFireUntilDays) {
          const double baseR = std::max(12000.0, st.radiusKm * 24.0);
          const double ang = std::fmod((double)(c.id % 100000) * 0.000173 + timeDays * 0.75, 6.283185307179586);
          const math::Vec3d offset = math::Vec3d{std::cos(ang), 0.14 * std::sin(ang * 0.7), std::sin(ang)} * baseR;

          chaseTarget(c.ship,
                     ai,
                     stPos + offset,
                     stationVelKmS(st, timeDays),
                     /*desiredDistKm=*/baseR * 0.65,
                     /*tangentialStrength=*/0.26,
                     /*maxAccelFrac=*/1.4);
          c.ship.step(dtSim, ai);
          continue;
        }
      }

      const bool hostile = c.hostileToPlayer || c.missionTarget || (c.attackTargetId != 0);

      // Resolve attack target: by default pirates harass the player, but mission-driven pirates
      // can be pointed at a specific contact (e.g. an escort convoy).
      const sim::Ship* targetShip = &ship;
      Contact* targetContact = nullptr;
      if (c.attackTargetId != 0) {
        for (auto& t : contacts) {
          if (!t.alive) continue;
          if (t.id == c.attackTargetId) {
            targetContact = &t;
            targetShip = &t.ship;
            break;
          }
        }
      }

      const math::Vec3d tgtPos = targetShip->positionKm();
      const math::Vec3d tgtVel = targetShip->velocityKmS();

      const bool fleeing = (timeDays < c.fleeUntilDays);
      if (fleeing) {
        // Break off: run away from the player.
        const math::Vec3d away = (c.ship.positionKm() - ship.positionKm());
        const math::Vec3d dir = (away.lengthSq() > 1e-6) ? away.normalized() : math::Vec3d{0,0,1};
        chaseTarget(c.ship, ai, c.ship.positionKm() + dir * 240000.0, ship.velocityKmS(), 0.0, 0.28, 1.3);
        c.ship.step(dtSim, ai);
        continue;
      }

      // Small squad spread so packs don't sit on top of each other.
      const math::Vec3d toP = tgtPos - c.ship.positionKm();
      const double distP = toP.length();
      const math::Vec3d toPN = (distP > 1e-6) ? (toP / distP) : math::Vec3d{0,0,1};
      math::Vec3d right = math::cross(toPN, math::Vec3d{0,1,0});
      if (right.lengthSq() < 1e-9) right = math::cross(toPN, math::Vec3d{1,0,0});
      right = right.normalized();
      const math::Vec3d up = math::cross(right, toPN).normalized();

      const double spreadOn = (c.groupId != 0) ? 1.0 : 0.0;
      const double s = spreadOn * (((int)(c.id % 7) - 3) * 2200.0);
      const double u = spreadOn * (((int)(c.id % 5) - 2) * 1500.0);
      const math::Vec3d aimPos = tgtPos + right * s + up * u;

      const double standOffKm = 35000.0 + ((c.leaderId != 0) ? 5000.0 : 0.0);
      const WeaponDef& w = weaponDef(c.weapon);

      // Aim at lead point for projectiles (so pirates with cannons/railguns actually hit).
      math::Vec3d aimPointKm = tgtPos;
      if (!w.beam && w.projSpeedKmS > 1e-6) {
        const double ttl = w.rangeKm / w.projSpeedKmS;
        if (auto sol = sim::solveProjectileLead(c.ship.positionKm(), c.ship.velocityKmS(), tgtPos, tgtVel, w.projSpeedKmS, ttl)) {
          aimPointKm = sol->leadPointKm;
        }
      }

      chaseTargetFace(c.ship, ai, aimPos, tgtVel, standOffKm, 0.22, 1.8, (aimPointKm - c.ship.positionKm()), /*allowBoost*/false);
      c.ship.step(dtSim, ai);

      // Fire if aligned and we have weapon capacitor.
      const math::Vec3d to = tgtPos - c.ship.positionKm();
      const double dist = to.length();
      if (hostile && c.fireCooldown <= 0.0 && dist < w.rangeKm) {
        // Recompute aim point (post-step) so the dot test matches current geometry.
        math::Vec3d fireAimPointKm = tgtPos;
        if (!w.beam && w.projSpeedKmS > 1e-6) {
          const double ttl = w.rangeKm / w.projSpeedKmS;
          if (auto sol = sim::solveProjectileLead(c.ship.positionKm(), c.ship.velocityKmS(), tgtPos, tgtVel, w.projSpeedKmS, ttl)) {
            fireAimPointKm = sol->leadPointKm;
          }
        }

        const math::Vec3d toAim = fireAimPointKm - c.ship.positionKm();
        const double distAim = toAim.length();
        if (distAim > 1e-6) {
          const math::Vec3d toAimN = toAim / distAim;
          const double aim = math::dot(c.ship.forward().normalized(), toAimN);

          // Lower skill => stricter aim requirement (so they don't hit with perfect turret tracking).
          // Missiles have a wider "acceptable" fire cone since they steer after launch.
          double baseThresh = 0.995;
          if (w.beam) baseThresh = 0.992;
          else if (w.guided) baseThresh = 0.975;
          const double aimThresh = std::clamp(baseThresh + (1.0 - c.aiSkill) * 0.010, 0.990, 0.999);

          if (aim > aimThresh) {
            const double capCost = sim::weaponCapacitorCost(w);
            if (c.distributorState.wep + 1e-9 >= capCost) {
              std::vector<sim::SphereTarget> fireTargets;
              if (w.beam) {
                fireTargets.reserve(1 + contacts.size() + asteroids.size());

                // Player target
                {
                  sim::SphereTarget t{};
                  t.kind = sim::CombatTargetKind::Player;
                  t.index = 0;
                  t.id = 0;
                  t.centerKm = ship.positionKm();
                  t.velKmS = ship.velocityKmS();
                  t.radiusKm = 900.0;
                  t.minAimCos = -1.0;
                  fireTargets.push_back(t);
                }

                for (std::size_t i = 0; i < contacts.size(); ++i) {
                  const auto& oc = contacts[i];
                  if (!oc.alive) continue;
                  sim::SphereTarget t{};
                  t.kind = sim::CombatTargetKind::Ship;
                  t.index = i;
                  t.id = oc.id;
                  t.centerKm = oc.ship.positionKm();
                  t.velKmS = oc.ship.velocityKmS();
                  t.radiusKm = 900.0;
                  t.minAimCos = -1.0;
                  fireTargets.push_back(t);
                }

                for (std::size_t i = 0; i < asteroids.size(); ++i) {
                  const auto& a = asteroids[i];
                  sim::SphereTarget t{};
                  t.kind = sim::CombatTargetKind::Asteroid;
                  t.index = i;
                  t.id = a.id;
                  t.centerKm = a.posKm;
                  t.velKmS = {0,0,0};
                  t.radiusKm = a.radiusKm;
                  t.minAimCos = -1.0;
                  fireTargets.push_back(t);
                }
              } else if (w.guided) {
                // Missiles: lock to the resolved attack target.
                fireTargets.reserve(1);
                if (targetContact) {
                  sim::SphereTarget t{};
                  t.kind = sim::CombatTargetKind::Ship;
                  t.index = (std::size_t)(targetContact - contacts.data());
                  t.id = targetContact->id;
                  t.centerKm = targetContact->ship.positionKm();
                  t.velKmS = targetContact->ship.velocityKmS();
                  t.radiusKm = 900.0;
                  t.minAimCos = 0.0;
                  fireTargets.push_back(t);
                } else {
                  sim::SphereTarget t{};
                  t.kind = sim::CombatTargetKind::Player;
                  t.index = 0;
                  t.id = 0;
                  t.centerKm = ship.positionKm();
                  t.velKmS = ship.velocityKmS();
                  t.radiusKm = 900.0;
                  t.minAimCos = 0.0;
                  fireTargets.push_back(t);
                }
              }

              const sim::FireResult fr = sim::tryFireWeapon(
                c.ship,
                c.weapon,
                c.fireCooldown,
                c.distributorMk,
                c.id,
                /*fromPlayer*/false,
                fireTargets.empty() ? nullptr : fireTargets.data(),
                fireTargets.size()
              );

              if (fr.fired) {
                c.distributorState.wep = std::max(0.0, c.distributorState.wep - capCost);
                c.fireCooldown = fr.newCooldownSimSec;

                if (fr.hasBeam) {
                  beams.push_back({toRenderU(fr.beam.aKm), toRenderU(fr.beam.bKm), fr.beam.r, fr.beam.g, fr.beam.b, 0.10});
                }
                if (fr.hasProjectile) {
                  projectiles.push_back(fr.projectile);
                }
                if (fr.hasMissile) {
                  missiles.push_back(fr.missile);
                }

                if (w.beam && fr.hit) {
                  if (fr.hitKind == sim::CombatTargetKind::Player) {
                    applyDamage(w.dmg, playerShield, playerHull);
                  } else if (fr.hitKind == sim::CombatTargetKind::Ship && fr.hitIndex < contacts.size()) {
                    auto& v = contacts[fr.hitIndex];
                    if (v.alive) {
                      applyDamage(w.dmg, v.shield, v.hull);
                      v.underFireUntilDays = timeDays + (18.0 / 86400.0);
                      if (v.hull <= 0.0) v.alive = false;
                    }
                  }
                }
              }
            }
          }
        }
      }
      continue;
    }

    // ---- TRADERS ----
    if (c.role == ContactRole::Trader) {
      sim::ShipInput ai{};
      ai.dampers = true;

      // Distress victims are "disabled": keep them stationary until the player helps.
      if (c.distressVictim && c.distressNeedUnits > 1e-6) {
        c.ship.step(dtSim, ai);
        continue;
      }

      double dtStep = dtSim;

      const bool fleeing = (timeDays < c.fleeUntilDays);
      if (fleeing) {
        // Run away from player.
        const math::Vec3d away = (c.ship.positionKm() - ship.positionKm());
        const math::Vec3d dir = (away.lengthSq() > 1e-6) ? away.normalized() : math::Vec3d{0,0,1};
        chaseTarget(c.ship, ai, c.ship.positionKm() + dir * 200000.0, ship.velocityKmS(), 0.0, 0.25, 1.4);
      } else if (currentSystem && !currentSystem->stations.empty()) {
        // Hauling between stations in the current system (if possible).
        if (currentSystem->stations.size() >= 2) {
          std::size_t destIdx = std::min(c.tradeDestStationIndex, currentSystem->stations.size() - 1);
          const auto& st = currentSystem->stations[destIdx];
          const math::Vec3d stPos = stationPosKm(st, timeDays);

          const math::Vec3d stVel = stationVelKmS(st, timeDays);

          // Approach to comms-range and "trade" (we simulate docking by adjusting inventories).
          const double desiredDist = std::clamp(st.commsRangeKm * 0.55, st.radiusKm * 6.0, st.radiusKm * 14.0);

          double distKm = (c.ship.positionKm() - stPos).length();
          const double playerDistKm = (c.ship.positionKm() - ship.positionKm()).length();

          // Trader "supercruise": when far away (and not in the player's immediate vicinity),
          // move quickly along the route so in-system economy actually changes during play.
          const bool allowSuper = docked || (playerDistKm > 160000.0);
          if (allowSuper && distKm > desiredDist + 60000.0) {
            const math::Vec3d to = (stPos - c.ship.positionKm());
            if (to.lengthSq() > 1e-6) {
              const math::Vec3d n = to.normalized();
              const double speed = std::max(1500.0, c.tradeSupercruiseSpeedKmS);
              const double stepKm = std::min(distKm - desiredDist, speed * dtSim);

              c.ship.setPositionKm(c.ship.positionKm() + n * stepKm);
              c.ship.setVelocityKmS(stVel);
              c.ship.setAngularVelocityRadS({0,0,0});
              c.ship.setOrientation(quatFromTo({0,0,1}, n));

              // We already moved this frame.
              dtStep = 0.0;
              distKm = (c.ship.positionKm() - stPos).length();
            }
          } else {
            // Normal flight when the player is near (keeps things readable/physical).
            chaseTarget(c.ship, ai, stPos, stVel, desiredDist, 6.0, 1.0);
          }

          const bool arrived = distKm < st.commsRangeKm * 0.65;

          if (arrived && timeDays >= c.tradeCooldownUntilDays) {
            auto& stEcon = universe.stationEconomy(st, timeDays);

            if (c.escortConvoy) {
              // Mission convoys don't directly mutate station inventories. We reserve goods when
              // the mission is accepted and settle inventories via MissionLogic upon completion.
              c.tradeUnits = 0.0;
              c.tradeCargoValueCr = 0.0;
              c.cargoValueCr = 0.0;
              c.tradeCooldownUntilDays = timeDays + (45.0 / 86400.0);
            } else {
              // Deliver current cargo (if any).
              if (c.tradeUnits > 0.0) {
                const double before = c.tradeUnits;
                const double delivered = econ::addInventory(stEcon, st.economyModel, c.tradeCommodity, c.tradeUnits);
                c.tradeUnits = std::max(0.0, c.tradeUnits - delivered);
                if (c.tradeUnits < 1e-6) c.tradeUnits = 0.0;

                // If storage was full and we couldn't unload everything, pick a different buyer.
                if (c.tradeUnits > 0.0 && delivered + 1e-6 < before) {
                  c.tradeDestStationIndex = chooseOtherStation(destIdx);
                }
              }

              // If empty, load a new haul from this station for the next leg.
              if (c.tradeUnits <= 0.0) {
                planTraderHaul(c, destIdx);
              }

              // Update loot value for piracy / UI.
              if (c.tradeUnits > 0.0) {
                const auto q = econ::quote(stEcon, st.economyModel, c.tradeCommodity, 0.10);
                c.cargoValueCr = std::max(0.0, c.tradeUnits * q.mid);
                c.tradeCargoValueCr = c.cargoValueCr;
              } else {
                c.cargoValueCr = 0.0;
                c.tradeCargoValueCr = 0.0;
              }

              // Cooldown to prevent multiple trades in one arrival.
              c.tradeCooldownUntilDays = timeDays + (20.0 / 86400.0);
            }
          }
        } else {
          // Single-station systems: lazy orbit/patrol around their home station.
          const auto& st = currentSystem->stations[std::min(c.homeStationIndex, currentSystem->stations.size()-1)];
          const math::Vec3d stPos = stationPosKm(st, timeDays);
          const double baseR = st.radiusKm * 26.0;
          const double ang = std::fmod((double)(c.id % 1000) * 0.01 + timeDays * 0.75, 2.0 * math::kPi);
          const math::Vec3d offset = math::Vec3d{std::cos(ang), 0.15*std::sin(ang*0.7), std::sin(ang)} * baseR;
          chaseTarget(c.ship, ai, stPos + offset, stationVelKmS(st, timeDays), baseR * 0.6, 0.18, 1.2);
        }
      }

      c.ship.step(dtStep, ai);
      continue;
    }

    // ---- POLICE ----
    if (c.role == ContactRole::Police) {
	      const double bountyHere = (c.factionId != 0) ? getBounty(c.factionId) : 0.0;
	      const bool wantedHere = bountyHere > 1e-6;
	      const bool repHostile = (c.factionId != 0) && (getRep(c.factionId) < -45.0);
	      const bool hostile = c.hostileToPlayer || repHostile || wantedHere;

	      // If you're wanted, police will try to scan and demand compliance before firing (unless you attacked them).
	      const bool demandExpired = wantedHere && (policeDemand.factionId == c.factionId) && !policeDemand.active
	                                && (policeDemand.untilDays > 0.0) && (timeDays >= policeDemand.untilDays);
	      const bool holdFireForDemand = wantedHere && !c.hostileToPlayer && !demandExpired;
	      const bool holdFireForScan = !c.hostileToPlayer && cargoScanActive && (cargoScanFactionId == c.factionId);
	      const bool holdFireForBribe = !c.hostileToPlayer && bribeOffer.active && (bribeOffer.factionId == c.factionId) && (timeDays < bribeOffer.untilDays);
	      const bool canFireAtPlayer = hostile && !(holdFireForDemand || holdFireForScan || holdFireForBribe);

	      // If not hostile to player, try to engage pirates near the player.
	      std::optional<std::size_t> pirateIdx{};
	      if (!hostile) {
	        pirateIdx = nearestPirateIndex(c.ship.positionKm(), 140000.0);
	      }

	      // Escort behavior: some police ships are spawned as convoy escorts and should follow a
	      // designated contact (usually a mission convoy) when not actively engaging pirates.
	      const Contact* followTarget = nullptr;
	      if (!hostile && c.followId != 0) {
	        for (const auto& other : contacts) {
	          if (other.alive && other.id == c.followId) {
	            followTarget = &other;
	            break;
	          }
	        }
	      }

	      // Decide who we're engaging this frame.
	      const sim::WeaponDef w = sim::weaponDef(c.weapon);
	      const sim::Ship* engageShip = nullptr;
	      const Contact* engageContact = nullptr;
	      const bool chasePlayer = hostile;
	      const bool shootPlayer = canFireAtPlayer;
	      if (chasePlayer) {
	        engageShip = &ship;
	      } else if (pirateIdx) {
	        engageContact = &contacts[*pirateIdx];
	        engageShip = &engageContact->ship;
	      }

	      sim::ShipInput ai{};
	      if (engageShip) {
	        const math::Vec3d tgtPos = engageShip->positionKm();
	        const math::Vec3d tgtVel = engageShip->velocityKmS();
	        const double standoff = 42000.0 + ((c.leaderId != 0) ? 7000.0 : 0.0) + (double)((c.id % 3) * 1500);

	        // Aim at either the current position (beams) or a predicted intercept point (projectiles).
	        math::Vec3d aimPointKm = tgtPos;
	        if (!w.beam && w.projSpeedKmS > 1e-6) {
	          const double ttl = w.rangeKm / w.projSpeedKmS;
	          if (auto lead = sim::solveProjectileLead(c.ship.positionKm(), c.ship.velocityKmS(), tgtPos, tgtVel, w.projSpeedKmS, ttl)) {
	            aimPointKm = lead->leadPointKm;
	          }
	        }
	        chaseTargetFace(c.ship, ai, tgtPos, tgtVel, standoff, 0.26, 2.0, aimPointKm - c.ship.positionKm(), false);
	      } else if (followTarget) {
	        const double standoff = 36000.0 + ((c.leaderId != 0) ? 5500.0 : 0.0) + (double)((c.id % 3) * 1200);
	        chaseTarget(c.ship, ai, followTarget->ship.positionKm(), followTarget->ship.velocityKmS(), standoff, 0.22, 1.8);
	      } else if (currentSystem && !currentSystem->stations.empty()) {
	        // Patrol around home station.
	        const auto& st = currentSystem->stations[std::min(c.homeStationIndex, currentSystem->stations.size()-1)];
	        const math::Vec3d stPos = stationPosKm(st, timeDays);
	        const double baseR = st.radiusKm * 22.0;
	        const double ang = std::fmod((double)(c.id % 1000) * 0.012 + timeDays * 1.10, 2.0 * math::kPi);
	        const math::Vec3d offset = math::Vec3d{std::cos(ang), 0.10*std::sin(ang*0.6), std::sin(ang)} * baseR;
	        chaseTarget(c.ship, ai, stPos + offset, stationVelKmS(st, timeDays), baseR * 0.6, 0.20, 1.6);
	      }

	      c.ship.step(dtSim, ai);

	      // Fire if aligned at the chosen target.
	      if (engageShip && (shootPlayer || engageContact)) {
	        const math::Vec3d tgtPos = engageShip->positionKm();
	        const math::Vec3d tgtVel = engageShip->velocityKmS();
	        const double dist = (tgtPos - c.ship.positionKm()).length();
	        if (c.fireCooldown <= 0.0 && dist < w.rangeKm) {
	          math::Vec3d aimPointKm = tgtPos;
	          if (!w.beam && w.projSpeedKmS > 1e-6) {
	            const double ttl = w.rangeKm / w.projSpeedKmS;
	            if (auto lead = sim::solveProjectileLead(c.ship.positionKm(), c.ship.velocityKmS(), tgtPos, tgtVel, w.projSpeedKmS, ttl)) {
	              aimPointKm = lead->leadPointKm;
	            }
	          }

	          const math::Vec3d toAim = aimPointKm - c.ship.positionKm();
	          const double distAim = toAim.length();
	          const math::Vec3d toAimN = (distAim > 1e-6) ? (toAim / distAim) : math::Vec3d{0,0,1};
	          const double aim = math::dot(c.ship.forward().normalized(), toAimN);
	          double baseThr = 0.995;
	          if (w.beam) baseThr = 0.993;
	          else if (w.guided) baseThr = 0.975;
	          const double aimThr = std::clamp(baseThr + (1.0 - c.aiSkill) * 0.010, 0.975, 0.9995);

	          if (aim > aimThr) {
	            const double capCost = sim::weaponCapacitorCost(w);
	            if (c.distributorState.wep + 1e-9 >= capCost) {
	              std::vector<sim::SphereTarget> fireTargets;
	              if (w.beam) {
	                fireTargets.reserve(1 + contacts.size() + asteroids.size());
	                if (shootPlayer) {
	                  sim::SphereTarget pt{};
	                  pt.kind = sim::CombatTargetKind::Player;
	                  pt.index = 0;
	                  pt.id = 0;
	                  pt.centerKm = ship.positionKm();
	                  pt.velKmS = ship.velocityKmS();
	                  pt.radiusKm = 900.0;
	                  pt.minAimCos = -1.0;
	                  fireTargets.push_back(pt);
	                }
	                for (std::size_t i = 0; i < contacts.size(); ++i) {
	                  const auto& t = contacts[i];
	                  if (!t.alive) continue;
	                  sim::SphereTarget st{};
	                  st.kind = sim::CombatTargetKind::Ship;
	                  st.index = i;
	                  st.id = t.id;
	                  st.centerKm = t.ship.positionKm();
	                  st.velKmS = t.ship.velocityKmS();
	                  st.radiusKm = 900.0;
	                  st.minAimCos = -1.0;
	                  fireTargets.push_back(st);
	                }
	                for (std::size_t i = 0; i < asteroids.size(); ++i) {
	                  const auto& a = asteroids[i];
	                  sim::SphereTarget at{};
	                  at.kind = sim::CombatTargetKind::Asteroid;
	                  at.index = i;
	                  at.id = a.id;
	                  at.centerKm = a.posKm;
	                  at.velKmS = {0,0,0};
	                  at.radiusKm = a.radiusKm;
	                  at.minAimCos = -1.0;
	                  fireTargets.push_back(at);
	                }
	              } else if (w.guided) {
	                // Missiles: lock the chosen engagement target.
	                fireTargets.reserve(1);
	                if (shootPlayer) {
	                  sim::SphereTarget pt{};
	                  pt.kind = sim::CombatTargetKind::Player;
	                  pt.index = 0;
	                  pt.id = 0;
	                  pt.centerKm = ship.positionKm();
	                  pt.velKmS = ship.velocityKmS();
	                  pt.radiusKm = 900.0;
	                  pt.minAimCos = 0.0;
	                  fireTargets.push_back(pt);
	                } else if (engageContact) {
	                  sim::SphereTarget st{};
	                  st.kind = sim::CombatTargetKind::Ship;
	                  st.index = (std::size_t)(engageContact - contacts.data());
	                  st.id = engageContact->id;
	                  st.centerKm = engageContact->ship.positionKm();
	                  st.velKmS = engageContact->ship.velocityKmS();
	                  st.radiusKm = 900.0;
	                  st.minAimCos = 0.0;
	                  fireTargets.push_back(st);
	                }
	              }

	              const sim::FireResult fr = sim::tryFireWeapon(c.ship,
	                                                          c.weapon,
	                                                          c.fireCooldown,
	                                                          c.distributorMk,
	                                                          c.id,
	                                                          /*fromPlayer=*/false,
	                                                          fireTargets.empty() ? nullptr : fireTargets.data(),
	                                                          (int)fireTargets.size());
	              if (fr.fired) {
	                c.distributorState.wep = std::max(0.0, c.distributorState.wep - capCost);
	                c.fireCooldown = fr.newCooldownSimSec;
	                if (fr.hasBeam) {
	                  beams.push_back({toRenderU(fr.beam.aKm), toRenderU(fr.beam.bKm), fr.beam.r, fr.beam.g, fr.beam.b, 0.10});
	                }
	                if (fr.hasProjectile) projectiles.push_back(fr.projectile);
	                if (fr.hasMissile) missiles.push_back(fr.missile);

	                if (fr.hit && w.beam) {
	                  if (fr.hitKind == sim::CombatTargetKind::Player) {
	                    if (shootPlayer) {
	                      applyDamage(w.dmg, playerShield, playerHull);
	                    }
	                  } else if (fr.hitKind == sim::CombatTargetKind::Ship) {
	                    if (fr.hitIndex >= 0 && (std::size_t)fr.hitIndex < contacts.size()) {
	                      auto& t = contacts[(std::size_t)fr.hitIndex];
	                      if (t.alive) {
	                        applyDamage(w.dmg, t.shield, t.hull);
	                        t.underFireUntilDays = timeDays + (14.0 / 86400.0);
	                        if (t.hull <= 0.0) {
	                          t.alive = false;
	                          if (t.role == ContactRole::Pirate) {
	                            credits += 180.0;
	                            toast(toasts, "Security destroyed a pirate (+180).", 2.0);
	                          }
	                        }
	                      }
	                    }
	                  }
	                }
	              }
	            }
	          }
	        }
	      }
      continue;
    }
  }

  // Station turrets:
  // - Help against pirates
  // - If you are WANTED with the station's faction, the station will also chip at you near the no-fire zone.
  for (const auto& st : currentSystem->stations) {
    const math::Vec3d stPos = stationPosKm(st, timeDays);
    const double zoneKm = st.radiusKm * 25.0;
    const double distShip = (ship.positionKm() - stPos).length();

    // Only if player is nearby.
    if (distShip > zoneKm) continue;

    // Shoot pirates in zone
    for (auto& c : contacts) {
      if (!c.alive || c.role != ContactRole::Pirate) continue;
      const double d = (c.ship.positionKm() - stPos).length();
      if (d < zoneKm) {
        applyDamage(4.0 * dtSim, c.shield, c.hull);
        if (c.hull <= 0.0) {
          c.alive = false;
          credits += 250.0;
          toast(toasts, "Station defenses destroyed a pirate (+250).", 2.5);
        }
      }
    }

	    // Station security vs wanted player
	    if (st.factionId != 0 && getBounty(st.factionId) > 0.0) {
	      const bool holdFire = (policeDemand.factionId == st.factionId && policeDemand.active && timeDays < policeDemand.untilDays);
	      if (!holdFire) {
	        // very light pressure - enough to create urgency without insta-kill
	        applyDamage(1.5 * dtSim, playerShield, playerHull);
	      }
	    }
  }
}

      // Collisions (player with station hull)
      if (!docked) {
        for (const auto& st : currentSystem->stations) {
          const math::Vec3d stPos = stationPosKm(st, timeDays);
          const math::Quatd stQ = stationOrient(st, stPos, timeDays);
          const math::Vec3d relLocal = stQ.conjugate().rotate(ship.positionKm() - stPos);

          if (sim::insideStationHullExceptSlot(st, relLocal)) {
            // Damage based on relative speed and push outward slightly.
            const math::Vec3d stV = stationVelKmS(st, timeDays);
            const double relSpeed = (ship.velocityKmS() - stV).length();
            applyDamage(relSpeed * 18.0, playerShield, playerHull);

            // Push out along local axis with largest penetration.
            const double wx = st.radiusKm * 0.70;
            const double wy = st.radiusKm * 0.70;
            const double wz = st.radiusKm * 1.10;

            double dx = wx - std::abs(relLocal.x);
            double dy = wy - std::abs(relLocal.y);
            double dz = wz - std::abs(relLocal.z);

            math::Vec3d pushLocal{0,0,0};
            if (dx <= dy && dx <= dz) pushLocal.x = (relLocal.x >= 0 ? 1 : -1) * (dx + 200.0);
            else if (dy <= dz) pushLocal.y = (relLocal.y >= 0 ? 1 : -1) * (dy + 200.0);
            else pushLocal.z = (relLocal.z >= 0 ? 1 : -1) * (dz + 200.0);

            ship.setPositionKm(ship.positionKm() + stQ.rotate(pushLocal));
            ship.setVelocityKmS(stV); // kill relative motion on impact
            toast(toasts, "Collision!", 1.2);
            break;
          }
        }
      }

    // --- Cargo scans / contraband ---
    {
      const core::u32 jurisdiction = currentSystem ? currentSystem->stub.factionId : 0;
      const bool canScan = (!docked && fsdState == FsdState::Idle && supercruiseState == SupercruiseState::Idle && jurisdiction != 0);

      auto cancelCargoScan = [&]() {
        cargoScanActive = false;
        cargoScanProgressSec = 0.0;
        cargoScanDurationSec = 0.0;
        cargoScanRangeKm = 0.0;
        cargoScanStationId = 0;
        cargoScanContactId = 0;
        cargoScanFactionId = 0;
        cargoScanSourceName.clear();
      };

      // Resolve any pending bribe/compliance window first (from a completed contraband scan).
      if (bribeOffer.active) {
        auto issueFleeBounty = [&](const std::string& why, bool escalateLocal) {
          const sim::LawProfile law = lawForFaction(bribeOffer.factionId);
          const auto res = sim::evadeContraband(law, bribeOffer.fineCr, bribeOffer.illegalValueCr, escalateLocal);

          addRep(bribeOffer.factionId, res.repPenalty);
          // Running turns the fine into a bounty. (You keep the cargo, but you're now WANTED.)
          addBounty(bribeOffer.factionId, res.bountyAddedCr);

          if (res.policeHeatDelta > 1e-9) {
            policeHeat = std::clamp(policeHeat + res.policeHeatDelta, 0.0, 6.0);
          }
          if (res.policeAlertSeconds > 1e-9) {
            policeAlertUntilDays = std::max(policeAlertUntilDays, timeDays + (res.policeAlertSeconds / 86400.0));
          }
          if (res.nextPoliceSpawnDelaySeconds > 1e-9) {
            encounterDirector.nextPoliceSpawnDays = std::min(encounterDirector.nextPoliceSpawnDays, timeDays + (res.nextPoliceSpawnDelaySeconds / 86400.0));
          }

          heatImpulse += res.shipHeatDelta;

          const std::string src = bribeOffer.sourceName.empty() ? std::string("Security") : bribeOffer.sourceName;
          toast(toasts,
                src + ": " + why + " Bounty issued (" + std::to_string((int)std::round(bribeOffer.fineCr)) + " cr).",
                3.6);

          bribeOffer = BribeOffer{};
        };

        // Left the issuing faction's jurisdiction: can't confiscate, but the bounty stands.
        if (jurisdiction != bribeOffer.factionId) {
          issueFleeBounty("you fled the scan.", false);
        } else if (docked) {
          // Docking while under a scan demand triggers immediate enforcement.
          enforceContraband(bribeOffer.factionId,
                            bribeOffer.sourceName,
                            bribeOffer.illegalValueCr,
                            bribeOffer.detail,
                            bribeOffer.scannedIllegal);
          bribeOffer = BribeOffer{};
        } else {
          // Distance-based escape: if the scanning contact is lost (or you run out of range),
          // the authorities issue a bounty but don't confiscate cargo.
          bool outOfRange = false;
          if (bribeOffer.scannerContactId != 0 && bribeOffer.scannerRangeKm > 1e-6) {
            bool found = false;
            for (const auto& c : contacts) {
              if (!c.alive) continue;
              if (c.id != bribeOffer.scannerContactId) continue;
              found = true;
              const double distKm = (c.ship.positionKm() - ship.positionKm()).length();
              outOfRange = distKm > (bribeOffer.scannerRangeKm * 1.25);
              break;
            }
            if (!found) outOfRange = true;
          }

          if (!canScan) {
            issueFleeBounty("no compliance detected.", true);
          } else if (outOfRange) {
            issueFleeBounty("scan contact lost.", true);
          } else if (timeDays > bribeOffer.untilDays) {
            // Time's up: enforce confiscation + fine.
            enforceContraband(bribeOffer.factionId,
                              bribeOffer.sourceName,
                              bribeOffer.illegalValueCr,
                              bribeOffer.detail,
                              bribeOffer.scannedIllegal);
            bribeOffer = BribeOffer{};
          }
        }
      }

      if (!canScan) {
        // Cancel scans when leaving normal space / docking.
        if (cargoScanActive) cancelCargoScan();
      } else {
        // Don't start a new scan while a bribe window is active.
        if (bribeOffer.active) {
          if (cargoScanActive) cancelCargoScan();
        } else {
          const bool contraband = hasIllegalCargo(jurisdiction);

          if (cargoScanActive) {
            // Scan breaks if you leave range (station comms-range or police scan sphere).
            bool inRange = true;

            if (cargoScanSourceKind == CargoScanSourceKind::Station) {
              // Station scan: player must remain within comms range.
              inRange = false;
              if (currentSystem) {
                for (const auto& st : currentSystem->stations) {
                  if (st.id == cargoScanStationId) {
                    const math::Vec3d stPos = stationPosKm(st, timeDays);
                    const double dist = (ship.positionKm() - stPos).length();
                    inRange = dist < cargoScanRangeKm;
                    break;
                  }
                }
              }
            } else {
              // Police scan: player must remain within police scanner range.
              inRange = false;
              for (const auto& c : contacts) {
                if (c.id == cargoScanContactId && c.role == ContactRole::Police && c.alive) {
                  const double dist = (ship.positionKm() - c.ship.positionKm()).length();
                  inRange = dist < cargoScanRangeKm;
                  break;
                }
              }
            }

            if (!inRange) {
              cancelCargoScan();
              cargoScanCooldownUntilDays = std::max(cargoScanCooldownUntilDays, timeDays + (20.0 / 86400.0));
            } else {
              cargoScanProgressSec += dtReal;
              if (cargoScanProgressSec >= cargoScanDurationSec) {
                // Determine reference prices from the station economy if possible (for fine sizing).
                const sim::Station* priceRef = nullptr;
                if (currentSystem && cargoScanSourceKind == CargoScanSourceKind::Station) {
                  for (const auto& st : currentSystem->stations) {
                    if (st.id == cargoScanStationId) { priceRef = &st; break; }
                  }
                }
                if (!priceRef && currentSystem && !currentSystem->stations.empty()) {
                  priceRef = &currentSystem->stations[0];
                }

                // Build a mid-price reference table for sizing fines in local market terms.
                std::array<double, econ::kCommodityCount> midPriceCr{};
                for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
                  const auto cid = (econ::CommodityId)i;
                  midPriceCr[i] = econ::commodityDef(cid).basePrice;
                }
                if (priceRef) {
                  auto& stEcon = universe.stationEconomy(*priceRef, timeDays);
                  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
                    const auto cid = (econ::CommodityId)i;
                    const auto q = econ::quote(stEcon, priceRef->economyModel, cid, 0.10);
                    midPriceCr[i] = q.mid;
                  }
                }

                const auto scan = sim::scanIllegalCargo(universe.seed(), jurisdiction, cargo, &midPriceCr);
                double illegalValueCr = scan.illegalValueCr;
                std::array<double, econ::kCommodityCount> scannedIllegal = scan.scannedIllegalUnits;

                // Build a short UI detail string (first N line items).
                std::string detail;
                int shown = 0;
                int total = 0;

                for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
                  const double units = scannedIllegal[i];
                  if (units <= 1e-6) continue;

                  ++total;
                  if (shown < 2) {
                    const auto cid = (econ::CommodityId)i;
                    if (!detail.empty()) detail += ", ";
                    detail += econ::commodityDef(cid).name;
                    detail += " x" + std::to_string((int)std::round(units));
                    ++shown;
                  }
                }
                if (total > shown && !detail.empty()) detail += " ...";

                if (illegalValueCr <= 1e-6) {
                  toast(toasts, (cargoScanSourceName.empty() ? "Security" : cargoScanSourceName) + std::string(": scan complete. Cargo is clean."), 2.4);
                  cancelCargoScan();
                  cargoScanCooldownUntilDays = std::max(cargoScanCooldownUntilDays, timeDays + (60.0 / 86400.0));

                  // Also check if you are WANTED here: a successful scan while wanted triggers a demand window.
                  const double bounty = getBounty(jurisdiction);
                  if (bounty > 1e-6) {
                    policeDemand.active = true;
                    policeDemand.factionId = jurisdiction;
                    policeDemand.amountCr = bounty;
                    policeDemand.untilDays = timeDays + (18.0 / 86400.0);
                    policeDemand.sourceName = cargoScanSourceName.empty() ? "Authorities" : cargoScanSourceName;
                    toast(toasts, policeDemand.sourceName + ": outstanding bounty detected. Submit (I) or face enforcement.", 3.2);
                  }
                } else {
                  const sim::LawProfile law = lawForFaction(jurisdiction);
                  const double fineCr = law.contrabandFineCr(illegalValueCr);

                  // Chance that the scanning police contact is corrupt and offers a bribe.
                  bool offerBribe = false;
                  if (cargoScanSourceKind == CargoScanSourceKind::Police) {
                    const double chance = sim::bribeOfferChance(law, getRep(jurisdiction), heat, illegalValueCr);
                    offerBribe = rng.nextUnit() < chance;
                  }

                  if (offerBribe) {
                    bribeOffer.active = true;
                    bribeOffer.factionId = jurisdiction;
                    bribeOffer.scannerContactId = cargoScanContactId;
                    bribeOffer.scannerRangeKm = cargoScanRangeKm;
                    bribeOffer.illegalValueCr = illegalValueCr;
                    bribeOffer.fineCr = fineCr;
                    bribeOffer.amountCr = sim::bribeAmountCrRounded(law, illegalValueCr, 10.0); // nicer UI numbers
                    bribeOffer.startDays = timeDays;
                    bribeOffer.untilDays = timeDays + (14.0 / 86400.0);
                    bribeOffer.sourceName = cargoScanSourceName.empty() ? "Security" : cargoScanSourceName;
                    bribeOffer.detail = detail;
                    bribeOffer.scannedIllegal = scannedIllegal;

                    // Keep some local pressure while you're deciding.
                    policeAlertUntilDays = std::max(policeAlertUntilDays, timeDays + (90.0 / 86400.0));
                    encounterDirector.nextPoliceSpawnDays = std::min(encounterDirector.nextPoliceSpawnDays, timeDays + (6.0 / 86400.0));

                    toast(toasts,
                          bribeOffer.sourceName + ": contraband detected (" + (detail.empty() ? std::string("illegal cargo") : detail)
                            + "). Offer: bribe " + std::to_string((int)std::round(bribeOffer.amountCr)) + " cr (C) to keep cargo, or comply (I).",
                          4.2);
                  } else {
                    enforceContraband(jurisdiction,
                                      cargoScanSourceName.empty() ? "Security" : cargoScanSourceName,
                                      illegalValueCr,
                                      detail,
                                      scannedIllegal);
                  }

                  cancelCargoScan();
                  cargoScanCooldownUntilDays = std::max(cargoScanCooldownUntilDays, timeDays + (90.0 / 86400.0));
                }
              }
            }
          } else if (timeDays >= cargoScanCooldownUntilDays) {
            // Start a scan sometimes when close to a station (or if police are nearby).
            // Contraband makes it more likely.
            const sim::LawProfile law = lawForFaction(jurisdiction);
            const double ratePerSec = sim::cargoScanStartRatePerSec(contraband, law, smuggleHoldMk);

            const double p = 1.0 - std::exp(-ratePerSec * dtReal);

            bool start = rng.nextUnit() < p;

            // Slightly bias scans to happen near stations.
            if (start && currentSystem && !currentSystem->stations.empty()) {
              // nearest station distance
              double nearestKm = 1e99;
              const sim::Station* nearest = nullptr;
              for (const auto& st : currentSystem->stations) {
                const math::Vec3d stPos = stationPosKm(st, timeDays);
                const double dist = (ship.positionKm() - stPos).length();
                if (dist < nearestKm) { nearestKm = dist; nearest = &st; }
              }

              if (nearest && nearestKm < nearest->commsRangeKm * 0.75) {
                // Station scan
                cargoScanActive = true;
                cargoScanSourceKind = CargoScanSourceKind::Station;
                cargoScanStationId = nearest->id;
                cargoScanFactionId = jurisdiction;
                cargoScanSourceName = nearest->name;
                cargoScanProgressSec = 0.0;

                // Duration depends on whether you're carrying contraband and your smuggle hold grade.
                cargoScanDurationSec = sim::cargoScanDurationSecStation(contraband, smuggleHoldMk);

                cargoScanRangeKm = nearest->commsRangeKm;
                toast(toasts, "Incoming cargo scan: " + nearest->name, 2.0);
              }
            }

            // Police scan: if there is a police contact in range, they may scan.
            if (start && !cargoScanActive) {
              for (const auto& c : contacts) {
                if (!c.alive || c.role != ContactRole::Police) continue;
                if (c.factionId != jurisdiction) continue;

                const double dist = (ship.positionKm() - c.ship.positionKm()).length();
                const double scanRange = 42000.0; // km

                if (dist < scanRange) {
                  cargoScanActive = true;
                  cargoScanSourceKind = CargoScanSourceKind::Police;
                  cargoScanContactId = c.id;
                  cargoScanFactionId = jurisdiction;
                  cargoScanSourceName = c.name;
                  cargoScanProgressSec = 0.0;

                  cargoScanDurationSec = sim::cargoScanDurationSecPolice(smuggleHoldMk);
                  cargoScanRangeKm = scanRange;

                  toast(toasts, "Police cargo scan: " + c.name, 2.0);
                  break;
                }
              }
            }
          }
        }
      }
    }

// --- Scanner progress (missions + exploration) ---
if (scanning && !docked && fsdState == FsdState::Idle && supercruiseState == SupercruiseState::Idle) {
  bool valid = false;

  // If the player changes target while scanning, cancel.
  if (target.kind != scanLockedTarget.kind || target.index != scanLockedTarget.index) {
    scanning = false;
    scanProgressSec = 0.0;
    scanLockedId = 0;
    scanLabel.clear();
  } else {
    auto completeScan = [&](const std::string& msg, double toastSec) {
      scanning = false;
      scanProgressSec = 0.0;
      scanLockedId = 0;
      scanLabel.clear();
      scanLockedTarget = Target{};
      toast(toasts, msg, toastSec);
    };

    // --- CONTACT SCAN (primarily for bounty-scan missions) ---
    if (scanLockedTarget.kind == TargetKind::Contact && scanLockedTarget.index < contacts.size()) {
      auto& c = contacts[scanLockedTarget.index];
      if (c.alive && c.id == scanLockedId) {
        const double dist = (c.ship.positionKm() - ship.positionKm()).length();
        if (dist <= scanRangeKm) {
          valid = true;
          scanProgressSec += dtReal;

          if (scanProgressSec >= scanDurationSec) {
	            // Distress victims: use a normal contact scan as a "comms handshake" and
	            // optionally transfer requested cargo for a small reward.
	            if (c.distressVictim && c.distressNeedUnits > 1e-6) {
	              const double kTransferRangeKm = 15000.0;
	              if (dist > kTransferRangeKm) {
	                completeScan("Too far for transfer. Close to within 15,000 km and scan again.", 2.8);
	                continue;
	              }

	              const econ::CommodityId needCid = c.distressNeedCommodity;
	              const auto def = econ::commodityDef(needCid);
	              const double have = cargo[(int)needCid];
	              const double need = std::max(0.0, c.distressNeedUnits);
	              const double moved = std::clamp(std::min(have, need), 0.0, have);

	              if (moved <= 1e-6) {
	                completeScan(std::string("Distress: no ") + def.name + " in cargo. Bring supplies and scan again.", 3.2);
	                continue;
	              }

	              cargo[(int)needCid] = std::max(0.0, cargo[(int)needCid] - moved);
	              c.distressNeedUnits = std::max(0.0, c.distressNeedUnits - moved);
	              const int movedI = std::max(1, (int)std::llround(moved));

	              if (c.distressNeedUnits <= 1e-6) {
	                // Complete rescue: pay out once.
	                const int rewardI = std::max(0, (int)std::llround(std::max(0.0, c.distressRewardCr)));
	                if (rewardI > 0) credits += (double)rewardI;
	
	                if (c.distressPayerFactionId != 0 && c.distressRepReward != 0.0) {
	                  repByFaction[c.distressPayerFactionId] = repByFaction[c.distressPayerFactionId] + c.distressRepReward;
	                }

	                // Mark the parent signal (if any) as completed for UI.
	                for (auto& sig : signals) {
	                  if (sig.id == c.distressSignalId) {
	                    sig.distressCompleted = true;
	                    break;
	                  }
	                }

	                c.distressVictim = false;
	                c.tradeCooldownUntilDays = timeDays; // allow normal trader behavior again
	
	                const std::string repStr = (c.distressPayerFactionId != 0 && c.distressRepReward != 0.0)
	                  ? (" | rep +" + std::to_string((int)std::llround(c.distressRepReward)))
	                  : std::string{};
	
	                completeScan(std::string("Rescue complete: transferred ") + std::to_string(movedI) + " " + def.name +
	                              " | +" + std::to_string(rewardI) + " cr" + repStr,
	                            4.0);
	              } else {
	                const int remainingI = std::max(1, (int)std::llround(c.distressNeedUnits));
	                completeScan(std::string("Aid delivered: transferred ") + std::to_string(movedI) + " " + def.name +
	                              " (remaining " + std::to_string(remainingI) + ").",
	                            3.2);
	              }
	              continue;
	            }

            // Complete bounty-scan missions via the shared mission logic (keeps prototype + tests consistent).
            sim::SaveGame tmp{};
            tmp.credits = credits;
            tmp.missions = missions;

            tmp.reputation.reserve(repByFaction.size());
            for (const auto& kv : repByFaction) {
              tmp.reputation.push_back(sim::FactionReputation{kv.first, kv.second});
            }

            const auto res = sim::tryCompleteBountyScan(tmp, currentSystem->stub.id, c.id, +2.0);
            if (res.completed > 0) {
              credits = tmp.credits;
              missions = std::move(tmp.missions);

              repByFaction.clear();
              for (const auto& r : tmp.reputation) {
                repByFaction[r.factionId] = r.rep;
              }

              toast(toasts, "Mission complete: bounty scan uploaded! +" + std::to_string((int)std::round(res.rewardCr)) + " cr", 3.0);

              scanning = false;
              scanProgressSec = 0.0;
              scanLockedId = 0;
              scanLabel.clear();
              scanLockedTarget = Target{};
            } else {
              completeScan("Scan complete.", 1.6);
            }
          }
        }
      }
    }

    // --- STATION SCAN ---
    if (!valid && scanLockedTarget.kind == TargetKind::Station && scanLockedTarget.index < currentSystem->stations.size()) {
      const auto& st = currentSystem->stations[scanLockedTarget.index];
      if (st.id == scanLockedId) {
        const double dist = (stationPosKm(st, timeDays) - ship.positionKm()).length();
        if (dist <= scanRangeKm) {
          valid = true;
          scanProgressSec += dtReal;

          if (scanProgressSec >= scanDurationSec) {
            const core::u64 key = scanKeyStation(st.id);
            if (scannedKeys.find(key) == scannedKeys.end()) {
              scannedKeys.insert(key);

              const double value = 180.0 + (double)st.type * 40.0;
              explorationDataCr += value;
              completeScan("Station scan logged (+data " + std::to_string((int)value) + " cr).", 2.5);
            } else {
              completeScan("Station already scanned.", 1.8);
            }
          }
        }
      }
    }

    // --- PLANET SCAN ---
    if (!valid && scanLockedTarget.kind == TargetKind::Planet && scanLockedTarget.index < currentSystem->planets.size()) {
      const auto& p = currentSystem->planets[scanLockedTarget.index];
      const math::Vec3d pPosKm = sim::orbitPosition3DAU(p.orbit, timeDays) * kAU_KM;
      const double dist = (pPosKm - ship.positionKm()).length();
      if (dist <= scanRangeKm) {
        valid = true;
        scanProgressSec += dtReal;

        if (scanProgressSec >= scanDurationSec) {
          const core::u64 key = scanKeyPlanet(currentSystem->stub.id, scanLockedTarget.index);
          if (scannedKeys.find(key) == scannedKeys.end()) {
            scannedKeys.insert(key);

            double typeMul = 1.0;
            if (p.type == sim::PlanetType::Ocean) typeMul = 1.25;
            if (p.type == sim::PlanetType::Ice) typeMul = 1.15;
            if (p.type == sim::PlanetType::GasGiant) typeMul = 1.45;

            const double radiusKm = p.radiusEarth * kEARTH_RADIUS_KM;
            const double base = 240.0 + radiusKm * 0.03;
            const double value = base * typeMul + (p.orbit.semiMajorAxisAU * 22.0);
            explorationDataCr += value;

            completeScan("Planet scan logged (+data " + std::to_string((int)value) + " cr).", 2.5);
          } else {
            completeScan("Planet already scanned.", 1.8);
          }
        }
      }
    }

	    // --- STAR SCAN ---
	    if (!valid && scanLockedTarget.kind == TargetKind::Star) {
	      valid = true;
	      scanProgressSec += dtReal;

	      if (scanProgressSec >= scanDurationSec) {
	        const core::u64 key = scanKeyStar(currentSystem->stub.id);
	        if (scannedKeys.find(key) == scannedKeys.end()) {
	          scannedKeys.insert(key);

	          const double value = 220.0 + (double)static_cast<int>(currentSystem->star.cls) * 80.0;
	          explorationDataCr += value;
	          completeScan("Star scan logged (+data " + std::to_string((int)value) + " cr).", 2.2);
	        } else {
	          completeScan("Star already scanned.", 1.8);
	        }
	      }
	    }

	    // --- SIGNAL SCAN (distress / derelicts / resource sites) ---
	    if (!valid && scanLockedTarget.kind == TargetKind::Signal && scanLockedTarget.index < signals.size()) {
	      auto& s = signals[scanLockedTarget.index];
	      if (s.id == scanLockedId) {
	        const double dist = (s.posKm - ship.positionKm()).length();
	        if (dist <= scanRangeKm) {
	          valid = true;
	          scanProgressSec += dtReal;

	          if (scanProgressSec >= scanDurationSec) {
	            const core::u64 key = scanKeySignal(s.id);
	            if (scannedKeys.find(key) == scannedKeys.end()) {
	              scannedKeys.insert(key);

	              double value = 140.0;
	              if (s.type == SignalType::Distress) value = 180.0;
	              if (s.type == SignalType::Derelict) value = 320.0;
	              explorationDataCr += value;

	              // Derelicts sometimes still hold an intact data core you can scoop.
	              if (s.type == SignalType::Derelict) {
	                spawnCargoPod(econ::CommodityId::Electronics, 1.0, s.posKm, {0,0,0}, 0.25);
	              }

	              completeScan(std::string("Signal scan logged (+data ") + std::to_string((int)value) + " cr).", 2.5);

	              // Distress scans can reveal what the call is likely asking for.
	              if (s.type == SignalType::Distress && currentSystem) {
	                if (!s.hasDistressPlan) {
	                  const core::u32 jurisdiction = currentSystem->stub.factionId;
	                  s.hasDistressPlan = true;
	                  s.distress = sim::planDistressEncounter(universe.seed(), currentSystem->stub.id, s.id, timeDays, jurisdiction);
	                }
	                if (s.hasDistressPlan) {
	                  if (s.distress.hasVictim && !s.distressCompleted) {
	                    const auto def = econ::commodityDef(s.distress.needCommodity);
	                    const int needUnits = std::max(1, (int)std::llround(s.distress.needUnits));
	                    toast(toasts,
	                          std::string("Distress analysis: request ") + std::to_string(needUnits) + " " + def.name +
	                            " | reward ~" + std::to_string((int)std::llround(s.distress.rewardCr)) + " cr",
	                          3.6);
	                  }
	                  if (s.distress.ambush) {
	                    toast(toasts, "Warning: hostile activity suspected.", 3.2);
	                  }
	                }
	              }

	              // Resource field scans can reveal composition and richness.
	              if (s.type == SignalType::Resource && s.hasResourcePlan) {
	                const auto& p = s.resource;
	                const auto& a = econ::commodityDef(p.primaryYield);
	                const auto& b = econ::commodityDef(p.secondaryYield);
	                const int pa = (int)std::lround(std::clamp(p.primaryChance, 0.0, 1.0) * 100.0);
	                const int pb = 100 - pa;
	                const int rich = (int)std::lround(std::clamp(p.richness01, 0.0, 1.0) * 100.0);
	                char buf[256];
	                std::snprintf(buf, sizeof(buf),
	                              "Resource analysis: %s | Richness %d%% | %s ~%d%% / %s ~%d%%",
	                              sim::resourceFieldKindName(p.kind),
	                              rich,
	                              a.name,
	                              pa,
	                              b.name,
	                              pb);
	                toast(toasts, buf, 3.8);
	              }
	            } else {
	              completeScan("Signal already scanned.", 1.8);
	            }
	          }
	        }
	      }
	    }

	    // --- ASTEROID PROSPECT ---
	    if (!valid && scanLockedTarget.kind == TargetKind::Asteroid && scanLockedTarget.index < asteroids.size()) {
	      const auto& a = asteroids[scanLockedTarget.index];
	      if (a.id == scanLockedId) {
	        const double dist = (a.posKm - ship.positionKm()).length();
	        if (dist <= scanRangeKm) {
	          valid = true;
	          scanProgressSec += dtReal;

	          if (scanProgressSec >= scanDurationSec) {
	            const core::u64 key = scanKeyAsteroid(a.id);
	            if (scannedKeys.find(key) == scannedKeys.end()) {
	              scannedKeys.insert(key);

	              const double value = 45.0;
	              explorationDataCr += value;
	              const std::string cname = econ::commodityDef(a.yield).name;
	              const int rem = (int)std::lround(std::max(0.0, a.remainingUnits));
	              completeScan(std::string("Asteroid prospected: ") + cname + " | Remaining ~" + std::to_string(rem) + " u (+data " + std::to_string((int)value) + " cr).", 2.8);
	            } else {
	              completeScan("Asteroid already prospected.", 1.8);
	            }
	          }
	        }
	      }
	    }

    // If we successfully scanned something, check "system completion" bonus once.
    if (!scanning) {
      const bool starDone = scannedKeys.find(scanKeyStar(currentSystem->stub.id)) != scannedKeys.end();
      bool planetsDone = true;
      for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
        if (scannedKeys.find(scanKeyPlanet(currentSystem->stub.id, i)) == scannedKeys.end()) {
          planetsDone = false;
          break;
        }
      }
      bool stationsDone = true;
      for (const auto& st : currentSystem->stations) {
        if (scannedKeys.find(scanKeyStation(st.id)) == scannedKeys.end()) {
          stationsDone = false;
          break;
        }
      }

      if (starDone && planetsDone && stationsDone) {
        const core::u64 compKey = scanKeySystemComplete(currentSystem->stub.id);
        if (scannedKeys.find(compKey) == scannedKeys.end()) {
          scannedKeys.insert(compKey);
          const double bonus = 600.0 + 90.0 * (double)currentSystem->planets.size() + 120.0 * (double)currentSystem->stations.size();
          explorationDataCr += bonus;
          const core::u32 lf = currentSystem ? currentSystem->stub.factionId : 0;
          if (lf != 0) addRep(lf, +1.0);
          toast(toasts, "System survey complete! Bonus data +" + std::to_string((int)bonus) + " cr.", 3.0);
        }
      }
    }

    if (!valid) {
      // Cancel silently if the target became invalid or we drifted out of range.
      scanning = false;
      scanProgressSec = 0.0;
      scanLockedId = 0;
      scanLabel.clear();
    }
  }
} else if (!scanning) {
  scanProgressSec = 0.0;
}

	      // --- Salvage / signals / mining objects ---
	      {
	        // Police "submit or fight" window expiry
	        if (policeDemand.active && timeDays > policeDemand.untilDays) {
	          toast(toasts, policeDemand.sourceName + ": no compliance detected. Lethal force authorized!", 3.0);
	          policeDemand.active = false;
	        }

	        // Pirate extortion: satisfy by jettisoning cargo value, or they open fire.
	        if (pirateDemand.active) {
	          const bool normalSpace = (!docked && fsdState == FsdState::Idle && supercruiseState == SupercruiseState::Idle);
	          const std::string src = pirateDemand.leaderName.empty() ? std::string("Pirates") : pirateDemand.leaderName;

	          if (!normalSpace) {
	            pirateDemand = PirateDemand{};
	          } else {
	            bool groupAlive = false;
	            bool groupUnderFire = false;
	            for (const auto& c : contacts) {
	              if (!c.alive) continue;
	              if (c.role != ContactRole::Pirate) continue;
	              if (c.groupId != pirateDemand.groupId) continue;
	              groupAlive = true;
	              if (c.underFireUntilDays > timeDays) groupUnderFire = true;
	            }

	            if (!groupAlive) {
	              pirateDemand = PirateDemand{};
	            } else if (groupUnderFire) {
	              resolvePirateDemand(false, src + ": betrayal! Open fire!");
	            } else if (pirateDemand.requiredValueCr > 0.0 && pirateDemand.deliveredValueCr + 1e-6 >= pirateDemand.requiredValueCr) {
	              resolvePirateDemand(true, src + ": good. We'll take it and leave.");
	            } else if (timeDays > pirateDemand.untilDays) {
	              resolvePirateDemand(false, src + ": time's up. Open fire!");
	            }
	          }
	        }

	        // Drift + despawn floating cargo pods
	        for (auto& pod : floatingCargo) {
	          pod.posKm += pod.velKmS * dtSim;
	          // light damping to keep the field readable
	          pod.velKmS *= std::pow(0.985, dtSim);
	        }
	        for (std::size_t i = 0; i < floatingCargo.size(); /*manual*/) {
	          const bool expired = (floatingCargo[i].expireDay > 0.0 && timeDays > floatingCargo[i].expireDay);
	          const bool empty = (floatingCargo[i].units <= 1e-4);
	          if (expired || empty) {
	            if (target.kind == TargetKind::Cargo) {
	              if (target.index == i) target = {};
	              else if (target.index > i) target.index--;
	            }
	            floatingCargo.erase(floatingCargo.begin() + (std::ptrdiff_t)i);
	            continue;
	          }
	          ++i;
	        }

	        // Despawn expired signal sources (and keep target index consistent)
	        for (std::size_t i = 0; i < signals.size(); /*manual*/) {
	          const bool expired = (signals[i].expireDay > 0.0 && timeDays > signals[i].expireDay);
	          if (expired) {
	            if (target.kind == TargetKind::Signal) {
	              if (target.index == i) target = {};
	              else if (target.index > i) target.index--;
	            }
	            signals.erase(signals.begin() + (std::ptrdiff_t)i);
	            continue;
	          }
	          ++i;
	        }

	        // Ensure mission-specific salvage signals exist in-system (so accepting a mission while docked
	        // still spawns the site without requiring a jump/reseed).
	        if (currentSystem) {
	          for (const auto& m : missions) {
	            if (m.completed || m.failed) continue;
	            if (m.type != sim::MissionType::Salvage) continue;
	            if (m.toSystem != currentSystem->stub.id) continue;
	            if (m.scanned) continue;
	            if (m.targetNpcId == 0) continue;

	            bool exists = false;
	            for (const auto& s : signals) {
	              if (s.id == m.targetNpcId) { exists = true; break; }
	            }
	            if (exists) continue;

	            const sim::Station* baseSt = nullptr;
	            for (const auto& st : currentSystem->stations) {
	              if (st.id == m.toStation) { baseSt = &st; break; }
	            }
	            const math::Vec3d basePos = baseSt ? stationPosKm(*baseSt, timeDays) : ship.positionKm();
	            const double comms = baseSt ? baseSt->commsRangeKm : 120000.0;

	            core::SplitMix64 srng(core::hashCombine((core::u64)currentSystem->stub.id, (core::u64)m.targetNpcId));
	            auto randUnitDet = [&]() -> math::Vec3d {
	              const double x = srng.range(-1.0, 1.0);
	              const double y = srng.range(-1.0, 1.0);
	              const double z = srng.range(-1.0, 1.0);
	              const double len = std::sqrt(std::max(1e-9, x*x + y*y + z*z));
	              return {x/len, y/len, z/len};
	            };

	            const double dist = comms * 1.8 + srng.range(150000.0, 260000.0);

	            SignalSource s{};
	            s.id = m.targetNpcId; // stable mission signal id
	            s.type = SignalType::Derelict;
	            s.posKm = basePos + randUnitDet() * dist;

	            const double ttl = (m.deadlineDay > timeDays) ? std::clamp(m.deadlineDay - timeDays, 0.25, 4.0) : 1.0;
	            s.expireDay = timeDays + ttl;
	            s.resolved = false;
	            signals.push_back(s);
	          }
	        }

	        // Resolve signal sites when you arrive in normal space
	        if (!docked && fsdState == FsdState::Idle && supercruiseState == SupercruiseState::Idle) {
	          const double kResolveRangeKm = 90000.0;
	          for (auto& s : signals) {
	            if (s.resolved) continue;
	            const double distKm = (s.posKm - ship.positionKm()).length();
	            if (distKm > kResolveRangeKm) continue;
	
	            s.resolved = true;

	            // Persist one-shot deterministic signals (e.g., system-entry derelicts) so players
	            // can't farm them by leaving/re-entering a system.
	            if ((s.id & kDeterministicWorldIdBit) && s.type != SignalType::Resource) {
	              resolvedSignalIds.insert(s.id);
	            }
	
	            if (s.type == SignalType::Resource) {
	              if (!s.fieldSpawned) {
	                s.fieldSpawned = true;
	                const int n = 24 + rng.range(0, 16);
	                for (int i = 0; i < n; ++i) {
	                  AsteroidNode a;
	                  a.id = allocWorldId();
	                  a.posKm = s.posKm + randUnit() * rng.range(20000.0, 120000.0);
	                  a.radiusKm = rng.range(1400.0, 4200.0);
	                  a.yield = (rng.nextUnit() < 0.18) ? econ::CommodityId::Metals : econ::CommodityId::Ore;
	                  a.remainingUnits = rng.range(70.0, 220.0) * (a.yield == econ::CommodityId::Metals ? 0.65 : 1.0);
	                  asteroids.push_back(a);
	                }
	              }
	              if (s.hasResourcePlan) {
	                toast(toasts,
	                      std::string("Arrived at ") + sim::resourceFieldKindName(s.resource.kind) + ": asteroid fragments detected.",
	                      3.0);
	              } else {
	                toast(toasts, "Arrived at Resource Site: asteroid fragments detected.", 3.0);
	              }
	            } else if (s.type == SignalType::Derelict) {
	              // Check if this derelict is a Salvage mission site (mission stores the signal id in targetNpcId).
	              sim::Mission* salvageM = nullptr;
	              if (currentSystem) {
	                for (auto& m : missions) {
	                  if (m.completed || m.failed) continue;
	                  if (m.type != sim::MissionType::Salvage) continue;
	                  if (m.toSystem != currentSystem->stub.id) continue;
	                  if (m.targetNpcId != s.id) continue;
	                  salvageM = &m;
	                  break;
	                }
	              }

	              if (salvageM && !salvageM->scanned) {
	                salvageM->scanned = true;

	                const auto def = econ::commodityDef(salvageM->commodity);
	                const int needUnits = std::max(1, (int)std::llround(salvageM->units));

	                toast(toasts,
	                      std::string("Mission site found: recover ") + std::to_string(needUnits) + " " + def.name + " and return to the station.",
	                      3.6);

	                // Guaranteed mission cargo burst (tagged with mission id so it can be highlighted in the scanner).
	                const int missionPods = 2 + rng.range(0, 2);
	                spawnCargoBurst(salvageM->commodity,
	                                (double)needUnits,
	                                s.posKm + randUnit() * rng.range(1500.0, 6000.0),
	                                {0,0,0},
	                                missionPods,
	                                false,
	                                salvageM->id);
	              } else {
	                toast(toasts, "Derelict located. Salvage pods drifting nearby.", 3.0);
	              }

	              // Spawn generic salvage pods once (reduced if this is a mission site so it doesn't overpay).
	              const int genericPods = salvageM ? (1 + rng.range(0, 2)) : (2 + rng.range(0, 3));
	              for (int i = 0; i < genericPods; ++i) {
	                const econ::CommodityId table[] = {econ::CommodityId::Machinery, econ::CommodityId::Electronics, econ::CommodityId::Metals, econ::CommodityId::Luxury};
	                const auto cid = table[rng.range(0, (int)std::size(table) - 1)];
	                const double units = rng.range(2.0, 10.0);
	                spawnCargoPod(cid, units, s.posKm + randUnit() * rng.range(1500.0, 9000.0), {0,0,0}, 0.35);
	              }
	            } else if (s.type == SignalType::Distress) {
	              // Ensure a deterministic plan exists (older save runs / modded spawners may not fill it).
	              if (!s.hasDistressPlan && currentSystem) {
	                const core::u32 jurisdiction = currentSystem->stub.factionId;
	                s.hasDistressPlan = true;
	                s.distress = sim::planDistressEncounter(universe.seed(), currentSystem->stub.id, s.id, timeDays, jurisdiction);
	              }

	              const sim::DistressPlan plan = s.distress;
	              toast(toasts, "Distress beacon acquired. Approach with caution.", 3.0);

	              // Some cargo may be drifting at the site (either spilled supplies or bait).
	              const econ::CommodityId driftCid = plan.hasVictim ? plan.needCommodity : econ::CommodityId::Food;
	              spawnCargoPod(driftCid,
	                            rng.range(2.0, 9.0),
	                            s.posKm + randUnit() * rng.range(1500.0, 8000.0),
	                            {0,0,0},
	                            0.35);

	              // Legit distress calls can spawn a stranded ship that requests supplies.
	              if (plan.hasVictim) {
	                Contact v{};
	                v.id = allocWorldId();
	                v.role = ContactRole::Trader;
	                v.factionId = currentSystem ? currentSystem->stub.factionId : 0;
	                v.name = "Distress Vessel";
	                v.hostileToPlayer = false;
	
	                v.ship = sim::Ship{};
	                v.ship.setPositionKm(s.posKm + randUnit() * rng.range(4000.0, 9000.0));
	                v.ship.setVelocityKmS({0,0,0});
	                v.ship.setOrientation(ship.orientation());
	                v.ship.setAngularVelocityRadS({0,0,0});

	                // A weak hauler archetype; the distress flag disables trader AI until rescued.
	                configureContactLoadout(v,
	                                        ShipHullClass::Hauler,
	                                        /*thrMk=*/1,
	                                        /*shieldMk=*/1,
	                                        /*distMk=*/1,
	                                        /*weapon=*/WeaponType::PulseLaser,
	                                        /*hullMul=*/0.80,
	                                        /*shieldMul=*/0.45,
	                                        /*regenMul=*/0.60,
	                                        /*accelMul=*/0.55,
	                                        /*aiSkill=*/0.35);
	                v.pips = sim::Pips{3,2,1};
	                sim::normalizePips(v.pips);

	                v.distressVictim = true;
	                v.distressSignalId = s.id;
	                v.distressNeedCommodity = plan.needCommodity;
	                v.distressNeedUnits = plan.needUnits;
	                v.distressRewardCr = plan.rewardCr;
	                v.distressRepReward = plan.repReward;
	                v.distressPayerFactionId = plan.payerFactionId;

	                // Keep them from affecting station inventories while disabled.
	                v.tradeUnits = 0.0;
	                v.tradeCooldownUntilDays = timeDays + 9999.0;
	                v.cargoValueCr = 0.0;
	                v.tradeCargoValueCr = 0.0;

	                contacts.push_back(v);
	                s.distressVictimId = v.id;

	                const auto def = econ::commodityDef(plan.needCommodity);
	                const int needUnits = std::max(1, (int)std::llround(plan.needUnits));
	                toast(toasts,
	                      std::string("Distress vessel requests ") + std::to_string(needUnits) + " " + def.name +
	                        " (reward " + std::to_string((int)std::llround(plan.rewardCr)) + " cr). Scan the ship to transfer supplies.",
	                      4.2);
	              }

	              // ...and sometimes it's an ambush.
	              if (plan.ambush && plan.pirateCount > 0) {
	                const int pirates = std::clamp(plan.pirateCount, 1, 4);
	                for (int i = 0; i < pirates; ++i) {
	                  Contact p;
	                  p.id = allocWorldId();
	                  p.role = ContactRole::Pirate;
	                  p.name = "Pirate " + std::to_string((int)contacts.size() + 1);
	                  p.factionId = 0;
	                  p.hostileToPlayer = true;
	                  p.ship = sim::Ship{};
	                  p.ship.setPositionKm(s.posKm + randUnit() * rng.range(80000.0, 140000.0));
	                  p.ship.setVelocityKmS({0,0,0});
	                  p.ship.setOrientation(ship.orientation());
	                  p.ship.setAngularVelocityRadS({0,0,0});
	                  const WeaponType w = (rng.nextUnit() < 0.55) ? WeaponType::Cannon : WeaponType::BeamLaser;
	                  configureContactLoadout(p,
	                                          ShipHullClass::Fighter,
	                                          1,
	                                          1,
	                                          1,
	                                          w,
	                                          /*hullMul=*/0.95,
	                                          /*shieldMul=*/0.67,
	                                          /*regenMul=*/0.63,
	                                          /*accelMul=*/0.70,
	                                          /*aiSkill=*/0.55);
	                  p.pips = sim::Pips{1, 3, 2};
	                  sim::normalizePips(p.pips);
	                  contacts.push_back(p);
	                }
	                toast(toasts, "Ambush! Pirate signatures inbound.", 3.0);
	              }
	            }
	          }
	        }

	        // Cargo scooping (ties combat/mining to trade)
	        if (!docked && fsdState == FsdState::Idle && supercruiseState == SupercruiseState::Idle && cargoScoopDeployed) {
	          const double kScoopRangeKm = 3500.0;
	          const double kMaxRelSpeedKmS = 18.0;
	          const double shipMassKg = cargoMassKg(cargo) + ship.massKg();
	          (void)shipMassKg; // reserved for future handling
	
	          const math::Vec3d shipPos = ship.positionKm();
	          const math::Vec3d shipVel = ship.velocityKmS();
	
	          for (std::size_t i = 0; i < floatingCargo.size(); ++i) {
	            auto& pod = floatingCargo[i];
	            const double distKm = (pod.posKm - shipPos).length();
	            if (distKm > kScoopRangeKm) continue;
	            const double relSpd = (pod.velKmS - shipVel).length();
	            if (relSpd > kMaxRelSpeedKmS) continue;
	
	            const auto def = econ::commodityDef(pod.commodity);
	            const double freeKg = cargoCapacityKg - cargoMassKg(cargo);
	            const double maxUnits = freeKg / std::max(0.001, def.massKg);
	            const double takeUnits = std::min(pod.units, std::floor(maxUnits + 1e-6));
	
	            if (takeUnits >= 1e-4) {
	              cargo[(int)pod.commodity] += takeUnits;
	              pod.units -= takeUnits;
	              toast(toasts, std::string("Scooped ") + def.name + " x" + std::to_string((int)takeUnits), 2.0);
	            } else {
	              if (timeDays > cargoFullToastCooldownUntilDays) {
	                cargoFullToastCooldownUntilDays = timeDays + (4.0 / 86400.0);
	                toast(toasts, "Cargo hold full.", 2.0);
	              }
	            }
	          }
	        }
	      }

	      // --- Heat model (real-time) ---
      {
        sim::ThermalInputs tin{};
        tin.dtReal = dtReal;
        tin.docked = docked;
        tin.supercruiseActive = (supercruiseState == SupercruiseState::Active);
        if (fsdState == FsdState::Charging) tin.fsd = sim::ThermalFsdState::Charging;
        else if (fsdState == FsdState::Jumping) tin.fsd = sim::ThermalFsdState::Jumping;
        else tin.fsd = sim::ThermalFsdState::Idle;
        tin.boostAppliedFrac = boostAppliedFrac;
        tin.heatImpulse = heatImpulse;
        tin.heatCoolRate = playerHeatCoolRate;
        tin.hullMax = playerHullMax;

        const auto tr = sim::stepThermal(heat, tin);
        heat = tr.heat;

        // Apply overheat hull damage.
        if (tr.hullDamage > 0.0) {
          playerHull = std::max(0.0, playerHull - tr.hullDamage);
        }

        // Consume per-frame impulse.
        heatImpulse = 0.0;
      }

      timeDays += dtSim / 86400.0;

      // --- Law/alert decay ---
      // policeHeat is local "security alert" pressure. It should cool off over time,
      // independent of ship thermal heat. Use an exponential decay with a time constant
      // so large timeScale jumps remain stable.
      {
        const double tauSec = 22.0 * 60.0; // ~22 min time constant (tunable)
        const double k = std::exp(-std::max(0.0, dtSim) / std::max(1.0, tauSec));
        policeHeat *= k;
        if (policeHeat < 0.0005) policeHeat = 0.0;
      }


      // --- Background NPC trade traffic (market nudging) ---
      if (currentSystem) {
        sim::simulateNpcTradeTraffic(universe, *currentSystem, timeDays, trafficDayStampBySystem);
      }

      // --- Escort mission runtime (convoy status + ambush events) ---
      if (currentSystem) {
        auto findContactById = [&](core::u64 id) -> Contact* {
          if (id == 0) return nullptr;
          for (auto& c : contacts) {
            if (c.id == id) return &c;
          }
          return nullptr;
        };

        auto findStationById = [&](sim::StationId id) -> const sim::Station* {
          if (id == 0) return nullptr;
          for (const auto& st : currentSystem->stations) {
            if (st.id == id) return &st;
          }
          return nullptr;
        };

        auto runtimeForEscort = [&](core::u64 missionId, core::u64 convoyId) -> EscortRuntime& {
          for (auto& er : escortRuntime) {
            if (er.missionId == missionId) {
              er.convoyId = convoyId;
              return er;
            }
          }
          escortRuntime.push_back(EscortRuntime{});
          escortRuntime.back().missionId = missionId;
          escortRuntime.back().convoyId = convoyId;
          escortRuntime.back().pirateGroupId = (convoyId != 0)
              ? (convoyId ^ 0x5a5a5a5a5a5a5a5aull)
              : std::max<core::u64>(1, rng.nextU64());
          escortRuntime.back().tooFarSec = 0.0;
          escortRuntime.back().ambushSpawned = false;
          escortRuntime.back().nextAmbushDays = timeDays + rng.range(70.0, 150.0) / 86400.0;
          return escortRuntime.back();
        };

        auto repDelta = [&](core::u32 fid, double delta) {
          auto& r = repByFaction[fid];
          r = std::clamp(r + delta, -100.0, 100.0);
        };

        auto cleanupEscort = [&](const sim::Mission& m) {
          // Clear convoy markers / escort follow links so these contacts revert to normal behavior.
          if (m.targetNpcId != 0) {
            for (auto& c : contacts) {
              if (c.id == m.targetNpcId) {
                c.escortConvoy = false;
              }
              if (c.followId == m.targetNpcId) {
                c.followId = 0;
              }
              if (c.attackTargetId == m.targetNpcId) {
                c.attackTargetId = 0;
              }
            }
          }

          // Drop runtime state.
          escortRuntime.erase(std::remove_if(escortRuntime.begin(), escortRuntime.end(),
                                             [&](const EscortRuntime& er) { return er.missionId == m.id; }),
                              escortRuntime.end());
        };

        auto spawnAmbush = [&](const sim::Mission& m, EscortRuntime& er, Contact& convoy) {
          // Keep the fight readable: a single small pack.
          int alivePiratesNow = 0;
          for (const auto& c : contacts) {
            if (c.alive && c.role == ContactRole::Pirate) ++alivePiratesNow;
          }

          if (alivePiratesNow >= 8) {
            er.ambushSpawned = true;
            return;
          }

          const double value = std::max(0.0, convoy.tradeCargoValueCr);
          int count = 2;
          if (rng.nextUnit() < 0.55) ++count;
          if (value > 9000.0 && rng.nextUnit() < 0.45) ++count;
          count = std::clamp(count, 2, 4);

          const core::u64 groupId = (er.pirateGroupId != 0) ? er.pirateGroupId : std::max<core::u64>(1, rng.nextU64());

          core::u64 leaderWorldId = 0;
          for (int i = 0; i < count; ++i) {
            Contact p{};
            p.alive = true;
            p.id = allocWorldId();
            if (i == 0) leaderWorldId = p.id;

            p.role = ContactRole::Pirate;
            p.groupId = groupId;
            p.leaderId = (i == 0) ? 0 : leaderWorldId;
            p.hostileToPlayer = true;
            p.attackTargetId = convoy.id;

            // Slightly tougher leader so the pack doesn't melt instantly.
            const double statMul = (i == 0) ? 1.22 : (0.90 + 0.18 * rng.nextUnit());
            const bool leader = (i == 0);
            const int tMk = leader ? 2 : 1;
            const int dMk = leader ? 2 : 1;

            const double r = rng.nextUnit();
            WeaponType w = WeaponType::Cannon;
            if (leader) {
              w = (r < 0.55) ? WeaponType::Railgun : WeaponType::Cannon;
            } else {
              w = (r < 0.55) ? WeaponType::Cannon : WeaponType::BeamLaser;
            }

            configureContactLoadout(p,
                                    ShipHullClass::Fighter,
                                    tMk,
                                    /*sMk=*/1,
                                    dMk,
                                    w,
                                    /*hullMul=*/0.95 * statMul,
                                    /*shieldMul=*/0.67 * statMul,
                                    /*regenMul=*/leader ? 0.70 : 0.63,
                                    /*accelMul=*/leader ? 0.72 : 0.70,
                                    /*aiSkill=*/leader ? 0.68 : 0.55);
            p.pips = leader ? sim::Pips{1, 3, 2} : sim::Pips{1, 3, 2};
            sim::normalizePips(p.pips);

            p.name = (i == 0) ? "Raid Leader" : "Raider";

            if (supercruiseState == SupercruiseState::Active) {
              // Attach as a shadow so it stays in the encounter cluster.
              p.supercruiseShadow = true;
              p.shadowOffsetLocalKm = {
                  rng.range(-28000.0, 28000.0),
                  rng.range(-9000.0, 9000.0),
                  -rng.range(200000.0, 260000.0),
              };

              const math::Quatd q = ship.orientation();
              const math::Vec3d worldOff = q * p.shadowOffsetLocalKm;
              p.ship.setPositionKm(ship.positionKm() + worldOff);
              p.ship.setVelocityKmS(ship.velocityKmS());
            } else {
              const math::Vec3d basePos = convoy.ship.positionKm();
              const math::Vec3d baseVel = convoy.ship.velocityKmS();
              const math::Vec3d dir = randUnit();
              p.ship.setPositionKm(basePos + dir * rng.range(45000.0, 82000.0));
              p.ship.setVelocityKmS(baseVel + randUnit() * rng.range(0.0, 1.2));
            }

            contacts.push_back(std::move(p));
          }

          er.ambushSpawned = true;
          toast(toasts, "Ambush! Pirates attacking the convoy.", 3.5);
        };

        for (auto& m : missions) {
          if (m.type != sim::MissionType::Escort) continue;

          // Cleanup once a mission resolves.
          if (m.failed || m.completed) {
            cleanupEscort(m);
            continue;
          }

          // If the player leaves the system after starting the escort, the convoy is abandoned.
          if (m.leg >= 1 && currentSystem->stub.id != m.fromSystem) {
            m.failed = true;
            repDelta(m.factionId, -4.0);
            toast(toasts, "Mission failed: you abandoned the convoy.", 3.5);
            if (trackedMissionId == m.id) trackedMissionId = 0;
            cleanupEscort(m);
            continue;
          }

          if (m.leg < 1) continue;  // not started yet

          EscortRuntime& er = runtimeForEscort(m.id, m.targetNpcId);

          Contact* convoy = findContactById(m.targetNpcId);
          if (!convoy || !convoy->alive || convoy->hull <= 0.0) {
            m.failed = true;
            repDelta(m.factionId, -4.0);
            toast(toasts, "Mission failed: convoy destroyed.", 3.5);
            if (trackedMissionId == m.id) trackedMissionId = 0;
            cleanupEscort(m);
            continue;
          }

          // Arrival: once the convoy is in range of the destination station, it "checks in".
          if (!m.scanned) {
            if (const sim::Station* destSt = findStationById(m.toStation)) {
              const math::Vec3d destPos = stationPosKm(*destSt, timeDays);
              const double d = (convoy->ship.positionKm() - destPos).length();
              if (d < destSt->commsRangeKm * 0.70) {
                m.scanned = true;
                toast(toasts, "Convoy arrived. Dock to claim your reward.", 3.5);
              }
            }
          }

          // Too-far failure: you must stay within a reasonable range (ignored during supercruise).
          if (!m.scanned) {
            if (supercruiseState == SupercruiseState::Active) {
              er.tooFarSec = 0.0;
            } else {
              const double distKm = (ship.positionKm() - convoy->ship.positionKm()).length();
              if (distKm > 180000.0) {
                er.tooFarSec += dtSim;
              } else {
                er.tooFarSec = 0.0;
              }

              if (er.tooFarSec > 35.0) {
                m.failed = true;
                repDelta(m.factionId, -4.0);
                toast(toasts, "Mission failed: you lost the convoy.", 3.5);
                if (trackedMissionId == m.id) trackedMissionId = 0;
                cleanupEscort(m);
                continue;
              }
            }
          }

          // Ambush event.
          if (!m.scanned && !er.ambushSpawned && timeDays >= er.nextAmbushDays) {
            spawnAmbush(m, er, *convoy);
          }
        }
      }

      // --- Mission deadlines / docked completion ---
      // Delegate state transitions to the shared headless module so gameplay + tooling/tests stay consistent.
      if (!missions.empty()) {
        sim::SaveGame tmp{};
        tmp.timeDays = timeDays;
        tmp.credits = credits;
        tmp.cargo = cargo;
        tmp.cargoCapacityKg = cargoCapacityKg;
        tmp.passengerSeats = passengerSeats;
        tmp.missions = missions;

        // Rep map -> SaveGame vector
        tmp.reputation.clear();
        tmp.reputation.reserve(repByFaction.size());
        for (const auto& kv : repByFaction) {
          tmp.reputation.push_back(sim::FactionReputation{kv.first, kv.second});
        }

        const auto prevMissions = tmp.missions;

        // Deadlines
        (void)sim::tickMissionDeadlines(tmp, timeDays, -4.0);

        // Dock completion
        if (docked && dockedStationId != 0 && currentSystem) {
          const sim::Station* dockSt = nullptr;
          for (const auto& st : currentSystem->stations) {
            if (st.id == dockedStationId) { dockSt = &st; break; }
          }
          if (dockSt) {
            (void)sim::tryCompleteMissionsAtDock(universe, *currentSystem, *dockSt, timeDays, tmp, +2.0);
          }
        }

        // Emit toasts for state changes (small vector, so O(n^2) is fine).
        auto findPrev = [&](core::u64 id) -> const sim::Mission* {
          for (const auto& pm : prevMissions) if (pm.id == id) return &pm;
          return nullptr;
        };

        for (const auto& m : tmp.missions) {
          const sim::Mission* pm = findPrev(m.id);
          if (!pm) continue;

          if (!pm->failed && m.failed) {
            toast(toasts, "Mission failed: deadline missed.", 3.0);
            if (trackedMissionId == m.id) trackedMissionId = 0;
          }

          if (!pm->completed && m.completed) {
            if (m.type == sim::MissionType::Courier) {
              toast(toasts, "Mission complete: courier delivery. +" + std::to_string((int)m.reward) + " cr", 3.0);
            } else if (m.type == sim::MissionType::Escort) {
              toast(toasts, "Mission complete: convoy escorted. +" + std::to_string((int)m.reward) + " cr", 3.0);
            } else if (m.type == sim::MissionType::Passenger) {
              toast(toasts, "Mission complete: passengers delivered. +" + std::to_string((int)m.reward) + " cr", 3.0);
            } else if (m.type == sim::MissionType::Smuggle) {
              toast(toasts, "Mission complete: contraband delivered. +" + std::to_string((int)m.reward) + " cr", 3.0);
            } else if (m.type == sim::MissionType::MultiDelivery) {
              toast(toasts, "Mission complete: multi-hop delivery. +" + std::to_string((int)m.reward) + " cr", 3.0);
            } else if (m.type == sim::MissionType::Delivery) {
              toast(toasts, "Mission complete: delivery. +" + std::to_string((int)m.reward) + " cr", 3.0);
            } else if (m.type == sim::MissionType::Salvage) {
              toast(toasts, "Mission complete: salvage recovered. +" + std::to_string((int)m.reward) + " cr", 3.0);
            }
            if (trackedMissionId == m.id) trackedMissionId = 0;
          }

          if (m.type == sim::MissionType::MultiDelivery && pm->leg == 0 && m.leg == 1 && !m.completed && !m.failed) {
            toast(toasts, "Multi-hop delivery: leg 1/2 complete.", 2.5);

            // If the player is tracking this mission, optionally auto-plot the route to the final destination.
            if (missionTrackerAutoPlotNextLeg && trackedMissionId == m.id) {
              if (plotRouteToSystem(m.toSystem, false)) {
                pendingArrivalTargetStationId = m.toStation;
                if (currentSystem && m.toSystem == currentSystem->stub.id && m.toStation != 0) tryTargetStationById(m.toStation);
                toast(toasts, "Route updated: final destination plotted (leg 2/2).", 2.6);
              }
            }
          }
        }

        // Write back state
        credits = tmp.credits;
        cargo = tmp.cargo;
        missions = std::move(tmp.missions);

        repByFaction.clear();
        for (const auto& r : tmp.reputation) {
          repByFaction[r.factionId] = r.rep;
        }
      }

      // Projectiles + missiles update (ballistic slugs + guided explosives)
      {
        // Build a light-weight target list for collision tests.
        // (We include the player so NPC projectiles can hit later; player shots ignore it.)
        std::vector<sim::SphereTarget> projTargets;
        projTargets.reserve(1 + contacts.size() + asteroids.size());

        sim::SphereTarget playerT{};
        playerT.kind = sim::CombatTargetKind::Player;
        playerT.index = 0;
        playerT.id = 0;
        playerT.centerKm = ship.positionKm();
        playerT.velKmS = ship.velocityKmS();
        playerT.radiusKm = 900.0;
        projTargets.push_back(playerT);

        for (std::size_t i = 0; i < contacts.size(); ++i) {
          const auto& c = contacts[i];
          if (!c.alive) continue;
          sim::SphereTarget t{};
          t.kind = sim::CombatTargetKind::Ship;
          t.index = i;
          t.id = c.id;
          t.centerKm = c.ship.positionKm();
          t.velKmS = c.ship.velocityKmS();
          t.radiusKm = 900.0;
          projTargets.push_back(t);
        }

        for (std::size_t i = 0; i < asteroids.size(); ++i) {
          const auto& a = asteroids[i];
          sim::SphereTarget t{};
          t.kind = sim::CombatTargetKind::Asteroid;
          t.index = i;
          t.id = a.id;
          t.centerKm = a.posKm;
          t.velKmS = {0,0,0};
          t.radiusKm = a.radiusKm;
          projTargets.push_back(t);
        }

        std::vector<sim::ProjectileHit> hits;
        sim::stepProjectiles(
          projectiles,
          dtSim,
          projTargets.empty() ? nullptr : projTargets.data(),
          projTargets.size(),
          hits
        );

        for (const auto& h : hits) {
          if (h.fromPlayer) {
            if (h.kind == sim::CombatTargetKind::Ship) {
              if (h.targetIndex < contacts.size()) {
                playerDamageContact((int)h.targetIndex, h.dmg);
              }
            } else if (h.kind == sim::CombatTargetKind::Asteroid) {
              // Projectiles don't mine, but they do spark (small feedback loop).
              if (vfxParticlesEnabled && vfxImpactsEnabled) {
                math::Vec3d n = h.pointKm - ship.positionKm();
                if (n.lengthSq() < 1e-12) n = ship.forward();
                n = n.normalized();
                const double energy = std::clamp((h.dmg / 18.0) * (double)vfxParticleIntensity, 0.10, 1.6);
                particles.spawnSparks(toRenderU(h.pointKm), n, toRenderU(ship.velocityKmS()), energy);
              }
            }
          } else {
            // NPC projectiles (not used heavily yet, but supported by the core module).
            if (!docked && h.kind == sim::CombatTargetKind::Player) {
              applyDamage(h.dmg, playerShield, playerHull);

              if (vfxParticlesEnabled && vfxImpactsEnabled) {
                // Approximate impact normal as projectile->player direction.
                math::Vec3d n = ship.positionKm() - h.pointKm;
                if (n.lengthSq() < 1e-12) n = ship.forward();
                n = n.normalized();
                const double energy = std::clamp((h.dmg / 16.0) * (double)vfxParticleIntensity, 0.12, 2.2);
                particles.spawnSparks(toRenderU(ship.positionKm()), n, toRenderU(ship.velocityKmS()), energy);
              }
            } else if (h.kind == sim::CombatTargetKind::Ship) {
              if (h.targetIndex < contacts.size()) {
                Contact& victim = contacts[h.targetIndex];
                if (victim.alive && victim.hull > 0.0) {
                  applyDamage(h.dmg, victim.shield, victim.hull);
                  victim.underFireUntilDays = timeDays + (14.0 / 86400.0);

                  if (victim.hull <= 0.0) {
                    victim.alive = false;

                    // Small reward if police clean up pirates (keeps the world feeling alive).
                    Contact* shooter = nullptr;
                    for (auto& s : contacts) {
                      if (s.id == h.shooterId) {
                        shooter = &s;
                        break;
                      }
                    }
                    if (shooter && shooter->role == ContactRole::Police && victim.role == ContactRole::Pirate) {
                      credits += 180.0;
                      toast(toasts, "Police destroyed a pirate. +180 cr", 2.2);
                    }
                  }
                }
              }
            }
          }
        }

        projectiles.erase(
          std::remove_if(projectiles.begin(), projectiles.end(), [](const sim::Projectile& p) { return p.ttlSimSec <= 0.0; }),
          projectiles.end());

        // Missiles update (homing guided explosives)
        {
          std::vector<sim::MissileDetonation> detonations;
          std::vector<sim::MissileHit> mHits;

          sim::stepMissiles(
            missiles,
            dtSim,
            projTargets.empty() ? nullptr : projTargets.data(),
            projTargets.size(),
            detonations,
            mHits
          );

          // Missile trails (cheap thruster emit)
          if (vfxParticlesEnabled && vfxThrustersEnabled) {
            const double intensity = std::clamp(0.55 * (double)vfxParticleIntensity, 0.20, 3.0);
            for (const auto& m : missiles) {
              if (m.ttlSimSec <= 0.0) continue;
              math::Vec3d dir = m.velKmS;
              if (dir.lengthSq() < 1e-12) continue;
              dir = dir.normalized();
              particles.emitThruster(toRenderU(m.posKm), -dir, intensity, dtReal, /*boost=*/true);
            }
          }

          // Throttled incoming missile warning.
          if (!docked && playerHull > 0.0 && timeDays > incomingMissileToastCooldownUntilDays) {
            bool incoming = false;
            double closestKm = 1e30;
            for (const auto& m : missiles) {
              if (m.ttlSimSec <= 0.0) continue;
              if (m.fromPlayer) continue;
              if (!m.hasTarget) continue;
              if (m.targetKind != sim::CombatTargetKind::Player) continue;
              incoming = true;
              closestKm = std::min(closestKm, (m.posKm - ship.positionKm()).length());
            }
            if (incoming && closestKm < 350000.0) {
              toast(toasts, "INCOMING MISSILE", 1.4);
              incomingMissileToastCooldownUntilDays = timeDays + (4.0 / 86400.0);
            }
          }

          // Detonation VFX (one per detonation)
          if (vfxParticlesEnabled && vfxExplosionsEnabled) {
            for (const auto& det : detonations) {
              const double eBase = std::clamp((det.baseDmg / 28.0) * (double)vfxParticleIntensity, 0.6, 4.2);
              particles.spawnExplosion(toRenderU(det.pointKm), {0,0,0}, eBase);
            }
          }

          // Apply splash hits.
          for (const auto& h : mHits) {
            if (h.fromPlayer) {
              if (h.kind == sim::CombatTargetKind::Ship) {
                if (h.targetIndex < contacts.size()) {
                  playerDamageContact((int)h.targetIndex, h.dmg);
                }
              }
              // Missiles don't mine asteroids.
            } else {
              if (!docked && h.kind == sim::CombatTargetKind::Player) {
                applyDamage(h.dmg, playerShield, playerHull);
              } else if (h.kind == sim::CombatTargetKind::Ship) {
                if (h.targetIndex < contacts.size()) {
                  Contact& victim = contacts[h.targetIndex];
                  if (victim.alive && victim.hull > 0.0) {
                    applyDamage(h.dmg, victim.shield, victim.hull);
                    victim.underFireUntilDays = timeDays + (16.0 / 86400.0);

                    if (victim.hull <= 0.0) {
                      victim.alive = false;

                      // Small reward if police clean up pirates.
                      Contact* shooter = nullptr;
                      for (auto& s : contacts) {
                        if (s.id == h.shooterId) { shooter = &s; break; }
                      }
                      if (shooter && shooter->role == ContactRole::Police && victim.role == ContactRole::Pirate) {
                        credits += 200.0;
                        toast(toasts, "Police destroyed a pirate. +200 cr", 2.2);
                      }
                    }
                  }
                }
              }
            }
          }

          missiles.erase(
            std::remove_if(missiles.begin(), missiles.end(), [](const sim::Missile& m) { return m.ttlSimSec <= 0.0; }),
            missiles.end());
        }
      }


      weaponPrimaryCooldown = std::max(0.0, weaponPrimaryCooldown - dtSim);
      weaponSecondaryCooldown = std::max(0.0, weaponSecondaryCooldown - dtSim);

      // Shield regen (slow)
      if (!paused && playerHull > 0.0 && playerShield < playerShieldMax) {
        const double regenMul = sim::shieldRegenMultiplierFromPips(distributorPips.sys);
        const double desiredRegen = playerShieldRegenPerSimMin * regenMul * (dtSim / 60.0);

        const double costPerPt = std::max(0.0, distributorCfg.shieldRegenCostPerPoint);
        const double affordableRegen = (costPerPt > 1e-12) ? (distributorState.sys / costPerPt) : desiredRegen;
        const double actualRegen = std::max(0.0, std::min(desiredRegen, affordableRegen));

        if (costPerPt > 1e-12) {
          distributorState.sys = std::max(0.0, distributorState.sys - actualRegen * costPerPt);
        }
        playerShield = std::min(playerShieldMax, playerShield + actualRegen);
      }

    }

    // Death / respawn
    if (playerHull <= 0.0) {
      toast(toasts, "Ship destroyed! Respawning (lost cargo, -10% credits).", 4.0);

      if (vfxParticlesEnabled && vfxExplosionsEnabled) {
        const double eBase = 2.2 * (double)vfxParticleIntensity;
        particles.spawnExplosion(toRenderU(ship.positionKm()), toRenderU(ship.velocityKmS()), eBase);
      }

      playerHull = playerHullMax;
      playerShield = playerShieldMax * 0.60;
      credits *= 0.90;
      cargo.fill(0.0);
      fuel = fuelMax;
      supercruiseState = SupercruiseState::Idle;
      supercruiseChargeRemainingSec = 0.0;
      supercruiseCooldownRemainingSec = 0.0;
      supercruiseDropRequested = false;
      interdiction = sim::InterdictionState{};
      interdictionSubmitRequested = false;
      interdictionPirateName.clear();
      interdictionPirateStrength = 1.0;
      scanning = false;
      scanProgressSec = 0.0;
      fsdState = FsdState::Idle;
      fsdTargetSystem = 0;
      fsdChargeRemainingSec = 0.0;
      fsdTravelRemainingSec = 0.0;
      fsdTravelTotalSec = 0.0;
      navAutoRun = false;
      incomingMissileToastCooldownUntilDays = 0.0;
      beams.clear();
      projectiles.clear();
      missiles.clear();
      contacts.clear();
      docked = true;
      if (!currentSystem->stations.empty()) {
        dockedStationId = currentSystem->stations.front().id;
        selectedStationIndex = 0;
      } else {
        dockedStationId = 0;
      }
      // Will snap to dock position next tick.
    }

    // Update beam TTL
    for (auto& b : beams) b.ttl -= dtReal;
    beams.erase(std::remove_if(beams.begin(), beams.end(), [](const Beam& b){ return b.ttl <= 0.0; }), beams.end());

    // Toast TTL
    for (auto& t : toasts) t.ttl -= dtReal;
    toasts.erase(std::remove_if(toasts.begin(), toasts.end(), [](const ToastMsg& t){ return t.ttl <= 0.0; }), toasts.end());

    // ---- Camera follow (third-person) ----
    render::Camera cam;
    int w = 1280, h = 720;
    SDL_GetWindowSize(window, &w, &h);
    const double aspect = (h > 0) ? (double)w / (double)h : 16.0/9.0;

    cam.setPerspective(math::degToRad(60.0), aspect, 0.01, 20000.0);

    const math::Vec3d shipPosU = toRenderU(ship.positionKm());
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
    // World lighting: treat the main star as a point light at the origin.
    // UI preview scenes can override this temporarily.
    meshRenderer.setLightPos(0.0f, 0.0f, 0.0f);
    atmosphereRenderer.setViewProj(viewF, projF);
    lineRenderer.setViewProj(viewF, projF);
    pointRenderer.setViewProj(viewF, projF);

    // Atmosphere needs the camera position in world/render units for the rim factor.
    {
      const math::Vec3d cp = cam.position();
      atmosphereRenderer.setCameraPos((float)cp.x, (float)cp.y, (float)cp.z);
    }

    // Update VFX render buffers.
    if (vfxStarfieldEnabled) {
      if ((int)starfield.starCount() != vfxStarCount) {
        starfield.regenerate(seed ^ 0xC5A0F1A9u, vfxStarCount);
      }
      if (std::abs(starfield.radius() - vfxStarRadiusU) > 1e-3) {
        starfield.setRadius(vfxStarRadiusU);
      }
      starfield.update(cam.position(), timeRealSec);
    }

    if (vfxNebulaEnabled) {
      const bool regen = ((int)nebula.puffCount() != vfxNebulaPuffCount)
                      || (vfxNebulaVariant != vfxNebulaVariantLast)
                      || (std::abs(vfxNebulaBandPower - vfxNebulaBandPowerLast) > 1e-3f);
      if (regen) {
        const core::u64 nSeed = core::hashCombine(seed ^ 0xBADC0FFEE0DDF00DULL, (core::u64)vfxNebulaVariant);
        nebula.regenerate(nSeed, vfxNebulaPuffCount, vfxNebulaBandPower);
        vfxNebulaVariantLast = vfxNebulaVariant;
        vfxNebulaBandPowerLast = vfxNebulaBandPower;
      }

      render::NebulaField::Settings ns{};
      ns.innerRadiusU = vfxNebulaInnerRadiusU;
      ns.outerRadiusU = vfxNebulaOuterRadiusU;
      ns.parallax = vfxNebulaParallax;
      ns.intensity = vfxNebulaIntensity;
      ns.opacity = vfxNebulaOpacity;
      ns.sizeMinPx = vfxNebulaSizeMinPx;
      ns.sizeMaxPx = vfxNebulaSizeMaxPx;
      ns.turbulence = vfxNebulaTurbulence;
      ns.turbulenceSpeed = vfxNebulaTurbulenceSpeed;

      nebula.update(cam.position(), timeRealSec, ns);
    }

    if (vfxParticlesEnabled) {
      particles.buildPoints(particleVerts);
    } else {
      particleVerts.clear();
    }

    // ---- Build instances (star + planets) ----
    // We keep a fallback 'spheres' vector (checker-textured) and also track per-body
    // procedural surface seeds for the optional procedural planet textures.
    std::vector<render::InstanceData> spheres;
    spheres.reserve(1 + currentSystem->planets.size());

    struct PlanetSurfaceDraw {
      render::InstanceData inst;
      render::SurfaceKind surface;
      core::u64 surfaceSeed;
    };
    std::vector<PlanetSurfaceDraw> planetSurfaceDraws;
    planetSurfaceDraws.reserve(currentSystem->planets.size());

    struct PlanetCloudDraw {
      render::InstanceData inst;
      core::u64 cloudSeed;
      float alphaMul;
    };
    std::vector<PlanetCloudDraw> planetCloudDraws;
    planetCloudDraws.reserve(currentSystem->planets.size());

    std::vector<render::InstanceData> planetAtmoDraws;
    planetAtmoDraws.reserve(currentSystem->planets.size());

    // Star at origin
    float starR = 1.0f, starG = 0.95f, starB = 0.75f;
    starClassRgb(currentSystem->star.cls, starR, starG, starB);
    const core::u64 starSurfaceSeed = core::hashCombine(
        core::hashCombine(core::fnv1a64("star_surface"), (core::u64)currentSystem->stub.seed),
        (core::u64)(int)currentSystem->star.cls);
    render::InstanceData starSphereInst{};
    {
      const double starRadiusKm = currentSystem->star.radiusSol * kSOLAR_RADIUS_KM;
      const float starScale = (float)std::max(0.8, (starRadiusKm / kRENDER_UNIT_KM) * 3.0);
      // HDR-friendly multiplier to make the star drive bloom.
      starSphereInst = makeInstUniform({0,0,0}, starScale,
                                      starR * worldStarIntensity,
                                      starG * worldStarIntensity,
                                      starB * worldStarIntensity);
      spheres.push_back(starSphereInst);
    }

    // Planets
    for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
      const auto& p = currentSystem->planets[i];

      const math::Vec3d posAU = sim::orbitPosition3DAU(p.orbit, timeDays);
      const math::Vec3d posKm = posAU * kAU_KM;
      const math::Vec3d posU = toRenderU(posKm);

      const double radiusKm = p.radiusEarth * kEARTH_RADIUS_KM;
      const float scale = (float)std::max(0.25, (radiusKm / kRENDER_UNIT_KM) * 200.0);

      // Fallback color palette by type (checker texture multiplied by this).
      float cr=0.6f, cg=0.6f, cb=0.6f;
      switch (p.type) {
        case sim::PlanetType::Rocky: cr=0.6f; cg=0.55f; cb=0.5f; break;
        case sim::PlanetType::Desert: cr=0.8f; cg=0.7f; cb=0.35f; break;
        case sim::PlanetType::Ocean: cr=0.25f; cg=0.45f; cb=0.85f; break;
        case sim::PlanetType::Ice: cr=0.7f; cg=0.85f; cb=0.95f; break;
        case sim::PlanetType::GasGiant: cr=0.7f; cg=0.55f; cb=0.35f; break;
        default: break;
      }

      spheres.push_back(makeInstUniform(posU, scale, cr,cg,cb));

      const core::u64 pSurfaceSeed = core::hashCombine(
          core::hashCombine(core::fnv1a64("planet_surface"), (core::u64)currentSystem->stub.seed),
          (core::u64)i);
      planetSurfaceDraws.push_back({makeInstUniform(posU, scale, 1.0f, 1.0f, 1.0f),
                                    planetSurfaceKind(p.type),
                                    pSurfaceSeed});

      // ---- Secondary layers: cloud shell + atmosphere rim ----
      // These are purely visual and intentionally lightweight.
      const float var = (float)((pSurfaceSeed >> 24) & 0xFFull) / 255.0f;
      float cloudness = 0.0f;
      float atmoStrength = 0.0f;
      float atmoR = 0.45f, atmoG = 0.70f, atmoB = 1.00f; // default Rayleigh-ish

      switch (p.type) {
        case sim::PlanetType::Rocky:
          cloudness = 0.15f;
          atmoStrength = 0.35f;
          atmoR = 0.45f; atmoG = 0.70f; atmoB = 1.00f;
          break;
        case sim::PlanetType::Desert:
          cloudness = 0.08f;
          atmoStrength = 0.25f;
          atmoR = 0.90f; atmoG = 0.78f; atmoB = 0.55f;
          break;
        case sim::PlanetType::Ocean:
          cloudness = 0.35f;
          atmoStrength = 0.55f;
          atmoR = 0.38f; atmoG = 0.67f; atmoB = 1.00f;
          break;
        case sim::PlanetType::Ice:
          cloudness = 0.20f;
          atmoStrength = 0.28f;
          atmoR = 0.65f; atmoG = 0.85f; atmoB = 1.00f;
          break;
        case sim::PlanetType::GasGiant:
          cloudness = 0.65f;
          atmoStrength = 0.75f;
          if (((pSurfaceSeed >> 5) & 1ull) != 0ull) {
            atmoR = 1.00f; atmoG = 0.65f; atmoB = 0.35f;
          } else {
            atmoR = 0.55f; atmoG = 0.75f; atmoB = 1.00f;
          }
          break;
        default: break;
      }

      cloudness *= (0.85f + 0.30f * var);
      atmoStrength *= (0.85f + 0.25f * var);

      if (worldAtmoTintWithStar) {
        const float tr = 0.35f + 0.65f * starR;
        const float tg = 0.35f + 0.65f * starG;
        const float tb = 0.35f + 0.65f * starB;
        atmoR *= tr;
        atmoG *= tg;
        atmoB *= tb;
      }

      if (worldUseProceduralSurfaces && worldCloudsEnabled && cloudness > 1e-4f) {
        core::u64 cloudSeed = core::hashCombine(pSurfaceSeed, core::fnv1a64("clouds"));
        cloudSeed = core::hashCombine(cloudSeed, (core::u64)(int)p.type);

        const double cloudScale = (double)scale * (double)worldCloudShellScale;
        const double phase = (double)((cloudSeed >> 10) & 0xFFFFull) / 65535.0 * (2.0 * math::kPi);
        const double speed = math::degToRad((double)worldCloudSpinDegPerSec) * (0.4 + 1.2 * (double)var);
        const double ang = phase + timeRealSec * speed;
        const math::Quatd q = math::Quatd::fromAxisAngle({0,1,0}, ang);

        PlanetCloudDraw cd{};
        cd.inst = makeInst(posU, {cloudScale, cloudScale, cloudScale}, q, 1.0f, 1.0f, 1.0f);
        cd.cloudSeed = cloudSeed;
        cd.alphaMul = worldCloudOpacity * cloudness;
        planetCloudDraws.push_back(cd);
      }

      if (worldAtmospheresEnabled && atmoStrength > 1e-4f) {
        const double atmoScale = (double)scale * (double)worldAtmoShellScale;
        planetAtmoDraws.push_back(makeInstUniform(posU, atmoScale,
                                                  atmoR * atmoStrength,
                                                  atmoG * atmoStrength,
                                                  atmoB * atmoStrength));
      }
    }

    // Orbit lines (planets)
    std::vector<render::LineVertex> lines;
    lines.reserve(currentSystem->planets.size() * 128 + currentSystem->stations.size() * 64 + beams.size() * 2);

    for (const auto& p : currentSystem->planets) {
      const int seg = 96;
      math::Vec3d prev{};
      for (int s = 0; s <= seg; ++s) {
        const double t = (double)s / (double)seg * p.orbit.periodDays;
        const math::Vec3d posAU = sim::orbitPosition3DAU(p.orbit, t);
        const math::Vec3d posU = toRenderU(posAU * kAU_KM);

        if (s > 0) {
          lines.push_back({(float)prev.x,(float)prev.y,(float)prev.z, 0.22f,0.22f,0.25f});
          lines.push_back({(float)posU.x,(float)posU.y,(float)posU.z, 0.22f,0.22f,0.25f});
        }
        prev = posU;
      }
    }

    // Station orbit lines + slot axis lines
    for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
      const auto& st = currentSystem->stations[i];

      // orbit (fewer segs, stations are secondary)
      const int seg = 48;
      math::Vec3d prev{};
      for (int s = 0; s <= seg; ++s) {
        const double t = (double)s / (double)seg * st.orbit.periodDays;
        const math::Vec3d posU = toRenderU(sim::orbitPosition3DAU(st.orbit, t) * kAU_KM);
        if (s > 0) {
          lines.push_back({(float)prev.x,(float)prev.y,(float)prev.z, 0.16f,0.18f,0.22f});
          lines.push_back({(float)posU.x,(float)posU.y,(float)posU.z, 0.16f,0.18f,0.22f});
        }
        prev = posU;
      }

      // docking corridor guidance for targeted station
      if (target.kind == TargetKind::Station && target.index == i) {
        const math::Vec3d stPos = stationPosKm(st, timeDays);
        const math::Quatd stQ = stationOrient(st, stPos, timeDays);
        const math::Vec3d axis = stQ.rotate({0,0,1});
        const math::Vec3d right = stQ.rotate({1,0,0});
        const math::Vec3d up = stQ.rotate({0,1,0});

        const double zEntrance = st.radiusKm * 1.10;
        const double zEnd = zEntrance + st.approachLengthKm;

        const math::Vec3d startKm = stPos + axis * zEntrance;
        const math::Vec3d endKm = stPos + axis * zEnd;

        bool clearanceGranted = false;
        if (auto it = clearances.find(st.id); it != clearances.end()) {
          clearanceGranted = sim::dockingClearanceValid(it->second, timeDays);
        }

        const float cr = clearanceGranted ? 0.35f : 0.85f;
        const float cg = clearanceGranted ? 0.85f : 0.55f;
        const float cb = clearanceGranted ? 0.45f : 0.25f;

        auto addLineKm = [&](const math::Vec3d& aKm, const math::Vec3d& bKm, float r, float g, float b) {
          const math::Vec3d aU = toRenderU(aKm);
          const math::Vec3d bU = toRenderU(bKm);
          lines.push_back({(float)aU.x,(float)aU.y,(float)aU.z, r,g,b});
          lines.push_back({(float)bU.x,(float)bU.y,(float)bU.z, r,g,b});
        };

        // Centerline
        addLineKm(startKm, endKm, cr,cg,cb);

        // Corridor wireframe box
        const double hx = st.approachRadiusKm;
        const double hy = st.approachRadiusKm;

        math::Vec3d c0[4] = {
          startKm + right*hx + up*hy,
          startKm + right*hx - up*hy,
          startKm - right*hx - up*hy,
          startKm - right*hx + up*hy,
        };
        math::Vec3d c1[4] = {
          endKm + right*hx + up*hy,
          endKm + right*hx - up*hy,
          endKm - right*hx - up*hy,
          endKm - right*hx + up*hy,
        };

        for (int k = 0; k < 4; ++k) {
          const int k2 = (k + 1) % 4;
          addLineKm(c0[k], c0[k2], cr,cg,cb);
          addLineKm(c1[k], c1[k2], cr,cg,cb);
          addLineKm(c0[k], c1[k], cr,cg,cb);
        }

        // Slot frame at the corridor entrance (yellow)
        const double sw = st.slotWidthKm * 0.5;
        const double sh = st.slotHeightKm * 0.5;
        math::Vec3d s[4] = {
          startKm + right*sw + up*sh,
          startKm + right*sw - up*sh,
          startKm - right*sw - up*sh,
          startKm - right*sw + up*sh,
        };
        for (int k = 0; k < 4; ++k) {
          addLineKm(s[k], s[(k + 1) % 4], 0.95f,0.9f,0.15f);
        }

        // Small staging cross at the corridor end (cyan)
        const double cross = std::max(250.0, st.approachRadiusKm * 0.15);
        addLineKm(endKm - right*cross, endKm + right*cross, 0.25f,0.7f,1.0f);
        addLineKm(endKm - up*cross, endKm + up*cross, 0.25f,0.7f,1.0f);
      }
    }

    // Laser beams
    for (const auto& b : beams) {
      lines.push_back({(float)b.aU.x,(float)b.aU.y,(float)b.aU.z, b.r,b.g,b.b});
      lines.push_back({(float)b.bU.x,(float)b.bU.y,(float)b.bU.z, b.r,b.g,b.b});
    }

    // Projectile tracers (draw a fixed-length tail so they're visible at astronomical scales)
    for (const auto& p : projectiles) {
      math::Vec3d tailKm = p.prevKm;
      if (p.velKmS.lengthSq() > 1e-12) {
        const double tracerLenKm = 15000.0;
        tailKm = p.posKm - p.velKmS.normalized() * tracerLenKm;
      }
      const math::Vec3d aU = toRenderU(tailKm);
      const math::Vec3d bU = toRenderU(p.posKm);
      lines.push_back({(float)aU.x,(float)aU.y,(float)aU.z, p.r,p.g,p.b});
      lines.push_back({(float)bU.x,(float)bU.y,(float)bU.z, p.r,p.g,p.b});
    }

    // Missile tracers (longer tail so you can visually track guidance)
    for (const auto& m : missiles) {
      math::Vec3d tailKm = m.prevKm;
      if (m.velKmS.lengthSq() > 1e-12) {
        const double tracerLenKm = 24000.0;
        tailKm = m.posKm - m.velKmS.normalized() * tracerLenKm;
      }
      const math::Vec3d aU = toRenderU(tailKm);
      const math::Vec3d bU = toRenderU(m.posKm);
      lines.push_back({(float)aU.x,(float)aU.y,(float)aU.z, m.r,m.g,m.b});
      lines.push_back({(float)bU.x,(float)bU.y,(float)bU.z, m.r,m.g,m.b});
    }

    // Station geometry (cubes)
    std::vector<render::InstanceData> cubes;
    std::vector<render::InstanceData> shipCubes;
    cubes.reserve(1 + currentSystem->stations.size() * 18 + contacts.size());
    shipCubes.reserve(1);

    for (const auto& st : currentSystem->stations) {
      const math::Vec3d stPos = stationPosKm(st, timeDays);
      const math::Quatd stQ = stationOrient(st, stPos, timeDays);
      emitStationGeometry(st, stPos, stQ, cubes);
    }

    // Ship instance (cube, rotated)
    {
      const auto inst = makeInst(toRenderU(ship.positionKm()),
                                 {0.35, 0.20, 0.60},
                                 ship.orientation(),
                                 1.00f, 1.00f, 1.00f);

      // If the livery system is enabled, draw the ship in a separate pass with
      // the procedurally generated texture.
      if (liveryCfg.applyInWorld && shipLiveryTex.handle() != 0) {
        shipCubes.push_back(inst);
      } else {
        cubes.push_back(inst);
      }
    }

    // Contacts (pirates)
    for (const auto& c : contacts) {
      if (!c.alive) continue;
      float r = 0.65f, g = 0.75f, b = 0.85f;
      if (c.role == ContactRole::Pirate) { r = 1.0f; g = 0.25f; b = 0.25f; }
      if (c.role == ContactRole::Police) { r = 0.35f; g = 0.75f; b = 1.0f; }
      if (c.role == ContactRole::Trader) { r = 0.45f; g = 0.95f; b = 0.45f; }
      cubes.push_back(makeInst(toRenderU(c.ship.positionKm()),
                               {0.25, 0.18, 0.45},
                               c.ship.orientation(),
                               r,g,b));
    }

	    // Floating salvage / cargo pods
	    for (const auto& pod : floatingCargo) {
	      if (pod.units <= 0.0) continue;
	      const float r = 1.00f, g = 0.85f, b = 0.25f;
	      cubes.push_back(makeInst(toRenderU(pod.posKm),
	                               {0.12, 0.12, 0.12},
	                               math::Quatd::identity(),
	                               r,g,b));
	    }

	    // Asteroid mining nodes
	    for (const auto& a : asteroids) {
	      const float r = 0.55f, g = 0.55f, b = 0.58f;
	      const double s = std::clamp(a.radiusKm / 3500.0, 0.35, 1.25);
	      cubes.push_back(makeInst(toRenderU(a.posKm),
	                               {0.55 * s, 0.55 * s, 0.55 * s},
	                               math::Quatd::identity(),
	                               r,g,b));
	    }

	    // Signal sources (distress / derelicts / resource sites)
	    for (const auto& s : signals) {
	      float r = 0.85f, g = 0.85f, b = 0.85f;
	      if (s.type == SignalType::Distress) { r = 1.0f; g = 0.6f; b = 0.2f; }
	      if (s.type == SignalType::Derelict) { r = 0.75f; g = 0.75f; b = 1.0f; }
	      if (s.type == SignalType::Resource) { r = 0.55f; g = 0.95f; b = 0.55f; }
	      cubes.push_back(makeInst(toRenderU(s.posKm),
	                               {0.25, 0.25, 0.25},
	                               math::Quatd::identity(),
	                               r,g,b));
	    }

    // ---- Render ---
    // World pass
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);

    if (postFxSettings.enabled) {
      postFx.ensureSize(w, h);
      postFx.beginScene(w, h);
    } else {
      render::gl::BindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    glViewport(0, 0, w, h);
    glClearColor(0.01f, 0.01f, 0.02f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Background nebula (large point sprites, additive blend). Don't write depth.
    if (vfxNebulaEnabled) {
      glDepthMask(GL_FALSE);
      pointRenderer.drawPointsSprite(nebula.points(), nebulaPointSpriteTex, render::PointBlendMode::Additive);
      glDepthMask(GL_TRUE);
    }

    // Background stars (point sprites, additive blend). Don't write depth.
    if (vfxStarfieldEnabled) {
      glDepthMask(GL_FALSE);
      if (vfxStarfieldTextured) {
        pointRenderer.drawPointsSprite(starfield.points(), starPointSpriteTex, render::PointBlendMode::Additive);
      } else {
        pointRenderer.drawPoints(starfield.points(), render::PointBlendMode::Additive);
      }
      glDepthMask(GL_TRUE);
    }

    // Lines (orbits, corridor, beams)
    lineRenderer.drawLines(lines);

    // Star + planets
    meshRenderer.setMesh(&sphere);
    meshRenderer.setAlphaFromTexture(false);
    meshRenderer.setAlphaMul(1.0f);
    std::vector<render::InstanceData> tmp;
    tmp.reserve(1);
    if (worldUseProceduralSurfaces) {
      // Star (procedural surface + optional unlit)
      {
        const auto& starSurf = surfaceTexCache.get(render::SurfaceKind::Star, starSurfaceSeed, worldSurfaceTexWidth);
        meshRenderer.setTexture(&starSurf);
        meshRenderer.setUnlit(worldStarUnlit);
        tmp.clear();
        tmp.push_back(starSphereInst);
        meshRenderer.drawInstances(tmp);
        meshRenderer.setUnlit(false);
      }

      // Planets (per-type procedural surfaces)
      for (const auto& pd : planetSurfaceDraws) {
        const auto& surf = surfaceTexCache.get(pd.surface, pd.surfaceSeed, worldSurfaceTexWidth);
        meshRenderer.setTexture(&surf);
        tmp.clear();
        tmp.push_back(pd.inst);
        meshRenderer.drawInstances(tmp);
      }

      // Cloud shell (alpha-blended) on top of the planet surface.
      if (worldCloudsEnabled && !planetCloudDraws.empty()) {
        glDepthMask(GL_FALSE);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        meshRenderer.setAlphaFromTexture(true);
        meshRenderer.setUnlit(false);

        for (const auto& cd : planetCloudDraws) {
          const auto& cloud = surfaceTexCache.get(render::SurfaceKind::Clouds, cd.cloudSeed, worldSurfaceTexWidth);
          meshRenderer.setTexture(&cloud);
          meshRenderer.setAlphaMul(cd.alphaMul);
          tmp.clear();
          tmp.push_back(cd.inst);
          meshRenderer.drawInstances(tmp);
        }

        meshRenderer.setAlphaMul(1.0f);
        meshRenderer.setAlphaFromTexture(false);
        glDepthMask(GL_TRUE);
      }

      // Restore default texture for other meshes.
      meshRenderer.setTexture(&checker);
    } else {
      // Fallback (single instanced draw, checker texture)
      meshRenderer.setTexture(&checker);
      meshRenderer.setUnlit(false);
      meshRenderer.drawInstances(spheres);
    }

    // Atmosphere limb glow (additive). This pass is independent of the surface texture mode.
    if (worldAtmospheresEnabled && !planetAtmoDraws.empty()) {
      glDepthMask(GL_FALSE);
      atmosphereRenderer.setIntensity(worldAtmoIntensity);
      atmosphereRenderer.setPower(worldAtmoPower);
      atmosphereRenderer.setSunLitBoost(worldAtmoSunLitBoost);
      atmosphereRenderer.setForwardScatter(worldAtmoForwardScatter);
      atmosphereRenderer.drawInstances(planetAtmoDraws);
      glDepthMask(GL_TRUE);
    }

    // Cubes (stations, contacts, salvage, asteroids, signal sources)
    meshRenderer.setMesh(&cube);
    meshRenderer.setTexture(&checker);
    meshRenderer.setUnlit(false);
    meshRenderer.setAlphaFromTexture(false);
    meshRenderer.drawInstances(cubes);

    // Player ship livery pass (separate texture).
    if (!shipCubes.empty()) {
      meshRenderer.setTexture(&shipLiveryTex);
      meshRenderer.setUnlit(false);
      meshRenderer.setAlphaFromTexture(false);
      meshRenderer.drawInstances(shipCubes);
      meshRenderer.setTexture(&checker);
    }

    // Particles (thrusters, impacts, explosions). Depth-tested but no depth writes.
    if (vfxParticlesEnabled) {
      glDepthMask(GL_FALSE);
      if (vfxParticlesTextured) {
        pointRenderer.drawPointsSprite(particleVerts, particlePointSpriteTex, render::PointBlendMode::Additive);
      } else {
        pointRenderer.drawPoints(particleVerts, render::PointBlendMode::Additive);
      }
      glDepthMask(GL_TRUE);
    }

    // PostFX pass (HDR tonemap / bloom) - draw scene to the backbuffer, then render ImGui on top.
    if (postFxSettings.enabled) {
      // Auto-warp from ship speed gives a subtle "speed lines" vibe in supercruise without touching gameplay.
      render::PostFXSettings s = postFxSettings;
      if (postFxAutoWarpFromSpeed) {
        const double spdKmS = ship.velocityKmS().length();
        const float a = (float)std::clamp((spdKmS - 25.0) / 450.0, 0.0, 1.0);
        s.warp = std::max(s.warp, a * 0.030f);
      }

      // FSD-driven hyperspace tunnel (purely visual).
      if (postFxAutoHyperspaceFromFsd && fsdState != FsdState::Idle) {
        auto smooth01 = [](float x) -> float {
          x = std::clamp(x, 0.0f, 1.0f);
          return x * x * (3.0f - 2.0f * x);
        };

        float fx = 0.0f;
        if (fsdState == FsdState::Charging) {
          const float t = 1.0f - (float)std::clamp(fsdChargeRemainingSec / kFsdChargeSec, 0.0, 1.0);
          fx = smooth01(t);
        } else if (fsdState == FsdState::Jumping) {
          const float total = (float)std::max(0.001, fsdTravelTotalSec);
          const float t = 1.0f - (float)std::clamp(fsdTravelRemainingSec / (double)total, 0.0, 1.0);
          const float fadeIn = smooth01(std::clamp(t / 0.12f, 0.0f, 1.0f));
          const float fadeOut = 1.0f - smooth01(std::clamp((t - 0.80f) / 0.20f, 0.0f, 1.0f));
          fx = std::clamp(fadeIn * fadeOut, 0.0f, 1.0f);
        }

        s.hyperspace = std::max(s.hyperspace, fx * postFxFsdHyperspaceBoost);
        s.warp = std::max(s.warp, fx * postFxFsdWarpBoost);
      }

      postFx.present(w, h, s, (float)timeRealSec);
    } else {
      // Ensure UI draws to the backbuffer even if PostFX was previously enabled.
      render::gl::BindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    // ---- UI ----
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplSDL2_NewFrame();
    ImGui::NewFrame();

    // HUD layout editor: allow dragging borderless overlay windows.
    // (ImGui defaults to moving windows from the title bar only.)
    io.ConfigWindowsMoveFromTitleBarOnly = !hudLayoutEditMode;

    // Keep persisted toggles synced with runtime settings so saving is always accurate.
    hudLayout.widget(ui::HudWidgetId::Radar).enabled = showRadarHud;
    hudLayout.widget(ui::HudWidgetId::Objective).enabled = objectiveHudEnabled;
    hudLayout.widget(ui::HudWidgetId::Threat).enabled = hudThreatOverlayEnabled;
    hudLayout.widget(ui::HudWidgetId::Jump).enabled = hudJumpOverlay;

    // Keep HUD settings synced so manual save / auto-save always matches runtime.
    syncHudSettingsFromRuntime();
    hudSettingsDirty = !hudSettingsEquivalent(hudSettings, hudSettingsSaved);

    // Menu bar + (optional) DockSpace layout.
    // When built against Dear ImGui's docking branch/tag, we create a fullscreen dockspace and
    // seed a sensible default layout. Without docking support, windows simply float as before.

    auto drawMenuBarContents = [&]() {
      auto withKey = [&](const char* label, const game::KeyChord& chord) -> std::string {
        const std::string k = game::chordLabel(chord);
        if (k == "(unbound)") return std::string(label);
        return std::string(label) + " (" + k + ")";
      };

      if (ImGui::BeginMenu("Windows")) {
        ImGui::MenuItem(withKey("Galaxy", controls.actions.toggleGalaxy).c_str(), nullptr, &showGalaxy);
        ImGui::MenuItem("Bookmarks", nullptr, &showBookmarksWindow);
        ImGui::MenuItem(withKey("Ship/Status", controls.actions.toggleShip).c_str(), nullptr, &showShip);
        ImGui::MenuItem(withKey("Market", controls.actions.toggleMarket).c_str(), nullptr, &showEconomy);
        ImGui::MenuItem(withKey("Contacts", controls.actions.toggleContacts).c_str(), nullptr, &showContacts);
        ImGui::MenuItem(withKey("Missions", controls.actions.toggleMissions).c_str(), nullptr, &showMissions);
        ImGui::MenuItem(withKey("Scanner", controls.actions.toggleScanner).c_str(), nullptr, &showScanner);
        ImGui::MenuItem(withKey("Trade Finder", controls.actions.toggleTrade).c_str(), nullptr, &showTrade);
        ImGui::MenuItem(withKey("Guide", controls.actions.toggleGuide).c_str(), nullptr, &showGuide);
        ImGui::MenuItem(withKey("Hangar", controls.actions.toggleHangar).c_str(), nullptr, &showHangar);
        ImGui::MenuItem(withKey("World Visuals", controls.actions.toggleWorldVisuals).c_str(), nullptr, &showWorldVisuals);
        ImGui::Separator();
        ImGui::MenuItem(withKey("Sprite Lab", controls.actions.toggleSpriteLab).c_str(), nullptr, &showSprites);
        ImGui::MenuItem(withKey("VFX Lab", controls.actions.toggleVfxLab).c_str(), nullptr, &showVfx);
        ImGui::MenuItem(withKey("Post FX", controls.actions.togglePostFx).c_str(), nullptr, &showPostFx);
        ImGui::Separator();
        ImGui::MenuItem(withKey("Controls", controls.actions.toggleControlsWindow).c_str(), nullptr, &controlsWindow.open);
        ImGui::MenuItem("Notifications", nullptr, &showNotifications);
        if (ImGui::MenuItem("Console", nullptr, &consoleWindow.open)) { if (consoleWindow.open) consoleWindow.focusInput = true; }
        ImGui::EndMenu();
      }

      if (ImGui::BeginMenu("HUD")) {
        ImGui::MenuItem(withKey("Radar", controls.actions.toggleRadarHud).c_str(), nullptr, &showRadarHud);
        ImGui::MenuItem("Objective", nullptr, &objectiveHudEnabled);
        ImGui::MenuItem("Threat", nullptr, &hudThreatOverlayEnabled);
        ImGui::MenuItem("Jump", nullptr, &hudJumpOverlay);
        ImGui::Separator();
        ImGui::MenuItem(withKey("Tactical overlay", controls.actions.toggleTacticalOverlay).c_str(), nullptr, &showTacticalOverlay);
        ImGui::Separator();
        ImGui::MenuItem(withKey("Edit HUD layout", controls.actions.hudLayoutToggleEdit).c_str(), nullptr, &hudLayoutEditMode);
        ImGui::MenuItem("HUD layout window", nullptr, &showHudLayoutWindow);
        ImGui::MenuItem("HUD settings window", nullptr, &showHudSettingsWindow);
        ImGui::Separator();

        const double minRadarKm = 25000.0;
        const double maxRadarKm = 1200000.0;
        ImGui::SliderScalar("Radar range (km)", ImGuiDataType_Double, &radarRangeKm,
                            &minRadarKm, &maxRadarKm, "%.0f", ImGuiSliderFlags_Logarithmic);

        float safePx = (float)hudLayout.safeMarginPx;
        if (ImGui::SliderFloat("Safe margin (px)", &safePx, 0.0f, 80.0f, "%.0f")) {
          hudLayout.safeMarginPx = safePx;
        }

        if (ImGui::MenuItem(withKey("Save HUD layout", controls.actions.hudLayoutSave).c_str())) {
          hudLayout.widget(ui::HudWidgetId::Radar).enabled = showRadarHud;
          hudLayout.widget(ui::HudWidgetId::Objective).enabled = objectiveHudEnabled;
          hudLayout.widget(ui::HudWidgetId::Threat).enabled = hudThreatOverlayEnabled;
          hudLayout.widget(ui::HudWidgetId::Jump).enabled = hudJumpOverlay;
          const bool ok = ui::saveToFile(hudLayout, hudLayoutPath);
          toast(toasts, ok ? "Saved HUD layout." : "Failed to save HUD layout.", 1.8);
        }
        if (ImGui::MenuItem(withKey("Load HUD layout", controls.actions.hudLayoutLoad).c_str())) {
          ui::HudLayout loaded = ui::makeDefaultHudLayout();
          if (ui::loadFromFile(hudLayoutPath, loaded)) {
            hudLayout = loaded;
            showRadarHud = hudLayout.widget(ui::HudWidgetId::Radar).enabled;
            objectiveHudEnabled = hudLayout.widget(ui::HudWidgetId::Objective).enabled;
            hudThreatOverlayEnabled = hudLayout.widget(ui::HudWidgetId::Threat).enabled;
            hudJumpOverlay = hudLayout.widget(ui::HudWidgetId::Jump).enabled;
            toast(toasts, "Loaded HUD layout.", 1.8);
          } else {
            toast(toasts, "Failed to load HUD layout.", 1.8);
          }
        }
        if (ImGui::MenuItem(withKey("Reset HUD layout", controls.actions.hudLayoutReset).c_str())) {
          hudLayout = ui::makeDefaultHudLayout();
          showRadarHud = hudLayout.widget(ui::HudWidgetId::Radar).enabled;
          objectiveHudEnabled = hudLayout.widget(ui::HudWidgetId::Objective).enabled;
          hudThreatOverlayEnabled = hudLayout.widget(ui::HudWidgetId::Threat).enabled;
          hudJumpOverlay = hudLayout.widget(ui::HudWidgetId::Jump).enabled;
          toast(toasts, "HUD layout reset.", 1.6);
        }
        ImGui::EndMenu();
      }

      if (ImGui::BeginMenu("Controls")) {
        ImGui::MenuItem(withKey("Controls window", controls.actions.toggleControlsWindow).c_str(), nullptr, &controlsWindow.open);
        ImGui::Separator();
        if (ImGui::MenuItem("Save controls")) {
          const bool ok = game::saveToFile(controls, controlsPath);
          if (ok) controlsDirty = false;
          toast(toasts, ok ? "Saved controls." : "Failed to save controls.", 1.8);
        }
        if (ImGui::MenuItem("Load controls")) {
          game::ControlsConfig loaded = game::makeDefaultControls();
          if (game::loadFromFile(controlsPath, loaded)) {
            controls = loaded;
            controlsDirty = false;
            toast(toasts, "Loaded controls.", 1.8);
          } else {
            toast(toasts, "No controls file found. Using defaults.", 2.2);
          }
        }
        if (ImGui::MenuItem("Restore defaults")) {
          controls = game::makeDefaultControls();
          controlsDirty = true;
          toast(toasts, "Controls reset to defaults.", 1.8);
        }
        ImGui::Separator();
        ImGui::MenuItem("Auto-save on exit", nullptr, &controlsAutoSaveOnExit);
        if (controlsDirty) {
          ImGui::TextDisabled("(unsaved changes)");
        }
        ImGui::EndMenu();
      }

      if (ImGui::BeginMenu("UI")) {
        if (ImGui::MenuItem(withKey("Command Palette...", controls.actions.commandPalette).c_str())) {
          game::openCommandPalette(commandPalette);
        }
        ImGui::MenuItem("Notifications", nullptr, &showNotifications);
        if (ImGui::MenuItem("Console", nullptr, &consoleWindow.open)) { if (consoleWindow.open) consoleWindow.focusInput = true; }
        ImGui::MenuItem("Bookmarks", nullptr, &showBookmarksWindow);
        ImGui::MenuItem("UI Settings...", nullptr, &showUiSettingsWindow);
        ImGui::Separator();

        ImGui::MenuItem("ImGui Demo Window", nullptr, &showImGuiDemo);
        ImGui::MenuItem("ImGui Metrics Window", nullptr, &showImGuiMetrics);

        ImGui::Separator();
        if (ImGui::BeginMenu("Theme")) {
          const auto setTheme = [&](ui::UiTheme t) {
            if (uiTheme == t) return;
            uiTheme = t;
            rebuildUiStyle(uiTheme, uiScale);
            io.FontGlobalScale = uiScale;
            uiScaleApplied = uiScale;
            uiSettingsDirty = true;
          };

          if (ImGui::MenuItem("Dark", nullptr, uiTheme == ui::UiTheme::Dark)) setTheme(ui::UiTheme::Dark);
          if (ImGui::MenuItem("Light", nullptr, uiTheme == ui::UiTheme::Light)) setTheme(ui::UiTheme::Light);
          if (ImGui::MenuItem("Classic", nullptr, uiTheme == ui::UiTheme::Classic)) setTheme(ui::UiTheme::Classic);
          if (ImGui::MenuItem("High Contrast", nullptr, uiTheme == ui::UiTheme::HighContrast)) setTheme(ui::UiTheme::HighContrast);
          ImGui::EndMenu();
        }

        ImGui::Separator();
        if (ImGui::Checkbox("Auto scale from DPI", &uiAutoScaleFromDpi)) {
          recomputeUiDpiScale();
          applyUiScaleNow();
          uiSettingsDirty = true;
        }

        ImGui::SetNextItemWidth(180.0f);
        if (ImGui::SliderFloat("User scale", &uiScaleUser, 0.50f, 2.50f, "%.2fx")) {
          applyUiScaleNow();
          uiSettingsDirty = true;
        }

        ImGui::TextDisabled("Effective: %.2fx", uiScale);

        ImGui::Separator();
        if (ImGui::MenuItem("Save UI settings")) {
          syncUiSettingsFromRuntime();
          const bool ok = ui::saveUiSettingsToFile(uiSettings, uiSettingsPath);
          uiSettingsDirty = uiSettingsDirty && !ok;
          if (ok) uiSettingsDirty = false;
          toast(toasts, ok ? "Saved UI settings." : "Failed to save UI settings.", 1.8);
        }
        if (ImGui::MenuItem("Load UI settings")) {
          ui::UiSettings loaded = ui::makeDefaultUiSettings();
          if (ui::loadUiSettingsFromFile(uiSettingsPath, loaded)) {
            uiSettings = loaded;
            applyUiSettingsToRuntime(uiSettings, /*loadImGuiIni=*/true);
            uiSettingsDirty = false;
            toast(toasts, "Loaded UI settings.", 1.8);
          } else {
            toast(toasts, "No UI settings file found.", 1.8);
          }
        }
        if (ImGui::MenuItem("Auto-save UI settings on exit", nullptr, &uiSettingsAutoSaveOnExit)) {
          uiSettingsDirty = true;
        }
        if (uiSettingsDirty) {
          ImGui::TextDisabled("(unsaved changes)");
        }

#ifdef IMGUI_HAS_DOCK
        ImGui::Separator();
        ImGui::TextDisabled("Docking");
        if (ImGui::Checkbox("Enable DockSpace", &uiDockingEnabled)) uiSettingsDirty = true;
        if (ImGui::MenuItem("Reset dock layout")) uiDockResetLayout = true;
        if (ImGui::Checkbox("Passthrough central view", &uiDockPassthruCentral)) { uiDockResetLayout = true; uiSettingsDirty = true; }
        if (ImGui::Checkbox("Lock central view", &uiDockLockCentralView)) { uiDockResetLayout = true; uiSettingsDirty = true; }

        ImGui::SetNextItemWidth(180.0f);
        if (ImGui::SliderFloat("Left width", &uiDockLeftRatio, 0.10f, 0.45f, "%.2f")) { uiDockResetLayout = true; uiSettingsDirty = true; }
        ImGui::SetNextItemWidth(180.0f);
        if (ImGui::SliderFloat("Right width", &uiDockRightRatio, 0.10f, 0.45f, "%.2f")) { uiDockResetLayout = true; uiSettingsDirty = true; }
        ImGui::SetNextItemWidth(180.0f);
        if (ImGui::SliderFloat("Bottom height", &uiDockBottomRatio, 0.10f, 0.45f, "%.2f")) { uiDockResetLayout = true; uiSettingsDirty = true; }

        ImGui::TextDisabled("Tip: drag window tab bars to rearrange panels.");
#endif
        ImGui::EndMenu();
      }

      // Right-side status text.
      std::string status;
      status.reserve(256);
      status += docked ? "Docked" : "In flight";
      status += " | System: ";
      status += currentStub.name;

      status += " | Credits: ";
      status += std::to_string((int)std::round(credits));
      status += " cr";

      if (currentSystem) {
        status += " | Jurisdiction: ";
        status += factionName(currentSystem->stub.factionId);
      }

      const int fpsNow = (dtReal > 1e-6) ? (int)std::round(1.0 / dtReal) : 0;
      status += " | ";
      status += std::to_string(fpsNow);
      status += " fps";

      if (controlsDirty) status += " | controls*";
	      if (uiSettingsDirty) status += " | ui*";
	      if (bookmarksDirty) status += " | bm*";

      const float wText = ImGui::CalcTextSize(status.c_str()).x;
      ImGui::SameLine(ImGui::GetWindowWidth() - wText - 14.0f);
      ImGui::TextUnformatted(status.c_str());
    };

#ifdef IMGUI_HAS_DOCK
    if (uiDockingEnabled) {
      ImGuiDockNodeFlags dockFlags = ImGuiDockNodeFlags_None;
      if (uiDockPassthruCentral) dockFlags |= ImGuiDockNodeFlags_PassthruCentralNode;
      if (uiDockLockCentralView) dockFlags |= ImGuiDockNodeFlags_NoDockingOverCentralNode;

      ImGuiViewport* viewport = ImGui::GetMainViewport();
      ImGui::SetNextWindowPos(viewport->Pos);
      ImGui::SetNextWindowSize(viewport->Size);
      ImGui::SetNextWindowViewport(viewport->ID);

      ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
      ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
      ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
      if (dockFlags & ImGuiDockNodeFlags_PassthruCentralNode) ImGui::SetNextWindowBgAlpha(0.0f);

      ImGuiWindowFlags hostFlags = ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar
                                   | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize
                                   | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoBringToFrontOnFocus
                                   | ImGuiWindowFlags_NoNavFocus | ImGuiWindowFlags_MenuBar;
      if (dockFlags & ImGuiDockNodeFlags_PassthruCentralNode) hostFlags |= ImGuiWindowFlags_NoBackground;

      ImGui::Begin("##DockSpaceHost", nullptr, hostFlags);
      ImGui::PopStyleVar(3);

      if (ImGui::BeginMenuBar()) {
        drawMenuBarContents();
        ImGui::EndMenuBar();
      }

      ImGuiID dockspaceId = ImGui::GetID("StellarForgeDockSpace");
      ImGui::DockSpace(dockspaceId, ImVec2(0.0f, 0.0f), dockFlags);

      if (uiDockResetLayout) {
        // Build a predictable default layout around a free central 3D view.
        ImGui::DockBuilderRemoveNode(dockspaceId);
        ImGui::DockBuilderAddNode(dockspaceId, dockFlags | ImGuiDockNodeFlags_DockSpace);
        ImGui::DockBuilderSetNodeSize(dockspaceId, viewport->Size);

        ImGuiID dockMain = dockspaceId;
        ImGuiID dockLeft = 0, dockRight = 0, dockBottom = 0;
        ImGui::DockBuilderSplitNode(dockMain, ImGuiDir_Left, uiDockLeftRatio, &dockLeft, &dockMain);
        ImGui::DockBuilderSplitNode(dockMain, ImGuiDir_Right, uiDockRightRatio, &dockRight, &dockMain);
        ImGui::DockBuilderSplitNode(dockMain, ImGuiDir_Down, uiDockBottomRatio, &dockBottom, &dockMain);

        ImGuiID dockLeftBottom = 0;
        ImGui::DockBuilderSplitNode(dockLeft, ImGuiDir_Down, 0.45f, &dockLeftBottom, &dockLeft);

        ImGuiID dockRightBottom = 0;
        ImGui::DockBuilderSplitNode(dockRight, ImGuiDir_Down, 0.52f, &dockRightBottom, &dockRight);

        // Left column
        ImGui::DockBuilderDockWindow("Galaxy / Streaming", dockLeft);
        ImGui::DockBuilderDockWindow("Contacts / Combat", dockLeftBottom);
        ImGui::DockBuilderDockWindow("Missions", dockLeftBottom);

        // Right column
        ImGui::DockBuilderDockWindow("Ship / Status", dockRight);
        ImGui::DockBuilderDockWindow("Trade Helper", dockRightBottom);
        ImGui::DockBuilderDockWindow("Market Details", dockRightBottom);
        ImGui::DockBuilderDockWindow("Dock / Station", dockRightBottom);

        // Bottom strip (tools)
        ImGui::DockBuilderDockWindow("System Scanner", dockBottom);
        ImGui::DockBuilderDockWindow("Pilot Guide", dockBottom);
        ImGui::DockBuilderDockWindow("HUD Layout", dockBottom);
        ImGui::DockBuilderDockWindow("World Visuals", dockBottom);
        ImGui::DockBuilderDockWindow("Sprite Lab", dockBottom);
        ImGui::DockBuilderDockWindow("VFX Lab", dockBottom);
        ImGui::DockBuilderDockWindow("Post FX", dockBottom);

        ImGui::DockBuilderFinish(dockspaceId);

        // After a reset, show the primary panels so the layout is immediately useful.
        showShip = true;
        showContacts = true;
        showMissions = true;

        uiDockResetLayout = false;
      }

      ImGui::End();
    } else
#endif
    {
      // Fallback (no docking / disabled dockspace): classic main menu bar + floating windows.
      if (ImGui::BeginMainMenuBar()) {
        drawMenuBarContents();
        ImGui::EndMainMenuBar();
      }
    }

    if (showImGuiDemo) ImGui::ShowDemoWindow(&showImGuiDemo);
    if (showImGuiMetrics) ImGui::ShowMetricsWindow(&showImGuiMetrics);

    // HUD Layout editor window
    if (showHudLayoutWindow) {
      ImGui::SetNextWindowSize(ImVec2(520.0f, 450.0f), ImGuiCond_FirstUseEver);
      ImGui::Begin("HUD Layout", &showHudLayoutWindow);

      ImGui::TextDisabled(
          "Reposition in-flight HUD overlays. Drag the HUD panels while Edit mode is enabled.");
      {
        std::string shortcuts;
        shortcuts.reserve(160);
        shortcuts += "Shortcuts: ";
        shortcuts += game::chordLabel(controls.actions.hudLayoutToggleEdit);
        shortcuts += " edit | ";
        shortcuts += game::chordLabel(controls.actions.hudLayoutSave);
        shortcuts += " save | ";
        shortcuts += game::chordLabel(controls.actions.hudLayoutLoad);
        shortcuts += " load | ";
        shortcuts += game::chordLabel(controls.actions.hudLayoutReset);
        shortcuts += " reset";
        ImGui::TextDisabled("%s", shortcuts.c_str());
      }

      ImGui::Checkbox("Edit mode (drag panels)", &hudLayoutEditMode);
      ImGui::SameLine();
      ImGui::Checkbox("Auto-save on exit", &hudLayoutAutoSaveOnExit);

      {
        float safePx = (float)hudLayout.safeMarginPx;
        ImGui::SetNextItemWidth(240.0f);
        if (ImGui::SliderFloat("Safe margin guide (px)", &safePx, 0.0f, 64.0f, "%.0f")) {
          hudLayout.safeMarginPx = (double)safePx;
        }
      }

      const std::string saveLabel = "Save (" + game::chordLabel(controls.actions.hudLayoutSave) + ")";
      if (ImGui::Button(saveLabel.c_str())) {
        // Sync enabled toggles into the layout before saving.
        hudLayout.widget(ui::HudWidgetId::Radar).enabled = showRadarHud;
        hudLayout.widget(ui::HudWidgetId::Objective).enabled = objectiveHudEnabled;
        hudLayout.widget(ui::HudWidgetId::Threat).enabled = hudThreatOverlayEnabled;
        hudLayout.widget(ui::HudWidgetId::Jump).enabled = hudJumpOverlay;
        if (ui::saveToFile(hudLayout, hudLayoutPath)) {
          toast(toasts, "Saved HUD layout to " + hudLayoutPath, 2.0);
        } else {
          toast(toasts, "Failed to save HUD layout.", 2.0);
        }
      }
      ImGui::SameLine();
      const std::string loadLabel = "Load (" + game::chordLabel(controls.actions.hudLayoutLoad) + ")";
      if (ImGui::Button(loadLabel.c_str())) {
        ui::HudLayout loaded = ui::makeDefaultHudLayout();
        if (ui::loadFromFile(hudLayoutPath, loaded)) {
          hudLayout = loaded;
          showRadarHud = hudLayout.widget(ui::HudWidgetId::Radar).enabled;
          objectiveHudEnabled = hudLayout.widget(ui::HudWidgetId::Objective).enabled;
          hudThreatOverlayEnabled = hudLayout.widget(ui::HudWidgetId::Threat).enabled;
          hudJumpOverlay = hudLayout.widget(ui::HudWidgetId::Jump).enabled;
          toast(toasts, "Loaded HUD layout from " + hudLayoutPath, 2.0);
        } else {
          toast(toasts, "HUD layout file not found.", 2.0);
        }
      }
      ImGui::SameLine();
      const std::string resetLabel = "Reset (" + game::chordLabel(controls.actions.hudLayoutReset) + ")";
      if (ImGui::Button(resetLabel.c_str())) {
        hudLayout = ui::makeDefaultHudLayout();
        showRadarHud = hudLayout.widget(ui::HudWidgetId::Radar).enabled;
        objectiveHudEnabled = hudLayout.widget(ui::HudWidgetId::Objective).enabled;
        hudThreatOverlayEnabled = hudLayout.widget(ui::HudWidgetId::Threat).enabled;
        hudJumpOverlay = hudLayout.widget(ui::HudWidgetId::Jump).enabled;
        toast(toasts, "HUD layout reset to defaults.", 2.0);
      }

      ImGui::Separator();

      // Anchor presets (pivot)
      struct AnchorPreset { const char* label; float px; float py; };
      static constexpr AnchorPreset kAnchors[] = {
        {"TopLeft", 0.0f, 0.0f},
        {"Top", 0.5f, 0.0f},
        {"TopRight", 1.0f, 0.0f},
        {"Left", 0.0f, 0.5f},
        {"Center", 0.5f, 0.5f},
        {"Right", 1.0f, 0.5f},
        {"BottomLeft", 0.0f, 1.0f},
        {"Bottom", 0.5f, 1.0f},
        {"BottomRight", 1.0f, 1.0f},
      };

      auto anchorIndexFor = [&](float px, float py) -> int {
        for (int i = 0; i < (int)std::size(kAnchors); ++i) {
          if (std::abs(px - kAnchors[i].px) < 1e-4f && std::abs(py - kAnchors[i].py) < 1e-4f) return i;
        }
        return 4; // Center
      };

      if (ImGui::BeginTable("##hudLayoutTable", 6,
                            ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit)) {
        ImGui::TableSetupColumn("Widget");
        ImGui::TableSetupColumn("Enabled");
        ImGui::TableSetupColumn("Scale");
        ImGui::TableSetupColumn("Anchor");
        ImGui::TableSetupColumn("X");
        ImGui::TableSetupColumn("Y");
        ImGui::TableHeadersRow();

        auto row = [&](ui::HudWidgetId id, bool* enabledPtr) {
          auto& wLay = hudLayout.widget(id);
          ImGui::TableNextRow();

          ImGui::TableSetColumnIndex(0);
          ImGui::TextUnformatted(ui::toString(id));

          ImGui::TableSetColumnIndex(1);
          if (enabledPtr) {
            bool en = *enabledPtr;
            if (ImGui::Checkbox((std::string("##en_") + ui::toString(id)).c_str(), &en)) {
              *enabledPtr = en;
              wLay.enabled = en;
            }
          } else {
            ImGui::TextDisabled("(auto)");
          }

          ImGui::TableSetColumnIndex(2);
          ImGui::SetNextItemWidth(90.0f);
          ImGui::SliderFloat((std::string("##sc_") + ui::toString(id)).c_str(), &wLay.scale, 0.50f, 2.50f, "%.2f");

          ImGui::TableSetColumnIndex(3);
          int aIdx = anchorIndexFor(wLay.pivotX, wLay.pivotY);
          ImGui::SetNextItemWidth(130.0f);
          if (ImGui::BeginCombo((std::string("##anc_") + ui::toString(id)).c_str(), kAnchors[aIdx].label)) {
            for (int i = 0; i < (int)std::size(kAnchors); ++i) {
              const bool sel = (i == aIdx);
              if (ImGui::Selectable(kAnchors[i].label, sel)) {
                wLay.pivotX = kAnchors[i].px;
                wLay.pivotY = kAnchors[i].py;
              }
              if (sel) ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
          }

          ImGui::TableSetColumnIndex(4);
          ImGui::SetNextItemWidth(90.0f);
          ImGui::DragFloat((std::string("##x_") + ui::toString(id)).c_str(), &wLay.posNormX, 0.0025f, 0.0f, 1.0f, "%.3f");

          ImGui::TableSetColumnIndex(5);
          ImGui::SetNextItemWidth(90.0f);
          ImGui::DragFloat((std::string("##y_") + ui::toString(id)).c_str(), &wLay.posNormY, 0.0025f, 0.0f, 1.0f, "%.3f");
        };

        row(ui::HudWidgetId::Radar, &showRadarHud);
        row(ui::HudWidgetId::Objective, &objectiveHudEnabled);
        row(ui::HudWidgetId::Threat, &hudThreatOverlayEnabled);
        row(ui::HudWidgetId::Jump, &hudJumpOverlay);

        ImGui::EndTable();
      }

      ImGui::TextDisabled("Tip: When Edit mode is ON, you can drag borderless HUD panels from anywhere.");
      ImGui::End();
    }

    // Controls (keybinds) window
    game::drawControlsWindow(controlsWindow, controls, controlsDefaults, controlsDirty, controlsAutoSaveOnExit, controlsPath,
                             [&](const std::string& msg, double ttlSec) { toast(toasts, msg, ttlSec); });

    // Console window (log + commands)
    game::drawConsoleWindow(consoleWindow);

    // HUD overlay: combat symbology + target marker + toasts
    {
      ImDrawList* draw = ImGui::GetForegroundDrawList();
      const ImVec2 center((float)w * 0.5f, (float)h * 0.5f);

      // HUD layout edit guides (safe area + pivot markers)
      if (hudLayoutEditMode) {
        ImGuiViewport* vp = ImGui::GetMainViewport();
        const float m = std::max(0.0f, (float)hudLayout.safeMarginPx);
        const ImVec2 a(vp->WorkPos.x + m, vp->WorkPos.y + m);
        const ImVec2 b(vp->WorkPos.x + vp->WorkSize.x - m, vp->WorkPos.y + vp->WorkSize.y - m);
        const ImU32 col = IM_COL32(255, 210, 120, 150);
        draw->AddRect(a, b, col, 0.0f, 0, 1.0f);
        draw->AddText(ImVec2(a.x + 6.0f, a.y + 6.0f), IM_COL32(255, 210, 120, 220), "HUD EDIT MODE");

        auto pivotPt = [&](ui::HudWidgetId id) -> ImVec2 {
          const auto& wl = hudLayout.widget(id);
          return ImVec2(vp->WorkPos.x + vp->WorkSize.x * wl.posNormX,
                        vp->WorkPos.y + vp->WorkSize.y * wl.posNormY);
        };
        auto drawPivot = [&](ui::HudWidgetId id, bool enabled) {
          if (!enabled) return;
          const ImVec2 p = pivotPt(id);
          const float r = 6.0f;
          draw->AddLine(ImVec2(p.x - r, p.y), ImVec2(p.x + r, p.y), col, 1.0f);
          draw->AddLine(ImVec2(p.x, p.y - r), ImVec2(p.x, p.y + r), col, 1.0f);
          draw->AddText(ImVec2(p.x + r + 3.0f, p.y - r - 2.0f), col, ui::toString(id));
        };

        drawPivot(ui::HudWidgetId::Radar, showRadarHud);
        drawPivot(ui::HudWidgetId::Objective, objectiveHudEnabled);
        drawPivot(ui::HudWidgetId::Threat, hudThreatOverlayEnabled);
        drawPivot(ui::HudWidgetId::Jump, hudJumpOverlay);
      }

      // Shared atlas texture for HUD icons.
      const auto& atlasTex = hudAtlas.texture();
      const ImTextureID atlasId = (ImTextureID)(intptr_t)atlasTex.handle();

      auto clampToEdge = [&](ImVec2 p, float margin) -> ImVec2 {
        ImVec2 dir(p.x - center.x, p.y - center.y);
        const float len = std::sqrt(dir.x * dir.x + dir.y * dir.y);
        if (len < 1.0e-3f) return p;
        dir.x /= len;
        dir.y /= len;

        const float halfW = center.x - margin;
        const float halfH = center.y - margin;
        const float ax = std::abs(dir.x);
        const float ay = std::abs(dir.y);
        const float tx = (ax > 1.0e-4f) ? (halfW / ax) : 1.0e9f;
        const float ty = (ay > 1.0e-4f) ? (halfH / ay) : 1.0e9f;
        const float t = std::min(tx, ty);
        return ImVec2(center.x + dir.x * t, center.y + dir.y * t);
      };

      auto drawArc = [&](float radius, float a0, float a1, ImU32 col, float thickness) {
        const int segments = std::max(12, (int)(radius * 0.35f));
        draw->PathClear();
        draw->PathArcTo(center, radius, a0, a1, segments);
        draw->PathStroke(col, false, thickness);
      };

      const float retA = std::clamp(hudReticleAlpha, 0.0f, 1.0f);

      // Center reticle (procedural sprite) + weapon rings.
      if (hudCombatHud) {
        if (hudUseProceduralReticle) {
          const auto uv = hudAtlas.get(render::SpriteKind::HudReticle,
                                       core::hashCombine(core::fnv1a64("hud_reticle"), seed));
          const float sz = std::max(12.0f, hudReticleSizePx);
          const float oSz = sz * 1.08f;
          const ImVec2 o0(center.x - oSz * 0.5f, center.y - oSz * 0.5f);
          const ImVec2 o1(center.x + oSz * 0.5f, center.y + oSz * 0.5f);
          const ImVec2 p0(center.x - sz * 0.5f, center.y - sz * 0.5f);
          const ImVec2 p1(center.x + sz * 0.5f, center.y + sz * 0.5f);

          draw->AddImage(atlasId, o0, o1, ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1),
                         IM_COL32(0, 0, 0, (int)(100.0f * retA)));
          draw->AddImage(atlasId, p0, p1, ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1),
                         IM_COL32(210, 230, 255, (int)(255.0f * retA)));
        } else {
          const ImU32 col = IM_COL32(160,160,170, (int)(200.0f * retA));
          draw->AddLine({center.x - 8, center.y}, {center.x + 8, center.y}, col, 1.0f);
          draw->AddLine({center.x, center.y - 8}, {center.x, center.y + 8}, col, 1.0f);
        }

        const float baseR = std::max(14.0f, (hudReticleSizePx * 0.5f) + 12.0f);
        const float start = -(float)math::kPi * 0.5f;
        const float full = (float)math::kPi * 2.0f;

        auto clampByte = [](double x) -> int { return (int)std::clamp(x, 0.0, 255.0); };

        auto drawWeaponRing = [&](WeaponType wt, double cd, float radius, float thickness) {
          const WeaponDef& wd = weaponDef(wt);
          if (wd.cooldownSimSec <= 1e-6) return;

          const float prog = (float)std::clamp(1.0 - (cd / wd.cooldownSimSec), 0.0, 1.0);
          const ImU32 bg = IM_COL32(40, 40, 55, (int)(105.0f * retA));
          drawArc(radius, start, start + full, bg, thickness);

          const int r = clampByte((double)wd.r * 255.0);
          const int g = clampByte((double)wd.g * 255.0);
          const int b = clampByte((double)wd.b * 255.0);
          const ImU32 fg = IM_COL32(r, g, b, (int)(210.0f * retA));
          drawArc(radius, start, start + full * prog, fg, thickness);
        };

        if (hudShowWeaponRings) {
          drawWeaponRing(weaponPrimary, weaponPrimaryCooldown, baseR + 4.0f, 3.0f);
          drawWeaponRing(weaponSecondary, weaponSecondaryCooldown, baseR - 2.0f, 3.0f);
        }

        if (hudShowHeatRing) {
          const float heatFrac = (float)std::clamp(heat / 100.0, 0.0, 1.0);
          const ImU32 bg = IM_COL32(35, 35, 45, (int)(95.0f * retA));
          drawArc(baseR + 10.0f, start, start + full, bg, 2.0f);

          const int r = (int)std::clamp(255.0f * heatFrac, 0.0f, 255.0f);
          const int g = (int)std::clamp(255.0f * (1.0f - heatFrac), 0.0f, 255.0f);
          const int b = 60;
          const ImU32 fg = IM_COL32(r, g, b, (int)(190.0f * retA));
          drawArc(baseR + 10.0f, start, start + full * heatFrac, fg, 2.0f);
        }

        // Capacitor rings (ENG/WEP/SYS)
        if (hudShowDistributorRings) {
          auto capFrac = [](double v, double cap) {
            return (float)std::clamp((cap > 1e-9) ? (v / cap) : 1.0, 0.0, 1.0);
          };

          const float engFrac = capFrac(distributorState.eng, distributorCfg.capEng);
          const float wepFrac = capFrac(distributorState.wep, distributorCfg.capWep);
          const float sysFrac = capFrac(distributorState.sys, distributorCfg.capSys);

          const float ringR = baseR + 15.0f;
          const float seg = full / 3.0f;
          const float th = 2.0f;

          const ImU32 bg = IM_COL32(35, 35, 45, (int)(85.0f * retA));

          auto drawSeg = [&](int idx, float frac, ImU32 fg) {
            const float a0 = start + seg * (float)idx;
            const float a1 = a0 + seg;
            drawArc(ringR, a0, a1, bg, th);
            drawArc(ringR, a0, a0 + seg * frac, fg, th);
          };

          drawSeg(0, engFrac, IM_COL32(80, 200, 255, (int)(190.0f * retA))); // ENG
          drawSeg(1, wepFrac, IM_COL32(255, 170, 80, (int)(190.0f * retA))); // WEP
          drawSeg(2, sysFrac, IM_COL32(110, 255, 150, (int)(190.0f * retA))); // SYS
        }

        // Interdiction HUD: escape vector marker + tug-of-war meter (supercruise only).
        if (supercruiseState == SupercruiseState::Active && sim::interdictionInProgress(interdiction)) {
          const float meter = (float)std::clamp(interdiction.meter, 0.0, 1.0);
          const float ringR = baseR + 22.0f;
          const ImU32 bg = IM_COL32(40, 32, 28, (int)(95.0f * retA));
          const ImU32 fg = IM_COL32(255, 150, 100, (int)(205.0f * retA));
          drawArc(ringR, start, start + full, bg, 3.0f);
          drawArc(ringR, start, start + full * meter, fg, 3.0f);

          // Escape marker projected in ship-local space (forward is +Z).
          const math::Vec3d escW = interdiction.escapeDir.normalized();
          const math::Vec3d escL = ship.orientation().conjugate().rotate(escW);
          const double z = std::max(0.18, escL.z);
          double px = escL.x / z;
          double py = escL.y / z;
          const double mag = std::sqrt(px * px + py * py);
          if (mag > 1.0) {
            px /= mag;
            py /= mag;
          }

          const float markerR = ringR - 6.0f;
          const ImVec2 mpos(center.x + (float)px * markerR, center.y - (float)py * markerR);
          const ImU32 stroke = IM_COL32(0, 0, 0, (int)(140.0f * retA));
          draw->AddCircleFilled(mpos, 3.2f, fg);
          draw->AddCircle(mpos, 7.0f, stroke, 0, 1.0f);

          // Brief label during warning to make it obvious *why* the ring appeared.
          if (interdiction.phase == sim::InterdictionPhase::Warning) {
            draw->AddText(ImVec2(center.x + ringR + 10.0f, center.y - 6.0f), fg, "WARNING");
          }
        }

      } else {
        // Minimal crosshair if the combat HUD is disabled.
        draw->AddLine({center.x - 8, center.y}, {center.x + 8, center.y}, IM_COL32(160,160,170,140), 1.0f);
        draw->AddLine({center.x, center.y - 8}, {center.x, center.y + 8}, IM_COL32(160,160,170,140), 1.0f);
      }

      // Flight path marker (velocity vector) - useful with Newtonian drift.
      if (hudCombatHud && hudShowFlightPathMarker && !docked) {
        math::Vec3d v = ship.velocityKmS();
        if (hudFlightMarkerUseLocalFrame) v = v - localFrameVelKmS;
        const double spd = v.length();
        if (spd > 0.05) {
          const math::Vec3d dir = v / spd;
          const math::Vec3d pKm = ship.positionKm() + dir * 120000.0; // direction only (projection stable)
          ImVec2 p{};
          bool off = false;
          if (projectToScreenAny(toRenderU(pKm), view, proj, w, h, p, off)) {
            ImVec2 drawP = p;
            if (off && hudFlightMarkerClampToEdge) {
              drawP = clampToEdge(p, 34.0f);
            }

            const auto uv = hudAtlas.get(render::SpriteKind::HudVelocity,
                                         core::hashCombine(core::fnv1a64("hud_vel"), seed));
            const float sz = std::max(8.0f, hudFlightMarkerSizePx);
            const float oSz = sz * 1.10f;
            const ImVec2 o0(drawP.x - oSz * 0.5f, drawP.y - oSz * 0.5f);
            const ImVec2 o1(drawP.x + oSz * 0.5f, drawP.y + oSz * 0.5f);
            const ImVec2 p0(drawP.x - sz * 0.5f, drawP.y - sz * 0.5f);
            const ImVec2 p1(drawP.x + sz * 0.5f, drawP.y + sz * 0.5f);

            draw->AddImage(atlasId, o0, o1, ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1), IM_COL32(0,0,0,120));
            draw->AddImage(atlasId, p0, p1, ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1), IM_COL32(120, 220, 255, 210));

            if (io.KeyShift) {
              char buf[64];
              std::snprintf(buf, sizeof(buf), "v %.2f km/s", spd);
              draw->AddText(ImVec2(drawP.x + sz * 0.58f, drawP.y - 7.0f), IM_COL32(160, 240, 255, 200), buf);
            }
          }
        }
      }

      // Target marker
      std::optional<math::Vec3d> tgtKm{};
      std::optional<math::Vec3d> tgtVelKmS{};
      std::string tgtLabel;
      std::optional<std::pair<render::SpriteKind, core::u64>> tgtIcon{};

      if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
        const auto& st = currentSystem->stations[target.index];
        tgtKm = stationPosKm(st, timeDays);
        tgtLabel = st.name;
        tgtIcon = {render::SpriteKind::Station, core::hashCombine(core::fnv1a64("station"), st.id)};
      } else if (target.kind == TargetKind::Planet && target.index < currentSystem->planets.size()) {
        const auto& p = currentSystem->planets[target.index];
        tgtKm = sim::orbitPosition3DAU(p.orbit, timeDays) * kAU_KM;
        tgtLabel = p.name;
        tgtIcon = {render::SpriteKind::Planet,
                   core::hashCombine(core::hashCombine(core::fnv1a64("planet"), currentSystem->stub.seed), (core::u64)target.index)};
      } else if (target.kind == TargetKind::Star) {
        tgtKm = math::Vec3d{0,0,0};
        tgtLabel = std::string("Star ") + starClassName(currentSystem->star.cls);
        tgtIcon = {render::SpriteKind::Star, core::hashCombine(core::fnv1a64("star"), currentSystem->stub.seed)};
      } else if (target.kind == TargetKind::Contact && target.index < contacts.size()) {
        const auto& c = contacts[target.index];
        if (c.alive) {
          tgtKm = c.ship.positionKm();
          tgtVelKmS = c.ship.velocityKmS();
          tgtLabel = c.name + std::string(" [") + contactRoleName(c.role) + "]";
          tgtIcon = {render::SpriteKind::Ship, core::hashCombine(core::fnv1a64("ship"), c.id)};
        }
	      } else if (target.kind == TargetKind::Signal && target.index < signals.size()) {
	        const auto& s = signals[target.index];
	        tgtKm = s.posKm;
	        const char* base = signalTypeName(s.type);
	        if (s.type == SignalType::Resource && s.hasResourcePlan) {
	          base = sim::resourceFieldKindName(s.resource.kind);
	        }
	        tgtLabel = std::string(base) + " Signal";
	        tgtIcon = {render::SpriteKind::Signal,
	                   core::hashCombine(core::hashCombine(core::fnv1a64("signal"), s.id), (core::u64)(int)s.type)};
	      } else if (target.kind == TargetKind::Cargo && target.index < floatingCargo.size()) {
	        const auto& pod = floatingCargo[target.index];
	        tgtKm = pod.posKm;
	        const auto& def = econ::commodityDef(pod.commodity);
	        tgtLabel = def.name + std::string(" Pod") + (pod.units >= 1.0 ? (" x" + std::to_string((int)std::round(pod.units))) : "");
	        tgtIcon = {render::SpriteKind::Cargo,
	                   core::hashCombine(core::hashCombine(core::fnv1a64("cargo"), pod.id), (core::u64)pod.commodity)};
	      } else if (target.kind == TargetKind::Asteroid && target.index < asteroids.size()) {
	        const auto& a = asteroids[target.index];
	        tgtKm = a.posKm;
	        const auto& def = econ::commodityDef(a.yield);
	        tgtLabel = std::string("Asteroid [") + def.name + "]";
	        tgtIcon = {render::SpriteKind::Asteroid,
	                   core::hashCombine(core::hashCombine(core::fnv1a64("asteroid"), a.id), (core::u64)a.yield)};
      }

      if (tgtKm) {
        ImVec2 px{};
        if (projectToScreen(toRenderU(*tgtKm), view, proj, w, h, px)) {
          const double distKm = (*tgtKm - ship.positionKm()).length();
          draw->AddCircle({px.x, px.y}, 14.0f, IM_COL32(255,170,80,190), 1.5f);

          const float iconSize = 18.0f;
          float textX = px.x + 18.0f;
          if (tgtIcon) {
            const auto uv = hudAtlas.get(tgtIcon->first, tgtIcon->second);
            const auto& tex = hudAtlas.texture();
            const ImVec2 p0(textX, px.y - iconSize * 0.55f);
            const ImVec2 p1(textX + iconSize, px.y - iconSize * 0.55f + iconSize);
            draw->AddImage((ImTextureID)(intptr_t)tex.handle(), p0, p1,
                           ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1),
                           IM_COL32(255,255,255,235));
            textX = p1.x + 6.0f;
          }

          std::string s = tgtLabel + "  " + std::to_string((int)distKm) + " km";
          draw->AddText({textX, px.y - 8}, IM_COL32(255,210,170,210), s.c_str());

          // Combat HUD: projectile lead indicator (contacts only).
          if (hudCombatHud && hudShowLeadIndicator && tgtVelKmS && target.kind == TargetKind::Contact && !docked) {
            WeaponType wLeadType = hudLeadUseLastFiredWeapon ? (hudLastFiredPrimary ? weaponPrimary : weaponSecondary) : weaponPrimary;

            auto pickLeadWeapon = [&](WeaponType preferred) -> std::optional<WeaponType> {
              const WeaponDef& w0 = weaponDef(preferred);
              if (!w0.beam && w0.projSpeedKmS > 1e-6) return preferred;
              const WeaponType other = (preferred == weaponPrimary) ? weaponSecondary : weaponPrimary;
              const WeaponDef& w1 = weaponDef(other);
              if (!w1.beam && w1.projSpeedKmS > 1e-6) return other;
              return std::nullopt;
            };

            if (auto wTypeOpt = pickLeadWeapon(wLeadType)) {
              wLeadType = *wTypeOpt;
              const WeaponDef& wLead = weaponDef(wLeadType);
              const double s = wLead.projSpeedKmS;

              const double ttl = wLead.rangeKm / std::max(1e-9, s);
              const double maxT = std::min(ttl, hudLeadMaxTimeSec);

              if (auto sol = sim::solveProjectileLead(ship.positionKm(), ship.velocityKmS(),
                                                      *tgtKm, *tgtVelKmS,
                                                      s, maxT)) {
                const double t = sol->tSec;
                const math::Vec3d leadKm = sol->leadPointKm;
                ImVec2 leadPx{};
                if (projectToScreen(toRenderU(leadKm), view, proj, w, h, leadPx)) {
                  const auto uv = hudAtlas.get(render::SpriteKind::HudLead,
                                               core::hashCombine(core::fnv1a64("hud_lead"), (core::u64)(int)wLeadType));

                  auto clampByte = [](double x) -> int { return (int)std::clamp(x, 0.0, 255.0); };
                  const int cr = clampByte((double)wLead.r * 255.0);
                  const int cg = clampByte((double)wLead.g * 255.0);
                  const int cb = clampByte((double)wLead.b * 255.0);

                  // Connector line from target to lead.
                  draw->AddLine(px, leadPx, IM_COL32(cr, cg, cb, 70), 1.0f);

                  const float sz = std::max(8.0f, hudLeadSizePx);
                  const float oSz = sz * 1.12f;
                  const ImVec2 o0(leadPx.x - oSz * 0.5f, leadPx.y - oSz * 0.5f);
                  const ImVec2 o1(leadPx.x + oSz * 0.5f, leadPx.y + oSz * 0.5f);
                  const ImVec2 p0(leadPx.x - sz * 0.5f, leadPx.y - sz * 0.5f);
                  const ImVec2 p1(leadPx.x + sz * 0.5f, leadPx.y + sz * 0.5f);

                  draw->AddImage(atlasId, o0, o1, ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1), IM_COL32(0,0,0,120));
                  draw->AddImage(atlasId, p0, p1, ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1), IM_COL32(cr, cg, cb, 230));

                  if (io.KeyShift) {
                    char buf[32];
                    std::snprintf(buf, sizeof(buf), "%.1fs", t);
                    draw->AddText(ImVec2(leadPx.x + sz * 0.58f, leadPx.y - 7.0f), IM_COL32(cr, cg, cb, 200), buf);
                  }
                }
              }
            }
          }
        } else if (hudOffscreenTargetIndicator) {
          // If the target is off-screen, draw an edge-of-screen indicator so the player
          // can still navigate toward it without opening the tactical overlay.
          bool off = false;
          ImVec2 pAny{};
          if (projectToScreenAny(toRenderU(*tgtKm), view, proj, w, h, pAny, off) && off) {
            const ImVec2 center((float)w * 0.5f, (float)h * 0.5f);
            ImVec2 dir(pAny.x - center.x, pAny.y - center.y);
            const float len = std::sqrt(dir.x * dir.x + dir.y * dir.y);
            if (len > 1.0e-3f) {
              dir.x /= len;
              dir.y /= len;
            } else {
              dir = ImVec2(1.0f, 0.0f);
            }

            const float margin = 34.0f;
            const float halfW = center.x - margin;
            const float halfH = center.y - margin;
            const float ax = std::abs(dir.x);
            const float ay = std::abs(dir.y);
            const float tx = (ax > 1.0e-3f) ? (halfW / ax) : 1.0e9f;
            const float ty = (ay > 1.0e-3f) ? (halfH / ay) : 1.0e9f;
            const float t = std::min(tx, ty);

            const ImVec2 tip(center.x + dir.x * t, center.y + dir.y * t);
            const float arrowLen = 18.0f;
            const float arrowW = 12.0f;
            const ImVec2 base(tip.x - dir.x * arrowLen, tip.y - dir.y * arrowLen);
            const ImVec2 perp(-dir.y, dir.x);
            const ImVec2 p1(base.x + perp.x * (arrowW * 0.5f), base.y + perp.y * (arrowW * 0.5f));
            const ImVec2 p2(base.x - perp.x * (arrowW * 0.5f), base.y - perp.y * (arrowW * 0.5f));

            draw->AddTriangleFilled(tip, p1, p2, IM_COL32(255,170,80,200));
            draw->AddTriangle(tip, p1, p2, IM_COL32(15, 15, 15, 220), 1.2f);

            // Small icon + distance label near the arrow.
            const double distKm = (*tgtKm - ship.positionKm()).length();
            const float iconSz = 20.0f;
            const ImVec2 iconC(base.x - dir.x * (iconSz * 0.55f), base.y - dir.y * (iconSz * 0.55f));
            const ImVec2 i0(iconC.x - iconSz * 0.5f, iconC.y - iconSz * 0.5f);
            const ImVec2 i1(iconC.x + iconSz * 0.5f, iconC.y + iconSz * 0.5f);

            if (tgtIcon) {
              const auto uv = hudAtlas.get(tgtIcon->first, tgtIcon->second);
              const auto& tex = hudAtlas.texture();
              draw->AddImage((ImTextureID)(intptr_t)tex.handle(), i0, i1,
                             ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1),
                             IM_COL32(255,255,255,235));
            }

            const std::string s = tgtLabel + "  " + std::to_string((int)distKm) + " km";
            draw->AddText({iconC.x + 14.0f, iconC.y - 7.0f}, IM_COL32(255,210,170,210), s.c_str());
          }
        }
      }

      // Docking corridor HUD for targeted station
      if (currentSystem && target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
        const auto& st = currentSystem->stations[target.index];

        const math::Vec3d stPos = stationPosKm(st, timeDays);
        const math::Quatd stQ = stationOrient(st, stPos, timeDays);
        const math::Vec3d stV = stationVelKmS(st, timeDays);

        const math::Vec3d relLocal = stQ.conjugate().rotate(ship.positionKm() - stPos);
        const math::Vec3d vRelLocal = stQ.conjugate().rotate(ship.velocityKmS() - stV);

        const double zEntrance = st.radiusKm * 1.10;
        const double zEnd = zEntrance + st.approachLengthKm;

        const double hx = st.approachRadiusKm;
        const double hy = st.approachRadiusKm;

        const bool inZ = (relLocal.z >= zEntrance && relLocal.z <= zEnd);
        const bool inXY = (std::abs(relLocal.x) <= hx && std::abs(relLocal.y) <= hy);
        const bool inCorridor = inZ && inXY;

        bool clearanceGranted = false;
        if (auto it = clearances.find(st.id); it != clearances.end()) {
          clearanceGranted = sim::dockingClearanceValid(it->second, timeDays);
        }

        const math::Vec3d axisOut = stQ.rotate({0,0,1});
        const math::Vec3d slotFwd = (-axisOut).normalized();
        const math::Vec3d slotUp = stQ.rotate({0,1,0}).normalized();

        const double fwdAlign = math::dot(ship.forward().normalized(), slotFwd);
        const double upAlign = math::dot(ship.up().normalized(), slotUp);

        const double relSpeed = (ship.velocityKmS() - stV).length();

        ImVec2 ds = io.DisplaySize;
        float y0 = ds.y - 78.0f;

        char buf[256];
        std::snprintf(buf, sizeof(buf),
                      "Corridor: %s | offX %.0f km offY %.0f km | z %.0f km | vrel %.2f km/s",
                      inCorridor ? "IN" : "OUT",
                      relLocal.x, relLocal.y, relLocal.z, relSpeed);

        const ImU32 col = clearanceGranted ? IM_COL32(160,255,180,220) : IM_COL32(255,220,140,220);
        draw->AddText({18.0f, y0}, col, buf);

        std::snprintf(buf, sizeof(buf),
                      "Align: fwd %.2f up %.2f | speed limit %.2f km/s | %s",
                      fwdAlign, upAlign, st.maxApproachSpeedKmS,
                      clearanceGranted ? "CLEARANCE GRANTED" : "NO CLEARANCE (press L to request)");
        draw->AddText({18.0f, y0 + 18.0f}, col, buf);
      }

      // Tactical overlay: draw projected world icons directly on the screen.
      // This intentionally does not create an ImGui window (so it won't steal clicks
      // from other UI) and uses MMB for targeting to avoid conflicting with weapons.
      if (showTacticalOverlay && currentSystem && !docked) {
        struct TacMarker {
          TargetKind kind;
          std::size_t index;
          ImVec2 px;
          double distKm;
          render::SpriteKind iconKind;
          core::u64 iconSeed;
          std::string label;
        };

        std::vector<TacMarker> markers;
        markers.reserve(128);

        const math::Vec3d shipPosKm = ship.positionKm();

        auto addMarker = [&](TargetKind kind,
                             std::size_t index,
                             const math::Vec3d& posKm,
                             render::SpriteKind iconKind,
                             core::u64 iconSeed,
                             std::string label) {
          const double distKm = (posKm - shipPosKm).length();
          if (distKm > tacticalRangeKm) return;
          ImVec2 px{};
          if (!projectToScreen(toRenderU(posKm), view, proj, w, h, px)) return;
          markers.push_back({kind, index, px, distKm, iconKind, iconSeed, std::move(label)});
        };

        if (tacticalShowStations) {
          for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
            const auto& st = currentSystem->stations[i];
            const math::Vec3d posKm = stationPosKm(st, timeDays);
            addMarker(TargetKind::Station,
                      i,
                      posKm,
                      render::SpriteKind::Station,
                      core::hashCombine(core::fnv1a64("station"), st.id),
                      st.name);
          }
        }

        if (tacticalShowPlanets) {
          for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
            const auto& p = currentSystem->planets[i];
            const math::Vec3d posKm = sim::orbitPosition3DAU(p.orbit, timeDays) * kAU_KM;
            const core::u64 seedP = core::hashCombine(core::hashCombine(core::fnv1a64("planet"), currentSystem->stub.seed), (core::u64)i);
            addMarker(TargetKind::Planet, i, posKm, render::SpriteKind::Planet, seedP, p.name);
          }
        }

        if (tacticalShowContacts) {
          for (std::size_t i = 0; i < contacts.size(); ++i) {
            const auto& c = contacts[i];
            if (!c.alive) continue;
            addMarker(TargetKind::Contact,
                      i,
                      c.ship.positionKm(),
                      render::SpriteKind::Ship,
                      core::hashCombine(core::fnv1a64("ship"), c.id),
                      c.name + std::string(" [") + contactRoleName(c.role) + "]");
          }
        }

        if (tacticalShowCargo) {
          for (std::size_t i = 0; i < floatingCargo.size(); ++i) {
            const auto& pod = floatingCargo[i];
            if (pod.units <= 0.0) continue;
            const auto& def = econ::commodityDef(pod.commodity);
            std::string label = def.name + std::string(" Pod");
            if (pod.units >= 1.0) label += " x" + std::to_string((int)std::round(pod.units));
            addMarker(TargetKind::Cargo,
                      i,
                      pod.posKm,
                      render::SpriteKind::Cargo,
                      core::hashCombine(core::hashCombine(core::fnv1a64("cargo"), pod.id), (core::u64)pod.commodity),
                      std::move(label));
          }
        }

        if (tacticalShowAsteroids) {
          for (std::size_t i = 0; i < asteroids.size(); ++i) {
            const auto& a = asteroids[i];
            if (a.remainingUnits <= 0.0) continue;
            const auto& def = econ::commodityDef(a.yield);
            addMarker(TargetKind::Asteroid,
                      i,
                      a.posKm,
                      render::SpriteKind::Asteroid,
                      core::hashCombine(core::hashCombine(core::fnv1a64("asteroid"), a.id), (core::u64)a.yield),
                      std::string("Asteroid [") + def.name + "]");
          }
        }

        if (tacticalShowSignals) {
          for (std::size_t i = 0; i < signals.size(); ++i) {
            const auto& s = signals[i];
            if (timeDays > s.expireDay) continue;
	        const char* name = signalTypeName(s.type);
	        if (s.type == SignalType::Resource && s.hasResourcePlan) {
	          name = sim::resourceFieldKindName(s.resource.kind);
	        }
            addMarker(TargetKind::Signal,
                      i,
                      s.posKm,
                      render::SpriteKind::Signal,
                      core::hashCombine(core::hashCombine(core::fnv1a64("signal"), s.id), (core::u64)(int)s.type),
	                  std::string(name) + " Signal");
          }
        }

        // Keep the overlay readable: render closest first and cap count.
        std::sort(markers.begin(), markers.end(), [](const TacMarker& a, const TacMarker& b) {
          return a.distKm < b.distKm;
        });
        if ((int)markers.size() > tacticalMaxMarkers) markers.resize((std::size_t)tacticalMaxMarkers);

        auto baseSizeFor = [&](TargetKind k) -> float {
          switch (k) {
            case TargetKind::Station: return 20.0f;
            case TargetKind::Planet: return 20.0f;
            case TargetKind::Contact: return 18.0f;
            case TargetKind::Cargo: return 16.0f;
            case TargetKind::Asteroid: return 16.0f;
            case TargetKind::Signal: return 16.0f;
            default: return 16.0f;
          }
        };

        auto quantSize = [&](float s) -> int {
          // Keep marker sizes stable: quantize to 4px steps to reduce jitter/shimmer.
          int q = (int)std::lround(s / 4.0f) * 4;
          return std::clamp(q, 12, 40);
        };

        // Hover/click selection (MMB) without an ImGui window.
        const ImVec2 mouse = ImGui::GetMousePos();
        int hovered = -1;
        double bestD2 = 1e30;

        if (!io.WantCaptureMouse) {
          for (int i = 0; i < (int)markers.size(); ++i) {
            const auto& m = markers[(std::size_t)i];
            const double frac = std::clamp(1.0 - (m.distKm / std::max(1.0, tacticalRangeKm)), 0.0, 1.0);
            float sz = baseSizeFor(m.kind) * (0.75f + 0.55f * (float)std::sqrt(frac));
            if (m.kind == target.kind && m.index == target.index) sz *= 1.25f;
            const int szPx = quantSize(sz);
            const float rad = (float)szPx * 0.60f;
            const double dx = (double)mouse.x - (double)m.px.x;
            const double dy = (double)mouse.y - (double)m.px.y;
            const double d2 = dx*dx + dy*dy;
            if (d2 <= (double)(rad * rad) && d2 < bestD2) {
              bestD2 = d2;
              hovered = i;
            }
          }
        }

        if (hovered >= 0 && !io.WantCaptureMouse && ImGui::IsMouseClicked(ImGuiMouseButton_Middle)) {
          target.kind = markers[(std::size_t)hovered].kind;
          target.index = markers[(std::size_t)hovered].index;
        }

        // Draw markers
        for (int i = 0; i < (int)markers.size(); ++i) {
          const auto& m = markers[(std::size_t)i];
          const double frac = std::clamp(1.0 - (m.distKm / std::max(1.0, tacticalRangeKm)), 0.0, 1.0);
          float sz = baseSizeFor(m.kind) * (0.75f + 0.55f * (float)std::sqrt(frac));
          if (m.kind == target.kind && m.index == target.index) sz *= 1.25f;
          if (i == hovered) sz *= 1.10f;
          const int szPx = quantSize(sz);
          const float szF = (float)szPx;

          const auto uv = hudAtlas.get(m.iconKind, m.iconSeed);
          const auto& tex = hudAtlas.texture();
          const ImVec2 p0(m.px.x - szF * 0.5f, m.px.y - szF * 0.5f);
          const ImVec2 p1(m.px.x + szF * 0.5f, m.px.y + szF * 0.5f);
          draw->AddImage((ImTextureID)(intptr_t)tex.handle(), p0, p1,
                         ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1),
                         IM_COL32(255,255,255,235));

          if (m.kind == target.kind && m.index == target.index) {
            draw->AddCircle(m.px, szF * 0.62f, IM_COL32(255,170,80,170), 16, 1.5f);
          }
          if (i == hovered && !io.WantCaptureMouse) {
            draw->AddCircle(m.px, szF * 0.62f, IM_COL32(255,255,255,110), 16, 1.25f);
          }

          if (tacticalShowLabels) {
            const std::string label = m.label + "  " + std::to_string((int)std::llround(m.distKm)) + " km";
            draw->AddText({m.px.x + szF * 0.55f + 4.0f, m.px.y - szF * 0.52f}, IM_COL32(230,230,240,190), label.c_str());
          }
        }

        // Tooltip (hover)
        if (hovered >= 0 && !io.WantCaptureMouse) {
          const auto& m = markers[(std::size_t)hovered];
          ImGui::BeginTooltip();
          const auto& big = spriteCache.get(m.iconKind, m.iconSeed, 96);
          ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(64, 64));
          ImGui::SameLine();
          ImGui::BeginGroup();
          ImGui::TextUnformatted(m.label.c_str());
          ImGui::TextDisabled("Distance: %.0f km", m.distKm);
          ImGui::TextDisabled("MMB: target");
          ImGui::EndGroup();
          ImGui::EndTooltip();
        }
      }

      // Toasts
      float y = 18.0f;
      const float x = 18.0f;
      const ImVec2 pad(8.0f, 4.0f);

      for (const auto& t : toasts) {
        const float denom = (float)std::max(0.001, t.ttlTotal);
        const float a = std::clamp((float)(t.ttl / denom), 0.0f, 1.0f);

        const ImVec2 textSz = ImGui::CalcTextSize(t.text.c_str());
        const ImVec2 p0(x, y);
        const ImVec2 r0(p0.x - pad.x, p0.y - pad.y);
        const ImVec2 r1(p0.x + textSz.x + pad.x, p0.y + textSz.y + pad.y);

        const int bgA = (int)std::llround(150.0f * a);
        const int fgA = (int)std::llround(230.0f * a);
        draw->AddRectFilled(r0, r1, IM_COL32(12, 12, 16, bgA), 4.0f);
        draw->AddRect(r0, r1, IM_COL32(255, 255, 255, (int)std::llround(42.0f * a)), 4.0f);
        draw->AddText(p0, IM_COL32(240, 240, 240, fgA), t.text.c_str());

        y += textSz.y + pad.y * 2.0f + 4.0f;
      }
    }

// --- Shared UI helpers (missions / tracker) ---
auto uiSystemNameById = [&](sim::SystemId sysId) -> std::string {
  if (sysId == 0) return "(none)";
  return universe.getSystem(sysId).stub.name;
};

auto uiStationNameById = [&](sim::SystemId sysId, sim::StationId stId) -> std::string {
  if (sysId == 0 || stId == 0) return "(none)";
  const auto& sys = universe.getSystem(sysId);
  for (const auto& st : sys.stations) {
    if (st.id == stId) return st.name;
  }
  return "Station " + std::to_string(stId);
};

auto uiMissionNextStop = [&](const sim::Mission& m) -> std::pair<sim::SystemId, sim::StationId> {
  // Multi-leg deliveries: go to the via stop first.
  if (m.type == sim::MissionType::MultiDelivery && m.viaSystem != 0 && m.viaStation != 0 && m.leg == 0) {
    return {m.viaSystem, m.viaStation};
  }

  // Salvage jobs: before the site is visited, the objective is the mission signal (not a station).
  if (m.type == sim::MissionType::Salvage && !m.scanned) {
    return {m.toSystem, 0};
  }

  // Bounty missions provide a best-effort station hint where the target tends to lurk.
  if (m.type == sim::MissionType::BountyScan || m.type == sim::MissionType::BountyKill) {
    return {m.toSystem, m.toStation};
  }

  return {m.toSystem, m.toStation};
};

auto uiDescribeMission = [&](const sim::Mission& m) -> std::string {
  std::string out;

  if (m.type == sim::MissionType::Courier) {
    out = "Courier: Deliver data to " + uiStationNameById(m.toSystem, m.toStation) + " (" + uiSystemNameById(m.toSystem) + ")";
  } else if (m.type == sim::MissionType::Delivery) {
    out = std::string("Delivery: ")
          + econ::commodityDef(m.commodity).name + " x" + std::to_string((int)std::round(m.units))
          + " to " + uiStationNameById(m.toSystem, m.toStation) + " (" + uiSystemNameById(m.toSystem) + ")";
  } else if (m.type == sim::MissionType::Salvage) {
    out = std::string("Salvage: ")
          + (m.scanned ? "deliver " : "recover ")
          + econ::commodityDef(m.commodity).name + " x" + std::to_string((int)std::round(m.units))
          + " -> " + uiStationNameById(m.toSystem, m.toStation) + " (" + uiSystemNameById(m.toSystem) + ")";
    if (!m.scanned) out += " [visit mission derelict]";
  } else if (m.type == sim::MissionType::Smuggle) {
    out = std::string("Smuggle: ")
          + econ::commodityDef(m.commodity).name + " x" + std::to_string((int)std::round(m.units))
          + " to " + uiStationNameById(m.toSystem, m.toStation) + " (" + uiSystemNameById(m.toSystem) + ") [CONTRABAND]";
  } else if (m.type == sim::MissionType::Escort) {
    out = std::string("Escort: ")
          + (m.scanned ? "report in" : "protect convoy")
          + " ";
    if (m.units > 0.0) {
      out += econ::commodityDef(m.commodity).name + " x" + std::to_string((int)std::round(m.units)) + " -> ";
    }
    out += uiStationNameById(m.toSystem, m.toStation) + " (" + uiSystemNameById(m.toSystem) + ")";
    if (!m.scanned) out += " [stay close]";
  } else if (m.type == sim::MissionType::MultiDelivery) {
    out = std::string("Multi-delivery: ")
          + econ::commodityDef(m.commodity).name + " x" + std::to_string((int)std::round(m.units));
    if (m.viaSystem != 0 && m.viaStation != 0 && m.leg == 0) {
      out += " via " + uiStationNameById(m.viaSystem, m.viaStation) + " (" + uiSystemNameById(m.viaSystem) + ")";
    }
    out += " -> " + uiStationNameById(m.toSystem, m.toStation) + " (" + uiSystemNameById(m.toSystem) + ")";
  } else if (m.type == sim::MissionType::Passenger) {
    out = "Passengers: deliver party of "
          + std::to_string((int)std::round(m.units)) + " to "
          + uiStationNameById(m.toSystem, m.toStation) + " (" + uiSystemNameById(m.toSystem) + ")";
  } else if (m.type == sim::MissionType::BountyScan) {
    out = "Bounty Scan: find & scan target in " + uiSystemNameById(m.toSystem);
    if (m.toStation != 0) {
      out += " (near " + uiStationNameById(m.toSystem, m.toStation) + ")";
    }
  } else if (m.type == sim::MissionType::BountyKill) {
    out = "Bounty Hunt: eliminate target in " + uiSystemNameById(m.toSystem);
    if (m.toStation != 0) {
      out += " (near " + uiStationNameById(m.toSystem, m.toStation) + ")";
    }
  } else {
    out = "Mission";
  }

  if (m.completed) out += " [COMPLETED]";
  if (m.failed) out += " [FAILED]";
  return out;
};



// HUD layout helpers (persistent position for overlay windows)
auto hudSetNextWindowPosFromLayout = [&](ui::HudWidgetId id) {
  ImGuiViewport* vp = ImGui::GetMainViewport();
  const auto& wl = hudLayout.widget(id);
  const ImVec2 pos(vp->WorkPos.x + vp->WorkSize.x * wl.posNormX,
                   vp->WorkPos.y + vp->WorkSize.y * wl.posNormY);
  ImGui::SetNextWindowPos(pos, ImGuiCond_Always, ImVec2(wl.pivotX, wl.pivotY));
};

auto hudCaptureWindowPosToLayout = [&](ui::HudWidgetId id) {
  ImGuiViewport* vp = ImGui::GetMainViewport();
  if (vp->WorkSize.x <= 1.0f || vp->WorkSize.y <= 1.0f) return;

  auto& wl = hudLayout.widget(id);
  const ImVec2 winPos = ImGui::GetWindowPos();
  const ImVec2 winSize = ImGui::GetWindowSize();

  // We store the location of the pivot point, not the top-left.
  const ImVec2 pivotPoint(winPos.x + winSize.x * wl.pivotX,
                          winPos.y + winSize.y * wl.pivotY);

  wl.posNormX = (pivotPoint.x - vp->WorkPos.x) / vp->WorkSize.x;
  wl.posNormY = (pivotPoint.y - vp->WorkPos.y) / vp->WorkSize.y;
  wl.posNormX = std::clamp(wl.posNormX, 0.0f, 1.0f);
  wl.posNormY = std::clamp(wl.posNormY, 0.0f, 1.0f);
};

// Radar HUD overlay (local space awareness)
if (showRadarHud && currentSystem) {
  ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration
                      | ImGuiWindowFlags_NoFocusOnAppearing
                      | ImGuiWindowFlags_NoNav
                      | ImGuiWindowFlags_NoSavedSettings;
#ifdef IMGUI_HAS_DOCK
  flags |= ImGuiWindowFlags_NoDocking;
#endif
  if (!hudLayoutEditMode) flags |= ImGuiWindowFlags_NoBringToFrontOnFocus;

  const auto& wLay = hudLayout.widget(ui::HudWidgetId::Radar);
  const float uiScale = std::max(0.50f, wLay.scale);
  const float baseRad = 94.0f;
  const float rad = baseRad * uiScale;
  const float winSide = (baseRad * 2.0f + 22.0f) * uiScale;
  const ImVec2 winSize(winSide, winSide);

  ImGui::SetNextWindowBgAlpha(hudLayoutEditMode ? 0.45f : 0.35f);
  hudSetNextWindowPosFromLayout(ui::HudWidgetId::Radar);
  ImGui::SetNextWindowSize(winSize, ImGuiCond_Always);

  if (hudLayoutEditMode) {
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 2.0f);
  }

  ImGui::Begin("Radar HUD##radar", nullptr, flags);
  if (hudLayoutEditMode) {
    hudCaptureWindowPosToLayout(ui::HudWidgetId::Radar);
  }

  // Context menu (right click)
  if (ImGui::BeginPopupContextWindow("RadarCtx", ImGuiPopupFlags_MouseButtonRight)) {
    ImGui::TextDisabled("Radar HUD");
    ImGui::Separator();
    ImGui::Checkbox("Enable", &showRadarHud);
    const double minKm = 25000.0;
    const double maxKm = 1200000.0;
    ImGui::DragScalar("Range (km)", ImGuiDataType_Double, &radarRangeKm, 5000.0, &minKm, &maxKm, "%.0f");
    ImGui::SliderInt("Max blips", &radarMaxBlips, 16, 160);
    ImGui::EndPopup();
  }

  ImDrawList* draw = ImGui::GetWindowDrawList();
  const ImVec2 p0 = ImGui::GetCursorScreenPos();
  const ImVec2 canvasSize(rad * 2.0f, rad * 2.0f);
  ImGui::Dummy(canvasSize);
  const ImVec2 p1(p0.x + canvasSize.x, p0.y + canvasSize.y);
  const ImVec2 c((p0.x + p1.x) * 0.5f, (p0.y + p1.y) * 0.5f);

  // Background + rings
  draw->AddCircleFilled(c, rad, IM_COL32(0, 0, 0, 120));
  draw->AddCircle(c, rad, IM_COL32(120, 140, 170, 160), 0, 1.2f * std::max(0.75f, uiScale));
  draw->AddCircle(c, rad * 0.50f, IM_COL32(120, 140, 170, 70), 0, 1.0f);
  draw->AddCircle(c, rad * 0.25f, IM_COL32(120, 140, 170, 50), 0, 1.0f);
  draw->AddLine(ImVec2(c.x - rad, c.y), ImVec2(c.x + rad, c.y), IM_COL32(120, 140, 170, 70), 1.0f * std::max(0.75f, uiScale));
  draw->AddLine(ImVec2(c.x, c.y - rad), ImVec2(c.x, c.y + rad), IM_COL32(120, 140, 170, 70), 1.0f * std::max(0.75f, uiScale));

  // Ship marker (always centered)
  draw->AddTriangleFilled(ImVec2(c.x, c.y - 6 * uiScale),
                          ImVec2(c.x - 5 * uiScale, c.y + 7 * uiScale),
                          ImVec2(c.x + 5 * uiScale, c.y + 7 * uiScale),
                          IM_COL32(230, 230, 240, 210));
  draw->AddTriangle(ImVec2(c.x, c.y - 6 * uiScale),
                    ImVec2(c.x - 5 * uiScale, c.y + 7 * uiScale),
                    ImVec2(c.x + 5 * uiScale, c.y + 7 * uiScale),
                    IM_COL32(20, 20, 20, 220),
                    1.0f * std::max(0.75f, uiScale));

  // Range label
  {
    const std::string rTxt = std::to_string((int)std::round(radarRangeKm)) + " km";
    draw->AddText(ImVec2(p0.x + 6.0f * uiScale, p1.y - 18.0f * uiScale),
                  IM_COL32(200, 210, 230, 170), rTxt.c_str());
  }

  struct Blip {
    TargetKind kind{TargetKind::None};
    std::size_t index{0};
    math::Vec3d posKm{0,0,0};
    double distKm{0.0};
    bool force{false};
  };

  const math::Vec3d shipPos = ship.positionKm();
  const math::Quatd shipInv = ship.orientation().conjugate();
  const double rangeKm = std::max(1.0, radarRangeKm);
  const float s = (float)(rad / rangeKm);

  std::vector<Blip> blips;
  blips.reserve(256);

  auto addBlip = [&](TargetKind kind, std::size_t index, const math::Vec3d& posKm, bool force) {
    const double d = (posKm - shipPos).length();
    if (force || d <= rangeKm) {
      blips.push_back(Blip{kind, index, posKm, d, force});
    }
  };

  // Always include current target (clamped to edge if out of range)
  if (target.kind != TargetKind::None) {
    if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
      addBlip(TargetKind::Station, target.index, stationPosKm(currentSystem->stations[target.index], timeDays), true);
    } else if (target.kind == TargetKind::Planet && target.index < currentSystem->planets.size()) {
      addBlip(TargetKind::Planet, target.index, sim::orbitPosition3DAU(currentSystem->planets[target.index].orbit, timeDays) * kAU_KM, true);
    } else if (target.kind == TargetKind::Contact && target.index < contacts.size() && contacts[target.index].alive) {
      addBlip(TargetKind::Contact, target.index, contacts[target.index].ship.positionKm(), true);
    } else if (target.kind == TargetKind::Signal && target.index < signals.size()) {
      addBlip(TargetKind::Signal, target.index, signals[target.index].posKm, true);
    } else if (target.kind == TargetKind::Cargo && target.index < floatingCargo.size()) {
      addBlip(TargetKind::Cargo, target.index, floatingCargo[target.index].posKm, true);
    } else if (target.kind == TargetKind::Asteroid && target.index < asteroids.size()) {
      addBlip(TargetKind::Asteroid, target.index, asteroids[target.index].posKm, true);
    }
  }

  // Contacts
  for (std::size_t i = 0; i < contacts.size(); ++i) {
    const auto& ctc = contacts[i];
    if (!ctc.alive) continue;
    addBlip(TargetKind::Contact, i, ctc.ship.positionKm(), false);
  }

  // Signals, floating cargo, asteroids
  for (std::size_t i = 0; i < signals.size(); ++i) addBlip(TargetKind::Signal, i, signals[i].posKm, false);
  for (std::size_t i = 0; i < floatingCargo.size(); ++i) addBlip(TargetKind::Cargo, i, floatingCargo[i].posKm, false);
  for (std::size_t i = 0; i < asteroids.size(); ++i) addBlip(TargetKind::Asteroid, i, asteroids[i].posKm, false);

  // Nearby stations/planets (only if they actually fall within the radar range).
  for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
    addBlip(TargetKind::Station, i, stationPosKm(currentSystem->stations[i], timeDays), false);
  }
  for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
    addBlip(TargetKind::Planet, i, sim::orbitPosition3DAU(currentSystem->planets[i].orbit, timeDays) * kAU_KM, false);
  }

  // Sort by (force, distance) and cap count.
  std::sort(blips.begin(), blips.end(), [](const Blip& a, const Blip& b) {
    if (a.force != b.force) return a.force > b.force;
    return a.distKm < b.distKm;
  });
  if ((int)blips.size() > radarMaxBlips) blips.resize((std::size_t)radarMaxBlips);

  // Draw + interaction.
  // Preserve cursor position after our absolute-position buttons.
  const ImVec2 cursorAfter = ImGui::GetCursorScreenPos();

  for (std::size_t bi = 0; bi < blips.size(); ++bi) {
    const Blip& b = blips[bi];
    const math::Vec3d relKm = b.posKm - shipPos;
    const math::Vec3d local = shipInv.rotate(relKm);

    float dx = (float)(local.x * (double)s);
    float dy = (float)(-local.z * (double)s);
    const float r2 = dx*dx + dy*dy;
    const float rMax = rad * 0.92f;
    if (r2 > rMax * rMax) {
      const float r = std::sqrt(r2);
      const float k = (r > 1e-5f) ? (rMax / r) : 1.0f;
      dx *= k;
      dy *= k;
    }

    const ImVec2 bp(c.x + dx, c.y + dy);

    // Pick icon kind/seed.
    render::SpriteKind iKind = render::SpriteKind::Signal;
    core::u64 iSeed = 0;
    float iconSz = (b.force ? 22.0f : 16.0f) * uiScale;
    switch (b.kind) {
      case TargetKind::Station: {
        const auto& st = currentSystem->stations[b.index];
        iKind = render::SpriteKind::Station;
        iSeed = core::hashCombine(core::fnv1a64("station"), st.id);
        iconSz = (b.force ? 24.0f : 18.0f) * uiScale;
      } break;
      case TargetKind::Planet: {
        iKind = render::SpriteKind::Planet;
        iSeed = core::hashCombine(core::hashCombine(core::fnv1a64("planet"), currentSystem->stub.seed), (core::u64)b.index);
        iconSz = (b.force ? 26.0f : 20.0f) * uiScale;
      } break;
      case TargetKind::Contact: {
        const auto& ctc = contacts[b.index];
        iKind = render::SpriteKind::Ship;
        iSeed = core::hashCombine(core::fnv1a64("ship"), ctc.id);
        iconSz = (b.force ? 22.0f : 16.0f) * uiScale;
      } break;
      case TargetKind::Cargo: {
        const auto& pod = floatingCargo[b.index];
        iKind = render::SpriteKind::Cargo;
        iSeed = core::hashCombine(core::hashCombine(core::fnv1a64("cargo"), pod.id), (core::u64)pod.commodity);
        iconSz = (b.force ? 22.0f : 15.0f) * uiScale;
      } break;
      case TargetKind::Asteroid: {
        const auto& a = asteroids[b.index];
        iKind = render::SpriteKind::Asteroid;
        iSeed = core::hashCombine(core::hashCombine(core::fnv1a64("asteroid"), a.id), (core::u64)a.yield);
        iconSz = (b.force ? 22.0f : 15.0f) * uiScale;
      } break;
      case TargetKind::Signal: {
        const auto& ssrc = signals[b.index];
        iKind = render::SpriteKind::Signal;
        iSeed = core::hashCombine(core::hashCombine(core::fnv1a64("signal"), ssrc.id), (core::u64)(int)ssrc.type);
        iconSz = (b.force ? 22.0f : 15.0f) * uiScale;
      } break;
      default: break;
    }

    const ImVec2 b0(bp.x - iconSz * 0.5f, bp.y - iconSz * 0.5f);
    const ImVec2 b1(bp.x + iconSz * 0.5f, bp.y + iconSz * 0.5f);

    // Interactive hitbox (uses ImGui's ID system for hover/click + tooltips)
    ImGui::SetCursorScreenPos(b0);
    ImGui::PushID((int)bi);
    ImGui::InvisibleButton("blip", ImVec2(iconSz, iconSz));

    const bool hovered = ImGui::IsItemHovered();
    if (hovered) {
      ImGui::BeginTooltip();
      // Show a larger icon and a short label.
      const auto& big = spriteCache.get(iKind, iSeed, 96);
      ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(64, 64));

      if (b.kind == TargetKind::Station && b.index < currentSystem->stations.size()) {
        ImGui::SameLine();
        ImGui::Text("%s", currentSystem->stations[b.index].name.c_str());
      } else if (b.kind == TargetKind::Planet && b.index < currentSystem->planets.size()) {
        ImGui::SameLine();
        ImGui::Text("%s", currentSystem->planets[b.index].name.c_str());
      } else if (b.kind == TargetKind::Contact && b.index < contacts.size()) {
        ImGui::SameLine();
        const auto& ctc = contacts[b.index];
        ImGui::Text("%s [%s]", ctc.name.c_str(), contactRoleName(ctc.role));
      } else if (b.kind == TargetKind::Signal && b.index < signals.size()) {
        ImGui::SameLine();
        ImGui::Text("%s Signal", signalTypeName(signals[b.index].type));
      } else if (b.kind == TargetKind::Cargo && b.index < floatingCargo.size()) {
        ImGui::SameLine();
        const auto& pod = floatingCargo[b.index];
        const auto& def = econ::commodityDef(pod.commodity);
        ImGui::Text("%s Pod", def.name);
      } else if (b.kind == TargetKind::Asteroid && b.index < asteroids.size()) {
        ImGui::SameLine();
        const auto& a = asteroids[b.index];
        const auto& def = econ::commodityDef(a.yield);
        ImGui::Text("Asteroid [%s]", def.name);
      }

      ImGui::TextDisabled("Dist %.0f km | RelY %.0f km", b.distKm, local.y);
      ImGui::EndTooltip();
    }

    if (ImGui::IsItemClicked(ImGuiMouseButton_Left) && b.kind != TargetKind::None) {
      target.kind = b.kind;
      target.index = b.index;
    }
    ImGui::PopID();

    // Draw icon (atlas-batched)
    const auto uv = hudAtlas.get(iKind, iSeed);
    const auto& tex = hudAtlas.texture();
    const ImU32 tint = b.force ? IM_COL32(255, 255, 255, 255) : IM_COL32(255, 255, 255, 215);
    draw->AddImage((ImTextureID)(intptr_t)tex.handle(), b0, b1,
                   ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1),
                   tint);

    if (b.force) {
      draw->AddCircle(bp,
                      iconSz * 0.65f + 5.0f * uiScale,
                      IM_COL32(255, 180, 90, 210),
                      0,
                      1.5f * std::max(0.75f, uiScale));
    }

    // Vertical hint for off-plane objects (up/down)
    if (std::abs(local.y) > 2500.0) {
      const char* sym = (local.y > 0.0) ? "^" : "v";
      draw->AddText(ImVec2(bp.x + iconSz * 0.30f, bp.y - iconSz * 0.65f), IM_COL32(230, 230, 240, 190), sym);
    }
  }

  ImGui::SetCursorScreenPos(cursorAfter);
  ImGui::End();
  if (hudLayoutEditMode) ImGui::PopStyleVar();
}


// Objective HUD overlay (tracked mission summary)
if (objectiveHudEnabled) {
  sim::Mission* tracked = nullptr;

  if (trackedMissionId != 0) {
    for (auto& m : missions) {
      if (m.id == trackedMissionId) { tracked = &m; break; }
    }
    if (tracked && (tracked->completed || tracked->failed)) tracked = nullptr;
  }

  if (!tracked) {
    for (auto& m : missions) {
      if (m.completed || m.failed) continue;
      trackedMissionId = m.id;
      tracked = &m;
      break;
    }
  }

  if (tracked && currentSystem) {
    const auto dest = uiMissionNextStop(*tracked);
    const sim::SystemId destSys = dest.first;
    const sim::StationId destSt = dest.second;

    ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration
                        | ImGuiWindowFlags_AlwaysAutoResize
                        | ImGuiWindowFlags_NoFocusOnAppearing
                        | ImGuiWindowFlags_NoNav
                        | ImGuiWindowFlags_NoSavedSettings;

    auto& wLay = hudLayout.widget(ui::HudWidgetId::Objective);
    const float uiScale = std::max(0.50f, wLay.scale);

    ImGui::SetNextWindowBgAlpha(hudLayoutEditMode ? 0.45f : 0.35f);
    hudSetNextWindowPosFromLayout(ui::HudWidgetId::Objective);
    if (hudLayoutEditMode) {
      ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 2.0f);
    }
    ImGui::Begin("Objective HUD##objective", nullptr, flags);
    ImGui::SetWindowFontScale(uiScale);
    if (hudLayoutEditMode) {
      hudCaptureWindowPosToLayout(ui::HudWidgetId::Objective);
    }

    ImGui::Text("Objective");
    ImGui::Separator();
    ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + 360.0f * uiScale);
    ImGui::TextWrapped("%s", uiDescribeMission(*tracked).c_str());
    ImGui::PopTextWrapPos();

    if (tracked->deadlineDay > 0.0 && !tracked->completed && !tracked->failed) {
      const double hrsLeft = (tracked->deadlineDay - timeDays) * 24.0;
      ImGui::TextDisabled("Reward %.0f cr | Deadline %.1f h", tracked->reward, hrsLeft);
    } else {
      ImGui::TextDisabled("Reward %.0f cr", tracked->reward);
    }

    // Escort missions: surface convoy proximity + status to reduce "where did it go?" friction.
    if (tracked->type == sim::MissionType::Escort && !tracked->completed && !tracked->failed) {
      const Contact* convoy = nullptr;
      for (const auto& c : contacts) {
        if (c.alive && c.id == tracked->targetNpcId) { convoy = &c; break; }
      }
      if (convoy) {
        const double distKm = (convoy->ship.positionKm() - ship.positionKm()).length();
        ImGui::Text("Convoy: %s (%.0f km)", convoy->name.c_str(), distKm);
        if (!tracked->scanned && supercruiseState == SupercruiseState::Idle) {
          if (distKm > 180000.0) {
            ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f), "WARNING: Too far from convoy!");
          } else {
            ImGui::TextDisabled("Stay within 180,000 km to avoid losing the convoy.");
          }
        }
      } else {
        ImGui::Text("Convoy: (not in sensors)");
      }

      if (tracked->scanned) {
        ImGui::TextDisabled("Convoy has arrived. Dock at the destination station to complete.");
      }
    }

    // Next stop hint
    if (destSys != 0) {
      if (tracked->type == sim::MissionType::Salvage && !tracked->scanned && destSys == currentSystem->stub.id) {
        // For salvage jobs, the next objective is a mission derelict signal (not a station).
        const SignalSource* ms = nullptr;
        for (const auto& s : signals) {
          if (s.id == tracked->targetNpcId) { ms = &s; break; }
        }
        if (ms) {
          const double distKm = (ms->posKm - ship.positionKm()).length();
          ImGui::Text("Next: Mission derelict (%.0f km)", distKm);
        } else {
          ImGui::Text("Next: Mission derelict (not detected)");
        }
      } else if (destSys == currentSystem->stub.id) {
        if (destSt != 0) {
          const sim::Station* st = nullptr;
          for (const auto& s : currentSystem->stations) {
            if (s.id == destSt) { st = &s; break; }
          }
          if (st) {
            const double distKm = (stationPosKm(*st, timeDays) - ship.positionKm()).length();
            ImGui::Text("Next: %s (%.0f km)", st->name.c_str(), distKm);
          } else {
            ImGui::Text("Next: %s (this system)", uiStationNameById(destSys, destSt).c_str());
          }
        } else {
          ImGui::Text("Next: %s (this system)", uiSystemNameById(destSys).c_str());
        }
      } else {
        const auto& destStub = universe.getSystem(destSys).stub;
        const double distLy = (destStub.posLy - currentSystem->stub.posLy).length();
        ImGui::Text("Next: %s (%.0f ly)", destStub.name.c_str(), distLy);
      }
    }

    // Route summary (if a route is plotted)
    if (navRoute.size() >= 2) {
      const int totalJumps = (int)navRoute.size() - 1;
      const int remaining = std::max(0, totalJumps - (int)navRouteHop);

      double remFuel = 0.0;
      for (std::size_t i = navRouteHop; i + 1 < navRoute.size(); ++i) {
        const auto& a = universe.getSystem(navRoute[i]).stub;
        const auto& b = universe.getSystem(navRoute[i + 1]).stub;
        const double d = (b.posLy - a.posLy).length();
        remFuel += fsdFuelCostFor(d);
      }

      std::string suffix;
      if (destSys != 0 && navRoute.back() == destSys) suffix = " (mission)";
      ImGui::TextDisabled("Route: %d jumps left | est fuel %.1f%s", remaining, remFuel, suffix.c_str());
      if (navAutoRun) ImGui::TextDisabled("Auto-run: ON");
    } else if (destSys != 0 && currentSystem && destSys != currentSystem->stub.id) {
      ImGui::TextDisabled("No route plotted (Plot Route in Ship/Status).");
    }

    ImGui::End();
    if (hudLayoutEditMode) ImGui::PopStyleVar();
  }
}

// Pirate demand HUD (threat/tribute)
if (hudThreatOverlayEnabled && pirateDemand.active && pirateDemand.requiredValueCr > 0.0
    && !docked && fsdState == FsdState::Idle && supercruiseState == SupercruiseState::Idle) {

  ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration
                      | ImGuiWindowFlags_AlwaysAutoResize
                      | ImGuiWindowFlags_NoFocusOnAppearing
                      | ImGuiWindowFlags_NoNav
                      | ImGuiWindowFlags_NoSavedSettings;

  auto& wLay = hudLayout.widget(ui::HudWidgetId::Threat);
  const float uiScale = std::max(0.50f, wLay.scale);

  ImGui::SetNextWindowBgAlpha(hudLayoutEditMode ? 0.45f : 0.35f);
  hudSetNextWindowPosFromLayout(ui::HudWidgetId::Threat);
  if (hudLayoutEditMode) {
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 2.0f);
  }
  ImGui::Begin("Threat HUD##pirate_demand", nullptr, flags);
  ImGui::SetWindowFontScale(uiScale);
  if (hudLayoutEditMode) {
    hudCaptureWindowPosToLayout(ui::HudWidgetId::Threat);
  }

  const std::string src = pirateDemand.leaderName.empty() ? std::string("Pirates") : pirateDemand.leaderName;
  const double need = std::max(0.0, pirateDemand.requiredValueCr);
  const double have = std::max(0.0, pirateDemand.deliveredValueCr);
  const double leftSec = std::max(0.0, (pirateDemand.untilDays - timeDays) * 86400.0);

  ImGui::Text("Threat");
  ImGui::Separator();
  ImGui::TextWrapped("%s demand tribute!", src.c_str());

  const float frac = (need > 1e-6) ? (float)std::clamp(have / need, 0.0, 1.0) : 0.0f;
  ImGui::ProgressBar(frac, ImVec2(240.0f * uiScale, 0.0f));
  ImGui::TextDisabled("Drop cargo value: %.0f / %.0f cr | %.0f s left", have, need, leftSec);
  ImGui::TextDisabled("Jettison: Ship/Status -> Cargo management");
  ImGui::End();
  if (hudLayoutEditMode) ImGui::PopStyleVar();
}


if (hudJumpOverlay && fsdState != FsdState::Idle) {
  ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration
                      | ImGuiWindowFlags_AlwaysAutoResize
                      | ImGuiWindowFlags_NoFocusOnAppearing
                      | ImGuiWindowFlags_NoNav
                      | ImGuiWindowFlags_NoSavedSettings;

  auto& wLay = hudLayout.widget(ui::HudWidgetId::Jump);
  const float uiScale = std::max(0.50f, wLay.scale);

  ImGui::SetNextWindowBgAlpha(hudLayoutEditMode ? 0.45f : 0.28f);
  hudSetNextWindowPosFromLayout(ui::HudWidgetId::Jump);
  if (hudLayoutEditMode) {
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 2.0f);
  }
  ImGui::Begin("FSD Jump##jump_hud", nullptr, flags);
  ImGui::SetWindowFontScale(uiScale);
  if (hudLayoutEditMode) {
    hudCaptureWindowPosToLayout(ui::HudWidgetId::Jump);
  }

  std::string destName = "(none)";
  if (fsdTargetSystem != 0) {
    destName = universe.getSystem(fsdTargetSystem).stub.name;
  }

  if (fsdState == FsdState::Charging) {
    const float t = (float)std::clamp(1.0 - (fsdChargeRemainingSec / std::max(0.001, kFsdChargeSec)), 0.0, 1.0);
    ImGui::Text("FSD CHARGING");
    ImGui::TextDisabled("Destination: %s | %.1f ly", destName.c_str(), fsdJumpDistanceLy);
    ImGui::ProgressBar(t, ImVec2(320.0f * uiScale, 0.0f));
    ImGui::TextDisabled("Press J again to cancel before fuel is consumed.");
  } else if (fsdState == FsdState::Jumping) {
    const float t = (float)std::clamp(1.0 - (fsdTravelRemainingSec / std::max(0.001, fsdTravelTotalSec)), 0.0, 1.0);
    const double eta = std::max(0.0, fsdTravelRemainingSec);
    ImGui::Text("HYPERSPACE");
    ImGui::TextDisabled("Destination: %s", destName.c_str());
    ImGui::ProgressBar(t, ImVec2(320.0f * uiScale, 0.0f));
    ImGui::TextDisabled("ETA: %.1fs", eta);
  }

  ImGui::End();
  if (hudLayoutEditMode) ImGui::PopStyleVar();
}


if (showVfx) {
  ImGui::Begin("VFX Lab");

  ImGui::TextDisabled(
      "%s toggles this panel. Background stars + particles are rendered as point sprites.",
      game::chordLabel(controls.actions.toggleVfxLab).c_str());

  if (ImGui::CollapsingHeader("Starfield", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Checkbox("Enabled##stars", &vfxStarfieldEnabled);

    if (vfxStarfieldEnabled) {
      ImGui::SameLine();
      ImGui::Checkbox("Textured", &vfxStarfieldTextured);
    }

    if (vfxStarfieldEnabled) {
      ImGui::SliderInt("Star count", &vfxStarCount, 0, 20000);
      float radiusF = (float)vfxStarRadiusU;
      if (ImGui::SliderFloat("Star radius (U)", &radiusF, 2000.0f, 20000.0f, "%.0f")) {
        vfxStarRadiusU = (double)radiusF;
      }

      ImGui::Text("Generated stars: %zu", starfield.starCount());
    } else {
      ImGui::TextDisabled("(starfield disabled)");
    }
  }



  if (ImGui::CollapsingHeader("Nebula", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Checkbox("Enabled##nebula", &vfxNebulaEnabled);

    if (vfxNebulaEnabled) {
      ImGui::SliderInt("Puff count", &vfxNebulaPuffCount, 0, 6000);
      ImGui::SliderInt("Variant##nebula", &vfxNebulaVariant, 0, 50);
      ImGui::SliderFloat("Band power##nebula", &vfxNebulaBandPower, 1.0f, 4.0f, "%.2f");

      float innerF = (float)vfxNebulaInnerRadiusU;
      float outerF = (float)vfxNebulaOuterRadiusU;
      if (ImGui::SliderFloat("Inner radius (U)", &innerF, 2000.0f, 40000.0f, "%.0f")) {
        vfxNebulaInnerRadiusU = (double)innerF;
      }
      if (ImGui::SliderFloat("Outer radius (U)", &outerF, 4000.0f, 80000.0f, "%.0f")) {
        vfxNebulaOuterRadiusU = (double)outerF;
      }

      float par = (float)vfxNebulaParallax;
      if (ImGui::SliderFloat("Parallax", &par, 0.0f, 1.0f, "%.2f")) {
        vfxNebulaParallax = (double)par;
      }

      ImGui::SliderFloat("Intensity##nebula", &vfxNebulaIntensity, 0.0f, 6.0f, "%.2f");
      ImGui::SliderFloat("Opacity##nebula", &vfxNebulaOpacity, 0.0f, 0.60f, "%.3f");

      ImGui::SliderFloat("Size min (px)##nebula", &vfxNebulaSizeMinPx, 8.0f, 400.0f, "%.0f");
      ImGui::SliderFloat("Size max (px)##nebula", &vfxNebulaSizeMaxPx, 8.0f, 700.0f, "%.0f");
      if (vfxNebulaSizeMaxPx < vfxNebulaSizeMinPx) vfxNebulaSizeMaxPx = vfxNebulaSizeMinPx;

      ImGui::SliderFloat("Turbulence##nebula", &vfxNebulaTurbulence, 0.0f, 1.0f, "%.2f");
      ImGui::SliderFloat("Turbulence speed##nebula", &vfxNebulaTurbulenceSpeed, 0.0f, 0.50f, "%.3f");

      ImGui::TextDisabled("Puffs generated: %zu", nebula.puffCount());
      ImGui::TextDisabled("Tip: Nebula is additive; tune PostFX bloom/warp for extra glow.");
    }
  }
  if (ImGui::CollapsingHeader("Particles", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Checkbox("Enabled##particles", &vfxParticlesEnabled);

    if (vfxParticlesEnabled) {
      ImGui::SameLine();
      ImGui::Checkbox("Textured", &vfxParticlesTextured);
    }

    if (vfxParticlesEnabled) {
      ImGui::Text("Alive: %zu / %zu", particles.aliveCount(), particles.maxParticles());

      float inten = vfxParticleIntensity;
      if (ImGui::SliderFloat("Intensity", &inten, 0.0f, 2.0f, "%.2f")) {
        vfxParticleIntensity = inten;
      }

      ImGui::Checkbox("Thrusters", &vfxThrustersEnabled);
      ImGui::Checkbox("Impacts", &vfxImpactsEnabled);
      ImGui::Checkbox("Explosions", &vfxExplosionsEnabled);

      static int maxP = 14000;
      if (maxP < 256) maxP = 256;
      if (ImGui::SliderInt("Max particles", &maxP, 256, 30000)) {
        particles.setMaxParticles((std::size_t)maxP);
      }

      if (ImGui::Button("Clear particles")) {
        particles.clear();
      }

      ImGui::SameLine();
      if (ImGui::Button("Test explosion")) {
        const double eBase = 1.75 * (double)vfxParticleIntensity;
        particles.spawnExplosion(toRenderU(ship.positionKm()), toRenderU(ship.velocityKmS()), eBase);
      }
    } else {
      ImGui::TextDisabled("(particles disabled)");
    }
  }

  ImGui::End();
}


if (showPostFx) {
  ImGui::Begin("Post FX");

  ImGui::TextDisabled(
      "%s toggles this panel. Scene renders to an HDR target and is post-processed before UI.",
      game::chordLabel(controls.actions.togglePostFx).c_str());

  ImGui::Checkbox("Enabled", &postFxSettings.enabled);

  if (postFxSettings.enabled) {
    if (ImGui::CollapsingHeader("Bloom", ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::Checkbox("Bloom enabled", &postFxSettings.bloomEnabled);
      ImGui::SliderFloat("Threshold", &postFxSettings.bloomThreshold, 0.10f, 3.50f, "%.2f");
      ImGui::SliderFloat("Soft knee", &postFxSettings.bloomKnee, 0.00f, 2.50f, "%.2f");
      ImGui::SliderFloat("Boost", &postFxSettings.bloomBoost, 0.00f, 3.00f, "%.2f");
      ImGui::SliderFloat("Intensity", &postFxSettings.bloomIntensity, 0.00f, 2.50f, "%.2f");
      ImGui::SliderInt("Blur passes", &postFxSettings.bloomPasses, 0, 16);
      ImGui::TextDisabled("Bloom is extracted at half resolution, then ping-pong blurred.");
    }

    if (ImGui::CollapsingHeader("Tonemap / Output", ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::SliderFloat("Exposure", &postFxSettings.exposure, 0.10f, 3.50f, "%.2f");
      ImGui::SliderFloat("Gamma", &postFxSettings.gamma, 1.20f, 2.80f, "%.2f");
    }

    if (ImGui::CollapsingHeader("Extras", ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::SliderFloat("Vignette", &postFxSettings.vignette, 0.00f, 0.75f, "%.2f");
      ImGui::SliderFloat("Film grain", &postFxSettings.grain, 0.00f, 0.20f, "%.3f");
      ImGui::SliderFloat("Chromatic aberration", &postFxSettings.chromaticAberration, 0.0f, 0.01f, "%.4f");
      ImGui::Separator();
      ImGui::Checkbox("Auto warp from speed", &postFxAutoWarpFromSpeed);
      ImGui::SliderFloat("Warp (manual)", &postFxSettings.warp, 0.0f, 0.08f, "%.4f");
    }

    if (ImGui::CollapsingHeader("Hyperspace / Jump FX", ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::Checkbox("Auto hyperspace from FSD", &postFxAutoHyperspaceFromFsd);
      ImGui::SliderFloat("FSD warp boost", &postFxFsdWarpBoost, 0.0f, 0.12f, "%.3f");
      ImGui::SliderFloat("FSD hyperspace boost", &postFxFsdHyperspaceBoost, 0.0f, 1.5f, "%.2f");
      ImGui::Checkbox("Jump overlay HUD", &hudJumpOverlay);
      ImGui::Separator();
      ImGui::SliderFloat("Hyperspace (manual)", &postFxSettings.hyperspace, 0.0f, 1.0f, "%.2f");
      ImGui::SliderFloat("Twist", &postFxSettings.hyperspaceTwist, 0.0f, 1.0f, "%.2f");
      ImGui::SliderFloat("Density", &postFxSettings.hyperspaceDensity, 0.0f, 1.0f, "%.2f");
      ImGui::SliderFloat("Noise", &postFxSettings.hyperspaceNoise, 0.0f, 1.0f, "%.2f");
      ImGui::SliderFloat("Intensity", &postFxSettings.hyperspaceIntensity, 0.0f, 3.0f, "%.2f");
      ImGui::TextDisabled("Tip: During FSD charge, press J again to cancel before fuel is consumed.");
    }
  } else {
    ImGui::TextDisabled("(PostFX disabled)");
  }

  ImGui::End();
}



if (showShip) {
  ImGui::Begin("Ship / Status");

  ImGui::Text("System: %s", currentSystem->stub.name.c_str());

  const core::u32 localFaction = currentSystem ? currentSystem->stub.factionId : 0;
  if (localFaction != 0) {
    const core::u64 fSeed = core::hashCombine(core::fnv1a64("faction"), (core::u64)localFaction);
    const auto& fIcon = spriteCache.get(render::SpriteKind::Faction, fSeed, 32);
    ImGui::Image((ImTextureID)(intptr_t)fIcon.handle(), ImVec2(18, 18));
    if (ImGui::IsItemHovered()) {
      const auto& big = spriteCache.get(render::SpriteKind::Faction, fSeed, 96);
      ImGui::BeginTooltip();
      ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(96, 96));
      ImGui::EndTooltip();
    }
    ImGui::SameLine();
    ImGui::Text("Local faction: %s | Rep %.1f | Bounty %.0f | Alert %.1f",
                factionName(localFaction).c_str(),
                getRep(localFaction),
                getBounty(localFaction),
                policeHeat);
  } else {
    ImGui::Text("Local faction: (none)");
  }

  // Mission tracker (quality-of-life HUD)
  {
    ImGui::Separator();
    ImGui::Text("Mission tracker");

    sim::Mission* tracked = nullptr;

    if (trackedMissionId != 0) {
      for (auto& m : missions) {
        if (m.id == trackedMissionId) { tracked = &m; break; }
      }
      if (tracked && (tracked->completed || tracked->failed)) tracked = nullptr;
    }

    if (!tracked) {
      for (auto& m : missions) {
        if (m.completed || m.failed) continue;
        trackedMissionId = m.id;
        tracked = &m;
        break;
      }
    }

    if (tracked) {
      ImGui::TextWrapped("%s", uiDescribeMission(*tracked).c_str());

      if (tracked->deadlineDay > 0.0 && !tracked->completed && !tracked->failed) {
        const double hrsLeft = (tracked->deadlineDay - timeDays) * 24.0;
        ImGui::TextDisabled("Reward %.0f cr | Deadline in %.1f h", tracked->reward, hrsLeft);
      } else {
        ImGui::TextDisabled("Reward %.0f cr", tracked->reward);
      }

      // Progress hints
      if (tracked->type == sim::MissionType::Salvage) {
        ImGui::TextDisabled("Site: %s", tracked->scanned ? "found" : "not found");

        const std::size_t ci = (std::size_t)tracked->commodity;
        if (ci < econ::kCommodityCount) {
          const double have = cargo[ci];
          ImGui::TextDisabled("Cargo: %.0f / %.0f units", have, tracked->units);
        }
      } else if (tracked->type == sim::MissionType::Delivery || tracked->type == sim::MissionType::MultiDelivery || tracked->type == sim::MissionType::Smuggle) {
        const std::size_t ci = (std::size_t)tracked->commodity;
        if (ci < econ::kCommodityCount) {
          const double have = cargo[ci];
          ImGui::TextDisabled("Cargo: %.0f / %.0f units", have, tracked->units);
        }
      } else if (tracked->type == sim::MissionType::Passenger) {
        const int party = std::max(0, (int)std::llround(tracked->units));
        int used = 0;
        for (const auto& m : missions) {
          if (m.completed || m.failed) continue;
          if (m.type != sim::MissionType::Passenger) continue;
          used += std::max(0, (int)std::llround(m.units));
        }
        ImGui::TextDisabled("Passengers: %d (seats %d/%d)", party, used, std::max(0, passengerSeats));
      } else if (tracked->type == sim::MissionType::BountyScan) {
        ImGui::TextDisabled("Status: %s", tracked->scanned ? "scanned" : "not scanned");
      }

      const auto dest = uiMissionNextStop(*tracked);
      const sim::SystemId destSys = dest.first;
      const sim::StationId destSt = dest.second;

      if (destSys != 0) {
        if (ImGui::SmallButton("Select")) {
          galaxySelectedSystemId = destSys;
          showGalaxy = true;
        }
        ImGui::SameLine();
        if (ImGui::SmallButton("Plot Route")) {
          if (plotRouteToSystem(destSys)) {
            pendingArrivalTargetStationId = destSt;
            if (destSys == currentSystem->stub.id && destSt != 0) tryTargetStationById(destSt);
            showGalaxy = true;
          }
        }
        ImGui::SameLine();
        if (tracked->type == sim::MissionType::Salvage && destSys == currentSystem->stub.id && !tracked->scanned) {
          if (ImGui::SmallButton("Target site")) {
            if (!tryTargetSignalById(tracked->targetNpcId)) {
              toast(toasts, "Mission salvage site not detected in this system.", 2.2);
            }
          }
          ImGui::SameLine();
        }
        ImGui::SameLine();
        if (destSys == currentSystem->stub.id && destSt != 0) {
          if (ImGui::SmallButton("Target")) {
            tryTargetStationById(destSt);
          }
          ImGui::SameLine();
        }
        if (ImGui::SmallButton("Untrack")) {
          trackedMissionId = 0;
        }

        ImGui::Checkbox("Objective HUD overlay", &objectiveHudEnabled);
        ImGui::SameLine();
        ImGui::TextDisabled("(shows tracked mission in-flight)");

        ImGui::Checkbox("Auto-plot next leg (multi-delivery)", &missionTrackerAutoPlotNextLeg);
      } else {
        if (ImGui::SmallButton("Untrack")) trackedMissionId = 0;
      }
    } else {
      ImGui::TextDisabled("No active missions.");
    }
  }

  ImGui::Text("Credits: %.0f | Exploration data: %.0f cr", credits, explorationDataCr);
  ImGui::Text("Fuel: %.1f | Ship heat: %.0f", fuel, heat);

  // Power distributor
  {
    ImGui::Separator();
    ImGui::Text("Power distributor");

    auto capFrac = [](double v, double cap) {
      return (cap > 1e-9) ? std::clamp(v / cap, 0.0, 1.0) : 1.0;
    };

    const double engFrac = capFrac(distributorState.eng, distributorCfg.capEng);
    const double wepFrac = capFrac(distributorState.wep, distributorCfg.capWep);
    const double sysFrac = capFrac(distributorState.sys, distributorCfg.capSys);
    ImGui::Text("Capacitors: ENG %.0f%% | WEP %.0f%% | SYS %.0f%%", engFrac * 100.0, wepFrac * 100.0, sysFrac * 100.0);

    bool pipsChanged = false;
    pipsChanged |= ImGui::SliderInt("ENG pips", &distributorPips.eng, 0, sim::kPipMax);
    pipsChanged |= ImGui::SliderInt("WEP pips", &distributorPips.wep, 0, sim::kPipMax);
    pipsChanged |= ImGui::SliderInt("SYS pips", &distributorPips.sys, 0, sim::kPipMax);
    if (pipsChanged) {
      sim::normalizePips(distributorPips);
    }

    ImGui::TextDisabled("Pips total: %d (target %d)", distributorPips.eng + distributorPips.wep + distributorPips.sys, sim::kPipTotal);
    if (ImGui::SmallButton("Reset pips")) {
      distributorPips = sim::Pips{2, 2, 2};
    }
  }

  ImGui::Text("Cargo: %.0f / %.0f kg", cargoMassKg(cargo), cargoCapacityKg);
  ImGui::Text("Passengers: %d seats", std::max(0, passengerSeats));
	  ImGui::Text(
	      "Cargo scoop (%s): %s | Floating pods: %d", game::chordLabel(controls.actions.toggleCargoScoop).c_str(),
	      cargoScoopDeployed ? "DEPLOYED" : "RETRACTED", (int)floatingCargo.size());

	  {
	    double totalVouchers = 0.0;
	    for (const auto& kv : bountyVoucherByFaction) totalVouchers += std::max(0.0, kv.second);
	    if (totalVouchers > 0.0) {
	      ImGui::Text("Bounty vouchers: %.0f cr (redeem at matching faction stations)", totalVouchers);
	    }
	  }
  

  // Cargo management (jettison / dump)
  {
    ImGui::Separator();
    ImGui::Text("Cargo management");

    if (docked) {
      ImGui::TextDisabled("Docked: use the Market to trade cargo. Undock to jettison.");
    } else if (fsdState != FsdState::Idle || supercruiseState != SupercruiseState::Idle) {
      ImGui::TextDisabled("Cargo bay locked during FSD / supercruise.");
    } else {
      // Reserve mission cargo (prevents accidental dumping like the Market sell guard).
      std::array<double, econ::kCommodityCount> reservedUnits{};
      reservedUnits.fill(0.0);
      for (const auto& m : missions) {
        if (m.completed || m.failed) continue;
        if (m.type != sim::MissionType::Delivery
         && m.type != sim::MissionType::MultiDelivery
         && m.type != sim::MissionType::Smuggle) continue;
        const std::size_t idx = (std::size_t)m.commodity;
        if (idx >= reservedUnits.size()) continue;
        reservedUnits[idx] += std::max(0.0, m.units);
      }

      static bool allowDumpMissionCargo = false;
      static int jettisonChoice = 0;
      static float jettisonQty = 5.0f;

      std::vector<int> have;
      have.reserve((std::size_t)econ::kCommodityCount);
      for (int i = 0; i < (int)econ::kCommodityCount; ++i) {
        if (cargo[i] > 1e-4) have.push_back(i);
      }

      if (have.empty()) {
        ImGui::TextDisabled("Cargo hold empty.");
      } else {
        jettisonChoice = std::clamp(jettisonChoice, 0, (int)have.size() - 1);
        const int ci = have[(std::size_t)jettisonChoice];
        const auto cid = (econ::CommodityId)ci;

        const std::string preview =
          std::string(econ::commodityName(cid)) + " (" + std::to_string((int)std::round(cargo[ci])) + "u)";

        if (ImGui::BeginCombo("Commodity", preview.c_str())) {
          for (int n = 0; n < (int)have.size(); ++n) {
            const int idx = have[(std::size_t)n];
            const auto c = (econ::CommodityId)idx;
            std::string label =
              std::string(econ::commodityName(c)) + " (" + std::to_string((int)std::round(cargo[idx])) + "u)";
            if (ImGui::Selectable(label.c_str(), n == jettisonChoice)) {
              jettisonChoice = n;
            }
            if (n == jettisonChoice) ImGui::SetItemDefaultFocus();
          }
          ImGui::EndCombo();
        }

        ImGui::SetNextItemWidth(90.0f);
        ImGui::InputFloat("Units##jettisonQty", &jettisonQty, 1.0f, 10.0f, "%.1f");

        const double reserved = reservedUnits[(std::size_t)ci];
        const double freeUnits = std::max(0.0, cargo[ci] - reserved);

        if (reserved > 1e-4 && !allowDumpMissionCargo) {
          ImGui::TextDisabled("Onboard %.1fu | Reserved %.1fu | Free %.1fu", cargo[ci], reserved, freeUnits);
        } else {
          ImGui::TextDisabled("Onboard %.1fu", cargo[ci]);
        }

        ImGui::Checkbox("Allow dumping mission cargo", &allowDumpMissionCargo);
        ImGui::SameLine();
        ImGui::TextDisabled("(dangerous: can strand deliveries/smuggles)");

        // Detect whether an authority is likely to witness the dump (ports/scanners/nearby police).
        auto authorityWitness = [&]() -> std::pair<bool, std::string> {
          if (!currentSystem) return {false, ""};

          if (cargoScanActive && cargoScanFactionId != 0) {
            const std::string src = cargoScanSourceName.empty() ? std::string("Authorities") : cargoScanSourceName;
            return {true, src};
          }

          for (const auto& st : currentSystem->stations) {
            const math::Vec3d stPos = stationPosKm(st, timeDays);
            const double d = (stPos - ship.positionKm()).length();
            if (d <= st.commsRangeKm * 1.05) return {true, st.name};
          }

          for (const auto& c : contacts) {
            if (!c.alive || c.role != ContactRole::Police) continue;
            const double d = (c.ship.positionKm() - ship.positionKm()).length();
            if (d <= 120000.0) {
              const std::string src = c.name.empty() ? std::string("Police") : c.name;
              return {true, src};
            }
          }

          return {false, ""};
        };

        auto doJettison = [&](econ::CommodityId what, double units) {
          if (units <= 1e-6) return;
          const math::Vec3d fwd = ship.forward().normalized();
          const math::Vec3d dropPos = ship.positionKm() - fwd * 2600.0 + randUnit() * 250.0;
          spawnCargoBurst(what, units, dropPos, ship.velocityKmS(), 2 + rng.range(0, 2), true);

          // Counts toward pirate extortion demands (value by base price).
          if (pirateDemand.active && pirateDemand.requiredValueCr > 0.0) {
            pirateDemand.deliveredValueCr += units * econ::commodityDef(what).basePrice;
          }
        };

        auto applyDumpCrime = [&](double dumpedValueCr, const std::string& witnessName, bool contraband) {
          if (!currentSystem) return;
          const core::u32 j = currentSystem->stub.factionId;
          if (j == 0) return;

          const double bountyAdd = (contraband ? 250.0 : 110.0) + dumpedValueCr * (contraband ? 0.22 : 0.12);
          const double repPenalty = -std::clamp(0.8 + dumpedValueCr / (contraband ? 4200.0 : 6500.0), 0.8, contraband ? 8.0 : 5.0);

          commitCrime(j, bountyAdd, repPenalty, (contraband ? "contraband dumping" : "cargo dumping") + std::string(" reported by ") + witnessName, false);
          toast(toasts,
                witnessName + ": illegal dumping detected. Bounty +" + std::to_string((int)std::round(bountyAdd)) + " cr",
                3.0);
        };

        // Jettison selected commodity
        {
          const double maxUnits = allowDumpMissionCargo ? cargo[ci] : freeUnits;
          const double req = std::max(0.0, (double)jettisonQty);
          const double take = std::min(maxUnits, req);
          const bool can = (take > 1e-6);

          ImGui::BeginDisabled(!can);
          if (ImGui::Button("Jettison")) {
            cargo[ci] = std::max(0.0, cargo[ci] - take);
            doJettison(cid, take);

            const bool contraband = (currentSystem && currentSystem->stub.factionId != 0)
                                      ? isIllegalCommodity(currentSystem->stub.factionId, cid)
                                      : false;

            const auto cNameSv = econ::commodityName(cid);
            std::string msg = "Jettisoned: ";
            msg.append(cNameSv.data(), cNameSv.size());
            msg += " x";
            msg += std::to_string((int)std::round(take));
            toast(toasts, msg, contraband ? 2.6 : 2.0);

            const auto wit = authorityWitness();
            if (wit.first) {
              const double v = take * econ::commodityDef(cid).basePrice;
              applyDumpCrime(v, wit.second, contraband);
            }
          }
          ImGui::EndDisabled();

          ImGui::SameLine();
          ImGui::TextDisabled("(drops cargo pods behind your ship)");
        }

        // Dump all contraband (quick smuggler panic button)
        if (currentSystem && currentSystem->stub.factionId != 0) {
          const core::u32 j = currentSystem->stub.factionId;

          bool anyIllegal = false;
          for (int i = 0; i < (int)econ::kCommodityCount; ++i) {
            if (cargo[i] <= 1e-6) continue;
            if (!isIllegalCommodity(j, (econ::CommodityId)i)) continue;
            const double free = allowDumpMissionCargo ? cargo[i] : std::max(0.0, cargo[i] - reservedUnits[(std::size_t)i]);
            if (free > 1e-6) { anyIllegal = true; break; }
          }

          ImGui::BeginDisabled(!anyIllegal);
          if (ImGui::Button("Dump contraband")) {
            double dumpedValue = 0.0;
            int dumpedTypes = 0;

            for (int i = 0; i < (int)econ::kCommodityCount; ++i) {
              const auto c = (econ::CommodityId)i;
              if (cargo[i] <= 1e-6) continue;
              if (!isIllegalCommodity(j, c)) continue;

              const double free = allowDumpMissionCargo ? cargo[i] : std::max(0.0, cargo[i] - reservedUnits[(std::size_t)i]);
              if (free <= 1e-6) continue;

              cargo[i] = std::max(0.0, cargo[i] - free);
              doJettison(c, free);
              dumpedValue += free * econ::commodityDef(c).basePrice;
              dumpedTypes++;
            }

            toast(toasts, "Dumped contraband pods (" + std::to_string(dumpedTypes) + " types).", 2.6);

            const auto wit = authorityWitness();
            if (wit.first && dumpedValue > 1e-6) {
              applyDumpCrime(dumpedValue, wit.second, true);
            }
          }
          ImGui::EndDisabled();
          ImGui::SameLine();
          ImGui::TextDisabled("(jettison all illegal goods in this jurisdiction)");
        }
      }
    }
  }

ImGui::Text("Shield: %.0f/%.0f | Hull: %.0f/%.0f", playerShield, playerShieldMax, playerHull, playerHullMax);

  {
    const auto& hd = kHullDefs[(int)shipHullClass];
    ImGui::Text("Hull: %s | Thrusters Mk%d | Shield Mk%d | Distributor Mk%d", hd.name, thrusterMk, shieldMk, distributorMk);
    if (smuggleHoldMk > 0) {
      ImGui::TextDisabled("Smuggling compartments: Mk%d (reduced scan chance)", smuggleHoldMk);
    } else {
      ImGui::TextDisabled("Smuggling compartments: none");
    }
    ImGui::Text("Weapons: LMB %s | RMB %s", weaponDef(weaponPrimary).name, weaponDef(weaponSecondary).name);
    ImGui::TextDisabled("Accel %.3f km/s^2 | Turn %.2f rad/s^2", playerBaseLinAccelKmS2, playerBaseAngAccelRadS2);
  }

  ImGui::Text("Ship pos: (%.0f,%.0f,%.0f) km", ship.positionKm().x, ship.positionKm().y, ship.positionKm().z);
  ImGui::Text("Vel: %.3f km/s", ship.velocityKmS().length());

  ImGui::Separator();
  ImGui::Text("Controls");
  ImGui::Checkbox("Mouse steer##mouseSteer", &mouseSteer);
  ImGui::SameLine();
  ImGui::TextDisabled("(%s)", game::chordLabel(controls.actions.toggleMouseSteer).c_str());
  ImGui::Checkbox("Invert mouse Y", &mouseInvertY);
  ImGui::SliderFloat("Mouse sensitivity", &mouseSensitivity, 0.0006f, 0.0080f, "%.4f");
  ImGui::TextDisabled("Mouse steer captures the cursor (relative mouse mode).");

  if (scanning) {
    ImGui::Separator();
    ImGui::TextColored(ImVec4(0.9f, 0.95f, 1.0f, 1.0f), "%s", scanLabel.empty() ? "Scanning..." : scanLabel.c_str());
    const float frac = (float)std::clamp(scanProgressSec / std::max(0.01, scanDurationSec), 0.0, 1.0);
    ImGui::ProgressBar(frac, ImVec2(-1, 0));
    ImGui::TextDisabled("%s cancels scan.", game::chordLabel(controls.actions.scannerAction).c_str());
  }

  if (cargoScanActive) {
    ImGui::Separator();
    ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.3f, 1.0f), "Cargo scan: %s",
                       cargoScanSourceName.empty() ? "Security" : cargoScanSourceName.c_str());
    const float frac = (float)std::clamp(cargoScanProgressSec / std::max(0.01, cargoScanDurationSec), 0.0, 1.0);
    ImGui::ProgressBar(frac, ImVec2(-1, 0));
    ImGui::TextDisabled("Leave scan range / dock to break scan.");
  }

  if (bribeOffer.active && timeDays < bribeOffer.untilDays) {
    ImGui::Separator();
    const std::string src = bribeOffer.sourceName.empty() ? std::string("Security") : bribeOffer.sourceName;
    ImGui::TextColored(ImVec4(1.0f, 0.65f, 0.35f, 1.0f), "Authority comms: %s", src.c_str());
    ImGui::Text("Bribe: %.0f cr (C) | Comply: (I) | Run: leave scan range", bribeOffer.amountCr);
    ImGui::TextDisabled("Fine: %.0f cr | Bounty if you run: %.0f cr", bribeOffer.fineCr, bribeOffer.fineCr);
    if (!bribeOffer.detail.empty()) ImGui::TextDisabled("Detected: %s", bribeOffer.detail.c_str());

    const double totalSec = std::max(0.001, (bribeOffer.untilDays - bribeOffer.startDays) * 86400.0);
    const double remainSec = std::max(0.0, (bribeOffer.untilDays - timeDays) * 86400.0);
    const float frac = (float)std::clamp(1.0 - (remainSec / totalSec), 0.0, 1.0);
    ImGui::ProgressBar(frac, ImVec2(-1, 0));
    ImGui::TextDisabled("Time remaining: %.1f s", remainSec);
  }

  if (currentSystem && currentSystem->stub.factionId != 0) {
    const core::u32 j = currentSystem->stub.factionId;
    ImGui::Separator();
    ImGui::Text("Jurisdiction: %s", factionName(j).c_str());
    ImGui::Text("Illegal goods: %s", illegalListString(j).c_str());
    if (hasIllegalCargo(j)) {
      ImGui::TextColored(ImVec4(1.0f, 0.55f, 0.55f, 1.0f), "Warning: CONTRABAND on board");
    }
  }

  // FSD status (Charging / Jumping / Cooling)
  {
    const bool cooling = (fsdState == FsdState::Idle && timeDays < fsdReadyDay);
    if (fsdState != FsdState::Idle || cooling) {
      ImGui::Separator();

      const char* label = "Idle";
      double timer = 0.0;
      double total = 1.0;

      if (fsdState == FsdState::Charging) {
        label = "Charging";
        total = kFsdChargeSec;
        timer = std::clamp(total - fsdChargeRemainingSec, 0.0, total);
      } else if (fsdState == FsdState::Jumping) {
        label = "Jumping";
        total = std::max(0.001, fsdTravelTotalSec);
        timer = std::clamp(total - fsdTravelRemainingSec, 0.0, total);
      } else {
        // Cooldown uses simulated time.
        label = "Cooling";
        total = kFsdCooldownSec;
        const double remain = std::max(0.0, (fsdReadyDay - timeDays) * 86400.0);
        timer = std::clamp(total - remain, 0.0, total);
      }

      ImGui::Text("FSD: %s", label);
      ImGui::ProgressBar((float)std::clamp(timer / total, 0.0, 1.0), ImVec2(-1, 0));
      if (fsdState == FsdState::Charging) ImGui::TextDisabled("J cancels during charge.");
    }
  }

  if (supercruiseState != SupercruiseState::Idle) {
    ImGui::Separator();

    if (supercruiseState == SupercruiseState::Charging) {
      ImGui::TextColored(ImVec4(0.6f, 0.9f, 1.0f, 1.0f), "Supercruise: CHARGING");
      const float t = (float)std::clamp(1.0 - (supercruiseChargeRemainingSec / std::max(0.001, kSupercruiseChargeSec)), 0.0, 1.0);
      ImGui::ProgressBar(t, ImVec2(-1, 0));
    } else if (supercruiseState == SupercruiseState::Active) {
      ImGui::TextColored(ImVec4(0.6f, 0.9f, 1.0f, 1.0f), "Supercruise: ACTIVE");
      if (supercruiseDistKm > 0.0) ImGui::Text("Range: %.0f km", supercruiseDistKm);
      if (supercruiseTtaSec > 0.0 && supercruiseTtaSec < 1e8) ImGui::Text("ETA: %.1f s", supercruiseTtaSec);

      if (supercruiseSafeDropReady) {
        ImGui::TextColored(ImVec4(0.2f, 1.0f, 0.2f, 1.0f), "SAFE DROP");
      } else {
        ImGui::TextColored(ImVec4(1.0f, 0.9f, 0.2f, 1.0f), "Keep ~7s ETA for safe drop");
      }

      if (!supercruiseAssist) {
        ImGui::Text("Manual drop: press H");
      } else {
        ImGui::Text("Nav assist: auto-drop when safe");
      }
    } else if (supercruiseState == SupercruiseState::Cooldown) {
      ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "Supercruise: COOLDOWN");
      const float t = (float)std::clamp(1.0 - (supercruiseCooldownRemainingSec / std::max(0.001, supercruiseCooldownTotalSec)), 0.0, 1.0);
      ImGui::ProgressBar(t, ImVec2(-1, 0));
    }
  }

  // Navigation / route status
  if (navRoute.size() >= 2) {
    ImGui::Separator();
    const int totalJumps = (int)navRoute.size() - 1;
    const int remaining = std::max(0, totalJumps - (int)navRouteHop);

    double remDist = 0.0;
    double remFuel = 0.0;
    for (std::size_t i = navRouteHop; i + 1 < navRoute.size(); ++i) {
      const auto& a = universe.getSystem(navRoute[i]).stub;
      const auto& b = universe.getSystem(navRoute[i + 1]).stub;
      const double d = (b.posLy - a.posLy).length();
      remDist += d;
      remFuel += fsdFuelCostFor(d);
    }

    ImGui::Text("Route: %d jumps (%d remaining) | est fuel %.1f", totalJumps, remaining, remFuel);
    ImGui::TextDisabled("Auto-run: %s", navAutoRun ? "ON" : "OFF");

    if (remFuel > fuel + 1e-6) {
      ImGui::TextColored(ImVec4(1.0f, 0.55f, 0.45f, 1.0f), "Fuel shortfall: +%.1f", remFuel - fuel);
    }
  }

  // Target summary
  ImGui::Separator();
  if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
    const auto& st = currentSystem->stations[target.index];
    ImGui::Text("Target: %s (%.0f km)", st.name.c_str(), (stationPosKm(st, timeDays) - ship.positionKm()).length());
  } else if (target.kind == TargetKind::Planet && target.index < currentSystem->planets.size()) {
    const auto& p = currentSystem->planets[target.index];
    const math::Vec3d pPos = sim::orbitPosition3DAU(p.orbit, timeDays) * kAU_KM;
    ImGui::Text("Target: %s (%.0f km)", p.name.c_str(), (pPos - ship.positionKm()).length());
  } else if (target.kind == TargetKind::Contact && target.index < contacts.size()) {
    const auto& c = contacts[target.index];
    ImGui::Text("Target: %s [%s] (%.0f km)", c.name.c_str(), contactRoleName(c.role), (c.ship.positionKm() - ship.positionKm()).length());
    if (c.distressVictim && c.distressNeedUnits > 1e-6) {
      const auto def = econ::commodityDef(c.distressNeedCommodity);
      const int needUnits = std::max(1, (int)std::llround(c.distressNeedUnits));
      ImGui::TextDisabled("Distress request: %d %s", needUnits, def.name);
      ImGui::TextDisabled("Reward: ~%d cr | Rep: +%.1f",
                          (int)std::llround(c.distressRewardCr),
                          c.distressRepReward);
      const std::string kScan = game::chordLabel(controls.actions.scannerAction);
      ImGui::TextDisabled("Scan (%s) to transfer supplies.", kScan.c_str());
    }
  } else if (target.kind == TargetKind::Signal && target.index < signals.size()) {
    const auto& s = signals[target.index];
    const double distKm = (s.posKm - ship.positionKm()).length();
    const char* base = signalTypeName(s.type);
    if (s.type == SignalType::Resource && s.hasResourcePlan) {
      base = sim::resourceFieldKindName(s.resource.kind);
    }
    ImGui::Text("Target: %s Signal (%.0f km)", base, distKm);
    if (s.type == SignalType::Resource && s.hasResourcePlan) {
      const auto& p = s.resource;
      const auto& a = econ::commodityDef(p.primaryYield);
      const auto& b = econ::commodityDef(p.secondaryYield);
      const int pa = (int)std::lround(std::clamp(p.primaryChance, 0.0, 1.0) * 100.0);
      const int pb = 100 - pa;
      const int rich = (int)std::lround(std::clamp(p.richness01, 0.0, 1.0) * 100.0);
      ImGui::TextDisabled("Richness: %d%% | %s ~%d%% / %s ~%d%%", rich, a.name, pa, b.name, pb);
    }
    if (s.type == SignalType::Distress && s.hasDistressPlan && !s.distressCompleted) {
      if (s.distress.hasVictim) {
        const auto def = econ::commodityDef(s.distress.needCommodity);
        const int needUnits = std::max(1, (int)std::llround(s.distress.needUnits));
        ImGui::TextDisabled("Request: %d %s | reward ~%d cr",
                            needUnits,
                            def.name,
                            (int)std::llround(s.distress.rewardCr));
      }
      if (s.distress.ambush) {
        ImGui::TextDisabled("Warning: hostile activity suspected.");
      }
    }
  } else if (target.kind == TargetKind::Star) {
    ImGui::Text("Target: Star (%s)", starClassName(currentSystem->star.cls));
  } else {
    const std::string kT = game::chordLabel(controls.actions.targetStationCycle);
    const std::string kB = game::chordLabel(controls.actions.targetPlanetCycle);
    const std::string kN = game::chordLabel(controls.actions.targetContactCycle);
    const std::string kU = game::chordLabel(controls.actions.targetStar);
    const std::string hint = "[" + kT + " station / " + kB + " planet / " + kN + " contact / " + kU + " star]";
    ImGui::Text("Target: (none)  %s", hint.c_str());
  }

  ImGui::Separator();
  ImGui::TextDisabled("Controls:");
  {
    auto axisLabel = [&](const game::AxisPair& a) {
      return game::scancodeLabel(a.positive) + "/" + game::scancodeLabel(a.negative);
    };

    const std::string trFwd = axisLabel(controls.axes.thrustForward);
    const std::string trStr = axisLabel(controls.axes.thrustRight);
    const std::string trVert = axisLabel(controls.axes.thrustUp);
    const std::string rotPitch = axisLabel(controls.axes.pitch);
    const std::string rotYaw = axisLabel(controls.axes.yaw);
    const std::string rotRoll = axisLabel(controls.axes.roll);

    const std::string boostKey = game::scancodeLabel(controls.holds.boost);
    const std::string brakeKey = game::scancodeLabel(controls.holds.brake);

    const std::string kSuper = game::chordLabel(controls.actions.supercruise);
    const std::string kJump = game::chordLabel(controls.actions.fsdJump);
    const std::string kAuto = game::chordLabel(controls.actions.toggleAutopilot);
    const std::string kT = game::chordLabel(controls.actions.targetStationCycle);
    const std::string kB = game::chordLabel(controls.actions.targetPlanetCycle);
    const std::string kN = game::chordLabel(controls.actions.targetContactCycle);
    const std::string kU = game::chordLabel(controls.actions.targetStar);
    const std::string kClear = game::chordLabel(controls.actions.clearTarget);
    const std::string kScan = game::chordLabel(controls.actions.scannerAction);
    const std::string kClearance = game::chordLabel(controls.actions.requestDockingClearance);
    const std::string kDock = game::chordLabel(controls.actions.dockOrUndock);

    const std::string kGalaxy = game::chordLabel(controls.actions.toggleGalaxy);
    const std::string kShip = game::chordLabel(controls.actions.toggleShip);
    const std::string kMarket = game::chordLabel(controls.actions.toggleMarket);
    const std::string kContacts = game::chordLabel(controls.actions.toggleContacts);
    const std::string kMissions = game::chordLabel(controls.actions.toggleMissions);
    const std::string kScannerW = game::chordLabel(controls.actions.toggleScanner);
    const std::string kTrade = game::chordLabel(controls.actions.toggleTrade);
    const std::string kGuide = game::chordLabel(controls.actions.toggleGuide);

    const std::string kSave = game::chordLabel(controls.actions.quicksave);
    const std::string kLoad = game::chordLabel(controls.actions.quickload);

    const std::string kSpriteLab = game::chordLabel(controls.actions.toggleSpriteLab);
    const std::string kVfxLab = game::chordLabel(controls.actions.toggleVfxLab);
    const std::string kPostFx = game::chordLabel(controls.actions.togglePostFx);

    ImGui::BulletText(
        "Translate: fwd/back %s | strafe %s | vertical %s", trFwd.c_str(), trStr.c_str(), trVert.c_str());
    ImGui::BulletText(
        "Rotate: pitch %s | yaw %s | roll %s", rotPitch.c_str(), rotYaw.c_str(), rotRoll.c_str());
    ImGui::BulletText(
        "%s boost | %s brake | LMB: %s | RMB: %s",
        boostKey.c_str(),
        brakeKey.c_str(),
        weaponDef(weaponPrimary).name,
        weaponDef(weaponSecondary).name);
    ImGui::BulletText(
        "%s supercruise (charge/engage). While active: %s drops (~7s safe). During interdiction: %s submits",
        kSuper.c_str(),
        kSuper.c_str(),
        kSuper.c_str());
    ImGui::BulletText("%s engage FSD jump (uses plotted route if present)", kJump.c_str());
    ImGui::BulletText("%s docking computer (auto-clearance + auto-dock)", kAuto.c_str());
    ImGui::BulletText(
        "%s/%s/%s/%s cycle targets, %s clear target",
        kT.c_str(),
        kB.c_str(),
        kN.c_str(),
        kU.c_str(),
        kClear.c_str());
    ImGui::BulletText("%s scan target (missions + exploration scans)", kScan.c_str());
    ImGui::BulletText("%s request docking clearance | %s dock / undock", kClearance.c_str(), kDock.c_str());
    ImGui::BulletText(
        "%s Galaxy map, %s Ship, %s Market, %s Contacts, %s Missions, %s Scanner, %s Trade, %s Guide",
        kGalaxy.c_str(),
        kShip.c_str(),
        kMarket.c_str(),
        kContacts.c_str(),
        kMissions.c_str(),
        kScannerW.c_str(),
        kTrade.c_str(),
        kGuide.c_str());
    ImGui::BulletText("%s quicksave, %s quickload", kSave.c_str(), kLoad.c_str());
    ImGui::BulletText("%s Sprite Lab, %s VFX Lab, %s Post FX", kSpriteLab.c_str(), kVfxLab.c_str(), kPostFx.c_str());
  }

  ImGui::Separator();
  float scMax = (float)supercruiseMaxSpeedKmS;
  if (ImGui::SliderFloat("Supercruise max speed (km/s)", &scMax, 4000.0f, 60000.0f, "%.0f")) {
    supercruiseMaxSpeedKmS = (double)scMax;
  }
  ImGui::Checkbox("Supercruise assist", &supercruiseAssist);

  ImGui::Separator();
  ImGui::TextDisabled("Flight Assist:");
  ImGui::Checkbox("Local reference frame", &localFrameEnabled);
  float lfTau = (float)localFrameBlendTauSec;
  if (ImGui::SliderFloat("Local frame blend (s)", &lfTau, 0.05f, 6.0f, "%.2f")) {
    localFrameBlendTauSec = (double)lfTau;
  }
  ImGui::Checkbox("Blue-zone turn assist", &blueZoneTurnAssist);

  ImGui::End();
}

if (showGuide) {
  ImGui::Begin("Pilot Guide");

  const std::string kStation = game::chordLabel(controls.actions.targetStationCycle);
  const std::string kClearance = game::chordLabel(controls.actions.requestDockingClearance);
  const std::string kAuto = game::chordLabel(controls.actions.toggleAutopilot);
  const std::string kDock = game::chordLabel(controls.actions.dockOrUndock);
  const std::string kMissions = game::chordLabel(controls.actions.toggleMissions);
  const std::string kJump = game::chordLabel(controls.actions.fsdJump);
  const std::string kSuper = game::chordLabel(controls.actions.supercruise);
  const std::string kTrade = game::chordLabel(controls.actions.toggleTrade);
  const std::string kScan = game::chordLabel(controls.actions.scannerAction);
  const std::string kSprite = game::chordLabel(controls.actions.toggleSpriteLab);
  const std::string kVfx = game::chordLabel(controls.actions.toggleVfxLab);
  const std::string kGalaxy = game::chordLabel(controls.actions.toggleGalaxy);
  const std::string kMarket = game::chordLabel(controls.actions.toggleMarket);

  ImGui::Text("Quick start checklist");
  ImGui::BulletText(
      "Dock: target station (%s), request clearance (%s) (or let docking computer do it), engage docking computer (%s), dock/undock (%s)",
      kStation.c_str(),
      kClearance.c_str(),
      kAuto.c_str(),
      kDock.c_str());
  ImGui::BulletText(
      "Missions: open Mission Board (%s). Accept/Plot auto-plot; track a mission for the Objective HUD",
      kMissions.c_str());
  ImGui::BulletText("Jump: %s (uses plotted route if present)", kJump.c_str());
  ImGui::BulletText(
      "Supercruise: %s engage, %s drop near ~7s ETA. Interdiction: align to escape vector; %s submits",
      kSuper.c_str(),
      kSuper.c_str(),
      kSuper.c_str());
  ImGui::BulletText("Trade helper: %s (best local routes)", kTrade.c_str());
  ImGui::BulletText("Scan: %s (stars/planets/contacts). Sell exploration data at stations", kScan.c_str());
  ImGui::BulletText("Sprite Lab: %s (procedural UI icons preview + cache)", kSprite.c_str());
  ImGui::BulletText("VFX Lab: %s (starfield + particles)", kVfx.c_str());

  ImGui::Separator();
  ImGui::Text("Quick actions");
  const std::string openGalaxyLabel = "Open Galaxy (" + kGalaxy + ")";
  const std::string openMarketLabel = "Open Market (" + kMarket + ")";
  const std::string openMissionsLabel = "Open Missions (" + kMissions + ")";
  const std::string openTradeLabel = "Open Trade (" + kTrade + ")";

  if (ImGui::Button(openGalaxyLabel.c_str())) showGalaxy = true;
  ImGui::SameLine();
  if (ImGui::Button(openMarketLabel.c_str())) showEconomy = true;
  ImGui::SameLine();
  if (ImGui::Button(openMissionsLabel.c_str())) showMissions = true;
  ImGui::SameLine();
  if (ImGui::Button(openTradeLabel.c_str())) showTrade = true;

  ImGui::TextDisabled("Tip: %s opens the Command Palette (search actions, systems, stations, missions).",
                      game::chordLabel(controls.actions.commandPalette).c_str());
  ImGui::TextDisabled("Tip: Windows > Notifications keeps a history of HUD toasts.");

  ImGui::Separator();
  ImGui::Text("Status");
  ImGui::Text("Docked: %s", docked ? "YES" : "NO");
  if (docked && dockedStationId != 0) {
    // Find the docked station type for a small hint.
    const sim::Station* dockSt = nullptr;
    for (const auto& st : currentSystem->stations) {
      if (st.id == dockedStationId) {
        dockSt = &st;
        break;
      }
    }
    if (dockSt) {
      if (dockSt->type == econ::StationType::Shipyard) {
        ImGui::TextDisabled("Shipyard available here: upgrade hull/modules/weapons under Market > Shipyard");
      } else {
        ImGui::TextDisabled("No shipyard at this station. Try another station for upgrades.");
      }
    }
  }

  if (supercruiseState == SupercruiseState::Active && sim::interdictionInProgress(interdiction)) {
    ImGui::Separator();
    ImGui::TextColored(ImVec4(1.0f, 0.55f, 0.45f, 1.0f), "Interdiction in progress!");
    ImGui::TextDisabled("Point your ship toward the escape vector and keep it aligned.");
    ImGui::TextDisabled("Press H to submit (safer, but drops you)."
    );
  }

  ImGui::End();
}


if (showSprites) {
  ImGui::Begin("Sprite Lab");

  ImGui::TextDisabled(
      "%s toggles this panel. Procedural icons are generated deterministically from (kind, seed).",
      game::chordLabel(controls.actions.toggleSpriteLab).c_str());
  ImGui::Text("Sprite cache: %zu / %zu", spriteCache.size(), spriteCache.maxEntries());
  if (ImGui::Button("Clear sprite cache")) spriteCache.clear();

  ImGui::Text("HUD icon atlas: %zu / %zu (grid cap %d)", hudAtlas.size(), hudAtlas.maxEntries(), hudAtlas.capacity());
  ImGui::SameLine();
  if (ImGui::Button("Clear HUD atlas")) hudAtlas.clear();
  ImGui::TextDisabled("Atlas %dx%d | Cell %d | Pad %d", hudAtlas.atlasSizePx(), hudAtlas.atlasSizePx(), hudAtlas.cellSizePx(), hudAtlas.paddingPx());
  ImGui::Image((ImTextureID)(intptr_t)hudAtlas.texture().handle(), ImVec2(256, 256), ImVec2(0,0), ImVec2(1,1));

  ImGui::Separator();

  static int kindIdx = 0;
  static int texSize = 64;
  static core::u64 seed = 0xC0FFEEULL;

  const char* kindNames[] = {"Commodity", "Faction", "Mission", "Ship", "Station", "Planet", "Star",
                             "Cargo", "Asteroid", "Signal",
                             "HudReticle", "HudLead", "HudVelocity"};
  ImGui::Combo("Kind", &kindIdx, kindNames, IM_ARRAYSIZE(kindNames));
  ImGui::SliderInt("Texture size", &texSize, 16, 128);
  ImGui::InputScalar("Seed (u64)", ImGuiDataType_U64, &seed);
  ImGui::SameLine();
  if (ImGui::Button("Randomize")) {
    seed = core::hashCombine((core::u64)SDL_GetPerformanceCounter(), core::fnv1a64("sprite_rand"));
  }

  const auto kind = (render::SpriteKind)std::clamp(kindIdx, 0, (int)IM_ARRAYSIZE(kindNames) - 1);
  const auto& tex = spriteCache.get(kind, seed, texSize);

  const float preview = std::max(128.0f, (float)texSize * 3.0f);
  ImGui::Image((ImTextureID)(intptr_t)tex.handle(), ImVec2(preview, preview));

  ImGui::Separator();

  if (ImGui::CollapsingHeader("Commodity icon atlas", ImGuiTreeNodeFlags_DefaultOpen)) {
    if (ImGui::BeginTable("commod_atlas", 6, ImGuiTableFlags_SizingFixedFit)) {
      for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
        const auto cid = (econ::CommodityId)i;
        const core::u64 cSeed = core::hashCombine(core::fnv1a64("commodity"), (core::u64)i);
        const auto& cTex = spriteCache.get(render::SpriteKind::Commodity, cSeed, 32);

        ImGui::TableNextColumn();
        ImGui::PushID((int)i);

        const auto nameSv = econ::commodityName(cid);
        const auto codeSv = econ::commodityCode(cid);

        ImGui::Image((ImTextureID)(intptr_t)cTex.handle(), ImVec2(26, 26));
        if (ImGui::IsItemHovered()) {
          const auto& big = spriteCache.get(render::SpriteKind::Commodity, cSeed, 96);
          ImGui::BeginTooltip();
          ImGui::Text("%.*s (%.*s)", (int)nameSv.size(), nameSv.data(), (int)codeSv.size(), codeSv.data());
          ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(96, 96));
          ImGui::EndTooltip();
        }

        ImGui::TextUnformatted(codeSv.data(), codeSv.data() + codeSv.size());

        ImGui::PopID();
      }
      ImGui::EndTable();
    }
  }

  if (ImGui::CollapsingHeader("System icon samples", ImGuiTreeNodeFlags_DefaultOpen)) {
    if (currentSystem) {
      const core::u64 sSeed = core::hashCombine(core::fnv1a64("star"), (core::u64)currentSystem->stub.seed);
      const auto& sTex = spriteCache.get(render::SpriteKind::Star, sSeed, 64);
      ImGui::Image((ImTextureID)(intptr_t)sTex.handle(), ImVec2(64, 64));
      ImGui::SameLine();
      ImGui::Text("Star: %s", starClassName(currentSystem->star.cls));

      ImGui::Separator();

      for (std::size_t i = 0; i < std::min<std::size_t>(currentSystem->planets.size(), 8); ++i) {
        const auto& p = currentSystem->planets[i];
        const core::u64 pSeed = core::hashCombine(core::hashCombine(core::fnv1a64("planet"), (core::u64)currentSystem->stub.seed), (core::u64)i);
        const auto& pTex = spriteCache.get(render::SpriteKind::Planet, pSeed, 64);
        ImGui::Image((ImTextureID)(intptr_t)pTex.handle(), ImVec2(48, 48));
        ImGui::SameLine();
        ImGui::Text("%s (%s)", p.name.c_str(), planetTypeName(p.type));
      }
    } else {
      ImGui::TextDisabled("No current system loaded.");
    }
  }

  ImGui::End();
}


if (showWorldVisuals) {
  ImGui::Begin("World Visuals", &showWorldVisuals);

  ImGui::TextDisabled("Procedural surface textures for stars/planets + star-relative lighting.");
  ImGui::Text("Surface texture cache: %zu / %zu", surfaceTexCache.size(), surfaceTexCache.maxEntries());
  ImGui::SameLine();
  if (ImGui::Button("Clear surface cache")) surfaceTexCache.clear();

  ImGui::Separator();

  ImGui::Checkbox("Procedural star/planet surfaces", &worldUseProceduralSurfaces);
  ImGui::SameLine();
  ImGui::Checkbox("Surface previews in UI", &worldUseSurfaceInUi);

  if (worldUseProceduralSurfaces) {
    ImGui::SetNextItemWidth(260.0f);
    ImGui::SliderInt("Surface tex width", &worldSurfaceTexWidth, 128, 1024, "%d");
    ImGui::SameLine();
    ImGui::TextDisabled("(height = width/2)");

    ImGui::SetNextItemWidth(260.0f);
    ImGui::SliderFloat("Star intensity", &worldStarIntensity, 0.5f, 8.0f, "%.2f");
    ImGui::SameLine();
    ImGui::Checkbox("Unlit star", &worldStarUnlit);
  } else {
    ImGui::TextDisabled("Meshes fall back to the checker texture.");
  }

  ImGui::Separator();
  ImGui::TextDisabled("Secondary layers (visual only)");

  ImGui::Checkbox("Atmospheres (limb glow)", &worldAtmospheresEnabled);
  if (worldAtmospheresEnabled) {
    ImGui::SetNextItemWidth(260.0f);
    ImGui::SliderFloat("Atmosphere intensity", &worldAtmoIntensity, 0.0f, 3.0f, "%.2f");
    ImGui::SetNextItemWidth(260.0f);
    ImGui::SliderFloat("Atmosphere rim power", &worldAtmoPower, 1.5f, 10.0f, "%.2f");
    ImGui::SetNextItemWidth(260.0f);
    ImGui::SliderFloat("Atmosphere shell scale", &worldAtmoShellScale, 1.005f, 1.10f, "%.3f");
    ImGui::SetNextItemWidth(260.0f);
    ImGui::SliderFloat("Day-side boost", &worldAtmoSunLitBoost, 0.0f, 1.0f, "%.2f");
    ImGui::SetNextItemWidth(260.0f);
    ImGui::SliderFloat("Forward scatter", &worldAtmoForwardScatter, 0.0f, 1.0f, "%.2f");
    ImGui::Checkbox("Tint with star color", &worldAtmoTintWithStar);
  }

  ImGui::BeginDisabled(!worldUseProceduralSurfaces);
  ImGui::Checkbox("Cloud shell (alpha)", &worldCloudsEnabled);
  if (worldCloudsEnabled) {
    ImGui::SetNextItemWidth(260.0f);
    ImGui::SliderFloat("Cloud opacity", &worldCloudOpacity, 0.0f, 1.0f, "%.2f");
    ImGui::SetNextItemWidth(260.0f);
    ImGui::SliderFloat("Cloud shell scale", &worldCloudShellScale, 1.002f, 1.06f, "%.3f");
    ImGui::SetNextItemWidth(260.0f);
    ImGui::SliderFloat("Cloud spin (deg/s)", &worldCloudSpinDegPerSec, 0.0f, 18.0f, "%.1f");
  }
  ImGui::EndDisabled();

  if (!worldUseProceduralSurfaces) {
    ImGui::TextDisabled("Cloud shells require procedural surfaces enabled.");
  }

  if (ImGui::CollapsingHeader("Surface generator preview", ImGuiTreeNodeFlags_DefaultOpen)) {
    const float thumbW = 192.0f;
    const float thumbH = thumbW * 0.5f;

    auto preview = [&](render::SurfaceKind kind, const char* name, core::u64 s) {
      const auto& t = surfaceTexCache.get(kind, s, worldSurfaceTexWidth);
      ImGui::TextUnformatted(name);
      ImGui::Image((ImTextureID)(intptr_t)t.handle(), ImVec2(thumbW, thumbH), ImVec2(0, 0), ImVec2(1, 1));
    };

    const core::u64 base = currentSystem ? (core::u64)currentSystem->stub.seed : seed;
    preview(render::SurfaceKind::Star, "Star", core::hashCombine(base, core::fnv1a64("surf_star")));
    ImGui::Separator();

    if (ImGui::BeginTable("surfprev", 2, ImGuiTableFlags_SizingFixedFit)) {
      ImGui::TableNextRow();
      ImGui::TableNextColumn();
      preview(render::SurfaceKind::Rocky, "Rocky", core::hashCombine(base, core::fnv1a64("surf_rocky")));
      ImGui::TableNextColumn();
      preview(render::SurfaceKind::Desert, "Desert", core::hashCombine(base, core::fnv1a64("surf_desert")));
      ImGui::TableNextRow();
      ImGui::TableNextColumn();
      preview(render::SurfaceKind::Ocean, "Ocean", core::hashCombine(base, core::fnv1a64("surf_ocean")));
      ImGui::TableNextColumn();
      preview(render::SurfaceKind::Ice, "Ice", core::hashCombine(base, core::fnv1a64("surf_ice")));
      ImGui::TableNextRow();
      ImGui::TableNextColumn();
      preview(render::SurfaceKind::GasGiant, "Gas Giant", core::hashCombine(base, core::fnv1a64("surf_gas")));
      ImGui::TableNextColumn();
      preview(render::SurfaceKind::Clouds, "Clouds", core::hashCombine(base, core::fnv1a64("surf_clouds")));
      ImGui::EndTable();
    }

    ImGui::TextDisabled("Tip: star/planet shading is now computed from a point-light at the origin (the star). ");
  }

  ImGui::End();
}


if (showHangar) {
  const bool wasOpen = showHangar;
  ImGui::Begin("Hangar / Livery", &showHangar);

  ImGui::TextDisabled("3D hangar preview is rendered offscreen (OpenGL FBO) and displayed via ImGui::Image().");
  ImGui::Separator();

  // Preview configuration
  {
    ImGui::SetNextItemWidth(240.0f);
    const int minPrev = 256;
    const int maxPrev = 1024;
    const int before = hangarPreviewSizePx;
    if (ImGui::SliderInt("Preview size (px)", &hangarPreviewSizePx, minPrev, maxPrev)) {
      hangarPreviewSizePx = std::clamp(hangarPreviewSizePx, minPrev, maxPrev);
    }
    if (hangarPreviewSizePx != before) {
      hangarTarget.ensureSize(hangarPreviewSizePx, hangarPreviewSizePx);
    }

    ImGui::SameLine();
    ImGui::Checkbox("Animate", &hangarAnimate);
    if (hangarAnimate) {
      ImGui::SameLine();
      ImGui::SetNextItemWidth(160.0f);
      ImGui::SliderFloat("Spin (rad/s)", &hangarSpinRadPerSec, 0.0f, 2.5f, "%.2f");
    }

    ImGui::SetNextItemWidth(240.0f);
    ImGui::SliderFloat("Zoom", &hangarZoom, 1.5f, 7.5f, "%.2f");
    ImGui::SameLine();
    ImGui::Checkbox("Auto-save on close", &liveryAutoSaveOnExit);
  }

  if (ImGui::BeginTable("hangar_table", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_SizingStretchProp)) {
    ImGui::TableNextRow();

    // --- Left: preview ---
    ImGui::TableNextColumn();
    {
      const float maxSize = 560.0f;
      const float avail = ImGui::GetContentRegionAvail().x;
      const float imgSize = std::max(140.0f, std::min(avail, maxSize));

      const GLuint prevTex = hangarTarget.color().handle();
      if (prevTex != 0) {
        ImGui::Image((ImTextureID)(intptr_t)prevTex, ImVec2(imgSize, imgSize), ImVec2(0, 1), ImVec2(1, 0));

        // Drag to rotate
        if (ImGui::IsItemActive() && ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
          const ImVec2 d = ImGui::GetMouseDragDelta(ImGuiMouseButton_Left);
          hangarYawDeg += d.x * 0.25f;
          hangarPitchDeg = std::clamp(hangarPitchDeg + d.y * 0.25f, -85.0f, 85.0f);
          ImGui::ResetMouseDragDelta(ImGuiMouseButton_Left);
          hangarAnimate = false;
        }

        // Scroll to zoom
        if (ImGui::IsItemHovered()) {
          const float wheel = ImGui::GetIO().MouseWheel;
          if (wheel != 0.0f) {
            hangarZoom = std::clamp(hangarZoom - wheel * 0.25f, 1.5f, 7.5f);
          }
        }

        ImGui::TextDisabled("Drag: rotate | Scroll: zoom");
      } else {
        ImGui::TextColored(ImVec4(1, 0.6f, 0.4f, 1), "Hangar preview unavailable (FBO init failed).");
      }
    }

    // --- Right: livery controls ---
    ImGui::TableNextColumn();
    {
      bool changed = false;

      ImGui::SeparatorText("Livery");
      {
        int pat = (int)liveryCfg.pattern;
        const char* items[] = {"Solid", "Stripes", "Hazard", "Camo", "Hex", "Digital"};
        if (ImGui::Combo("Pattern", &pat, items, (int)(sizeof(items) / sizeof(items[0])))) {
          pat = std::clamp(pat, 0, (int)render::LiveryPattern::Count - 1);
          liveryCfg.pattern = (render::LiveryPattern)pat;
          changed = true;
        }

        ImGui::SetNextItemWidth(240.0f);
        changed |= ImGui::InputScalar("Seed", ImGuiDataType_U64, &liveryCfg.seed);
        ImGui::SameLine();
        if (ImGui::SmallButton("Random")) {
          liveryCfg.seed = liveryRng.nextU64();
          changed = true;
        }

        if (ImGui::SmallButton("Random palette")) {
          core::SplitMix64 tmp(core::hashCombine(liveryCfg.seed, core::fnv1a64("palette")));
          auto rnd = [&]() -> float { return (float)tmp.range(0.08, 0.98); };
          liveryCfg.base[0] = rnd(); liveryCfg.base[1] = rnd(); liveryCfg.base[2] = rnd();
          liveryCfg.accent1[0] = rnd(); liveryCfg.accent1[1] = rnd(); liveryCfg.accent1[2] = rnd();
          liveryCfg.accent2[0] = rnd(); liveryCfg.accent2[1] = rnd(); liveryCfg.accent2[2] = rnd();
          changed = true;
        }
      }

      changed |= ImGui::ColorEdit3("Base", liveryCfg.base);
      changed |= ImGui::ColorEdit3("Accent A", liveryCfg.accent1);
      changed |= ImGui::ColorEdit3("Accent B", liveryCfg.accent2);

      ImGui::SeparatorText("Pattern params");
      ImGui::SetNextItemWidth(240.0f);
      changed |= ImGui::SliderFloat("Scale", &liveryCfg.scale, 0.25f, 4.0f, "%.2f");
      ImGui::SetNextItemWidth(240.0f);
      changed |= ImGui::SliderFloat("Angle (deg)", &liveryCfg.angleDeg, -90.0f, 90.0f, "%.1f");
      ImGui::SetNextItemWidth(240.0f);
      changed |= ImGui::SliderFloat("Detail", &liveryCfg.detail, 0.0f, 1.0f, "%.2f");
      ImGui::SetNextItemWidth(240.0f);
      changed |= ImGui::SliderFloat("Wear", &liveryCfg.wear, 0.0f, 1.0f, "%.2f");
      ImGui::SetNextItemWidth(240.0f);
      changed |= ImGui::SliderFloat("Contrast", &liveryCfg.contrast, 0.0f, 1.0f, "%.2f");
      changed |= ImGui::Checkbox("Decal", &liveryCfg.decal);

      ImGui::SeparatorText("Texture");
      {
        ImGui::SetNextItemWidth(240.0f);
        changed |= ImGui::SliderInt("Livery tex size", &liveryCfg.textureSizePx, 128, 1024);
        ImGui::SetNextItemWidth(240.0f);
        ImGui::SliderFloat("Regen throttle (s)", &shipLiveryRegenMinIntervalSec, 0.0f, 0.25f, "%.3f");
      }
      changed |= ImGui::Checkbox("Apply in world", &liveryCfg.applyInWorld);

      ImGui::Separator();
      if (ImGui::Button("Save livery")) {
        ui::saveLiveryToFile(liveryCfg, liveryPath);
      }
      ImGui::SameLine();
      if (ImGui::Button("Load livery")) {
        ui::LiveryConfig loaded = ui::makeDefaultLivery();
        if (ui::loadLiveryFromFile(liveryPath, loaded)) {
          liveryCfg = loaded;
          changed = true;
        }
      }
      ImGui::SameLine();
      if (ImGui::Button("Reset")) {
        liveryCfg = ui::makeDefaultLivery();
        changed = true;
      }
      ImGui::TextDisabled("File: %s", liveryPath.c_str());

      // Mini 2D preview of the generated albedo.
      if (shipLiveryTex.handle() != 0) {
        ImGui::SeparatorText("Albedo preview");
        ImGui::Image((ImTextureID)(intptr_t)shipLiveryTex.handle(), ImVec2(192, 192), ImVec2(0, 1), ImVec2(1, 0));
      }

      if (changed) shipLiveryDirty = true;

      // Regenerate the livery texture with a simple throttle so dragging sliders
      // doesn't re-run the generator at unbounded rates.
      if (shipLiveryDirty) {
        const bool allow = (!ImGui::IsAnyItemActive()) || ((timeRealSec - shipLiveryLastRegenTimeSec) > (double)shipLiveryRegenMinIntervalSec);
        if (allow) {
          rebuildShipLivery();
          shipLiveryLastRegenTimeSec = timeRealSec;
        }
      }
    }

    ImGui::EndTable();
  }

  ImGui::End();

  if (wasOpen && !showHangar && liveryAutoSaveOnExit) {
    ui::saveLiveryToFile(liveryCfg, liveryPath);
  }
}


if (showScanner) {
  ImGui::Begin("System Scanner");

  ImGui::Text("System: %s", currentSystem->stub.name.c_str());
  {
    const core::u64 sSeed = core::hashCombine(core::fnv1a64("star"), (core::u64)currentSystem->stub.seed);
    const auto& sIcon = spriteCache.get(render::SpriteKind::Star, sSeed, 48);
    ImGui::Image((ImTextureID)(intptr_t)sIcon.handle(), ImVec2(22, 22));
    if (ImGui::IsItemHovered()) {
      const auto& big = spriteCache.get(render::SpriteKind::Star, sSeed, 96);
      ImGui::BeginTooltip();
      ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(96, 96));
      if (worldUseProceduralSurfaces && worldUseSurfaceInUi) {
        const core::u64 surfSeed = core::hashCombine(
            core::hashCombine(core::fnv1a64("star_surface"), (core::u64)currentSystem->stub.seed),
            (core::u64)(int)currentSystem->star.cls);
        const auto& surf = surfaceTexCache.get(render::SurfaceKind::Star, surfSeed, worldSurfaceTexWidth);
        ImGui::Separator();
        ImGui::TextDisabled("Surface");
        ImGui::Image((ImTextureID)(intptr_t)surf.handle(), ImVec2(192, 96));
      }
      ImGui::EndTooltip();
    }
    ImGui::SameLine();
    ImGui::Text("Star: %s | Mass %.2f | Radius %.2f | Temp %.0fK",
                starClassName(currentSystem->star.cls),
                currentSystem->star.massSol,
                currentSystem->star.radiusSol,
                currentSystem->star.temperatureK);
  }

  const bool starScanned = (scannedKeys.find(scanKeyStar(currentSystem->stub.id)) != scannedKeys.end());

  // Quick UI helper: start a scan for a given target.
  auto uiStartScan = [&](TargetKind kind, std::size_t idx) {
    if (scanning) {
      toast(toasts, "Already scanning. (K cancels)", 2.0);
      return;
    }
    if (docked) {
      toast(toasts, "Undock to scan.", 1.8);
      return;
    }
    if (supercruiseState != SupercruiseState::Idle || fsdState != FsdState::Idle) {
      toast(toasts, "Scanning unavailable in supercruise / hyperspace.", 2.0);
      return;
    }

    target.kind = kind;
    target.index = idx;

    scanLockedTarget = target;
    scanLockedId = 0;
    scanLabel.clear();
    scanProgressSec = 0.0;

    bool ok = false;

    if (kind == TargetKind::Star) {
      scanLockedId = currentSystem->stub.id;
      scanDurationSec = 2.5;
      scanRangeKm = 1.0e18;
      scanLabel = std::string("Star scan: ") + starClassName(currentSystem->star.cls);
      ok = true;
    } else if (kind == TargetKind::Planet && idx < currentSystem->planets.size()) {
      const auto& p = currentSystem->planets[idx];
      scanLockedId = (core::u64)idx;
      scanDurationSec = 5.0;
      const double rKm = p.radiusEarth * kEARTH_RADIUS_KM;
      scanRangeKm = std::max(200000.0, rKm * 45.0);
      scanLabel = "Planet scan: " + p.name;
      ok = true;
    } else if (kind == TargetKind::Station && idx < currentSystem->stations.size()) {
      const auto& st = currentSystem->stations[idx];
      scanLockedId = st.id;
      scanDurationSec = 3.0;
      scanRangeKm = std::max(25000.0, st.commsRangeKm * 0.9);
      scanLabel = "Station scan: " + st.name;
      ok = true;
	    } else if (kind == TargetKind::Contact && idx < contacts.size()) {
	      const auto& c = contacts[idx];
	      if (c.alive) {
	        scanLockedId = c.id;
	        scanDurationSec = 3.5;
	        scanRangeKm = 140000.0;
	        scanLabel = "Ship scan: " + c.name;
	        ok = true;
	      }
	    } else if (kind == TargetKind::Signal && idx < signals.size()) {
	      const auto& s = signals[idx];
	      scanLockedId = s.id;
	      scanDurationSec = (s.type == SignalType::Resource) ? 3.5 : 4.0;
	      scanRangeKm = 200000.0;
	      const char* base = signalTypeName(s.type);
	      if (s.type == SignalType::Resource && s.hasResourcePlan) {
	        base = sim::resourceFieldKindName(s.resource.kind);
	      }
	      scanLabel = std::string("Signal scan: ") + base;
	      ok = true;
	    } else if (kind == TargetKind::Asteroid && idx < asteroids.size()) {
	      const auto& a = asteroids[idx];
	      scanLockedId = a.id;
	      scanDurationSec = 4.5;
	      scanRangeKm = 150000.0;
	      scanLabel = std::string("Asteroid prospect: ") + econ::commodityDef(a.yield).name;
	      ok = true;
    }

    if (ok) {
      scanning = true;
      toast(toasts, scanLabel + " (hold steady)...", 2.0);
    } else {
      toast(toasts, "Invalid scan target.", 2.0);
    }
  };

  ImGui::Separator();

  // Survey status
  int scannedPlanets = 0;
  for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
    if (scannedKeys.find(scanKeyPlanet(currentSystem->stub.id, i)) != scannedKeys.end()) ++scannedPlanets;
  }
  int scannedStations = 0;
  for (const auto& st : currentSystem->stations) {
    if (scannedKeys.find(scanKeyStation(st.id)) != scannedKeys.end()) ++scannedStations;
  }

  ImGui::Text("Survey: Star %s | Planets %d/%d | Stations %d/%d",
              starScanned ? "OK" : "UNSCANNED",
              scannedPlanets, (int)currentSystem->planets.size(),
              scannedStations, (int)currentSystem->stations.size());

  if (!starScanned) {
    if (ImGui::Button("Target Star")) {
      target.kind = TargetKind::Star;
      target.index = 0;
    }
    ImGui::SameLine();
    if (ImGui::Button("Scan Star")) {
      uiStartScan(TargetKind::Star, 0);
    }
  }

  if (ImGui::BeginTabBar("ScannerTabs")) {
    if (ImGui::BeginTabItem("Catalog")) {
      ImGui::Separator();
      ImGui::Text("Planets");

  for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
    const auto& p = currentSystem->planets[i];
    const bool scanned = (scannedKeys.find(scanKeyPlanet(currentSystem->stub.id, i)) != scannedKeys.end());

    ImGui::PushID((int)i);

    {
      const core::u64 pSeed = core::hashCombine(core::hashCombine(core::fnv1a64("planet"), (core::u64)currentSystem->stub.seed), (core::u64)i);
      const auto& pIcon = spriteCache.get(render::SpriteKind::Planet, pSeed, 48);
      ImGui::Image((ImTextureID)(intptr_t)pIcon.handle(), ImVec2(18, 18));
      if (ImGui::IsItemHovered()) {
        const auto& big = spriteCache.get(render::SpriteKind::Planet, pSeed, 96);
        ImGui::BeginTooltip();
        ImGui::TextUnformatted(p.name.c_str());
        ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(96, 96));
        if (worldUseProceduralSurfaces && worldUseSurfaceInUi) {
          const core::u64 surfSeed = core::hashCombine(
              core::hashCombine(core::fnv1a64("planet_surface"), (core::u64)currentSystem->stub.seed),
              (core::u64)i);
          const auto& surf = surfaceTexCache.get(planetSurfaceKind(p.type), surfSeed, worldSurfaceTexWidth);
          ImGui::Separator();
          ImGui::TextDisabled("Surface");
          ImGui::Image((ImTextureID)(intptr_t)surf.handle(), ImVec2(192, 96));

          // Optional cloud layer preview (alpha mask)
          if (worldCloudsEnabled) {
            core::u64 cloudSeed = core::hashCombine(surfSeed, core::fnv1a64("clouds"));
            cloudSeed = core::hashCombine(cloudSeed, (core::u64)(int)p.type);
            const auto& cloud = surfaceTexCache.get(render::SurfaceKind::Clouds, cloudSeed, worldSurfaceTexWidth);
            ImGui::TextDisabled("Clouds");
            ImGui::Image((ImTextureID)(intptr_t)cloud.handle(), ImVec2(192, 96));
          }
        }
        ImGui::EndTooltip();
      }
      ImGui::SameLine();
      ImGui::Text("%s %s", p.name.c_str(), scanned ? "" : "(unscanned)");
    }
    ImGui::SameLine();
    if (ImGui::SmallButton("Target##p")) {
      target.kind = TargetKind::Planet;
      target.index = i;
    }
    ImGui::SameLine();
    if (ImGui::SmallButton("Scan##p")) {
      uiStartScan(TargetKind::Planet, i);
    }

    if (scanned) {
      const double radiusKm = p.radiusEarth * kEARTH_RADIUS_KM;
      const double periodDays = p.orbit.periodDays;
      ImGui::TextDisabled("Type: %s | Radius: %.0f km | a: %.2f AU | Period: %.1f d",
                          planetTypeName(p.type), radiusKm, p.orbit.semiMajorAxisAU, periodDays);
    } else {
      ImGui::TextDisabled("Type: ??? | Radius: ??? | Orbit: ???");
    }

    ImGui::Separator();
    ImGui::PopID();
  }

  ImGui::Text("Stations");

  for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
    const auto& st = currentSystem->stations[i];
    const bool scanned = (scannedKeys.find(scanKeyStation(st.id)) != scannedKeys.end());

    ImGui::PushID((int)(1000 + i));

    {
      const core::u64 stSeed = core::hashCombine(core::fnv1a64("station"), (core::u64)st.id);
      const auto& stIcon = spriteCache.get(render::SpriteKind::Station, stSeed, 32);
      ImGui::Image((ImTextureID)(intptr_t)stIcon.handle(), ImVec2(18, 18));
      if (ImGui::IsItemHovered()) {
        const auto& big = spriteCache.get(render::SpriteKind::Station, stSeed, 96);
        ImGui::BeginTooltip();
        ImGui::TextUnformatted(st.name.c_str());
        ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(96, 96));
        ImGui::EndTooltip();
      }
      ImGui::SameLine();
      ImGui::Text("%s %s", st.name.c_str(), scanned ? "" : "(unscanned)");
    }
    ImGui::SameLine();
    if (ImGui::SmallButton("Target##s")) {
      target.kind = TargetKind::Station;
      target.index = i;
    }
    ImGui::SameLine();
    if (ImGui::SmallButton("Scan##s")) {
      uiStartScan(TargetKind::Station, i);
    }

    if (scanned) {
      ImGui::TextDisabled("Type: %s | Faction: %s | Fee: %.1f%%",
                          stationTypeName(st.type),
                          factionName(st.factionId).c_str(),
                          st.feeRate * 100.0);
    } else {
      ImGui::TextDisabled("Type: ??? | Faction: ??? | Services: ???");
    }

    ImGui::Separator();
    ImGui::PopID();
  }

	  // Signals / salvage / mining
	  ImGui::Separator();
	  ImGui::Text("Signal Sources");
	  if (signals.empty()) {
	    ImGui::TextDisabled("(none detected)");
	  } else {
	    for (std::size_t i = 0; i < signals.size(); ++i) {
	      const auto& s = signals[i];
	      const double distKm = (s.posKm - ship.positionKm()).length();
	      const double minsLeft = (s.expireDay > 0.0) ? (s.expireDay - timeDays) * 1440.0 : 0.0;
	      ImGui::PushID((int)i + 9000);
      const char* mtag = "";
      for (const auto& m : missions) {
        if (m.completed || m.failed) continue;
        if (m.type != sim::MissionType::Salvage) continue;
        if (m.targetNpcId != s.id) continue;
        mtag = " [MISSION]";
        break;
      }

      // Procedural signal icon
      {
        const core::u64 sigSeed = core::hashCombine(core::hashCombine(core::fnv1a64("signal"), (core::u64)s.id), (core::u64)(int)s.type);
        const auto& sigIcon = spriteCache.get(render::SpriteKind::Signal, sigSeed, 32);
        ImGui::Image((ImTextureID)(intptr_t)sigIcon.handle(), ImVec2(18, 18));
        if (ImGui::IsItemHovered()) {
          const auto& big = spriteCache.get(render::SpriteKind::Signal, sigSeed, 96);
          ImGui::BeginTooltip();
	          const char* name = signalTypeName(s.type);
	          if (s.type == SignalType::Resource && s.hasResourcePlan) {
	            name = sim::resourceFieldKindName(s.resource.kind);
	          }
	          ImGui::Text("%s Signal", name);
          ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(96, 96));
	          if (s.type == SignalType::Resource && s.hasResourcePlan) {
	            const auto& p = s.resource;
	            const auto& a = econ::commodityDef(p.primaryYield);
	            const auto& b = econ::commodityDef(p.secondaryYield);
	            const int pa = (int)std::lround(std::clamp(p.primaryChance, 0.0, 1.0) * 100.0);
	            const int pb = 100 - pa;
	            const int rich = (int)std::lround(std::clamp(p.richness01, 0.0, 1.0) * 100.0);
	            ImGui::Separator();
	            ImGui::TextDisabled("Richness %d%%", rich);
	            ImGui::TextDisabled("%s ~%d%% / %s ~%d%%", a.name, pa, b.name, pb);
	          }
          ImGui::EndTooltip();
        }
        ImGui::SameLine();
      }

	      const char* name = signalTypeName(s.type);
	      if (s.type == SignalType::Resource && s.hasResourcePlan) {
	        name = sim::resourceFieldKindName(s.resource.kind);
	      }
	      ImGui::Text("%s%s%s | %.0f km%s", name, s.resolved ? " (resolved)" : "", mtag, distKm,
	                  (s.expireDay > 0.0) ? (std::string(" | ") + std::to_string((int)std::round(std::max(0.0, minsLeft)))
	                                        + " min").c_str()
	                                     : "");
	      ImGui::SameLine();
	      if (ImGui::SmallButton("Target##sig")) {
	        target.kind = TargetKind::Signal;
	        target.index = i;
	      }
	      ImGui::SameLine();
	      if (ImGui::SmallButton("Scan##sig")) {
	        uiStartScan(TargetKind::Signal, i);
	      }
	      ImGui::PopID();
	    }
	  }

	  ImGui::Separator();
	  ImGui::Text("Asteroids");
	  if (asteroids.empty()) {
	    ImGui::TextDisabled("(no nearby belts)");
	  } else {
	    for (std::size_t i = 0; i < asteroids.size(); ++i) {
	      const auto& a = asteroids[i];
	      if (a.remainingUnits <= 1e-3) continue;
	      const double distKm = (a.posKm - ship.positionKm()).length();
	      ImGui::PushID((int)i + 12000);

      // Procedural asteroid icon
      {
        const core::u64 aSeed = core::hashCombine(core::hashCombine(core::fnv1a64("asteroid"), (core::u64)a.id), (core::u64)a.yield);
        const auto& aIcon = spriteCache.get(render::SpriteKind::Asteroid, aSeed, 32);
        ImGui::Image((ImTextureID)(intptr_t)aIcon.handle(), ImVec2(18, 18));
        if (ImGui::IsItemHovered()) {
          const auto& big = spriteCache.get(render::SpriteKind::Asteroid, aSeed, 96);
          ImGui::BeginTooltip();
          ImGui::Text("Asteroid [%s]", econ::commodityDef(a.yield).name);
          ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(96, 96));
          ImGui::EndTooltip();
        }
        ImGui::SameLine();
      }
	      ImGui::Text("%s | %.0f u | %.0f km", econ::commodityDef(a.yield).name,
	                  std::round(std::max(0.0, a.remainingUnits)), distKm);
	      ImGui::SameLine();
	      if (ImGui::SmallButton("Target##ast")) {
	        target.kind = TargetKind::Asteroid;
	        target.index = i;
	      }
	      ImGui::SameLine();
	      if (ImGui::SmallButton("Prospect##ast")) {
	        uiStartScan(TargetKind::Asteroid, i);
	      }
	      ImGui::PopID();
	    }
	  }

	  ImGui::Separator();
	  ImGui::Text("Floating Cargo");
	  if (floatingCargo.empty()) {
	    ImGui::TextDisabled("(none)");
	  } else {
	    for (std::size_t i = 0; i < floatingCargo.size(); ++i) {
	      const auto& pod = floatingCargo[i];
	      const double distKm = (pod.posKm - ship.positionKm()).length();
	      const double minsLeft = (pod.expireDay > 0.0) ? (pod.expireDay - timeDays) * 1440.0 : 0.0;
	      ImGui::PushID((int)i + 15000);

      // Procedural cargo pod icon
      {
        const core::u64 cSeed = core::hashCombine(core::hashCombine(core::fnv1a64("cargo"), (core::u64)pod.id), (core::u64)pod.commodity);
        const auto& cIcon = spriteCache.get(render::SpriteKind::Cargo, cSeed, 32);
        ImGui::Image((ImTextureID)(intptr_t)cIcon.handle(), ImVec2(18, 18));
        if (ImGui::IsItemHovered()) {
          const auto& big = spriteCache.get(render::SpriteKind::Cargo, cSeed, 96);
          ImGui::BeginTooltip();
          ImGui::Text("%s Pod", econ::commodityDef(pod.commodity).name);
          ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(96, 96));
          ImGui::EndTooltip();
        }
        ImGui::SameLine();
      }
	      const char* ptag = (pod.missionId != 0) ? " [MISSION]" : "";
	      ImGui::Text("%s%s x%.0f | %.0f km | %d min", econ::commodityDef(pod.commodity).name, ptag,
	                  std::round(std::max(0.0, pod.units)), distKm, (int)std::round(std::max(0.0, minsLeft)));
	      ImGui::SameLine();
	      if (ImGui::SmallButton("Target##pod")) {
	        target.kind = TargetKind::Cargo;
	        target.index = i;
	      }
	      ImGui::PopID();
	    }
	  }

	  ImGui::Separator();
	  ImGui::Text("Exploration data bank: %.0f cr", explorationDataCr);
  if (docked && selectedStationIndex >= 0 && selectedStationIndex < (int)currentSystem->stations.size()) {
    if (explorationDataCr > 0.0) {
      if (ImGui::Button("Sell exploration data here")) {
        credits += explorationDataCr;
        toast(toasts, "Sold exploration data +" + std::to_string((int)explorationDataCr) + " cr.", 2.5);
        explorationDataCr = 0.0;
      }
    } else {
      ImGui::TextDisabled("No data to sell.");
    }
  } else {
    ImGui::TextDisabled("Dock at a station to sell exploration data.");
  }

      // End Catalog tab.
      ImGui::EndTabItem();
    }

    // System Map tab (Orrery 2.0: infinite-canvas pan/zoom + time scrub + atlas icons)
    if (ImGui::BeginTabItem("System Map")) {
      struct Vec2AU {
        double x{0.0}, y{0.0};
      };
      struct SystemMapState {
        core::u64 systemId{0};
        Vec2AU centerAU{0.0, 0.0};
        float zoom{1.0f};
        bool followShip{false};

        bool showGrid{true};
        bool showOrbits{true};
        bool showStations{true};
        bool showContacts{true};
        bool showSignals{true};
        bool showAsteroids{true};
        bool showCargo{true};

        double previewOffsetDays{0.0};
        bool previewAnimate{false};
        float previewSpeedDaysPerSec{0.35f};
      };
      static SystemMapState map{};

      // Reset view when changing star system.
      if (map.systemId != currentSystem->stub.id) {
        map.systemId = currentSystem->stub.id;
        map.centerAU = {0.0, 0.0};
        map.zoom = 1.0f;
        map.followShip = false;
        map.previewOffsetDays = 0.0;
        map.previewAnimate = false;
      }

      if (map.previewAnimate && !paused) {
        map.previewOffsetDays += (double)map.previewSpeedDaysPerSec * dtReal;
        map.previewOffsetDays = std::clamp(map.previewOffsetDays, -2.0, 120.0);
      }

      const double previewTimeDays = timeDays + map.previewOffsetDays;

      // Shared atlas texture for all map icons.
      const auto& atlasTex = hudAtlas.texture();
      const ImTextureID atlasId = (ImTextureID)(intptr_t)atlasTex.handle();

      // Mirror the keyboard 'H' supercruise toggle as a UI action.
      auto uiToggleSupercruise = [&]() {
        if (docked) { toast(toasts, "Undock to use supercruise.", 1.8); return; }
        if (fsdState != FsdState::Idle) { toast(toasts, "Supercruise unavailable in hyperspace.", 2.0); return; }

        if (supercruiseState == SupercruiseState::Idle) {
          if (target.kind == TargetKind::None) {
            toast(toasts, "No nav target set.", 1.8);
            return;
          }

          double nearestPirateKm = 1e99;
          for (const auto& c : contacts) {
            if (!c.alive || c.role != ContactRole::Pirate) continue;
            nearestPirateKm = std::min(nearestPirateKm, (c.ship.positionKm() - ship.positionKm()).length());
          }
          const bool pirateNearby = (nearestPirateKm < 160000.0);

          supercruiseState = SupercruiseState::Charging;
          supercruiseChargeRemainingSec = kSupercruiseChargeSec;
          supercruiseDropRequested = false;
          interdiction = sim::InterdictionState{};
          interdictionSubmitRequested = false;
          interdictionPirateName.clear();
          interdictionPirateStrength = 1.0;
          autopilot = false;
          dockingComputer.reset();
          scanning = false;
          scanProgressSec = 0.0;
          navAutoRun = false;

          toast(toasts,
                pirateNearby ? "Supercruise charging... (pirate signals nearby)" : "Supercruise charging...",
                1.8);
        } else if (supercruiseState == SupercruiseState::Charging) {
          supercruiseState = SupercruiseState::Idle;
          supercruiseChargeRemainingSec = 0.0;
          interdiction = sim::InterdictionState{};
          interdictionSubmitRequested = false;
          interdictionPirateName.clear();
          interdictionPirateStrength = 1.0;
          toast(toasts, "Supercruise canceled.", 1.4);
        } else if (supercruiseState == SupercruiseState::Active) {
          if (sim::interdictionInProgress(interdiction)) {
            interdictionSubmitRequested = true;
            toast(toasts, "Submitted to interdiction.", 1.4);
          } else {
            supercruiseDropRequested = true;
            toast(toasts, "Drop requested.", 1.2);
          }
        } else {
          toast(toasts, "Supercruise cooling down...", 1.6);
        }
      };

      auto centerOnStar = [&]() {
        map.followShip = false;
        map.centerAU = {0.0, 0.0};
      };

      auto shipPosAU2 = [&]() -> Vec2AU {
        const math::Vec3d au = ship.positionKm() / kAU_KM;
        return {au.x, au.z};
      };

      ImGui::TextDisabled("Pan: RMB/MMB drag | Zoom: wheel (cursor anchored) | Double-click: center | Right-click: context menu");

      if (ImGui::Button("Center Star")) { centerOnStar(); }
      ImGui::SameLine();
      if (ImGui::Button("Center Ship")) { map.centerAU = shipPosAU2(); map.followShip = true; }
      ImGui::SameLine();
      ImGui::Checkbox("Follow ship", &map.followShip);
      ImGui::SameLine();
      ImGui::SetNextItemWidth(180.0f);
      ImGui::SliderFloat("Zoom", &map.zoom, 0.12f, 14.0f, "%.2fx", ImGuiSliderFlags_Logarithmic);
      ImGui::SameLine();
      ImGui::Checkbox("Grid", &map.showGrid);

      ImGui::Checkbox("Orbits", &map.showOrbits);
      ImGui::SameLine();
      ImGui::Checkbox("Stations", &map.showStations);
      ImGui::SameLine();
      ImGui::Checkbox("Contacts", &map.showContacts);
      ImGui::SameLine();
      ImGui::Checkbox("Signals", &map.showSignals);
      ImGui::SameLine();
      ImGui::Checkbox("Asteroids", &map.showAsteroids);
      ImGui::SameLine();
      ImGui::Checkbox("Cargo", &map.showCargo);

      ImGui::Separator();

      // Orbit preview time scrub (does not change sim-time; purely for map visualization).
      {
        ImGui::TextDisabled("Orbit preview:");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(240.0f);
        double tOff = map.previewOffsetDays;
        const double tMin = -2.0;
        const double tMax = 120.0;
        if (ImGui::SliderScalar("##preview_days", ImGuiDataType_Double, &tOff, &tMin, &tMax, "%+.2f days", ImGuiSliderFlags_Logarithmic)) {
          map.previewOffsetDays = tOff;
        }
        ImGui::SameLine();
        ImGui::Checkbox("Animate##preview", &map.previewAnimate);
        ImGui::SameLine();
        ImGui::SetNextItemWidth(120.0f);
        ImGui::SliderFloat("Speed (d/s)##preview", &map.previewSpeedDaysPerSec, 0.02f, 5.0f, "%.2f", ImGuiSliderFlags_Logarithmic);
        ImGui::SameLine();
        if (ImGui::Button("Reset##preview")) { map.previewOffsetDays = 0.0; map.previewAnimate = false; }
        ImGui::TextDisabled("Preview time: t %+.2f days", (float)map.previewOffsetDays);
      }

      ImGui::Separator();

      // Layout: map canvas + details panel.
      if (ImGui::BeginTable("sysmap2_layout", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("Map", ImGuiTableColumnFlags_WidthStretch, 0.68f);
        ImGui::TableSetupColumn("Details", ImGuiTableColumnFlags_WidthStretch, 0.32f);
        ImGui::TableNextRow();

        // ---- Map column ----
        ImGui::TableNextColumn();
        ImGui::BeginChild("sysmap2_canvas", ImVec2(0, 0), true, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);

        ImDrawList* draw = ImGui::GetWindowDrawList();
        const ImVec2 p0 = ImGui::GetCursorScreenPos();
        ImVec2 canvasSize = ImGui::GetContentRegionAvail();
        canvasSize.x = std::max(canvasSize.x, 60.0f);
        canvasSize.y = std::max(canvasSize.y, 60.0f);
        const ImVec2 p1(p0.x + canvasSize.x, p0.y + canvasSize.y);
        const ImVec2 centerScreen(p0.x + canvasSize.x * 0.5f, p0.y + canvasSize.y * 0.5f);

        ImGui::InvisibleButton("sysmap2_invis", canvasSize,
                               ImGuiButtonFlags_MouseButtonLeft |
                               ImGuiButtonFlags_MouseButtonRight |
                               ImGuiButtonFlags_MouseButtonMiddle);
        const bool canvasHovered = ImGui::IsItemHovered();
        const bool canvasActive = ImGui::IsItemActive();
        const ImVec2 mouse = io.MousePos;

        // Determine system extent (AU) for scale.
        double maxAU = 0.25;
        for (const auto& p : currentSystem->planets) {
          maxAU = std::max(maxAU, p.orbit.semiMajorAxisAU * (1.0 + p.orbit.eccentricity));
        }
        for (const auto& st : currentSystem->stations) {
          maxAU = std::max(maxAU, st.orbit.semiMajorAxisAU * (1.0 + st.orbit.eccentricity));
        }
        // Include ship radius so you can always find yourself if you're far out.
        {
          const math::Vec3d shipAU = ship.positionKm() / kAU_KM;
          maxAU = std::max(maxAU, std::sqrt(shipAU.x * shipAU.x + shipAU.z * shipAU.z));
        }

        const float pad = 16.0f;
        const float rad = std::max(20.0f, 0.5f * std::min(canvasSize.x, canvasSize.y) - pad);
        const float baseScale = (float)(rad / std::max(0.01, maxAU));

        auto clampZoom = [](float z) { return std::clamp(z, 0.08f, 22.0f); };
        map.zoom = clampZoom(map.zoom);
        float pxPerAU = baseScale * map.zoom;

        // Follow ship unless the user is actively panning.
        if (map.followShip && !canvasActive) {
          map.centerAU = shipPosAU2();
        }

        auto worldToScreen = [&](const Vec2AU& au) -> ImVec2 {
          const float x = (float)((au.x - map.centerAU.x) * (double)pxPerAU);
          const float y = (float)((au.y - map.centerAU.y) * (double)pxPerAU);
          return ImVec2(centerScreen.x + x, centerScreen.y - y);
        };

        auto screenToWorld = [&](const ImVec2& sp, float scale) -> Vec2AU {
          const double wx = map.centerAU.x + (double)(sp.x - centerScreen.x) / (double)scale;
          const double wy = map.centerAU.y - (double)(sp.y - centerScreen.y) / (double)scale;
          return {wx, wy};
        };

        // Zoom (anchored under cursor).
        if (canvasHovered && std::abs(io.MouseWheel) > 1e-6f) {
          const Vec2AU wUnder = screenToWorld(mouse, pxPerAU);
          const float k = (io.MouseWheel > 0.0f) ? 1.18f : 0.85f;
          map.zoom = clampZoom(map.zoom * k);

          const float newPxPerAU = baseScale * map.zoom;
          map.centerAU.x = wUnder.x - (double)(mouse.x - centerScreen.x) / (double)newPxPerAU;
          map.centerAU.y = wUnder.y + (double)(mouse.y - centerScreen.y) / (double)newPxPerAU;
          pxPerAU = newPxPerAU;
        }

        // Pan (drag) - disables follow-ship.
        if (canvasHovered && (ImGui::IsMouseDragging(ImGuiMouseButton_Right, 0.0f) ||
                              ImGui::IsMouseDragging(ImGuiMouseButton_Middle, 0.0f))) {
          map.followShip = false;
          map.centerAU.x -= (double)io.MouseDelta.x / (double)pxPerAU;
          map.centerAU.y += (double)io.MouseDelta.y / (double)pxPerAU;
        }

        // Background
        draw->AddRectFilled(p0, p1, IM_COL32(0, 0, 0, 90), 6.0f);
        draw->AddRect(p0, p1, IM_COL32(120, 140, 170, 120), 6.0f, 0, 1.0f);

        // Nice-number helper for grid / scale.
        auto niceNum = [](double x) -> double {
          if (x <= 0.0) return 1.0;
          const double expv = std::floor(std::log10(x));
          const double f = x / std::pow(10.0, expv);
          double nf = 1.0;
          if (f < 1.5) nf = 1.0;
          else if (f < 3.0) nf = 2.0;
          else if (f < 7.0) nf = 5.0;
          else nf = 10.0;
          return nf * std::pow(10.0, expv);
        };

        // Grid (anchored on star at (0,0) AU)
        if (map.showGrid) {
          const float halfAUx = (float)(canvasSize.x * 0.5f / pxPerAU);
          const float halfAUy = (float)(canvasSize.y * 0.5f / pxPerAU);
          const float halfAU = std::max(halfAUx, halfAUy);
          const double stepAU = niceNum((double)halfAU / 5.5);
          const ImVec2 starP = worldToScreen({0.0, 0.0});

          // Axes
          draw->AddLine(ImVec2(p0.x, starP.y), ImVec2(p1.x, starP.y), IM_COL32(120, 140, 170, 50), 1.0f);
          draw->AddLine(ImVec2(starP.x, p0.y), ImVec2(starP.x, p1.y), IM_COL32(120, 140, 170, 50), 1.0f);

          // Concentric circles
          for (int i = 1; i <= 12; ++i) {
            const double rAU = stepAU * (double)i;
            const float rPx = (float)(rAU * (double)pxPerAU);
            if (rPx > std::max(canvasSize.x, canvasSize.y)) break;
            const ImU32 col = (i % 2 == 0) ? IM_COL32(120, 140, 170, 38) : IM_COL32(120, 140, 170, 26);
            draw->AddCircle(starP, rPx, col, 0, 1.0f);

            if (i == 1 || (i % 2 == 0)) {
              char buf[32];
              std::snprintf(buf, sizeof(buf), "%.3g AU", rAU);
              draw->AddText(ImVec2(starP.x + rPx + 4.0f, starP.y + 2.0f), IM_COL32(160, 190, 220, 120), buf);
            }
          }

          // Scale bar (bottom-left)
          double barAU = stepAU;
          float barPx = (float)(barAU * (double)pxPerAU);
          while (barPx < 70.0f) { barAU *= 2.0; barPx = (float)(barAU * (double)pxPerAU); }
          while (barPx > 180.0f) { barAU *= 0.5; barPx = (float)(barAU * (double)pxPerAU); }

          const ImVec2 b0(p0.x + 12.0f, p1.y - 18.0f);
          const ImVec2 b1(b0.x + barPx, b0.y);
          draw->AddLine(b0, b1, IM_COL32(200, 220, 255, 180), 2.0f);
          draw->AddLine(ImVec2(b0.x, b0.y - 4.0f), ImVec2(b0.x, b0.y + 4.0f), IM_COL32(200, 220, 255, 180), 2.0f);
          draw->AddLine(ImVec2(b1.x, b1.y - 4.0f), ImVec2(b1.x, b1.y + 4.0f), IM_COL32(200, 220, 255, 180), 2.0f);
          char sb[32];
          std::snprintf(sb, sizeof(sb), "%.3g AU", barAU);
          draw->AddText(ImVec2(b0.x, b0.y - 16.0f), IM_COL32(200, 220, 255, 180), sb);
        }

        // Hover bookkeeping
        struct Hover {
          bool has{false};
          TargetKind kind{TargetKind::None};
          std::size_t idx{0};
          render::SpriteKind sprite{render::SpriteKind::Star};
          core::u64 seed{0};
          std::string label{};
          ImVec2 screenPos{};
          float iconSz{0.0f};
          float d2{1e30f};
        } hover;

        auto considerHover = [&](TargetKind kind, std::size_t idx, const ImVec2& pos, float iconSz,
                                 render::SpriteKind spriteKind, core::u64 spriteSeed,
                                 const std::string& label) {
          const float hx = iconSz * 0.5f;
          if (mouse.x < pos.x - hx || mouse.x > pos.x + hx || mouse.y < pos.y - hx || mouse.y > pos.y + hx) return;
          const float dx = mouse.x - pos.x;
          const float dy = mouse.y - pos.y;
          const float d2 = dx * dx + dy * dy;
          if (d2 < hover.d2) {
            hover.has = true;
            hover.kind = kind;
            hover.idx = idx;
            hover.sprite = spriteKind;
            hover.seed = spriteSeed;
            hover.label = label;
            hover.screenPos = pos;
            hover.iconSz = iconSz;
            hover.d2 = d2;
          }
        };

        auto drawIcon = [&](const ImVec2& pos, float iconSz,
                            render::SpriteKind kind, core::u64 seed,
                            ImU32 tint, bool selected) {
          const auto uv = hudAtlas.get(kind, seed);
          const ImVec2 b0(pos.x - iconSz * 0.5f, pos.y - iconSz * 0.5f);
          const ImVec2 b1(pos.x + iconSz * 0.5f, pos.y + iconSz * 0.5f);

          if (selected) {
            const float ring = iconSz * 0.72f;
            const float pulse = 0.5f + 0.5f * std::sin((float)timeRealSec * 4.0f);
            draw->AddCircle(pos, ring, IM_COL32(255, 210, 120, (int)(120.0f + 80.0f * pulse)), 0, 2.0f);
          }

          // Subtle outline for contrast
          draw->AddImage(atlasId,
                         ImVec2(b0.x - 1, b0.y - 1), ImVec2(b1.x + 1, b1.y + 1),
                         ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1),
                         IM_COL32(0, 0, 0, 140));
          draw->AddImage(atlasId, b0, b1,
                         ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1),
                         tint);
        };

        // Orbits (planets + stations)
        if (map.showOrbits) {
          for (const auto& p : currentSystem->planets) {
            const int N = 96;
            std::array<ImVec2, N + 1> pts;
            for (int i = 0; i <= N; ++i) {
              const double t = previewTimeDays + (double)i / (double)N * p.orbit.periodDays;
              const math::Vec3d posAU3 = sim::orbitPosition3DAU(p.orbit, t);
              pts[(std::size_t)i] = worldToScreen({posAU3.x, posAU3.z});
            }
            draw->AddPolyline(pts.data(), (int)pts.size(), IM_COL32(120, 140, 170, 55), false, 1.0f);
          }
          for (const auto& st : currentSystem->stations) {
            const int N = 72;
            std::array<ImVec2, N + 1> pts;
            for (int i = 0; i <= N; ++i) {
              const double t = previewTimeDays + (double)i / (double)N * st.orbit.periodDays;
              const math::Vec3d posAU3 = sim::orbitPosition3DAU(st.orbit, t);
              pts[(std::size_t)i] = worldToScreen({posAU3.x, posAU3.z});
            }
            draw->AddPolyline(pts.data(), (int)pts.size(), IM_COL32(100, 130, 170, 38), false, 1.0f);
          }
        }

        const float iconScale = std::clamp(0.95f + 0.10f * std::log2(std::max(0.20f, map.zoom)), 0.70f, 1.35f);

        // Star at origin
        {
          const core::u64 sSeed = core::hashCombine(core::fnv1a64("star"), (core::u64)currentSystem->stub.seed);
          const ImVec2 sp = worldToScreen({0.0, 0.0});
          const bool sel = (target.kind == TargetKind::Star);
          drawIcon(sp, 30.0f * iconScale, render::SpriteKind::Star, sSeed, IM_COL32(255, 255, 255, 230), sel);
          considerHover(TargetKind::Star, 0, sp, 30.0f * iconScale, render::SpriteKind::Star, sSeed,
                        std::string("Star (") + starClassName(currentSystem->star.cls) + ")");
        }

        // Ship marker (Ship sprite + forward tick)
        {
          const Vec2AU shipAU = shipPosAU2();
          const ImVec2 sp = worldToScreen(shipAU);
          const core::u64 sSeed = core::hashCombine(core::fnv1a64("player_ship"), seed);
          drawIcon(sp, 18.0f * iconScale, render::SpriteKind::Ship, sSeed, IM_COL32(210, 230, 255, 230), false);

          const math::Vec3d f = ship.forward().normalized();
          const ImVec2 dir((float)f.x, (float)-f.z);
          const float len = 16.0f;
          draw->AddLine(sp, ImVec2(sp.x + dir.x * len, sp.y + dir.y * len), IM_COL32(210, 230, 255, 160), 1.5f);
          draw->AddCircleFilled(sp, 2.0f, IM_COL32(210, 230, 255, 190));
        }

        // Planets
        for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
          const auto& p = currentSystem->planets[i];
          const math::Vec3d posAU3 = sim::orbitPosition3DAU(p.orbit, previewTimeDays);
          const ImVec2 sp = worldToScreen({posAU3.x, posAU3.z});
          const core::u64 pSeed = core::hashCombine(core::hashCombine(core::fnv1a64("planet"), (core::u64)currentSystem->stub.seed), (core::u64)i);
          const bool sel = (target.kind == TargetKind::Planet && target.index == i);

          ImU32 tint = IM_COL32(220, 230, 255, 220);
          if (p.type == sim::PlanetType::Ocean) tint = IM_COL32(120, 190, 255, 230);
          else if (p.type == sim::PlanetType::Desert) tint = IM_COL32(240, 205, 135, 230);
          else if (p.type == sim::PlanetType::Ice) tint = IM_COL32(210, 245, 255, 230);
          else if (p.type == sim::PlanetType::GasGiant) tint = IM_COL32(195, 210, 255, 230);

          drawIcon(sp, 20.0f * iconScale, render::SpriteKind::Planet, pSeed, tint, sel);
          considerHover(TargetKind::Planet, i, sp, 20.0f * iconScale, render::SpriteKind::Planet, pSeed, p.name);
        }

        // Stations
        if (map.showStations) {
          for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
            const auto& st = currentSystem->stations[i];
            const math::Vec3d posAU3 = sim::orbitPosition3DAU(st.orbit, previewTimeDays);
            const ImVec2 sp = worldToScreen({posAU3.x, posAU3.z});
            const core::u64 stSeed = core::hashCombine(core::fnv1a64("station"), (core::u64)st.id);
            const bool sel = (target.kind == TargetKind::Station && target.index == i);
            drawIcon(sp, 17.0f * iconScale, render::SpriteKind::Station, stSeed, IM_COL32(210, 210, 220, 225), sel);
            considerHover(TargetKind::Station, i, sp, 17.0f * iconScale, render::SpriteKind::Station, stSeed, st.name);
          }
        }

        // Contacts / signals / asteroids / cargo
        if (map.showContacts) {
          for (std::size_t i = 0; i < contacts.size(); ++i) {
            const auto& ctc = contacts[i];
            if (!ctc.alive) continue;
            const math::Vec3d posAU3 = ctc.ship.positionKm() / kAU_KM;
            const ImVec2 sp = worldToScreen({posAU3.x, posAU3.z});
            const core::u64 cSeed = core::hashCombine(core::fnv1a64("ship"), ctc.id);
            const bool sel = (target.kind == TargetKind::Contact && target.index == i);

            ImU32 tint = IM_COL32(180, 210, 235, 210);
            if (ctc.role == ContactRole::Pirate) tint = IM_COL32(255, 80, 80, 230);
            if (ctc.role == ContactRole::Police) tint = IM_COL32(90, 190, 255, 230);
            if (ctc.role == ContactRole::Trader) tint = IM_COL32(120, 245, 140, 230);

            drawIcon(sp, 14.0f * iconScale, render::SpriteKind::Ship, cSeed, tint, sel);
            considerHover(TargetKind::Contact, i, sp, 14.0f * iconScale, render::SpriteKind::Ship, cSeed,
                          ctc.name + " [" + contactRoleName(ctc.role) + "]");
          }
        }

        if (map.showSignals) {
          for (std::size_t i = 0; i < signals.size(); ++i) {
            const auto& ssrc = signals[i];
            if (timeDays > ssrc.expireDay) continue;
            const math::Vec3d posAU3 = ssrc.posKm / kAU_KM;
            const ImVec2 sp = worldToScreen({posAU3.x, posAU3.z});
            const core::u64 sSeed = core::hashCombine(core::hashCombine(core::fnv1a64("signal"), ssrc.id), (core::u64)(int)ssrc.type);
            const bool sel = (target.kind == TargetKind::Signal && target.index == i);
            drawIcon(sp, 13.0f * iconScale, render::SpriteKind::Signal, sSeed, IM_COL32(255, 210, 120, 230), sel);
            considerHover(TargetKind::Signal, i, sp, 13.0f * iconScale, render::SpriteKind::Signal, sSeed,
                          std::string(signalTypeName(ssrc.type)) + " Signal");
          }
        }

        if (map.showAsteroids) {
          for (std::size_t i = 0; i < asteroids.size(); ++i) {
            const auto& a = asteroids[i];
            if (a.remainingUnits <= 1e-3) continue;
            const math::Vec3d posAU3 = a.posKm / kAU_KM;
            const ImVec2 sp = worldToScreen({posAU3.x, posAU3.z});
            const core::u64 aSeed = core::hashCombine(core::hashCombine(core::fnv1a64("asteroid"), a.id), (core::u64)a.yield);
            const bool sel = (target.kind == TargetKind::Asteroid && target.index == i);
            drawIcon(sp, 13.0f * iconScale, render::SpriteKind::Asteroid, aSeed, IM_COL32(190, 190, 200, 220), sel);
            considerHover(TargetKind::Asteroid, i, sp, 13.0f * iconScale, render::SpriteKind::Asteroid, aSeed,
                          std::string("Asteroid [") + econ::commodityDef(a.yield).name + "]");
          }
        }

        if (map.showCargo) {
          for (std::size_t i = 0; i < floatingCargo.size(); ++i) {
            const auto& pod = floatingCargo[i];
            if (pod.units <= 0.0) continue;
            const math::Vec3d posAU3 = pod.posKm / kAU_KM;
            const ImVec2 sp = worldToScreen({posAU3.x, posAU3.z});
            const core::u64 cSeed = core::hashCombine(core::hashCombine(core::fnv1a64("cargo"), pod.id), (core::u64)pod.commodity);
            const bool sel = (target.kind == TargetKind::Cargo && target.index == i);
            drawIcon(sp, 13.0f * iconScale, render::SpriteKind::Cargo, cSeed, IM_COL32(255, 230, 140, 230), sel);
            considerHover(TargetKind::Cargo, i, sp, 13.0f * iconScale, render::SpriteKind::Cargo, cSeed,
                          std::string(econ::commodityDef(pod.commodity).name) + " Pod");
          }
        }

        // Show supercruise drop ring for the current target (helps the "7-second rule").
        if (target.kind != TargetKind::None && target.kind != TargetKind::Star) {
          double dropKm = 0.0;
          math::Vec3d destKm{0,0,0};
          bool ok = false;

          if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
            const auto& st = currentSystem->stations[target.index];
            destKm = stationPosKm(st, timeDays);
            dropKm = std::max(15000.0, st.radiusKm * 3.0);
            ok = true;
          } else if (target.kind == TargetKind::Planet && target.index < currentSystem->planets.size()) {
            const auto& p = currentSystem->planets[target.index];
            destKm = planetPosKm(p, timeDays);
            const double rKm = p.radiusEarth * 6371.0;
            dropKm = std::max(60000.0, rKm * 12.0);
            ok = true;
          } else if (target.kind == TargetKind::Signal && target.index < signals.size()) {
            destKm = signals[target.index].posKm;
            dropKm = 70000.0;
            ok = true;
          }

          if (ok && dropKm > 1.0) {
            const ImVec2 sp = worldToScreen({destKm.x / kAU_KM, destKm.z / kAU_KM});
            const float rPx = (float)((dropKm / kAU_KM) * (double)pxPerAU);
            if (rPx > 4.0f && rPx < 50000.0f) {
              draw->AddCircle(sp, rPx, IM_COL32(255, 210, 120, 80), 0, 1.5f);
            }
          }
        }

        // Tooltip
        if (canvasHovered && hover.has) {
          ImGui::BeginTooltip();
          ImGui::TextUnformatted(hover.label.c_str());

          // Distance from ship (now)
          math::Vec3d wPosKm{0,0,0};
          bool hasPos = false;
          if (hover.kind == TargetKind::Star) {
            wPosKm = {0,0,0};
            hasPos = true;
          } else if (hover.kind == TargetKind::Planet && hover.idx < currentSystem->planets.size()) {
            wPosKm = planetPosKm(currentSystem->planets[hover.idx], timeDays);
            hasPos = true;
          } else if (hover.kind == TargetKind::Station && hover.idx < currentSystem->stations.size()) {
            wPosKm = stationPosKm(currentSystem->stations[hover.idx], timeDays);
            hasPos = true;
          } else if (hover.kind == TargetKind::Contact && hover.idx < contacts.size()) {
            wPosKm = contacts[hover.idx].ship.positionKm();
            hasPos = true;
          } else if (hover.kind == TargetKind::Signal && hover.idx < signals.size()) {
            wPosKm = signals[hover.idx].posKm;
            hasPos = true;
          } else if (hover.kind == TargetKind::Asteroid && hover.idx < asteroids.size()) {
            wPosKm = asteroids[hover.idx].posKm;
            hasPos = true;
          } else if (hover.kind == TargetKind::Cargo && hover.idx < floatingCargo.size()) {
            wPosKm = floatingCargo[hover.idx].posKm;
            hasPos = true;
          }

          if (hasPos) {
            const double dKm = (wPosKm - ship.positionKm()).length();
            ImGui::TextDisabled("Distance: %.0f km (%.4f AU)", dKm, dKm / kAU_KM);
          }

          const auto& big = spriteCache.get(hover.sprite, hover.seed, 96);
          ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(96, 96));
          ImGui::TextDisabled("Preview: t %+.2f days", (float)map.previewOffsetDays);
          ImGui::EndTooltip();
        }

        // Click interaction (target / context menu)
        static TargetKind ctxKind = TargetKind::None;
        static std::size_t ctxIdx = 0;
        static std::string ctxLabel;
        static render::SpriteKind ctxSpriteKind = render::SpriteKind::Star;
        static core::u64 ctxSpriteSeed = 0;

        if (canvasHovered && hover.has) {
          if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            target.kind = hover.kind;
            target.index = hover.idx;
          }

          if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
            map.followShip = false;
            if (hover.kind == TargetKind::Star) map.centerAU = {0.0, 0.0};
            else if (hover.kind == TargetKind::Planet && hover.idx < currentSystem->planets.size()) {
              const auto& p = currentSystem->planets[hover.idx];
              const math::Vec3d posAU3 = sim::orbitPosition3DAU(p.orbit, previewTimeDays);
              map.centerAU = {posAU3.x, posAU3.z};
            } else if (hover.kind == TargetKind::Station && hover.idx < currentSystem->stations.size()) {
              const auto& st = currentSystem->stations[hover.idx];
              const math::Vec3d posAU3 = sim::orbitPosition3DAU(st.orbit, previewTimeDays);
              map.centerAU = {posAU3.x, posAU3.z};
            } else if (hover.kind == TargetKind::Contact && hover.idx < contacts.size()) {
              const math::Vec3d au = contacts[hover.idx].ship.positionKm() / kAU_KM;
              map.centerAU = {au.x, au.z};
            } else if (hover.kind == TargetKind::Signal && hover.idx < signals.size()) {
              const math::Vec3d au = signals[hover.idx].posKm / kAU_KM;
              map.centerAU = {au.x, au.z};
            } else if (hover.kind == TargetKind::Asteroid && hover.idx < asteroids.size()) {
              const math::Vec3d au = asteroids[hover.idx].posKm / kAU_KM;
              map.centerAU = {au.x, au.z};
            } else if (hover.kind == TargetKind::Cargo && hover.idx < floatingCargo.size()) {
              const math::Vec3d au = floatingCargo[hover.idx].posKm / kAU_KM;
              map.centerAU = {au.x, au.z};
            }
          }

          if (ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
            ctxKind = hover.kind;
            ctxIdx = hover.idx;
            ctxLabel = hover.label;
            ctxSpriteKind = hover.sprite;
            ctxSpriteSeed = hover.seed;
            ImGui::OpenPopup("sysmap2_ctx");
          }
        } else if (canvasHovered && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
          ctxKind = TargetKind::None;
          ctxIdx = 0;
          ctxLabel.clear();
          ImGui::OpenPopup("sysmap2_ctx");
        }

        if (ImGui::BeginPopup("sysmap2_ctx")) {
          if (!ctxLabel.empty()) {
            ImGui::TextUnformatted(ctxLabel.c_str());
            ImGui::Separator();

            const auto& big = spriteCache.get(ctxSpriteKind, ctxSpriteSeed, 64);
            ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(48, 48));
            ImGui::SameLine();
            ImGui::TextDisabled("Actions");
            ImGui::Separator();

            if (ImGui::MenuItem("Target")) {
              target.kind = ctxKind;
              target.index = ctxIdx;
            }
            const std::string scanLabel =
                std::string("Scan (") + game::chordLabel(controls.actions.scannerAction) + ")";
            if (ImGui::MenuItem(scanLabel.c_str())) {
              uiStartScan(ctxKind, ctxIdx);
            }

            // Supercruise only makes sense for station/planet/signal.
            const bool canSc = (ctxKind == TargetKind::Station || ctxKind == TargetKind::Planet || ctxKind == TargetKind::Signal);
            const std::string scLabel =
                std::string("Engage supercruise (") + game::chordLabel(controls.actions.supercruise) + ")";
            if (ImGui::MenuItem(scLabel.c_str(), nullptr, false, canSc)) {
              target.kind = ctxKind;
              target.index = ctxIdx;
              uiToggleSupercruise();
            }

            if (ctxKind == TargetKind::Station) {
              const std::string apLabel = std::string("Docking computer (auto-dock) (") +
                                         game::chordLabel(controls.actions.toggleAutopilot) + ")";
              if (ImGui::MenuItem(apLabel.c_str())) {
                target.kind = ctxKind;
                target.index = ctxIdx;
                autopilot = true;
                dockingComputer.reset();
              }
            }

            if (ImGui::MenuItem("Center view")) {
              map.followShip = false;
              if (ctxKind == TargetKind::Star) map.centerAU = {0.0, 0.0};
              else if (ctxKind == TargetKind::Planet && ctxIdx < currentSystem->planets.size()) {
                const auto& p = currentSystem->planets[ctxIdx];
                const math::Vec3d posAU3 = sim::orbitPosition3DAU(p.orbit, previewTimeDays);
                map.centerAU = {posAU3.x, posAU3.z};
              } else if (ctxKind == TargetKind::Station && ctxIdx < currentSystem->stations.size()) {
                const auto& st = currentSystem->stations[ctxIdx];
                const math::Vec3d posAU3 = sim::orbitPosition3DAU(st.orbit, previewTimeDays);
                map.centerAU = {posAU3.x, posAU3.z};
              } else if (ctxKind == TargetKind::Contact && ctxIdx < contacts.size()) {
                const math::Vec3d au = contacts[ctxIdx].ship.positionKm() / kAU_KM;
                map.centerAU = {au.x, au.z};
              } else if (ctxKind == TargetKind::Signal && ctxIdx < signals.size()) {
                const math::Vec3d au = signals[ctxIdx].posKm / kAU_KM;
                map.centerAU = {au.x, au.z};
              } else if (ctxKind == TargetKind::Asteroid && ctxIdx < asteroids.size()) {
                const math::Vec3d au = asteroids[ctxIdx].posKm / kAU_KM;
                map.centerAU = {au.x, au.z};
              } else if (ctxKind == TargetKind::Cargo && ctxIdx < floatingCargo.size()) {
                const math::Vec3d au = floatingCargo[ctxIdx].posKm / kAU_KM;
                map.centerAU = {au.x, au.z};
              }
            }
          } else {
            ImGui::TextDisabled("System Map");
            ImGui::Separator();
            if (ImGui::MenuItem("Center on star")) centerOnStar();
            if (ImGui::MenuItem("Center on ship")) { map.centerAU = shipPosAU2(); map.followShip = true; }
            ImGui::MenuItem("Follow ship", nullptr, &map.followShip);
            ImGui::MenuItem("Grid", nullptr, &map.showGrid);
            ImGui::MenuItem("Orbits", nullptr, &map.showOrbits);
            if (ImGui::MenuItem("Reset zoom")) map.zoom = 1.0f;
            if (ImGui::MenuItem("Reset preview time")) { map.previewOffsetDays = 0.0; map.previewAnimate = false; }
          }
          ImGui::EndPopup();
        }

        // Double-click empty canvas centers on star.
        if (canvasHovered && !hover.has && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
          centerOnStar();
        }

        ImGui::EndChild();

        // ---- Details column ----
        ImGui::TableNextColumn();
        ImGui::BeginChild("sysmap2_details", ImVec2(0, 0), true);

        ImGui::Text("Selection");
        ImGui::Separator();

        auto drawAtlasIconInline = [&](render::SpriteKind kind, core::u64 seed, float sizePx, ImVec4 tint) {
          const auto uv = hudAtlas.get(kind, seed);
          ImGui::Image(atlasId, ImVec2(sizePx, sizePx), ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1), tint);
        };

        auto showTargetHeader = [&](const std::string& title, render::SpriteKind sk, core::u64 sseed, ImVec4 tint) {
          drawAtlasIconInline(sk, sseed, 22.0f, tint);
          ImGui::SameLine();
          ImGui::TextUnformatted(title.c_str());
        };

        if (target.kind == TargetKind::None) {
          ImGui::TextDisabled("No target selected.");
        } else if (target.kind == TargetKind::Star) {
          const core::u64 sSeed = core::hashCombine(core::fnv1a64("star"), (core::u64)currentSystem->stub.seed);
          showTargetHeader(std::string("Star (") + starClassName(currentSystem->star.cls) + ")",
                           render::SpriteKind::Star, sSeed, ImVec4(1,1,1,1));
        } else if (target.kind == TargetKind::Planet && target.index < currentSystem->planets.size()) {
          const auto& p = currentSystem->planets[target.index];
          const core::u64 pSeed = core::hashCombine(core::hashCombine(core::fnv1a64("planet"), (core::u64)currentSystem->stub.seed), (core::u64)target.index);
          showTargetHeader(p.name, render::SpriteKind::Planet, pSeed, ImVec4(1,1,1,1));
          const math::Vec3d pPos = planetPosKm(p, timeDays);
          ImGui::TextDisabled("Distance: %.0f km", (pPos - ship.positionKm()).length());
        } else if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
          const auto& st = currentSystem->stations[target.index];
          const core::u64 stSeed = core::hashCombine(core::fnv1a64("station"), (core::u64)st.id);
          showTargetHeader(st.name, render::SpriteKind::Station, stSeed, ImVec4(1,1,1,1));
          const math::Vec3d stPos = stationPosKm(st, timeDays);
          ImGui::TextDisabled("Distance: %.0f km", (stPos - ship.positionKm()).length());
          ImGui::TextDisabled("Type: %s", stationTypeName(st.type));
        } else if (target.kind == TargetKind::Contact && target.index < contacts.size()) {
          const auto& ctc = contacts[target.index];
          const core::u64 cSeed = core::hashCombine(core::fnv1a64("ship"), ctc.id);
          showTargetHeader(ctc.name + " [" + contactRoleName(ctc.role) + "]", render::SpriteKind::Ship, cSeed, ImVec4(1,1,1,1));
          ImGui::TextDisabled("Distance: %.0f km", (ctc.ship.positionKm() - ship.positionKm()).length());
        } else if (target.kind == TargetKind::Signal && target.index < signals.size()) {
          const auto& ssrc = signals[target.index];
          const core::u64 sSeed = core::hashCombine(core::hashCombine(core::fnv1a64("signal"), ssrc.id), (core::u64)(int)ssrc.type);
          showTargetHeader(std::string(signalTypeName(ssrc.type)) + " Signal", render::SpriteKind::Signal, sSeed, ImVec4(1,1,1,1));
          ImGui::TextDisabled("Distance: %.0f km", (ssrc.posKm - ship.positionKm()).length());
        } else if (target.kind == TargetKind::Asteroid && target.index < asteroids.size()) {
          const auto& a = asteroids[target.index];
          const core::u64 aSeed = core::hashCombine(core::hashCombine(core::fnv1a64("asteroid"), a.id), (core::u64)a.yield);
          showTargetHeader("Asteroid", render::SpriteKind::Asteroid, aSeed, ImVec4(1,1,1,1));
          ImGui::TextDisabled("Yield: %s", econ::commodityDef(a.yield).name);
          ImGui::TextDisabled("Remaining: %.0f", a.remainingUnits);
          ImGui::TextDisabled("Distance: %.0f km", (a.posKm - ship.positionKm()).length());
        } else if (target.kind == TargetKind::Cargo && target.index < floatingCargo.size()) {
          const auto& pod = floatingCargo[target.index];
          const core::u64 cSeed = core::hashCombine(core::hashCombine(core::fnv1a64("cargo"), pod.id), (core::u64)pod.commodity);
          showTargetHeader("Cargo Pod", render::SpriteKind::Cargo, cSeed, ImVec4(1,1,1,1));
          ImGui::TextDisabled("Commodity: %s", econ::commodityDef(pod.commodity).name);
          ImGui::TextDisabled("Units: %.0f", pod.units);
          ImGui::TextDisabled("Distance: %.0f km", (pod.posKm - ship.positionKm()).length());
        } else {
          ImGui::TextDisabled("Invalid selection.");
        }

        ImGui::Separator();

        // Actions
        ImGui::BeginDisabled(target.kind == TargetKind::None);
        const std::string scanBtnLabel =
            std::string("Scan target (") + game::chordLabel(controls.actions.scannerAction) + ")";
        if (ImGui::Button(scanBtnLabel.c_str())) uiStartScan(target.kind, target.index);
        ImGui::EndDisabled();

        const bool canScTarget = (target.kind == TargetKind::Station || target.kind == TargetKind::Planet || target.kind == TargetKind::Signal);
        ImGui::BeginDisabled(!canScTarget);
        const std::string supercruiseBtnLabel =
            std::string("Supercruise (") + game::chordLabel(controls.actions.supercruise) + ")";
        if (ImGui::Button(supercruiseBtnLabel.c_str())) uiToggleSupercruise();
        ImGui::EndDisabled();

        if (target.kind == TargetKind::Station) {
          ImGui::BeginDisabled(docked);
          const std::string autopilotBtnLabel =
              std::string("Docking computer (auto-dock) (") + game::chordLabel(controls.actions.toggleAutopilot) + ")";
          if (ImGui::Button(autopilotBtnLabel.c_str())) {
            autopilot = true;
            dockingComputer.reset();
          }
          ImGui::EndDisabled();
        }

        if (ImGui::Button("Center map on selection")) {
          map.followShip = false;
          if (target.kind == TargetKind::Star) map.centerAU = {0.0, 0.0};
          else if (target.kind == TargetKind::Planet && target.index < currentSystem->planets.size()) {
            const auto& p = currentSystem->planets[target.index];
            const math::Vec3d posAU3 = sim::orbitPosition3DAU(p.orbit, previewTimeDays);
            map.centerAU = {posAU3.x, posAU3.z};
          } else if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
            const auto& st = currentSystem->stations[target.index];
            const math::Vec3d posAU3 = sim::orbitPosition3DAU(st.orbit, previewTimeDays);
            map.centerAU = {posAU3.x, posAU3.z};
          } else if (target.kind == TargetKind::Contact && target.index < contacts.size()) {
            const math::Vec3d au = contacts[target.index].ship.positionKm() / kAU_KM;
            map.centerAU = {au.x, au.z};
          } else if (target.kind == TargetKind::Signal && target.index < signals.size()) {
            const math::Vec3d au = signals[target.index].posKm / kAU_KM;
            map.centerAU = {au.x, au.z};
          } else if (target.kind == TargetKind::Asteroid && target.index < asteroids.size()) {
            const math::Vec3d au = asteroids[target.index].posKm / kAU_KM;
            map.centerAU = {au.x, au.z};
          } else if (target.kind == TargetKind::Cargo && target.index < floatingCargo.size()) {
            const math::Vec3d au = floatingCargo[target.index].posKm / kAU_KM;
            map.centerAU = {au.x, au.z};
          }
        }

        ImGui::Separator();
        ImGui::TextDisabled("Quick targets");

        if (ImGui::BeginTabBar("sysmap2_quicktabs", ImGuiTabBarFlags_NoCloseWithMiddleMouseButton)) {
          if (ImGui::BeginTabItem("Planets")) {
            for (std::size_t i = 0; i < currentSystem->planets.size(); ++i) {
              const auto& p = currentSystem->planets[i];
              const core::u64 pSeed = core::hashCombine(core::hashCombine(core::fnv1a64("planet"), (core::u64)currentSystem->stub.seed), (core::u64)i);
              drawAtlasIconInline(render::SpriteKind::Planet, pSeed, 18.0f, ImVec4(1,1,1,1));
              ImGui::SameLine();
              const bool sel = (target.kind == TargetKind::Planet && target.index == i);
              if (ImGui::Selectable(p.name.c_str(), sel)) {
                target.kind = TargetKind::Planet;
                target.index = i;
              }
            }
            ImGui::EndTabItem();
          }
          if (ImGui::BeginTabItem("Stations")) {
            for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
              const auto& st = currentSystem->stations[i];
              const core::u64 stSeed = core::hashCombine(core::fnv1a64("station"), (core::u64)st.id);
              drawAtlasIconInline(render::SpriteKind::Station, stSeed, 18.0f, ImVec4(1,1,1,1));
              ImGui::SameLine();
              const bool sel = (target.kind == TargetKind::Station && target.index == i);
              if (ImGui::Selectable(st.name.c_str(), sel)) {
                target.kind = TargetKind::Station;
                target.index = i;
              }
            }
            ImGui::EndTabItem();
          }
          if (ImGui::BeginTabItem("Signals")) {
            for (std::size_t i = 0; i < signals.size(); ++i) {
              const auto& ssrc = signals[i];
              if (timeDays > ssrc.expireDay) continue;
              const core::u64 sSeed = core::hashCombine(core::hashCombine(core::fnv1a64("signal"), ssrc.id), (core::u64)(int)ssrc.type);
              drawAtlasIconInline(render::SpriteKind::Signal, sSeed, 18.0f, ImVec4(1,1,1,1));
              ImGui::SameLine();
              const std::string name = std::string(signalTypeName(ssrc.type)) + " #" + std::to_string((int)i);
              const bool sel = (target.kind == TargetKind::Signal && target.index == i);
              if (ImGui::Selectable(name.c_str(), sel)) {
                target.kind = TargetKind::Signal;
                target.index = i;
              }
            }
            ImGui::EndTabItem();
          }
          ImGui::EndTabBar();
        }

        ImGui::EndChild();
        ImGui::EndTable();
      }

      ImGui::EndTabItem();
    }

    ImGui::EndTabBar();
  }

  ImGui::End();
}

    if (showEconomy) {
      beginStationSelectorHUD(*currentSystem, selectedStationIndex, docked, dockedStationId);

      if (!currentSystem->stations.empty()) {
        const auto& station = currentSystem->stations[(std::size_t)selectedStationIndex];
        auto& stEcon = universe.stationEconomy(station, timeDays);
        const double rep = getRep(station.factionId);
        const double feeEff = effectiveFeeRate(station);
        double cargoKgNow = cargoMassKg(cargo);

        ImGui::Begin("Market Details");
        ImGui::Text("Credits: %.2f", credits);

        ImGui::Text("Station: %s (%s)", station.name.c_str(), stationTypeName(station.type));
        ImGui::Text("Faction: %s   Rep: %.0f", factionName(station.factionId).c_str(), rep);
        ImGui::Text("Fees: base %.2f%%  effective %.2f%%", station.feeRate * 100.0, feeEff * 100.0);
        ImGui::Text("Cargo: %.0f / %.0f kg", cargoKgNow, cargoCapacityKg);

        const bool canTrade = docked && (station.id == dockedStationId);
        if (!canTrade) {
          ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.35f, 1.0f), "Trade disabled: dock at this station to buy/sell.");
        }

        // Simple docked services
        if (canTrade) {
          const double hullMissing = std::max(0.0, playerHullMax - playerHull);
          const double repairBase = hullMissing * 12.0;
          const double repairCost = repairBase * (1.0 + feeEff);

          if (ImGui::Button("Repair hull")) {
            if (hullMissing <= 0.01) {
              toast(toasts, "Hull already full.", 2.0);
            } else if (credits >= repairCost) {
              credits -= repairCost;
              playerHull = playerHullMax;
              toast(toasts, "Ship repaired.", 2.0);
            } else {
              toast(toasts, "Not enough credits for repair.", 2.0);
            }
          }
          ImGui::SameLine();
          ImGui::TextDisabled("(%.0f cr)", repairCost);

          // Refuel service: buys Fuel commodity into the ship tank.
          const auto fuelQuote = econ::quote(stEcon, station.economyModel, econ::CommodityId::Fuel, 0.10);
          const double fuelNeed = std::max(0.0, fuelMax - fuel);
          const double fuelAvail = std::max(0.0, fuelQuote.inventory);
          const double fuelBuy = std::min(fuelNeed, fuelAvail);
          const double fuelCost = fuelBuy * fuelQuote.ask * (1.0 + feeEff);

          if (ImGui::Button("Refuel")) {
            if (fuelNeed <= 0.01) {
              toast(toasts, "Fuel tank already full.", 2.0);
            } else if (fuelBuy <= 0.01) {
              toast(toasts, "Station is out of fuel.", 2.0);
            } else if (credits >= fuelCost) {
              // Take fuel from the station's commodity inventory (clamped). We charge for the
              // amount actually transferred (useful if inventories are changing due to NPC traffic).
              const double taken = econ::takeInventory(stEcon, station.economyModel, econ::CommodityId::Fuel, fuelBuy);
              if (taken <= 0.01) {
                toast(toasts, "Station is out of fuel.", 2.0);
              } else {
                const double cost = taken * fuelQuote.ask * (1.0 + feeEff);
                credits -= cost;
                fuel += taken;
                toast(toasts, "Refueled.", 2.0);
              }
            } else {
              toast(toasts, "Not enough credits to refuel.", 2.0);
            }
          }
          ImGui::SameLine();
          ImGui::TextDisabled("(%.1f units, %.0f cr)", fuelBuy, fuelCost);

          // Industry / fabrication orders
          {
            const auto recipes = sim::availableIndustryRecipes(station.type);
            if (!recipes.empty()) {
              ImGui::Separator();
              ImGui::Text("Industry (processing)");

              // Summary counts
              int activeHere = 0;
              int readyHere = 0;
              int activeElsewhere = 0;
              for (const auto& o : industryOrders) {
                if (o.claimed) continue;
                if (o.stationId == station.id) {
                  ++activeHere;
                  if (timeDays >= o.readyDay - 1e-6) ++readyHere;
                } else {
                  ++activeElsewhere;
                }
              }

              if (activeHere == 0) {
                ImGui::TextDisabled("No active orders at this station.");
              } else {
                ImGui::TextDisabled("Active here: %d  ready: %d", activeHere, readyHere);
              }
              if (activeElsewhere > 0) {
                ImGui::TextDisabled("Orders elsewhere: %d", activeElsewhere);
              }

              // Claim ready orders at this station
              for (auto& o : industryOrders) {
                if (o.claimed) continue;
                if (o.stationId != station.id) continue;

                const auto* def = sim::findIndustryRecipe(o.recipe);
                const char* code = def ? def->code : "ORDER";
                const bool ready = (timeDays >= o.readyDay - 1e-6);

                ImGui::PushID((int)o.id);
                ImGui::Bullet();
                ImGui::SameLine();
                if (o.inputBUnits > 0.0) {
                  ImGui::Text("%s #%llu: %s %.1f + %s %.1f -> %s %.1f",
                             code,
                             (unsigned long long)o.id,
                             econ::commodityCode(o.inputA), o.inputAUnits,
                             econ::commodityCode(o.inputB), o.inputBUnits,
                             econ::commodityCode(o.output), o.outputUnits);
                } else {
                  ImGui::Text("%s #%llu: %s %.1f -> %s %.1f",
                             code,
                             (unsigned long long)o.id,
                             econ::commodityCode(o.inputA), o.inputAUnits,
                             econ::commodityCode(o.output), o.outputUnits);
                }

                ImGui::SameLine();
                if (ready) {
                  ImGui::BeginDisabled(o.outputUnits <= 0.01);
                  if (ImGui::SmallButton("Claim")) {
                    const auto r = sim::claimIndustryOrderToCargo(o, station, timeDays, cargo, cargoCapacityKg);
                    if (!r.ok) {
                      const std::string why = r.reason ? r.reason : "";
                      if (why == "cargo_full") {
                        toast(toasts, "Not enough free cargo capacity.", 2.0);
                      } else {
                        toast(toasts, "Unable to claim order.", 2.0);
                      }
                    } else {
                      if (r.completed) {
                        toast(toasts, "Order claimed.", 2.0);
                      } else {
                        toast(toasts, "Claimed " + std::to_string((int)std::round(r.unitsMoved)) + " units (cargo full).", 2.0);
                      }
                    }
                  }

                  ImGui::SameLine();
                  if (ImGui::SmallButton("Store")) {
                    const auto r = sim::moveIndustryOrderOutputToWarehouse(o, station, timeDays, rep, stationStorage);
                    if (!r.ok) {
                      toast(toasts, "Unable to move output to warehouse.", 2.0);
                    } else {
                      toast(toasts, "Order output moved to warehouse storage.", 2.0);
                    }
                  }
                  ImGui::EndDisabled();
                } else {
                  const double hoursLeft = std::max(0.0, (o.readyDay - timeDays) * 24.0);
                  ImGui::TextDisabled("Ready in %.1f h", hoursLeft);
                }
                ImGui::PopID();
              }

              // Cleanup
              if (activeHere + activeElsewhere > 0) {
                if (ImGui::SmallButton("Clear completed")) {
                  sim::pruneClaimedIndustryOrders(industryOrders);
                  toast(toasts, "Cleared completed orders.", 1.6);
                }
              }

              // New order
              ImGui::Separator();
              ImGui::Text("New order");

              static int recipeChoice = 0;
              recipeChoice = std::clamp(recipeChoice, 0, (int)recipes.size() - 1);

              auto recipeLabel = [&](const sim::IndustryRecipeDef& r) {
                return std::string(r.code) + " - " + r.name;
              };

              const auto* r = recipes[recipeChoice];
              std::string preview = recipeLabel(*r);
              if (ImGui::BeginCombo("Recipe", preview.c_str())) {
                for (int i = 0; i < (int)recipes.size(); ++i) {
                  const bool selected = (i == recipeChoice);
                  std::string label = recipeLabel(*recipes[i]);
                  if (ImGui::Selectable(label.c_str(), selected)) {
                    recipeChoice = i;
                  }
                  if (selected) ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
              }

              r = recipes[recipeChoice];
              ImGui::TextDisabled("%s", r->desc);

              // Compute max batches from cargo.
              const double haveA = cargo[(int)r->inputA];
              const double maxA = (r->inputAUnits > 1e-9) ? (haveA / r->inputAUnits) : 0.0;
              double maxB = 1e9;
              if (r->inputBUnits > 1e-9) {
                const double haveB = cargo[(int)r->inputB];
                maxB = haveB / r->inputBUnits;
              }
              const int maxBatches = (int)std::floor(std::max(0.0, std::min(maxA, maxB)) + 1e-9);

              static int batchCount = 1;
              batchCount = std::clamp(batchCount, 1, std::max(1, maxBatches));
              ImGui::BeginDisabled(maxBatches <= 0);
              ImGui::SliderInt("Batches", &batchCount, 1, std::max(1, maxBatches));
              ImGui::EndDisabled();

              const double batches = (maxBatches <= 0) ? 0.0 : (double)batchCount;
              const auto q = sim::quoteIndustryOrder(*r, station.id, station.type, batches, feeEff, rep);

              ImGui::Text("Input: %s %.1f", econ::commodityCode(q.inputA), q.inputAUnits);
              if (q.inputBUnits > 0.0) {
                ImGui::Text("       %s %.1f", econ::commodityCode(q.inputB), q.inputBUnits);
              }
              ImGui::Text("Output: %s %.1f", econ::commodityCode(q.output), q.outputUnits);
              ImGui::TextDisabled("Time: %.1f h   Fee: %.0f cr", q.timeDays * 24.0, q.serviceFeeCr);
              ImGui::TextDisabled("Station mods: yield x%.3f  speed x%.3f  fee x%.3f",
                                  q.mods.yieldMul, q.mods.speedMul, q.mods.feeMul);

              const bool hasInputs =
                (batches > 0.0) &&
                (cargo[(int)q.inputA] + 1e-6 >= q.inputAUnits) &&
                ((q.inputBUnits <= 0.0) || (cargo[(int)q.inputB] + 1e-6 >= q.inputBUnits));
              const bool canAfford = (credits + 1e-6 >= q.serviceFeeCr);
              ImGui::BeginDisabled(!(hasInputs && canAfford));
              if (ImGui::Button("Submit order")) {
                const auto res = sim::submitIndustryOrder(nextIndustryOrderId,
                                                          industryOrders,
                                                          cargo,
                                                          credits,
                                                          station,
                                                          timeDays,
                                                          rep,
                                                          *r,
                                                          batchCount,
                                                          feeEff);
                if (res.ok) {
                  toast(toasts, "Industry order submitted.", 2.0);
                } else {
                  toast(toasts, "Unable to submit industry order.", 2.0);
                }
              }
              ImGui::EndDisabled();

              if (!hasInputs) {
                ImGui::SameLine();
                ImGui::TextDisabled("(need inputs)");
              } else if (!canAfford) {
                ImGui::SameLine();
                ImGui::TextDisabled("(need credits)");
              }
            }
          }

          // Shipyard upgrades
          if (station.type == econ::StationType::Shipyard) {
            ImGui::Separator();
            ImGui::Text("Shipyard Upgrades");

            const double cargoUpgradeKg = 200.0;
            const double cargoUpgradeCost = 2000.0 * (1.0 + feeEff);
            if (ImGui::Button("Cargo racks +200kg")) {
              if (credits >= cargoUpgradeCost) {
                credits -= cargoUpgradeCost;
                cargoCapacityKg += cargoUpgradeKg;
                toast(toasts, "Cargo capacity upgraded.", 2.0);
              } else {
                toast(toasts, "Not enough credits.", 2.0);
              }
            }
            ImGui::SameLine();
            ImGui::TextDisabled("(%.0f cr)", cargoUpgradeCost);

            const int cabinUpgradeSeats = 2;
            const double cabinUpgradeCost = 3000.0 * (1.0 + feeEff);
            if (ImGui::Button("Passenger cabins +2 seats")) {
              if (credits >= cabinUpgradeCost) {
                credits -= cabinUpgradeCost;
                passengerSeats += cabinUpgradeSeats;
                toast(toasts, "Passenger cabins upgraded.", 2.0);
              } else {
                toast(toasts, "Not enough credits.", 2.0);
              }
            }
            ImGui::SameLine();
            ImGui::TextDisabled("(%.0f cr)", cabinUpgradeCost);

            // Smuggling / stealth upgrade (reduces cargo scan chance when carrying contraband).
            {
              static const double kSmuggleUpgradeCost[4] = {0.0, 4500.0, 8500.0, 14000.0}; // Mk1..Mk3
              if (smuggleHoldMk >= 3) {
                ImGui::TextDisabled("Smuggling compartments: Mk3 (max)");
              } else {
                const int nextMk = smuggleHoldMk + 1;
                const double cost = kSmuggleUpgradeCost[nextMk] * (1.0 + feeEff);
                std::string label = std::string("Smuggling compartments Mk") + std::to_string(nextMk);

                ImGui::BeginDisabled(credits + 1e-6 < cost);
                if (ImGui::Button(label.c_str())) {
                  credits -= cost;
                  smuggleHoldMk = nextMk;
                  toast(toasts, "Hidden compartments installed.", 2.0);
                }
                ImGui::EndDisabled();

                ImGui::SameLine();
                ImGui::TextDisabled("(%.0f cr) Reduce cargo scan chance for contraband.", cost);
              }
            }

            const double fuelUpgrade = 10.0;
            const double fuelUpgradeCost = 2500.0 * (1.0 + feeEff);
            if (ImGui::Button("Fuel tank +10")) {
              if (credits >= fuelUpgradeCost) {
                credits -= fuelUpgradeCost;
                fuelMax += fuelUpgrade;
                toast(toasts, "Fuel tank upgraded.", 2.0);
              } else {
                toast(toasts, "Not enough credits.", 2.0);
              }
            }
            ImGui::SameLine();
            ImGui::TextDisabled("(%.0f cr)", fuelUpgradeCost);

            const double rangeUpgrade = 2.0;
            const double rangeUpgradeCost = 9000.0 * (1.0 + feeEff);
            if (ImGui::Button("FSD tuning +2ly")) {
              if (credits >= rangeUpgradeCost) {
                credits -= rangeUpgradeCost;
                fsdRangeLy += rangeUpgrade;
                toast(toasts, "FSD range improved.", 2.0);
              } else {
                toast(toasts, "Not enough credits.", 2.0);
              }
            }
            ImGui::SameLine();
            ImGui::TextDisabled("(%.0f cr)", rangeUpgradeCost);

            ImGui::Separator();
            ImGui::Text("Core Loadout");

            // Hull selection
            {
              static int hullChoice = 0;
              static int hullChoiceSync = -1;
              if (hullChoiceSync != (int)shipHullClass) {
                hullChoice = (int)shipHullClass;
                hullChoiceSync = (int)shipHullClass;
              }

              if (ImGui::BeginCombo("Hull class", kHullDefs[hullChoice].name)) {
                for (int i = 0; i < (int)(sizeof(kHullDefs) / sizeof(kHullDefs[0])); ++i) {
                  const bool selected = (hullChoice == i);
                  if (ImGui::Selectable(kHullDefs[i].name, selected)) {
                    hullChoice = i;
                  }
                  if (selected) ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
              }

              if (hullChoice != (int)shipHullClass) {
                const auto& curHull = kHullDefs[(int)shipHullClass];
                const auto& newHull = kHullDefs[hullChoice];
                const double newCargoCapKg = cargoCapacityKg * (newHull.cargoMult / std::max(1e-6, curHull.cargoMult));
                const double newFuelMax = fuelMax * (newHull.fuelMult / std::max(1e-6, curHull.fuelMult));
                const double hullCost = newHull.priceCr * (1.0 + feeEff);

                ImGui::TextDisabled("New caps: Cargo %.0f kg, Fuel %.1f", newCargoCapKg, newFuelMax);
                const bool cargoOk = (cargoMassKg(cargo) <= newCargoCapKg + 1e-6);
                if (!cargoOk) {
                  ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.6f, 1.0f), "Sell cargo first (new hull has less capacity).");
                }

                ImGui::BeginDisabled(!cargoOk);
                if (ImGui::Button("Buy selected hull")) {
                  if (credits >= hullCost) {
                    const double hullFrac = std::clamp(playerHull / std::max(1e-6, playerHullMax), 0.0, 1.0);
                    const double shieldFrac = std::clamp(playerShield / std::max(1e-6, playerShieldMax), 0.0, 1.0);
                    credits -= hullCost;

                    cargoCapacityKg = newCargoCapKg;
                    fuelMax = newFuelMax;
                    fuel = std::min(fuel, fuelMax);

                    shipHullClass = (ShipHullClass)hullChoice;
                    recalcPlayerStats();

                    playerHull = hullFrac * playerHullMax;
                    playerShield = shieldFrac * playerShieldMax;

                    toast(toasts, "Hull purchased: " + std::string(kHullDefs[hullChoice].name), 2.2);
                  } else {
                    toast(toasts, "Not enough credits.", 2.0);
                  }
                }
                ImGui::SameLine();
                ImGui::TextDisabled("(%.0f cr)", hullCost);
                ImGui::EndDisabled();
              }
            }

            // Core modules (thrusters / shield gen / distributor)
            {
              auto installMk = [&](const char* label, int& mkVar, const MkDef (&defs)[4]) {
                static int choice[3] = {1, 1, 1};
                static int sync[3] = {-1, -1, -1};
                int slot = 0;
                if (std::string(label).find("Thrusters") != std::string::npos) slot = 0;
                else if (std::string(label).find("Shield") != std::string::npos) slot = 1;
                else slot = 2;

                mkVar = std::clamp(mkVar, 1, 3);

                if (sync[slot] != mkVar) {
                  choice[slot] = mkVar;
                  sync[slot] = mkVar;
                }

                std::string comboId = std::string(label) + "##mk";
                if (ImGui::BeginCombo(comboId.c_str(), defs[choice[slot]].name)) {
                  for (int mk = 1; mk <= 3; ++mk) {
                    const bool selected = (choice[slot] == mk);
                    if (ImGui::Selectable(defs[mk].name, selected)) {
                      choice[slot] = mk;
                    }
                    if (selected) ImGui::SetItemDefaultFocus();
                  }
                  ImGui::EndCombo();
                }

                const int newMk = choice[slot];
                const double deltaCost = std::max(0.0, defs[newMk].priceCr - defs[mkVar].priceCr);
                const double cost = deltaCost * (1.0 + feeEff);

                ImGui::SameLine();
                ImGui::BeginDisabled(newMk == mkVar);
                if (ImGui::SmallButton((std::string("Install##") + label).c_str())) {
                  if (credits >= cost) {
                    const double hullFrac = std::clamp(playerHull / std::max(1e-6, playerHullMax), 0.0, 1.0);
                    const double shieldFrac = std::clamp(playerShield / std::max(1e-6, playerShieldMax), 0.0, 1.0);
                    credits -= cost;
                    mkVar = newMk;
                    recalcPlayerStats();
                    playerHull = hullFrac * playerHullMax;
                    playerShield = shieldFrac * playerShieldMax;
                    toast(toasts, std::string(label) + " installed.", 1.8);
                  } else {
                    toast(toasts, "Not enough credits.", 2.0);
                  }
                }
                ImGui::EndDisabled();
                ImGui::SameLine();
                ImGui::TextDisabled("(%.0f cr)", cost);
              };

              installMk("Thrusters", thrusterMk, kThrusters);
              installMk("Shield Gen", shieldMk, kShields);
              installMk("Distributor", distributorMk, kDistributors);
            }

            // Weapons / hardpoints
            {
              ImGui::Separator();
              ImGui::Text("Hardpoints");

              auto weaponCombo = [&](const char* label, WeaponType& wVar) {
                static int choiceP = 0;
                static int choiceS = 2;
                static int syncP = -1;
                static int syncS = -1;

                const bool isPrimary = (std::string(label).find("Primary") != std::string::npos);
                int& choice = isPrimary ? choiceP : choiceS;
                int& sync = isPrimary ? syncP : syncS;
                const int cur = (int)wVar;

                if (sync != cur) {
                  choice = cur;
                  sync = cur;
                }

                if (ImGui::BeginCombo(label, weaponDef((WeaponType)choice).name)) {
                  for (int i = 0; i < (int)(sizeof(kWeaponDefs) / sizeof(kWeaponDefs[0])); ++i) {
                    const bool selected = (choice == i);
                    if (ImGui::Selectable(kWeaponDefs[i].name, selected)) {
                      choice = i;
                    }
                    if (selected) ImGui::SetItemDefaultFocus();
                  }
                  ImGui::EndCombo();
                }

                const WeaponType newW = (WeaponType)choice;
                const double cost = (newW == wVar) ? 0.0 : (weaponDef(newW).priceCr * (1.0 + feeEff));
                ImGui::SameLine();
                ImGui::BeginDisabled(newW == wVar);
                if (ImGui::SmallButton((std::string("Buy##") + label).c_str())) {
                  if (credits >= cost) {
                    credits -= cost;
                    wVar = newW;
                    toast(toasts, "Weapon installed: " + std::string(weaponDef(newW).name), 1.8);
                  } else {
                    toast(toasts, "Not enough credits.", 2.0);
                  }
                }
                ImGui::EndDisabled();
                ImGui::SameLine();
                ImGui::TextDisabled("(%.0f cr)", cost);
              };

              weaponCombo("Primary (LMB)", weaponPrimary);
              weaponCombo("Secondary (RMB)", weaponSecondary);
            }
          }
        }

// Exploration / legal services
if (canTrade) {
  ImGui::Separator();
  ImGui::Text("Exploration Data");
  ImGui::Text("Bank: %.0f cr", explorationDataCr);
  if (explorationDataCr > 0.0) {
    if (ImGui::Button("Sell all exploration data")) {
      credits += explorationDataCr;
      toast(toasts, "Sold exploration data +" + std::to_string((int)explorationDataCr) + " cr.", 2.5);
      explorationDataCr = 0.0;
      addRep(station.factionId, +0.5);
    }
  } else {
    ImGui::TextDisabled("No scan data to sell.");
  }

  ImGui::Separator();
  ImGui::Text("Legal");
  const double bounty = getBounty(station.factionId);
  ImGui::Text("Outstanding bounty (this faction): %.0f cr", bounty);
  if (bounty > 0.0) {
    if (ImGui::Button("Pay bounty")) {
      if (credits >= bounty) {
        credits -= bounty;
        clearBounty(station.factionId);
        toast(toasts, "Bounty cleared.", 2.0);
      } else {
        toast(toasts, "Not enough credits to pay bounty.", 2.0);
      }
    }
  } else {
    ImGui::TextDisabled("No bounty.");
  }

	  const double vouchers = getVoucher(station.factionId);
	  ImGui::Text("Bounty vouchers (this faction): %.0f cr", vouchers);
	  if (vouchers > 0.0) {
	    if (ImGui::Button("Redeem vouchers")) {
	      credits += vouchers;
	      clearVoucher(station.factionId);
	      toast(toasts, "Redeemed bounty vouchers +" + std::to_string((int)vouchers) + " cr.", 2.5);
	      addRep(station.factionId, +0.6);
	    }
	  } else {
	    ImGui::TextDisabled("No vouchers.");
	  }
	  // Small hint if you have vouchers elsewhere
	  {
	    double other = 0.0;
	    for (const auto& kv : bountyVoucherByFaction) {
	      if (kv.first == station.factionId) continue;
	      other += kv.second;
	    }
	    if (other > 0.0) {
	      ImGui::TextDisabled("Other jurisdictions: %.0f cr", other);
	    }
	  }

  ImGui::Separator();
  ImGui::Text("Security Services");
  {
    const bool isAlert = (policeHeat > 0.15) || (policeAlertUntilDays > timeDays);
    const double fee = 150.0 + 250.0 * std::max(0.0, policeHeat);
    if (isAlert) {
      ImGui::Text("Local alert level: %.2f", policeHeat);
      if (ImGui::Button("Bribe port authority (reduce alert)")) {
        if (credits >= fee) {
          credits -= fee;
          policeHeat = std::max(0.0, policeHeat - 2.0);
          policeAlertUntilDays = std::min(policeAlertUntilDays, timeDays);
          toast(toasts, "Port authority bribed. Local alert reduced.", 2.2);
        } else {
          toast(toasts, "Not enough credits to bribe port authority.", 2.2);
        }
      }
      ImGui::SameLine();
      ImGui::TextDisabled("(%.0f cr)", fee);
      ImGui::TextDisabled("Reduces patrol response; does not clear bounties.");
    } else {
      ImGui::TextDisabled("No active alert.");
    }
  
  ImGui::Separator();
  ImGui::Text("Lay Low");
  {
    // Pay to spend time docked and reduce local alert. This trades off against mission deadlines.
    const bool canLayLow = (policeHeat > 0.05) || ((policeAlertUntilDays > timeDays) && (policeAlertUntilDays - timeDays) > 1e-6);
    const double hours = 3.0;
    const double layLowDays = hours / 24.0;

    if (canLayLow) {
      // Fee scales with current alert.
      const double fee = std::clamp(180.0 + policeHeat * 900.0 + (policeAlertUntilDays > timeDays ? 220.0 : 0.0), 120.0, 2500.0);

      if (ImGui::Button("Lay low (3h)")) {
        if (credits >= fee) {
          credits -= fee;

          // Advance time and cool alert. Decay is applied using the same curve as passive decay.
          timeDays += layLowDays;
          {
            const double tauSec = 22.0 * 60.0;
            const double k = std::exp(- (layLowDays * 86400.0) / std::max(1.0, tauSec));
            policeHeat *= k;
            if (policeHeat < 0.0005) policeHeat = 0.0;
          }

          // Paying to lay low also shortens the "active alert" window.
          if (policeAlertUntilDays > timeDays) {
            policeAlertUntilDays = std::max(timeDays, policeAlertUntilDays - layLowDays * 1.25);
          }

          toast(toasts, "You lay low for 3 hours. Local alert reduced.", 2.5);
        } else {
          toast(toasts, "Not enough credits to lay low.", 2.2);
        }
      }
      ImGui::SameLine();
      ImGui::TextDisabled("(%.0f cr)", fee);
      ImGui::TextDisabled("Advances time; reduces local alert/response.");
    } else {
      ImGui::TextDisabled("No need to lay low.");
    }
  }

}

}

        // Prevent accidentally selling mission-critical cargo. Delivery / multi-hop / smuggle
        // contracts often require specific commodities; reserve those units so the Sell buttons
        // only operate on your "free" cargo.
        std::array<double, econ::kCommodityCount> reservedMissionUnits{};
        reservedMissionUnits.fill(0.0);
        for (const auto& m : missions) {
          if (m.completed || m.failed) continue;
          if (m.type != sim::MissionType::Delivery
           && m.type != sim::MissionType::MultiDelivery
           && m.type != sim::MissionType::Smuggle) continue;
          const std::size_t idx = (std::size_t)m.commodity;
          if (idx >= reservedMissionUnits.size()) continue;
          reservedMissionUnits[idx] += std::max(0.0, m.units);
        }

        // Warehouse / station storage
        {
          ImGui::Separator();
          ImGui::Text("Warehouse (station storage)");

          // Accrue fees for this station's entry (if it exists) before displaying numbers.
          if (auto* e = sim::findStorage(stationStorage, station.id)) {
            sim::accrueStorageFees(*e, timeDays, rep);
          }

          const auto* entry = sim::findStorage(stationStorage, station.id);
          const double storedKg = entry ? sim::storageMassKg(*entry) : 0.0;
          const double feesDue = entry ? std::max(0.0, entry->feesDueCr) : 0.0;
          const double dailyFee = entry ? sim::estimateStorageDailyFeeCr(*entry, rep) : 0.0;

          if (!entry || storedKg <= 1e-6) {
            ImGui::TextDisabled("No cargo stored at this station.");
          } else {
            ImGui::Text("Stored mass: %.0f kg", storedKg);
          }
          ImGui::TextDisabled("Fees due: %.0f cr  (est. %.0f cr/day)", feesDue, dailyFee);

          if (entry && feesDue > 1e-6) {
            const bool canPay = (credits + 1e-6 >= feesDue);
            ImGui::BeginDisabled(!canPay);
            if (ImGui::SmallButton("Pay all fees")) {
              const double before = credits;
              const auto r = sim::payStorageFees(stationStorage, station, timeDays, rep, credits, 1e30);
              const double paid = std::max(0.0, before - credits);
              if (paid > 0.01) {
                toast(toasts, "Paid " + std::to_string((int)std::round(paid)) + " cr in storage fees.", 2.0);
              } else {
                toast(toasts, "No storage fees due.", 1.5);
              }
              sim::pruneEmptyStorage(stationStorage);
            }
            ImGui::EndDisabled();
            if (!canPay && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
              ImGui::SetTooltip("Not enough credits to pay fees (need %.0f cr).", feesDue);
            }
          }

          // Per-commodity deposit / withdraw table.
          if (ImGui::BeginTable("warehouse", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
            ImGui::TableSetupColumn("Commodity");
            ImGui::TableSetupColumn("Stored");
            ImGui::TableSetupColumn("Ship");
            ImGui::TableSetupColumn("Qty");
            ImGui::TableSetupColumn("Actions");
            ImGui::TableHeadersRow();

            static float whQty[(int)econ::kCommodityCount] = {};

            for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
              const auto cid = (econ::CommodityId)i;
              const auto* eNow = sim::findStorage(stationStorage, station.id);
              const double storedUnits = eNow ? eNow->cargo[i] : 0.0;
              const double dueNow = eNow ? std::max(0.0, eNow->feesDueCr) : 0.0;

              if (storedUnits <= 1e-6 && cargo[i] <= 1e-6) continue;

              ImGui::TableNextRow();
              ImGui::PushID((int)i);

              ImGui::TableSetColumnIndex(0);
              ImGui::Text("%s", std::string(econ::commodityName(cid)).c_str());

              ImGui::TableSetColumnIndex(1);
              ImGui::Text("%.0f", storedUnits);

              ImGui::TableSetColumnIndex(2);
              {
                const double reserved = reservedMissionUnits[i];
                if (reserved > 1e-6) {
                  const double keep = std::min(reserved, cargo[i]);
                  ImGui::Text("%.0f (%.0f res)", cargo[i], keep);
                } else {
                  ImGui::Text("%.0f", cargo[i]);
                }
              }

              ImGui::TableSetColumnIndex(3);
              if (whQty[i] <= 0.0f) whQty[i] = 10.0f;
              ImGui::SetNextItemWidth(70);
              ImGui::InputFloat("##wqty", &whQty[i], 1.0f, 10.0f, "%.0f");

              ImGui::TableSetColumnIndex(4);

              const double reqUnits = std::max(0.0, (double)whQty[i]);
              const double reserved = reservedMissionUnits[i];
              const double freeForStorage = std::max(0.0, cargo[i] - reserved);

              // Store (deposit)
              ImGui::BeginDisabled(freeForStorage <= 1e-6 || reqUnits <= 0.0);
              if (ImGui::SmallButton("Store")) {
                const double want = std::min(reqUnits, freeForStorage);
                const double beforeUnits = cargo[i];
                const auto r = sim::depositToStorage(stationStorage, station, timeDays, rep, cargo, cid, want);
                if (r.ok) {
                  const double moved = std::max(0.0, beforeUnits - cargo[i]);
                  cargoKgNow = std::max(0.0, cargoKgNow - moved * econ::commodityDef(cid).massKg);
                  toast(toasts, "Stored " + std::to_string((int)std::round(moved)) + " units.", 1.7);
                } else {
                  toast(toasts, std::string("Storage failed: ") + (r.reason ? r.reason : "error"), 2.0);
                }
              }
              ImGui::EndDisabled();
              if (freeForStorage <= 1e-6 && reserved > 1e-6 && cargo[i] > 1e-6
                  && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
                ImGui::SetTooltip("Reserved for active mission(s): %.0f units", std::min(reserved, cargo[i]));
              }

              ImGui::SameLine();

              // Withdraw (retrieve)
              const bool withdrawBlockedByFees = (dueNow > 1e-6) && (credits + 1e-6 < dueNow);
              const bool withdrawDisabled = (storedUnits <= 1e-6) || (reqUnits <= 0.0) || withdrawBlockedByFees;
              ImGui::BeginDisabled(withdrawDisabled);
              if (ImGui::SmallButton("Withdraw")) {
                const auto* eBefore = sim::findStorage(stationStorage, station.id);
                const double dueBefore = eBefore ? std::max(0.0, eBefore->feesDueCr) : 0.0;
                const double beforeUnits = cargo[i];
                const auto r = sim::withdrawFromStorage(stationStorage, station, timeDays, rep, cargo, credits, cargoCapacityKg, cid, reqUnits);
                if (r.ok) {
                  const double moved = std::max(0.0, cargo[i] - beforeUnits);
                  cargoKgNow += moved * econ::commodityDef(cid).massKg;
                  std::string msg = "Withdrew " + std::to_string((int)std::round(moved)) + " units.";
                  if (dueBefore > 1e-6) {
                    msg += " (paid fees)";
                  }
                  toast(toasts, msg, 1.9);
                  sim::pruneEmptyStorage(stationStorage);
                } else {
                  const char* why = r.reason ? r.reason : "error";
                  if (std::string(why) == "fees_due") {
                    toast(toasts, "Fees due: pay warehouse fees before withdrawing.", 2.2);
                  } else if (std::string(why) == "cargo_full") {
                    toast(toasts, "Cargo hold full (mass limit).", 2.2);
                  } else {
                    toast(toasts, std::string("Withdraw failed: ") + why, 2.0);
                  }
                }
              }
              ImGui::EndDisabled();
              if (withdrawBlockedByFees && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
                ImGui::SetTooltip("Need %.0f cr to pay storage fees before withdrawing.", dueNow);
              }

              ImGui::PopID();
            }

            ImGui::EndTable();
          }
        }

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
            const bool illegalHere = isIllegalCommodity(station.factionId, cid);
            const auto q = econ::quote(stEcon, station.economyModel, cid, 0.10);

            ImGui::TableNextRow();

            ImGui::TableSetColumnIndex(0);
            {
              const core::u64 iconSeed = core::hashCombine(core::fnv1a64("commodity"), (core::u64)i);
              const auto& iconTex = spriteCache.get(render::SpriteKind::Commodity, iconSeed, 32);
              ImGui::Image((ImTextureID)(intptr_t)iconTex.handle(), ImVec2(18, 18));
              if (ImGui::IsItemHovered()) {
                const auto& big = spriteCache.get(render::SpriteKind::Commodity, iconSeed, 96);
                ImGui::BeginTooltip();
                ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(96, 96));
                ImGui::EndTooltip();
              }
              ImGui::SameLine();
              ImGui::Text("%s%s", std::string(econ::commodityName(cid)).c_str(), illegalHere ? " (ILLEGAL)" : "");
            }

            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%.0f", q.inventory);

            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%.2f", q.ask);

            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%.2f", q.bid);

            ImGui::TableSetColumnIndex(4);
	            {
	              const double reserved = reservedMissionUnits[i];
	              if (reserved > 1e-6) {
	                const double keep = std::min(reserved, cargo[i]);
	                ImGui::Text("%.0f (%.0f res)", cargo[i], keep);
	              } else {
	                ImGui::Text("%.0f", cargo[i]);
	              }
	            }

            ImGui::TableSetColumnIndex(5);
            ImGui::PushID((int)i);

            static float qty[ (int)econ::kCommodityCount ] = {};
            if (qty[i] <= 0.0f) qty[i] = 10.0f;
            ImGui::SetNextItemWidth(70);
            ImGui::InputFloat("##qty", &qty[i], 1.0f, 10.0f, "%.0f");

            ImGui::SameLine();

            // Buying illegal goods isn't available on the open market; you can still *sell* them (black market).
            ImGui::BeginDisabled(!canTrade || illegalHere);
            if (ImGui::SmallButton("Buy")) {
              const double buyUnits = std::max(0.0, (double)qty[i]);
              const double addKg = buyUnits * econ::commodityDef(cid).massKg;
              if (cargoKgNow + addKg > cargoCapacityKg + 1e-6) {
                toast(toasts, "Cargo hold full (mass limit).", 2.0);
              } else {
                auto tr = econ::buy(stEcon, station.economyModel, cid, buyUnits, credits, 0.10, feeEff);
                if (tr.ok) {
                  cargo[i] += buyUnits;
                  cargoKgNow += addKg;
                }
              }
            }
            ImGui::EndDisabled();

            ImGui::SameLine();
	            const double reserved = reservedMissionUnits[i];
	            const double freeForSale = std::max(0.0, cargo[i] - reserved);
	            const bool sellDisabled = (!canTrade) || (freeForSale <= 1e-6);
	            ImGui::BeginDisabled(sellDisabled);
            const char* sellLabel = illegalHere ? "Sell (BM)" : "Sell";
            if (ImGui::SmallButton(sellLabel)) {
	              const double sellUnitsReq = std::min<double>(qty[i], freeForSale);
              if (sellUnitsReq > 0.0) {
                if (illegalHere) {
                  // Black market: pay a premium, but storage capacity still applies.
                  constexpr double kBlackMarketMult = 1.35;
                  const double moved = econ::addInventory(stEcon, station.economyModel, cid, sellUnitsReq);
                  if (moved > 0.0) {
                    credits += q.bid * moved * kBlackMarketMult;
                    cargo[i] -= moved;
                    cargoKgNow = std::max(0.0, cargoKgNow - moved * econ::commodityDef(cid).massKg);

                    // Dealing contraband attracts attention over time.
                    policeHeat = std::clamp(policeHeat + std::min(0.75, 0.10 + moved * 0.01), 0.0, 6.0);
                    addRep(station.factionId, -std::min(1.2, moved * 0.02));
                  } else {
                    toast(toasts, "No space to sell (storage full).", 2.0);
                  }
                } else {
                  auto tr = econ::sell(stEcon, station.economyModel, cid, sellUnitsReq, credits, 0.10, feeEff);
                  if (tr.ok) {
                    cargo[i] -= sellUnitsReq;
                    cargoKgNow = std::max(0.0, cargoKgNow - sellUnitsReq * econ::commodityDef(cid).massKg);
                  }
                }
              }
            }
            ImGui::EndDisabled();
	            if (sellDisabled && reserved > 1e-6 && cargo[i] > 1e-6
	                && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
	              ImGui::SetTooltip("Reserved for active mission(s): %.0f units", std::min(reserved, cargo[i]));
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

        ImGui::End();
      }
    }

    if (showTrade) {
      ImGui::Begin("Trade Helper");

      if (!currentSystem) {
        ImGui::TextDisabled("No current system.");
      } else {
        // Choose a 'from' station: prefer docked station; otherwise use the currently targeted station.
        sim::StationId fromStationId = 0;
        if (docked && dockedStationId != 0) {
          fromStationId = dockedStationId;
        } else if (target.kind == TargetKind::Station && target.index >= 0 && (std::size_t)target.index < currentSystem->stations.size()) {
          fromStationId = currentSystem->stations[(std::size_t)target.index].id;
        }

        if (fromStationId == 0) {
          ImGui::TextDisabled("Dock at (or target) a station to see suggested trade routes.");
        } else {
          const sim::Station* fromSt = nullptr;
          for (const auto& st : currentSystem->stations) {
            if (st.id == fromStationId) {
              fromSt = &st;
              break;
            }
          }

          if (!fromSt) {
            ImGui::TextDisabled("Station not found in current system.");
          } else {
            const bool fromDocked = docked && (dockedStationId == fromStationId);
            const double feeEff = effectiveFeeRate(*fromSt);
            ImGui::Text("From: %s%s", fromSt->name.c_str(), fromDocked ? " (docked)" : "");
            ImGui::TextDisabled("Fees (rep-adjusted): %.1f%%", feeEff * 100.0);

            float radiusLy = (float)tradeSearchRadiusLy;
            if (ImGui::SliderFloat("Search radius (ly)", &radiusLy, 80.0f, 500.0f, "%.0f")) {
              tradeSearchRadiusLy = (double)radiusLy;
              tradeIdeasDayStamp = -1; // force refresh
              tradeMixDayStamp = -1;
            }
            ImGui::SameLine();
            if (ImGui::SmallButton("Refresh")) {
              tradeIdeasDayStamp = -1;
              tradeMixDayStamp = -1;
            }

            {
              bool useFree = tradeUseFreeHold;
              if (ImGui::Checkbox("Use free hold", &useFree)) {
                tradeUseFreeHold = useFree;
                tradeIdeasDayStamp = -1;
                tradeMixDayStamp = -1;
              }
              ImGui::SameLine();
              bool includeSame = tradeIncludeSameSystem;
              if (ImGui::Checkbox("Include same system", &includeSame)) {
                tradeIncludeSameSystem = includeSame;
                tradeIdeasDayStamp = -1;
                tradeMixDayStamp = -1;
              }
              ImGui::SameLine();
              bool useManifest = tradeUseManifest;
              if (ImGui::Checkbox("Cargo mix", &useManifest)) {
                tradeUseManifest = useManifest;
                tradeIdeasDayStamp = -1;
                tradeMixDayStamp = -1;
              }

              float minProfit = (float)tradeMinNetProfit;
              if (ImGui::SliderFloat("Min net profit / trip", &minProfit, 0.0f, 25000.0f, "%.0f")) {
                tradeMinNetProfit = (double)minProfit;
                tradeIdeasDayStamp = -1;
                tradeMixDayStamp = -1;
              }

              if (tradeUseManifest) {
                float stepKg = (float)tradeMixStepKg;
                if (ImGui::SliderFloat("Mix step (kg)", &stepKg, 0.2f, 5.0f, "%.1f")) {
                  tradeMixStepKg = std::max(0.05, (double)stepKg);
                  tradeMixDayStamp = -1;
                }
                ImGui::SameLine();
                bool impact = tradeMixSimulateImpact;
                if (ImGui::Checkbox("Price impact", &impact)) {
                  tradeMixSimulateImpact = impact;
                  tradeMixDayStamp = -1;
                }

                int linesShown = std::max(1, tradeMixLinesShown);
                if (ImGui::SliderInt("Lines shown", &linesShown, 1, 6)) {
                  tradeMixLinesShown = linesShown;
                }
              } else {
                int perStation = std::max(1, tradeIdeasPerStation);
                if (ImGui::SliderInt("Ideas / station", &perStation, 1, 3)) {
                  tradeIdeasPerStation = perStation;
                  tradeIdeasDayStamp = -1;
                }

                const char* curFilter = "Any";
                if (tradeCommodityFilter >= 0) {
                  const auto id = (econ::CommodityId)tradeCommodityFilter;
                  curFilter = econ::commodityCode(id);
                }

                if (ImGui::BeginCombo("Commodity filter", curFilter)) {
                  const bool anySelected = tradeCommodityFilter < 0;
                  if (ImGui::Selectable("Any", anySelected)) {
                    tradeCommodityFilter = -1;
                    tradeIdeasDayStamp = -1;
                  }
                  if (anySelected) ImGui::SetItemDefaultFocus();

                  for (std::size_t ci = 0; ci < econ::kCommodityCount; ++ci) {
                    const auto id = (econ::CommodityId)ci;
                    const bool selected = tradeCommodityFilter == (int)ci;
                    std::string label = std::string(econ::commodityCode(id)) + " - " + econ::commodityDef(id).name;
                    if (ImGui::Selectable(label.c_str(), selected)) {
                      tradeCommodityFilter = (int)ci;
                      tradeIdeasDayStamp = -1;
                    }
                    if (selected) ImGui::SetItemDefaultFocus();
                  }

                  ImGui::EndCombo();
                }
              }
            }

            const int dayStamp = (int)std::floor(timeDays);

            if (!tradeUseManifest) {
              if (tradeFromStationId != fromStationId || tradeIdeasDayStamp != dayStamp) {
                tradeIdeas.clear();
                tradeFromStationId = fromStationId;
                tradeIdeasDayStamp = dayStamp;

                sim::TradeScanParams scan;
                scan.radiusLy = tradeSearchRadiusLy;
                scan.maxSystems = 256;
                scan.maxResults = 12;
                scan.perStationLimit = (std::size_t)tradeIdeasPerStation;

                scan.cargoCapacityKg = cargoCapacityKg;
                scan.cargoUsedKg = cargoMassKg(cargo);
                scan.useFreeHold = tradeUseFreeHold;

                scan.bidAskSpread = 0.10;
                scan.minNetProfit = tradeMinNetProfit;
                scan.includeSameSystem = tradeIncludeSameSystem;

                if (tradeCommodityFilter >= 0) {
                  scan.hasCommodityFilter = true;
                  scan.commodityFilter = (econ::CommodityId)tradeCommodityFilter;
                }

                const auto results = sim::scanTradeOpportunities(universe,
                                                                 currentSystem->stub,
                                                                 *fromSt,
                                                                 timeDays,
                                                                 scan,
                                                                 effectiveFeeRate);

                for (const auto& r : results) {
                  TradeIdea t;
                  t.toSystem = r.toSystem;
                  t.toStation = r.toStation;
                  t.toSystemName = r.toSystemName;
                  t.toStationName = r.toStationName;
                  t.commodity = r.commodity;
                  t.buy = r.buy;
                  t.sell = r.sell;
                  t.units = r.units;
                  t.massKg = r.massKg;
                  t.netProfit = r.netProfit;
                  t.netProfitPerKg = r.netProfitPerKg;
                  t.netTripProfit = r.netTripProfit;
                  t.distanceLy = r.distanceLy;
                  tradeIdeas.push_back(std::move(t));
                }
              }
            } else {
              if (tradeMixFromStationId != fromStationId || tradeMixDayStamp != dayStamp) {
                tradeMixIdeas.clear();
                tradeMixFromStationId = fromStationId;
                tradeMixDayStamp = dayStamp;

                sim::TradeManifestScanParams scan;
                scan.radiusLy = tradeSearchRadiusLy;
                scan.maxSystems = 256;
                scan.maxResults = 12;

                scan.cargoCapacityKg = cargoCapacityKg;
                scan.cargoUsedKg = cargoMassKg(cargo);
                scan.useFreeHold = tradeUseFreeHold;

                scan.bidAskSpread = 0.10;

                scan.stepKg = tradeMixStepKg;
                scan.maxBuyCreditsCr = 0.0; // display mode: ignore player credits
                scan.simulatePriceImpact = tradeMixSimulateImpact;

                scan.minNetProfit = tradeMinNetProfit;
                scan.includeSameSystem = tradeIncludeSameSystem;

                tradeMixIdeas = sim::scanTradeManifests(universe,
                                                        currentSystem->stub,
                                                        *fromSt,
                                                        timeDays,
                                                        scan,
                                                        effectiveFeeRate);
              }
            }

            const double jr = std::max(1e-6, fsdCurrentRangeLy());
            const double usedKg = cargoMassKg(cargo);
            const double freeKg = std::max(0.0, cargoCapacityKg - usedKg);

            if (!tradeUseManifest) {
              if (tradeUseFreeHold) {
                ImGui::TextDisabled("Top trades (net, using free hold: %.0f kg free)", freeKg);
              } else {
                ImGui::TextDisabled("Top trades (net, using full hold: %.0f kg)", cargoCapacityKg);
              }

              if (tradeIdeas.empty()) {
                ImGui::TextDisabled("No profitable trades found in the selected radius.");
                ImGui::TextDisabled("Tip: expand the radius or wait a day for prices to move.");
              } else {
                ImGui::TextDisabled("Profit assumes a full fill and current prices. Use 'Fill' docked to buy.");

                if (ImGui::BeginTable("trade_table", 6,
                                      ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable |
                                          ImGuiTableFlags_ScrollY,
                                      ImVec2(0, 260))) {
                  ImGui::TableSetupColumn("To");
                  ImGui::TableSetupColumn("Commodity");
                  ImGui::TableSetupColumn("Buy");
                  ImGui::TableSetupColumn("Sell");
                  ImGui::TableSetupColumn("Net profit");
                  ImGui::TableSetupColumn("Actions");
                  ImGui::TableHeadersRow();

                  for (std::size_t i = 0; i < tradeIdeas.size(); ++i) {
                    const auto& t = tradeIdeas[i];

                    ImGui::TableNextRow();

                    ImGui::TableSetColumnIndex(0);
                    std::string toLabel = t.toSystemName + " / " + t.toStationName;
                    ImGui::TextUnformatted(toLabel.c_str());
                    const double jumps = t.distanceLy / jr;
                    ImGui::TextDisabled("%.1f ly (%.1f jumps)", t.distanceLy, jumps);

                    ImGui::TableSetColumnIndex(1);
                    ImGui::Text("%s", econ::commodityDef(t.commodity).name.c_str());
                    ImGui::TextDisabled("%s", econ::commodityCode(t.commodity));

                    ImGui::TableSetColumnIndex(2);
                    ImGui::Text("%.1f", t.buy);
                    ImGui::TextDisabled("x%.1f", t.units);

                    ImGui::TableSetColumnIndex(3);
                    ImGui::Text("%.1f", t.sell);

                    ImGui::TableSetColumnIndex(4);
                    ImGui::Text("+%.0f", t.netTripProfit);
                    ImGui::TextDisabled("+%.1f/kg", t.netProfitPerKg);

                    ImGui::TableSetColumnIndex(5);
                    ImGui::PushID((int)i);

                    if (ImGui::SmallButton("Plot")) {
                      plotRouteToSystem(t.toSystem);
                    }

                    if (fromDocked) {
                      ImGui::SameLine();
                      if (ImGui::SmallButton("Fill")) {
                        // Buy enough of the commodity to fill the hold (or free hold), limited by inventory & credits.
                        auto& fromEcon = universe.stationEconomy(*fromSt, timeDays);

                        const double unitMass = std::max(1e-6, econ::commodityDef(t.commodity).massKg);
                        const double freeNowKg = std::max(0.0, cargoCapacityKg - cargoMassKg(cargo));
                        const double maxUnitsHold = freeNowKg / unitMass;

                        // unitsPossible already accounts for fees/spread and station caps, but may exceed free hold.
                        const double unitsToBuy = std::min({t.units, maxUnitsHold, fromEcon.inventory[(std::size_t)t.commodity]});

                        const auto tr = econ::buy(fromEcon,
                                                  fromSt->economyModel,
                                                  t.commodity,
                                                  unitsToBuy,
                                                  credits,
                                                  0.10,
                                                  feeEff);

                        if (tr.ok && tr.unitsDelta > 1e-9) {
                          cargo[(std::size_t)t.commodity] += tr.unitsDelta;
                          toast(toasts, ("Bought " + std::to_string((int)std::round(tr.unitsDelta)) + " " +
                                         econ::commodityDef(t.commodity).name + " (" +
                                         std::to_string((int)std::round(-tr.creditsDelta)) + " cr)"));
                          tradeIdeasDayStamp = -1;
                          tradeMixDayStamp = -1;
                        } else {
                          toast(toasts, ("Can't buy: " + tr.reason));
                        }
                      }
                    }

                    ImGui::PopID();
                  }

                  ImGui::EndTable();
                }
              }
            } else {
              const double capShownKg = tradeUseFreeHold ? freeKg : cargoCapacityKg;

              ImGui::TextDisabled("Top cargo mixes (net, step=%.1f kg, impact=%s)", tradeMixStepKg,
                                  tradeMixSimulateImpact ? "yes" : "no");

              if (tradeMixIdeas.empty()) {
                ImGui::TextDisabled("No profitable cargo mixes found in the selected radius.");
                ImGui::TextDisabled("Tip: expand the radius or wait a day for prices to move.");
              } else {
                ImGui::TextDisabled("Profit assumes a full fill and current prices. Use 'Fill mix' docked to buy.");

                if (ImGui::BeginTable("trade_mix_table", 5,
                                      ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable |
                                          ImGuiTableFlags_ScrollY,
                                      ImVec2(0, 260))) {
                  ImGui::TableSetupColumn("To");
                  ImGui::TableSetupColumn("Mix");
                  ImGui::TableSetupColumn("Fill");
                  ImGui::TableSetupColumn("Net profit");
                  ImGui::TableSetupColumn("Actions");
                  ImGui::TableHeadersRow();

                  for (std::size_t i = 0; i < tradeMixIdeas.size(); ++i) {
                    const auto& t = tradeMixIdeas[i];

                    ImGui::TableNextRow();

                    ImGui::TableSetColumnIndex(0);
                    std::string toLabel = t.toSystemName + " / " + t.toStationName;
                    ImGui::TextUnformatted(toLabel.c_str());
                    const double jumps = t.distanceLy / jr;
                    ImGui::TextDisabled("%.1f ly (%.1f jumps)", t.distanceLy, jumps);

                    ImGui::TableSetColumnIndex(1);
                    ImGui::BeginGroup();
                    const int take = std::min<int>(std::max(1, tradeMixLinesShown), (int)t.lines.size());
                    for (int li = 0; li < take; ++li) {
                      const auto& ln = t.lines[(std::size_t)li];
                      ImGui::Text("%s  %.1fu", econ::commodityCode(ln.commodity), ln.units);
                      ImGui::TextDisabled("+%.0f (%.0f/kg)", ln.netProfitCr, ln.netProfitPerKg);
                    }
                    if ((int)t.lines.size() > take) {
                      ImGui::TextDisabled("... +%d more", (int)t.lines.size() - take);
                    }
                    ImGui::EndGroup();

                    if (ImGui::IsItemHovered()) {
                      ImGui::BeginTooltip();
                      ImGui::TextUnformatted("Manifest details:");
                      ImGui::Separator();
                      for (const auto& ln : t.lines) {
                        ImGui::Text("%s  units=%.2f  mass=%.1fkg  net=%.0f",
                                    econ::commodityCode(ln.commodity),
                                    ln.units,
                                    ln.massKg,
                                    ln.netProfitCr);
                        ImGui::TextDisabled("  buy=%.1f  sell=%.1f  net/kg=%.0f",
                                            ln.avgNetBuyPrice,
                                            ln.avgNetSellPrice,
                                            ln.netProfitPerKg);
                      }
                      ImGui::EndTooltip();
                    }

                    ImGui::TableSetColumnIndex(2);
                    ImGui::Text("%.0f / %.0f kg", t.cargoFilledKg, capShownKg);
                    ImGui::TextDisabled("buy %.0f", t.netBuyCr);

                    ImGui::TableSetColumnIndex(3);
                    const double tripPerKg = (t.cargoFilledKg > 1e-6) ? (t.netProfitCr / t.cargoFilledKg) : 0.0;
                    ImGui::Text("+%.0f", t.netProfitCr);
                    ImGui::TextDisabled("+%.0f/kg", tripPerKg);

                    ImGui::TableSetColumnIndex(4);
                    ImGui::PushID((int)i);

                    if (ImGui::SmallButton("Plot")) {
                      plotRouteToSystem(t.toSystem);
                    }

                    if (fromDocked) {
                      ImGui::SameLine();
                      if (ImGui::SmallButton("Fill mix")) {
                        auto& fromEcon = universe.stationEconomy(*fromSt, timeDays);

                        const auto& toSys = universe.getSystem(t.toSystem);
                        const sim::Station* toStPtr = nullptr;
                        for (const auto& st : toSys.stations) {
                          if (st.id == t.toStation) {
                            toStPtr = &st;
                            break;
                          }
                        }

                        if (!toStPtr) {
                          toast(toasts, "Destination station not found.");
                        } else {
                          auto& toEcon = universe.stationEconomy(*toStPtr, timeDays);

                          const double freeNowKg = std::max(0.0, cargoCapacityKg - cargoMassKg(cargo));
                          econ::CargoManifestParams mp;
                          mp.cargoCapacityKg = freeNowKg;
                          mp.bidAskSpread = 0.10;
                          mp.fromFeeRate = feeEff;
                          mp.toFeeRate = effectiveFeeRate(*toStPtr);
                          mp.stepKg = tradeMixStepKg;
                          mp.maxBuyCreditsCr = credits;
                          mp.simulatePriceImpact = tradeMixSimulateImpact;

                          const auto plan = econ::bestManifestForCargo(fromEcon,
                                                                       fromSt->economyModel,
                                                                       toEcon,
                                                                       toStPtr->economyModel,
                                                                       mp);

                          if (plan.lines.empty()) {
                            toast(toasts, "No profitable cargo mix available to buy right now.");
                          } else {
                            double boughtKg = 0.0;
                            double spentCr = 0.0;

                            for (const auto& ln : plan.lines) {
                              const double unitMass = std::max(1e-6, econ::commodityDef(ln.commodity).massKg);
                              const double freeKgNow = std::max(0.0, cargoCapacityKg - cargoMassKg(cargo));
                              if (freeKgNow <= 1e-6) break;

                              const double maxUnitsHold = freeKgNow / unitMass;
                              const double unitsToBuy = std::min(ln.units, maxUnitsHold);
                              if (unitsToBuy <= 1e-9) continue;

                              const auto tr = econ::buy(fromEcon,
                                                        fromSt->economyModel,
                                                        ln.commodity,
                                                        unitsToBuy,
                                                        credits,
                                                        0.10,
                                                        feeEff);

                              if (tr.ok && tr.unitsDelta > 1e-9) {
                                cargo[(std::size_t)ln.commodity] += tr.unitsDelta;
                                boughtKg += tr.unitsDelta * unitMass;
                                spentCr += -tr.creditsDelta;
                              }
                            }

                            if (boughtKg > 1e-6) {
                              toast(toasts, ("Bought mix (" + std::to_string((int)std::round(boughtKg)) + " kg, " +
                                             std::to_string((int)std::round(spentCr)) + " cr)"));
                              tradeIdeasDayStamp = -1;
                              tradeMixDayStamp = -1;
                            } else {
                              toast(toasts, "Couldn't buy any of the suggested mix (no credits / no stock / no hold).");
                            }
                          }
                        }
                      }
                    }

                    ImGui::PopID();
                  }

                  ImGui::EndTable();
                }
              }
            }

          }
        }
      }

      ImGui::End();
    }

    if (showMissions) {
      ImGui::Begin("Missions");

      auto stationNameById = [&](sim::SystemId sysId, sim::StationId stId) -> std::string {
        if (stId == 0) return "";
        const auto& sys = universe.getSystem(sysId);
        for (const auto& st : sys.stations) {
          if (st.id == stId) return st.name;
        }
        return "Station #" + std::to_string((std::uint64_t)stId);
      };

      auto systemNameById = [&](sim::SystemId sysId) -> std::string {
        if (sysId == 0) return "";
        return universe.getSystem(sysId).stub.name;
      };

      auto describeMission = [&](const sim::Mission& m) -> std::string {
        std::string out;
        switch (m.type) {
          case sim::MissionType::Courier: {
            out = "Courier: Deliver data to " + stationNameById(m.toSystem, m.toStation) + " (" + systemNameById(m.toSystem) + ")";
          } break;
          case sim::MissionType::Delivery: {
            const econ::CommodityId cid = m.commodity;
            out = "Delivery: Deliver " + std::to_string((int)m.units) + " units of " + std::string(econ::commodityDef(cid).name)
                + " to " + stationNameById(m.toSystem, m.toStation) + " (" + systemNameById(m.toSystem) + ")";
          } break;
          case sim::MissionType::Salvage: {
            const econ::CommodityId cid = m.commodity;
            const std::string dst = stationNameById(m.toSystem, m.toStation) + " (" + systemNameById(m.toSystem) + ")";
            out = std::string("Salvage: ") + (m.scanned ? "Deliver " : "Recover ")
                + std::to_string((int)m.units) + " units of " + std::string(econ::commodityDef(cid).name)
                + " -> " + dst;
            if (!m.scanned) out += " [visit mission derelict]";
          } break;
          case sim::MissionType::Smuggle: {
            const econ::CommodityId cid = m.commodity;
            out = "Smuggle: Deliver " + std::to_string((int)m.units) + " units of " + std::string(econ::commodityDef(cid).name)
                + " to " + stationNameById(m.toSystem, m.toStation) + " (" + systemNameById(m.toSystem) + ") [CONTRABAND]";
          } break;
          case sim::MissionType::Escort: {
            const econ::CommodityId cid = m.commodity;
            const std::string dst = stationNameById(m.toSystem, m.toStation) + " (" + systemNameById(m.toSystem) + ")";
            out = std::string("Escort: ") + (m.scanned ? "Report in " : "Protect convoy ");
            if (m.units > 0.0) {
              out += std::to_string((int)m.units) + " units of " + std::string(econ::commodityDef(cid).name) + " -> " + dst;
            } else {
              out += "to " + dst;
            }
            if (!m.scanned) out += " [stay close]";
          } break;
          case sim::MissionType::MultiDelivery: {
            const econ::CommodityId cid = m.commodity;
            out = "Multi-hop: Deliver " + std::to_string((int)m.units) + " units of " + std::string(econ::commodityDef(cid).name);
            if (m.viaSystem != 0 && m.viaStation != 0) {
              out += " via " + stationNameById(m.viaSystem, m.viaStation) + " (" + systemNameById(m.viaSystem) + ")";
            }
            out += " -> " + stationNameById(m.toSystem, m.toStation) + " (" + systemNameById(m.toSystem) + ")";
          } break;
          case sim::MissionType::Passenger: {
            out = "Passenger: Transport " + std::to_string((int)std::llround(m.units)) + " pax to "
                + stationNameById(m.toSystem, m.toStation) + " (" + systemNameById(m.toSystem) + ")";
          } break;
          case sim::MissionType::BountyScan: {
            out = "Bounty: Scan wanted pirate in " + systemNameById(m.toSystem);
          } break;
          case sim::MissionType::BountyKill: {
            out = "Bounty: Eliminate wanted pirate in " + systemNameById(m.toSystem);
          } break;
          default:
            out = "Mission";
            break;
        }
        return out;
      };

      if (ImGui::BeginTabBar("missions_tabs")) {
        if (ImGui::BeginTabItem("Active")) {
          if (missions.empty()) {
            ImGui::TextDisabled("No missions accepted.");
          }

          for (auto& m : missions) {
            ImGui::PushID((int)m.id);
            const bool active = !(m.completed || m.failed);
            const char* status = m.completed ? "COMPLETED" : (m.failed ? "FAILED" : "ACTIVE");
            // Procedural mission icon(s)
            {
              core::u64 mSeed = core::hashCombine(core::fnv1a64("mission"), (core::u64)m.id);
              mSeed = core::hashCombine(mSeed, (core::u64)(core::i64)m.type);
              mSeed = core::hashCombine(mSeed, (core::u64)m.factionId);
              const auto& mIcon = spriteCache.get(render::SpriteKind::Mission, mSeed, 32);
              ImGui::Image((ImTextureID)(intptr_t)mIcon.handle(), ImVec2(18, 18));
              if (ImGui::IsItemHovered()) {
                const auto& big = spriteCache.get(render::SpriteKind::Mission, mSeed, 96);
                ImGui::BeginTooltip();
                ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(96, 96));
                ImGui::EndTooltip();
              }

              if (m.type == sim::MissionType::Delivery || m.type == sim::MissionType::MultiDelivery ||
                  m.type == sim::MissionType::Smuggle || m.type == sim::MissionType::Salvage) {
                const core::u64 cSeed = core::hashCombine(core::fnv1a64("commodity"), (core::u64)(core::i64)m.commodity);
                const auto& cIcon = spriteCache.get(render::SpriteKind::Commodity, cSeed, 32);
                ImGui::SameLine();
                ImGui::Image((ImTextureID)(intptr_t)cIcon.handle(), ImVec2(18, 18));
                if (ImGui::IsItemHovered()) {
                  const auto& cBig = spriteCache.get(render::SpriteKind::Commodity, cSeed, 96);
                  ImGui::BeginTooltip();
                  ImGui::Image((ImTextureID)(intptr_t)cBig.handle(), ImVec2(96, 96));
                  ImGui::EndTooltip();
                }
              }

              ImGui::SameLine();
              ImGui::Text("[%s]%s %s", status, (trackedMissionId == m.id ? " [TRACKED]" : ""), describeMission(m).c_str());
            }
            ImGui::TextDisabled("Reward %.0f cr | Faction: %s (rep %.1f)", m.reward, factionName(m.factionId).c_str(), getRep(m.factionId));

            if (m.deadlineDay > 0.0) {
              const double hrsLeft = (m.deadlineDay - timeDays) * 24.0;
              ImGui::TextDisabled("Deadline: day %.2f (%.1f h left)", m.deadlineDay, hrsLeft);
            }

            if (active && (m.type == sim::MissionType::Delivery || m.type == sim::MissionType::MultiDelivery || m.type == sim::MissionType::Smuggle)) {
              const std::size_t ci = (std::size_t)m.commodity;
              const double have = (ci < econ::kCommodityCount) ? cargo[ci] : 0.0;
              const double need = m.units;
              if (have + 1e-6 < need) {
                ImGui::TextColored(ImVec4(1.0f, 0.55f, 0.55f, 1.0f), "Cargo: %.0f / %.0f (missing %.0f)", have, need, need - have);
              } else {
                ImGui::TextDisabled("Cargo: %.0f / %.0f", have, need);
              }
            }

            if (active && m.type == sim::MissionType::MultiDelivery && m.viaSystem != 0 && m.viaStation != 0) {
              ImGui::TextDisabled("Progress: leg %d / 2", (int)m.leg + 1);
            }

            if (active) {
              if (ImGui::SmallButton("Set destination")) {
                sim::SystemId destSys = m.toSystem;
                if (m.type == sim::MissionType::MultiDelivery && m.viaSystem != 0 && m.leg == 0) {
                  destSys = m.viaSystem;
                }
                galaxySelectedSystemId = destSys;
                showGalaxy = true;
                toast(toasts, "Galaxy: mission destination selected.", 2.0);
              }
              ImGui::SameLine();
              if (ImGui::SmallButton("Plot route")) {
                sim::SystemId destSys = m.toSystem;
                sim::StationId destSt = m.toStation;
                if (m.type == sim::MissionType::MultiDelivery && m.viaSystem != 0 && m.viaStation != 0 && m.leg == 0) {
                  destSys = m.viaSystem;
                  destSt = m.viaStation;
                }
                if (plotRouteToSystem(destSys)) {
                  pendingArrivalTargetStationId = destSt;
                  if (currentSystem && destSys == currentSystem->stub.id && destSt != 0) tryTargetStationById(destSt);
                  showGalaxy = true;
                }
              }
              ImGui::SameLine();
              if (ImGui::SmallButton(trackedMissionId == m.id ? "Untrack" : "Track")) {
                if (trackedMissionId == m.id) {
                  trackedMissionId = 0;
                  toast(toasts, "Mission untracked.", 1.6);
                } else {
                  trackedMissionId = m.id;
                  toast(toasts, "Mission tracked (shown in Ship/Status).", 2.0);
                }
              }
              ImGui::SameLine();
              if (ImGui::SmallButton("Abandon")) {
                m.failed = true;
                addRep(m.factionId, -2.0);
                toast(toasts, "Mission abandoned.", 2.0);
              }
            }

            ImGui::Separator();
            ImGui::PopID();
          }

          ImGui::EndTabItem();
        }

        if (ImGui::BeginTabItem("Mission Board")) {
          if (!docked || dockedStationId == 0) {
            ImGui::TextDisabled("Dock at a station to browse missions.");
            ImGui::EndTabItem();
          } else {
            // Resolve docked station
            int dockedIdx = -1;
            for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
              if (currentSystem->stations[i].id == dockedStationId) {
                dockedIdx = (int)i;
                break;
              }
            }

            if (dockedIdx < 0) {
              ImGui::TextDisabled("Docked station not found.");
              ImGui::EndTabItem();
            } else {
              const auto& st = currentSystem->stations[(std::size_t)dockedIdx];
              const double rep = getRep(st.factionId);
              ImGui::Text("Station: %s", st.name.c_str());
              ImGui::Text("Faction: %s (rep %.1f)", factionName(st.factionId).c_str(), rep);

              // Mission board offers are generated by the shared headless mission module.
              // This keeps offer weights, accept rules, and future mission types consistent
              // between the prototype UI and unit tests.
              {
                sim::SaveGame board{};
                board.missionOffersStationId = missionOffersStationId;
                board.missionOffersDayStamp = missionOffersDayStamp;
                board.missionOffers = missionOffers;
                sim::refreshMissionOffers(universe, *currentSystem, st, timeDays, rep, board);
                missionOffersStationId = board.missionOffersStationId;
                missionOffersDayStamp = board.missionOffersDayStamp;
                missionOffers = std::move(board.missionOffers);
              }

              if (missionOffers.empty()) {
                ImGui::TextDisabled("No offers available right now.");
              }

              const bool dockedHere = docked && (dockedStationId == st.id);
              if (!dockedHere) {
                ImGui::TextDisabled("Dock at this station to accept missions (G to dock).");
              }

              const double feeEff = effectiveFeeRate(st);
              auto& stEcon = universe.stationEconomy(st, timeDays);

              static bool autoTrackAccepted = true;
              ImGui::Checkbox("Auto-track accepted missions", &autoTrackAccepted);
              ImGui::SameLine();
              ImGui::Checkbox("Show route preview (A*)", &missionOfferRoutePreviewEnabled);
              ImGui::SameLine();
              ImGui::TextDisabled("(jumps / distance / fuel to next stop)");

              // Recompute cached route previews when offers or planner settings change.
              const double jrPreviewLy = navConstrainToCurrentFuelRange ? fsdCurrentRangeLy() : fsdBaseRangeLy();
              const bool previewKeyChanged = (missionOfferRoutePreviewStationId != missionOffersStationId)
                                        || (missionOfferRoutePreviewDayStamp != missionOffersDayStamp)
                                        || (missionOfferRoutePreviewFromSystem != currentSystem->stub.id)
                                        || (missionOfferRoutePreviewMode != navRouteMode)
                                        || (missionOfferRoutePreviewConstrainToCurrentFuelRange != navConstrainToCurrentFuelRange)
                                        || (std::abs(missionOfferRoutePreviewMaxJumpLy - jrPreviewLy) > 1e-6)
                                        || (missionOfferRoutePreview.size() != missionOffers.size());

              if (missionOfferRoutePreviewEnabled && previewKeyChanged) {
                missionOfferRoutePreviewStationId = missionOffersStationId;
                missionOfferRoutePreviewDayStamp = missionOffersDayStamp;
                missionOfferRoutePreviewFromSystem = currentSystem->stub.id;
                missionOfferRoutePreviewMode = navRouteMode;
                missionOfferRoutePreviewConstrainToCurrentFuelRange = navConstrainToCurrentFuelRange;
                missionOfferRoutePreviewMaxJumpLy = jrPreviewLy;

                missionOfferRoutePreview.clear();
                missionOfferRoutePreview.resize(missionOffers.size());

                for (std::size_t pi = 0; pi < missionOffers.size(); ++pi) {
                  const auto& offer = missionOffers[pi];
                  auto& out = missionOfferRoutePreview[pi];

                  const bool hasVia = (offer.type == sim::MissionType::MultiDelivery
                                    && offer.viaSystem != 0 && offer.viaStation != 0
                                    && offer.leg == 0);
                  const sim::SystemId destSys = hasVia ? offer.viaSystem : offer.toSystem;

                  if (!currentSystem || destSys == 0) {
                    out.ok = false;
                    continue;
                  }

                  if (destSys == currentSystem->stub.id) {
                    out.ok = true;
                    out.jumps = 0;
                    out.distanceLy = 0.0;
                    out.fuel = 0.0;
                    continue;
                  }

                  if (jrPreviewLy <= 0.0) {
                    out.ok = false;
                    continue;
                  }

                  const auto& destStub = universe.getSystem(destSys).stub;
                  const double directLy = (destStub.posLy - currentSystem->stub.posLy).length();
                  double radius = std::clamp(directLy * 1.20 + 60.0, 180.0, 1400.0);

                  for (int attempt = 0; attempt < 6; ++attempt) {
                    auto nearby = universe.queryNearby(currentSystem->stub.posLy, radius);
                    sim::RoutePlanStats stats{};
                    std::vector<sim::SystemId> route;

                    if (navRouteMode == NavRouteMode::Hops) {
                      route = sim::plotRouteAStarHops(nearby, currentSystem->stub.id, destSys, jrPreviewLy, &stats);
                    } else if (navRouteMode == NavRouteMode::Distance) {
                      route = sim::plotRouteAStarCost(nearby, currentSystem->stub.id, destSys, jrPreviewLy, 0.0, 1.0, &stats);
                    } else {
                      route = sim::plotRouteAStarCost(nearby, currentSystem->stub.id, destSys, jrPreviewLy, kFsdFuelBase, kFsdFuelPerLy, &stats);
                    }

                    if (!route.empty()) {
                      out.ok = true;
                      out.jumps = (int)route.size() - 1;
                      out.distanceLy = sim::routeDistanceLy(nearby, route);
                      out.fuel = sim::routeCost(nearby, route, kFsdFuelBase, kFsdFuelPerLy);
                      break;
                    }

                    radius *= 1.35;
                  }
                }
              }

              auto acceptOffer = [&](std::size_t offerIdx, bool autoPlot) -> bool {
                if (!dockedHere) {
                  toast(toasts, "Dock at this station to accept missions.", 2.2);
                  return false;
                }
                if (offerIdx >= missionOffers.size()) return false;

                // Keep a copy for route plotting (accept may erase the offer).
                const sim::Mission chosen = missionOffers[offerIdx];

                // Delegate accept checks (cargo + inventory + credits + passenger seats) to the
                // shared headless module so game/UI stays consistent with tests.
                sim::SaveGame tmp{};
                tmp.credits = credits;
                tmp.cargo = cargo;
                tmp.cargoCapacityKg = cargoCapacityKg;
                tmp.passengerSeats = passengerSeats;
                tmp.nextMissionId = nextMissionId;
                tmp.missions = missions;
                tmp.missionOffers = missionOffers;

                std::string err;
                if (!sim::acceptMissionOffer(universe, st, timeDays, rep, tmp, offerIdx, &err)) {
                  toast(toasts, err.empty() ? "Couldn't accept mission." : err, 3.0);
                  return false;
                }

                credits = tmp.credits;
                cargo = tmp.cargo;
                nextMissionId = tmp.nextMissionId;
                missions = std::move(tmp.missions);
                missionOffers = std::move(tmp.missionOffers);
                if (autoTrackAccepted && !missions.empty()) {
                  trackedMissionId = missions.back().id;
                }
                toast(toasts, "Mission accepted.", 2.0);

                if (autoPlot) {
                  const bool hasVia = (chosen.viaSystem != 0 && chosen.viaStation != 0 && chosen.leg == 0);
                  const sim::SystemId destSys = hasVia ? chosen.viaSystem : chosen.toSystem;
                  const sim::StationId destSt = hasVia ? chosen.viaStation : chosen.toStation;
                  plotRouteToSystem(destSys);
                  pendingArrivalTargetStationId = destSt;
                  if (destSys == currentSystem->stub.id) {
                    tryTargetStationById(destSt);
                  }
                  showGalaxy = true;
                }

                return true;
              };

              // Starter helper: pick a simple Courier/Delivery we can complete from here.
              int starterIdx = -1;
              double starterScore = 1e99;
              const double jrLy = std::max(1e-6, fsdCurrentRangeLy());
              bool acceptedOfferThisFrame = false;
              for (std::size_t i = 0; i < missionOffers.size(); ++i) {
                const auto& o = missionOffers[i];
                if (o.type != sim::MissionType::Courier && o.type != sim::MissionType::Delivery) continue;

                // Avoid complex / multi-leg missions in the starter picker.
                if (o.viaSystem != 0 || o.viaStation != 0) continue;

                if (o.type == sim::MissionType::Delivery) {
                  const double massKg = econ::commodityDef(o.commodity).massKg;
                  const double needKg = massKg * (double)o.units;
                  if (needKg > cargoCapacityKg + 1e-6) continue;

                  if (o.cargoProvided) {
                    if (stEcon.inventory[(int)o.commodity] + 1e-6 < (double)o.units) continue;
                  } else {
                    const auto q = econ::quote(stEcon, st.economyModel, o.commodity, 0.10);
                    const double totalCost = q.ask * (double)o.units * (1.0 + feeEff);
                    if (q.inventory + 1e-6 < (double)o.units) continue;
                    if (credits + 1e-6 < totalCost) continue;
                  }
                }

                const sim::SystemId destSys = o.toSystem;
                const auto& destStub = universe.getSystem(destSys).stub;
                const double distLy = (destStub.posLy - currentSystem->stub.posLy).length();
                const int estJumps = (int)std::ceil(distLy / jrLy);

                // Heuristic: prefer fewer jumps and shorter distance; mildly prefer higher reward.
                const double score = (double)estJumps * 100000.0 + distLy * 200.0 - (double)o.reward;
                if (score < starterScore) {
                  starterScore = score;
                  starterIdx = (int)i;
                }
              }

              if (starterIdx >= 0) {
                ImGui::Separator();
                ImGui::TextDisabled("Starter");
                ImGui::BeginDisabled(!dockedHere);
                if (ImGui::Button("Accept starter mission & auto-plot")) {
                  acceptOffer((std::size_t)starterIdx, true);
                }
                ImGui::EndDisabled();
                ImGui::SameLine();
                ImGui::TextDisabled("(one click: accepts a simple Courier/Delivery and plots a route)");
              }

              const int paxUsed = [&] {
                int used = 0;
                for (const auto& m : missions) {
                  if (m.completed || m.failed) continue;
                  if (m.type != sim::MissionType::Passenger) continue;
                  used += std::max(0, (int)std::llround(m.units));
                }
                return used;
              }();
              const int paxCap = std::max(0, passengerSeats);

              for (std::size_t i = 0; i < missionOffers.size(); ++i) {
                const auto offer = missionOffers[i];
                ImGui::PushID((int)i);

                ImGui::TextWrapped("%s", describeMission(offer).c_str());
                if (offer.deadlineDay > 0.0) {
                  const double hrsLeft = (offer.deadlineDay - timeDays) * 24.0;
                  ImGui::TextDisabled("Reward %.0f cr | Deadline in %.1f h", offer.reward, hrsLeft);
                } else {
                  ImGui::TextDisabled("Reward %.0f cr", offer.reward);
                }

                if (missionOfferRoutePreviewEnabled && i < missionOfferRoutePreview.size()) {
                  const auto& rp = missionOfferRoutePreview[i];
                  if (rp.ok) {
                    if (rp.jumps <= 0) {
                      ImGui::TextDisabled("Route: local system");
                    } else {
                      ImGui::TextDisabled("Route: ~%d jumps | %.0f ly | est fuel %.1f", rp.jumps, rp.distanceLy, rp.fuel);
                    }
                  } else {
                    ImGui::TextDisabled("Route: (no route w/ current planner settings)");
                  }
                }

                bool canAccept = (missions.size() < 16);

                // Basic pre-flight validation for cargo missions (so the "Accept" button behaviour matches the toast checks).
                const bool isCargoJob = (offer.type == sim::MissionType::Delivery
                                      || offer.type == sim::MissionType::MultiDelivery
                                      || offer.type == sim::MissionType::Smuggle);
                if (isCargoJob && offer.units > 0.0) {
                  const econ::CommodityId cid = offer.commodity;
                  const double massKg = econ::commodityDef(cid).massKg;
                  const double addKg = (double)offer.units * massKg;
                  canAccept = canAccept && (cargoMassKg(cargo) + addKg <= cargoCapacityKg + 1e-6);

                  // Smuggling cargo is provided by the contact; it doesn't require legal market inventory or credits.
                  if (offer.type != sim::MissionType::Smuggle) {
                    if (offer.cargoProvided) {
                      // Cargo-provided offers must be backed by real station inventory.
                      canAccept = canAccept && (stEcon.inventory[(int)cid] + 1e-6 >= (double)offer.units);
                    } else {
                      const auto q = econ::quote(stEcon, st.economyModel, cid, 0.10);
                      const double totalCost = q.ask * (double)offer.units * (1.0 + feeEff);
                      canAccept = canAccept && (q.inventory + 1e-6 >= (double)offer.units) && (credits + 1e-6 >= totalCost);
                    }
                  }
                }

                // Passenger jobs require enough free seats.
                if (offer.type == sim::MissionType::Passenger) {
                  const int party = std::max(0, (int)std::llround(offer.units));
                  canAccept = canAccept && (party > 0) && (paxUsed + party <= paxCap);
                }

                const bool acceptDisabled = (!dockedHere) || (!canAccept);
                if (acceptDisabled) ImGui::BeginDisabled();
                if (ImGui::SmallButton("Accept & Plot")) {
                  acceptedOfferThisFrame = acceptOffer(i, true);
                }
                if (acceptDisabled) ImGui::EndDisabled();
                ImGui::SameLine();
                if (ImGui::SmallButton("Plot")) {
                  const sim::SystemId legSys = (offer.viaSystem != 0 && offer.viaStation != 0 && offer.leg == 0) ? offer.viaSystem : offer.toSystem;
                  const sim::StationId legSt = (offer.viaSystem != 0 && offer.viaStation != 0 && offer.leg == 0) ? offer.viaStation : offer.toStation;
                  if (plotRouteToSystem(legSys)) {
                    pendingArrivalTargetStationId = legSt;
                    if (legSys == currentSystem->stub.id) tryTargetStationById(legSt);
                    showGalaxy = true;
                  }
                }

                ImGui::Separator();
                ImGui::PopID();
              }

              ImGui::EndTabItem();
            }
          }
        }

        ImGui::EndTabBar();
      }

      ImGui::End();
    }

if (showContacts) {
  ImGui::Begin("Contacts / Combat");

  int alivePirates = 0;
  int aliveTraders = 0;
  int alivePolice = 0;
  for (const auto& c : contacts) {
    if (!c.alive) continue;
    if (c.role == ContactRole::Pirate) ++alivePirates;
    if (c.role == ContactRole::Trader) ++aliveTraders;
    if (c.role == ContactRole::Police) ++alivePolice;
  }

  const core::u32 localFaction = currentSystem ? currentSystem->stub.factionId : 0;
  if (localFaction != 0) {
    const core::u64 fSeed = core::hashCombine(core::fnv1a64("faction"), (core::u64)localFaction);
    const auto& fIcon = spriteCache.get(render::SpriteKind::Faction, fSeed, 32);
    ImGui::Image((ImTextureID)(intptr_t)fIcon.handle(), ImVec2(18, 18));
    if (ImGui::IsItemHovered()) {
      const auto& big = spriteCache.get(render::SpriteKind::Faction, fSeed, 96);
      ImGui::BeginTooltip();
      ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(96, 96));
      ImGui::EndTooltip();
    }
    ImGui::SameLine();
    ImGui::Text("Local faction: %s | Rep %.1f | Bounty %.0f | Alert %.1f",
                factionName(localFaction).c_str(),
                getRep(localFaction),
                getBounty(localFaction),
                policeHeat);
  } else {
    ImGui::Text("Local faction: (none)");
  }

  ImGui::Text("Pirates: %d   Traders: %d   Police: %d", alivePirates, aliveTraders, alivePolice);

  if (ImGui::Button("Panic: clear contacts")) {
    for (auto& c : contacts) c.alive = false;
    toast(toasts, "Contacts cleared.", 2.0);
  }

  ImGui::Separator();

  for (std::size_t i = 0; i < contacts.size(); ++i) {
    const auto& c = contacts[i];
    if (!c.alive) continue;

    const double distKm = (c.ship.positionKm() - ship.positionKm()).length();

    ImGui::PushID((int)i);

    std::string tag = std::string("[") + contactRoleName(c.role) + "]";
    if (c.missionTarget) tag += " [BOUNTY]";
    if (c.escortConvoy) tag += " [CONVOY]";
    if (c.followId != 0 && c.role == ContactRole::Police) tag += " [ESCORT]";
    if (c.groupId != 0) {
      tag += (c.leaderId == 0) ? " [LEAD]" : " [WING]";
    }
    if (c.role == ContactRole::Police) {
      const bool hostile = c.hostileToPlayer || (c.factionId != 0 && getBounty(c.factionId) > 0.0);
      if (hostile) tag += " [HOSTILE]";
    }

    if (c.role == ContactRole::Pirate) {
      if (pirateDemand.active && c.groupId != 0 && c.groupId == pirateDemand.groupId && !c.hostileToPlayer) {
        tag += " [DEMAND]";
      } else if (c.hostileToPlayer || c.missionTarget) {
        tag += " [HOSTILE]";
      }
    }

    ImGui::Text("%s %s", c.name.c_str(), tag.c_str());
    if (c.factionId != 0) {
      ImGui::SameLine();
      ImGui::TextDisabled("(%s)", factionName(c.factionId).c_str());
    }

    ImGui::SameLine();
    if (ImGui::SmallButton("Target")) {
      target.kind = TargetKind::Contact;
      target.index = i;
    }

    ImGui::TextDisabled("Dist %.0f km | Hull %.0f | Shield %.0f", distKm, c.hull, c.shield);

    if (c.role == ContactRole::Trader) {
      std::string destName = "?";
      if (currentSystem && c.tradeDestStationIndex < currentSystem->stations.size()) {
        destName = currentSystem->stations[c.tradeDestStationIndex].name;
      }

      const char* cName = econ::commodityDef(c.tradeCommodity).name;
      if (c.tradeUnits > 0.0) {
        ImGui::TextDisabled("Hauling %.0fu %s -> %s | Cargo value ~%.0f cr (piracy is illegal).",
                            c.tradeUnits,
                            cName,
                            destName.c_str(),
                            c.cargoValueCr);
      } else {
        ImGui::TextDisabled("En route -> %s | Empty hold (piracy is illegal).", destName.c_str());
      }
    }

    ImGui::Separator();
    ImGui::PopID();
  }

  ImGui::End();
}
    if (showGalaxy) {
      ImGui::Begin("Galaxy / Streaming");

      const auto center = currentSystem->stub.posLy;
      static float radius = 200.0f;
      ImGui::SliderFloat("Query radius (ly)", &radius, 20.0f, 2000.0f);

      auto nearby = universe.queryNearby(center, radius, 256);

      const double galaxyJrMaxLy = fsdBaseRangeLy();
      const double galaxyJrNowLy = fsdCurrentRangeLy();
      const double galaxyPlanMaxLy = navConstrainToCurrentFuelRange ? galaxyJrNowLy : galaxyJrMaxLy;

      ImGui::Text("Nearby systems: %d", (int)nearby.size());

      // Plot-route helper so both the detail panel and the map context menu can trigger it.
      auto plotRouteTo = [&](sim::SystemId dstId) {
        if (dstId == 0) return;

        const double jrMaxLy = fsdBaseRangeLy();
        const double jrNowLy = fsdCurrentRangeLy();
        const double planMaxLy = navConstrainToCurrentFuelRange ? jrNowLy : jrMaxLy;

        if (planMaxLy <= 0.0) {
          toast(toasts, "No jump range available (refuel or reduce cargo).", 2.6);
          return;
        }

        sim::RoutePlanStats stats{};
        double costPerJump = 1.0;
        double costPerLy = 0.0;

        if (navRouteMode == NavRouteMode::Hops) {
          costPerJump = 1.0;
          costPerLy = 0.0;
          navRoute = sim::plotRouteAStarHops(nearby, currentSystem->stub.id, dstId, planMaxLy, &stats);
        } else if (navRouteMode == NavRouteMode::Distance) {
          costPerJump = 0.0;
          costPerLy = 1.0;
          navRoute = sim::plotRouteAStarCost(nearby, currentSystem->stub.id, dstId, planMaxLy, costPerJump, costPerLy, &stats);
        } else {
          costPerJump = kFsdFuelBase;
          costPerLy = kFsdFuelPerLy;
          navRoute = sim::plotRouteAStarCost(nearby, currentSystem->stub.id, dstId, planMaxLy, costPerJump, costPerLy, &stats);
        }

        navRouteHop = 0;
        navAutoRun = false;

        if (navRoute.empty()) {
          navRoutePlanStatsValid = false;
          toast(toasts, "No route found. Try increasing the scan radius or switching route mode.", 3.0);
          return;
        }

        navRoutePlanStats = stats;
        navRoutePlanStatsValid = true;
        navRoutePlannedMode = navRouteMode;
        navRoutePlannedUsedCurrentFuelRange = navConstrainToCurrentFuelRange;
        navRoutePlanMaxJumpLy = planMaxLy;
        navRoutePlanCostPerJump = costPerJump;
        navRoutePlanCostPerLy = costPerLy;

        const double totalDist = sim::routeDistanceLy(nearby, navRoute);
        const double totalFuel = sim::routeCost(nearby, navRoute, kFsdFuelBase, kFsdFuelPerLy);
        toast(toasts,
              "Route plotted (" + std::to_string((int)navRoute.size() - 1)
              + " jumps, " + std::to_string((int)std::round(totalDist))
              + " ly, est fuel " + std::to_string((int)std::round(totalFuel)) + ").",
              2.6);
      };

      // Star-class tinting for the map.
      auto starTint = [&](sim::StarClass cls, float alpha) -> ImU32 {
        alpha = std::clamp(alpha, 0.0f, 1.0f);
        auto pack = [&](int r, int g, int b) -> ImU32 {
          return IM_COL32(r, g, b, (int)(alpha * 255.0f));
        };
        switch (cls) {
          case sim::StarClass::O: return pack(120, 180, 255);
          case sim::StarClass::B: return pack(155, 205, 255);
          case sim::StarClass::A: return pack(240, 240, 255);
          case sim::StarClass::F: return pack(255, 248, 220);
          case sim::StarClass::G: return pack(255, 236, 175);
          case sim::StarClass::K: return pack(255, 200, 140);
          case sim::StarClass::M: return pack(255, 165, 140);
          default: return pack(200, 200, 220);
        }
      };

      // Persistent "infinite canvas" style map view state.
      static sim::SystemId mapLastSystemId = 0;
      static math::Vec3d mapViewCenterLy{0,0,0};
      static float mapZoom = 1.0f;
      static float mapHeight = 520.0f;

      static bool mapShowLanes = true;
      static bool mapShowGrid = true;
      static bool mapShowLabels = false;
      static bool mapShowFactions = true;
      static bool mapDepthFade = true;
      static bool mapShowJumpRings = true;
      static bool mapShowRoute = true;

      if (mapLastSystemId != currentSystem->stub.id) {
        mapLastSystemId = currentSystem->stub.id;
        mapViewCenterLy = currentSystem->stub.posLy;
      }

      if (ImGui::BeginTable("galaxy_layout", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV)) {
        ImGui::TableSetupColumn("Map", ImGuiTableColumnFlags_WidthStretch, 0.62f);
        ImGui::TableSetupColumn("Details", ImGuiTableColumnFlags_WidthStretch, 0.38f);
        ImGui::TableNextRow();

        // ---------------- Left column (map + list) ----------------
        ImGui::TableSetColumnIndex(0);

        if (ImGui::CollapsingHeader("Map options", ImGuiTreeNodeFlags_DefaultOpen)) {
          ImGui::Checkbox("Star lanes", &mapShowLanes); ImGui::SameLine();
          ImGui::Checkbox("Grid", &mapShowGrid); ImGui::SameLine();
          ImGui::Checkbox("Labels", &mapShowLabels);

          ImGui::Checkbox("Factions", &mapShowFactions); ImGui::SameLine();
          ImGui::Checkbox("Depth fade (Z)", &mapDepthFade); ImGui::SameLine();
          ImGui::Checkbox("Jump rings", &mapShowJumpRings);

          ImGui::Checkbox("Route overlay", &mapShowRoute);

          ImGui::SetNextItemWidth(180.0f);
          ImGui::SliderFloat("Zoom", &mapZoom, 0.25f, 8.0f, "%.2fx");

          ImGui::SameLine();
          if (ImGui::Button("Center current")) {
            mapViewCenterLy = currentSystem->stub.posLy;
          }
          ImGui::SameLine();
          if (ImGui::Button("Center selection") && galaxySelectedSystemId != 0) {
            mapViewCenterLy = universe.getSystem(galaxySelectedSystemId).stub.posLy;
          }

          ImGui::SetNextItemWidth(220.0f);
          ImGui::SliderFloat("Map height", &mapHeight, 260.0f, 900.0f, "%.0f px");

          ImGui::TextDisabled("Wheel: zoom | MMB/RMB drag: pan | LMB: select | RMB: context");
        }

        // Interactive galaxy map canvas.
        {
          const ImVec2 childSize(0, mapHeight);
          ImGui::BeginChild("map", childSize, false, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);

          ImDrawList* draw = ImGui::GetWindowDrawList();
          const ImVec2 canvasP0 = ImGui::GetCursorScreenPos();
          ImVec2 canvasSz = ImGui::GetContentRegionAvail();
          canvasSz.x = std::max(canvasSz.x, 60.0f);
          canvasSz.y = std::max(canvasSz.y, 60.0f);

          const ImVec2 canvasP1(canvasP0.x + canvasSz.x, canvasP0.y + canvasSz.y);
          const ImVec2 canvasCenter((canvasP0.x + canvasP1.x) * 0.5f, (canvasP0.y + canvasP1.y) * 0.5f);

          // Background
          draw->AddRectFilled(canvasP0, canvasP1, IM_COL32(10, 10, 14, 255));
          draw->AddRect(canvasP0, canvasP1, IM_COL32(80, 80, 95, 255));

          // Interaction surface
          ImGui::InvisibleButton("galaxy_map_btn", canvasSz,
                                 ImGuiButtonFlags_MouseButtonLeft |
                                 ImGuiButtonFlags_MouseButtonRight |
                                 ImGuiButtonFlags_MouseButtonMiddle);
          const bool isHovered = ImGui::IsItemHovered();
          const bool isActive = ImGui::IsItemActive();
          const ImVec2 mouse = ImGui::GetIO().MousePos;

          const float halfMin = 0.5f * std::min(canvasSz.x, canvasSz.y);
          float pxPerLy = std::max(halfMin / std::max(1.0f, radius) * mapZoom, 1e-6f);

          // Zoom (mouse wheel) - keep the world point under the cursor stable.
          if (isHovered) {
            const float wheel = ImGui::GetIO().MouseWheel;
            if (wheel != 0.0f) {
              const float zoomFactor = std::pow(1.2f, wheel);
              const float newZoom = std::clamp(mapZoom * zoomFactor, 0.25f, 12.0f);

              const ImVec2 rel(mouse.x - canvasCenter.x, mouse.y - canvasCenter.y);
              const math::Vec3d worldUnderMouse = mapViewCenterLy + math::Vec3d((double)rel.x / (double)pxPerLy,
                                                                               (double)rel.y / (double)pxPerLy, 0.0);

              const float newPxPerLy = std::max(halfMin / std::max(1.0f, radius) * newZoom, 1e-6f);

              mapViewCenterLy = worldUnderMouse - math::Vec3d((double)rel.x / (double)newPxPerLy,
                                                             (double)rel.y / (double)newPxPerLy, 0.0);
              mapZoom = newZoom;
              pxPerLy = newPxPerLy;
            }
          }

          // Pan (drag with middle or right mouse button).
          if (isActive && (ImGui::IsMouseDragging(ImGuiMouseButton_Middle, 0.0f) ||
                           ImGui::IsMouseDragging(ImGuiMouseButton_Right, 0.0f))) {
            const ImVec2 d = ImGui::GetIO().MouseDelta;
            mapViewCenterLy.x -= (double)d.x / (double)pxPerLy;
            mapViewCenterLy.y -= (double)d.y / (double)pxPerLy;
          }

          auto toPx = [&](const math::Vec3d& posLy) -> ImVec2 {
            const math::Vec3d d = posLy - mapViewCenterLy;
            return ImVec2(canvasCenter.x + (float)(d.x * (double)pxPerLy),
                          canvasCenter.y + (float)(d.y * (double)pxPerLy));
          };

          draw->PushClipRect(canvasP0, canvasP1, true);

          // Grid (nice step size) in view space.
          if (mapShowGrid) {
            const double halfWLy = (double)(canvasSz.x * 0.5f) / (double)pxPerLy;
            const double halfHLy = (double)(canvasSz.y * 0.5f) / (double)pxPerLy;
            const double minX = mapViewCenterLy.x - halfWLy;
            const double maxX = mapViewCenterLy.x + halfWLy;
            const double minY = mapViewCenterLy.y - halfHLy;
            const double maxY = mapViewCenterLy.y + halfHLy;

            const double targetPx = 70.0;
            const double rawStep = targetPx / (double)pxPerLy;
            const double p10 = std::pow(10.0, std::floor(std::log10(std::max(1e-6, rawStep))));
            const double m = rawStep / p10;
            const double step = (m < 2.0) ? 2.0 * p10 : (m < 5.0) ? 5.0 * p10 : 10.0 * p10;

            const double gx0 = std::floor(minX / step) * step;
            const double gy0 = std::floor(minY / step) * step;

            const ImU32 gridCol = IM_COL32(28, 30, 40, 140);
            for (double gx = gx0; gx <= maxX + 1e-6; gx += step) {
              const float x = canvasCenter.x + (float)((gx - mapViewCenterLy.x) * (double)pxPerLy);
              draw->AddLine(ImVec2(x, canvasP0.y), ImVec2(x, canvasP1.y), gridCol, 1.0f);
            }
            for (double gy = gy0; gy <= maxY + 1e-6; gy += step) {
              const float y = canvasCenter.y + (float)((gy - mapViewCenterLy.y) * (double)pxPerLy);
              draw->AddLine(ImVec2(canvasP0.x, y), ImVec2(canvasP1.x, y), gridCol, 1.0f);
            }

            // Center crosshair (view-center).
            draw->AddLine(ImVec2(canvasCenter.x - 8, canvasCenter.y), ImVec2(canvasCenter.x + 8, canvasCenter.y), IM_COL32(60, 60, 80, 200), 1.0f);
            draw->AddLine(ImVec2(canvasCenter.x, canvasCenter.y - 8), ImVec2(canvasCenter.x, canvasCenter.y + 8), IM_COL32(60, 60, 80, 200), 1.0f);
          }

          // Build a lookup of stub positions for route drawing.
          std::unordered_map<sim::SystemId, math::Vec3d> stubPosById;
          stubPosById.reserve(nearby.size());
          for (const auto& s : nearby) {
            stubPosById[s.id] = s.posLy;
          }

          // Jump range feedback rings around the CURRENT system.
          const double jrMaxLy = fsdBaseRangeLy();
          const double jrNowLy = fsdCurrentRangeLy();

          if (mapShowJumpRings) {
            const float jrMaxPx = (float)(jrMaxLy * (double)pxPerLy);
            const float jrNowPx = (float)(jrNowLy * (double)pxPerLy);
            const ImVec2 curPx = toPx(currentSystem->stub.posLy);

            if (jrMaxPx > 2.0f && jrMaxPx < 50000.0f) {
              draw->AddCircle(curPx, jrMaxPx, IM_COL32(90, 150, 240, 120), 96, 1.5f);
            }
            if (jrNowPx > 2.0f && jrNowPx < 50000.0f) {
              draw->AddCircle(curPx, jrNowPx, IM_COL32(90, 220, 150, 120), 96, 1.5f);
            }
          }

          // Plotted route overlay.
          if (mapShowRoute && navRoute.size() >= 2) {
            for (std::size_t i = 0; i + 1 < navRoute.size(); ++i) {
              auto itA = stubPosById.find(navRoute[i]);
              auto itB = stubPosById.find(navRoute[i + 1]);
              if (itA == stubPosById.end() || itB == stubPosById.end()) continue;
              const ImVec2 a = toPx(itA->second);
              const ImVec2 b = toPx(itB->second);
              draw->AddLine(a, b, IM_COL32(255, 140, 80, 210), 2.5f);
            }
          }

          // Star lanes: connect each system to 3 nearest neighbors in 3D (draw projected in XY).
          if (mapShowLanes) {
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
                draw->AddLine(a, b, IM_COL32(50, 80, 120, 110), 1.0f);
              }
            }
          }

          // Hover detection for click + tooltip.
          sim::SystemId hoveredId = 0;
          float hoveredD2 = 1.0e30f;
          if (isHovered) {
            for (const auto& s : nearby) {
              const ImVec2 p = toPx(s.posLy);
              const bool isCurrent = (s.id == currentSystem->stub.id);
              const float iconSz = isCurrent ? 26.0f : 18.0f;
              const float r = iconSz * 0.6f;
              const float dx = mouse.x - p.x;
              const float dy = mouse.y - p.y;
              const float d2 = dx*dx + dy*dy;
              if (d2 <= r*r && d2 < hoveredD2) {
                hoveredD2 = d2;
                hoveredId = s.id;
              }
            }
          }

          // Click to select.
          if (isHovered && hoveredId != 0 && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            galaxySelectedSystemId = hoveredId;
          }

          // Double-click to center view on the hovered system.
          if (isHovered && hoveredId != 0 && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
            mapViewCenterLy = universe.getSystem(hoveredId).stub.posLy;
          }

          // Context menu (right click release without drag).
          static sim::SystemId mapContextId = 0;
          if (isHovered && ImGui::IsMouseReleased(ImGuiMouseButton_Right) &&
              !ImGui::IsMouseDragging(ImGuiMouseButton_Right, 4.0f)) {
            mapContextId = hoveredId;
            ImGui::OpenPopup("galaxy_map_context");
          }

          if (ImGui::BeginPopup("galaxy_map_context")) {
            if (mapContextId != 0) {
              const auto& sys = universe.getSystem(mapContextId);
              ImGui::TextUnformatted(sys.stub.name.c_str());
              ImGui::Separator();

              if (ImGui::MenuItem("Select")) {
                galaxySelectedSystemId = mapContextId;
              }
              if (ImGui::MenuItem("Center view")) {
                mapViewCenterLy = sys.stub.posLy;
              }

              ImGui::SeparatorText("Bookmarks");
              {
                const bool hasSysBm = ui::findSystemBookmarkIndex(bookmarks, mapContextId).has_value();
                if (ImGui::MenuItem(hasSysBm ? "Unbookmark system" : "Bookmark system")) {
                  if (hasSysBm) {
                    if (ui::removeSystemBookmark(bookmarks, mapContextId)) {
                      bookmarksDirty = true;
                      toast(toasts, "Removed bookmark.", 1.4);
                    }
                  } else {
                    ui::addOrUpdateSystemBookmark(bookmarks, mapContextId, sys.stub.name);
                    bookmarksDirty = true;
                    toast(toasts, "Bookmarked system.", 1.4);
                  }
                }

                if (ImGui::BeginMenu("Stations")) {
                  for (const auto& st : sys.stations) {
                    const bool hasSt = ui::findStationBookmarkIndex(bookmarks, mapContextId, st.id).has_value();
                    std::string label = st.name;
                    if (hasSt) label += "  (bookmarked)";
                    if (ImGui::MenuItem(label.c_str())) {
                      if (hasSt) {
                        if (ui::removeStationBookmark(bookmarks, mapContextId, st.id)) {
                          bookmarksDirty = true;
                          toast(toasts, "Removed station bookmark.", 1.4);
                        }
                      } else {
                        ui::addOrUpdateStationBookmark(bookmarks, mapContextId, st.id, sys.stub.name + " / " + st.name);
                        bookmarksDirty = true;
                        toast(toasts, "Bookmarked station.", 1.4);
                      }
                    }
                  }
                  ImGui::EndMenu();
                }

                if (ui::hasAnyInSystem(bookmarks, mapContextId)) {
                  if (ImGui::MenuItem("Remove all in system")) {
                    if (ui::removeAllInSystem(bookmarks, mapContextId)) {
                      bookmarksDirty = true;
                      toast(toasts, "Removed all bookmarks in system.", 1.6);
                    }
                  }
                }
              }

              ImGui::Separator();

              if (ImGui::MenuItem("Plot route")) {
                plotRouteTo(mapContextId);
              }
              if (ImGui::MenuItem("Engage FSD (direct)")) {
                startFsdJumpTo(mapContextId);
              }
            } else {
              ImGui::TextUnformatted("Map");
            }

            ImGui::Separator();
            if (ImGui::MenuItem("Center current")) {
              mapViewCenterLy = currentSystem->stub.posLy;
            }
            if (ImGui::MenuItem("Reset zoom")) {
              mapZoom = 1.0f;
            }
            ImGui::EndPopup();
          }

          // Draw systems using the HUD atlas star icon (tinted by star class).
          const auto& atlasTex = hudAtlas.texture();
          const ImTextureID atlasId = (ImTextureID)(intptr_t)atlasTex.handle();

          const sim::SystemId nextHopId = (!navRoute.empty() && navRouteHop + 1 < navRoute.size()) ? navRoute[navRouteHop + 1] : 0;
          const float pulse = 0.55f + 0.45f * (float)std::sin((float)timeRealSec * 3.0f);

          for (const auto& s : nearby) {
            const ImVec2 p = toPx(s.posLy);
            const bool isCurrent = (s.id == currentSystem->stub.id);
            const bool isSel = (s.id == galaxySelectedSystemId);
            const bool isHover = (s.id == hoveredId);
            const bool isNextHop = (s.id == nextHopId);

            const double distLy = (s.posLy - currentSystem->stub.posLy).length();
            const bool inJumpRange = (distLy <= jrMaxLy + 1e-9);

            float alpha = 1.0f;
            if (mapDepthFade) {
              const double dz = s.posLy.z - center.z;
              const double zn = std::min(1.0, std::abs(dz) / std::max(1.0, (double)radius));
              alpha *= 1.0f - 0.55f * (float)zn;
            }

            float iconSz = isCurrent ? 26.0f : 18.0f;
            if (isSel || isHover) iconSz = std::max(iconSz, 22.0f);

            // Star icon
            const auto uv = hudAtlas.get(render::SpriteKind::Star, core::hashCombine(core::fnv1a64("galaxy_star"), s.seed));
            const ImVec2 i0(p.x - iconSz * 0.5f, p.y - iconSz * 0.5f);
            const ImVec2 i1(p.x + iconSz * 0.5f, p.y + iconSz * 0.5f);
            draw->AddImage(atlasId, i0, i1, ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1), starTint(s.primaryClass, alpha));

            // Bookmark marker: a small warm dot offset from the system icon.
            if (ui::hasAnyInSystem(bookmarks, s.id)) {
              const float rB = std::max(2.5f, iconSz * 0.18f);
              const ImVec2 bp(p.x - iconSz * 0.55f, p.y - iconSz * 0.55f);
              draw->AddCircleFilled(bp, rB, IM_COL32(255, 210, 90, (int)(alpha * 220.0f)), 12);
              draw->AddCircle(bp, rB, IM_COL32(20, 18, 10, (int)(alpha * 180.0f)), 12, 1.0f);
            }

            // Faction marker (small overlay glyph)
            if (mapShowFactions && s.factionId != 0) {
              const float fSz = std::max(8.0f, iconSz * 0.55f);
              const auto uvF = hudAtlas.get(render::SpriteKind::Faction, core::hashCombine(core::fnv1a64("faction"), (core::u64)s.factionId));
              const ImVec2 f0(p.x + iconSz * 0.12f, p.y - iconSz * 0.55f);
              const ImVec2 f1(f0.x + fSz, f0.y + fSz);
              draw->AddImage(atlasId, f0, f1, ImVec2(uvF.u0, uvF.v0), ImVec2(uvF.u1, uvF.v1), IM_COL32(200, 255, 210, (int)(alpha * 235.0f)));
            }

            // Jump range ring
            if (!isCurrent && inJumpRange) {
              draw->AddCircle(p, iconSz * 0.62f, IM_COL32(80, 140, 220, (int)(alpha * 150.0f)), 24, 1.0f);
            }

            // Selected / hovered highlight
            if (isSel) {
              draw->AddCircle(p, iconSz * 0.82f, IM_COL32(255, 110, 110, 220), 32, 2.0f);
            } else if (isHover) {
              draw->AddCircle(p, iconSz * 0.82f, IM_COL32(255, 255, 255, 160), 32, 2.0f);
            } else if (isCurrent) {
              draw->AddCircle(p, iconSz * 0.82f, IM_COL32(255, 240, 160, 210), 32, 2.0f);
            }

            // Next hop pulse highlight
            if (isNextHop) {
              draw->AddCircle(p, iconSz * (0.95f + 0.12f * pulse), IM_COL32(255, 170, 90, (int)(170.0f * pulse)), 48, 2.0f);
            }

            // Labels
            if (mapShowLabels || isCurrent || isSel || isHover) {
              const std::string label = s.name;
              draw->AddText(ImVec2(p.x + iconSz * 0.65f, p.y - iconSz * 0.35f),
                            IM_COL32(220, 220, 235, (int)(alpha * 200.0f)), label.c_str());
            }
          }

          // Tooltip for hovered system.
          if (isHovered && hoveredId != 0) {
            const auto& sys = universe.getSystem(hoveredId);
            const double distLy = (sys.stub.posLy - currentSystem->stub.posLy).length();
            const double fuelCost = fsdFuelCostFor(distLy);
            const bool inRange = (distLy <= jrMaxLy + 1e-6);
            const bool hasFuel = (fuel + 1e-6 >= fuelCost);

            ImGui::BeginTooltip();

            // Icon + name line
            {
              const auto uv = hudAtlas.get(render::SpriteKind::Star, core::hashCombine(core::fnv1a64("galaxy_star"), sys.stub.seed));
              ImGui::Image((ImTextureID)(intptr_t)hudAtlas.texture().handle(), ImVec2(22, 22),
                           ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1),
                           ImVec4(1, 1, 1, 1));
              ImGui::SameLine();
              ImGui::Text("%s", sys.stub.name.c_str());
            }

            ImGui::Text("Class: %s", starClassName(sys.stub.primaryClass));
            ImGui::Text("Distance: %.1f ly", distLy);
            ImGui::Text("Stations: %d  |  Planets: %d", sys.stub.stationCount, sys.stub.planetCount);
            ImGui::Text("Faction: %s", factionName(sys.stub.factionId).c_str());
            ImGui::Separator();
            if (!inRange) {
              ImGui::TextColored(ImVec4(1.0f, 0.45f, 0.45f, 1.0f), "Direct jump: OUT OF RANGE");
            } else if (!hasFuel) {
              ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.35f, 1.0f), "Direct jump: IN RANGE (LOW FUEL)");
            } else {
              ImGui::TextColored(ImVec4(0.70f, 0.95f, 0.70f, 1.0f), "Direct jump: IN RANGE");
            }
            ImGui::TextDisabled("Fuel cost %.1f | fuel %.1f", fuelCost, fuel);
            ImGui::TextDisabled("LMB select | Double-click center | RMB menu");
            ImGui::EndTooltip();
          }

          // Map legend (top-left).
          draw->AddText(ImVec2(canvasP0.x + 8.0f, canvasP0.y + 8.0f),
                        IM_COL32(170, 170, 190, 170),
                        "Galaxy Map");

          draw->PopClipRect();
          ImGui::EndChild();
        }

        // System list (search + sort) for fast selection (useful once the map gets dense).
        if (ImGui::CollapsingHeader("System list", ImGuiTreeNodeFlags_DefaultOpen)) {
          static char sysFilter[64] = "";
          ImGui::InputText("Filter", sysFilter, sizeof(sysFilter));
          ImGui::SameLine();
          static int sortMode = 0;
          const char* sortNames[] = {"Distance", "Name"};
          ImGui::Combo("Sort", &sortMode, sortNames, 2);

          auto toLower = [](std::string s) {
            for (auto& ch : s) ch = (char)std::tolower((unsigned char)ch);
            return s;
          };

          const std::string f = toLower(std::string(sysFilter));

          struct Row {
            const sim::SystemStub* s;
            double dist;
          };

          std::vector<Row> rows;
          rows.reserve(nearby.size());

          for (const auto& s : nearby) {
            if (!f.empty()) {
              const std::string n = toLower(s.name);
              if (n.find(f) == std::string::npos) continue;
            }
            const double d = (s.posLy - currentSystem->stub.posLy).length();
            rows.push_back(Row{&s, d});
          }

          if (sortMode == 0) {
            std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.dist < b.dist; });
          } else {
            std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.s->name < b.s->name; });
          }

          ImGui::BeginChild("system_list", ImVec2(0, 240), true);
          if (ImGui::BeginTable("sys_table", 5, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_SizingFixedFit)) {
            ImGui::TableSetupColumn("Name");
            ImGui::TableSetupColumn("Dist (ly)");
            ImGui::TableSetupColumn("St");
            ImGui::TableSetupColumn("Faction");
            ImGui::TableSetupColumn("Range");
            ImGui::TableHeadersRow();

            for (const auto& r : rows) {
              ImGui::TableNextRow();
              ImGui::TableSetColumnIndex(0);
              const bool isSel = (r.s->id == galaxySelectedSystemId);
              std::string rowLabel = r.s->name;
              if (ui::hasAnyInSystem(bookmarks, r.s->id)) rowLabel += "  [B]";
              if (ImGui::Selectable(rowLabel.c_str(), isSel)) {
                galaxySelectedSystemId = r.s->id;
              }

              ImGui::TableSetColumnIndex(1);
              ImGui::Text("%.1f", r.dist);

              ImGui::TableSetColumnIndex(2);
              ImGui::Text("%d", r.s->stationCount);

              ImGui::TableSetColumnIndex(3);
              ImGui::TextUnformatted(factionName(r.s->factionId).c_str());

              ImGui::TableSetColumnIndex(4);
              const bool inRange = (r.s->id == currentSystem->stub.id) || (r.dist <= galaxyPlanMaxLy + 1e-6);
              if (inRange) ImGui::TextColored(ImVec4(0.70f, 0.95f, 0.70f, 1.0f), "IN");
              else ImGui::TextColored(ImVec4(1.0f, 0.45f, 0.45f, 1.0f), "OUT");
            }

            ImGui::EndTable();
          }
          ImGui::EndChild();
        }

        // ---------------- Right column (details + planner) ----------------
        ImGui::TableSetColumnIndex(1);

        ImGui::SeparatorText("Selection");

        if (galaxySelectedSystemId == 0) {
          ImGui::Text("Select a system on the map.");
        } else {
          const auto& selSys = universe.getSystem(galaxySelectedSystemId);

          // Icon + title row.
          {
            const auto uv = hudAtlas.get(render::SpriteKind::Star, core::hashCombine(core::fnv1a64("galaxy_star"), selSys.stub.seed));
            ImGui::Image((ImTextureID)(intptr_t)hudAtlas.texture().handle(), ImVec2(28, 28),
                         ImVec2(uv.u0, uv.v0), ImVec2(uv.u1, uv.v1), ImVec4(1, 1, 1, 1));
            ImGui::SameLine();
            ImGui::Text("%s", selSys.stub.name.c_str());
          }

          const double distLy = (selSys.stub.posLy - currentSystem->stub.posLy).length();
          const double jrMaxLySel = fsdBaseRangeLy();
          const double jrNowLySel = fsdCurrentRangeLy();

          ImGui::Text("Class: %s", starClassName(selSys.stub.primaryClass));
          ImGui::Text("Distance: %.1f ly", distLy);
          ImGui::Text("Stations: %d  |  Planets: %d", selSys.stub.stationCount, selSys.stub.planetCount);
          ImGui::Text("Faction: %s", factionName(selSys.stub.factionId).c_str());

          ImGui::SeparatorText("Bookmarks");
          {
            const bool hasSysBm = ui::findSystemBookmarkIndex(bookmarks, galaxySelectedSystemId).has_value();
            if (ImGui::Button(hasSysBm ? "Unbookmark system" : "Bookmark system")) {
              if (hasSysBm) {
                if (ui::removeSystemBookmark(bookmarks, galaxySelectedSystemId)) {
                  bookmarksDirty = true;
                  toast(toasts, "Removed bookmark.", 1.4);
                }
              } else {
                ui::addOrUpdateSystemBookmark(bookmarks, galaxySelectedSystemId, selSys.stub.name);
                bookmarksDirty = true;
                toast(toasts, "Bookmarked system.", 1.4);
              }
            }
            ImGui::SameLine();
            if (ImGui::Button("Open Bookmarks")) {
              showBookmarksWindow = true;
            }

            if (!selSys.stations.empty()) {
              // Station bookmark quick-toggle.
              static sim::SystemId lastSysForBm = 0;
              static int stationPick = 0;
              if (lastSysForBm != galaxySelectedSystemId) {
                lastSysForBm = galaxySelectedSystemId;
                stationPick = 0;
              }

              if (stationPick < 0) stationPick = 0;
              if (stationPick >= (int)selSys.stations.size()) stationPick = (int)selSys.stations.size() - 1;

              ImGui::SetNextItemWidth(-FLT_MIN);
              const char* preview = selSys.stations[(std::size_t)stationPick].name.c_str();
              if (ImGui::BeginCombo("Station", preview)) {
                for (int i = 0; i < (int)selSys.stations.size(); ++i) {
                  const auto& st = selSys.stations[(std::size_t)i];
                  const bool isSel = (stationPick == i);
                  const bool isBm = ui::findStationBookmarkIndex(bookmarks, galaxySelectedSystemId, st.id).has_value();
                  std::string label = st.name;
                  if (isBm) label += "  (bookmarked)";
                  if (ImGui::Selectable(label.c_str(), isSel)) stationPick = i;
                }
                ImGui::EndCombo();
              }

              const auto& st = selSys.stations[(std::size_t)stationPick];
              const bool hasSt = ui::findStationBookmarkIndex(bookmarks, galaxySelectedSystemId, st.id).has_value();
              if (ImGui::Button(hasSt ? "Unbookmark station" : "Bookmark station")) {
                if (hasSt) {
                  if (ui::removeStationBookmark(bookmarks, galaxySelectedSystemId, st.id)) {
                    bookmarksDirty = true;
                    toast(toasts, "Removed station bookmark.", 1.4);
                  }
                } else {
                  ui::addOrUpdateStationBookmark(bookmarks, galaxySelectedSystemId, st.id, selSys.stub.name + " / " + st.name);
                  bookmarksDirty = true;
                  toast(toasts, "Bookmarked station.", 1.4);
                }
              }
            }
          }

          ImGui::Separator();

          ImGui::Text("Jump range: %.1f ly max, %.1f ly current-fuel", jrMaxLySel, jrNowLySel);

          if (galaxySelectedSystemId == currentSystem->stub.id) {
            ImGui::TextDisabled("This is your current system.");
          } else {
            const double fuelCost = fsdFuelCostFor(distLy);
            const bool inRange = (distLy <= jrMaxLySel + 1e-6);
            const bool hasFuel = (fuel + 1e-6 >= fuelCost);
            const char* rangeText = !inRange ? "OUT OF RANGE" : (hasFuel ? "IN RANGE" : "IN RANGE (LOW FUEL)");
            const ImVec4 rangeCol = !inRange ? ImVec4(1.0f, 0.45f, 0.45f, 1.0f)
                                            : (hasFuel ? ImVec4(0.70f, 0.95f, 0.70f, 1.0f)
                                                      : ImVec4(1.0f, 0.75f, 0.35f, 1.0f));

            ImGui::TextColored(rangeCol, "Direct jump: %s", rangeText);
            ImGui::SameLine();
            ImGui::TextDisabled("(fuel cost %.1f, fuel %.1f)", fuelCost, fuel);

            ImGui::SeparatorText("Route planner");

            {
              static const char* kRouteModes[] = {"Min jumps", "Min distance", "Min fuel"};
              int modeIdx = (int)navRouteMode;
              if (ImGui::Combo("Mode", &modeIdx, kRouteModes, (int)(sizeof(kRouteModes) / sizeof(kRouteModes[0])))) {
                navRouteMode = (NavRouteMode)modeIdx;
              }
            }

            ImGui::Checkbox("Constrain edges to current-fuel range", &navConstrainToCurrentFuelRange);
            const double planMaxLy = navConstrainToCurrentFuelRange ? jrNowLySel : jrMaxLySel;
            ImGui::TextDisabled("Planner max hop: %.1f ly", planMaxLy);

            if (ImGui::Button("Plot route")) {
              plotRouteTo(galaxySelectedSystemId);
            }
            ImGui::SameLine();
            if (ImGui::Button("Clear route")) {
              navRoute.clear();
              navRouteHop = 0;
              navAutoRun = false;
              navRoutePlanStatsValid = false;
            }

            ImGui::Checkbox("Auto-run route", &navAutoRun);

            if (!navRoute.empty()) {
              const int totalJumps = (int)navRoute.size() - 1;
              const int remaining = std::max(0, totalJumps - (int)navRouteHop);

              double totalDist = 0.0;
              double totalFuel = 0.0;
              double maxHopDist = 0.0;
              for (std::size_t i = 0; i + 1 < navRoute.size(); ++i) {
                const auto& a = universe.getSystem(navRoute[i]).stub;
                const auto& b = universe.getSystem(navRoute[i + 1]).stub;
                const double d = (b.posLy - a.posLy).length();
                totalDist += d;
                totalFuel += fsdFuelCostFor(d);
                if (d > maxHopDist) maxHopDist = d;
              }

              ImGui::Text("Route: %d jumps (%d remaining)", totalJumps, remaining);
              ImGui::Text("Total: %.1f ly, est fuel %.1f (no refuel)", totalDist, totalFuel);
              ImGui::TextDisabled("Max hop: %.1f ly", maxHopDist);

              if (totalFuel > fuel + 1e-6) {
                ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.35f, 1.0f), "Warning: route needs %.1f fuel, you have %.1f.", totalFuel, fuel);
              }

              if (navRoutePlanStatsValid) {
                const char* planName = (navRoutePlannedMode == NavRouteMode::Hops) ? "Min jumps"
                                     : (navRoutePlannedMode == NavRouteMode::Distance) ? "Min distance"
                                     : "Min fuel";
                ImGui::TextDisabled("Planner: %s | max hop %.1f ly | expansions %d", planName, navRoutePlanMaxJumpLy, navRoutePlanStats.expansions);
                ImGui::TextDisabled("Cost model: %.2f/jump + %.2f/ly", navRoutePlanCostPerJump, navRoutePlanCostPerLy);
                ImGui::TextDisabled("Constraint: %s", navRoutePlannedUsedCurrentFuelRange ? "current-fuel range" : "max range");
              }

              if (navRouteHop + 1 < navRoute.size()) {
                const auto& nextSys = universe.getSystem(navRoute[navRouteHop + 1]);
                const double hopDist = (nextSys.stub.posLy - currentSystem->stub.posLy).length();
                const double hopFuel = fsdFuelCostFor(hopDist);
                if (fuel + 1e-6 < hopFuel) {
                  ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.35f, 1.0f),
                                     "Next hop: %s (%.1f ly, fuel %.1f > %.1f)",
                                     nextSys.stub.name.c_str(), hopDist, hopFuel, fuel);
                } else {
                  ImGui::Text("Next hop: %s (%.1f ly, fuel %.1f)", nextSys.stub.name.c_str(), hopDist, hopFuel);
                }
              }
            }

            const std::string jumpLabel =
                std::string("Engage FSD (") + game::chordLabel(controls.actions.fsdJump) + ")";
            if (ImGui::Button(jumpLabel.c_str())) {
              if (!navRoute.empty() && navRouteHop + 1 < navRoute.size()) {
                startFsdJumpTo(navRoute[navRouteHop + 1]);
              } else {
                startFsdJumpTo(galaxySelectedSystemId);
              }
            }
          }
        }

        ImGui::EndTable();
      }

      {
        std::string hint;
        hint.reserve(160);
        hint += game::chordLabel(controls.actions.toggleGalaxy);
        hint += " toggles this window, ";
        hint += game::chordLabel(controls.actions.toggleShip);
        hint += " Ship, ";
        hint += game::chordLabel(controls.actions.toggleMarket);
        hint += " Market, ";
        hint += game::chordLabel(controls.actions.toggleContacts);
        hint += " Contacts, ";
        hint += game::chordLabel(controls.actions.togglePostFx);
        hint += " PostFX";
        ImGui::TextDisabled("%s", hint.c_str());
      }
      ImGui::End();
    }

    // ---- Bookmarks ----
    if (showBookmarksWindow) {
      ImGui::SetNextWindowSize(ImVec2(760, 520), ImGuiCond_FirstUseEver);
      if (ImGui::Begin("Bookmarks", &showBookmarksWindow)) {
        ImGui::TextDisabled("Saved to %s", bookmarksPath.c_str());
        ImGui::SameLine();
        if (bookmarksDirty) {
          ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.35f, 1.0f), "(unsaved changes)");
        }

        if (ImGui::Button("Save")) {
          const bool ok = ui::saveBookmarksToFile(bookmarks, bookmarksPath);
          if (ok) {
            bookmarksDirty = false;
            toast(toasts, "Saved bookmarks.", 1.4);
          } else {
            toast(toasts, "Failed to save bookmarks.", 1.8);
          }
        }
        ImGui::SameLine();
        if (ImGui::Button("Load")) {
          ui::Bookmarks loaded = ui::makeDefaultBookmarks();
          if (ui::loadBookmarksFromFile(bookmarksPath, loaded)) {
            bookmarks = std::move(loaded);
            bookmarksDirty = false;
            toast(toasts, "Loaded bookmarks.", 1.4);
          } else {
            toast(toasts, "No bookmarks file found.", 1.6);
          }
        }
        ImGui::SameLine();
        if (ImGui::Button("Clear all")) {
          bookmarks.items.clear();
          bookmarksDirty = true;
          toast(toasts, "Cleared all bookmarks.", 1.4);
        }

        ImGui::SameLine();
        ImGui::Checkbox("Auto-save on exit", &bookmarksAutoSaveOnExit);

        ImGui::Separator();

        static char bmFilter[96] = "";
        ImGui::InputTextWithHint("##bm_filter", "Filter (name)...", bmFilter, sizeof(bmFilter));

        auto anyNonSpaceLocal = [&](const char* s) {
          if (!s) return false;
          for (const unsigned char* p = (const unsigned char*)s; *p; ++p) {
            if (std::isspace(*p) == 0) return true;
          }
          return false;
        };

        auto drawHighlighted = [&](std::string_view text, const std::vector<int>& positions, ImU32 colText) {
          if (positions.empty()) {
            ImGui::TextUnformatted(text.data(), text.data() + text.size());
            return;
          }

          std::vector<unsigned char> mark(text.size(), 0);
          for (int p : positions) {
            if (p >= 0 && (std::size_t)p < text.size()) mark[(std::size_t)p] = 1;
          }

          ImDrawList* dl = ImGui::GetWindowDrawList();
          ImVec2 cur = ImGui::GetCursorScreenPos();
          const ImU32 colHi = ImGui::GetColorU32(ImGuiCol_PlotHistogram);
          const float lineH = ImGui::GetTextLineHeight();

          const char* base = text.data();
          std::size_t i = 0;
          while (i < text.size()) {
            const bool hi = mark[i] != 0;
            std::size_t j = i + 1;
            while (j < text.size() && (mark[j] != 0) == hi) ++j;
            const char* segB = base + i;
            const char* segE = base + j;
            dl->AddText(cur, hi ? colHi : colText, segB, segE);
            const ImVec2 sz = ImGui::CalcTextSize(segB, segE);
            cur.x += sz.x;
            i = j;
          }

          ImGui::Dummy(ImVec2(cur.x - ImGui::GetCursorScreenPos().x, lineH));
        };

        // Build a filtered+sorted list with fuzzy scores + highlight positions.
        struct BmRow {
          std::size_t idx;
          int score;
          std::vector<int> labelPos;
          std::vector<int> sysPos;
        };
        std::vector<BmRow> rows;
        rows.reserve(bookmarks.items.size());

        const bool hasFilter = anyNonSpaceLocal(bmFilter);
        for (std::size_t i = 0; i < bookmarks.items.size(); ++i) {
          const auto& bm = bookmarks.items[i];
          const std::string sysName = uiSystemNameById(bm.systemId);

          if (!hasFilter) {
            rows.push_back({i, 0, {}, {}});
            continue;
          }

          std::string hay;
          hay.reserve(bm.label.size() + sysName.size() + 2);
          hay += bm.label;
          hay.push_back(' ');
          hay += sysName;

          const auto r = ui::fuzzyMatch(bmFilter, hay);
          if (r.score < 0) continue;

          BmRow row{ i, r.score, {}, {} };
          const int split = (int)bm.label.size();
          for (int p : r.positions) {
            if (p < 0) continue;
            if (p < split) {
              row.labelPos.push_back(p);
            } else if (p > split) {
              row.sysPos.push_back(p - split - 1);
            }
          }
          rows.push_back(std::move(row));
        }

        std::sort(rows.begin(), rows.end(), [&](const BmRow& a, const BmRow& b) {
          if (hasFilter && a.score != b.score) return a.score > b.score;
          const auto& A = bookmarks.items[a.idx].label;
          const auto& B = bookmarks.items[b.idx].label;
          return A < B;
        });

        if (ImGui::BeginTable("##bm_table", 4,
                              ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY,
                              ImVec2(0.0f, 320.0f))) {
          ImGui::TableSetupColumn("Type", ImGuiTableColumnFlags_WidthFixed, 58.0f);
          ImGui::TableSetupColumn("Label");
          ImGui::TableSetupColumn("System", ImGuiTableColumnFlags_WidthFixed, 200.0f);
          ImGui::TableSetupColumn("Actions", ImGuiTableColumnFlags_WidthFixed, 190.0f);
          ImGui::TableHeadersRow();

          for (std::size_t n = 0; n < rows.size(); ++n) {
            const BmRow& row = rows[n];
            const std::size_t i = row.idx;
            const auto& bm = bookmarks.items[i];
            ImGui::TableNextRow();

            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(bm.kind == ui::BookmarkKind::Station ? "Station" : "System");

            ImGui::TableSetColumnIndex(1);
            if (hasFilter) {
              drawHighlighted(bm.label, row.labelPos, ImGui::GetColorU32(ImGuiCol_Text));
            } else {
              ImGui::TextUnformatted(bm.label.c_str());
            }

            ImGui::TableSetColumnIndex(2);
            {
              const std::string sysName = uiSystemNameById(bm.systemId);
              if (hasFilter) {
                drawHighlighted(sysName, row.sysPos, ImGui::GetColorU32(ImGuiCol_Text));
              } else {
                ImGui::TextUnformatted(sysName.c_str());
              }
            }

            ImGui::TableSetColumnIndex(3);
            const std::string bGo = "Go##" + std::to_string((unsigned long long)i);
            const std::string bDel = "Remove##" + std::to_string((unsigned long long)i);
            if (ImGui::SmallButton(bGo.c_str())) {
              galaxySelectedSystemId = bm.systemId;
              showGalaxy = true;
              if (currentSystem) {
                plotRouteToSystem(bm.systemId);
                if (bm.kind == ui::BookmarkKind::Station) {
                  pendingArrivalTargetStationId = bm.stationId;
                }
              }
            }
            ImGui::SameLine();
            if (ImGui::SmallButton(bDel.c_str())) {
              bool removed = false;
              if (bm.kind == ui::BookmarkKind::System) {
                removed = ui::removeSystemBookmark(bookmarks, bm.systemId);
              } else {
                removed = ui::removeStationBookmark(bookmarks, bm.systemId, bm.stationId);
              }
              if (removed) {
                bookmarksDirty = true;
                toast(toasts, "Removed bookmark.", 1.2);
              }
            }
          }

          ImGui::EndTable();
        }

        ImGui::Separator();

        // Quick add helpers.
        if (currentSystem) {
          const bool hasCur = ui::findSystemBookmarkIndex(bookmarks, currentSystem->stub.id).has_value();
          if (ImGui::Button(hasCur ? "Unbookmark current system" : "Bookmark current system")) {
            if (hasCur) {
              ui::removeSystemBookmark(bookmarks, currentSystem->stub.id);
              toast(toasts, "Removed bookmark.", 1.2);
            } else {
              ui::addOrUpdateSystemBookmark(bookmarks, currentSystem->stub.id, currentSystem->stub.name);
              toast(toasts, "Bookmarked current system.", 1.2);
            }
            bookmarksDirty = true;
          }
        }
        ImGui::SameLine();
        if (galaxySelectedSystemId != 0) {
          const bool hasSel = ui::findSystemBookmarkIndex(bookmarks, galaxySelectedSystemId).has_value();
          const auto& sys = universe.getSystem(galaxySelectedSystemId);
          if (ImGui::Button(hasSel ? "Unbookmark selected system" : "Bookmark selected system")) {
            if (hasSel) {
              ui::removeSystemBookmark(bookmarks, galaxySelectedSystemId);
              toast(toasts, "Removed bookmark.", 1.2);
            } else {
              ui::addOrUpdateSystemBookmark(bookmarks, galaxySelectedSystemId, sys.stub.name);
              toast(toasts, "Bookmarked selected system.", 1.2);
            }
            bookmarksDirty = true;
          }
        } else {
          ImGui::BeginDisabled();
          ImGui::Button("Bookmark selected system");
          ImGui::EndDisabled();
        }
      }
      ImGui::End();
    }

    // ---- Notifications (toast history) ----
    if (showNotifications) {
      ImGui::SetNextWindowSize(ImVec2(820, 420), ImGuiCond_FirstUseEver);
      if (ImGui::Begin("Notifications", &showNotifications)) {
        ImGui::TextUnformatted("On-screen toasts are automatically recorded here.");
        ImGui::SameLine();
        ImGui::TextDisabled("(%zu entries)", toastHistory.size());

        ImGui::InputTextWithHint("##notif_filter", "Filter...", notificationsFilter, sizeof(notificationsFilter));

        if (ImGui::Button("Clear")) {
          toastHistory.clear();
          notificationsSelected = -1;
          toast(toasts, "Notifications cleared.", 1.4);
        }
        ImGui::SameLine();
        const bool hasSelection = (notificationsSelected >= 0 && (std::size_t)notificationsSelected < toastHistory.size());
        if (!hasSelection) {
          ImGui::BeginDisabled();
        }
        if (ImGui::Button("Copy selected") && hasSelection) {
          ImGui::SetClipboardText(toastHistory[(std::size_t)notificationsSelected].text.c_str());
          toast(toasts, "Copied notification to clipboard.", 1.2);
        }
        if (!hasSelection) {
          ImGui::EndDisabled();
        }

        ImGui::Separator();

        auto anyNonSpaceLocal = [&](const char* s) {
          if (!s) return false;
          for (const unsigned char* p = (const unsigned char*)s; *p; ++p) {
            if (std::isspace(*p) == 0) return true;
          }
          return false;
        };

        auto drawHighlighted = [&](std::string_view text, const std::vector<int>& positions, ImU32 colText) {
          if (positions.empty()) {
            ImGui::TextUnformatted(text.data(), text.data() + text.size());
            return;
          }

          std::vector<unsigned char> mark(text.size(), 0);
          for (int p : positions) {
            if (p >= 0 && (std::size_t)p < text.size()) mark[(std::size_t)p] = 1;
          }

          ImDrawList* dl = ImGui::GetWindowDrawList();
          ImVec2 cur = ImGui::GetCursorScreenPos();
          const ImU32 colHi = ImGui::GetColorU32(ImGuiCol_PlotHistogram);
          const float lineH = ImGui::GetTextLineHeight();

          const char* base = text.data();
          std::size_t i = 0;
          while (i < text.size()) {
            const bool hi = mark[i] != 0;
            std::size_t j = i + 1;
            while (j < text.size() && (mark[j] != 0) == hi) ++j;
            const char* segB = base + i;
            const char* segE = base + j;
            dl->AddText(cur, hi ? colHi : colText, segB, segE);
            const ImVec2 sz = ImGui::CalcTextSize(segB, segE);
            cur.x += sz.x;
            i = j;
          }

          ImGui::Dummy(ImVec2(cur.x - ImGui::GetCursorScreenPos().x, lineH));
        };

        const bool hasFilter = anyNonSpaceLocal(notificationsFilter);
        static bool notifSortByRelevance = false;
        if (hasFilter) {
          ImGui::SameLine();
          ImGui::Checkbox("Sort by relevance", &notifSortByRelevance);
        }

        struct NotifRow {
          int idx;
          int score;
          std::vector<int> pos;
        };
        std::vector<NotifRow> notifRows;
        notifRows.reserve(toastHistory.size());

        for (int i = 0; i < (int)toastHistory.size(); ++i) {
          const auto& e = toastHistory[(std::size_t)i];
          if (!hasFilter) {
            notifRows.push_back({i, 0, {}});
            continue;
          }
          const auto r = ui::fuzzyMatch(notificationsFilter, e.text);
          if (r.score < 0) continue;
          notifRows.push_back({i, r.score, r.positions});
        }

        if (hasFilter && notifSortByRelevance) {
          std::stable_sort(notifRows.begin(), notifRows.end(), [&](const NotifRow& a, const NotifRow& b) {
            if (a.score != b.score) return a.score > b.score;
            return a.idx > b.idx; // prefer newer when scores tie
          });
        }

        auto fmtSimTime = [&](double tDays) -> std::string {
          const long long totalSec = (long long)std::llround(tDays * 86400.0);
          const long long day = totalSec / 86400;
          const int secInDay = (int)(totalSec - day * 86400);
          const int hh = secInDay / 3600;
          const int mm = (secInDay / 60) % 60;
          const int ss = secInDay % 60;
          char buf[64];
          if (day > 0) {
            std::snprintf(buf, sizeof(buf), "D%lld %02d:%02d:%02d", day, hh, mm, ss);
          } else {
            std::snprintf(buf, sizeof(buf), "%02d:%02d:%02d", hh, mm, ss);
          }
          return std::string(buf);
        };

        if (ImGui::BeginTable("##notif_table", 2,
                              ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY
                              | ImGuiTableFlags_SizingStretchProp)) {
          ImGui::TableSetupColumn("Time", ImGuiTableColumnFlags_WidthFixed, 120.0f);
          ImGui::TableSetupColumn("Message", ImGuiTableColumnFlags_WidthStretch);
          ImGui::TableHeadersRow();

          for (int n = 0; n < (int)notifRows.size(); ++n) {
            const NotifRow& row = notifRows[(std::size_t)n];
            const int i = row.idx;
            const auto& e = toastHistory[(std::size_t)i];

            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextDisabled("%s", fmtSimTime(e.timeDays).c_str());
            ImGui::TableSetColumnIndex(1);

            const bool sel = (notificationsSelected == i);
            const std::string selId = "##notif_sel_" + std::to_string(i);
            const ImVec2 startPos = ImGui::GetCursorScreenPos();
            const bool clicked = ImGui::Selectable(selId.c_str(), sel,
                                                  ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowItemOverlap);
            const ImVec2 afterPos = ImGui::GetCursorScreenPos();

            ImGui::SetCursorScreenPos(startPos);
            if (hasFilter) {
              drawHighlighted(e.text, row.pos, ImGui::GetColorU32(ImGuiCol_Text));
            } else {
              ImGui::TextUnformatted(e.text.c_str());
            }
            ImGui::SetCursorScreenPos(afterPos);

            if (clicked) notificationsSelected = i;
          }

          ImGui::EndTable();
        }
      }
      ImGui::End();
    }

    // ---- UI Settings ----
    if (showUiSettingsWindow) {
      ImGui::SetNextWindowSize(ImVec2(760, 560), ImGuiCond_FirstUseEver);
      if (ImGui::Begin("UI Settings", &showUiSettingsWindow)) {
        ImGui::TextDisabled("Saved to %s", uiSettingsPath.c_str());
        ImGui::SameLine();
        ImGui::TextDisabled("| Layout: %s", (io.IniFilename ? io.IniFilename : "(none)"));

        if (uiSettingsDirty) {
          ImGui::SameLine();
          ImGui::TextColored(ImVec4(1.0f, 0.85f, 0.25f, 1.0f), "*unsaved*");
        }

        ImGui::Separator();

        // Top actions.
        if (ImGui::Button("Save")) {
          syncUiSettingsFromRuntime();
          const bool ok = ui::saveUiSettingsToFile(uiSettings, uiSettingsPath);
          if (ok) {
            uiSettingsDirty = false;
          }
          toast(toasts, ok ? "Saved UI settings." : "Failed to save UI settings.", 1.8);
        }
        ImGui::SameLine();
        if (ImGui::Button("Load")) {
          ui::UiSettings loaded = ui::makeDefaultUiSettings();
          if (ui::loadUiSettingsFromFile(uiSettingsPath, loaded)) {
            uiSettings = loaded;
            applyUiSettingsToRuntime(uiSettings, /*loadImGuiIni=*/true);
            uiSettingsDirty = false;
            toast(toasts, "Loaded UI settings.", 1.8);
          } else {
            toast(toasts, "No UI settings file found.", 1.8);
          }
        }
        ImGui::SameLine();
        if (ImGui::Button("Restore defaults")) {
          uiSettings = ui::makeDefaultUiSettings();
          applyUiSettingsToRuntime(uiSettings, /*loadImGuiIni=*/false);
          uiSettingsDirty = true;
          toast(toasts, "Restored default UI settings (not saved).", 2.2);
        }
        ImGui::SameLine();
        if (ImGui::Checkbox("Auto-save on exit", &uiSettingsAutoSaveOnExit)) {
          uiSettingsDirty = true;
        }

        ImGui::Separator();

        if (ImGui::CollapsingHeader("Appearance", ImGuiTreeNodeFlags_DefaultOpen)) {
          const char* themeNames[] = {"Dark", "Light", "Classic", "High Contrast"};
          int themeIdx = (int)uiTheme;
          if (ImGui::Combo("Theme", &themeIdx, themeNames, IM_ARRAYSIZE(themeNames))) {
            uiTheme = (ui::UiTheme)themeIdx;
            rebuildUiStyle(uiTheme, uiScale);
            io.FontGlobalScale = uiScale;
            uiScaleApplied = uiScale;
            uiSettingsDirty = true;
          }

          if (uiTheme == ui::UiTheme::HighContrast) {
            ImGui::TextDisabled("High contrast mode targets readability on dark backgrounds.");
          }
        }

        if (ImGui::CollapsingHeader("Scaling", ImGuiTreeNodeFlags_DefaultOpen)) {
          if (ImGui::Checkbox("Auto scale from DPI", &uiAutoScaleFromDpi)) {
            recomputeUiDpiScale();
            applyUiScaleNow();
            uiSettingsDirty = true;
          }
          ImGui::SameLine();
          if (ImGui::Button("Recompute DPI")) {
            recomputeUiDpiScale();
            applyUiScaleNow();
            uiSettingsDirty = true;
          }

          ImGui::SetNextItemWidth(260.0f);
          if (ImGui::SliderFloat("User scale", &uiScaleUser, 0.50f, 2.50f, "%.2fx")) {
            applyUiScaleNow();
            uiSettingsDirty = true;
          }
          ImGui::TextDisabled("Effective: %.2fx (DPI %.2fx × User %.2fx)", uiScale, uiScaleDpi, uiScaleUser);
        }

        if (ImGui::CollapsingHeader("Layout profiles", ImGuiTreeNodeFlags_DefaultOpen)) {
          ImGui::TextUnformatted("Dear ImGui stores window positions/docking in an ini file.");
          ImGui::TextUnformatted("Profiles let you keep multiple layouts (e.g. combat vs. trading)."
                                " Use Save/Load to persist high-level UI preferences.");

          namespace fs = std::filesystem;
          const fs::path dir("ui_layouts");

          struct LayoutProfile { std::string name; std::string path; bool isDefault; };
          std::vector<LayoutProfile> profiles;
          profiles.push_back(LayoutProfile{"Default", "imgui.ini", true});
          if (fs::exists(dir) && fs::is_directory(dir)) {
            for (const auto& ent : fs::directory_iterator(dir)) {
              if (!ent.is_regular_file()) continue;
              const fs::path p = ent.path();
              if (p.extension() != ".ini") continue;
              const std::string name = p.stem().string();
              profiles.push_back(LayoutProfile{name, p.string(), false});
            }
          }
          std::sort(profiles.begin() + 1, profiles.end(), [](const LayoutProfile& a, const LayoutProfile& b) {
            return a.name < b.name;
          });

          // Pick current index by matching the active ini filename.
          int currentProfile = 0;
          for (int i = 0; i < (int)profiles.size(); ++i) {
            if (profiles[(std::size_t)i].path == uiIniFilename) {
              currentProfile = i;
              break;
            }
          }

          static int selectedProfile = 0;
          if (selectedProfile < 0 || selectedProfile >= (int)profiles.size()) selectedProfile = currentProfile;

          if (ImGui::BeginListBox("Profiles", ImVec2(-FLT_MIN, 160.0f))) {
            for (int i = 0; i < (int)profiles.size(); ++i) {
              const bool sel = (selectedProfile == i);
              std::string label = profiles[(std::size_t)i].name;
              if (profiles[(std::size_t)i].path == uiIniFilename) label += "  (active)";
              if (ImGui::Selectable(label.c_str(), sel)) selectedProfile = i;
            }
            ImGui::EndListBox();
          }

          const bool canActivate = (selectedProfile >= 0 && selectedProfile < (int)profiles.size());
          if (!canActivate) ImGui::BeginDisabled();
          if (ImGui::Button("Activate selected") && canActivate) {
            uiIniFilename = profiles[(std::size_t)selectedProfile].path;
            io.IniFilename = uiIniFilename.c_str();
            if (io.IniFilename && io.IniFilename[0] != '\0') {
              ImGui::LoadIniSettingsFromDisk(io.IniFilename);
            }
            uiSettingsDirty = true;
            toast(toasts, "Activated layout profile.", 1.6);
          }
          ImGui::SameLine();
          if (ImGui::Button("Save layout now")) {
            if (io.IniFilename && io.IniFilename[0] != '\0') {
              ImGui::SaveIniSettingsToDisk(io.IniFilename);
              toast(toasts, "Saved layout ini.", 1.4);
            }
          }
          if (!canActivate) ImGui::EndDisabled();

          // Save-as workflow.
          static char newProfileName[64] = "combat";
          ImGui::Separator();
          ImGui::InputTextWithHint("##new_profile", "New profile name (letters/numbers/-/_)", newProfileName, sizeof(newProfileName));
          ImGui::SameLine();
          if (ImGui::Button("Save As")) {
            auto sanitize = [](const char* in) {
              std::string out;
              for (const unsigned char* p = (const unsigned char*)in; *p; ++p) {
                const char c = (char)*p;
                if (std::isalnum(*p) || c == '_' || c == '-') out.push_back(c);
              }
              if (out.empty()) out = "profile";
              return out;
            };

            const std::string safe = sanitize(newProfileName);
            fs::create_directories(dir);
            const fs::path p = dir / (safe + ".ini");
            ImGui::SaveIniSettingsToDisk(p.string().c_str());
            uiIniFilename = p.string();
            io.IniFilename = uiIniFilename.c_str();
            uiSettingsDirty = true;
            toast(toasts, "Saved new layout profile.", 1.8);
          }

          // Deleting is optional but handy.
          const bool selectedIsDefault = canActivate ? profiles[(std::size_t)selectedProfile].isDefault : true;
          if (selectedIsDefault) ImGui::BeginDisabled();
          if (ImGui::Button("Delete selected") && canActivate && !selectedIsDefault) {
            const std::string path = profiles[(std::size_t)selectedProfile].path;
            // Safety: only delete files inside ui_layouts.
            if (path.rfind("ui_layouts", 0) == 0) {
              std::error_code ec;
              fs::remove(fs::path(path), ec);
              if (!ec) {
                toast(toasts, "Deleted layout profile.", 1.6);
                if (uiIniFilename == path) {
                  uiIniFilename = "imgui.ini";
                  io.IniFilename = uiIniFilename.c_str();
                  uiSettingsDirty = true;
                }
              } else {
                toast(toasts, "Failed to delete profile.", 1.6);
              }
            } else {
              toast(toasts, "Refusing to delete outside ui_layouts/", 2.0);
            }
          }
          if (selectedIsDefault) ImGui::EndDisabled();
        }

#ifdef IMGUI_HAS_DOCK
        if (ImGui::CollapsingHeader("Docking", ImGuiTreeNodeFlags_DefaultOpen)) {
          if (ImGui::Checkbox("Enable DockSpace", &uiDockingEnabled)) uiSettingsDirty = true;
          if (ImGui::Checkbox("Passthrough central view", &uiDockPassthruCentral)) { uiDockResetLayout = true; uiSettingsDirty = true; }
          if (ImGui::Checkbox("Lock central view", &uiDockLockCentralView)) { uiDockResetLayout = true; uiSettingsDirty = true; }

          ImGui::SetNextItemWidth(260.0f);
          if (ImGui::SliderFloat("Left width", &uiDockLeftRatio, 0.10f, 0.45f, "%.2f")) { uiDockResetLayout = true; uiSettingsDirty = true; }
          ImGui::SetNextItemWidth(260.0f);
          if (ImGui::SliderFloat("Right width", &uiDockRightRatio, 0.10f, 0.45f, "%.2f")) { uiDockResetLayout = true; uiSettingsDirty = true; }
          ImGui::SetNextItemWidth(260.0f);
          if (ImGui::SliderFloat("Bottom height", &uiDockBottomRatio, 0.10f, 0.45f, "%.2f")) { uiDockResetLayout = true; uiSettingsDirty = true; }

          if (ImGui::Button("Reset dock layout")) {
            uiDockResetLayout = true;
            toast(toasts, "Rebuilding default dock layout...", 1.4);
          }
          ImGui::TextDisabled("Tip: ImGui also autosaves window positions/layout to the active .ini.");
        }
#endif
      }
      ImGui::End();
    }

    // ---- HUD Settings ----
    if (showHudSettingsWindow) {
      ImGui::SetNextWindowSize(ImVec2(760, 560), ImGuiCond_FirstUseEver);
      if (ImGui::Begin("HUD Settings", &showHudSettingsWindow)) {
        ImGui::TextDisabled("Saved to %s", hudSettingsPath.c_str());

        if (hudSettingsDirty) {
          ImGui::SameLine();
          ImGui::TextColored(ImVec4(1.0f, 0.85f, 0.25f, 1.0f), "*unsaved*");
        }

        ImGui::Separator();

        // Top actions.
        if (ImGui::Button("Save")) {
          syncHudSettingsFromRuntime();
          const bool ok = ui::saveHudSettingsToFile(hudSettings, hudSettingsPath);
          if (ok) {
            hudSettingsSaved = hudSettings;
            hudSettingsDirty = false;
          }
          toast(toasts, ok ? "Saved HUD settings." : "Failed to save HUD settings.", 1.8);
        }
        ImGui::SameLine();
        if (ImGui::Button("Load")) {
          ui::HudSettings loaded = ui::makeDefaultHudSettings();
          if (ui::loadHudSettingsFromFile(hudSettingsPath, loaded)) {
            hudSettings = loaded;
            applyHudSettingsToRuntime(hudSettings);
            hudSettingsSaved = hudSettings;
            hudSettingsDirty = false;
            toast(toasts, "Loaded HUD settings.", 1.8);
          } else {
            toast(toasts, "No HUD settings file found.", 1.8);
          }
        }
        ImGui::SameLine();
        if (ImGui::Button("Restore defaults")) {
          hudSettings = ui::makeDefaultHudSettings();
          applyHudSettingsToRuntime(hudSettings);
          // Keep 'saved' snapshot so this shows as unsaved until the user saves.
          toast(toasts, "Restored default HUD settings (not saved).", 2.2);
        }
        ImGui::SameLine();
        ImGui::Checkbox("Auto-save on exit", &hudSettingsAutoSaveOnExit);

        ImGui::Separator();

        // Presets
        ImGui::TextDisabled("Presets (applied immediately; click Save to persist):");
        if (ImGui::Button("Combat")) {
          showRadarHud = true;
          objectiveHudEnabled = true;
          hudThreatOverlayEnabled = true;
          hudJumpOverlay = true;
          showTacticalOverlay = true;
          tacticalShowLabels = true;
          hudOffscreenTargetIndicator = true;
          hudCombatHud = true;
          hudUseProceduralReticle = true;
          hudShowWeaponRings = true;
          hudShowHeatRing = true;
          hudShowDistributorRings = true;
          hudShowLeadIndicator = true;
          hudShowFlightPathMarker = true;
          toast(toasts, "Applied Combat HUD preset.", 1.6);
        }
        ImGui::SameLine();
        if (ImGui::Button("Exploration")) {
          showRadarHud = true;
          objectiveHudEnabled = true;
          hudThreatOverlayEnabled = true;
          hudJumpOverlay = true;
          showTacticalOverlay = true;
          tacticalShowLabels = true;
          hudOffscreenTargetIndicator = true;
          hudCombatHud = true;
          hudUseProceduralReticle = true;
          hudShowWeaponRings = false;
          hudShowHeatRing = false;
          hudShowDistributorRings = false;
          hudShowLeadIndicator = false;
          hudShowFlightPathMarker = true;
          toast(toasts, "Applied Exploration HUD preset.", 1.6);
        }
        ImGui::SameLine();
        if (ImGui::Button("Cinematic")) {
          showRadarHud = false;
          objectiveHudEnabled = false;
          hudThreatOverlayEnabled = false;
          hudJumpOverlay = false;
          showTacticalOverlay = false;
          hudOffscreenTargetIndicator = false;
          hudCombatHud = false;
          toast(toasts, "Applied Cinematic HUD preset.", 1.6);
        }

        ImGui::Separator();

        if (ImGui::CollapsingHeader("Overlays", ImGuiTreeNodeFlags_DefaultOpen)) {
          ImGui::Checkbox("Radar HUD", &showRadarHud);
          ImGui::SameLine();
          ImGui::Checkbox("Objective", &objectiveHudEnabled);
          ImGui::SameLine();
          ImGui::Checkbox("Threat", &hudThreatOverlayEnabled);
          ImGui::SameLine();
          ImGui::Checkbox("Jump", &hudJumpOverlay);

          ImGui::Checkbox("Tactical overlay", &showTacticalOverlay);
          ImGui::SameLine();
          ImGui::Checkbox("Show tactical labels", &tacticalShowLabels);

          ImGui::Separator();
          ImGui::Checkbox("Offscreen target indicator", &hudOffscreenTargetIndicator);

          if (ImGui::Button("Open HUD layout editor")) {
            showHudLayoutWindow = true;
            hudLayoutEditMode = true;
          }
          ImGui::SameLine();
          ImGui::TextDisabled("(drag overlays; right-click radar for quick tweaks)");
        }

        if (ImGui::CollapsingHeader("Radar", ImGuiTreeNodeFlags_DefaultOpen)) {
          const double minKm = 25000.0;
          const double maxKm = 1200000.0;
          ImGui::SliderScalar("Range (km)", ImGuiDataType_Double, &radarRangeKm, &minKm, &maxKm, "%.0f", ImGuiSliderFlags_Logarithmic);
          ImGui::SliderInt("Max blips", &radarMaxBlips, 16, 160);
        }

        if (ImGui::CollapsingHeader("Combat symbology", ImGuiTreeNodeFlags_DefaultOpen)) {
          ImGui::Checkbox("Enable combat HUD", &hudCombatHud);
          if (!hudCombatHud) ImGui::BeginDisabled();

          ImGui::Checkbox("Procedural reticle", &hudUseProceduralReticle);
          ImGui::Checkbox("Weapon rings", &hudShowWeaponRings);
          ImGui::Checkbox("Heat ring", &hudShowHeatRing);
          ImGui::Checkbox("Distributor rings", &hudShowDistributorRings);

          ImGui::SetNextItemWidth(260.0f);
          ImGui::SliderFloat("Reticle size (px)", &hudReticleSizePx, 10.0f, 120.0f, "%.0f");
          ImGui::SetNextItemWidth(260.0f);
          ImGui::SliderFloat("Reticle alpha", &hudReticleAlpha, 0.05f, 1.0f, "%.2f");

          ImGui::Separator();
          ImGui::Checkbox("Lead indicator", &hudShowLeadIndicator);
          if (hudShowLeadIndicator) {
            ImGui::Indent();
            ImGui::Checkbox("Use last fired weapon", &hudLeadUseLastFiredWeapon);
            ImGui::SetNextItemWidth(260.0f);
            ImGui::SliderFloat("Lead size (px)", &hudLeadSizePx, 8.0f, 80.0f, "%.0f");
            ImGui::SetNextItemWidth(260.0f);
            const double minSec = 1.0;
            const double maxSec = 120.0;
            ImGui::SliderScalar("Max lead time (sec)", ImGuiDataType_Double, &hudLeadMaxTimeSec, &minSec, &maxSec, "%.1f");
            ImGui::Unindent();
          }

          ImGui::Separator();
          ImGui::Checkbox("Flight path marker", &hudShowFlightPathMarker);
          if (hudShowFlightPathMarker) {
            ImGui::Indent();
            ImGui::Checkbox("Local frame", &hudFlightMarkerUseLocalFrame);
            ImGui::Checkbox("Clamp to screen", &hudFlightMarkerClampToEdge);
            ImGui::SetNextItemWidth(260.0f);
            ImGui::SliderFloat("Marker size (px)", &hudFlightMarkerSizePx, 8.0f, 80.0f, "%.0f");
            ImGui::Unindent();
          }

          if (!hudCombatHud) ImGui::EndDisabled();
        }

        if (ImGui::CollapsingHeader("Tactical overlay", ImGuiTreeNodeFlags_DefaultOpen)) {
          const double minKm = 20000.0;
          const double maxKm = 2000000.0;
          ImGui::SliderScalar("Range (km)", ImGuiDataType_Double, &tacticalRangeKm, &minKm, &maxKm, "%.0f", ImGuiSliderFlags_Logarithmic);
          ImGui::SliderInt("Max markers", &tacticalMaxMarkers, 16, 256);

          ImGui::Separator();
          ImGui::TextDisabled("Filters");
          ImGui::Checkbox("Stations", &tacticalShowStations);
          ImGui::SameLine();
          ImGui::Checkbox("Planets", &tacticalShowPlanets);
          ImGui::SameLine();
          ImGui::Checkbox("Contacts", &tacticalShowContacts);

          ImGui::Checkbox("Cargo", &tacticalShowCargo);
          ImGui::SameLine();
          ImGui::Checkbox("Asteroids", &tacticalShowAsteroids);
          ImGui::SameLine();
          ImGui::Checkbox("Signals", &tacticalShowSignals);
        }
      }
      ImGui::End();
    }

    // ---- Command Palette ----
    if (commandPalette.open) {
      // Rebuild the item list only while the palette is open.
      commandPaletteItems.clear();
      commandPaletteItems.reserve(512);

      auto onOff = [](bool v) { return v ? "ON" : "OFF"; };
      auto addItem = [&](std::string label, std::string detail, std::string shortcut, int priority,
                         std::function<void()> fn) {
        commandPaletteItems.push_back(game::PaletteItem{std::move(label), std::move(detail), std::move(shortcut), priority, std::move(fn)});
      };

      // UI / windows
      addItem("Toggle Galaxy window", std::string("Window • ") + onOff(showGalaxy), game::chordLabel(controls.actions.toggleGalaxy), 100, [&]() { showGalaxy = !showGalaxy; });
      addItem("Toggle Ship/Status window", std::string("Window • ") + onOff(showShip), game::chordLabel(controls.actions.toggleShip), 100, [&]() { showShip = !showShip; });
      addItem("Toggle Market window", std::string("Window • ") + onOff(showEconomy), game::chordLabel(controls.actions.toggleMarket), 100, [&]() { showEconomy = !showEconomy; });
      addItem("Toggle Contacts window", std::string("Window • ") + onOff(showContacts), game::chordLabel(controls.actions.toggleContacts), 100, [&]() { showContacts = !showContacts; });
      addItem("Toggle Missions window", std::string("Window • ") + onOff(showMissions), game::chordLabel(controls.actions.toggleMissions), 100, [&]() { showMissions = !showMissions; });
      addItem("Toggle Scanner window", std::string("Window • ") + onOff(showScanner), game::chordLabel(controls.actions.toggleScanner), 100, [&]() { showScanner = !showScanner; });
      addItem("Toggle Trade Finder", std::string("Window • ") + onOff(showTrade), game::chordLabel(controls.actions.toggleTrade), 100, [&]() { showTrade = !showTrade; });
      addItem("Toggle Pilot Guide", std::string("Window • ") + onOff(showGuide), game::chordLabel(controls.actions.toggleGuide), 90, [&]() { showGuide = !showGuide; });
      addItem("Toggle Hangar", std::string("Window • ") + onOff(showHangar), game::chordLabel(controls.actions.toggleHangar), 90, [&]() { showHangar = !showHangar; });
      addItem("Toggle World Visuals", std::string("Window • ") + onOff(showWorldVisuals), game::chordLabel(controls.actions.toggleWorldVisuals), 80, [&]() { showWorldVisuals = !showWorldVisuals; });
      addItem("Toggle Controls window", std::string("Window • ") + onOff(controlsWindow.open), game::chordLabel(controls.actions.toggleControlsWindow), 80,
              [&]() {
                controlsWindow.open = !controlsWindow.open;
                if (controlsWindow.open) controlsWindow.focusFilter = true;
              });
      addItem("Toggle Notifications window", std::string("Window • ") + onOff(showNotifications), std::string(), 80, [&]() { showNotifications = !showNotifications; });
      addItem("Toggle Console window", std::string("Window • ") + onOff(consoleWindow.open), std::string(), 80, [&]() { consoleWindow.open = !consoleWindow.open; if (consoleWindow.open) consoleWindow.focusInput = true; });
      addItem("Toggle Bookmarks window", std::string("Window • ") + onOff(showBookmarksWindow), std::string(), 80, [&]() { showBookmarksWindow = !showBookmarksWindow; });
      addItem("Toggle UI Settings window", std::string("Window • ") + onOff(showUiSettingsWindow), std::string(), 80, [&]() { showUiSettingsWindow = !showUiSettingsWindow; });
      addItem("Toggle HUD Settings window", std::string("Window • ") + onOff(showHudSettingsWindow), std::string(), 80, [&]() { showHudSettingsWindow = !showHudSettingsWindow; });

      // HUD
      addItem("Toggle Radar HUD", std::string("HUD • ") + onOff(showRadarHud), game::chordLabel(controls.actions.toggleRadarHud), 70, [&]() { showRadarHud = !showRadarHud; });
      addItem("Toggle Tactical overlay", std::string("HUD • ") + onOff(showTacticalOverlay), game::chordLabel(controls.actions.toggleTacticalOverlay), 70, [&]() { showTacticalOverlay = !showTacticalOverlay; });
      addItem("Toggle HUD edit mode", std::string("HUD • ") + onOff(hudLayoutEditMode), game::chordLabel(controls.actions.hudLayoutToggleEdit), 60, [&]() { hudLayoutEditMode = !hudLayoutEditMode; });

      // Gameplay
      addItem(paused ? "Unpause" : "Pause", "Gameplay", game::chordLabel(controls.actions.pause), 60, [&]() { paused = !paused; });
      addItem("Set time scale: 1x", "Time", std::string(), 40, [&]() { timeScale = 1.0; paused = false; });
      addItem("Set time scale: 10x", "Time", std::string(), 40, [&]() { timeScale = 10.0; paused = false; });
      addItem("Set time scale: 60x", "Time", std::string(), 40, [&]() { timeScale = 60.0; paused = false; });
      addItem("Set time scale: 600x", "Time", std::string(), 40, [&]() { timeScale = 600.0; paused = false; });

      // Missions (quick tracking)
      if (trackedMissionId != 0) {
        addItem("Untrack current mission", "Missions", std::string(), 55, [&]() {
          trackedMissionId = 0;
          toast(toasts, "Mission tracking cleared.", 1.6);
        });
      }
      for (const auto& m : missions) {
        if (m.failed || m.completed) continue;
        const core::u64 id = m.id;
        const auto [nextSys, nextSt] = uiMissionNextStop(m);
        const std::string label = uiDescribeMission(m);

        addItem(label, "Missions • Track", std::string(), 20, [&, id, nextSys, nextSt]() {
          trackedMissionId = id;
          objectiveHudEnabled = true;
          showMissions = true;
          if (nextSys != 0) {
            plotRouteToSystem(nextSys);
            pendingArrivalTargetStationId = nextSt;
          }
          toast(toasts, "Tracking mission.", 1.6);
        });
      }

      // Bookmarks
      if (currentSystem) {
        const bool hasCur = ui::findSystemBookmarkIndex(bookmarks, currentSystem->stub.id).has_value();
        addItem(hasCur ? "Remove bookmark: current system" : "Bookmark current system",
                "Bookmarks", std::string(), 58, [&, hasCur]() {
                  if (!currentSystem) return;
                  if (hasCur) {
                    if (ui::removeSystemBookmark(bookmarks, currentSystem->stub.id)) {
                      bookmarksDirty = true;
                      toast(toasts, "Removed bookmark.", 1.4);
                    }
                  } else {
                    ui::addOrUpdateSystemBookmark(bookmarks, currentSystem->stub.id, currentSystem->stub.name);
                    bookmarksDirty = true;
                    toast(toasts, "Bookmarked current system.", 1.4);
                  }
                });
      }
      if (galaxySelectedSystemId != 0) {
        const bool hasSel = ui::findSystemBookmarkIndex(bookmarks, galaxySelectedSystemId).has_value();
        addItem(hasSel ? "Remove bookmark: selected system" : "Bookmark selected system",
                "Bookmarks", std::string(), 48, [&, hasSel]() {
                  if (galaxySelectedSystemId == 0) return;
                  const auto& sys = universe.getSystem(galaxySelectedSystemId);
                  if (hasSel) {
                    if (ui::removeSystemBookmark(bookmarks, galaxySelectedSystemId)) {
                      bookmarksDirty = true;
                      toast(toasts, "Removed bookmark.", 1.4);
                    }
                  } else {
                    ui::addOrUpdateSystemBookmark(bookmarks, galaxySelectedSystemId, sys.stub.name);
                    bookmarksDirty = true;
                    toast(toasts, "Bookmarked selected system.", 1.4);
                  }
                });
      }

      for (const auto& bm : bookmarks.items) {
        if (bm.systemId == 0) continue;
        const sim::SystemId sysId = bm.systemId;
        const sim::StationId stId = bm.stationId;

        std::string label = bm.label;
        if (label.empty()) label = (bm.kind == ui::BookmarkKind::Station) ? "Station" : "System";

        if (bm.kind == ui::BookmarkKind::System) {
          addItem(label, "Bookmarks • Route", std::string(), 30, [&, sysId]() {
            galaxySelectedSystemId = sysId;
            showGalaxy = true;
            plotRouteToSystem(sysId);
          });
        } else {
          addItem(label, "Bookmarks • Route (station)", std::string(), 30, [&, sysId, stId]() {
            galaxySelectedSystemId = sysId;
            showGalaxy = true;
            plotRouteToSystem(sysId);
            pendingArrivalTargetStationId = stId;
          });
        }
      }

      // Systems (nearby search radius scales with jump range)
      if (currentSystem) {
        const double radiusLy = std::clamp(fsdRangeLy * 90.0, 600.0, 3500.0);
        const auto stubs = universe.queryNearby(currentSystem->stub.posLy, radiusLy, 512);
        for (const auto& s : stubs) {
          if (s.id == currentSystem->stub.id) continue;
          const double d = (s.posLy - currentSystem->stub.posLy).length();
          const core::u64 sysId = s.id;
          std::string detail = std::string("System • ") + starClassName(s.primaryClass) + " • " + std::to_string((int)std::llround(d)) + " ly";
          if (s.stationCount > 0) {
            detail += " • ";
            detail += std::to_string(s.stationCount);
            detail += " st";
          }

          addItem(s.name, std::move(detail), std::string(), 10, [&, sysId]() {
            galaxySelectedSystemId = (sim::SystemId)sysId;
            showGalaxy = true;
            plotRouteToSystem((sim::SystemId)sysId);
          });
        }
      }

      // Stations (in current system)
      if (currentSystem) {
        for (std::size_t i = 0; i < currentSystem->stations.size(); ++i) {
          const auto& st = currentSystem->stations[i];
          const std::size_t idx = i;
          addItem(st.name, "Station • Target", std::string(), 15, [&, idx]() {
            target.kind = TargetKind::Station;
            target.index = idx;
            selectedStationIndex = (int)idx;
            showShip = true;
          });
        }
      }

      // Contacts (local space)
      for (std::size_t i = 0; i < contacts.size(); ++i) {
        const std::size_t idx = i;
        const auto& c = contacts[i];
        if (!c.alive) continue;
        const double dKm = (c.ship.positionKm() - ship.positionKm()).length();
        std::string detail = std::string("Contact • ") + contactRoleName(c.role) + " • " + std::to_string((int)std::llround(dKm)) + " km";
        addItem(c.name, std::move(detail), std::string(), 5, [&, idx]() {
          target.kind = TargetKind::Contact;
          target.index = idx;
          toast(toasts, "Target set.", 1.0);
        });
      }
    }

    // Draw palette last so it appears above other windows.
    (void)game::drawCommandPalette(commandPalette, commandPaletteItems);

    // ---- Offscreen hangar preview render (render-to-texture -> ImGui) ----
    // Render AFTER building UI so the preview uses the latest yaw/pitch/zoom (drag UI)
    // but BEFORE ImGui::Render, so ImGui draws the up-to-date texture.
    if (showHangar && hangarTarget.color().handle() != 0) {
      if (hangarAnimate) {
        hangarYawDeg += (float)(dtReal * (double)hangarSpinRadPerSec * (180.0 / math::kPi));
      }

      // Render into the hangar target.
      hangarTarget.begin();
      glViewport(0, 0, hangarTarget.width(), hangarTarget.height());
      glEnable(GL_DEPTH_TEST);
      glDepthMask(GL_TRUE);
      glDisable(GL_BLEND);
      glClearColor(0.035f, 0.040f, 0.050f, 1.0f);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      // Simple camera looking down -Z with the ship at the origin.
      render::Camera hangarCam;
      const double aspect = (double)hangarTarget.width() / (double)std::max(1, hangarTarget.height());
      hangarCam.setPerspective(math::degToRad(45.0), aspect, 0.01, 50.0);
      hangarCam.setPosition({0.0, 0.0, (double)hangarZoom});
      hangarCam.setOrientation(math::Quatd::identity());

      float viewF[16], projF[16];
      matToFloat(hangarCam.viewMatrix(), viewF);
      matToFloat(hangarCam.projectionMatrix(), projF);
      meshRenderer.setViewProj(viewF, projF);

      // Preview light: offset point-light so the model reads nicely.
      meshRenderer.setLightPos(2.5f, 1.8f, 3.5f);
      meshRenderer.setMesh(&cube);
      meshRenderer.setUnlit(false);
      meshRenderer.setAlphaFromTexture(false);
      meshRenderer.setTexture((shipLiveryTex.handle() != 0) ? &shipLiveryTex : &checker);

      const math::Quatd q = math::Quatd::fromAxisAngle({0, 1, 0}, math::degToRad((double)hangarYawDeg))
                          * math::Quatd::fromAxisAngle({1, 0, 0}, math::degToRad((double)hangarPitchDeg));

      std::vector<render::InstanceData> inst;
      inst.reserve(1);
      inst.push_back(makeInst({0.0, 0.0, 0.0},
                              {0.35, 0.20, 0.60},
                              q,
                              1.0f, 1.0f, 1.0f));
      meshRenderer.drawInstances(inst);

      // Restore backbuffer for UI.
      hangarTarget.end();
      render::gl::BindFramebuffer(GL_FRAMEBUFFER, 0);
      glViewport(0, 0, w, h);

      // Restore default light so next frame's world render is correct even if
      // something draws before the per-frame setLightPos call.
      meshRenderer.setLightPos(0.0f, 0.0f, 0.0f);
    }

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    SDL_GL_SwapWindow(window);
  }

  // Persist HUD layout on exit (optional).
  if (hudLayoutAutoSaveOnExit) {
    hudLayout.widget(ui::HudWidgetId::Radar).enabled = showRadarHud;
    hudLayout.widget(ui::HudWidgetId::Objective).enabled = objectiveHudEnabled;
    hudLayout.widget(ui::HudWidgetId::Threat).enabled = hudThreatOverlayEnabled;
    hudLayout.widget(ui::HudWidgetId::Jump).enabled = hudJumpOverlay;
    ui::saveToFile(hudLayout, hudLayoutPath);
  }

  // Persist HUD settings on exit.
  {
    syncHudSettingsFromRuntime();
    const bool autoSavePrefChanged = (hudSettingsSaved.autoSaveOnExit != hudSettingsAutoSaveOnExit);
    if ((hudSettingsAutoSaveOnExit || autoSavePrefChanged) && !hudSettingsEquivalent(hudSettings, hudSettingsSaved)) {
      (void)ui::saveHudSettingsToFile(hudSettings, hudSettingsPath);
    }
  }

  // Persist controls on exit (optional).
  if (controlsAutoSaveOnExit && controlsDirty) {
    game::saveToFile(controls, controlsPath);
  }

  // Persist bookmarks on exit (optional).
  if (bookmarksAutoSaveOnExit && bookmarksDirty) {
    (void)ui::saveBookmarksToFile(bookmarks, bookmarksPath);
  }

  // Persist UI settings on exit (optional).
  if (uiSettingsAutoSaveOnExit && uiSettingsDirty) {
    syncUiSettingsFromRuntime();
    (void)ui::saveUiSettingsToFile(uiSettings, uiSettingsPath);
  }

  // Unhook core log forwarding before shutdown.
  game::consoleRemoveCoreLogSink(consoleWindow);

  // Cleanup
  spriteCache.clear();
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplSDL2_Shutdown();
  ImGui::DestroyContext();

  SDL_GL_DeleteContext(glContext);
  SDL_DestroyWindow(window);
  SDL_Quit();

  return 0;
}
