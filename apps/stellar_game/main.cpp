// On Windows, SDL may redefine main() to SDL_main unless SDL2main is linked.
// We provide our own entry point, so prevent SDL from overriding it.
#ifdef _WIN32
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
#include "stellar/ui/Livery.h"
#include "stellar/sim/Orbit.h"
#include "stellar/sim/MissionLogic.h"
#include "stellar/sim/Contraband.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/Ship.h"
#include "stellar/sim/Traffic.h"
#include "stellar/sim/NavRoute.h"
#include "stellar/sim/Universe.h"

#include <SDL.h>
#include <SDL_opengl.h>

// Some SDL configurations still define main -> SDL_main; ensure our main symbol remains intact.
#ifdef main
#undef main
#endif

#include <imgui.h>
#include "stellar/ui/ImGuiCompat.h"
#include <backends/imgui_impl_opengl3.h>
#include <backends/imgui_impl_sdl2.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cctype>
#include <cstdint>
#include <limits>
#include <optional>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace stellar;

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
  const math::Vec3d ab = bKm - aKm;
  const double abLenSq = ab.lengthSq();
  if (abLenSq < 1e-12) {
    return (aKm - centerKm).lengthSq() <= radiusKm * radiusKm;
  }

  const double t = std::clamp(math::dot(centerKm - aKm, ab) / abLenSq, 0.0, 1.0);
  const math::Vec3d closest = aKm + ab * t;
  return (closest - centerKm).lengthSq() <= radiusKm * radiusKm;
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
  double ttl{3.0}; // seconds remaining
};

static void toast(std::vector<ToastMsg>& toasts, std::string msg, double ttlSec=3.0) {
  toasts.push_back({std::move(msg), ttlSec});
}

struct ClearanceState {
  bool granted{false};
  double expiresDays{0.0};
  double cooldownUntilDays{0.0};
};

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

  // Short "under fire" window for basic behavior tuning (e.g., morale / evasion).
  double underFireUntilDays{0.0};

  // Behavior flags
  bool hostileToPlayer{false}; // police become hostile when you're wanted / shoot them
  double fleeUntilDays{0.0};   // traders flee after being attacked

  // Combat stats (very lightweight for now, but enough for early combat/mining loops)
  double shield{60.0};
  double shieldMax{60.0};
  double shieldRegenPerSec{0.0}; // points per simulated second

  double hull{70.0};
  double hullMax{70.0};

  // Traders can have an approximate loot value (paid on destruction for now).
  double cargoValueCr{0.0};

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

  double fireCooldown{0.0}; // seconds
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

  // Whether the player has "resolved" this site (i.e., we have fired its one-shot content).
  bool resolved{false};

  // Some sites (resource fields) need a second one-shot flag even after 'resolved' for clarity / future expansion.
  bool fieldSpawned{false};
};

static void applyDamage(double dmg, double& shield, double& hull) {
  if (shield > 0.0) {
    const double s = std::min(shield, dmg);
    shield -= s;
    dmg -= s;
  }
  if (dmg > 0.0) {
    hull -= dmg;
  }
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

static bool insideStationHullExceptSlot(const sim::Station& st,
                                       const math::Vec3d& relLocalKm) {
  // Local station hull approximation: box (wx,wy,wz), with a rectangular slot tunnel cut out.
  const double wx = st.radiusKm * 0.70;
  const double wy = st.radiusKm * 0.70;
  const double wz = st.radiusKm * 1.10;

  const bool insideBox = (std::abs(relLocalKm.x) < wx) && (std::abs(relLocalKm.y) < wy) && (std::abs(relLocalKm.z) < wz);
  if (!insideBox) return false;

  // Slot tunnel cutout near +Z face (entrance at +wz)
  const double slotHalfW = st.slotWidthKm * 0.5;
  const double slotHalfH = st.slotHeightKm * 0.5;

  const double zEntrance = wz;
  const double zMin = zEntrance - st.slotDepthKm;

  const bool insideTunnel =
    (std::abs(relLocalKm.x) < slotHalfW) &&
    (std::abs(relLocalKm.y) < slotHalfH) &&
    (relLocalKm.z <= zEntrance) &&
    (relLocalKm.z >= zMin);

  return !insideTunnel;
}

static bool dockingSlotConditions(const sim::Station& st,
                                 const math::Vec3d& relLocalKm,
                                 const math::Vec3d& shipVelRelLocalKmS,
                                 const math::Vec3d& shipForwardLocal,
                                 const math::Vec3d& shipUpLocal,
                                 bool clearanceGranted) {
  if (!clearanceGranted) return false;

  const double wx = st.radiusKm * 0.70;
  const double wy = st.radiusKm * 0.70;
  const double wz = st.radiusKm * 1.10;
  const double zEntrance = wz;
  const double zMin = zEntrance - st.slotDepthKm;

  // Must be inside tunnel volume (i.e. have flown into the mail-slot)
  const double slotHalfW = st.slotWidthKm * 0.5;
  const double slotHalfH = st.slotHeightKm * 0.5;

  const bool insideTunnel =
    (std::abs(relLocalKm.x) < slotHalfW) &&
    (std::abs(relLocalKm.y) < slotHalfH) &&
    (relLocalKm.z <= zEntrance - 0.05 * st.radiusKm) &&
    (relLocalKm.z >= zMin + 0.10 * st.radiusKm);

  if (!insideTunnel) return false;

  const double relSpeed = shipVelRelLocalKmS.length();
  if (relSpeed > st.maxApproachSpeedKmS) return false;

  // Orientation: ship forward should point into the station (-Z local),
  // and ship up should be roughly +Y local (roll alignment).
  const double fwdAlign = math::dot(shipForwardLocal.normalized(), math::Vec3d{0,0,-1});
  const double upAlign = math::dot(shipUpLocal.normalized(), math::Vec3d{0,1,0});
  return (fwdAlign > 0.92 && upAlign > 0.70);
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
  // Note: Docking is only available in Dear ImGui's "docking" branch.
  // The project pins ImGui's mainline tag by default, so keep docking optional.
  // (If you later switch to the docking branch, feel free to re-enable this.)
  // io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
  io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
  ImGui::StyleColorsDark();

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

  // --- Ship loadout / progression (early, lightweight) ---
  enum class ShipHullClass : int {
    Scout = 0,
    Hauler,
    Fighter,
  };

  enum class WeaponType : int {
    BeamLaser = 0,
    PulseLaser,
    Cannon,
    Railgun,
    MiningLaser,
  };

  struct HullDef {
    const char* name;
    double priceCr;
    double hullMax;
    double shieldBase;
    double linAccelKmS2;
    double angAccelRadS2;
    double cargoMult;
    double fuelMult;
  };

  static const HullDef kHullDefs[] = {
    // name,  price, hull, shield, linAccel, angAccel, cargoMult, fuelMult
    {"Scout",  0.0,   100.0, 100.0, 0.080,  1.20,     1.00,      1.00},
    {"Hauler", 14000.0, 125.0,  90.0, 0.070,  1.05,     1.35,      1.25},
    {"Fighter", 22000.0,  95.0, 120.0, 0.095,  1.35,     0.90,      1.05},
  };

  struct MkDef {
    const char* name;
    double priceCr;
    double mult;
  };

  static const MkDef kThrusters[] = {
    {"(invalid)", 0.0, 1.0},
    {"Thrusters Mk1", 0.0, 1.00},
    {"Thrusters Mk2", 8500.0, 1.18},
    {"Thrusters Mk3", 17500.0, 1.35},
  };

  static const MkDef kShields[] = {
    {"(invalid)", 0.0, 1.0},
    {"Shields Mk1", 0.0, 1.00},
    {"Shields Mk2", 9000.0, 1.25},
    {"Shields Mk3", 18500.0, 1.55},
  };

  static const MkDef kDistributors[] = {
    {"(invalid)", 0.0, 1.0},
    {"Distributor Mk1", 0.0, 1.00},
    {"Distributor Mk2", 7500.0, 1.22},
    {"Distributor Mk3", 16000.0, 1.45},
  };

  struct WeaponDef {
    WeaponType type;
    const char* name;
    double priceCr;
    double cooldownSimSec;
    double heatPerShot;
    double dmg;
    double rangeKm;
    double projSpeedKmS; // ignored for beam
    bool beam;
    float r, g, b;
  };

  static const WeaponDef kWeaponDefs[] = {
    {WeaponType::BeamLaser, "Beam Laser", 0.0, 0.18, 2.2, 7.5, 210000.0, 0.0, true, 1.00f, 0.35f, 0.15f},
    {WeaponType::PulseLaser, "Pulse Laser", 5200.0, 0.28, 1.9, 6.0, 230000.0, 0.0, true, 1.00f, 0.80f, 0.20f},
    {WeaponType::Cannon, "Cannon", 0.0, 0.90, 4.5, 22.0, 260000.0, 120.0, false, 1.00f, 1.00f, 0.90f},
    {WeaponType::Railgun, "Railgun", 9800.0, 1.65, 7.5, 45.0, 320000.0, 240.0, false, 0.60f, 0.90f, 1.00f},
    {WeaponType::MiningLaser, "Mining Laser", 6500.0, 0.22, 1.6, 3.0, 180000.0, 0.0, true, 0.30f, 1.00f, 0.35f},
  };

  auto weaponDef = [&](WeaponType t) -> const WeaponDef& {
    return kWeaponDefs[(int)t];
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

  // Baseline NPC combat tuning (used for random encounters / ambushes)
  const double npcShieldMax = 80.0;
  const double npcHullMax = 90.0;
  const double npcShieldRegenPerSec = 0.035; // points per simulated second (~2.1 per sim minute)

  auto recalcPlayerStats = [&]() {
    const auto& hull = kHullDefs[(int)shipHullClass];

    const int tMk = std::clamp(thrusterMk, 1, 3);
    const int sMk = std::clamp(shieldMk, 1, 3);
    const int dMk = std::clamp(distributorMk, 1, 3);

    const double thrMult = kThrusters[tMk].mult;
    const double shMult  = kShields[sMk].mult;
    const double distMult = kDistributors[dMk].mult;

    playerHullMax = hull.hullMax;
    playerShieldMax = hull.shieldBase * shMult;

    playerBaseLinAccelKmS2 = hull.linAccelKmS2 * thrMult;
    playerBaseAngAccelRadS2 = hull.angAccelRadS2 * (0.92 + 0.08 * thrMult);

    // Regen & cooling scale mostly with distributor. Shields help slightly.
    playerShieldRegenPerSimMin = 2.5 * (0.85 + 0.15 * (double)sMk) * (0.85 + 0.15 * (double)dMk);
    playerHeatCoolRate = 10.0 * distMult;

    ship.setMaxLinearAccelKmS2(playerBaseLinAccelKmS2);
    ship.setMaxAngularAccelRadS2(playerBaseAngAccelRadS2);

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

  // Missions + reputation
  core::u64 nextMissionId = 1;
  std::vector<sim::Mission> missions;
  // Mission UI convenience: which mission is currently tracked in the HUD.
  core::u64 trackedMissionId = 0;
  bool missionTrackerAutoPlotNextLeg = true;
  bool objectiveHudEnabled = true; // overlay showing tracked mission objective in-flight
  std::unordered_map<core::u32, double> repByFaction;

  // Smuggling / contraband legality (per faction, deterministic bitmask of illegal commodities)
  std::unordered_map<core::u32, core::u32> illegalMaskByFaction;

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
    econ::CommodityId commodity{econ::CommodityId::Food};
    double buyPrice{0.0};
    double sellPrice{0.0};
    double profitPerUnit{0.0};
    double profitPerKg{0.0};
    double estTripProfit{0.0};
    double distanceLy{0.0};
  };

  std::vector<TradeIdea> tradeIdeas;
  sim::StationId tradeFromStationId = 0;
  int tradeIdeasDayStamp = -1;
  double tradeSearchRadiusLy = 200.0;

  // Docking state
  bool docked = false;
  sim::StationId dockedStationId = 0;
  std::unordered_map<sim::StationId, ClearanceState> clearances;

  // Contacts (pirates etc.)
  std::vector<Contact> contacts;
  core::SplitMix64 rng(seed ^ 0xC0FFEEu);
  double nextPirateSpawnDays = 0.01; // soon after start
  double nextTraderSpawnDays = 0.008;
  double nextPoliceSpawnDays = 0.006;

  // Salvage / mining / signal sources
  bool cargoScoopDeployed = true;
  double cargoFullToastCooldownUntilDays = 0.0;
  std::vector<FloatingCargo> floatingCargo;
  std::vector<AsteroidNode> asteroids;
  std::vector<SignalSource> signals;
  core::u64 nextWorldObjectId = 1;
  double nextSignalSpawnDays = 0.01;
  double nextSupercruiseShadowSpawnDays = 0.01;

  // Beams (for laser visuals)
  struct Beam { math::Vec3d aU, bU; float r,g,b; double ttl; };
  std::vector<Beam> beams;

  // Projectiles (kinetic cannons / slugs)
  struct Projectile {
    math::Vec3d prevKm;
    math::Vec3d posKm;
    math::Vec3d velKmS;
    float r{1}, g{1}, b{1};
    double ttlSim{0.0};      // simulated seconds remaining
    double radiusKm{450.0};  // collision radius
    double dmg{0.0};
    bool fromPlayer{false};
    core::u64 shooterId{0};
  };
  std::vector<Projectile> projectiles;

  // Save/load
  const std::string savePath = "savegame.txt";

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

  // Draw an edge-of-screen arrow for the current target when it is off-screen.
  bool hudOffscreenTargetIndicator = true;

  // Combat HUD symbology (procedural reticle + weapon rings + lead indicator + flight path marker).
  bool hudCombatHud = true;
  bool hudUseProceduralReticle = true;
  bool hudShowWeaponRings = true;
  bool hudShowHeatRing = true;
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

  // Optional mouse steering (relative mouse mode). Toggle with M.
  bool mouseSteer = false;
  float mouseSensitivity = 0.0025f; // torque intent per pixel
  bool mouseInvertY = false;

  // Flight assistance
  bool autopilot = false;
  int autopilotPhase = 0; // 0=staging, 1=corridor

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

  // Interdiction (very early minigame while in supercruise)
  enum class InterdictionState { None, Warning, Active };
  InterdictionState interdictionState = InterdictionState::None;
  double interdictionWarningRemainingSec = 0.0;
  double interdictionActiveRemainingSec = 0.0;
  double interdictionEscapeMeter = 0.0; // 0..1
  math::Vec3d interdictionEscapeDir{0,0,1};
  core::u64 interdictionPirateId = 0;
  std::string interdictionPirateName;
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

    // Always seed at least one resource field near a "useful" station so mining is easy to find.
    if (sys.stations.empty()) return;
    std::size_t anchorIdx = 0;
    for (std::size_t i = 0; i < sys.stations.size(); ++i) {
      const auto t = sys.stations[i].type;
      if (t == econ::StationType::Mining || t == econ::StationType::Refinery) { anchorIdx = i; break; }
    }

    const auto& anchor = sys.stations[anchorIdx];
    const math::Vec3d anchorPos = stationPosKm(anchor, timeDays);

    auto spawnSignal = [&](SignalType type, const math::Vec3d& posKm, double ttlDays) {
      SignalSource s{};
      s.id = allocWorldId();
      s.type = type;
      s.posKm = posKm;
      s.expireDay = timeDays + ttlDays;
      s.resolved = false;
      signals.push_back(s);
    };

    auto spawnField = [&](const math::Vec3d& centerKm, int count, econ::CommodityId yield) {
      for (int i = 0; i < count; ++i) {
        AsteroidNode a{};
        a.id = allocWorldId();
        a.posKm = centerKm + randUnit() * rng.range(15000.0, 75000.0);
        a.radiusKm = rng.range(2500.0, 7500.0);
        a.yield = yield;
        a.remainingUnits = rng.range(90.0, 260.0) * (a.radiusKm / 5000.0);
        a.chunkAccumulator = 0.0;
        asteroids.push_back(a);
      }
    };

    // A permanent-ish resource field.
    const math::Vec3d resourcePos = anchorPos + randUnit() * (anchor.commsRangeKm * 1.3 + 120000.0);
    spawnSignal(SignalType::Resource, resourcePos, 5.0); // lasts several days
    spawnField(resourcePos, 28, econ::CommodityId::Ore);

    // A derelict somewhere nearby for early salvage.
    const math::Vec3d derelictPos = anchorPos + randUnit() * (anchor.commsRangeKm * 1.6 + 190000.0);
    spawnSignal(SignalType::Derelict, derelictPos, 1.0);

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
  nextPoliceSpawnDays = std::min(nextPoliceSpawnDays, timeDays + (respSec / 86400.0));
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

// Apply confiscation + fine + reputation/bounty/alert effects when contraband is enforced.
// `scannedIllegal` is what security saw during the scan (used for messaging + smuggle mission attribution).
auto enforceContraband = [&](core::u32 jurisdiction,
                             const std::string& sourceName,
                             double illegalValueCr,
                             const std::string& detail,
                             const std::array<double, econ::kCommodityCount>& scannedIllegal) {
  if (jurisdiction == 0) return;

  // Confiscate the scanned illegal goods (or whatever remains of them).
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const double seen = scannedIllegal[i];
    if (seen <= 1e-6) continue;
    const double take = std::min(cargo[i], seen);
    if (take > 1e-6) {
      cargo[i] = std::max(0.0, cargo[i] - take);
    }
  }

  const double fineCr = 200.0 + illegalValueCr * 0.65;
  const double paidCr = std::min(credits, fineCr);
  credits -= paidCr;
  const double unpaidCr = fineCr - paidCr;

  const double repPenalty = -std::clamp(2.0 + illegalValueCr / 3000.0, 2.0, 10.0);
  addRep(jurisdiction, repPenalty);

  if (unpaidCr > 1e-6) {
    addBounty(jurisdiction, unpaidCr);
  }

  // Heightened police attention for a short window.
  policeAlertUntilDays = std::max(policeAlertUntilDays, timeDays + (90.0 / 86400.0));
  nextPoliceSpawnDays = std::min(nextPoliceSpawnDays, timeDays + (6.0 / 86400.0));

  // Raise local security alert (affects patrol spawn rate / response).
  policeHeat = std::clamp(policeHeat + std::min(1.8, 0.5 + illegalValueCr / 5000.0), 0.0, 6.0);

  std::string msg = "Contraband detected! Confiscated: "
                    + (detail.empty() ? std::string("illegal cargo") : detail)
                    + ". Fine: " + std::to_string((int)std::round(fineCr)) + " cr"
                    + (unpaidCr > 1e-6 ? " (unpaid -> bounty)" : "")
                    + ". Rep " + std::to_string((int)std::round(repPenalty));

  if (unpaidCr > 1e-6 && !docked) {
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
        if (event.key.keysym.sym == SDLK_ESCAPE) running = false;

        const bool ctrlDown = (event.key.keysym.mod & KMOD_CTRL) != 0;

        if (event.key.keysym.sym == SDLK_F5) {
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
          s.fsdReadyDay = fsdReadyDay;

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
          s.trackedMissionId = trackedMissionId;

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

            nextMissionId = s.nextMissionId;
            missions = s.missions;
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

	            bountyByFaction.clear();
	            for (const auto& b : s.bounties) bountyByFaction[b.factionId] = b.bountyCr;
	            bountyVoucherByFaction.clear();
	            for (const auto& v : s.bountyVouchers) bountyVoucherByFaction[v.factionId] = v.bountyCr;

                // Background traffic stamps
                trafficDayStampBySystem.clear();
                for (const auto& t : s.trafficStamps) {
                  trafficDayStampBySystem[t.systemId] = t.dayStamp;
                }
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
            beams.clear();
            projectiles.clear();

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

            nextPirateSpawnDays = timeDays + (rng.range(25.0, 55.0) / 86400.0);
            nextTraderSpawnDays = timeDays + (rng.range(15.0, 35.0) / 86400.0);
            nextPoliceSpawnDays = timeDays + (rng.range(10.0, 25.0) / 86400.0);
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
            navRoute.clear();
            navRouteHop = 0;
            navAutoRun = false;
            scanning = false;
            scanProgressSec = 0.0;
            scanLockedId = 0;
            scanLabel.clear();
            scanLockedTarget = Target{};
            clearances.clear();
            target = Target{};

            galaxySelectedSystemId = currentSystem->stub.id;

            toast(toasts, "Loaded " + savePath, 2.5);
          }
        }

        if (event.key.keysym.sym == SDLK_TAB) showGalaxy = !showGalaxy;
        if (event.key.keysym.sym == SDLK_F1) showShip = !showShip;
        if (event.key.keysym.sym == SDLK_F2) showEconomy = !showEconomy;
        if (event.key.keysym.sym == SDLK_F4) showMissions = !showMissions;
        if (event.key.keysym.sym == SDLK_F3) showContacts = !showContacts;
        if (event.key.keysym.sym == SDLK_F6) showScanner = !showScanner;
        if (event.key.keysym.sym == SDLK_F7) showTrade = !showTrade;
        if (event.key.keysym.sym == SDLK_F8) showGuide = !showGuide;
        if (event.key.keysym.sym == SDLK_F10) showSprites = !showSprites;
        if (event.key.keysym.sym == SDLK_F11) showVfx = !showVfx;
        if (event.key.keysym.sym == SDLK_F12) showPostFx = !showPostFx;

        // HUD layout editor shortcuts
        // Ctrl+H : toggle edit mode
        // Ctrl+S : save HUD layout
        // Ctrl+L : load HUD layout
        // Ctrl+R : reset HUD layout to defaults
        if (ctrlDown && !io.WantCaptureKeyboard) {
          if (event.key.keysym.sym == SDLK_h) {
            hudLayoutEditMode = !hudLayoutEditMode;
            showHudLayoutWindow = true;
            toast(toasts,
                  std::string("HUD layout edit ") + (hudLayoutEditMode ? "ON" : "OFF") + " (Ctrl+H)",
                  1.8);
          }

          if (event.key.keysym.sym == SDLK_s) {
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

          if (event.key.keysym.sym == SDLK_l) {
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

          if (event.key.keysym.sym == SDLK_r) {
            hudLayout = ui::makeDefaultHudLayout();
            showRadarHud = hudLayout.widget(ui::HudWidgetId::Radar).enabled;
            objectiveHudEnabled = hudLayout.widget(ui::HudWidgetId::Objective).enabled;
            hudThreatOverlayEnabled = hudLayout.widget(ui::HudWidgetId::Threat).enabled;
            hudJumpOverlay = hudLayout.widget(ui::HudWidgetId::Jump).enabled;
            toast(toasts, "HUD layout reset to defaults.", 2.0);
            showHudLayoutWindow = true;
          }
        }

        if (event.key.keysym.sym == SDLK_r) {
          if (!ctrlDown && !io.WantCaptureKeyboard) {
            showRadarHud = !showRadarHud;
            toast(toasts, std::string("Radar HUD ") + (showRadarHud ? "ON" : "OFF") + " (R)", 1.6);
          }
        }

        if (event.key.keysym.sym == SDLK_BACKQUOTE) {
          if (!io.WantCaptureKeyboard) {
            showTacticalOverlay = !showTacticalOverlay;
            toast(toasts, std::string("Tactical overlay ") + (showTacticalOverlay ? "ON" : "OFF") + " (`)", 1.6);
          }
        }

        if (event.key.keysym.sym == SDLK_SPACE) paused = !paused;

        if (event.key.keysym.sym == SDLK_p) {
          autopilot = !autopilot;
          autopilotPhase = 0;
        }

	        if (event.key.keysym.sym == SDLK_o) {
	          if (!io.WantCaptureKeyboard) {
	            cargoScoopDeployed = !cargoScoopDeployed;
	            toast(toasts,
	                  std::string("Cargo scoop ") + (cargoScoopDeployed ? "DEPLOYED" : "RETRACTED") + " (O)",
	                  1.8);
	          }
	        }

	        if (event.key.keysym.sym == SDLK_v) {
	          if (!io.WantCaptureKeyboard && !docked && currentSystem) {
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

	        if (event.key.keysym.sym == SDLK_i) {
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

        if (event.key.keysym.sym == SDLK_c) {
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

if (event.key.keysym.sym == SDLK_m) {
          if (!io.WantCaptureKeyboard) {
            mouseSteer = !mouseSteer;
            SDL_SetRelativeMouseMode(mouseSteer ? SDL_TRUE : SDL_FALSE);
            SDL_GetRelativeMouseState(nullptr, nullptr); // flush deltas
            toast(toasts, std::string("Mouse steer ") + (mouseSteer ? "ON" : "OFF") + " (M)", 1.8);
          }
        }

        if (event.key.keysym.sym == SDLK_h) {
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
                interdictionState = InterdictionState::None;
                interdictionWarningRemainingSec = 0.0;
                interdictionActiveRemainingSec = 0.0;
                interdictionEscapeMeter = 0.0;
                interdictionSubmitRequested = false;
                autopilot = false;
                autopilotPhase = 0;
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
              interdictionState = InterdictionState::None;
              interdictionSubmitRequested = false;
              toast(toasts, "Supercruise canceled.", 1.4);
            } else if (supercruiseState == SupercruiseState::Active) {
              if (interdictionState != InterdictionState::None) {
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

        if (event.key.keysym.sym == SDLK_k) {
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
                toast(toasts, "No scan target selected. (Use T/B/N/U or System Scanner)", 2.5);
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
	                  scanDurationSec = 4.0;
	                  scanRangeKm = 120000.0;
	                  scanLabel = std::string("Signal scan: ") + signalTypeName(s.type);
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

        if (event.key.keysym.sym == SDLK_j) {
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
if (event.key.keysym.sym == SDLK_t) {
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

if (event.key.keysym.sym == SDLK_b) {
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

if (event.key.keysym.sym == SDLK_n) {
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

if (event.key.keysym.sym == SDLK_u) {
  // target system primary star (for scanning / nav reference)
  if (scanning) { scanning = false; scanProgressSec = 0.0; scanLockedId = 0; scanLabel.clear(); }
  target.kind = TargetKind::Star;
  target.index = 0;
}

if (event.key.keysym.sym == SDLK_y) {
  if (scanning) { scanning = false; scanProgressSec = 0.0; scanLockedId = 0; scanLabel.clear(); }
  target = Target{};
}

        if (event.key.keysym.sym == SDLK_l) {
          // Request docking clearance from targeted station
          if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
            const auto& st = currentSystem->stations[target.index];
            const math::Vec3d stPos = stationPosKm(st, timeDays);
            const double dist = (ship.positionKm() - stPos).length();

            auto& cs = clearances[st.id];

            if (dist > st.commsRangeKm) {
              toast(toasts, "Out of comms range for clearance.", 2.5);
            } else if (timeDays < cs.cooldownUntilDays) {
              toast(toasts, "Clearance channel busy. Try again soon.", 2.5);
            } else {
              // Simple logic: usually granted; sometimes denied (simulate traffic / capacity).
              const double pGrant = 0.82;
              cs.granted = rng.chance(pGrant);
              if (cs.granted) {
                cs.expiresDays = timeDays + (12.0 * 60.0) / 86400.0; // 12 minutes
                toast(toasts, "Docking clearance GRANTED.", 3.0);
              } else {
                cs.cooldownUntilDays = timeDays + (90.0) / 86400.0; // 90 seconds
                toast(toasts, "Docking clearance DENIED. (traffic)", 3.0);
              }
            }
          } else {
            toast(toasts, "No station targeted for clearance.", 2.0);
          }
        }

        if (event.key.keysym.sym == SDLK_g) {
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

              docked = false;
              dockedStationId = 0;
              toast(toasts, "Undocked.", 2.0);
            } else {
              docked = false;
              dockedStationId = 0;
            }
          } else {
            // Dock attempt: must be inside slot tunnel, aligned, and have clearance.
            if (target.kind == TargetKind::Station && target.index < currentSystem->stations.size()) {
              const auto& st = currentSystem->stations[target.index];
              auto& cs = clearances[st.id];
              const bool clearanceValid = cs.granted && (timeDays <= cs.expiresDays);

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
              } else if (!dockingSlotConditions(st, relLocalKm, relVelLocal, fwdLocal, upLocal, clearanceValid)) {
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
            // NOTE: cooldown is in *sim* seconds, so it naturally scales with timeScale.
            cd = w.cooldownSimSec;

            // Track last fired weapon for HUD lead indicator.
            hudLastFiredPrimary = primary;

            // Better distributors run cooler (a small benefit, not a hard gate).
            const int dMk = std::clamp(distributorMk, 1, 3);
            const double heatFactor = std::clamp(1.08 - 0.08 * (double)dMk, 0.80, 1.10);
            heat = std::min(120.0, heat + w.heatPerShot * heatFactor);

            if (w.beam) {
              const double rangeKm = w.rangeKm;
              const double dmg = w.dmg;

              const math::Vec3d aKm = ship.positionKm();
              const math::Vec3d dir = ship.forward().normalized();

	              // Find best hit (contacts: simple cone + nearest along ray, asteroids: sphere intersection).
	              int bestIdx = -1;
	              double bestT = rangeKm;
	              const double cone = (wType == WeaponType::PulseLaser) ? 0.993 : 0.995;

              for (int i = 0; i < (int)contacts.size(); ++i) {
                auto& c = contacts[(std::size_t)i];
                if (!c.alive) continue;

                const math::Vec3d to = c.ship.positionKm() - aKm;
                const double dist = to.length();
                if (dist > rangeKm) continue;

                const math::Vec3d toN = to / std::max(1e-9, dist);
                const double aim = math::dot(dir, toN);
                if (aim < cone) continue;

                const double t = dist * aim;
                if (t < bestT) { bestT = t; bestIdx = i; }
              }

	              // Asteroid hit test
	              int bestAstIdx = -1;
	              double bestAstT = rangeKm;
	              for (int i = 0; i < (int)asteroids.size(); ++i) {
	                const auto& a = asteroids[(std::size_t)i];
	                const math::Vec3d toCenter = a.posKm - aKm;
	                const double tProj = math::dot(toCenter, dir);
	                if (tProj < 0.0 || tProj > rangeKm) continue;
	                const math::Vec3d closest = aKm + dir * tProj;
	                const double d2 = (a.posKm - closest).lengthSq();
	                const double r2 = a.radiusKm * a.radiusKm;
	                if (d2 > r2) continue;
	                const double thc = std::sqrt(std::max(0.0, r2 - d2));
	                const double tHit = std::clamp(tProj - thc, 0.0, rangeKm);
	                if (tHit < bestAstT) { bestAstT = tHit; bestAstIdx = i; }
	              }

	              // Choose nearest along the ray.
	              const bool hitAsteroidFirst = (bestAstIdx >= 0 && bestAstT <= bestT);
	              const double finalT = hitAsteroidFirst ? bestAstT : bestT;
	              const math::Vec3d bKm = aKm + dir * finalT;
	              beams.push_back({toRenderU(aKm), toRenderU(bKm), w.r, w.g, w.b, 0.10});

	              if (hitAsteroidFirst) {
	                if (wType == WeaponType::MiningLaser) {
	                  auto& a = asteroids[(std::size_t)bestAstIdx];
	                  if (a.remainingUnits > 1e-6) {
	                    const double chunk = std::min(a.remainingUnits, rng.range(6.0, 14.0));
	                    a.remainingUnits -= chunk;
	                    spawnCargoPod(a.yield, chunk, bKm, math::Vec3d{0,0,0}, 0.45);
	                    toast(toasts,
	                          std::string("Mined ") + econ::commodityDef(a.yield).name + " x" + std::to_string((int)std::round(chunk)),
	                          2.0);
	                  } else {
	                    toast(toasts, "Asteroid depleted.", 1.5);
	                  }
	                }
	              } else {
	                if (bestIdx >= 0) playerDamageContact(bestIdx, dmg);
	              }
            } else {
              const double muzzleSpeedKmS = std::max(1e-6, w.projSpeedKmS);
              const double rangeKm = w.rangeKm;
              const double dmg = w.dmg;

              const double ttlSim = rangeKm / muzzleSpeedKmS;

              const math::Vec3d fwd = ship.forward().normalized();
              const math::Vec3d spawnKm = ship.positionKm() + fwd * 400.0;

              Projectile p{};
              p.prevKm = spawnKm;
              p.posKm = spawnKm;
              p.velKmS = ship.velocityKmS() + fwd * muzzleSpeedKmS;
              p.r = w.r; p.g = w.g; p.b = w.b;
              p.ttlSim = ttlSim;
              p.radiusKm = (wType == WeaponType::Railgun) ? 520.0 : 700.0;
              p.dmg = dmg;
              p.fromPlayer = true;
              p.shooterId = 0;
              projectiles.push_back(p);

              // Tiny recoil impulse
              ship.setVelocityKmS(ship.velocityKmS() - fwd * (wType == WeaponType::Railgun ? 0.003 : 0.002));
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

      input.boost = keys[SDL_SCANCODE_LSHIFT] != 0;
      input.brake = keys[SDL_SCANCODE_X] != 0;

      static bool dampers = true;
      if (keys[SDL_SCANCODE_Z]) dampers = true;
      if (keys[SDL_SCANCODE_C]) dampers = false;
      input.dampers = dampers;
    }

    // Autopilot: station approach assist (staging + docking corridor guidance).
    if (autopilot && !docked && currentSystem && target.kind == TargetKind::Station && target.index < currentSystem->stations.size() && !captureKeys) {
      const auto& st = currentSystem->stations[target.index];

      const math::Vec3d stPos = stationPosKm(st, timeDays);
      const math::Quatd stQ = stationOrient(st, stPos, timeDays);
      const math::Vec3d stV = stationVelKmS(st, timeDays);

      const double zEntrance = st.radiusKm * 1.10;
      const double zCorridorEnd = zEntrance + st.approachLengthKm;
      const double zStage = zCorridorEnd + st.radiusKm * 0.8;

      const math::Vec3d relLocal = stQ.conjugate().rotate(ship.positionKm() - stPos);
      const double lateralDist = math::Vec3d{relLocal.x, relLocal.y, 0.0}.length();

      // Phase switching: start by staging at the corridor end, then commit through the corridor.
      if (autopilotPhase == 0) {
        if (relLocal.z < zCorridorEnd * 1.05 && lateralDist < st.approachRadiusKm * 0.65) autopilotPhase = 1;
      } else {
        // If we drift far off-axis, go back to staging.
        if (relLocal.z > zStage || lateralDist > st.approachRadiusKm * 1.25) autopilotPhase = 0;
      }

      math::Vec3d desiredLocal{0,0,0};
      if (autopilotPhase == 0) desiredLocal = {0,0, zStage};
      else desiredLocal = {0,0, zEntrance + st.radiusKm * 0.7};

      const math::Vec3d desiredPoint = stPos + stQ.rotate(desiredLocal);
      const math::Vec3d rel = desiredPoint - ship.positionKm();
      const double dist = rel.length();
      const math::Vec3d dir = (dist > 1e-6) ? (rel / dist) : math::Vec3d{0,0,0};

      const double maxV = (autopilotPhase == 0) ? std::max(st.maxApproachSpeedKmS * 2.5, 0.45)
                                                : (st.maxApproachSpeedKmS * 0.85);
      double vMag = std::min(maxV, 0.004 * dist);

      // Slow down more if we're not centered when committing.
      if (autopilotPhase == 1) {
        const double frac = std::clamp(lateralDist / std::max(1.0, st.approachRadiusKm), 0.0, 1.0);
        vMag *= (1.0 - 0.55 * frac);
      }

      const math::Vec3d desiredVel = stV + dir * vMag;
      const math::Vec3d dv = desiredVel - ship.velocityKmS();

      if (dv.lengthSq() > 1e-12) input.thrustLocal = ship.orientation().conjugate().rotate(dv.normalized());
      else input.thrustLocal = {0,0,0};

      input.dampers = true;

      // Align forward into station (towards -axisOut).
      const math::Vec3d axisOut = stQ.rotate({0,0,1});
      const math::Vec3d desiredFwdWorld = -axisOut;

      const math::Vec3d desiredFwdLocal = ship.orientation().conjugate().rotate(desiredFwdWorld);
      const double yawErr = std::atan2(desiredFwdLocal.x, desiredFwdLocal.z);
      const double pitchErr = -std::atan2(desiredFwdLocal.y, desiredFwdLocal.z);

      // Roll: keep ship "up" aligned to station +Y so you fly the slot level.
      const math::Vec3d desiredUpWorld = stQ.rotate({0,1,0});
      const math::Vec3d desiredUpLocal = ship.orientation().conjugate().rotate(desiredUpWorld);
      const double rollErr = std::atan2(desiredUpLocal.x, desiredUpLocal.y);

      input.torqueLocal.x = (float)std::clamp(pitchErr * 1.8, -1.0, 1.0);
      input.torqueLocal.y = (float)std::clamp(yawErr * 1.8, -1.0, 1.0);
      input.torqueLocal.z = (float)std::clamp(rollErr * 1.6, -1.0, 1.0);
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
        interdictionState = InterdictionState::None;
        interdictionWarningRemainingSec = 0.0;
        interdictionActiveRemainingSec = 0.0;
        interdictionEscapeMeter = 0.0;
        interdictionSubmitRequested = false;
        interdictionPirateId = 0;
        interdictionPirateName.clear();
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
      interdictionState = InterdictionState::None;
      interdictionWarningRemainingSec = 0.0;
      interdictionActiveRemainingSec = 0.0;
      interdictionEscapeMeter = 0.0;
      interdictionSubmitRequested = false;
      interdictionPirateId = 0;
      interdictionPirateName.clear();

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
        heat = std::min(120.0, heat + 25.0);
        playerShield = std::max(0.0, playerShield - playerShieldMax * 0.12);
        playerHull = std::max(0.0, playerHull - playerHullMax * 0.04);
        ship.setAngularVelocityRadS({rng.range(-0.6, 0.6), rng.range(-0.6, 0.6), rng.range(-0.4, 0.4)});
      } else {
        ship.setAngularVelocityRadS({0,0,0});
      }

      toast(toasts, msg, 2.2);
    };

    if (supercruiseState == SupercruiseState::Active && !docked && !captureKeys && fsdState == FsdState::Idle) {
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
	            p.ship.setMaxLinearAccelKmS2(0.06);
	            p.ship.setMaxAngularAccelRadS2(0.9);
	            p.hostileToPlayer = true;
	            p.supercruiseShadow = true;
	            // Keep them behind you so they can pressure an interdiction.
	            p.shadowOffsetLocalKm = {rng.range(-25000.0, 25000.0), rng.range(-12000.0, 12000.0), -rng.range(165000.0, 220000.0)};
	            p.ship.setPositionKm(ship.positionKm() + ship.orientation().rotate(p.shadowOffsetLocalKm));
	            p.ship.setVelocityKmS(ship.velocityKmS());
	            p.ship.setOrientation(ship.orientation());
	            p.ship.setAngularVelocityRadS({0,0,0});
	            p.shield = npcShieldMax;
	            p.shieldMax = npcShieldMax;
	            p.shieldRegenPerSec = npcShieldRegenPerSec;
	            p.hull = npcHullMax;
	            p.hullMax = npcHullMax;
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
	            pc.ship.setMaxLinearAccelKmS2(0.07);
	            pc.ship.setMaxAngularAccelRadS2(1.0);
	            pc.supercruiseShadow = true;
	            pc.shadowOffsetLocalKm = {rng.range(-22000.0, 22000.0), rng.range(-10000.0, 10000.0), -rng.range(140000.0, 200000.0)};
	            pc.ship.setPositionKm(ship.positionKm() + ship.orientation().rotate(pc.shadowOffsetLocalKm));
	            pc.ship.setVelocityKmS(ship.velocityKmS());
	            pc.ship.setOrientation(ship.orientation());
	            pc.ship.setAngularVelocityRadS({0,0,0});
	            pc.shield = npcShieldMax;
	            pc.shieldMax = npcShieldMax;
	            pc.shieldRegenPerSec = npcShieldRegenPerSec;
	            pc.hull = npcHullMax;
	            pc.hullMax = npcHullMax;
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
	        signals.push_back(s);

	        toast(toasts, std::string("Signal detected: ") + signalTypeName(t), 3.0);
	      }

	      // Interdiction (very early, but clearer + playable than the old "instant drop")
      double nearestPirateKm = 1e99;
      core::u64 nearestPirateId = 0;
      std::string nearestPirateName;
      for (const auto& c : contacts) {
        if (!c.alive || c.role != ContactRole::Pirate) continue;
        const double d = (c.ship.positionKm() - ship.positionKm()).length();
        if (d < nearestPirateKm) {
          nearestPirateKm = d;
          nearestPirateId = c.id;
          nearestPirateName = c.name;
        }
      }
      const bool pirateThreat = (nearestPirateKm < 240000.0);

      // "Submit" (H) while the interdiction is ongoing: safe-ish drop with short cooldown.
      if (interdictionSubmitRequested && interdictionState != InterdictionState::None) {
        interdictionSubmitRequested = false;
        interdictionCooldownUntilDays = timeDays + (90.0 / 86400.0);
        endSupercruise(false, destVelKmS, dirToDest, "Submitted: dropped from supercruise.");
      } else if (!hasDest) {
        endSupercruise(true, localFrameVelKmS, dirToDest, "Supercruise lost target. Emergency drop.");
      } else {
        // Start interdiction warning (random-ish cadence based on proximity)
        if (interdictionState == InterdictionState::None && pirateThreat && timeDays >= interdictionCooldownUntilDays) {
          const double closeness = std::clamp(1.0 - (nearestPirateKm / 240000.0), 0.0, 1.0);
          const double ratePerSec = 0.12 + 0.75 * (closeness * closeness); // events / real second
          const double chance = 1.0 - std::exp(-ratePerSec * dtReal);
          if (rng.chance(std::clamp(chance, 0.0, 0.85))) {
            interdictionState = InterdictionState::Warning;
            interdictionWarningRemainingSec = 2.2;
            interdictionActiveRemainingSec = 0.0;
            interdictionEscapeMeter = 0.42;
            interdictionPirateId = nearestPirateId;
            interdictionPirateName = nearestPirateName;

            // Pick an escape direction that isn't trivially "forward".
            const math::Vec3d fwd = ship.forward().normalized();
            math::Vec3d jitter{rng.range(-1.0, 1.0), rng.range(-0.35, 0.35), rng.range(-1.0, 1.0)};
            if (jitter.lengthSq() < 1e-9) jitter = {0.2, 0.0, 0.8};
            jitter = jitter.normalized();
            interdictionEscapeDir = (fwd * 0.35 + jitter * 0.65).normalized();

            supercruiseDropRequested = false; // disable manual drop while interdicted
            toast(toasts, "INTERDICTION WARNING! Align and hold to evade (H to submit).", 2.6);
          }
        }

        // Update interdiction states
        if (interdictionState == InterdictionState::Warning) {
          interdictionWarningRemainingSec = std::max(0.0, interdictionWarningRemainingSec - dtReal);
          if (interdictionWarningRemainingSec <= 0.0) {
            interdictionState = InterdictionState::Active;
            interdictionActiveRemainingSec = 11.0;
            toast(toasts, "Interdiction engaged! Keep the escape vector aligned.", 2.2);
          }
        }

        if (interdictionState == InterdictionState::Active) {
          interdictionActiveRemainingSec = std::max(0.0, interdictionActiveRemainingSec - dtReal);

          const math::Vec3d fwd = ship.forward().normalized();
          const math::Vec3d esc = interdictionEscapeDir.normalized();
          const double align = math::dot(fwd, esc); // -1..1

          const double gain = std::clamp((align - 0.72) / 0.28, 0.0, 1.0);
          const double pull = 0.26 + (pirateThreat ? std::clamp((240000.0 - nearestPirateKm) / 240000.0, 0.0, 1.0) * 0.22 : 0.0);

          interdictionEscapeMeter = std::clamp(interdictionEscapeMeter + (gain * 0.78 - pull) * dtReal, 0.0, 1.0);

          if (interdictionEscapeMeter >= 1.0) {
            interdictionState = InterdictionState::None;
            interdictionCooldownUntilDays = timeDays + (150.0 / 86400.0);
            toast(toasts, "Interdiction evaded.", 2.0);
          } else if (interdictionEscapeMeter <= 0.0 || interdictionActiveRemainingSec <= 0.0) {
            interdictionCooldownUntilDays = timeDays + (180.0 / 86400.0);
            endSupercruise(true, destVelKmS, dirToDest, "Interdicted! Emergency drop.");
          }
        }

        // Compute safe-drop window ("7-second rule")
        const math::Vec3d rel = destPosKm - ship.positionKm();
        const double dist = rel.length();
        supercruiseDistKm = dist;

        if (dist < 1e-6) {
          endSupercruise(false, destVelKmS, dirToDest, "Supercruise drop.");
        } else {
          const math::Vec3d dir = rel / dist;
          const math::Vec3d vRel = ship.velocityKmS() - destVelKmS;
          const double closing = math::dot(vRel, dir);
          supercruiseClosingKmS = closing;

          const double tta = (closing > 1e-3) ? (dist / closing) : 1e9;
          supercruiseTtaSec = tta;

          const double safeMin = kSupercruiseSafeTtaSec - 2.0;
          const double safeMax = kSupercruiseSafeTtaSec + 2.0;
          const bool safeWindow = (dist < dropKm) && (tta > safeMin) && (tta < safeMax) && (closing > 0.05);
          supercruiseSafeDropReady = safeWindow;

          const bool interdicted = (interdictionState != InterdictionState::None);

          // Manual drop request (H while in supercruise)
          if (!interdicted && supercruiseDropRequested) {
            if (safeWindow) endSupercruise(false, destVelKmS, dir, "Supercruise drop.");
            else endSupercruise(true, destVelKmS, dir, "EMERGENCY DROP!");
          } else if (!interdicted && supercruiseAssist && safeWindow) {
            // Nav assist auto-drop when safe
            endSupercruise(false, destVelKmS, dir, "Supercruise drop.");
          } else {
            // Follow a speed profile toward the target
            const double desiredSpeed =
              supercruiseAssist ? std::clamp(dist / kSupercruiseSafeTtaSec, 60.0, supercruiseMaxSpeedKmS)
                               : std::clamp(dist * 0.0008, 90.0, supercruiseMaxSpeedKmS);

            const math::Vec3d desiredVel = destVelKmS + dir * desiredSpeed;
            const math::Vec3d dv = desiredVel - ship.velocityKmS();

            if (dv.lengthSq() > 1e-12) {
              const math::Vec3d accelDir = dv.normalized();
              input.thrustLocal = ship.orientation().conjugate().rotate(accelDir);
            } else {
              input.thrustLocal = {0,0,0};
            }

            // Face travel direction in normal supercruise. During interdiction, we let the player steer.
            if (!interdicted) {
              const math::Vec3d desiredFwdLocal = ship.orientation().conjugate().rotate(dir);
              const double yawErr = std::atan2(desiredFwdLocal.x, desiredFwdLocal.z);
              const double pitchErr = -std::atan2(desiredFwdLocal.y, desiredFwdLocal.z);

              input.torqueLocal.x = (float)std::clamp(pitchErr * 1.6, -1.0, 1.0);
              input.torqueLocal.y = (float)std::clamp(yawErr * 1.6, -1.0, 1.0);
              input.torqueLocal.z = 0.0f;
            }

            input.dampers = true;
            input.brake = false;
            input.boost = false;

            ship.setMaxLinearAccelKmS2(6.0);
            ship.setMaxAngularAccelRadS2(1.2);
          }
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
      // Also: avoid spamming toasts if the next hop is currently impossible.
      if (navAutoRun && fsdState == FsdState::Idle && timeDays >= fsdReadyDay && !docked && supercruiseState == SupercruiseState::Idle) {
        if (!navRoute.empty() && navRouteHop + 1 < navRoute.size() && !isMassLocked()) {
          const auto& nextSys = universe.getSystem(navRoute[navRouteHop + 1]);
          const double hopDistLy = (nextSys.stub.posLy - currentSystem->stub.posLy).length();

          const double rangeLy = fsdBaseRangeLy();
          if (hopDistLy > rangeLy + 1e-9) {
            navAutoRun = false;
            toast(toasts, "Auto-run stopped: next hop is out of range (replot route).", 2.8);
          } else {
            const double hopFuel = fsdFuelCostFor(hopDistLy);
            if (fuel + 1e-6 < hopFuel) {
              navAutoRun = false;
              toast(toasts, "Auto-run stopped: not enough fuel for next hop.", 2.8);
            } else {
              startFsdJumpTo(navRoute[navRouteHop + 1]);
            }
          }
        }
      }

      // --- Ship physics ---
      if (!docked) {
        if (fsdState != FsdState::Idle || fsdJustArrived) {
          sim::ShipInput hold{};
          hold.dampers = true;
          hold.brake = true;
          ship.step(dtSim, hold);
        } else {
          ship.step(dtSim, input);
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

          heat = std::min(120.0, heat + 25.0);
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


    const double statMul = (k == 0) ? 1.15 : 0.92;
    p.shield = p.shield * statMul;
    p.shieldMax = p.shield;
    p.hull = p.hull * statMul;
    p.hullMax = p.hull;
    p.shieldRegenPerSec = npcShieldRegenPerSec * ((k == 0) ? 0.85 : 0.70);

    p.ship.setMaxLinearAccelKmS2((k == 0) ? 0.065 : 0.058);
    p.ship.setMaxAngularAccelRadS2((k == 0) ? 0.95 : 0.85);

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
    t.shield = 35.0;
    t.shieldMax = t.shield;
    t.hull = 55.0;
    t.hullMax = t.hull;

    // Cargo capacity (kg) varies per trader.
    t.tradeCapacityKg = rng.range(120.0, 420.0);

    // Traders have a high-speed shortcut when far from the player (keeps economy moving).
    t.tradeSupercruiseSpeedKmS = rng.range(7000.0, 14000.0);
    t.tradeSupercruiseDropDistKm = rng.range(70000.0, 140000.0);

    // Assign an initial hauling job (if possible) and load from source inventory.
    planTraderHaul(t, t.homeStationIndex);

    t.ship.setMaxLinearAccelKmS2(0.05);
    t.ship.setMaxAngularAccelRadS2(0.6);
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

    p.shield = 90.0 * tierMul;
    p.shieldMax = p.shield;
    p.shieldRegenPerSec = npcShieldRegenPerSec * tierMul;
    p.hull = 90.0 * tierMul;
    p.hullMax = p.hull;

    p.ship.setMaxLinearAccelKmS2(0.075 * tierMul);
    p.ship.setMaxAngularAccelRadS2(1.05 * tierMul);

    // Wingmen spawn slightly further out to avoid immediate overlaps.
    spawnNearStation(p.homeStationIndex, 16.0 + 1.5 * (double)k, 22.0 + 2.0 * (double)k, p.ship);

    contacts.push_back(std::move(p));
  }
  return count;
};

  // Pirates occasionally (baseline threat).
  if (timeDays >= nextPirateSpawnDays && alivePirates < 4 && aliveTotal < 14) {
    nextPirateSpawnDays = timeDays + (rng.range(120.0, 220.0) / 86400.0); // every ~2-4 minutes
    const int cap = std::max(1, std::min({14 - aliveTotal, 4 - alivePirates, 3}));
    const int spawned = spawnPiratePack(cap);
    alivePirates += spawned;
    aliveTotal += spawned;
  }

  // Traders / traffic: gives you something to pirate (but doing so triggers police).
  if (timeDays >= nextTraderSpawnDays && aliveTraders < 3 && currentSystem && !currentSystem->stations.empty() && aliveTotal < 14) {
    nextTraderSpawnDays = timeDays + (rng.range(70.0, 140.0) / 86400.0);
    spawnTrader();
    ++aliveTraders;
    ++aliveTotal;
  }

  // Police / patrols: scale with threat + your legal status.
  if (localFaction != 0) {
    int desiredPolice = 1;
    if (playerWantedHere) desiredPolice += 2;
    if (localBounty > 1200.0) desiredPolice += 1;
    if (localBounty > 4500.0) desiredPolice += 1;
    if (policeHeat > 3.0) desiredPolice += 1;
    if (alivePirates > 0) desiredPolice += 1;
    if (localRep < -25.0) desiredPolice += 1;
    desiredPolice = std::clamp(desiredPolice, 0, 7);

    // If recently alerted by a crime, tighten the spawn interval.
    const bool alert = (policeAlertUntilDays > timeDays);
    const double heat = policeHeat;
    const double baseMinSec = alert ? 12.0 : 55.0;
    const double baseMaxSec = alert ? 24.0 : 95.0;

    // As heat rises, response gets a bit quicker (but stays bounded).
    const double spawnMinSec = std::clamp(baseMinSec - heat * (alert ? 1.1 : 1.8), 4.0, baseMinSec);
    const double spawnMaxSec = std::clamp(baseMaxSec - heat * (alert ? 1.6 : 2.4), 10.0, baseMaxSec);

    if (timeDays >= nextPoliceSpawnDays && alivePolice < desiredPolice && aliveTotal < 16) {
      nextPoliceSpawnDays = timeDays + (rng.range(spawnMinSec, spawnMaxSec) / 86400.0);
      const int cap = std::max(1, std::min({16 - aliveTotal, desiredPolice - alivePolice, 3}));
      const int spawned = spawnPolicePack(cap);
      alivePolice += spawned;
      aliveTotal += spawned;
    }
  }

  // Ensure mission bounty targets exist in their target system.
  if (!docked && !missions.empty()) {
    for (const auto& m : missions) {
      if (m.completed || m.failed) continue;
      if (!((m.type == sim::MissionType::BountyScan) || (m.type == sim::MissionType::BountyKill))) continue;
      if (m.toSystem != currentSystem->stub.id) continue;
      if (m.targetNpcId == 0) continue;

      bool present = false;
      for (auto& c : contacts) {
        if (c.alive && c.id == m.targetNpcId) {
          c.missionTarget = true;
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
      tgt.shield = 80.0;
      tgt.shieldMax = tgt.shield;
      tgt.shieldRegenPerSec = npcShieldRegenPerSec * 0.9;
      tgt.hull = 90.0;
      tgt.hullMax = tgt.hull;
      tgt.ship.setMaxLinearAccelKmS2(0.07);
      tgt.ship.setMaxAngularAccelRadS2(1.0);

      const math::Vec3d d = randDir();
      const double distKm = rng.range(65000.0, 150000.0);
      tgt.ship.setPositionKm(ship.positionKm() + d * distKm);
      tgt.ship.setVelocityKmS(ship.velocityKmS());
      tgt.ship.setOrientation(quatFromTo({0,0,1}, (-d).normalized()));

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
    ai.dampers = true;

    const math::Vec3d to = targetPosKm - selfShip.positionKm();
    const double dist = to.length();
    const math::Vec3d toN = (dist > 1e-6) ? (to / dist) : math::Vec3d{0,0,1};

    // Face target
    const math::Vec3d desiredFwdWorld = toN;
    const math::Vec3d desiredFwdLocal = selfShip.orientation().conjugate().rotate(desiredFwdWorld);
    const double yawErr = std::atan2(desiredFwdLocal.x, desiredFwdLocal.z);
    const double pitchErr = -std::atan2(desiredFwdLocal.y, desiredFwdLocal.z);
    ai.torqueLocal.x = std::clamp(pitchErr * faceGain, -1.0, 1.0);
    ai.torqueLocal.y = std::clamp(yawErr * faceGain, -1.0, 1.0);
    ai.torqueLocal.z = 0.0;

    // Speed control
    double vAim = 0.0;
    if (dist > desiredDistKm) vAim = std::min(maxSpeedKmS, 0.000004 * (dist - desiredDistKm));
    if (dist < desiredDistKm * 0.60) vAim = -0.10;

    const math::Vec3d desiredVel = targetVelKmS + toN * vAim;
    const math::Vec3d dv = desiredVel - selfShip.velocityKmS();
    const math::Vec3d accelWorldDir = (dv.lengthSq() > 1e-9) ? dv.normalized() : math::Vec3d{0,0,0};
    ai.thrustLocal = selfShip.orientation().conjugate().rotate(accelWorldDir);
  };

  // Contacts AI + combat
  for (auto& c : contacts) {
    if (!c.alive) continue;

    // Keep contact dampers in the same local reference frame as the player.
    c.ship.setDampingFrameVelocityKmS(localFrameVelKmS);

    // cooldowns
    c.fireCooldown = std::max(0.0, c.fireCooldown - dtSim);

    // ---- PIRATES ----
    if (c.role == ContactRole::Pirate) {
      sim::ShipInput ai{};
      const bool hostile = c.hostileToPlayer || c.missionTarget;

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
      const math::Vec3d toP = ship.positionKm() - c.ship.positionKm();
      const double distP = toP.length();
      const math::Vec3d toPN = (distP > 1e-6) ? (toP / distP) : math::Vec3d{0,0,1};
      math::Vec3d right = math::cross(toPN, math::Vec3d{0,1,0});
      if (right.lengthSq() < 1e-9) right = math::cross(toPN, math::Vec3d{1,0,0});
      right = right.normalized();
      const math::Vec3d up = math::cross(right, toPN).normalized();

      const double spreadOn = (c.groupId != 0) ? 1.0 : 0.0;
      const double s = spreadOn * (((int)(c.id % 7) - 3) * 2200.0);
      const double u = spreadOn * (((int)(c.id % 5) - 2) * 1500.0);
      const math::Vec3d aimPos = ship.positionKm() + right * s + up * u;

      const double standOffKm = 35000.0 + ((c.leaderId != 0) ? 5000.0 : 0.0);
      chaseTarget(c.ship, ai, aimPos, ship.velocityKmS(), standOffKm, 0.22, 1.8);
      c.ship.step(dtSim, ai);

      // Fire if aligned
      const math::Vec3d to = ship.positionKm() - c.ship.positionKm();
      const double dist = to.length();
      if (hostile && c.fireCooldown <= 0.0 && dist < 90000.0) {
        const math::Vec3d toN = (dist > 1e-6) ? (to / dist) : math::Vec3d{0,0,1};
        const double aim = math::dot(c.ship.forward().normalized(), toN);
        if (aim > 0.992) {
          c.fireCooldown = 0.35;
          const double dmg = 11.0;
          applyDamage(dmg, playerShield, playerHull);

          const math::Vec3d aKm = c.ship.positionKm();
          const math::Vec3d bKm = ship.positionKm();
          beams.push_back({toRenderU(aKm), toRenderU(bKm), 0.95f, 0.45f, 0.10f, 0.08});
        }
      }
      continue;
    }

    // ---- TRADERS ----
    if (c.role == ContactRole::Trader) {
      sim::ShipInput ai{};
      ai.dampers = true;

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
            // Deliver current cargo (if any).
            auto& stEcon = universe.stationEconomy(st, timeDays);
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
            } else {
              c.cargoValueCr = 0.0;
            }

            // Cooldown to prevent multiple trades in one arrival.
            c.tradeCooldownUntilDays = timeDays + (20.0 / 86400.0);
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

	      sim::ShipInput ai{};
	      if (hostile) {
	        const double standoff = 42000.0 + ((c.leaderId != 0) ? 7000.0 : 0.0) + (double)((c.id % 3) * 1500);
	        chaseTarget(c.ship, ai, ship.positionKm(), ship.velocityKmS(), standoff, 0.26, 2.0);
	      } else if (pirateIdx) {
	        const auto& p = contacts[*pirateIdx];
	        const double standoff = 42000.0 + ((c.leaderId != 0) ? 6500.0 : 0.0) + (double)((c.id % 3) * 1400);
	        chaseTarget(c.ship, ai, p.ship.positionKm(), p.ship.velocityKmS(), standoff, 0.26, 2.0);
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

	      // Fire if aligned (at player OR pirate target).
	      const math::Vec3d targetPos = canFireAtPlayer ? ship.positionKm()
	                                                    : (pirateIdx ? contacts[*pirateIdx].ship.positionKm() : math::Vec3d{0,0,0});
	      const bool hasTarget = canFireAtPlayer || (bool)pirateIdx;
	      if (hasTarget) {
	        const math::Vec3d to = targetPos - c.ship.positionKm();
	        const double dist = to.length();
	        if (c.fireCooldown <= 0.0 && dist < 95000.0) {
	          const math::Vec3d toN = (dist > 1e-6) ? (to / dist) : math::Vec3d{0,0,1};
	          const double aim = math::dot(c.ship.forward().normalized(), toN);
	          if (aim > 0.993) {
	            c.fireCooldown = 0.28;
	            const double dmg = 12.0;

	            if (canFireAtPlayer) {
	              applyDamage(dmg, playerShield, playerHull);
	              const math::Vec3d aKm = c.ship.positionKm();
	              const math::Vec3d bKm = ship.positionKm();
	              beams.push_back({toRenderU(aKm), toRenderU(bKm), 0.35f, 0.75f, 1.0f, 0.07});
	            } else if (pirateIdx) {
	              auto& p = contacts[*pirateIdx];
	              applyDamage(dmg, p.shield, p.hull);
	              const math::Vec3d aKm = c.ship.positionKm();
	              const math::Vec3d bKm = p.ship.positionKm();
	              beams.push_back({toRenderU(aKm), toRenderU(bKm), 0.35f, 0.75f, 1.0f, 0.07});
	              if (p.hull <= 0.0) {
	                p.alive = false;
	                credits += 180.0;
	                toast(toasts, "Security destroyed a pirate (+180).", 2.0);
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

          if (insideStationHullExceptSlot(st, relLocal)) {
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
          const double repPenalty = -std::clamp(3.0 + bribeOffer.illegalValueCr / 2800.0, 3.0, 12.0);
          addRep(bribeOffer.factionId, repPenalty);

          // Running turns the fine into a bounty. (You keep the cargo, but you're now WANTED.)
          addBounty(bribeOffer.factionId, bribeOffer.fineCr);

          if (escalateLocal) {
            policeHeat = std::clamp(policeHeat + 1.2 + std::min(2.0, bribeOffer.fineCr / 2500.0), 0.0, 6.0);
            policeAlertUntilDays = std::max(policeAlertUntilDays, timeDays + (140.0 / 86400.0));
            nextPoliceSpawnDays = std::min(nextPoliceSpawnDays, timeDays + (4.0 / 86400.0));
          }

          heat = std::min(100.0, heat + 5.0);

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

                double illegalValueCr = 0.0;
                std::array<double, econ::kCommodityCount> scannedIllegal{};
                scannedIllegal.fill(0.0);

                std::string detail;
                int shown = 0;

                for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
                  const auto cid = (econ::CommodityId)i;
                  if (!isIllegalCommodity(jurisdiction, cid)) continue;

                  const double units = cargo[i];
                  if (units <= 1e-6) continue;

                  double mid = econ::commodityDef(cid).basePrice;
                  if (priceRef) {
                    auto& stEcon = universe.stationEconomy(*priceRef, timeDays);
                    const auto q = econ::quote(stEcon, priceRef->economyModel, cid, 0.10);
                    mid = q.mid;
                  }

                  illegalValueCr += units * mid;
                  scannedIllegal[i] = units;

                  if (shown < 2) {
                    if (!detail.empty()) detail += ", ";
                    detail += econ::commodityDef(cid).name;
                    detail += " x" + std::to_string((int)std::round(units));
                    ++shown;
                  }
                }
                if ((int)shown >= 2) detail += " ...";

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
                  const double fineCr = 200.0 + illegalValueCr * 0.65;

                  // Chance that the scanning police contact is corrupt and offers a bribe.
                  bool offerBribe = false;
                  if (cargoScanSourceKind == CargoScanSourceKind::Police) {
                    double chance = 0.33;

                    // Less bribable when your heat is already high.
                    chance *= std::clamp(1.10 - (heat / 110.0), 0.25, 1.10);

                    // Reputation influences corruption willingness.
                    const double rep = getRep(jurisdiction);
                    if (rep < -30.0) chance *= 0.35;
                    if (rep > 55.0) chance *= 0.75;

                    // Huge hauls are riskier to "look the other way" on.
                    chance *= std::clamp(1.05 - (illegalValueCr / 20000.0), 0.35, 1.05);

                    offerBribe = rng.nextUnit() < chance;
                  }

                  if (offerBribe) {
                    bribeOffer.active = true;
                    bribeOffer.factionId = jurisdiction;
                    bribeOffer.scannerContactId = cargoScanContactId;
                    bribeOffer.scannerRangeKm = cargoScanRangeKm;
                    bribeOffer.illegalValueCr = illegalValueCr;
                    bribeOffer.fineCr = fineCr;
                    bribeOffer.amountCr = std::max(200.0, 250.0 + illegalValueCr * 0.55);
                    bribeOffer.amountCr = std::round(bribeOffer.amountCr / 10.0) * 10.0; // nicer UI numbers
                    bribeOffer.startDays = timeDays;
                    bribeOffer.untilDays = timeDays + (14.0 / 86400.0);
                    bribeOffer.sourceName = cargoScanSourceName.empty() ? "Security" : cargoScanSourceName;
                    bribeOffer.detail = detail;
                    bribeOffer.scannedIllegal = scannedIllegal;

                    // Keep some local pressure while you're deciding.
                    policeAlertUntilDays = std::max(policeAlertUntilDays, timeDays + (90.0 / 86400.0));
                    nextPoliceSpawnDays = std::min(nextPoliceSpawnDays, timeDays + (6.0 / 86400.0));

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
            const double baseRatePerSec = contraband ? 0.06 : 0.012;

            // Smuggling compartments reduce scan frequency significantly.
            const double holdMult = (smuggleHoldMk <= 0) ? 1.0
                              : (smuggleHoldMk == 1 ? 0.72
                              : (smuggleHoldMk == 2 ? 0.50
                              : 0.35));
            const double ratePerSec = baseRatePerSec * holdMult;

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
                const double baseDur = contraband ? 7.0 : 4.0;
                const double durMult = (smuggleHoldMk <= 0) ? 1.0
                                    : (smuggleHoldMk == 1 ? 1.15
                                    : (smuggleHoldMk == 2 ? 1.35
                                    : 1.60));
                cargoScanDurationSec = baseDur * durMult;

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

                  const double baseDur = 6.0;
                  const double durMult = (smuggleHoldMk <= 0) ? 1.0
                                      : (smuggleHoldMk == 1 ? 1.15
                                      : (smuggleHoldMk == 2 ? 1.35
                                      : 1.60));
                  cargoScanDurationSec = baseDur * durMult;
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
	      const auto& s = signals[scanLockedTarget.index];
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
	              completeScan(std::string("Asteroid prospected: ") + cname + " (+data " + std::to_string((int)value) + " cr).", 2.5);
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
	              toast(toasts, "Arrived at Resource Site: asteroid fragments detected.", 3.0);
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
	              toast(toasts, "Distress beacon acquired. Approach with caution.", 3.0);
	              // Some supplies may be drifting...
	              spawnCargoPod(econ::CommodityId::Food, rng.range(3.0, 12.0), s.posKm + randUnit() * rng.range(1500.0, 7000.0), {0,0,0}, 0.35);
	              // ...and sometimes it's an ambush.
	              if (rng.nextUnit() < 0.55) {
	                const int pirates = 1 + rng.range(0, 2);
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
	                  p.ship.setMaxLinearAccelKmS2(0.06);
	                  p.ship.setMaxAngularAccelRadS2(0.9);
	                  p.hullMax = npcHullMax;
	                  p.hull = npcHullMax;
	                  p.shieldMax = npcShieldMax;
	                  p.shield = npcShieldMax;
	                  p.shieldRegenPerSec = npcShieldRegenPerSec;
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
        double heatIn = 0.0;
        if (!docked) {
          if (input.boost) heatIn += 18.0;
          if (supercruiseState == SupercruiseState::Active) heatIn += 6.0;
          if (fsdState == FsdState::Charging) heatIn += 12.0;
          if (fsdState == FsdState::Jumping) heatIn += 16.0;
        }

        const double baseCool = docked ? 20.0 : 10.0;
        const double coolRate = baseCool * (playerHeatCoolRate / 10.0);
        heat += (heatIn - coolRate) * dtReal;
        heat = std::clamp(heat, 0.0, 120.0);

        // Very simple overheat consequence for now: take slow hull damage.
        if (heat > 100.0) {
          const double over = heat - 100.0;
          playerHull = std::max(0.0, playerHull - over * 0.0008 * playerHullMax * dtReal);
        }
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

      // Projectiles update (simple ballistic slugs)
      for (auto& p : projectiles) {
        if (p.ttlSim <= 0.0) continue;

        const math::Vec3d a = p.posKm;
        const math::Vec3d b = p.posKm + p.velKmS * dtSim;

        p.prevKm = a;
        p.posKm = b;
        p.ttlSim -= dtSim;

        // Collisions
        if (p.fromPlayer) {
          for (int i = 0; i < (int)contacts.size(); ++i) {
            const auto& c = contacts[(std::size_t)i];
            if (!c.alive) continue;
            if (p.shooterId != 0 && c.id == p.shooterId) continue;

            const double hitRadiusKm = 900.0 + p.radiusKm;
            if (segmentHitsSphere(a, b, c.ship.positionKm(), hitRadiusKm)) {
              playerDamageContact(i, p.dmg);
              p.ttlSim = 0.0;
              break;
            }
          }
        } else {
          if (!docked) {
            const double playerHitRadiusKm = 900.0 + p.radiusKm;
            if (segmentHitsSphere(a, b, ship.positionKm(), playerHitRadiusKm)) {
              applyDamage(p.dmg, playerShield, playerHull);

              if (vfxParticlesEnabled && vfxImpactsEnabled) {
                // Approximate impact normal by projectile travel direction.
                math::Vec3d n = -p.velKmS;
                if (n.lengthSq() < 1e-12) n = ship.forward();
                n = n.normalized();

                const double energy = std::clamp((p.dmg / 16.0) * (double)vfxParticleIntensity, 0.12, 2.2);
                particles.spawnSparks(toRenderU(ship.positionKm()), n, toRenderU(ship.velocityKmS()), energy);
              }

              p.ttlSim = 0.0;
            }
          }
        }
      }

      projectiles.erase(
        std::remove_if(projectiles.begin(), projectiles.end(), [](const Projectile& p) { return p.ttlSim <= 0.0; }),
        projectiles.end());


      weaponPrimaryCooldown = std::max(0.0, weaponPrimaryCooldown - dtSim);
      weaponSecondaryCooldown = std::max(0.0, weaponSecondaryCooldown - dtSim);

      // Shield regen (slow)
      if (!paused && playerHull > 0.0 && playerShield < playerShieldMax) {
        playerShield = std::min(playerShieldMax, playerShield + playerShieldRegenPerSimMin * (dtSim / 60.0));
      }

      // NPC shield regen (simple; gives fights a bit of pacing without heavy AI)
      for (auto& c : contacts) {
        if (!c.alive) continue;
        if (c.hull <= 0.0) continue;
        if (c.shieldRegenPerSec <= 0.0) continue;
        if (c.shieldMax <= 0.0) continue;
        if (c.shield >= c.shieldMax) continue;
        c.shield = std::min(c.shieldMax, c.shield + c.shieldRegenPerSec * dtSim);
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
      interdictionState = InterdictionState::None;
      interdictionWarningRemainingSec = 0.0;
      interdictionActiveRemainingSec = 0.0;
      interdictionEscapeMeter = 0.0;
      interdictionSubmitRequested = false;
      interdictionPirateId = 0;
      interdictionPirateName.clear();
      scanning = false;
      scanProgressSec = 0.0;
      fsdState = FsdState::Idle;
      fsdTargetSystem = 0;
      fsdChargeRemainingSec = 0.0;
      fsdTravelRemainingSec = 0.0;
      fsdTravelTotalSec = 0.0;
      navAutoRun = false;
      beams.clear();
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
          clearanceGranted = it->second.granted && (timeDays <= it->second.expiresDays);
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

    // DockSpace is only available in Dear ImGui's docking branch; keep UI
    // functional without it (windows will simply be floating).

    // Main menu bar (quick window toggles + HUD settings)
    if (ImGui::BeginMainMenuBar()) {
      if (ImGui::BeginMenu("Windows")) {
        ImGui::MenuItem("Galaxy (TAB)", nullptr, &showGalaxy);
        ImGui::MenuItem("Ship/Status (F1)", nullptr, &showShip);
        ImGui::MenuItem("Contacts (F3)", nullptr, &showContacts);
        ImGui::MenuItem("Market (F2)", nullptr, &showEconomy);
        ImGui::MenuItem("Missions (F4)", nullptr, &showMissions);
        ImGui::MenuItem("Scanner (F6)", nullptr, &showScanner);
        ImGui::MenuItem("Trade Finder (F7)", nullptr, &showTrade);
        ImGui::Separator();
        ImGui::MenuItem("Guide (F8)", nullptr, &showGuide);
        ImGui::MenuItem("Sprite Lab (F10)", nullptr, &showSprites);
        ImGui::MenuItem("World Visuals", nullptr, &showWorldVisuals);
        ImGui::MenuItem("Hangar / Livery", nullptr, &showHangar);
        ImGui::MenuItem("VFX Lab (F11)", nullptr, &showVfx);
        ImGui::MenuItem("PostFX (F12)", nullptr, &showPostFx);
        ImGui::EndMenu();
      }

      if (ImGui::BeginMenu("HUD")) {
        ImGui::MenuItem("Radar (R)", nullptr, &showRadarHud);
        if (showRadarHud) {
          const double minKm = 25000.0;
          const double maxKm = 1200000.0;
          ImGui::SetNextItemWidth(220.0f);
          ImGui::DragScalar("Range (km)", ImGuiDataType_Double, &radarRangeKm, 5000.0, &minKm, &maxKm, "%.0f");
          ImGui::SetNextItemWidth(220.0f);
          ImGui::SliderInt("Max blips", &radarMaxBlips, 16, 160);
        }

        ImGui::MenuItem("Objective HUD overlay", nullptr, &objectiveHudEnabled);
        ImGui::MenuItem("Threat overlay HUD", nullptr, &hudThreatOverlayEnabled);
        ImGui::MenuItem("Jump overlay HUD", nullptr, &hudJumpOverlay);

        ImGui::MenuItem("Offscreen target indicator", nullptr, &hudOffscreenTargetIndicator);

        ImGui::SeparatorText("Combat");
        ImGui::Checkbox("Combat HUD", &hudCombatHud);
        if (hudCombatHud) {
          ImGui::Checkbox("Procedural reticle", &hudUseProceduralReticle);
          ImGui::SameLine();
          ImGui::Checkbox("Weapon rings", &hudShowWeaponRings);
          ImGui::SameLine();
          ImGui::Checkbox("Heat ring", &hudShowHeatRing);

          ImGui::SetNextItemWidth(220.0f);
          ImGui::SliderFloat("Reticle size (px)", &hudReticleSizePx, 24.0f, 96.0f, "%.0f");
          ImGui::SetNextItemWidth(220.0f);
          ImGui::SliderFloat("Reticle alpha", &hudReticleAlpha, 0.10f, 1.0f, "%.2f");

          ImGui::Separator();
          ImGui::Checkbox("Lead indicator (projectiles)", &hudShowLeadIndicator);
          ImGui::SameLine();
          ImGui::Checkbox("Use last fired weapon", &hudLeadUseLastFiredWeapon);
          ImGui::SetNextItemWidth(220.0f);
          ImGui::SliderFloat("Lead size (px)", &hudLeadSizePx, 10.0f, 64.0f, "%.0f");
          const double leadMinT = 0.5;
          const double leadMaxT = 60.0;
          ImGui::SetNextItemWidth(220.0f);
          ImGui::DragScalar("Lead max time (s)", ImGuiDataType_Double, &hudLeadMaxTimeSec, 0.5, &leadMinT, &leadMaxT, "%.1f");

          ImGui::Separator();
          ImGui::Checkbox("Flight path marker", &hudShowFlightPathMarker);
          ImGui::SameLine();
          ImGui::Checkbox("Local frame", &hudFlightMarkerUseLocalFrame);
          ImGui::SameLine();
          ImGui::Checkbox("Clamp to edge", &hudFlightMarkerClampToEdge);
          ImGui::SetNextItemWidth(220.0f);
          ImGui::SliderFloat("Marker size (px)", &hudFlightMarkerSizePx, 10.0f, 64.0f, "%.0f");

          ImGui::TextDisabled("Reticle rings show cooldown + heat; lead marker predicts projectile intercept.");
        }

        ImGui::Separator();
        ImGui::MenuItem("Tactical overlay (`)", nullptr, &showTacticalOverlay);
        if (showTacticalOverlay) {
          ImGui::Checkbox("Labels", &tacticalShowLabels);
          const double minKm = 25000.0;
          const double maxKm = 2000000.0;
          ImGui::SetNextItemWidth(220.0f);
          ImGui::DragScalar("Range (km)##tac", ImGuiDataType_Double, &tacticalRangeKm, 5000.0, &minKm, &maxKm, "%.0f");
          ImGui::SetNextItemWidth(220.0f);
          ImGui::SliderInt("Max markers##tac", &tacticalMaxMarkers, 16, 200);

          ImGui::SeparatorText("Types");
          ImGui::Checkbox("Stations", &tacticalShowStations);
          ImGui::SameLine();
          ImGui::Checkbox("Planets", &tacticalShowPlanets);
          ImGui::SameLine();
          ImGui::Checkbox("Starships", &tacticalShowContacts);
          ImGui::Checkbox("Cargo", &tacticalShowCargo);
          ImGui::SameLine();
          ImGui::Checkbox("Asteroids", &tacticalShowAsteroids);
          ImGui::SameLine();
          ImGui::Checkbox("Signals", &tacticalShowSignals);

          ImGui::TextDisabled("Hover markers for tooltips; MMB targets.");
        }

        ImGui::SeparatorText("Layout");
        ImGui::Checkbox("Edit HUD layout (Ctrl+H)", &hudLayoutEditMode);
        ImGui::SameLine();
        ImGui::Checkbox("Layout window", &showHudLayoutWindow);
        ImGui::SetNextItemWidth(220.0f);
        ImGui::SliderFloat("Safe margin (px)", &hudLayout.safeMarginPx, 0.0f, 64.0f, "%.0f");

        if (ImGui::MenuItem("Save HUD layout (Ctrl+S)")) {
          hudLayout.widget(ui::HudWidgetId::Radar).enabled = showRadarHud;
          hudLayout.widget(ui::HudWidgetId::Objective).enabled = objectiveHudEnabled;
          hudLayout.widget(ui::HudWidgetId::Threat).enabled = hudThreatOverlayEnabled;
          hudLayout.widget(ui::HudWidgetId::Jump).enabled = hudJumpOverlay;
          ui::saveToFile(hudLayout, hudLayoutPath);
        }
        if (ImGui::MenuItem("Load HUD layout (Ctrl+L)")) {
          ui::HudLayout loaded = ui::makeDefaultHudLayout();
          if (ui::loadFromFile(hudLayoutPath, loaded)) {
            hudLayout = loaded;
            showRadarHud = hudLayout.widget(ui::HudWidgetId::Radar).enabled;
            objectiveHudEnabled = hudLayout.widget(ui::HudWidgetId::Objective).enabled;
            hudThreatOverlayEnabled = hudLayout.widget(ui::HudWidgetId::Threat).enabled;
            hudJumpOverlay = hudLayout.widget(ui::HudWidgetId::Jump).enabled;
          }
        }
        if (ImGui::MenuItem("Reset HUD layout (Ctrl+R)")) {
          hudLayout = ui::makeDefaultHudLayout();
          showRadarHud = hudLayout.widget(ui::HudWidgetId::Radar).enabled;
          objectiveHudEnabled = hudLayout.widget(ui::HudWidgetId::Objective).enabled;
          hudThreatOverlayEnabled = hudLayout.widget(ui::HudWidgetId::Threat).enabled;
          hudJumpOverlay = hudLayout.widget(ui::HudWidgetId::Jump).enabled;
        }

        ImGui::EndMenu();
      }

      // Right-side status text.
      {
        const std::string sysName = currentSystem ? currentSystem->stub.name : "(no system)";
        const std::string status = " | " + sysName + " | " + (docked ? "DOCKED" : "IN FLIGHT");
        const float wText = ImGui::CalcTextSize(status.c_str()).x;
        ImGui::SameLine(ImGui::GetWindowWidth() - wText - 14.0f);
        ImGui::TextDisabled("%s", status.c_str());
      }

      ImGui::EndMainMenuBar();
    }

    // HUD Layout editor window
    if (showHudLayoutWindow) {
      ImGui::SetNextWindowSize(ImVec2(520.0f, 450.0f), ImGuiCond_FirstUseEver);
      ImGui::Begin("HUD Layout", &showHudLayoutWindow);

      ImGui::TextDisabled(
          "Reposition in-flight HUD overlays. Drag the HUD panels while Edit mode is enabled.");
      ImGui::TextDisabled("Shortcuts: Ctrl+H edit | Ctrl+S save | Ctrl+L load | Ctrl+R reset");

      ImGui::Checkbox("Edit mode (drag panels)", &hudLayoutEditMode);
      ImGui::SameLine();
      ImGui::Checkbox("Auto-save on exit", &hudLayoutAutoSaveOnExit);

      ImGui::SetNextItemWidth(240.0f);
      ImGui::SliderFloat("Safe margin guide (px)", &hudLayout.safeMarginPx, 0.0f, 64.0f, "%.0f");

      if (ImGui::Button("Save (Ctrl+S)")) {
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
      if (ImGui::Button("Load (Ctrl+L)")) {
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
      if (ImGui::Button("Reset (Ctrl+R)")) {
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

    // HUD overlay: combat symbology + target marker + toasts
    {
      ImDrawList* draw = ImGui::GetForegroundDrawList();
      const ImVec2 center((float)w * 0.5f, (float)h * 0.5f);

      // HUD layout edit guides (safe area + pivot markers)
      if (hudLayoutEditMode) {
        ImGuiViewport* vp = ImGui::GetMainViewport();
        const float m = std::max(0.0f, hudLayout.safeMarginPx);
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
	        tgtLabel = std::string(signalTypeName(s.type)) + " Signal";
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

              // Solve |r + v t| = s t for t>0 (classic lead/intercept quadratic).
              const math::Vec3d r = (*tgtKm - ship.positionKm());
              const math::Vec3d v = (*tgtVelKmS - ship.velocityKmS());
              const double a = v.lengthSq() - s * s;
              const double b = 2.0 * math::dot(r, v);
              const double c = r.lengthSq();

              double t = -1.0;
              const double eps = 1e-9;
              if (std::abs(a) < eps) {
                if (std::abs(b) > eps) {
                  t = -c / b;
                }
              } else {
                const double disc = b * b - 4.0 * a * c;
                if (disc >= 0.0) {
                  const double sd = std::sqrt(disc);
                  const double t1 = (-b - sd) / (2.0 * a);
                  const double t2 = (-b + sd) / (2.0 * a);
                  const double tMin = std::min(t1, t2);
                  const double tMax = std::max(t1, t2);
                  if (tMin > 1.0e-5) t = tMin;
                  else if (tMax > 1.0e-5) t = tMax;
                }
              }

              const double ttl = wLead.rangeKm / std::max(1e-9, s);
              if (t > 0.0 && t <= ttl && t <= hudLeadMaxTimeSec) {
                const math::Vec3d leadKm = *tgtKm + (*tgtVelKmS) * t;
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
          clearanceGranted = it->second.granted && (timeDays <= it->second.expiresDays);
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
            addMarker(TargetKind::Signal,
                      i,
                      s.posKm,
                      render::SpriteKind::Signal,
                      core::hashCombine(core::hashCombine(core::fnv1a64("signal"), s.id), (core::u64)(int)s.type),
                      std::string(signalTypeName(s.type)) + " Signal");
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
      for (const auto& t : toasts) {
        draw->AddText({18.0f, y}, IM_COL32(240,240,240,220), t.text.c_str());
        y += 18.0f;
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

  // Bounty missions only care about the target system (station is not applicable).
  if (m.type == sim::MissionType::BountyScan || m.type == sim::MissionType::BountyKill) {
    return {m.toSystem, 0};
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
  } else if (m.type == sim::MissionType::BountyKill) {
    out = "Bounty Hunt: eliminate target in " + uiSystemNameById(m.toSystem);
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

  ImGui::TextDisabled("F11 toggles this panel. Background stars + particles are rendered as point sprites.");

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

  ImGui::TextDisabled("F12 toggles this panel. Scene renders to an HDR target and is post-processed before UI.");

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
  ImGui::Text("Cargo: %.0f / %.0f kg", cargoMassKg(cargo), cargoCapacityKg);
  ImGui::Text("Passengers: %d seats", std::max(0, passengerSeats));
	  ImGui::Text("Cargo scoop (O): %s | Floating pods: %d", cargoScoopDeployed ? "DEPLOYED" : "RETRACTED",
	              (int)floatingCargo.size());

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
  ImGui::Checkbox("Mouse steer (M)", &mouseSteer);
  ImGui::Checkbox("Invert mouse Y", &mouseInvertY);
  ImGui::SliderFloat("Mouse sensitivity", &mouseSensitivity, 0.0006f, 0.0080f, "%.4f");
  ImGui::TextDisabled("Mouse steer captures the cursor (relative mouse mode).");

  if (scanning) {
    ImGui::Separator();
    ImGui::TextColored(ImVec4(0.9f, 0.95f, 1.0f, 1.0f), "%s", scanLabel.empty() ? "Scanning..." : scanLabel.c_str());
    const float frac = (float)std::clamp(scanProgressSec / std::max(0.01, scanDurationSec), 0.0, 1.0);
    ImGui::ProgressBar(frac, ImVec2(-1, 0));
    ImGui::TextDisabled("K cancels scan.");
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
  } else if (target.kind == TargetKind::Star) {
    ImGui::Text("Target: Star (%s)", starClassName(currentSystem->star.cls));
  } else {
    ImGui::Text("Target: (none)  [T station / B planet / N contact / U star]");
  }

  ImGui::Separator();
  ImGui::TextDisabled("Controls:");
  ImGui::BulletText("WASD / Space/Ctrl: translate | Arrows: pitch/yaw | Q/E: roll");
  ImGui::BulletText("Shift: boost | X: brake | LMB: %s | RMB: %s", weaponDef(weaponPrimary).name, weaponDef(weaponSecondary).name);
  ImGui::BulletText("H: supercruise (charge/engage). While active: H drops (~7s safe). During interdiction: H submits");
  ImGui::BulletText("J engage FSD jump (uses plotted route if present)");
  ImGui::BulletText("P: autopilot to station (staging + corridor guidance)");
  ImGui::BulletText("T/B/N/U cycle targets, Y clear target");
  ImGui::BulletText("K scan target (missions + exploration scans)");
  ImGui::BulletText("L request docking clearance | G dock / undock");
  ImGui::BulletText("TAB Galaxy map, F1 Ship, F2 Market, F3 Contacts, F4 Missions, F6 Scanner, F7 Trade, F8 Guide");
  ImGui::BulletText("F5 quicksave, F9 quickload");
  ImGui::BulletText("F10 Sprite Lab, F11 VFX Lab");

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

  ImGui::Text("Quick start checklist");
  ImGui::BulletText("Dock: target station (T), request clearance (L), fly corridor (P), dock/undock (G)");
  ImGui::BulletText("Missions: open Mission Board (F4). Accept/Plot auto-plot; track a mission for the Objective HUD");
  ImGui::BulletText("Jump: J (uses plotted route if present)");
  ImGui::BulletText("Supercruise: H engage, H drop near ~7s ETA. Interdiction: align to escape vector; H submits");
  ImGui::BulletText("Trade helper: F7 (best local routes)");
  ImGui::BulletText("Scan: K (stars/planets/contacts). Sell exploration data at stations");
  ImGui::BulletText("Sprite Lab: F10 (procedural UI icons preview + cache)");
  ImGui::BulletText("VFX Lab: F11 (starfield + particles)");

  ImGui::Separator();
  ImGui::Text("Quick actions");
  if (ImGui::Button("Open Galaxy (TAB)")) showGalaxy = true;
  ImGui::SameLine();
  if (ImGui::Button("Open Market (F2)")) showEconomy = true;
  ImGui::SameLine();
  if (ImGui::Button("Open Missions (F4)")) showMissions = true;
  ImGui::SameLine();
  if (ImGui::Button("Open Trade (F7)")) showTrade = true;

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

  if (supercruiseState == SupercruiseState::Active && interdictionState != InterdictionState::None) {
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

  ImGui::TextDisabled("F10 toggles this panel. Procedural icons are generated deterministically from (kind, seed).");
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
	      scanDurationSec = 4.0;
	      scanRangeKm = 200000.0;
	      scanLabel = std::string("Signal scan: ") + signalTypeName(s.type);
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
          ImGui::Text("%s Signal", signalTypeName(s.type));
          ImGui::Image((ImTextureID)(intptr_t)big.handle(), ImVec2(96, 96));
          ImGui::EndTooltip();
        }
        ImGui::SameLine();
      }

	      ImGui::Text("%s%s%s | %.0f km%s", signalTypeName(s.type), s.resolved ? " (resolved)" : "", mtag, distKm,
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
          interdictionState = InterdictionState::None;
          interdictionWarningRemainingSec = 0.0;
          interdictionActiveRemainingSec = 0.0;
          interdictionEscapeMeter = 0.0;
          interdictionSubmitRequested = false;
          autopilot = false;
          autopilotPhase = 0;
          scanning = false;
          scanProgressSec = 0.0;
          navAutoRun = false;

          toast(toasts,
                pirateNearby ? "Supercruise charging... (pirate signals nearby)" : "Supercruise charging...",
                1.8);
        } else if (supercruiseState == SupercruiseState::Charging) {
          supercruiseState = SupercruiseState::Idle;
          supercruiseChargeRemainingSec = 0.0;
          interdictionState = InterdictionState::None;
          interdictionSubmitRequested = false;
          toast(toasts, "Supercruise canceled.", 1.4);
        } else if (supercruiseState == SupercruiseState::Active) {
          if (interdictionState != InterdictionState::None) {
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
            if (ImGui::MenuItem("Scan (K)")) {
              uiStartScan(ctxKind, ctxIdx);
            }

            // Supercruise only makes sense for station/planet/signal.
            const bool canSc = (ctxKind == TargetKind::Station || ctxKind == TargetKind::Planet || ctxKind == TargetKind::Signal);
            if (ImGui::MenuItem("Engage supercruise (H)", nullptr, false, canSc)) {
              target.kind = ctxKind;
              target.index = ctxIdx;
              uiToggleSupercruise();
            }

            if (ctxKind == TargetKind::Station) {
              if (ImGui::MenuItem("Autopilot to station (P)")) {
                target.kind = ctxKind;
                target.index = ctxIdx;
                autopilot = true;
                autopilotPhase = 0;
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
          ImGui::TextDisabled("Yield: %s", econ::commodityDef(a.yield).name.c_str());
          ImGui::TextDisabled("Remaining: %.0f", a.remainingUnits);
          ImGui::TextDisabled("Distance: %.0f km", (a.posKm - ship.positionKm()).length());
        } else if (target.kind == TargetKind::Cargo && target.index < floatingCargo.size()) {
          const auto& pod = floatingCargo[target.index];
          const core::u64 cSeed = core::hashCombine(core::hashCombine(core::fnv1a64("cargo"), pod.id), (core::u64)pod.commodity);
          showTargetHeader("Cargo Pod", render::SpriteKind::Cargo, cSeed, ImVec4(1,1,1,1));
          ImGui::TextDisabled("Commodity: %s", econ::commodityDef(pod.commodity).name.c_str());
          ImGui::TextDisabled("Units: %.0f", pod.units);
          ImGui::TextDisabled("Distance: %.0f km", (pod.posKm - ship.positionKm()).length());
        } else {
          ImGui::TextDisabled("Invalid selection.");
        }

        ImGui::Separator();

        // Actions
        ImGui::BeginDisabled(target.kind == TargetKind::None);
        if (ImGui::Button("Scan target (K)")) uiStartScan(target.kind, target.index);
        ImGui::EndDisabled();

        const bool canScTarget = (target.kind == TargetKind::Station || target.kind == TargetKind::Planet || target.kind == TargetKind::Signal);
        ImGui::BeginDisabled(!canScTarget);
        if (ImGui::Button("Supercruise (H)")) uiToggleSupercruise();
        ImGui::EndDisabled();

        if (target.kind == TargetKind::Station) {
          ImGui::BeginDisabled(docked);
          if (ImGui::Button("Autopilot to station (P)")) {
            autopilot = true;
            autopilotPhase = 0;
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
            }
            ImGui::SameLine();
            if (ImGui::SmallButton("Refresh")) {
              tradeIdeasDayStamp = -1;
            }

            const int dayStamp = (int)std::floor(timeDays);
            if (tradeFromStationId != fromStationId || tradeIdeasDayStamp != dayStamp) {
              tradeIdeas.clear();
              tradeFromStationId = fromStationId;
              tradeIdeasDayStamp = dayStamp;

              auto& fromEcon = universe.stationEconomy(*fromSt, timeDays);

              const auto nearby = universe.queryNearby(currentSystem->stub.posLy, tradeSearchRadiusLy);
              for (const auto& stub : nearby) {
                const auto& sys = universe.getSystem(stub.id, &stub);
                for (const auto& toSt : sys.stations) {
                  if (sys.stub.id == currentSystem->stub.id && toSt.id == fromStationId) continue;

                  auto& toEcon = universe.stationEconomy(toSt, timeDays);
                  const auto routes = econ::bestRoutes(fromEcon, fromSt->economyModel, toEcon, toSt.economyModel, 0.10, 2);
                  for (const auto& r : routes) {
                    const double massKg = std::max(1e-6, econ::commodityDef(r.commodity).massKg);

                    // Approximate *net* profit per unit (rep-adjusted fees at both stations).
                    const double feeTo = effectiveFeeRate(toSt);
                    const double netBuy = r.buyPrice * (1.0 + feeEff);
                    const double netSell = r.sellPrice * (1.0 - feeTo);
                    const double netProfit = netSell - netBuy;
                    if (netProfit <= 0.0) continue;

                    TradeIdea t;
                    t.toSystem = sys.stub.id;
                    t.toStation = toSt.id;
                    t.commodity = r.commodity;
                    t.buyPrice = r.buyPrice;
                    t.sellPrice = r.sellPrice;
                    t.profitPerUnit = netProfit;
                    t.profitPerKg = netProfit / massKg;
                    t.distanceLy = (sys.stub.posLy - currentSystem->stub.posLy).length();

                    const double unitsFull = std::floor(cargoCapacityKg / massKg);
                    t.estTripProfit = unitsFull * netProfit;
                    tradeIdeas.push_back(t);
                  }
                }
              }

              std::sort(tradeIdeas.begin(), tradeIdeas.end(), [](const TradeIdea& a, const TradeIdea& b) {
                return a.estTripProfit > b.estTripProfit;
              });
              if (tradeIdeas.size() > 12) tradeIdeas.resize(12);
            }

            if (tradeIdeas.empty()) {
              ImGui::TextDisabled("No profitable routes found in the selected radius.");
              ImGui::TextDisabled("Tip: expand the radius or wait a day for prices to move.");
            } else {
              auto stationNameInSystem = [&](sim::SystemId sysId, sim::StationId stId) -> std::string {
                const auto& sys = universe.getSystem(sysId);
                for (const auto& st : sys.stations) {
                  if (st.id == stId) return st.name;
                }
                return "Station #" + std::to_string((std::uint64_t)stId);
              };

              const double jr = std::max(1e-6, fsdCurrentRangeLy());

              ImGui::TextDisabled("Top routes (estimated with a full hold):");
              if (ImGui::BeginTable("trade_table", 6, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_SizingStretchSame)) {
                ImGui::TableSetupColumn("To");
                ImGui::TableSetupColumn("Commodity");
                ImGui::TableSetupColumn("Buy/Sell");
                ImGui::TableSetupColumn("Profit / unit");
                ImGui::TableSetupColumn("Profit / trip");
                ImGui::TableSetupColumn("Actions");
                ImGui::TableHeadersRow();

                for (std::size_t i = 0; i < tradeIdeas.size(); ++i) {
                  const auto& t = tradeIdeas[i];
                  const auto& sys = universe.getSystem(t.toSystem);
                  const std::string toStName = stationNameInSystem(t.toSystem, t.toStation);
                  const std::string toLabel = sys.stub.name + " / " + toStName;

                  ImGui::TableNextRow();
                  ImGui::TableNextColumn();
                  ImGui::Text("%s", toLabel.c_str());
                  ImGui::TextDisabled("%.0f ly (~%d jumps)", t.distanceLy, (int)std::ceil(t.distanceLy / jr));

                  ImGui::TableNextColumn();
                  // CommodityDef::name is a const char* (not std::string)
                  ImGui::Text("%s", econ::commodityDef(t.commodity).name);
                  ImGui::TextDisabled("%.0f kg/u", econ::commodityDef(t.commodity).massKg);

                  ImGui::TableNextColumn();
                  ImGui::Text("%.0f / %.0f", t.buyPrice, t.sellPrice);

                  ImGui::TableNextColumn();
                  ImGui::Text("+%.0f", t.profitPerUnit);
                  ImGui::TextDisabled("(+%.0f / kg)", t.profitPerKg);

                  ImGui::TableNextColumn();
                  ImGui::Text("+%.0f", t.estTripProfit);

                  ImGui::TableNextColumn();
                  ImGui::PushID((int)i);

                  if (ImGui::SmallButton("Plot")) {
                    if (plotRouteToSystem(t.toSystem)) {
                      pendingArrivalTargetStationId = t.toStation;
                      if (currentSystem && t.toSystem == currentSystem->stub.id) {
                        tryTargetStationById(t.toStation);
                      }
                      showGalaxy = true;
                    }
                  }

                  if (fromDocked) {
                    ImGui::SameLine();
                    if (ImGui::SmallButton("Fill")) {
                      auto& fromEcon = universe.stationEconomy(*fromSt, timeDays);
                      const double massKg = std::max(1e-6, econ::commodityDef(t.commodity).massKg);
                      const double freeKg = std::max(0.0, cargoCapacityKg - cargoMassKg(cargo));
                      const double maxUnitsHold = std::floor(freeKg / massKg);
                      const double maxUnitsCredits = std::floor(credits / std::max(1e-6, t.buyPrice * (1.0 + feeEff)));
                      const double maxUnits = std::max(0.0, std::min({maxUnitsHold, maxUnitsCredits, fromEcon.inventory[(int)t.commodity]}));

                      if (maxUnits <= 0.0) {
                        toast(toasts, "Can't buy: insufficient hold / credits / inventory.", 2.2);
                      } else {
                        auto tr = econ::buy(fromEcon, fromSt->economyModel, t.commodity, maxUnits, credits, 0.10, feeEff);
                        if (tr.ok) {
                          cargo[(int)t.commodity] += maxUnits;
                          toast(toasts, ("Bought " + std::to_string((int)maxUnits) + "u of " + econ::commodityDef(t.commodity).name).c_str(), 2.2);
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
              if (ImGui::Selectable(r.s->name.c_str(), isSel)) {
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

            if (ImGui::Button("Engage FSD (J)")) {
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

      ImGui::TextDisabled("TAB toggles this window, F1 Flight, F2 Market, F3 Contacts, F12 PostFX");
      ImGui::End();
    }

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

  // Cleanup
  spriteCache.clear();
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplSDL2_Shutdown();
  ImGui::DestroyContext();
  }

  SDL_GL_DeleteContext(glContext);
  SDL_DestroyWindow(window);
  SDL_Quit();

  return 0;
}
