#include "stellar/render/ProceduralPlanet.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/proc/Noise.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace stellar::render {

namespace {

constexpr double kPi = 3.1415926535897932384626433832795;
constexpr double kTwoPi = 2.0 * kPi;

struct Vec3 {
  double x{0.0};
  double y{0.0};
  double z{0.0};
};

struct Color {
  double r{0.0};
  double g{0.0};
  double b{0.0};
};

static inline double clamp01(double x) { return std::clamp(x, 0.0, 1.0); }

static inline double lerp(double a, double b, double t) { return a + (b - a) * t; }

static inline Color lerp(Color a, Color b, double t) {
  return {lerp(a.r, b.r, t), lerp(a.g, b.g, t), lerp(a.b, b.b, t)};
}

static inline Color mul(Color c, double k) { return {c.r * k, c.g * k, c.b * k}; }
static inline Color add(Color a, Color b) { return {a.r + b.r, a.g + b.g, a.b + b.b}; }

static inline Color clampColor(Color c) {
  c.r = clamp01(c.r);
  c.g = clamp01(c.g);
  c.b = clamp01(c.b);
  return c;
}

static inline void writeRGBA(std::vector<std::uint8_t>& rgba, std::size_t idx, Color c, double a = 1.0) {
  c = clampColor(c);
  a = clamp01(a);
  rgba[idx + 0] = static_cast<std::uint8_t>(std::lround(c.r * 255.0));
  rgba[idx + 1] = static_cast<std::uint8_t>(std::lround(c.g * 255.0));
  rgba[idx + 2] = static_cast<std::uint8_t>(std::lround(c.b * 255.0));
  rgba[idx + 3] = static_cast<std::uint8_t>(std::lround(a * 255.0));
}

static inline double smoothstep(double e0, double e1, double x) {
  if (e0 == e1) return (x < e0) ? 0.0 : 1.0;
  const double t = clamp01((x - e0) / (e1 - e0));
  return t * t * (3.0 - 2.0 * t);
}

static inline double fade(double t) {
  // Perlin fade polynomial: 6t^5 - 15t^4 + 10t^3
  return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
}

static inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

static inline Vec3 add3(const Vec3& a, const Vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }

static inline Vec3 mul3(const Vec3& a, double k) { return {a.x * k, a.y * k, a.z * k}; }

static inline double lengthSq(const Vec3& v) { return dot(v, v); }

static inline Vec3 normalize(const Vec3& v) {
  const double lsq = lengthSq(v);
  if (lsq <= 1e-18) return {0.0, 1.0, 0.0};
  const double inv = 1.0 / std::sqrt(lsq);
  return {v.x * inv, v.y * inv, v.z * inv};
}

static inline core::u64 tagSeed(core::u64 seed, const char* tag) {
  return core::hashCombine(seed, core::fnv1a64(tag));
}

// Smooth 3D value-noise sampled at fractional coordinates (trilinear with fade).
static double smoothNoise3D(core::u64 seed, double x, double y, double z) {
  const int x0 = static_cast<int>(std::floor(x));
  const int y0 = static_cast<int>(std::floor(y));
  const int z0 = static_cast<int>(std::floor(z));
  const int x1 = x0 + 1;
  const int y1 = y0 + 1;
  const int z1 = z0 + 1;

  const double tx = x - (double)x0;
  const double ty = y - (double)y0;
  const double tz = z - (double)z0;

  const double n000 = proc::valueNoise3D(seed, x0, y0, z0);
  const double n100 = proc::valueNoise3D(seed, x1, y0, z0);
  const double n010 = proc::valueNoise3D(seed, x0, y1, z0);
  const double n110 = proc::valueNoise3D(seed, x1, y1, z0);
  const double n001 = proc::valueNoise3D(seed, x0, y0, z1);
  const double n101 = proc::valueNoise3D(seed, x1, y0, z1);
  const double n011 = proc::valueNoise3D(seed, x0, y1, z1);
  const double n111 = proc::valueNoise3D(seed, x1, y1, z1);

  const double u = fade(tx);
  const double v = fade(ty);
  const double w = fade(tz);

  const double x00 = lerp(n000, n100, u);
  const double x10 = lerp(n010, n110, u);
  const double x01 = lerp(n001, n101, u);
  const double x11 = lerp(n011, n111, u);

  const double y0v = lerp(x00, x10, v);
  const double y1v = lerp(x01, x11, v);
  return lerp(y0v, y1v, w);
}

static double fbm3D(core::u64 seed,
                    double x, double y, double z,
                    int octaves,
                    double lacunarity,
                    double gain) {
  double amp = 0.5;
  double freq = 1.0;
  double sum = 0.0;
  for (int i = 0; i < octaves; ++i) {
    sum += amp * smoothNoise3D(seed + static_cast<core::u64>(i) * 1013ull, x * freq, y * freq, z * freq);
    freq *= lacunarity;
    amp *= gain;
  }
  return sum;
}

static inline double ampSum(int octaves, double gain) {
  double amp = 0.5;
  double s = 0.0;
  for (int i = 0; i < octaves; ++i) {
    s += amp;
    amp *= gain;
  }
  return s;
}

static inline double fbm01(core::u64 seed,
                          double x, double y, double z,
                          int octaves,
                          double lacunarity,
                          double gain) {
  if (octaves <= 0) return 0.5;
  const double s = std::max(1e-9, ampSum(octaves, gain));
  return clamp01(fbm3D(seed, x, y, z, octaves, lacunarity, gain) / s);
}

static inline double ridged01(double n01) {
  // n01 in [0,1] -> ridged in [0,1]
  return 1.0 - std::abs(n01 * 2.0 - 1.0);
}

static inline Vec3 noiseVec3(core::u64 seed, const Vec3& p, double freq, int octaves) {
  // 3 correlated FBMs with different seed offsets.
  const double x = fbm01(seed ^ 0xA11CE5u, p.x * freq, p.y * freq, p.z * freq, octaves, 2.0, 0.5) - 0.5;
  const double y = fbm01(seed ^ 0xBADC0FFEu, p.x * freq, p.y * freq, p.z * freq, octaves, 2.0, 0.5) - 0.5;
  const double z = fbm01(seed ^ 0xC0FFEEu, p.x * freq, p.y * freq, p.z * freq, octaves, 2.0, 0.5) - 0.5;
  return {x, y, z};
}

static inline Vec3 warpDir(core::u64 seed, const Vec3& p, double freq, double amp) {
  const Vec3 w = noiseVec3(seed, p, freq, 4);
  return normalize(add3(p, mul3(w, amp)));
}

static inline Vec3 randUnitVec(core::SplitMix64& rng) {
  // Uniform sphere distribution via (z, phi).
  const double z = rng.nextDouble() * 2.0 - 1.0;
  const double phi = rng.nextDouble() * kTwoPi;
  const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
  return {r * std::cos(phi), z, r * std::sin(phi)};
}

struct Plate {
  Vec3 dir{};
  double heightBias{0.5};
  double tempBias{0.0};
  double moistBias{0.0};
};

struct Crater {
  Vec3 dir{};
  double cosRadius{0.0};
  double invSpan{1.0};
  double depth{0.0};
  double rim{0.0};
};

struct Vortex {
  Vec3 dir{};
  double cosRadius{0.0};
  double invSpan{1.0};
  double strength{0.0};
  Color tint{1.0, 1.0, 1.0};
};

struct SurfaceContext {
  SurfaceKind kind{SurfaceKind::Rocky};
  core::u64 seed{0};

  int paletteIx{0};
  double seaLevel{0.54};

  std::vector<Plate> plates;
  std::vector<Crater> craters;
  std::vector<Vortex> vortices; // gas giants + clouds
};

static std::vector<Plate> makePlates(SurfaceKind kind, core::u64 seed) {
  std::vector<Plate> out;

  // No plates needed for stars.
  if (kind == SurfaceKind::Star) return out;

  core::SplitMix64 rng(tagSeed(seed, "plates"));

  int base = 20;
  switch (kind) {
    case SurfaceKind::Rocky: base = 26; break;
    case SurfaceKind::Desert: base = 22; break;
    case SurfaceKind::Ocean: base = 18; break;
    case SurfaceKind::Ice: base = 20; break;
    case SurfaceKind::GasGiant: base = 12; break;
    case SurfaceKind::Clouds: base = 12; break;
    default: base = 20; break;
  }

  const int count = std::clamp(base + rng.range(-4, 6), 8, 48);
  out.reserve((std::size_t)count);

  for (int i = 0; i < count; ++i) {
    Plate p;
    p.dir = randUnitVec(rng);

    if (kind == SurfaceKind::Ocean) {
      // Bias towards lower values so oceans dominate.
      const double t = rng.nextDouble();
      p.heightBias = std::pow(t, 1.35);
      // Climate biases to create continent-scale patterns.
      p.tempBias = rng.range(-0.12, 0.12);
      p.moistBias = rng.range(-0.20, 0.20);
    } else if (kind == SurfaceKind::GasGiant || kind == SurfaceKind::Clouds) {
      p.heightBias = rng.range(0.35, 0.65);
      p.tempBias = rng.range(-0.08, 0.08);
      p.moistBias = rng.range(-0.08, 0.08);
    } else {
      p.heightBias = rng.range(0.25, 0.85);
      p.tempBias = rng.range(-0.10, 0.10);
      p.moistBias = rng.range(-0.10, 0.10);
    }

    out.push_back(p);
  }

  return out;
}

static std::vector<Crater> makeCraters(SurfaceKind kind, core::u64 seed) {
  std::vector<Crater> out;

  // Oceans and gas giants don't get a visible crater field here.
  if (kind == SurfaceKind::Ocean || kind == SurfaceKind::GasGiant || kind == SurfaceKind::Star || kind == SurfaceKind::Clouds) {
    return out;
  }

  core::SplitMix64 rng(tagSeed(seed, "craters"));

  int baseCount = 48;
  double rMinDeg = 4.0;
  double rMaxDeg = 18.0;
  double depthBase = 0.085;
  double rim = 0.25;

  if (kind == SurfaceKind::Desert) {
    baseCount = 38;
    rMinDeg = 5.0;
    rMaxDeg = 20.0;
    depthBase = 0.070;
    rim = 0.22;
  } else if (kind == SurfaceKind::Ice) {
    baseCount = 52;
    rMinDeg = 6.0;
    rMaxDeg = 22.0;
    depthBase = 0.080;
    rim = 0.28;
  }

  const int count = std::clamp(baseCount + rng.range(-10, 16), 0, 96);
  if (count <= 0 || depthBase <= 0.0) return out;

  out.reserve((std::size_t)count);

  const double rMin = (rMinDeg * kPi) / 180.0;
  const double rMax = (rMaxDeg * kPi) / 180.0;

  for (int i = 0; i < count; ++i) {
    const Vec3 d = randUnitVec(rng);
    const double rad = rng.range(rMin, rMax);
    const double cosR = std::cos(rad);
    const double invSpan = 1.0 / std::max(1e-6, (1.0 - cosR));

    Crater c;
    c.dir = d;
    c.cosRadius = cosR;
    c.invSpan = invSpan;
    c.depth = depthBase * rng.range(0.65, 1.15);
    c.rim = std::clamp(rim * rng.range(0.85, 1.15), 0.0, 1.0);
    out.push_back(c);
  }

  return out;
}

static std::vector<Vortex> makeVortices(SurfaceKind kind, core::u64 seed) {
  std::vector<Vortex> out;
  if (kind != SurfaceKind::GasGiant && kind != SurfaceKind::Clouds) return out;

  core::SplitMix64 rng(tagSeed(seed, kind == SurfaceKind::GasGiant ? "gas_vortices" : "cloud_vortices"));

  const int count = std::clamp(3 + rng.range(0, 5), 0, 10);
  out.reserve((std::size_t)count);

  auto pickTint = [&](int i) -> Color {
    // A small palette of storm tints.
    const int v = (i + (int)(seed & 7ull)) % 5;
    switch (v) {
      case 0: return {0.95, 0.60, 0.45}; // red-ish
      case 1: return {0.85, 0.92, 1.00}; // pale blue
      case 2: return {1.00, 0.92, 0.70}; // creamy
      case 3: return {0.72, 0.78, 0.86}; // cool gray
      default: return {0.92, 0.82, 0.62};
    }
  };

  for (int i = 0; i < count; ++i) {
    const Vec3 d = randUnitVec(rng);

    // Keep most vortices away from the poles (aesthetic).
    const double y = std::clamp(d.y, -0.92, 0.92);
    const double phi = std::atan2(d.z, d.x);
    const double r = std::sqrt(std::max(0.0, 1.0 - y * y));
    Vec3 dir{r * std::cos(phi), y, r * std::sin(phi)};

    const double radDeg = rng.range(8.0, 28.0);
    const double rad = (radDeg * kPi) / 180.0;

    Vortex v;
    v.dir = normalize(dir);
    v.cosRadius = std::cos(rad);
    v.invSpan = 1.0 / std::max(1e-6, (1.0 - v.cosRadius));
    v.strength = rng.range(0.35, 0.85);
    v.tint = pickTint(i);
    out.push_back(v);
  }

  return out;
}

struct PlateSample {
  int idx{0};
  double boundary01{0.0};
  double bestDot{-2.0};
  double secondDot{-2.0};
};

static PlateSample samplePlates(const SurfaceContext& ctx, const Vec3& p) {
  PlateSample s;
  if (ctx.plates.empty()) return s;

  double best = -2.0;
  double second = -2.0;
  int bestIdx = 0;

  for (int i = 0; i < (int)ctx.plates.size(); ++i) {
    const double d = dot(p, ctx.plates[(std::size_t)i].dir);
    if (d > best) {
      second = best;
      best = d;
      bestIdx = i;
    } else if (d > second) {
      second = d;
    }
  }

  const double diff = best - second;
  // Near the boundary, best and second-best are close.
  const double t = clamp01((0.11 - diff) / 0.11);
  const double boundary = smoothstep(0.0, 1.0, t);

  s.idx = bestIdx;
  s.boundary01 = boundary;
  s.bestDot = best;
  s.secondDot = second;
  return s;
}

static double applyCraters(const SurfaceContext& ctx,
                           const Vec3& p,
                           double h,
                           double strengthMul,
                           double* outCraterInterior01,
                           double* outCraterRim01) {
  if (outCraterInterior01) *outCraterInterior01 = 0.0;
  if (outCraterRim01) *outCraterRim01 = 0.0;
  if (ctx.craters.empty() || strengthMul <= 0.0) return h;

  double interior = 0.0;
  double rim = 0.0;

  for (const auto& c : ctx.craters) {
    const double d = dot(p, c.dir); // cos(angle)
    if (d <= c.cosRadius) continue;

    const double t = clamp01((d - c.cosRadius) * c.invSpan); // 0=edge -> 1=center
    const double w = smoothstep(0.0, 1.0, t);

    // Depression.
    h -= c.depth * strengthMul * w;
    interior = std::max(interior, w);

    // Ridge ring near the edge.
    if (c.rim > 0.0) {
      const double ring = smoothstep(0.06, 0.22, t) * (1.0 - smoothstep(0.22, 0.58, t));
      h += c.depth * strengthMul * c.rim * ring;
      rim = std::max(rim, ring);
    }
  }

  if (outCraterInterior01) *outCraterInterior01 = clamp01(interior);
  if (outCraterRim01) *outCraterRim01 = clamp01(rim);
  return clamp01(h);
}

static Color biomeColor(double temp01, double moist01) {
  temp01 = clamp01(temp01);
  moist01 = clamp01(moist01);

  // Coarse biome palette (stylized but plausible).
  const Color desert{0.86, 0.78, 0.45};
  const Color savanna{0.62, 0.70, 0.32};
  const Color grass{0.34, 0.60, 0.26};
  const Color forest{0.16, 0.42, 0.18};
  const Color jungle{0.12, 0.50, 0.22};
  const Color steppe{0.55, 0.50, 0.26};
  const Color tundra{0.62, 0.62, 0.54};

  if (temp01 < 0.22) {
    // Polar/tundra.
    return lerp(tundra, {0.90, 0.95, 1.00}, smoothstep(0.00, 0.20, (0.22 - temp01) / 0.22));
  }

  if (temp01 < 0.42) {
    // Cool climates.
    if (moist01 > 0.55) return lerp(steppe, forest, (moist01 - 0.55) / 0.45);
    return lerp(tundra, steppe, moist01 / 0.55);
  }

  if (temp01 < 0.72) {
    // Temperate.
    if (moist01 > 0.70) return lerp(forest, jungle, (moist01 - 0.70) / 0.30);
    if (moist01 > 0.45) return lerp(grass, forest, (moist01 - 0.45) / 0.25);
    if (moist01 > 0.25) return lerp(savanna, grass, (moist01 - 0.25) / 0.20);
    return lerp(desert, savanna, moist01 / 0.25);
  }

  // Tropical.
  if (moist01 > 0.70) return jungle;
  if (moist01 > 0.48) return lerp(savanna, jungle, (moist01 - 0.48) / 0.22);
  if (moist01 > 0.30) return lerp(desert, savanna, (moist01 - 0.30) / 0.18);
  return desert;
}

// --- Height + albedo functions ------------------------------------------------

static double surfaceHeightRocky(const SurfaceContext& ctx, const Vec3& p, double lat, double lon) {
  (void)lat;
  (void)lon;
  const Vec3 q = warpDir(tagSeed(ctx.seed, "rock_warp"), p, 2.1, 0.22);
  const PlateSample ps = samplePlates(ctx, q);
  const double plateH = ctx.plates.empty() ? 0.55 : ctx.plates[(std::size_t)ps.idx].heightBias;

  const double macro = fbm01(tagSeed(ctx.seed, "rock_macro"), q.x * 2.8, q.y * 2.8, q.z * 2.8, 6, 2.0, 0.5);
  const double rN = fbm01(tagSeed(ctx.seed, "rock_ridged"), q.x * 7.5, q.y * 7.5, q.z * 7.5, 4, 2.1, 0.55);
  const double ridged = ridged01(rN);
  const double dust = fbm01(tagSeed(ctx.seed, "rock_dust"), q.x * 15.0, q.y * 15.0, q.z * 15.0, 3, 2.4, 0.55);

  const double mountain = std::pow(ps.boundary01, 1.65) * (0.18 + 0.22 * dust);

  double h = clamp01(0.52 * plateH + 0.30 * macro + 0.18 * ridged + mountain);
  h = clamp01(h * (0.92 + 0.10 * dust));

  h = applyCraters(ctx, p, h, 1.0, nullptr, nullptr);
  return h;
}

static Color surfaceRocky(const SurfaceContext& ctx, const Vec3& p, double lat, double lon) {
  (void)lat;
  (void)lon;

  const Vec3 q = warpDir(tagSeed(ctx.seed, "rock_warp"), p, 2.1, 0.22);
  const PlateSample ps = samplePlates(ctx, q);

  double craterIn = 0.0;
  double craterRim = 0.0;

  // Recompute height with crater masks (keeps albedo + normal maps visually aligned).
  const double plateH = ctx.plates.empty() ? 0.55 : ctx.plates[(std::size_t)ps.idx].heightBias;
  const double macro = fbm01(tagSeed(ctx.seed, "rock_macro"), q.x * 2.8, q.y * 2.8, q.z * 2.8, 6, 2.0, 0.5);
  const double rN = fbm01(tagSeed(ctx.seed, "rock_ridged"), q.x * 7.5, q.y * 7.5, q.z * 7.5, 4, 2.1, 0.55);
  const double ridged = ridged01(rN);
  const double dust = fbm01(tagSeed(ctx.seed, "rock_dust"), q.x * 15.0, q.y * 15.0, q.z * 15.0, 3, 2.4, 0.55);
  const double mountain = std::pow(ps.boundary01, 1.65) * (0.18 + 0.22 * dust);

  double h = clamp01(0.52 * plateH + 0.30 * macro + 0.18 * ridged + mountain);
  h = clamp01(h * (0.92 + 0.10 * dust));
  h = applyCraters(ctx, p, h, 1.0, &craterIn, &craterRim);

  Color baseLo{0.18, 0.16, 0.15};
  Color baseHi{0.62, 0.56, 0.48};

  if ((ctx.paletteIx % 3) == 1) {
    baseLo = {0.20, 0.12, 0.10};
    baseHi = {0.70, 0.46, 0.32};
  } else if ((ctx.paletteIx % 3) == 2) {
    baseLo = {0.10, 0.10, 0.12};
    baseHi = {0.46, 0.46, 0.50};
  }

  Color col = lerp(baseLo, baseHi, h);

  // Shadow valleys a bit.
  col = lerp(col, mul(col, 0.65), (1.0 - ridged) * 0.35);
  // Dusty variation.
  col = mul(col, (0.85 + 0.28 * dust));

  // Plate boundaries -> bright mountains.
  col = lerp(col, {0.74, 0.72, 0.70}, ps.boundary01 * 0.18);

  // Craters: dark interior + bright rim.
  col = lerp(col, mul(col, 0.55), craterIn * 0.55);
  col = lerp(col, {0.95, 0.93, 0.90}, craterRim * 0.28);

  // Speckles.
  const double speck = fbm01(tagSeed(ctx.seed, "rock_speck"), q.x * 26.0, q.y * 26.0, q.z * 26.0, 2, 2.5, 0.6);
  const double speckMask = smoothstep(0.84, 0.97, speck);
  col = lerp(col, {0.88, 0.84, 0.78}, speckMask * 0.32);

  return col;
}

static double surfaceHeightDesert(const SurfaceContext& ctx, const Vec3& p, double lat, double lon) {
  const Vec3 q = warpDir(tagSeed(ctx.seed, "des_warp"), p, 1.8, 0.18);
  const PlateSample ps = samplePlates(ctx, q);
  const double plateH = ctx.plates.empty() ? 0.55 : ctx.plates[(std::size_t)ps.idx].heightBias;

  const double macro = fbm01(tagSeed(ctx.seed, "des_macro"), q.x * 2.2, q.y * 2.2, q.z * 2.2, 6, 2.0, 0.5);
  const double warp = fbm01(tagSeed(ctx.seed, "des_warp2"), q.x * 3.4, q.y * 3.4, q.z * 3.4, 4, 2.1, 0.55);

  const double band = 0.5 + 0.5 * std::sin((lat * 16.0) + (lon * 2.2) + (warp - 0.5) * 3.2);

  const double rock = fbm01(tagSeed(ctx.seed, "des_rock"), q.x * 6.0, q.y * 6.0, q.z * 6.0, 4, 2.1, 0.52);
  const double rockMask = smoothstep(0.70, 0.88, rock);

  const double mountain = std::pow(ps.boundary01, 1.55) * (0.12 + 0.16 * rock);

  double h = clamp01(0.55 * plateH + 0.22 * macro + 0.20 * band + 0.08 * rockMask + mountain);
  h = applyCraters(ctx, p, h, 0.75, nullptr, nullptr);
  return h;
}

static Color surfaceDesert(const SurfaceContext& ctx, const Vec3& p, double lat, double lon) {
  const Vec3 q = warpDir(tagSeed(ctx.seed, "des_warp"), p, 1.8, 0.18);
  const PlateSample ps = samplePlates(ctx, q);

  double craterIn = 0.0;
  double craterRim = 0.0;

  const double plateH = ctx.plates.empty() ? 0.55 : ctx.plates[(std::size_t)ps.idx].heightBias;
  const double macro = fbm01(tagSeed(ctx.seed, "des_macro"), q.x * 2.2, q.y * 2.2, q.z * 2.2, 6, 2.0, 0.5);
  const double warp = fbm01(tagSeed(ctx.seed, "des_warp2"), q.x * 3.4, q.y * 3.4, q.z * 3.4, 4, 2.1, 0.55);
  const double band = 0.5 + 0.5 * std::sin((lat * 16.0) + (lon * 2.2) + (warp - 0.5) * 3.2);
  const double rock = fbm01(tagSeed(ctx.seed, "des_rock"), q.x * 6.0, q.y * 6.0, q.z * 6.0, 4, 2.1, 0.52);
  const double rockMask = smoothstep(0.70, 0.88, rock);
  const double mountain = std::pow(ps.boundary01, 1.55) * (0.12 + 0.16 * rock);

  double h = clamp01(0.55 * plateH + 0.22 * macro + 0.20 * band + 0.08 * rockMask + mountain);
  h = applyCraters(ctx, p, h, 0.75, &craterIn, &craterRim);

  Color sandLo{0.62, 0.52, 0.28};
  Color sandHi{0.92, 0.86, 0.50};
  if ((ctx.paletteIx & 1) != 0) {
    sandLo = {0.56, 0.34, 0.22};
    sandHi = {0.88, 0.62, 0.40};
  }

  Color col = lerp(sandLo, sandHi, clamp01(0.35 + 0.65 * macro));
  col = mul(col, (0.90 + 0.18 * band));

  // Rock outcrops.
  col = lerp(col, {0.34, 0.28, 0.22}, rockMask * 0.65);

  // Plate boundary mountains.
  col = lerp(col, {0.70, 0.64, 0.52}, ps.boundary01 * 0.16);

  // Craters.
  col = lerp(col, mul(col, 0.62), craterIn * 0.50);
  col = lerp(col, {0.98, 0.94, 0.86}, craterRim * 0.18);

  // Subtle darker dunes at low "height".
  col = lerp(col, mul(col, 0.92), (1.0 - h) * 0.25);

  return col;
}

static double surfaceHeightOcean(const SurfaceContext& ctx, const Vec3& p, double lat, double lon) {
  (void)lon;

  const Vec3 q = warpDir(tagSeed(ctx.seed, "ocean_warp"), p, 1.6, 0.20);
  const PlateSample ps = samplePlates(ctx, q);
  const double plateH = ctx.plates.empty() ? 0.50 : ctx.plates[(std::size_t)ps.idx].heightBias;

  const double macro = fbm01(tagSeed(ctx.seed, "ocean_macro"), q.x * 2.1, q.y * 2.1, q.z * 2.1, 6, 2.0, 0.5);
  const double detail = fbm01(tagSeed(ctx.seed, "ocean_detail"), q.x * 6.5, q.y * 6.5, q.z * 6.5, 5, 2.1, 0.55);

  // Tectonic uplift near plate boundaries.
  const double mountain = std::pow(ps.boundary01, 1.85) * (0.18 + 0.26 * detail);

  double elevRaw = clamp01(0.68 * plateH + 0.32 * macro + mountain);

  const double sea = ctx.seaLevel;
  const double land = smoothstep(sea - 0.02, sea + 0.02, elevRaw);

  const double landElev = clamp01((elevRaw - sea) / std::max(1e-6, (1.0 - sea)));

  // Waves (subtle; we don't want oceans to look like terrain).
  const double wave = fbm01(tagSeed(ctx.seed, "ocean_wave"), q.x * 18.0, q.y * 18.0, q.z * 18.0, 3, 2.2, 0.55);
  const double hOcean = 0.03 * (wave - 0.5) * (0.85 + 0.15 * std::abs(std::sin(lat)));

  const double hLand = 0.18 + 0.82 * landElev;

  return clamp01(lerp(hOcean, hLand, land));
}

static Color surfaceOcean(const SurfaceContext& ctx, const Vec3& p, double lat, double lon) {
  const Vec3 q = warpDir(tagSeed(ctx.seed, "ocean_warp"), p, 1.6, 0.20);
  const PlateSample ps = samplePlates(ctx, q);
  const Plate plate = ctx.plates.empty() ? Plate{} : ctx.plates[(std::size_t)ps.idx];

  const double plateH = ctx.plates.empty() ? 0.50 : plate.heightBias;
  const double macro = fbm01(tagSeed(ctx.seed, "ocean_macro"), q.x * 2.1, q.y * 2.1, q.z * 2.1, 6, 2.0, 0.5);
  const double detail = fbm01(tagSeed(ctx.seed, "ocean_detail"), q.x * 6.5, q.y * 6.5, q.z * 6.5, 5, 2.1, 0.55);

  const double mountain = std::pow(ps.boundary01, 1.85) * (0.18 + 0.26 * detail);
  double elevRaw = clamp01(0.68 * plateH + 0.32 * macro + mountain);

  const double sea = ctx.seaLevel;
  const double land = smoothstep(sea - 0.02, sea + 0.02, elevRaw);
  const double landElev = clamp01((elevRaw - sea) / std::max(1e-6, (1.0 - sea)));

  // Coastline band.
  const double coast = clamp01(1.0 - std::abs(elevRaw - sea) * 22.0);

  Color waterDeep{0.02, 0.07, 0.20};
  Color waterShallow{0.07, 0.26, 0.52};
  if ((ctx.paletteIx & 2) != 0) {
    // Slightly greener oceans.
    waterDeep = {0.03, 0.09, 0.18};
    waterShallow = {0.08, 0.30, 0.40};
  }

  Color water = lerp(waterDeep, waterShallow, coast);

  // Climate model (stylized): temperature from latitude + altitude; moisture from noise.
  const double lat01 = std::abs(std::sin(lat)); // 0 at equator, 1 at poles

  double temp = 1.0 - std::pow(lat01, 0.95);
  temp += plate.tempBias;
  temp -= landElev * 0.35;
  temp = clamp01(temp);

  double moist = fbm01(tagSeed(ctx.seed, "ocean_moist"), q.x * 3.5, q.y * 3.5, q.z * 3.5, 5, 2.0, 0.5);
  moist = clamp01((moist - 0.35) * 1.45);
  // Hadley-ish belt: moister near equator and mid-latitudes.
  moist *= (0.62 + 0.38 * std::cos(lat * 2.0));
  moist = clamp01(moist + plate.moistBias);

  Color landCol = biomeColor(temp, moist);

  // Beaches.
  if (coast > 0.0) {
    const double beach = coast * (1.0 - smoothstep(sea + 0.01, sea + 0.04, elevRaw));
    if (beach > 1e-6) {
      landCol = lerp(landCol, {0.92, 0.86, 0.62}, beach);
    }
  }

  // Mountains: gray rock.
  const double mountainMask = std::max(smoothstep(0.72, 0.92, landElev), ps.boundary01 * 0.55);
  landCol = lerp(landCol, {0.52, 0.50, 0.48}, mountainMask * 0.55);

  // Snow/ice caps.
  const double snowLat = smoothstep(0.62, 0.88, lat01);
  const double snowElev = smoothstep(0.80, 0.98, landElev);
  const double snow = clamp01(std::max(snowLat, snowElev) * (0.55 + 0.45 * (1.0 - temp)));
  landCol = lerp(landCol, {0.90, 0.95, 1.00}, snow);

  // Combine land + water.
  Color col = lerp(water, landCol, land);

  // Subtle cloud reflection on ocean.
  const double cloud = fbm01(tagSeed(ctx.seed, "ocean_cloud_hint"), q.x * 4.0, q.y * 4.0, q.z * 4.0, 5, 2.0, 0.5);
  const double cloudMask = smoothstep(0.74, 0.93, cloud) * (1.0 - land * 0.35);
  col = lerp(col, {1.0, 1.0, 1.0}, cloudMask * 0.08);

  return col;
}

static double surfaceHeightIce(const SurfaceContext& ctx, const Vec3& p, double lat, double lon) {
  (void)lon;
  const Vec3 q = warpDir(tagSeed(ctx.seed, "ice_warp"), p, 1.9, 0.16);
  const PlateSample ps = samplePlates(ctx, q);
  const double plateH = ctx.plates.empty() ? 0.55 : ctx.plates[(std::size_t)ps.idx].heightBias;

  const double n = fbm01(tagSeed(ctx.seed, "ice_macro"), q.x * 3.3, q.y * 3.3, q.z * 3.3, 6, 2.0, 0.5);
  const double crack = fbm01(tagSeed(ctx.seed, "ice_crack"), q.x * 20.0, q.y * 20.0, q.z * 20.0, 3, 2.4, 0.55);
  const double crackMask = smoothstep(0.78, 0.93, crack);

  const double ridge = std::pow(ps.boundary01, 1.35) * 0.28;

  double h = clamp01(0.52 * plateH + 0.32 * n + 0.16 * crackMask + ridge);

  // Shallow waves/ripples near the equator.
  const double eq = 1.0 - std::abs(std::sin(lat));
  const double ripple = fbm01(tagSeed(ctx.seed, "ice_ripple"), q.x * 9.0, q.y * 9.0, q.z * 9.0, 4, 2.1, 0.55);
  h = clamp01(h + (ripple - 0.5) * 0.05 * eq);

  h = applyCraters(ctx, p, h, 0.85, nullptr, nullptr);
  return h;
}

static Color surfaceIce(const SurfaceContext& ctx, const Vec3& p, double lat, double lon) {
  (void)lon;
  const Vec3 q = warpDir(tagSeed(ctx.seed, "ice_warp"), p, 1.9, 0.16);
  const PlateSample ps = samplePlates(ctx, q);

  double craterIn = 0.0;
  double craterRim = 0.0;

  const double plateH = ctx.plates.empty() ? 0.55 : ctx.plates[(std::size_t)ps.idx].heightBias;
  const double n = fbm01(tagSeed(ctx.seed, "ice_macro"), q.x * 3.3, q.y * 3.3, q.z * 3.3, 6, 2.0, 0.5);
  const double crack = fbm01(tagSeed(ctx.seed, "ice_crack"), q.x * 20.0, q.y * 20.0, q.z * 20.0, 3, 2.4, 0.55);
  const double crackMask = smoothstep(0.78, 0.93, crack);
  const double ridge = std::pow(ps.boundary01, 1.35) * 0.28;

  double h = clamp01(0.52 * plateH + 0.32 * n + 0.16 * crackMask + ridge);
  h = applyCraters(ctx, p, h, 0.85, &craterIn, &craterRim);

  const Color iceLo{0.70, 0.82, 0.95};
  const Color iceHi{0.93, 0.98, 1.00};
  Color col = lerp(iceLo, iceHi, clamp01(0.15 + 0.85 * n));

  // Cracks.
  col = lerp(col, {0.28, 0.32, 0.38}, crackMask * 0.35);

  // Dust near equator.
  const double dust = fbm01(tagSeed(ctx.seed, "ice_dust"), q.x * 10.0, q.y * 10.0, q.z * 10.0, 4, 2.1, 0.55);
  const double eq = 1.0 - std::abs(std::sin(lat));
  const double dustMask = clamp01((dust - 0.45) * 1.45) * eq;
  col = lerp(col, {0.62, 0.56, 0.44}, dustMask * 0.35);

  // Plate boundary ridges.
  col = lerp(col, {0.86, 0.90, 0.98}, ps.boundary01 * 0.10);

  // Crater shading.
  col = lerp(col, mul(col, 0.75), craterIn * 0.55);
  col = lerp(col, {1.00, 1.00, 1.00}, craterRim * 0.18);

  // Poles are cleaner.
  const double pole = smoothstep(0.62, 0.92, std::abs(std::sin(lat)));
  col = lerp(col, {0.92, 0.98, 1.00}, pole * 0.15);

  return col;
}

static double surfaceHeightGasGiant(const SurfaceContext& ctx, const Vec3& p, double lat, double lon) {
  // Domain warp gives "fluid" swirl.
  const Vec3 q = warpDir(tagSeed(ctx.seed, "gas_warp"), p, 2.6, 0.25);
  const double warp = fbm01(tagSeed(ctx.seed, "gas_warpn"), q.x * 2.6, q.y * 2.6, q.z * 2.6, 5, 2.0, 0.55);

  const double bandFreq = 20.0 + 10.0 * warp;
  const double band = 0.5 + 0.5 * std::sin((lat * bandFreq) + (warp - 0.5) * 4.2);
  const double turb = fbm01(tagSeed(ctx.seed, "gas_turb"), q.x * 7.0, q.y * 7.0, q.z * 7.0, 4, 2.1, 0.55);

  double h = 0.5 + 0.10 * (band - 0.5) + 0.08 * (turb - 0.5);

  // Vortex height bumps.
  for (const auto& v : ctx.vortices) {
    const double d = dot(p, v.dir);
    if (d <= v.cosRadius) continue;
    const double t = clamp01((d - v.cosRadius) * v.invSpan);
    const double w = smoothstep(0.0, 1.0, t);
    const double swirl = 0.5 + 0.5 * std::sin(lon * 3.2 + lat * 5.1 + (double)(ctx.paletteIx) * 0.7);
    h += (swirl - 0.5) * 0.08 * w * v.strength;
  }

  return clamp01(h);
}

static Color surfaceGasGiant(const SurfaceContext& ctx, const Vec3& p, double lat, double lon) {
  const Vec3 q = warpDir(tagSeed(ctx.seed, "gas_warp"), p, 2.6, 0.25);
  const double warp = fbm01(tagSeed(ctx.seed, "gas_warpn"), q.x * 2.6, q.y * 2.6, q.z * 2.6, 5, 2.0, 0.55);

  // Compress along Y to emphasize horizontal bands.
  const double bandFreq = 22.0 + 12.0 * warp;
  double band = 0.5 + 0.5 * std::sin((lat * bandFreq) + (warp - 0.5) * 4.2);
  const double fine = fbm01(tagSeed(ctx.seed, "gas_fine"), q.x * 12.0, q.y * 12.0, q.z * 12.0, 3, 2.3, 0.55);
  const double turb = fbm01(tagSeed(ctx.seed, "gas_turb"), q.x * 7.0, q.y * 7.0, q.z * 7.0, 4, 2.1, 0.55);

  double t = clamp01(0.60 * band + 0.25 * turb + 0.15 * fine);

  // Base palette varies slightly.
  Color lo{0.84, 0.68, 0.42};
  Color hi{0.96, 0.88, 0.66};
  if ((ctx.paletteIx & 1) != 0) {
    lo = {0.70, 0.76, 0.86};
    hi = {0.95, 0.94, 0.88};
  }

  Color col = lerp(lo, hi, t);

  // Polar cooling.
  const double pole = smoothstep(0.55, 0.92, std::abs(std::sin(lat)));
  col = lerp(col, {0.68, 0.76, 0.88}, pole * 0.18);

  // Vortices (spots).
  for (const auto& v : ctx.vortices) {
    const double d = dot(p, v.dir);
    if (d <= v.cosRadius) continue;
    const double u = clamp01((d - v.cosRadius) * v.invSpan);
    const double w = smoothstep(0.0, 1.0, u);

    // Add some swirl to the band intensity.
    const double swirl = 0.5 + 0.5 * std::sin(lon * 2.8 + lat * 7.5 + (double)(ctx.paletteIx) * 0.9);
    const double local = clamp01(t + (swirl - 0.5) * 0.25 * w * v.strength);
    col = lerp(col, lerp(lo, hi, local), 0.35 * w * v.strength);

    // Tint the vortex.
    col = lerp(col, v.tint, w * v.strength * 0.65);
  }

  return col;
}

static double surfaceHeightStar(const SurfaceContext& ctx, const Vec3& p, double lat, double lon) {
  (void)lat;
  (void)lon;
  const double coarse = fbm01(tagSeed(ctx.seed, "star_coarse"), p.x * 5.0, p.y * 5.0, p.z * 5.0, 5, 2.0, 0.55);
  const double fine = fbm01(tagSeed(ctx.seed, "star_fine"), p.x * 18.0, p.y * 18.0, p.z * 18.0, 3, 2.4, 0.55);
  const double h = 0.55 + 0.55 * coarse + 0.20 * (fine - 0.5);
  return clamp01(h);
}

static Color surfaceStar(const SurfaceContext& ctx, const Vec3& p, double lat, double lon) {
  (void)lat;
  (void)lon;

  // Star surface: hotter/brighter granulation + subtle limb variation encoded in the texture.
  const double coarse = fbm01(tagSeed(ctx.seed, "star_coarse"), p.x * 5.0, p.y * 5.0, p.z * 5.0, 5, 2.0, 0.55);
  const double fine = fbm01(tagSeed(ctx.seed, "star_fine"), p.x * 18.0, p.y * 18.0, p.z * 18.0, 3, 2.4, 0.55);

  const double t = clamp01(0.65 * coarse + 0.35 * fine);

  // Warm palette.
  Color c0{1.00, 0.70, 0.18};
  Color c1{1.00, 0.92, 0.58};
  if ((ctx.paletteIx & 1) != 0) {
    // Slightly whiter star.
    c0 = {1.00, 0.84, 0.52};
    c1 = {1.00, 0.98, 0.86};
  }

  Color col = lerp(c0, c1, t);

  // Subtle darker spots.
  const double spot = fbm01(tagSeed(ctx.seed, "star_spot"), p.x * 9.0, p.y * 9.0, p.z * 9.0, 4, 2.0, 0.5);
  const double spotMask = smoothstep(0.78, 0.92, spot);
  col = lerp(col, mul(col, 0.82), spotMask * 0.25);

  return col;
}

static double surfaceCloudAlpha(const SurfaceContext& ctx, const Vec3& p, double lat, double lon) {
  const Vec3 q = warpDir(tagSeed(ctx.seed, "cloud_warp"), p, 2.2, 0.22);

  // Cloud density based on layered FBM.
  const double n1 = fbm01(tagSeed(ctx.seed, "cloud_n1"), q.x * 3.6, q.y * 3.6, q.z * 3.6, 6, 2.0, 0.5);
  const double n2 = fbm01(tagSeed(ctx.seed, "cloud_n2"), q.x * 10.5 + lat * 0.35, q.y * 10.5, q.z * 10.5 + lon * 0.15, 4, 2.1, 0.55);
  const double n = 0.65 * n1 + 0.35 * n2;
  double d = smoothstep(0.55, 0.82, n);

  // Latitude banding (simple circulation cells).
  const double bandWarp = fbm01(tagSeed(ctx.seed, "cloud_band"), q.x * 2.8, q.y * 2.8, q.z * 2.8, 4, 2.0, 0.55);
  const double band = 0.5 + 0.5 * std::sin(lat * (10.0 + 8.0 * bandWarp) + (n2 - 0.5) * 3.0);
  d *= 0.78 + 0.22 * band;

  // Add a few storm systems (vortices).
  for (const auto& v : ctx.vortices) {
    const double dd = dot(p, v.dir);
    if (dd <= v.cosRadius) continue;
    const double u = clamp01((dd - v.cosRadius) * v.invSpan);
    const double w = smoothstep(0.0, 1.0, u);
    const double swirl = 0.5 + 0.5 * std::sin(lon * 4.5 + lat * 7.0 + (double)ctx.paletteIx);
    d += (swirl * 0.55 + 0.45) * w * v.strength * 0.30;
  }

  return clamp01(d);
}

static double surfaceHeight(const SurfaceContext& ctx, const Vec3& p, double lat, double lon) {
  switch (ctx.kind) {
    case SurfaceKind::Rocky: return surfaceHeightRocky(ctx, p, lat, lon);
    case SurfaceKind::Desert: return surfaceHeightDesert(ctx, p, lat, lon);
    case SurfaceKind::Ocean: return surfaceHeightOcean(ctx, p, lat, lon);
    case SurfaceKind::Ice: return surfaceHeightIce(ctx, p, lat, lon);
    case SurfaceKind::GasGiant: return surfaceHeightGasGiant(ctx, p, lat, lon);
    case SurfaceKind::Star: return surfaceHeightStar(ctx, p, lat, lon);
    case SurfaceKind::Clouds: return surfaceCloudAlpha(ctx, p, lat, lon);
    default: return 0.5;
  }
}

static Color surfaceAlbedo(const SurfaceContext& ctx, const Vec3& p, double lat, double lon, double* outAlpha) {
  if (outAlpha) *outAlpha = 1.0;

  switch (ctx.kind) {
    case SurfaceKind::Rocky:
      return surfaceRocky(ctx, p, lat, lon);
    case SurfaceKind::Desert:
      return surfaceDesert(ctx, p, lat, lon);
    case SurfaceKind::Ocean:
      return surfaceOcean(ctx, p, lat, lon);
    case SurfaceKind::Ice:
      return surfaceIce(ctx, p, lat, lon);
    case SurfaceKind::GasGiant:
      return surfaceGasGiant(ctx, p, lat, lon);
    case SurfaceKind::Star:
      return surfaceStar(ctx, p, lat, lon);
    case SurfaceKind::Clouds: {
      if (outAlpha) *outAlpha = surfaceCloudAlpha(ctx, p, lat, lon);
      return {0.94, 0.97, 1.00};
    }
    default:
      return {0.6, 0.6, 0.6};
  }
}

static SurfaceContext makeContext(SurfaceKind kind, core::u64 seed) {
  SurfaceContext ctx;
  ctx.kind = kind;
  ctx.seed = seed;

  // Small palette selector. (Used for rocky/desert/ocean/gas/star variants.)
  ctx.paletteIx = (int)((seed >> 21) & 7ull);

  // Ocean sea level per seed.
  if (kind == SurfaceKind::Ocean) {
    const double t = (double)((seed >> 11) & 0xFFull) / 255.0;
    ctx.seaLevel = 0.50 + t * 0.09; // 0.50 .. 0.59
  }

  ctx.plates = makePlates(kind, seed);
  ctx.craters = makeCraters(kind, seed);
  ctx.vortices = makeVortices(kind, seed);

  return ctx;
}

} // namespace

SurfaceImage generateSurfaceTexture(SurfaceKind kind, core::u64 seed, int widthPx) {
  SurfaceImage img{};
  if (widthPx <= 0) return img;

  img.w = std::max(4, widthPx);
  img.h = std::max(2, img.w / 2);
  img.rgba.assign(static_cast<std::size_t>(img.w * img.h * 4), 255);

  const SurfaceContext ctx = makeContext(kind, seed);

  // Deterministic rotation so systems don't all align on the seam.
  const double rot = (double)((seed >> 8) & 0xFFFFull) / 65535.0 * kTwoPi;

  for (int y = 0; y < img.h; ++y) {
    const double v = (double)(y + 0.5) / (double)img.h;
    const double lat = (v - 0.5) * kPi; // -pi/2..pi/2
    const double cosLat = std::cos(lat);
    const double sinLat = std::sin(lat);

    for (int x = 0; x < img.w; ++x) {
      const double u = (double)(x + 0.5) / (double)img.w;
      const double lon = u * kTwoPi - kPi + rot;

      const double cosLon = std::cos(lon);
      const double sinLon = std::sin(lon);

      // Unit sphere sample point.
      const Vec3 p{cosLat * cosLon, sinLat, cosLat * sinLon};

      double a = 1.0;
      const Color c = surfaceAlbedo(ctx, p, lat, lon, &a);

      const std::size_t idx = (static_cast<std::size_t>(y) * static_cast<std::size_t>(img.w) + static_cast<std::size_t>(x)) * 4;
      writeRGBA(img.rgba, idx, c, a);
    }
  }

  return img;
}

SurfaceImage generateSurfaceNormalMap(SurfaceKind kind, core::u64 seed, int widthPx) {
  SurfaceImage img{};
  if (widthPx <= 0) return img;

  img.w = std::max(4, widthPx);
  img.h = std::max(2, img.w / 2);
  img.rgba.assign(static_cast<std::size_t>(img.w * img.h * 4), 255);

  const SurfaceContext ctx = makeContext(kind, seed);

  // Keep the same deterministic seam rotation as the albedo so features line up.
  const double rot = (double)((seed >> 8) & 0xFFFFull) / 65535.0 * kTwoPi;

  // Precompute height map in one pass. This makes normal-map generation ~4x cheaper
  // for heavier procedural surfaces (plates/craters).
  std::vector<double> heights;
  heights.resize((std::size_t)img.w * (std::size_t)img.h);

  for (int y = 0; y < img.h; ++y) {
    const double v = (double)(y + 0.5) / (double)img.h;
    const double lat = (v - 0.5) * kPi;
    const double cosLat = std::cos(lat);
    const double sinLat = std::sin(lat);

    for (int x = 0; x < img.w; ++x) {
      const double u = (double)(x + 0.5) / (double)img.w;
      const double lon = u * kTwoPi - kPi + rot;

      const double cosLon = std::cos(lon);
      const double sinLon = std::sin(lon);

      const Vec3 p{cosLat * cosLon, sinLat, cosLat * sinLon};
      heights[(std::size_t)y * (std::size_t)img.w + (std::size_t)x] = surfaceHeight(ctx, p, lat, lon);
    }
  }

  const auto idxWrapX = [&](int x) {
    while (x < 0) x += img.w;
    while (x >= img.w) x -= img.w;
    return x;
  };
  const auto idxClampY = [&](int y) { return std::clamp(y, 0, img.h - 1); };

  const double du = 1.0 / (double)img.w;
  const double dv = 1.0 / (double)img.h;

  // Per-kind slope scaling (the shader also has a global normal strength multiplier).
  double kindScale = 1.0;
  switch (kind) {
    case SurfaceKind::Rocky: kindScale = 1.25; break;
    case SurfaceKind::Desert: kindScale = 1.05; break;
    case SurfaceKind::Ocean: kindScale = 0.85; break;
    case SurfaceKind::Ice: kindScale = 1.10; break;
    case SurfaceKind::GasGiant: kindScale = 0.22; break;
    case SurfaceKind::Star: kindScale = 0.34; break;
    case SurfaceKind::Clouds: kindScale = 0.60; break;
    default: kindScale = 1.0; break;
  }

  for (int y = 0; y < img.h; ++y) {
    const double v = (double)(y + 0.5) / (double)img.h;
    const double lat = (v - 0.5) * kPi;
    const double cosLat = std::cos(lat);

    // Physical distances on a unit sphere for one texel step.
    const double dLon = du * kTwoPi;
    const double dLat = dv * kPi;

    const double dx = std::max(1.0e-6, std::abs(cosLat) * dLon);
    const double dy = std::max(1.0e-6, dLat);

    for (int x = 0; x < img.w; ++x) {
      const int xL = idxWrapX(x - 1);
      const int xR = idxWrapX(x + 1);
      const int yD = idxClampY(y - 1);
      const int yU = idxClampY(y + 1);

      const double hL = heights[(std::size_t)y * (std::size_t)img.w + (std::size_t)xL];
      const double hR = heights[(std::size_t)y * (std::size_t)img.w + (std::size_t)xR];
      const double hD = heights[(std::size_t)yD * (std::size_t)img.w + (std::size_t)x];
      const double hU = heights[(std::size_t)yU * (std::size_t)img.w + (std::size_t)x];

      const double dHdx = (hR - hL) / (2.0 * dx);
      const double dHdy = (hU - hD) / (2.0 * dy);

      // Tangent-space normal: +Z points "out of the surface".
      double nx = -dHdx * kindScale;
      double ny = -dHdy * kindScale;
      double nz = 1.0;

      const double invLen = 1.0 / std::sqrt(nx * nx + ny * ny + nz * nz + 1e-12);
      nx *= invLen;
      ny *= invLen;
      nz *= invLen;

      const std::size_t idx = (static_cast<std::size_t>(y) * static_cast<std::size_t>(img.w) + static_cast<std::size_t>(x)) * 4;
      img.rgba[idx + 0] = static_cast<std::uint8_t>(std::lround((nx * 0.5 + 0.5) * 255.0));
      img.rgba[idx + 1] = static_cast<std::uint8_t>(std::lround((ny * 0.5 + 0.5) * 255.0));
      img.rgba[idx + 2] = static_cast<std::uint8_t>(std::lround((nz * 0.5 + 0.5) * 255.0));
      img.rgba[idx + 3] = 255;
    }
  }

  return img;
}

void SurfaceTextureCache::clear() {
  cache_.clear();
  tick_ = 0;
}

core::u64 SurfaceTextureCache::makeKey(SurfaceKind kind, core::u64 seed, int widthPx) const {
  core::u64 h = core::fnv1a64("surface_tex");
  h = core::hashCombine(h, (core::u64)(core::u8)kind);
  h = core::hashCombine(h, seed);
  h = core::hashCombine(h, (core::u64)(core::i64)widthPx);
  return h;
}

void SurfaceTextureCache::evictIfNeeded() {
  if (maxEntries_ == 0) {
    cache_.clear();
    return;
  }
  while (cache_.size() > maxEntries_) {
    core::u64 oldestKey = 0;
    core::u64 oldestTick = (core::u64)-1;
    for (const auto& kv : cache_) {
      if (kv.second.lastUseTick < oldestTick) {
        oldestTick = kv.second.lastUseTick;
        oldestKey = kv.first;
      }
    }
    if (oldestKey != 0) cache_.erase(oldestKey);
    else break;
  }
}

const Texture2D& SurfaceTextureCache::get(SurfaceKind kind, core::u64 seed, int widthPx) {
  ++tick_;
  widthPx = std::clamp(widthPx, 64, 2048);

  const core::u64 k = makeKey(kind, seed, widthPx);
  if (auto it = cache_.find(k); it != cache_.end()) {
    it->second.lastUseTick = tick_;
    return it->second.tex;
  }

  // Generate.
  SurfaceImage img = generateSurfaceTexture(kind, seed, widthPx);

  Entry e{};
  // For planets/stars we want smooth sampling + mipmaps to avoid shimmering.
  // We repeat on both axes; the noise is sampled on the unit sphere so U=0/1 is seamless.
  e.tex.createRGBA(img.w, img.h, img.rgba.data(),
                   /*generateMips=*/true,
                   /*nearestFilter=*/false,
                   /*clampToEdge=*/false);
  e.lastUseTick = tick_;

  auto [it, inserted] = cache_.emplace(k, std::move(e));
  (void)inserted;

  evictIfNeeded();
  return it->second.tex;
}

void SurfaceNormalMapCache::clear() {
  cache_.clear();
  tick_ = 0;
}

core::u64 SurfaceNormalMapCache::makeKey(SurfaceKind kind, core::u64 seed, int widthPx) const {
  core::u64 h = core::fnv1a64("surface_nrm");
  h = core::hashCombine(h, (core::u64)(core::u8)kind);
  h = core::hashCombine(h, seed);
  h = core::hashCombine(h, (core::u64)(core::i64)widthPx);
  return h;
}

void SurfaceNormalMapCache::evictIfNeeded() {
  if (maxEntries_ == 0) {
    cache_.clear();
    return;
  }
  while (cache_.size() > maxEntries_) {
    core::u64 oldestKey = 0;
    core::u64 oldestTick = (core::u64)-1;
    for (const auto& kv : cache_) {
      if (kv.second.lastUseTick < oldestTick) {
        oldestTick = kv.second.lastUseTick;
        oldestKey = kv.first;
      }
    }
    if (oldestKey != 0) cache_.erase(oldestKey);
    else break;
  }
}

const Texture2D& SurfaceNormalMapCache::get(SurfaceKind kind, core::u64 seed, int widthPx) {
  ++tick_;
  widthPx = std::clamp(widthPx, 64, 2048);

  const core::u64 k = makeKey(kind, seed, widthPx);
  if (auto it = cache_.find(k); it != cache_.end()) {
    it->second.lastUseTick = tick_;
    return it->second.tex;
  }

  SurfaceImage img = generateSurfaceNormalMap(kind, seed, widthPx);

  Entry e{};
  // Normal maps want linear sampling + mips; repeat on both axes.
  e.tex.createRGBA(img.w, img.h, img.rgba.data(),
                   /*generateMips=*/true,
                   /*nearestFilter=*/false,
                   /*clampToEdge=*/false);
  e.lastUseTick = tick_;

  auto [it, inserted] = cache_.emplace(k, std::move(e));
  (void)inserted;

  evictIfNeeded();
  return it->second.tex;
}

} // namespace stellar::render
