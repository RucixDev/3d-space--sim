#include "stellar/render/ProceduralRings.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/proc/Noise.h"

#include <algorithm>
#include <cmath>

namespace stellar::render {

static inline double clamp01(double v) { return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v); }
static inline double lerp(double a, double b, double t) { return a + (b - a) * t; }

static inline double smoothstep(double edge0, double edge1, double x) {
  if (edge0 == edge1) return (x < edge0) ? 0.0 : 1.0;
  const double t = clamp01((x - edge0) / (edge1 - edge0));
  return t * t * (3.0 - 2.0 * t);
}

struct Col3 {
  double r{1}, g{1}, b{1};
};

static inline Col3 lerpCol(const Col3& a, const Col3& b, double t) {
  return { lerp(a.r, b.r, t), lerp(a.g, b.g, t), lerp(a.b, b.b, t) };
}

static inline std::uint8_t toU8(double v01) {
  const double v = clamp01(v01);
  return static_cast<std::uint8_t>(std::lround(v * 255.0));
}

RingImage generateRingTexture(core::u64 seed, int widthPx, int heightPx) {
  RingImage img;
  widthPx = std::clamp(widthPx, 32, 2048);
  heightPx = std::clamp(heightPx, 16, 1024);
  img.w = widthPx;
  img.h = heightPx;
  img.rgba.assign(static_cast<std::size_t>(widthPx * heightPx * 4), 255);

  // Parameter RNG (kept deterministic by seed).
  core::SplitMix64 rng(core::hashCombine(seed, core::fnv1a64("rings_params")));

  // Palette: warm dusty ... cool icy.
  const bool warm = rng.chance(0.55);
  Col3 baseLo = warm ? Col3{0.55, 0.52, 0.48} : Col3{0.50, 0.56, 0.66};
  Col3 baseHi = warm ? Col3{0.90, 0.84, 0.74} : Col3{0.78, 0.84, 0.92};

  // Subtle tint variation per system.
  const double tint = rng.range(0.85, 1.15);
  baseLo = { clamp01(baseLo.r * tint), clamp01(baseLo.g * tint), clamp01(baseLo.b * tint) };
  baseHi = { clamp01(baseHi.r * tint), clamp01(baseHi.g * tint), clamp01(baseHi.b * tint) };

  // Banding + noise parameters.
  const double bandFreq = rng.range(12.0, 28.0);
  const double bandSharp = rng.range(1.8, 3.8);
  const double edgeWidth = rng.range(0.018, 0.042);

  // One or two big gaps (Cassini-division-ish).
  struct Gap { double c, w, depth; };
  const int gapCount = rng.range(1, 2);
  Gap gaps[2]{};
  for (int gi = 0; gi < gapCount; ++gi) {
    gaps[gi].c = rng.range(0.30, 0.86);
    gaps[gi].w = rng.range(0.012, 0.045);
    gaps[gi].depth = rng.range(0.35, 0.85);
  }

  const double twoPi = 2.0 * 3.14159265358979323846;

  // Noise seeds for different channels.
  const core::u64 sWarp = core::hashCombine(seed, core::fnv1a64("rings_warp"));
  const core::u64 sGrain = core::hashCombine(seed, core::fnv1a64("rings_grain"));
  const core::u64 sClump = core::hashCombine(seed, core::fnv1a64("rings_clump"));
  const core::u64 sHue = core::hashCombine(seed, core::fnv1a64("rings_hue"));

  for (int y = 0; y < heightPx; ++y) {
    const double v = (heightPx <= 1) ? 0.0 : (double)y / (double)(heightPx - 1); // radial [0..1]

    // Smooth fade on both radial edges.
    double edge = smoothstep(0.0, edgeWidth, v) * (1.0 - smoothstep(1.0 - edgeWidth, 1.0, v));

    for (int x = 0; x < widthPx; ++x) {
      const double u = (widthPx <= 1) ? 0.0 : (double)x / (double)(widthPx - 1); // angle [0..1]
      const double th = u * twoPi;
      const double ca = std::cos(th);
      const double sa = std::sin(th);

      // Periodic warp so the texture is seam-free at u=0/1.
      const double warp = (proc::fbm2D(sWarp + 0x11u, v * 6.0, ca * 2.6, 4) +
                           proc::fbm2D(sWarp + 0x22u, v * 6.0, sa * 2.6, 4) - 0.5) * 0.9;

      // Base radial band pattern.
      double bands = 0.5 + 0.5 * std::sin((v * bandFreq + warp) * twoPi);
      bands = std::pow(clamp01(bands), bandSharp);

      // Fine grain + azimuthal clumping.
      const double grain = proc::fbm2D(sGrain, v * 64.0, ca * 10.0 + sa * 10.0, 3);
      const double clump = proc::fbm2D(sClump, v * 2.4, ca * 1.8 + sa * 1.8, 5);

      // Alpha density: bands dominate, clumps modulate.
      double a = edge;
      a *= (0.18 + 0.82 * bands);
      a *= (0.55 + 0.75 * clump);
      a *= (0.75 + 0.55 * grain);

      // Apply a couple of major gaps.
      for (int gi = 0; gi < gapCount; ++gi) {
        const double d = (v - gaps[gi].c) / std::max(1e-6, gaps[gi].w);
        const double g = std::exp(-0.5 * d * d); // gaussian
        a *= (1.0 - gaps[gi].depth * g);
      }

      a = clamp01(a);

      // Color: blend palette by band intensity + subtle hue noise.
      const double hueN = proc::fbm2D(sHue, v * 4.0, ca * 1.0 + sa * 1.0, 4);
      double t = clamp01(0.15 + 0.85 * (0.65 * bands + 0.35 * hueN));
      Col3 c = lerpCol(baseLo, baseHi, t);

      // Darken sparse parts to keep the rings from looking like a solid disk.
      const double lum = lerp(0.30, 1.0, a);
      c.r *= lum;
      c.g *= lum;
      c.b *= lum;

      const std::size_t idx = (static_cast<std::size_t>(y) * static_cast<std::size_t>(widthPx) +
                               static_cast<std::size_t>(x)) * 4u;
      img.rgba[idx + 0] = toU8(c.r);
      img.rgba[idx + 1] = toU8(c.g);
      img.rgba[idx + 2] = toU8(c.b);
      img.rgba[idx + 3] = toU8(a);
    }
  }

  // Hard guarantee the seam is clean: copy column 0 -> last column.
  for (int y = 0; y < heightPx; ++y) {
    const std::size_t row = static_cast<std::size_t>(y) * static_cast<std::size_t>(widthPx) * 4u;
    const std::size_t first = row + 0u;
    const std::size_t last = row + static_cast<std::size_t>(widthPx - 1) * 4u;
    img.rgba[last + 0] = img.rgba[first + 0];
    img.rgba[last + 1] = img.rgba[first + 1];
    img.rgba[last + 2] = img.rgba[first + 2];
    img.rgba[last + 3] = img.rgba[first + 3];
  }

  return img;
}

void RingTextureCache::clear() {
  cache_.clear();
  tick_ = 0;
}

void RingTextureCache::setMaxEntries(std::size_t m) {
  maxEntries_ = std::max<std::size_t>(1u, m);
  evictIfNeeded();
}

core::u64 RingTextureCache::makeKey(core::u64 seed, int w, int h) {
  core::u64 k = core::hashCombine(seed, core::fnv1a64("ring_tex"));
  k = core::hashCombine(k, static_cast<core::u64>(static_cast<std::uint32_t>(w)));
  k = core::hashCombine(k, static_cast<core::u64>(static_cast<std::uint32_t>(h)));
  return k;
}

void RingTextureCache::evictIfNeeded() {
  while (cache_.size() > maxEntries_) {
    // Find least-recently-used entry.
    auto itLRU = cache_.begin();
    for (auto it = cache_.begin(); it != cache_.end(); ++it) {
      if (it->second.lastUseTick < itLRU->second.lastUseTick) itLRU = it;
    }
    cache_.erase(itLRU);
  }
}

const Texture2D& RingTextureCache::get(core::u64 seed, int widthPx, int heightPx) {
  widthPx = std::clamp(widthPx, 32, 2048);
  heightPx = std::clamp(heightPx, 16, 1024);
  const core::u64 key = makeKey(seed, widthPx, heightPx);

  auto it = cache_.find(key);
  if (it != cache_.end()) {
    it->second.lastUseTick = ++tick_;
    return it->second.tex;
  }

  RingImage img = generateRingTexture(seed, widthPx, heightPx);
  Entry e;
  // Linear filtering + mips help the thin banding hold up at distance.
  e.tex.createRGBA(img.w, img.h, img.rgba.data(), true, false, true);
  e.lastUseTick = ++tick_;

  auto [insIt, _] = cache_.emplace(key, std::move(e));
  evictIfNeeded();
  return insIt->second.tex;
}

} // namespace stellar::render
