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
  //
  // NOTE:
  // We intentionally keep the seed string stable so rings remain deterministic
  // across releases. The *look* evolves as the shader-less generator improves,
  // but the seed-space stays coherent.
  core::SplitMix64 rng(core::hashCombine(seed, core::fnv1a64("rings_params")));

  constexpr double twoPi = 2.0 * 3.14159265358979323846;

  // ---- Palette selection ----
  // Rings trend either "warm dusty" or "cold icy". We blend between a low and high
  // albedo tint and then add subtle radial tinting.
  const bool warm = rng.chance(0.52);

  Col3 baseLo = warm ? Col3{0.55, 0.52, 0.48} : Col3{0.50, 0.56, 0.66};
  Col3 baseHi = warm ? Col3{0.92, 0.86, 0.76} : Col3{0.80, 0.86, 0.94};

  // Subtle per-system tint shift.
  const double tint = rng.range(0.85, 1.18);
  baseLo = { clamp01(baseLo.r * tint), clamp01(baseLo.g * tint), clamp01(baseLo.b * tint) };
  baseHi = { clamp01(baseHi.r * tint), clamp01(baseHi.g * tint), clamp01(baseHi.b * tint) };

  const Col3 outerTint = warm ? Col3{0.92, 0.90, 0.86} : Col3{0.78, 0.86, 0.98};
  const Col3 ringletTint = warm ? Col3{1.00, 0.96, 0.90} : Col3{0.92, 0.98, 1.00};

  // ---- Multi-scale radial structure knobs ----
  const double edgeWidth = rng.range(0.016, 0.040);

  const double bandFreqA = rng.range(8.0, 16.0);
  const double bandFreqB = rng.range(22.0, 48.0);
  const double bandSharpA = rng.range(1.6, 3.2);
  const double bandSharpB = rng.range(1.8, 4.6);
  const double bandMix = rng.range(0.55, 0.75);

  // Ringlets: narrow high-density stripes.
  struct Ringlet {
    double c{0.5};
    double w{0.006};
    double strength{0.15};
    double arcChance{0.0};  // chance of becoming an "arc" (azimuth-localized)
    double arcWidthRad{0.5};
    double arcDepth{0.0};   // how much the arc disappears outside its azimuth window
    double arcAngleRad{0.0};
  };

  const int ringletCount = rng.range(18, 44);
  std::vector<Ringlet> ringlets;
  ringlets.reserve((std::size_t)ringletCount);

  for (int i = 0; i < ringletCount; ++i) {
    Ringlet r{};
    // Bias a little toward the mid/outer region (stylized Saturn-ish).
    const double u = rng.nextUnit();
    r.c = std::clamp(0.06 + 0.92 * std::pow(u, rng.range(0.80, 1.35)), 0.02, 0.98);
    r.w = rng.range(0.0018, 0.013);
    r.strength = rng.range(0.06, 0.28);

    // Some ringlets become arcs (localized azimuthal segments).
    r.arcChance = warm ? 0.18 : 0.28;
    if (rng.chance(r.arcChance)) {
      const double arcDeg = rng.range(18.0, 120.0);
      r.arcWidthRad = arcDeg / 180.0 * 3.14159265358979323846;
      r.arcDepth = rng.range(0.45, 0.92);
      r.arcAngleRad = rng.range(0.0, twoPi);
    }
    ringlets.push_back(r);
  }

  // Major divisions + many micro-gaps.
  struct Gap { double c, w, depth; };
  const int majorGapCount = rng.range(1, 2);
  Gap majorGaps[2]{};
  for (int gi = 0; gi < majorGapCount; ++gi) {
    majorGaps[gi].c = rng.range(0.28, 0.90);
    majorGaps[gi].w = rng.range(0.010, 0.055);
    majorGaps[gi].depth = rng.range(0.35, 0.92);
  }

  const int microGapCount = rng.range(10, 34);
  std::vector<Gap> microGaps;
  microGaps.reserve((std::size_t)microGapCount);
  for (int gi = 0; gi < microGapCount; ++gi) {
    Gap g{};
    g.c = rng.range(0.08, 0.96);
    g.w = rng.range(0.0012, 0.010);
    g.depth = rng.range(0.06, 0.38);
    microGaps.push_back(g);
  }

  // Resonance-ish spiral density waves: a small phase ripple that varies with both
  // radius and azimuth. This is a stylized approximation of density waves and
  // wakes observed in planetary rings.
  struct Wave {
    double r0{0.6};   // resonance radius in v-space
    double w{0.06};   // radial envelope width
    int m{5};         // azimuthal mode count (integer -> seamless at u=0/1)
    double amp{0.10}; // modulation amplitude
    double pitch{6.0};
    double phase{0.0};
  };

  const int waveCount = rng.range(1, 3);
  std::vector<Wave> waves;
  waves.reserve((std::size_t)waveCount);
  for (int wi = 0; wi < waveCount; ++wi) {
    Wave w{};
    w.r0 = rng.range(0.18, 0.92);
    w.w = rng.range(0.025, 0.090);
    w.m = rng.range(2, 12);
    w.amp = rng.range(0.05, 0.18) * (rng.chance(0.5) ? 1.0 : -1.0);
    w.pitch = rng.range(2.0, 10.0);
    w.phase = rng.range(0.0, twoPi);
    waves.push_back(w);
  }

  // Occasional spokes: radial dark features (rare).
  struct Spoke {
    double ang{0.0};
    double widthRad{0.03};
    double strength{0.35};
    double rMid{0.55};
    double rWidth{0.22};
  };

  std::vector<Spoke> spokes;
  if (rng.chance(0.22)) {
    const int spokeCount = rng.range(1, 3);
    spokes.reserve((std::size_t)spokeCount);
    for (int si = 0; si < spokeCount; ++si) {
      Spoke s{};
      s.ang = rng.range(0.0, twoPi);
      s.widthRad = (rng.range(0.6, 3.2) / 180.0) * 3.14159265358979323846;
      s.strength = rng.range(0.18, 0.70);
      s.rMid = rng.range(0.35, 0.75);
      s.rWidth = rng.range(0.12, 0.32);
      spokes.push_back(s);
    }
  }

  // Noise seeds for different channels.
  const core::u64 sWarp = core::hashCombine(seed, core::fnv1a64("rings_warp"));
  const core::u64 sGrain = core::hashCombine(seed, core::fnv1a64("rings_grain"));
  const core::u64 sClump = core::hashCombine(seed, core::fnv1a64("rings_clump"));
  const core::u64 sHue = core::hashCombine(seed, core::fnv1a64("rings_hue"));
  const core::u64 sWake = core::hashCombine(seed, core::fnv1a64("rings_wake"));
  const core::u64 sWave = core::hashCombine(seed, core::fnv1a64("rings_wave"));

  // Wake anisotropy basis.
  const double wakePhi = rng.range(0.0, twoPi);
  const double wakeC = std::cos(wakePhi);
  const double wakeS = std::sin(wakePhi);
  const double wakeFreqU = rng.range(120.0, 220.0);
  const double wakeFreqV = rng.range(20.0, 60.0);
  const double wakeAmp = rng.range(0.08, 0.22);

  for (int y = 0; y < heightPx; ++y) {
    const double v = (heightPx <= 1) ? 0.0 : (double)y / (double)(heightPx - 1); // radial [0..1]

    // Smooth fade on both radial edges.
    const double edge = smoothstep(0.0, edgeWidth, v) * (1.0 - smoothstep(1.0 - edgeWidth, 1.0, v));

    // Outer rings trend slightly brighter/colder (visual interest).
    const double outer01 = std::pow(clamp01(v), 0.85);

    for (int x = 0; x < widthPx; ++x) {
      const double u = (widthPx <= 1) ? 0.0 : (double)x / (double)(widthPx - 1); // angle [0..1]
      const double th = u * twoPi;
      const double ca = std::cos(th);
      const double sa = std::sin(th);

      // Periodic warp so the texture is seam-free at u=0/1.
      const double warp = (proc::fbm2D(sWarp, v * 3.2, ca * 1.2 + sa * 1.2, 4) - 0.5) * 0.40;

      // Multi-scale radial banding.
      double bA = 0.5 + 0.5 * std::sin((v * bandFreqA + warp) * twoPi);
      double bB = 0.5 + 0.5 * std::sin((v * bandFreqB + warp * 0.7 + 0.17) * twoPi);
      bA = std::pow(clamp01(bA), bandSharpA);
      bB = std::pow(clamp01(bB), bandSharpB);
      const double bands = clamp01(bandMix * bA + (1.0 - bandMix) * bB);

      // Fine grain + azimuthal clumping.
      const double grain = proc::fbm2D(sGrain, v * 86.0, ca * 12.0 + sa * 12.0, 3);
      const double clump = proc::fbm2D(sClump, v * 2.9, ca * 2.1 + sa * 2.1, 5);

      // Anisotropic "wake" noise (elongated streaks).
      const double px = ca * (v * 2.1 + 0.15);
      const double py = sa * (v * 2.1 + 0.15);
      const double rx = (px * wakeC - py * wakeS);
      const double ry = (px * wakeS + py * wakeC);
      const double wakes = proc::fbm2D(sWake, rx * wakeFreqU, ry * wakeFreqV, 4);

      // ---- Alpha density ----
      double a = edge;
      a *= (0.10 + 0.90 * bands);
      a *= (0.55 + 0.75 * clump);
      a *= (0.78 + 0.48 * grain);
      a *= (1.0 + wakeAmp * (wakes - 0.5) * 2.0);

      // Ringlets: multiplicative boosts (kept gentle).
      double ringletBoost = 0.0;
      double ringletHighlight = 0.0;

      for (const auto& r : ringlets) {
        const double d = (v - r.c) / std::max(1e-6, r.w);
        const double g = std::exp(-0.5 * d * d);

        // Optional "arc" window (localized along azimuth).
        double arcMul = 1.0;
        if (r.arcDepth > 1e-4) {
          const double dotDir = ca * std::cos(r.arcAngleRad) + sa * std::sin(r.arcAngleRad);
          const double cosW = std::cos(std::clamp(r.arcWidthRad, 1e-4, 3.14159265358979323846));
          double t = 0.0;
          if (dotDir > cosW) {
            t = (dotDir - cosW) / std::max(1e-6, (1.0 - cosW));
          }
          // Inside arc => t=1, outside => t=0
          const double w = smoothstep(0.0, 1.0, t);
          arcMul = lerp(1.0 - r.arcDepth, 1.0, w);
        }

        const double s = r.strength * g * arcMul;
        ringletBoost += s;
        ringletHighlight = std::max(ringletHighlight, s);
      }
      a *= (1.0 + ringletBoost);

      // Spiral density waves.
      for (const auto& w : waves) {
        const double dr = (v - w.r0) / std::max(1e-6, w.w);
        const double env = std::exp(-0.5 * dr * dr);
        const double jitter = (proc::fbm2D(sWave, v * 5.0, ca * 1.5 + sa * 1.5, 2) - 0.5) * 0.9;
        const double phase = (double)w.m * th + w.pitch * std::log(std::max(0.02, v)) + w.phase + jitter;
        const double s = std::sin(phase);
        a *= (1.0 + w.amp * s * env);
      }

      // Apply major gaps.
      for (int gi = 0; gi < majorGapCount; ++gi) {
        const double d = (v - majorGaps[gi].c) / std::max(1e-6, majorGaps[gi].w);
        const double g = std::exp(-0.5 * d * d); // gaussian
        a *= (1.0 - majorGaps[gi].depth * g);
      }

      // Apply micro gaps (ringlets being carved out).
      for (const auto& g0 : microGaps) {
        const double d = (v - g0.c) / std::max(1e-6, g0.w);
        const double g = std::exp(-0.5 * d * d);
        a *= (1.0 - g0.depth * g);
      }

      // Spokes: radial darkening in a limited radial band.
      double spokeMask = 0.0;
      for (const auto& s : spokes) {
        const double dotDir = ca * std::cos(s.ang) + sa * std::sin(s.ang);
        const double cosW = std::cos(std::clamp(s.widthRad, 1e-4, 3.14159265358979323846));
        double t = 0.0;
        if (dotDir > cosW) {
          t = (dotDir - cosW) / std::max(1e-6, (1.0 - cosW));
        }
        const double wAng = smoothstep(0.0, 1.0, t);

        const double dr = (v - s.rMid) / std::max(1e-6, s.rWidth);
        const double wRad = std::exp(-0.5 * dr * dr);

        spokeMask = std::max(spokeMask, wAng * wRad);
      }
      if (spokeMask > 1e-6) {
        a *= (1.0 - std::clamp(spokeMask, 0.0, 1.0) * 0.65);
      }

      a = clamp01(a);

      // ---- Color ----
      // Blend palette by band intensity + subtle hue noise.
      const double hueN = proc::fbm2D(sHue, v * 5.0, ca * 1.1 + sa * 1.1, 4);
      double t = clamp01(0.10 + 0.90 * (0.55 * bands + 0.25 * clump + 0.20 * hueN));
      Col3 c = lerpCol(baseLo, baseHi, t);

      // Outer tint (icy edge).
      c = lerpCol(c, outerTint, 0.12 * outer01);

      // Ringlet highlights are slightly brighter/whiter.
      c = lerpCol(c, ringletTint, std::clamp(ringletHighlight * 2.4, 0.0, 0.35));

      // Darken sparse parts to keep rings from looking like a solid disk.
      double lum = lerp(0.28, 1.00, a);
      if (spokeMask > 1e-6) {
        lum *= (1.0 - 0.45 * std::clamp(spokeMask, 0.0, 1.0));
      }
      c.r = clamp01(c.r * lum);
      c.g = clamp01(c.g * lum);
      c.b = clamp01(c.b * lum);

      const std::size_t idx = static_cast<std::size_t>((y * widthPx + x) * 4);
      img.rgba[idx + 0] = toU8(c.r);
      img.rgba[idx + 1] = toU8(c.g);
      img.rgba[idx + 2] = toU8(c.b);
      img.rgba[idx + 3] = toU8(a);
    }
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
