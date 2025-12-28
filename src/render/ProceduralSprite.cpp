#include "stellar/render/ProceduralSprite.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/proc/Noise.h"

#include <algorithm>
#include <cmath>
#include <queue>

namespace stellar::render {

namespace {

struct FColor {
  float r{0}, g{0}, b{0}, a{1};
};

static inline float clamp01(float v) { return std::clamp(v, 0.0f, 1.0f); }

static inline std::uint8_t toByte(float v) {
  const float c = clamp01(v);
  return (std::uint8_t)std::lround(c * 255.0f);
}

static inline FColor mul(const FColor& c, float s) {
  return {c.r * s, c.g * s, c.b * s, c.a};
}

static inline FColor lerp(const FColor& a, const FColor& b, float t) {
  t = clamp01(t);
  return {a.r + (b.r - a.r) * t,
          a.g + (b.g - a.g) * t,
          a.b + (b.b - a.b) * t,
          a.a + (b.a - a.a) * t};
}

// h,s,v in [0,1]
static inline FColor hsv(float h, float s, float v, float a = 1.0f) {
  h = h - std::floor(h);
  s = clamp01(s);
  v = clamp01(v);

  const float c = v * s;
  const float hh = h * 6.0f;
  const float x = c * (1.0f - std::fabs(std::fmod(hh, 2.0f) - 1.0f));
  float r1=0,g1=0,b1=0;
  if (0.0f <= hh && hh < 1.0f) { r1 = c; g1 = x; b1 = 0; }
  else if (1.0f <= hh && hh < 2.0f) { r1 = x; g1 = c; b1 = 0; }
  else if (2.0f <= hh && hh < 3.0f) { r1 = 0; g1 = c; b1 = x; }
  else if (3.0f <= hh && hh < 4.0f) { r1 = 0; g1 = x; b1 = c; }
  else if (4.0f <= hh && hh < 5.0f) { r1 = x; g1 = 0; b1 = c; }
  else { r1 = c; g1 = 0; b1 = x; }
  const float m = v - c;
  return {r1 + m, g1 + m, b1 + m, a};
}

static inline float smoothstep(float e0, float e1, float x) {
  const float t = clamp01((x - e0) / (e1 - e0));
  return t * t * (3.0f - 2.0f * t);
}

static core::u64 spriteKey(SpriteKind kind, core::u64 seed, int size) {
  core::u64 h = core::hashCombine(core::fnv1a64("sprite"), (core::u64)kind);
  h = core::hashCombine(h, seed);
  h = core::hashCombine(h, (core::u64)(core::i64)size);
  return h;
}

static int neigh8(const std::vector<std::uint8_t>& g, int w, int h, int x, int y) {
  int n = 0;
  for (int oy = -1; oy <= 1; ++oy) {
    for (int ox = -1; ox <= 1; ++ox) {
      if (ox == 0 && oy == 0) continue;
      const int nx = x + ox;
      const int ny = y + oy;
      if (nx < 0 || ny < 0 || nx >= w || ny >= h) continue;
      n += (g[(std::size_t)ny * (std::size_t)w + (std::size_t)nx] != 0) ? 1 : 0;
    }
  }
  return n;
}

static std::vector<std::uint8_t> genSymMask(core::SplitMix64& rng,
                                            int w, int h,
                                            bool mirrorX, bool mirrorY,
                                            float fillProb) {
  std::vector<std::uint8_t> g((std::size_t)w * (std::size_t)h, 0);

  const int halfW = mirrorX ? (w + 1) / 2 : w;
  const int halfH = mirrorY ? (h + 1) / 2 : h;

  for (int y = 0; y < halfH; ++y) {
    for (int x = 0; x < halfW; ++x) {
      const bool filled = rng.chance(fillProb);
      const int xs[2] = {x, w - 1 - x};
      const int ys[2] = {y, h - 1 - y};

      for (int yi = 0; yi < (mirrorY ? 2 : 1); ++yi) {
        for (int xi = 0; xi < (mirrorX ? 2 : 1); ++xi) {
          const int px = xs[xi];
          const int py = ys[yi];
          g[(std::size_t)py * (std::size_t)w + (std::size_t)px] = filled ? 1 : 0;
        }
      }
    }
  }

  // Ensure we have something in the center.
  const int cx = w / 2;
  const int cy = h / 2;
  g[(std::size_t)cy * (std::size_t)w + (std::size_t)cx] = 1;

  return g;
}

static void smoothMask(std::vector<std::uint8_t>& g, int w, int h, int iters) {
  std::vector<std::uint8_t> tmp = g;
  for (int it = 0; it < iters; ++it) {
    tmp = g;
    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        const int n = neigh8(g, w, h, x, y);
        const std::size_t idx = (std::size_t)y * (std::size_t)w + (std::size_t)x;
        if (g[idx]) {
          // "death" rule
          tmp[idx] = (n >= 3) ? 1 : 0;
        } else {
          // "birth" rule
          tmp[idx] = (n >= 5) ? 1 : 0;
        }
      }
    }
    g.swap(tmp);
  }
}

static void keepMainComponent(std::vector<std::uint8_t>& g, int w, int h) {
  const int cx = w / 2;
  const int cy = h / 2;

  // Find a start cell closest to the center.
  int sx = -1, sy = -1;
  int bestD2 = 1 << 30;
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      if (!g[(std::size_t)y * (std::size_t)w + (std::size_t)x]) continue;
      const int dx = x - cx;
      const int dy = y - cy;
      const int d2 = dx * dx + dy * dy;
      if (d2 < bestD2) {
        bestD2 = d2;
        sx = x;
        sy = y;
      }
    }
  }
  if (sx < 0) return;

  std::vector<std::uint8_t> keep((std::size_t)w * (std::size_t)h, 0);

  std::queue<std::pair<int,int>> q;
  q.push({sx, sy});
  keep[(std::size_t)sy * (std::size_t)w + (std::size_t)sx] = 1;

  auto push = [&](int x, int y) {
    if (x < 0 || y < 0 || x >= w || y >= h) return;
    const std::size_t idx = (std::size_t)y * (std::size_t)w + (std::size_t)x;
    if (!g[idx] || keep[idx]) return;
    keep[idx] = 1;
    q.push({x, y});
  };

  while (!q.empty()) {
    auto [x, y] = q.front();
    q.pop();
    push(x + 1, y);
    push(x - 1, y);
    push(x, y + 1);
    push(x, y - 1);
  }

  // Keep only the chosen component.
  for (std::size_t i = 0; i < g.size(); ++i) g[i] = keep[i];
}

static std::vector<std::uint8_t> addOutline(const std::vector<std::uint8_t>& g, int w, int h) {
  std::vector<std::uint8_t> m = g;
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      const std::size_t idx = (std::size_t)y * (std::size_t)w + (std::size_t)x;
      if (!g[idx]) continue;

      bool edge = false;
      const int ox[4] = {1, -1, 0, 0};
      const int oy[4] = {0, 0, 1, -1};
      for (int k = 0; k < 4; ++k) {
        const int nx = x + ox[k];
        const int ny = y + oy[k];
        if (nx < 0 || ny < 0 || nx >= w || ny >= h) { edge = true; break; }
        if (!g[(std::size_t)ny * (std::size_t)w + (std::size_t)nx]) { edge = true; break; }
      }
      if (edge) m[idx] = 2; // outline pixel
    }
  }
  return m;
}

static SpriteImage renderMask(const std::vector<std::uint8_t>& mask,
                              int gridW, int gridH,
                              core::u64 seed, int size,
                              const FColor& base,
                              const FColor& accent,
                              const FColor& outline) {
  SpriteImage img;
  img.w = size;
  img.h = size;
  img.rgba.resize((std::size_t)size * (std::size_t)size * 4, 0);

  const core::u64 shadeSeed = core::hashCombine(seed, core::fnv1a64("shade"));
  const core::u64 hiSeed    = core::hashCombine(seed, core::fnv1a64("highlight"));

  for (int y = 0; y < size; ++y) {
    for (int x = 0; x < size; ++x) {
      const int gx = (x * gridW) / size;
      const int gy = (y * gridH) / size;
      const std::uint8_t m = mask[(std::size_t)gy * (std::size_t)gridW + (std::size_t)gx];

      if (!m) continue;

      FColor c = base;

      // Simple "pixel art" lighting: top-left is brighter.
      const float fx = (size <= 1) ? 0.0f : (float)x / (float)(size - 1);
      const float fy = (size <= 1) ? 0.0f : (float)y / (float)(size - 1);
      const float light = 0.72f + 0.28f * (1.0f - (fx * 0.65f + fy * 0.85f));

      const double n = stellar::proc::smoothNoise2D(shadeSeed, (double)x * 0.17, (double)y * 0.17);
      const float noise = (float)(n - 0.5) * 0.22f;

      const float shade = clamp01(light + noise);

      if (m == 2) {
        c = outline;
      } else {
        c = mul(c, shade);

        // Dither-ish highlights.
        const double hn = stellar::proc::smoothNoise2D(hiSeed, (double)x * 0.31, (double)y * 0.31);
        if (hn > 0.72) {
          c = lerp(c, accent, 0.55f);
        }
      }

      const std::size_t out = ((std::size_t)y * (std::size_t)size + (std::size_t)x) * 4;
      img.rgba[out + 0] = toByte(c.r);
      img.rgba[out + 1] = toByte(c.g);
      img.rgba[out + 2] = toByte(c.b);
      img.rgba[out + 3] = toByte(c.a);
    }
  }

  return img;
}

static SpriteImage genPixelArt(SpriteKind kind, core::u64 seed, int size) {
  core::SplitMix64 rng(seed);

  // Parameterize by kind.
  int gw = 16, gh = 16;
  bool mirrorX = true;
  bool mirrorY = false;
  float fill = 0.42f;
  int smooth = 2;

  switch (kind) {
    case SpriteKind::Commodity:
      gw = 12; gh = 12; mirrorX = true; mirrorY = false; fill = 0.46f; smooth = 2;
      break;
    case SpriteKind::Faction:
      gw = 16; gh = 16; mirrorX = true; mirrorY = true; fill = 0.34f; smooth = 1;
      break;
    case SpriteKind::Mission:
      gw = 14; gh = 14; mirrorX = true; mirrorY = false; fill = 0.40f; smooth = 1;
      break;
    case SpriteKind::Ship:
      gw = 12; gh = 18; mirrorX = true; mirrorY = false; fill = 0.36f; smooth = 1;
      break;
    default:
      break;
  }

  // Palette.
  const float h0 = (float)rng.nextDouble();
  const float s0 = 0.55f + (float)rng.nextDouble() * 0.35f;
  const float v0 = 0.70f + (float)rng.nextDouble() * 0.25f;

  FColor base = hsv(h0, s0, v0, 1.0f);
  FColor accent = hsv(h0 + 0.08f + (float)rng.nextDouble() * 0.12f, std::max(0.35f, s0 - 0.15f), 0.95f, 1.0f);
  FColor outline = mul(base, 0.18f);

  // Certain kinds skew toward "industrial" palettes.
  if (kind == SpriteKind::Station) {
    base = hsv(h0, 0.15f + (float)rng.nextDouble() * 0.15f, 0.78f, 1.0f);
    accent = hsv(h0 + 0.55f, 0.18f, 0.92f, 1.0f);
    outline = mul(base, 0.12f);
  }
  if (kind == SpriteKind::Mission) {
    accent = hsv(h0 + 0.5f, 0.25f + (float)rng.nextDouble() * 0.25f, 1.0f, 1.0f);
  }
  if (kind == SpriteKind::Faction) {
    outline = mul(base, 0.10f);
  }

  // Mask generation.
  std::vector<std::uint8_t> g = genSymMask(rng, gw, gh, mirrorX, mirrorY, fill);

  // Ship silhouettes benefit from a clear spine.
  if (kind == SpriteKind::Ship) {
    const int cx = gw / 2;
    for (int y = gh / 4; y < gh; ++y) {
      g[(std::size_t)y * (std::size_t)gw + (std::size_t)cx] = 1;
    }
  }

  smoothMask(g, gw, gh, smooth);
  keepMainComponent(g, gw, gh);

  // Station: override with a more "mechanical" silhouette.
  if (kind == SpriteKind::Station) {
    g.assign((std::size_t)gw * (std::size_t)gh, 0);
    const int cx = gw / 2;
    const int cy = gh / 2;
    const int hubR = 2 + (int)rng.range(0.0, 2.0);
    for (int y = -hubR; y <= hubR; ++y) {
      for (int x = -hubR; x <= hubR; ++x) {
        const int px = cx + x;
        const int py = cy + y;
        if (px < 0 || py < 0 || px >= gw || py >= gh) continue;
        if (x*x + y*y <= hubR*hubR) g[(std::size_t)py * (std::size_t)gw + (std::size_t)px] = 1;
      }
    }

    const int armLenX = 3 + (int)rng.range(0.0, 5.0);
    const int armLenY = 3 + (int)rng.range(0.0, 5.0);
    const int thick = 1 + (int)rng.range(0.0, 2.0);

    auto drawRect = [&](int x0, int y0, int x1, int y1) {
      if (x0 > x1) std::swap(x0, x1);
      if (y0 > y1) std::swap(y0, y1);
      for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
          if (x < 0 || y < 0 || x >= gw || y >= gh) continue;
          g[(std::size_t)y * (std::size_t)gw + (std::size_t)x] = 1;
        }
      }
    };

    drawRect(cx - armLenX, cy - thick, cx + armLenX, cy + thick);
    drawRect(cx - thick, cy - armLenY, cx + thick, cy + armLenY);

    // End caps.
    drawRect(cx - armLenX - 1, cy - 1, cx - armLenX + 1, cy + 1);
    drawRect(cx + armLenX - 1, cy - 1, cx + armLenX + 1, cy + 1);
    drawRect(cx - 1, cy - armLenY - 1, cx + 1, cy - armLenY + 1);
    drawRect(cx - 1, cy + armLenY - 1, cx + 1, cy + armLenY + 1);

    // Optional ring.
    if (rng.chance(0.55)) {
      const int ringR = 5 + (int)rng.range(0.0, 3.0);
      for (int y = 0; y < gh; ++y) {
        for (int x = 0; x < gw; ++x) {
          const int dx = x - cx;
          const int dy = y - cy;
          const int d2 = dx*dx + dy*dy;
          if (d2 >= ringR*ringR - ringR && d2 <= ringR*ringR + ringR) {
            g[(std::size_t)y * (std::size_t)gw + (std::size_t)x] = 1;
          }
        }
      }
    }

    smoothMask(g, gw, gh, 1);
    keepMainComponent(g, gw, gh);
  }

  std::vector<std::uint8_t> mask = addOutline(g, gw, gh);
  return renderMask(mask, gw, gh, seed, size, base, accent, outline);
}

static SpriteImage genPlanet(core::u64 seed, int size) {
  core::SplitMix64 rng(seed);

  const float h0 = (float)rng.nextDouble();
  const float s0 = 0.35f + (float)rng.nextDouble() * 0.35f;
  const float v0 = 0.65f + (float)rng.nextDouble() * 0.25f;

  const FColor ocean = hsv(h0 + 0.55f, 0.55f, 0.55f, 1.0f);
  const FColor landA = hsv(h0 + 0.08f, s0, v0, 1.0f);
  const FColor landB = hsv(h0 + 0.15f, std::max(0.15f, s0 - 0.15f), std::min(1.0f, v0 + 0.15f), 1.0f);

  const core::u64 nSeed = core::hashCombine(seed, core::fnv1a64("planet"));
  const core::u64 cSeed = core::hashCombine(seed, core::fnv1a64("cloud"));

  SpriteImage img;
  img.w = size;
  img.h = size;
  img.rgba.resize((std::size_t)size * (std::size_t)size * 4, 0);

  // Light direction (normalized).
  float lx = -0.55f, ly = -0.35f, lz = 0.75f;
  const float ll = std::sqrt(lx*lx + ly*ly + lz*lz);
  lx /= ll; ly /= ll; lz /= ll;

  for (int y = 0; y < size; ++y) {
    for (int x = 0; x < size; ++x) {
      const float u = ((float)x + 0.5f) / (float)size * 2.0f - 1.0f;
      const float v = ((float)y + 0.5f) / (float)size * 2.0f - 1.0f;

      const float r2 = u*u + v*v;
      if (r2 > 1.0f) continue;

      const float z = std::sqrt(std::max(0.0f, 1.0f - r2));

      // Terrain noise in [0,1).
      const double n = stellar::proc::fbm2D(nSeed, (double)u * 2.2 + 11.7, (double)v * 2.2 - 3.1, 5);
      const float nn = (float)n;

      // Simple biome mix.
      const float waterLine = 0.48f + 0.10f * (float)rng.nextDouble();
      FColor col = (nn < waterLine) ? ocean : lerp(landA, landB, (nn - waterLine) / std::max(1e-6f, 1.0f - waterLine));

      // Lighting.
      const float ndotl = std::max(0.0f, u*lx + v*ly + z*lz);
      const float ambient = 0.22f;
      const float lit = ambient + (1.0f - ambient) * ndotl;

      col = mul(col, lit);

      // Cloud layer (very cheap).
      const double cn = stellar::proc::fbm2D(cSeed, (double)u * 5.5, (double)v * 5.5, 4);
      if (cn > 0.62) {
        const float ca = (float)((cn - 0.62) / 0.38);
        col = lerp(col, {1, 1, 1, 1}, 0.35f * clamp01(ca));
      }

      // Rim atmosphere.
      const float r = std::sqrt(r2);
      const float rim = smoothstep(0.78f, 1.00f, r);
      const FColor atm = {0.35f, 0.60f, 0.95f, 1.0f};
      col = lerp(col, atm, 0.22f * rim);

      // Alpha falloff at edge.
      const float edge = smoothstep(0.98f, 1.00f, r);
      col.a = 1.0f - edge;

      const std::size_t out = ((std::size_t)y * (std::size_t)size + (std::size_t)x) * 4;
      img.rgba[out + 0] = toByte(col.r);
      img.rgba[out + 1] = toByte(col.g);
      img.rgba[out + 2] = toByte(col.b);
      img.rgba[out + 3] = toByte(col.a);
    }
  }

  return img;
}

static SpriteImage genStar(core::u64 seed, int size) {
  core::SplitMix64 rng(seed);

  const float h0 = (float)rng.nextDouble();
  const FColor coreCol = hsv(h0, 0.25f + (float)rng.nextDouble() * 0.25f, 1.0f, 1.0f);
  const FColor glowCol = hsv(h0, 0.40f, 1.0f, 1.0f);

  const core::u64 raySeed = core::hashCombine(seed, core::fnv1a64("rays"));

  SpriteImage img;
  img.w = size;
  img.h = size;
  img.rgba.resize((std::size_t)size * (std::size_t)size * 4, 0);

  const int rays = 4 + (int)rng.range(0.0, 5.0);

  for (int y = 0; y < size; ++y) {
    for (int x = 0; x < size; ++x) {
      const float u = ((float)x + 0.5f) / (float)size * 2.0f - 1.0f;
      const float v = ((float)y + 0.5f) / (float)size * 2.0f - 1.0f;
      const float r = std::sqrt(u*u + v*v);

      // Core + outer glow.
      const float core = std::exp(-r * 4.5f);
      const float glow = std::exp(-r * 2.2f);

      if (glow < 0.02f) continue;

      // Rays modulate the glow.
      const float ang = std::atan2(v, u);
      const double rn = stellar::proc::smoothNoise2D(raySeed, ang * 2.0, r * 3.0);
      const float ray = 0.55f + 0.45f * std::sin(ang * (float)rays + (float)rn * 2.0f);

      const float a = clamp01(glow * (0.35f + 0.65f * ray));
      FColor col = lerp(glowCol, coreCol, clamp01(core * 1.6f));
      col = mul(col, 0.65f + 0.45f * core);
      col.a = a;

      const std::size_t out = ((std::size_t)y * (std::size_t)size + (std::size_t)x) * 4;
      img.rgba[out + 0] = toByte(col.r);
      img.rgba[out + 1] = toByte(col.g);
      img.rgba[out + 2] = toByte(col.b);
      img.rgba[out + 3] = toByte(col.a);
    }
  }

  return img;
}

} // namespace

SpriteImage generateSprite(SpriteKind kind, core::u64 seed, int size) {
  size = std::clamp(size, 8, 256);

  switch (kind) {
    case SpriteKind::Planet:
      return genPlanet(seed, size);
    case SpriteKind::Star:
      return genStar(seed, size);
    case SpriteKind::Station:
      return genPixelArt(kind, seed, size);
    case SpriteKind::Ship:
    case SpriteKind::Mission:
    case SpriteKind::Faction:
    case SpriteKind::Commodity:
    default:
      return genPixelArt(kind, seed, size);
  }
}

void SpriteCache::clear() { cache_.clear(); }

void SpriteCache::evictIfNeeded() {
  if (maxEntries_ == 0) {
    cache_.clear();
    return;
  }
  if (cache_.size() <= maxEntries_) return;

  // Evict least-recently-used entries. (O(n) scan; cache sizes are small.)
  while (cache_.size() > maxEntries_) {
    auto oldest = cache_.begin();
    for (auto it = cache_.begin(); it != cache_.end(); ++it) {
      if (it->second.lastUseTick < oldest->second.lastUseTick) oldest = it;
    }
    cache_.erase(oldest);
  }
}

const Texture2D& SpriteCache::get(SpriteKind kind, core::u64 seed, int size) {
  const core::u64 k = spriteKey(kind, seed, size);
  ++tick_;

  auto it = cache_.find(k);
  if (it != cache_.end()) {
    it->second.lastUseTick = tick_;
    return it->second.tex;
  }

  SpriteImage img = generateSprite(kind, seed, size);

  Entry e;
  e.lastUseTick = tick_;
  e.tex.createRGBA(img.w, img.h, img.rgba.data(),
                   /*generateMips=*/false,
                   /*nearestFilter=*/true,
                   /*clampToEdge=*/true);

  auto [insIt, ok] = cache_.emplace(k, std::move(e));
  (void)ok;

  evictIfNeeded();
  return insIt->second.tex;
}

} // namespace stellar::render
