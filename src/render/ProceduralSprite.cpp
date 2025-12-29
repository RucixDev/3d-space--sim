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

static SpriteImage genCargo(core::u64 seed, int size) {
  // Sci-fi cargo pod / container icon.
  core::SplitMix64 rng(seed);

  // Palette: generally neutral greys with a seed-based accent.
  const float hAcc = (float)rng.nextDouble();
  const float sAcc = 0.45f + (float)rng.nextDouble() * 0.40f;

  FColor base = hsv(hAcc, 0.08f + (float)rng.nextDouble() * 0.10f, 0.72f + (float)rng.nextDouble() * 0.16f, 1.0f);
  FColor accent = hsv(hAcc + 0.52f, sAcc, 0.98f, 1.0f);
  FColor outline = mul(base, 0.10f);

  const bool hazard = rng.chance(0.30);
  const bool band = rng.chance(0.65);
  const bool label = rng.chance(0.40);

  const core::u64 nSeed = core::hashCombine(seed, core::fnv1a64("cargo_noise"));

  SpriteImage img;
  img.w = size;
  img.h = size;
  img.rgba.resize((std::size_t)size * (std::size_t)size * 4, 0);

  // Rounded-rect body in normalized coordinates.
  const float hx = 0.46f; // half-width
  const float hy = 0.24f; // half-height
  const float rad = 0.16f;

  for (int y = 0; y < size; ++y) {
    for (int x = 0; x < size; ++x) {
      const float u = ((float)x + 0.5f) / (float)size * 2.0f - 1.0f;
      const float v = ((float)y + 0.5f) / (float)size * 2.0f - 1.0f;

      // Signed distance to rounded rectangle.
      const float qx = std::abs(u) - hx;
      const float qy = std::abs(v) - hy;
      const float ax = std::max(qx, 0.0f);
      const float ay = std::max(qy, 0.0f);
      const float outside = std::sqrt(ax * ax + ay * ay);
      const float inside = std::min(std::max(qx, qy), 0.0f);
      const float dist = outside + inside - rad;

      const float distPx = dist * (float)size * 0.5f;
      const float alpha = clamp01(0.5f - distPx);
      if (alpha <= 0.0f) continue;

      // Edge / outline.
      const bool isEdge = (distPx > -0.9f && distPx < 0.9f);

      // Simple shading + noise.
      const float fx = (size <= 1) ? 0.0f : (float)x / (float)(size - 1);
      const float fy = (size <= 1) ? 0.0f : (float)y / (float)(size - 1);
      const float light = 0.72f + 0.28f * (1.0f - (fx * 0.45f + fy * 0.85f));
      const double n = stellar::proc::fbm2D(nSeed, (double)u * 3.1, (double)v * 3.1, 4, 2.0, 0.5);
      const float noise = (float)(n - 0.5) * 0.10f;
      const float shade = clamp01(light + noise);

      FColor c = isEdge ? outline : mul(base, shade);
      c.a = alpha;

      // Central band / stripes.
      if (!isEdge && band && std::abs(v) < 0.08f) {
        if (hazard) {
          // Diagonal-ish stripes in UV space.
          const float t = u * 9.0f + v * 13.0f;
          const int k = (int)std::floor(t);
          const float m = (k % 2 == 0) ? 0.65f : 0.15f;
          c = lerp(c, accent, m);
        } else {
          c = lerp(c, accent, 0.55f);
        }
      }

      // Panel seams.
      if (!isEdge) {
        const float seam = 0.012f;
        if (std::abs(u) < seam || std::abs(u - 0.20f) < seam || std::abs(u + 0.20f) < seam) {
          c = lerp(c, outline, 0.65f);
        }
      }

      // Tiny label dot (like a serial/marking).
      if (!isEdge && label) {
        const float dx = u - 0.22f;
        const float dy = v + 0.11f;
        if (dx * dx + dy * dy < 0.0028f) {
          c = lerp(c, accent, 0.90f);
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

static SpriteImage genAsteroid(core::u64 seed, int size) {
  // Irregular rocky asteroid icon with sphere-ish shading.
  core::SplitMix64 rng(seed);

  const float hBase = 0.06f + 0.08f * (float)rng.nextDouble(); // brown/ochre range
  const float sBase = 0.10f + 0.12f * (float)rng.nextDouble();
  const float vBase = 0.55f + 0.18f * (float)rng.nextDouble();

  FColor base = hsv(hBase, sBase, vBase, 1.0f);
  FColor accent = hsv(hBase + 0.02f, std::min(0.35f, sBase + 0.12f), 0.92f, 1.0f);
  FColor outline = mul(base, 0.10f);

  const core::u64 edgeSeed = core::hashCombine(seed, core::fnv1a64("ast_edge"));
  const core::u64 bumpSeed = core::hashCombine(seed, core::fnv1a64("ast_bump"));
  const core::u64 shadeSeed = core::hashCombine(seed, core::fnv1a64("ast_shade"));

  SpriteImage img;
  img.w = size;
  img.h = size;
  img.rgba.resize((std::size_t)size * (std::size_t)size * 4, 0);

  const float rBase = 0.62f;
  // Light from top-left.
  float lx = -0.55f, ly = -0.45f, lz = 0.70f;
  {
    const float inv = 1.0f / std::sqrt(lx*lx + ly*ly + lz*lz);
    lx *= inv; ly *= inv; lz *= inv;
  }

  for (int y = 0; y < size; ++y) {
    for (int x = 0; x < size; ++x) {
      const float u = ((float)x + 0.5f) / (float)size * 2.0f - 1.0f;
      const float v = ((float)y + 0.5f) / (float)size * 2.0f - 1.0f;
      const float rr = std::sqrt(u*u + v*v);
      if (rr > 1.05f) continue;

      const float ang = std::atan2(v, u);
      const double nAng = stellar::proc::smoothNoise2D(edgeSeed, std::cos(ang) * 1.7, std::sin(ang) * 1.7);
      const double nXY  = stellar::proc::smoothNoise2D(edgeSeed, (double)u * 2.2, (double)v * 2.2);
      const float rEdge = rBase * (0.80f + 0.38f * (float)nAng) * (0.92f + 0.18f * (float)(nXY - 0.5));

      const float dist = rr - rEdge;
      const float distPx = dist * (float)size * 0.5f;
      const float alpha = clamp01(0.5f - distPx);
      if (alpha <= 0.0f) continue;

      const bool isEdge = (distPx > -1.1f && distPx < 0.9f);

      // Normal for a deformed sphere.
      const float z = std::sqrt(std::max(0.0f, rEdge*rEdge - rr*rr));
      float nx = (rEdge > 1e-6f) ? (u / rEdge) : 0.0f;
      float ny = (rEdge > 1e-6f) ? (v / rEdge) : 0.0f;
      float nz = (rEdge > 1e-6f) ? (z / rEdge) : 1.0f;

      // Bumpy lighting.
      const double bump = stellar::proc::fbm2D(bumpSeed, (double)u * 3.4, (double)v * 3.4, 4, 2.0, 0.5);
      nz = clamp01(nz + (float)(bump - 0.5) * 0.55f);

      const float ndotl = clamp01(nx * lx + ny * ly + nz * lz);
      const float rim = clamp01((rr / std::max(1e-6f, rEdge)));
      const float ao = 0.58f + 0.42f * (1.0f - rim);

      const double grain = stellar::proc::smoothNoise2D(shadeSeed, (double)u * 6.0, (double)v * 6.0);
      const float g = 0.92f + 0.16f * (float)(grain - 0.5);

      FColor c = isEdge ? outline : mul(base, (0.35f + 0.65f * ndotl) * ao * g);
      c.a = alpha;

      // Occasional bright "chip" highlight.
      if (!isEdge && ndotl > 0.85f && grain > 0.72) {
        c = lerp(c, accent, 0.65f);
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

static SpriteImage genSignal(core::u64 seed, int size) {
  // Holo "radar ping" / signal icon.
  core::SplitMix64 rng(seed);

  const float h = 0.42f + 0.20f * (float)rng.nextDouble(); // cyan/blue range
  const float s = 0.55f + 0.35f * (float)rng.nextDouble();
  const float v0 = 0.92f;

  FColor col = hsv(h, s, v0, 1.0f);
  FColor col2 = hsv(h + 0.08f, std::max(0.35f, s - 0.20f), 1.0f, 1.0f);
  FColor outline = mul(col, 0.20f);

  const float phase = (float)rng.nextDouble() * 6.2831853f;
  const float arcGate = 0.45f + 0.35f * (float)rng.nextDouble();
  const bool sweep = rng.chance(0.55);

  const core::u64 nSeed = core::hashCombine(seed, core::fnv1a64("sig_noise"));

  SpriteImage img;
  img.w = size;
  img.h = size;
  img.rgba.resize((std::size_t)size * (std::size_t)size * 4, 0);

  const float R = 0.78f;

  for (int y = 0; y < size; ++y) {
    for (int x = 0; x < size; ++x) {
      const float u = ((float)x + 0.5f) / (float)size * 2.0f - 1.0f;
      const float v = ((float)y + 0.5f) / (float)size * 2.0f - 1.0f;
      const float rr = std::sqrt(u*u + v*v);
      if (rr > 1.05f) continue;

      const float dist = rr - R;
      const float distPx = dist * (float)size * 0.5f;
      const float outer = clamp01(0.5f - distPx);
      if (outer <= 0.0f) continue;

      // Rings.
      float ring = 0.0f;
      auto addRing = [&](float r0, float w) {
        const float d = rr - r0;
        ring += std::exp(-(d*d) * w);
      };
      addRing(0.05f, 220.0f);
      addRing(0.24f, 70.0f);
      addRing(0.42f, 70.0f);
      addRing(0.60f, 70.0f);
      addRing(0.76f, 100.0f);

      // Arc gating to avoid perfect circles.
      float gate = 1.0f;
      const float ang = std::atan2(v, u) + phase;
      if (rr > 0.18f) {
        const float s1 = std::sin(ang * 2.0f);
        const float s2 = std::sin(ang * 5.0f + 0.7f);
        const float g = 0.55f + 0.25f * s1 + 0.20f * s2;
        gate = (g > arcGate) ? 1.0f : 0.12f;
      }

      // Optional sweeping beam.
      if (sweep && rr < 0.82f) {
        const float beam = std::exp(-(ang*ang) * 1.8f) * (1.0f - rr);
        ring += beam * 0.9f;
      }

      const double n = stellar::proc::smoothNoise2D(nSeed, (double)u * 4.0, (double)v * 4.0);
      const float noise = 0.88f + 0.24f * (float)(n - 0.5);

      float a = clamp01(ring * 0.18f) * gate * outer * noise;
      if (a <= 0.01f) continue;

      FColor c = lerp(col, col2, clamp01(ring * 0.65f));
      c.a = a;

      // Subtle outline near the outer edge.
      if (distPx > -1.0f && distPx < 0.9f) {
        c = lerp(c, outline, 0.35f);
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



static void writeHudPixel(SpriteImage& img, int x, int y, float a) {
  if (a <= 0.0f) return;
  a = clamp01(a);
  const std::size_t out = ((std::size_t)y * (std::size_t)img.w + (std::size_t)x) * 4;
  img.rgba[out + 0] = 255;
  img.rgba[out + 1] = 255;
  img.rgba[out + 2] = 255;
  img.rgba[out + 3] = toByte(a);
}

static inline float wrapPi(float a) {
  // Wrap to [-pi, pi]
  constexpr float kTwoPi = 6.28318530718f;
  constexpr float kPi = 3.14159265359f;
  a = std::fmod(a + kPi, kTwoPi);
  if (a < 0.0f) a += kTwoPi;
  return a - kPi;
}

static SpriteImage genHudReticle(core::u64 seed, int size) {
  // White alpha-mask reticle. Tint it in the HUD via vertex colors.
  // The seed subtly changes notch layout so different ships/factions can feel distinct.
  core::SplitMix64 rng(seed);

  SpriteImage img;
  img.w = size;
  img.h = size;
  img.rgba.resize((std::size_t)size * (std::size_t)size * 4, 0);

  const float invR = 2.0f / (float)size;
  const float aa = invR * 1.35f;

  const float rOuter = 0.80f;
  const float rInner = 0.66f;

  const float crossLen = 0.42f;
  const float crossW = 0.020f + 0.008f * (float)rng.nextDouble();

  const float gapW = 0.17f + 0.06f * (float)rng.nextDouble();

  const core::u64 nSeed = core::hashCombine(seed, core::fnv1a64("hud_ret_noise"));

  for (int y = 0; y < size; ++y) {
    for (int x = 0; x < size; ++x) {
      const float u = ((float)x + 0.5f) / (float)size * 2.0f - 1.0f;
      const float v = ((float)y + 0.5f) / (float)size * 2.0f - 1.0f;
      const float rr = std::sqrt(u*u + v*v);
      if (rr > 1.06f) continue;

      const float ang = std::atan2(v, u);

      // Ring between rInner..rOuter with soft AA.
      const float outerGate = smoothstep(rOuter + aa, rOuter - aa, rr);
      const float innerGate = smoothstep(rInner - aa, rInner + aa, rr);
      float ring = outerGate * innerGate;

      // Gaps at the 4 cardinal directions.
      float minD = std::fabs(wrapPi(ang - 0.0f));
      minD = std::min(minD, std::fabs(wrapPi(ang - 1.57079632679f)));
      minD = std::min(minD, std::fabs(wrapPi(ang - 3.14159265359f)));
      minD = std::min(minD, std::fabs(wrapPi(ang + 1.57079632679f)));
      const float gap = smoothstep(gapW * 0.55f, gapW, minD);
      ring *= gap;

      // Crosshair lines (kept inside the ring).
      float cross = 0.0f;
      if (std::fabs(u) < crossW && std::fabs(v) < crossLen) {
        cross = std::max(cross, 1.0f - smoothstep(crossW, crossW + aa, std::fabs(u)));
      }
      if (std::fabs(v) < crossW && std::fabs(u) < crossLen) {
        cross = std::max(cross, 1.0f - smoothstep(crossW, crossW + aa, std::fabs(v)));
      }

      // Center dot.
      const float dot = 1.0f - smoothstep(0.045f, 0.065f, rr);

      // Gentle scanline/noise variation (keeps it from being perfectly sterile).
      const double n = stellar::proc::smoothNoise2D(nSeed, (double)u * 5.0, (double)v * 5.0);
      const float noise = 0.92f + 0.16f * (float)(n - 0.5);

      float a = std::max(ring * 0.85f, std::max(cross * 0.95f, dot * 0.90f));
      a *= noise;

      if (a > 0.01f) writeHudPixel(img, x, y, a);
    }
  }

  return img;
}

static SpriteImage genHudLead(core::u64 seed, int size) {
  // Lead/aim indicator: ring + dot (white alpha-mask).
  core::SplitMix64 rng(seed);

  SpriteImage img;
  img.w = size;
  img.h = size;
  img.rgba.resize((std::size_t)size * (std::size_t)size * 4, 0);

  const float invR = 2.0f / (float)size;
  const float aa = invR * 1.35f;

  const float rOuter = 0.70f;
  const float rInner = 0.54f;
  const float dotR0 = 0.06f + 0.02f * (float)rng.nextDouble();

  const float notchW = 0.11f + 0.06f * (float)rng.nextDouble();

  for (int y = 0; y < size; ++y) {
    for (int x = 0; x < size; ++x) {
      const float u = ((float)x + 0.5f) / (float)size * 2.0f - 1.0f;
      const float v = ((float)y + 0.5f) / (float)size * 2.0f - 1.0f;
      const float rr = std::sqrt(u*u + v*v);
      if (rr > 1.06f) continue;

      const float ang = std::atan2(v, u);

      const float outerGate = smoothstep(rOuter + aa, rOuter - aa, rr);
      const float innerGate = smoothstep(rInner - aa, rInner + aa, rr);
      float ring = outerGate * innerGate;

      // Small notches at 45-degree diagonals.
      float minD = std::fabs(wrapPi(ang - 0.78539816339f));
      minD = std::min(minD, std::fabs(wrapPi(ang - 2.35619449019f)));
      minD = std::min(minD, std::fabs(wrapPi(ang + 2.35619449019f)));
      minD = std::min(minD, std::fabs(wrapPi(ang + 0.78539816339f)));
      const float notch = smoothstep(notchW * 0.55f, notchW, minD);
      ring *= notch;

      const float dot = 1.0f - smoothstep(dotR0, dotR0 + 0.03f, rr);

      float a = std::max(ring * 0.90f, dot * 0.95f);

      if (a > 0.01f) writeHudPixel(img, x, y, a);
    }
  }

  return img;
}

static SpriteImage genHudVelocity(core::u64 seed, int size) {
  // Flight-path marker: circle with 3 ticks (120deg), white alpha-mask.
  core::SplitMix64 rng(seed);

  SpriteImage img;
  img.w = size;
  img.h = size;
  img.rgba.resize((std::size_t)size * (std::size_t)size * 4, 0);

  const float invR = 2.0f / (float)size;
  const float aa = invR * 1.35f;

  const float rOuter = 0.74f;
  const float rInner = 0.60f;

  const float tickAngW = 0.16f;
  const float tickRad0 = 0.40f + 0.06f * (float)rng.nextDouble();
  const float tickRad1 = tickRad0 + 0.22f;

  for (int y = 0; y < size; ++y) {
    for (int x = 0; x < size; ++x) {
      const float u = ((float)x + 0.5f) / (float)size * 2.0f - 1.0f;
      const float v = ((float)y + 0.5f) / (float)size * 2.0f - 1.0f;
      const float rr = std::sqrt(u*u + v*v);
      if (rr > 1.06f) continue;

      const float ang = std::atan2(v, u);

      const float outerGate = smoothstep(rOuter + aa, rOuter - aa, rr);
      const float innerGate = smoothstep(rInner - aa, rInner + aa, rr);
      float ring = outerGate * innerGate;

      // Three ticks, evenly spaced.
      float tick = 0.0f;
      for (int k = 0; k < 3; ++k) {
        const float a0 = (float)k * 2.09439510239f; // 2*pi/3
        const float dA = std::fabs(wrapPi(ang - a0));
        if (dA < tickAngW && rr > tickRad0 && rr < tickRad1) {
          tick = std::max(tick, 1.0f - smoothstep(tickAngW * 0.6f, tickAngW, dA));
        }
      }

      // Small center dot.
      const float dot = 1.0f - smoothstep(0.045f, 0.075f, rr);

      float a = std::max(ring * 0.75f, std::max(tick * 0.95f, dot * 0.65f));
      if (a > 0.01f) writeHudPixel(img, x, y, a);
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
    case SpriteKind::Cargo:
      return genCargo(seed, size);
    case SpriteKind::Asteroid:
      return genAsteroid(seed, size);
    case SpriteKind::Signal:
      return genSignal(seed, size);
    case SpriteKind::HudReticle:
      return genHudReticle(seed, size);
    case SpriteKind::HudLead:
      return genHudLead(seed, size);
    case SpriteKind::HudVelocity:
      return genHudVelocity(seed, size);
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

// ---------------------------- SpriteAtlas ----------------------------

core::u64 SpriteAtlas::atlasKey(SpriteKind kind, core::u64 seed) const {
  core::u64 k = core::fnv1a64("sprite_atlas");
  k = core::hashCombine(k, (core::u64)(core::u8)kind);
  k = core::hashCombine(k, seed);
  k = core::hashCombine(k, (core::u64)cellSizePx_);
  k = core::hashCombine(k, (core::u64)paddingPx_);
  return k;
}

void SpriteAtlas::clear() {
  entries_.clear();
  freeCells_.clear();
  tick_ = 0;
  // Keep texture allocation around; it will be overwritten as new cells are uploaded.
  for (int i = 0; i < capacity_; ++i) freeCells_.push_back(i);
}

void SpriteAtlas::init(int atlasSizePx, int cellSizePx, int paddingPx, bool nearestFilter) {
  atlasSizePx_ = std::clamp(atlasSizePx, 128, 4096);
  cellSizePx_ = std::clamp(cellSizePx, 8, 256);
  paddingPx_ = std::clamp(paddingPx, 0, 16);

  stridePx_ = cellSizePx_ + paddingPx_ * 2;
  cols_ = atlasSizePx_ / stridePx_;
  rows_ = atlasSizePx_ / stridePx_;
  capacity_ = std::max(0, cols_ * rows_);

  // Clamp max entries to capacity (but allow 0 to mean "no caching" behavior).
  if (maxEntries_ > (std::size_t)capacity_) {
    maxEntries_ = (std::size_t)capacity_;
  }

  entries_.clear();
  freeCells_.clear();
  freeCells_.reserve((std::size_t)capacity_);
  for (int i = 0; i < capacity_; ++i) freeCells_.push_back(i);
  tick_ = 0;

  // Clear the atlas texture to transparent.
  std::vector<std::uint8_t> zeros((std::size_t)atlasSizePx_ * (std::size_t)atlasSizePx_ * 4u, 0);
  tex_.createRGBA(atlasSizePx_, atlasSizePx_, zeros.data(),
                  /*generateMips=*/false,
                  /*nearestFilter=*/nearestFilter,
                  /*clampToEdge=*/true);

  inited_ = true;
}

int SpriteAtlas::allocCell() {
  if (!freeCells_.empty()) {
    int idx = freeCells_.back();
    freeCells_.pop_back();
    return idx;
  }

  // No free cells: evict LRU.
  if (entries_.empty()) return 0;

  auto oldest = entries_.begin();
  for (auto it = entries_.begin(); it != entries_.end(); ++it) {
    if (it->second.lastUseTick < oldest->second.lastUseTick) oldest = it;
  }
  int idx = oldest->second.cellIndex;
  entries_.erase(oldest);
  return idx;
}

static void blitWithPadding(std::vector<std::uint8_t>& dst, int dstW, int dstH,
                            const SpriteImage& src, int pad) {
  // dst is expected to be (src.w + 2*pad) square-ish.
  const int stride = dstW;

  // Copy center.
  for (int y = 0; y < src.h; ++y) {
    for (int x = 0; x < src.w; ++x) {
      const std::size_t s = ((std::size_t)y * (std::size_t)src.w + (std::size_t)x) * 4u;
      const int dx = pad + x;
      const int dy = pad + y;
      const std::size_t d = ((std::size_t)dy * (std::size_t)stride + (std::size_t)dx) * 4u;
      dst[d + 0] = src.rgba[s + 0];
      dst[d + 1] = src.rgba[s + 1];
      dst[d + 2] = src.rgba[s + 2];
      dst[d + 3] = src.rgba[s + 3];
    }
  }

  if (pad <= 0) return;

  // Left/right pad by extending edge texels.
  for (int y = 0; y < src.h; ++y) {
    const int row = pad + y;
    const int leftX = pad;
    const int rightX = pad + src.w - 1;

    for (int p = 1; p <= pad; ++p) {
      const int xL = leftX - p;
      const int xR = rightX + p;
      const std::size_t srcL = ((std::size_t)row * (std::size_t)stride + (std::size_t)leftX) * 4u;
      const std::size_t srcR = ((std::size_t)row * (std::size_t)stride + (std::size_t)rightX) * 4u;
      const std::size_t dstL = ((std::size_t)row * (std::size_t)stride + (std::size_t)xL) * 4u;
      const std::size_t dstR = ((std::size_t)row * (std::size_t)stride + (std::size_t)xR) * 4u;
      for (int c = 0; c < 4; ++c) {
        dst[dstL + (std::size_t)c] = dst[srcL + (std::size_t)c];
        dst[dstR + (std::size_t)c] = dst[srcR + (std::size_t)c];
      }
    }
  }

  // Top/bottom pad by extending edge rows (after left/right fill so corners propagate).
  const int topY = pad;
  const int botY = pad + src.h - 1;
  for (int p = 1; p <= pad; ++p) {
    const int yT = topY - p;
    const int yB = botY + p;
    for (int x = 0; x < dstW; ++x) {
      const std::size_t srcT = ((std::size_t)topY * (std::size_t)stride + (std::size_t)x) * 4u;
      const std::size_t srcB = ((std::size_t)botY * (std::size_t)stride + (std::size_t)x) * 4u;
      const std::size_t dstT = ((std::size_t)yT * (std::size_t)stride + (std::size_t)x) * 4u;
      const std::size_t dstB = ((std::size_t)yB * (std::size_t)stride + (std::size_t)x) * 4u;
      for (int c = 0; c < 4; ++c) {
        dst[dstT + (std::size_t)c] = dst[srcT + (std::size_t)c];
        dst[dstB + (std::size_t)c] = dst[srcB + (std::size_t)c];
      }
    }
  }
}

void SpriteAtlas::uploadCell(int cellIndex, SpriteKind kind, core::u64 seed) {
  if (capacity_ <= 0 || stridePx_ <= 0) return;

  const int cx = cellIndex % std::max(1, cols_);
  const int cy = cellIndex / std::max(1, cols_);

  const int px0 = cx * stridePx_;
  const int py0 = cy * stridePx_;

  SpriteImage img = generateSprite(kind, seed, cellSizePx_);
  if (img.w != cellSizePx_ || img.h != cellSizePx_) {
    // generateSprite clamps to [8,256], but keep defensive.
    img.w = cellSizePx_;
    img.h = cellSizePx_;
  }

  std::vector<std::uint8_t> cell((std::size_t)stridePx_ * (std::size_t)stridePx_ * 4u, 0);
  blitWithPadding(cell, stridePx_, stridePx_, img, paddingPx_);

  tex_.updateRGBA(px0, py0, stridePx_, stridePx_, cell.data());
}

void SpriteAtlas::evictIfNeeded() {
  if (maxEntries_ == 0) {
    entries_.clear();
    freeCells_.clear();
    for (int i = 0; i < capacity_; ++i) freeCells_.push_back(i);
    return;
  }
  if (entries_.size() <= maxEntries_) return;

  while (entries_.size() > maxEntries_) {
    auto oldest = entries_.begin();
    for (auto it = entries_.begin(); it != entries_.end(); ++it) {
      if (it->second.lastUseTick < oldest->second.lastUseTick) oldest = it;
    }
    freeCells_.push_back(oldest->second.cellIndex);
    entries_.erase(oldest);
  }
}

SpriteUvRect SpriteAtlas::get(SpriteKind kind, core::u64 seed) {
  if (!inited_) {
    init();
  }

  const core::u64 k = atlasKey(kind, seed);
  ++tick_;

  if (auto it = entries_.find(k); it != entries_.end()) {
    it->second.lastUseTick = tick_;
    return it->second.uv;
  }

  // Evict LRU if needed (based on maxEntries_ as well as physical capacity).
  if ((int)entries_.size() >= capacity_ && capacity_ > 0) {
    // Force a free cell by evicting one entry.
    auto oldest = entries_.begin();
    for (auto it = entries_.begin(); it != entries_.end(); ++it) {
      if (it->second.lastUseTick < oldest->second.lastUseTick) oldest = it;
    }
    freeCells_.push_back(oldest->second.cellIndex);
    entries_.erase(oldest);
  }

  int cellIndex = allocCell();

  // Compute UV rect for the icon area (excluding padding).
  const int cx = cellIndex % std::max(1, cols_);
  const int cy = cellIndex / std::max(1, cols_);
  const int px0 = cx * stridePx_ + paddingPx_;
  const int py0 = cy * stridePx_ + paddingPx_;

  SpriteUvRect uv;
  uv.u0 = (float)px0 / (float)atlasSizePx_;
  uv.v0 = (float)py0 / (float)atlasSizePx_;
  uv.u1 = (float)(px0 + cellSizePx_) / (float)atlasSizePx_;
  uv.v1 = (float)(py0 + cellSizePx_) / (float)atlasSizePx_;

  uploadCell(cellIndex, kind, seed);

  Entry e;
  e.uv = uv;
  e.cellIndex = cellIndex;
  e.lastUseTick = tick_;
  entries_.emplace(k, e);

  evictIfNeeded();
  return uv;
}

} // namespace stellar::render
