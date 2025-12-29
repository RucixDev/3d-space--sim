#include "stellar/render/ProceduralLivery.h"

#include "stellar/core/Hash.h"

#include <algorithm>
#include <cmath>
#include <cctype>

namespace stellar::render {

static std::string lowerAscii(std::string s) {
  for (char& c : s) c = (char)std::tolower((unsigned char)c);
  return s;
}

const char* toString(LiveryPattern p) {
  switch (p) {
    case LiveryPattern::Solid:   return "Solid";
    case LiveryPattern::Stripes: return "Stripes";
    case LiveryPattern::Hazard:  return "Hazard";
    case LiveryPattern::Camo:    return "Camo";
    case LiveryPattern::Hex:     return "Hex";
    case LiveryPattern::Digital: return "Digital";
    default: return "Unknown";
  }
}

std::optional<LiveryPattern> patternFromString(const std::string& s) {
  const std::string k = lowerAscii(s);
  if (k == "solid") return LiveryPattern::Solid;
  if (k == "stripes" || k == "stripe") return LiveryPattern::Stripes;
  if (k == "hazard" || k == "caution") return LiveryPattern::Hazard;
  if (k == "camo" || k == "camouflage") return LiveryPattern::Camo;
  if (k == "hex" || k == "honeycomb") return LiveryPattern::Hex;
  if (k == "digital" || k == "digitalcamo") return LiveryPattern::Digital;
  return std::nullopt;
}

// --- Tiny deterministic noise helpers (no dependencies on proc/Noise) ---

static inline core::u32 hash32(core::u32 x) {
  // Wang hash / mix.
  x = (x ^ 61u) ^ (x >> 16);
  x *= 9u;
  x = x ^ (x >> 4);
  x *= 0x27d4eb2du;
  x = x ^ (x >> 15);
  return x;
}

static inline float u32ToUnit(core::u32 x) {
  // 24 bits of precision.
  return (float)(x & 0x00FFFFFFu) / (float)0x01000000u;
}

static inline float lerp(float a, float b, float t) { return a + (b - a) * t; }

static inline float smooth(float t) { return t * t * (3.0f - 2.0f * t); }

static float valueNoise2D(core::u64 seed, float x, float y) {
  const int xi = (int)std::floor(x);
  const int yi = (int)std::floor(y);
  const float xf = x - (float)xi;
  const float yf = y - (float)yi;

  const float u = smooth(xf);
  const float v = smooth(yf);

  auto h = [&](int ix, int iy) -> float {
    core::u64 k = seed;
    k = core::hashCombine(k, (core::u64)(core::u32)ix);
    k = core::hashCombine(k, (core::u64)(core::u32)iy);
    return u32ToUnit(hash32((core::u32)(k ^ (k >> 32))));
  };

  const float a = h(xi + 0, yi + 0);
  const float b = h(xi + 1, yi + 0);
  const float c = h(xi + 0, yi + 1);
  const float d = h(xi + 1, yi + 1);

  const float x0 = lerp(a, b, u);
  const float x1 = lerp(c, d, u);
  return lerp(x0, x1, v);
}

static float fbm2D(core::u64 seed, float x, float y, int octaves, float lacunarity, float gain) {
  float sum = 0.0f;
  float amp = 0.5f;
  float freq = 1.0f;
  for (int i = 0; i < octaves; ++i) {
    sum += amp * valueNoise2D(seed + (core::u64)i * 1013ull, x * freq, y * freq);
    freq *= lacunarity;
    amp *= gain;
  }
  return sum;
}

static inline void clamp01(float& x) { x = std::clamp(x, 0.0f, 1.0f); }

static void mix3(const float a[3], const float b[3], float t, float out[3]) {
  out[0] = lerp(a[0], b[0], t);
  out[1] = lerp(a[1], b[1], t);
  out[2] = lerp(a[2], b[2], t);
}

static void mul3(float c[3], float s) {
  c[0] *= s; c[1] *= s; c[2] *= s;
}

static void add3(float c[3], float s) {
  c[0] += s; c[1] += s; c[2] += s;
}

static void clamp3(float c[3]) {
  clamp01(c[0]);
  clamp01(c[1]);
  clamp01(c[2]);
}

static void putPixel(std::vector<std::uint8_t>& rgba, int w, int x, int y, const float c[3], std::uint8_t a = 255) {
  const std::size_t i = (std::size_t)(y * w + x) * 4;
  rgba[i + 0] = (std::uint8_t)std::clamp((int)std::lround(c[0] * 255.0f), 0, 255);
  rgba[i + 1] = (std::uint8_t)std::clamp((int)std::lround(c[1] * 255.0f), 0, 255);
  rgba[i + 2] = (std::uint8_t)std::clamp((int)std::lround(c[2] * 255.0f), 0, 255);
  rgba[i + 3] = a;
}

// Tiny 7-seg-like digit glyph for decals.
static bool digitMask(int d, float u, float v) {
  // u,v in [0..1]
  const float t = 0.12f; // thickness
  const bool top    = (v > 1.0f - t) && (u > t) && (u < 1.0f - t);
  const bool mid    = (v > 0.5f - t * 0.5f) && (v < 0.5f + t * 0.5f) && (u > t) && (u < 1.0f - t);
  const bool bot    = (v < t) && (u > t) && (u < 1.0f - t);
  const bool tl     = (u < t) && (v > 0.5f) && (v < 1.0f - t);
  const bool bl     = (u < t) && (v < 0.5f) && (v > t);
  const bool tr     = (u > 1.0f - t) && (v > 0.5f) && (v < 1.0f - t);
  const bool br     = (u > 1.0f - t) && (v < 0.5f) && (v > t);

  auto seg = [&](int s) -> bool {
    switch (s) {
      case 0: return top;
      case 1: return tr;
      case 2: return br;
      case 3: return bot;
      case 4: return bl;
      case 5: return tl;
      case 6: return mid;
      default: return false;
    }
  };

  // Segment table (0..6) for digits.
  static const int mask[10] = {
    0b0111111, // 0
    0b0000110, // 1
    0b1011011, // 2
    0b1001111, // 3
    0b1100110, // 4
    0b1101101, // 5
    0b1111101, // 6
    0b0000111, // 7
    0b1111111, // 8
    0b1101111  // 9
  };
  if (d < 0 || d > 9) return false;
  const int m = mask[d];
  bool on = false;
  for (int s = 0; s < 7; ++s) {
    if ((m >> s) & 1) on = on || seg(s);
  }
  return on;
}

static void drawDecal(std::vector<std::uint8_t>& rgba, int w, int h,
                      core::u64 seed, const float ink[3]) {
  // A subtle registration stencil placed in the upper-left quadrant.
  // The UV placement is intentionally simple: it's just a rectangle in texture space.
  const int boxW = std::max(24, w / 4);
  const int boxH = std::max(18, h / 8);
  const int x0 = w / 14;
  const int y0 = h - boxH - h / 14;

  // Deterministic 3-digit id.
  core::u32 r = hash32((core::u32)(seed ^ (seed >> 32)));
  const int d0 = (int)(r % 10u);
  const int d1 = (int)((r / 10u) % 10u);
  const int d2 = (int)((r / 100u) % 10u);

  auto blend = [&](int x, int y, float a) {
    const std::size_t i = (std::size_t)(y * w + x) * 4;
    const float dst[3] = {
      rgba[i + 0] / 255.0f,
      rgba[i + 1] / 255.0f,
      rgba[i + 2] / 255.0f,
    };
    float out[3] = {
      lerp(dst[0], ink[0], a),
      lerp(dst[1], ink[1], a),
      lerp(dst[2], ink[2], a),
    };
    putPixel(rgba, w, x, y, out);
  };

  // Border
  for (int y = 0; y < boxH; ++y) {
    for (int x = 0; x < boxW; ++x) {
      const bool edge = (x == 0 || y == 0 || x == boxW - 1 || y == boxH - 1);
      if (edge) blend(x0 + x, y0 + y, 0.55f);
    }
  }

  // Digits
  const int pad = boxW / 10;
  const int digitW = (boxW - pad * 4) / 3;
  const int digitH = boxH - pad * 2;
  const int digits[3] = {d0, d1, d2};

  for (int k = 0; k < 3; ++k) {
    const int dx0 = x0 + pad + k * (digitW + pad);
    const int dy0 = y0 + pad;
    for (int y = 0; y < digitH; ++y) {
      for (int x = 0; x < digitW; ++x) {
        const float u = (digitW > 1) ? (float)x / (float)(digitW - 1) : 0.0f;
        const float v = (digitH > 1) ? (float)y / (float)(digitH - 1) : 0.0f;
        if (digitMask(digits[k], u, v)) {
          blend(dx0 + x, dy0 + y, 0.70f);
        }
      }
    }
  }
}

// Distance to nearest edge of a unit hex cell (approx), used to draw outlines.
static float hexEdgeDist(float x, float y) {
  // Convert to axial-ish coords; this is a cheap approximation for UI art.
  // Based on projecting to 60-degree axes.
  const float q = (2.0f / 3.0f) * x;
  const float r = (-1.0f / 3.0f) * x + (std::sqrt(3.0f) / 3.0f) * y;
  const float s = -q - r;

  // fractional part around the nearest integer cube coords
  const float rq = std::round(q);
  const float rr = std::round(r);
  const float rs = std::round(s);

  float dq = std::fabs(q - rq);
  float dr = std::fabs(r - rr);
  float ds = std::fabs(s - rs);
  return std::min({dq, dr, ds});
}

LiveryImage generateLiveryTexture(const LiveryDesc& in, int sizePx) {
  LiveryDesc d = in;
  sizePx = std::clamp(sizePx, 32, 2048);
  d.scale = std::clamp(d.scale, 0.25f, 4.0f);
  d.detail = std::clamp(d.detail, 0.0f, 1.0f);
  d.wear = std::clamp(d.wear, 0.0f, 1.0f);
  d.contrast = std::clamp(d.contrast, 0.0f, 1.0f);

  LiveryImage img;
  img.w = sizePx;
  img.h = sizePx;
  img.rgba.resize((std::size_t)img.w * (std::size_t)img.h * 4u, 255u);

  const float rad = d.angleDeg * 3.14159265f / 180.0f;
  const float ca = std::cos(rad);
  const float sa = std::sin(rad);

  const float base[3]    = {d.base[0], d.base[1], d.base[2]};
  const float accent1[3] = {d.accent1[0], d.accent1[1], d.accent1[2]};
  const float accent2[3] = {d.accent2[0], d.accent2[1], d.accent2[2]};

  const core::u64 seed = d.seed;

  for (int y = 0; y < img.h; ++y) {
    for (int x = 0; x < img.w; ++x) {
      // Normalized UV in [-0.5..0.5] for pattern math.
      const float u = ((float)x + 0.5f) / (float)img.w;
      const float v = ((float)y + 0.5f) / (float)img.h;
      const float px = (u - 0.5f) * d.scale * 2.0f;
      const float py = (v - 0.5f) * d.scale * 2.0f;

      float c[3] = {base[0], base[1], base[2]};

      const float n0 = fbm2D(seed ^ 0xC0FFEEull, px * 1.6f, py * 1.6f, 4, 2.15f, 0.55f);
      const float n1 = fbm2D(seed ^ 0xBADC0DEull, px * 3.0f + 17.0f, py * 3.0f - 9.0f, 3, 2.2f, 0.55f);

      if (d.pattern == LiveryPattern::Solid) {
        const float t = (n0 - 0.5f) * 0.18f * d.contrast;
        add3(c, t);
      } else if (d.pattern == LiveryPattern::Stripes) {
        // Rotate space for angled stripes.
        const float rx = px * ca - py * sa;
        const float stripe = std::sin((rx * 3.14159265f) * (2.0f + 10.0f * d.detail));
        const float m = (stripe > 0.0f) ? 1.0f : 0.0f;
        float sc[3];
        mix3(base, accent1, 0.92f, sc);
        mix3(c, sc, m * 0.85f, c);
        // Add small panel variation.
        add3(c, (n1 - 0.5f) * 0.12f);
      } else if (d.pattern == LiveryPattern::Hazard) {
        // Classic caution stripes: accent1 = yellow-ish, base = black-ish.
        const float rx = px * ca - py * sa;
        const float stripe = std::sin((rx * 3.14159265f) * (4.0f + 10.0f * d.detail));
        const float m = (stripe > 0.0f) ? 1.0f : 0.0f;
        float yellow[3];
        mix3(accent1, accent2, 0.25f, yellow);
        mix3(base, yellow, m, c);
        // Darken edges for depth.
        const float edge = std::pow(std::clamp(std::fabs(rx) * 0.7f, 0.0f, 1.0f), 2.0f);
        mul3(c, 0.80f + 0.20f * (1.0f - edge));
      } else if (d.pattern == LiveryPattern::Camo) {
        // Thresholded FBM gives soft blobs.
        const float t = fbm2D(seed ^ 0xA11CEull, px * 2.0f, py * 2.0f, 5, 2.0f, 0.55f);
        const float a = fbm2D(seed ^ 0xDEADC0DEull, px * 1.2f + 11.0f, py * 1.2f - 7.0f, 3, 2.2f, 0.60f);
        float p0 = (t > 0.55f) ? 1.0f : 0.0f;
        float p1 = (a > 0.58f) ? 1.0f : 0.0f;
        float ca0[3];
        float ca1[3];
        mix3(base, accent1, 0.65f, ca0);
        mix3(base, accent2, 0.55f, ca1);
        if (p0 > 0.5f) mix3(c, ca0, 0.95f, c);
        if (p1 > 0.5f) mix3(c, ca1, 0.85f, c);
        // Edge soften
        add3(c, (n1 - 0.5f) * 0.10f);
      } else if (d.pattern == LiveryPattern::Hex) {
        // Hex tiling + outlines.
        const float sx = px * (3.0f + 10.0f * d.detail);
        const float sy = py * (3.0f + 10.0f * d.detail);
        const float dist = hexEdgeDist(sx, sy);
        const float outline = std::clamp((0.15f - dist) * 5.0f, 0.0f, 1.0f);
        const float cell = valueNoise2D(seed ^ 0x1234ull, std::floor(sx), std::floor(sy));
        float cc[3];
        mix3(accent1, accent2, cell, cc);
        mix3(c, cc, 0.65f, c);
        // Outline dark
        mul3(c, 1.0f - 0.55f * outline);
      } else if (d.pattern == LiveryPattern::Digital) {
        // Pixelated camo blocks.
        const float grid = 6.0f + 44.0f * d.detail;
        const float gx = std::floor((u * grid));
        const float gy = std::floor((v * grid));
        const float r0 = valueNoise2D(seed ^ 0xFACEB00Cull, gx, gy);
        float cc[3];
        if (r0 < 0.35f) mix3(base, accent1, 0.85f, cc);
        else if (r0 < 0.65f) mix3(base, accent2, 0.80f, cc);
        else mix3(accent1, accent2, 0.5f, cc);
        mix3(c, cc, 0.88f, c);
        add3(c, (n1 - 0.5f) * 0.08f);
      }

      // Wear pass: scratches + soot. (A subtle darkening with bright edge specks.)
      if (d.wear > 1e-4f) {
        const float scratches = fbm2D(seed ^ 0x5157ull, px * 14.0f, py * 14.0f, 2, 2.0f, 0.5f);
        const float pits = fbm2D(seed ^ 0xD17Dull, px * 6.0f + 3.0f, py * 6.0f - 4.0f, 3, 2.3f, 0.55f);
        const float s = std::clamp((scratches - 0.72f) * 8.0f, 0.0f, 1.0f);
        const float p = std::clamp((pits - 0.55f) * 2.2f, 0.0f, 1.0f);
        const float grime = d.wear * (0.35f * p + 0.65f * (1.0f - n0));

        // Darken overall.
        mul3(c, 1.0f - grime * 0.55f);

        // Scratch highlights.
        add3(c, s * d.wear * 0.12f);
      }

      // Final contrast tweak (simple power curve).
      if (d.contrast > 1e-4f) {
        const float p = 0.85f + 0.65f * (1.0f - d.contrast);
        c[0] = std::pow(std::clamp(c[0], 0.0f, 1.0f), p);
        c[1] = std::pow(std::clamp(c[1], 0.0f, 1.0f), p);
        c[2] = std::pow(std::clamp(c[2], 0.0f, 1.0f), p);
      }

      clamp3(c);
      putPixel(img.rgba, img.w, x, y, c, 255);
    }
  }

  if (d.decal) {
    // ink color derived from accent2 but darker.
    float ink[3] = {accent2[0], accent2[1], accent2[2]};
    mul3(ink, 0.65f);
    // Stable seed salt ("decal") - avoids non-hex literal typos like 0xDECAL.
    const core::u64 decalSeed = core::hashCombine(seed, core::fnv1a64("decal"));
    drawDecal(img.rgba, img.w, img.h, decalSeed, ink);
  }

  return img;
}

} // namespace stellar::render
