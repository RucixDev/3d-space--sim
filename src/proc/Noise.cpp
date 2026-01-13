#include "stellar/proc/Noise.h"

#include "stellar/core/Hash.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace stellar::proc {

static inline double hash01(stellar::core::u64 h) {
  // map to [0,1)
  return (double)((h >> 11) & ((1ull<<53)-1)) / (double)(1ull<<53);
}

double valueNoise2D(core::u64 seed, int x, int y) {
  core::u64 h = seed;
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(x)));
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(y)));
  return hash01(h);
}

double valueNoise3D(core::u64 seed, int x, int y, int z) {
  core::u64 h = seed;
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(x)));
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(y)));
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(z)));
  return hash01(h);
}

static inline double fade(double t) {
  // Smoothstep-like (Perlin fade)
  return t*t*t*(t*(t*6 - 15) + 10);
}

static inline double lerp(double a, double b, double t) { return a + (b-a)*t; }

double smoothNoise2D(core::u64 seed, double x, double y) {
  const int x0 = static_cast<int>(std::floor(x));
  const int y0 = static_cast<int>(std::floor(y));
  const int x1 = x0 + 1;
  const int y1 = y0 + 1;

  const double tx = x - x0;
  const double ty = y - y0;

  const double n00 = valueNoise2D(seed, x0, y0);
  const double n10 = valueNoise2D(seed, x1, y0);
  const double n01 = valueNoise2D(seed, x0, y1);
  const double n11 = valueNoise2D(seed, x1, y1);

  const double u = fade(tx);
  const double v = fade(ty);

  const double a = lerp(n00, n10, u);
  const double b = lerp(n01, n11, u);
  return lerp(a, b, v);
}

double smoothNoise3D(core::u64 seed, double x, double y, double z) {
  const int x0 = static_cast<int>(std::floor(x));
  const int y0 = static_cast<int>(std::floor(y));
  const int z0 = static_cast<int>(std::floor(z));
  const int x1 = x0 + 1;
  const int y1 = y0 + 1;
  const int z1 = z0 + 1;

  const double tx = x - x0;
  const double ty = y - y0;
  const double tz = z - z0;

  const double n000 = valueNoise3D(seed, x0, y0, z0);
  const double n100 = valueNoise3D(seed, x1, y0, z0);
  const double n010 = valueNoise3D(seed, x0, y1, z0);
  const double n110 = valueNoise3D(seed, x1, y1, z0);
  const double n001 = valueNoise3D(seed, x0, y0, z1);
  const double n101 = valueNoise3D(seed, x1, y0, z1);
  const double n011 = valueNoise3D(seed, x0, y1, z1);
  const double n111 = valueNoise3D(seed, x1, y1, z1);

  const double u = fade(tx);
  const double v = fade(ty);
  const double w = fade(tz);

  const double a0 = lerp(n000, n100, u);
  const double b0 = lerp(n010, n110, u);
  const double a1 = lerp(n001, n101, u);
  const double b1 = lerp(n011, n111, u);

  const double c0 = lerp(a0, b0, v);
  const double c1 = lerp(a1, b1, v);
  return lerp(c0, c1, w);
}

double fbm2D(core::u64 seed, double x, double y, int octaves, double lacunarity, double gain) {
  double amp = 0.5;
  double freq = 1.0;
  double sum = 0.0;
  for (int i = 0; i < octaves; ++i) {
    sum += amp * smoothNoise2D(seed + static_cast<core::u64>(i)*1013ull, x * freq, y * freq);
    freq *= lacunarity;
    amp *= gain;
  }
  return sum;
}

double fbm3D(core::u64 seed, double x, double y, double z, int octaves, double lacunarity, double gain) {
  double amp = 0.5;
  double freq = 1.0;
  double sum = 0.0;
  for (int i = 0; i < octaves; ++i) {
    sum += amp * smoothNoise3D(seed + static_cast<core::u64>(i)*1013ull, x * freq, y * freq, z * freq);
    freq *= lacunarity;
    amp *= gain;
  }
  return sum;
}

// -----------------------------------------------------------------------------
// Gradient noise (Perlin-style)

namespace {

struct Grad2 {
  double x;
  double y;
};

struct Grad3 {
  double x;
  double y;
  double z;
};

constexpr double kInvSqrt2 = 0.70710678118654752440;
constexpr double kInvSqrt3 = 0.57735026918962576450;

constexpr std::array<Grad2, 8> kGrad2 = {
    Grad2{ 1.0, 0.0},
    Grad2{-1.0, 0.0},
    Grad2{ 0.0, 1.0},
    Grad2{ 0.0,-1.0},
    Grad2{ kInvSqrt2,  kInvSqrt2},
    Grad2{-kInvSqrt2,  kInvSqrt2},
    Grad2{ kInvSqrt2, -kInvSqrt2},
    Grad2{-kInvSqrt2, -kInvSqrt2},
};

constexpr std::array<Grad3, 12> kGrad3 = {
    // Classic Perlin gradient set (normalized).
    Grad3{ kInvSqrt2,  kInvSqrt2, 0.0},
    Grad3{-kInvSqrt2,  kInvSqrt2, 0.0},
    Grad3{ kInvSqrt2, -kInvSqrt2, 0.0},
    Grad3{-kInvSqrt2, -kInvSqrt2, 0.0},
    Grad3{ kInvSqrt2, 0.0,  kInvSqrt2},
    Grad3{-kInvSqrt2, 0.0,  kInvSqrt2},
    Grad3{ kInvSqrt2, 0.0, -kInvSqrt2},
    Grad3{-kInvSqrt2, 0.0, -kInvSqrt2},
    Grad3{0.0,  kInvSqrt2,  kInvSqrt2},
    Grad3{0.0, -kInvSqrt2,  kInvSqrt2},
    Grad3{0.0,  kInvSqrt2, -kInvSqrt2},
    Grad3{0.0, -kInvSqrt2, -kInvSqrt2},
};

static inline core::u64 hashCorner2(core::u64 seed, int x, int y) {
  core::u64 h = seed;
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(x)));
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(y)));
  return h;
}

static inline core::u64 hashCorner3(core::u64 seed, int x, int y, int z) {
  core::u64 h = seed;
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(x)));
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(y)));
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(z)));
  return h;
}

static inline Grad2 grad2(core::u64 seed, int x, int y) {
  const core::u64 h = hashCorner2(seed, x, y);
  const std::size_t idx = (std::size_t)(h & 7ull);
  return kGrad2[idx];
}

static inline Grad3 grad3(core::u64 seed, int x, int y, int z) {
  const core::u64 h = hashCorner3(seed, x, y, z);
  const std::size_t idx = (std::size_t)((h >> 24) % (core::u64)kGrad3.size());
  return kGrad3[idx];
}

static inline double clamp01(double v) {
  if (v < 0.0) return 0.0;
  if (v > 1.0) return 1.0;
  return v;
}

} // namespace

double perlin2D(core::u64 seed, double x, double y) {
  const int x0 = static_cast<int>(std::floor(x));
  const int y0 = static_cast<int>(std::floor(y));
  const int x1 = x0 + 1;
  const int y1 = y0 + 1;

  const double fx = x - (double)x0;
  const double fy = y - (double)y0;

  const Grad2 g00 = grad2(seed, x0, y0);
  const Grad2 g10 = grad2(seed, x1, y0);
  const Grad2 g01 = grad2(seed, x0, y1);
  const Grad2 g11 = grad2(seed, x1, y1);

  const double dx0 = fx;
  const double dy0 = fy;
  const double dx1 = fx - 1.0;
  const double dy1 = fy - 1.0;

  const double d00 = g00.x * dx0 + g00.y * dy0;
  const double d10 = g10.x * dx1 + g10.y * dy0;
  const double d01 = g01.x * dx0 + g01.y * dy1;
  const double d11 = g11.x * dx1 + g11.y * dy1;

  const double u = fade(fx);
  const double v = fade(fy);

  const double a = lerp(d00, d10, u);
  const double b = lerp(d01, d11, u);
  double n = lerp(a, b, v);

  // Scale so the common case stays inside [-1,1].
  n *= kInvSqrt2;
  return clamp01(0.5 + 0.5 * n);
}

double perlin3D(core::u64 seed, double x, double y, double z) {
  const int x0 = static_cast<int>(std::floor(x));
  const int y0 = static_cast<int>(std::floor(y));
  const int z0 = static_cast<int>(std::floor(z));
  const int x1 = x0 + 1;
  const int y1 = y0 + 1;
  const int z1 = z0 + 1;

  const double fx = x - (double)x0;
  const double fy = y - (double)y0;
  const double fz = z - (double)z0;

  const Grad3 g000 = grad3(seed, x0, y0, z0);
  const Grad3 g100 = grad3(seed, x1, y0, z0);
  const Grad3 g010 = grad3(seed, x0, y1, z0);
  const Grad3 g110 = grad3(seed, x1, y1, z0);
  const Grad3 g001 = grad3(seed, x0, y0, z1);
  const Grad3 g101 = grad3(seed, x1, y0, z1);
  const Grad3 g011 = grad3(seed, x0, y1, z1);
  const Grad3 g111 = grad3(seed, x1, y1, z1);

  const double dx0 = fx;
  const double dy0 = fy;
  const double dz0 = fz;
  const double dx1 = fx - 1.0;
  const double dy1 = fy - 1.0;
  const double dz1 = fz - 1.0;

  const double d000 = g000.x * dx0 + g000.y * dy0 + g000.z * dz0;
  const double d100 = g100.x * dx1 + g100.y * dy0 + g100.z * dz0;
  const double d010 = g010.x * dx0 + g010.y * dy1 + g010.z * dz0;
  const double d110 = g110.x * dx1 + g110.y * dy1 + g110.z * dz0;
  const double d001 = g001.x * dx0 + g001.y * dy0 + g001.z * dz1;
  const double d101 = g101.x * dx1 + g101.y * dy0 + g101.z * dz1;
  const double d011 = g011.x * dx0 + g011.y * dy1 + g011.z * dz1;
  const double d111 = g111.x * dx1 + g111.y * dy1 + g111.z * dz1;

  const double u = fade(fx);
  const double v = fade(fy);
  const double w = fade(fz);

  const double a0 = lerp(d000, d100, u);
  const double b0 = lerp(d010, d110, u);
  const double a1 = lerp(d001, d101, u);
  const double b1 = lerp(d011, d111, u);
  const double c0 = lerp(a0, b0, v);
  const double c1 = lerp(a1, b1, v);
  double n = lerp(c0, c1, w);

  n *= kInvSqrt3;
  return clamp01(0.5 + 0.5 * n);
}

double fbmPerlin2D(core::u64 seed, double x, double y, int octaves, double lacunarity, double gain) {
  double amp = 0.5;
  double freq = 1.0;
  double sum = 0.0;
  for (int i = 0; i < octaves; ++i) {
    sum += amp * perlin2D(seed + static_cast<core::u64>(i) * 1013ull, x * freq, y * freq);
    freq *= lacunarity;
    amp *= gain;
  }
  return sum;
}

double fbmPerlin3D(core::u64 seed, double x, double y, double z, int octaves, double lacunarity, double gain) {
  double amp = 0.5;
  double freq = 1.0;
  double sum = 0.0;
  for (int i = 0; i < octaves; ++i) {
    sum += amp * perlin3D(seed + static_cast<core::u64>(i) * 1013ull, x * freq, y * freq, z * freq);
    freq *= lacunarity;
    amp *= gain;
  }
  return sum;
}

double ridgedFbmPerlin2D(core::u64 seed, double x, double y, int octaves, double lacunarity, double gain) {
  double amp = 0.5;
  double freq = 1.0;
  double sum = 0.0;
  for (int i = 0; i < octaves; ++i) {
    const double n = perlin2D(seed + static_cast<core::u64>(i) * 1013ull, x * freq, y * freq);
    double ridge = 1.0 - std::fabs(2.0 * n - 1.0); // [0,1]
    ridge = ridge * ridge; // sharpen
    sum += amp * ridge;
    freq *= lacunarity;
    amp *= gain;
  }
  return sum;
}

double ridgedFbmPerlin3D(core::u64 seed, double x, double y, double z, int octaves, double lacunarity, double gain) {
  double amp = 0.5;
  double freq = 1.0;
  double sum = 0.0;
  for (int i = 0; i < octaves; ++i) {
    const double n = perlin3D(seed + static_cast<core::u64>(i) * 1013ull, x * freq, y * freq, z * freq);
    double ridge = 1.0 - std::fabs(2.0 * n - 1.0);
    ridge = ridge * ridge;
    sum += amp * ridge;
    freq *= lacunarity;
    amp *= gain;
  }
  return sum;
}

} // namespace stellar::proc
