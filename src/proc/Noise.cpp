#include "stellar/proc/Noise.h"

#include "stellar/core/Hash.h"

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

} // namespace stellar::proc
