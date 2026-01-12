#include "stellar/render/ProceduralPlanet.h"

#include "stellar/core/Hash.h"
#include "stellar/proc/Noise.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace stellar::render {

namespace {

constexpr double kPi = 3.1415926535897932384626433832795;

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
  return {
      lerp(a.r, b.r, t),
      lerp(a.g, b.g, t),
      lerp(a.b, b.b, t),
  };
}

static inline double smoothstep(double e0, double e1, double x) {
  const double t = clamp01((x - e0) / (e1 - e0));
  return t * t * (3.0 - 2.0 * t);
}

static inline double fade(double t) {
  // Perlin fade polynomial
  return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
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

static double fbm3D(core::u64 seed, double x, double y, double z, int octaves, double lacunarity, double gain) {
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

static Color surfaceRocky(core::u64 seed, const Vec3& p, double lat, double lon) {
  (void)lat;
  (void)lon;

  const double e0 = fbm3D(seed ^ 0xA11CE5u, p.x * 2.8, p.y * 2.8, p.z * 2.8, 6, 2.0, 0.5);
  const double e = clamp01((e0 - 0.25) * 1.35);

  // Ridged detail.
  const double r0 = fbm3D(seed ^ 0xBADC0FFEu, p.x * 7.5, p.y * 7.5, p.z * 7.5, 4, 2.1, 0.55);
  const double ridged = 1.0 - std::abs(r0 * 2.0 - 1.0);

  // Micro-variation.
  const double dust = fbm3D(seed ^ 0xC0FFEEu, p.x * 15.0, p.y * 15.0, p.z * 15.0, 3, 2.4, 0.55);

  Color baseLo{0.18, 0.16, 0.15};
  Color baseHi{0.62, 0.56, 0.48};
  Color col = lerp(baseLo, baseHi, e);

  // Darken creases.
  col = lerp(col, mul(col, 0.65), (1.0 - ridged) * 0.35);

  // Light dusting.
  col = mul(col, 0.85 + 0.28 * dust);

  // Occasional bright regolith speckles.
  const double speck = fbm3D(seed ^ 0xFEEDBEEFu, p.x * 26.0, p.y * 26.0, p.z * 26.0, 2, 2.5, 0.6);
  const double speckMask = smoothstep(0.82, 0.96, speck);
  col = lerp(col, Color{0.88, 0.84, 0.78}, speckMask * 0.32);

  return col;
}

static Color surfaceDesert(core::u64 seed, const Vec3& p, double lat, double lon) {
  // Big dunes and warm albedo.
  const double e0 = fbm3D(seed ^ 0xD35E7u, p.x * 2.2, p.y * 2.2, p.z * 2.2, 6, 2.0, 0.5);
  const double e = clamp01((e0 - 0.22) * 1.30);

  Color sandLo{0.62, 0.52, 0.28};
  Color sandHi{0.92, 0.86, 0.50};
  Color col = lerp(sandLo, sandHi, e);

  // Dune bands mostly aligned along latitude with some longitudinal warping.
  const double warp = fbm3D(seed ^ 0x51DEu, p.x * 3.4, p.y * 3.4, p.z * 3.4, 4, 2.1, 0.55);
  const double band = 0.5 + 0.5 * std::sin((lat * 16.0) + (lon * 2.2) + (warp - 0.5) * 3.2);
  col = mul(col, 0.90 + 0.18 * band);

  // Rocky outcrops.
  const double rock = fbm3D(seed ^ 0x0BADC0DEu, p.x * 6.0, p.y * 6.0, p.z * 6.0, 4, 2.1, 0.52);
  const double rockMask = smoothstep(0.70, 0.88, rock);
  col = lerp(col, Color{0.36, 0.30, 0.22}, rockMask * 0.55);

  return col;
}

static Color surfaceOcean(core::u64 seed, const Vec3& p, double lat, double lon) {
  (void)lon;

  // Continents: low-frequency field.
  const double n = fbm3D(seed ^ 0x0CE4u, p.x * 2.1, p.y * 2.1, p.z * 2.1, 6, 2.0, 0.5);
  const double land = smoothstep(0.48, 0.55, n);

  // Land elevation.
  const double e0 = fbm3D(seed ^ 0xE1E7u, p.x * 6.5, p.y * 6.5, p.z * 6.5, 5, 2.1, 0.55);
  const double elev = clamp01((e0 - 0.30) * 1.55);

  // Coastline factor for shallow water tint.
  const double coast = clamp01(1.0 - std::abs(n - 0.52) * 18.0);

  Color waterDeep{0.03, 0.08, 0.22};
  Color waterShallow{0.07, 0.22, 0.45};
  Color water = lerp(waterDeep, waterShallow, coast);

  Color landLo{0.12, 0.30, 0.16};
  Color landHi{0.52, 0.46, 0.26};
  Color landCol = lerp(landLo, landHi, elev);

  Color col = lerp(water, landCol, land);

  // Polar ice caps.
  const double cap = smoothstep(0.62, 0.90, std::abs(std::sin(lat)));
  if (cap > 1.0e-6) {
    const double capNoise = fbm3D(seed ^ 0xC4F5u, p.x * 9.0, p.y * 9.0, p.z * 9.0, 3, 2.2, 0.55);
    const double capMask = cap * (0.65 + 0.35 * capNoise);
    col = lerp(col, Color{0.90, 0.95, 1.00}, capMask);
  }

  // Subtle clouds to break up the flatness.
  const double cloudN = fbm3D(seed ^ 0xC10Du, p.x * 4.0, p.y * 4.0, p.z * 4.0, 5, 2.0, 0.5);
  const double cloud = smoothstep(0.70, 0.92, cloudN) * (1.0 - land * 0.25);
  col = lerp(col, Color{1.0, 1.0, 1.0}, cloud * 0.10);

  return col;
}

static Color surfaceIce(core::u64 seed, const Vec3& p, double lat, double lon) {
  (void)lon;

  const double n = fbm3D(seed ^ 0x1CEu, p.x * 3.3, p.y * 3.3, p.z * 3.3, 6, 2.0, 0.5);

  Color iceLo{0.70, 0.84, 0.94};
  Color iceHi{0.95, 0.99, 1.00};
  Color col = lerp(iceLo, iceHi, clamp01((n - 0.22) * 1.35));

  // Cracks / ridges.
  const double crack = fbm3D(seed ^ 0xC24Cu, p.x * 20.0, p.y * 20.0, p.z * 20.0, 3, 2.4, 0.55);
  const double crackMask = smoothstep(0.76, 0.93, crack);
  col = lerp(col, Color{0.55, 0.70, 0.82}, crackMask * 0.35);

  // Slight polar brightening.
  const double cap = smoothstep(0.55, 0.88, std::abs(std::sin(lat)));
  col = lerp(col, Color{0.95, 0.99, 1.00}, cap * 0.30);

  return col;
}

static Color surfaceGasGiant(core::u64 seed, const Vec3& p, double lat, double lon) {
  (void)lon;

  const bool warm = ((seed >> 5) & 1ull) != 0ull;
  Color a0 = warm ? Color{0.62, 0.46, 0.24} : Color{0.22, 0.42, 0.62};
  Color a1 = warm ? Color{0.88, 0.74, 0.42} : Color{0.52, 0.70, 0.88};
  Color a2 = warm ? Color{0.38, 0.28, 0.20} : Color{0.14, 0.22, 0.38};

  const double warp = fbm3D(seed ^ 0x6A59u, p.x * 2.6, p.y * 2.6, p.z * 2.6, 5, 2.0, 0.55);
  const double band = 0.5 + 0.5 * std::sin((lat * 22.0) + (warp - 0.5) * 4.2);

  Color col = lerp(a0, a1, band);
  col = lerp(col, a2, clamp01((warp - 0.55) * 1.8) * 0.25);

  const double turb = fbm3D(seed ^ 0x7B3Bu, p.x * 7.0, p.y * 7.0, p.z * 7.0, 4, 2.1, 0.55);
  col = mul(col, 0.85 + 0.35 * turb);

  // Storms.
  const double stormN = fbm3D(seed ^ 0x570A5u, p.x * 11.0, p.y * 11.0, p.z * 11.0, 4, 2.2, 0.55);
  const double storm = smoothstep(0.78, 0.92, stormN);
  if (storm > 1.0e-6) {
    const double swirl = 0.5 + 0.5 * std::sin((lat * 8.0) + (turb - 0.5) * 6.0);
    Color sc = warm ? Color{0.95, 0.86, 0.72} : Color{0.88, 0.92, 0.98};
    col = lerp(col, sc, storm * (0.35 + 0.25 * swirl));
  }

  return col;
}

static Color surfaceStar(core::u64 seed, const Vec3& p, double lat, double lon) {
  (void)lat;
  (void)lon;

  // Greyscale granulation field; intended to be tinted by instance color.
  const double coarse = fbm3D(seed ^ 0x57A4u, p.x * 5.0, p.y * 5.0, p.z * 5.0, 5, 2.0, 0.55);
  const double fine = fbm3D(seed ^ 0xF1AEu, p.x * 18.0, p.y * 18.0, p.z * 18.0, 3, 2.4, 0.55);

  double v = 0.60 + 0.70 * coarse + 0.25 * (fine - 0.5);

  // Star spots (darken).
  const double spotN = fbm3D(seed ^ 0x5B07u, p.x * 8.0, p.y * 8.0, p.z * 8.0, 4, 2.1, 0.55);
  const double spot = smoothstep(0.74, 0.92, spotN);
  v *= (1.0 - 0.22 * spot);

  // Clamp to [0,1] since tint can push brightness beyond 1 in HDR.
  v = clamp01(v);
  return {v, v, v};
}

static double surfaceCloudAlpha(core::u64 seed, const Vec3& p, double lat, double lon) {
  // A lightweight, deterministic cloud mask. Intended to be rendered as a slightly
  // larger shell around the planet with alpha blending.
  //
  // We keep it seam-free by sampling noise on the unit sphere (p) rather than on UV.

  // Base coverage and large features.
  const double base = fbm3D(seed ^ 0xC10DF00Du, p.x * 2.6, p.y * 2.6, p.z * 2.6, 6, 2.0, 0.5);

  // Detail wisps.
  const double detail = fbm3D(seed ^ 0xBADC0FFEu, p.x * 9.0, p.y * 9.0, p.z * 9.0, 3, 2.35, 0.55);

  // Mild banding/warping so the field doesn't look like static Perlin.
  const double warp = fbm3D(seed ^ 0xFEEDu, p.x * 1.3, p.y * 1.3, p.z * 1.3, 4, 2.0, 0.55);
  const double band = 0.5 + 0.5 * std::sin((lat * (12.0 + 6.0 * warp)) + (lon * 0.7) + (warp - 0.5) * 1.9);

  // More clouds toward the equator, fewer toward the poles.
  const double equ = std::pow(std::abs(std::cos(lat)), 0.35);

  double v = base + 0.35 * (detail - 0.5) + 0.18 * (band - 0.5);

  // Threshold varies with latitude to create a believable band of cloud activity.
  const double thresh = 0.56 - 0.18 * equ;
  double dens = smoothstep(thresh, thresh + 0.20, v);

  // Carve out holes (breaks up the mask so it isn't a full blanket).
  const double hole = fbm3D(seed ^ 0xA5A5u, p.x * 4.2, p.y * 4.2, p.z * 4.2, 4, 2.0, 0.5);
  dens *= (1.0 - 0.45 * smoothstep(0.52, 0.78, hole));

  // Sharpen slightly.
  dens = std::pow(clamp01(dens), 1.15);
  return clamp01(dens);
}



static double surfaceHeight(SurfaceKind kind, core::u64 seed, const Vec3& p, double lat, double lon) {
  // Height fields are intentionally approximate; they exist to generate tangent-space
  // normal maps that add micro-relief to the procedural albedo.
  //
  // Return value is in [0,1] (roughly), but only relative variation matters.

  switch (kind) {
    case SurfaceKind::Rocky: {
      const double e0 = fbm3D(seed ^ 0xA11CE5u, p.x * 2.8, p.y * 2.8, p.z * 2.8, 6, 2.0, 0.5);
      const double e = clamp01((e0 - 0.25) * 1.35);

      const double r0 = fbm3D(seed ^ 0xBADC0FFEu, p.x * 7.5, p.y * 7.5, p.z * 7.5, 4, 2.1, 0.55);
      const double ridged = 1.0 - std::abs(r0 * 2.0 - 1.0);

      const double dust = fbm3D(seed ^ 0xC0FFEEu, p.x * 15.0, p.y * 15.0, p.z * 15.0, 3, 2.4, 0.55);

      double h = 0.70 * e + 0.30 * ridged;
      h = clamp01(h);

      // Dust slightly fills cavities.
      h = clamp01(h * (0.92 + 0.10 * dust));
      return h;
    }
    case SurfaceKind::Desert: {
      const double e0 = fbm3D(seed ^ 0xD35E7u, p.x * 2.2, p.y * 2.2, p.z * 2.2, 6, 2.0, 0.5);
      const double e = clamp01((e0 - 0.22) * 1.30);

      const double warp = fbm3D(seed ^ 0x51DEu, p.x * 3.4, p.y * 3.4, p.z * 3.4, 4, 2.1, 0.55);
      const double band = 0.5 + 0.5 * std::sin((lat * 16.0) + (lon * 2.2) + (warp - 0.5) * 3.2);

      const double rock = fbm3D(seed ^ 0x0BADC0DEu, p.x * 6.0, p.y * 6.0, p.z * 6.0, 4, 2.1, 0.52);
      const double rockMask = smoothstep(0.70, 0.88, rock);

      double h = 0.78 * e + 0.22 * band;
      h += rockMask * 0.12; // outcrops
      return clamp01(h);
    }
    case SurfaceKind::Ocean: {
      const double n = fbm3D(seed ^ 0x0CE4u, p.x * 2.1, p.y * 2.1, p.z * 2.1, 6, 2.0, 0.5);
      const double land = smoothstep(0.48, 0.55, n);

      const double e0 = fbm3D(seed ^ 0xE1E7u, p.x * 6.5, p.y * 6.5, p.z * 6.5, 5, 2.1, 0.55);
      const double elev = clamp01((e0 - 0.30) * 1.55);

      // Mild ocean waves (very small amplitude).
      const double wave = fbm3D(seed ^ 0xA57EA11Cu, p.x * 18.0, p.y * 18.0, p.z * 18.0, 3, 2.2, 0.55);

      const double hOcean = 0.03 * (wave - 0.5);
      const double hLand = 0.25 + 0.75 * elev;
      const double h = lerp(hOcean, hLand, land);
      return clamp01(h);
    }
    case SurfaceKind::Ice: {
      const double n = fbm3D(seed ^ 0x1CEu, p.x * 3.3, p.y * 3.3, p.z * 3.3, 6, 2.0, 0.5);
      const double crack = fbm3D(seed ^ 0xC24Cu, p.x * 20.0, p.y * 20.0, p.z * 20.0, 3, 2.4, 0.55);
      const double crackMask = smoothstep(0.76, 0.93, crack);

      double h = clamp01((n - 0.22) * 1.35);
      h = clamp01(h + 0.22 * crackMask);
      return h;
    }
    case SurfaceKind::GasGiant: {
      const double warp = fbm3D(seed ^ 0x6A59u, p.x * 2.6, p.y * 2.6, p.z * 2.6, 5, 2.0, 0.55);
      const double band = 0.5 + 0.5 * std::sin((lat * 22.0) + (warp - 0.5) * 4.2);
      const double turb = fbm3D(seed ^ 0x7B3Bu, p.x * 7.0, p.y * 7.0, p.z * 7.0, 4, 2.1, 0.55);

      // Keep very subtle to avoid looking like solid rock.
      double h = 0.5 + 0.10 * (band - 0.5) + 0.08 * (turb - 0.5);
      return clamp01(h);
    }
    case SurfaceKind::Star: {
      const double coarse = fbm3D(seed ^ 0x57A4u, p.x * 5.0, p.y * 5.0, p.z * 5.0, 5, 2.0, 0.55);
      const double fine = fbm3D(seed ^ 0xF1AEu, p.x * 18.0, p.y * 18.0, p.z * 18.0, 3, 2.4, 0.55);
      double h = 0.55 + 0.55 * coarse + 0.20 * (fine - 0.5);
      return clamp01(h);
    }
    case SurfaceKind::Clouds: {
      // Clouds are rendered as an alpha shell; use the density field as a height proxy.
      return surfaceCloudAlpha(seed, p, lat, lon);
    }
    default:
      return 0.5;
  }
}

} // namespace

SurfaceImage generateSurfaceTexture(SurfaceKind kind, core::u64 seed, int widthPx) {
  SurfaceImage img{};
  if (widthPx <= 0) return img;

  img.w = std::max(4, widthPx);
  img.h = std::max(2, img.w / 2);
  img.rgba.assign(static_cast<std::size_t>(img.w * img.h * 4), 255);

  // Deterministic rotation so systems don't all align on the seam.
  const double rot = (double)((seed >> 8) & 0xFFFFull) / 65535.0 * (2.0 * kPi);

  for (int y = 0; y < img.h; ++y) {
    const double v = (double)(y + 0.5) / (double)img.h;
    const double lat = (v - 0.5) * kPi; // -pi/2..pi/2
    const double cosLat = std::cos(lat);
    const double sinLat = std::sin(lat);

    for (int x = 0; x < img.w; ++x) {
      const double u = (double)(x + 0.5) / (double)img.w;
      const double lon = u * (2.0 * kPi) - kPi + rot;

      const double cosLon = std::cos(lon);
      const double sinLon = std::sin(lon);

      // Unit sphere sample point.
      Vec3 p{cosLat * cosLon, sinLat, cosLat * sinLon};

      Color c{};
      double a = 1.0;
      switch (kind) {
        case SurfaceKind::Rocky: c = surfaceRocky(seed, p, lat, lon); break;
        case SurfaceKind::Desert: c = surfaceDesert(seed, p, lat, lon); break;
        case SurfaceKind::Ocean: c = surfaceOcean(seed, p, lat, lon); break;
        case SurfaceKind::Ice: c = surfaceIce(seed, p, lat, lon); break;
        case SurfaceKind::GasGiant: c = surfaceGasGiant(seed, p, lat, lon); break;
        case SurfaceKind::Star: c = surfaceStar(seed, p, lat, lon); break;
        case SurfaceKind::Clouds:
          a = surfaceCloudAlpha(seed, p, lat, lon);
          // Slight blue tint gives a bit of depth when lit by the star.
          c = {0.94, 0.97, 1.00};
          break;
        default: c = {0.6, 0.6, 0.6}; break;
      }

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

  // Keep the same deterministic seam rotation as the albedo so features line up.
  const double rot = (double)((seed >> 8) & 0xFFFFull) / 65535.0 * (2.0 * kPi);

  auto wrap01 = [](double x) {
    x = std::fmod(x, 1.0);
    if (x < 0.0) x += 1.0;
    return x;
  };

  auto clamp01d = [](double x) { return std::clamp(x, 0.0, 1.0); };

  auto sampleH = [&](double u, double v) -> double {
    u = wrap01(u);
    v = clamp01d(v);

    const double lat = (v - 0.5) * kPi;
    const double cosLat = std::cos(lat);
    const double sinLat = std::sin(lat);

    const double lon = u * (2.0 * kPi) - kPi + rot;
    const double cosLon = std::cos(lon);
    const double sinLon = std::sin(lon);

    Vec3 p{cosLat * cosLon, sinLat, cosLat * sinLon};
    return surfaceHeight(kind, seed, p, lat, lon);
  };

  const double du = 1.0 / (double)img.w;
  const double dv = 1.0 / (double)img.h;

  // Per-kind slope scaling (the shader also has a global normal strength multiplier).
  double kindScale = 1.0;
  switch (kind) {
    case SurfaceKind::Rocky: kindScale = 1.35; break;
    case SurfaceKind::Desert: kindScale = 1.05; break;
    case SurfaceKind::Ocean: kindScale = 0.85; break;
    case SurfaceKind::Ice: kindScale = 1.05; break;
    case SurfaceKind::GasGiant: kindScale = 0.20; break;
    case SurfaceKind::Star: kindScale = 0.30; break;
    case SurfaceKind::Clouds: kindScale = 0.55; break;
    default: kindScale = 1.0; break;
  }

  for (int y = 0; y < img.h; ++y) {
    const double v = (double)(y + 0.5) / (double)img.h;
    const double lat = (v - 0.5) * kPi;
    const double cosLat = std::cos(lat);

    // Physical distances on a unit sphere for one texel step.
    const double dLon = du * (2.0 * kPi);
    const double dLat = dv * kPi;

    const double dx = std::max(1.0e-6, std::abs(cosLat) * dLon);
    const double dy = std::max(1.0e-6, dLat);

    for (int x = 0; x < img.w; ++x) {
      const double u = (double)(x + 0.5) / (double)img.w;

      const double hL = sampleH(u - du, v);
      const double hR = sampleH(u + du, v);
      const double hD = sampleH(u, v - dv);
      const double hU = sampleH(u, v + dv);

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
