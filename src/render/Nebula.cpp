#include "stellar/render/Nebula.h"

#include "stellar/proc/Noise.h"

#include <algorithm>
#include <cmath>

namespace stellar::render {

static constexpr double kTwoPi = 6.283185307179586476925286766559;

static inline void hsvToRgb(float hDeg, float s, float v, float& r, float& g, float& b) {
  // h: [0,360)
  float h = std::fmod(hDeg, 360.0f);
  if (h < 0.0f) h += 360.0f;

  const float c = v * s;
  const float x = c * (1.0f - std::fabs(std::fmod(h / 60.0f, 2.0f) - 1.0f));
  const float m = v - c;

  float rr = 0.0f, gg = 0.0f, bb = 0.0f;
  if (h < 60.0f) {
    rr = c; gg = x; bb = 0.0f;
  } else if (h < 120.0f) {
    rr = x; gg = c; bb = 0.0f;
  } else if (h < 180.0f) {
    rr = 0.0f; gg = c; bb = x;
  } else if (h < 240.0f) {
    rr = 0.0f; gg = x; bb = c;
  } else if (h < 300.0f) {
    rr = x; gg = 0.0f; bb = c;
  } else {
    rr = c; gg = 0.0f; bb = x;
  }

  r = rr + m;
  g = gg + m;
  b = bb + m;
}

void NebulaField::regenerate(core::u64 seed, int puffCount, float bandPower) {
  seed_ = seed;
  bandPower_ = std::max(1.0f, bandPower);

  if (puffCount < 0) puffCount = 0;

  puffs_.clear();
  puffs_.reserve((std::size_t)puffCount);

  core::SplitMix64 rng(seed);

  auto randDir = [&]() -> math::Vec3d {
    const double ang = rng.nextDouble() * kTwoPi;

    // Start with uniform y in [-1,1], then squeeze toward 0 to create a "galactic plane" band.
    const double y0 = rng.nextDouble() * 2.0 - 1.0;
    const double ySign = (y0 >= 0.0) ? 1.0 : -1.0;
    const double y = ySign * std::pow(std::abs(y0), (double)bandPower_);

    const double r = std::sqrt(std::max(0.0, 1.0 - y * y));
    const double x = r * std::cos(ang);
    const double z = r * std::sin(ang);
    return {x, y, z};
  };

  // Seeded palettes in HSV space (nebulae tend to live in blue/purple/teal ranges)
  const core::u64 hueSeed = core::hashCombine(seed, core::fnv1a64("nebula_hue"));

  for (int i = 0; i < puffCount; ++i) {
    Puff p{};
    p.dir = randDir();

    // Depth selects where between inner/outer radius a puff lives.
    // Bias toward deeper (larger) clouds slightly.
    const double d0 = rng.nextDouble();
    p.depth01 = (float)std::pow(d0, 0.75);

    // Color: sample a smooth-ish hue field on the XZ plane so adjacent puffs tend to cluster by hue.
    const double hx = p.dir.x * 3.0 + 12.7;
    const double hz = p.dir.z * 3.0 - 8.3;
    const double hn = proc::fbm2D(hueSeed, hx, hz, 4, 2.0, 0.55);

    const float hue = (float)(200.0 + 145.0 * hn + 25.0 * rng.nextDouble());
    const float sat = (float)(0.40 + 0.55 * rng.nextDouble());
    const float val = (float)(0.50 + 0.45 * rng.nextDouble());
    hsvToRgb(hue, sat, val, p.r, p.g, p.b);

    // Per-puff base alpha and size scalars.
    p.alpha = (float)(0.35 + 0.65 * rng.nextDouble());
    p.size01 = (float)rng.nextDouble();

    // Turbulence / twinkle
    p.twinkleSpeed = (float)(0.15 + 0.85 * rng.nextDouble());
    p.phase = (float)(rng.nextDouble() * kTwoPi);

    puffs_.push_back(p);
  }

  points_.clear();
  points_.reserve(puffs_.size());
}

void NebulaField::update(const math::Vec3d& cameraPosU, double timeSeconds, const Settings& s) {
  points_.clear();
  points_.reserve(puffs_.size());

  const double par = std::clamp(s.parallax, 0.0, 1.0);
  const math::Vec3d anchor = cameraPosU * par;

  const double inner = std::max(0.0, std::min(s.innerRadiusU, s.outerRadiusU));
  const double outer = std::max(s.innerRadiusU, s.outerRadiusU);

  const float sizeMin = std::max(1.0f, s.sizeMinPx);
  const float sizeMax = std::max(sizeMin, s.sizeMaxPx);

  const float opacity = std::clamp(s.opacity, 0.0f, 1.0f);
  const float intensity = std::max(0.0f, s.intensity);

  const float turb = std::clamp(s.turbulence, 0.0f, 1.0f);
  const float turbSpeed = std::max(0.0f, s.turbulenceSpeed);

  const core::u64 noiseSeed = core::hashCombine(seed_, core::fnv1a64("nebula_update"));

  for (const Puff& p : puffs_) {
    const double r = inner + (outer - inner) * (double)p.depth01;

    // Mild per-puff noise to break up uniform shells.
    const double n = proc::fbm2D(noiseSeed, p.dir.x * 4.0 + timeSeconds * 0.025, p.dir.z * 4.0 - timeSeconds * 0.020, 3, 2.1, 0.55);
    const double rr = r * (0.82 + 0.36 * n);

    const math::Vec3d pos = anchor + p.dir * rr;

    // Alpha wobble.
    float wobble = 1.0f;
    if (turb > 1e-4f && turbSpeed > 1e-4f) {
      const float w = std::sin((float)(timeSeconds * (double)turbSpeed) * (float)(kTwoPi * (double)p.twinkleSpeed) + p.phase);
      wobble = 1.0f + turb * 0.45f * w;
    }

    float a = opacity * p.alpha * wobble;
    // Use noise as additional breakup.
    a *= (float)(0.72 + 0.55 * n);
    a = std::clamp(a, 0.0f, 1.0f);

    const float sz = sizeMin + (sizeMax - sizeMin) * std::clamp(p.size01, 0.0f, 1.0f);

    PointVertex v{};
    v.px = (float)pos.x;
    v.py = (float)pos.y;
    v.pz = (float)pos.z;
    v.cr = p.r * intensity;
    v.cg = p.g * intensity;
    v.cb = p.b * intensity;
    v.a = a;
    v.size = sz;

    points_.push_back(v);
  }
}

} // namespace stellar::render
