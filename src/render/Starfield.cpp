#include "stellar/render/Starfield.h"

#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>

namespace stellar::render {

static math::Vec3d unitVecOnSphere(core::SplitMix64& rng) {
  // Uniform sampling on a sphere using z in [-1,1] and phi in [0,2pi).
  const double z = rng.range(-1.0, 1.0);
  const double phi = rng.range(0.0, 2.0 * math::kPi);
  const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
  return {r * std::cos(phi), r * std::sin(phi), z};
}

static void pickStarColor(core::SplitMix64& rng, float& r, float& g, float& b) {
  // Very rough "blackbody-ish" palette: interpolate between warm and cool tints.
  const double t = std::clamp(rng.nextUnit(), 0.0, 1.0);

  // Warm (K/M-ish) -> Cool (A/B-ish)
  const float wr = 1.00f, wg = 0.88f, wb = 0.70f;
  const float cr = 0.68f, cg = 0.80f, cb = 1.00f;

  auto lerp = [](float a, float b, float t) { return a + (b - a) * t; };
  const float tt = (float)t;
  r = lerp(wr, cr, tt);
  g = lerp(wg, cg, tt);
  b = lerp(wb, cb, tt);

  // A bit of per-star variance so the field isn't too uniform.
  const float v = (float)rng.range(0.90, 1.10);
  r *= v; g *= v; b *= v;
}

void Starfield::regenerate(core::u64 seed, int starCount) {
  seed_ = seed;
  stars_.clear();
  points_.clear();

  starCount = std::clamp(starCount, 0, 200000);
  stars_.reserve((std::size_t)starCount);

  core::SplitMix64 rng(seed);

  for (int i = 0; i < starCount; ++i) {
    Star s;
    s.dir = unitVecOnSphere(rng);

    pickStarColor(rng, s.r, s.g, s.b);

    // Magnitude-ish distribution: many dim, few bright.
    // baseAlpha around [0.08 .. 1.0]
    const double m = std::pow(rng.nextUnit(), 2.8); // skew toward 0
    s.baseAlpha = (float)std::clamp(0.08 + (1.0 - m) * 0.92, 0.02, 1.0);

    // Slight size correlation with brightness.
    const double baseSize = 1.0 + (1.0 - m) * 2.6;
    s.sizePx = (float)std::clamp(baseSize + rng.range(-0.25, 0.35), 0.75, 4.25);

    // Twinkle: slow for dim, slightly faster for bright.
    s.twinkleSpeed = (float)std::clamp(0.35 + (1.0 - m) * 1.8, 0.25, 2.5);
    s.phase = (float)rng.range(0.0, 2.0 * math::kPi);

    stars_.push_back(s);
  }
}

void Starfield::update(const math::Vec3d& cameraPosU, double timeSeconds) {
  points_.clear();
  points_.reserve(stars_.size());

  const double t = timeSeconds;

  for (const auto& s : stars_) {
    const math::Vec3d posU = cameraPosU + s.dir * radiusU_;

    const float tw = 0.85f + 0.15f * (float)std::sin(t * (double)s.twinkleSpeed + (double)s.phase);
    const float a = std::clamp(s.baseAlpha * tw, 0.0f, 1.0f);

    PointVertex v;
    v.px = (float)posU.x;
    v.py = (float)posU.y;
    v.pz = (float)posU.z;
    v.cr = s.r;
    v.cg = s.g;
    v.cb = s.b;
    v.a = a;
    v.size = s.sizePx;
    points_.push_back(v);
  }
}

} // namespace stellar::render
