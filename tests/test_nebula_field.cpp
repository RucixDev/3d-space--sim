#include "test_harness.h"

// NebulaField is CPU-only and can be tested even in headless builds.
#include "stellar/render/Nebula.h"

#include "stellar/core/Hash.h"

#include <cmath>
#include <cstdint>
#include <vector>

using namespace stellar;

static inline long long q(double v, double s) {
  return (long long)std::llround(v * s);
}

static core::u64 hashPoints(const std::vector<render::PointVertex>& pts) {
  core::u64 h = core::hashCombine(core::fnv1a64("nebula_points"), (core::u64)pts.size());

  // Quantize before hashing so the test is robust to tiny floating differences.
  for (const auto& p : pts) {
    h = core::hashCombine(h, (core::u64)q((double)p.px, 16.0));
    h = core::hashCombine(h, (core::u64)q((double)p.py, 16.0));
    h = core::hashCombine(h, (core::u64)q((double)p.pz, 16.0));

    h = core::hashCombine(h, (core::u64)q((double)p.cr, 4096.0));
    h = core::hashCombine(h, (core::u64)q((double)p.cg, 4096.0));
    h = core::hashCombine(h, (core::u64)q((double)p.cb, 4096.0));

    h = core::hashCombine(h, (core::u64)q((double)p.a, 65535.0));
    h = core::hashCombine(h, (core::u64)q((double)p.size, 16.0));
  }

  return h;
}

static double dist3(double ax, double ay, double az, double bx, double by, double bz) {
  const double dx = ax - bx;
  const double dy = ay - by;
  const double dz = az - bz;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

int test_nebula_field() {
  int failures = 0;

  const core::u64 seed = 0x123456789abcdef0ull;
  const int puffCount = 512;
  const float bandPower = 1.8f;

  render::NebulaField a;
  render::NebulaField b;
  a.regenerate(seed, puffCount, bandPower);
  b.regenerate(seed, puffCount, bandPower);

  render::NebulaField::Settings s{};
  s.innerRadiusU = 9000.0;
  s.outerRadiusU = 22000.0;
  s.parallax = 0.25;
  s.intensity = 1.4f;
  s.opacity = 0.18f;
  s.sizeMinPx = 90.0f;
  s.sizeMaxPx = 320.0f;
  s.turbulence = 0.35f;
  s.turbulenceSpeed = 0.10f;

  const math::Vec3d camPosU{1234.5, -678.9, 420.0};
  const double t = 12.345;

  a.update(camPosU, t, s);
  b.update(camPosU, t, s);

  CHECK(a.points().size() == (std::size_t)puffCount);
  CHECK(b.points().size() == (std::size_t)puffCount);

  const core::u64 ha = hashPoints(a.points());
  const core::u64 hb = hashPoints(b.points());
  CHECK(ha == hb);

  // Changing the seed should change the field.
  render::NebulaField c;
  c.regenerate(seed + 1u, puffCount, bandPower);
  c.update(camPosU, t, s);
  const core::u64 hc = hashPoints(c.points());
  CHECK(hc != ha);

  // Sanity: alpha/size ranges and broad radius bounds.
  const math::Vec3d anchor = camPosU * s.parallax;
  const double rMin = s.innerRadiusU * 0.60;
  const double rMax = s.outerRadiusU * 1.25;

  for (const auto& p : a.points()) {
    CHECK(p.a >= -1e-4f && p.a <= 1.0001f);
    CHECK(p.size >= s.sizeMinPx - 1e-3f && p.size <= s.sizeMaxPx + 1e-3f);

    const double rr = dist3((double)p.px, (double)p.py, (double)p.pz, anchor.x, anchor.y, anchor.z);
    CHECK(rr >= rMin && rr <= rMax);
  }

  return failures;
}
