#include "stellar/sim/Occlusion.h"

#include "test_harness.h"

#include <cmath>
#include <vector>

using namespace stellar;

namespace {

static bool approx(double a, double b, double eps = 1e-12) {
  return std::abs(a - b) <= eps;
}

} // namespace

int test_occlusion_field() {
  int failures = 0;

  // ---- Basic composition semantics (independent occluders). ----
  {
    sim::OcclusionFieldKm f;
    std::vector<sim::SphereOccluderKm> o;
    o.push_back(sim::SphereOccluderKm{{0.0, 0.0, 0.0}, 1.0, 0.25});
    o.push_back(sim::SphereOccluderKm{{5.0, 0.0, 0.0}, 1.0, 0.40});
    f.build(o, /*leafSize=*/4);

    CHECK(!f.empty());
    CHECK(f.occluderCount() == 2);

    const double occ1 = f.occlusion01Segment({-2.0, 0.0, 0.0}, {2.0, 0.0, 0.0});
    CHECK(approx(occ1, 0.25, 1e-12));

    const double occ2 = f.occlusion01Segment({-2.0, 0.0, 0.0}, {7.0, 0.0, 0.0});
    const double expected = 1.0 - (1.0 - 0.25) * (1.0 - 0.40);
    CHECK(approx(occ2, expected, 1e-12));

    const double occ0 = f.occlusion01Segment({-2.0, 3.0, 0.0}, {7.0, 3.0, 0.0});
    CHECK(approx(occ0, 0.0, 1e-12));
  }

  // ---- Endpoint-inside skip (avoid always-occluded artifacts). ----
  {
    sim::OcclusionFieldKm f;
    std::vector<sim::SphereOccluderKm> o;
    o.push_back(sim::SphereOccluderKm{{0.0, 0.0, 0.0}, 2.0, 0.75});
    f.build(o);
    const double occ = f.occlusion01Segment({0.0, 0.0, 0.0}, {10.0, 0.0, 0.0});
    CHECK(approx(occ, 0.0, 1e-12));
  }

  // ---- Max refinement cap reduces work (can underestimate occlusion). ----
  {
    sim::OcclusionFieldKm f;
    std::vector<sim::SphereOccluderKm> o;
    o.reserve(128);
    for (int i = 0; i < 100; ++i) {
      o.push_back(sim::SphereOccluderKm{{(double)i, 0.0, 0.0}, 0.6, 0.10});
    }
    f.build(o, /*leafSize=*/6);

    const math::Vec3d a{-1.0, 0.0, 0.0};
    const math::Vec3d b{101.0, 0.0, 0.0};

    const double occSmall = f.occlusion01Segment(a, b, /*maxSphereTests=*/2);
    const double occLarge = f.occlusion01Segment(a, b, /*maxSphereTests=*/256);

    CHECK(occLarge + 1e-12 >= occSmall);
    CHECK(occLarge <= 1.0 + 1e-12);
  }

  return failures;
}
