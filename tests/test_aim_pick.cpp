#include "test_harness.h"

#include "stellar/sim/AimPick.h"

#include <vector>

int test_aim_pick() {
  int failures = 0;

  using stellar::math::Vec3d;
  using stellar::sim::AimPickHit;
  using stellar::sim::AimPickSphere;

  const Vec3d origin{0, 0, 0};
  const Vec3d dir{0, 0, 1};

  // Two hits along the +Z ray, plus one miss to the side.
  AimPickSphere spheres[3];
  spheres[0] = AimPickSphere{0, Vec3d{0, 0, 10}, 1.0, -1.0};
  spheres[1] = AimPickSphere{1, Vec3d{0, 0, 5}, 1.0, -1.0};
  spheres[2] = AimPickSphere{2, Vec3d{2, 0, 5}, 1.0, -1.0};

  std::vector<AimPickHit> hits;
  stellar::sim::aimPickSphereRayHitsKm(hits, origin, dir, 100.0, spheres, 3);

  CHECK_EQ(hits.size(), (std::size_t)2);
  if (hits.size() >= 2) {
    CHECK_EQ(hits[0].tag, (std::size_t)1);
    CHECK_EQ(hits[1].tag, (std::size_t)0);
    CHECK(hits[0].tEnterKm <= hits[1].tEnterKm);
  }

  // Pad expands the selection radius (the side sphere becomes a hit).
  hits.clear();
  stellar::sim::aimPickSphereRayHitsKm(hits, origin, dir, 100.0, spheres, 3, /*padKm=*/1.2);
  CHECK_EQ(hits.size(), (std::size_t)3);
  if (hits.size() == 3) {
    CHECK_EQ(hits[0].tag, (std::size_t)1);
    CHECK_EQ(hits[1].tag, (std::size_t)2);
    CHECK_EQ(hits[2].tag, (std::size_t)0);
  }

  // Aim cone filter rejects off-axis targets.
  AimPickSphere cone[2];
  cone[0] = AimPickSphere{0, Vec3d{0, 0, 10}, 1.0, 0.99};
  cone[1] = AimPickSphere{1, Vec3d{10, 0, 10}, 1.0, 0.99}; // ~45deg off-axis

  hits.clear();
  stellar::sim::aimPickSphereRayHitsKm(hits, origin, dir, 100.0, cone, 2);
  CHECK_EQ(hits.size(), (std::size_t)1);
  if (!hits.empty()) {
    CHECK_EQ(hits[0].tag, (std::size_t)0);
  }

  // Nearest helper should match the first hit.
  stellar::sim::AimPickHit nearest;
  const bool ok = stellar::sim::aimPickSphereRayNearestKm(nearest, origin, dir, 100.0, spheres, 3);
  CHECK(ok);
  if (ok) {
    CHECK_EQ(nearest.tag, (std::size_t)1);
  }

  return failures;
}
