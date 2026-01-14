#include "stellar/proc/GalaxyRegions.h"

#include "test_harness.h"

#include <iostream>

using namespace stellar;

int test_galaxy_regions() {
  int failures = 0;

  const core::u64 seed = 1337;
  const double cellSize = 800.0;

  const math::Vec3d pos{1234.5, -987.6, 42.0};

  const auto a = proc::sampleGalaxyRegion(seed, pos, cellSize);
  const auto b = proc::sampleGalaxyRegion(seed, pos, cellSize);

  CHECK(a.regionId == b.regionId);
  CHECK(a.regionSeed == b.regionSeed);
  CHECK(a.kind == b.kind);
  CHECK(a.name == b.name);
  CHECK(a.distanceToCenterLy >= 0.0);
  CHECK(a.approxRadiusLy > 0.0);
  CHECK(a.edge01 >= 0.0 && a.edge01 <= 1.0);
  CHECK(proc::galaxyRegionKindName(a.kind) != nullptr);

  // Spatial coherence: a small nudge should usually remain in the same region.
  {
    const auto c = proc::sampleGalaxyRegion(seed, pos + math::Vec3d{25.0, -10.0, 5.0}, cellSize);
    CHECK(c.regionId == a.regionId);
  }

  // A large move should very likely land in a different region.
  {
    const auto d = proc::sampleGalaxyRegion(seed, pos + math::Vec3d{5000.0, 0.0, 0.0}, cellSize);
    CHECK(d.regionId != a.regionId);
  }

  if (failures == 0) {
    std::cout << "[test_galaxy_regions] PASS\n";
  }
  return failures;
}
