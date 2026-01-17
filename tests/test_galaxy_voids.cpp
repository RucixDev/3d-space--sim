#include "stellar/proc/GalaxyVoids.h"

#include "test_harness.h"

#include <iostream>

using namespace stellar;

int test_galaxy_voids() {
  int failures = 0;

  const core::u64 seed = 424242;

  // Deterministic setup that guarantees coverage: chance=1, radius > cell diagonal.
  proc::GalaxyVoidsParams p{};
  p.cellSizeLy = 1000.0;
  p.chancePerCell = 1.0;
  p.radiusLy = 2000.0;
  p.radiusJitter01 = 0.0;
  p.strengthJitter01 = 0.0;
  p.falloffPower = 2.0;

  const math::Vec3d pos{1234.5, -987.6, 42.0};

  const auto a = proc::sampleGalaxyVoids(seed, pos, p);
  const auto b = proc::sampleGalaxyVoids(seed, pos, p);

  CHECK(a.hasVoid);
  CHECK(b.hasVoid);
  CHECK(a.voidId == b.voidId);
  CHECK(a.voidSeed == b.voidSeed);
  CHECK(a.centerLy.x == b.centerLy.x);
  CHECK(a.centerLy.y == b.centerLy.y);
  CHECK(a.centerLy.z == b.centerLy.z);
  CHECK(a.radiusLy == b.radiusLy);
  CHECK(a.void01 == b.void01);

  CHECK(a.radiusLy > 0.0);
  CHECK(a.void01 >= 0.0 && a.void01 <= 1.0);

  // A very large move should land under a different strongest void.
  {
    const auto c = proc::sampleGalaxyVoids(seed, pos + math::Vec3d{25'000.0, 0.0, 0.0}, p);
    CHECK(c.hasVoid);
    CHECK(c.voidId != a.voidId);
  }

  // Disabling voids should always return "no void".
  {
    proc::GalaxyVoidsParams off = p;
    off.chancePerCell = 0.0;

    const auto d = proc::sampleGalaxyVoids(seed, pos, off);
    CHECK(!d.hasVoid);
    CHECK(d.void01 == 0.0);
    CHECK(d.voidId == 0);
  }

  if (failures == 0) {
    std::cout << "[test_galaxy_voids] PASS\n";
  }
  return failures;
}
