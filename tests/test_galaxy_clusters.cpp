#include "stellar/proc/GalaxyClusters.h"

#include "test_harness.h"

#include <iostream>

using namespace stellar;

int test_galaxy_clusters() {
  int failures = 0;

  const core::u64 seed = 424242;

  // A deterministic setup that guarantees coverage: chance=1, radius > cell diagonal.
  proc::GalaxyClustersParams p{};
  p.cellSizeLy = 1000.0;
  p.chancePerCell = 1.0;
  p.radiusLy = 2000.0;
  p.radiusJitter01 = 0.0;
  p.strengthJitter01 = 0.0;
  p.falloffPower = 2.0;

  const math::Vec3d pos{1234.5, -987.6, 42.0};

  const auto a = proc::sampleGalaxyClusters(seed, pos, p);
  const auto b = proc::sampleGalaxyClusters(seed, pos, p);

  CHECK(a.hasCluster);
  CHECK(b.hasCluster);
  CHECK(a.clusterId == b.clusterId);
  CHECK(a.clusterSeed == b.clusterSeed);
  CHECK(a.centerLy.x == b.centerLy.x);
  CHECK(a.centerLy.y == b.centerLy.y);
  CHECK(a.centerLy.z == b.centerLy.z);
  CHECK(a.radiusLy == b.radiusLy);
  CHECK(a.cluster01 == b.cluster01);

  CHECK(a.radiusLy > 0.0);
  CHECK(a.cluster01 >= 0.0 && a.cluster01 <= 1.0);

  // A very large move should land under a different strongest cluster.
  {
    const auto c = proc::sampleGalaxyClusters(seed, pos + math::Vec3d{25'000.0, 0.0, 0.0}, p);
    CHECK(c.hasCluster);
    CHECK(c.clusterId != a.clusterId);
  }

  // Disabling clusters should always return "no cluster".
  {
    proc::GalaxyClustersParams off = p;
    off.chancePerCell = 0.0;

    const auto d = proc::sampleGalaxyClusters(seed, pos, off);
    CHECK(!d.hasCluster);
    CHECK(d.cluster01 == 0.0);
    CHECK(d.clusterId == 0);
  }

  if (failures == 0) {
    std::cout << "[test_galaxy_clusters] PASS\n";
  }
  return failures;
}
