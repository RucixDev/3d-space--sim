#include "stellar/proc/GalaxyHazards.h"

#include "test_harness.h"

#include <cmath>
#include <iostream>

using namespace stellar;

int test_galaxy_hazards() {
  int failures = 0;

  const core::u64 seed = 1337;

  proc::GalaxyHazardsParams params{};
  params.nebulaScaleLy = 700.0;
  params.stormScaleLy = 210.0;
  params.regionCellSizeLy = 900.0;
  params.regionInfluence = 0.35;
  params.timeDays = 0.0;
  params.driftLyPerDay = 0.65;

  const math::Vec3d pos{1234.5, -987.6, 42.0};

  const auto a = proc::sampleGalaxyHazards(seed, pos, params);
  const auto b = proc::sampleGalaxyHazards(seed, pos, params);

  // Determinism.
  CHECK(a.kind == b.kind);
  CHECK(a.nebula01 == b.nebula01);
  CHECK(a.storm01 == b.storm01);
  CHECK(a.hazard01 == b.hazard01);

  // Ranges.
  CHECK(a.nebula01 >= 0.0 && a.nebula01 <= 1.0);
  CHECK(a.storm01 >= 0.0 && a.storm01 <= 1.0);
  CHECK(a.hazard01 >= 0.0 && a.hazard01 <= 1.0);
  CHECK(a.sensorOcclusion01 == a.nebula01);
  CHECK(a.navDisruption01 == a.storm01);
  CHECK(proc::galaxyHazardKindName(a.kind) != nullptr);

  // Spatial coherence: a small nudge should not completely decorrelate.
  {
    const auto c = proc::sampleGalaxyHazards(seed, pos + math::Vec3d{25.0, -10.0, 5.0}, params);
    CHECK(std::abs(c.hazard01 - a.hazard01) < 0.60);
  }

  // Time drift: should generally change with time.
  {
    proc::GalaxyHazardsParams p2 = params;
    p2.timeDays = 250.0;
    const auto t = proc::sampleGalaxyHazards(seed, pos, p2);
    CHECK(std::abs(t.hazard01 - a.hazard01) > 1e-9);
  }

  // Segment sampling helper.
  {
    const math::Vec3d bpos = pos + math::Vec3d{120.0, 50.0, -25.0};
    const double avg = proc::sampleGalaxyHazardAvgOnSegment(seed, pos, bpos, params, 5);
    const double avg2 = proc::sampleGalaxyHazardAvgOnSegment(seed, pos, bpos, params, 5);
    CHECK(avg == avg2);
    CHECK(avg >= 0.0 && avg <= 1.0);
  }

  if (failures == 0) {
    std::cout << "[test_galaxy_hazards] PASS\n";
  }
  return failures;
}
