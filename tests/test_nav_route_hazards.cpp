#include "test_harness.h"

#include "stellar/proc/GalaxyHazards.h"
#include "stellar/sim/NavRoute.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

stellar::sim::SystemStub makeStub(stellar::sim::SystemId id, const stellar::math::Vec3d& posLy) {
  stellar::sim::SystemStub s{};
  s.id = id;
  s.name = "T" + std::to_string((unsigned long long)id);
  s.posLy = posLy;
  s.population = 0;
  s.factionId = 0;
  s.econId = 0;
  s.govId = 0;
  s.securityId = 0;
  s.techLevel = 0;
  s.habitable = false;
  s.radiusLy = 0.0;
  return s;
}

double avgNavDisruption01(stellar::core::u64 seed,
                          const stellar::math::Vec3d& aLy,
                          const stellar::math::Vec3d& bLy,
                          const stellar::proc::GalaxyHazardsParams& hp) {
  const int samples = 5;
  double sum = 0.0;
  for (int i = 0; i < samples; ++i) {
    const double t = (i + 0.5) / (double)samples;
    const stellar::math::Vec3d p = aLy + (bLy - aLy) * t;
    sum += stellar::proc::sampleGalaxyHazards(seed, p, hp).navDisruption01;
  }
  return std::clamp(sum / (double)samples, 0.0, 1.0);
}

double navIntegral(stellar::core::u64 seed,
                   const stellar::math::Vec3d& aLy,
                   const stellar::math::Vec3d& bLy,
                   const stellar::proc::GalaxyHazardsParams& hp) {
  const double d = (bLy - aLy).length();
  return avgNavDisruption01(seed, aLy, bLy, hp) * d;
}

} // namespace

int test_nav_route_hazards() {
  int failures = 0;

  using namespace stellar;

  const core::u64 seed = 0xC0FFEEULL;
  const double timeDays = 123.0;
  const double maxJumpLy = 8.0;

  // Build a symmetric 2-hop geometry where start->mid and mid->goal are exactly maxJumpLy,
  // but the two mids are farther apart than maxJumpLy. This yields two disjoint 2-hop routes.
  const double spanLy = 12.0;
  const double halfSpan = spanLy / 2.0;
  const double dy = std::sqrt(std::max(0.0, maxJumpLy * maxJumpLy - halfSpan * halfSpan));

  proc::GalaxyHazardsParams hp{};
  hp.timeDays = timeDays;

  struct Candidate {
    double tx{0.0};
    double ty{0.0};
    double hazTop{0.0};
    double hazBottom{0.0};
    double diff{0.0};
  };

  Candidate best{};
  best.diff = -1.0;

  // Search a small grid of translations to find a region where the hazard field makes
  // the two routes meaningfully different.
  for (int ix = -10; ix <= 10; ++ix) {
    for (int iy = -10; iy <= 10; ++iy) {
      const double tx = ix * 24.0;
      const double ty = iy * 24.0;

      const math::Vec3d start{tx, ty, 0.0};
      const math::Vec3d goal{tx + spanLy, ty, 0.0};
      const math::Vec3d midTop{tx + halfSpan, ty + dy, 0.0};
      const math::Vec3d midBottom{tx + halfSpan, ty - dy, 0.0};

      const double hazTop = navIntegral(seed, start, midTop, hp) + navIntegral(seed, midTop, goal, hp);
      const double hazBottom =
        navIntegral(seed, start, midBottom, hp) + navIntegral(seed, midBottom, goal, hp);
      const double diff = std::abs(hazTop - hazBottom);

      if (diff > best.diff) {
        best = Candidate{tx, ty, hazTop, hazBottom, diff};
      }
    }
  }

  // Keep the threshold loose to avoid flakiness across platforms.
  CHECK(best.diff > 1e-6);

  const math::Vec3d start{best.tx, best.ty, 0.0};
  const math::Vec3d goal{best.tx + spanLy, best.ty, 0.0};
  const math::Vec3d midTop{best.tx + halfSpan, best.ty + dy, 0.0};
  const math::Vec3d midBottom{best.tx + halfSpan, best.ty - dy, 0.0};

  std::vector<sim::SystemStub> nodes;
  nodes.reserve(4);
  nodes.push_back(makeStub(1, start));
  nodes.push_back(makeStub(2, midTop));
  nodes.push_back(makeStub(3, midBottom));
  nodes.push_back(makeStub(4, goal));

  const sim::SystemId expectedMid = (best.hazTop < best.hazBottom) ? 2 : 3;

  sim::RoutePlanStats stats{};
  const auto route = sim::plotRouteAStarCostHazards(nodes,
                                                    /*startId=*/1,
                                                    /*goalId=*/4,
                                                    maxJumpLy,
                                                    /*costPerJump=*/0.0,
                                                    /*costPerLy=*/0.0,
                                                    /*hazardWeightPerLy=*/5.0,
                                                    seed,
                                                    timeDays,
                                                    &stats);

  CHECK(stats.found);
  CHECK(route.size() == 3);
  CHECK(route[0] == 1);
  CHECK(route[2] == 4);
  CHECK(route[1] == expectedMid);

  return failures;
}
