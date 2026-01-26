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
  s.seed = (stellar::core::u64)id * 1337ULL;
  s.name = "T" + std::to_string((unsigned long long)id);
  s.posLy = posLy;
  s.primaryClass = stellar::sim::StarClass::G;
  s.planetCount = 0;
  s.stationCount = 0;
  s.factionId = 0;
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

  CHECK(stats.reached);
  CHECK(route.size() == 3);
  CHECK(route[0] == 1);
  CHECK(route[2] == 4);
  CHECK(route[1] == expectedMid);

  // Validate hazard summary helpers against the same sampling logic used above.
  const int samplesPerLeg = 5;
  const auto hz = sim::routeHazardStats(nodes, route, seed, timeDays, samplesPerLeg);

  const std::vector<sim::SystemId> expectTop{1, 2, 4};
  const std::vector<sim::SystemId> expectBottom{1, 3, 4};
  const double expectedIntegral = (expectedMid == 2) ? navIntegral(nodes, expectTop, seed, hp, samplesPerLeg)
                                                     : navIntegral(nodes, expectBottom, seed, hp, samplesPerLeg);
  CHECK(std::abs(hz.integralLy - expectedIntegral) < 1e-9);

  const double expectedDist = sim::routeDistanceLy(nodes, route);
  CHECK(expectedDist > 0.0);
  CHECK(std::abs(hz.average01 - (expectedIntegral / expectedDist)) < 1e-9);

  const auto& a = nodes.at(0);
  const auto& b = nodes.at(expectedMid == 2 ? 1 : 2);
  const auto& c = nodes.at(3);
  const double leg0Avg = avgNavDisruption01(seed, a.posLy, b.posLy, hp, samplesPerLeg);
  const double leg1Avg = avgNavDisruption01(seed, b.posLy, c.posLy, hp, samplesPerLeg);
  const double expectedMax = std::max(leg0Avg, leg1Avg);
  CHECK(std::abs(hz.max01 - expectedMax) < 1e-9);
  if (std::abs(leg0Avg - leg1Avg) > 1e-12) {
    CHECK(hz.maxLegIndex == (leg0Avg > leg1Avg ? 0 : 1));
  } else {
    CHECK(hz.maxLegIndex == 0 || hz.maxLegIndex == 1);
  }

  // routeHazardIntegralLy should match the integral reported by routeHazardStats.
  CHECK(std::abs(sim::routeHazardIntegralLy(nodes, route, seed, timeDays, samplesPerLeg) - hz.integralLy) < 1e-9);

  return failures;
}
