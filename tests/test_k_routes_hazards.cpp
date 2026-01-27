#include "test_harness.h"

#include "stellar/sim/Celestial.h"
#include "stellar/sim/NavRoute.h"

#include <algorithm>
#include <cmath>

using namespace stellar;

static double distLy(const sim::SystemStub& a, const sim::SystemStub& b) {
  const auto d = a.posLy - b.posLy;
  return std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
}

int test_k_routes_hazards() {
  int failures = 0;

  // A small graph with two distinct viable routes A->C.
  //
  // A (0,0) -4ly- B (4,0) -4ly- C (8,0)
  //                 |
  //               4ly
  //                 |
  //               D (4,4) -4ly- E (8,4) -4ly- C (8,0)
  std::vector<sim::SystemStub> systems;
  systems.push_back(sim::SystemStub{sim::SystemId{1}, 0, "A", math::Vec3d{0, 0, 0}});
  systems.push_back(sim::SystemStub{sim::SystemId{2}, 0, "B", math::Vec3d{4, 0, 0}});
  systems.push_back(sim::SystemStub{sim::SystemId{3}, 0, "C", math::Vec3d{8, 0, 0}});
  systems.push_back(sim::SystemStub{sim::SystemId{4}, 0, "D", math::Vec3d{4, 4, 0}});
  systems.push_back(sim::SystemStub{sim::SystemId{5}, 0, "E", math::Vec3d{8, 4, 0}});

  const double maxJumpLy = 4.1;
  const int k = 2;
  const double costPerLy = 1.0;
  const double costPerJump = 1.0;
  const auto routes = sim::plotKRoutesAStarCost(
      systems,
      systems[0].id,
      systems[2].id,
      maxJumpLy,
      costPerJump,
      costPerLy,
      k);

  CHECK(!routes.empty());
  CHECK(routes.size() >= 1);

  // Best route should be A->B->C.
  CHECK(routes[0].path.size() == 3u);
  CHECK(routes[0].path[0] == systems[0].id);
  CHECK(routes[0].path[1] == systems[1].id);
  CHECK(routes[0].path[2] == systems[2].id);

  // If a 2nd route is returned, it should be distinct and valid.
  if (routes.size() >= 2) {
    CHECK(routes[1].path != routes[0].path);

    // Validate jump range constraints for the 2nd route.
    for (size_t i = 1; i < routes[1].path.size(); ++i) {
      const auto idA = routes[1].path[i - 1];
      const auto idB = routes[1].path[i];
      const auto itA = std::find_if(systems.begin(), systems.end(), [&](const sim::SystemStub& s) { return s.id == idA; });
      const auto itB = std::find_if(systems.begin(), systems.end(), [&](const sim::SystemStub& s) { return s.id == idB; });
      CHECK(itA != systems.end());
      CHECK(itB != systems.end());
      CHECK(distLy(*itA, *itB) <= maxJumpLy + 1e-6);
    }
  }

  return failures;
}
