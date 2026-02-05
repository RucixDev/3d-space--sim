#include "stellar/sim/NavRoute.h"

#include <cmath>
#include <iostream>
#include <vector>

int test_nav_route_similarity() {
  int fails = 0;

  using namespace stellar;

  auto expectNear = [&](double got, double want, double eps, const char* msg) {
    if (!(std::isfinite(got) && std::isfinite(want)) || std::abs(got - want) > eps) {
      std::cerr << "[test_nav_route_similarity] " << msg << " got=" << got << " want=" << want << "\n";
      ++fails;
    }
  };

  // Identical routes -> 1.
  {
    std::vector<sim::SystemId> a{1, 2, 3, 4};
    std::vector<sim::SystemId> b{1, 2, 3, 4};
    expectNear(sim::routeNodeJaccard01(a, b, true), 1.0, 1e-12, "identical (ignore endpoints)");
    expectNear(sim::routeNodeJaccard01(a, b, false), 1.0, 1e-12, "identical (include endpoints)");
  }

  // Completely disjoint internal nodes (same endpoints) -> 0 when ignoring endpoints.
  {
    std::vector<sim::SystemId> a{1, 2, 3, 9};
    std::vector<sim::SystemId> b{1, 5, 6, 9};
    expectNear(sim::routeNodeJaccard01(a, b, true), 0.0, 1e-12, "disjoint internals");

    // Including endpoints, they share {1,9}, union {1,2,3,5,6,9} => 2/6.
    expectNear(sim::routeNodeJaccard01(a, b, false), 2.0 / 6.0, 1e-12, "disjoint internals (include endpoints)");
  }

  // Partial overlap.
  {
    std::vector<sim::SystemId> a{1, 2, 3, 4};
    std::vector<sim::SystemId> b{1, 5, 3, 4};

    // Ignore endpoints: A={2,3}, B={5,3} => inter=1, union=3 => 1/3.
    expectNear(sim::routeNodeJaccard01(a, b, true), 1.0 / 3.0, 1e-12, "partial overlap (ignore endpoints)");

    // Include endpoints: A={1,2,3,4}, B={1,5,3,4} => inter={1,3,4}=3, union={1,2,3,4,5}=5 => 3/5.
    expectNear(sim::routeNodeJaccard01(a, b, false), 3.0 / 5.0, 1e-12, "partial overlap (include endpoints)");
  }

  // Degenerate routes.
  {
    std::vector<sim::SystemId> a{7};
    std::vector<sim::SystemId> b{7};
    expectNear(sim::routeNodeJaccard01(a, b, true), 1.0, 1e-12, "single node routes");

    std::vector<sim::SystemId> c{};
    std::vector<sim::SystemId> d{};
    expectNear(sim::routeNodeJaccard01(c, d, true), 1.0, 1e-12, "empty routes");

    std::vector<sim::SystemId> e{1};
    std::vector<sim::SystemId> f{2};
    // Both sets empty when ignoring endpoints (size<2) -> treat as identical empty => 1.
    expectNear(sim::routeNodeJaccard01(e, f, true), 1.0, 1e-12, "single node different ids (ignore endpoints)");

    // Including endpoints: inter=0 union=2 => 0.
    expectNear(sim::routeNodeJaccard01(e, f, false), 0.0, 1e-12, "single node different ids (include endpoints)");
  }

  return fails;
}
