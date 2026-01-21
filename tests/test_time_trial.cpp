#include "test_harness.h"

#include "stellar/sim/TimeTrial.h"

#include "stellar/math/Quat.h"

#include <cmath>

using namespace stellar;

int test_time_trial() {
  int failures = 0;

  // --- Gate crossing ---
  {
    sim::TimeTrialGate g{};
    g.posKm = {0, 0, 0};
    g.normal = {0, 0, 1};
    g.radiusKm = 1.0;

    CHECK(sim::timeTrialGatePassed(g, {0, 0, -2}, {0, 0, +2}, {0, 0, +1}));
    CHECK(!sim::timeTrialGatePassed(g, {0, 0, +2}, {0, 0, -2}, {0, 0, -1}));
    CHECK(!sim::timeTrialGatePassed(g, {2.5, 0, -2}, {2.5, 0, +2}, {0, 0, +1}));
    CHECK(!sim::timeTrialGatePassed(g, {0, 0, -2}, {0, 0, +2}, {0, 0, -1}));
  }

  // --- Deterministic generator invariants ---
  {
    sim::TimeTrialCourseParams p{};
    p.gateCount = 16;
    p.gateRadiusKm = 2500.0;
    p.courseRadiusKm = 60000.0;
    p.courseHeightKm = 12000.0;
    p.jitterKm = 8000.0;
    p.loops = 2;
    p.closedLoop = false;

    const auto c0 = sim::generateTimeTrialCourseStationSlalomKm({0, 0, 0}, math::Quatd::identity(), 12345u, p);
    const auto c1 = sim::generateTimeTrialCourseStationSlalomKm({0, 0, 0}, math::Quatd::identity(), 12345u, p);
    const auto c2 = sim::generateTimeTrialCourseStationSlalomKm({0, 0, 0}, math::Quatd::identity(), 12346u, p);

    CHECK(c0.gates.size() == (std::size_t)p.gateCount);
    CHECK(c1.gates.size() == (std::size_t)p.gateCount);
    CHECK(c0.key == c1.key);
    CHECK(c0.key != c2.key);

    for (std::size_t i = 0; i < c0.gates.size(); ++i) {
      const auto& a = c0.gates[i];
      const auto& b = c1.gates[i];
      CHECK((a.posKm - b.posKm).length() < 1e-6);
      CHECK(std::abs(a.normal.length() - 1.0) < 1e-6);
      CHECK(std::abs(a.radiusKm - p.gateRadiusKm) < 1e-12);
    }
  }

  return failures;
}
