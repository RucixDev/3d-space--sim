#include "stellar/sim/FieldRepair.h"

#include "test_harness.h"

#include <cmath>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

int test_field_repair() {
  int failures = 0;

  // --- Quote: full repair to cap when stock is sufficient ---
  {
    sim::FieldRepairParams p{};
    p.maxRepairFrac = 0.75;
    p.hullPerRealSec = 2.5;
    p.metalsPerHull = 0.30;
    p.machineryPerHull = 0.22;

    sim::FieldRepairInventory inv{};
    inv.metals = 30.0;
    inv.machinery = 30.0;

    const double hullCur = 20.0;
    const double hullMax = 100.0;

    const auto q = sim::quoteFieldRepairToCap(inv, hullCur, hullMax, p);
    CHECK(q.ok);
    CHECK(approx(q.capHull, 75.0));
    CHECK(approx(q.hullMissingToCap, 55.0));
    CHECK(approx(q.hullToRepair, 55.0));
    CHECK(q.limitedByCap);
    CHECK(!q.limitedByStock);
    CHECK(approx(q.metalsNeeded, 55.0 * p.metalsPerHull));
    CHECK(approx(q.machineryNeeded, 55.0 * p.machineryPerHull));
    CHECK(approx(q.timeSec, 55.0 / p.hullPerRealSec));
  }

  // --- Quote: limited by stock ---
  {
    sim::FieldRepairParams p{};
    p.maxRepairFrac = 0.75;
    p.hullPerRealSec = 5.0;
    p.metalsPerHull = 0.30;
    p.machineryPerHull = 0.22;

    sim::FieldRepairInventory inv{};
    inv.metals = 3.0;    // limiting
    inv.machinery = 99.0;

    const auto q = sim::quoteFieldRepairToCap(inv, /*hullCurrent=*/10.0, /*hullMax=*/100.0, p);
    CHECK(q.ok);
    // cap hull = 75, missing=65. metals allow 3/0.3=10 hull.
    CHECK(approx(q.hullToRepair, 10.0));
    CHECK(q.limitedByStock);
    CHECK(approx(q.metalsNeeded, 3.0));
    CHECK(approx(q.machineryNeeded, 10.0 * p.machineryPerHull));
    CHECK(approx(q.timeSec, 10.0 / p.hullPerRealSec));
  }

  // --- Step: consumes materials and repairs, clamps to cap ---
  {
    sim::FieldRepairParams p{};
    p.maxRepairFrac = 0.75;
    p.hullPerRealSec = 10.0;
    p.metalsPerHull = 0.50;
    p.machineryPerHull = 0.25;
    p.heatPerHull = 2.0;

    sim::FieldRepairInventory inv{};
    inv.metals = 10.0;
    inv.machinery = 10.0;

    double hullCur = 70.0;
    const double hullMax = 100.0;

    // missing to cap = 5 hull; at 10 hull/s, dt=1 would want 10 but should clamp to 5.
    const auto r = sim::stepFieldRepair(/*dtRealSec=*/1.0, hullMax, inv, hullCur, p);
    CHECK(r.progressed);
    CHECK(r.done);
    CHECK(r.reachedCap);
    CHECK(approx(r.hullRepaired, 5.0));
    CHECK(approx(hullCur, 75.0));

    CHECK(approx(r.metalsTaken, 5.0 * p.metalsPerHull));
    CHECK(approx(r.machineryTaken, 5.0 * p.machineryPerHull));
    CHECK(approx(r.heatAdded, 5.0 * p.heatPerHull));
  }

  // --- Step: out of supplies ---
  {
    sim::FieldRepairParams p{};
    p.maxRepairFrac = 1.0;
    p.hullPerRealSec = 10.0;
    p.metalsPerHull = 1.0;
    p.machineryPerHull = 1.0;

    sim::FieldRepairInventory inv{};
    inv.metals = 0.0;
    inv.machinery = 5.0;

    double hullCur = 50.0;
    const auto r = sim::stepFieldRepair(/*dtRealSec=*/1.0, /*hullMax=*/100.0, inv, hullCur, p);
    CHECK(r.done);
    CHECK(r.outOfSupplies);
    CHECK(!r.progressed);
    CHECK(approx(hullCur, 50.0));
  }

  return failures;
}
