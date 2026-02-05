#include "stellar/sim/RadarPing.h"

#include "test_harness.h"

#include <cmath>

using namespace stellar;

namespace {

static bool approx(double a, double b, double eps = 1e-12) {
  return std::abs(a - b) <= eps;
}

} // namespace

int test_radar_ping() {
  int failures = 0;

  // ---- Sweep fraction shape (out-and-back triangle wave). ----
  {
    sim::RadarPingParams p;
    p.outAndBack = true;

    const double start = 10.0;
    const double dur = 2.0;

    CHECK(approx(sim::pingFrac01(start, start, dur, p), 0.0, 1e-12));
    CHECK(approx(sim::pingFrac01(start + dur * 0.25, start, dur, p), 0.5, 1e-12));
    CHECK(approx(sim::pingFrac01(start + dur * 0.50, start, dur, p), 1.0, 1e-12));
    CHECK(approx(sim::pingFrac01(start + dur * 0.75, start, dur, p), 0.5, 1e-12));
    CHECK(approx(sim::pingFrac01(start + dur * 1.00, start, dur, p), 0.0, 1e-12));
  }

  // ---- Sweep fraction shape (outbound-only). ----
  {
    sim::RadarPingParams p;
    p.outAndBack = false;

    const double start = 0.0;
    const double dur = 4.0;

    CHECK(approx(sim::pingFrac01(start, start, dur, p), 0.0, 1e-12));
    CHECK(approx(sim::pingFrac01(start + dur * 0.25, start, dur, p), 0.25, 1e-12));
    CHECK(approx(sim::pingFrac01(start + dur * 1.00, start, dur, p), 1.0, 1e-12));
  }

  // ---- Ring boost behavior ----
  {
    sim::RadarPingParams p;
    p.outAndBack = false;
    p.ringThicknessFrac = 0.10; // 10% of range
    p.ringFeatherFrac = 1.0;    // peak-like (no plateau)

    const double rangeKm = 1000.0;
    const double frac = 0.50;
    const double radiusKm = frac * rangeKm;

    const double peak = sim::pingRingBoost01(radiusKm, frac, rangeKm, p);
    CHECK(peak > 0.95);

    // Far outside the ring band should be near zero.
    const double off = sim::pingRingBoost01(radiusKm + 300.0, frac, rangeKm, p);
    CHECK(off < 0.05);

    // Symmetry around the ring radius.
    const double a = sim::pingRingBoost01(radiusKm - 30.0, frac, rangeKm, p);
    const double b = sim::pingRingBoost01(radiusKm + 30.0, frac, rangeKm, p);
    CHECK(approx(a, b, 1e-9));
  }

  // ---- Invalid ranges should return 0. ----
  {
    sim::RadarPingParams p;
    CHECK(approx(sim::pingRingBoost01(10.0, 0.5, 0.0, p), 0.0, 1e-12));
    CHECK(approx(sim::pingRingBoost01(10.0, 0.5, -1.0, p), 0.0, 1e-12));
  }

  return failures;
}
