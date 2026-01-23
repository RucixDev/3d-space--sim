#include "stellar/sim/ElectronicWarfare.h"

#include "test_harness.h"

#include <cmath>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

int test_electronic_warfare() {
  int failures = 0;

  // --- Half-range sanity: halfRange -> jamming01 == 0.5 for jammerPower=1 ---
  {
    sim::EwJammerParams p{};
    p.halfRangeKm = 1000.0;

    const double snr = sim::computeJammingSnr(/*distKm=*/1000.0, /*jammerPower=*/1.0, p);
    CHECK(approx(snr, 1.0, 1e-12));

    const double j = sim::jamming01FromSnr(snr);
    CHECK(approx(j, 0.5, 1e-12));
  }

  // --- Inverse-square: double range => snr=0.25 => j=0.2 ---
  {
    sim::EwJammerParams p{};
    p.halfRangeKm = 1000.0;

    const double j = sim::computeJamming01(/*distKm=*/2000.0, /*jammerPower=*/1.0, p);
    CHECK(approx(j, 0.2, 1e-12));
  }

  // --- applyJammingToSensorPower: full jam reduces power; ping punches through ---
  {
    sim::EwJammerParams p{};
    p.suppressionGain = 1.0;
    p.pingJammingMult = 0.25;

    const double base = 1.0;
    const double jam = 1.0;

    const double eff = sim::applyJammingToSensorPower(base, jam, /*pingActive=*/false, p);
    // base / (1 + 1*1) = 0.5
    CHECK(approx(eff, 0.5, 1e-12));

    const double effPing = sim::applyJammingToSensorPower(base, jam, /*pingActive=*/true, p);
    // base / (1 + 1*(0.25)) = 0.8
    CHECK(approx(effPing, 0.8, 1e-12));
  }

  // --- Ghosts: jamming==0 yields empty ---
  {
    sim::EwGhostParams gp{};
    const auto g = sim::generateGhostBlips(/*seed=*/123, /*timeSec=*/0.0, /*rangeKm=*/500.0, /*jamming01=*/0.0, gp);
    CHECK(g.empty());
  }

  // --- Ghosts: deterministic for same (seed, time bucket) and bounded in range ---
  {
    sim::EwGhostParams gp{};
    gp.maxBlips = 10;
    gp.reseedHz = 0.20; // 5s buckets

    const double range = 1000.0;

    const auto g1 = sim::generateGhostBlips(/*seed=*/123, /*timeSec=*/10.0, range, /*jamming01=*/0.85, gp);
    CHECK(!g1.empty());
    CHECK((int)g1.size() <= gp.maxBlips);

    for (const auto& b : g1) {
      CHECK(b.strength01 >= 0.0 && b.strength01 <= 1.0);
      CHECK((b.xKm * b.xKm + b.zKm * b.zKm) <= range * range + 1e-6);
      CHECK(b.id != 0);
    }

    const auto g2 = sim::generateGhostBlips(/*seed=*/123, /*timeSec=*/10.0, range, /*jamming01=*/0.85, gp);
    CHECK(g2.size() == g1.size());
    for (std::size_t i = 0; i < g1.size(); ++i) {
      CHECK(approx(g2[i].xKm, g1[i].xKm, 1e-12));
      CHECK(approx(g2[i].zKm, g1[i].zKm, 1e-12));
      CHECK(g2[i].id == g1[i].id);
    }

    // Different reseed bucket: should differ in pattern (not necessarily size).
    const auto g3 = sim::generateGhostBlips(/*seed=*/123, /*timeSec=*/22.0, range, /*jamming01=*/0.85, gp);
    bool different = (g3.size() != g1.size());
    if (!different && !g3.empty()) {
      for (std::size_t i = 0; i < g1.size() && i < g3.size(); ++i) {
        if (!approx(g3[i].xKm, g1[i].xKm, 1e-6)
            || !approx(g3[i].zKm, g1[i].zKm, 1e-6)
            || g3[i].id != g1[i].id) {
          different = true;
          break;
        }
      }
    }
    CHECK(different);
  }

  return failures;
}
