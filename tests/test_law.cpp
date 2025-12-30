#include "stellar/sim/Law.h"

#include <cmath>
#include <iostream>

static bool nearly(double a, double b, double eps = 1e-12) {
  return std::abs(a - b) <= eps;
}

int test_law() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;

  const core::u64 seed = 1337ull;
  const core::u32 fA = 7;
  const core::u32 fB = 9;

  // Determinism: same inputs => bitwise-stable doubles should match exactly in practice,
  // but we compare with an epsilon to be safe across platforms.
  {
    const LawProfile a = lawProfile(seed, fA);
    const LawProfile b = lawProfile(seed, fA);

    if (!nearly(a.scanStrictness, b.scanStrictness)) {
      std::cerr << "[test_law] scanStrictness not deterministic\n";
      ++fails;
    }
    if (!nearly(a.corruption, b.corruption)) {
      std::cerr << "[test_law] corruption not deterministic\n";
      ++fails;
    }
    if (!nearly(a.fineRate, b.fineRate) || !nearly(a.fineBaseCr, b.fineBaseCr)) {
      std::cerr << "[test_law] fine schedule not deterministic\n";
      ++fails;
    }
  }

  // Basic bounds and monotonicity.
  {
    const LawProfile a = lawProfile(seed, fA);

    if (a.scanStrictness < 0.50 || a.scanStrictness > 2.00) {
      std::cerr << "[test_law] scanStrictness out of expected bounds\n";
      ++fails;
    }
    if (a.corruption < 0.0 || a.corruption > 1.0) {
      std::cerr << "[test_law] corruption out of expected bounds\n";
      ++fails;
    }

    const double v0 = 0.0;
    const double v1 = 500.0;
    const double v2 = 20000.0;

    const double f0 = a.contrabandFineCr(v0);
    const double f1 = a.contrabandFineCr(v1);
    const double f2 = a.contrabandFineCr(v2);

    if (!(f0 >= 0.0 && f1 >= f0 && f2 >= f1)) {
      std::cerr << "[test_law] fine not monotonic\n";
      ++fails;
    }

    const double b0 = a.contrabandBribeCr(v0);
    const double b1 = a.contrabandBribeCr(v1);
    const double b2 = a.contrabandBribeCr(v2);

    if (!(b0 >= 0.0 && b1 >= b0 && b2 >= b1)) {
      std::cerr << "[test_law] bribe not monotonic\n";
      ++fails;
    }

    const double r1 = a.contrabandRepPenalty(v1);
    const double r2 = a.contrabandRepPenalty(v2);
    if (!(r1 <= 0.0 && r2 <= r1 + 1e-9)) {
      std::cerr << "[test_law] rep penalty should be negative and become more negative for larger values\n";
      ++fails;
    }

    const double e1 = a.contrabandEvadeRepPenalty(v1);
    if (!(e1 <= r1 + 1e-9)) {
      std::cerr << "[test_law] evade rep penalty should be >= severity of comply penalty\n";
      ++fails;
    }
  }

  // Different factions should usually produce a different profile.
  // We do not require all fields to differ (some collisions are possible), but at least one should.
  {
    const LawProfile a = lawProfile(seed, fA);
    const LawProfile b = lawProfile(seed, fB);

    const bool anyDifferent =
        !nearly(a.scanStrictness, b.scanStrictness) ||
        !nearly(a.corruption, b.corruption) ||
        !nearly(a.fineRate, b.fineRate) ||
        !nearly(a.fineBaseCr, b.fineBaseCr);

    if (!anyDifferent) {
      std::cerr << "[test_law] expected different factions to produce different law profiles\n";
      ++fails;
    }
  }

  return fails;
}
