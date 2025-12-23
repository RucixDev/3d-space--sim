#include "stellar/sim/Universe.h"

#include <iostream>

int test_proc() {
  int fails = 0;

  stellar::sim::UniverseConfig cfg;
  cfg.seed = 12345;
  cfg.systemCountHint = 10;

  stellar::sim::Universe a(cfg);
  stellar::sim::Universe b(cfg);

  // Determinism check: first N systems should match.
  for (std::size_t i = 0; i < cfg.systemCountHint; ++i) {
    const auto sa = a.system(i);
    const auto sb = b.system(i);

    if (sa->primary.name != sb->primary.name) {
      std::cerr << "[test_proc] star name mismatch at " << i << ": "
                << sa->primary.name << " vs " << sb->primary.name << "\n";
      ++fails;
    }
    if (sa->planets.size() != sb->planets.size()) {
      std::cerr << "[test_proc] planet count mismatch at " << i << "\n";
      ++fails;
    }
    if (sa->id != sb->id) {
      std::cerr << "[test_proc] id mismatch at " << i << "\n";
      ++fails;
    }
  }

  // Different seed should produce different output (very likely).
  stellar::sim::UniverseConfig cfg2 = cfg;
  cfg2.seed = 54321;
  stellar::sim::Universe c(cfg2);
  const auto a0 = a.system(0);
  const auto c0 = c.system(0);
  if (a0->primary.name == c0->primary.name && a0->id == c0->id) {
    std::cerr << "[test_proc] different seed produced identical first system (unexpected)\n";
    ++fails;
  }

  return fails;
}
