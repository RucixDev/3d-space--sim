#include "stellar/sim/Universe.h"

#include <iostream>

int test_streaming() {
  int fails = 0;

  const stellar::core::u64 seed = 42;
  stellar::sim::Universe a(seed);
  stellar::sim::Universe b(seed);

  const stellar::math::Vec3d pos{0,0,0};
  const auto la = a.queryNearby(pos, 60.0, 32);
  const auto lb = b.queryNearby(pos, 60.0, 32);

  if (la.size() != lb.size()) {
    std::cerr << "[test_streaming] size mismatch " << la.size() << " vs " << lb.size() << "\n";
    ++fails;
  }

  const std::size_t n = std::min<std::size_t>({la.size(), lb.size(), 8});
  for (std::size_t i = 0; i < n; ++i) {
    if (la[i].id != lb[i].id || la[i].name != lb[i].name) {
      std::cerr << "[test_streaming] mismatch at i=" << i << " id/name\n";
      ++fails;
      break;
    }
  }

  if (!la.empty()) {
    // getSystem without hint should still work because id encodes sector
    const auto& sys = a.getSystem(la.front().id, &la.front());
    const auto& sys2 = a.getSystem(la.front().id); // cached
    if (sys.planets.size() != sys2.planets.size()) {
      std::cerr << "[test_streaming] cached system mismatch\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_streaming] pass\n";
  return fails;
}
