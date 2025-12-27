#include "stellar/sim/Universe.h"
#include "stellar/core/Log.h"

#include <iostream>

int test_universe_cache() {
  int fails = 0;

  // Keep unit-test output clean even if a fallback path logs a warning.
  stellar::core::setLogLevel(stellar::core::LogLevel::Off);

  const stellar::core::u64 seed = 42;

  // ---------------------------------------------------------------------------
  // Cache caps: clamped to at least 1 (so getSystem/sector/stationEconomy can
  // safely return references).
  // ---------------------------------------------------------------------------
  {
    stellar::sim::Universe u(seed);
    u.setCacheCaps(0, 0, 0);

    const auto s = u.cacheStats();
    if (s.sectors.capacity < 1 || s.systems.capacity < 1 || s.stationEconomies.capacity < 1) {
      std::cerr << "[test_universe_cache] cache caps not clamped ("
                << "sector=" << s.sectors.capacity
                << " system=" << s.systems.capacity
                << " station=" << s.stationEconomies.capacity << ")\n";
      ++fails;
    }
  }

  // ---------------------------------------------------------------------------
  // Sector cache stats: deterministic hit/miss/eviction behavior for single-sector queries.
  // ---------------------------------------------------------------------------
  {
    stellar::sim::Universe u(seed);
    u.setCacheCaps(2, 2, 2);
    u.resetCacheStats();

    const stellar::math::Vec3d p0{5.0, 5.0, 5.0};   // sector (0,0,0)
    const stellar::math::Vec3d p1{15.0, 5.0, 5.0};  // sector (1,0,0)
    const stellar::math::Vec3d p2{25.0, 5.0, 5.0};  // sector (2,0,0)

    // Each of these radius values should touch exactly one sector.
    (void)u.queryNearby(p0, 4.0, 8);
    (void)u.queryNearby(p0, 4.0, 8); // hit
    (void)u.queryNearby(p1, 4.0, 8);
    (void)u.queryNearby(p2, 4.0, 8); // should evict one

    const auto s = u.cacheStats();

    if (s.sectors.capacity != 2 || s.sectors.size != 2) {
      std::cerr << "[test_universe_cache] sector cache size/cap mismatch "
                << "size=" << s.sectors.size << " cap=" << s.sectors.capacity << "\n";
      ++fails;
    }

    if (s.sectors.hits != 1 || s.sectors.misses != 3 || s.sectors.puts != 3) {
      std::cerr << "[test_universe_cache] unexpected sector cache stats "
                << "(hits=" << s.sectors.hits
                << " misses=" << s.sectors.misses
                << " puts=" << s.sectors.puts << ")\n";
      ++fails;
    }

    if (s.sectors.evictions < 1) {
      std::cerr << "[test_universe_cache] expected sector cache eviction(s)\n";
      ++fails;
    }
  }

  // ---------------------------------------------------------------------------
  // System cache stats: verify hit/miss + eviction when capacity is exceeded.
  // ---------------------------------------------------------------------------
  {
    // Pick real ids from a separate Universe so our cache-stat assertions remain crisp.
    stellar::sim::Universe picker(seed);

    auto stubs = picker.queryNearby(stellar::math::Vec3d{0.0, 0.0, 0.0}, 25.0, 32);
    if (stubs.size() < 3) {
      stubs = picker.queryNearby(stellar::math::Vec3d{0.0, 0.0, 0.0}, 40.0, 64);
    }

    if (stubs.size() < 3) {
      std::cerr << "[test_universe_cache] not enough nearby stubs to run system-cache test\n";
      ++fails;
    } else {
      const stellar::sim::SystemId a = stubs[0].id;
      const stellar::sim::SystemId b = stubs[1].id;
      const stellar::sim::SystemId c = stubs[2].id;

      stellar::sim::Universe u(seed);
      u.setCacheCaps(64, 2, 2); // keep sector cache roomy to avoid noisy evictions
      u.resetCacheStats();

      (void)u.getSystem(a); // miss + put
      (void)u.getSystem(a); // hit
      (void)u.getSystem(b); // miss + put
      (void)u.getSystem(c); // miss + put (evict one)

      const auto s = u.cacheStats();

      if (s.systems.capacity != 2 || s.systems.size != 2) {
        std::cerr << "[test_universe_cache] system cache size/cap mismatch "
                  << "size=" << s.systems.size << " cap=" << s.systems.capacity << "\n";
        ++fails;
      }

      if (s.systems.hits != 1 || s.systems.misses != 3 || s.systems.puts != 3) {
        std::cerr << "[test_universe_cache] unexpected system cache stats "
                  << "(hits=" << s.systems.hits
                  << " misses=" << s.systems.misses
                  << " puts=" << s.systems.puts << ")\n";
        ++fails;
      }

      if (s.systems.evictions < 1) {
        std::cerr << "[test_universe_cache] expected system cache eviction(s)\n";
        ++fails;
      }
    }
  }

  if (fails == 0) std::cout << "[test_universe_cache] pass\n";
  return fails;
}
