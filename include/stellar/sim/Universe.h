#pragma once

#include "stellar/core/Types.h"
#include "stellar/proc/GalaxyGenerator.h"
#include "stellar/proc/SystemGenerator.h"
#include "stellar/sim/Celestial.h"

#include <cstddef>
#include <list>
#include <memory>
#include <unordered_map>

namespace stellar::sim {

struct UniverseConfig {
  stellar::core::u64 seed = 1;

  // Optional hint for tooling/sandbox apps that want to sample the universe.
  // The Universe itself is effectively infinite: any system index can be requested.
  std::size_t systemCountHint = 1000;

  // LRU cache size for streamed star systems.
  // (Systems are deterministic, so they can be regenerated if evicted.)
  std::size_t systemCacheSize = 512;

  // Simple disc galaxy shape in light-years
  double radiusLy = 50000.0;
  double thicknessLy = 2000.0;
};

class Universe {
public:
  explicit Universe(UniverseConfig cfg);

  using SystemIndex = std::size_t;
  using SystemHandle = std::shared_ptr<const StarSystem>;

  // Clears internal caches. Universe data is generated on-demand.
  void reset();

  // Returns a deterministic star system for the given index.
  // Note: the returned shared_ptr keeps the system alive even if it is evicted from the cache.
  SystemHandle system(SystemIndex index);

  const UniverseConfig& config() const { return m_cfg; }

  std::size_t cachedSystemCount() const { return m_cache.size(); }

  // Convenience: warm up the cache with the first N systems.
  void prewarm(std::size_t count);

private:
  struct CacheEntry {
    SystemHandle sys;
    std::list<SystemIndex>::iterator lruIt;
  };

  UniverseConfig m_cfg{};

  // Stateless procedural generators (deterministic per index).
  proc::GalaxyGenerator m_galaxy;
  proc::SystemGenerator m_sysgen;

  // LRU cache
  std::list<SystemIndex> m_lru; // front = most recently used
  std::unordered_map<SystemIndex, CacheEntry> m_cache;
};

} // namespace stellar::sim
