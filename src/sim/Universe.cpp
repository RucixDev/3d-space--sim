#include "stellar/sim/Universe.h"

#include "stellar/core/Log.h"

#include <algorithm>
#include <sstream>

namespace stellar::sim {

Universe::Universe(UniverseConfig cfg)
  : m_cfg(cfg)
  , m_galaxy(proc::GalaxyGenConfig{
      .seed = m_cfg.seed,
      .systemCount = m_cfg.systemCountHint,
      .radiusLy = m_cfg.radiusLy,
      .thicknessLy = m_cfg.thicknessLy,
    })
  , m_sysgen(proc::SystemGenConfig{ .galaxySeed = m_cfg.seed })
{
  // Avoid pathological configs.
  if (m_cfg.systemCacheSize == 0) {
    m_cfg.systemCacheSize = 1;
  }
}

void Universe::reset() {
  m_lru.clear();
  m_cache.clear();
}

Universe::SystemHandle Universe::system(SystemIndex index) {
  // Hit cache
  if (auto it = m_cache.find(index); it != m_cache.end()) {
    // Move to front (MRU)
    m_lru.splice(m_lru.begin(), m_lru, it->second.lruIt);
    it->second.lruIt = m_lru.begin();
    return it->second.sys;
  }

  // Miss: generate deterministically.
  const auto stub = m_galaxy.stubAt(index);
  auto sys = std::make_shared<StarSystem>(m_sysgen.generate(stub));

  // Insert into LRU
  m_lru.push_front(index);
  m_cache.emplace(index, CacheEntry{sys, m_lru.begin()});

  // Evict
  while (m_cache.size() > m_cfg.systemCacheSize && !m_lru.empty()) {
    const auto victim = m_lru.back();
    m_lru.pop_back();
    m_cache.erase(victim);
  }

  return sys;
}

void Universe::prewarm(std::size_t count) {
  std::ostringstream oss;
  oss << "Prewarming universe cache with " << count << " systems...";
  STELLAR_LOG_INFO(oss.str());

  for (std::size_t i = 0; i < count; ++i) {
    (void)system(i);
  }
}

} // namespace stellar::sim
