#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Economy.h"
#include "stellar/proc/GalaxyGenerator.h"
#include "stellar/sim/Faction.h"
#include "stellar/sim/System.h"
#include "stellar/sim/SaveGame.h"

#include <list>
#include <optional>
#include <unordered_map>

namespace stellar::sim {

// A streaming universe:
// - galaxy is generated in fixed-size sectors on-demand
// - star systems are generated on-demand from system stubs
class Universe {
public:
  explicit Universe(core::u64 seed, proc::GalaxyParams params = {});

  core::u64 seed() const { return seed_; }
  const proc::GalaxyParams& galaxyParams() const { return galaxyParams_; }

  const std::vector<Faction>& factions() const { return factions_; }

  // Return stubs for systems within `radiusLy` of `posLy`.
  // Deterministic order by distance, then by id.
  std::vector<SystemStub> queryNearby(const math::Vec3d& posLy,
                                      double radiusLy,
                                      std::size_t maxResults = 256);

  // Find closest system stub within `maxRadiusLy` (returns nullopt if none found).
  std::optional<SystemStub> findClosestSystem(const math::Vec3d& posLy, double maxRadiusLy);

  // Generate (or fetch cached) full system data.
  const StarSystem& getSystem(SystemId id, const SystemStub* hintStub = nullptr);

  // Economy state access for a station. Will create deterministic initial state if missing,
  // then advance it to `timeDays`.
  econ::StationEconomyState& stationEconomy(const Station& station, double timeDays);

  // Persist/restore known station economy states (only what you have in cache / visited).
  std::vector<StationEconomyOverride> exportStationOverrides() const;
  void importStationOverrides(const std::vector<StationEconomyOverride>& overrides);

  void setCacheCaps(std::size_t sectorCap, std::size_t systemCap, std::size_t stationCap);

  // Lightweight cache diagnostics (useful for tests/tools and for tuning cache caps).
  struct CacheStats {
    struct One {
      std::size_t capacity{0};
      std::size_t size{0};
      std::size_t hits{0};
      std::size_t misses{0};
      std::size_t puts{0};
      std::size_t evictions{0};
    };

    One sectors{};
    One systems{};
    One stationEconomies{};
  };

  CacheStats cacheStats() const;
  void resetCacheStats();

private:
  template <class Key, class Value, class Hash = std::hash<Key>>
  class LruCache {
  public:
    explicit LruCache(std::size_t cap = 128) : cap_(cap ? cap : 1) {}

    void setCapacity(std::size_t cap) {
      cap_ = cap ? cap : 1;
      evictIfNeeded();
    }
    std::size_t capacity() const { return cap_; }
    std::size_t size() const { return map_.size(); }

    struct Stats {
      std::size_t capacity{0};
      std::size_t size{0};
      std::size_t hits{0};
      std::size_t misses{0};
      std::size_t puts{0};
      std::size_t evictions{0};
    };

    Stats stats() const {
      Stats s;
      s.capacity = cap_;
      s.size = map_.size();
      s.hits = hits_;
      s.misses = misses_;
      s.puts = puts_;
      s.evictions = evictions_;
      return s;
    }

    void resetStats() { hits_ = misses_ = puts_ = evictions_ = 0; }

    Value* get(const Key& key) {
      auto it = map_.find(key);
      if (it == map_.end()) {
        ++misses_;
        return nullptr;
      }
      ++hits_;
      touch(it);
      return &it->second.value;
    }

    const Value* get(const Key& key) const {
      auto it = map_.find(key);
      if (it == map_.end()) {
        ++misses_;
        return nullptr;
      }
      ++hits_;
      // const cache doesn't update LRU; ok for now
      return &it->second.value;
    }

    Value& put(const Key& key, Value value) {
      ++puts_;
      auto it = map_.find(key);
      if (it != map_.end()) {
        it->second.value = std::move(value);
        touch(it);
        return it->second.value;
      }

      order_.push_front(key);
      Entry e;
      e.value = std::move(value);
      e.it = order_.begin();
      auto itIns = map_.emplace(key, std::move(e)).first;

      evictIfNeeded();
      return itIns->second.value;
    }

    std::unordered_map<Key, Value, Hash> snapshot() const {
      std::unordered_map<Key, Value, Hash> r;
      r.reserve(map_.size());
      for (const auto& kv : map_) {
        r.emplace(kv.first, kv.second.value);
      }
      return r;
    }

  private:
    struct Entry {
      Value value{};
      typename std::list<Key>::iterator it{};
    };

    void touch(typename std::unordered_map<Key, Entry, Hash>::iterator it) {
      order_.erase(it->second.it);
      order_.push_front(it->first);
      it->second.it = order_.begin();
    }

    void evictIfNeeded() {
      while (map_.size() > cap_) {
        const Key& key = order_.back();
        map_.erase(key);
        order_.pop_back();
        ++evictions_;
      }
    }

    std::size_t cap_{128};
    std::list<Key> order_{};
    std::unordered_map<Key, Entry, Hash> map_{};

    // These are intentionally mutable so callers can query cache stats on const caches
    // (e.g. via const Universe) without changing functional behavior.
    mutable std::size_t hits_{0};
    mutable std::size_t misses_{0};
    mutable std::size_t puts_{0};
    mutable std::size_t evictions_{0};
  };

  // Return a cached sector (generated on-demand). Returning by reference avoids
  // copying potentially large stub vectors during frequent queries.
  const proc::Sector& sector(const proc::SectorCoord& coord);

  core::u64 seed_{0};
  proc::GalaxyParams galaxyParams_{};
  proc::GalaxyGenerator galaxyGen_{0, {}};

  std::vector<Faction> factions_{};

  LruCache<proc::SectorCoord, proc::Sector, proc::SectorCoordHash> sectorCache_{256};
  LruCache<SystemId, StarSystem> systemCache_{256};
  LruCache<StationId, econ::StationEconomyState> stationEconomyCache_{512};
};

} // namespace stellar::sim
