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

private:
  template <class Key, class Value, class Hash = std::hash<Key>>
  class LruCache {
  public:
    explicit LruCache(std::size_t cap = 128) : cap_(cap) {}

    void setCapacity(std::size_t cap) { cap_ = cap; evictIfNeeded(); }
    std::size_t capacity() const { return cap_; }
    std::size_t size() const { return map_.size(); }

    Value* get(const Key& key) {
      auto it = map_.find(key);
      if (it == map_.end()) return nullptr;
      touch(it);
      return &it->second.value;
    }

    const Value* get(const Key& key) const {
      auto it = map_.find(key);
      if (it == map_.end()) return nullptr;
      // const cache doesn't update LRU; ok for now
      return &it->second.value;
    }

    Value& put(const Key& key, Value value) {
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
      map_.emplace(key, std::move(e));

      evictIfNeeded();
      return map_.find(key)->second.value;
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
      }
    }

    std::size_t cap_{128};
    std::list<Key> order_{};
    std::unordered_map<Key, Entry, Hash> map_{};
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
