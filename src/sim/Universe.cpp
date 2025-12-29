#include "stellar/sim/Universe.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Log.h"
#include "stellar/core/Random.h"
#include "stellar/proc/NameGenerator.h"
#include "stellar/proc/SystemGenerator.h"
#include "stellar/econ/Economy.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <unordered_set>
#include <sstream>

namespace stellar::sim {

static proc::SectorCoord decodeSector(SystemId id, core::u32& outLocalIndex) {
  const auto unbias = [](core::u16 b) -> core::i32 { return static_cast<core::i32>(b) - 32768; };

  const core::u16 bx = static_cast<core::u16>((id >> 48) & 0xFFFFu);
  const core::u16 by = static_cast<core::u16>((id >> 32) & 0xFFFFu);
  const core::u16 bz = static_cast<core::u16>((id >> 16) & 0xFFFFu);
  const core::u16 bi = static_cast<core::u16>((id >> 0)  & 0xFFFFu);

  outLocalIndex = static_cast<core::u32>(bi);
  return proc::SectorCoord{ unbias(bx), unbias(by), unbias(bz) };
}

static StarClass pickStarClass(core::SplitMix64& rng) {
  const double r = rng.nextDouble();
  // Very rough main-sequence-ish distribution.
  if (r < 0.0003) return StarClass::O;
  if (r < 0.0016) return StarClass::B;
  if (r < 0.006)  return StarClass::A;
  if (r < 0.03)   return StarClass::F;
  if (r < 0.10)   return StarClass::G;
  if (r < 0.30)   return StarClass::K;
  return StarClass::M;
}

Universe::Universe(core::u64 seed, proc::GalaxyParams params)
: seed_(seed),
  galaxyParams_(params),
  galaxyGen_(seed, params) {
  factions_ = generateFactions(seed_, 8);
}

void Universe::setCacheCaps(std::size_t sectorCap, std::size_t systemCap, std::size_t stationCap) {
  sectorCache_.setCapacity(sectorCap);
  systemCache_.setCapacity(systemCap);
  stationEconomyCache_.setCapacity(stationCap);
}

Universe::CacheStats Universe::cacheStats() const {
  CacheStats s;

  {
    const auto st = sectorCache_.stats();
    s.sectors.capacity = st.capacity;
    s.sectors.size = st.size;
    s.sectors.hits = st.hits;
    s.sectors.misses = st.misses;
    s.sectors.puts = st.puts;
    s.sectors.evictions = st.evictions;
  }

  {
    const auto st = systemCache_.stats();
    s.systems.capacity = st.capacity;
    s.systems.size = st.size;
    s.systems.hits = st.hits;
    s.systems.misses = st.misses;
    s.systems.puts = st.puts;
    s.systems.evictions = st.evictions;
  }

  {
    const auto st = stationEconomyCache_.stats();
    s.stationEconomies.capacity = st.capacity;
    s.stationEconomies.size = st.size;
    s.stationEconomies.hits = st.hits;
    s.stationEconomies.misses = st.misses;
    s.stationEconomies.puts = st.puts;
    s.stationEconomies.evictions = st.evictions;
  }

  return s;
}

void Universe::resetCacheStats() {
  sectorCache_.resetStats();
  systemCache_.resetStats();
  stationEconomyCache_.resetStats();
}


const proc::Sector& Universe::sector(const proc::SectorCoord& coord) {
  if (auto* cached = sectorCache_.get(coord)) return *cached;
  auto sec = galaxyGen_.generateSector(coord, factions_);
  return sectorCache_.put(coord, std::move(sec));
}

std::vector<SystemStub> Universe::queryNearby(const math::Vec3d& posLy,
                                              double radiusLy,
                                              std::size_t maxResults) {
  std::vector<SystemStub> out;
  if (radiusLy <= 0.0 || maxResults == 0) return out;

  const double r2 = radiusLy * radiusLy;
  const double s = galaxyParams_.sectorSizeLy;

  // Conservative bounds on sector coords that could intersect the query sphere.
  const proc::SectorCoord minC{
    static_cast<core::i32>(std::floor((posLy.x - radiusLy) / s)),
    static_cast<core::i32>(std::floor((posLy.y - radiusLy) / s)),
    static_cast<core::i32>(std::floor((posLy.z - radiusLy) / s)),
  };
  const proc::SectorCoord maxC{
    static_cast<core::i32>(std::floor((posLy.x + radiusLy) / s)),
    static_cast<core::i32>(std::floor((posLy.y + radiusLy) / s)),
    static_cast<core::i32>(std::floor((posLy.z + radiusLy) / s)),
  };

  const auto inBounds = [&](const proc::SectorCoord& c) {
    return c.x >= minC.x && c.x <= maxC.x &&
           c.y >= minC.y && c.y <= maxC.y &&
           c.z >= minC.z && c.z <= maxC.z;
  };

  const auto axisDist = [](double v, double lo, double hi) {
    if (v < lo) return lo - v;
    if (v > hi) return v - hi;
    return 0.0;
  };

  // Lower bound on the distance from posLy to *any point* in a sector cube.
  // Used for best-first sector expansion + safe early-out when maxResults is reached.
  const auto minDist2ToSector = [&](const proc::SectorCoord& c) -> double {
    const double x0 = static_cast<double>(c.x) * s;
    const double y0 = static_cast<double>(c.y) * s;
    const double z0 = static_cast<double>(c.z) * s;

    const double x1 = x0 + s;
    const double y1 = y0 + s;
    const double z1 = z0 + s;

    const double dx = axisDist(posLy.x, x0, x1);
    const double dy = axisDist(posLy.y, y0, y1);
    const double dz = axisDist(posLy.z, z0, z1);
    return dx*dx + dy*dy + dz*dz;
  };

  struct SectorItem {
    proc::SectorCoord c{};
    double minD2{0.0};
  };

  // Min-heap by minD2 (then coord) so we visit the nearest possible sectors first.
  struct SectorItemGreater {
    bool operator()(const SectorItem& a, const SectorItem& b) const {
      if (a.minD2 != b.minD2) return a.minD2 > b.minD2;
      if (a.c.x != b.c.x) return a.c.x > b.c.x;
      if (a.c.y != b.c.y) return a.c.y > b.c.y;
      return a.c.z > b.c.z;
    }
  };

  std::priority_queue<SectorItem, std::vector<SectorItem>, SectorItemGreater> open;
  std::unordered_set<proc::SectorCoord, proc::SectorCoordHash> seen;
  seen.reserve(1024);

  // Keep the best maxResults stubs seen so far ("best" = smallest distance, then id).
  struct Candidate {
    SystemStub stub{};
    double d2{0.0};
  };

  const auto better = [](const Candidate& a, const Candidate& b) {
    if (a.d2 != b.d2) return a.d2 < b.d2;
    return a.stub.id < b.stub.id;
  };

  // Max-heap by (d2,id) so top() is the current *worst* of the kept candidates.
  struct CandidateLess {
    bool operator()(const Candidate& a, const Candidate& b) const {
      if (a.d2 != b.d2) return a.d2 < b.d2;
      return a.stub.id < b.stub.id;
    }
  };

  std::priority_queue<Candidate, std::vector<Candidate>, CandidateLess> best;

  // Dynamic cutoff: before we have maxResults candidates we must consider the full radius.
  // Once we have maxResults, we can safely ignore sectors whose AABB is farther than our
  // current worst candidate (because they cannot contain a closer stub).
  double stopD2 = r2;

  const auto tryPush = [&](const proc::SectorCoord& c) {
    if (!inBounds(c)) return;
    if (!seen.insert(c).second) return;
    const double md2 = minDist2ToSector(c);
    if (md2 <= stopD2) open.push(SectorItem{c, md2});
  };

  const proc::SectorCoord startC = galaxyGen_.sectorOf(posLy);
  tryPush(startC);

  while (!open.empty()) {
    // Refresh the cutoff.
    stopD2 = r2;
    if (best.size() >= maxResults) stopD2 = std::min(stopD2, best.top().d2);

    const SectorItem cur = open.top();
    if (cur.minD2 > stopD2) break; // remaining sectors cannot improve the top-N
    open.pop();

    const proc::Sector& sec = sector(cur.c);
    for (const auto& stub : sec.systems) {
      const double dd = (stub.posLy - posLy).lengthSq();
      if (dd > r2) continue;

      Candidate cand{stub, dd};

      if (best.size() < maxResults) {
        best.push(std::move(cand));
      } else {
        const Candidate& worst = best.top();
        if (better(cand, worst)) {
          best.pop();
          best.push(std::move(cand));
        }
      }
    }

    // Update cutoff after processing this sector (it may have improved).
    stopD2 = r2;
    if (best.size() >= maxResults) stopD2 = std::min(stopD2, best.top().d2);

    // Expand in the 6 face-adjacent directions. The set of sectors that intersect
    // a sphere (or any smaller cutoff sphere) is connected under face adjacency.
    tryPush(proc::SectorCoord{cur.c.x + 1, cur.c.y, cur.c.z});
    tryPush(proc::SectorCoord{cur.c.x - 1, cur.c.y, cur.c.z});
    tryPush(proc::SectorCoord{cur.c.x, cur.c.y + 1, cur.c.z});
    tryPush(proc::SectorCoord{cur.c.x, cur.c.y - 1, cur.c.z});
    tryPush(proc::SectorCoord{cur.c.x, cur.c.y, cur.c.z + 1});
    tryPush(proc::SectorCoord{cur.c.x, cur.c.y, cur.c.z - 1});
  }

  // Unheap + sort by the required stable order (distance, then id).
  std::vector<Candidate> items;
  items.reserve(best.size());
  while (!best.empty()) {
    items.push_back(best.top());
    best.pop();
  }

  std::sort(items.begin(), items.end(), [&](const Candidate& a, const Candidate& b) {
    return better(a, b);
  });

  out.reserve(items.size());
  for (auto& it : items) out.push_back(std::move(it.stub));
  return out;
}

std::optional<SystemStub> Universe::findClosestSystem(const math::Vec3d& posLy, double maxRadiusLy) {
  auto list = queryNearby(posLy, maxRadiusLy, 64);
  if (list.empty()) return std::nullopt;

  // queryNearby is distance-sorted already
  return list.front();
}

const StarSystem& Universe::getSystem(SystemId id, const SystemStub* hintStub) {
  if (auto* cached = systemCache_.get(id)) return *cached;

  SystemStub stub{};
  bool haveStub = false;

  if (hintStub && hintStub->id == id) {
    stub = *hintStub;
    haveStub = true;
  } else {
    // Decode sector from id and search it.
    core::u32 localIndex = 0;
    const proc::SectorCoord c = decodeSector(id, localIndex);

    const proc::Sector& sec = sector(c);

    auto it = std::lower_bound(sec.systems.begin(), sec.systems.end(), id,
                               [](const SystemStub& s, SystemId idv) { return s.id < idv; });
    if (it != sec.systems.end() && it->id == id) {
      stub = *it;
      haveStub = true;
    } else {
      std::ostringstream oss;
      oss << "Universe::getSystem: stub not found for id=" << id
          << " (sector " << c.x << "," << c.y << "," << c.z
          << " localIndex=" << localIndex << "). Generating fallback stub.";
      stellar::core::log(stellar::core::LogLevel::Warn, oss.str());

      // Fallback: deterministic from id + seed.
      //
      // We also synthesize a plausible disc position so legacy/unresolvable ids don't
      // all collapse to the origin (which can make debugging/nav UI confusing).
      stub.id = id;
      stub.seed = core::hashCombine(seed_, static_cast<core::u64>(id));
      proc::NameGenerator ng(stub.seed);
      stub.name = ng.systemName();

      core::SplitMix64 frng(core::hashCombine(stub.seed, core::seedFromText("fallback_stub")));

      const double r = galaxyParams_.radiusLy * std::sqrt(frng.nextDouble());
      const double a = frng.nextDouble() * 6.283185307179586;
      const double z = (frng.nextDouble() - 0.5) * galaxyParams_.thicknessLy;
      stub.posLy = { r * std::cos(a), r * std::sin(a), z };

      stub.primaryClass = pickStarClass(frng);
      stub.planetCount = frng.range(0, 12);
      stub.stationCount = std::max(1, frng.range(0, 3));
      stub.factionId = 0;
      haveStub = true;


    }
  }

  (void)haveStub;

  StarSystem sys = proc::generateSystem(stub, factions_);
  return systemCache_.put(id, std::move(sys));
}

econ::StationEconomyState& Universe::stationEconomy(const Station& station, double timeDays) {
  econ::StationEconomyState* st = stationEconomyCache_.get(station.id);
  if (!st) {
    core::SplitMix64 rng(core::hashCombine(seed_, static_cast<core::u64>(station.id)));
    econ::StationEconomyState init = econ::makeInitialState(station.economyModel, rng);
    st = &stationEconomyCache_.put(station.id, std::move(init));
  }

  core::SplitMix64 rng(core::hashCombine(seed_, static_cast<core::u64>(station.id)));
  econ::updateEconomyTo(*st, station.economyModel, timeDays, rng);
  return *st;
}

std::vector<StationEconomyOverride> Universe::exportStationOverrides() const {
  std::vector<StationEconomyOverride> out;

  auto snap = stationEconomyCache_.snapshot();
  out.reserve(snap.size());
  for (auto& kv : snap) {
    StationEconomyOverride ov{};
    ov.stationId = kv.first;
    ov.state = std::move(kv.second);
    out.push_back(std::move(ov));
  }

  return out;
}

void Universe::importStationOverrides(const std::vector<StationEconomyOverride>& overrides) {
  for (const auto& ov : overrides) {
    stationEconomyCache_.put(ov.stationId, ov.state);
  }
}

} // namespace stellar::sim
