#include "stellar/sim/TrafficLanes.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/sim/Units.h"
#include "stellar/sim/WorldIds.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <unordered_set>

namespace stellar::sim {

static constexpr double kPi = std::numbers::pi_v<double>;

// Stable tags for behavior continuity across versions.
static constexpr std::string_view kConvoyTag     = "traffic_convoy_v1";
static constexpr std::string_view kLaneTagV1     = "traffic_lane_v1";
static constexpr std::string_view kLaneTagV2     = "traffic_lane_corridor_v2";
static constexpr std::string_view kLaneDirTag    = "traffic_lane_dir_v1";
static constexpr std::string_view kLaneArcTag    = "traffic_lane_arc_v1";
static constexpr std::string_view kLaneJitterTag = "traffic_lane_arcjitter_v1";

static const Station* findStation(const StarSystem& sys, StationId id) {
  for (const auto& st : sys.stations) {
    if (st.id == id) return &st;
  }
  return nullptr;
}

static int chooseWeighted(core::SplitMix64& rng, const std::array<double, econ::kCommodityCount>& w) {
  double sum = 0.0;
  for (double x : w) sum += std::max(0.0, x);
  if (sum <= 0.0) return -1;
  double r = rng.nextUnit() * sum;
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const double wi = std::max(0.0, w[i]);
    r -= wi;
    if (r <= 0.0) return (int)i;
  }
  return (int)econ::kCommodityCount - 1;
}

static int chooseWeightedStations(core::SplitMix64& rng, const std::vector<double>& w) {
  double sum = 0.0;
  for (double x : w) sum += std::max(0.0, x);
  if (sum <= 0.0) return -1;
  double r = rng.nextUnit() * sum;
  for (std::size_t i = 0; i < w.size(); ++i) {
    const double wi = std::max(0.0, w[i]);
    r -= wi;
    if (r <= 0.0) return (int)i;
  }
  return (int)w.size() - 1;
}

// Deterministic unit-vector sampler with fixed RNG consumption.
//
// This avoids rejection sampling so arc magnitude draws remain stable even if
// code changes the number of samples taken in the direction generator.
static math::Vec3d randomUnitVec(core::SplitMix64& rng) {
  // z in [-1,1], phi in [0,2pi)
  const double z = rng.range(-1.0, 1.0);
  const double phi = rng.range(0.0, 2.0 * kPi);

  const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
  const double x = r * std::cos(phi);
  const double y = z;
  const double zz = r * std::sin(phi);
  return {x, y, zz};
}

static core::u64 convoySeed(core::u64 universeSeed, SystemId systemId, int dayStamp) {
  core::u64 s = universeSeed;
  s ^= core::seedFromText(kConvoyTag);
  s = core::hashCombine(s, (core::u64)systemId);
  s = core::hashCombine(s, (core::u64)(core::u32)dayStamp);
  return s;
}

static core::u64 corridorLaneKey(SystemId systemId, StationId a, StationId b) {
  const StationId lo = std::min(a, b);
  const StationId hi = std::max(a, b);

  // Note: StationId/SystemId are already deterministic world ids derived from
  // the universe seed, so hashing them is sufficient to disambiguate lanes
  // across different universes.
  const core::u64 tag = core::seedFromText(kLaneTagV2);
  const core::u64 h = core::hashCombine((core::u64)systemId,
                                       core::hashCombine((core::u64)lo, (core::u64)hi));
  return makeDeterministicWorldId(tag, h);
}

static core::u64 legacyLaneKey(const TrafficConvoy& convoy, const StarSystem& system) {
  // Legacy: lane seed derived from convoy id.
  const core::u64 tag = core::seedFromText(kLaneTagV1);
  const core::u64 h = core::hashCombine((core::u64)system.stub.id, convoy.id);
  return makeDeterministicWorldId(tag, h);
}

std::vector<TrafficConvoy> generateTrafficConvoysForDay(core::u64 universeSeed,
                                                        const StarSystem& system,
                                                        int dayStamp,
                                                        const TrafficLaneParams& params) {
  std::vector<TrafficConvoy> out;

  if (system.stations.size() < 2) return out;

  core::SplitMix64 rng(convoySeed(universeSeed, system.stub.id, dayStamp));

  int count = params.convoysPerDayBase + (int)system.stations.size() * params.convoysPerStation;
  count += rng.range(-1, 2);
  count = std::clamp(count, 0, std::max(0, params.maxConvoysPerDay));

  // Avoid duplicate (from,to,commodity) triples inside the day stamp.
  std::unordered_set<core::u64> used;
  used.reserve((std::size_t)count * 2);

  out.reserve((std::size_t)count);

  for (int ci = 0; ci < count; ++ci) {
    const int nStations = (int)system.stations.size();

    int fromIdx = 0;
    int toIdx = 1;
    int commodityIdx = 0;
    core::u64 key = 0;

    // Pick a (from,to,commodity) triple with light de-duplication.
    // We intentionally keep this inexpensive — station counts are small.
    for (int tries = 0; tries < 8; ++tries) {
      fromIdx = rng.range(0, nStations - 1);
      toIdx = rng.range(0, nStations - 2);
      if (toIdx >= fromIdx) toIdx += 1;

      const Station& fromTry = system.stations[(std::size_t)fromIdx];

      // Choose a commodity biased by what the origin produces.
      std::array<double, econ::kCommodityCount> cWeights{};
      for (std::size_t k = 0; k < econ::kCommodityCount; ++k) {
        const double prod = fromTry.economyModel.productionPerDay[k];
        const double cons = fromTry.economyModel.consumptionPerDay[k];
        cWeights[k] = std::max(0.0, prod - 0.5 * cons);
      }

      commodityIdx = chooseWeighted(rng, cWeights);
      if (commodityIdx < 0) commodityIdx = rng.range(0, (int)econ::kCommodityCount - 1);

      // Choose a destination biased by who consumes that commodity.
      std::vector<double> dWeights;
      dWeights.resize(system.stations.size(), 0.0);
      for (std::size_t si = 0; si < system.stations.size(); ++si) {
        if ((int)si == fromIdx) continue;
        const Station& st = system.stations[si];
        const double prod = st.economyModel.productionPerDay[(std::size_t)commodityIdx];
        const double cons = st.economyModel.consumptionPerDay[(std::size_t)commodityIdx];
        dWeights[si] = std::max(0.0, cons - 0.5 * prod) + 0.05;
      }
      const int chosenTo = chooseWeightedStations(rng, dWeights);
      if (chosenTo >= 0 && chosenTo != fromIdx) {
        toIdx = chosenTo;
      }

      const Station& toTry = system.stations[(std::size_t)toIdx];

      key = core::hashCombine((core::u64)fromTry.id,
                             core::hashCombine((core::u64)toTry.id, (core::u64)(core::u32)commodityIdx));
      if (used.insert(key).second) {
        break;
      }
    }

    const Station& from = system.stations[(std::size_t)fromIdx];
    const Station& to = system.stations[(std::size_t)toIdx];
    const econ::CommodityId commodity = (econ::CommodityId)commodityIdx;

    // Schedule.
    const double departDay = (double)dayStamp + rng.nextUnit();
    const math::Vec3d p0 = stationPosKm(from, departDay);

    // Pick a target duration; clamp by speed bounds.
    const double minDur = std::max(1e-6, params.minDurationDays);
    const double maxDur = std::max(minDur, params.maxDurationDays);
    double durationDays = rng.range(minDur, maxDur);
    durationDays = std::max(1e-6, durationDays);

    const double vMin = std::max(0.0, params.speedMinKmS);
    const double vMax = std::max(vMin, params.speedMaxKmS);

    auto clampDurationDays = [&](double distKm, double durDays) {
      if (distKm <= 1e-6) return std::max(1e-6, durDays);
      double durSec = std::max(1e-6, durDays) * kSecondsPerDay;
      double speed = (durSec > 0.0) ? (distKm / durSec) : 0.0;

      if (speed < vMin && vMin > 1e-9) {
        speed = vMin;
        durSec = distKm / speed;
      }
      if (speed > vMax && vMax > 1e-9) {
        speed = vMax;
        durSec = distKm / speed;
      }
      return std::max(1e-6, durSec / kSecondsPerDay);
    };

    // Initial estimate based on destination position at departure.
    double distKm = (stationPosKm(to, departDay) - p0).length();
    durationDays = clampDurationDays(distKm, durationDays);

    // One refinement step based on destination position at arrival.
    double arriveDay = departDay + durationDays;
    distKm = (stationPosKm(to, arriveDay) - p0).length();
    durationDays = clampDurationDays(distKm, durationDays);
    arriveDay = departDay + durationDays;

    // Units: roughly a few hours of net export, with some noise.
    const double prod = from.economyModel.productionPerDay[(std::size_t)commodityIdx];
    const double cons = from.economyModel.consumptionPerDay[(std::size_t)commodityIdx];
    const double net = std::max(0.0, prod - cons);
    const double baseUnits = std::clamp(net * durationDays * 0.75, 6.0, 240.0);
    const double units = std::max(1.0, baseUnits * (0.75 + 0.55 * rng.nextUnit()));

    // Deterministic ID.
    core::u64 a = universeSeed ^ core::seedFromText(kConvoyTag);
    core::u64 b = (core::u64)system.stub.id;
    b = core::hashCombine(b, (core::u64)(core::u32)dayStamp);
    b = core::hashCombine(b, (core::u64)(core::u32)ci);
    b = core::hashCombine(b, (core::u64)from.id);
    b = core::hashCombine(b, (core::u64)to.id);
    b = core::hashCombine(b, (core::u64)(core::u32)commodityIdx);
    const core::u64 id = makeDeterministicWorldId(a, b);

    TrafficConvoy c;
    c.id = id;
    c.systemId = system.stub.id;
    c.fromStation = from.id;
    c.toStation = to.id;
    c.factionId = from.factionId;
    c.commodity = commodity;
    c.units = units;
    c.departDay = departDay;
    c.arriveDay = arriveDay;
    out.push_back(c);
  }

  // Stable ordering for deterministic output.
  std::sort(out.begin(), out.end(), [](const TrafficConvoy& a, const TrafficConvoy& b) {
    if (a.departDay != b.departDay) return a.departDay < b.departDay;
    return a.id < b.id;
  });

  return out;
}

TrafficLaneGeometry computeTrafficLaneGeometry(const TrafficConvoy& convoy,
                                               const StarSystem& system,
                                               const TrafficLaneParams& params) {
  TrafficLaneGeometry lane{};

  lane.systemId = convoy.systemId;
  lane.fromStation = convoy.fromStation;
  lane.toStation = convoy.toStation;
  lane.departDay = convoy.departDay;
  lane.arriveDay = convoy.arriveDay;

  const Station* from = findStation(system, convoy.fromStation);
  const Station* to = findStation(system, convoy.toStation);
  if (!from || !to) return lane;

  // Endpoints.
  const double depart = convoy.departDay;
  const double arrive = convoy.arriveDay;
  const double durDays = std::max(1e-9, arrive - depart);
  const double durSec = durDays * kSecondsPerDay;

  lane.p0Km = stationPosKm(*from, depart);
  lane.p1Km = stationPosKm(*to, arrive);
  const math::Vec3d delta = lane.p1Km - lane.p0Km;
  lane.distKm = delta.length();
  lane.durSec = durSec;

  if (lane.distKm > 1e-6) {
    lane.dir = delta * (1.0 / lane.distKm);
  } else {
    lane.dir = {0, 0, 1};
  }

  // Direction sign for dual-carriageway lanes.
  lane.directionSign = 1;
  if (params.dualCarriageway && convoy.fromStation != convoy.toStation) {
    lane.directionSign = (convoy.fromStation < convoy.toStation) ? 1 : -1;
  }

  // Stable key used for corridor grouping and sub-stream RNG.
  if (params.bundleByStationPair) {
    lane.laneKey = corridorLaneKey(system.stub.id, convoy.fromStation, convoy.toStation);
  } else {
    lane.laneKey = legacyLaneKey(convoy, system);
  }

  // Side direction.
  // Use independent sub-streams so future tweaks to the direction sampler don't
  // perturb arc magnitudes.
  core::SplitMix64 rngDir(core::hashCombine(lane.laneKey, core::seedFromText(kLaneDirTag)));
  const math::Vec3d rv = randomUnitVec(rngDir);

  math::Vec3d side = math::cross(lane.dir, rv);
  double sideLen = side.length();
  if (sideLen < 1e-6) {
    side = math::cross(lane.dir, {0, 1, 0});
    sideLen = side.length();
  }
  if (sideLen < 1e-6) {
    side = {1, 0, 0};
    sideLen = 1.0;
  }
  side = side * (1.0 / sideLen);
  side = side * (double)lane.directionSign;
  lane.side = side;

  // Up axis completes the frame. (Useful for optional future lane ribbon jitter.)
  math::Vec3d up = math::cross(lane.side, lane.dir);
  const double upLen = up.length();
  if (upLen > 1e-6) up = up * (1.0 / upLen);
  lane.up = up;

  // Arc magnitude: clamp by distance so short hops don't produce huge arcs.
  const double arcCapByDist = std::max(0.0, params.arcMaxFracOfDistance) * lane.distKm;
  const double arcCap = std::max(0.0, std::min(params.arcMaxKm, arcCapByDist));

  // arcMinKm is a "preferred" minimum, but we never exceed the cap.
  const double arcMin = std::clamp(params.arcMinKm, 0.0, arcCap);
  const double arcMax = std::max(arcMin, arcCap);

  core::SplitMix64 rngArc(core::hashCombine(lane.laneKey, core::seedFromText(kLaneArcTag)));
  lane.arcMagKm = (arcCap > 0.0) ? rngArc.range(arcMin, arcMax) : 0.0;

  // Optional per-convoy arc jitter (keeps overlapping convoys distinct).
  lane.arcJitterKm = 0.0;
  if (lane.arcMagKm > 0.0 && params.arcJitterFrac > 0.0) {
    const double jf = std::max(0.0, params.arcJitterFrac);
    core::SplitMix64 rngJ(core::hashCombine(convoy.id ^ core::seedFromText(kLaneJitterTag), lane.laneKey));
    const double s = (rngJ.nextUnit() * 2.0 - 1.0); // [-1,1]
    lane.arcJitterKm = lane.arcMagKm * (s * jf);
  }

  return lane;
}

TrafficConvoyState evaluateTrafficConvoy(const TrafficLaneGeometry& lane,
                                         double timeDays,
                                         const TrafficLaneParams&) {
  TrafficConvoyState st;

  const double depart = lane.departDay;
  const double arrive = lane.arriveDay;
  const double durDays = std::max(1e-9, arrive - depart);
  const double durSec = durDays * kSecondsPerDay;

  double t = 0.0;
  if (durDays > 0.0) t = (timeDays - depart) / durDays;

  st.active = (timeDays >= depart && timeDays <= arrive);
  st.progress01 = std::clamp(t, 0.0, 1.0);

  const math::Vec3d delta = lane.p1Km - lane.p0Km;
  st.distKm = lane.distKm;
  st.speedKmS = (durSec > 0.0) ? (st.distKm / durSec) : 0.0;

  st.dir = lane.dir;

  // Base linear interpolation along the chord.
  const math::Vec3d basePos = lane.p0Km + delta * st.progress01;
  const math::Vec3d baseVel = (durSec > 0.0) ? (delta * (1.0 / durSec)) : math::Vec3d{0, 0, 0};

  // Use an ease-in/out arc that has *zero* lateral velocity at the endpoints.
  // offset(phase) = sin^2(pi*phase) * arcMag
  const double phase = st.progress01;
  const double sp = std::sin(kPi * phase);
  const double cp = std::cos(kPi * phase);
  const double s2 = sp * sp;

  const double arcMag = std::max(0.0, lane.arcMagKm + lane.arcJitterKm);
  const math::Vec3d arcOffset = lane.side * (s2 * arcMag);

  // d/dt[sin^2(pi*phase)] = 2*sin(pi*phase)*cos(pi*phase) * pi / durSec
  const double ds2_dt = (durSec > 0.0) ? (2.0 * sp * cp * kPi / durSec) : 0.0;
  const math::Vec3d arcVel = lane.side * (ds2_dt * arcMag);

  st.posKm = basePos + arcOffset;
  st.velKmS = baseVel + arcVel;

  // When the convoy is outside its active window, treat it as stationary at its endpoint.
  // (This avoids leaking a non-zero velocity into UI/navigation when includeInactive=true.)
  if (!st.active) {
    st.velKmS = {0.0, 0.0, 0.0};
  }

  return st;
}

TrafficConvoyState evaluateTrafficConvoy(const TrafficConvoy& convoy,
                                         const StarSystem& system,
                                         double timeDays,
                                         const TrafficLaneParams& params) {
  const TrafficLaneGeometry lane = computeTrafficLaneGeometry(convoy, system, params);
  return evaluateTrafficConvoy(lane, timeDays, params);
}

std::vector<math::Vec3d> sampleTrafficLanePathKm(const TrafficLaneGeometry& lane,
                                                 int segments,
                                                 const TrafficLaneParams& params) {
  std::vector<math::Vec3d> pts;
  if (segments < 1) return pts;

  const double depart = lane.departDay;
  const double arrive = lane.arriveDay;
  const double dur = std::max(1e-12, arrive - depart);

  pts.reserve((std::size_t)segments + 1);
  for (int i = 0; i <= segments; ++i) {
    const double u = (double)i / (double)segments;
    const double t = depart + u * dur;
    const auto st = evaluateTrafficConvoy(lane, t, params);
    pts.push_back(st.posKm);
  }
  return pts;
}

std::vector<TrafficConvoyView> generateTrafficConvoys(core::u64 universeSeed,
                                                      const StarSystem& system,
                                                      double timeDays,
                                                      const TrafficLaneParams& params) {
  std::vector<TrafficConvoyView> out;
  if (system.stations.size() < 2) return out;

  const int day = (int)std::floor(timeDays);
  const int w = std::max(0, params.genWindowDays);

  for (int d = day - w; d <= day + w; ++d) {
    auto sched = generateTrafficConvoysForDay(universeSeed, system, d, params);
    for (const auto& c : sched) {
      TrafficConvoyState st = evaluateTrafficConvoy(c, system, timeDays, params);
      if (!params.includeInactive && !st.active) continue;
      TrafficConvoyView v;
      v.convoy = c;
      v.state = st;
      out.push_back(v);
    }
  }

  std::sort(out.begin(), out.end(), [](const TrafficConvoyView& a, const TrafficConvoyView& b) {
    if (a.convoy.departDay != b.convoy.departDay) return a.convoy.departDay < b.convoy.departDay;
    return a.convoy.id < b.convoy.id;
  });

  return out;
}

} // namespace stellar::sim
