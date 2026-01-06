#include "stellar/sim/TrafficLedger.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/sim/Units.h"
#include "stellar/sim/WorldIds.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

// Stable tags for behavior continuity across versions.
static constexpr std::string_view kShipmentTag  = "npc_trade_shipment_v1";
static constexpr std::string_view kScheduleTag  = "npc_trade_schedule_v1";

static core::u64 shipmentId(core::u64 universeSeed,
                            SystemId systemId,
                            int dayStamp,
                            int runIndex,
                            StationId from,
                            StationId to,
                            econ::CommodityId commodity) {
  core::u64 a = universeSeed ^ core::seedFromText(kShipmentTag);

  core::u64 b = (core::u64)systemId;
  b = core::hashCombine(b, (core::u64)(core::u32)dayStamp);
  b = core::hashCombine(b, (core::u64)(core::u32)runIndex);
  b = core::hashCombine(b, (core::u64)from);
  b = core::hashCombine(b, (core::u64)to);
  b = core::hashCombine(b, (core::u64)(core::u32)commodity);

  return makeDeterministicWorldId(a, b);
}

TrafficShipment makeNpcTradeShipment(core::u64 universeSeed,
                                    const StarSystem& system,
                                    int dayStamp,
                                    int runIndex,
                                    const Station& from,
                                    const Station& to,
                                    econ::CommodityId commodity,
                                    double units,
                                    const TrafficLedgerParams& params) {
  TrafficShipment s{};
  s.systemId = system.stub.id;
  s.dayStamp = dayStamp;
  s.fromStation = from.id;
  s.toStation = to.id;
  s.factionId = from.factionId;
  s.commodity = commodity;
  s.units = units;

  s.id = shipmentId(universeSeed, system.stub.id, dayStamp, runIndex, from.id, to.id, commodity);

  // Schedule metadata is derived from the shipment id so that recording shipments
  // does NOT perturb the traffic simulation RNG sequence.
  core::SplitMix64 rng(core::hashCombine(s.id ^ core::seedFromText(kScheduleTag), (core::u64)system.stub.id));

  s.departDay = (double)dayStamp + rng.nextUnit();

  const math::Vec3d p0 = stationPosKm(from, s.departDay);

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
  double distKm = (stationPosKm(to, s.departDay) - p0).length();
  durationDays = clampDurationDays(distKm, durationDays);

  // One refinement step based on destination position at arrival.
  double arriveDay = s.departDay + durationDays;
  distKm = (stationPosKm(to, arriveDay) - p0).length();
  durationDays = clampDurationDays(distKm, durationDays);
  arriveDay = s.departDay + durationDays;

  // Final metrics captured at the refined arrival time.
  distKm = (stationPosKm(to, arriveDay) - p0).length();
  const double durSec = durationDays * kSecondsPerDay;

  s.distKm = distKm;
  s.speedKmS = (durSec > 0.0) ? (distKm / durSec) : 0.0;
  s.arriveDay = arriveDay;
  return s;
}


static const Station* findStationById(const StarSystem& sys, StationId id) {
  for (const auto& st : sys.stations) {
    if (st.id == id) return &st;
  }
  return nullptr;
}

bool shipmentScheduleLooksValid(const TrafficShipment& s) {
  if (!std::isfinite(s.departDay) || !std::isfinite(s.arriveDay)) return false;
  if (s.arriveDay <= s.departDay + 1e-9) return false;

  // makeNpcTradeShipment() schedules departDay within [dayStamp, dayStamp+1).
  // If the recorded value is wildly outside that range, treat it as missing/corrupt.
  const double day = (double)s.dayStamp;
  if (std::abs(s.departDay - day) > 2.0) return false;

  return true;
}

bool hydrateShipmentScheduleFromId(TrafficShipment& s, const StarSystem& system, const TrafficLedgerParams& params) {
  if (shipmentScheduleLooksValid(s)) return true;

  if (s.id == 0) return false;
  if (s.systemId != 0 && s.systemId != system.stub.id) return false;

  const Station* from = findStationById(system, s.fromStation);
  const Station* to = findStationById(system, s.toStation);
  if (!from || !to) return false;

  // Schedule metadata is derived from the shipment id so that recording shipments
  // does NOT perturb the traffic simulation RNG sequence.
  core::SplitMix64 rng(core::hashCombine(s.id ^ core::seedFromText(kScheduleTag), (core::u64)system.stub.id));

  s.departDay = (double)s.dayStamp + rng.nextUnit();

  const math::Vec3d p0 = stationPosKm(*from, s.departDay);

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
  double distKm = (stationPosKm(*to, s.departDay) - p0).length();
  durationDays = clampDurationDays(distKm, durationDays);

  // One refinement step based on destination position at arrival.
  double arriveDay = s.departDay + durationDays;
  distKm = (stationPosKm(*to, arriveDay) - p0).length();
  durationDays = clampDurationDays(distKm, durationDays);
  arriveDay = s.departDay + durationDays;

  // Final metrics captured at the refined arrival time.
  distKm = (stationPosKm(*to, arriveDay) - p0).length();
  const double durSec = durationDays * kSecondsPerDay;

  s.distKm = distKm;
  s.speedKmS = (durSec > 0.0) ? (distKm / durSec) : 0.0;
  s.arriveDay = arriveDay;

  return true;
}

void TrafficLedger::clear() {
  shipments.clear();
}

void TrafficLedger::record(const TrafficShipment& s) {
  // Guard against accidental duplication if the caller re-simulates a day.
  // (Counts are small, so an O(n) check is fine.)
  for (const auto& e : shipments) {
    if (e.id == s.id) return;
  }
  shipments.push_back(s);
}

void TrafficLedger::record(TrafficShipment&& s) {
  // Guard against accidental duplication if the caller re-simulates a day.
  // (Counts are small, so an O(n) check is fine.)
  for (const auto& e : shipments) {
    if (e.id == s.id) return;
  }
  shipments.push_back(std::move(s));
}

void TrafficLedger::prune(double timeDays) {
  if (params.keepDays <= 0) {
    shipments.clear();
    return;
  }
  if (timeDays < 0.0) return;

  const int day = (int)std::floor(timeDays);
  const int minDay = day - params.keepDays;

  shipments.erase(std::remove_if(shipments.begin(), shipments.end(), [&](const TrafficShipment& s) {
    return s.dayStamp < minDay;
  }),
                 shipments.end());
}

std::vector<TrafficShipment> TrafficLedger::query(SystemId systemId,
                                                  double timeDays,
                                                  int windowDays,
                                                  bool includeInactive) const {
  std::vector<TrafficShipment> out;
  if (shipments.empty()) return out;
  if (timeDays < 0.0) return out;

  const int day = (int)std::floor(timeDays);
  const int w = std::max(0, windowDays);
  const int minDay = day - w;
  const int maxDay = day + w;

  for (const auto& s : shipments) {
    if (s.systemId != systemId) continue;
    if (s.dayStamp < minDay || s.dayStamp > maxDay) continue;
    if (!includeInactive) {
      if (!(timeDays >= s.departDay && timeDays <= s.arriveDay)) continue;
    }
    out.push_back(s);
  }

  std::sort(out.begin(), out.end(), [](const TrafficShipment& a, const TrafficShipment& b) {
    if (a.departDay != b.departDay) return a.departDay < b.departDay;
    return a.id < b.id;
  });

  return out;
}

} // namespace stellar::sim
