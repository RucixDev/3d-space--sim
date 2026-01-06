#include "stellar/sim/TrafficConvoyLayer.h"

#include "stellar/sim/Units.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static const Station* findStationById(const StarSystem& sys, StationId id) {
  for (const auto& st : sys.stations) {
    if (st.id == id) return &st;
  }
  return nullptr;
}

TrafficConvoy convoyFromShipment(const TrafficShipment& s) {
  TrafficConvoy c{};
  c.id = s.id;
  c.systemId = s.systemId;
  c.fromStation = s.fromStation;
  c.toStation = s.toStation;
  c.factionId = s.factionId;
  c.commodity = s.commodity;
  c.units = s.units;
  c.departDay = s.departDay;
  c.arriveDay = s.arriveDay;
  return c;
}

std::vector<TrafficConvoyView> generateTrafficConvoysFromLedger(const TrafficLedger& ledger,
                                                                const StarSystem& system,
                                                                double timeDays,
                                                                int windowDays,
                                                                bool includeInactive,
                                                                const TrafficLaneParams& laneParams) {
  std::vector<TrafficConvoyView> out;
  if (timeDays < 0.0) return out;
  if (system.stations.size() < 2) return out;

  // Query with includeInactive=true so we can recover shipments that may be missing
  // schedule metadata (older/corrupted save files). We'll apply the active filter
  // ourselves after hydrating the schedule.
  const auto shipments = ledger.query(system.stub.id, timeDays, windowDays, /*includeInactive=*/true);
  out.reserve(shipments.size());

  for (auto s : shipments) {
    // Defensive: ignore broken references.
    const Station* from = findStationById(system, s.fromStation);
    const Station* to = findStationById(system, s.toStation);
    if (!from || !to) continue;
    if (s.fromStation == s.toStation) continue;

    // Ensure we have a usable schedule (required for replay/visualization).
    if (!hydrateShipmentScheduleFromId(s, system, ledger.params)) continue;

    TrafficConvoyView v;
    v.convoy = convoyFromShipment(s);
    v.state = evaluateTrafficConvoy(v.convoy, system, timeDays, laneParams);
    if (!includeInactive && !v.state.active) continue;
    out.push_back(std::move(v));
  }

  std::sort(out.begin(), out.end(), [](const TrafficConvoyView& a, const TrafficConvoyView& b) {
    if (a.convoy.departDay != b.convoy.departDay) return a.convoy.departDay < b.convoy.departDay;
    return a.convoy.id < b.convoy.id;
  });

  return out;
}

std::vector<math::Vec3d> sampleTrafficConvoyPathKm(const TrafficConvoy& convoy,
                                                   const StarSystem& system,
                                                   int segments,
                                                   const TrafficLaneParams& laneParams) {
  std::vector<math::Vec3d> pts;
  if (segments < 1) return pts;

  const double depart = convoy.departDay;
  const double arrive = convoy.arriveDay;
  const double dur = std::max(1e-12, arrive - depart);

  pts.reserve((std::size_t)segments + 1);
  for (int i = 0; i <= segments; ++i) {
    const double u = (double)i / (double)segments;
    const double t = depart + u * dur;
    const auto st = evaluateTrafficConvoy(convoy, system, t, laneParams);
    pts.push_back(st.posKm);
  }
  return pts;
}

} // namespace stellar::sim