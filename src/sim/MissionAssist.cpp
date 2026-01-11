#include "stellar/sim/MissionAssist.h"

#include "stellar/econ/Market.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static bool commodityValid(stellar::econ::CommodityId id) {
  const auto idx = (std::size_t)id;
  return idx < stellar::econ::kCommodityCount;
}

MissionCargoSourcePlan planMissionCargoSourcing(Universe& universe,
                                                const StarSystem& originSystem,
                                                double timeDays,
                                                econ::CommodityId commodity,
                                                double missingUnits,
                                                const MissionCargoSourceParams& params) {
  MissionCargoSourcePlan out{};

  if (!std::isfinite(timeDays)) timeDays = 0.0;
  if (!std::isfinite(missingUnits)) missingUnits = 0.0;
  missingUnits = std::max(0.0, missingUnits);

  out.commodity = commodity;
  out.missingUnits = missingUnits;
  if (commodityValid(commodity)) {
    out.missingMassKg = missingUnits * econ::commodityDef(commodity).massKg;
  }

  if (missingUnits <= 1e-9) {
    out.ok = true;
    return out;
  }
  if (!commodityValid(commodity)) {
    out.ok = false;
    return out;
  }

  const double rLy = std::max(0.0, params.searchRadiusLy);
  const std::size_t maxSys = std::max<std::size_t>(1, params.maxSystems);
  const auto stubs = universe.queryNearby(originSystem.stub.posLy, rLy, maxSys);

  out.candidates.reserve(64);

  for (const auto& stub : stubs) {
    if (!params.includeCurrentSystem && stub.id == originSystem.stub.id) continue;

    const auto& sys = universe.getSystem(stub.id, &stub);
    if (sys.stations.empty()) continue;

    const double distLy = (stub.posLy - originSystem.stub.posLy).length();

    for (const auto& st : sys.stations) {
      auto& econState = universe.stationEconomy(st, timeDays);
      const auto q = econ::quote(econState, st.economyModel, commodity, params.bidAskSpread);
      const double inv = std::max(0.0, q.inventory);
      if (inv <= 1e-6) continue;
      if (params.requireEnoughInventory && inv + 1e-6 < missingUnits) continue;

      MissionCargoSourceCandidate c;
      c.systemId = sys.stub.id;
      c.stationId = st.id;
      c.stationType = st.type;
      c.distanceLy = distLy;
      c.inventoryUnits = inv;
      c.askCr = std::max(0.0, q.ask);
      c.feeRate = std::clamp(st.feeRate, 0.0, 1.0);
      c.askEffCr = c.askCr * (1.0 + c.feeRate);

      out.candidates.push_back(std::move(c));
    }
  }

  // Stable, deterministic ordering.
  std::stable_sort(out.candidates.begin(), out.candidates.end(), [&](const auto& a, const auto& b) {
    const bool aEnough = (a.inventoryUnits + 1e-6 >= missingUnits);
    const bool bEnough = (b.inventoryUnits + 1e-6 >= missingUnits);
    if (aEnough != bEnough) return aEnough > bEnough;

    if (a.distanceLy != b.distanceLy) return a.distanceLy < b.distanceLy;
    if (a.askEffCr != b.askEffCr) return a.askEffCr < b.askEffCr;
    if (a.systemId != b.systemId) return a.systemId < b.systemId;
    return a.stationId < b.stationId;
  });

  if (out.candidates.size() > params.maxResults) {
    out.candidates.resize(params.maxResults);
  }

  out.ok = true;
  return out;
}

MissionCargoSourcePlan planMissionCargoSourcingForMission(Universe& universe,
                                                          const StarSystem& originSystem,
                                                          double timeDays,
                                                          const Mission& mission,
                                                          const std::array<double, econ::kCommodityCount>& cargo,
                                                          const MissionCargoSourceParams& params) {
  const bool cargoMission =
    (mission.type == MissionType::Delivery) ||
    (mission.type == MissionType::MultiDelivery) ||
    (mission.type == MissionType::Smuggle);

  if (!cargoMission) {
    return {};
  }

  const auto idx = (std::size_t)mission.commodity;
  if (idx >= econ::kCommodityCount) {
    return {};
  }

  const double have = std::max(0.0, cargo[idx]);
  const double need = std::max(0.0, mission.units);
  const double missing = std::max(0.0, need - have);

  return planMissionCargoSourcing(universe, originSystem, timeDays, mission.commodity, missing, params);
}

} // namespace stellar::sim
