#include "test_harness.h"

#include "stellar/sim/MissionBriefing.h"

#include "stellar/econ/Commodity.h"
#include "stellar/sim/Contraband.h"
#include "stellar/sim/Universe.h"

#include <cmath>
#include <string>
#include <vector>

using namespace stellar;

static bool in01(double x) {
  return std::isfinite(x) && x >= -1e-6 && x <= 1.0 + 1e-6;
}

int test_mission_briefing() {
  int failures = 0;

  sim::Universe u(1337ull);

  const auto stubs = u.queryNearby(math::Vec3d{0.0, 0.0, 0.0}, 420.0, 220);
  CHECK(!stubs.empty());
  if (stubs.empty()) return failures;

  const sim::SystemStub originStub = stubs.front();
  const auto& originSys = u.getSystem(originStub.id, &originStub);
  CHECK(!originSys.stations.empty());
  if (originSys.stations.empty()) return failures;

  const sim::Station originStation = originSys.stations.front();
  const core::u64 seed = u.seed();

  // Find a destination station with at least one illegal commodity that is legal at origin.
  bool found = false;
  sim::SystemStub destStub{};
  sim::Station destStation{};
  econ::CommodityId chosenCid = econ::CommodityId::Food;

  for (int pass = 0; pass < 2 && !found; ++pass) {
    for (std::size_t si = 0; si < stubs.size() && !found; ++si) {
      const auto& stub = stubs[si];
      if (stub.id == 0 || stub.id == originStub.id) continue;

      const auto& sys = u.getSystem(stub.id, &stub);
      if (sys.stations.empty()) continue;
      if (pass == 0 && sys.stations.size() != 1) continue;

      for (const auto& st : sys.stations) {
        if (st.id == 0) continue;

        const core::u32 illegalMask = sim::illegalCommodityMaskForStation(seed, st.factionId, st.id, st.type);
        if (illegalMask == 0u) continue;

        const std::size_t maxBits = std::min<std::size_t>(econ::kCommodityCount, 32);
        for (std::size_t i = 0; i < maxBits; ++i) {
          if ((illegalMask & ((core::u32)1u << (core::u32)i)) == 0u) continue;
          const auto cid = static_cast<econ::CommodityId>(i);

          // Must be legal at origin so the mission is a true "smuggle" case.
          if (sim::isIllegalCommodityAtStation(seed,
                                              originStation.factionId,
                                              originStation.id,
                                              originStation.type,
                                              cid)) {
            continue;
          }

          found = true;
          destStub = stub;
          destStation = st;
          chosenCid = cid;
          break;
        }

        if (found) break;
      }
    }
  }

  CHECK(found);
  if (!found) return failures;

  const double timeDays = 5.25;
  const double rep = 0.0;

  sim::Mission m{};
  m.id = 42;
  m.type = sim::MissionType::Smuggle;
  m.factionId = originStation.factionId;
  m.toSystem = destStub.id;
  m.toStation = destStation.id;
  m.commodity = chosenCid;
  m.units = 12.0;
  m.reward = 25000.0;
  m.deadlineDay = timeDays + 2.0;
  m.cargoProvided = true;

  sim::MissionBriefingParams p;
  p.applySecurityDeltas = false; // keep this test independent of dynamic delta setup
  p.useMarkup = true;

  const auto b1 = sim::generateMissionBriefing(u, originSys, originStation, timeDays, rep, m, {}, p);
  const auto b2 = sim::generateMissionBriefing(u, originSys, originStation, timeDays, rep, m, {}, p);

  CHECK(!b1.titleMarkup.empty());
  CHECK(!b1.synopsisMarkup.empty());
  CHECK(!b1.contractCode.empty());
  CHECK(!b1.contactName.empty());

  CHECK(b1.contractCode == b2.contractCode);
  CHECK(b1.contactName == b2.contactName);
  CHECK(b1.titleMarkup == b2.titleMarkup);
  CHECK(b1.synopsisMarkup == b2.synopsisMarkup);

  CHECK(in01(b1.risk.overall01));
  CHECK(in01(b1.risk.danger01));
  CHECK(in01(b1.risk.lawRisk01));
  CHECK(in01(b1.risk.blackMarketAccess01));
  CHECK(in01(b1.risk.stingChance01));

  // For smuggling the law risk should be present and the commodity should be illegal at destination.
  CHECK(b1.risk.lawRisk01 > 1e-4);
  CHECK(b1.risk.contrabandIllegalAtDestination);

  bool sawContraband = false;
  for (const auto& line : b1.bulletsMarkup) {
    if (line.find("CONTRABAND") != std::string::npos) {
      sawContraband = true;
      break;
    }
  }
  CHECK(sawContraband);

  return failures;
}
