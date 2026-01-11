#include "stellar/sim/MissionAssist.h"

#include "stellar/sim/MissionLogic.h"
#include "stellar/sim/Universe.h"

#include <cmath>
#include <iostream>

int test_mission_assist() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;

  Universe u(424242);

  // Find a nearby system with at least one station.
  const auto stubs = u.queryNearby({0, 0, 0}, 200.0, 64);
  if (stubs.empty()) {
    std::cerr << "[test_mission_assist] no systems returned\n";
    return 1;
  }

  const SystemStub* chosenStub = nullptr;
  for (const auto& s : stubs) {
    if (s.stationCount > 0) { chosenStub = &s; break; }
  }
  if (!chosenStub) {
    std::cerr << "[test_mission_assist] no systems with stations in query\n";
    return 1;
  }

  const auto& sys = u.getSystem(chosenStub->id, chosenStub);
  if (sys.stations.empty()) {
    std::cerr << "[test_mission_assist] stub said stations but system had none\n";
    return 1;
  }

  Mission m{};
  m.id = 1;
  m.type = MissionType::Delivery;
  m.commodity = econ::CommodityId::Food;
  m.units = 22.0;
  m.toSystem = sys.stub.id;
  m.toStation = sys.stations.front().id;

  std::array<double, econ::kCommodityCount> cargo{};

  // Basic deterministic scan.
  const double t = 12.25;
  const auto a = planMissionCargoSourcingForMission(u, sys, t, m, cargo);
  const auto b = planMissionCargoSourcingForMission(u, sys, t, m, cargo);

  if (!a.ok || !b.ok) {
    std::cerr << "[test_mission_assist] expected ok scan\n";
    ++fails;
  }

  if (std::abs(a.missingUnits - 22.0) > 1e-9) {
    std::cerr << "[test_mission_assist] missingUnits mismatch\n";
    ++fails;
  }

  if (a.candidates.size() != b.candidates.size()) {
    std::cerr << "[test_mission_assist] deterministic candidates size mismatch\n";
    ++fails;
  } else if (!a.candidates.empty()) {
    const auto& ca = a.candidates.front();
    const auto& cb = b.candidates.front();
    if (ca.systemId != cb.systemId || ca.stationId != cb.stationId) {
      std::cerr << "[test_mission_assist] deterministic best candidate mismatch\n";
      ++fails;
    }

    if (ca.inventoryUnits <= 0.0 || ca.askCr < 0.0 || ca.askEffCr < 0.0) {
      std::cerr << "[test_mission_assist] invalid quote fields\n";
      ++fails;
    }
    if (ca.askEffCr + 1e-9 < ca.askCr) {
      std::cerr << "[test_mission_assist] askEff should be >= ask\n";
      ++fails;
    }
  }

  // If you already have the cargo, the plan should be empty but still ok.
  cargo[(std::size_t)m.commodity] = 100.0;
  const auto done = planMissionCargoSourcingForMission(u, sys, t, m, cargo);
  if (!done.ok || done.missingUnits > 1e-9 || !done.candidates.empty()) {
    std::cerr << "[test_mission_assist] expected empty plan when not missing cargo\n";
    ++fails;
  }

  return fails;
}
