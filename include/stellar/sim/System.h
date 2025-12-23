#pragma once

#include "stellar/econ/Economy.h"
#include "stellar/sim/Celestial.h"

#include <vector>

namespace stellar::sim {

struct Station {
  StationId id{0};
  std::string name;
  econ::StationType type{econ::StationType::Outpost};
  core::u32 factionId{0};       // 0 = independent
  double feeRate{0.0};          // market fee rate (0..1)
  econ::StationEconomyModel economyModel{};
};

struct StarSystem {
  SystemStub stub{};
  Star star{};
  std::vector<Planet> planets{};
  std::vector<Station> stations{};
};

} // namespace stellar::sim
