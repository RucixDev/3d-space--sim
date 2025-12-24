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

  // In-system physical placement (simple Keplerian orbit around the primary star).
  // Units:
  //  - orbit.semiMajorAxisAU: AU
  //  - angles: radians
  //  - orbit.periodDays: days
  OrbitElements orbit{};

  // Physical size (used for simple collision / docking distances).
  double radiusKm{6.0};

  // Docking / approach corridor parameters (local frame defined by station position
  // and the star at the origin; corridor axis points away from the star).
  double corridorLengthKm{120.0};
  double corridorRadiusKm{25.0};
  double corridorSpeedLimitKmS{0.12};
  double corridorAlignCos{0.90};

  // Comms range required to request docking.
  double commsRangeKm{2500.0};
};

struct StarSystem {
  SystemStub stub{};
  Star star{};
  std::vector<Planet> planets{};
  std::vector<Station> stations{};
};

} // namespace stellar::sim
