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

  // Physical placement (currently: heliocentric orbit).
  OrbitElements orbit{};

  // Physical scale for rendering + navigation. This is NOT (yet) used for
  // collision.
  double radiusKm{12.0};

  // Docking parameters.
  // - You must be within dockingRangeKm
  // - Approach must be within dockingCorridorHalfAngleRad of the corridor axis
  // - Relative speed must be <= dockingSpeedLimitKmS
  double dockingRangeKm{60.0};
  double dockingCorridorLengthKm{3000.0};
  double dockingCorridorHalfAngleRad{0.25}; // ~14 degrees
  double dockingSpeedLimitKmS{0.10};        // ~100 m/s
};

struct StarSystem {
  SystemStub stub{};
  Star star{};
  std::vector<Planet> planets{};
  std::vector<Station> stations{};
};

} // namespace stellar::sim
