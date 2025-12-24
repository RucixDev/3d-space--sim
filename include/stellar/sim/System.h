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

  // Physical placement (for in-system flight / docking).
  // Orbits the system primary at the given Keplerian elements.
  OrbitElements orbit{};

  // Approximate physical radius of the station (km). Used for docking checks.
  double radiusKm{25.0};

  // Simple docking approach "corridor" (cylinder) extending out from the station.
  // The corridor axis is defined at runtime (see game prototype) and is intended
  // to enforce safe approach speed and alignment.
  double corridorLengthKm{50.0};
  double corridorRadiusKm{15.0};

  // Speed limit inside the corridor (km/s). 0.10 km/s = 100 m/s.
  double corridorSpeedLimitKmS{0.10};

  // Alignment requirement inside corridor: ship must face the entrance within
  // a cone of this half-angle (degrees).
  double corridorAlignHalfAngleDeg{20.0};
};

struct StarSystem {
  SystemStub stub{};
  Star star{};
  std::vector<Planet> planets{};
  std::vector<Station> stations{};
};

} // namespace stellar::sim
