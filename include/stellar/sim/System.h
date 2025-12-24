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

  // Physical/orbital placement (around the primary star).
  // NOTE: this is intentionally light-weight / "good enough" for prototype gameplay.
  // Stations are placed on their own Keplerian orbits so the player can fly to them and dock.
  OrbitElements orbit{};

  // Approximate physical size (km). Used for docking range hints / soft collision.
  double radiusKm{12.0};
};

struct StarSystem {
  SystemStub stub{};
  Star star{};
  std::vector<Planet> planets{};
  std::vector<Station> stations{};
};

} // namespace stellar::sim
