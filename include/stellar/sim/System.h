#pragma once

#include "stellar/econ/Economy.h"
#include "stellar/sim/Celestial.h"

#include <vector>

namespace stellar::sim {

// Simple docking/approach parameters for gameplay.
// These are intentionally "game-y" and can be iterated on as we add real station meshes.
struct DockingParams {
  // Physical station radius (used for docking offset / visuals).
  double radiusKm{6.0};

  // Comms range to request docking.
  double commsRangeKm{50000.0};

  // Approach corridor (cylinder) in the station's local approach axis.
  // The axis is defined in gameplay code; currently it is the outward radial from the star.
  double corridorRadiusKm{80.0};
  double corridorLengthKm{1800.0};

  // Maximum relative speed allowed inside the corridor for docking.
  double speedLimitKmS{0.12};

  // Alignment requirement: ship must face within acos(alignCos) of the corridor axis (inbound).
  double alignCos{0.9659258263}; // cos(15 deg)
};

struct Station {
  StationId id{0};
  std::string name;
  econ::StationType type{econ::StationType::Outpost};
  core::u32 factionId{0};       // 0 = independent
  double feeRate{0.0};          // market fee rate (0..1)
  econ::StationEconomyModel economyModel{};

  // Orbital parameters (around the system's primary star).
  OrbitElements orbit{};

  // Docking/approach parameters.
  DockingParams docking{};
};

struct StarSystem {
  SystemStub stub{};
  Star star{};
  std::vector<Planet> planets{};
  std::vector<Station> stations{};
};

} // namespace stellar::sim
