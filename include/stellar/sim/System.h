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

  // Physical orbit around the system primary (origin).
  // Units: AU for orbit elements, days for period/epoch.
  OrbitElements orbit{};

  // Visual scale hint (km). Gameplay currently treats stations as points, but this
  // can be used for rendering and safe-zone logic.
  double radiusKm{5.0};

  // Docking/traffic control parameters.
  struct Docking {
    // Range in which docking clearance can be requested.
    double commsRangeKm{25.0};

    // Approach corridor: cylinder from the docking port outward.
    // Corridor axis points away from the star (radial) at the station's position.
    double corridorLengthKm{50.0};
    double corridorRadiusKm{8.0};

    // Must be within this distance (km) to complete docking.
    double dockRangeKm{1.5};

    // Relative speed cap (km/s) to dock.
    double speedLimitKmS{0.25};

    // Minimum alignment with docking axis to dock. (1.0 = perfect)
    double alignCosMin{0.95};

    // Clearance timer (minutes) once granted.
    double clearanceDurationMin{10.0};

    // Cooldown after a denial (seconds).
    double deniedCooldownSec{60.0};
  } docking{};
};

struct StarSystem {
  SystemStub stub{};
  Star star{};
  std::vector<Planet> planets{};
  std::vector<Station> stations{};
};

} // namespace stellar::sim
