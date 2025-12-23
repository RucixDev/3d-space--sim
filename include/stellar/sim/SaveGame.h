#pragma once

#include "stellar/econ/Economy.h"
#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Celestial.h"

#include <array>
#include <string>
#include <vector>

namespace stellar::sim {

struct StationEconomyOverride {
  StationId stationId{0};
  econ::StationEconomyState state{};
};

struct SaveGame {
  int version{1};

  core::u64 seed{0};
  double timeDays{0.0};

  SystemId currentSystem{0};
  StationId dockedStation{0};

  // Player ship
  math::Vec3d shipPosKm{0,0,0};
  math::Vec3d shipVelKmS{0,0,0};
  math::Quatd shipOrient{1,0,0,0};
  math::Vec3d shipAngVelRadS{0,0,0};

  // Economy
  double credits{1000.0};
  std::array<double, econ::kCommodityCount> cargo{}; // units

  std::vector<StationEconomyOverride> stationOverrides{};
};

bool saveToFile(const SaveGame& s, const std::string& path);
bool loadFromFile(const std::string& path, SaveGame& out);

} // namespace stellar::sim
