#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"

#include <string>
#include <vector>

namespace stellar::sim {

struct Faction {
  core::u32 id{0};              // 0 reserved for "Independent"
  std::string name{"Independent"};
  math::Vec3d homePosLy{0,0,0}; // galaxy-space
  math::Vec3d color{0.8,0.8,0.8};
  double influenceRadiusLy{0.0}; // how far this faction tends to own stations
  double taxRate{0.02};          // market fee baseline
  double industryBias{0.0};      // -1..+1 (agri <-> industrial)
};

std::vector<Faction> generateFactions(core::u64 seed, int count);

} // namespace stellar::sim
