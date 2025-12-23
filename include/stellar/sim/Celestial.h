#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"

#include <string>
#include <vector>

namespace stellar::sim {

using SystemId  = core::u64;
using StationId = core::u64;

enum class StarClass : core::u8 {
  O, B, A, F, G, K, M,
  Count
};

enum class PlanetType : core::u8 {
  Rocky,
  Desert,
  Ocean,
  Ice,
  GasGiant,
  Count
};

struct Star {
  StarClass cls{StarClass::G};
  double massSol{1.0};
  double radiusSol{1.0};
  double luminositySol{1.0};
  double temperatureK{5778.0};
};

struct OrbitElements {
  // Units:
  //  - semiMajorAxisAU: AU
  //  - angles: radians
  //  - epochDays: days
  //  - periodDays: days
  double semiMajorAxisAU{1.0};
  double eccentricity{0.0};
  double inclinationRad{0.0};
  double ascendingNodeRad{0.0};
  double argPeriapsisRad{0.0};
  double meanAnomalyAtEpochRad{0.0};
  double epochDays{0.0};
  double periodDays{365.25};
};

struct Planet {
  std::string name;
  PlanetType type{PlanetType::Rocky};
  double radiusEarth{1.0};
  double massEarth{1.0};
  OrbitElements orbit{};
};

struct SystemStub {
  SystemId id{0};
  core::u64 seed{0};
  std::string name;
  math::Vec3d posLy{0,0,0}; // galaxy-space (light-years)
  StarClass primaryClass{StarClass::G};
  int planetCount{0};
  int stationCount{0};
  core::u32 factionId{0}; // 0 = independent
};

} // namespace stellar::sim
