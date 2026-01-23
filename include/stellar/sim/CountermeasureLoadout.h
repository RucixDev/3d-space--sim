#pragma once

#include "stellar/sim/ShipLoadout.h"

#include <algorithm>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Countermeasure loadout caps (flares / chaff / heat sinks)
// -----------------------------------------------------------------------------
//
// The prototype game historically kept player countermeasure caps inside the
// UI layer. We centralize them here so gameplay services (station restock,
// loadout validation, future tools) share the same deterministic rules.
//
// NOTE: Values here implicitly define save-file and gameplay semantics.
// Keep them stable unless you also migrate saves / UI assumptions.

struct CountermeasureCaps {
  int flares{0};
  int chaff{0};
  int heatSinks{0};
};

// Capacity by hull. Larger hulls carry a bit more.
inline constexpr CountermeasureCaps countermeasureCapsForHull(ShipHullClass h) {
  switch (h) {
    case ShipHullClass::Scout:   return {6, 4, 1};
    case ShipHullClass::Hauler:  return {8, 6, 2};
    case ShipHullClass::Fighter: return {10, 8, 2};
  }
  return {6, 4, 1};
}

inline void clampCountermeasureAmmo(int& ioFlares, int& ioChaff, int& ioHeatSinks, ShipHullClass hull) {
  const auto cap = countermeasureCapsForHull(hull);
  ioFlares = std::clamp(ioFlares, 0, cap.flares);
  ioChaff = std::clamp(ioChaff, 0, cap.chaff);
  ioHeatSinks = std::clamp(ioHeatSinks, 0, cap.heatSinks);
}

inline void restockCountermeasureAmmo(int& ioFlares, int& ioChaff, int& ioHeatSinks, ShipHullClass hull) {
  const auto cap = countermeasureCapsForHull(hull);
  ioFlares = cap.flares;
  ioChaff = cap.chaff;
  ioHeatSinks = cap.heatSinks;
}

} // namespace stellar::sim
