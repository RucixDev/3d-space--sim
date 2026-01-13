#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/proc/GalaxyGenerator.h"
#include "stellar/sim/Celestial.h"
#include "stellar/sim/Faction.h"

#include <functional>
#include <string>
#include <vector>

namespace stellar::game {

struct ProceduralGalaxyLabWindowState {
  bool open{false};

  core::u64 seed{1337ull};
  int factionCount{8};

  proc::GalaxyParams params{};

  // View / preview controls (top-down XY projection).
  math::Vec3d centerLy{0.0, 0.0, 0.0};
  double viewRadiusLy{250.0};
  double zHalfLy{200.0};
  std::size_t maxStubs{20000};

  bool autoRegenerate{true};
  bool colorByFaction{false};
  bool showArmGuides{true};
  bool showLegend{true};

  // Cached preview data.
  bool dirty{true};
  std::vector<sim::Faction> factions{};
  std::vector<sim::SystemStub> stubs{};
  double lastGenMs{0.0};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

void drawProceduralGalaxyLabWindow(ProceduralGalaxyLabWindowState& state, float timeSec, const ToastFn& toast);

} // namespace stellar::game
