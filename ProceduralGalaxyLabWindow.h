#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/proc/GalaxyGenerator.h"
#include "stellar/proc/GalaxyRegions.h"
#include "stellar/proc/Hyperlanes.h"
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
  bool colorByRegion{false};
  bool showArmGuides{true};
  bool showLegend{true};
  // Hyperlane overlay (procedural sparse navigation/trade corridor graph).
  bool showHyperlanes{false};
  proc::HyperlaneParams hyperlaneParams{};


  // Galaxy regions (Worley/Voronoi). Used for visualization only.
  double regionCellSizeLy{900.0};

  // Cached preview data.
  bool dirty{true};
  std::vector<sim::Faction> factions{};
  std::vector<sim::SystemStub> stubs{};

  // Cached per-stub region kind/id for preview rendering.
  std::vector<proc::GalaxyRegionKind> stubRegionKind{};
  std::vector<core::u64> stubRegionId{};
  // Cached hyperlane overlay.
  std::vector<proc::HyperlaneEdge> hyperlanes{};
  double lastLaneMs{0.0};

  double lastGenMs{0.0};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

void drawProceduralGalaxyLabWindow(ProceduralGalaxyLabWindowState& state, float timeSec, const ToastFn& toast);

} // namespace stellar::game
