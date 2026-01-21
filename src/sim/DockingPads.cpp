#include "stellar/sim/DockingPads.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/core/Clamp.h"

#include "stellar/math/Vec2.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace stellar::sim {

namespace {

static core::u64 padsSeed(core::u64 universeSeed, StationId stId) {
  // Use a domain-separated deterministic seed so pad layouts are stable and
  // independent from other procedural systems.
  const core::u64 salt = core::fnv1a64("DockingPads.v1");
  return core::hashCombine(core::hashCombine(universeSeed, salt), (core::u64)stId);
}

static void generatePads2D(core::SplitMix64& rng,
                           int count,
                           const math::Vec2d& halfExt,
                           std::vector<math::Vec2d>& outPts) {
  outPts.clear();
  outPts.reserve((std::size_t)std::max(0, count));

  auto randPt = [&]() -> math::Vec2d {
    return {
      rng.range(-halfExt.x, halfExt.x),
      rng.range(-halfExt.y, halfExt.y)
    };
  };

  if (count <= 0) return;

  // Farthest-point sampling ("best candidate") gives a good spread while staying
  // deterministic and lightweight.
  constexpr int kCandidates = 128;

  outPts.push_back(randPt());

  for (int i = 1; i < count; ++i) {
    math::Vec2d best = randPt();
    double bestMinD2 = -1.0;

    for (int c = 0; c < kCandidates; ++c) {
      const math::Vec2d cand = randPt();

      double minD2 = std::numeric_limits<double>::infinity();
      for (const auto& p : outPts) {
        const double dx = cand.x - p.x;
        const double dy = cand.y - p.y;
        const double d2 = dx * dx + dy * dy;
        if (d2 < minD2) minD2 = d2;
      }

      if (minD2 > bestMinD2) {
        bestMinD2 = minD2;
        best = cand;
      }
    }

    outPts.push_back(best);
  }
}

} // namespace

int stationDockingPadCount(const Station& st) {
  // Default a "reasonable" radius if the caller didn't populate one (common in tests).
  const double radiusKm = (st.radiusKm > 1.0) ? st.radiusKm : 6000.0;

  // Base counts tuned for gameplay feel, not realism.
  int base = 12;
  switch (st.type) {
    case econ::StationType::Outpost:     base = 8;  break;
    case econ::StationType::Mining:      base = 12; break;
    case econ::StationType::Refinery:    base = 14; break;
    case econ::StationType::Industrial:  base = 16; break;
    case econ::StationType::Research:    base = 12; break;
    case econ::StationType::TradeHub:    base = 22; break;
    case econ::StationType::Shipyard:    base = 28; break;
    default:                              base = 12; break;
  }

  // Scale gently with radius.
  const double rNorm = std::clamp(radiusKm / 6000.0, 0.35, 2.8);
  const double scale = std::pow(rNorm, 1.10);

  const int count = (int)std::lround((double)base * scale);
  return core::clamp(count, 4, 48);
}

DockingPad stationDockingPad(core::u64 universeSeed, const Station& st, core::u16 index1Based) {
  DockingPad pad{};

  const int padCount = stationDockingPadCount(st);
  if (padCount <= 0) {
    pad.index = 1;
    pad.localPosKm = {0,0,0};
    pad.localOrient = math::Quatd::identity();
    pad.radiusKm = 250.0;
    return pad;
  }

  const int idx0 = (int)core::clamp<int>((int)index1Based - 1, 0, padCount - 1);

  core::SplitMix64 rng(padsSeed(universeSeed, st.id));

  // Hangar cross-section region where pads are placed.
  // We keep pads well inside the slot bounds so the resulting docked positions
  // do not appear "clipped" into the tunnel wall.
  const double slotW = (st.slotWidthKm > 1.0) ? st.slotWidthKm : (st.radiusKm > 1.0 ? st.radiusKm * 0.85 : 5000.0);
  const double slotH = (st.slotHeightKm > 1.0) ? st.slotHeightKm : (st.radiusKm > 1.0 ? st.radiusKm * 0.75 : 4000.0);

  math::Vec2d halfExt{ slotW * 0.5 * 0.78, slotH * 0.5 * 0.78 };
  halfExt.x = std::max(halfExt.x, 120.0);
  halfExt.y = std::max(halfExt.y, 120.0);

  std::vector<math::Vec2d> pts;
  generatePads2D(rng, padCount, halfExt, pts);

  // Deterministic depth inside the station.
  const double radiusKm = (st.radiusKm > 1.0) ? st.radiusKm : 6000.0;
  const double wz = radiusKm * 1.10;
  const double dockZ = wz - st.slotDepthKm - radiusKm * 0.25;
  const double layerSpacing = std::max(160.0, radiusKm * 0.035);
  const int layer = (idx0 % 3) - 1; // -1,0,1
  double z = dockZ - (double)layer * layerSpacing;

  // Keep within the station hull bounds (box approximation used by collision code).
  const double zMax = wz * 0.95;
  z = std::clamp(z, -zMax, zMax);

  pad.index = (core::u16)(idx0 + 1);
  if (idx0 < (int)pts.size()) {
    pad.localPosKm = { pts[(std::size_t)idx0].x, pts[(std::size_t)idx0].y, z };
  } else {
    pad.localPosKm = { 0, 0, z };
  }

  // By default, a docked ship faces +Z in station-local space (outward from the slot).
  pad.localOrient = math::Quatd::identity();

  // Approximate pad radius as a fraction of the hangar width.
  const double r = std::max(160.0, std::min(halfExt.x, halfExt.y) * 0.18);
  pad.radiusKm = r;

  return pad;
}

core::u16 nearestDockingPadIndex(core::u64 universeSeed, const Station& st, const math::Vec3d& relLocalKm) {
  const int padCount = stationDockingPadCount(st);
  if (padCount <= 0) return 1;

  core::u16 bestIdx = 1;
  double bestD2 = std::numeric_limits<double>::infinity();

  for (int i = 1; i <= padCount; ++i) {
    const DockingPad p = stationDockingPad(universeSeed, st, (core::u16)i);
    const math::Vec3d d = relLocalKm - p.localPosKm;
    const double d2 = d.lengthSq();
    if (d2 < bestD2) {
      bestD2 = d2;
      bestIdx = (core::u16)i;
    }
  }

  return bestIdx;
}

} // namespace stellar::sim
