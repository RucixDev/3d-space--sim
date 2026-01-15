#include "stellar/ui/ProceduralShipHud.h"

#include <algorithm>
#include <cmath>
#include <vector>

using namespace stellar;

static float area(const ui::ShipHudPanel& p) {
  return std::max(0.0f, p.x1 - p.x0) * std::max(0.0f, p.y1 - p.y0);
}

static float overlapArea(const ui::ShipHudPanel& a, const ui::ShipHudPanel& b) {
  const float x0 = std::max(a.x0, b.x0);
  const float y0 = std::max(a.y0, b.y0);
  const float x1 = std::min(a.x1, b.x1);
  const float y1 = std::min(a.y1, b.y1);
  return std::max(0.0f, x1 - x0) * std::max(0.0f, y1 - y0);
}

int test_ship_hud_plan() {
  // Determinism: same inputs -> identical plan.
  const core::u64 seed = 0xC0FFEEu;
  const int detail = 3;
  const ui::ProceduralShipHudPlan a = ui::makeProceduralShipHudPlan(seed, detail);
  const ui::ProceduralShipHudPlan b = ui::makeProceduralShipHudPlan(seed, detail);

  if (a.seed != b.seed) return 1;
  if (a.detailLevel != b.detailLevel) return 2;
  if (a.themeVariant != b.themeVariant) return 3;
  if (a.panels.size() != b.panels.size()) return 4;

  for (std::size_t i = 0; i < a.panels.size(); ++i) {
    const ui::ShipHudPanel& pa = a.panels[i];
    const ui::ShipHudPanel& pb = b.panels[i];
    if (pa.instrument != pb.instrument) return 10;
    if (pa.style != pb.style) return 11;
    if (std::abs(pa.x0 - pb.x0) > 1e-6f) return 12;
    if (std::abs(pa.y0 - pb.y0) > 1e-6f) return 13;
    if (std::abs(pa.x1 - pb.x1) > 1e-6f) return 14;
    if (std::abs(pa.y1 - pb.y1) > 1e-6f) return 15;
    if (pa.sparkline != pb.sparkline) return 16;
  }

  // Instrument set: Attitude appears at detail>=2 (and not at detail<2).
  {
    bool hasAtt = false;
    for (const ui::ShipHudPanel& p : a.panels) {
      if (p.instrument == ui::ShipHudInstrument::Attitude) hasAtt = true;
    }
    if (detail >= 2 && !hasAtt) return 17;
    if (detail < 2 && hasAtt) return 18;
  }

  // Sanity: panels are within [0,1], have positive area, and do not overlap.
  float sumArea = 0.0f;
  for (const ui::ShipHudPanel& p : a.panels) {
    if (!(p.x0 >= 0.0f && p.y0 >= 0.0f && p.x1 <= 1.0f && p.y1 <= 1.0f)) return 20;
    if (area(p) <= 1e-4f) return 21;
    sumArea += area(p);
  }

  // Pairwise overlap should be ~0.
  for (std::size_t i = 0; i < a.panels.size(); ++i) {
    for (std::size_t j = i + 1; j < a.panels.size(); ++j) {
      if (overlapArea(a.panels[i], a.panels[j]) > 1e-4f) return 22;
    }
  }

  // The treemap packing should cover the full unit square (within tolerance).
  if (std::abs(sumArea - 1.0f) > 1e-2f) return 23;

  // Detail levels should produce expected instrument counts.
  const int n = (int)a.panels.size();
  if (detail == 0 && n < 5) return 30;
  if (detail == 3 && n < 10) return 31;
  return 0;
}
