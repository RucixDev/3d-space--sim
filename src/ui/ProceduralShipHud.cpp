#include "stellar/ui/ProceduralShipHud.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace stellar::ui {

const char* toString(ShipHudInstrument i) {
  switch (i) {
    case ShipHudInstrument::Speed: return "Speed";
    case ShipHudInstrument::Shield: return "Shield";
    case ShipHudInstrument::Hull: return "Hull";
    case ShipHudInstrument::Heat: return "Heat";
    case ShipHudInstrument::Fuel: return "Fuel";
    case ShipHudInstrument::Pips: return "Pips";
    case ShipHudInstrument::Fsd: return "FSD";
    case ShipHudInstrument::Throttle: return "Throttle";
    case ShipHudInstrument::Target: return "Target";
    case ShipHudInstrument::Cargo: return "Cargo";
    case ShipHudInstrument::GForce: return "G";
    case ShipHudInstrument::Attitude: return "Attitude";
    default: return "Unknown";
  }
}

const char* toString(ShipHudGaugeStyle s) {
  switch (s) {
    case ShipHudGaugeStyle::Bar: return "Bar";
    case ShipHudGaugeStyle::Dial: return "Dial";
    case ShipHudGaugeStyle::Arc: return "Arc";
    case ShipHudGaugeStyle::Segmented: return "Segmented";
    default: return "Unknown";
  }
}

namespace {

// Normalized unit-square packing rectangle.
struct Rect {
  double x0{0.0}, y0{0.0}, x1{1.0}, y1{1.0};
  double w() const { return x1 - x0; }
  double h() const { return y1 - y0; }
  double area() const { return std::max(0.0, w()) * std::max(0.0, h()); }
};

struct Item {
  ShipHudInstrument inst{ShipHudInstrument::Speed};
  double weight{1.0};
  ShipHudGaugeStyle style{ShipHudGaugeStyle::Bar};
  bool spark{false};
  core::u8 accent{0};
};

struct AreaItem {
  Item it{};
  double area{0.0}; // normalized area in [0..1], sum(items.area) ~= 1
};

static double importance(ShipHudInstrument inst) {
  // Higher = larger panel area in the treemap.
  // Keep this stable so seeds preserve their "feel" across versions.
  switch (inst) {
    case ShipHudInstrument::Speed: return 3.00;
    case ShipHudInstrument::Shield: return 2.65;
    case ShipHudInstrument::Attitude: return 2.45;
    case ShipHudInstrument::Hull: return 2.25;
    case ShipHudInstrument::Fuel: return 1.95;
    case ShipHudInstrument::Heat: return 1.85;
    case ShipHudInstrument::Pips: return 1.65;
    case ShipHudInstrument::Fsd: return 1.55;
    case ShipHudInstrument::Throttle: return 1.45;
    case ShipHudInstrument::Target: return 1.35;
    case ShipHudInstrument::Cargo: return 1.15;
    case ShipHudInstrument::GForce: return 1.05;
    default: return 1.0;
  }
}

static bool defaultSparkline(ShipHudInstrument inst) {
  switch (inst) {
    case ShipHudInstrument::Speed:
    case ShipHudInstrument::Shield:
    case ShipHudInstrument::Hull:
    case ShipHudInstrument::Heat:
    case ShipHudInstrument::Fuel:
      return true;
    default:
      return false;
  }
}

static ShipHudGaugeStyle chooseStyle(ShipHudInstrument inst, core::SplitMix64& rng) {
  auto pick = [&](ShipHudGaugeStyle a, ShipHudGaugeStyle b) {
    return rng.chance(0.5) ? a : b;
  };

  switch (inst) {
    case ShipHudInstrument::Speed: return pick(ShipHudGaugeStyle::Dial, ShipHudGaugeStyle::Arc);
    case ShipHudInstrument::Shield: return pick(ShipHudGaugeStyle::Arc, ShipHudGaugeStyle::Segmented);
    case ShipHudInstrument::Hull: return pick(ShipHudGaugeStyle::Bar, ShipHudGaugeStyle::Segmented);
    case ShipHudInstrument::Heat: return ShipHudGaugeStyle::Arc;
    case ShipHudInstrument::Fuel: return ShipHudGaugeStyle::Bar;
    case ShipHudInstrument::Pips: return ShipHudGaugeStyle::Segmented;
    case ShipHudInstrument::Fsd: return pick(ShipHudGaugeStyle::Bar, ShipHudGaugeStyle::Segmented);
    case ShipHudInstrument::Throttle: return ShipHudGaugeStyle::Bar;
    case ShipHudInstrument::Target: return pick(ShipHudGaugeStyle::Bar, ShipHudGaugeStyle::Segmented);
    case ShipHudInstrument::Cargo: return ShipHudGaugeStyle::Bar;
    case ShipHudInstrument::GForce: return pick(ShipHudGaugeStyle::Dial, ShipHudGaugeStyle::Bar);
    case ShipHudInstrument::Attitude: return ShipHudGaugeStyle::Dial; // rendered as a custom instrument
    default: return ShipHudGaugeStyle::Bar;
  }
}

static double sumArea(const std::vector<AreaItem>& row) {
  double s = 0.0;
  for (const auto& it : row) s += it.area;
  return s;
}

static double worstAspect(const std::vector<AreaItem>& row, double shortSide) {
  if (row.empty()) return std::numeric_limits<double>::infinity();
  shortSide = std::max(1e-9, shortSide);

  double s = 0.0;
  double mn = std::numeric_limits<double>::infinity();
  double mx = 0.0;
  for (const auto& it : row) {
    s += it.area;
    mn = std::min(mn, it.area);
    mx = std::max(mx, it.area);
  }
  s = std::max(1e-12, s);
  mn = std::max(1e-12, mn);
  mx = std::max(1e-12, mx);

  const double w2 = shortSide * shortSide;
  const double s2 = s * s;

  const double r1 = (w2 * mx) / s2;
  const double r2 = s2 / (w2 * mn);
  return std::max(r1, r2);
}

static void pushPanel(std::vector<ShipHudPanel>& out, const Item& it, const Rect& r) {
  ShipHudPanel p{};
  p.instrument = it.inst;
  p.style = it.style;
  p.sparkline = it.spark;
  p.accentVariant = it.accent;

  p.x0 = (float)r.x0;
  p.y0 = (float)r.y0;
  p.x1 = (float)r.x1;
  p.y1 = (float)r.y1;

  out.push_back(p);
}

// Squarified treemap row layout.
// This yields much nicer aspect ratios than the old greedy BSP split, which helps
// readability (especially for circular instruments like Attitude).
static void layoutRowAndAdvance(const std::vector<AreaItem>& row, Rect& r, core::SplitMix64& rng, std::vector<ShipHudPanel>& out) {
  if (row.empty()) return;
  const double rw = std::max(1e-9, r.w());
  const double rh = std::max(1e-9, r.h());
  const double aRow = std::max(0.0, sumArea(row));
  if (aRow <= 1e-12) return;

  // Orientation: lay rows along the shorter edge.
  const bool horizontal = (rw >= rh);

  if (horizontal) {
    const double hRow = std::clamp(aRow / rw, 0.0, rh);
    const bool placeTop = rng.chance(0.5);
    const double y0 = placeTop ? r.y0 : (r.y1 - hRow);
    const double y1 = y0 + hRow;

    bool leftToRight = rng.chance(0.5);
    double x = leftToRight ? r.x0 : r.x1;

    // Optionally reverse item order within the row for more variety.
    const bool reverse = rng.chance(0.45);
    const int n = (int)row.size();
    for (int i = 0; i < n; ++i) {
      const int idx = reverse ? (n - 1 - i) : i;
      const AreaItem& ai = row[(std::size_t)idx];
      const double wItem = std::clamp(ai.area / std::max(1e-9, hRow), 0.0, rw);

      Rect pr{};
      pr.y0 = y0;
      pr.y1 = y1;
      if (leftToRight) {
        pr.x0 = x;
        pr.x1 = x + wItem;
        x = pr.x1;
      } else {
        pr.x1 = x;
        pr.x0 = x - wItem;
        x = pr.x0;
      }

      // Clamp to be safe.
      pr.x0 = std::clamp(pr.x0, r.x0, r.x1);
      pr.x1 = std::clamp(pr.x1, r.x0, r.x1);

      pushPanel(out, ai.it, pr);
    }

    if (placeTop) r.y0 = std::min(r.y0 + hRow, r.y1);
    else r.y1 = std::max(r.y1 - hRow, r.y0);

  } else {
    const double wRow = std::clamp(aRow / rh, 0.0, rw);
    const bool placeLeft = rng.chance(0.5);
    const double x0 = placeLeft ? r.x0 : (r.x1 - wRow);
    const double x1 = x0 + wRow;

    bool topToBottom = rng.chance(0.5);
    double y = topToBottom ? r.y0 : r.y1;

    const bool reverse = rng.chance(0.45);
    const int n = (int)row.size();
    for (int i = 0; i < n; ++i) {
      const int idx = reverse ? (n - 1 - i) : i;
      const AreaItem& ai = row[(std::size_t)idx];
      const double hItem = std::clamp(ai.area / std::max(1e-9, wRow), 0.0, rh);

      Rect pr{};
      pr.x0 = x0;
      pr.x1 = x1;
      if (topToBottom) {
        pr.y0 = y;
        pr.y1 = y + hItem;
        y = pr.y1;
      } else {
        pr.y1 = y;
        pr.y0 = y - hItem;
        y = pr.y0;
      }

      pr.y0 = std::clamp(pr.y0, r.y0, r.y1);
      pr.y1 = std::clamp(pr.y1, r.y0, r.y1);

      pushPanel(out, ai.it, pr);
    }

    if (placeLeft) r.x0 = std::min(r.x0 + wRow, r.x1);
    else r.x1 = std::max(r.x1 - wRow, r.x0);
  }
}

static void squarifyTreemap(const std::vector<AreaItem>& items, core::SplitMix64& rng, std::vector<ShipHudPanel>& out) {
  out.clear();
  out.reserve(items.size());

  Rect r{};
  std::vector<AreaItem> row;
  row.reserve(12);

  const std::size_t n = items.size();
  std::size_t i = 0;
  while (i < n) {
    const AreaItem& next = items[i];

    const double shortSide = std::max(1e-9, std::min(r.w(), r.h()));

    if (row.empty()) {
      row.push_back(next);
      ++i;
      continue;
    }

    // Compare the worst aspect ratio with/without 'next' to decide whether to append.
    std::vector<AreaItem> rowPlus = row;
    rowPlus.push_back(next);

    const double w0 = worstAspect(row, shortSide);
    const double w1 = worstAspect(rowPlus, shortSide);

    if (w1 <= w0) {
      row.push_back(next);
      ++i;
    } else {
      layoutRowAndAdvance(row, r, rng, out);
      row.clear();
    }
  }

  if (!row.empty()) {
    layoutRowAndAdvance(row, r, rng, out);
  }

  // Final clamp and ordering stability.
  for (ShipHudPanel& p : out) {
    p.x0 = std::clamp(p.x0, 0.0f, 1.0f);
    p.y0 = std::clamp(p.y0, 0.0f, 1.0f);
    p.x1 = std::clamp(p.x1, 0.0f, 1.0f);
    p.y1 = std::clamp(p.y1, 0.0f, 1.0f);
    if (p.x1 < p.x0) std::swap(p.x0, p.x1);
    if (p.y1 < p.y0) std::swap(p.y0, p.y1);
  }
}

} // namespace

ProceduralShipHudPlan makeProceduralShipHudPlan(core::u64 seed, int detailLevel) {
  ProceduralShipHudPlan plan{};
  plan.seed = seed;
  plan.detailLevel = std::clamp(detailLevel, 0, 3);

  core::SplitMix64 rng(core::hashCombine(seed, core::fnv1a64("ship_hud_plan")));

  // A stable visual theme selector for renderers. This is derived from the same
  // plan RNG so it changes when the user rerolls the layout seed.
  plan.themeVariant = (core::u8)rng.range<int>(0, 3);

  // Base size scales with detail. With the squarified treemap layout we can
  // pack panels more efficiently, but we still grow slightly with instrument count.
  plan.baseWidthPx = (340.0f + (float)plan.detailLevel * 60.0f) + rng.range<float>(-8.0f, 8.0f);
  plan.baseHeightPx = (210.0f + (float)plan.detailLevel * 44.0f) + rng.range<float>(-6.0f, 6.0f);
  plan.cellPaddingNorm = rng.range<float>(0.010f, 0.016f);

  std::vector<ShipHudInstrument> insts;
  insts.reserve(14);

  // Always-present core instruments.
  insts.push_back(ShipHudInstrument::Speed);
  insts.push_back(ShipHudInstrument::Shield);
  insts.push_back(ShipHudInstrument::Hull);
  insts.push_back(ShipHudInstrument::Fuel);
  insts.push_back(ShipHudInstrument::Heat);

  if (plan.detailLevel >= 1) {
    insts.push_back(ShipHudInstrument::Pips);
    insts.push_back(ShipHudInstrument::Fsd);
  }
  if (plan.detailLevel >= 2) {
    insts.push_back(ShipHudInstrument::Attitude);
    insts.push_back(ShipHudInstrument::Throttle);
    insts.push_back(ShipHudInstrument::Target);
  }
  if (plan.detailLevel >= 3) {
    insts.push_back(ShipHudInstrument::Cargo);
    insts.push_back(ShipHudInstrument::GForce);
  }

  std::vector<Item> items;
  items.reserve(insts.size());

  const core::u8 accentBase = (core::u8)rng.range<int>(0, 3);
  for (std::size_t idx = 0; idx < insts.size(); ++idx) {
    const ShipHudInstrument inst = insts[idx];
    Item it{};
    it.inst = inst;

    // Tiny per-instrument jitter keeps layouts feeling organic while staying seed-stable.
    it.weight = std::max(0.18, importance(inst) + (double)rng.range<float>(-0.18f, 0.18f));

    // Attitude is a special renderer path; keep its "style" deterministic anyway.
    it.style = chooseStyle(inst, rng);
    it.spark = defaultSparkline(inst);

    // Accent variant is a small renderer hint.
    it.accent = (core::u8)((accentBase + (core::u8)idx + (core::u8)rng.range<int>(0, 2)) & 3u);
    items.push_back(it);
  }

  // Sort by descending weight so the treemap's row-building is stable and tends to reserve
  // good aspect ratio real estate for important instruments.
  std::stable_sort(items.begin(), items.end(), [&](const Item& a, const Item& b) {
    if (a.weight == b.weight) return (int)a.inst < (int)b.inst;
    return a.weight > b.weight;
  });

  // Convert weights to normalized areas.
  double sumW = 0.0;
  for (const Item& it : items) sumW += std::max(0.0, it.weight);
  if (sumW <= 1e-9) sumW = 1.0;

  std::vector<AreaItem> aitems;
  aitems.reserve(items.size());
  for (const Item& it : items) {
    AreaItem ai{};
    ai.it = it;
    ai.area = std::max(1e-8, it.weight / sumW);
    aitems.push_back(ai);
  }

  // Normalize precisely so the unit square is covered (within floating tolerance).
  double sumA = 0.0;
  for (const auto& ai : aitems) sumA += ai.area;
  if (sumA > 1e-12) {
    for (auto& ai : aitems) ai.area /= sumA;
  }

  // Generate panels.
  squarifyTreemap(aitems, rng, plan.panels);

  return plan;
}

} // namespace stellar::ui
