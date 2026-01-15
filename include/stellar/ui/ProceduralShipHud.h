#pragma once

#include "stellar/core/Types.h"

#include <vector>

namespace stellar::ui {

// Deterministic procedural "instrument cluster" layout for the in-flight Ship HUD.
//
// This module is intentionally ImGui-free so it can be generated/validated in tests
// and reused by future tooling.

enum class ShipHudInstrument : core::u8 {
  Speed = 0,
  Shield,
  Hull,
  Heat,
  Fuel,
  Pips,
  Fsd,
  Throttle,
  Target,
  Cargo,
  GForce,

  // Attitude / reference frame instrument (horizon + vectors).
  Attitude,

  Count
};

enum class ShipHudGaugeStyle : core::u8 {
  Bar = 0,
  Dial,
  Arc,
  Segmented
};

struct ShipHudPanel {
  ShipHudInstrument instrument{ShipHudInstrument::Speed};
  ShipHudGaugeStyle style{ShipHudGaugeStyle::Bar};

  // Normalized panel rectangle within the Ship HUD client area (0..1).
  float x0{0.0f};
  float y0{0.0f};
  float x1{1.0f};
  float y1{1.0f};

  // Whether this instrument benefits from a small telemetry sparkline.
  bool sparkline{false};

  // Small integer used by the renderer to pick subtle styling variations.
  core::u8 accentVariant{0};
};

struct ProceduralShipHudPlan {
  core::u64 seed{0};
  int detailLevel{2}; // 0=minimal .. 3=max

  // Stable, seed-derived theme variant for renderers.
  //
  // NOTE: This is intentionally *not* a color palette. Colors are sourced from
  // the global HUD palette (HudSettings). The theme variant represents
  // structural styling such as corner rounding, tick density, and background
  // decor patterns.
  core::u8 themeVariant{0};

  // Suggested base size before applying per-widget scale.
  float baseWidthPx{420.0f};
  float baseHeightPx{260.0f};

  // Normalized padding inset applied to each panel by renderers.
  float cellPaddingNorm{0.012f};

  std::vector<ShipHudPanel> panels;
};

const char* toString(ShipHudInstrument i);
const char* toString(ShipHudGaugeStyle s);

// Generate a deterministic plan. Same seed+detail => identical plan.
ProceduralShipHudPlan makeProceduralShipHudPlan(core::u64 seed, int detailLevel);

} // namespace stellar::ui
