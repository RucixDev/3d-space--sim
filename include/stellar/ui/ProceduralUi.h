#pragma once

#include "stellar/core/Types.h"

#include "stellar/ui/HudLayout.h"
#include "stellar/ui/HudSettings.h"

namespace stellar::ui {

// A small, deterministic "skin" for HUD overlays (palette + cosmetics), derived from a seed.
//
// This is intentionally ImGui-free so it can be used by headless tests and future tooling.
struct ProceduralHudSkin {
  float overlayBgAlpha{0.35f};
  float overlayBgAlphaEdit{0.45f};
  bool tintRadarIcons{false};
  bool tintTacticalIcons{false};

  Color4f colorPrimary{0.82f, 0.90f, 1.00f, 1.00f};
  Color4f colorAccent{1.00f, 0.70f, 0.35f, 1.00f};
  Color4f colorDanger{1.00f, 0.25f, 0.25f, 1.00f};
  Color4f colorGrid{0.47f, 0.55f, 0.67f, 1.00f};
  Color4f colorText{0.90f, 0.90f, 0.94f, 1.00f};
  Color4f colorBackground{0.00f, 0.00f, 0.00f, 1.00f};
};

// Generate a HUD skin from a seed. The same seed always yields the same skin.
ProceduralHudSkin makeProceduralHudSkin(core::u64 seed);

// Apply a generated skin to an existing HudSettings object (preserves non-style fields).
void applyHudSkinToSettings(const ProceduralHudSkin& skin, HudSettings& inOut);

// Apply a deterministic procedural HUD layout to an existing layout.
//
// Design intent:
// - preserve per-widget scale/enable flags (users often tune those manually)
// - vary the *placement* in a safe, predictable way (mirrors/corners, small jitters)
void applyProceduralHudLayout(core::u64 seed, HudLayout& inOut);

// Convenience helper: generate a new layout from defaults.
HudLayout makeProceduralHudLayout(core::u64 seed);

// Convenience helper: derive an RGB accent color from a seed.
// (Currently matches makeProceduralHudSkin(seed).colorPrimary).
void makeProceduralAccentRgb(core::u64 seed, float outRgb[3]);

} // namespace stellar::ui
