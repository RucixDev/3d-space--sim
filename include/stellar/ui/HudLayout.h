#pragma once

#include "stellar/core/Types.h"

#include <array>
#include <optional>
#include <string>

namespace stellar::ui {

// Identifiers for HUD overlays that can be positioned/scaled by the layout system.
//
// Note: this is intentionally small and focused on in-flight overlays that are drawn
// with hard-coded positions. "Regular" ImGui windows already have built-in persistence
// via imgui.ini.
enum class HudWidgetId : core::u8 {
  Radar = 0,
  Objective,
  Threat,
  Jump,

  Count
};

struct HudWidgetLayout {
  // Normalized position (0..1) within the main viewport's WorkRect.
  // This represents the location of the widget's pivot point.
  float posNormX{0.5f};
  float posNormY{0.5f};

  // Pivot point inside the widget rect (0..1). (0,0)=top-left, (1,1)=bottom-right.
  float pivotX{0.5f};
  float pivotY{0.5f};

  // Per-widget scale multiplier. Interpreted by the UI code.
  float scale{1.0f};

  // Optional visibility switch for widgets that have a persistent toggle.
  // (Example: Radar, Objective HUD, Jump overlay HUD.)
  bool enabled{true};
};

struct HudLayout {
  int version{1};

  // Cosmetic "safe area" in pixels (used by the editor as a guide).
  float safeMarginPx{14.0f};

  std::array<HudWidgetLayout, (std::size_t)HudWidgetId::Count> widgets{};

  HudWidgetLayout& widget(HudWidgetId id) { return widgets[(std::size_t)id]; }
  const HudWidgetLayout& widget(HudWidgetId id) const { return widgets[(std::size_t)id]; }
};

// Default config file path (relative to the working directory).
std::string defaultHudLayoutPath();

// Default layout used when no file exists or parsing fails.
HudLayout makeDefaultHudLayout();

// Human-readable widget name (used in the layout file + UI).
const char* toString(HudWidgetId id);

// Parse widget name (case-insensitive, supports a few aliases).
std::optional<HudWidgetId> widgetFromString(const std::string& s);

bool saveToFile(const HudLayout& layout, const std::string& path);
bool loadFromFile(const std::string& path, HudLayout& out);

} // namespace stellar::ui
