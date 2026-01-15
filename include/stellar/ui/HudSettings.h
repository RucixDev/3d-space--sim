#pragma once

#include <string>
#include <cstdint>

namespace stellar::ui {

// Lightweight, dependency-free color type used by HUD settings (floats in [0,1]).
//
// This is intentionally separate from Dear ImGui types so HudSettings remains usable
// in headless builds and tests.
struct Color4f {
  float r{1.0f};
  float g{1.0f};
  float b{1.0f};
  float a{1.0f};
};

// Persistent configuration for in-flight HUD and tactical overlays.
//
// This intentionally does NOT depend on Dear ImGui or SDL. It is a small, human-editable
// text file that stores player preferences such as radar range, symbology toggles, and
// tactical overlay filters.
//
// Saved separately from:
//  - imgui.ini (window positions / docking)
//  - hud_layout.txt (overlay widget positions / per-widget scale)

struct HudSettings {
  // File format version.
  // 1: core toggles/tuning (radar/combat/tactical)
  // 2: adds HUD palette + overlay cosmetics (colors, background alpha, tint flags)
  // 3: adds procedural Ship HUD settings
  // 4: adds procedural Ship HUD cosmetics (decor/glyphs/microtext/dropouts)
  // 5: adds Ship HUD segment-display readouts (procedural vector text)
  int version{5};

  // Persist this file on exit.
  bool autoSaveOnExit{true};

  // --- Primary overlays (these map to the toggleable HUD widgets) ---
  bool showRadarHud{true};
  bool objectiveHudEnabled{true};
  bool threatHudEnabled{true};
  bool jumpHudEnabled{true};

  // Procedural Ship HUD (instrument cluster)
  bool shipHudEnabled{true};
  int shipHudDetailLevel{2}; // 0=minimal .. 3=max
  bool shipHudShowSparklines{true};
  bool shipHudGlitchFx{true};
  bool shipHudPanelDropouts{true};
  bool shipHudShowDecor{true};
  float shipHudDecorAlpha{0.35f}; // 0..1 (multiplies grid color)
  bool shipHudShowGlyphs{true};
  bool shipHudShowMicrotext{true};
  bool shipHudSegmentText{false}; // procedural 14-seg segment-display readouts for values
  bool shipHudSeedFromLoadout{true};
  std::uint64_t shipHudSeed{0};       // used when seedFromLoadout=false
  std::uint64_t shipHudSeedNonce{0};  // increments when user "rerolls" the procedural layout

  // --- Radar ---
  double radarRangeKm{220000.0};
  int radarMaxBlips{72};

  // --- Misc indicators ---
  bool offscreenTargetIndicator{true};

  // --- Combat symbology ---
  bool combatHudEnabled{true};
  bool useProceduralReticle{true};
  bool showWeaponRings{true};
  bool showHeatRing{true};
  bool showDistributorRings{true};
  float reticleSizePx{44.0f};
  float reticleAlpha{0.80f};

  bool showLeadIndicator{true};
  bool leadUseLastFiredWeapon{true};
  float leadSizePx{22.0f};
  double leadMaxTimeSec{18.0};

  bool showFlightPathMarker{true};
  bool flightMarkerUseLocalFrame{true};
  bool flightMarkerClampToEdge{true};
  float flightMarkerSizePx{22.0f};

  // --- Tactical overlay ---
  bool tacticalOverlayEnabled{true};
  bool tacticalShowLabels{true};
  double tacticalRangeKm{450000.0};
  int tacticalMaxMarkers{96};
  bool tacticalShowStations{true};
  bool tacticalShowPlanets{true};
  bool tacticalShowContacts{true};
  bool tacticalShowCargo{true};
  bool tacticalShowAsteroids{true};
  bool tacticalShowSignals{true};

  // --- Style / Colors (v2+) ---
  // Background opacity for borderless HUD overlay windows.
  float overlayBgAlpha{0.35f};
  float overlayBgAlphaEdit{0.45f};

  // Optional tinting of atlas icons.
  bool tintRadarIcons{false};
  bool tintTacticalIcons{false};

  // HUD palette.
  Color4f colorPrimary{0.82f, 0.90f, 1.00f, 1.00f};
  Color4f colorAccent{1.00f, 0.70f, 0.35f, 1.00f};
  Color4f colorDanger{1.00f, 0.25f, 0.25f, 1.00f};
  Color4f colorGrid{0.47f, 0.55f, 0.67f, 1.00f};
  Color4f colorText{0.90f, 0.90f, 0.94f, 1.00f};
  Color4f colorBackground{0.00f, 0.00f, 0.00f, 1.00f};
};

// Default config file path (relative to the working directory).
std::string defaultHudSettingsPath();

HudSettings makeDefaultHudSettings();

bool saveHudSettingsToFile(const HudSettings& s, const std::string& path);
bool loadHudSettingsFromFile(const std::string& path, HudSettings& out);

} // namespace stellar::ui
