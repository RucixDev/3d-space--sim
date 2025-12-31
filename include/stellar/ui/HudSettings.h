#pragma once

#include <string>

namespace stellar::ui {

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
  int version{1};

  // Persist this file on exit.
  bool autoSaveOnExit{true};

  // --- Primary overlays (these map to the toggleable HUD widgets) ---
  bool showRadarHud{true};
  bool objectiveHudEnabled{true};
  bool threatHudEnabled{true};
  bool jumpHudEnabled{true};

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
};

// Default config file path (relative to the working directory).
std::string defaultHudSettingsPath();

HudSettings makeDefaultHudSettings();

bool saveHudSettingsToFile(const HudSettings& s, const std::string& path);
bool loadHudSettingsFromFile(const std::string& path, HudSettings& out);

} // namespace stellar::ui
