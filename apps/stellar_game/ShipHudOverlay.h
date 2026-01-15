#pragma once

#include "stellar/ui/HudSettings.h"
#include "stellar/ui/ProceduralShipHud.h"

#include <array>

#include <imgui.h>

namespace stellar::ui {

// Lightweight ring-buffer used for Ship HUD telemetry sparklines.
struct ShipHudSparkline {
  std::array<float, 128> samples{};
  int head{0};
  int count{0};

  void push(float v);

  // Oldest-to-newest sample access.
  float get(int i) const;
};

struct ShipHudHistory {
  double acc{0.0};
  ShipHudSparkline speed{};
  ShipHudSparkline shield{};
  ShipHudSparkline hull{};
  ShipHudSparkline heat{};
  ShipHudSparkline fuel{};
};

// Update sparkline history at a fixed sampling rate.
void shipHudHistoryUpdate(ShipHudHistory& h,
                          double dtRealSec,
                          double speedKmS,
                          float shieldFrac,
                          float hullFrac,
                          float heatFrac,
                          float fuelFrac);

struct ShipHudOverlayConfig {
  bool showSparklines{true};
  bool glitchFx{true};
  bool panelDropouts{true};
  bool showDecor{true};
  float decorAlpha{0.35f}; // 0..1 multiplier on the grid color
  bool showGlyphs{true};
  bool showMicrotext{true};
  bool segmentText{false}; // render value readouts using a procedural segment-display font
};

struct ShipHudOverlayPalette {
  Color4f primary{0.82f, 0.90f, 1.00f, 1.00f};
  Color4f accent{1.00f, 0.70f, 0.35f, 1.00f};
  Color4f danger{1.00f, 0.25f, 0.25f, 1.00f};
  Color4f grid{0.47f, 0.55f, 0.67f, 1.00f};
  Color4f text{0.90f, 0.90f, 0.94f, 1.00f};
  Color4f background{0.00f, 0.00f, 0.00f, 1.00f};
};

struct ShipHudOverlayTelemetry {
  double speedKmS{0.0};
  float shieldFrac{0.0f};
  float hullFrac{0.0f};
  float heatFrac{0.0f};
  float fuelFrac{0.0f};

  int pipSys{0};
  int pipEng{0};
  int pipWep{0};

  float fsdFrac{0.0f};
  char fsdText[32]{}; // short status string

  float throttleFrac{0.0f};

  bool hasNearest{false};
  double nearestKm{0.0};

  double cargoCr{0.0};
  double gForce{0.0};

  // --- Attitude / navigation instrument ---
  // These values are computed by the game layer (main.cpp) so the renderer stays lightweight.
  bool attitudeValid{false};
  bool attitudeUsesGravity{false}; // true: reference frame is gravity up; false: orbit-plane up
  float attitudePitchDeg{0.0f};    // nose up (+) / down (-) relative to reference "horizon"
  float attitudeRollDeg{0.0f};     // right wing down (+) / left wing down (-) relative to reference up
  float attitudeRefG{0.0f};        // reference magnitude (in G) when using gravity

  // Ship-local velocity direction (prograde marker) expressed as yaw/pitch (degrees).
  bool progradeValid{false};
  float progradeYawDeg{0.0f};
  float progradePitchDeg{0.0f};

  // Ship-local retrograde direction (opposite of velocity) expressed as yaw/pitch (degrees).
  bool retrogradeValid{false};
  float retrogradeYawDeg{0.0f};
  float retrogradePitchDeg{0.0f};

  // Ship-local direction to the nearest contact expressed as yaw/pitch (degrees).
  bool targetDirValid{false};
  float targetYawDeg{0.0f};
  float targetPitchDeg{0.0f};

  // Ship-local gravity direction (down vector) expressed as yaw/pitch (degrees).
  bool gravityValid{false};
  float gravityYawDeg{0.0f};
  float gravityPitchDeg{0.0f};


  double timeRealSec{0.0};
};

// Draw the full Ship HUD overlay contents into a borderless ImGui window.
//
// The caller owns the Dear ImGui window creation and only supplies the window
// rect + draw list.
void drawShipHudOverlay(ImDrawList* draw,
                        ImVec2 winPos,
                        ImVec2 winSize,
                        float uiScale,
                        const ProceduralShipHudPlan& plan,
                        const ShipHudHistory& hist,
                        const ShipHudOverlayTelemetry& t,
                        const ShipHudOverlayConfig& cfg,
                        const ShipHudOverlayPalette& pal);

} // namespace stellar::ui
