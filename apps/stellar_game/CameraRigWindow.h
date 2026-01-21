#pragma once

#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"

#include <SDL.h>

#include <functional>
#include <string>

struct ImGuiIO;
// SDL provides Uint8 and SDL_Event.

namespace stellar::sim {
struct GravityParams;
struct StarSystem;
class Ship;
} // namespace stellar::sim

namespace stellar::game {

// Camera modes for the 3D camera rig.
//
// Chase: legacy behavior (ship-local offset + optional look-at).
// Orbit: editor-style orbit around a pivot with horizon locking, look-ahead, and collision avoidance.
enum class CameraRigMode : int {
  Chase = 0,
  Orbit = 1,
};

struct CameraRigWindowState {
  bool open{false};

  // Active mode.
  CameraRigMode mode{CameraRigMode::Orbit};

  // Last camera preset applied (for UI/debug). -1 means none yet.
  int lastPresetId{-1};

  // If true and the cinematic camera rig is enabled, use it as the *driver* (offset/FOV),
  // but still apply smoothing + collision avoidance here.
  bool followCinematicWhenEnabled{true};

  // ----- Input (Alt+mouse orbit controls) -----
  bool enableAltMouseOrbit{true};
  float mouseSensitivityDegPerPx{0.15f};
  bool invertY{false};
  float scrollZoomSpeed{0.25f};

  // ----- Orbit camera params -----
  float orbitYawDeg{0.0f};
  float orbitPitchDeg{15.0f};
  float orbitDistanceU{6.0f};

  // Pivot offsets in the camera reference frame (right/up/forward-flat).
  float orbitPivotHeightU{1.2f};
  float orbitPivotShoulderU{0.0f};
  math::Vec3d orbitPanU{0.0, 0.0, 0.0};

  // Horizon stabilization.
  bool orbitUseReferenceUp{true};

  // Look-ahead uses ship velocity to bias the look target.
  bool orbitLookAhead{true};
  float orbitLookAheadSec{1.15f};
  float orbitLookAheadMaxU{2.0f};

  // ----- Chase camera params -----
  // Ship-local offset (X=right, Y=up, Z=forward) in render-units.
  math::Vec3d chaseOffsetLocalU{0.0, 2.0, -6.0};
  bool chaseLookAtShip{false};
  bool chaseUseReferenceUp{false};

  // ----- Smoothing -----
  // Exponential half-life controls (seconds). 0 disables smoothing.
  float posHalfLifeSec{0.18f};
  float rotHalfLifeSec{0.14f};
  float fovHalfLifeSec{0.20f};

  // ----- Field of view -----
  float baseFovDeg{60.0f};
  bool speedFov{true};
  float speedFovGainDeg{14.0f};
  float speedFovMaxAddDeg{18.0f};

  // ----- Distance & FOV tools -----
  bool showDistanceTelemetry{true};

  // Dolly zoom / framing lock: keep apparent framing constant while changing camera distance.
  // When enabled, FOV is computed so that the projection scale at the look-at depth remains constant.
  bool dollyZoomLock{false};
  bool dollyAutoCapture{true};
  bool dollyApplySpeedFov{false};
  float dollyRefDistanceU{6.0f};
  float dollyRefFovDeg{60.0f};
  float dollyMinFovDeg{12.0f};
  float dollyMaxFovDeg{120.0f};

  // Physical camera readout (mm) used only for UI equivalence (does not affect rendering).
  float sensorWidthMm{36.0f};
  float sensorHeightMm{24.0f};

  // Focal length (mm) used for the physical camera equivalence UI.
  float focalLengthMm{35.0f};

  // When framing the dominant body, target what fraction of the vertical view it should fill.
  float frameBodyFill{0.85f};

  // ----- Telemetry (filled by computeCameraRigFrame) -----
  double lastCamDistU{0.0};
  double lastCamDistKm{0.0};
  double lastFinalFovDeg{60.0};
  math::Vec3d lastCamOffsetU{0.0, 0.0, 0.0};

  // ----- Collision / safety -----
  bool avoidBodies{true};
  float bodyPaddingU{0.08f};
  float minCamDistU{0.75f};

  // ----- Runtime state (input) -----
  bool draggingOrbit{false};
  bool draggingPan{false};
  int lastMouseX{0};
  int lastMouseY{0};

  // ----- Runtime state (smoothing) -----
  bool smoothValid{false};
  math::Vec3d smoothPosU{0.0, 0.0, 0.0};
  math::Quatd smoothOrient{math::Quatd::identity()};
  double smoothFovDeg{60.0};
};

struct CameraRigFrame {
  math::Vec3d posU{0.0, 0.0, 0.0};
  math::Quatd orient{math::Quatd::identity()};
  double fovDeg{60.0};
  double nearPlane{0.01};
  double farPlane{20000.0};

  // Debug/visualization: reference up vector actually used for horizon locking.
  math::Vec3d refUp{0.0, 1.0, 0.0};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

// Camera presets: quick ways to configure camera mode + tuning for common contexts.
//
// The numeric ids are intentionally stable because they can be referenced from
// Integration Hub automation actions (GameActionKind::SetCameraRigPreset).
enum class CameraRigPreset : int {
  DefaultOrbit = 0,
  Travel      = 1,
  Docking     = 2,
  Combat      = 3,
  Cinematic   = 4,
};

inline const char* cameraRigPresetName(CameraRigPreset p) {
  switch (p) {
    case CameraRigPreset::DefaultOrbit: return "DefaultOrbit";
    case CameraRigPreset::Travel:      return "Travel";
    case CameraRigPreset::Docking:     return "Docking";
    case CameraRigPreset::Combat:      return "Combat";
    case CameraRigPreset::Cinematic:   return "Cinematic";
    default:                           return "?";
  }
}

// Apply a camera preset. This is a cross-system helper used by main.cpp to
// implement the Integration Hub camera preset action.
void applyCameraRigPreset(CameraRigWindowState& st, CameraRigPreset preset, const ToastFn& toast = {});

// Handle a single SDL event for camera controls (Alt+mouse orbit + hotkeys).
// Returns true if the event was consumed for camera input.
bool handleCameraRigEvent(CameraRigWindowState& st,
                          const SDL_Event& event,
                          const ImGuiIO& io,
                          bool mouseSteerActive,
                          bool actionWheelOpen,
                          const ToastFn& toast);

// Advance any continuous camera state (currently small; reserved for future features).
void tickCameraRig(CameraRigWindowState& st, double dtRealSec, const Uint8* keys, const ImGuiIO& io);

// Compute the final camera pose/FOV for this frame.
CameraRigFrame computeCameraRigFrame(CameraRigWindowState& st,
                                     const sim::Ship& ship,
                                     const math::Vec3d& shipPosU,
                                     const sim::StarSystem* system,
                                     double timeDays,
                                     const sim::GravityParams& gravityParams,
                                     double dtRealSec,
                                     bool cinematicEnabled,
                                     const math::Vec3d& cinematicOffsetLocalU,
                                     bool cinematicLookAtShip,
                                     double cinematicFovDeg,
                                     bool cinematicOverrideFov);

// Utility: retarget the orbit camera pivot to a specific world-space point.
//
// This is a *cross-system* helper used by tools like the Integration Hub to
// quickly focus the camera on interesting events (collisions, gates, docking).
//
// Inputs are expressed in the current render frame (U) to avoid duplicating
// origin-shifting logic in multiple places.
void cameraRigFocusOrbitPivot(CameraRigWindowState& st,
                              const sim::Ship& ship,
                              const sim::StarSystem* system,
                              double timeDays,
                              const sim::GravityParams& gravityParams,
                              const math::Vec3d& shipPosU,
                              const math::Vec3d& targetPosU);

void drawCameraRigWindow(CameraRigWindowState& st,
                         const sim::Ship& ship,
                         const sim::StarSystem* system,
                         double timeDays,
                         const sim::GravityParams& gravityParams,
                         const ToastFn& toast);

} // namespace stellar::game
