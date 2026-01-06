#pragma once

#include "stellar/math/Vec3.h"

#include <functional>
#include <string>
#include <vector>

namespace stellar::game {

// Keyframe for the cinematic camera rig.
//
// offsetLocalU is expressed in ship-local axes (X=right, Y=up, Z=forward) in render-units.
struct CinematicCameraKeyframe {
  double tSec{0.0};
  math::Vec3d offsetLocalU{0.0, 2.0, -6.0};
  double fovDeg{60.0};
};

struct CinematicCameraSample {
  math::Vec3d offsetLocalU{0.0, 2.0, -6.0};
  double fovDeg{60.0};
  bool hasFov{false};
};

struct CinematicCameraWindowState {
  bool open{false};

  // If true, the rig overrides the default chase-camera offsets.
  bool enabled{false};

  // Playback controls
  bool playing{false};
  bool loop{true};
  bool playUsesSimTime{false}; // advance by dtSim (affected by timeScale)
  bool playWhilePaused{false}; // if false, playback freezes while paused

  // Orientation helper: if true, the main loop can aim the camera at the ship.
  bool lookAtShip{false};

  // If true, apply keyframed FOV. If false, keep the default 60°.
  bool overrideFov{false};

  // Interpolation for offsets between keyframes.
  bool catmullRom{true};

  // Playback state
  double playheadSec{0.0};
  double playbackRate{1.0}; // 1.0 = realtime, 2.0 = double speed, etc.

  // Editor state
  int selectedIndex{-1};
  double editorTimeSec{0.0};
  math::Vec3d editorOffsetU{0.0, 2.0, -6.0};
  double editorFovDeg{60.0};

  // Import/export
  char csvPath[256]{"camera_path.csv"};

  // Keyframes
  std::vector<CinematicCameraKeyframe> keys;
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

// Advance playhead and keep it within the clip bounds.
void tickCinematicCamera(CinematicCameraWindowState& st,
                         double dtRealSec,
                         double dtSimSec,
                         bool paused);

// Sample offset + (optional) FOV at the current playhead.
CinematicCameraSample sampleCinematicCamera(const CinematicCameraWindowState& st);

// Draw UI for the cinematic camera window (safe to call even when st.open == false).
void drawCinematicCameraWindow(CinematicCameraWindowState& st, const ToastFn& toast);

} // namespace stellar::game
