#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"

#include <deque>
#include <functional>
#include <string>

namespace stellar::sim {
class Ship;
}

namespace stellar::game {

// Lightweight flight telemetry recorder.
//
// Records the player's ship state at a configurable rate and allows exporting
// it for offline analysis. Intended as a debugging + content-creation tool.

struct FlightRecorderSample {
  // Seconds since the start of the current recording session.
  double tRealSec{0.0};

  // Absolute simulation time (days).
  double tSimDays{0.0};

  // Physics state (km / km/s / rad/s).
  math::Vec3d posKm{};
  math::Vec3d velKmS{};

  // Full attitude (for ghost replays / debugging / trace exports).
  math::Quatd orient{1,0,0,0};
  math::Vec3d angVelRadS{};
};

struct FlightRecorderWindowState {
  bool open{false};

  bool recording{false};
  bool recordWhilePaused{false};

  // Sample rate (Hz). Actual sampling is frame-based with an accumulator.
  double sampleHz{20.0};

  // Ring buffer size (oldest samples dropped).
  int maxSamples{20000};

  bool showPlot{true};

  // Export paths (relative to working directory unless absolute).
  char csvPath[256]{"flight_recorder.csv"};
  char tracePath[256]{"flight_recorder_trace.json"};

  // Trace export toggles (counter tracks).
  bool traceIncludeSpeed{true};
  bool traceIncludePosition{true};
  bool traceIncludeVelocity{false};
  bool traceIncludeSimTime{true};
  bool traceIncludeOrientation{false};
  bool traceIncludeAngularVelocity{false};
  bool tracePretty{false};

  // --- Replay / Ghost ---
  // This is intentionally simple and deterministic:
  // - Replay playhead is in recording-relative seconds.
  // - Sampling is linear (nlerp for orientation).
  bool playing{false};
  bool playWhilePaused{false};
  bool playUsesSimTime{true};
  bool loop{true};
  double playbackRate{1.0};
  double playheadSec{0.0};

  // Render a debug "ghost" ship at the current playhead.
  bool ghostEnabled{false};
  bool ghostWireframe{true};
  bool ghostDrawTrail{true};

  // Internal: sampled pose at the current playhead.
  bool ghostSampleValid{false};
  FlightRecorderSample ghostSample{};

  // Internal session state.
  double startRealSec{0.0};
  double startSimDays{0.0};
  double sampleAccumulatorSec{0.0};

  // Internal tick bookkeeping to derive simulation delta-time for replay.
  bool hasPrevTimes{false};
  double prevRealSec{0.0};
  double prevSimDays{0.0};

  std::deque<FlightRecorderSample> samples;
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

// Record telemetry (call once per frame).
void tickFlightRecorder(FlightRecorderWindowState& st,
                        double dtRealSec,
                        double timeRealSec,
                        double simTimeDays,
                        const sim::Ship& ship,
                        bool paused);

// Draw UI for the window (safe to call even when st.open == false).
void drawFlightRecorderWindow(FlightRecorderWindowState& st,
                              const sim::Ship& ship,
                              double timeRealSec,
                              double simTimeDays,
                              bool paused,
                              const ToastFn& toast);

} // namespace stellar::game
