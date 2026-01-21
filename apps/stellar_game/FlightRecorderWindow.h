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

struct GameEvent;
enum class GameEventKind : core::u8;

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

// Discrete marker captured from Integration Hub events.
//
// Markers are stored with absolute timestamps and are mapped to the current
// recording session when displayed / exported.
struct FlightRecorderMarker {
  double tRealSec{0.0};
  double tSimDays{0.0};

  GameEventKind kind{};
  std::string tag;
  std::string msg;

  core::u64 u64a{0};
  core::u64 u64b{0};
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

  // Trace: include Integration Hub markers as instant events.
  bool traceIncludeMarkers{true};

  // --- Markers (Integration Hub event stream) ---
  // When enabled, Integration Hub events are mirrored into the flight recorder
  // as discrete markers (useful for correlation in replays/exports).
  bool captureMarkers{true};
  bool captureMarkersIncludeDebug{false};
  bool captureMarkersWhileNotRecording{false};
  int maxMarkers{2048};

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

  // Marker ring buffer (oldest markers dropped).
  std::deque<FlightRecorderMarker> markers;
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

// Programmatic control helpers (used by cross-system gameplay/devtools).
void startFlightRecorderSession(FlightRecorderWindowState& st,
                               double timeRealSec,
                               double simTimeDays,
                               const sim::Ship& ship);

void stopFlightRecorderSession(FlightRecorderWindowState& st);

// Jump the replay playhead to an absolute wall-clock timestamp (timeRealSec).
//
// This is used by cross-system tooling (Integration Hub) to quickly scrub the
// flight recorder to an interesting event time.
//
// If the recorder has no samples, this is a no-op.
void flightRecorderSeekToRealTime(FlightRecorderWindowState& st, double timeRealSec);

// Marker helpers (cross-system integration with the Integration Hub).
void flightRecorderCaptureMarkerFromEvent(FlightRecorderWindowState& st, const GameEvent& ev);
void flightRecorderClearMarkers(FlightRecorderWindowState& st);

// Export helpers (also used by Integration Hub actions).
bool exportFlightRecorderCsv(const FlightRecorderWindowState& st, const char* path, std::string* err = nullptr);
bool exportFlightRecorderTraceJson(const FlightRecorderWindowState& st, const char* path, std::string* err = nullptr);

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
