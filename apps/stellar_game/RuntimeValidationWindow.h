#pragma once

#include <functional>
#include <string>
#include <vector>

namespace stellar::sim {
class Ship;
} // namespace stellar::sim

namespace stellar::game {

struct FlightRecorderWindowState;
struct GameEvent;
struct GameAction;
struct GameEventLog;

// Runtime Validation window
//
// Purpose:
//  - Catch state corruption early (NaNs/Infs, absurd magnitudes).
//  - Provide a "copy/paste" report for bug reports.
//  - Offer a tiny suite of deterministic smoke checks for recently added gameplay systems.
//
// This is intentionally a pure ImGui tool (only compiled in the SDL/OpenGL app).

struct ValidationCheck {
  std::string name;
  bool ok{true};
  std::string details;
};

struct RuntimeValidationWindowState {
  bool open{false};

  // Continuous watchdog
  bool watchdogEnabled{true};
  bool pauseOnFailure{true};
  bool logOnFailure{true};
  bool suppressDuplicateReports{true};

  // Magnitude clamps (purely for *detecting* issues, not changing simulation).
  bool magnitudeWatchdog{false};
  double maxAbsPosKm{1.0e12};
  double maxAbsVelKmS{1.0e6};
  double maxAbsAngVelRadS{1.0e6};

  // Latching state
  bool watchdogLatched{false};
  int watchdogHits{0};
  std::string watchdogLastMessage;

  // Optional screenshot request when the watchdog trips (uses Integration Hub action queue).
  bool screenshotOnFailure{true};
  bool screenshotIncludeUi{true};
  bool screenshotIncludeWorld{false};
  bool screenshotCopyToClipboard{false};
  bool screenshotPauseForCapture{true};
  bool screenshotTimestamp{true};
  char screenshotOutDir[128]{"screenshots"};
  char screenshotBaseName[64]{"validation"};
  bool screenshotLatchedThisFailure{false};

  // Repro pack export (JSON)
  // When the watchdog trips, write a small JSON bundle to help reproduce/debug.
  bool dumpReproOnFailure{true};
  bool dumpIncludeTelemetry{true};
  bool dumpIncludeEvents{true};
  int dumpMaxEvents{200};
  double dumpEventWindowSec{30.0};
  bool dumpPretty{true};
  bool dumpUniquePerHit{true};
  double telemetryWindowSec{8.0};
  int telemetryMaxSamples{400};
  char dumpPath[256]{"repro_pack.json"};

  int dumpsWritten{0};
  bool dumpLatchedThisFailure{false};
  std::string lastDumpPath;
  std::string lastDumpError;

  // Deterministic smoke checks
  std::vector<ValidationCheck> lastChecks;
  int lastFailCount{0};
  std::string lastReport;
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;
using EventSinkFn = std::function<void(const GameEvent& ev)>;
using ActionSinkFn = std::function<void(const GameAction& a)>;

// Lightweight per-frame validation (can pause simulation on failure).
void tickRuntimeValidation(RuntimeValidationWindowState& st,
                         const sim::Ship& ship,
                         const FlightRecorderWindowState* flightRecorder,
                         double timeRealSec,
                         double simTimeDays,
                         bool& paused,
                         const ToastFn& toast,
                         const EventSinkFn& emitEvent = {},
                         const ActionSinkFn& emitAction = {},
                         const GameEventLog* eventLogForDump = nullptr);


// ImGui window for running smoke checks + toggling watchdog behavior.
void drawRuntimeValidationWindow(RuntimeValidationWindowState& st,
                               const sim::Ship& ship,
                               const FlightRecorderWindowState* flightRecorder,
                               double timeRealSec,
                               double simTimeDays,
                               bool paused,
                               const ToastFn& toast,
                               const EventSinkFn& emitEvent = {},
                               const ActionSinkFn& emitAction = {},
                               const GameEventLog* eventLogForDump = nullptr);

// Export a repro pack snapshot to an explicit path (without mutating the window state).
// This is used by higher-level capture actions (e.g. CaptureBundle) to gather artifacts
// into a single folder.
bool exportReproPackJsonToPath(const RuntimeValidationWindowState& st,
                              const sim::Ship& ship,
                              const FlightRecorderWindowState* flightRecorder,
                              const GameEventLog* eventLogForDump,
                              double timeRealSec,
                              double simTimeDays,
                              bool paused,
                              const char* path,
                              std::string* outErr = nullptr);

} // namespace stellar::game
