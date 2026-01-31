#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Celestial.h"
#include "stellar/sim/TimeTrial.h"

// Cross-integration: Time Trials can replay a "ghost" of the best run using the
// same lightweight sample representation as the Flight Recorder.
#include "FlightRecorderWindow.h"

#include "GameSignals.h"

#include <deque>
#include <functional>
#include <string>
#include <unordered_map>

namespace stellar::sim {
struct StarSystem;
class Ship;
}

namespace stellar::game {

// Lightweight 3D gameplay loop: "Time Trials" (gate courses).
//
// - UI: generate a deterministic course around a station.
// - Runtime: detect gate passes via segment-plane intersection.
// - HUD: main.cpp draws markers and an objective panel from this state.

enum class TimeTrialPhase : int {
  Inactive = 0,
  Ready,
  Running,
  Docking,
  Finished
};

enum class TimeTrialFinishMode : int {
  GatesOnly = 0,
  DockAtAnchorStation
};

struct TimeTrialWindowState {
  bool open{false};

  // HUD / marker toggles.
  bool hudEnabled{true};
  bool showGateMarker{true};
  bool showAllGates{false};
  bool clampOffscreen{true};
  float markerAlpha{0.9f};
  float markerThickness{2.0f};

  // Course generation settings.
  int stationIndex{0};
  sim::TimeTrialCourseParams params{};

  // Finish mode: either stop at last gate, or require docking at the anchor station
  // after clearing the final gate.
  TimeTrialFinishMode finishMode{TimeTrialFinishMode::GatesOnly};

  // Optional cross-system integrations (handled by main.cpp).
  bool autoTargetAnchorOnDocking{true};
  bool autoEngageDockingComputerOnDocking{true};
  bool autoCameraChaseOnStart{true};
  bool autoCameraOrbitOnDocking{true};

  // Optional run capture (cross-integrates gameplay with devtools).
  // When enabled, a time-trial run can automatically record ship telemetry and
  // export traces on finish for easy repro / sharing.
  bool autoStartFlightRecorderOnStart{true};
  bool autoStopFlightRecorderOnFinish{true};
  bool autoExportTracesOnFinish{false};

  // Clear Integration Hub history when (re)arming a course so exports contain a clean run.
  bool autoClearIntegrationHubOnArm{true};

  // Seed controls.
  bool seedFromSystem{true};
  core::u64 userSeed{0x1234u};

  // Runtime session.
  TimeTrialPhase phase{TimeTrialPhase::Inactive};
  bool hasCourse{false};
  sim::TimeTrialCourse course{};

  // Anchor metadata for cross-system integration.
  sim::StationId anchorStationId{0};
  int anchorStationIndex{-1};

  int nextGate{0};
  double timeSec{0.0};
  double finishTimeSec{0.0};

  // Cached best time for the currently loaded course (course.key).
  double bestTimeSec{0.0};

  bool hasPrevShipPos{false};
  math::Vec3d prevShipPosKm{};

  // Best times (by course key).
  std::unordered_map<core::u64, double> bestTimesSec;

  // --- Ghost racing (best-run replay) ---
  // Records a lightweight ship telemetry trace during runs and replays the
  // best run as an in-world ghost ship + trail so you can race yourself.
  //
  // Best-times + best-ghost data are persisted via SaveGame; live run state is in-memory only.
  bool ghostEnabled{true};
  bool ghostDrawShip{true};
  bool ghostDrawTrail{true};
  bool ghostShowSplitHud{true};
  bool ghostRecordRun{true};
  bool ghostSaveBestRun{true};
  float ghostOpacity{0.55f};
  double ghostSampleHz{30.0};
  int ghostMaxSamples{18000};
  // Optional time offset applied when sampling the best-run ghost (seconds).
  // Positive values push the ghost ahead; negative values lag it behind.
  double ghostLeadSec{0.0};

  // Internal: recording buffer for the active run.
  double ghostSampleAccumulatorSec{0.0};
  std::deque<FlightRecorderSample> ghostRunSamples;

  // Best ghost samples (by course key).
  std::unordered_map<core::u64, std::deque<FlightRecorderSample>> bestGhostSamples;

  // Sampled pose for rendering (consumed by main.cpp).
  bool ghostSampleValid{false};
  FlightRecorderSample ghostSample{};

  // Internal: latch set when a course is armed/re-armed (for one-shot actions).
  bool justArmed{false};

  // One-shot integration requests (set by tickTimeTrial, consumed by main.cpp).
  bool requestTargetAnchor{false};
  bool requestEngageDockingComputer{false};
  bool requestCameraChase{false};
  bool requestCameraOrbit{false};


  // UI feedback.
  int lastGatePassed{-1};
  double gatePulseSec{0.0};

  void clearCourse();
  void armCourse();
  void cancelRun();
  void armFromCourse(const sim::TimeTrialCourse& c);
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

using ActionSinkFn = std::function<void(const GameAction& action)>;
using EventSinkFn  = std::function<void(const GameEvent& event)>;

void tickTimeTrial(TimeTrialWindowState& st,
                   double dtRealSec,
                   double timeRealSec,
                   double simTimeDays,
                   const sim::Ship& ship,
                   bool paused,
                   bool docked,
                   sim::StationId dockedStationId,
                   const ToastFn& toast,
                   const ActionSinkFn& emitAction = {},
                   const EventSinkFn& emitEvent = {});


void drawTimeTrialWindow(TimeTrialWindowState& st,
                         const sim::StarSystem* sys,
                         double timeDays,
                         const sim::Ship& ship,
                         const ToastFn& toast);

} // namespace stellar::game
