#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"

#include <deque>
#include <string>
#include <utility>

namespace stellar::game {

// --- Cross-system integration primitives ---
//
// Many gameplay/devtools systems in stellar_game want to:
//  - request actions that are executed by the central game loop (main.cpp), and
//  - emit structured events for debugging/repro packs.
//
// These lightweight structs provide a single vocabulary for that integration.

enum class GameEventKind : core::u8 {
  // NOTE: These numeric IDs are persisted in SaveGame (Integration Hub automations).
  // Do not reorder existing values.
  Info           = 0,
  TimeTrial      = 1,
  Docking        = 2,
  Camera         = 3,
  NavAssist      = 4,
  Validation     = 5,
  FlightRecorder = 6,
  Gameplay       = 7,

  // Engine/devtools level events.
  Debug          = 8,
};

inline const char* gameEventKindName(GameEventKind k) {
  switch (k) {
    case GameEventKind::Info:           return "Info";
    case GameEventKind::TimeTrial:      return "TimeTrial";
    case GameEventKind::Docking:        return "Docking";
    case GameEventKind::Camera:         return "Camera";
    case GameEventKind::NavAssist:      return "NavAssist";
    case GameEventKind::Validation:     return "Validation";
    case GameEventKind::FlightRecorder: return "FlightRecorder";
    case GameEventKind::Gameplay:       return "Gameplay";
    case GameEventKind::Debug:          return "Debug";
    default:                            return "?";
  }
}

struct GameEvent {
  double tRealSec{0.0};
  double tSimDays{0.0};

  GameEventKind kind{GameEventKind::Info};

  // Short type-like tag (e.g. "GatePassed", "ClearanceGranted").
  std::string tag;

  // Human-readable summary.
  std::string msg;

  // Optional spatial context.
  bool hasPos{false};
  math::Vec3d posKm{};

  // Optional numeric fields for cross-referencing IDs or small payloads.
  core::u64 u64a{0};
  core::u64 u64b{0};
};

// Helper for concise event emission without fragile aggregate ordering.
inline GameEvent makeEvent(double tRealSec,
                           double tSimDays,
                           GameEventKind kind,
                           std::string tag,
                           std::string msg,
                           bool hasPos = false,
                           math::Vec3d posKm = {},
                           core::u64 u64a = 0,
                           core::u64 u64b = 0) {
  GameEvent ev;
  ev.tRealSec = tRealSec;
  ev.tSimDays = tSimDays;
  ev.kind = kind;
  ev.tag = std::move(tag);
  ev.msg = std::move(msg);
  ev.hasPos = hasPos;
  ev.posKm = posKm;
  ev.u64a = u64a;
  ev.u64b = u64b;
  return ev;
}

struct GameEventLog {
  int maxEvents{2000};
  std::deque<GameEvent> events;

  void clear() { events.clear(); }
};

inline void pushGameEvent(GameEventLog& log, GameEvent ev) {
  if (log.maxEvents <= 0) return;
  while ((int)log.events.size() >= log.maxEvents) {
    log.events.pop_front();
  }
  log.events.push_back(std::move(ev));
}


enum class GameActionKind : core::u8 {
  // Stable numeric IDs (persisted in SaveGame via Integration Hub automation rules).
  // Do not reorder existing values. New values must use an explicit number.

  // Targeting/automation
  SetTargetStationId    = 0,
  EngageDockingComputer = 1,

  // Mission tracker / navigation bridge
  // Uses GameAction::u64a = MissionId, and GameAction::b = arm auto-run (optional)
  // to cross-integrate missions with the route planner and auto-run runner.
  SyncNavToMission      = 2,

  // Mission tracker / HUD bridge
  // Uses GameAction::u64a = MissionId (0 clears the tracked mission).
  SetTrackedMissionId   = 3,

  // Camera
  SetCameraRigMode      = 4,

  // Generic (UI)
  Toast                 = 5,

  // ---- Devtools / capture / exports ----
  StartFlightRecorder        = 6,
  StopFlightRecorder         = 7,

  // Export paths are carried in GameAction::msg.
  ExportFlightRecorderTrace  = 8,
  ExportIntegrationTrace     = 9,

  // Comms / diegetic notifications
  // The Comms payload is carried in GameAction::msg. Use `|` to separate subject and body.
  TransmitComms              = 10,

  // Clears the Integration Hub event + action history (useful for capturing a clean run).
  ClearIntegrationHub        = 11,

  // Screenshot request (bridges Integration Hub / automation to Photo Mode capture pipeline)
  // Payload:
  //   b    : include UI (capture after UI pass)
  //   i32a : ScreenshotActionFlags bitmask
  //   msg  : "outDir|baseName" OR "baseName" (dir defaults to "screenshots")
  RequestScreenshot          = 12,

  // Capture a multi-artifact diagnostics bundle (screenshots + traces + optional repro pack).
  // Payload:
  //   i32a : CaptureBundleFlags bitmask
  //   msg  : "outDir|label" OR "label" (outDir defaults to "captures")
  CaptureBundle              = 13,

  // High-level travel command.
  // Payload:
  //   u64b : SystemId (0 => current system)
  //   u64a : StationId (0 => system-only travel)
  //   b   : arm auto-run (auto jump/supercruise + auto-dock when possible)
  GoToStation                = 14,

  // Camera preset (mode + tuning).
  // Payload:
  //   i32a : preset id (0=DefaultOrbit, 1=Travel, 2=Docking, 3=Combat, 4=Cinematic)
  SetCameraRigPreset         = 15,
};

inline const char* gameActionKindName(GameActionKind k) {
  switch (k) {
    case GameActionKind::SetTargetStationId:       return "SetTargetStationId";
    case GameActionKind::EngageDockingComputer:    return "EngageDockingComputer";
    case GameActionKind::SyncNavToMission:       return "SyncNavToMission";
    case GameActionKind::SetTrackedMissionId:     return "SetTrackedMissionId";
    case GameActionKind::SetCameraRigMode:         return "SetCameraRigMode";
    case GameActionKind::Toast:                    return "Toast";
    case GameActionKind::StartFlightRecorder:      return "StartFlightRecorder";
    case GameActionKind::StopFlightRecorder:       return "StopFlightRecorder";
    case GameActionKind::ExportFlightRecorderTrace:return "ExportFlightRecorderTrace";
    case GameActionKind::ExportIntegrationTrace:   return "ExportIntegrationTrace";
    case GameActionKind::TransmitComms:            return "TransmitComms";
    case GameActionKind::ClearIntegrationHub:      return "ClearIntegrationHub";
    case GameActionKind::RequestScreenshot:        return "RequestScreenshot";
    case GameActionKind::CaptureBundle:            return "CaptureBundle";
    case GameActionKind::GoToStation:              return "GoToStation";
    case GameActionKind::SetCameraRigPreset:        return "SetCameraRigPreset";
    default:                                       return "?";
  }
}

// Flags for GameActionKind::RequestScreenshot (stored in GameAction::i32a).
enum ScreenshotActionFlags : int {
  ScreenshotFlag_Timestamp       = 1 << 0,
  ScreenshotFlag_CopyToClipboard = 1 << 1,
  ScreenshotFlag_PauseForCapture = 1 << 2,
};

// Flags for GameActionKind::CaptureBundle (stored in GameAction::i32a).
//
// CaptureBundle is a high-level *meta action* that orchestrates multiple
// lower-level capture/export actions into a single folder for bug reports
// and content sharing.
//
// The main loop expands CaptureBundle into concrete actions (screenshots,
// trace exports) so the Integration Hub remains dependency-light.
enum CaptureBundleFlags : int {
  CaptureBundle_WorldScreenshot      = 1 << 0,
  CaptureBundle_UiScreenshot         = 1 << 1,
  CaptureBundle_IntegrationTrace     = 1 << 2,
  CaptureBundle_FlightRecorderTrace  = 1 << 3,
  CaptureBundle_ReproPack            = 1 << 4,
  CaptureBundle_StopFlightRecorder   = 1 << 5,
  CaptureBundle_PauseForScreenshots  = 1 << 6,
  CaptureBundle_CopyDirToClipboard   = 1 << 7,

  CaptureBundle_Default =
      CaptureBundle_WorldScreenshot |
      CaptureBundle_UiScreenshot |
      CaptureBundle_IntegrationTrace |
      CaptureBundle_FlightRecorderTrace |
      CaptureBundle_ReproPack |
      CaptureBundle_CopyDirToClipboard,
};

struct GameAction {
  double tRealSec{0.0};
  double tSimDays{0.0};

  // Optional human-friendly source tag (e.g. "TimeTrial").
  std::string origin;

  GameActionKind kind{GameActionKind::Toast};

  // Generic payload fields.
  core::u64 u64a{0};
  core::u64 u64b{0};
  int i32a{0};
  bool b{false};

  // Free-form payload. Used for Toast, export paths, Comms transmissions, and screenshot requests.
  // For TransmitComms: `subject|body` (if `|` missing, subject is derived from origin and body=msg).
  // For RequestScreenshot: `outDir|baseName` (if `|` missing, defaults outDir to "screenshots").
  std::string msg;
};

struct GameActionQueue {
  int maxPending{256};
  int maxHistory{512};

  std::deque<GameAction> pending;
  std::deque<GameAction> history;

  void clear() {
    pending.clear();
    history.clear();
  }
};

inline void pushGameAction(GameActionQueue& q, GameAction a) {
  if (q.maxPending <= 0) return;
  while ((int)q.pending.size() >= q.maxPending) {
    q.pending.pop_front();
  }
  q.pending.push_back(std::move(a));
}

// Pop the next pending action and also store it in history (for debugging/export).
inline bool popGameAction(GameActionQueue& q, GameAction* outAction) {
  if (q.pending.empty()) return false;

  GameAction a = std::move(q.pending.front());
  q.pending.pop_front();

  // Record history.
  if (q.maxHistory > 0) {
    while ((int)q.history.size() >= q.maxHistory) {
      q.history.pop_front();
    }
    q.history.push_back(a); // copy for history
  }

  if (outAction) *outAction = std::move(a);
  return true;
}


// Convenience constructors.
inline GameAction makeActionToast(double tRealSec, double tSimDays, std::string origin, std::string msg) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::Toast;
  a.msg = std::move(msg);
  return a;
}

inline GameAction makeActionSetTargetStation(double tRealSec, double tSimDays, std::string origin, core::u64 stationId) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::SetTargetStationId;
  a.u64a = stationId;
  return a;
}

inline GameAction makeActionEngageDockingComputer(double tRealSec, double tSimDays, std::string origin, bool engage) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::EngageDockingComputer;
  a.b = engage;
  return a;
}

inline GameAction makeActionGoToStation(double tRealSec, double tSimDays, std::string origin,
                                        core::u64 systemId, core::u64 stationId, bool armAutoRun) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::GoToStation;
  a.u64b = systemId;
  a.u64a = stationId;
  a.b = armAutoRun;
  return a;
}

inline GameAction makeActionSetTrackedMission(double tRealSec, double tSimDays, std::string origin, core::u64 missionId) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::SetTrackedMissionId;
  a.u64a = missionId;
  return a;
}

inline GameAction makeActionSyncNavToMission(double tRealSec,
                                             double tSimDays,
                                             std::string origin,
                                             core::u64 missionId,
                                             bool armAutoRun,
                                             std::string msg = {}) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::SyncNavToMission;
  a.u64a = missionId;
  a.b = armAutoRun;
  a.msg = std::move(msg);
  return a;
}

inline GameAction makeActionSetCameraRigPreset(double tRealSec, double tSimDays, std::string origin, int presetId) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::SetCameraRigPreset;
  a.i32a = presetId;
  return a;
}


inline GameAction makeActionSetCameraRigMode(double tRealSec, double tSimDays, std::string origin, int mode) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::SetCameraRigMode;
  a.i32a = mode;
  return a;
}

inline GameAction makeActionStartFlightRecorder(double tRealSec, double tSimDays, std::string origin) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::StartFlightRecorder;
  return a;
}

inline GameAction makeActionStopFlightRecorder(double tRealSec, double tSimDays, std::string origin) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::StopFlightRecorder;
  return a;
}

inline GameAction makeActionExportFlightRecorderTrace(double tRealSec,
                                                     double tSimDays,
                                                     std::string origin,
                                                     std::string path) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::ExportFlightRecorderTrace;
  a.msg = std::move(path);
  return a;
}

inline GameAction makeActionExportIntegrationTrace(double tRealSec,
                                                  double tSimDays,
                                                  std::string origin,
                                                  std::string path) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::ExportIntegrationTrace;
  a.msg = std::move(path);
  return a;
}


inline GameAction makeActionTransmitComms(double tRealSec,
                                         double tSimDays,
                                         std::string origin,
                                         int channel,
                                         core::u64 stationId,
                                         bool showOverlay,
                                         std::string msg) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::TransmitComms;
  a.i32a = channel;
  a.u64a = stationId;
  a.b = showOverlay;
  a.msg = std::move(msg);
  return a;
}

inline GameAction makeActionClearIntegrationHub(double tRealSec, double tSimDays, std::string origin) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::ClearIntegrationHub;
  return a;
}

inline GameAction makeActionRequestScreenshot(double tRealSec,
                                             double tSimDays,
                                             std::string origin,
                                             bool includeUi,
                                             int flags,
                                             std::string outDir,
                                             std::string baseName) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::RequestScreenshot;
  a.b = includeUi;
  a.i32a = flags;

  if (outDir.empty()) outDir = "screenshots";
  if (baseName.empty()) baseName = "shot";

  a.msg = std::move(outDir);
  a.msg.push_back('|');
  a.msg += baseName;
  return a;
}

inline GameAction makeActionCaptureBundle(double tRealSec,
                                         double tSimDays,
                                         std::string origin,
                                         int flags,
                                         std::string outDir,
                                         std::string label) {
  GameAction a;
  a.tRealSec = tRealSec;
  a.tSimDays = tSimDays;
  a.origin = std::move(origin);
  a.kind = GameActionKind::CaptureBundle;
  a.i32a = flags;

  if (outDir.empty()) outDir = "captures";
  if (label.empty()) label = "bundle";

  a.msg = std::move(outDir);
  a.msg.push_back('|');
  a.msg += label;
  return a;
}

} // namespace stellar::game
