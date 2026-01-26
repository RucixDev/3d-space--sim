#pragma once

#include "GameSignals.h"

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"

#include <deque>
#include <functional>
#include <string>
#include <vector>

namespace stellar::game {

// Integration Hub
//
// A cross-system "glue" window that:
//  - collects structured events (GameEventLog)
//  - queues cross-system actions (GameActionQueue)
//  - exports a compact JSON trace for bug reports
//  - (NEW) can run simple automation rules (IFTTT-style) that translate events -> actions
//
// This is intentionally tiny and dependency-light, so other windows can emit
// actions/events without including main.cpp.

// -----------------------------------------------------------------------------
// Automations
// -----------------------------------------------------------------------------
//
// The goal is to let you stitch existing systems together *without* hard-coding
// that glue into every window:
//
//   "When a TimeTrial finishes -> stop flight recorder -> export traces -> open FlightRecorder"
//
// Rules are persisted in SaveGame (quicksave/quickload) and are also exported
// in integration traces so bug repro packs can include the exact automation setup.

enum class AutomationTagMatch : int {
  Any = 0,
  Equals,
  Contains,
  Prefix,
  Suffix,
};

inline const char* automationTagMatchName(AutomationTagMatch m) {
  switch (m) {
    case AutomationTagMatch::Any:      return "Any";
    case AutomationTagMatch::Equals:   return "Equals";
    case AutomationTagMatch::Contains: return "Contains";
    case AutomationTagMatch::Prefix:   return "Prefix";
    case AutomationTagMatch::Suffix:   return "Suffix";
    default:                           return "?";
  }
}

enum class AutomationValueSource : int {
  Constant = 0,
  EventU64a,
  EventU64b,
};

inline const char* automationValueSourceName(AutomationValueSource s) {
  switch (s) {
    case AutomationValueSource::Constant: return "Constant";
    case AutomationValueSource::EventU64a: return "Event.u64a";
    case AutomationValueSource::EventU64b: return "Event.u64b";
    default: return "?";
  }
}

struct AutomationActionTemplate {
  GameActionKind kind{GameActionKind::Toast};

  // For actions that need numeric ids.
  //
  // Most actions historically used only u64a, but newer higher-level actions
  // (e.g. GoToStation) can use both u64a and u64b. The automation system lets
  // each field be sourced independently from the incoming event.
  AutomationValueSource u64aSource{AutomationValueSource::Constant};
  core::u64 u64aConst{0};

  AutomationValueSource u64bSource{AutomationValueSource::Constant};
  core::u64 u64bConst{0};

  // For actions that need a small int (e.g., SetCameraRigMode).
  int i32Const{0};

  // For actions that need a boolean (e.g., EngageDockingComputer).
  bool bConst{false};

  // Optional delay for this action relative to the triggering event (seconds).
  // Useful for sequencing multi-step workflows (stop recorder -> export, etc.).
  float delaySec{0.0f};

  // For actions that carry strings (Toast message, export path, ...).
  // Supports lightweight {placeholders}; see IntegrationHubWindow.cpp.
  char msgTemplate[192]{0};
};

struct AutomationRule {
  bool enabled{false};
  char name[64]{"Rule"};

  GameEventKind kind{GameEventKind::Info};
  AutomationTagMatch tagMatch{AutomationTagMatch::Any};
  char tagText[64]{0};

  // Cooldown in seconds to prevent accidental per-frame spam.
  double cooldownSec{0.25};
  double lastFiredRealSec{-1e9};

  std::vector<AutomationActionTemplate> actions;
};

// -----------------------------------------------------------------------------
// Hub state
// -----------------------------------------------------------------------------

struct IntegrationHubWindowState {
  bool open{false};

  // Central queues.
  GameActionQueue actions;
  GameEventLog events;

  // Scheduled actions are held here until their tRealSec is reached, then
  // moved into the main pending action queue.
  std::deque<GameAction> scheduledActions;
  int maxScheduledActions{256};

  // UI
  bool autoScroll{true};
  bool showPendingActions{true};
  bool showScheduledActions{true};
  bool showActionHistory{true};

  // Event filters
  bool showInfo{true};
  bool showTimeTrial{true};
  bool showDocking{true};
  bool showCamera{true};
  bool showNavAssist{true};
  bool showValidation{true};
  bool showFlightRecorder{true};
  bool showGameplay{true};
  bool showDebug{true};

  // Automations
  bool automationsEnabled{true};
  bool automationsInitialized{false};
  int automationsMaxActionsPerFrame{32};
  int automationsActionsFiredThisFrame{0};
  int automationsRecursionDepth{0};
  std::vector<AutomationRule> automationRules;

  // Export
  bool exportPretty{true};
  bool exportIncludePendingActions{true};
  bool exportIncludeScheduledActions{true};
  bool exportIncludeActionHistory{true};
  bool exportIncludeEvents{true};
  bool exportIncludeAutomationRules{true};
  char exportPath[256]{"integration_trace.json"};

  // Lightweight stats
  int totalEventsPushed{0};
  int totalActionsPushed{0};
  int totalActionsScheduled{0};

  // Current frame times (populated by the main loop).
  // Useful for UI-triggered actions that should be timestamped consistently.
  double nowRealSec{0.0};
  double nowSimDays{0.0};

  // Optional sink: observe every incoming event (used to mirror events into
  // other tools such as the Flight Recorder marker track).
  //
  // The sink is invoked *before* the event is moved into the log.
  std::function<void(const GameEvent&)> onEvent;

  // --- UI state ---
  // Selected event index in events.events (or -1 for none).
  int selectedEventIndex{-1};
  // Splitter for the Events tab.
  float eventListWidth{520.0f};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

// Optional UI hooks for cross-system navigation.
//
// The Integration Hub is intentionally dependency-light; these callbacks allow
// the main game loop to wire the timeline to other systems (camera rig, flight
// recorder, etc.) without hard dependencies.
struct IntegrationHubUiHooks {
  // Focus the camera on a world-space point (km).
  std::function<void(const math::Vec3d& posKm)> focusCameraAtPosKm;

  // Scrub the flight recorder to a wall-clock timestamp (timeRealSec).
  std::function<void(double tRealSec)> scrubFlightRecorderToRealTimeSec;
};

// Initialize a set of helpful (disabled-by-default) starter rules.
void initDefaultAutomationRules(IntegrationHubWindowState& st);

// Evaluate rules for a single event and enqueue any resulting actions.
void applyAutomationRules(IntegrationHubWindowState& st, const GameEvent& ev);

inline void hubPushEvent(IntegrationHubWindowState& st, GameEvent ev) {
  st.totalEventsPushed++;

  // Optional sink (e.g. Flight Recorder markers) observes the raw incoming event.
  if (st.onEvent) {
    st.onEvent(ev);
  }

  // Automations are evaluated on the *incoming* event (even if the event log is disabled).
  if (!st.automationsInitialized) initDefaultAutomationRules(st);
  if (st.automationsEnabled) {
    applyAutomationRules(st, ev);
  }

  // Always record the event for debugging (subject to log capacity).
  pushGameEvent(st.events, std::move(ev));
}

inline void hubPushAction(IntegrationHubWindowState& st, GameAction a) {
  st.totalActionsPushed++;
  pushGameAction(st.actions, std::move(a));
}

// Schedule an action for future execution. The action's tRealSec is treated
// as its execution time.
void hubScheduleAction(IntegrationHubWindowState& st, GameAction a);

// Move any due scheduled actions into the pending queue.
void hubTickScheduledActions(IntegrationHubWindowState& st, double nowRealSec);

// Write a JSON trace to disk.
bool writeIntegrationTraceJson(const IntegrationHubWindowState& st, const char* path, std::string* outErr);

// ImGui window.
void drawIntegrationHubWindow(IntegrationHubWindowState& st,
                              const ToastFn& toast,
                              const IntegrationHubUiHooks* hooks = nullptr);

} // namespace stellar::game
