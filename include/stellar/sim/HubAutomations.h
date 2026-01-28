#pragma once

#include "stellar/core/Types.h"

#include <deque>
#include <string>
#include <string_view>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Integration Hub Automations (headless-friendly)
// -----------------------------------------------------------------------------
//
// The playable prototype has an "Integration Hub" window (stellar_game) that
// allows IFTTT-style rules: IF event THEN actions.
//
// SaveGame persists these rules in a dependency-free representation so:
//  - quicksave/quickload retains automation setup
//  - headless builds/tests can validate the serialization format
//
// IMPORTANT:
// The numeric values used here intentionally mirror the enums in
// apps/stellar_game/IntegrationHubWindow.h:
//   - AutomationTagMatch
//   - AutomationValueSource
//   - GameEventKind
//   - GameActionKind
//
// This header is part of the core library and must not depend on the game app.

// Numeric mirror of stellar::game::AutomationTagMatch (stored in SaveGame).
enum class HubAutomationTagMatch : core::u8 {
  Any      = 0,
  Equals   = 1,
  Contains = 2,
  Prefix   = 3,
  Suffix   = 4,
};

// Numeric mirror of stellar::game::AutomationValueSource (stored in SaveGame).
enum class HubAutomationValueSource : core::u8 {
  Constant = 0,
  EventU64a = 1,
  EventU64b = 2,
};

// Persisted "action template" for hub automations.
//
// The engine can resolve ids from the triggering event, apply a delay, and
// optionally expand the msgTemplate (placeholders).
struct IntegrationHubAutomationActionState {
  // Numeric value of stellar::game::GameActionKind.
  core::u8 kind{0};

  // Numeric value of HubAutomationValueSource.
  core::u8 u64aSource{0};
  core::u64 u64aConst{0};

  core::u8 u64bSource{0};
  core::u64 u64bConst{0};

  int i32Const{0};
  bool bConst{false};

  // Optional delay relative to triggering event time.
  double delaySec{0.0};

  // Stored as a normal string; the save format base64-encodes it into a single token.
  std::string msgTemplate{};
};

// Persisted "rule" for hub automations.
struct IntegrationHubAutomationRuleState {
  bool enabled{false};

  // Human-friendly name.
  std::string name{"Rule"};

  // Numeric value of stellar::game::GameEventKind.
  core::u8 eventKind{0};

  // Numeric value of HubAutomationTagMatch.
  core::u8 tagMatch{0};
  std::string tagText{};

  // Cooldown in seconds.
  double cooldownSec{0.25};

  std::vector<IntegrationHubAutomationActionState> actions{};
};

// Lightweight headless event used by the sim-side automation engine.
//
// This intentionally uses simple POD-ish fields so it can be populated by
// higher-level layers (stellar_game) without pulling in UI dependencies.
struct HubAutomationEvent {
  core::u8 kind{0};
  // Optional human-friendly kind name (used by template expansion).
  // Game-side layers can populate this with gameEventKindName(kind).
  std::string kindName{};
  std::string tag{};
  std::string msg{};

  core::u64 u64a{0};
  core::u64 u64b{0};

  // Wall-clock timestamp.
  double tRealSec{0.0};

  // Simulation timestamp (days).
  double tSimDays{0.0};
};

// Generic scheduled action produced by the sim-side automation engine.
//
// Interpretation of the numeric kind and payload fields is owned by the
// consumer (typically stellar_game).
struct HubAutomationAction {
  core::u8 kind{0};
  core::u64 u64a{0};
  core::u64 u64b{0};
  int i32a{0};
  bool b{false};

  double tRealSec{0.0};
  double tSimDays{0.0};

  std::string msg{};
  std::string origin{};
};

// Runtime bookkeeping for applying a persistent rule set.
//
// The SaveGame representation is intentionally "pure data". Cooldowns and
// scheduling queues are runtime-only and stored here.
struct HubAutomationRuntime {
  bool enabled{true};

  // Guardrails against runaway recursion / action spam.
  int maxActionsPerEvent{32};
  int maxScheduledActions{512};

  // Recursion depth counter (incremented per apply call).
  int recursionDepth{0};

  // Per-rule last fired time (real seconds). Size should match rules size.
  std::vector<double> lastFiredRealSec{};

  // Scheduled actions sorted by tRealSec (ascending).
  std::deque<HubAutomationAction> scheduled{};
};

// Returns true if `tag` matches (matchKind, pattern).
bool hubAutomationTagMatches(core::u8 matchKind, std::string_view pattern, std::string_view tag);

// Resolve a u64 parameter from an event or a constant based on sourceKind.
core::u64 hubAutomationResolveU64(core::u8 sourceKind, core::u64 constant, const HubAutomationEvent& e);

// Expand a message template using event fields.
//
// Supported placeholders:
//   {tag}, {tag_safe}, {kind}, {msg}, {msg_safe}, {u64a}, {u64b}, {tRealMs}, {tRealSec}, {tSimDays}, {rule}
std::string hubAutomationExpandTemplate(std::string_view tmpl, const HubAutomationEvent& e, std::string_view ruleName);

// Initialize runtime storage for a given rule set.
// Safe to call repeatedly; will preserve existing last-fired times when possible.
void hubAutomationInitRuntime(HubAutomationRuntime& rt, const std::vector<IntegrationHubAutomationRuleState>& rules);

// Evaluate all matching rules for a single event and schedule resulting actions.
void hubAutomationApplyRules(
    HubAutomationRuntime& rt,
    const std::vector<IntegrationHubAutomationRuleState>& rules,
    const HubAutomationEvent& ev);

// Pop any scheduled actions whose time is <= nowRealSec.
// If outDue is null, actions are simply dropped.
void hubAutomationTickScheduled(
    HubAutomationRuntime& rt,
    double nowRealSec,
    std::vector<HubAutomationAction>* outDue);

} // namespace stellar::sim
