#include "IntegrationHubWindow.h"

#include "stellar/core/JsonWriter.h"
#include "stellar/core/AtomicWriteFile.h"
#include "stellar/core/Types.h"

#include <imgui.h>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string_view>

namespace stellar::game {


static void jsonVec3(stellar::core::JsonWriter& w, const stellar::math::Vec3d& v) {
  w.beginArray();
  w.value(v.x);
  w.value(v.y);
  w.value(v.z);
  w.endArray();
}

static void jsonAction(stellar::core::JsonWriter& w, const GameAction& a) {
  w.beginObject();
  w.key("tRealSec"); w.value(a.tRealSec);
  w.key("tSimDays"); w.value(a.tSimDays);
  w.key("origin"); w.value(a.origin);
  w.key("kind"); w.value(gameActionKindName(a.kind));

  // Payload (compact but explicit).
  switch (a.kind) {
    case GameActionKind::SetTargetStationId:
      w.key("stationId"); w.value((long long)a.u64a);
      break;
    case GameActionKind::EngageDockingComputer:
      w.key("engage"); w.value(a.b);
      break;
    case GameActionKind::SyncNavToMission:
      w.key("missionId"); w.value((long long)a.u64a);
      w.key("armAutoRun"); w.value(a.b);
      break;
    case GameActionKind::SetTrackedMissionId:
      w.key("missionId"); w.value((long long)a.u64a);
      break;
    case GameActionKind::GoToStation:
      w.key("systemId"); w.value((long long)a.u64b);
      w.key("stationId"); w.value((long long)a.u64a);
      w.key("armAutoRun"); w.value(a.b);
      break;
    case GameActionKind::SetCameraRigMode:
      w.key("mode"); w.value(a.i32a);
      break;
    case GameActionKind::SetCameraRigPreset:
      w.key("presetId"); w.value(a.i32a);
      break;
    case GameActionKind::Toast:
      w.key("msg"); w.value(a.msg);
      break;

    case GameActionKind::TransmitComms:
      w.key("channel"); w.value(a.i32a);
      w.key("showOverlay"); w.value(a.b);
      if (a.u64a) { w.key("stationId"); w.value((long long)a.u64a); }
      w.key("msg"); w.value(a.msg);
      break;

    case GameActionKind::StartFlightRecorder:
    case GameActionKind::StopFlightRecorder:
    case GameActionKind::ClearIntegrationHub:
      // no payload
      break;

    case GameActionKind::ExportFlightRecorderTrace:
    case GameActionKind::ExportIntegrationTrace:
      w.key("path"); w.value(a.msg);
      break;

    case GameActionKind::RequestScreenshot:
      w.key("includeUi"); w.value(a.b);
      w.key("flags"); w.value(a.i32a);
      w.key("spec"); w.value(a.msg);
      break;

    case GameActionKind::CaptureBundle:
      w.key("flags"); w.value(a.i32a);
      w.key("spec"); w.value(a.msg);
      break;

    default:
      w.key("u64a"); w.value((long long)a.u64a);
      if (a.u64a) {
        const std::string s = std::to_string((unsigned long long)a.u64a);
        w.key("u64a_str"); w.value(std::string_view{s});
      }
      w.key("u64b"); w.value((long long)a.u64b);
      if (a.u64b) {
        const std::string s = std::to_string((unsigned long long)a.u64b);
        w.key("u64b_str"); w.value(std::string_view{s});
      }
      w.key("i32a"); w.value(a.i32a);
      w.key("b"); w.value(a.b);
      w.key("msg"); w.value(a.msg);
      break;
  }

  w.endObject();
}

static void jsonEvent(stellar::core::JsonWriter& w, const GameEvent& e) {
  w.beginObject();
  w.key("tRealSec"); w.value(e.tRealSec);
  w.key("tSimDays"); w.value(e.tSimDays);
  w.key("kind"); w.value(gameEventKindName(e.kind));
  w.key("tag"); w.value(e.tag);
  w.key("msg"); w.value(e.msg);
  if (e.hasPos) {
    w.key("posKm"); jsonVec3(w, e.posKm);
  }
  if (e.u64a) {
    w.key("u64a"); w.value((long long)e.u64a);
    const std::string s = std::to_string((unsigned long long)e.u64a);
    w.key("u64a_str"); w.value(std::string_view{s});
  }
  if (e.u64b) {
    w.key("u64b"); w.value((long long)e.u64b);
    const std::string s = std::to_string((unsigned long long)e.u64b);
    w.key("u64b_str"); w.value(std::string_view{s});
  }
  w.endObject();
}

static void jsonAutomationAction(stellar::core::JsonWriter& w, const AutomationActionTemplate& a) {
  w.beginObject();
  w.key("kind"); w.value(gameActionKindName(a.kind));
  w.key("u64aSource"); w.value(automationValueSourceName(a.u64aSource));
  w.key("u64aConst"); w.value((long long)a.u64aConst);
  w.key("u64bSource"); w.value(automationValueSourceName(a.u64bSource));
  w.key("u64bConst"); w.value((long long)a.u64bConst);
  w.key("i32Const"); w.value(a.i32Const);
  w.key("bConst"); w.value(a.bConst);
  w.key("delaySec"); w.value(a.delaySec);
  w.key("msgTemplate"); w.value(a.msgTemplate);
  w.endObject();
}

static void jsonAutomationRule(stellar::core::JsonWriter& w, const AutomationRule& r) {
  w.beginObject();
  w.key("enabled"); w.value(r.enabled);
  w.key("name"); w.value(r.name);
  w.key("eventKind"); w.value(gameEventKindName(r.kind));
  w.key("tagMatch"); w.value(automationTagMatchName(r.tagMatch));
  w.key("tagText"); w.value(r.tagText);
  w.key("cooldownSec"); w.value(r.cooldownSec);

  w.key("actions");
  w.beginArray();
  for (const auto& a : r.actions) jsonAutomationAction(w, a);
  w.endArray();

  w.endObject();
}

static bool passesFilter(const IntegrationHubWindowState& st, const GameEvent& e) {
  switch (e.kind) {
    case GameEventKind::Info: return st.showInfo;
    case GameEventKind::TimeTrial: return st.showTimeTrial;
    case GameEventKind::Docking: return st.showDocking;
    case GameEventKind::Camera: return st.showCamera;
    case GameEventKind::NavAssist: return st.showNavAssist;
    case GameEventKind::Validation: return st.showValidation;
    case GameEventKind::FlightRecorder: return st.showFlightRecorder;
    case GameEventKind::Gameplay: return st.showGameplay;
    case GameEventKind::Debug: return st.showDebug;
    default: return true;
  }
}

static std::string formatActionLine(const GameAction& a) {
  std::ostringstream ss;
  ss.setf(std::ios::fixed);
  ss.precision(3);
  ss << a.tRealSec << "s | " << gameActionKindName(a.kind);
  if (!a.origin.empty()) ss << " (" << a.origin << ")";

  switch (a.kind) {
    case GameActionKind::SetTargetStationId:
      ss << " stationId=" << a.u64a;
      break;
    case GameActionKind::SyncNavToMission:
      ss << " missionId=" << a.u64a;
      if (a.b) ss << " autorun";
      break;
    case GameActionKind::GoToStation:
      ss << " sysId=" << a.u64b;
      if (a.u64a) ss << " stationId=" << a.u64a;
      if (a.b) ss << " autorun";
      break;
    case GameActionKind::EngageDockingComputer:
      ss << (a.b ? " engage" : " disengage");
      break;
    case GameActionKind::SetCameraRigMode:
      ss << " mode=" << a.i32a;
      break;
    case GameActionKind::SetCameraRigPreset:
      ss << " preset=" << a.i32a;
      break;
    case GameActionKind::Toast:
      ss << " \"" << a.msg << "\"";
      break;

    case GameActionKind::StartFlightRecorder:
      ss << " start";
      break;
    case GameActionKind::StopFlightRecorder:
      ss << " stop";
      break;
    case GameActionKind::ExportFlightRecorderTrace:
      ss << " path=" << a.msg;
      break;
    case GameActionKind::ExportIntegrationTrace:
      ss << " path=" << a.msg;
      break;
    case GameActionKind::ClearIntegrationHub:
      ss << " clear";
      break;

    case GameActionKind::RequestScreenshot:
      ss << (a.b ? " ui" : " world");
      ss << " flags=" << a.i32a;
      if (!a.msg.empty()) ss << " spec=" << a.msg;
      break;

    case GameActionKind::CaptureBundle:
      ss << " flags=" << a.i32a;
      if (!a.msg.empty()) ss << " spec=" << a.msg;
      break;

    default:
      break;
  }

  return ss.str();
}
void hubScheduleAction(IntegrationHubWindowState& st, GameAction a) {
  st.totalActionsScheduled++;

  // Enforce a hard cap to avoid unbounded growth if the user creates a runaway
  // automation. Oldest scheduled actions are dropped first.
  if (st.maxScheduledActions > 0 && (int)st.scheduledActions.size() >= st.maxScheduledActions) {
    st.scheduledActions.pop_front();
  }

  // Keep scheduled actions ordered by execution time.
  const double t = a.tRealSec;
  auto it = std::lower_bound(st.scheduledActions.begin(), st.scheduledActions.end(), t,
                             [](const GameAction& lhs, double rhsT) { return lhs.tRealSec < rhsT; });
  st.scheduledActions.insert(it, std::move(a));
}

void hubTickScheduledActions(IntegrationHubWindowState& st, double nowRealSec) {
  // Move any due actions into the pending queue.
  while (!st.scheduledActions.empty()) {
    const GameAction& a = st.scheduledActions.front();
    if (a.tRealSec > nowRealSec + 1e-9) break;

    // NOTE: scheduled actions are pushed into the normal action queue so they
    // appear in the history and can be exported as part of the trace.
    hubPushAction(st, st.scheduledActions.front());
    st.scheduledActions.pop_front();
  }
}

static std::string formatEventLine(const GameEvent& e) {
  std::ostringstream ss;
  ss.setf(std::ios::fixed);
  ss.precision(3);
  ss << e.tRealSec << "s | " << gameEventKindName(e.kind);
  if (!e.tag.empty()) ss << ":" << e.tag;
  if (!e.msg.empty()) ss << " — " << e.msg;
  return ss.str();
}

static void copyTextToClipboard(const std::string& s) {
  ImGui::SetClipboardText(s.c_str());
}

static bool tagMatches(const AutomationRule& r, std::string_view tag) {
  const std::string_view pat(r.tagText);
  switch (r.tagMatch) {
    case AutomationTagMatch::Any:
      return true;
    case AutomationTagMatch::Equals:
      return (!pat.empty()) && (tag == pat);
    case AutomationTagMatch::Contains:
      return (!pat.empty()) && (tag.find(pat) != std::string_view::npos);
    case AutomationTagMatch::Prefix:
      return (!pat.empty()) && (tag.size() >= pat.size()) && (tag.substr(0, pat.size()) == pat);
    case AutomationTagMatch::Suffix:
      return (!pat.empty()) && (tag.size() >= pat.size()) && (tag.substr(tag.size() - pat.size()) == pat);
    default:
      return false;
  }
}

static stellar::core::u64 resolveU64(
    AutomationValueSource src,
    stellar::core::u64 constant,
    const GameEvent& e) {
  switch (src) {
    case AutomationValueSource::Constant: return constant;
    case AutomationValueSource::EventU64a: return e.u64a;
    case AutomationValueSource::EventU64b: return e.u64b;
    default: return constant;
  }
}

static std::string sanitizeForFilename(std::string s) {
  // Keep it conservative: avoid platform-specific illegal chars.
  for (char& c : s) {
    const unsigned char uc = (unsigned char)c;
    if (std::isalnum(uc) || c == '_' || c == '-' || c == '.' ) continue;
    c = '_';
  }
  return s;
}

static std::string expandTemplate(std::string_view tmpl, const GameEvent& e, std::string_view ruleName) {
  // Small, dependency-free template replacement:
  //   {tag}, {kind}, {msg}, {u64a}, {u64b}, {tRealMs}, {tRealSec}, {tSimDays}, {rule}
  std::string out;
  out.reserve(tmpl.size() + 32);

  auto appendNum = [&](stellar::core::u64 v) {
    out += std::to_string((unsigned long long)v);
  };

  for (std::size_t i = 0; i < tmpl.size(); ) {
    const char c = tmpl[i];
    if (c != '{') {
      out.push_back(c);
      ++i;
      continue;
    }

    const std::size_t j = tmpl.find('}', i + 1);
    if (j == std::string_view::npos) {
      out.append(tmpl.substr(i));
      break;
    }

    const std::string_view key = tmpl.substr(i + 1, j - (i + 1));

    if (key == "tag") {
      out.append(e.tag);
    } else if (key == "tag_safe") {
      out.append(sanitizeForFilename(e.tag));
    } else if (key == "kind") {
      out.append(gameEventKindName(e.kind));
    } else if (key == "msg") {
      out.append(e.msg);
    } else if (key == "msg_safe") {
      out.append(sanitizeForFilename(e.msg));
    } else if (key == "u64a") {
      appendNum(e.u64a);
    } else if (key == "u64b") {
      appendNum(e.u64b);
    } else if (key == "tRealMs") {
      const long long ms = (long long)std::llround(e.tRealSec * 1000.0);
      out += std::to_string(ms);
    } else if (key == "tRealSec") {
      std::ostringstream ss;
      ss.setf(std::ios::fixed);
      ss.precision(3);
      ss << e.tRealSec;
      out += ss.str();
    } else if (key == "tSimDays") {
      std::ostringstream ss;
      ss.setf(std::ios::fixed);
      ss.precision(6);
      ss << e.tSimDays;
      out += ss.str();
    } else if (key == "rule") {
      out.append(ruleName);
    } else {
      // Unknown key: keep literal.
      out.push_back('{');
      out.append(key);
      out.push_back('}');
    }

    i = j + 1;
  }

  return out;
}

static void setCStr(char* dst, std::size_t dstSz, std::string_view s) {
  if (!dst || dstSz == 0) return;
  std::snprintf(dst, dstSz, "%.*s", (int)std::min<std::size_t>(dstSz - 1, s.size()), s.data());
}

void initDefaultAutomationRules(IntegrationHubWindowState& st) {
  if (st.automationsInitialized) return;
  st.automationsInitialized = true;

  st.automationRules.clear();
  st.automationRules.reserve(12);

  // Starter rules are intentionally disabled by default.
  // They are meant as examples + quick toggles.

  {
    AutomationRule r;
    r.enabled = false;
    setCStr(r.name, sizeof(r.name), "Validation Watchdog -> Capture bundle");
    r.kind = GameEventKind::Validation;
    r.tagMatch = AutomationTagMatch::Equals;
    setCStr(r.tagText, sizeof(r.tagText), "Watchdog");
    r.cooldownSec = 0.25;

        // One-click bundle: screenshots + traces + repro JSON in a single folder.
    {
      AutomationActionTemplate a;
      a.kind = GameActionKind::CaptureBundle;
      a.i32Const = CaptureBundle_Default | CaptureBundle_PauseForScreenshots | CaptureBundle_StopFlightRecorder;
      setCStr(a.msgTemplate, sizeof(a.msgTemplate), "captures|watchdog_{tRealMs}_{tag_safe}");
      r.actions.push_back(a);
    }

    st.automationRules.push_back(std::move(r));
  }

  {
    AutomationRule r;
    r.enabled = false;
    setCStr(r.name, sizeof(r.name), "TimeTrial Finish -> Stop recorder");
    r.kind = GameEventKind::TimeTrial;
    r.tagMatch = AutomationTagMatch::Equals;
    setCStr(r.tagText, sizeof(r.tagText), "Finish");
    r.cooldownSec = 0.25;

    AutomationActionTemplate a;
    a.kind = GameActionKind::StopFlightRecorder;
    r.actions.push_back(a);

    st.automationRules.push_back(std::move(r));
  }

  {
    AutomationRule r;
    r.enabled = false;
    setCStr(r.name, sizeof(r.name), "TimeTrial Finish -> Export traces");
    r.kind = GameEventKind::TimeTrial;
    r.tagMatch = AutomationTagMatch::Equals;
    setCStr(r.tagText, sizeof(r.tagText), "Finish");
    r.cooldownSec = 0.25;

    AutomationActionTemplate a1;
    a1.kind = GameActionKind::ExportFlightRecorderTrace;
    a1.delaySec = 0.15f;
    setCStr(a1.msgTemplate, sizeof(a1.msgTemplate), "flight_trace_{tRealMs}_{tag_safe}.json");
    r.actions.push_back(a1);

    AutomationActionTemplate a2;
    a2.kind = GameActionKind::ExportIntegrationTrace;
    a2.delaySec = 0.20f;
    setCStr(a2.msgTemplate, sizeof(a2.msgTemplate), "integration_trace_{tRealMs}_{tag_safe}.json");
    r.actions.push_back(a2);

    st.automationRules.push_back(std::move(r));
  }

  {
    AutomationRule r;
    r.enabled = false;
    setCStr(r.name, sizeof(r.name), "Dock success -> Toast");
    r.kind = GameEventKind::Docking;
    r.tagMatch = AutomationTagMatch::Equals;
    setCStr(r.tagText, sizeof(r.tagText), "Dock");
    r.cooldownSec = 1.0;

    AutomationActionTemplate a;
    a.kind = GameActionKind::Toast;
    setCStr(a.msgTemplate, sizeof(a.msgTemplate), "Dock event: station={u64a} pad={u64b}");
    r.actions.push_back(a);

    st.automationRules.push_back(std::move(r));
  }

  {
    AutomationRule r;
    r.enabled = false;
    setCStr(r.name, sizeof(r.name), "Mission Complete -> Export integration trace");
    r.kind = GameEventKind::Gameplay;
    r.tagMatch = AutomationTagMatch::Equals;
    setCStr(r.tagText, sizeof(r.tagText), "MissionComplete");
    r.cooldownSec = 0.25;

    AutomationActionTemplate a;
    a.kind = GameActionKind::ExportIntegrationTrace;
    a.delaySec = 0.10f;
    setCStr(a.msgTemplate, sizeof(a.msgTemplate), "integration_mission_{u64a}_{tRealMs}.json");
    r.actions.push_back(a);

    st.automationRules.push_back(std::move(r));
  }

  {
    AutomationRule r;
    r.enabled = false;
    setCStr(r.name, sizeof(r.name), "Docking Clearance -> Comms");
    r.kind = GameEventKind::Docking;
    r.tagMatch = AutomationTagMatch::Equals;
    setCStr(r.tagText, sizeof(r.tagText), "Clearance");
    r.cooldownSec = 0.35;

    AutomationActionTemplate a;
    a.kind = GameActionKind::TransmitComms;
    a.i32Const = 0; // System
    a.bConst = true; // show overlay
    a.u64aSource = AutomationValueSource::EventU64a; // station id
    setCStr(a.msgTemplate, sizeof(a.msgTemplate),
            "[color #67B7FF]DOCKING CLEARANCE[/color]|[type cps=60 fade=0.04]{msg}[/type]");
    r.actions.push_back(a);

    st.automationRules.push_back(std::move(r));
  }

  {
    AutomationRule r;
    r.enabled = false;
    setCStr(r.name, sizeof(r.name), "Mission Accept -> Sync Nav");
    r.kind = GameEventKind::Gameplay;
    r.tagMatch = AutomationTagMatch::Equals;
    setCStr(r.tagText, sizeof(r.tagText), "MissionAccepted");
    r.cooldownSec = 0.5;

    AutomationActionTemplate a;
    a.kind = GameActionKind::SyncNavToMission;
    a.u64aSource = AutomationValueSource::EventU64a; // mission id
    a.u64aConst = 0;
    a.i32Const = 0;
    a.bConst = false; // arm auto-run
    setCStr(a.msgTemplate, sizeof(a.msgTemplate), "");
    r.actions.push_back(a);

    st.automationRules.push_back(std::move(r));
  }

  {
    AutomationRule r;
    r.enabled = false;
    setCStr(r.name, sizeof(r.name), "Mission Accept -> Track mission");
    r.kind = GameEventKind::Gameplay;
    r.tagMatch = AutomationTagMatch::Equals;
    setCStr(r.tagText, sizeof(r.tagText), "MissionAccepted");
    r.cooldownSec = 0.5;

    AutomationActionTemplate a;
    a.kind = GameActionKind::SetTrackedMissionId;
    a.u64aSource = AutomationValueSource::EventU64a; // mission id
    a.u64aConst = 0;
    r.actions.push_back(a);

    st.automationRules.push_back(std::move(r));
  }

  {
    AutomationRule r;
    r.enabled = false;
    setCStr(r.name, sizeof(r.name), "Mission Complete -> Comms");
    r.kind = GameEventKind::Gameplay;
    r.tagMatch = AutomationTagMatch::Equals;
    setCStr(r.tagText, sizeof(r.tagText), "MissionComplete");
    r.cooldownSec = 0.5;

    AutomationActionTemplate a;
    a.kind = GameActionKind::TransmitComms;
    a.i32Const = 1; // Mission
    a.bConst = true;
    a.u64aSource = AutomationValueSource::Constant;
    a.u64aConst = 0;
    setCStr(a.msgTemplate, sizeof(a.msgTemplate),
            "[color #67B7FF]MISSION COMPLETE[/color]|[type cps=72 fade=0.04]{msg}[/type]");
    r.actions.push_back(a);

    st.automationRules.push_back(std::move(r));
  }

  {
    AutomationRule r;
    r.enabled = false;
    setCStr(r.name, sizeof(r.name), "TimeTrial Finish -> Comms");
    r.kind = GameEventKind::TimeTrial;
    r.tagMatch = AutomationTagMatch::Equals;
    setCStr(r.tagText, sizeof(r.tagText), "Finish");
    r.cooldownSec = 0.5;

    AutomationActionTemplate a;
    a.kind = GameActionKind::TransmitComms;
    a.i32Const = 0; // System
    a.bConst = true;
    a.u64aSource = AutomationValueSource::Constant;
    a.u64aConst = 0;
    setCStr(a.msgTemplate, sizeof(a.msgTemplate),
            "[color #9BA6B4]RACE CONTROL[/color]|[type cps=68 fade=0.04]{msg}[/type]");
    r.actions.push_back(a);

    st.automationRules.push_back(std::move(r));
  }

  // --- Nav auto-run phase hooks (disabled by default) ---
  // These rules demonstrate how the nav auto-run state machine can drive other systems.
  {
    AutomationRule r;
    r.enabled = false;
    setCStr(r.name, sizeof(r.name), "AutoRun Supercruise -> CameraPreset Travel");
    r.kind = GameEventKind::NavAssist;
    r.tagMatch = AutomationTagMatch::Equals;
    setCStr(r.tagText, sizeof(r.tagText), "AutoRunSupercruise");
    r.cooldownSec = 0.15;

    AutomationActionTemplate a;
    a.kind = GameActionKind::SetCameraRigPreset;
    a.i32Const = 1; // Travel preset
    r.actions.push_back(a);

    st.automationRules.push_back(std::move(r));
  }

  {
    AutomationRule r;
    r.enabled = false;
    setCStr(r.name, sizeof(r.name), "AutoRun DockingComputer -> CameraPreset Docking");
    r.kind = GameEventKind::NavAssist;
    r.tagMatch = AutomationTagMatch::Equals;
    setCStr(r.tagText, sizeof(r.tagText), "AutoRunDockingComputer");
    r.cooldownSec = 0.15;

    AutomationActionTemplate a;
    a.kind = GameActionKind::SetCameraRigPreset;
    a.i32Const = 2; // Docking preset
    r.actions.push_back(a);

    st.automationRules.push_back(std::move(r));
  }

  {
    AutomationRule r;
    r.enabled = false;
    setCStr(r.name, sizeof(r.name), "AutoRun Stop -> Comms");
    r.kind = GameEventKind::NavAssist;
    r.tagMatch = AutomationTagMatch::Prefix;
    setCStr(r.tagText, sizeof(r.tagText), "AutoRunStop");
    r.cooldownSec = 0.25;

    AutomationActionTemplate a;
    a.kind = GameActionKind::TransmitComms;
    a.i32Const = 0; // System
    a.bConst = true;
    a.u64aSource = AutomationValueSource::Constant;
    a.u64aConst = 0;
    a.delaySec = 0.05f;
    setCStr(a.msgTemplate, sizeof(a.msgTemplate),
            "[color #9BA6B4]NAV ASSIST[/color]|[type cps=72 fade=0.04]{msg}[/type]");
    r.actions.push_back(a);

    st.automationRules.push_back(std::move(r));
  }

  // --- Camera presets (disabled by default) ---
  // These are intentionally *examples* of cross-system integration:
  // docking computer events drive the camera rig preset via the automation engine.
  {
    AutomationRule r;
    r.enabled = false;
    setCStr(r.name, sizeof(r.name), "DockingComputer Engaged -> CameraPreset Docking");
    r.kind = GameEventKind::Docking;
    r.tagMatch = AutomationTagMatch::Equals;
    setCStr(r.tagText, sizeof(r.tagText), "DockingComputerEngaged");
    r.cooldownSec = 0.2;

    AutomationActionTemplate a;
    a.kind = GameActionKind::SetCameraRigPreset;
    a.i32Const = 2; // Docking preset
    r.actions.push_back(a);

    st.automationRules.push_back(std::move(r));
  }

  {
    AutomationRule r;
    r.enabled = false;
    setCStr(r.name, sizeof(r.name), "DockingComputer Disengaged -> CameraPreset DefaultOrbit");
    r.kind = GameEventKind::Docking;
    r.tagMatch = AutomationTagMatch::Equals;
    setCStr(r.tagText, sizeof(r.tagText), "DockingComputerDisengaged");
    r.cooldownSec = 0.2;

    AutomationActionTemplate a;
    a.kind = GameActionKind::SetCameraRigPreset;
    a.i32Const = 0; // DefaultOrbit preset
    r.actions.push_back(a);

    st.automationRules.push_back(std::move(r));
  }


}

void applyAutomationRules(IntegrationHubWindowState& st, const GameEvent& ev) {
  if (!st.automationsEnabled) return;
  if (st.automationRules.empty()) return;

  if (st.automationsRecursionDepth > 3) return;
  struct DepthGuard {
    IntegrationHubWindowState& st;
    explicit DepthGuard(IntegrationHubWindowState& s) : st(s) { st.automationsRecursionDepth++; }
    ~DepthGuard() { st.automationsRecursionDepth--; }
  } guard(st);

  for (auto& r : st.automationRules) {
    if (!r.enabled) continue;
    if (r.kind != ev.kind) continue;
    if (!tagMatches(r, ev.tag)) continue;

    const double dt = ev.tRealSec - r.lastFiredRealSec;
    if (dt < r.cooldownSec) continue;

    r.lastFiredRealSec = ev.tRealSec;

    const std::string origin = std::string("Automation:") + r.name;

    for (const auto& t : r.actions) {
      if (st.automationsActionsFiredThisFrame >= st.automationsMaxActionsPerFrame) return;

      const double delaySec = std::max(0.0, (double)t.delaySec);

      GameAction a;
      a.tRealSec = ev.tRealSec + delaySec;
      a.tSimDays = ev.tSimDays + (delaySec / 86400.0);
      a.origin = origin;
      a.kind = t.kind;

      switch (t.kind) {
        case GameActionKind::SetTargetStationId:
          a.u64a = resolveU64(t.u64aSource, t.u64aConst, ev);
          break;
        case GameActionKind::EngageDockingComputer:
          a.b = t.bConst;
          break;
        case GameActionKind::SyncNavToMission:
          a.u64a = resolveU64(t.u64aSource, t.u64aConst, ev);
          a.b = t.bConst;
          break;
        case GameActionKind::SetTrackedMissionId:
          a.u64a = resolveU64(t.u64aSource, t.u64aConst, ev);
          break;
        case GameActionKind::GoToStation:
          a.u64a = resolveU64(t.u64aSource, t.u64aConst, ev);
          a.u64b = resolveU64(t.u64bSource, t.u64bConst, ev);
          a.b = t.bConst;
          break;
        case GameActionKind::SetCameraRigMode:
          a.i32a = t.i32Const;
          break;
        case GameActionKind::SetCameraRigPreset:
          a.i32a = t.i32Const;
          break;

        case GameActionKind::TransmitComms: {
          a.i32a = t.i32Const;
          a.b = t.bConst;
          a.u64a = resolveU64(t.u64aSource, t.u64aConst, ev);
          a.msg = expandTemplate(t.msgTemplate, ev, r.name);
          break;
        }
        case GameActionKind::RequestScreenshot:
          a.b = t.bConst;
          a.i32a = t.i32Const;
          a.msg = expandTemplate(t.msgTemplate, ev, r.name);
          break;

        case GameActionKind::CaptureBundle:
          a.i32a = t.i32Const;
          a.msg = expandTemplate(t.msgTemplate, ev, r.name);
          break;

        case GameActionKind::Toast:
        case GameActionKind::ExportFlightRecorderTrace:
        case GameActionKind::ExportIntegrationTrace:
          a.msg = expandTemplate(t.msgTemplate, ev, r.name);
          break;
        case GameActionKind::StartFlightRecorder:
        case GameActionKind::StopFlightRecorder:
        case GameActionKind::ClearIntegrationHub:
        default:
          break;
      }

      hubScheduleAction(st, std::move(a));
      st.automationsActionsFiredThisFrame++;
    }
  }
}

static stellar::core::u64 nowEpochMs() {
  using namespace std::chrono;
  return (stellar::core::u64)duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

bool writeIntegrationTraceJson(const IntegrationHubWindowState& st, const char* path, std::string* outErr) {
  // Use atomicWriteFile to avoid partially-written JSON being consumed by external tools.
  const std::filesystem::path outPath(path);

  stellar::core::JsonWriter::Options opt;
  opt.pretty = st.exportPretty;

  stellar::core::JsonWriter w(opt);

  w.beginObject();
    w.key("type"); w.value("stellar_integration_trace");
    w.key("version"); w.value(4);

    // Small metadata block to make traces self-describing even when consumers only
    // stream-read the JSON.
    w.object("metadata", [&]() {
      w.field("generatedAtEpochMs", (double)nowEpochMs());
      w.field("exportPretty", st.exportPretty);
      w.field("pendingCount", (double)st.actions.pending.size());
      w.field("scheduledCount", (double)st.scheduledActions.size());
      w.field("historyCount", (double)st.actions.history.size());
      w.field("eventCount", (double)st.events.events.size());
      w.field("automationRuleCount", (double)st.automationRules.size());

      // Lightweight session stats (useful for judging if a trace was truncated by ring buffers).
      w.field("totalEventsPushed", (double)st.totalEventsPushed);
      w.field("totalActionsPushed", (double)st.totalActionsPushed);
      w.field("totalActionsScheduled", (double)st.totalActionsScheduled);

      // Embed current timestamps so actions/events can be correlated to wall clock.
      w.field("nowRealSec", st.nowRealSec);
      w.field("nowSimDays", st.nowSimDays);
    });

    if (st.exportIncludePendingActions) {
      w.key("pendingActions");
      w.beginArray();
      for (const auto& a : st.actions.pending) jsonAction(w, a);
      w.endArray();
    }

    if (st.exportIncludeScheduledActions) {
      w.key("scheduledActions");
      w.beginArray();
      for (const auto& a : st.scheduledActions) jsonAction(w, a);
      w.endArray();
    }

    if (st.exportIncludeActionHistory) {
      w.key("actionHistory");
      w.beginArray();
      for (const auto& a : st.actions.history) jsonAction(w, a);
      w.endArray();
    }

    if (st.exportIncludeEvents) {
      w.key("events");
      w.beginArray();
      for (const auto& e : st.events.events) jsonEvent(w, e);
      w.endArray();
    }

    if (st.exportIncludeAutomationRules) {
      w.key("automationRules");
      w.beginArray();
      for (const auto& r : st.automationRules) jsonAutomationRule(w, r);
      w.endArray();
    }

  w.endObject();

  const auto json = w.takeString();
  return stellar::core::atomicWriteFile(
      outPath,
      [&](std::ostream& out, std::string* writeErr) -> bool {
        (void)writeErr;
        out.write(json.data(), static_cast<std::streamsize>(json.size()));
        return out.good();
      },
      outErr);
}

void drawIntegrationHubWindow(IntegrationHubWindowState& st, const ToastFn& toast, const IntegrationHubUiHooks* hooks) {
  if (!st.open) return;

  if (!st.automationsInitialized) initDefaultAutomationRules(st);

  if (!ImGui::Begin("Integration Hub", &st.open)) {
    ImGui::End();
    return;
  }

  ImGui::TextUnformatted(
      "Cross-system actions + event timeline (for debugging and repro packs).\n"
      "Now includes a tiny rule engine to translate events -> actions (automation tab)."
  );

  ImGui::Separator();

  if (ImGui::BeginTabBar("##IntegrationHubTabs")) {
    if (ImGui::BeginTabItem("Actions")) {
      ImGui::Checkbox("Show scheduled", &st.showScheduledActions);
      ImGui::SameLine();
      ImGui::Checkbox("Show pending", &st.showPendingActions);
      ImGui::SameLine();
      ImGui::Checkbox("Show history", &st.showActionHistory);
      ImGui::SameLine();
      ImGui::Checkbox("Auto-scroll", &st.autoScroll);

      ImGui::Text("Scheduled: %d  |  Pending: %d  |  History: %d", (int)st.scheduledActions.size(), (int)st.actions.pending.size(), (int)st.actions.history.size());

      if (ImGui::Button("Clear scheduled")) {
        st.scheduledActions.clear();
      }
      ImGui::SameLine();
      if (ImGui::Button("Clear pending")) {
        st.actions.pending.clear();
      }
      ImGui::SameLine();
      if (ImGui::Button("Clear history")) {
        st.actions.history.clear();
      }
      ImGui::SameLine();
      if (ImGui::Button("Copy scheduled as text")) {
        std::ostringstream ss;
        for (const auto& a : st.scheduledActions) ss << formatActionLine(a) << "\n";
        copyTextToClipboard(ss.str());
        if (toast) toast("Copied scheduled actions to clipboard.", 1.6);
      }
      ImGui::SameLine();
      if (ImGui::Button("Copy pending as text")) {
        std::ostringstream ss;
        for (const auto& a : st.actions.pending) ss << formatActionLine(a) << "\n";
        copyTextToClipboard(ss.str());
        if (toast) toast("Copied pending actions to clipboard.", 1.6);
      }

      ImGui::Separator();

      const float childH = ImGui::GetContentRegionAvail().y;
      if (ImGui::BeginChild("##ActionList", ImVec2(0, childH), true, ImGuiWindowFlags_HorizontalScrollbar)) {
        auto drawList = [&](const std::deque<GameAction>& list, const char* header, bool scheduled) {
          if (!list.empty()) {
            ImGui::TextDisabled("%s", header);
            ImGui::Separator();
          }
          for (const auto& a : list) {
            std::string line = formatActionLine(a);
            if (scheduled) {
              const double dt = a.tRealSec - st.nowRealSec;
              char buf[64];
              if (dt > 0.001) {
                std::snprintf(buf, sizeof(buf), "  (in +%.2fs)", dt);
              } else {
                std::snprintf(buf, sizeof(buf), "  (due)");
              }
              line += buf;
            }
            ImGui::TextUnformatted(line.c_str());
          }
        };

        if (st.showScheduledActions) drawList(st.scheduledActions, "Scheduled", true);
        if (st.showPendingActions) {
          if (st.showScheduledActions && !st.actions.pending.empty()) ImGui::Separator();
          drawList(st.actions.pending, "Pending", false);
        }
        if (st.showActionHistory) {
          if (st.showPendingActions && !st.actions.history.empty()) ImGui::Separator();
          drawList(st.actions.history, "History", false);
        }

        if (st.autoScroll && ImGui::GetScrollY() >= ImGui::GetScrollMaxY() - 2.0f) {
          ImGui::SetScrollHereY(1.0f);
        }
      }
      ImGui::EndChild();

      ImGui::EndTabItem();
    }

    if (ImGui::BeginTabItem("Events")) {
      ImGui::Checkbox("Auto-scroll", &st.autoScroll);
      ImGui::SameLine();
      if (ImGui::Button("Clear")) {
        st.events.clear();
        st.selectedEventIndex = -1;
      }
      ImGui::SameLine();
      if (ImGui::Button("Copy filtered as text")) {
        std::ostringstream ss;
        for (const auto& e : st.events.events) {
          if (!passesFilter(st, e)) continue;
          ss << formatEventLine(e) << "\n";
        }
        copyTextToClipboard(ss.str());
        if (toast) toast("Copied events to clipboard.", 1.6);
      }

      ImGui::Separator();

      ImGui::TextDisabled("Filters:");
      ImGui::Checkbox("Info", &st.showInfo); ImGui::SameLine();
      ImGui::Checkbox("TimeTrial", &st.showTimeTrial); ImGui::SameLine();
      ImGui::Checkbox("Docking", &st.showDocking); ImGui::SameLine();
      ImGui::Checkbox("Camera", &st.showCamera);

      ImGui::Checkbox("NavAssist", &st.showNavAssist); ImGui::SameLine();
      ImGui::Checkbox("Validation", &st.showValidation); ImGui::SameLine();
      ImGui::Checkbox("FlightRecorder", &st.showFlightRecorder); ImGui::SameLine();
      ImGui::Checkbox("Gameplay", &st.showGameplay); ImGui::SameLine();
      ImGui::Checkbox("Debug", &st.showDebug);

      ImGui::Separator();

      // Two-pane layout: selectable timeline list + event inspector.
      const float childH = ImGui::GetContentRegionAvail().y;
      const float availW = ImGui::GetContentRegionAvail().x;
      const float kMinPaneW = 240.0f;
      const float kSplitterW = 6.0f;
      st.eventListWidth = std::clamp(st.eventListWidth, kMinPaneW, std::max(kMinPaneW, availW - kMinPaneW - kSplitterW));

      // Keep the selection valid if the ring buffer drops old events.
      if (st.selectedEventIndex >= (int)st.events.events.size()) st.selectedEventIndex = -1;

      // Left: list
      if (ImGui::BeginChild("##EventList", ImVec2(st.eventListWidth, childH), true, ImGuiWindowFlags_HorizontalScrollbar)) {
        int idx = 0;
        for (const auto& e : st.events.events) {
          const int thisIdx = idx;
          idx++;

          if (!passesFilter(st, e)) continue;

          const std::string line = formatEventLine(e);
          const bool selected = (st.selectedEventIndex == thisIdx);
          if (ImGui::Selectable(line.c_str(), selected, ImGuiSelectableFlags_AllowDoubleClick)) {
            st.selectedEventIndex = thisIdx;

            // Double-click is a quick "focus" gesture.
            if (ImGui::IsMouseDoubleClicked(0) && hooks && hooks->focusCameraAtPosKm && e.hasPos) {
              hooks->focusCameraAtPosKm(e.posKm);
              if (toast) toast("Camera focused on event position.", 1.4);
            }
          }
        }
        if (st.autoScroll && ImGui::GetScrollY() >= ImGui::GetScrollMaxY() - 2.0f) {
          ImGui::SetScrollHereY(1.0f);
        }
      }
      ImGui::EndChild();

      ImGui::SameLine();
      ImGui::InvisibleButton("##EventSplitter", ImVec2(kSplitterW, childH));
      if (ImGui::IsItemActive()) {
        st.eventListWidth += ImGui::GetIO().MouseDelta.x;
      }
      if (ImGui::IsItemHovered()) {
        ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeEW);
      }
      ImGui::SameLine();

      // Right: inspector
      if (ImGui::BeginChild("##EventInspector", ImVec2(0, childH), true)) {
        ImGui::TextDisabled("Event Inspector");
        ImGui::Separator();

        if (st.selectedEventIndex < 0 || st.selectedEventIndex >= (int)st.events.events.size()) {
          ImGui::TextWrapped("Select an event from the timeline.");
          ImGui::TextDisabled("Tip: double-click an event with a position to focus the camera.");
        } else {
          const GameEvent& e = st.events.events[st.selectedEventIndex];

          ImGui::Text("tReal: %.3fs", e.tRealSec);
          ImGui::Text("tSim : %.6f days", e.tSimDays);
          ImGui::Text("Kind : %s", gameEventKindName(e.kind));
          ImGui::Text("Tag  : %s", e.tag.c_str());

          if (!e.msg.empty()) {
            ImGui::Separator();
            ImGui::TextWrapped("%s", e.msg.c_str());
          }

          if (e.hasPos) {
            ImGui::Separator();
            ImGui::Text("posKm: [%.3f, %.3f, %.3f]", e.posKm.x, e.posKm.y, e.posKm.z);
          }

          if (e.u64a || e.u64b) {
            ImGui::Text("u64a : %llu", (unsigned long long)e.u64a);
            ImGui::Text("u64b : %llu", (unsigned long long)e.u64b);
          }

          ImGui::Separator();

          if (ImGui::Button("Copy line")) {
            copyTextToClipboard(formatEventLine(e));
            if (toast) toast("Copied event line.", 1.3);
          }
          ImGui::SameLine();
          if (ImGui::Button("Copy JSON")) {
            stellar::core::JsonWriter::Options opt;
            opt.pretty = true;
            stellar::core::JsonWriter w(opt);
            jsonEvent(w, e);
            auto json = w.takeString();
            copyTextToClipboard(json);
            if (toast) toast("Copied event JSON.", 1.3);
          }

          // Cross-system hooks.
          if (hooks && hooks->focusCameraAtPosKm && e.hasPos) {
            if (ImGui::Button("Focus camera")) {
              hooks->focusCameraAtPosKm(e.posKm);
              if (toast) toast("Camera focused on event position.", 1.4);
            }
          } else {
            ImGui::BeginDisabled();
            ImGui::Button("Focus camera");
            ImGui::EndDisabled();
          }

          ImGui::SameLine();
          if (hooks && hooks->scrubFlightRecorderToRealTimeSec) {
            if (ImGui::Button("Scrub flight recorder")) {
              hooks->scrubFlightRecorderToRealTimeSec(e.tRealSec);
              if (toast) toast("Flight recorder playhead moved.", 1.4);
            }
          } else {
            ImGui::BeginDisabled();
            ImGui::Button("Scrub flight recorder");
            ImGui::EndDisabled();
          }

          // Quick capture: let the timeline drive the screenshot pipeline via the Integration Hub action queue.
          {
            const std::string base = sanitizeForFilename(std::string("hub_") + gameEventKindName(e.kind) + (e.tag.empty() ? std::string() : (std::string("_") + e.tag)));
            const int flags = ScreenshotFlag_Timestamp | ScreenshotFlag_PauseForCapture;

            ImGui::SameLine();
            if (ImGui::SmallButton("Shot UI")) {
              hubPushAction(st, makeActionRequestScreenshot(st.nowRealSec, st.nowSimDays, "IntegrationHub", /*includeUi=*/true, flags, "screenshots", base));
              if (toast) toast("Screenshot requested (UI).", 1.3);
            }
            ImGui::SameLine();
            if (ImGui::SmallButton("Shot World")) {
              hubPushAction(st, makeActionRequestScreenshot(st.nowRealSec, st.nowSimDays, "IntegrationHub", /*includeUi=*/false, flags, "screenshots", base));
              if (toast) toast("Screenshot requested (world).", 1.3);
            }

            ImGui::SameLine();
            if (ImGui::SmallButton("Bundle")) {
              const int bFlags = CaptureBundle_Default | CaptureBundle_PauseForScreenshots;
              hubPushAction(st, makeActionCaptureBundle(st.nowRealSec, st.nowSimDays, "IntegrationHub", bFlags, "captures", base));
              if (toast) toast("Capture bundle requested.", 1.5);
            }
          }

          ImGui::Separator();

          if (ImGui::Button("Create automation rule from this event")) {
            AutomationRule r;
            r.enabled = true;
            r.kind = e.kind;
            r.tagMatch = AutomationTagMatch::Equals;
            r.cooldownSec = 0.25;

            std::string name = std::string("On ") + gameEventKindName(e.kind);
            if (!e.tag.empty()) name += std::string(":") + e.tag;
            setCStr(r.name, sizeof(r.name), name);
            setCStr(r.tagText, sizeof(r.tagText), e.tag);

            AutomationActionTemplate a;
            a.kind = GameActionKind::Toast;
            setCStr(a.msgTemplate, sizeof(a.msgTemplate), "Automation: {kind}:{tag} {msg}");
            r.actions.push_back(a);

            st.automationRules.push_back(std::move(r));
            if (toast) toast("Automation rule created (see Automation tab).", 2.0);
          }
        }
      }
      ImGui::EndChild();

      ImGui::EndTabItem();
    }

    if (ImGui::BeginTabItem("Automation")) {
      ImGui::Checkbox("Enable automations", &st.automationsEnabled);
      ImGui::SameLine();
      ImGui::TextDisabled("Fired this frame: %d / %d", st.automationsActionsFiredThisFrame, st.automationsMaxActionsPerFrame);

      ImGui::SliderInt("Max actions/frame", &st.automationsMaxActionsPerFrame, 0, 256);

      if (ImGui::Button("Reset starter rules")) {
        st.automationsInitialized = false;
        initDefaultAutomationRules(st);
        if (toast) toast("Starter automation rules reset.", 1.8);
      }
      ImGui::SameLine();
      if (ImGui::Button("Add blank rule")) {
        AutomationRule r;
        r.enabled = true;
        setCStr(r.name, sizeof(r.name), "New rule");
        r.kind = GameEventKind::Info;
        r.tagMatch = AutomationTagMatch::Any;
        r.cooldownSec = 0.25;
        st.automationRules.push_back(std::move(r));
      }

      ImGui::Separator();

      static const GameEventKind kEventKinds[] = {
          GameEventKind::Info,
          GameEventKind::TimeTrial,
          GameEventKind::Docking,
          GameEventKind::Camera,
          GameEventKind::NavAssist,
          GameEventKind::Validation,
          GameEventKind::FlightRecorder,
          GameEventKind::Gameplay,
          GameEventKind::Debug,
      };
      static const char* kEventKindLabels[] = {
          "Info", "TimeTrial", "Docking", "Camera", "NavAssist", "Validation", "FlightRecorder", "Gameplay", "Debug",
      };

      static const AutomationTagMatch kTagMatches[] = {
          AutomationTagMatch::Any,
          AutomationTagMatch::Equals,
          AutomationTagMatch::Contains,
          AutomationTagMatch::Prefix,
          AutomationTagMatch::Suffix,
      };
      static const char* kTagMatchLabels[] = {
          "Any", "Equals", "Contains", "Prefix", "Suffix",
      };

      static const AutomationValueSource kU64Sources[] = {
          AutomationValueSource::Constant,
          AutomationValueSource::EventU64a,
          AutomationValueSource::EventU64b,
      };
      static const char* kU64SourceLabels[] = {
          "Constant", "Event.u64a", "Event.u64b",
      };

      static const GameActionKind kActionKinds[] = {
          GameActionKind::Toast,
          GameActionKind::SetTargetStationId,
          GameActionKind::EngageDockingComputer,
          GameActionKind::SyncNavToMission,
          GameActionKind::SetTrackedMissionId,
          GameActionKind::GoToStation,
          GameActionKind::SetCameraRigMode,
          GameActionKind::SetCameraRigPreset,
          GameActionKind::RequestScreenshot,
          GameActionKind::CaptureBundle,
          GameActionKind::StartFlightRecorder,
          GameActionKind::StopFlightRecorder,
          GameActionKind::ExportFlightRecorderTrace,
          GameActionKind::ExportIntegrationTrace,
          GameActionKind::TransmitComms,
          GameActionKind::ClearIntegrationHub,
      };
      static const char* kActionKindLabels[] = {
          "Toast",
          "SetTargetStationId",
          "EngageDockingComputer",
          "SyncNavToMission",
          "SetTrackedMissionId",
          "GoToStation",
          "SetCameraRigMode",
          "SetCameraRigPreset",
          "RequestScreenshot",
          "CaptureBundle",
          "StartFlightRecorder",
          "StopFlightRecorder",
          "ExportFlightRecorderTrace",
          "ExportIntegrationTrace",
          "TransmitComms",
          "ClearIntegrationHub",
      };

      for (std::size_t ri = 0; ri < st.automationRules.size(); ++ri) {
        AutomationRule& r = st.automationRules[ri];

        ImGui::PushID((int)ri);

        // Header
        const bool open = ImGui::TreeNodeEx("Rule", ImGuiTreeNodeFlags_DefaultOpen, "%s%s", r.enabled ? "[on] " : "[off] ", r.name);
        ImGui::SameLine();
        if (ImGui::SmallButton("Delete")) {
          st.automationRules.erase(st.automationRules.begin() + (long long)ri);
          ImGui::PopID();
          break;
        }

        if (open) {
          ImGui::Checkbox("Enabled", &r.enabled);
          ImGui::InputText("Name", r.name, IM_ARRAYSIZE(r.name));

          int kindIdx = 0;
          for (int i = 0; i < (int)(sizeof(kEventKinds) / sizeof(kEventKinds[0])); ++i) {
            if (kEventKinds[i] == r.kind) { kindIdx = i; break; }
          }
          if (ImGui::Combo("Event kind", &kindIdx, kEventKindLabels, IM_ARRAYSIZE(kEventKindLabels))) {
            r.kind = kEventKinds[kindIdx];
          }

          int matchIdx = 0;
          for (int i = 0; i < (int)(sizeof(kTagMatches) / sizeof(kTagMatches[0])); ++i) {
            if (kTagMatches[i] == r.tagMatch) { matchIdx = i; break; }
          }
          if (ImGui::Combo("Tag match", &matchIdx, kTagMatchLabels, IM_ARRAYSIZE(kTagMatchLabels))) {
            r.tagMatch = kTagMatches[matchIdx];
          }

          ImGui::InputText("Tag text", r.tagText, IM_ARRAYSIZE(r.tagText));
          ImGui::SliderFloat("Cooldown (sec)", (float*)&r.cooldownSec, 0.0f, 10.0f, "%.2f");

          ImGui::Separator();
          ImGui::TextDisabled("Actions:");

          for (std::size_t ai = 0; ai < r.actions.size(); ++ai) {
            AutomationActionTemplate& a = r.actions[ai];
            ImGui::PushID((int)ai);

            int actIdx = 0;
            for (int i = 0; i < (int)(sizeof(kActionKinds) / sizeof(kActionKinds[0])); ++i) {
              if (kActionKinds[i] == a.kind) { actIdx = i; break; }
            }
            const GameActionKind prevKind = a.kind;
            if (ImGui::Combo("Kind", &actIdx, kActionKindLabels, IM_ARRAYSIZE(kActionKindLabels))) {
              const GameActionKind newKind = kActionKinds[actIdx];
              if (newKind != prevKind) {
                a.kind = newKind;

                // Helpful defaults for newly-selected kinds.
                if (a.kind == GameActionKind::RequestScreenshot) {
                  a.bConst = true; // include UI by default
                  a.i32Const = ScreenshotFlag_Timestamp | ScreenshotFlag_PauseForCapture;
                  if (a.msgTemplate[0] == 0) {
                    setCStr(a.msgTemplate, sizeof(a.msgTemplate), "screenshots|shot_{tag_safe}_{tRealMs}");
                  }
                } else if (a.kind == GameActionKind::CaptureBundle) {
                  a.i32Const = CaptureBundle_Default | CaptureBundle_PauseForScreenshots;
                  if (a.msgTemplate[0] == 0) {
                    setCStr(a.msgTemplate, sizeof(a.msgTemplate), "captures|bundle_{tag_safe}_{tRealMs}");
                  }
                }
              }
            }

            // Delay (sec): schedule this action relative to the triggering event.
            ImGui::SliderFloat("Delay (sec)", &a.delaySec, 0.0f, 10.0f, "%.2f");

            if (a.kind == GameActionKind::SetTargetStationId) {
              int srcIdx = 0;
              for (int i = 0; i < (int)(sizeof(kU64Sources) / sizeof(kU64Sources[0])); ++i) {
                if (kU64Sources[i] == a.u64aSource) { srcIdx = i; break; }
              }
              if (ImGui::Combo("Station id source", &srcIdx, kU64SourceLabels, IM_ARRAYSIZE(kU64SourceLabels))) {
                a.u64aSource = kU64Sources[srcIdx];
              }
              if (a.u64aSource == AutomationValueSource::Constant) {
                unsigned long long v = (unsigned long long)a.u64aConst;
                if (ImGui::InputScalar("Station id", ImGuiDataType_U64, &v)) {
                  a.u64aConst = (stellar::core::u64)v;
                }
              } else {
                ImGui::TextDisabled("Using %s", automationValueSourceName(a.u64aSource));
              }
            } else if (a.kind == GameActionKind::EngageDockingComputer) {
              ImGui::Checkbox("Engage", &a.bConst);
            } else if (a.kind == GameActionKind::SyncNavToMission) {
              int srcIdx = 0;
              for (int i = 0; i < (int)(sizeof(kU64Sources) / sizeof(kU64Sources[0])); ++i) {
                if (kU64Sources[i] == a.u64aSource) { srcIdx = i; break; }
              }
              if (ImGui::Combo("Mission id source", &srcIdx, kU64SourceLabels, IM_ARRAYSIZE(kU64SourceLabels))) {
                a.u64aSource = kU64Sources[srcIdx];
              }
              if (a.u64aSource == AutomationValueSource::Constant) {
                unsigned long long v = (unsigned long long)a.u64aConst;
                if (ImGui::InputScalar("Mission id", ImGuiDataType_U64, &v)) {
                  a.u64aConst = (stellar::core::u64)v;
                }
              } else {
                ImGui::TextDisabled("Using %s", automationValueSourceName(a.u64aSource));
              }

              ImGui::Checkbox("Arm auto-run", &a.bConst);
              ImGui::TextDisabled("Syncs the nav route + arrival target to the mission's next stop.");
            } else if (a.kind == GameActionKind::SetTrackedMissionId) {
              int srcIdx = 0;
              for (int i = 0; i < (int)(sizeof(kU64Sources) / sizeof(kU64Sources[0])); ++i) {
                if (kU64Sources[i] == a.u64aSource) { srcIdx = i; break; }
              }
              if (ImGui::Combo("Mission id source", &srcIdx, kU64SourceLabels, IM_ARRAYSIZE(kU64SourceLabels))) {
                a.u64aSource = kU64Sources[srcIdx];
              }
              if (a.u64aSource == AutomationValueSource::Constant) {
                unsigned long long v = (unsigned long long)a.u64aConst;
                if (ImGui::InputScalar("Mission id", ImGuiDataType_U64, &v)) {
                  a.u64aConst = (stellar::core::u64)v;
                }
              } else {
                ImGui::TextDisabled("Using %s", automationValueSourceName(a.u64aSource));
              }
              ImGui::TextDisabled("Sets the tracked mission shown in the Objective HUD (0 clears).");
            } else if (a.kind == GameActionKind::GoToStation) {
              // u64a: station id
              int aSrcIdx = 0;
              for (int i = 0; i < (int)(sizeof(kU64Sources) / sizeof(kU64Sources[0])); ++i) {
                if (kU64Sources[i] == a.u64aSource) { aSrcIdx = i; break; }
              }
              if (ImGui::Combo("Station id source", &aSrcIdx, kU64SourceLabels, IM_ARRAYSIZE(kU64SourceLabels))) {
                a.u64aSource = kU64Sources[aSrcIdx];
              }
              if (a.u64aSource == AutomationValueSource::Constant) {
                unsigned long long v = (unsigned long long)a.u64aConst;
                if (ImGui::InputScalar("Station id", ImGuiDataType_U64, &v)) {
                  a.u64aConst = (stellar::core::u64)v;
                }
              } else {
                ImGui::TextDisabled("Station id uses %s", automationValueSourceName(a.u64aSource));
              }

              // u64b: system id
              int bSrcIdx = 0;
              for (int i = 0; i < (int)(sizeof(kU64Sources) / sizeof(kU64Sources[0])); ++i) {
                if (kU64Sources[i] == a.u64bSource) { bSrcIdx = i; break; }
              }
              if (ImGui::Combo("System id source", &bSrcIdx, kU64SourceLabels, IM_ARRAYSIZE(kU64SourceLabels))) {
                a.u64bSource = kU64Sources[bSrcIdx];
              }
              if (a.u64bSource == AutomationValueSource::Constant) {
                unsigned long long v = (unsigned long long)a.u64bConst;
                if (ImGui::InputScalar("System id", ImGuiDataType_U64, &v)) {
                  a.u64bConst = (stellar::core::u64)v;
                }
              } else {
                ImGui::TextDisabled("System id uses %s", automationValueSourceName(a.u64bSource));
              }

              ImGui::Checkbox("Arm auto-run", &a.bConst);
              ImGui::TextDisabled("Plots route / targets the station, then optionally auto-jumps/supercruises + docks.");
            } else if (a.kind == GameActionKind::SetCameraRigMode) {
              ImGui::SliderInt("Mode (0=Chase,1=Orbit)", &a.i32Const, 0, 1);
            } else if (a.kind == GameActionKind::SetCameraRigPreset) {
              static const char* kPresetLabels[] = {"DefaultOrbit", "Travel", "Docking", "Combat", "Cinematic"};
              int p = std::clamp(a.i32Const, 0, (int)IM_ARRAYSIZE(kPresetLabels) - 1);
              if (ImGui::Combo("Preset", &p, kPresetLabels, IM_ARRAYSIZE(kPresetLabels))) {
                a.i32Const = p;
              }
              ImGui::TextDisabled("Applies a camera rig preset (mode + tuning) via main.cpp.");
            } else if (a.kind == GameActionKind::RequestScreenshot) {
              ImGui::Checkbox("Include UI (after UI pass)", &a.bConst);

              int flags = a.i32Const;
              bool ts = (flags & ScreenshotFlag_Timestamp) != 0;
              bool cb = (flags & ScreenshotFlag_CopyToClipboard) != 0;
              bool pause = (flags & ScreenshotFlag_PauseForCapture) != 0;

              if (ImGui::Checkbox("Timestamp filename", &ts)) {
                flags = ts ? (flags | ScreenshotFlag_Timestamp) : (flags & ~ScreenshotFlag_Timestamp);
              }
              if (ImGui::Checkbox("Copy path to clipboard", &cb)) {
                flags = cb ? (flags | ScreenshotFlag_CopyToClipboard) : (flags & ~ScreenshotFlag_CopyToClipboard);
              }
              if (ImGui::Checkbox("Pause for capture", &pause)) {
                flags = pause ? (flags | ScreenshotFlag_PauseForCapture) : (flags & ~ScreenshotFlag_PauseForCapture);
              }
              a.i32Const = flags;

              ImGui::InputText("Spec template", a.msgTemplate, IM_ARRAYSIZE(a.msgTemplate));
              ImGui::TextDisabled("Format: outDir|baseName (if no '|' then baseName only).");
              ImGui::TextDisabled("Tokens: {tag} {kind} {msg} {u64a} {u64b} {tRealMs} {rule} (+ *_safe)");
            } else if (a.kind == GameActionKind::CaptureBundle) {
              int flags = a.i32Const;
              bool world = (flags & CaptureBundle_WorldScreenshot) != 0;
              bool ui = (flags & CaptureBundle_UiScreenshot) != 0;
              bool it = (flags & CaptureBundle_IntegrationTrace) != 0;
              bool ft = (flags & CaptureBundle_FlightRecorderTrace) != 0;
              bool rp = (flags & CaptureBundle_ReproPack) != 0;
              bool stop = (flags & CaptureBundle_StopFlightRecorder) != 0;
              bool pause = (flags & CaptureBundle_PauseForScreenshots) != 0;
              bool clip = (flags & CaptureBundle_CopyDirToClipboard) != 0;

              if (ImGui::Checkbox("World screenshot", &world)) {
                flags = world ? (flags | CaptureBundle_WorldScreenshot) : (flags & ~CaptureBundle_WorldScreenshot);
              }
              if (ImGui::Checkbox("UI screenshot", &ui)) {
                flags = ui ? (flags | CaptureBundle_UiScreenshot) : (flags & ~CaptureBundle_UiScreenshot);
              }
              if (ImGui::Checkbox("Integration trace (JSON)", &it)) {
                flags = it ? (flags | CaptureBundle_IntegrationTrace) : (flags & ~CaptureBundle_IntegrationTrace);
              }
              if (ImGui::Checkbox("Flight recorder trace (JSON)", &ft)) {
                flags = ft ? (flags | CaptureBundle_FlightRecorderTrace) : (flags & ~CaptureBundle_FlightRecorderTrace);
              }
              if (ImGui::Checkbox("Repro pack snapshot (JSON)", &rp)) {
                flags = rp ? (flags | CaptureBundle_ReproPack) : (flags & ~CaptureBundle_ReproPack);
              }

              ImGui::Separator();
              if (ImGui::Checkbox("Stop flight recorder before export", &stop)) {
                flags = stop ? (flags | CaptureBundle_StopFlightRecorder) : (flags & ~CaptureBundle_StopFlightRecorder);
              }
              if (ImGui::Checkbox("Pause for screenshots", &pause)) {
                flags = pause ? (flags | CaptureBundle_PauseForScreenshots) : (flags & ~CaptureBundle_PauseForScreenshots);
              }
              if (ImGui::Checkbox("Copy bundle dir to clipboard", &clip)) {
                flags = clip ? (flags | CaptureBundle_CopyDirToClipboard) : (flags & ~CaptureBundle_CopyDirToClipboard);
              }

              a.i32Const = flags;

              ImGui::InputText("Spec template", a.msgTemplate, IM_ARRAYSIZE(a.msgTemplate));
              ImGui::TextDisabled("Format: outDir|label (label is sanitized + timestamped into a folder).");
              ImGui::TextDisabled("Tokens: {tag} {kind} {msg} {u64a} {u64b} {tRealMs} {rule} (+ *_safe)");
            } else if (a.kind == GameActionKind::Toast) {
              ImGui::InputText("Message template", a.msgTemplate, IM_ARRAYSIZE(a.msgTemplate));
              ImGui::TextDisabled("Tokens: {tag} {kind} {msg} {u64a} {u64b} {tRealMs} {rule} (+ *_safe variants)");
            } else if (a.kind == GameActionKind::ExportFlightRecorderTrace || a.kind == GameActionKind::ExportIntegrationTrace) {
              ImGui::InputText("Path template", a.msgTemplate, IM_ARRAYSIZE(a.msgTemplate));
              ImGui::TextDisabled("Tokens: {tRealMs} {tag_safe} {u64a} {u64b} {rule}");
            } else if (a.kind == GameActionKind::TransmitComms) {
              // i32Const: channel (matches sim::CommsChannel values)
              static const char* kChanLabels[] = {"System", "Mission", "Security", "Pirate", "Trade", "Distress", "Custom"};
              int chan = std::clamp(a.i32Const, 0, (int)IM_ARRAYSIZE(kChanLabels) - 1);
              if (ImGui::Combo("Channel", &chan, kChanLabels, IM_ARRAYSIZE(kChanLabels))) {
                a.i32Const = chan;
              }

              int srcIdx = 0;
              for (int i = 0; i < (int)(sizeof(kU64Sources) / sizeof(kU64Sources[0])); ++i) {
                if (kU64Sources[i] == a.u64aSource) { srcIdx = i; break; }
              }
              if (ImGui::Combo("Station id source", &srcIdx, kU64SourceLabels, IM_ARRAYSIZE(kU64SourceLabels))) {
                a.u64aSource = kU64Sources[srcIdx];
              }
              if (a.u64aSource == AutomationValueSource::Constant) {
                unsigned long long v = (unsigned long long)a.u64aConst;
                if (ImGui::InputScalar("Station id", ImGuiDataType_U64, &v)) {
                  a.u64aConst = (stellar::core::u64)v;
                }
              } else {
                ImGui::TextDisabled("Using %s", automationValueSourceName(a.u64aSource));
              }

              ImGui::Checkbox("Show overlay", &a.bConst);
              ImGui::InputText("Message template", a.msgTemplate, IM_ARRAYSIZE(a.msgTemplate));
              ImGui::TextDisabled("Tokens: {tag} {kind} {msg} {u64a} {u64b} {tRealMs} {rule} (+ *_safe variants)");
              ImGui::TextDisabled("Format: subject|body (if no '|' then subject is derived from origin)");
            } else {
              ImGui::TextDisabled("(no parameters)");
            }

            if (ImGui::SmallButton("Remove action")) {
              r.actions.erase(r.actions.begin() + (long long)ai);
              ImGui::PopID();
              break;
            }

            ImGui::Separator();
            ImGui::PopID();
          }

          if (ImGui::Button("Add action")) {
            AutomationActionTemplate a;
            a.kind = GameActionKind::Toast;
            setCStr(a.msgTemplate, sizeof(a.msgTemplate), "Automation fired: {kind}:{tag}");
            r.actions.push_back(a);
          }

          ImGui::TreePop();
        }

        ImGui::PopID();
      }

      ImGui::EndTabItem();
    }

    if (ImGui::BeginTabItem("Export")) {
      ImGui::Checkbox("Pretty JSON", &st.exportPretty);
      ImGui::Checkbox("Include pending actions", &st.exportIncludePendingActions);
      ImGui::Checkbox("Include scheduled actions", &st.exportIncludeScheduledActions);
      ImGui::Checkbox("Include action history", &st.exportIncludeActionHistory);
      ImGui::Checkbox("Include events", &st.exportIncludeEvents);
      ImGui::Checkbox("Include automation rules", &st.exportIncludeAutomationRules);

      ImGui::InputText("Path", st.exportPath, IM_ARRAYSIZE(st.exportPath));

      if (ImGui::Button("Write trace JSON")) {
        std::string err;
        if (writeIntegrationTraceJson(st, st.exportPath, &err)) {
          if (toast) toast(std::string("Wrote trace: ") + st.exportPath, 2.2);
        } else {
          if (toast) toast(std::string("Trace export failed: ") + err, 3.0);
        }
      }
      ImGui::SameLine();
      if (ImGui::Button("Copy JSON to clipboard")) {
        stellar::core::JsonWriter::Options opt;
        opt.pretty = st.exportPretty;
        stellar::core::JsonWriter w(opt);

        w.beginObject();
          w.key("type"); w.value("stellar_integration_trace");
          w.key("version"); w.value(4);

          if (st.exportIncludePendingActions) {
            w.key("pendingActions");
            w.beginArray();
            for (const auto& a : st.actions.pending) jsonAction(w, a);
            w.endArray();
          }
          if (st.exportIncludeScheduledActions) {
            w.key("scheduledActions");
            w.beginArray();
            for (const auto& a : st.scheduledActions) jsonAction(w, a);
            w.endArray();
          }
          if (st.exportIncludeActionHistory) {
            w.key("actionHistory");
            w.beginArray();
            for (const auto& a : st.actions.history) jsonAction(w, a);
            w.endArray();
          }
          if (st.exportIncludeEvents) {
            w.key("events");
            w.beginArray();
            for (const auto& e : st.events.events) jsonEvent(w, e);
            w.endArray();
          }
          if (st.exportIncludeAutomationRules) {
            w.key("automationRules");
            w.beginArray();
            for (const auto& r : st.automationRules) jsonAutomationRule(w, r);
            w.endArray();
          }
          w.endObject();
        auto json = w.takeString();
        copyTextToClipboard(json);
        if (toast) toast("Copied trace JSON to clipboard.", 1.8);
      }

      ImGui::Separator();
      ImGui::TextDisabled(
          "Tip: attach integration_trace.json to bug reports.\n"
          "It captures cross-system requests (actions), resulting state changes (events),\n"
          "and the automation rules in effect."
      );

      ImGui::EndTabItem();
    }

    ImGui::EndTabBar();
  }

  ImGui::End();
}

} // namespace stellar::game
