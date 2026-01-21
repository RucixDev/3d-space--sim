#include "RuntimeValidationWindow.h"

#include "FlightRecorderWindow.h"

#include "GameSignals.h"

#include "stellar/core/Log.h"
#include "stellar/core/JsonWriter.h"

#include "stellar/math/Math.h"
#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"

#include "stellar/sim/Aerodynamics.h"
#include "stellar/sim/Atmosphere.h"
#include "stellar/sim/TimeTrial.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

namespace stellar::game {

namespace {

static bool finite(double v) {
  return std::isfinite(v);
}

static bool finiteVec(const math::Vec3d& v) {
  return finite(v.x) && finite(v.y) && finite(v.z);
}

static bool finiteQuat(const math::Quatd& q) {
  return finite(q.w) && finite(q.x) && finite(q.y) && finite(q.z);
}

static std::string fmtVec3(const math::Vec3d& v, const char* unit) {
  char buf[256];
  std::snprintf(buf, sizeof(buf), "(%.6g, %.6g, %.6g)%s", v.x, v.y, v.z, unit ? unit : "");
  return std::string(buf);
}

static void jsonVec3(core::JsonWriter& w, const math::Vec3d& v) {
  w.beginArray();
  w.value(v.x);
  w.value(v.y);
  w.value(v.z);
  w.endArray();
}

static void jsonQuat(core::JsonWriter& w, const math::Quatd& q) {
  w.beginArray();
  w.value(q.w);
  w.value(q.x);
  w.value(q.y);
  w.value(q.z);
  w.endArray();
}

static std::string uniqueDumpPath(const RuntimeValidationWindowState& st) {
  std::string path = st.dumpPath;
  if (!st.dumpUniquePerHit) return path;

  // Insert _<hit> before extension (if any).
  const std::size_t sep = path.find_last_of("/\\\\");
  std::size_t dot = path.find_last_of('.');
  if (dot == std::string::npos || (sep != std::string::npos && dot < sep)) {
    dot = path.size();
  }

  const std::string prefix = path.substr(0, dot);
  const std::string ext = path.substr(dot);
  return prefix + "_" + std::to_string(st.watchdogHits) + ext;
}

static bool writeReproPackJson(const RuntimeValidationWindowState& st,
                              const sim::Ship& ship,
                              const FlightRecorderWindowState* flightRecorder,
                              const GameEventLog* eventLog,
                              double timeRealSec,
                              double simTimeDays,
                              bool paused,
                              std::string* outPath,
                              std::string* outErr) {
  const std::string path = uniqueDumpPath(st);

  std::ofstream out(path, std::ios::binary);
  if (!out) {
    if (outErr) *outErr = "Could not open file for write: " + path;
    return false;
  }

  core::JsonWriter w(out, st.dumpPretty);

  w.beginObject();
  w.key("type"); w.value("stellar_repro_pack");
  w.key("version"); w.value(1);

  w.key("timeRealSec"); w.value(timeRealSec);
  w.key("simTimeDays"); w.value(simTimeDays);
  w.key("paused"); w.value(paused);

  w.key("watchdog");
  w.beginObject();
  w.key("hits"); w.value(st.watchdogHits);
  w.key("message"); w.value(st.watchdogLastMessage);
  w.endObject();

  w.key("build");
  w.beginObject();
  w.key("date"); w.value(__DATE__);
  w.key("time"); w.value(__TIME__);
  w.endObject();

  w.key("ship");
  w.beginObject();
  w.key("posKm"); jsonVec3(w, ship.positionKm());
  w.key("velKmS"); jsonVec3(w, ship.velocityKmS());
  w.key("angVelRadS"); jsonVec3(w, ship.angularVelocityRadS());
  w.key("orient"); jsonQuat(w, ship.orientation());
  w.endObject();

  if (st.dumpIncludeTelemetry && flightRecorder) {
    const auto& samples = flightRecorder->samples;

    w.key("flightRecorder");
    w.beginObject();
    w.key("recording"); w.value(flightRecorder->recording);
    w.key("sampleHz"); w.value(flightRecorder->sampleHz);
    w.key("sampleCount"); w.value((long long)samples.size());
    w.key("windowSec"); w.value(st.telemetryWindowSec);
    w.key("maxSamples"); w.value(st.telemetryMaxSamples);

    w.key("samples");
    w.beginArray();

    if (!samples.empty()) {
      const double lastT = samples.back().tRealSec;
      const double minT = lastT - std::max(0.0, st.telemetryWindowSec);

      // Collect from the back so we can stop early.
      std::vector<const FlightRecorderSample*> slice;
      slice.reserve((std::size_t)std::max(0, st.telemetryMaxSamples));

      for (auto it = samples.rbegin(); it != samples.rend(); ++it) {
        if ((int)slice.size() >= st.telemetryMaxSamples) break;
        if (it->tRealSec < minT) break;
        slice.push_back(&*it);
      }

      // Emit in chronological order.
      for (auto it = slice.rbegin(); it != slice.rend(); ++it) {
        const FlightRecorderSample& s = **it;
        w.beginObject();
        w.key("tRealSec"); w.value(s.tRealSec);
        w.key("tSimDays"); w.value(s.tSimDays);
        w.key("posKm"); jsonVec3(w, s.posKm);
        w.key("velKmS"); jsonVec3(w, s.velKmS);
        w.key("orient"); jsonQuat(w, s.orient);
        w.key("angVelRadS"); jsonVec3(w, s.angVelRadS);
        w.endObject();
      }
    }

    w.endArray();
    w.endObject();
  }

  if (st.dumpIncludeEvents && eventLog) {
    const auto& evs = eventLog->events;

    w.key("eventLog");
    w.beginObject();
    w.key("count"); w.value((long long)evs.size());
    w.key("windowSec"); w.value(st.dumpEventWindowSec);
    w.key("maxEvents"); w.value(st.dumpMaxEvents);

    w.key("events");
    w.beginArray();

    if (!evs.empty()) {
      const double minT = timeRealSec - std::max(0.0, st.dumpEventWindowSec);

      std::vector<const GameEvent*> slice;
      slice.reserve((std::size_t)std::max(0, st.dumpMaxEvents));

      for (auto it = evs.rbegin(); it != evs.rend(); ++it) {
        if ((int)slice.size() >= st.dumpMaxEvents) break;
        if (it->tRealSec < minT) break;
        slice.push_back(&*it);
      }

      for (auto it = slice.rbegin(); it != slice.rend(); ++it) {
        const GameEvent& e = **it;
        w.beginObject();
        w.key("tRealSec"); w.value(e.tRealSec);
        w.key("tSimDays"); w.value(e.tSimDays);
        w.key("kind"); w.value(gameEventKindName(e.kind));
        w.key("tag"); w.value(e.tag);
        w.key("msg"); w.value(e.msg);
        if (e.hasPos) {
          w.key("posKm"); jsonVec3(w, e.posKm);
        }
        if (e.u64a != 0) {
          w.key("u64a"); w.value((long long)e.u64a);
        }
        if (e.u64b != 0) {
          w.key("u64b"); w.value((long long)e.u64b);
        }
        w.endObject();
      }
    }

    w.endArray();
    w.endObject();
  }

  w.endObject();

  if (!out) {
    if (outErr) *outErr = "Write failed: " + path;
    return false;
  }

  if (outPath) *outPath = path;
  return true;
}
static void addCheck(std::vector<ValidationCheck>& out, const std::string& name, bool ok, const std::string& details = {}) {
  out.push_back(ValidationCheck{name, ok, details});
}

static std::string buildReport(const std::vector<ValidationCheck>& checks, int failCount) {
  std::ostringstream ss;
  ss << "Runtime Validation Report\n";
  ss << "--------------------------\n";
  ss << "Checks: " << checks.size() << ", Failures: " << failCount << "\n\n";
  for (const auto& c : checks) {
    ss << (c.ok ? "[PASS] " : "[FAIL] ") << c.name;
    if (!c.details.empty()) ss << " — " << c.details;
    ss << "\n";
  }
  return ss.str();
}

static void runSmokeChecks(RuntimeValidationWindowState& st) {
  std::vector<ValidationCheck> checks;

  // --- TimeTrial gate crossing ---
  {
    sim::TimeTrialGate g{};
    g.posKm = {0, 0, 0};
    g.normal = {0, 0, 1};
    g.radiusKm = 1.0;

    const bool passForward = sim::timeTrialGatePassed(g, {0, 0, -2}, {0, 0, +2}, {0, 0, +1});
    const bool rejectBack = !sim::timeTrialGatePassed(g, {0, 0, +2}, {0, 0, -2}, {0, 0, -1});
    const bool rejectOutside = !sim::timeTrialGatePassed(g, {2.5, 0, -2}, {2.5, 0, +2}, {0, 0, +1});

    addCheck(checks, "TimeTrial: gate crossing direction", passForward && rejectBack,
             (passForward && rejectBack) ? "ok" : "direction checks failed");
    addCheck(checks, "TimeTrial: gate radius reject", rejectOutside,
             rejectOutside ? "ok" : "outside-radius crossing was accepted");
  }

  // --- TimeTrial generator invariants ---
  {
    sim::TimeTrialCourseParams p{};
    p.gateCount = 16;
    p.gateRadiusKm = 2500.0;
    p.courseRadiusKm = 60000.0;
    p.courseHeightKm = 12000.0;
    p.jitterKm = 8000.0;
    p.loops = 2;
    p.closedLoop = false;

    const auto c0 = sim::generateTimeTrialCourseStationSlalomKm({0, 0, 0}, math::Quatd::identity(), 12345u, p);
    const auto c1 = sim::generateTimeTrialCourseStationSlalomKm({0, 0, 0}, math::Quatd::identity(), 12345u, p);

    bool ok = true;
    ok = ok && (c0.gates.size() == (std::size_t)p.gateCount);
    ok = ok && (c1.gates.size() == (std::size_t)p.gateCount);
    ok = ok && (c0.key == c1.key);

    for (std::size_t i = 0; ok && i < c0.gates.size(); ++i) {
      const auto& a = c0.gates[i];
      const auto& b = c1.gates[i];
      ok = ok && ((a.posKm - b.posKm).length() < 1e-6);
      ok = ok && (std::abs(a.radiusKm - p.gateRadiusKm) < 1e-12);
      ok = ok && (std::abs(a.normal.length() - 1.0) < 1e-6);
    }

    addCheck(checks, "TimeTrial: generator determinism", ok, ok ? "stable" : "course mismatch / invariant failure");
  }

  // --- Aerodynamics basic invariants ---
  {
    sim::AtmosphereSample atmo{};
    atmo.inAtmosphere = true;
    atmo.dynamicPressurePa = 5000.0;
    atmo.relVelKmS = {0, 0, 10};

    sim::AerodynamicsParams ap{};
    ap.enabled = true;
    ap.minDynamicPressurePa = 25.0;
    ap.wingAreaM2 = 120.0;
    ap.liftSlopePerRad = 5.0;
    ap.clMax = 1.25;
    ap.stallAngleDeg = 18.0;
    ap.controlSurfaces = false;

    const auto s0 = sim::computeAerodynamics(atmo, math::Quatd::identity(), {0, 0, 0}, 10000.0, ap);
    bool ok0 = s0.active;
    ok0 = ok0 && (std::abs(s0.alphaRad) < 1e-6);
    ok0 = ok0 && (std::abs(s0.cl) < 1e-6);
    ok0 = ok0 && (s0.liftAccelKmS2.length() < 1e-9);
    ok0 = ok0 && finiteVec(s0.extraDragAccelKmS2);
    addCheck(checks, "Aerodynamics: neutral AoA yields ~zero lift", ok0, ok0 ? "ok" : "unexpected lift at alpha=0");

    // Small positive AoA (downward velocity component produces wind from below).
    atmo.relVelKmS = {0, -1, 10};
    const auto s1 = sim::computeAerodynamics(atmo, math::Quatd::identity(), {0, 0, 0}, 10000.0, ap);
    bool ok1 = s1.active;
    ok1 = ok1 && (s1.alphaRad > 0.0);
    ok1 = ok1 && (s1.cl > 0.0);
    ok1 = ok1 && (s1.liftAccelKmS2.y > 0.0);
    ok1 = ok1 && finiteVec(s1.liftAccelKmS2);
    addCheck(checks, "Aerodynamics: positive AoA produces positive lift", ok1, ok1 ? "ok" : "lift sign/invariants failed");
  }

  // Summarize
  int fails = 0;
  for (const auto& c : checks) {
    if (!c.ok) ++fails;
  }

  st.lastChecks = std::move(checks);
  st.lastFailCount = fails;
  st.lastReport = buildReport(st.lastChecks, fails);
}

} // namespace

bool exportReproPackJsonToPath(const RuntimeValidationWindowState& st,
                              const sim::Ship& ship,
                              const FlightRecorderWindowState* flightRecorder,
                              const GameEventLog* eventLogForDump,
                              double timeRealSec,
                              double simTimeDays,
                              bool paused,
                              const char* path,
                              std::string* outErr) {
  if (!path || !path[0]) {
    if (outErr) *outErr = "Empty path";
    return false;
  }

  // Write using the existing repro pack machinery, but keep the UI state immutable.
  RuntimeValidationWindowState tmp = st;
  tmp.dumpUniquePerHit = false;
  std::snprintf(tmp.dumpPath, sizeof(tmp.dumpPath), "%s", path);

  std::string outPath;
  std::string err;
  const bool ok = writeReproPackJson(tmp, ship, flightRecorder, eventLogForDump, timeRealSec, simTimeDays, paused, &outPath, &err);
  if (!ok) {
    if (outErr) *outErr = err;
    return false;
  }
  return true;
}

void tickRuntimeValidation(RuntimeValidationWindowState& st,
                           const sim::Ship& ship,
                           const FlightRecorderWindowState* flightRecorder,
                           double timeRealSec,
                           double simTimeDays,
                           bool& paused,
                           const ToastFn& toast,
                           const EventSinkFn& emitEvent,
                           const ActionSinkFn& emitAction,
                           const GameEventLog* eventLogForDump) {
  if (!st.watchdogEnabled) {
    st.watchdogLatched = false;
    return;
  }

  const math::Vec3d p = ship.positionKm();
  const math::Vec3d v = ship.velocityKmS();
  const math::Vec3d w = ship.angularVelocityRadS();
  const math::Quatd q = ship.orientation();

  bool ok = true;
  ok = ok && finiteVec(p) && finiteVec(v) && finiteVec(w) && finiteQuat(q);

  if (ok && st.magnitudeWatchdog) {
    const auto absMax = [](const math::Vec3d& a) {
      return std::max({std::abs(a.x), std::abs(a.y), std::abs(a.z)});
    };
    ok = ok && (absMax(p) <= st.maxAbsPosKm);
    ok = ok && (absMax(v) <= st.maxAbsVelKmS);
    ok = ok && (absMax(w) <= st.maxAbsAngVelRadS);
  }

  if (ok) {
    st.watchdogLatched = false;
    st.dumpLatchedThisFailure = false;
    st.screenshotLatchedThisFailure = false;
    return;
  }

  // Latch so we don't spam logs/toasts every frame.
  if (st.watchdogLatched && st.suppressDuplicateReports) {
    return;
  }

  st.watchdogLatched = true;
  st.watchdogHits++;

  std::string msg = "Ship state invalid: ";
  if (!finiteVec(p)) msg += "pos=" + fmtVec3(p, " km") + " ";
  if (!finiteVec(v)) msg += "vel=" + fmtVec3(v, " km/s") + " ";
  if (!finiteVec(w)) msg += "angVel=" + fmtVec3(w, " rad/s") + " ";
  if (!finiteQuat(q)) msg += "orient=(non-finite quaternion) ";

  if (st.magnitudeWatchdog) {
    const auto absMax = [](const math::Vec3d& a) {
      return std::max({std::abs(a.x), std::abs(a.y), std::abs(a.z)});
    };
    if (absMax(p) > st.maxAbsPosKm) msg += "|pos| too large ";
    if (absMax(v) > st.maxAbsVelKmS) msg += "|vel| too large ";
    if (absMax(w) > st.maxAbsAngVelRadS) msg += "|angVel| too large ";
  }

  st.watchdogLastMessage = msg;

  {
    GameEvent ev{};
    ev.tRealSec = timeRealSec;
    ev.tSimDays = simTimeDays;
    ev.kind = GameEventKind::Validation;
    ev.tag = "Watchdog";
    ev.msg = msg;
    ev.hasPos = true;
    ev.posKm = p;
    ev.u64a = (core::u64)st.watchdogHits;
    if (emitEvent) emitEvent(std::move(ev));
  }

  // Optional: request screenshots via the Integration Hub action queue.
  if (st.screenshotOnFailure && !st.screenshotLatchedThisFailure) {
    int flags = 0;
    if (st.screenshotTimestamp) flags |= ScreenshotFlag_Timestamp;
    if (st.screenshotCopyToClipboard) flags |= ScreenshotFlag_CopyToClipboard;
    if (st.screenshotPauseForCapture) flags |= ScreenshotFlag_PauseForCapture;

    const std::string outDir = st.screenshotOutDir[0] ? st.screenshotOutDir : "screenshots";
    const std::string baseName = st.screenshotBaseName[0] ? st.screenshotBaseName : "validation";

    // Capture UI pass (fast, same-frame if invoked early enough).
    if (emitAction && st.screenshotIncludeUi) {
      emitAction(makeActionRequestScreenshot(timeRealSec, simTimeDays, "RuntimeValidation", /*includeUi=*/true, flags, outDir, baseName));
    }
    // Capture world pass (may land next frame depending on scheduling point).
    if (emitAction && st.screenshotIncludeWorld) {
      emitAction(makeActionRequestScreenshot(timeRealSec, simTimeDays, "RuntimeValidation", /*includeUi=*/false, flags, outDir, baseName));
    }

    st.screenshotLatchedThisFailure = true;
  }

  // Optional: dump a tiny repro pack JSON for bug reports / regression tracking.
  if (st.dumpReproOnFailure && !st.dumpLatchedThisFailure) {
    std::string outPath;
    std::string err;
    if (writeReproPackJson(st, ship, flightRecorder, eventLogForDump, timeRealSec, simTimeDays, paused, &outPath, &err)) {
      st.dumpsWritten++;
      st.lastDumpPath = outPath;
      st.lastDumpError.clear();
      if (toast) toast("Wrote repro pack: " + outPath, 2.0);
      if (emitEvent) {
        emitEvent(GameEvent{timeRealSec, simTimeDays, GameEventKind::Validation, "ReproPack", "Wrote repro pack: " + outPath});
      }
    } else {
      st.lastDumpError = err;
      if (toast) toast("Repro pack write failed.", 2.2);
      if (emitEvent) {
        emitEvent(GameEvent{timeRealSec, simTimeDays, GameEventKind::Validation, "ReproPack", "Repro pack write failed: " + err});
      }
    }
    st.dumpLatchedThisFailure = true;
  }

  if (st.logOnFailure) {
    core::log(core::LogLevel::Error, "[Validation] " + msg);
  }
  if (toast) {
    toast("Validation failure detected. See Runtime Validation window.", 2.6);
  }
  if (st.pauseOnFailure) {
    paused = true;
  }
}

void drawRuntimeValidationWindow(RuntimeValidationWindowState& st,
                                 const sim::Ship& ship,
                                 const FlightRecorderWindowState* flightRecorder,
                                 double timeRealSec,
                                 double simTimeDays,
                                 bool paused,
                                 const ToastFn& toast,
                                 const EventSinkFn& emitEvent,
                                 const ActionSinkFn& emitAction,
                                 const GameEventLog* eventLogForDump) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(560, 520), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Runtime Validation", &st.open)) {
    ImGui::End();
    return;
  }

  ImGui::TextDisabled("Catch simulation corruption early and run deterministic smoke checks.");
  ImGui::Separator();

  // ---- Watchdog ----
  ImGui::Text("Watchdog");
  ImGui::Indent();
  ImGui::Checkbox("Enable watchdog", &st.watchdogEnabled);
  ImGui::Checkbox("Pause on failure", &st.pauseOnFailure);
  ImGui::Checkbox("Log on failure", &st.logOnFailure);
  ImGui::Checkbox("Suppress duplicate reports", &st.suppressDuplicateReports);
  ImGui::Checkbox("Magnitude watchdog", &st.magnitudeWatchdog);
  // Screenshot hook
  ImGui::Separator();
  ImGui::Text("On watchdog failure");
  ImGui::Indent();
  ImGui::Checkbox("Request screenshot(s)", &st.screenshotOnFailure);
  if (st.screenshotOnFailure) {
    ImGui::Checkbox("Include UI", &st.screenshotIncludeUi);
    ImGui::Checkbox("Include world", &st.screenshotIncludeWorld);
    ImGui::Checkbox("Pause for capture", &st.screenshotPauseForCapture);
    ImGui::Checkbox("Timestamp filename", &st.screenshotTimestamp);
    ImGui::Checkbox("Copy path to clipboard", &st.screenshotCopyToClipboard);
    ImGui::InputText("Out dir", st.screenshotOutDir, IM_ARRAYSIZE(st.screenshotOutDir));
    ImGui::InputText("Base name", st.screenshotBaseName, IM_ARRAYSIZE(st.screenshotBaseName));
  }
  ImGui::Unindent();

  if (st.magnitudeWatchdog) {
    ImGui::SetNextItemWidth(180);
    ImGui::InputDouble("Max |pos| (km)", &st.maxAbsPosKm, 0, 0, "%.3g");
    ImGui::SetNextItemWidth(180);
    ImGui::InputDouble("Max |vel| (km/s)", &st.maxAbsVelKmS, 0, 0, "%.3g");
    ImGui::SetNextItemWidth(180);
    ImGui::InputDouble("Max |angVel| (rad/s)", &st.maxAbsAngVelRadS, 0, 0, "%.3g");
  }

  ImGui::TextDisabled("Hits: %d", st.watchdogHits);
  if (!st.watchdogLastMessage.empty()) {
    const ImVec4 col = st.watchdogLatched ? ImVec4(1.0f, 0.35f, 0.35f, 1.0f) : ImVec4(0.7f, 0.7f, 0.7f, 1.0f);
    ImGui::TextColored(col, "%s", st.watchdogLastMessage.c_str());
  }
  ImGui::Unindent();

  ImGui::Separator();

  // ---- Repro pack ----
  ImGui::Text("Repro pack (JSON)");
  ImGui::Indent();
  ImGui::Checkbox("Dump on watchdog failure", &st.dumpReproOnFailure);
  ImGui::SameLine();
  ImGui::Checkbox("Unique name per hit", &st.dumpUniquePerHit);
  ImGui::Checkbox("Pretty JSON", &st.dumpPretty);
  ImGui::Checkbox("Include Flight Recorder telemetry", &st.dumpIncludeTelemetry);
  if (st.dumpIncludeTelemetry) {
    ImGui::SetNextItemWidth(180);
    ImGui::InputDouble("Telemetry window (sec)", &st.telemetryWindowSec, 0, 0, "%.3g");
    ImGui::SetNextItemWidth(180);
    ImGui::InputInt("Telemetry max samples", &st.telemetryMaxSamples);
    if (st.telemetryMaxSamples < 0) st.telemetryMaxSamples = 0;
  }

  ImGui::Checkbox("Include Integration Hub events", &st.dumpIncludeEvents);
  if (st.dumpIncludeEvents) {
    ImGui::SetNextItemWidth(180);
    ImGui::InputDouble("Event window (sec)", &st.dumpEventWindowSec, 0, 0, "%.3g");
    ImGui::SetNextItemWidth(180);
    ImGui::InputInt("Event max entries", &st.dumpMaxEvents);
    if (st.dumpMaxEvents < 0) st.dumpMaxEvents = 0;
    if (st.dumpEventWindowSec < 0.0) st.dumpEventWindowSec = 0.0;
  }

  ImGui::SetNextItemWidth(420);
  ImGui::InputText("Path", st.dumpPath, sizeof(st.dumpPath));

  if (ImGui::Button("Dump now")) {
    std::string outPath;
    std::string err;
    if (writeReproPackJson(st, ship, flightRecorder, eventLogForDump, timeRealSec, simTimeDays, paused, &outPath, &err)) {
      st.dumpsWritten++;
      st.lastDumpPath = outPath;
      st.lastDumpError.clear();
      if (toast) toast("Wrote repro pack: " + outPath, 2.0);
      if (emitEvent) {
        emitEvent(GameEvent{timeRealSec, simTimeDays, GameEventKind::Validation, "ReproPack", "Wrote repro pack: " + outPath});
      }
    } else {
      st.lastDumpError = err;
      if (toast) toast("Repro pack write failed.", 2.2);
      if (emitEvent) {
        emitEvent(GameEvent{timeRealSec, simTimeDays, GameEventKind::Validation, "ReproPack", "Repro pack write failed: " + err});
      }
    }
  }

  ImGui::SameLine();
  if (emitAction) {
    if (ImGui::Button("Capture bundle")) {
      const int flags = CaptureBundle_Default | CaptureBundle_PauseForScreenshots | CaptureBundle_StopFlightRecorder;
      emitAction(makeActionCaptureBundle(timeRealSec, simTimeDays, "RuntimeValidation", flags, "captures", "validation"));
      if (toast) toast("Capture bundle requested.", 1.6);
    }
  } else {
    ImGui::BeginDisabled();
    ImGui::Button("Capture bundle");
    ImGui::EndDisabled();
  }

  ImGui::TextDisabled("Dumps written: %d", st.dumpsWritten);
  if (!st.lastDumpPath.empty()) {
    ImGui::TextWrapped("Last dump: %s", st.lastDumpPath.c_str());
  }
  if (!st.lastDumpError.empty()) {
    ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f), "Last dump error: %s", st.lastDumpError.c_str());
  }
  ImGui::Unindent();

  ImGui::Separator();

  // ---- Smoke tests ----
  ImGui::Text("Deterministic smoke checks");
  ImGui::Indent();
  if (ImGui::Button("Run checks")) {
    runSmokeChecks(st);
    if (toast) {
      if (st.lastFailCount == 0) toast("Validation checks: all PASS", 1.6);
      else toast("Validation checks: failures found", 2.2);
    }
  }
  ImGui::SameLine();
  if (ImGui::Button("Copy report")) {
    if (st.lastReport.empty()) runSmokeChecks(st);
    ImGui::SetClipboardText(st.lastReport.c_str());
    if (toast) toast("Validation report copied.", 1.4);
  }
  ImGui::SameLine();
  if (ImGui::Button("Clear")) {
    st.lastChecks.clear();
    st.lastReport.clear();
    st.lastFailCount = 0;
  }

  if (!st.lastChecks.empty()) {
    const ImVec4 passCol = ImVec4(0.35f, 0.95f, 0.55f, 1.0f);
    const ImVec4 failCol = ImVec4(1.0f, 0.35f, 0.35f, 1.0f);
    ImGui::Text("Failures: %d", st.lastFailCount);
    ImGui::Separator();
    for (const auto& c : st.lastChecks) {
      ImGui::PushID(&c);
      ImGui::TextColored(c.ok ? passCol : failCol, "%s", c.ok ? "PASS" : "FAIL");
      ImGui::SameLine();
      ImGui::TextUnformatted(c.name.c_str());
      if (!c.details.empty()) {
        ImGui::Indent();
        ImGui::TextWrapped("%s", c.details.c_str());
        ImGui::Unindent();
      }
      ImGui::PopID();
    }
  } else {
    ImGui::TextDisabled("No checks run yet.");
  }
  ImGui::Unindent();

  ImGui::End();
}

} // namespace stellar::game
