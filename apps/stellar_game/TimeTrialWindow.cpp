#include "TimeTrialWindow.h"

#include "stellar/core/Hash.h"
#include "stellar/sim/System.h"
#include "stellar/sim/Units.h"
#include "stellar/sim/Ship.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

namespace stellar::game {

namespace {

static std::string formatTime(double sec) {
  if (!(sec > 0.0)) return "--:--.---";
  const int minutes = (int)std::floor(sec / 60.0);
  const double s = sec - (double)minutes * 60.0;
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%02d:%06.3f", minutes, s);
  return std::string(buf);
}

static core::u64 computeCourseSeed(const sim::StarSystem& sys,
                                  const sim::Station& st,
                                  core::u64 userSeed,
                                  bool includeSystemSeed) {
  core::u64 h = core::fnv1a64("TimeTrialSeed");
  if (includeSystemSeed) {
    h = core::hashCombine(h, sys.stub.seed);
    h = core::hashCombine(h, static_cast<core::u64>(sys.stub.id));
  }
  h = core::hashCombine(h, static_cast<core::u64>(st.id));
  h = core::hashCombine(h, userSeed);
  return h;
}

static std::string courseCode(const sim::StarSystem& sys,
                              const sim::Station& st,
                              core::u64 seed,
                              const sim::TimeTrialCourseParams& p) {
  // Shareable string that uniquely captures the generator knobs.
  // (Not a formal serialization; just a nice "human" fingerprint.)
  char buf[128];
  std::snprintf(buf, sizeof(buf),
                "TT-%llu-%llu-%llu-G%d-R%.0f-CR%.0f-H%.0f",
                (unsigned long long)sys.stub.id,
                (unsigned long long)st.id,
                (unsigned long long)seed,
                p.gateCount,
                p.gateRadiusKm,
                p.courseRadiusKm,
                p.courseHeightKm);
  return std::string(buf);
}

// --- Ghost replay helpers ---
// We reuse the FlightRecorderSample struct for a tiny, interpolation-friendly
// representation of ship motion.


static bool sampleGhostAtTime(const std::deque<FlightRecorderSample>& samples,
                             double tSec,
                             FlightRecorderSample* out) {
  if (!out) return false;
  if (samples.empty()) return false;

  if (samples.size() == 1) {
    *out = samples.front();
    return true;
  }

  tSec = std::max(0.0, tSec);

  if (tSec <= samples.front().tRealSec) {
    *out = samples.front();
    return true;
  }
  if (tSec >= samples.back().tRealSec) {
    *out = samples.back();
    return true;
  }

  // Binary search for bracketing samples.
  std::size_t lo = 0;
  std::size_t hi = samples.size() - 1;
  while (hi - lo > 1) {
    const std::size_t mid = (lo + hi) / 2;
    if (samples[mid].tRealSec <= tSec) lo = mid;
    else hi = mid;
  }

  const auto& a = samples[lo];
  const auto& b = samples[hi];
  const double denom = std::max(1e-9, b.tRealSec - a.tRealSec);
  const double u = std::clamp((tSec - a.tRealSec) / denom, 0.0, 1.0);

  FlightRecorderSample s{};
  s.tRealSec = tSec;
  s.tSimDays = a.tSimDays * (1.0 - u) + b.tSimDays * u;
  s.posKm = a.posKm * (1.0 - u) + b.posKm * u;
  s.velKmS = a.velKmS * (1.0 - u) + b.velKmS * u;
  s.orient = math::Quatd::nlerp(a.orient, b.orient, u);
  s.angVelRadS = a.angVelRadS * (1.0 - u) + b.angVelRadS * u;

  *out = s;
  return true;
}

static void pushGhostSample(TimeTrialWindowState& st,
                            double tRunSec,
                            double simTimeDays,
                            const sim::Ship& ship,
                            const math::Vec3d* overridePosKm = nullptr) {
  FlightRecorderSample s;
  s.tRealSec = std::max(0.0, tRunSec);
  s.tSimDays = simTimeDays;
  s.posKm = overridePosKm ? *overridePosKm : ship.positionKm();
  s.velKmS = ship.velocityKmS();
  s.orient = ship.orientation();
  s.angVelRadS = ship.angularVelocityRadS();

  // Enforce ring buffer.
  if (st.ghostMaxSamples <= 0) st.ghostMaxSamples = 1;
  while ((int)st.ghostRunSamples.size() >= st.ghostMaxSamples) {
    st.ghostRunSamples.pop_front();
  }
  st.ghostRunSamples.push_back(s);
}

} // namespace

void TimeTrialWindowState::clearCourse() {
  phase = TimeTrialPhase::Inactive;
  hasCourse = false;
  course = sim::TimeTrialCourse{};
  anchorStationId = 0;
  anchorStationIndex = -1;
  nextGate = 0;
  timeSec = 0.0;
  finishTimeSec = 0.0;
  bestTimeSec = 0.0;
  hasPrevShipPos = false;
  prevShipPosKm = {};
  requestTargetAnchor = false;
  requestEngageDockingComputer = false;
  requestCameraChase = false;
  requestCameraOrbit = false;
  lastGatePassed = -1;
  gatePulseSec = 0.0;

  ghostSampleValid = false;
  ghostSampleAccumulatorSec = 0.0;
  ghostRunSamples.clear();

  justArmed = false;
}

void TimeTrialWindowState::armCourse() {
  if (!hasCourse) {
    phase = TimeTrialPhase::Inactive;
    return;
  }
  phase = TimeTrialPhase::Ready;
  nextGate = 0;
  timeSec = 0.0;
  finishTimeSec = 0.0;
  requestTargetAnchor = false;
  requestEngageDockingComputer = false;
  requestCameraChase = false;
  requestCameraOrbit = false;
  hasPrevShipPos = false;
  lastGatePassed = -1;
  gatePulseSec = 0.0;

  ghostSampleValid = false;
  ghostSampleAccumulatorSec = 0.0;
  ghostRunSamples.clear();

  justArmed = true;
}

void TimeTrialWindowState::cancelRun() {
  // Cancel the active run, but keep the course loaded so the player can re-arm
  // quickly from the UI/HUD.
  if (!hasCourse) {
    clearCourse();
    return;
  }
  phase = TimeTrialPhase::Inactive;
  nextGate = 0;
  timeSec = 0.0;
  finishTimeSec = 0.0;
  hasPrevShipPos = false;
  lastGatePassed = -1;
  gatePulseSec = 0.0;
  requestTargetAnchor = false;
  requestEngageDockingComputer = false;
  requestCameraChase = false;
  requestCameraOrbit = false;

  ghostSampleValid = false;
  ghostSampleAccumulatorSec = 0.0;
  ghostRunSamples.clear();

  justArmed = false;
}

void TimeTrialWindowState::armFromCourse(const sim::TimeTrialCourse& c) {
  course = c;
  hasCourse = !course.gates.empty();
  // Refresh cached best time.
  bestTimeSec = (bestTimesSec.count(course.key) > 0) ? bestTimesSec[course.key] : 0.0;
  armCourse();
}


void tickTimeTrial(TimeTrialWindowState& st,
                   double dtRealSec,
                   double timeRealSec,
                   double simTimeDays,
                   const sim::Ship& ship,
                   bool paused,
                   bool docked,
                   sim::StationId dockedStationId,
                   const ToastFn& toast,
                   const ActionSinkFn& emitAction,
                   const EventSinkFn& emitEvent) {
  if (!st.hasCourse || st.course.gates.empty()) {
    st.hasPrevShipPos = false;
    st.ghostSampleValid = false;
    return;
  }

  const math::Vec3d posKm = ship.positionKm();

  auto emitTTEvent = [&](const std::string& tag, const std::string& msg, core::u64 u64b = 0) {
    if (!emitEvent) return;
    GameEvent e;
    e.tRealSec = timeRealSec;
    e.tSimDays = simTimeDays;
    e.kind = GameEventKind::TimeTrial;
    e.tag = tag;
    e.msg = msg;
    e.hasPos = true;
    e.posKm = posKm;
    e.u64a = st.course.key;
    e.u64b = u64b;
    emitEvent(e);
  };

  auto updateGhostSample = [&](double tNowSec) {
    st.ghostSampleValid = false;
    if (!st.ghostEnabled) return;
    auto it = st.bestGhostSamples.find(st.course.key);
    if (it == st.bestGhostSamples.end() || it->second.empty()) return;
    const double tGhost = tNowSec + st.ghostLeadSec;
    st.ghostSampleValid = sampleGhostAtTime(it->second, tGhost, &st.ghostSample);
  };

  // Pulse decay.
  if (st.gatePulseSec > 0.0) {
    st.gatePulseSec = std::max(0.0, st.gatePulseSec - dtRealSec);
  }

  // When (re)arming a course, optionally clear the Integration Hub log so
  // exported traces contain a clean run.
  if (st.justArmed) {
    st.justArmed = false;
    if (st.autoClearIntegrationHubOnArm && emitAction) {
      emitAction(makeActionClearIntegrationHub(timeRealSec, simTimeDays, "TimeTrial"));
    }
  }

  if (paused) {
    st.prevShipPosKm = posKm;
    st.hasPrevShipPos = true;
    const double tNow = (st.phase == TimeTrialPhase::Finished) ? st.finishTimeSec : st.timeSec;
    updateGhostSample(tNow);
    return;
  }

  // One-shot integration requests are edge-triggered and consumed by main.cpp.
  st.requestTargetAnchor = false;
  st.requestEngageDockingComputer = false;
  st.requestCameraChase = false;
  st.requestCameraOrbit = false;

  bool startedThisFrame = false;
  double startFrac = 0.0;

  // Gate progression (only while in the gate phases).
  if (st.hasPrevShipPos && (st.phase == TimeTrialPhase::Ready || st.phase == TimeTrialPhase::Running)) {
    const int idx = std::clamp(st.nextGate, 0, (int)st.course.gates.size() - 1);
    const auto& g = st.course.gates[idx];

    double tHit = 0.0;
    if (sim::timeTrialGatePassed(g, st.prevShipPosKm, posKm, ship.velocityKmS(), &tHit)) {
      st.lastGatePassed = idx;
      st.gatePulseSec = 0.35;

      if (st.phase == TimeTrialPhase::Ready) {
        // Crossing the first gate starts the timer.
        st.phase = TimeTrialPhase::Running;
        st.timeSec = 0.0;
        startedThisFrame = true;
        startFrac = std::clamp(tHit, 0.0, 1.0);

        // Ghost recording starts on the start gate crossing.
        st.ghostRunSamples.clear();
        st.ghostSampleAccumulatorSec = 0.0;
        if (st.ghostRecordRun) {
          const math::Vec3d startPosKm = st.prevShipPosKm * (1.0 - startFrac) + posKm * startFrac;
          pushGhostSample(st, 0.0, simTimeDays, ship, &startPosKm);
        }

        if (st.autoCameraChaseOnStart) {
          st.requestCameraChase = true;
          if (emitAction) {
            // CameraRigMode::Chase == 0
            emitAction(makeActionSetCameraRigPreset(timeRealSec, simTimeDays, "TimeTrial", 1)); // Travel preset
          }
        }

        // Run capture: start telemetry recording on the start gate.
        if (st.autoStartFlightRecorderOnStart && emitAction) {
          emitAction(makeActionStartFlightRecorder(timeRealSec, simTimeDays, "TimeTrial"));
        }

        emitTTEvent("Start", "Time Trial started", (core::u64)idx);
        if (toast) toast("Time Trial started!", 1.6);
      } else {
        emitTTEvent("Gate", "Gate passed", (core::u64)idx);
        if (toast) toast("Gate passed", 0.8);
      }

      st.nextGate = idx + 1;

      // Last gate?
      if (st.nextGate >= (int)st.course.gates.size()) {
        if (st.finishMode == TimeTrialFinishMode::GatesOnly) {
          // Finish at the instant we crossed the final gate (sub-frame accurate).
          const double finishAt = st.timeSec + dtRealSec * std::clamp(tHit, 0.0, 1.0);
          st.phase = TimeTrialPhase::Finished;
          st.finishTimeSec = finishAt;
          st.timeSec = finishAt;

          // Ghost recording: ensure we capture a terminal sample for the run.
          if (st.ghostRecordRun) {
            pushGhostSample(st, st.finishTimeSec, simTimeDays, ship);
          }

          const double prevBest = (st.bestTimesSec.count(st.course.key) > 0)
                                    ? st.bestTimesSec[st.course.key]
                                    : 0.0;
          const bool isBest = !(prevBest > 0.0) || (st.finishTimeSec < prevBest);
          if (isBest) {
            st.bestTimesSec[st.course.key] = st.finishTimeSec;
            st.bestTimeSec = st.finishTimeSec;

            if (st.ghostSaveBestRun && st.ghostRecordRun && !st.ghostRunSamples.empty()) {
              st.bestGhostSamples[st.course.key] = st.ghostRunSamples;
            }

            if (toast) toast("Finished! New best time: " + formatTime(st.finishTimeSec), 2.5);
            emitTTEvent("Finish", "Finished (new best)", (core::u64)(st.finishTimeSec * 1000.0));
          } else {
            st.bestTimeSec = prevBest;
            if (toast) toast("Finished: " + formatTime(st.finishTimeSec) + " (best " + formatTime(prevBest) + ")", 2.4);
            emitTTEvent("Finish", "Finished", (core::u64)(st.finishTimeSec * 1000.0));
          }

          // Run capture: stop recording and optionally export traces.
          if (emitAction) {
            if (st.autoStopFlightRecorderOnFinish) {
              emitAction(makeActionStopFlightRecorder(timeRealSec, simTimeDays, "TimeTrial"));
            }
            if (st.autoExportTracesOnFinish) {
              const core::u64 stampMs = (core::u64)std::llround(timeRealSec * 1000.0);
              char flightPath[256];
              char integPath[256];
              std::snprintf(flightPath, sizeof(flightPath), "time_trial_%llu_%llu_flight_trace.json",
                            (unsigned long long)st.course.key, (unsigned long long)stampMs);
              std::snprintf(integPath, sizeof(integPath), "time_trial_%llu_%llu_integration_trace.json",
                            (unsigned long long)st.course.key, (unsigned long long)stampMs);
              emitAction(makeActionExportFlightRecorderTrace(timeRealSec, simTimeDays, "TimeTrial", flightPath));
              emitAction(makeActionExportIntegrationTrace(timeRealSec, simTimeDays, "TimeTrial", integPath));
            }
          }
        } else {
          // Cross-system finish: after the last gate, require docking at the anchor station.
          st.phase = TimeTrialPhase::Docking;

          if (st.autoTargetAnchorOnDocking) {
            st.requestTargetAnchor = true;
            if (emitAction) emitAction(makeActionSetTargetStation(timeRealSec, simTimeDays, "TimeTrial", (core::u64)st.anchorStationId));
          }

          if (st.autoEngageDockingComputerOnDocking) {
            st.requestEngageDockingComputer = true;
            if (emitAction) emitAction(makeActionEngageDockingComputer(timeRealSec, simTimeDays, "TimeTrial", true));
          }

          if (st.autoCameraOrbitOnDocking) {
            st.requestCameraOrbit = true;
            if (emitAction) {
              // CameraRigMode::Orbit == 1
              emitAction(makeActionSetCameraRigPreset(timeRealSec, simTimeDays, "TimeTrial", 2)); // Docking preset
            }
          }

          emitTTEvent("Docking", "Gates complete; dock at anchor to finish", (core::u64)st.anchorStationId);
          if (toast) toast("Gates complete. Dock at the anchor station to finish.", 2.4);
        }
      }
    }
  }

  // Docking completion.
  if (st.phase == TimeTrialPhase::Docking && st.finishMode == TimeTrialFinishMode::DockAtAnchorStation) {
    if (docked && st.anchorStationId != 0 && dockedStationId == st.anchorStationId) {
      // Docking is detected from the *post-physics* ship state, so treat the
      // finish time as the end of this tick.
      const double finishAt = st.timeSec + dtRealSec;
      st.phase = TimeTrialPhase::Finished;
      st.finishTimeSec = finishAt;
      st.timeSec = finishAt;

      // Ghost recording: ensure we capture a terminal sample for the run.
      if (st.ghostRecordRun) {
        pushGhostSample(st, st.finishTimeSec, simTimeDays, ship);
      }

      const double prevBest = (st.bestTimesSec.count(st.course.key) > 0)
                                ? st.bestTimesSec[st.course.key]
                                : 0.0;
      const bool isBest = !(prevBest > 0.0) || (st.finishTimeSec < prevBest);
      if (isBest) {
        st.bestTimesSec[st.course.key] = st.finishTimeSec;
        st.bestTimeSec = st.finishTimeSec;

        if (st.ghostSaveBestRun && st.ghostRecordRun && !st.ghostRunSamples.empty()) {
          st.bestGhostSamples[st.course.key] = st.ghostRunSamples;
        }
        if (toast) toast("Docked! New best time: " + formatTime(st.finishTimeSec), 2.8);
        emitTTEvent("DockFinish", "Docked at anchor (new best)", (core::u64)(st.finishTimeSec * 1000.0));
      } else {
        st.bestTimeSec = prevBest;
        if (toast) toast("Docked: " + formatTime(st.finishTimeSec) + " (best " + formatTime(prevBest) + ")", 2.6);
        emitTTEvent("DockFinish", "Docked at anchor", (core::u64)(st.finishTimeSec * 1000.0));
      }

      // Run capture: stop recording and optionally export traces.
      if (emitAction) {
        if (st.autoStopFlightRecorderOnFinish) {
          emitAction(makeActionStopFlightRecorder(timeRealSec, simTimeDays, "TimeTrial"));
        }
        if (st.autoExportTracesOnFinish) {
          const core::u64 stampMs = (core::u64)std::llround(timeRealSec * 1000.0);
          char flightPath[256];
          char integPath[256];
          std::snprintf(flightPath, sizeof(flightPath), "time_trial_%llu_%llu_flight_trace.json",
                        (unsigned long long)st.course.key, (unsigned long long)stampMs);
          std::snprintf(integPath, sizeof(integPath), "time_trial_%llu_%llu_integration_trace.json",
                        (unsigned long long)st.course.key, (unsigned long long)stampMs);
          emitAction(makeActionExportFlightRecorderTrace(timeRealSec, simTimeDays, "TimeTrial", flightPath));
          emitAction(makeActionExportIntegrationTrace(timeRealSec, simTimeDays, "TimeTrial", integPath));
        }
      }
    }
  }

  // Timer advance.
  if (st.phase == TimeTrialPhase::Running || st.phase == TimeTrialPhase::Docking) {
    if (startedThisFrame) {
      // Only count the sub-frame portion after the start gate crossing.
      st.timeSec += dtRealSec * (1.0 - std::clamp(startFrac, 0.0, 1.0));
    } else {
      st.timeSec += dtRealSec;
    }
  }

  // Ghost recording: sample ship pose while the timer is running.
  // This buffer can be promoted to a "best run" ghost on finish.
  if (st.ghostRecordRun && (st.phase == TimeTrialPhase::Running || st.phase == TimeTrialPhase::Docking)) {
    const double hz = std::clamp(st.ghostSampleHz, 5.0, 240.0);
    const double step = 1.0 / hz;
    const double dtGhost = startedThisFrame
                          ? (dtRealSec * (1.0 - std::clamp(startFrac, 0.0, 1.0)))
                          : dtRealSec;
    st.ghostSampleAccumulatorSec += std::max(0.0, dtGhost);

    // If the buffer is empty (e.g., user toggled recording mid-run), seed it.
    if (st.ghostRunSamples.empty()) {
      pushGhostSample(st, st.timeSec, simTimeDays, ship);
    } else if (st.ghostSampleAccumulatorSec >= step) {
      // Keep the remainder so we don't drift too hard under varying frame times.
      st.ghostSampleAccumulatorSec = std::fmod(st.ghostSampleAccumulatorSec, step);
      pushGhostSample(st, st.timeSec, simTimeDays, ship);
    }
  }

  // Update best-run ghost pose for rendering / HUD.
  const double tNow = (st.phase == TimeTrialPhase::Finished) ? st.finishTimeSec : st.timeSec;
  updateGhostSample(tNow);


  st.prevShipPosKm = posKm;
  st.hasPrevShipPos = true;
}

void drawTimeTrialWindow(TimeTrialWindowState& st,
                         const sim::StarSystem* sys,
                         double timeDays,
                         const sim::Ship& ship,
                         const ToastFn& toast) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(520, 520), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Time Trials", &st.open)) {
    ImGui::End();
    return;
  }

  ImGui::TextDisabled(
    "Deterministic gate courses for hands-on 3D gameplay (racing / training).\n"
    "Fly through Gate 1 to start. Finish mode can be last-gate, or last-gate + dock at station.");
  ImGui::Separator();

  const bool hasSys = (sys != nullptr);
  const bool hasStations = hasSys && !sys->stations.empty();

  if (!hasStations) {
    ImGui::TextColored(ImVec4(1,0.8f,0.2f,1), "No stations available in the current system.");
    ImGui::TextDisabled("Tip: open a system with stations, then generate a course.");
    ImGui::End();
    return;
  }

  // ---- Generation ----
  ImGui::Text("Course Generation");
  ImGui::Indent();

  // Station selector.
  st.stationIndex = std::clamp(st.stationIndex, 0, (int)sys->stations.size() - 1);
  if (ImGui::BeginCombo("Anchor station", sys->stations[st.stationIndex].name.c_str())) {
    for (int i = 0; i < (int)sys->stations.size(); ++i) {
      const bool selected = (i == st.stationIndex);
      if (ImGui::Selectable(sys->stations[i].name.c_str(), selected)) {
        st.stationIndex = i;
      }
      if (selected) ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }

  // Generator knobs.
  ImGui::SliderInt("Gate count", &st.params.gateCount, 4, 64);
  const double kMinGateRadiusKm = 250.0;
  const double kMaxGateRadiusKm = 15000.0;
  ImGui::SliderScalar("Gate radius (km)", ImGuiDataType_Double, &st.params.gateRadiusKm,
                      &kMinGateRadiusKm, &kMaxGateRadiusKm, "%.0f");

  const double kMinCourseRadiusKm = 5000.0;
  const double kMaxCourseRadiusKm = 250000.0;
  ImGui::SliderScalar("Course radius (km)", ImGuiDataType_Double, &st.params.courseRadiusKm,
                      &kMinCourseRadiusKm, &kMaxCourseRadiusKm, "%.0f");

  const double kMinCourseHeightKm = 0.0;
  const double kMaxCourseHeightKm = 150000.0;
  ImGui::SliderScalar("Course height (km)", ImGuiDataType_Double, &st.params.courseHeightKm,
                      &kMinCourseHeightKm, &kMaxCourseHeightKm, "%.0f");

  const double kMinJitterKm = 0.0;
  const double kMaxJitterKm = 80000.0;
  ImGui::SliderScalar("Jitter (km)", ImGuiDataType_Double, &st.params.jitterKm,
                      &kMinJitterKm, &kMaxJitterKm, "%.0f");
  ImGui::SliderInt("Loops", &st.params.loops, 1, 8);
  ImGui::Checkbox("Closed loop", &st.params.closedLoop);

  ImGui::Separator();
  ImGui::Checkbox("Seed from system", &st.seedFromSystem);
  ImGui::InputScalar("User seed", ImGuiDataType_U64, &st.userSeed);

  ImGui::Separator();
  ImGui::Text("Finish condition");
  const int fm = (int)st.finishMode;
  if (ImGui::RadioButton("Finish at last gate", fm == (int)TimeTrialFinishMode::GatesOnly)) {
    st.finishMode = TimeTrialFinishMode::GatesOnly;
  }
  if (ImGui::RadioButton("Finish by docking at anchor station", fm == (int)TimeTrialFinishMode::DockAtAnchorStation)) {
    st.finishMode = TimeTrialFinishMode::DockAtAnchorStation;
  }
  if (st.finishMode == TimeTrialFinishMode::DockAtAnchorStation) {
    ImGui::Indent();
    ImGui::Checkbox("Auto-target anchor when docking", &st.autoTargetAnchorOnDocking);
    ImGui::Checkbox("Auto-engage docking computer", &st.autoEngageDockingComputerOnDocking);
    ImGui::Checkbox("Auto camera preset: Travel on start", &st.autoCameraChaseOnStart);
    ImGui::Checkbox("Auto camera preset: Docking on docking", &st.autoCameraOrbitOnDocking);
    ImGui::TextDisabled("Timer continues after the last gate until you dock.");
    ImGui::Unindent();
  }

  ImGui::Separator();
  ImGui::Text("Run capture");
  ImGui::Indent();
  ImGui::Checkbox("Clear Integration Hub log on arm", &st.autoClearIntegrationHubOnArm);
  ImGui::Checkbox("Auto start flight recorder on start gate", &st.autoStartFlightRecorderOnStart);
  ImGui::Checkbox("Auto stop flight recorder on finish", &st.autoStopFlightRecorderOnFinish);
  ImGui::Checkbox("Auto export traces on finish", &st.autoExportTracesOnFinish);
  ImGui::TextDisabled("Exports: time_trial_<courseKey>_<timestamp>_{flight,integration}_trace.json");
  ImGui::Unindent();

  ImGui::Separator();
  ImGui::Text("Ghost (best run)");
  ImGui::Indent();
  ImGui::Checkbox("Enable ghost", &st.ghostEnabled);
  ImGui::SameLine();
  ImGui::Checkbox("Draw ship", &st.ghostDrawShip);
  ImGui::SameLine();
  ImGui::Checkbox("Trail", &st.ghostDrawTrail);
  ImGui::Checkbox("Show split in Objective HUD", &st.ghostShowSplitHud);
  ImGui::Checkbox("Record runs", &st.ghostRecordRun);
  ImGui::SameLine();
  ImGui::Checkbox("Save on new best", &st.ghostSaveBestRun);

  // Visual dial (implemented as color scaling in the renderer).
  ImGui::SliderFloat("Ghost opacity", &st.ghostOpacity, 0.05f, 1.0f, "%.2f");

  float hz = (float)st.ghostSampleHz;
  if (ImGui::DragFloat("Record Hz", &hz, 1.0f, 5.0f, 120.0f, "%.0f")) {
    st.ghostSampleHz = std::clamp((double)hz, 5.0, 240.0);
  }

  float lead = (float)st.ghostLeadSec;
  if (ImGui::DragFloat("Lead / lag (s)", &lead, 0.05f, -10.0f, 10.0f, "%.2f")) {
    st.ghostLeadSec = (double)lead;
  }

  if (st.hasCourse) {
    const auto itG = st.bestGhostSamples.find(st.course.key);
    const int ghostCount = (itG != st.bestGhostSamples.end()) ? (int)itG->second.size() : 0;
    ImGui::TextDisabled("Best ghost samples: %d", ghostCount);

    if (ImGui::SmallButton("Reset best time + ghost (this course)")) {
      st.bestTimesSec.erase(st.course.key);
      st.bestTimeSec = 0.0;
      st.bestGhostSamples.erase(st.course.key);
      st.ghostSampleValid = false;
      if (toast) toast("Best time + ghost reset for this course.", 1.6);
    }
  } else {
    ImGui::TextDisabled("Finish at least one run to record a ghost.");
  }

  ImGui::TextDisabled("Tip: the ghost is the best time you've set for this exact course key.");
  ImGui::Unindent();

  const sim::Station& anchor = sys->stations[st.stationIndex];
  const core::u64 seed = computeCourseSeed(*sys, anchor, st.userSeed, st.seedFromSystem);
  const std::string code = courseCode(*sys, anchor, seed, st.params);

  ImGui::TextDisabled("Course code:");
  ImGui::SameLine();
  ImGui::Text("%s", code.c_str());
  ImGui::SameLine();
  if (ImGui::SmallButton("Copy")) {
    ImGui::SetClipboardText(code.c_str());
    if (toast) toast("Course code copied.", 1.4);
  }

  const bool canGenerate = hasSys;
  if (ImGui::Button("Generate course")) {
    if (!canGenerate) {
      if (toast) toast("No current system.", 1.8);
    } else {
      const math::Vec3d anchorPosKm = sim::stationPosKm(anchor, timeDays);
      // Keep the anchor orientation neutral for now (world axes). The slalom generator uses it
      // only as a local frame for building offsets.
      const math::Quatd anchorOrient = math::Quatd::identity();

      st.course = sim::generateTimeTrialCourseStationSlalomKm(anchorPosKm, anchorOrient, seed, st.params);
      st.hasCourse = !st.course.gates.empty();
      st.anchorStationId = anchor.id;
      st.anchorStationIndex = st.stationIndex;
      st.bestTimeSec = (st.bestTimesSec.count(st.course.key) > 0) ? st.bestTimesSec[st.course.key] : 0.0;
      st.armCourse();

      if (toast) toast("Generated course: " + anchor.name + " (" + std::to_string((int)st.course.gates.size()) + " gates)", 2.0);
    }
  }

  ImGui::SameLine();
  if (ImGui::Button("Clear course")) {
    st.clearCourse();
    if (toast) toast("Time Trial cancelled.", 1.3);
  }

  ImGui::Unindent();

  // ---- Status ----
  ImGui::Separator();
  ImGui::Text("Run Status");
  ImGui::Indent();

  auto phaseLabel = [&]() -> const char* {
    switch (st.phase) {
      case TimeTrialPhase::Inactive: return "Inactive";
      case TimeTrialPhase::Ready:    return "Ready (cross Gate 1 to start)";
      case TimeTrialPhase::Running:  return "Running";
      case TimeTrialPhase::Docking:  return "Docking (dock at anchor station)";
      case TimeTrialPhase::Finished: return "Finished";
      default: return "?";
    }
  };

  ImGui::Text("Phase: %s", phaseLabel());

  const int totalGates = (int)st.course.gates.size();
  const int next = std::clamp(st.nextGate, 0, std::max(0, totalGates));
  ImGui::Text("Progress: %d / %d", next, totalGates);

  const double best = st.bestTimeSec;
  ImGui::Text("Time: %s", formatTime((st.phase == TimeTrialPhase::Finished) ? st.finishTimeSec : st.timeSec).c_str());
  ImGui::Text("Best: %s", formatTime(best).c_str());

  if (st.phase == TimeTrialPhase::Finished) {
    if (ImGui::Button("Re-arm course")) {
      st.armCourse();
      if (toast) toast("Ready: fly through Gate 1 to start.", 1.8);
    }
  } else if (st.phase != TimeTrialPhase::Inactive) {
    ImGui::SameLine();
    if (ImGui::Button("Cancel run")) {
      st.cancelRun();
      if (toast) toast("Run cancelled.", 1.4);
    }
  }

  // Helpful distance readout.
  if (st.hasCourse) {
    if ((st.phase == TimeTrialPhase::Ready || st.phase == TimeTrialPhase::Running) && st.nextGate < (int)st.course.gates.size()) {
      const auto& g = st.course.gates[st.nextGate];
      const double distKm = (g.posKm - ship.positionKm()).length();
      ImGui::TextDisabled("Next gate distance: %.0f km", distKm);
    } else if (st.phase == TimeTrialPhase::Docking && sys && st.anchorStationId != 0) {
      const sim::Station* anchorSt = nullptr;
      for (const auto& s : sys->stations) {
        if (s.id == st.anchorStationId) { anchorSt = &s; break; }
      }
      if (anchorSt) {
        const math::Vec3d anchorPosKm = sim::stationPosKm(*anchorSt, timeDays);
        const double distKm = (anchorPosKm - ship.positionKm()).length();
        ImGui::TextDisabled("Anchor station distance: %.0f km", distKm);
      }
    }
  }

  ImGui::Unindent();

  // ---- HUD toggles ----
  ImGui::Separator();
  ImGui::Text("HUD");
  ImGui::Indent();
  ImGui::Checkbox("Enable Time Trial HUD", &st.hudEnabled);
  ImGui::Checkbox("Show gate marker", &st.showGateMarker);
  ImGui::Checkbox("Show all gates (debug)", &st.showAllGates);
  ImGui::Checkbox("Clamp offscreen", &st.clampOffscreen);
  ImGui::SliderFloat("Marker alpha", &st.markerAlpha, 0.1f, 1.0f, "%.2f");
  ImGui::SliderFloat("Marker thickness", &st.markerThickness, 1.0f, 5.0f, "%.1f");
  ImGui::Unindent();

  ImGui::End();
}

} // namespace stellar::game
