#include "FlightRecorderWindow.h"

#include "stellar/core/JsonWriter.h"

#include "GameSignals.h"
#include "stellar/sim/Ship.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace stellar::game {

static double speedKmS(const FlightRecorderSample& s) {
  return s.velKmS.length();
}

static double recordingDurationSec(const FlightRecorderWindowState& st) {
  if (st.samples.empty()) return 0.0;
  return std::max(0.0, st.samples.back().tRealSec);
}


static bool sampleAtTime(const std::deque<FlightRecorderSample>& samples,
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

static void pushSample(FlightRecorderWindowState& st,
                       double timeRealSec,
                       double simTimeDays,
                       const sim::Ship& ship) {
  FlightRecorderSample s;
  s.tRealSec = timeRealSec - st.startRealSec;
  s.tSimDays = simTimeDays;
  s.posKm = ship.positionKm();
  s.velKmS = ship.velocityKmS();
  s.orient = ship.orientation();
  s.angVelRadS = ship.angularVelocityRadS();

  // Enforce ring buffer.
  if (st.maxSamples <= 0) st.maxSamples = 1;
  while ((int)st.samples.size() >= st.maxSamples) {
    st.samples.pop_front();
  }
  st.samples.push_back(s);
}



static void pushMarker(FlightRecorderWindowState& st, const FlightRecorderMarker& mIn) {
  FlightRecorderMarker m = mIn;

  // Enforce ring buffer.
  if (st.maxMarkers <= 0) st.maxMarkers = 1;
  while ((int)st.markers.size() >= st.maxMarkers) {
    st.markers.pop_front();
  }
  st.markers.push_back(std::move(m));
}

void flightRecorderCaptureMarkerFromEvent(FlightRecorderWindowState& st, const GameEvent& ev) {
  if (!st.captureMarkers) return;
  if (!st.recording && !st.captureMarkersWhileNotRecording) return;

  // Avoid flooding the marker track with ultra-verbose debug events unless requested.
  if (ev.kind == GameEventKind::Debug && !st.captureMarkersIncludeDebug) return;

  FlightRecorderMarker m;
  m.tRealSec = ev.tRealSec;
  m.tSimDays = ev.tSimDays;
  m.kind = ev.kind;
  m.tag = ev.tag;
  m.msg = ev.msg;
  m.u64a = ev.u64a;
  m.u64b = ev.u64b;

  pushMarker(st, m);
}

void flightRecorderClearMarkers(FlightRecorderWindowState& st) {
  st.markers.clear();
}

void startFlightRecorderSession(FlightRecorderWindowState& st,
                               double timeRealSec,
                               double simTimeDays,
                               const sim::Ship& ship) {
  st.recording = true;
  st.sampleAccumulatorSec = 0.0;
  st.startRealSec = timeRealSec;
  st.startSimDays = simTimeDays;
  st.samples.clear();
  st.markers.clear();

  // Reset replay to a sane baseline (but keep ghost toggles).
  st.playing = false;
  st.playheadSec = 0.0;
  st.ghostSampleValid = false;

  // Reset time bookkeeping so replay dt doesn't spike on the first tick.
  st.hasPrevTimes = false;
  st.prevRealSec = timeRealSec;
  st.prevSimDays = simTimeDays;

  // Capture an initial sample immediately.
  pushSample(st, timeRealSec, simTimeDays, ship);
  st.ghostSampleValid = sampleAtTime(st.samples, 0.0, &st.ghostSample);
}

void stopFlightRecorderSession(FlightRecorderWindowState& st) {
  st.recording = false;

  // When stopping, snap playhead to the end so the ghost reflects the last recorded pose.
  st.playing = false;
  st.playheadSec = recordingDurationSec(st);
  st.ghostSampleValid = sampleAtTime(st.samples, st.playheadSec, &st.ghostSample);
}

void flightRecorderSeekToRealTime(FlightRecorderWindowState& st, double timeRealSec) {
  if (st.samples.size() < 2) return;

  const double dur = recordingDurationSec(st);
  if (!(dur > 0.0)) return;

  // Map absolute -> session-relative time.
  double t = timeRealSec - st.startRealSec;
  if (!std::isfinite(t)) t = 0.0;
  t = std::clamp(t, 0.0, dur);

  st.playheadSec = t;
  st.playing = false;
  st.ghostSampleValid = sampleAtTime(st.samples, st.playheadSec, &st.ghostSample);
}

static bool exportCsv(const char* path,
                      const std::deque<FlightRecorderSample>& samples,
                      std::string* err) {
  if (!path || !*path) {
    if (err) *err = "Invalid CSV path.";
    return false;
  }
  if (samples.empty()) {
    if (err) *err = "No samples to export.";
    return false;
  }

  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f.is_open()) {
    if (err) *err = "Failed to open CSV file for writing.";
    return false;
  }

  f.setf(std::ios::fixed);
  f.precision(9);

  // Note: new columns appended at the end for backward compatibility with existing
  // analysis scripts that parse the initial set.
  f << "t_real_sec,t_sim_days,"
       "pos_x_km,pos_y_km,pos_z_km,"
       "vel_x_kms,vel_y_kms,vel_z_kms,"
       "speed_kms,"
       "orient_w,orient_x,orient_y,orient_z,"
       "ang_x_rads,ang_y_rads,ang_z_rads\n";

  for (const auto& s : samples) {
    f << s.tRealSec << ','
      << s.tSimDays << ','
      << s.posKm.x << ',' << s.posKm.y << ',' << s.posKm.z << ','
      << s.velKmS.x << ',' << s.velKmS.y << ',' << s.velKmS.z << ','
      << speedKmS(s) << ','
      << s.orient.w << ',' << s.orient.x << ',' << s.orient.y << ',' << s.orient.z << ','
      << s.angVelRadS.x << ',' << s.angVelRadS.y << ',' << s.angVelRadS.z
      << "\n";
  }

  return true;
}


static void writeProcessNameMetadata(core::JsonWriter& w, int pid, const char* processName) {
  w.beginObject();
  w.key("name"); w.value("process_name");
  w.key("ph"); w.value("M");
  w.key("pid"); w.value(pid);
  w.key("tid"); w.value(0);
  w.key("args");
  w.beginObject();
  w.key("name"); w.value(processName);
  w.endObject();
  w.endObject();
}

static void writeThreadNameMetadata(core::JsonWriter& w, int pid, int tid, const char* threadName) {
  w.beginObject();
  w.key("name"); w.value("thread_name");
  w.key("ph"); w.value("M");
  w.key("pid"); w.value(pid);
  w.key("tid"); w.value(tid);
  w.key("args");
  w.beginObject();
  w.key("name"); w.value(threadName);
  w.endObject();
  w.endObject();
}

static void writeCounterEvent(core::JsonWriter& w,
                              const std::vector<std::string>& keys,
                              const double* values,
                              std::uint64_t tsUs,
                              int pid,
                              int tid) {
  w.beginObject();
  w.key("name"); w.value("ship");
  w.key("cat"); w.value("flight");
  w.key("ph"); w.value("C");
  w.key("ts"); w.value((unsigned long long)tsUs);
  w.key("pid"); w.value(pid);
  w.key("tid"); w.value(tid);
  w.key("args");
  w.beginObject();
  for (std::size_t i = 0; i < keys.size(); ++i) {
    w.key(keys[i]);
    w.value(values ? values[i] : 0.0);
  }
  w.endObject();
  w.endObject();
}

static void writeInstantMarkerEvent(core::JsonWriter& w,
                                   const FlightRecorderMarker& m,
                                   std::uint64_t tsUs,
                                   int pid,
                                   int tid) {
  w.beginObject();
  w.key("name");
  w.value(m.tag.empty() ? std::string_view{"marker"} : std::string_view{m.tag});
  w.key("cat"); w.value("marker");
  w.key("ph"); w.value("i");
  w.key("s"); w.value("t"); // thread scope
  w.key("ts"); w.value((unsigned long long)tsUs);
  w.key("pid"); w.value(pid);
  w.key("tid"); w.value(tid);
  w.key("args");
  w.beginObject();
  w.key("kind"); w.value(gameEventKindName(m.kind));
  if (!m.msg.empty()) { w.key("msg"); w.value(m.msg); }
  if (m.u64a != 0) { w.key("u64a"); w.value((unsigned long long)m.u64a); }
  if (m.u64b != 0) { w.key("u64b"); w.value((unsigned long long)m.u64b); }
  w.endObject();
  w.endObject();
}

static bool exportTrace(const char* path,
                        const FlightRecorderWindowState& st,
                        std::string* err) {
  if (!path || !*path) {
    if (err) *err = "Invalid trace path.";
    return false;
  }
  if (st.samples.empty()) {
    if (err) *err = "No samples to export.";
    return false;
  }

  // Build counter keys.
  std::vector<std::string> keys;
  keys.reserve(24);

  auto add = [&](const char* k) { keys.emplace_back(k); };

  if (st.traceIncludeSimTime) add("sim_days");
  if (st.traceIncludeSpeed) add("speed_kms");
  if (st.traceIncludePosition) {
    add("pos_x_km");
    add("pos_y_km");
    add("pos_z_km");
  }
  if (st.traceIncludeVelocity) {
    add("vel_x_kms");
    add("vel_y_kms");
    add("vel_z_kms");
  }
  if (st.traceIncludeOrientation) {
    add("q_w");
    add("q_x");
    add("q_y");
    add("q_z");
  }
  if (st.traceIncludeAngularVelocity) {
    add("ang_x_rads");
    add("ang_y_rads");
    add("ang_z_rads");
  }

  if (keys.empty()) {
    if (err) *err = "No trace tracks selected.";
    return false;
  }

  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f.is_open()) {
    if (err) *err = "Failed to open trace file for writing.";
    return false;
  }

  const int pid = 1;
  const int tidCounter = 1;
  const int tidMarkers = 2;

  core::JsonWriter w(f, st.tracePretty);
  w.beginObject();
  w.key("displayTimeUnit");
  w.value("ms");
  w.key("traceEvents");
  w.beginArray();

  // Metadata.
  writeProcessNameMetadata(w, pid, "stellar");
  writeThreadNameMetadata(w, pid, tidCounter, "Flight");
  if (st.traceIncludeMarkers) {
    writeThreadNameMetadata(w, pid, tidMarkers, "Markers");
  }

  // Write counter events.
  std::vector<double> row;
  row.resize(keys.size());

  for (const auto& s : st.samples) {
    const std::uint64_t tsUs = (std::uint64_t)std::llround(s.tRealSec * 1e6);

    std::size_t col = 0;
    if (st.traceIncludeSimTime) row[col++] = s.tSimDays;
    if (st.traceIncludeSpeed) row[col++] = speedKmS(s);
    if (st.traceIncludePosition) {
      row[col++] = s.posKm.x;
      row[col++] = s.posKm.y;
      row[col++] = s.posKm.z;
    }
    if (st.traceIncludeVelocity) {
      row[col++] = s.velKmS.x;
      row[col++] = s.velKmS.y;
      row[col++] = s.velKmS.z;
    }
    if (st.traceIncludeOrientation) {
      row[col++] = s.orient.w;
      row[col++] = s.orient.x;
      row[col++] = s.orient.y;
      row[col++] = s.orient.z;
    }
    if (st.traceIncludeAngularVelocity) {
      row[col++] = s.angVelRadS.x;
      row[col++] = s.angVelRadS.y;
      row[col++] = s.angVelRadS.z;
    }

    writeCounterEvent(w, keys, row.data(), tsUs, pid, tidCounter);
  }

  // Write markers as instant events (optional).
  if (st.traceIncludeMarkers && !st.markers.empty()) {
    for (const auto& m : st.markers) {
      double rel = m.tRealSec - st.startRealSec;
      if (!std::isfinite(rel)) rel = 0.0;
      if (rel < 0.0) rel = 0.0;
      const std::uint64_t tsUs = (std::uint64_t)std::llround(rel * 1e6);
      writeInstantMarkerEvent(w, m, tsUs, pid, tidMarkers);
    }
  }

  w.endArray();
  w.endObject();
  return true;
}


bool exportFlightRecorderCsv(const FlightRecorderWindowState& st, const char* path, std::string* err) {
  return exportCsv(path, st.samples, err);
}

bool exportFlightRecorderTraceJson(const FlightRecorderWindowState& st, const char* path, std::string* err) {
  return exportTrace(path, st, err);
}

static void tickReplay(FlightRecorderWindowState& st,
                       double dtRealSec,
                       double dtSimSec,
                       bool paused) {
  const double dur = recordingDurationSec(st);
  if (dur <= 1e-9 || st.samples.empty()) {
    st.playing = false;
    st.playheadSec = 0.0;
    st.ghostSampleValid = false;
    return;
  }

  // Advance playhead.
  if (st.playing) {
    if (!(paused && !st.playWhilePaused)) {
      double dt = st.playUsesSimTime ? dtSimSec : dtRealSec;
      if (paused && st.playWhilePaused) {
        // When paused, sim dt is typically 0. Use real dt so the ghost can still scrub/animate.
        dt = dtRealSec;
      }
      dt *= st.playbackRate;

      if (dt != 0.0) {
        st.playheadSec += dt;
      }
    }
  }

  st.playheadSec = std::max(0.0, st.playheadSec);

  if (st.loop) {
    st.playheadSec = std::fmod(st.playheadSec, dur);
    if (st.playheadSec < 0.0) st.playheadSec += dur;
  } else {
    if (st.playheadSec >= dur) {
      st.playheadSec = dur;
      st.playing = false;
    }
  }

  st.ghostSampleValid = sampleAtTime(st.samples, st.playheadSec, &st.ghostSample);
}

void tickFlightRecorder(FlightRecorderWindowState& st,
                        double dtRealSec,
                        double timeRealSec,
                        double simTimeDays,
                        const sim::Ship& ship,
                        bool paused) {
  // Derive simulation dt (seconds) for replay, even if recording is off.
  double dtSimSec = 0.0;
  if (st.hasPrevTimes) {
    dtSimSec = (simTimeDays - st.prevSimDays) * 86400.0;
    if (!std::isfinite(dtSimSec)) dtSimSec = 0.0;
    // Guard against load/warp discontinuities.
    if (dtSimSec < 0.0) dtSimSec = 0.0;
    if (dtSimSec > 60.0) dtSimSec = 60.0;
  }
  st.prevRealSec = timeRealSec;
  st.prevSimDays = simTimeDays;
  st.hasPrevTimes = true;

  // --- Recording ---
  if (st.recording) {
    if (!(paused && !st.recordWhilePaused)) {
      // Guard against invalid sample rates.
      if (st.sampleHz <= 0.0) st.sampleHz = 20.0;

      const double intervalSec = 1.0 / st.sampleHz;
      st.sampleAccumulatorSec += dtRealSec;

      // We only sample at most once per frame (we don't attempt to reconstruct
      // intermediate ship states).
      if (st.sampleAccumulatorSec >= intervalSec) {
        st.sampleAccumulatorSec = 0.0;
        pushSample(st, timeRealSec, simTimeDays, ship);
      }
    }
  }

  // --- Replay ---
  tickReplay(st, dtRealSec, dtSimSec, paused);
}

void drawFlightRecorderWindow(FlightRecorderWindowState& st,
                              const sim::Ship& ship,
                              double timeRealSec,
                              double simTimeDays,
                              bool paused,
                              const ToastFn& toast) {
  if (!st.open) return;

  if (!ImGui::Begin("Flight Recorder", &st.open)) {
    ImGui::End();
    return;
  }

  ImGui::TextUnformatted("Record ship telemetry and export to CSV or Perfetto/Chrome trace.");
  ImGui::TextDisabled("Tip: Perfetto UI can open Chrome JSON traces (ui.perfetto.dev). ");
  ImGui::TextDisabled("New: Replay a recording as a debug ghost (pose + angular velocity).");

  ImGui::Separator();

  // ---- Controls ----
  if (!st.recording) {
    if (ImGui::Button("Start recording")) {
      startFlightRecorderSession(st, timeRealSec, simTimeDays, ship);
      toast("Flight recorder started.", 1.4);
    }
  } else {
    if (ImGui::Button("Stop recording")) {
      stopFlightRecorderSession(st);
      toast("Flight recorder stopped.", 1.4);
    }
  }
  ImGui::SameLine();
  if (ImGui::Button("Clear")) {
    const bool wasRecording = st.recording;

    st.samples.clear();
    st.markers.clear();
    st.sampleAccumulatorSec = 0.0;
    st.startRealSec = timeRealSec;
    st.startSimDays = simTimeDays;

    st.playing = false;
    st.playheadSec = 0.0;
    st.ghostSampleValid = false;

    // Reset time bookkeeping so replay dt doesn't spike after a clear.
    st.hasPrevTimes = false;
    st.prevRealSec = timeRealSec;
    st.prevSimDays = simTimeDays;

    if (wasRecording) {
      pushSample(st, timeRealSec, simTimeDays, ship);
      st.ghostSampleValid = sampleAtTime(st.samples, 0.0, &st.ghostSample);
    }
    toast("Flight recorder cleared.", 1.2);
  }

  ImGui::SameLine();
  ImGui::Checkbox("Record while paused", &st.recordWhilePaused);
  if (paused && !st.recordWhilePaused) {
    ImGui::SameLine();
    ImGui::TextDisabled("(paused)");
  }

  const double minHz = 1.0;
  const double maxHz = 120.0;
  ImGui::SliderScalar("Sample rate (Hz)", ImGuiDataType_Double, &st.sampleHz, &minHz, &maxHz, "%.1f", ImGuiSliderFlags_Logarithmic);

  ImGui::SliderInt("Max samples", &st.maxSamples, 500, 200000, "%d");
  ImGui::Checkbox("Show plot", &st.showPlot);

  ImGui::Separator();

  // ---- Status ----
  const math::Vec3d curPos = ship.positionKm();
  const math::Vec3d curVel = ship.velocityKmS();
  const double curSpeed = curVel.length();

  ImGui::Text("Current speed: %.3f km/s", curSpeed);
  ImGui::Text("Current pos (km): (%.1f, %.1f, %.1f)", curPos.x, curPos.y, curPos.z);

  if (!st.samples.empty()) {
    const FlightRecorderSample& last = st.samples.back();
    ImGui::Text("Recorded: %zu samples | duration %.2f s", st.samples.size(), last.tRealSec);

    if (st.samples.size() >= 2) {
      const FlightRecorderSample& prev = *(++st.samples.rbegin());
      const double dt = std::max(1e-9, last.tRealSec - prev.tRealSec);
      const math::Vec3d dv = last.velKmS - prev.velKmS;
      const double acc = dv.length() / dt;
      ImGui::Text("Approx accel: %.3f km/s^2", acc);
    }
  } else {
    ImGui::TextDisabled("No samples recorded yet.");
  }

  // ---- Plot ----
  if (st.showPlot) {
    if (st.samples.size() >= 2) {
      static std::vector<float> values;
      values.clear();
      values.reserve(st.samples.size());

      float maxV = 0.0f;
      for (const auto& s : st.samples) {
        const float v = (float)speedKmS(s);
        values.push_back(v);
        maxV = std::max(maxV, v);
      }

      ImGui::PlotLines("Speed (km/s)", values.data(), (int)values.size(), 0, nullptr, 0.0f, maxV * 1.10f, ImVec2(0, 90));
    } else {
      ImGui::TextDisabled("Need at least 2 samples to plot.");
    }
  }

  // ---- Markers ----
  if (ImGui::CollapsingHeader("Markers (Integration Hub)", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Checkbox("Capture markers", &st.captureMarkers);
    ImGui::SameLine();
    ImGui::Checkbox("Include Debug", &st.captureMarkersIncludeDebug);
    ImGui::SameLine();
    ImGui::Checkbox("Capture while not recording", &st.captureMarkersWhileNotRecording);

    ImGui::PushItemWidth(160.0f);
    ImGui::DragInt("Max markers", &st.maxMarkers, 1.0f, 64, 20000);
    ImGui::PopItemWidth();

    ImGui::SameLine();
    if (ImGui::Button("Clear markers")) {
      st.markers.clear();
      toast("Markers cleared.", 1.2);
    }

    ImGui::Text("Markers: %d", (int)st.markers.size());

    if (!st.markers.empty()) {
      const float height = 160.0f;
      if (ImGui::BeginChild("##flight_markers", ImVec2(0, height), true)) {
        for (int i = 0; i < (int)st.markers.size(); ++i) {
          const auto& m = st.markers[(std::size_t)i];
          const double rel = m.tRealSec - st.startRealSec;

          // Build a compact label.
          std::string label;
          label.reserve(128);
          label += "[";
          char buf[64];
          std::snprintf(buf, sizeof(buf), "%.3fs", rel);
          label += buf;
          label += "] ";
          label += gameEventKindName(m.kind);
          if (!m.tag.empty()) {
            label += " / ";
            label += m.tag;
          }
          if (!m.msg.empty()) {
            label += " : ";
            // Trim to avoid giant lines.
            const std::size_t kMax = 80;
            if (m.msg.size() <= kMax) label += m.msg;
            else { label += m.msg.substr(0, kMax); label += "..."; }
          }

          if (ImGui::Selectable(label.c_str(), false)) {
            flightRecorderSeekToRealTime(st, m.tRealSec);
          }
        }
        ImGui::EndChild();
      }

      ImGui::TextDisabled("Click a marker to scrub the playhead.");
    }
  }

  ImGui::Separator();

  // ---- Replay / Ghost ----
  ImGui::TextUnformatted("Replay / Ghost");

  const double dur = recordingDurationSec(st);
  if (st.samples.size() < 2 || dur <= 1e-9) {
    ImGui::TextDisabled("Record at least 2 samples to enable replay.");
  } else {
    ImGui::Checkbox("Show ghost", &st.ghostEnabled);
    ImGui::SameLine();
    ImGui::Checkbox("Wireframe ghost", &st.ghostWireframe);
    ImGui::SameLine();
    ImGui::Checkbox("Draw trail", &st.ghostDrawTrail);

    ImGui::Checkbox("Loop", &st.loop);
    ImGui::SameLine();
    ImGui::Checkbox("Advance using sim dt", &st.playUsesSimTime);
    ImGui::SameLine();
    ImGui::Checkbox("Play while paused", &st.playWhilePaused);

    ImGui::PushItemWidth(160.0f);
    float rate = (float)st.playbackRate;
    if (ImGui::DragFloat("Playback rate", &rate, 0.01f, 0.0f, 10.0f, "%.2fx")) {
      st.playbackRate = std::max(0.0, (double)rate);
    }
    ImGui::PopItemWidth();

    // Buttons.
    if (!st.playing) {
      if (ImGui::Button("Play")) st.playing = true;
    } else {
      if (ImGui::Button("Pause")) st.playing = false;
    }
    ImGui::SameLine();
    if (ImGui::Button("Rewind")) {
      st.playheadSec = 0.0;
      st.playing = false;
    }
    ImGui::SameLine();
    if (ImGui::Button("To end")) {
      st.playheadSec = dur;
      st.playing = false;
    }

    double tMin = 0.0;
    double tMax = std::max(0.001, dur);
    ImGui::SliderScalar("Playhead (s)", ImGuiDataType_Double, &st.playheadSec, &tMin, &tMax, "%.3f");

    if (st.ghostSampleValid) {
      const auto& gs = st.ghostSample;
      ImGui::Text("Ghost t=%.3fs | speed=%.3f km/s", gs.tRealSec, gs.velKmS.length());
      ImGui::Text("Ghost pos (km): (%.1f, %.1f, %.1f)", gs.posKm.x, gs.posKm.y, gs.posKm.z);
    }
  }

  ImGui::Separator();

  // ---- Export ----
  ImGui::TextUnformatted("Export");

  ImGui::InputText("CSV path", st.csvPath, (int)sizeof(st.csvPath));
  if (ImGui::Button("Export CSV")) {
    std::string err;
    const bool ok = exportFlightRecorderCsv(st, st.csvPath, &err);
    toast(ok ? (std::string("Wrote CSV: ") + st.csvPath) : (err.empty() ? "CSV export failed." : err), 2.2);
  }

  ImGui::Spacing();
  ImGui::InputText("Trace path", st.tracePath, (int)sizeof(st.tracePath));

  ImGui::Checkbox("Trace: include sim time", &st.traceIncludeSimTime);
  ImGui::Checkbox("Trace: include speed", &st.traceIncludeSpeed);
  ImGui::Checkbox("Trace: include position", &st.traceIncludePosition);
  ImGui::Checkbox("Trace: include velocity", &st.traceIncludeVelocity);
  ImGui::Checkbox("Trace: include orientation", &st.traceIncludeOrientation);
  ImGui::Checkbox("Trace: include angular velocity", &st.traceIncludeAngularVelocity);
  ImGui::Checkbox("Trace: include markers", &st.traceIncludeMarkers);
  ImGui::Checkbox("Trace: pretty JSON", &st.tracePretty);

  if (ImGui::Button("Export trace")) {
    std::string err;
    const bool ok = exportFlightRecorderTraceJson(st, st.tracePath, &err);
    toast(ok ? (std::string("Wrote trace: ") + st.tracePath) : (err.empty() ? "Trace export failed." : err), 2.4);
  }

  ImGui::End();
}

} // namespace stellar::game
