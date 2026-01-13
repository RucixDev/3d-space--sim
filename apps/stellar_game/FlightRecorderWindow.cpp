#include "FlightRecorderWindow.h"

#include "stellar/core/ChromeTrace.h"
#include "stellar/sim/Ship.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
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

static double quatDot(const math::Quatd& a, const math::Quatd& b) {
  return a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z;
}

// Normalized linear interpolation (nlerp) with shortest-arc handling.
// Good enough for smooth debug ghost playback without pulling in a full slerp.
static math::Quatd quatNlerp(const math::Quatd& a, const math::Quatd& bIn, double t) {
  math::Quatd b = bIn;
  if (quatDot(a, b) < 0.0) {
    b.w = -b.w; b.x = -b.x; b.y = -b.y; b.z = -b.z;
  }

  math::Quatd q{
    a.w * (1.0 - t) + b.w * t,
    a.x * (1.0 - t) + b.x * t,
    a.y * (1.0 - t) + b.y * t,
    a.z * (1.0 - t) + b.z * t
  };
  return q.normalized();
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
  s.orient = quatNlerp(a.orient, b.orient, u);
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

  core::ChromeTraceCounterTable table;
  table.name = "ship";
  table.category = "flight";
  table.keys = keys;

  const std::size_t n = st.samples.size();
  const std::size_t k = table.keys.size();

  table.tsUs.resize(n);
  table.values.resize(n * k);

  // Fill a row-major value table.
  std::size_t i = 0;
  for (const auto& s : st.samples) {
    const std::uint64_t ts = (std::uint64_t)std::llround(s.tRealSec * 1e6);
    table.tsUs[i] = ts;

    std::size_t col = 0;
    if (st.traceIncludeSimTime) table.values[i * k + col++] = s.tSimDays;
    if (st.traceIncludeSpeed) table.values[i * k + col++] = speedKmS(s);
    if (st.traceIncludePosition) {
      table.values[i * k + col++] = s.posKm.x;
      table.values[i * k + col++] = s.posKm.y;
      table.values[i * k + col++] = s.posKm.z;
    }
    if (st.traceIncludeVelocity) {
      table.values[i * k + col++] = s.velKmS.x;
      table.values[i * k + col++] = s.velKmS.y;
      table.values[i * k + col++] = s.velKmS.z;
    }
    if (st.traceIncludeOrientation) {
      table.values[i * k + col++] = s.orient.w;
      table.values[i * k + col++] = s.orient.x;
      table.values[i * k + col++] = s.orient.y;
      table.values[i * k + col++] = s.orient.z;
    }
    if (st.traceIncludeAngularVelocity) {
      table.values[i * k + col++] = s.angVelRadS.x;
      table.values[i * k + col++] = s.angVelRadS.y;
      table.values[i * k + col++] = s.angVelRadS.z;
    }

    ++i;
  }

  core::ChromeTraceWriteOptions opt;
  opt.includeFrameEvents = false;
  opt.pretty = st.tracePretty;
  opt.pid = 1;
  opt.tid = 1;

  return core::writeCounterChromeTraceJson(path, table, opt, err);
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
      st.recording = true;
      st.sampleAccumulatorSec = 0.0;
      st.startRealSec = timeRealSec;
      st.startSimDays = simTimeDays;
      st.samples.clear();

      // Reset replay to a sane baseline (but keep ghost toggles).
      st.playing = false;
      st.playheadSec = 0.0;

      // Capture an initial sample immediately.
      pushSample(st, timeRealSec, simTimeDays, ship);
      toast("Flight recorder started.", 1.4);
    }
  } else {
    if (ImGui::Button("Stop recording")) {
      st.recording = false;

      // When stopping, snap playhead to the end so the ghost reflects the last recorded pose.
      st.playing = false;
      st.playheadSec = recordingDurationSec(st);

      toast("Flight recorder stopped.", 1.4);
    }
  }
  ImGui::SameLine();
  if (ImGui::Button("Clear")) {
    st.samples.clear();
    st.sampleAccumulatorSec = 0.0;
    st.startRealSec = timeRealSec;
    st.startSimDays = simTimeDays;

    st.playing = false;
    st.playheadSec = 0.0;
    st.ghostSampleValid = false;

    if (st.recording) {
      pushSample(st, timeRealSec, simTimeDays, ship);
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
    const bool ok = exportCsv(st.csvPath, st.samples, &err);
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
  ImGui::Checkbox("Trace: pretty JSON", &st.tracePretty);

  if (ImGui::Button("Export trace")) {
    std::string err;
    const bool ok = exportTrace(st.tracePath, st, &err);
    toast(ok ? (std::string("Wrote trace: ") + st.tracePath) : (err.empty() ? "Trace export failed." : err), 2.4);
  }

  ImGui::End();
}

} // namespace stellar::game
