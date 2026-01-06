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
#include <string>
#include <vector>

namespace stellar::game {

static double speedKmS(const FlightRecorderSample& s) {
  return s.velKmS.length();
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

  f << "t_real_sec,t_sim_days,pos_x_km,pos_y_km,pos_z_km,vel_x_kms,vel_y_kms,vel_z_kms,speed_kms\n";
  for (const auto& s : samples) {
    f << s.tRealSec << ','
      << s.tSimDays << ','
      << s.posKm.x << ',' << s.posKm.y << ',' << s.posKm.z << ','
      << s.velKmS.x << ',' << s.velKmS.y << ',' << s.velKmS.z << ','
      << speedKmS(s)
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
  keys.reserve(16);

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

    ++i;
  }

  core::ChromeTraceWriteOptions opt;
  opt.includeFrameEvents = false;
  opt.pretty = st.tracePretty;
  opt.pid = 1;
  opt.tid = 1;

  return core::writeCounterChromeTraceJson(path, table, opt, err);
}

void tickFlightRecorder(FlightRecorderWindowState& st,
                        double dtRealSec,
                        double timeRealSec,
                        double simTimeDays,
                        const sim::Ship& ship,
                        bool paused) {
  if (!st.recording) return;
  if (paused && !st.recordWhilePaused) return;

  // Guard against invalid sample rates.
  if (st.sampleHz <= 0.0) st.sampleHz = 20.0;

  const double intervalSec = 1.0 / st.sampleHz;
  st.sampleAccumulatorSec += dtRealSec;

  // We only sample at most once per frame (we don't attempt to reconstruct
  // intermediate ship states).
  if (st.sampleAccumulatorSec < intervalSec) return;
  st.sampleAccumulatorSec = 0.0;

  pushSample(st, timeRealSec, simTimeDays, ship);
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

  ImGui::Separator();

  // ---- Controls ----
  if (!st.recording) {
    if (ImGui::Button("Start recording")) {
      st.recording = true;
      st.sampleAccumulatorSec = 0.0;
      st.startRealSec = timeRealSec;
      st.startSimDays = simTimeDays;
      st.samples.clear();
      // Capture an initial sample immediately.
      pushSample(st, timeRealSec, simTimeDays, ship);
      toast("Flight recorder started.", 1.4);
    }
  } else {
    if (ImGui::Button("Stop recording")) {
      st.recording = false;
      toast("Flight recorder stopped.", 1.4);
    }
  }
  ImGui::SameLine();
  if (ImGui::Button("Clear")) {
    st.samples.clear();
    st.sampleAccumulatorSec = 0.0;
    st.startRealSec = timeRealSec;
    st.startSimDays = simTimeDays;
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
  ImGui::Checkbox("Trace: pretty JSON", &st.tracePretty);

  if (ImGui::Button("Export trace")) {
    std::string err;
    const bool ok = exportTrace(st.tracePath, st, &err);
    toast(ok ? (std::string("Wrote trace: ") + st.tracePath) : (err.empty() ? "Trace export failed." : err), 2.4);
  }

  ImGui::End();
}

} // namespace stellar::game
