#include "CinematicCameraWindow.h"

#include <imgui.h>

#include <algorithm>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

namespace stellar::game {

static double clipDurationSec(const CinematicCameraWindowState& st) {
  double dur = 0.0;
  for (const auto& k : st.keys) dur = std::max(dur, k.tSec);
  return std::max(0.0, dur);
}

static void sortKeys(std::vector<CinematicCameraKeyframe>& keys) {
  for (auto& k : keys) k.tSec = std::max(0.0, k.tSec);
  std::sort(keys.begin(), keys.end(), [](const auto& a, const auto& b) { return a.tSec < b.tSec; });
}

static math::Vec3d lerp(const math::Vec3d& a, const math::Vec3d& b, double t) {
  return a * (1.0 - t) + b * t;
}

// Uniform Catmull-Rom spline (t in [0,1]).
static math::Vec3d catmullRom(const math::Vec3d& p0, const math::Vec3d& p1, const math::Vec3d& p2, const math::Vec3d& p3, double t) {
  const double t2 = t * t;
  const double t3 = t2 * t;
  return (p1 * 2.0 + (p2 - p0) * t + (p0 * 2.0 - p1 * 5.0 + p2 * 4.0 - p3) * t2 + (-p0 + p1 * 3.0 - p2 * 3.0 + p3) * t3) * 0.5;
}

static bool parseCsvLine(const std::string& line, double* out5) {
  // Very small parser: split by commas, trim whitespace, parse with stod.
  std::vector<std::string> cols;
  cols.reserve(5);
  {
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
      // trim
      auto begin = tok.begin();
      while (begin != tok.end() && std::isspace(static_cast<unsigned char>(*begin))) ++begin;
      auto end = tok.end();
      while (end != begin && std::isspace(static_cast<unsigned char>(*(end - 1)))) --end;
      cols.emplace_back(begin, end);
    }
  }
  if (cols.size() < 4) return false;

  try {
    out5[0] = std::stod(cols[0]);
    out5[1] = std::stod(cols[1]);
    out5[2] = std::stod(cols[2]);
    out5[3] = std::stod(cols[3]);
    out5[4] = (cols.size() >= 5) ? std::stod(cols[4]) : 60.0;
  } catch (...) {
    return false;
  }

  return true;
}

void tickCinematicCamera(CinematicCameraWindowState& st, double dtRealSec, double dtSimSec, bool paused) {
  if (!st.playing) return;

  if (paused && !st.playWhilePaused) return;

  double dt = st.playUsesSimTime ? dtSimSec : dtRealSec;
  dt *= st.playbackRate;
  if (dt == 0.0) return;

  st.playheadSec += dt;
  st.playheadSec = std::max(0.0, st.playheadSec);

  const double dur = clipDurationSec(st);
  if (dur <= 1e-9) return;

  if (st.loop) {
    st.playheadSec = std::fmod(st.playheadSec, dur);
    if (st.playheadSec < 0.0) st.playheadSec += dur;
  } else {
    if (st.playheadSec >= dur) {
      st.playheadSec = dur;
      st.playing = false;
    }
  }
}

CinematicCameraSample sampleCinematicCamera(const CinematicCameraWindowState& st) {
  CinematicCameraSample sample{};
  sample.offsetLocalU = st.editorOffsetU;
  sample.fovDeg = st.editorFovDeg;
  sample.hasFov = st.overrideFov;

  if (st.keys.empty()) return sample;

  const double dur = clipDurationSec(st);
  double t = st.playheadSec;
  if (dur > 1e-9 && st.loop) {
    t = std::fmod(t, dur);
    if (t < 0.0) t += dur;
  }

  // Clamp into key domain.
  t = std::max(0.0, t);
  if (t <= st.keys.front().tSec) {
    sample.offsetLocalU = st.keys.front().offsetLocalU;
    sample.fovDeg = st.keys.front().fovDeg;
    return sample;
  }
  if (t >= st.keys.back().tSec) {
    sample.offsetLocalU = st.keys.back().offsetLocalU;
    sample.fovDeg = st.keys.back().fovDeg;
    return sample;
  }

  // Find segment.
  int seg = 0;
  for (int i = 0; i + 1 < static_cast<int>(st.keys.size()); ++i) {
    if (t < st.keys[i + 1].tSec) {
      seg = i;
      break;
    }
  }

  const auto& k1 = st.keys[seg];
  const auto& k2 = st.keys[seg + 1];
  const double denom = std::max(1e-9, k2.tSec - k1.tSec);
  const double u = std::clamp((t - k1.tSec) / denom, 0.0, 1.0);

  // FOV: linear.
  sample.fovDeg = k1.fovDeg * (1.0 - u) + k2.fovDeg * u;

  // Offset: Catmull-Rom or linear.
  if (!st.catmullRom) {
    sample.offsetLocalU = lerp(k1.offsetLocalU, k2.offsetLocalU, u);
    return sample;
  }

  const math::Vec3d p0 = (seg - 1 >= 0) ? st.keys[seg - 1].offsetLocalU : k1.offsetLocalU;
  const math::Vec3d p1 = k1.offsetLocalU;
  const math::Vec3d p2 = k2.offsetLocalU;
  const math::Vec3d p3 = (seg + 2 < static_cast<int>(st.keys.size())) ? st.keys[seg + 2].offsetLocalU : k2.offsetLocalU;

  sample.offsetLocalU = catmullRom(p0, p1, p2, p3, u);
  return sample;
}

void drawCinematicCameraWindow(CinematicCameraWindowState& st, const ToastFn& toast) {
  if (!st.open) return;

  if (!ImGui::Begin("Cinematic Camera", &st.open)) {
    ImGui::End();
    return;
  }

  ImGui::TextUnformatted("Keyframed camera offsets for smooth fly-bys / chase-cam moves.");
  ImGui::Spacing();

  ImGui::Checkbox("Enable override", &st.enabled);
  ImGui::SameLine();
  ImGui::Checkbox("Look at ship", &st.lookAtShip);
  ImGui::SameLine();
  ImGui::Checkbox("Override FOV", &st.overrideFov);

  ImGui::Checkbox("Catmull-Rom (smooth)", &st.catmullRom);
  ImGui::SameLine();
  ImGui::Checkbox("Loop", &st.loop);

  ImGui::Checkbox("Advance using simulation dt", &st.playUsesSimTime);
  ImGui::SameLine();
  ImGui::Checkbox("Play while paused", &st.playWhilePaused);

  ImGui::PushItemWidth(160.0f);
  {
    float rate = static_cast<float>(st.playbackRate);
    if (ImGui::DragFloat("Playback rate", &rate, 0.01f, 0.0f, 10.0f, "%.2fx")) {
      st.playbackRate = std::max(0.0, static_cast<double>(rate));
    }
  }
  ImGui::PopItemWidth();

  const double dur = clipDurationSec(st);
  double tMin = 0.0;
  double tMax = std::max(0.001, dur);
  if (dur <= 1e-9) {
    ImGui::BeginDisabled();
  }
  ImGui::SliderScalar("Playhead (s)", ImGuiDataType_Double, &st.playheadSec, &tMin, &tMax, "%.3f");
  if (dur <= 1e-9) {
    ImGui::EndDisabled();
  }

  ImGui::SameLine();
  if (!st.playing) {
    if (ImGui::Button("Play")) {
      st.playing = true;
      if (st.playheadSec >= dur && !st.loop) st.playheadSec = 0.0;
    }
  } else {
    if (ImGui::Button("Pause")) st.playing = false;
  }
  ImGui::SameLine();
  if (ImGui::Button("Stop")) {
    st.playing = false;
    st.playheadSec = 0.0;
  }

  ImGui::Text("Keys: %d    Duration: %.3fs", static_cast<int>(st.keys.size()), dur);

  ImGui::Separator();

  if (ImGui::BeginTable("CineCamSplit", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV)) {
    ImGui::TableNextColumn();

    ImGui::TextUnformatted("Keyframes");
    ImGui::BeginChild("##cinecam_keys", ImVec2(0, 240), true);
    for (int i = 0; i < static_cast<int>(st.keys.size()); ++i) {
      char label[96];
      std::snprintf(label, sizeof(label), "%02d  t=%.3f  off=(%.2f,%.2f,%.2f)", i, st.keys[i].tSec,
                    st.keys[i].offsetLocalU.x, st.keys[i].offsetLocalU.y, st.keys[i].offsetLocalU.z);
      const bool selected = (st.selectedIndex == i);
      if (ImGui::Selectable(label, selected)) {
        st.selectedIndex = i;
        st.editorTimeSec = st.keys[i].tSec;
        st.editorOffsetU = st.keys[i].offsetLocalU;
        st.editorFovDeg = st.keys[i].fovDeg;
      }
    }
    ImGui::EndChild();

    if (ImGui::Button("Add @ playhead")) {
      st.editorTimeSec = st.playheadSec;
      st.keys.push_back(CinematicCameraKeyframe{st.editorTimeSec, st.editorOffsetU, st.editorFovDeg});
      sortKeys(st.keys);
      if (toast) toast("Added keyframe", 2.0);
    }
    ImGui::SameLine();
    if (ImGui::Button("Add @ editor time")) {
      st.keys.push_back(CinematicCameraKeyframe{st.editorTimeSec, st.editorOffsetU, st.editorFovDeg});
      sortKeys(st.keys);
      if (toast) toast("Added keyframe", 2.0);
    }

    ImGui::SameLine();
    if (ImGui::Button("Clear")) {
      st.keys.clear();
      st.selectedIndex = -1;
      st.playheadSec = 0.0;
      st.playing = false;
      if (toast) toast("Cleared keyframes", 2.0);
    }

    ImGui::TableNextColumn();

    ImGui::TextUnformatted("Editor");

    const bool hasSel = (st.selectedIndex >= 0 && st.selectedIndex < static_cast<int>(st.keys.size()));
    if (!hasSel) {
      ImGui::TextUnformatted("(No selection)\nUse the editor to set a preview offset and then add a keyframe.");
    }

    ImGui::PushItemWidth(220.0f);
    ImGui::DragScalar("Time (s)", ImGuiDataType_Double, &st.editorTimeSec, 0.05, &tMin, nullptr, "%.3f");

    float off[3] = {
      static_cast<float>(st.editorOffsetU.x),
      static_cast<float>(st.editorOffsetU.y),
      static_cast<float>(st.editorOffsetU.z)
    };
    if (ImGui::DragFloat3("Offset (R,U,F)", off, 0.05f, -200.0f, 200.0f, "%.2f")) {
      st.editorOffsetU = math::Vec3d{static_cast<double>(off[0]), static_cast<double>(off[1]), static_cast<double>(off[2])};
    }

    float fov = static_cast<float>(st.editorFovDeg);
    if (ImGui::SliderFloat("FOV (deg)", &fov, 15.0f, 120.0f, "%.1f")) {
      st.editorFovDeg = static_cast<double>(fov);
    }
    ImGui::PopItemWidth();

    if (ImGui::Button("Reset chase cam")) {
      st.editorOffsetU = math::Vec3d{0.0, 2.0, -6.0};
      st.editorFovDeg = 60.0;
    }
    ImGui::SameLine();
    if (ImGui::Button("Snap editor <- sample")) {
      const auto s = sampleCinematicCamera(st);
      st.editorOffsetU = s.offsetLocalU;
      st.editorFovDeg = s.fovDeg;
    }

    if (hasSel) {
      ImGui::Separator();
      if (ImGui::Button("Apply to selected")) {
        auto& k = st.keys[st.selectedIndex];
        k.tSec = st.editorTimeSec;
        k.offsetLocalU = st.editorOffsetU;
        k.fovDeg = st.editorFovDeg;
        sortKeys(st.keys);
        if (toast) toast("Updated keyframe", 2.0);
      }
      ImGui::SameLine();
      if (ImGui::Button("Delete selected")) {
        st.keys.erase(st.keys.begin() + st.selectedIndex);
        st.selectedIndex = -1;
        if (toast) toast("Deleted keyframe", 2.0);
      }
      ImGui::SameLine();
      if (ImGui::Button("Set playhead")) {
        st.playheadSec = st.editorTimeSec;
      }
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Import / Export");
    ImGui::InputText("CSV path", st.csvPath, sizeof(st.csvPath));

    if (ImGui::Button("Export CSV")) {
      std::ofstream f(st.csvPath);
      if (!f) {
        if (toast) toast("Failed to open CSV for writing", 3.0);
      } else {
        f << std::fixed << std::setprecision(6);
        f << "tSec,offX,offY,offZ,fovDeg\n";
        for (const auto& k : st.keys) {
          f << k.tSec << "," << k.offsetLocalU.x << "," << k.offsetLocalU.y << "," << k.offsetLocalU.z << "," << k.fovDeg << "\n";
        }
        if (toast) toast("Exported CSV", 2.0);
      }
    }

    ImGui::SameLine();
    if (ImGui::Button("Import CSV")) {
      std::ifstream f(st.csvPath);
      if (!f) {
        if (toast) toast("Failed to open CSV for reading", 3.0);
      } else {
        std::vector<CinematicCameraKeyframe> imported;
        std::string line;
        while (std::getline(f, line)) {
          // skip empty/comment
          bool onlyWs = true;
          for (char c : line) {
            if (!std::isspace(static_cast<unsigned char>(c))) { onlyWs = false; break; }
          }
          if (onlyWs) continue;
          if (!line.empty() && line[0] == '#') continue;

          double vals[5] = {0,0,0,0,60};
          if (!parseCsvLine(line, vals)) continue;
          CinematicCameraKeyframe k;
          k.tSec = std::max(0.0, vals[0]);
          k.offsetLocalU = math::Vec3d{vals[1], vals[2], vals[3]};
          k.fovDeg = vals[4];
          imported.push_back(k);
        }

        sortKeys(imported);
        st.keys = std::move(imported);
        st.selectedIndex = -1;
        st.playheadSec = 0.0;
        st.playing = false;
        if (toast) toast("Imported CSV", 2.0);
      }
    }

    ImGui::EndTable();
  }

  if (!st.enabled) {
    ImGui::Spacing();
    ImGui::TextUnformatted("Tip: enable the override to drive the in-game camera using this rig.");
  }

  ImGui::End();
}

} // namespace stellar::game
