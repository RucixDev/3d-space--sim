#include "AudioAnalyzerWindow.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

namespace stellar::game {
namespace {

static void writeLe16(std::ofstream& f, std::uint16_t v) {
  char b[2];
  b[0] = (char)(v & 0xFF);
  b[1] = (char)((v >> 8) & 0xFF);
  f.write(b, 2);
}

static void writeLe32(std::ofstream& f, std::uint32_t v) {
  char b[4];
  b[0] = (char)(v & 0xFF);
  b[1] = (char)((v >> 8) & 0xFF);
  b[2] = (char)((v >> 16) & 0xFF);
  b[3] = (char)((v >> 24) & 0xFF);
  f.write(b, 4);
}

static bool writeWavMono16(const char* path, const float* samples, int count, int sampleRate) {
  if (!path || !samples || count <= 0 || sampleRate <= 0) return false;

  try {
    const std::filesystem::path p(path);
    if (p.has_parent_path()) {
      std::filesystem::create_directories(p.parent_path());
    }
  } catch (...) {
    // If directory creation fails, we still attempt to write and let the stream
    // decide.
  }

  std::ofstream f(path, std::ios::binary);
  if (!f.is_open()) return false;

  const std::uint16_t channels = 1;
  const std::uint16_t bitsPerSample = 16;
  const std::uint32_t byteRate = (std::uint32_t)sampleRate * channels * (bitsPerSample / 8);
  const std::uint16_t blockAlign = channels * (bitsPerSample / 8);

  const std::uint32_t dataBytes = (std::uint32_t)count * channels * (bitsPerSample / 8);
  const std::uint32_t riffSize = 36u + dataBytes;

  // RIFF header
  f.write("RIFF", 4);
  writeLe32(f, riffSize);
  f.write("WAVE", 4);

  // fmt chunk
  f.write("fmt ", 4);
  writeLe32(f, 16);
  writeLe16(f, 1); // PCM
  writeLe16(f, channels);
  writeLe32(f, (std::uint32_t)sampleRate);
  writeLe32(f, byteRate);
  writeLe16(f, blockAlign);
  writeLe16(f, bitsPerSample);

  // data chunk
  f.write("data", 4);
  writeLe32(f, dataBytes);

  // Samples
  for (int i = 0; i < count; ++i) {
    const float x = std::clamp(samples[i], -1.0f, 1.0f);
    int v = (int)std::lround(x * 32767.0f);
    v = std::clamp(v, -32767, 32767);
    writeLe16(f, (std::uint16_t)(std::int16_t)v);
  }

  return true;
}

static int sanitizePow2(int v) {
  v = std::clamp(v, 64, 16384);
  // Round to nearest power-of-two choice if a user typed something custom.
  int p = 64;
  while (p < v) p <<= 1;
  return p;
}

} // namespace

void drawAudioAnalyzerWindow(AudioAnalyzerWindowState& st,
                             AudioEngine& audio,
                             const ToastFn& toast) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(820, 640), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Audio Analyzer", &st.open)) {
    ImGui::End();
    return;
  }

  const bool active = audio.active();
  const float sr = audio.sampleRate();

  ImGui::Text("Device: %s", active ? "Active" : "Inactive");
  ImGui::SameLine();
  ImGui::TextDisabled("(sampleRate: %.0f Hz)", sr);

  ImGui::Separator();

  // ---- Controls ----
  ImGui::Checkbox("Pause", &st.pause);
  ImGui::SameLine();
  ImGui::Checkbox("Waveform", &st.showWaveform);
  ImGui::SameLine();
  ImGui::Checkbox("Spectrum", &st.showSpectrum);

  // FFT size combo
  {
    const int sizes[] = {256, 512, 1024, 2048, 4096, 8192, 16384};
    int current = st.params.fftSize;
    int curIdx = 0;
    for (int i = 0; i < (int)std::size(sizes); ++i) {
      if (sizes[i] == current) curIdx = i;
    }

    ImGui::SetNextItemWidth(140.0f);
    if (ImGui::BeginCombo("FFT Size", std::to_string(current).c_str())) {
      for (int i = 0; i < (int)std::size(sizes); ++i) {
        const bool sel = (i == curIdx);
        if (ImGui::Selectable(std::to_string(sizes[i]).c_str(), sel)) {
          st.params.fftSize = sizes[i];
        }
        if (sel) ImGui::SetItemDefaultFocus();
      }
      ImGui::EndCombo();
    }
  }

  ImGui::SameLine();
  ImGui::SetNextItemWidth(120.0f);
  ImGui::DragInt("Bands", &st.params.bands, 1.0f, 8, 256);

  ImGui::SetNextItemWidth(200.0f);
  ImGui::SliderFloat("Smoothing", &st.params.smoothing, 0.0f, 0.95f, "%.2f");
  ImGui::SameLine();
  ImGui::SetNextItemWidth(200.0f);
  ImGui::SliderFloat("Floor dB", &st.params.floorDb, -120.0f, -20.0f, "%.0f dB");

  ImGui::SetNextItemWidth(200.0f);
  ImGui::DragInt("Capture Samples", &st.captureSamples, 256.0f, 2048, 32768);

  st.params.fftSize = sanitizePow2(st.params.fftSize);
  st.captureSamples = std::clamp(st.captureSamples, 512, 32768);

  // Apply parameter changes.
  st.analyzer.setParams(st.params);

  // Beat tracker params are independent from the FFT/banding params.
  st.beat.setParams(st.beatParams);

  // Pull capture + analyze.
  if (!st.pause) {
    st.capture.resize((std::size_t)st.captureSamples);
    const int n = audio.copyRecentMono(st.capture.data(), st.captureSamples);
    st.capture.resize((std::size_t)std::max(0, n));

    if (!st.capture.empty()) {
      st.analyzer.analyzeWindow(st.capture.data(), (int)st.capture.size(), sr);
    } else {
      st.analyzer.analyzeWindow(nullptr, 0, sr);
    }
  }

  // Reset beat tracker when toggled so it doesn't "remember" stale spectra.
  if (st.beatEnabled != st.beatPrevEnabled) {
    st.beat.reset();
    st.onsetHistory.clear();
    st.thresholdHistory.clear();
    st.beatPrevEnabled = st.beatEnabled;
  }

  // Beat/onset analysis runs on the most recent spectrum (even when paused).
  // This makes the beat pulse decay feel responsive in the UI.
  const float dtSec = ImGui::GetIO().DeltaTime;
  const auto& spec = st.analyzer.spectrumMag();
  const dsp::BeatTrackerOutput beatOut = (st.beatEnabled && !spec.empty())
    ? st.beat.processSpectrum(spec.data(), (int)spec.size(), dtSec)
    : st.beat.processSpectrum(nullptr, 0, dtSec);

  // Maintain a small onset history for plotting.
  st.onsetHistoryMax = std::clamp(st.onsetHistoryMax, 60, 4096);
  if (st.beatEnabled) {
    st.onsetHistory.push_back(beatOut.onset);
    st.thresholdHistory.push_back(beatOut.threshold);
    while ((int)st.onsetHistory.size() > st.onsetHistoryMax) st.onsetHistory.erase(st.onsetHistory.begin());
    while ((int)st.thresholdHistory.size() > st.onsetHistoryMax) st.thresholdHistory.erase(st.thresholdHistory.begin());
  }

  ImGui::Separator();

  // ---- Stats ----
  ImGui::Text("RMS: %.4f", st.analyzer.rms());
  ImGui::SameLine();
  ImGui::Text("Peak: %.3f", st.analyzer.peak());

  if (st.beatEnabled) {
    ImGui::SameLine();
    ImGui::Text("| BPM: %s%.1f", (beatOut.bpmConfidence01 >= 0.35f) ? "" : "(~)", beatOut.bpm);

    ImGui::SameLine();
    if (beatOut.beat) {
      ImGui::TextColored(ImVec4(0.95f, 0.90f, 0.20f, 1.0f), "BEAT");
    } else {
      ImGui::TextDisabled("beatPulse=%.2f", beatOut.beatPulse01);
    }
  }

  // ---- Plots ----
  if (st.showWaveform) {
    const auto& w = st.analyzer.waveform();
    if (!w.empty()) {
      ImGui::PlotLines("Waveform", w.data(), (int)w.size(), 0, nullptr, -1.0f, 1.0f, ImVec2(0, 140));
    }
  }

  if (st.showSpectrum) {
    ImGui::Checkbox("Spectrum in dB", &st.spectrumUseDb);

    const auto& bands = st.spectrumUseDb ? st.analyzer.bandsDb01() : st.analyzer.bandsMag();
    if (!bands.empty()) {
      float yMax = 1.0f;
      float yMin = 0.0f;
      if (!st.spectrumUseDb) {
        yMax = *std::max_element(bands.begin(), bands.end());
        yMax = std::max(1.0e-5f, yMax);
      }
      ImGui::PlotHistogram("Spectrum", bands.data(), (int)bands.size(), 0, nullptr, yMin, yMax, ImVec2(0, 160));
    }

    ImGui::TextDisabled("Tip: A Hann window + FFT reduces leakage vs raw DFT; smoothing stabilizes UI.");
  }

  // ---- Beat / Onset ----
  if (ImGui::CollapsingHeader("Beat / Onset", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Checkbox("Enable", &st.beatEnabled);
    ImGui::SameLine();
    ImGui::Checkbox("Plot", &st.beatShowPlot);
    ImGui::SameLine();
    ImGui::Checkbox("Debug", &st.beatShowDebug);

    ImGui::SetNextItemWidth(200.0f);
    ImGui::SliderFloat("Sensitivity", &st.beatParams.sensitivity, 0.0f, 4.0f, "%.2f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(200.0f);
    ImGui::SliderFloat("Min beat interval", &st.beatParams.minBeatIntervalSec, 0.05f, 1.0f, "%.2f s");

    ImGui::SetNextItemWidth(200.0f);
    ImGui::SliderFloat("Mean tau", &st.beatParams.meanTauSec, 0.02f, 2.0f, "%.2f s");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(200.0f);
    ImGui::SliderFloat("Dev tau", &st.beatParams.devTauSec, 0.02f, 2.0f, "%.2f s");

    ImGui::SetNextItemWidth(200.0f);
    ImGui::SliderFloat("Pulse decay", &st.beatParams.beatPulseDecaySec, 0.05f, 1.5f, "%.2f s");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::DragInt("History", &st.beatParams.bpmHistory, 1.0f, 4, 48);

    ImGui::SetNextItemWidth(120.0f);
    ImGui::DragFloat("BPM min", &st.beatParams.bpmMin, 1.0f, 20.0f, 240.0f, "%.0f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::DragFloat("BPM max", &st.beatParams.bpmMax, 1.0f, st.beatParams.bpmMin + 1.0f, 420.0f, "%.0f");

    ImGui::SetNextItemWidth(120.0f);
    ImGui::DragInt("Plot samples", &st.onsetHistoryMax, 2.0f, 64, 2048);

    if (!st.beatEnabled) {
      ImGui::TextDisabled("Beat tracking is disabled.");
    } else {
      ImGui::Text("Onset: %.4f | Thr: %.4f | mean: %.4f | dev: %.4f",
                  beatOut.onset, beatOut.threshold, beatOut.mean, beatOut.deviation);
      ImGui::Text("BPM: %.1f (conf %.2f)", beatOut.bpm, beatOut.bpmConfidence01);
    }

    if (st.beatEnabled && st.beatShowPlot && !st.onsetHistory.empty() && st.onsetHistory.size() == st.thresholdHistory.size()) {
      // Plot onset and threshold as two stacked plots so you can tune sensitivity.
      const float onsetMax = *std::max_element(st.onsetHistory.begin(), st.onsetHistory.end());
      const float thrMax = *std::max_element(st.thresholdHistory.begin(), st.thresholdHistory.end());
      const float yMax = std::max(1.0e-5f, std::max(onsetMax, thrMax) * 1.10f);
      ImGui::PlotLines("Onset", st.onsetHistory.data(), (int)st.onsetHistory.size(), 0, nullptr, 0.0f, yMax, ImVec2(0, 90));
      ImGui::PlotLines("Threshold", st.thresholdHistory.data(), (int)st.thresholdHistory.size(), 0, nullptr, 0.0f, yMax, ImVec2(0, 70));
    }

    if (st.beatEnabled && st.beatShowDebug) {
      ImGui::TextDisabled("Notes: This is a tiny spectral-flux onset tracker with adaptive threshold and peak-picking.");
      ImGui::TextDisabled("If you see double-triggers, increase Min beat interval; if you see misses, lower Sensitivity.");
    }
  }

  ImGui::Separator();

  // ---- SFX preview ----
  if (ImGui::CollapsingHeader("SFX Preview", ImGuiTreeNodeFlags_DefaultOpen)) {
    const struct { AudioEngine::Sfx sfx; const char* label; } kButtons[] = {
        {AudioEngine::Sfx::UiClick, "UI Click"},
        {AudioEngine::Sfx::UiConfirm, "UI Confirm"},
        {AudioEngine::Sfx::UiError, "UI Error"},
        {AudioEngine::Sfx::WeaponLaser, "Laser"},
        {AudioEngine::Sfx::WeaponPulse, "Pulse"},
        {AudioEngine::Sfx::WeaponCannon, "Cannon"},
        {AudioEngine::Sfx::WeaponRailgun, "Railgun"},
        {AudioEngine::Sfx::WeaponMining, "Mining"},
        {AudioEngine::Sfx::WeaponMissile, "Missile"},
        {AudioEngine::Sfx::Explosion, "Explosion"},
        {AudioEngine::Sfx::FsdCharge, "FSD Charge"},
        {AudioEngine::Sfx::FsdJump, "FSD Jump"},
        {AudioEngine::Sfx::FsdArrive, "FSD Arrive"},
    };

    if (ImGui::BeginTable("##sfx", 4, ImGuiTableFlags_SizingStretchProp)) {
      for (int i = 0; i < (int)std::size(kButtons); ++i) {
        ImGui::TableNextColumn();
        if (ImGui::Button(kButtons[i].label)) {
          audio.play(kButtons[i].sfx, 1.0f, 0.0f);
        }
      }
      ImGui::EndTable();
    }
  }

  // ---- Export ----
  if (ImGui::CollapsingHeader("Export", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::InputText("Path", st.exportPath, (int)sizeof(st.exportPath));

    ImGui::BeginDisabled(st.capture.empty());
    if (ImGui::Button("Export WAV (mono)")) {
      const bool ok = writeWavMono16(st.exportPath,
                                    st.capture.data(),
                                    (int)st.capture.size(),
                                    (int)std::lround(sr));
      if (ok) {
        toast(std::string("Exported: ") + st.exportPath, 2.5);
      } else {
        toast(std::string("Failed to export: ") + st.exportPath, 2.5);
      }
    }
    ImGui::EndDisabled();

    ImGui::TextDisabled("Exports the most recent capture window (mono 16-bit PCM).\n"
                        "If audio is disabled, the capture will export silence.");
  }

  ImGui::End();
}

} // namespace stellar::game
