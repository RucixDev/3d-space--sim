#include "SpectralMieLabWindow.h"

#include "Screenshot.h"

#include <imgui.h>

#include <algorithm>
#include <chrono>
#include <exception>
#include <cmath>

namespace stellar::game {

using Clock = std::chrono::high_resolution_clock;

static bool dragDouble(const char* label,
                       double& v,
                       double step,
                       double minV,
                       double maxV,
                       const char* fmt = "%.3f") {
  double tmp = v;
  const bool changed = ImGui::DragScalar(label, ImGuiDataType_Double, &tmp, step, &minV, &maxV, fmt);
  if (changed) v = tmp;
  return changed;
}

static SpectralMieLabJobOutput runJob(proc::SpectralMiePhaseSettings settings) {
  SpectralMieLabJobOutput out;

  const auto t0 = Clock::now();
  try {
    out.result = proc::generateSpectralMiePhase(settings);
    out.rgba = proc::spectralPhaseToRgba8(out.result);
  } catch (const std::exception& e) {
    out.error = std::string("Mie job exception: ") + e.what();
  } catch (...) {
    out.error = "Mie job exception (unknown)";
  }
  const auto t1 = Clock::now();
  out.ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
  return out;
}

static void rebuildMultipleScattering(SpectralMieLabWindowState& st, const ToastFn& toast, bool notify) {
  st.lastMsError.clear();
  st.lastMsGenMs = 0.0;

  if (!st.enableMultipleScattering || st.result.angleSamples <= 0 || st.result.mu.empty()) {
    st.msResult = {};
    st.msRgba.clear();
    return;
  }

  // Ensure the MS builder knows whether the single-scatter LUT was peak-normalized.
  st.msSettings.inputWasPeakNormalized = st.settings.peakNormalize;

  const auto t0 = Clock::now();
  try {
    st.msResult = proc::generateMultipleScatteringPhase(st.result, st.msSettings);
    st.msRgba = proc::spectralPhaseToRgba8(st.msResult);

    if (!st.msRgba.empty() && st.msResult.angleSamples > 0) {
      st.msLutTex.createRGBA(st.msResult.angleSamples, 1, st.msRgba.data(),
                             /*generateMips=*/false,
                             /*nearestFilter=*/false,
                             /*clampToEdge=*/true);
    }
  } catch (const std::exception& e) {
    st.lastMsError = std::string("MS phase exception: ") + e.what();
  } catch (...) {
    st.lastMsError = "MS phase exception (unknown)";
  }
  const auto t1 = Clock::now();
  st.lastMsGenMs = std::chrono::duration<double, std::milli>(t1 - t0).count();

  if (notify && st.lastMsError.empty() && st.msResult.angleSamples > 0) {
    toast("Analytic multiple-scatter LUT updated", 1.25);
  }
}

static void requestGenerate(SpectralMieLabWindowState& st) {
  st.lastError.clear();

  if (st.jobRunning) {
    st.jobQueued = true;
    st.queuedSettings = st.settings;
    return;
  }

  st.dirty = false;
  st.jobRunning = true;
  const proc::SpectralMiePhaseSettings settingsCopy = st.settings;
  st.jobFuture = std::async(std::launch::async, [settingsCopy]() {
    return runJob(settingsCopy);
  });
}

static void pollJob(SpectralMieLabWindowState& st, const ToastFn& toast) {
  if (!st.jobRunning) return;
  if (!st.jobFuture.valid()) {
    st.jobRunning = false;
    return;
  }

  using namespace std::chrono_literals;
  if (st.jobFuture.wait_for(0ms) != std::future_status::ready) return;

  SpectralMieLabJobOutput res = st.jobFuture.get();
  st.jobRunning = false;

  st.lastGenMs = res.ms;
  st.lastError = res.error;

  if (res.error.empty()) {
    st.result = std::move(res.result);
    st.rgba = std::move(res.rgba);

    if (!st.rgba.empty() && st.result.angleSamples > 0) {
      // Upload (height=1). Use linear filtering + clamp so sampling is smooth.
      st.lutTex.createRGBA(st.result.angleSamples, 1, st.rgba.data(),
                           /*generateMips=*/false,
                           /*nearestFilter=*/false,
                           /*clampToEdge=*/true);
    }

    toast("Spectral Mie LUT generated (" + std::to_string((int)st.result.angleSamples) + " samples)", 1.25);

    // Rebuild derived multiple-scattering LUT (fast, synchronous).
    rebuildMultipleScattering(st, toast, /*notify=*/false);
  }

  // If new settings were queued while we were computing, run again.
  if (st.jobQueued) {
    st.jobQueued = false;
    st.settings = st.queuedSettings;
    requestGenerate(st);
  }
}

static const char* distLabel(proc::MieRadiusDistribution d) {
  switch (d) {
    case proc::MieRadiusDistribution::Mono: return "Mono (single radius)";
    case proc::MieRadiusDistribution::LogNormal: return "Log-normal (approx)";
    default: return "(unknown)";
  }
}

void drawSpectralMieLabWindow(SpectralMieLabWindowState& st,
                             float /*timeSec*/,
                             const ToastFn& toast) {
  if (!st.open) return;

  pollJob(st, toast);

  ImGui::SetNextWindowSize(ImVec2(1140, 720), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Spectral Mie Lab", &st.open)) {
    ImGui::End();
    return;
  }

  bool dirty = false;
  bool dirtyMs = false;

  if (ImGui::BeginTable("##mie_layout", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV)) {
    ImGui::TableNextColumn();

    // ----- Controls -----
    {
      ImGui::TextUnformatted("Generate a spectral Mie phase function LUT p(μ) where μ = cos(θ).");
      ImGui::TextDisabled("Left = μ=-1 (backscatter), Right = μ=+1 (forward scatter).\nGenerate peak-normalized curves for real-time shading.");

      ImGui::SeparatorText("Generation");

      dirty |= ImGui::Checkbox("Auto regenerate", &st.autoRegenerate);
      ImGui::SameLine();
      if (ImGui::Button(st.jobRunning ? "Generating..." : "Generate")) {
        requestGenerate(st);
      }

      if (st.jobRunning) {
        ImGui::SameLine();
        ImGui::TextDisabled("(running)");
      }

      int samples = st.settings.angleSamples;
      if (ImGui::InputInt("LUT samples", &samples)) {
        st.settings.angleSamples = std::clamp(samples, 8, 4096);
        dirty = true;
      }

      dirty |= ImGui::Checkbox("Peak normalize (max=1)", &st.settings.peakNormalize);

      ImGui::SeparatorText("Spectrum (nm)");
      dirty |= dragDouble("λ Red",   st.settings.wavelengthsNm[0], 1.0, 300.0, 900.0, "%.0f");
      dirty |= dragDouble("λ Green", st.settings.wavelengthsNm[1], 1.0, 300.0, 900.0, "%.0f");
      dirty |= dragDouble("λ Blue",  st.settings.wavelengthsNm[2], 1.0, 300.0, 900.0, "%.0f");

      ImGui::SeparatorText("Particle");
      dirty |= dragDouble("Radius (µm)", st.settings.radiusUm, 0.005, 0.001, 50.0, "%.4f");

      // Distribution selector.
      {
        int dist = (int)st.settings.radiusDist;
        const char* preview = distLabel(st.settings.radiusDist);
        if (ImGui::BeginCombo("Radius distribution", preview)) {
          for (int i = 0; i < 2; ++i) {
            const auto d = (proc::MieRadiusDistribution)i;
            const bool selected = (i == dist);
            if (ImGui::Selectable(distLabel(d), selected)) {
              dist = i;
              st.settings.radiusDist = d;
              dirty = true;
            }
            if (selected) ImGui::SetItemDefaultFocus();
          }
          ImGui::EndCombo();
        }
      }

      if (st.settings.radiusDist == proc::MieRadiusDistribution::LogNormal) {
        dirty |= dragDouble("Log-normal σ", st.settings.lognormalSigma, 0.01, 0.01, 1.50, "%.3f");
        int k = st.settings.radiusSampleCount;
        if (ImGui::InputInt("Radius samples", &k)) {
          st.settings.radiusSampleCount = std::clamp(k, 2, 64);
          dirty = true;
        }
      }

      ImGui::SeparatorText("Refractive index");
      dirty |= dragDouble("n", st.settings.refractiveIndexN, 0.005, 1.00, 2.50, "%.4f");
      dirty |= dragDouble("k (absorption)", st.settings.refractiveIndexK, 0.001, 0.00, 1.00, "%.4f");

      ImGui::SeparatorText("Multiple scattering (analytic)");
      dirtyMs |= ImGui::Checkbox("Enable analytic MS LUT", &st.enableMultipleScattering);
      ImGui::BeginDisabled(!st.enableMultipleScattering);
      {
        int L = st.msSettings.legendreOrder;
        if (ImGui::SliderInt("Legendre order", &L, 4, 128, "%d")) {
          st.msSettings.legendreOrder = std::clamp(L, 0, 256);
          dirtyMs = true;
        }

        float w = (float)st.msSettings.scatteringAlbedo;
        if (ImGui::SliderFloat("Albedo ω", &w, 0.0f, 0.99f, "%.3f")) {
          st.msSettings.scatteringAlbedo = (double)std::clamp(w, 0.0f, 0.999f);
          dirtyMs = true;
        }

        // Mode selector.
        {
          int mode = (int)st.msSettings.mode;
          const char* preview = (mode == (int)proc::MultipleScatteringMode::TotalOrders)
                                  ? "Total (orders 1+)"
                                  : "Multiple only (orders 2+)";
          if (ImGui::BeginCombo("MS mode", preview)) {
            for (int i = 0; i < 2; ++i) {
              const auto m = (proc::MultipleScatteringMode)i;
              const bool selected = (i == mode);
              const char* label = (m == proc::MultipleScatteringMode::TotalOrders)
                                    ? "Total (orders 1+)"
                                    : "Multiple only (orders 2+)";
              if (ImGui::Selectable(label, selected)) {
                mode = i;
                st.msSettings.mode = m;
                dirtyMs = true;
              }
              if (selected) ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
          }
        }

        dirtyMs |= ImGui::Checkbox("Peak normalize MS (max=1)", &st.msSettings.peakNormalize);

        if (ImGui::Button("Rebuild MS LUT")) {
          dirtyMs = true;
        }

        ImGui::TextDisabled("MS LUT is derived from the current single-scatter LUT\n"
                            "using Legendre moments + a closed-form series.");
      }
      ImGui::EndDisabled();

      ImGui::SeparatorText("Runtime hookup");
      dirty |= ImGui::Checkbox("Apply to atmosphere renderer", &st.applyToAtmospheres);
      ImGui::BeginDisabled(!st.applyToAtmospheres);
      dirty |= ImGui::SliderFloat("Mie strength", &st.atmosphereMieStrength, 0.0f, 1.0f, "%.2f");
      ImGui::EndDisabled();

      dirtyMs |= ImGui::Checkbox("Apply MS LUT to volumetric", &st.applyMsToVolumetric);
      ImGui::BeginDisabled(!st.applyMsToVolumetric);
      dirtyMs |= ImGui::SliderFloat("MS phase strength", &st.atmosphereMsPhaseStrength, 0.0f, 1.0f, "%.2f");
      ImGui::EndDisabled();

      ImGui::SeparatorText("Stats");
      if (!st.lastError.empty()) {
        ImGui::TextColored(ImVec4(1, 0.35f, 0.35f, 1), "%s", st.lastError.c_str());
      }

      if (!st.lastMsError.empty()) {
        ImGui::TextColored(ImVec4(1, 0.55f, 0.35f, 1), "%s", st.lastMsError.c_str());
      }

      if (st.lastGenMs > 0.0) {
        ImGui::Text("Last gen: %.2f ms", st.lastGenMs);
      }

      if (st.lastMsGenMs > 0.0) {
        ImGui::Text("MS rebuild: %.2f ms", st.lastMsGenMs);
      }

      if (st.result.angleSamples > 0) {
        ImGui::Text("g (R,G,B): %.3f  %.3f  %.3f",
                    st.result.asymmetryG[0], st.result.asymmetryG[1], st.result.asymmetryG[2]);
        ImGui::Text("Qsca (R,G,B): %.3f  %.3f  %.3f",
                    st.result.qSca[0], st.result.qSca[1], st.result.qSca[2]);
        ImGui::Text("pMax (R,G,B): %.3g  %.3g  %.3g",
                    st.result.pMax[0], st.result.pMax[1], st.result.pMax[2]);
      }

      if (st.enableMultipleScattering && st.msResult.angleSamples > 0) {
        ImGui::SeparatorText("MS stats");
        ImGui::Text("g_ms (R,G,B): %.3f  %.3f  %.3f",
                    st.msResult.asymmetryG[0], st.msResult.asymmetryG[1], st.msResult.asymmetryG[2]);
        ImGui::Text("pMax_ms (R,G,B): %.3g  %.3g  %.3g",
                    st.msResult.pMax[0], st.msResult.pMax[1], st.msResult.pMax[2]);
      }

      ImGui::SeparatorText("Export");
      ImGui::InputText("Dir", st.exportDir, sizeof(st.exportDir));
      ImGui::InputText("Base name", st.exportBaseName, sizeof(st.exportBaseName));
      ImGui::Checkbox("Timestamp", &st.exportTimestamp);

      if (ImGui::Button("Export LUT as PNG")) {
        if (st.rgba.empty() || st.result.angleSamples <= 0) {
          toast("No LUT to export (generate first)", 1.6);
        } else {
          ScreenshotRequest req;
          req.includeUi = false;
          req.outDir = st.exportDir;
          req.baseName = st.exportBaseName;
          req.timestamp = st.exportTimestamp;
          req.extension = "png";

          std::string err;
          const std::string path = buildScreenshotPath(req, &err);
          if (path.empty()) {
            toast("Export failed: " + err, 2.5);
          } else {
            std::string werr;
            const bool ok = writePixelsToPng(path,
                                            st.result.angleSamples,
                                            1,
                                            4,
                                            st.rgba.data(),
                                            st.result.angleSamples * 4,
                                            /*flipY=*/false,
                                            &werr);
            if (ok) toast("Exported: " + path, 3.0);
            else toast("Export failed: " + werr, 3.0);
          }
        }
      }

      ImGui::SameLine();
      if (ImGui::Button("Export MS LUT")) {
        if (!st.enableMultipleScattering || st.msRgba.empty() || st.msResult.angleSamples <= 0) {
          toast("No MS LUT to export (generate first)", 1.6);
        } else {
          ScreenshotRequest req;
          req.includeUi = false;
          req.outDir = st.exportDir;
          req.baseName = std::string(st.exportBaseName) + "_ms";
          req.timestamp = st.exportTimestamp;
          req.extension = "png";

          std::string err;
          const std::string path = buildScreenshotPath(req, &err);
          if (path.empty()) {
            toast("Export failed: " + err, 2.5);
          } else {
            std::string werr;
            const bool ok = writePixelsToPng(path,
                                             st.msResult.angleSamples,
                                             1,
                                             4,
                                             st.msRgba.data(),
                                             st.msResult.angleSamples * 4,
                                             /*flipY=*/false,
                                             &werr);
            if (ok) toast("Exported: " + path, 3.0);
            else toast("Export failed: " + werr, 3.0);
          }
        }
      }
    }

    ImGui::TableNextColumn();

    // ----- Preview -----
    {
      ImGui::TextUnformatted("Preview");
      if (st.lutTex.handle() != 0) {
        // Draw the 1px-high texture stretched to a visible strip.
        ImGui::Image((ImTextureID)(intptr_t)st.lutTex.handle(), ImVec2(640, 56), ImVec2(0, 0), ImVec2(1, 1));
      } else {
        ImGui::TextDisabled("(no LUT yet)");
      }

      if (st.enableMultipleScattering) {
        ImGui::SeparatorText("Multiple scattering LUT");
        if (st.msLutTex.handle() != 0) {
          ImGui::Image((ImTextureID)(intptr_t)st.msLutTex.handle(), ImVec2(640, 56), ImVec2(0, 0), ImVec2(1, 1));
        } else {
          ImGui::TextDisabled("(no MS LUT yet)");
        }
      }

      ImGui::Spacing();
      ImGui::SeparatorText("Phase curves");

      if (!st.result.phaseR.empty()) {
        ImGui::PlotLines("Red", st.result.phaseR.data(), (int)st.result.phaseR.size(), 0, nullptr, 0.0f, 1.0f, ImVec2(0, 80));
      }
      if (!st.result.phaseG.empty()) {
        ImGui::PlotLines("Green", st.result.phaseG.data(), (int)st.result.phaseG.size(), 0, nullptr, 0.0f, 1.0f, ImVec2(0, 80));
      }
      if (!st.result.phaseB.empty()) {
        ImGui::PlotLines("Blue", st.result.phaseB.data(), (int)st.result.phaseB.size(), 0, nullptr, 0.0f, 1.0f, ImVec2(0, 80));
      }

      if (st.enableMultipleScattering && !st.msResult.phaseR.empty()) {
        ImGui::SeparatorText("MS phase curves");
        ImGui::PlotLines("MS Red", st.msResult.phaseR.data(), (int)st.msResult.phaseR.size(), 0, nullptr, 0.0f, 1.0f, ImVec2(0, 80));
        ImGui::PlotLines("MS Green", st.msResult.phaseG.data(), (int)st.msResult.phaseG.size(), 0, nullptr, 0.0f, 1.0f, ImVec2(0, 80));
        ImGui::PlotLines("MS Blue", st.msResult.phaseB.data(), (int)st.msResult.phaseB.size(), 0, nullptr, 0.0f, 1.0f, ImVec2(0, 80));
      }

      ImGui::TextDisabled("Note: PlotLines is clamped to [0,1] for readability. Use physical normalization for pMax/g stats.");
    }

    ImGui::EndTable();
  }

  // Trigger regeneration if needed.
  if (dirty) st.dirty = true;
  if (st.autoRegenerate && st.dirty && !st.jobRunning) {
    requestGenerate(st);
  }

  // MS LUT rebuilds are cheap; do them synchronously when params change.
  if (dirtyMs && !st.jobRunning) {
    rebuildMultipleScattering(st, toast, /*notify=*/true);
  }

  ImGui::End();
}

} // namespace stellar::game
