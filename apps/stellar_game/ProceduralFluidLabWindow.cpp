#include "ProceduralFluidLabWindow.h"

#include "Screenshot.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <string>

namespace stellar::game {

namespace fs = std::filesystem;

namespace {

void ensureSimSize(ProceduralFluidLabWindowState& st) {
  const int n = std::clamp(st.gridSize, 32, 512);
  st.gridSize = n;
  if (st.sim.gridSize() != n) {
    st.sim.resize(n);
  }
}

void ensureTexture(ProceduralFluidLabWindowState& st) {
  const int n = st.sim.gridSize();
  if (n <= 0) return;

  if (st.tex.handle() == 0 || st.tex.width() != n || st.tex.height() != n) {
    st.tex.allocateRGBA(n, n, /*generateMips=*/false, /*nearestFilter=*/false, /*clampToEdge=*/true);
    st.rgba.assign((std::size_t)n * (std::size_t)n * 4u, 0);
    st.tex.updateRGBA(0, 0, n, n, st.rgba.data());
  }
}

void updateTextureFromSim(ProceduralFluidLabWindowState& st) {
  const int n = st.sim.gridSize();
  if (n <= 0 || st.tex.handle() == 0) return;

  st.rgba.resize((std::size_t)n * (std::size_t)n * 4u);

  const int stride = st.sim.paddedStride();
  const auto& rF = st.sim.dyeR();
  const auto& gF = st.sim.dyeG();
  const auto& bF = st.sim.dyeB();

  const float exposure = std::max(0.0001f, st.displayExposure);

  const bool showVel = st.showVelocity;
  const float velScale = st.velocityVizScale;
  const auto& uF = st.sim.u();
  const auto& vF = st.sim.v();

  auto tonemap = [&](float x) -> float {
    // Simple filmic-ish curve that behaves well for unbounded dye values.
    x = std::max(0.0f, x);
    return 1.0f - std::exp(-exposure * x);
  };

  for (int j = 1; j <= n; ++j) {
    const int outY = j - 1;
    for (int i = 1; i <= n; ++i) {
      const int outX = i - 1;
      const std::size_t outIdx = ((std::size_t)outY * (std::size_t)n + (std::size_t)outX) * 4u;
      const int id = i + stride * j;

      float rf = tonemap(rF[(std::size_t)id]);
      float gf = tonemap(gF[(std::size_t)id]);
      float bf = tonemap(bF[(std::size_t)id]);

      if (showVel) {
        const float vx = uF[(std::size_t)id] * velScale;
        const float vy = vF[(std::size_t)id] * velScale;
        const float sp = std::sqrt(vx * vx + vy * vy);
        rf = std::clamp(0.5f + 0.5f * vx, 0.0f, 1.0f);
        gf = std::clamp(0.5f + 0.5f * vy, 0.0f, 1.0f);
        bf = std::clamp(0.25f + 0.75f * sp, 0.0f, 1.0f);
      }

      st.rgba[outIdx + 0] = (unsigned char)std::clamp((int)std::lround(rf * 255.0f), 0, 255);
      st.rgba[outIdx + 1] = (unsigned char)std::clamp((int)std::lround(gf * 255.0f), 0, 255);
      st.rgba[outIdx + 2] = (unsigned char)std::clamp((int)std::lround(bf * 255.0f), 0, 255);
      st.rgba[outIdx + 3] = 255;
    }
  }

  st.tex.updateRGBA(0, 0, n, n, st.rgba.data());
}

void handlePreviewMouse(ProceduralFluidLabWindowState& st) {
  if (!ImGui::IsItemHovered()) {
    st.mouseDown = false;
    return;
  }

  const ImVec2 min = ImGui::GetItemRectMin();
  const ImVec2 max = ImGui::GetItemRectMax();
  const ImVec2 mp = ImGui::GetMousePos();

  const float w = std::max(1.0f, max.x - min.x);
  const float h = std::max(1.0f, max.y - min.y);

  float u = (mp.x - min.x) / w;
  float v = 1.0f - (mp.y - min.y) / h;
  u = std::clamp(u, 0.0f, 1.0f);
  v = std::clamp(v, 0.0f, 1.0f);

  const bool lDown = ImGui::IsMouseDown(0);
  const bool rDown = ImGui::IsMouseDown(1);
  const bool active = lDown || rDown;
  if (!active) {
    st.mouseDown = false;
    return;
  }

  const bool start = (ImGui::IsMouseClicked(0) || ImGui::IsMouseClicked(1) || !st.mouseDown);
  float du = 0.0f;
  float dv = 0.0f;
  if (!start) {
    du = u - st.lastMouseU;
    dv = v - st.lastMouseV;
  }

  st.lastMouseU = u;
  st.lastMouseV = v;
  st.mouseDown = true;

  float dyeR = st.brushDye * st.brushColor[0];
  float dyeG = st.brushDye * st.brushColor[1];
  float dyeB = st.brushDye * st.brushColor[2];

  if (rDown && st.rightClickErases) {
    dyeR = -dyeR;
    dyeG = -dyeG;
    dyeB = -dyeB;
  }

  const float velX = du * st.brushForce;
  const float velY = dv * st.brushForce;

  st.sim.splat(u, v, st.brushRadius01, velX, velY, dyeR, dyeG, dyeB);
}

float computeDt(ProceduralFluidLabWindowState& st, float timeSec) {
  if (st.autoDt) {
    const float dt = std::max(0.0f, timeSec - st.lastTimeSec);
    return std::min(dt, st.maxDt);
  }
  return st.fixedDt;
}

} // namespace

void drawProceduralFluidLabWindow(ProceduralFluidLabWindowState& st, float timeSec, const ToastFn& toast) {
  if (!st.open) return;

  ensureSimSize(st);
  ensureTexture(st);

  // Simulation params are edited live.
  auto& p = st.sim.params();

  ImGui::SetNextWindowSize(ImVec2(1200, 820), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Procedural Fluid Lab", &st.open)) {
    ImGui::End();
    return;
  }

  if (ImGui::BeginTable("##fluidlab", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV)) {
    ImGui::TableSetupColumn("Controls", ImGuiTableColumnFlags_WidthStretch);
    ImGui::TableSetupColumn("Preview", ImGuiTableColumnFlags_WidthFixed, (float)st.previewPixels + 36.0f);

    // ---------------- Controls ----------------
    ImGui::TableNextColumn();

    if (ImGui::CollapsingHeader("Simulation", ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::Checkbox("Paused", &st.paused);
      ImGui::SameLine();
      if (ImGui::Button("Step")) {
        st.singleStep = true;
      }

      ImGui::SetNextItemWidth(160.0f);
      ImGui::SliderInt("Grid", &st.gridSize, 32, 512);

      ImGui::SetNextItemWidth(160.0f);
      ImGui::SliderInt("Iterations", &st.iterations, 5, 60);

      // Live solver diagnostics (from the last simulation step).
      const auto& diag = st.sim.stats();
      ImGui::SeparatorText("Diagnostics");
      if (diag.usedMultigrid) {
        ImGui::Text("Projection: Multigrid (%d levels, %d V-cycles)", diag.multigridLevels, diag.multigridVCycles);
      } else {
        ImGui::Text("Projection: Gauss-Seidel (%d iters)", diag.gaussSeidelIterations);
      }
      ImGui::Text("Div max: %.4g  rms: %.4g", diag.maxAbsDivergence, diag.rmsDivergence);
      ImGui::Text("Pressure residual rms: %.4g", diag.pressureResidualRms);

      ImGui::Checkbox("Auto dt", &st.autoDt);
      if (!st.autoDt) {
        ImGui::SetNextItemWidth(160.0f);
        ImGui::SliderFloat("Fixed dt", &st.fixedDt, 1.0f / 240.0f, 1.0f / 20.0f, "%.4f");
      }
      ImGui::SetNextItemWidth(160.0f);
      ImGui::SliderFloat("Max dt", &st.maxDt, 1.0f / 120.0f, 1.0f / 10.0f, "%.4f");

      if (ImGui::Button("Clear")) {
        st.sim.clear();
      }
      ImGui::SameLine();
      if (ImGui::Button("Randomize Noise Seed")) {
        // Cheap deterministic-ish permutation.
        st.noiseSeed ^= (st.noiseSeed << 13);
        st.noiseSeed ^= (st.noiseSeed >> 7);
        st.noiseSeed ^= (st.noiseSeed << 17);
      }
    }

    if (ImGui::CollapsingHeader("Solver Params", ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::SetNextItemWidth(220.0f);
      ImGui::DragFloat("Viscosity", &p.viscosity, 0.00005f, 0.0f, 0.05f, "%.5f");
      ImGui::SetNextItemWidth(220.0f);
      ImGui::DragFloat("Dye Diffusion", &p.diffusion, 0.00005f, 0.0f, 0.05f, "%.5f");
      ImGui::SetNextItemWidth(220.0f);
      ImGui::DragFloat("Dye Dissipation", &p.dyeDissipation, 0.01f, 0.0f, 5.0f, "%.3f /s");
      ImGui::SetNextItemWidth(220.0f);
      ImGui::DragFloat("Velocity Damping", &p.velocityDamping, 0.01f, 0.0f, 5.0f, "%.3f /s");

      ImGui::SeparatorText("Projection");
      ImGui::Checkbox("Multigrid projection", &p.multigridProjection);
      if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("Accelerates the pressure Poisson solve using multigrid V-cycles.\n"
                          "Keeps the Stable Fluids boundary conditions, but converges much faster at higher grids.");
      }
      if (p.multigridProjection) {
        ImGui::SetNextItemWidth(220.0f);
        ImGui::SliderInt("MG Pre-smooth", &p.multigridPreSmooth, 1, 4);
        ImGui::SetNextItemWidth(220.0f);
        ImGui::SliderInt("MG Post-smooth", &p.multigridPostSmooth, 1, 4);
        ImGui::SetNextItemWidth(220.0f);
        ImGui::SliderInt("MG Coarse iters", &p.multigridCoarseIters, 4, 200);
      }

      ImGui::SeparatorText("Advection Quality");
      ImGui::SetNextItemWidth(220.0f);
      ImGui::SliderFloat("Dye BFECC", &p.dyeAdvectionCorrection, 0.0f, 1.0f, "%.2f");
      if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("0 = plain semi-Lagrangian (more diffusive)\n1 = BFECC/MacCormack correction (sharper dye)");
      }
      ImGui::Checkbox("Clamp BFECC extrema", &p.dyeAdvectionClamp);
      if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("Clamps corrected dye to the source stencil min/max to avoid ringing/overshoot.");
      }

      ImGui::SeparatorText("Detail");
      ImGui::SetNextItemWidth(220.0f);
      ImGui::DragFloat("Vorticity Confinement", &p.vorticityConfinement, 0.10f, 0.0f, 50.0f, "%.2f");
      ImGui::SetNextItemWidth(220.0f);
      ImGui::DragFloat("Curl Noise Strength", &p.curlNoiseStrength, 0.10f, 0.0f, 50.0f, "%.2f");
      ImGui::SetNextItemWidth(220.0f);
      ImGui::DragFloat("Curl Noise Frequency", &p.curlNoiseFrequency, 0.05f, 0.05f, 20.0f, "%.2f");
      ImGui::SetNextItemWidth(220.0f);
      ImGui::DragFloat("Curl Noise Time Scale", &p.curlNoiseTimeScale, 0.01f, 0.0f, 5.0f, "%.2f");

      ImGui::SeparatorText("Safety Clamps");
      ImGui::SetNextItemWidth(220.0f);
      ImGui::DragFloat("Max Speed", &p.maxSpeed, 0.5f, 0.0f, 500.0f, "%.1f");
      ImGui::SetNextItemWidth(220.0f);
      ImGui::DragFloat("Max Dye", &p.maxDye, 0.5f, 0.0f, 500.0f, "%.1f");
    }

    if (ImGui::CollapsingHeader("Brush", ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::TextDisabled("Paint into the sim by dragging on the preview.");
      ImGui::SetNextItemWidth(220.0f);
      ImGui::SliderFloat("Radius", &st.brushRadius01, 0.005f, 0.20f, "%.3f");
      ImGui::SetNextItemWidth(220.0f);
      ImGui::SliderFloat("Dye", &st.brushDye, 0.0f, 30.0f, "%.2f");
      ImGui::SetNextItemWidth(220.0f);
      ImGui::SliderFloat("Force", &st.brushForce, 0.0f, 200.0f, "%.1f");
      ImGui::ColorEdit3("Color", st.brushColor);
      ImGui::Checkbox("Right-click erases", &st.rightClickErases);
    }

    if (ImGui::CollapsingHeader("Display / Export", ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::Checkbox("Show velocity", &st.showVelocity);
      ImGui::SetNextItemWidth(220.0f);
      ImGui::SliderFloat("Exposure", &st.displayExposure, 0.005f, 0.25f, "%.3f");
      ImGui::SetNextItemWidth(220.0f);
      ImGui::SliderFloat("Vel Viz Scale", &st.velocityVizScale, 0.001f, 0.10f, "%.4f");
      ImGui::SetNextItemWidth(220.0f);
      ImGui::SliderInt("Preview px", &st.previewPixels, 256, 900);

      ImGui::InputText("Export Path", st.exportPath, sizeof(st.exportPath));
      ImGui::Checkbox("Flip Y", &st.exportFlipY);
      if (ImGui::Button("Export PNG")) {
        std::string err;
        try {
          const fs::path outPath(st.exportPath);
          if (!outPath.parent_path().empty()) {
            fs::create_directories(outPath.parent_path());
          }
        } catch (const std::exception& e) {
          err = e.what();
        }

        if (err.empty()) {
          if (!writePixelsToPng(st.exportPath, st.sim.gridSize(), st.sim.gridSize(), 4,
                                st.rgba.data(), st.sim.gridSize() * 4, st.exportFlipY, &err)) {
            toast(std::string("Export failed: ") + err, 4.0);
          } else {
            toast(std::string("Exported: ") + st.exportPath, 2.5);
          }
        } else {
          toast(std::string("Export failed: ") + err, 4.0);
        }
      }
    }

    // ---------------- Preview ----------------
    ImGui::TableNextColumn();

    const int n = st.sim.gridSize();
    if (n <= 0 || st.tex.handle() == 0) {
      ImGui::TextDisabled("Preview unavailable.");
    } else {
      const float size = (float)st.previewPixels;
      ImGui::TextDisabled("%dx%d", n, n);

      ImGui::Image((ImTextureID)(intptr_t)st.tex.handle(), ImVec2(size, size), ImVec2(0, 1), ImVec2(1, 0));
      handlePreviewMouse(st);

      ImGui::Spacing();
      ImGui::TextDisabled("Left drag: dye + push. Right drag: erase (optional).");
    }

    ImGui::EndTable();
  }

  // Step the sim and refresh the texture *after* the Image command has been
  // queued. Since ImGui renders later, the updated GL texture will be visible
  // in this same frame.
  const float dt = computeDt(st, timeSec);
  st.lastTimeSec = timeSec;

  const bool run = (!st.paused) || st.singleStep;
  if (run) {
    const float stepDt = st.singleStep ? (st.autoDt ? std::min(st.fixedDt, st.maxDt) : st.fixedDt) : dt;
    if (stepDt > 0.0f) {
      st.sim.step(stepDt, st.iterations, st.noiseSeed, timeSec);
    }
    st.singleStep = false;
  }

  updateTextureFromSim(st);

  ImGui::End();
}

} // namespace stellar::game
