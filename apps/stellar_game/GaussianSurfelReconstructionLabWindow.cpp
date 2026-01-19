#include "GaussianSurfelReconstructionLabWindow.h"

#include "stellar/core/Log.h"
#include "stellar/math/Mat4.h"

#include "stellar/ui/ImGuiCompat.h"

#include <SDL_opengl.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <future>
#include <string>
#include <utility>

namespace stellar::game {
namespace {

using Clock = std::chrono::high_resolution_clock;

template <class T>
bool isFutureReady(std::future<T>& f) {
  if (!f.valid()) return false;
  return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

void matToFloat(const math::Mat4d& m, float out[16]) {
  for (int i = 0; i < 16; ++i) out[i] = (float)m.m[i];
}

math::Vec3d orbitCameraPos(double yawDeg, double pitchDeg, double dist) {
  const double yaw = math::degToRad(yawDeg);
  const double pitch = math::degToRad(pitchDeg);

  const double cy = std::cos(yaw);
  const double sy = std::sin(yaw);
  const double cp = std::cos(pitch);
  const double sp = std::sin(pitch);

  // Y is up. Yaw rotates around Y; pitch around X' (orbit-style).
  const double x = dist * (cp * sy);
  const double y = dist * sp;
  const double z = dist * (cp * cy);
  return {x, y, z};
}

void resetReconstruction(GaussianSurfelReconstructionLabWindowState& s) {
  s.triObsCount.assign(s.triangles.size(), 0u);
  s.triToSurfel.assign(s.triangles.size(), -1);
  s.surfels.clear();
  s.currentViewIdx = -1;
  s.capturedViews.clear();
  s.planned = {};
  s.coveredOnceArea = 0.0;
  s.completeArea = 0.0;
  s.autoCapture = false;
  s.autoCaptureRemaining = 0;
}

void updateCoverageMetrics(GaussianSurfelReconstructionLabWindowState& s) {
  double coveredOnce = 0.0;
  double complete = 0.0;
  const std::uint8_t req = (std::uint8_t)std::max(1, s.requiredObservations);

  for (std::size_t i = 0; i < s.triangles.size(); ++i) {
    const double a = s.triangles[i].area;
    const std::uint8_t c = s.triObsCount.empty() ? 0u : s.triObsCount[i];

    if (c > 0u) coveredOnce += a;
    if (c >= req) complete += a;
  }

  s.coveredOnceArea = coveredOnce;
  s.completeArea = complete;
}

void replan(GaussianSurfelReconstructionLabWindowState& s) {
  if (!s.candidatesReady) return;
  if (s.triangles.empty()) return;

  const std::uint8_t req = (std::uint8_t)std::max(1, s.requiredObservations);
  s.planned = render::planNextBestPathBeam(s.candidateVisibility,
                                          s.candidateViewPositions,
                                          s.triangles,
                                          s.triObsCount,
                                          req,
                                          s.currentViewIdx,
                                          s.planner,
                                          &s.capturedViews);
}

void captureView(GaussianSurfelReconstructionLabWindowState& s, int viewIdx) {
  if (!s.candidatesReady) return;
  if (viewIdx < 0 || viewIdx >= (int)s.candidateVisibility.size()) return;

  const std::uint8_t req = (std::uint8_t)std::max(1, s.requiredObservations);

  const auto& vv = s.candidateVisibility[(std::size_t)viewIdx].visibleTriangles;
  for (int triIdx : vv) {
    if (triIdx < 0 || triIdx >= (int)s.triangles.size()) continue;

    std::uint8_t& c = s.triObsCount[(std::size_t)triIdx];
    if (c < req) c += 1;

    const int existing = s.triToSurfel[(std::size_t)triIdx];
    if (existing < 0) {
      render::GaussianSurfel surf = render::makeGaussianSurfelFromTriangle(
          s.triangles[(std::size_t)triIdx],
          (std::uint16_t)c,
          (std::uint16_t)req,
          s.surfelAlphaMul);

      surf.triIndex = triIdx;
      surf.obsCount = c;
      s.triToSurfel[(std::size_t)triIdx] = (int)s.surfels.size();
      s.surfels.push_back(surf);
    } else {
      render::GaussianSurfel& surf = s.surfels[(std::size_t)existing];
      surf.obsCount = c;
      const float conf = (float)math::clamp((double)c / (double)req, 0.0, 1.0);
      surf.opacity = math::clamp(s.surfelAlphaMul * conf, 0.0f, 1.0f);
    }
  }

  s.currentViewIdx = viewIdx;
  s.capturedViews.push_back(viewIdx);

  updateCoverageMetrics(s);
  if (s.autoPlan) replan(s);
}

void ensurePreviewInit(GaussianSurfelReconstructionLabWindowState& s) {
  if (s.previewReady) return;

  s.previewError.clear();

  if (!s.preview.init(s.previewRes, s.previewRes, &s.previewError)) {
    s.previewReady = false;
    return;
  }

  // Texture2D::createChecker is a convenience helper and does not currently
  // surface failure through a return value.
  s.checker.createChecker(256, 256, 16);

  if (!s.meshRenderer.init(&s.previewError)) {
    s.previewReady = false;
    return;
  }

  if (!s.surfelRenderer.init(&s.previewError)) {
    s.previewReady = false;
    return;
  }

  s.previewReady = true;
}

void renderPreview(GaussianSurfelReconstructionLabWindowState& s, float timeSec) {
  if (!s.meshGpuValid && s.surfels.empty()) return;
  ensurePreviewInit(s);
  if (!s.previewReady) return;

  const int w = s.previewRes;
  const int h = s.previewRes;

  // Simple auto-spin (non-destructive).
  const double yaw = s.autoSpin ? (double)s.camYawDeg + (double)s.spinDegPerSec * (double)timeSec : (double)s.camYawDeg;
  const math::Vec3d eye = orbitCameraPos(yaw, s.camPitchDeg, s.camDist);
  const math::Vec3d target{0, 0, 0};
  const math::Vec3d up{0, 1, 0};

  const math::Mat4d viewM = math::Mat4d::lookAt(eye, target, up);
  const math::Mat4d projM = math::Mat4d::perspective(math::degToRad((double)s.fovDeg), 1.0, 0.02, 50.0);

  float viewF[16];
  float projF[16];
  matToFloat(viewM, viewF);
  matToFloat(projM, projF);

  s.preview.begin();
  glViewport(0, 0, w, h);

  glEnable(GL_DEPTH_TEST);
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  glClearColor(0.06f, 0.07f, 0.08f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Draw ground-truth mesh (optional).
  if (s.showMesh && s.meshGpuValid) {
    s.meshRenderer.setMesh(&s.meshGpu);
    s.meshRenderer.setTexture(&s.checker);
    s.meshRenderer.setNormalTexture(nullptr);
    s.meshRenderer.setUnlit(false);
    s.meshRenderer.setLightPos(2.5f, 3.5f, 2.0f);
    s.meshRenderer.setCameraPos((float)eye.x, (float)eye.y, (float)eye.z);
    s.meshRenderer.setViewProj(viewF, projF);

    std::vector<render::InstanceData> inst(1);
    inst[0] = render::InstanceData{
        0.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 1.0f,
        0.0f, 0.0f, 0.0f, 1.0f,
        1.0f, 1.0f, 1.0f
    };

    if (s.meshWireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    s.meshRenderer.drawInstances(inst);
    if (s.meshWireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }

  // Draw surfels (optional).
  if (s.showSurfels && !s.surfels.empty()) {
    glDisable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);

    s.surfelRenderer.setViewProj(viewF, projF);
    // If ru/rv are treated as ~3σ extents, ~9 gives a very soft edge.
    s.surfelRenderer.setSharpness(9.0f);
    s.surfelRenderer.draw(s.surfels);

    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
    glEnable(GL_CULL_FACE);
  }

  s.preview.end();
}

void startMeshJob(GaussianSurfelReconstructionLabWindowState& s, const ToastFn& toast) {
  if (s.meshJobRunning) {
    s.meshJobQueued = true;
    s.queuedPreset = s.preset;
    s.queuedSeed = s.seed;
    s.queuedMesher = s.mesher;
    return;
  }

  s.meshDirty = false;
  s.meshJobRunning = true;

  const render::SdfGraphPreset preset = s.preset;
  const std::uint64_t seed = s.seed;
  const render::SdfMesherParams params = s.mesher;

  s.meshFuture = std::async(std::launch::async, [preset, seed, params]() -> SurfelReconMeshJobOutput {
    SurfelReconMeshJobOutput out{};
    const auto t0 = Clock::now();

    // Build a procedural graph and mesh it.
    const render::SdfGraph g = render::makeSdfGraphPreset(preset, (core::u64)seed);
    const render::ScalarField3D field = render::makeSdfField(g);

    out.mesh = render::meshIsosurfaceMarchingTetrahedra(field, params);

    const auto t1 = Clock::now();
    out.ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return out;
  });

  if (toast) toast("Meshing started…", 1.5);
}

void pollMeshJob(GaussianSurfelReconstructionLabWindowState& s, const ToastFn& toast) {
  if (!s.meshJobRunning) return;
  if (!isFutureReady(s.meshFuture)) return;

  SurfelReconMeshJobOutput out = s.meshFuture.get();
  s.meshJobRunning = false;

  if (!out.error.empty()) {
    s.meshGpuValid = false;
    if (toast) toast("Meshing failed: " + out.error, 4.0);
  } else {
    s.meshCpu = std::move(out.mesh);

    // Upload GPU mesh (must happen on the main thread / GL context).
    s.meshGpu.upload(s.meshCpu.vertices, s.meshCpu.indices);
    s.meshGpuValid = true;

    // Build triangle list for visibility/planning.
    s.triangles.clear();
    s.totalArea = 0.0;

    if (s.meshCpu.indices.size() >= 3) {
      s.triangles.reserve(s.meshCpu.indices.size() / 3);
      for (std::size_t i = 0; i + 2 < s.meshCpu.indices.size(); i += 3) {
        const std::uint32_t i0 = s.meshCpu.indices[i + 0];
        const std::uint32_t i1 = s.meshCpu.indices[i + 1];
        const std::uint32_t i2 = s.meshCpu.indices[i + 2];

        if (i0 >= s.meshCpu.vertices.size() ||
            i1 >= s.meshCpu.vertices.size() ||
            i2 >= s.meshCpu.vertices.size()) {
          continue;
        }

        const auto& v0 = s.meshCpu.vertices[i0];
        const auto& v1 = s.meshCpu.vertices[i1];
        const auto& v2 = s.meshCpu.vertices[i2];

        const math::Vec3d p0{v0.px, v0.py, v0.pz};
        const math::Vec3d p1{v1.px, v1.py, v1.pz};
        const math::Vec3d p2{v2.px, v2.py, v2.pz};

        render::TriangleGeom tri = render::makeTriangleGeom(p0, p1, p2);
        s.totalArea += tri.area;
        s.triangles.push_back(tri);
      }
    }

    resetReconstruction(s);
    s.visDirty = true;
    s.candidatesReady = false;

    if (toast) toast("Meshing done (" + std::to_string(out.ms) + " ms)", 2.0);
  }

  // If we were queued, start another mesh build immediately.
  if (s.meshJobQueued) {
    s.meshJobQueued = false;
    s.preset = s.queuedPreset;
    s.seed = s.queuedSeed;
    s.mesher = s.queuedMesher;
    s.meshDirty = true;
  }
}

void startVisibilityJob(GaussianSurfelReconstructionLabWindowState& s, const ToastFn& toast) {
  if (s.visJobRunning) {
    s.visJobQueued = true;
    return;
  }
  if (s.triangles.empty()) return;

  s.visDirty = false;
  s.visJobRunning = true;

  const int count = std::max(8, s.candidateCount);
  const double radius = std::max(0.25, (double)s.viewRadius);

  const render::VisibilitySettings vs{
      math::degToRad((double)s.fovDeg),
      std::max(16, s.visibilityRasterRes),
      s.occlusionAware,
      0.02,
      50.0};

  const std::vector<render::TriangleGeom> tris = s.triangles;

  s.visFuture = std::async(std::launch::async, [count, radius, vs, tris]() -> SurfelReconVisibilityJobOutput {
    SurfelReconVisibilityJobOutput out{};
    const auto t0 = Clock::now();

    // Fibonacci sphere / golden-angle distribution.
    const double goldenAngle = math::kPi * (3.0 - std::sqrt(5.0));

    out.viewPositions.reserve((std::size_t)count);
    out.visibility.reserve((std::size_t)count);

    for (int i = 0; i < count; ++i) {
      const double t = ((double)i + 0.5) / (double)count;
      const double y = 1.0 - 2.0 * t;
      const double r = std::sqrt(std::max(0.0, 1.0 - y * y));
      const double theta = goldenAngle * (double)i;

      const double x = std::cos(theta) * r;
      const double z = std::sin(theta) * r;

      const math::Vec3d pos = math::Vec3d{x, y, z} * radius;
      out.viewPositions.push_back(pos);

      const render::Viewpoint view = render::makeLookAtView(pos, math::Vec3d{0, 0, 0}, math::Vec3d{0, 1, 0});
      out.visibility.push_back(render::computeViewVisibility(tris, view, vs));
    }

    const auto t1 = Clock::now();
    out.ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return out;
  });

  if (toast) toast("Visibility precompute started…", 1.5);
}

void pollVisibilityJob(GaussianSurfelReconstructionLabWindowState& s, const ToastFn& toast) {
  if (!s.visJobRunning) return;
  if (!isFutureReady(s.visFuture)) return;

  SurfelReconVisibilityJobOutput out = s.visFuture.get();
  s.visJobRunning = false;

  if (!out.error.empty()) {
    s.lastVisibilityError = out.error;
    s.candidatesReady = false;
    if (toast) toast("Visibility failed: " + out.error, 4.0);
  } else {
    s.candidateViewPositions = std::move(out.viewPositions);
    s.candidateVisibility = std::move(out.visibility);
    s.candidatesReady = true;
    s.lastVisibilityError.clear();
    s.lastVisibilityMs = out.ms;

    if (toast) toast("Visibility done (" + std::to_string(out.ms) + " ms)", 2.0);

    if (s.autoPlan) replan(s);
  }

  // If queued, recompute again next frame (picks up latest settings).
  if (s.visJobQueued) {
    s.visJobQueued = false;
    s.visDirty = true;
  }
}

} // namespace

void drawGaussianSurfelReconstructionLabWindow(GaussianSurfelReconstructionLabWindowState& s,
                                               float timeSec,
                                               const ToastFn& toast) {
  if (!s.open) return;

  // Poll any async jobs.
  pollMeshJob(s, toast);
  pollVisibilityJob(s, toast);

  // Auto-capture: one planned view per frame.
  if (s.autoCapture && s.autoCaptureRemaining > 0 && s.candidatesReady) {
    if (s.planned.views.empty()) replan(s);
    if (!s.planned.views.empty()) {
      captureView(s, s.planned.views[0]);
      s.autoCaptureRemaining -= 1;
    } else {
      s.autoCapture = false;
      s.autoCaptureRemaining = 0;
    }
  } else if (s.autoCaptureRemaining <= 0) {
    s.autoCapture = false;
  }

  if (ImGui::Begin("Gaussian Surfel Reconstruction Lab", &s.open, ImGuiWindowFlags_NoCollapse)) {
    // Two-panel layout: controls on the left, preview on the right.
    ImGui::BeginChild("##gsurfel_controls", ImVec2(360, 0), true);

    ImGui::SeparatorText("Procedural Object");

    // Preset selector.
    const char* presetName = render::sdfGraphPresetName(s.preset);
    if (ImGui::BeginCombo("Preset", presetName)) {
      for (int i = 0; i <= (int)render::SdfGraphPreset::BooleanDemo; ++i) {
        const render::SdfGraphPreset p = (render::SdfGraphPreset)i;
        const bool selected = p == s.preset;
        if (ImGui::Selectable(render::sdfGraphPresetName(p), selected)) {
          s.preset = p;
          s.meshDirty = true;
        }
        if (selected) ImGui::SetItemDefaultFocus();
      }
      ImGui::EndCombo();
    }

    if (ImGui::InputScalar("Seed", ImGuiDataType_U64, &s.seed)) {
      s.meshDirty = true;
    }

    ImGui::Checkbox("Auto rebuild mesh", &s.autoRebuildMesh);

    ImGui::SeparatorText("Meshing (SDF → Triangles)");
    bool mesherChanged = false;
    mesherChanged |= ImGui::DragInt("Resolution", &s.mesher.resolution, 1.0f, 8, 256);
    mesherChanged |= ImGui::DragFloat("Bounds", &s.mesher.bounds, 0.05f, 0.25f, 10.0f);
    mesherChanged |= ImGui::DragFloat("Iso", &s.mesher.iso, 0.005f, -1.0f, 1.0f);
    mesherChanged |= ImGui::Checkbox("Compute normals", &s.mesher.computeNormalsFromField);
    if (mesherChanged) s.meshDirty = true;

    if (ImGui::Button("Rebuild Mesh")) {
      s.meshDirty = true;
    }
    ImGui::SameLine();
    if (s.meshJobRunning) {
      ImGui::TextUnformatted("(building…)"); 
    } else if (s.meshGpuValid) {
      ImGui::Text("Triangles: %zu", s.triangles.size());
    }

    if (s.meshDirty && s.autoRebuildMesh) {
      startMeshJob(s, toast);
    }

    ImGui::SeparatorText("Candidate Views + Visibility");
    bool visChanged = false;
    visChanged |= ImGui::DragInt("Candidates", &s.candidateCount, 1.0f, 8, 256);
    visChanged |= ImGui::DragFloat("View radius", &s.viewRadius, 0.05f, 0.25f, 20.0f);
    visChanged |= ImGui::DragFloat("FOV (deg)", &s.fovDeg, 0.25f, 15.0f, 120.0f);
    visChanged |= ImGui::Checkbox("Occlusion aware", &s.occlusionAware);
    if (s.occlusionAware) {
      visChanged |= ImGui::DragInt("Raster res", &s.visibilityRasterRes, 1.0f, 16, 256);
    }
    if (visChanged && s.meshGpuValid) {
      s.visDirty = true;
    }

    if (ImGui::Button("Recompute Visibility") && s.meshGpuValid) {
      s.visDirty = true;
    }
    ImGui::SameLine();
    ImGui::Checkbox("Auto plan", &s.autoPlan);

    if (s.visDirty && s.meshGpuValid) {
      startVisibilityJob(s, toast);
    }

    if (!s.lastVisibilityError.empty()) {
      ImGui::TextColored(ImVec4(1, 0.25f, 0.25f, 1), "Vis error: %s", s.lastVisibilityError.c_str());
    } else if (s.candidatesReady) {
      ImGui::Text("Vis: %.2f ms  (%d candidates)", s.lastVisibilityMs, (int)s.candidateViewPositions.size());
    } else {
      ImGui::TextUnformatted("(visibility not ready)");
    }

    ImGui::SeparatorText("Reconstruction Model");
    ImGui::DragInt("Required obs", &s.requiredObservations, 1.0f, 1, 6);
    ImGui::DragFloat("Surfel alpha", &s.surfelAlphaMul, 0.01f, 0.0f, 1.0f);

    ImGui::SeparatorText("Planner (Beam Search)");
    bool plannerChanged = false;
    plannerChanged |= ImGui::DragInt("Horizon", &s.planner.horizon, 1.0f, 1, 8);
    plannerChanged |= ImGui::DragInt("Beam width", &s.planner.beamWidth, 1.0f, 1, 32);
    plannerChanged |= ImGui::DragInt("Branching", &s.planner.localBranching, 1.0f, 1, 32);
    plannerChanged |= ImGui::DragDouble("Move cost", &s.planner.moveCostWeight, 0.01f, 0.0, 2.0, "%.3f");
    plannerChanged |= ImGui::DragDouble("Revisit penalty", &s.planner.revisitPenalty, 0.01f, 0.0, 5.0, "%.3f");
    if (plannerChanged && s.autoPlan && s.candidatesReady) {
      replan(s);
    }

    if (ImGui::Button("Plan Next Path") && s.candidatesReady) {
      replan(s);
    }
    ImGui::SameLine();
    if (ImGui::Button("Capture Next") && s.candidatesReady) {
      if (s.planned.views.empty()) replan(s);
      if (!s.planned.views.empty()) {
        captureView(s, s.planned.views[0]);
      }
    }

    ImGui::SameLine();
    if (ImGui::Button("Reset Recon") && s.meshGpuValid) {
      resetReconstruction(s);
      if (s.autoPlan) replan(s);
    }

    if (s.candidatesReady) {
      ImGui::Checkbox("Auto capture", &s.autoCapture);
      ImGui::SameLine();
      ImGui::DragInt("Steps", &s.autoCaptureRemaining, 1.0f, 0, 256);
    }

    // Progress metrics.
    if (s.totalArea > 0.0) {
      const float p1 = (float)(s.coveredOnceArea / s.totalArea);
      const float p2 = (float)(s.completeArea / s.totalArea);

      ImGui::Text("Coverage: %.1f%%  Complete: %.1f%%", 100.0f * p1, 100.0f * p2);
      ImGui::ProgressBar(p2, ImVec2(-1, 0), "complete");
      ImGui::ProgressBar(p1, ImVec2(-1, 0), "covered");
    }

    // Planned path view.
    if (!s.planned.views.empty()) {
      ImGui::Text("Planned score: %.3f", s.planned.score);
      ImGui::BeginChild("##planned_path", ImVec2(0, 88), true);
      for (std::size_t i = 0; i < s.planned.views.size(); ++i) {
        ImGui::Text("%zu: view %d", i, s.planned.views[i]);
      }
      ImGui::EndChild();
    }

    ImGui::SeparatorText("Preview");
    ImGui::Checkbox("Show mesh", &s.showMesh);
    ImGui::SameLine();
    ImGui::Checkbox("Wireframe", &s.meshWireframe);
    ImGui::Checkbox("Show surfels", &s.showSurfels);
    ImGui::Checkbox("Auto spin", &s.autoSpin);
    if (s.autoSpin) {
      ImGui::SameLine();
      ImGui::DragFloat("Spin deg/s", &s.spinDegPerSec, 0.25f, 0.0f, 90.0f);
    }
    ImGui::DragFloat("Yaw", &s.camYawDeg, 0.25f, -180.0f, 180.0f);
    ImGui::DragFloat("Pitch", &s.camPitchDeg, 0.25f, -89.0f, 89.0f);
    ImGui::DragFloat("Distance", &s.camDist, 0.05f, 0.1f, 25.0f);

    ImGui::EndChild();

    ImGui::SameLine();

    ImGui::BeginChild("##gsurfel_preview", ImVec2(0, 0), true);

    if (!s.previewError.empty()) {
      ImGui::TextColored(ImVec4(1, 0.25f, 0.25f, 1), "Preview init failed: %s", s.previewError.c_str());
    }

    if (s.meshGpuValid || !s.surfels.empty()) {
      renderPreview(s, timeSec);

      if (s.previewReady) {
        const ImVec2 size((float)s.preview.width(), (float)s.preview.height());
        ImGui::Image((ImTextureID)(intptr_t)s.preview.color().handle(), size, ImVec2(0, 1), ImVec2(1, 0));
      } else {
        ImGui::TextUnformatted("(preview not ready)");
      }
    } else {
      ImGui::TextUnformatted("(generate a mesh first)");
    }

    ImGui::EndChild();
  }

  ImGui::End();
}

} // namespace stellar::game
