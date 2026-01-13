#pragma once

#include "stellar/render/GaussianSurfels.h"
#include "stellar/render/Mesh.h"
#include "stellar/render/MeshRenderer.h"
#include "stellar/render/RenderTarget.h"
#include "stellar/render/SdfGraph.h"
#include "stellar/render/SdfMesher.h"
#include "stellar/render/Texture.h"

#include <cstdint>
#include <functional>
#include <future>
#include <string>
#include <vector>

namespace stellar::game {

struct SurfelReconMeshJobOutput {
  render::SdfMeshData mesh{};
  double ms{0.0};
  std::string error{};
};

struct SurfelReconVisibilityJobOutput {
  std::vector<math::Vec3d> viewPositions{};
  std::vector<render::ViewVisibility> visibility{};
  double ms{0.0};
  std::string error{};
};

struct GaussianSurfelReconstructionLabWindowState {
  bool open{false};

  // Procedural object source (SDF graph presets).
  render::SdfGraphPreset preset{render::SdfGraphPreset::Crystal};
  std::uint64_t seed{1337};

  // Mesher settings (marching tetrahedra via SdfMesher).
  render::SdfMesherParams mesher{};
  bool autoRebuildMesh{true};

  // Candidate viewpoints + visibility evaluation.
  int candidateCount{64};
  float viewRadius{3.0f};
  float fovDeg{60.0f};
  bool occlusionAware{true};
  int visibilityRasterRes{96};

  // Reconstruction model.
  int requiredObservations{2};
  float surfelAlphaMul{0.85f};

  // Next-best-path planner.
  render::PathPlannerSettings planner{};
  bool autoPlan{true};
  // Optional: capture one planned view per frame.
  bool autoCapture{false};
  int autoCaptureRemaining{0};


  // Preview camera (orbit).
  int previewRes{512};
  float camYawDeg{35.0f};
  float camPitchDeg{25.0f};
  float camDist{3.2f};
  bool autoSpin{false};
  float spinDegPerSec{10.0f};

  bool showMesh{true};
  bool meshWireframe{false};
  bool showSurfels{true};

  // Async mesh job bookkeeping.
  bool meshDirty{true};
  bool meshJobRunning{false};
  bool meshJobQueued{false};
  render::SdfGraphPreset queuedPreset{render::SdfGraphPreset::Crystal};
  std::uint64_t queuedSeed{1337};
  render::SdfMesherParams queuedMesher{};
  std::future<SurfelReconMeshJobOutput> meshFuture{};

  // Async visibility job bookkeeping.
  bool visDirty{true};
  bool visJobRunning{false};
  bool visJobQueued{false};
  std::future<SurfelReconVisibilityJobOutput> visFuture{};

  // Built mesh + derived triangle list (CPU).
  render::SdfMeshData meshCpu{};
  render::Mesh meshGpu{};
  bool meshGpuValid{false};

  std::vector<render::TriangleGeom> triangles{};
  double totalArea{0.0};

  // Candidate views + precomputed visibility sets.
  std::vector<math::Vec3d> candidateViewPositions{};
  std::vector<render::ViewVisibility> candidateVisibility{};
  bool candidatesReady{false};
  double lastVisibilityMs{0.0};
  std::string lastVisibilityError{};

  // Reconstruction state.
  std::vector<std::uint8_t> triObsCount{};
  std::vector<int> triToSurfel{};
  std::vector<render::GaussianSurfel> surfels{};

  int currentViewIdx{-1};
  std::vector<int> capturedViews{};
  render::PlannedViewPath planned{};

  double coveredOnceArea{0.0};
  double completeArea{0.0};

  // GPU preview resources.
  render::RenderTarget2D preview{};
  render::MeshRenderer meshRenderer{};
  render::GaussianSurfelRenderer surfelRenderer{};
  render::Texture2D checker{};

  bool previewReady{false};
  std::string previewError{};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

void drawGaussianSurfelReconstructionLabWindow(GaussianSurfelReconstructionLabWindowState& state,
                                               float timeSec,
                                               const ToastFn& toast);

} // namespace stellar::game
