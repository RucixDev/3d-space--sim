#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/Mesh.h"
#include "stellar/render/MeshRenderer.h"
#include "stellar/render/MeshSimplify.h"
#include "stellar/render/ProceduralGraph.h"
#include "stellar/render/RenderTarget.h"
#include "stellar/render/SdfRaymarcher.h"
#include "stellar/render/SdfGraph.h"
#include "stellar/render/SdfMesher.h"
#include "stellar/render/SdfDualContouring.h"
#include "stellar/render/Texture.h"

#include <functional>
#include <future>
#include <string>
#include <vector>

namespace stellar::game {

struct ProceduralMeshLabWindowState {
  bool open{false};

  // Quickstart: start from a preset, then tweak nodes.
  render::SdfGraphPreset preset{render::SdfGraphPreset::Asteroid};
  bool lockToPreset{true};

  // Graph seed
  core::u64 seed{0xC0FFEEULL};

  // Mesher settings
  enum class MesherType : core::u8 {
    MarchingTetrahedra = 0,
    DualContouring = 1,
  };

  MesherType mesherType{MesherType::MarchingTetrahedra};
  render::SdfMesherParams mesher{};

  // Dual Contouring extras (used when mesherType == DualContouring).
  float dcQefRegularization{1.0e-6f};
  bool dcClampToCell{true};
  bool dcProjectToIso{true};
  int dcProjectIterations{2};
  bool autoRemesh{true};
  bool asyncRemesh{true};

  // ---- LOD / mesh simplification ----
  // Builds an LOD chain from the generated mesh using QEM simplification.
  bool buildLods{true};
  bool autoBuildLods{true};
  int lodLevels{3};
  float lodRatioPerLevel{0.5f};
  int previewLod{0};
  int exportLod{0};
  bool exportAllLods{false};

  // Preview settings
  int previewResolution{512};
  float yawDeg{25.0f};
  float pitchDeg{-20.0f};
  float distance{2.5f};
  bool autoSpin{true};
  float spinDegPerSec{18.0f};
  bool wireframe{false};

  // Preview mode: mesh rasterization (triangles) OR direct SDF raymarch.
  bool previewRaymarch{true};
  bool raymarchAutoRender{true};
  float raymarchFovDeg{55.0f};
  bool showRaymarchShader{false};
  render::SdfRaymarchSettings raymarchSettings{};
  std::string raymarchError{};

  // Export
  char exportPath[256]{"screenshots/procedural_mesh.obj"};
  char exportGltfPath[256]{"screenshots/procedural_mesh.gltf"};
  bool exportWithMaterial{true};
  bool exportTextureFlipY{true};

  // glTF export extras
  bool exportGltfTangents{true};
  bool exportGltfPackLodsMsft{false};

  // glTF PBR factors (used for export only)
  float exportPbrMetallic{0.0f};
  float exportPbrRoughness{1.0f};

  // Graph I/O
  char graphPath[256]{"proc_graphs/asteroid.sdfgraph"};

  // Optional procedural surface texture (2D ProcGraph baked to a Texture2D)
  bool useProceduralTexture{true};
  bool applyTextureToRaymarch{true};
  render::ProcGraphPreset texPreset{render::ProcGraphPreset::Rocky};
  bool texLockToPreset{true};
  // Track to avoid reapplying preset every frame.
  render::ProcGraphPreset texAppliedPreset{render::ProcGraphPreset::Nebula};
  core::u64 texAppliedSeed{0};
  core::u64 texSeed{0xDEC0ADDEULL};
  int texResolution{512};
  bool texAutoBake{true};
  bool texGenerateMips{true};
  float texDitherStrength{1.0f};
  bool texPackHeightInAlpha{false};
  bool showTexShader{false};
  char texGraphPath[256]{"proc_graphs/rocky.procgraph"};
  render::ProcGraph texGraph{render::ProcGraph::makeDefault()};
  bool texDirty{true};
  render::ProceduralGraphBaker texBaker{};
  std::string texError{};
  render::SdfRaymarchMaterial raymarchMat{};

  // Graph data
  render::SdfGraph graph{render::SdfGraph::makeDefault()};
  bool dirty{true};

  // ---- Variation Studio (SDF idea generator) ----
  //
  // Generates mutated variants of the current SDF graph so you can quickly
  // explore interesting silhouettes without hand-editing every parameter.
  //
  // This is deliberately "tooling-grade": variants are cheap to generate and
  // come with quick sampled stats, but you still do a full remesh for final
  // export/LOD building.
  struct SdfVariantCandidate {
    render::SdfGraph graph{};
    core::u64 key{0}; // fingerprint for dedupe
    float insideFraction{0.0f}; // [0,1]
    float boundingRadius{0.0f};
    float score{0.0f};
    float minDistance{0.0f};
    float maxDistance{0.0f};
    int samples{0};
  };

  bool showVariationStudio{false};
  bool variationAutoRegenerate{false};
  core::u64 variationSeed{0xA11CEBEEF1234567ULL};
  int variationCount{8};
  float variationParamJitter{0.18f};
  float variationRewireChance{0.06f};
  float variationOpChance{0.04f};
  float variationWrapChance{0.15f};
  bool variationMutateNodeSeeds{true};
  int variationSampleRes{12};
  float variationTargetFill{0.25f};
  float variationTargetRadius{0.85f};
  bool variationSortByScore{true};
  std::vector<SdfVariantCandidate> variants{};
  int selectedVariant{-1};
  std::string variantsError{};


  // Last CPU mesh (kept for export)
  render::SdfMeshData cpuMesh{};
  render::SdfMeshStats cpuStats{};
  double lastMeshMs{0.0};
  bool hasMesh{false};
  std::string lastError{};

  // LOD chain (LOD0 is cpuMesh; lodMeshes stores LOD1..LOD(N)).
  std::vector<render::SdfMeshData> lodMeshes{};
  std::vector<render::SdfMeshStats> lodStats{};
  std::vector<render::MeshSimplifyStats> lodSimplifyStats{};
  double lastLodMs{0.0};
  bool hasLods{false};
  std::string lodError{};

  // GPU preview resources
  bool previewInited{false};
  std::string previewInitError{};
  render::RenderTarget2D previewTarget{};
  render::Mesh previewMesh{};
  render::MeshRenderer previewRenderer{};
  render::Texture2D checker{};

  // GPU raymarch preview resources
  render::SdfRaymarcher raymarcher{};

  // Async job
  struct MeshJobResult {
    int id{0};
    render::SdfMeshData mesh;
    render::SdfMeshStats stats;
    double ms{0.0};
    std::string error;
  };

  int nextJobId{1};
  int latestJobId{0};
  std::future<MeshJobResult> job{};
  bool jobRunning{false};

  // Async LOD job (separate from the mesher job)
  struct LodJobResult {
    int id{0};
    std::vector<render::SdfMeshData> meshes;
    std::vector<render::SdfMeshStats> stats;
    std::vector<render::MeshSimplifyStats> simplifyStats;
    double ms{0.0};
    std::string error;
  };

  int nextLodJobId{1};
  int latestLodJobId{0};
  std::future<LodJobResult> lodJob{};
  bool lodJobRunning{false};

  // ---- Undo / redo history (graph + key settings) ----
  //
  // This is a lightweight, CPU-only history intended for tooling iteration.
  // Snapshots are small value types (graphs + editable settings) so undo/redo is
  // deterministic and safe to use with async meshing.
  struct HistorySnapshot {
    // SDF graph + preset context.
    render::SdfGraph graph{};
    render::SdfGraphPreset preset{render::SdfGraphPreset::Asteroid};
    bool lockToPreset{true};
    core::u64 seed{0};

    // Mesher knobs.
    MesherType mesherType{MesherType::MarchingTetrahedra};
    render::SdfMesherParams mesher{};

    // Dual Contouring extras.
    float dcQefRegularization{1.0e-6f};
    bool dcClampToCell{true};
    bool dcProjectToIso{true};
    int dcProjectIterations{2};

    // LOD chain knobs.
    bool buildLods{true};
    int lodLevels{3};
    float lodRatioPerLevel{0.5f};
    int previewLod{0};
    int exportLod{0};
    bool exportAllLods{false};

    // Preview mode knobs.
    bool previewRaymarch{true};
    render::SdfRaymarchSettings raymarchSettings{};

    // Procedural surface texture knobs.
    bool useProceduralTexture{true};
    bool applyTextureToRaymarch{true};
    render::ProcGraphPreset texPreset{render::ProcGraphPreset::Rocky};
    bool texLockToPreset{true};
    core::u64 texSeed{0};
    int texResolution{512};
    bool texGenerateMips{true};
    float texDitherStrength{1.0f};
    bool texPackHeightInAlpha{false};
    render::ProcGraph texGraph{};
  };

  bool historyEnabled{true};
  int historyMax{64};
  float historyCoalesceSec{0.35f};

  bool historyBaselineValid{false};
  bool historyPending{false};
  float historyPendingSince{0.0f};

  HistorySnapshot historyBaseline{};
  std::vector<HistorySnapshot> historyUndo{};
  std::vector<HistorySnapshot> historyRedo{};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

void drawProceduralMeshLabWindow(ProceduralMeshLabWindowState& st, float timeSec, const ToastFn& toast);

} // namespace stellar::game
