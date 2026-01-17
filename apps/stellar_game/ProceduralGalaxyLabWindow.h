#pragma once

#include "stellar/core/LruCache.h"
#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/proc/GalaxyGenerator.h"
#include "stellar/proc/GalaxyRegions.h"
#include "stellar/proc/Hyperlanes.h"
#include "stellar/proc/HyperlaneRouter.h"
#include "stellar/proc/HyperlaneKRoutes.h"
#include "stellar/proc/HyperlaneVulnerability.h"
#include "stellar/render/RenderTarget.h"
#include "stellar/render/ShaderToy.h"
#include "stellar/sim/Celestial.h"
#include "stellar/sim/Faction.h"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace stellar::game {

struct ProceduralGalaxyLabWindowState {
  bool open{false};

  core::u64 seed{1337ull};
  int factionCount{8};

  proc::GalaxyParams params{};

  // View / preview controls (top-down XY projection).
  math::Vec3d centerLy{0.0, 0.0, 0.0};
  double viewRadiusLy{250.0};
  double zHalfLy{200.0};
  std::size_t maxStubs{20000};

  // Preview performance:
  // Cache generated sector stubs so small view changes (pan/zoom) don't
  // constantly re-run deterministic sector generation.
  bool previewUseSectorCache{true};
  int previewSectorCacheMaxEntries{20000};

  // Parallel preview generation (CPU). Uses the core::JobSystem thread pool to
  // generate and filter sampled sectors concurrently.
  bool previewParallel{true};
  // 0 = automatic (JobSystem default / hardware_concurrency); otherwise clamp to [1..64].
  int previewParallelThreads{0};
  // Diagnostics: number of workers used in the last rebuild.
  int previewParallelWorkersUsed{1};

  // Diagnostics (reset each regenerate).
  int previewSectorCacheHits{0};
  int previewSectorCacheMisses{0};
  int previewSectorCacheEvictions{0};

  // Preview sampling diagnostics (auto LOD).
  int previewSectorStrideXY{1};
  int previewSectorsGenerated{0};
  int previewCandidateStubs{0};

  bool autoRegenerate{true};
  bool colorByFaction{false};
  bool colorByRegion{false};
  bool colorByCluster{false};
  bool colorByVoid{false};
  bool colorByMorphology{false};
  bool showArmGuides{true};
  bool showMorphGuides{true};
  bool showLegend{true};

  // Hyperlane overlay (procedural navigation graph) for the current preview subset.
  bool showHyperlanes{false};

  // Hyperlane generation params (see proc::Hyperlanes).
  proc::HyperlaneParams hyperlaneParams{};
  // Cap the number of nodes used for lane generation (keeps O(n log n) work bounded).
  std::size_t hyperlaneMaxNodes{6000};
  // Cap edges drawn for readability. 0 = auto.
  std::size_t hyperlaneMaxEdgesDraw{12000};
  float hyperlaneOpacity{0.35f};
  float hyperlaneMinLenPx{2.0f};
  float hyperlaneWidthPx{1.0f};
  // 0: Risk, 1: Bandwidth, 2: Traffic (betweenness), 3: Criticality (bridges/articulation)
  int hyperlaneColorMode{0};

  // Traffic / chokepoint analytics (approx. betweenness centrality).
  // 0 = exact (all sources), higher = faster approximation.
  int hyperlaneTrafficSamples{96};
  bool hyperlaneTrafficHighlightChokepoints{true};
  int hyperlaneTrafficTopEdges{12};
  float hyperlaneTrafficWidthBoost{1.0f};
  float hyperlaneTrafficAlphaBoost{1.0f};
  bool hyperlaneHighlightHovered{true};

  // Structural criticality analytics (bridges + articulation points).
  // This highlights fragile connections where losing a single lane or system splits regions.
  bool hyperlaneCriticalHighlightBridges{true};
  bool hyperlaneCriticalHighlightArticulation{true};
  int hyperlaneCriticalTopBridges{12};
  int hyperlaneCriticalTopArticulation{12};
  float hyperlaneCriticalWidthBoost{1.0f};
  float hyperlaneCriticalAlphaBoost{1.0f};
  float hyperlaneCriticalNodeRingRadiusPx{6.0f};
  float hyperlaneCriticalNodeRingThicknessPx{2.0f};

  // Cached hyperlane network for the current stub set.
  bool hyperlanesDirty{true};
  core::u64 previewStubSetHash{0};
  core::u64 hyperlaneLastHash{core::u64(-1)};
  proc::HyperlaneNetwork hyperlaneNet{};
  // The exact node subset used to generate the cached hyperlane network (bounded by hyperlaneMaxNodes).
  std::vector<sim::SystemStub> hyperlaneNodes{};
  int hyperlaneNodesUsed{0};
  int hyperlaneEdgesUsed{0};
  double hyperlaneGenMs{0.0};

  // Hyperlane routing inspector (pathfinding over the cached lane network).
  //
  // Usage in the UI:
  //   Shift+Click = set route start
  //   Ctrl+Click  = set route end
  //   Right-click = clear
  bool hyperlaneRouteEnabled{true};
  bool hyperlaneRouteDraw{true};
  float hyperlaneRouteOpacity{0.90f};
  float hyperlaneRouteWidthPx{2.75f};

  sim::SystemId hyperlaneRouteFrom{0};
  sim::SystemId hyperlaneRouteTo{0};

  proc::HyperlaneTravelParams hyperlaneTravel{};

  // K-shortest (loopless) alternate routes between Route Start and Route End.
  //
  //  - k = 1 behaves like a classic single best path.
  //  - k > 1 enumerates alternatives using Yen's algorithm (with A* as the sub-solver).
  int hyperlaneRouteK{3};
  int hyperlaneRouteSelected{0};
  bool hyperlaneRouteDrawAlternatives{false};
  float hyperlaneRouteAltOpacity{0.25f};
  float hyperlaneRouteAltWidthScale{0.75f};

  // Cached router adjacency for the current hyperlane network (rebuilt when hyperlaneLastHash changes).
  bool hyperlaneRouterInited{false};
  core::u64 hyperlaneRouterHash{core::u64(-1)};
  proc::HyperlaneRouter hyperlaneRouter{};

  // Cached route results.
  bool hyperlaneRouteDirty{true};
  core::u64 hyperlaneRouteLastHash{core::u64(-1)};
  std::vector<sim::SystemId> hyperlaneRoutePath{};
  std::vector<proc::HyperlaneKRoute> hyperlaneKRoutes{};
  proc::HyperlanePathMetrics hyperlaneRouteMetrics{};
  double hyperlaneRouteMs{0.0};
  std::string hyperlaneRouteError{};

  // Hyperlane traffic analytics (betweenness centrality).
  bool hyperlaneCentralityDirty{true};
  core::u64 hyperlaneCentralityLastHash{core::u64(-1)};
  std::vector<double> hyperlaneEdgeBetweenness{};
  std::vector<double> hyperlaneNodeBetweenness{};
  std::vector<sim::SystemId> hyperlaneCentralityNodeIds{};
  double hyperlaneCentralityMaxEdge{0.0};
  double hyperlaneCentralityMaxNode{0.0};
  double hyperlaneCentralityMs{0.0};
  std::vector<int> hyperlaneChokepointEdgeIdx{};

  // Hyperlane structural vulnerability (bridges + articulation points).
  bool hyperlaneVulnerabilityDirty{true};
  core::u64 hyperlaneVulnerabilityLastHash{core::u64(-1)};
  proc::HyperlaneVulnerabilityResult hyperlaneVulnerability{};
  double hyperlaneVulnerabilityMs{0.0};

  // Cached top-N lists for rendering (rebuilt with the vulnerability analysis).
  std::vector<int> hyperlaneCriticalBridgeEdgeIdx{};
  std::vector<sim::SystemId> hyperlaneCriticalArticulationNodeIds{};
  std::vector<float> hyperlaneCriticalArticulationNodeImpact01{};

  // GPU heatmap background (custom rendering). This is an analytic / approximate view
  // that helps visualize global density structure quickly.
  bool heatmapEnabled{true};
  // 1 = full, 2 = half, 4 = quarter resolution relative to the UI canvas.
  int heatmapResolutionDiv{2};
  // 0: Density, 1: Spiral, 2: Clusters, 3: Voids, 4: Morphology
  int heatmapMode{0};
  float heatmapExposure{0.25f};
  float heatmapGamma{1.15f};
  bool heatmapContours{true};
  int heatmapContourCount{14};
  float heatmapContourWidth{0.030f};

  // Dynamic LOD / filtering for the heatmap shader. Uses screen-space derivatives
  // (dFdx/dFdy/fwidth) to band-limit high-frequency procedural detail when it
  // becomes sub-pixel during zoom-out. This reduces shimmer and improves performance.
  bool heatmapDynamicLod{true};
  // When a feature is smaller than heatmapLodPxLo pixels it is mostly filtered out;
  // at heatmapLodPxHi pixels and above it is fully visible.
  float heatmapLodPxLo{0.65f};
  float heatmapLodPxHi{2.25f};
  // 0 = let energy naturally fall off (smoother); 1 = preserve energy (crisper).
  float heatmapLodEnergy{0.0f};
  // If enabled, skip expensive blob scans (clusters/voids) when the feature radius is sub-pixel.
  bool heatmapLodSkipTinyBlobs{true};

  // If enabled, the heatmap re-renders every frame using iTime (useful for animated shaders).
  // When disabled, we cache the rendered texture and only re-render when inputs change.
  bool heatmapAnimate{false};

  // Cached key for the last rendered heatmap (so we can skip re-rendering when nothing changed).
  core::u64 heatmapLastHash{core::u64(-1)};
  int heatmapLastW{0};
  int heatmapLastH{0};


  // Galaxy regions (Worley/Voronoi). Used for visualization only.
  double regionCellSizeLy{900.0};

  // Cached preview data.
  bool dirty{true};
  std::vector<sim::Faction> factions{};
  std::vector<sim::SystemStub> stubs{};

  // --- internal: preview sector cache (LRU) ---
  struct SectorCoord {
    int x{0};
    int y{0};
    int z{0};
    bool operator==(const SectorCoord&) const = default;
  };

  struct SectorCoordHash {
    std::size_t operator()(const SectorCoord& c) const noexcept {
      // A small hash-combine (64-bit friendly). We don't need stable hashing
      // across platforms here; this is an in-memory cache key.
      std::size_t h = std::hash<int>{}(c.x);
      h ^= std::hash<int>{}(c.y) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
      h ^= std::hash<int>{}(c.z) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
      return h;
    }
  };

  using SectorStubVec = std::vector<sim::SystemStub>;
  using SectorStubVecPtr = std::shared_ptr<const SectorStubVec>;

  core::u64 previewSectorCacheCtxHash{0};
  std::size_t previewSectorCacheStubTotal{0};
  core::LruCache<SectorCoord, SectorStubVecPtr, SectorCoordHash> previewSectorCache{};

  // Cached per-stub region kind/id for preview rendering.
  std::vector<proc::GalaxyRegionKind> stubRegionKind{};
  std::vector<core::u64> stubRegionId{};

  // Cached per-stub cluster influence (0..1) for preview rendering.
  std::vector<float> stubCluster01{};

  // Cached per-stub void influence (0..1) for preview rendering.
  std::vector<float> stubVoid01{};

  // Cached per-stub morphology influence (0..1) for preview rendering.
  std::vector<float> stubMorph01{};
  double lastGenMs{0.0};

  // --- internal: heatmap render resources ---
  bool heatmapInited{false};
  stellar::render::RenderTarget2D heatmapTarget{};
  stellar::render::ShaderToy heatmapToy{};
  stellar::render::ShaderToyParamSet heatmapParams{};

  // Cached parameter indices for internal heatmap uniforms.
  // This avoids repeated string lookups when updating dozens of uniforms per draw.
  struct HeatmapParamHandles {
    bool inited{false};

    // View
    int viewCenter{-1};
    int viewCenterZLy{-1};
    int viewRadiusLy{-1};

    // Base density
    int baseMean{-1};
    int radiusLy{-1};
    int radialScaleLy{-1};
    int thicknessLy{-1};
    int zHalfLy{-1};

    // Warp / flare
    int warpAmplitudeLy{-1};
    int warpStartRadiusLy{-1};
    int warpPower{-1};
    int warpLobes{-1};
    int warpPhaseDeg{-1};
    int warpNoiseStrength{-1};
    int warpNoiseFreq{-1};
    int flareStrength{-1};
    int flarePower{-1};

    // Spiral
    int spiralArmCount{-1};
    int spiralArmStrength{-1};
    int spiralPitchDeg{-1};
    int spiralArmWidthDeg{-1};
    int spiralArmPhaseDeg{-1};
    int spiralArmNoiseStrength{-1};
    int spiralArmNoiseFreq{-1};

    // Density noise
    int densityNoiseStrength{-1};
    int densityNoiseFreq{-1};

    // Clusters
    int clusterStrength{-1};
    int clusterCellSizeLy{-1};
    int clusterChancePerCell{-1};
    int clusterRadiusLy{-1};
    int clusterRadiusJitter{-1};
    int clusterStrengthJitter{-1};
    int clusterFalloffPower{-1};

    // Voids
    int voidStrength{-1};
    int voidCellSizeLy{-1};
    int voidChancePerCell{-1};
    int voidRadiusLy{-1};
    int voidRadiusJitter{-1};
    int voidStrengthJitter{-1};
    int voidFalloffPower{-1};

    // Morphology
    int barStrength{-1};
    int barAngleDeg{-1};
    int barLengthLy{-1};
    int barWidthLy{-1};
    int barPower{-1};

    int ringStrength{-1};
    int ringRadiusLy{-1};
    int ringWidthLy{-1};
    int ringPower{-1};

    // Display
    int mode{-1};
    int exposure{-1};
    int gamma{-1};
    int contours{-1};
    int contourCount{-1};
    int contourWidth{-1};

    // Dynamic LOD
    int lodEnabled{-1};
    int lodPxLo{-1};
    int lodPxHi{-1};
    int lodEnergy{-1};
    int lodSkipTinyBlobs{-1};
  };
  HeatmapParamHandles heatmapHandles{};
  std::string heatmapError{};
  int heatmapFrame{0};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

void drawProceduralGalaxyLabWindow(ProceduralGalaxyLabWindowState& state, float timeSec, const ToastFn& toast);

} // namespace stellar::game
