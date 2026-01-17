#include "ProceduralGalaxyLabWindow.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/core/JobSystem.h"
#include "stellar/proc/GalaxyClusters.h"
#include "stellar/proc/GalaxyVoids.h"
#include "stellar/proc/HyperlaneCentrality.h"

#include "stellar/math/Math.h"

#include <imgui.h>

#include "stellar/render/Gl.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <mutex>
#include <future>
#include <atomic>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace stellar::game {

using Clock = std::chrono::high_resolution_clock;


// -----------------------------
// Preview generation helpers
// -----------------------------

static core::JobSystem& previewJobs(int requestedThreads, int* outThreadsUsed) {
  // Persist a single pool across regenerations so UI tweaks don't constantly
  // create/join threads.
  static std::unique_ptr<core::JobSystem> gPool{};
  static int gRequested{-1};

  int req = requestedThreads;
  if (req < 0) req = 0;
  if (req > 64) req = 64;

  if (!gPool || gRequested != req) {
    // JobSystem treats 0 as 'auto'.
    gPool = std::make_unique<core::JobSystem>(static_cast<std::size_t>(req));
    gRequested = req;
  }

  if (outThreadsUsed) *outThreadsUsed = static_cast<int>(gPool->threadCount());
  return *gPool;
}


// -----------------------------
// GPU heatmap background (custom rendering)
// -----------------------------

static int seed32(core::u64 s) {
  // Stable-ish mapping to 32-bit for GLSL.
  const core::u64 h = core::hashCombine(s, core::fnv1a64("galaxyHeatmap"));
  return (int)(h & 0x7FFFFFFFu);

}


// Hash the heatmap inputs into a stable key so we can skip expensive re-renders
// when nothing that affects the heatmap has changed.
//
// This is essentially a dirty-flag for an offscreen RT: when inputs change,
// the key changes and we redraw once; otherwise we reuse the previous texture.
// (This saves a lot of GPU time when the heatmap is purely analytic/static.)

template <typename T>
static core::u64 hashMix(core::u64 h, const T& v) {
  return core::hashCombine(h, core::hashBytes(&v, sizeof(T)));
}

static core::u64 hashGalaxyParams(core::u64 h0, const proc::GalaxyParams& p) {
  h0 = hashMix(h0, p.sectorSizeLy);
  h0 = hashMix(h0, p.radiusLy);
  h0 = hashMix(h0, p.thicknessLy);
  h0 = hashMix(h0, p.radialScaleLengthLy);
  h0 = hashMix(h0, p.baseMeanSystemsPerSector);

  h0 = hashMix(h0, p.minSystemSeparationLy);

  h0 = hashMix(h0, p.spiralArmCount);
  h0 = hashMix(h0, p.spiralArmStrength);
  h0 = hashMix(h0, p.spiralPitchDeg);
  h0 = hashMix(h0, p.spiralArmWidthDeg);
  h0 = hashMix(h0, p.spiralArmPhaseDeg);
  h0 = hashMix(h0, p.spiralArmNoiseStrength);
  h0 = hashMix(h0, p.spiralArmNoiseFreq);

  h0 = hashMix(h0, p.densityNoiseStrength);
  h0 = hashMix(h0, p.densityNoiseFreq);

  h0 = hashMix(h0, p.clusterStrength);
  h0 = hashMix(h0, p.clusterCellSizeLy);
  h0 = hashMix(h0, p.clusterChancePerCell);
  h0 = hashMix(h0, p.clusterRadiusLy);
  h0 = hashMix(h0, p.clusterRadiusJitter);
  h0 = hashMix(h0, p.clusterStrengthJitter);
  h0 = hashMix(h0, p.clusterFalloffPower);

  h0 = hashMix(h0, p.voidStrength);
  h0 = hashMix(h0, p.voidCellSizeLy);
  h0 = hashMix(h0, p.voidChancePerCell);
  h0 = hashMix(h0, p.voidRadiusLy);
  h0 = hashMix(h0, p.voidRadiusJitter);
  h0 = hashMix(h0, p.voidStrengthJitter);
  h0 = hashMix(h0, p.voidFalloffPower);

  h0 = hashMix(h0, p.barStrength);
  h0 = hashMix(h0, p.barAngleDeg);
  h0 = hashMix(h0, p.barLengthLy);
  h0 = hashMix(h0, p.barWidthLy);
  h0 = hashMix(h0, p.barPower);

  h0 = hashMix(h0, p.ringStrength);
  h0 = hashMix(h0, p.ringRadiusLy);
  h0 = hashMix(h0, p.ringWidthLy);
  h0 = hashMix(h0, p.ringPower);

  h0 = hashMix(h0, p.warpAmplitudeLy);
  h0 = hashMix(h0, p.warpStartRadiusLy);
  h0 = hashMix(h0, p.warpPower);
  h0 = hashMix(h0, p.warpLobes);
  h0 = hashMix(h0, p.warpPhaseDeg);
  h0 = hashMix(h0, p.warpNoiseStrength);
  h0 = hashMix(h0, p.warpNoiseFreq);

  h0 = hashMix(h0, p.flareStrength);
  h0 = hashMix(h0, p.flarePower);

  return h0;
}

static core::u64 hashHyperlaneParams(core::u64 h0, const proc::HyperlaneParams& p) {
  h0 = hashMix(h0, p.maxNeighborDistanceLy);
  h0 = hashMix(h0, p.neighborK);
  h0 = hashMix(h0, p.forceConnected);
  h0 = hashMix(h0, p.minDegree);
  h0 = hashMix(h0, p.extraEdgeChance);
  h0 = hashMix(h0, p.regionCellSizeLy);

  // Optional galaxy-scale hazard modulation (nebulae / storms).
  h0 = hashMix(h0, p.hazards.enabled);
  h0 = hashMix(h0, p.hazards.timeDays);
  h0 = hashMix(h0, p.hazards.driftLyPerDay);
  h0 = hashMix(h0, p.hazards.nebulaScaleLy);
  h0 = hashMix(h0, p.hazards.stormScaleLy);
  h0 = hashMix(h0, p.hazards.edgeSamples);
  h0 = hashMix(h0, p.hazards.hazardToCost);
  h0 = hashMix(h0, p.hazards.hazardToRisk);
  h0 = hashMix(h0, p.hazards.hazardToBandwidth);

  // Cap to keep generation bounded. (size_t -> u64 for stable hashing)
  h0 = hashMix(h0, (core::u64)p.maxEdges);

  return h0;
}

static core::u64 hashHyperlaneTravelParams(core::u64 h0, const proc::HyperlaneTravelParams& p) {
  h0 = hashMix(h0, p.riskWeight);
  h0 = hashMix(h0, p.bandwidthBias);
  h0 = hashMix(h0, p.minBandwidthFactor);
  return h0;
}


static core::u64 heatmapRenderKey(const ProceduralGalaxyLabWindowState& st, int w, int h) {
  core::u64 h0 = core::fnv1a64("galaxyHeatmapRenderKey/v1");

  // Render target dimensions.
  h0 = hashMix(h0, w);
  h0 = hashMix(h0, h);

  // View.
  h0 = hashMix(h0, st.seed);
  h0 = hashMix(h0, st.centerLy.x);
  h0 = hashMix(h0, st.centerLy.y);
  h0 = hashMix(h0, st.centerLy.z);
  h0 = hashMix(h0, st.viewRadiusLy);
  h0 = hashMix(h0, st.zHalfLy);

  // Heatmap display + shader LOD controls.
  h0 = hashMix(h0, st.heatmapResolutionDiv);
  h0 = hashMix(h0, st.heatmapMode);
  h0 = hashMix(h0, st.heatmapExposure);
  h0 = hashMix(h0, st.heatmapGamma);
  h0 = hashMix(h0, st.heatmapContours);
  h0 = hashMix(h0, st.heatmapContourCount);
  h0 = hashMix(h0, st.heatmapContourWidth);

  h0 = hashMix(h0, st.heatmapDynamicLod);
  h0 = hashMix(h0, st.heatmapLodPxLo);
  h0 = hashMix(h0, st.heatmapLodPxHi);
  h0 = hashMix(h0, st.heatmapLodEnergy);
  h0 = hashMix(h0, st.heatmapLodSkipTinyBlobs);

  // Procedural galaxy parameters that affect the analytic density.
  h0 = hashGalaxyParams(h0, st.params);

  return h0;
}

static void addHeatmapParam(stellar::render::ShaderToyParamSet& ps,
                           stellar::render::ShaderToyParamType t,
                           const char* name) {
  stellar::render::ShaderToyParamDef d{};
  d.type = t;
  d.name = name;
  // These heatmap uniforms are internal (synced from ProceduralGalaxyLabWindowState).
  // They are not edited through the ShaderLab param UI.
  //
  // We update them via ShaderToyParamSet::setRawValue so values like radii,
  // cell sizes, and world-space centers are never clamped to any UI range.
  d.defaultValue = {0, 0, 0, 0};
  ps.defs.push_back(d);
}

static void buildHeatmapParamSchema(stellar::render::ShaderToyParamSet& ps) {
  ps.clear();

  using stellar::render::ShaderToyParamType;

  // View
  addHeatmapParam(ps, ShaderToyParamType::Vec2,  "gViewCenter");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gViewCenterZLy");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gViewRadiusLy");

  // Base density
  addHeatmapParam(ps, ShaderToyParamType::Float, "gBaseMean");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gRadiusLy");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gRadialScaleLy");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gThicknessLy");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gZHalfLy");

  // Warp / flare (affects midplane + local thickness)
  addHeatmapParam(ps, ShaderToyParamType::Float, "gWarpAmplitudeLy");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gWarpStartRadiusLy");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gWarpPower");
  addHeatmapParam(ps, ShaderToyParamType::Int,   "gWarpLobes");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gWarpPhaseDeg");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gWarpNoiseStrength");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gWarpNoiseFreq");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gFlareStrength");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gFlarePower");

  // Spiral
  addHeatmapParam(ps, ShaderToyParamType::Int,   "gSpiralArmCount");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gSpiralArmStrength");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gSpiralPitchDeg");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gSpiralArmWidthDeg");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gSpiralArmPhaseDeg");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gSpiralArmNoiseStrength");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gSpiralArmNoiseFreq");

  // Density noise
  addHeatmapParam(ps, ShaderToyParamType::Float, "gDensityNoiseStrength");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gDensityNoiseFreq");

  // Clusters
  addHeatmapParam(ps, ShaderToyParamType::Float, "gClusterStrength");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gClusterCellSizeLy");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gClusterChancePerCell");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gClusterRadiusLy");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gClusterRadiusJitter");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gClusterStrengthJitter");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gClusterFalloffPower");

  // Voids
  addHeatmapParam(ps, ShaderToyParamType::Float, "gVoidStrength");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gVoidCellSizeLy");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gVoidChancePerCell");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gVoidRadiusLy");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gVoidRadiusJitter");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gVoidStrengthJitter");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gVoidFalloffPower");

  // Morphology (bar / ring)
  addHeatmapParam(ps, ShaderToyParamType::Float, "gBarStrength");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gBarAngleDeg");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gBarLengthLy");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gBarWidthLy");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gBarPower");

  addHeatmapParam(ps, ShaderToyParamType::Float, "gRingStrength");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gRingRadiusLy");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gRingWidthLy");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gRingPower");

  // Display
  addHeatmapParam(ps, ShaderToyParamType::Int,   "gMode");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gExposure");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gGamma");
  addHeatmapParam(ps, ShaderToyParamType::Bool,  "gContours");
  addHeatmapParam(ps, ShaderToyParamType::Int,   "gContourCount");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gContourWidth");

  // Dynamic LOD / band-limiting (shader-side).
  addHeatmapParam(ps, ShaderToyParamType::Bool,  "gLodEnabled");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gLodPxLo");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gLodPxHi");
  addHeatmapParam(ps, ShaderToyParamType::Float, "gLodEnergy");
  addHeatmapParam(ps, ShaderToyParamType::Bool,  "gLodSkipTinyBlobs");

  ps.resetToDefaults();
  ps.rebuildIndex();
}


static bool ensureHeatmapParamHandles(ProceduralGalaxyLabWindowState& st) {
  auto& h = st.heatmapHandles;
  if (h.inited) return st.heatmapError.empty();

  h = ProceduralGalaxyLabWindowState::HeatmapParamHandles{};
  h.inited = true;

  auto require = [&](const char* name) -> int {
    const int idx = st.heatmapParams.findIndex(name);
    if (idx < 0 && st.heatmapError.empty()) {
      st.heatmapError = std::string("Heatmap internal param schema mismatch: missing '") + name + "'";
    }
    return idx;
  };

  // View
  h.viewCenter = require("gViewCenter");
  h.viewCenterZLy = require("gViewCenterZLy");
  h.viewRadiusLy = require("gViewRadiusLy");

  // Base density
  h.baseMean = require("gBaseMean");
  h.radiusLy = require("gRadiusLy");
  h.radialScaleLy = require("gRadialScaleLy");
  h.thicknessLy = require("gThicknessLy");
  h.zHalfLy = require("gZHalfLy");

  // Warp / flare
  h.warpAmplitudeLy = require("gWarpAmplitudeLy");
  h.warpStartRadiusLy = require("gWarpStartRadiusLy");
  h.warpPower = require("gWarpPower");
  h.warpLobes = require("gWarpLobes");
  h.warpPhaseDeg = require("gWarpPhaseDeg");
  h.warpNoiseStrength = require("gWarpNoiseStrength");
  h.warpNoiseFreq = require("gWarpNoiseFreq");
  h.flareStrength = require("gFlareStrength");
  h.flarePower = require("gFlarePower");

  // Spiral
  h.spiralArmCount = require("gSpiralArmCount");
  h.spiralArmStrength = require("gSpiralArmStrength");
  h.spiralPitchDeg = require("gSpiralPitchDeg");
  h.spiralArmWidthDeg = require("gSpiralArmWidthDeg");
  h.spiralArmPhaseDeg = require("gSpiralArmPhaseDeg");
  h.spiralArmNoiseStrength = require("gSpiralArmNoiseStrength");
  h.spiralArmNoiseFreq = require("gSpiralArmNoiseFreq");

  // Density noise
  h.densityNoiseStrength = require("gDensityNoiseStrength");
  h.densityNoiseFreq = require("gDensityNoiseFreq");

  // Clusters
  h.clusterStrength = require("gClusterStrength");
  h.clusterCellSizeLy = require("gClusterCellSizeLy");
  h.clusterChancePerCell = require("gClusterChancePerCell");
  h.clusterRadiusLy = require("gClusterRadiusLy");
  h.clusterRadiusJitter = require("gClusterRadiusJitter");
  h.clusterStrengthJitter = require("gClusterStrengthJitter");
  h.clusterFalloffPower = require("gClusterFalloffPower");

  // Voids
  h.voidStrength = require("gVoidStrength");
  h.voidCellSizeLy = require("gVoidCellSizeLy");
  h.voidChancePerCell = require("gVoidChancePerCell");
  h.voidRadiusLy = require("gVoidRadiusLy");
  h.voidRadiusJitter = require("gVoidRadiusJitter");
  h.voidStrengthJitter = require("gVoidStrengthJitter");
  h.voidFalloffPower = require("gVoidFalloffPower");

  // Morphology
  h.barStrength = require("gBarStrength");
  h.barAngleDeg = require("gBarAngleDeg");
  h.barLengthLy = require("gBarLengthLy");
  h.barWidthLy = require("gBarWidthLy");
  h.barPower = require("gBarPower");

  h.ringStrength = require("gRingStrength");
  h.ringRadiusLy = require("gRingRadiusLy");
  h.ringWidthLy = require("gRingWidthLy");
  h.ringPower = require("gRingPower");

  // Display
  h.mode = require("gMode");
  h.exposure = require("gExposure");
  h.gamma = require("gGamma");
  h.contours = require("gContours");
  h.contourCount = require("gContourCount");
  h.contourWidth = require("gContourWidth");

  // Dynamic LOD
  h.lodEnabled = require("gLodEnabled");
  h.lodPxLo = require("gLodPxLo");
  h.lodPxHi = require("gLodPxHi");
  h.lodEnergy = require("gLodEnergy");
  h.lodSkipTinyBlobs = require("gLodSkipTinyBlobs");

  return st.heatmapError.empty();
}

static void syncHeatmapParamsFromState(ProceduralGalaxyLabWindowState& st) {
  // This is called from the heatmap render pass (potentially every frame).
  // Use cached parameter indices to avoid repeated string lookups.
  if (!st.heatmapHandles.inited) {
    if (!ensureHeatmapParamHandles(st)) return;
  }

  const auto& h = st.heatmapHandles;

  auto setF = [&](int idx, double v) {
    st.heatmapParams.setRawValue(idx, { (float)v, 0, 0, 0 });
  };
  auto setI = [&](int idx, int v) {
    st.heatmapParams.setRawValue(idx, { (float)v, 0, 0, 0 });
  };
  auto setB = [&](int idx, bool v) {
    st.heatmapParams.setRawValue(idx, { v ? 1.0f : 0.0f, 0, 0, 0 });
  };

  st.heatmapParams.setRawValue(h.viewCenter, { (float)st.centerLy.x, (float)st.centerLy.y, 0, 0 });
  setF(h.viewCenterZLy, st.centerLy.z);
  setF(h.viewRadiusLy, st.viewRadiusLy);

  setF(h.baseMean, st.params.baseMeanSystemsPerSector);
  setF(h.radiusLy, st.params.radiusLy);
  setF(h.radialScaleLy, st.params.radialScaleLengthLy);
  setF(h.thicknessLy, st.params.thicknessLy);
  setF(h.zHalfLy, st.zHalfLy);

  // Warp / flare
  setF(h.warpAmplitudeLy, st.params.warpAmplitudeLy);
  setF(h.warpStartRadiusLy, st.params.warpStartRadiusLy);
  setF(h.warpPower, st.params.warpPower);
  setI(h.warpLobes, st.params.warpLobes);
  setF(h.warpPhaseDeg, st.params.warpPhaseDeg);
  setF(h.warpNoiseStrength, st.params.warpNoiseStrength);
  setF(h.warpNoiseFreq, st.params.warpNoiseFreq);
  setF(h.flareStrength, st.params.flareStrength);
  setF(h.flarePower, st.params.flarePower);

  // Spiral
  setI(h.spiralArmCount, st.params.spiralArmCount);
  setF(h.spiralArmStrength, st.params.spiralArmStrength);
  setF(h.spiralPitchDeg, st.params.spiralPitchDeg);
  setF(h.spiralArmWidthDeg, st.params.spiralArmWidthDeg);
  setF(h.spiralArmPhaseDeg, st.params.spiralArmPhaseDeg);
  setF(h.spiralArmNoiseStrength, st.params.spiralArmNoiseStrength);
  setF(h.spiralArmNoiseFreq, st.params.spiralArmNoiseFreq);

  // Density noise
  setF(h.densityNoiseStrength, st.params.densityNoiseStrength);
  setF(h.densityNoiseFreq, st.params.densityNoiseFreq);

  // Clusters
  setF(h.clusterStrength, st.params.clusterStrength);
  setF(h.clusterCellSizeLy, st.params.clusterCellSizeLy);
  setF(h.clusterChancePerCell, st.params.clusterChancePerCell);
  setF(h.clusterRadiusLy, st.params.clusterRadiusLy);
  setF(h.clusterRadiusJitter, st.params.clusterRadiusJitter);
  setF(h.clusterStrengthJitter, st.params.clusterStrengthJitter);
  setF(h.clusterFalloffPower, st.params.clusterFalloffPower);

  // Voids
  setF(h.voidStrength, st.params.voidStrength);
  setF(h.voidCellSizeLy, st.params.voidCellSizeLy);
  setF(h.voidChancePerCell, st.params.voidChancePerCell);
  setF(h.voidRadiusLy, st.params.voidRadiusLy);
  setF(h.voidRadiusJitter, st.params.voidRadiusJitter);
  setF(h.voidStrengthJitter, st.params.voidStrengthJitter);
  setF(h.voidFalloffPower, st.params.voidFalloffPower);

  // Morphology (bar / ring)
  setF(h.barStrength, st.params.barStrength);
  setF(h.barAngleDeg, st.params.barAngleDeg);
  setF(h.barLengthLy, st.params.barLengthLy);
  setF(h.barWidthLy, st.params.barWidthLy);
  setF(h.barPower, st.params.barPower);

  setF(h.ringStrength, st.params.ringStrength);
  setF(h.ringRadiusLy, st.params.ringRadiusLy);
  setF(h.ringWidthLy, st.params.ringWidthLy);
  setF(h.ringPower, st.params.ringPower);

  // Display
  setI(h.mode, st.heatmapMode);
  setF(h.exposure, st.heatmapExposure);
  setF(h.gamma, st.heatmapGamma);
  setB(h.contours, st.heatmapContours);
  setI(h.contourCount, st.heatmapContourCount);
  setF(h.contourWidth, st.heatmapContourWidth);

  // Dynamic LOD controls.
  setB(h.lodEnabled, st.heatmapDynamicLod);

  // Keep the thresholds sane even if the UI temporarily allows overlap.
  const float pxLo = (float)std::max(0.05, (double)st.heatmapLodPxLo);
  const float pxHi = (float)std::max((double)pxLo + 1.0e-3, (double)st.heatmapLodPxHi);
  const float energy = (float)std::clamp((double)st.heatmapLodEnergy, 0.0, 1.0);

  st.heatmapParams.setRawValue(h.lodPxLo, { pxLo, 0, 0, 0 });
  st.heatmapParams.setRawValue(h.lodPxHi, { pxHi, 0, 0, 0 });
  st.heatmapParams.setRawValue(h.lodEnergy, { energy, 0, 0, 0 });
  setB(h.lodSkipTinyBlobs, st.heatmapLodSkipTinyBlobs);
}

static const char* kGalaxyHeatmapShader = R"GLSL(

float sat(float x) { return clamp(x, 0.0, 1.0); }

float wrapAnglePi(float a) {
  float twoPi = 2.0 * PI;
  a = mod(a + PI, twoPi);
  if (a < 0.0) a += twoPi;
  return a - PI;
}

float smooth01(float t) {
  t = sat(t);
  return t * t * (3.0 - 2.0 * t);
}

// Deterministic cell RNG with a salt. Uses the wrapper's integer hash library.
float rand1(ivec3 c, int salt) {
  return hash31(c + ivec3(salt, salt * 31, salt * 97));
}

vec3 rand3(ivec3 c, int salt) {
  return hash33(c + ivec3(salt, salt * 31, salt * 97));
}

// Streaming-safe blob influence field: scan a small neighborhood of cells for blob centers.
float blobInfluence(vec3 pos,
                    float cellSize,
                    float chance,
                    float radiusBase,
                    float radiusJitter,
                    float strengthJitter,
                    float falloffPower,
                    int salt) {
  if (cellSize <= 1.0e-6 || radiusBase <= 1.0e-6 || chance <= 0.0) return 0.0;

  chance = sat(chance);
  radiusJitter = sat(radiusJitter);
  strengthJitter = sat(strengthJitter);
  falloffPower = max(0.05, falloffPower);

  ivec3 cell = ivec3(floor(pos / cellSize));

  float maxR = radiusBase * (1.0 + radiusJitter);
  int search = int(clamp(ceil(maxR / cellSize) + 1.0, 1.0, 4.0));

  float best = 0.0;

  for (int dz = -4; dz <= 4; ++dz) {
    if (abs(dz) > search) continue;
    for (int dy = -4; dy <= 4; ++dy) {
      if (abs(dy) > search) continue;
      for (int dx = -4; dx <= 4; ++dx) {
        if (abs(dx) > search) continue;

        ivec3 c = cell + ivec3(dx, dy, dz);

        float r0 = rand1(c, salt);
        if (r0 >= chance) continue;

        vec3 jitter = rand3(c, salt + 11);
        vec3 center = (vec3(c) + jitter) * cellSize;

        float rMul = mix(1.0 - radiusJitter, 1.0 + radiusJitter, rand1(c, salt + 23));
        float rad = max(1.0, radiusBase * rMul);

        float sMul = mix(1.0 - strengthJitter, 1.0 + strengthJitter, rand1(c, salt + 37));

        float d = length(pos - center);
        if (d > rad) continue;

        float t = 1.0 - d / rad;
        float s = smooth01(t);
        float f = pow(s, falloffPower);
        float influ = f * sMul;
        best = max(best, influ);
      }
    }
  }

  return sat(best);
}

float spiralMask(vec2 p, float Rref) {
  if (gSpiralArmCount <= 0 || gSpiralArmStrength <= 0.0) return 0.0;

  int arms = clamp(gSpiralArmCount, 1, 8);
  float r = length(p);

  float pitchDeg = clamp(gSpiralPitchDeg, 1.0, 89.0);
  float pitch = radians(pitchDeg);
  float k = 1.0 / tan(pitch);

  float phase = radians(gSpiralArmPhaseDeg);
  float theta = atan(p.y, p.x);

  float lnTerm = log(max(1.0e-6, r / max(1.0e-6, Rref)));
  float twoPi = 2.0 * PI;

  float dMin = 1.0e9;
  for (int i = 0; i < 8; ++i) {
    if (i >= arms) break;
    float armTheta = k * lnTerm + phase + twoPi * (float(i) / float(arms));
    float d = abs(wrapAnglePi(theta - armTheta));
    dMin = min(dMin, d);
  }

  float width = radians(max(0.1, gSpiralArmWidthDeg));
  float x = dMin / max(1.0e-6, width);
  float mask = exp(-0.5 * x * x);

  // Fade arms near the center.
  float fade = sat((r - Rref * 0.25) / max(1.0, Rref * 0.35));
  mask *= fade;

  // Optional noise modulation to avoid perfectly smooth arms.
  if (gSpiralArmNoiseStrength > 0.0 && gSpiralArmNoiseFreq > 0.0) {
    vec2 np = p * gSpiralArmNoiseFreq + vec2(23.1, 91.7);
    float n = gLodEnabled ? fbm2_lod(np, gLodPxLo, gLodPxHi, gLodEnergy) : fbm2(np);
    n = sat(n);
    float mul = mix(1.0 - gSpiralArmNoiseStrength, 1.0 + gSpiralArmNoiseStrength, n);
    mask *= mul;
  }

  return sat(mask);
}

float barMask(vec2 p) {
  if (gBarStrength <= 0.0) return 0.0;

  float a = max(1.0e-6, abs(gBarLengthLy));
  float b = max(1.0e-6, abs(gBarWidthLy));
  float ang = radians(gBarAngleDeg);

  float ca = cos(ang);
  float sa = sin(ang);

  float xr = ca * p.x + sa * p.y;
  float yr = -sa * p.x + ca * p.y;

  float q = length(vec2(xr / a, yr / b));
  if (q >= 1.0) return 0.0;

  float power = max(1.0, gBarPower);
  float t = 1.0 - q;
  return sat(pow(t, power));
}

float ringMask(float r) {
  if (gRingStrength <= 0.0) return 0.0;

  float sigma = max(1.0e-6, abs(gRingWidthLy));
  float d = (r - gRingRadiusLy) / sigma;
  float m = exp(-0.5 * d * d);

  float power = max(1.0, gRingPower);
  if (abs(power - 1.0) > 1.0e-5) m = pow(m, power);

  return sat(m);
}

float thicknessHalf(float rxy) {
  float base = max(1.0, 0.5 * gThicknessLy);
  float fs = max(0.0, gFlareStrength);
  if (fs <= 0.0) return base;

  float R = max(1.0, gRadiusLy);
  float t0 = clamp(max(0.0, rxy) / R, 0.0, 1.0);
  float p = max(0.05, gFlarePower);
  float t = pow(t0, p);
  return base * (1.0 + fs * t);
}

// Warp midplane offset in ly (stylized; approximates CPU morphology helper).
float warpZLy(vec2 p, float rxy) {
  float amp0 = gWarpAmplitudeLy;
  if (abs(amp0) <= 1.0e-6) return 0.0;

  float start = max(0.0, gWarpStartRadiusLy);
  if (rxy <= start) return 0.0;

  float R = max(start + 1.0, gRadiusLy);
  float t0 = clamp((rxy - start) / max(1.0, (R - start)), 0.0, 1.0);
  float pw = max(0.05, gWarpPower);
  float amp = amp0 * pow(t0, pw);

  int lobes = gWarpLobes;
  if (lobes == 0) lobes = 2;
  lobes = clamp(abs(lobes), 1, 8);

  float theta = atan(p.y, p.x);
  float phase = radians(gWarpPhaseDeg);

  float w = amp * sin(float(lobes) * theta + phase);

  // Optional low-frequency modulation so the warp isn't perfectly sinusoidal.
  float ns = sat(gWarpNoiseStrength);
  float nf = gWarpNoiseFreq;
  if (ns > 0.0 && nf > 0.0) {
    vec2 np = p * nf + vec2(101.7, 37.2);
    float n = gLodEnabled ? fbm2_lod(np, gLodPxLo, gLodPxHi, gLodEnergy) : fbm2(np);
    n = sat(n);
    float mul = (1.0 - ns) + (2.0 * ns) * n; // [1-ns, 1+ns]
    w *= mul;
  }

  return w;
}

// Integral of exp(-|u|/halfT) over [u1,u2] where u1<=u2.
// halfT must be > 0.
float expAbsIntegral(float u1, float u2, float halfT) {
  // Both negative.
  if (u2 <= 0.0) {
    return halfT * (exp(u2 / halfT) - exp(u1 / halfT));
  }
  // Both positive.
  if (u1 >= 0.0) {
    return halfT * (exp(-u1 / halfT) - exp(-u2 / halfT));
  }
  // Crosses zero.
  return halfT * (2.0 - exp(u1 / halfT) - exp(-u2 / halfT));
}

// Fraction of the truncated vertical density profile captured by the current Z slice.
// The CPU generator hard-clips the disc to |z-warpZ| <= halfT.
float sliceVerticalFrac(float warpZ, float halfT) {
  float zHalf = max(0.0, gZHalfLy);
  float z0 = gViewCenterZLy;

  float zMin = z0 - zHalf;
  float zMax = z0 + zHalf;

  float discMin = warpZ - halfT;
  float discMax = warpZ + halfT;

  float a = max(zMin, discMin);
  float b = min(zMax, discMax);
  if (b <= a) return 0.0;

  float u1 = a - warpZ;
  float u2 = b - warpZ;

  float ht = max(1.0e-6, halfT);
  float integral = expAbsIntegral(u1, u2, ht);
  float total = 2.0 * halfT * (1.0 - exp(-1.0)); // integral over [-halfT, halfT]
  return sat(integral / max(1.0e-6, total));
}

vec2 toWorld(vec2 uv) {
  // Match the CPU worldToScreen mapping:
  // scale = min(w,h) / (2*viewRadius)
  // screen = center + (dx*scale, -dy*scale)
  vec2 fragPx = uv * iResolution;
  vec2 centerPx = 0.5 * iResolution;
  vec2 dPx = fragPx - centerPx;

  float minDim = min(iResolution.x, iResolution.y);
  float scale = minDim / (2.0 * max(1.0e-3, gViewRadiusLy));

  vec2 dLy = vec2(dPx.x, -dPx.y) / scale;
  return gViewCenter + dLy;
}

vec4 shaderMain(vec2 uv) {
  vec2 w = toWorld(uv);
  float r = length(w);

  // Pixel footprint in world units (ly per pixel). Used for dynamic LOD decisions.
  vec2 dw_dx = dFdx(w);
  vec2 dw_dy = dFdy(w);
  float pixelLy = max(length(dw_dx), length(dw_dy));

  // LOD thresholds are expressed in *pixels per feature*.
  float lodLo = min(gLodPxLo, gLodPxHi);
  float lodHi = max(gLodPxLo, gLodPxHi);
  lodHi = max(lodHi, lodLo + 1.0e-3);

  float inside = step(r, max(1.0, gRadiusLy));
  float radial = exp(-r / max(1.0, gRadialScaleLy));

  // Local disc midplane warp + flared half-thickness (matches CPU morphology helpers).
  float halfT = thicknessHalf(r);
  float warpZ = warpZLy(w, r);

  // Fraction of the truncated vertical density profile captured by the current Z slice.
  float verticalFrac = sliceVerticalFrac(warpZ, halfT);

  float mean = gBaseMean * radial * verticalFrac;
  mean *= inside;

  float rRef = max(1.0, gRadiusLy * 0.02);
  float arm = spiralMask(w, rRef);
  mean *= (1.0 + gSpiralArmStrength * arm);

  // Optional clumpy noise.
  if (gDensityNoiseStrength > 0.0 && gDensityNoiseFreq > 0.0) {
    vec2 np = w * gDensityNoiseFreq + vec2(17.13, 3.91);
    float n = gLodEnabled ? fbm2_lod(np, gLodPxLo, gLodPxHi, gLodEnergy) : fbm2(np);
    n = sat(n);
    float mul = mix(1.0 - gDensityNoiseStrength, 1.0 + gDensityNoiseStrength, n);
    mean *= mul;
  }

  // Cluster / void influence fields. Note: these are GPU-side approximations that are
  // deterministic and streaming-safe, but not guaranteed to match the CPU generator.

  float clusterLod = 1.0;
  float voidLod = 1.0;
  if (gLodEnabled) {
    clusterLod = smoothstep(lodLo, lodHi, max(0.0, gClusterRadiusLy) / max(1.0e-6, pixelLy));
    voidLod = smoothstep(lodLo, lodHi, max(0.0, gVoidRadiusLy) / max(1.0e-6, pixelLy));
  }

  float cluster01 = 0.0;
  if (gClusterStrength > 0.0) {
    bool doClusters = (!gLodEnabled) || (!gLodSkipTinyBlobs) || (clusterLod > 0.001);
    if (doClusters) {
      cluster01 = blobInfluence(vec3(w, gViewCenterZLy),
                               max(1.0, gClusterCellSizeLy),
                               gClusterChancePerCell,
                               max(1.0, gClusterRadiusLy),
                               gClusterRadiusJitter,
                               gClusterStrengthJitter,
                               gClusterFalloffPower,
                               1337);
      if (gLodEnabled) cluster01 *= clusterLod;
    }
  }

  float void01 = 0.0;
  if (gVoidStrength > 0.0) {
    bool doVoids = (!gLodEnabled) || (!gLodSkipTinyBlobs) || (voidLod > 0.001);
    if (doVoids) {
      void01 = blobInfluence(vec3(w, gViewCenterZLy),
                            max(1.0, gVoidCellSizeLy),
                            gVoidChancePerCell,
                            max(1.0, gVoidRadiusLy),
                            gVoidRadiusJitter,
                            gVoidStrengthJitter,
                            gVoidFalloffPower,
                            4242);
      if (gLodEnabled) void01 *= voidLod;
    }
  }

  mean *= (1.0 + gClusterStrength * cluster01);
  mean *= sat(1.0 - gVoidStrength * void01);

  float bar01 = barMask(w);
  float ring01 = ringMask(r);

  mean *= (1.0 + gBarStrength * bar01) * (1.0 + gRingStrength * ring01);

  float v = 0.0;
  if (gMode == 0) {
    // Density (tone-mapped).
    v = 1.0 - exp(-mean * max(0.001, gExposure));
  } else if (gMode == 1) {
    v = arm;
  } else if (gMode == 2) {
    v = cluster01;
  } else if (gMode == 3) {
    v = void01;
  } else if (gMode == 4) {
    v = max(bar01, ring01);
  } else {
    v = 1.0 - exp(-mean * max(0.001, gExposure));
  }

  v = pow(sat(v), 1.0 / max(0.05, gGamma));

  vec3 col;
  if (gMode == 3) {
    // Voids: cooler palette.
    col = palette(v,
                  vec3(0.04, 0.05, 0.07),
                  vec3(0.35, 0.40, 0.55),
                  vec3(1.0, 1.0, 1.0),
                  vec3(0.00, 0.15, 0.35));
  } else if (gMode == 2) {
    // Clusters: warm highlights.
    col = palette(v,
                  vec3(0.04, 0.05, 0.07),
                  vec3(0.55, 0.35, 0.20),
                  vec3(1.0, 1.0, 1.0),
                  vec3(0.10, 0.18, 0.32));
  } else {
    // Default density palette.
    col = palette(v,
                  vec3(0.03, 0.04, 0.06),
                  vec3(0.40, 0.40, 0.45),
                  vec3(1.0, 1.0, 1.0),
                  vec3(0.10, 0.20, 0.40));
  }

  // Contour lines (derivative-aware anti-aliasing via fwidth()).
  if (gContours) {
    float cc = float(max(2, gContourCount));
    float x = v * cc;
    float f = fract(x);
    float d = min(f, 1.0 - f);

    // fwidth(x) is roughly the change in x across one pixel.
    // Ensure the smoothing width never goes below that, preventing shimmer.
    float dv = fwidth(x);
    float w = max(max(1.0e-4, gContourWidth), 0.5 * dv);

    float line = 1.0 - smoothstep(0.0, w, d);
    col = mix(col, col + vec3(0.16), line);
  }

  // Slight vignette.
  vec2 q = (uv - vec2(0.5)) * 2.0;
  float vign = 1.0 - 0.22 * dot(q, q);
  col *= vign;

  // Fade outside disc.
  float edge = sat((max(1.0, gRadiusLy) - r) / max(1.0, max(1.0, gRadiusLy) * 0.05));
  col *= 0.30 + 0.70 * (inside * edge);

  return vec4(col, 1.0);
}

)GLSL";

static bool ensureHeatmapInit(ProceduralGalaxyLabWindowState& st) {
  if (st.heatmapInited) return st.heatmapError.empty();

  // Schema + shader compilation.
  if (st.heatmapParams.empty()) {
    buildHeatmapParamSchema(st.heatmapParams);
    st.heatmapHandles.inited = false;
  }

  // Cache param indices once and fail fast if the schema ever drifts.
  if (!ensureHeatmapParamHandles(st)) {
    st.heatmapInited = true;
    return false;
  }

  std::string err;
  if (!st.heatmapToy.init(&err)) {
    st.heatmapError = "Heatmap init failed: " + err;
    st.heatmapInited = true;
    return false;
  }

  const std::string header = st.heatmapParams.buildUniformDecls();
  if (!st.heatmapToy.buildWithHeader(kGalaxyHeatmapShader, header, &err)) {
    st.heatmapError = "Heatmap shader build failed: " + err;
    st.heatmapInited = true;
    return false;
  }

  // Dummy RT init (will be resized to match the canvas).
  if (!st.heatmapTarget.init(64, 64, &err)) {
    st.heatmapError = "Heatmap target init failed: " + err;
    st.heatmapInited = true;
    return false;
  }

  st.heatmapError.clear();
  st.heatmapFrame = 0;
  st.heatmapLastHash = core::u64(-1);
  st.heatmapLastW = 0;
  st.heatmapLastH = 0;
  st.heatmapInited = true;
  return true;
}

static void renderGalaxyHeatmap(ProceduralGalaxyLabWindowState& st, int canvasW, int canvasH, float timeSec) {
  if (!st.heatmapEnabled) return;
  if (canvasW <= 8 || canvasH <= 8) return;

  if (!ensureHeatmapInit(st)) return;
  if (!st.heatmapError.empty()) return;

  int div = st.heatmapResolutionDiv;
  if (div <= 0) div = 1;
  // Clamp to {1,2,4}.
  if (div != 1 && div != 2 && div != 4) {
    div = (div < 2) ? 1 : (div < 4 ? 2 : 4);
    st.heatmapResolutionDiv = div;
  }

  const int w = std::clamp(canvasW / div, 64, 2048);
  const int h = std::clamp(canvasH / div, 64, 2048);

  if (!st.heatmapTarget.isInited() || st.heatmapTarget.width() != w || st.heatmapTarget.height() != h) {
    std::string err;
    if (!st.heatmapTarget.init(w, h, &err)) {
      st.heatmapError = "Heatmap resize failed: " + err;
      return;
    }
  }

  const core::u64 key = heatmapRenderKey(st, w, h);
  if (!st.heatmapAnimate && st.heatmapFrame > 0 &&
      key == st.heatmapLastHash && st.heatmapLastW == w && st.heatmapLastH == h) {
    // Nothing affecting the heatmap has changed; reuse the previous texture.
    return;
  }

  // Save/restore a wider slice of GL state so this offscreen pass can't
  // accidentally perturb the main renderer (texture bindings, clear color, etc.).
  struct GlHeatmapStateGuard {
    GLint prevFbo{0};
    GLint prevViewport[4]{0, 0, 0, 0};
    GLint prevProg{0};
    GLint prevVao{0};
    GLint prevActiveTex{GL_TEXTURE0};
    GLint prevTex2D[4]{0, 0, 0, 0};
    GLfloat prevClearColor[4]{0.0f, 0.0f, 0.0f, 0.0f};
    GLboolean prevDepth{GL_FALSE};
    GLboolean prevCull{GL_FALSE};
    GLboolean prevBlend{GL_FALSE};
    GLboolean prevSciss{GL_FALSE};

    GlHeatmapStateGuard() {
      glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
      glGetIntegerv(GL_VIEWPORT, prevViewport);
      glGetIntegerv(GL_CURRENT_PROGRAM, &prevProg);
      glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &prevVao);
      glGetIntegerv(GL_ACTIVE_TEXTURE, &prevActiveTex);
      glGetFloatv(GL_COLOR_CLEAR_VALUE, prevClearColor);

      prevDepth = glIsEnabled(GL_DEPTH_TEST);
      prevCull  = glIsEnabled(GL_CULL_FACE);
      prevBlend = glIsEnabled(GL_BLEND);
      prevSciss = glIsEnabled(GL_SCISSOR_TEST);

      // Capture the 2D texture bindings for the texture units we touch (0..3).
      // NOTE: Do NOT call ::glActiveTexture directly (not exported on Windows).
      if (stellar::render::gl::ActiveTexture) {
        for (int i = 0; i < 4; ++i) {
          stellar::render::gl::ActiveTexture((GLenum)(GL_TEXTURE0 + i));
          glGetIntegerv(GL_TEXTURE_BINDING_2D, &prevTex2D[i]);
        }
        stellar::render::gl::ActiveTexture((GLenum)prevActiveTex);
      }
    }

    ~GlHeatmapStateGuard() {
      // Restore the texture bindings we modified.
      if (stellar::render::gl::ActiveTexture) {
        for (int i = 0; i < 4; ++i) {
          stellar::render::gl::ActiveTexture((GLenum)(GL_TEXTURE0 + i));
          glBindTexture(GL_TEXTURE_2D, (GLuint)prevTex2D[i]);
        }
        stellar::render::gl::ActiveTexture((GLenum)prevActiveTex);
      }

      if (prevDepth) glEnable(GL_DEPTH_TEST); else glDisable(GL_DEPTH_TEST);
      if (prevCull) glEnable(GL_CULL_FACE); else glDisable(GL_CULL_FACE);
      if (prevBlend) glEnable(GL_BLEND); else glDisable(GL_BLEND);
      if (prevSciss) glEnable(GL_SCISSOR_TEST); else glDisable(GL_SCISSOR_TEST);

      glUseProgram((GLuint)prevProg);
      glBindVertexArray((GLuint)prevVao);
      glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
      glViewport(prevViewport[0], prevViewport[1], prevViewport[2], prevViewport[3]);
      glClearColor(prevClearColor[0], prevClearColor[1], prevClearColor[2], prevClearColor[3]);
    }
  } guard;

  // Bind RT and render.
  st.heatmapTarget.begin();
  glViewport(0, 0, w, h);

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_CULL_FACE);
  glDisable(GL_BLEND);
  glDisable(GL_SCISSOR_TEST);

  glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  syncHeatmapParamsFromState(st);

  stellar::render::ShaderToyInputs in{};
  in.width = w;
  in.height = h;
  in.timeSec = timeSec;
  in.frame = st.heatmapFrame++;
  in.seed = seed32(st.seed);
  in.userParams = &st.heatmapParams;

  st.heatmapToy.draw(in);

  st.heatmapTarget.end();

  // Update cache key after a successful render.
  st.heatmapLastHash = key;
  st.heatmapLastW = w;
  st.heatmapLastH = h;
}


static ImU32 rgba(float r, float g, float b, float a = 1.0f) {
  return ImGui::GetColorU32(ImVec4(r, g, b, a));
}

static ImU32 colorForStarClass(sim::StarClass c) {
  // Rough temperature-ish palette; intentionally stylized.
  switch (c) {
    case sim::StarClass::O: return rgba(0.55f, 0.70f, 1.00f);
    case sim::StarClass::B: return rgba(0.65f, 0.75f, 1.00f);
    case sim::StarClass::A: return rgba(0.80f, 0.85f, 1.00f);
    case sim::StarClass::F: return rgba(1.00f, 0.95f, 0.80f);
    case sim::StarClass::G: return rgba(1.00f, 0.90f, 0.60f);
    case sim::StarClass::K: return rgba(1.00f, 0.75f, 0.45f);
    case sim::StarClass::M: return rgba(1.00f, 0.55f, 0.40f);
    default: return rgba(1.0f, 1.0f, 1.0f);
  }
}

static ImU32 colorForFaction(core::u32 id, const std::vector<sim::Faction>& factions) {
  if (id == 0) return rgba(0.75f, 0.75f, 0.75f);
  for (const auto& f : factions) {
    if (f.id == id) {
      return rgba(static_cast<float>(f.colorRgb.x),
                  static_cast<float>(f.colorRgb.y),
                  static_cast<float>(f.colorRgb.z));
    }
  }
  return rgba(0.75f, 0.75f, 0.75f);
}

static ImU32 colorForRegionKind(proc::GalaxyRegionKind k) {
  // Intentionally stylized palette (not physically based).
  switch (k) {
    case proc::GalaxyRegionKind::Core: return rgba(1.00f, 0.86f, 0.40f);
    case proc::GalaxyRegionKind::InnerDisc: return rgba(0.55f, 0.90f, 0.55f);
    case proc::GalaxyRegionKind::OuterRim: return rgba(0.55f, 0.70f, 1.00f);
    case proc::GalaxyRegionKind::Nebula: return rgba(0.95f, 0.55f, 0.95f);
    case proc::GalaxyRegionKind::Cluster: return rgba(0.85f, 0.85f, 0.85f);
    case proc::GalaxyRegionKind::Rift: return rgba(0.90f, 0.30f, 0.30f);
    default: return rgba(0.75f, 0.75f, 0.75f);
  }
}

 static ImU32 colorForClusterInfluence(float t) {
  t = std::clamp(t, 0.0f, 1.0f);
  // Simple ramp from dark -> bright; intentionally stylized.
  const float r = 0.10f + 0.90f * t;
  const float g = 0.20f + 0.80f * t;
  const float b = 0.35f + 0.65f * t;
  return rgba(r, g, b);
}

static ImU32 colorForVoidInfluence(float t) {
  t = std::clamp(t, 0.0f, 1.0f);
  // Highlight void "carving" with a distinct purple-ish ramp.
  // 0 -> muted grey, 1 -> vivid magenta.
  const float r = 0.25f + 0.65f * t;
  const float g = 0.25f + 0.20f * (1.0f - t);
  const float b = 0.30f + 0.70f * t;
  return rgba(r, g, b);
}

static ImU32 colorForMorphologyInfluence(float t) {
  t = std::clamp(t, 0.0f, 1.0f);
  // Warm ramp: low -> cool; high -> hot.
  const float r = 0.20f + 0.80f * t;
  const float g = 0.18f + 0.52f * (1.0f - std::abs(2.0f * t - 1.0f));
  const float b = 0.12f + 0.70f * (1.0f - t);
  return rgba(r, g, b);
}

static bool dragDouble(const char* label, double& v, double step, double minV, double maxV, const char* fmt = "%.3f") {
  double tmp = v;
  const bool changed = ImGui::DragScalar(label, ImGuiDataType_Double, &tmp, step, &minV, &maxV, fmt);
  if (changed) v = tmp;
  return changed;
}

static void rebuildPreview(ProceduralGalaxyLabWindowState& st) {
  const auto t0 = Clock::now();

  const int fc = std::max(1, st.factionCount);
  st.factions = sim::generateFactions(st.seed, static_cast<core::u32>(fc));

  // Reset diagnostics per regenerate.
  st.previewSectorCacheHits = 0;
  st.previewSectorCacheMisses = 0;
  st.previewSectorCacheEvictions = 0;
  st.previewParallelWorkersUsed = 1;

  // Preview sector cache is keyed by generator *context*.
  // Any changes to (seed, params, faction count) invalidate cached sectors.
  const core::u64 ctx = hashGalaxyParams(hashMix(hashMix(core::fnv1a64("galaxy_preview_cache/v1"), st.seed), fc), st.params);
  if (ctx != st.previewSectorCacheCtxHash) {
    st.previewSectorCache.clear();
    st.previewSectorCacheStubTotal = 0;
    st.previewSectorCacheCtxHash = ctx;
  }

  // Apply capacity / enable toggles.
  if (!st.previewUseSectorCache) {
    if (!st.previewSectorCache.empty()) {
      st.previewSectorCache.clear();
      st.previewSectorCacheStubTotal = 0;
    }
  } else {
    const std::size_t cap = st.previewSectorCacheMaxEntries < 0 ? 0u : static_cast<std::size_t>(st.previewSectorCacheMaxEntries);
    auto onEvict = [&](const ProceduralGalaxyLabWindowState::SectorCoord&, ProceduralGalaxyLabWindowState::SectorStubVecPtr& v) {
      const std::size_t n = v ? v->size() : 0u;
      if (st.previewSectorCacheStubTotal >= n) st.previewSectorCacheStubTotal -= n;
      else st.previewSectorCacheStubTotal = 0;
      if (st.previewSectorCacheEvictions < 0x7FFFFFFF) ++st.previewSectorCacheEvictions;
    };

    if (cap != st.previewSectorCache.capacity()) {
      st.previewSectorCache.setCapacity(cap, onEvict);
      if (cap > 0) st.previewSectorCache.reserve(std::min<std::size_t>(cap + cap / 4 + 64, 400000u));
    } else {
      st.previewSectorCache.pruneToCapacity(onEvict);
    }
  }

  proc::GalaxyGenerator gen(st.seed, st.params);

  st.stubs.clear();
  st.stubs.reserve(4096);

  // Clear cached per-stub metadata; rebuild only what is needed for the current visualization mode.
  st.stubRegionKind.clear();
  st.stubRegionId.clear();
  if (st.colorByRegion) {
    st.stubRegionKind.reserve(4096);
    st.stubRegionId.reserve(4096);
  }

  st.stubCluster01.clear();
  if (st.colorByCluster) st.stubCluster01.reserve(4096);

  st.stubVoid01.clear();
  if (st.colorByVoid) st.stubVoid01.reserve(4096);

  st.stubMorph01.clear();
  if (st.colorByMorphology) st.stubMorph01.reserve(4096);

  // Preview sampling diagnostics.
  st.previewSectorStrideXY = 1;
  st.previewSectorsGenerated = 0;
  st.previewCandidateStubs = 0;

  // Recompute a stable hash for the currently materialized preview stub set.
  // This is used to invalidate dependent overlays (hyperlanes, etc) without
  // forcing a full preview regeneration.
  auto updateStubSetHash = [&]() {
    core::u64 h = core::fnv1a64("galaxy_preview/stubset/v1");
    h = core::hashCombine(h, st.seed);
    h = core::hashCombine(h, (core::u64)st.stubs.size());
    for (const auto& s0 : st.stubs) {
      h = core::hashCombine(h, (core::u64)s0.id);
    }
    st.previewStubSetHash = h;
    st.hyperlanesDirty = true;
  };

  // Early out: base mean <= 0 implies no systems anywhere (all other modifiers only scale this mean).
  if (st.params.baseMeanSystemsPerSector <= 0.0) {
    updateStubSetHash();


    st.lastGenMs = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    st.dirty = false;
    return;
  }

  proc::GalaxyClustersParams clusterParams{};
  if (st.colorByCluster) {
    clusterParams.cellSizeLy = st.params.clusterCellSizeLy;
    clusterParams.chancePerCell = st.params.clusterChancePerCell;
    clusterParams.radiusLy = st.params.clusterRadiusLy;
    clusterParams.radiusJitter01 = st.params.clusterRadiusJitter;
    clusterParams.strengthJitter01 = st.params.clusterStrengthJitter;
    clusterParams.falloffPower = st.params.clusterFalloffPower;
  }

  proc::GalaxyVoidsParams voidParams{};
  if (st.colorByVoid) {
    voidParams.cellSizeLy = st.params.voidCellSizeLy;
    voidParams.chancePerCell = st.params.voidChancePerCell;
    voidParams.radiusLy = st.params.voidRadiusLy;
    voidParams.radiusJitter01 = st.params.voidRadiusJitter;
    voidParams.strengthJitter01 = st.params.voidStrengthJitter;
    voidParams.falloffPower = st.params.voidFalloffPower;
  }

  const double s = std::max(1.0, st.params.sectorSizeLy);
  const double r = std::max(1.0, st.viewRadiusLy);
  const double zHalf = std::max(0.0, st.zHalfLy);

  const math::Vec3d c = st.centerLy;

  const int minX = static_cast<int>(std::floor((c.x - r) / s));
  const int maxX = static_cast<int>(std::floor((c.x + r) / s));
  const int minY = static_cast<int>(std::floor((c.y - r) / s));
  const int maxY = static_cast<int>(std::floor((c.y + r) / s));
  const int minZ = static_cast<int>(std::floor((c.z - zHalf) / s));
  const int maxZ = static_cast<int>(std::floor((c.z + zHalf) / s));

  const int nX = std::max(1, maxX - minX + 1);
  const int nY = std::max(1, maxY - minY + 1);
  const int nZ = std::max(1, maxZ - minZ + 1);

  const double r2 = r * r;

  struct PreviewItem {
    core::u64 key{0};
    sim::SystemStub stub{};

    proc::GalaxyRegionKind regionKind{proc::GalaxyRegionKind::OuterRim};
    core::u64 regionId{0};

    float cluster01{0.0f};
    float void01{0.0f};
    float morph01{0.0f};
  };

  // Max-heap by key: top() is the *worst* item (largest key) so we can pop/replace.
  struct PreviewItemLess {
    bool operator()(const PreviewItem& a, const PreviewItem& b) const { return a.key < b.key; }
  };

  const std::size_t K = std::max<std::size_t>(1, st.maxStubs);

  // Keep track of sectors already generated (multi-pass sampling may revisit blocks).
  std::unordered_set<core::u64> visited;
  visited.reserve(static_cast<std::size_t>(std::min(400000.0, (double)nX * (double)nY * (double)nZ)));

  std::vector<ProceduralGalaxyLabWindowState::SectorCoord> sectors;

  auto axisDist = [](double v, double lo, double hi) {
    if (v < lo) return lo - v;
    if (v > hi) return v - hi;
    return 0.0;
  };

  auto floorDiv = [](int a, int b) {
    // b>0. C++ truncates toward 0; for negative a we need floor division.
    int q = a / b;
    int r = a % b;
    if (r != 0 && a < 0) --q;
    return q;
  };

  auto sectorVisitedKey = [&](int x, int y, int z) -> core::u64 {
    struct C {
      int x;
      int y;
      int z;
    };
    const C c0{x, y, z};
    return core::hashCombine(core::fnv1a64("galaxy_preview/sector"),
                             core::hashCombine(st.seed, core::hashBytes(&c0, sizeof(C))));
  };

  auto queueSector = [&](int x, int y, int z) {
    const core::u64 vk = sectorVisitedKey(x, y, z);
    if (!visited.insert(vk).second) return;
    sectors.push_back({x, y, z});
  };

  // Always sample a small neighborhood around the view center at full resolution.
  const int cx = static_cast<int>(std::floor(c.x / s));
  const int cy = static_cast<int>(std::floor(c.y / s));
  for (int z = minZ; z <= maxZ; ++z) {
    for (int y = cy - 1; y <= cy + 1; ++y) {
      if (y < minY || y > maxY) continue;
      for (int x = cx - 1; x <= cx + 1; ++x) {
        if (x < minX || x > maxX) continue;
        queueSector(x, y, z);
      }
    }
  }

  // Adaptive "preview LOD": sample ~O(maxStubs) sectors using a jittered grid.
  // This avoids scanning enormous bounding boxes when zoomed out, while keeping the sample
  // deterministic in world space (stable as you pan).
  const double mean = std::max(0.05, st.params.baseMeanSystemsPerSector);
  double targetSectorsD = (double)K / mean;
  targetSectorsD *= 1.75; // headroom for circle+z filtering
  targetSectorsD = std::clamp(targetSectorsD, 256.0, 200000.0);

  const double totalSectorsD = (double)nX * (double)nY * (double)nZ;
  const double ratio = totalSectorsD / targetSectorsD;

  int strideXY = 1;
  if (ratio > 1.0) strideXY = (int)std::ceil(std::sqrt(ratio));
  strideXY = std::clamp(strideXY, 1, std::max(nX, nY));
  st.previewSectorStrideXY = strideXY;

  const std::size_t targetSectors = static_cast<std::size_t>(targetSectorsD);
  const std::size_t reserveSectors = std::min<std::size_t>(
      std::max<std::size_t>(targetSectors + 256, 512),
      static_cast<std::size_t>(std::min<double>(400000.0, totalSectorsD)));
  sectors.reserve(reserveSectors);

  auto samplePass = [&](int stride, int passIndex) {
    if (stride <= 0) return;

    // Anchor blocks to the world grid to keep the sampling stable under panning.
    const int startX = floorDiv(minX, stride) * stride;
    const int startY = floorDiv(minY, stride) * stride;

    for (int z = minZ; z <= maxZ; ++z) {
      for (int by = startY; by <= maxY; by += stride) {
        const int y0 = std::max(minY, by);
        const int y1 = std::min(maxY, by + stride - 1);
        const int bh = y1 - y0 + 1;
        if (bh <= 0) continue;

        for (int bx = startX; bx <= maxX; bx += stride) {
          const int x0 = std::max(minX, bx);
          const int x1 = std::min(maxX, bx + stride - 1);
          const int bw = x1 - x0 + 1;
          if (bw <= 0) continue;

          // Skip blocks entirely outside the view circle.
          const double xLo = (double)x0 * s;
          const double xHi = (double)(x1 + 1) * s;
          const double yLo = (double)y0 * s;
          const double yHi = (double)(y1 + 1) * s;
          const double dx = axisDist(c.x, xLo, xHi);
          const double dy = axisDist(c.y, yLo, yHi);
          if (dx * dx + dy * dy > r2) continue;

          struct BlockKey {
            int bx;
            int by;
            int z;
            int stride;
            int pass;
          };
          const BlockKey bk{bx, by, z, stride, passIndex};
          const core::u64 h = core::hashCombine(st.seed, core::hashBytes(&bk, sizeof(BlockKey)));
          core::SplitMix64 rng(h);

          const int sx = x0 + (int)(rng.nextU32() % (core::u32)bw);
          const int sy = y0 + (int)(rng.nextU32() % (core::u32)bh);

          queueSector(sx, sy, z);
        }
      }
    }
  };

  // Multi-pass refinement: if we didn't sample enough sectors, halve the stride
  // a couple of times (higher quality) while keeping everything deterministic.
  for (int pass = 0; pass < 3; ++pass) {
    samplePass(strideXY, pass);
    if (strideXY <= 1) break;
    if (sectors.size() >= targetSectors) break;

    if (sectors.size() < (targetSectors * 3) / 4) {
      strideXY = std::max(1, strideXY / 2);
      st.previewSectorStrideXY = strideXY;
    } else {
      break;
    }
  }

  st.previewSectorsGenerated = (int)std::min<std::size_t>(sectors.size(), 0x7FFFFFFF);

  // -----------------
  // Process sampled sectors (optionally parallel)
  // -----------------

  auto stubSampleKey = [&](const sim::SystemStub& stub) -> core::u64 {
    return core::hashCombine(core::fnv1a64("galaxy_preview/stub"),
                             core::hashCombine(st.seed, (core::u64)stub.id));
  };

  struct WorkerResult {
    std::vector<PreviewItem> items;
    int candidateStubs{0};
    int cacheHits{0};
    int cacheMisses{0};
  };

  using StubVec = ProceduralGalaxyLabWindowState::SectorStubVec;
  using StubVecPtr = ProceduralGalaxyLabWindowState::SectorStubVecPtr;

  std::mutex cacheMx;

  auto cacheEvictCb = [&](const ProceduralGalaxyLabWindowState::SectorCoord&, StubVecPtr& v) {
    const std::size_t n = v ? v->size() : 0u;
    if (st.previewSectorCacheStubTotal >= n) st.previewSectorCacheStubTotal -= n;
    else st.previewSectorCacheStubTotal = 0;
    if (st.previewSectorCacheEvictions < 0x7FFFFFFF) ++st.previewSectorCacheEvictions;
  };

  const bool wantRegion = st.colorByRegion;
  const bool wantCluster = st.colorByCluster;
  const bool wantVoid = st.colorByVoid;
  const bool wantMorph = st.colorByMorphology;
  const double morphDenom = std::max(1.0e-9, std::abs(st.params.barStrength) + std::abs(st.params.ringStrength));

  auto workerFn = [&](std::size_t begin, std::size_t end) -> WorkerResult {
    WorkerResult out;

    std::priority_queue<PreviewItem, std::vector<PreviewItem>, PreviewItemLess> heap;

    auto maybeKeepStub = [&](const sim::SystemStub& stub) {
      if (std::abs(stub.posLy.z - c.z) > zHalf) return;

      const double dx = stub.posLy.x - c.x;
      const double dy = stub.posLy.y - c.y;
      if (dx * dx + dy * dy > r2) return;

      if (out.candidateStubs < 0x7FFFFFFF) ++out.candidateStubs;

      const core::u64 key = stubSampleKey(stub);
      if (heap.size() >= K && key >= heap.top().key) return;

      PreviewItem item;
      item.key = key;
      item.stub = stub;

      if (wantRegion) {
        const auto reg = proc::sampleGalaxyRegion(st.seed, stub.posLy, st.regionCellSizeLy);
        item.regionKind = reg.kind;
        item.regionId = reg.regionId;
      }

      if (wantCluster) {
        const auto cs = proc::sampleGalaxyClusters(st.seed, stub.posLy, clusterParams);
        item.cluster01 = static_cast<float>(cs.cluster01);
      }

      if (wantVoid) {
        const auto vs = proc::sampleGalaxyVoids(st.seed, stub.posLy, voidParams);
        item.void01 = static_cast<float>(vs.void01);
      }

      if (wantMorph) {
        const auto ms = proc::sampleGalaxyMorphology(st.seed, st.params, stub.posLy);
        double t = 0.0;
        if (morphDenom > 0.0) t = (ms.densityMul - 1.0) / morphDenom;
        t = std::clamp(t, 0.0, 1.0);
        item.morph01 = static_cast<float>(t);
      }

      if (heap.size() < K) {
        heap.push(std::move(item));
      } else {
        heap.pop();
        heap.push(std::move(item));
      }
    };

    for (std::size_t i = begin; i < end; ++i) {
      const auto sc = sectors[i];

      if (st.previewUseSectorCache) {
        StubVecPtr systems{};
        {
          std::lock_guard<std::mutex> lock(cacheMx);
          if (auto* p = st.previewSectorCache.find(sc)) {
            systems = *p;
          }
        }

        if (systems) {
          if (out.cacheHits < 0x7FFFFFFF) ++out.cacheHits;
        } else {
          proc::Sector sec = gen.generateSector({sc.x, sc.y, sc.z}, st.factions);
          const auto sp = std::make_shared<StubVec>(std::move(sec.systems));
          StubVecPtr ins = sp;

          {
            std::lock_guard<std::mutex> lock(cacheMx);
            bool inserted = false;
            auto& ref = st.previewSectorCache.getOrInsert(sc, std::move(ins), &inserted, cacheEvictCb);
            systems = ref;
            if (inserted) {
              const std::size_t n = systems ? systems->size() : 0u;
              st.previewSectorCacheStubTotal += n;
              if (out.cacheMisses < 0x7FFFFFFF) ++out.cacheMisses;
            } else {
              if (out.cacheHits < 0x7FFFFFFF) ++out.cacheHits;
            }
          }
        }

        if (systems) {
          for (const auto& stub : *systems) {
            maybeKeepStub(stub);
          }
        }
      } else {
        proc::Sector sec = gen.generateSector({sc.x, sc.y, sc.z}, st.factions);
        for (const auto& stub : sec.systems) {
          maybeKeepStub(stub);
        }
      }
    }

    out.items.reserve(heap.size());
    while (!heap.empty()) {
      out.items.push_back(heap.top());
      heap.pop();
    }
    return out;
  };

  std::vector<WorkerResult> results;
  results.reserve(8);

  const bool wantParallel = st.previewParallel && st.previewParallelThreads != 1 && sectors.size() >= 64;

  if (!wantParallel) {
    results.push_back(workerFn(0, sectors.size()));
    st.previewParallelWorkersUsed = 1;
  } else {
    int poolThreads = 1;
    core::JobSystem& jobs = previewJobs(st.previewParallelThreads, &poolThreads);
    st.previewParallelWorkersUsed = std::max(1, poolThreads);

    const std::size_t pool = static_cast<std::size_t>(std::max(1, poolThreads));
    const std::size_t taskCount = std::max<std::size_t>(1, std::min<std::size_t>(sectors.size(), pool * 4));
    const std::size_t chunkSize = (sectors.size() + taskCount - 1) / taskCount;

    std::vector<std::future<WorkerResult>> futs;
    futs.reserve(taskCount);

    for (std::size_t t = 0; t < taskCount; ++t) {
      const std::size_t b = t * chunkSize;
      if (b >= sectors.size()) break;
      const std::size_t e = std::min<std::size_t>(sectors.size(), b + chunkSize);
      futs.push_back(jobs.submit(workerFn, b, e));
    }

    results.reserve(futs.size());
    for (auto& f : futs) {
      results.push_back(f.get());
    }
  }

  // Merge: take the global bottom-K by key across worker local bottom-K sets.
  std::priority_queue<PreviewItem, std::vector<PreviewItem>, PreviewItemLess> global;

  int totalCandidates = 0;
  int totalHits = 0;
  int totalMisses = 0;

  for (auto& r0 : results) {
    totalCandidates = std::min(0x7FFFFFFF, totalCandidates + r0.candidateStubs);
    totalHits = std::min(0x7FFFFFFF, totalHits + r0.cacheHits);
    totalMisses = std::min(0x7FFFFFFF, totalMisses + r0.cacheMisses);

    for (auto& it : r0.items) {
      if (global.size() < K) {
        global.push(std::move(it));
      } else if (it.key < global.top().key) {
        global.pop();
        global.push(std::move(it));
      }
    }
  }

  st.previewCandidateStubs = totalCandidates;
  st.previewSectorCacheHits = totalHits;
  st.previewSectorCacheMisses = totalMisses;

  // Materialize sample set in a stable order.
  std::vector<PreviewItem> items;
  items.reserve(global.size());
  while (!global.empty()) {
    items.push_back(global.top());
    global.pop();
  }
  std::sort(items.begin(), items.end(), [](const PreviewItem& a, const PreviewItem& b) { return a.key < b.key; });

  st.stubs.clear();
  st.stubs.reserve(items.size());

  if (st.colorByRegion) {
    st.stubRegionKind.clear();
    st.stubRegionId.clear();
    st.stubRegionKind.reserve(items.size());
    st.stubRegionId.reserve(items.size());
  }
  if (st.colorByCluster) {
    st.stubCluster01.clear();
    st.stubCluster01.reserve(items.size());
  }
  if (st.colorByVoid) {
    st.stubVoid01.clear();
    st.stubVoid01.reserve(items.size());
  }
  if (st.colorByMorphology) {
    st.stubMorph01.clear();
    st.stubMorph01.reserve(items.size());
  }

  for (const auto& it : items) {
    st.stubs.push_back(it.stub);
    if (st.colorByRegion) {
      st.stubRegionKind.push_back(it.regionKind);
      st.stubRegionId.push_back(it.regionId);
    }
    if (st.colorByCluster) st.stubCluster01.push_back(it.cluster01);
    if (st.colorByVoid) st.stubVoid01.push_back(it.void01);
    if (st.colorByMorphology) st.stubMorph01.push_back(it.morph01);
  }

  updateStubSetHash();
  st.lastGenMs = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
  st.dirty = false;
}

// -----------------------------
// Hyperlane overlay (procedural navigation graph)
// -----------------------------

static void rebuildHyperlanes(ProceduralGalaxyLabWindowState& st) {
  if (!st.showHyperlanes) {
    // Keep cached data around, but don't spend time rebuilding when disabled.
    st.hyperlanesDirty = false;
    return;
  }

  if (st.stubs.size() < 2) {
    st.hyperlaneNet.edges.clear();
    st.hyperlaneNodes.clear();
    st.hyperlaneNodesUsed = 0;
    st.hyperlaneEdgesUsed = 0;
    st.hyperlaneGenMs = 0.0;
    st.hyperlanesDirty = false;
    st.hyperlaneLastHash = core::u64(-1);

    // Reset traffic / chokepoint analytics cache.
    st.hyperlaneCentralityDirty = false;
    st.hyperlaneCentralityLastHash = core::u64(-1);
    st.hyperlaneEdgeBetweenness.clear();
    st.hyperlaneNodeBetweenness.clear();
    st.hyperlaneCentralityNodeIds.clear();
    st.hyperlaneCentralityMaxEdge = 0.0;
    st.hyperlaneCentralityMaxNode = 0.0;
    st.hyperlaneCentralityMs = 0.0;
    st.hyperlaneChokepointEdgeIdx.clear();

    // Reset structural criticality analytics cache.
    st.hyperlaneVulnerabilityDirty = false;
    st.hyperlaneVulnerabilityLastHash = core::u64(-1);
    st.hyperlaneVulnerability = proc::HyperlaneVulnerabilityResult{};
    st.hyperlaneVulnerabilityMs = 0.0;
    st.hyperlaneCriticalBridgeEdgeIdx.clear();
    st.hyperlaneCriticalArticulationNodeIds.clear();
    st.hyperlaneCriticalArticulationNodeImpact01.clear();

    return;
  }

  std::size_t maxNodes = st.hyperlaneMaxNodes;
  if (maxNodes == 0) maxNodes = st.stubs.size();
  maxNodes = std::clamp<std::size_t>(maxNodes, 2, st.stubs.size());

  // Hash all generation inputs. If the hash is unchanged, we can reuse the cached
  // network without recomputing expensive O(n log n) topology.
  core::u64 h = core::fnv1a64("galaxy_preview/hyperlanes/v1");
  h = core::hashCombine(h, st.seed);
  h = core::hashCombine(h, st.previewStubSetHash);
  h = hashHyperlaneParams(h, st.hyperlaneParams);
  h = hashMix(h, (core::u64)maxNodes);

  // Edge generation cap: tie it loosely to the draw cap so the cached network is
  // bounded even if the neighborhood is dense.
  std::size_t genCap = st.hyperlaneParams.maxEdges;
  if (genCap == 0 && st.hyperlaneMaxEdgesDraw > 0) {
    genCap = std::min<std::size_t>(std::max<std::size_t>(1024, st.hyperlaneMaxEdgesDraw * 4), 250000);
  }
  h = hashMix(h, (core::u64)genCap);

  if (!st.hyperlanesDirty && h == st.hyperlaneLastHash) return;

  const auto t0 = Clock::now();

  // Deterministic node subset selection (bottom-K by stable hash key).
  st.hyperlaneNodes.clear();
  st.hyperlaneNodes.reserve(maxNodes);

  if (st.stubs.size() <= maxNodes) {
    st.hyperlaneNodes = st.stubs;
  } else {
    struct Cand {
      core::u64 key{0};
      const sim::SystemStub* stub{nullptr};
    };
    struct CandLess {
      // For a max-heap: top() has the *largest* key (worst of the bottom-K).
      bool operator()(const Cand& a, const Cand& b) const { return a.key < b.key; }
    };

    std::priority_queue<Cand, std::vector<Cand>, CandLess> heap;
    heap = {};

    const core::u64 base = core::hashCombine(st.seed, core::fnv1a64("galaxy_preview/hyperlane_nodes/v1"));

    for (const auto& s : st.stubs) {
      const core::u64 key = core::hashCombine(base, (core::u64)s.id);
      Cand c{key, &s};
      if (heap.size() < maxNodes) {
        heap.push(c);
      } else if (c.key < heap.top().key) {
        heap.pop();
        heap.push(c);
      }
    }

    std::vector<Cand> picks;
    picks.reserve(heap.size());
    while (!heap.empty()) {
      picks.push_back(heap.top());
      heap.pop();
    }

    std::sort(picks.begin(), picks.end(), [](const Cand& a, const Cand& b) {
      if (a.key != b.key) return a.key < b.key;
      return a.stub->id < b.stub->id;
    });

    for (const auto& p : picks) {
      if (p.stub) st.hyperlaneNodes.push_back(*p.stub);
    }
  }

  proc::HyperlaneParams hp = st.hyperlaneParams;
  // Use the derived generation cap (genCap), but never override an explicit user cap.
  if (hp.maxEdges == 0) hp.maxEdges = genCap;

  st.hyperlaneNet = proc::generateHyperlaneNetwork(st.seed, st.hyperlaneNodes, hp);

  // Topology changed; traffic analytics must be recomputed if/when requested.
  st.hyperlaneCentralityDirty = true;
  st.hyperlaneVulnerabilityDirty = true;

  st.hyperlaneNodesUsed = (int)st.hyperlaneNodes.size();
  st.hyperlaneEdgesUsed = (int)st.hyperlaneNet.edges.size();
  st.hyperlaneGenMs = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();

  st.hyperlaneLastHash = h;
  st.hyperlanesDirty = false;
}


// -----------------------------------------------------------------------------
// Hyperlane routing (pathfinding inspector)
// -----------------------------------------------------------------------------

static const sim::SystemStub* findStubById(const std::vector<sim::SystemStub>& stubs, sim::SystemId id) {
  if (id == 0) return nullptr;
  for (const auto& s : stubs) {
    if (s.id == id) return &s;
  }
  return nullptr;
}

static void ensureHyperlaneRouter(ProceduralGalaxyLabWindowState& st) {
  if (!st.showHyperlanes) {
    st.hyperlaneRouterInited = false;
    st.hyperlaneRouterHash = core::u64(-1);
    return;
  }

  if (st.hyperlaneLastHash == core::u64(-1) || st.hyperlaneNodes.empty() || st.hyperlaneNet.edges.empty()) {
    st.hyperlaneRouterInited = false;
    st.hyperlaneRouterHash = core::u64(-1);
    return;
  }

  if (st.hyperlaneRouterInited && st.hyperlaneRouterHash == st.hyperlaneLastHash) return;

  st.hyperlaneRouter.reset(st.hyperlaneNet, st.hyperlaneNodes);
  st.hyperlaneRouterInited = true;
  st.hyperlaneRouterHash = st.hyperlaneLastHash;

  // Any cached route is invalid once the topology changes.
  st.hyperlaneRouteDirty = true;
}

static void rebuildHyperlaneRoute(ProceduralGalaxyLabWindowState& st) {
  if (!st.showHyperlanes || !st.hyperlaneRouteEnabled) {
    st.hyperlaneRouteError.clear();
    st.hyperlaneRoutePath.clear();
    st.hyperlaneKRoutes.clear();
    st.hyperlaneRouteMetrics = proc::HyperlanePathMetrics{};
    st.hyperlaneRouteMs = 0.0;
    st.hyperlaneRouteDirty = false;
    st.hyperlaneRouteLastHash = core::u64(-1);
    return;
  }

  // Keep the router adjacency cache warm (useful for other tooling), and use its
  // hash invalidation to invalidate our cached K-routes when topology changes.
  ensureHyperlaneRouter(st);

  if (st.hyperlaneLastHash == core::u64(-1) || st.hyperlaneNodes.empty() || st.hyperlaneNet.edges.empty()) {
    st.hyperlaneRouteError = "Router not available (hyperlane graph empty).";
    st.hyperlaneRoutePath.clear();
    st.hyperlaneKRoutes.clear();
    st.hyperlaneRouteMetrics = proc::HyperlanePathMetrics{};
    st.hyperlaneRouteMs = 0.0;
    st.hyperlaneRouteDirty = false;
    return;
  }

  // Nothing selected yet.
  if (st.hyperlaneRouteFrom == 0 || st.hyperlaneRouteTo == 0) {
    st.hyperlaneRouteError.clear();
    st.hyperlaneRoutePath.clear();
    st.hyperlaneKRoutes.clear();
    st.hyperlaneRouteMetrics = proc::HyperlanePathMetrics{};
    st.hyperlaneRouteMs = 0.0;
    st.hyperlaneRouteDirty = false;
    return;
  }

  // Degenerate: start == goal.
  if (st.hyperlaneRouteFrom == st.hyperlaneRouteTo) {
    st.hyperlaneRouteError.clear();
    st.hyperlaneKRoutes.clear();
    st.hyperlaneRoutePath = {st.hyperlaneRouteFrom};
    st.hyperlaneRouteMetrics = proc::HyperlanePathMetrics{};
    st.hyperlaneRouteMetrics.reachable = true;
    st.hyperlaneRouteMetrics.hops = 0;
    st.hyperlaneRouteMetrics.costLy = 0.0;
    st.hyperlaneRouteMetrics.distanceLy = 0.0;
    st.hyperlaneRouteMetrics.risk01 = 0.0;
    st.hyperlaneRouteMetrics.bottleneckBandwidth01 = 1.0;
    st.hyperlaneRouteMs = 0.0;
    st.hyperlaneRouteDirty = false;
    return;
  }

  int k = st.hyperlaneRouteK;
  if (k < 1) k = 1;
  if (k > 12) k = 12;

  core::u64 h = core::fnv1a64("galaxy_preview/hyperlane_k_route/v1");
  h = hashMix(h, st.hyperlaneLastHash);
  h = hashMix(h, (core::u64)st.hyperlaneRouteFrom);
  h = hashMix(h, (core::u64)st.hyperlaneRouteTo);
  h = hashHyperlaneTravelParams(h, st.hyperlaneTravel);
  h = hashMix(h, (core::u64)k);

  const bool needSolve = st.hyperlaneRouteDirty || (h != st.hyperlaneRouteLastHash);

  if (needSolve) {
    const auto t0 = Clock::now();

    st.hyperlaneRouteError.clear();
    st.hyperlaneRoutePath.clear();
    st.hyperlaneRouteMetrics = proc::HyperlanePathMetrics{};
    st.hyperlaneKRoutes.clear();

    st.hyperlaneRouteSelected = 0;

    st.hyperlaneKRoutes = proc::plotKHyperlaneRoutesAStarCost(
      st.hyperlaneNodes,
      st.hyperlaneNet,
      st.hyperlaneRouteFrom,
      st.hyperlaneRouteTo,
      st.hyperlaneTravel,
      (std::size_t)k,
      /*maxExpansionsPerSolve=*/250000);

    st.hyperlaneRouteMs = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();

    if (st.hyperlaneKRoutes.empty()) {
      const bool hasStart = (findStubById(st.hyperlaneNodes, st.hyperlaneRouteFrom) != nullptr);
      const bool hasGoal  = (findStubById(st.hyperlaneNodes, st.hyperlaneRouteTo) != nullptr);
      if (!hasStart || !hasGoal) {
        st.hyperlaneRouteError = "Start/End system is not in the hyperlane node subset (increase Max Nodes).";
      } else {
        st.hyperlaneRouteError = "Destination unreachable within the current hyperlane subset.";
      }
    }

    st.hyperlaneRouteLastHash = h;
    st.hyperlaneRouteDirty = false;
  }

  // Even if we didn't re-solve, selection changes should update the drawn route.
  if (!st.hyperlaneKRoutes.empty()) {
    int sel = st.hyperlaneRouteSelected;
    if (sel < 0) sel = 0;
    if (sel >= (int)st.hyperlaneKRoutes.size()) sel = (int)st.hyperlaneKRoutes.size() - 1;
    st.hyperlaneRouteSelected = sel;

    st.hyperlaneRoutePath = st.hyperlaneKRoutes[(std::size_t)sel].path;
    st.hyperlaneRouteMetrics = st.hyperlaneKRoutes[(std::size_t)sel].metrics;
    if (!st.hyperlaneRouteMetrics.reachable) {
      st.hyperlaneRouteError = "Destination unreachable within the current hyperlane subset.";
    } else {
      st.hyperlaneRouteError.clear();
    }
  }
}

// -----------------------------------------------------------------------------
// Hyperlane traffic analytics (betweenness centrality)
// -----------------------------------------------------------------------------

static void rebuildHyperlaneCentrality(ProceduralGalaxyLabWindowState& st) {
  const bool wantTraffic = st.showHyperlanes && (st.hyperlaneColorMode == 2);
  if (!wantTraffic) {
    // Keep cached values around, but don't spend time recomputing unless the
    // visualization is actively requested.
    return;
  }

  if (st.hyperlaneLastHash == core::u64(-1) || st.hyperlaneNodes.empty() || st.hyperlaneNet.edges.empty()) {
    st.hyperlaneCentralityDirty = false;
    st.hyperlaneCentralityLastHash = core::u64(-1);
    st.hyperlaneEdgeBetweenness.clear();
    st.hyperlaneNodeBetweenness.clear();
    st.hyperlaneCentralityNodeIds.clear();
    st.hyperlaneCentralityMaxEdge = 0.0;
    st.hyperlaneCentralityMaxNode = 0.0;
    st.hyperlaneCentralityMs = 0.0;
    st.hyperlaneChokepointEdgeIdx.clear();
    return;
  }

  int samples = st.hyperlaneTrafficSamples;
  if (samples < 0) samples = 0;
  if (samples > 8192) samples = 8192;

  core::u64 h = core::fnv1a64("galaxy_preview/hyperlane_centrality/v1");
  h = hashMix(h, st.hyperlaneLastHash);
  h = hashHyperlaneTravelParams(h, st.hyperlaneTravel);
  h = hashMix(h, (core::u64)samples);
  h = hashMix(h, st.seed);

  if (!st.hyperlaneCentralityDirty && h == st.hyperlaneCentralityLastHash) return;

  const auto t0 = Clock::now();

  proc::HyperlaneBetweennessParams p{};
  p.travel = st.hyperlaneTravel;
  p.sampleSources = (std::size_t)samples;
  p.sampleSeed = core::hashCombine(st.seed, core::fnv1a64("galaxy_preview/hyperlane_centrality_seed/v1"));
  p.maxExpansionsPerSource = 0;
  p.tieEps = 1e-9;

  proc::HyperlaneBetweennessResult r = proc::estimateHyperlaneBetweennessCentrality(st.hyperlaneNodes, st.hyperlaneNet, p);

  st.hyperlaneEdgeBetweenness = std::move(r.edgeBetweenness);
  st.hyperlaneNodeBetweenness = std::move(r.nodeBetweenness);
  st.hyperlaneCentralityNodeIds = std::move(r.nodeIds);
  st.hyperlaneCentralityMaxEdge = r.maxEdge;
  st.hyperlaneCentralityMaxNode = r.maxNode;

  st.hyperlaneCentralityMs = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();

  // Cache top chokepoints (edges with highest betweenness).
  st.hyperlaneChokepointEdgeIdx.clear();
  if (st.hyperlaneTrafficHighlightChokepoints && st.hyperlaneTrafficTopEdges > 0 && !st.hyperlaneEdgeBetweenness.empty()) {
    int topN = st.hyperlaneTrafficTopEdges;
    if (topN < 1) topN = 1;
    if (topN > (int)st.hyperlaneEdgeBetweenness.size()) topN = (int)st.hyperlaneEdgeBetweenness.size();

    std::vector<int> idx;
    idx.reserve(st.hyperlaneEdgeBetweenness.size());
    for (int i = 0; i < (int)st.hyperlaneEdgeBetweenness.size(); ++i) idx.push_back(i);

    std::partial_sort(idx.begin(), idx.begin() + topN, idx.end(), [&](int a, int b) {
      const double va = st.hyperlaneEdgeBetweenness[(std::size_t)a];
      const double vb = st.hyperlaneEdgeBetweenness[(std::size_t)b];
      if (va != vb) return va > vb;
      return a < b;
    });

    idx.resize((std::size_t)topN);
    st.hyperlaneChokepointEdgeIdx = std::move(idx);
  }

  st.hyperlaneCentralityLastHash = h;
  st.hyperlaneCentralityDirty = false;
}


// -----------------------------------------------------------------------------
// Hyperlane structural criticality (bridges + articulation points)
// -----------------------------------------------------------------------------

static void rebuildHyperlaneVulnerability(ProceduralGalaxyLabWindowState& st) {
  const bool wantCritical = st.showHyperlanes && (st.hyperlaneColorMode == 3);
  if (!wantCritical) {
    // Keep cached values around, but don't spend time recomputing unless the
    // visualization is actively requested.
    return;
  }

  if (st.hyperlaneLastHash == core::u64(-1) || st.hyperlaneNodes.empty() || st.hyperlaneNet.edges.empty()) {
    st.hyperlaneVulnerabilityDirty = false;
    st.hyperlaneVulnerabilityLastHash = core::u64(-1);
    st.hyperlaneVulnerability = proc::HyperlaneVulnerabilityResult{};
    st.hyperlaneVulnerabilityMs = 0.0;
    st.hyperlaneCriticalBridgeEdgeIdx.clear();
    st.hyperlaneCriticalArticulationNodeIds.clear();
    st.hyperlaneCriticalArticulationNodeImpact01.clear();
    return;
  }

  int topEdges = st.hyperlaneCriticalTopBridges;
  if (topEdges < 1) topEdges = 1;
  if (topEdges > 128) topEdges = 128;

  int topNodes = st.hyperlaneCriticalTopArticulation;
  if (topNodes < 1) topNodes = 1;
  if (topNodes > 128) topNodes = 128;

  core::u64 h = core::fnv1a64("galaxy_preview/hyperlane_vulnerability/v1");
  h = hashMix(h, st.hyperlaneLastHash);
  h = hashMix(h, (core::u64)topEdges);
  h = hashMix(h, (core::u64)topNodes);
  h = hashMix(h, (core::u64)(st.hyperlaneCriticalHighlightBridges ? 1 : 0));
  h = hashMix(h, (core::u64)(st.hyperlaneCriticalHighlightArticulation ? 1 : 0));

  if (!st.hyperlaneVulnerabilityDirty && h == st.hyperlaneVulnerabilityLastHash) return;

  const auto t0 = Clock::now();

  st.hyperlaneVulnerability = proc::analyzeHyperlaneVulnerability(st.hyperlaneNodes, st.hyperlaneNet);
  st.hyperlaneVulnerabilityMs = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();

  // Cache top bridges by cut ratio (tie-break with absolute cut size).
  st.hyperlaneCriticalBridgeEdgeIdx.clear();
  if (st.hyperlaneCriticalHighlightBridges && topEdges > 0 && !st.hyperlaneVulnerability.edgeIsBridge.empty()) {
    std::vector<int> idx;
    idx.reserve(st.hyperlaneVulnerability.edgeIsBridge.size());
    for (int i = 0; i < (int)st.hyperlaneVulnerability.edgeIsBridge.size(); ++i) {
      if (st.hyperlaneVulnerability.edgeIsBridge[(std::size_t)i]) idx.push_back(i);
    }

    const int want = std::min<int>(topEdges, (int)idx.size());
    if (want > 0) {
      std::partial_sort(idx.begin(), idx.begin() + want, idx.end(), [&](int a, int b) {
        const float ca = st.hyperlaneVulnerability.edgeCut01[(std::size_t)a];
        const float cb = st.hyperlaneVulnerability.edgeCut01[(std::size_t)b];
        if (ca != cb) return ca > cb;
        const int sa = st.hyperlaneVulnerability.edgeCutSize[(std::size_t)a];
        const int sb = st.hyperlaneVulnerability.edgeCutSize[(std::size_t)b];
        if (sa != sb) return sa > sb;
        return a < b;
      });
      idx.resize((std::size_t)want);
      st.hyperlaneCriticalBridgeEdgeIdx = std::move(idx);
    }
  }

  // Cache top articulation systems by impact.
  st.hyperlaneCriticalArticulationNodeIds.clear();
  st.hyperlaneCriticalArticulationNodeImpact01.clear();
  if (st.hyperlaneCriticalHighlightArticulation && topNodes > 0 && !st.hyperlaneVulnerability.nodeIds.empty()) {
    std::vector<int> idx;
    idx.reserve(st.hyperlaneVulnerability.nodeIds.size());
    for (int i = 0; i < (int)st.hyperlaneVulnerability.nodeIds.size(); ++i) {
      if (st.hyperlaneVulnerability.nodeIsArticulation[(std::size_t)i]) idx.push_back(i);
    }

    const int want = std::min<int>(topNodes, (int)idx.size());
    if (want > 0) {
      std::partial_sort(idx.begin(), idx.begin() + want, idx.end(), [&](int a, int b) {
        const float ia = st.hyperlaneVulnerability.nodeImpact01[(std::size_t)a];
        const float ib = st.hyperlaneVulnerability.nodeImpact01[(std::size_t)b];
        if (ia != ib) return ia > ib;
        const int sa = st.hyperlaneVulnerability.nodeImpactSize[(std::size_t)a];
        const int sb = st.hyperlaneVulnerability.nodeImpactSize[(std::size_t)b];
        if (sa != sb) return sa > sb;
        const auto ida = st.hyperlaneVulnerability.nodeIds[(std::size_t)a];
        const auto idb = st.hyperlaneVulnerability.nodeIds[(std::size_t)b];
        return ida < idb;
      });
      idx.resize((std::size_t)want);

      st.hyperlaneCriticalArticulationNodeIds.reserve((std::size_t)want);
      st.hyperlaneCriticalArticulationNodeImpact01.reserve((std::size_t)want);
      for (int ni : idx) {
        st.hyperlaneCriticalArticulationNodeIds.push_back(st.hyperlaneVulnerability.nodeIds[(std::size_t)ni]);
        st.hyperlaneCriticalArticulationNodeImpact01.push_back(st.hyperlaneVulnerability.nodeImpact01[(std::size_t)ni]);
      }
    }
  }

  st.hyperlaneVulnerabilityLastHash = h;
  st.hyperlaneVulnerabilityDirty = false;
}


void drawProceduralGalaxyLabWindow(ProceduralGalaxyLabWindowState& st, float timeSec, const ToastFn& toast) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(1180, 760), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Procedural Galaxy Lab", &st.open)) {
    ImGui::End();
    return;
  }

  bool dirty = false;

  if (ImGui::BeginTable("##galaxy_lab_layout", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV)) {
    ImGui::TableNextColumn();

    // ----- Controls -----
    {
      ImGui::TextUnformatted("Preview generator params in a top-down XY slice.");
      ImGui::TextDisabled("Tip: set Spiral Arms > 0 and Arm Strength > 0 to enable.");
      ImGui::Separator();

      dirty |= ImGui::InputScalar("Seed", ImGuiDataType_U64, &st.seed);

      ImGui::SameLine();
      if (ImGui::Button("Randomize")) {
        const core::u64 t = static_cast<core::u64>(timeSec * 1'000'000.0f);
        const core::u64 h = core::hashCombine(st.seed, core::hashCombine(core::fnv1a64("galaxy_lab"), t));
        core::SplitMix64 rng(h);
        st.seed = rng.nextU64();
        dirty = true;
        if (toast) toast("Galaxy seed randomized", 2.0);
      }

      dirty |= ImGui::SliderInt("Factions", &st.factionCount, 1, 24);

      ImGui::Separator();
      ImGui::TextUnformatted("Galaxy Base");

      dirty |= dragDouble("Sector Size (ly)", st.params.sectorSizeLy, 0.1, 1.0, 200.0, "%.1f");
      dirty |= dragDouble("Disc Radius (ly)", st.params.radiusLy, 50.0, 500.0, 250'000.0, "%.0f");
      dirty |= dragDouble("Disc Thickness (ly)", st.params.thicknessLy, 10.0, 10.0, 25'000.0, "%.0f");
      dirty |= dragDouble("Radial Scale (ly)", st.params.radialScaleLengthLy, 50.0, 100.0, 250'000.0, "%.0f");
      dirty |= dragDouble("Mean Systems/Sector", st.params.baseMeanSystemsPerSector, 0.05, 0.0, 100.0, "%.2f");

      ImGui::Separator();
      ImGui::TextUnformatted("Star Placement");

      dirty |= dragDouble("Min System Separation (ly)", st.params.minSystemSeparationLy, 0.05, 0.0, 10.0, "%.2f");
      ImGui::TextDisabled("0 = legacy Poisson. >0 enables streaming-safe blue-noise-like placement (clamped to >= 0.25).");

      ImGui::Separator();
      ImGui::TextUnformatted("Spiral Arms (log spiral)");

      dirty |= ImGui::SliderInt("Arm Count", &st.params.spiralArmCount, 0, 8);
      dirty |= dragDouble("Arm Strength", st.params.spiralArmStrength, 0.02, 0.0, 5.0, "%.2f");
      dirty |= dragDouble("Pitch (deg)", st.params.spiralPitchDeg, 0.1, 1.0, 45.0, "%.1f");
      dirty |= dragDouble("Arm Width (deg)", st.params.spiralArmWidthDeg, 0.1, 1.0, 90.0, "%.1f");
      dirty |= dragDouble("Arm Phase (deg)", st.params.spiralArmPhaseDeg, 0.25, -360.0, 360.0, "%.2f");
      dirty |= dragDouble("Arm Noise Strength", st.params.spiralArmNoiseStrength, 0.01, 0.0, 1.0, "%.2f");
      dirty |= dragDouble("Arm Noise Freq", st.params.spiralArmNoiseFreq, 0.0001, 0.0, 0.05, "%.4f");

      ImGui::Separator();
      ImGui::TextUnformatted("Density Noise (clumpiness)");

      dirty |= dragDouble("Density Noise Strength", st.params.densityNoiseStrength, 0.01, 0.0, 1.0, "%.2f");
      dirty |= dragDouble("Density Noise Freq", st.params.densityNoiseFreq, 0.0001, 0.0, 0.05, "%.4f");

      ImGui::Separator();
      ImGui::TextUnformatted("Star Clusters (sparse blobs)");

      dirty |= dragDouble("Cluster Strength", st.params.clusterStrength, 0.02, 0.0, 10.0, "%.2f");
      dirty |= dragDouble("Cluster Cell Size (ly)", st.params.clusterCellSizeLy, 10.0, 50.0, 50'000.0, "%.0f");
      dirty |= dragDouble("Cluster Chance/Cell", st.params.clusterChancePerCell, 0.01, 0.0, 1.0, "%.2f");
      dirty |= dragDouble("Cluster Radius (ly)", st.params.clusterRadiusLy, 10.0, 10.0, 50'000.0, "%.0f");
      dirty |= dragDouble("Radius Jitter", st.params.clusterRadiusJitter, 0.01, 0.0, 1.0, "%.2f");
      dirty |= dragDouble("Strength Jitter", st.params.clusterStrengthJitter, 0.01, 0.0, 1.0, "%.2f");
      dirty |= dragDouble("Falloff Power", st.params.clusterFalloffPower, 0.05, 0.25, 8.0, "%.2f");
      ImGui::TextDisabled("Density multiplier: mean *= (1 + Strength * cluster01). Strength=0 disables.");

      ImGui::Separator();
      ImGui::TextUnformatted("Cosmic Voids (bubble cavities)");

      dirty |= dragDouble("Void Strength", st.params.voidStrength, 0.02, 0.0, 10.0, "%.2f");
      dirty |= dragDouble("Void Cell Size (ly)", st.params.voidCellSizeLy, 10.0, 50.0, 100'000.0, "%.0f");
      dirty |= dragDouble("Void Chance/Cell", st.params.voidChancePerCell, 0.01, 0.0, 1.0, "%.2f");
      dirty |= dragDouble("Void Radius (ly)", st.params.voidRadiusLy, 10.0, 10.0, 100'000.0, "%.0f");
      dirty |= dragDouble("Radius Jitter", st.params.voidRadiusJitter, 0.01, 0.0, 1.0, "%.2f");
      dirty |= dragDouble("Strength Jitter", st.params.voidStrengthJitter, 0.01, 0.0, 1.0, "%.2f");
      dirty |= dragDouble("Falloff Power", st.params.voidFalloffPower, 0.05, 0.25, 8.0, "%.2f");
      ImGui::TextDisabled("Density multiplier: mean *= clamp01(1 - Strength * void01). Strength=0 disables.");

      ImGui::Separator();
      ImGui::TextUnformatted("Galaxy Morphology (bar / ring / warp / flare)");

      ImGui::TextUnformatted("Bar");
      dirty |= dragDouble("Bar Strength", st.params.barStrength, 0.02, 0.0, 10.0, "%.2f");
      dirty |= dragDouble("Bar Angle (deg)", st.params.barAngleDeg, 0.5, -180.0, 180.0, "%.1f");
      dirty |= dragDouble("Bar Length (ly)", st.params.barLengthLy, 50.0, 10.0, 200'000.0, "%.0f");
      dirty |= dragDouble("Bar Width (ly)", st.params.barWidthLy, 25.0, 10.0, 200'000.0, "%.0f");
      dirty |= dragDouble("Bar Power", st.params.barPower, 0.05, 1.0, 8.0, "%.2f");

      ImGui::TextUnformatted("Ring");
      dirty |= dragDouble("Ring Strength", st.params.ringStrength, 0.02, 0.0, 10.0, "%.2f");
      dirty |= dragDouble("Ring Radius (ly)", st.params.ringRadiusLy, 50.0, 0.0, 200'000.0, "%.0f");
      dirty |= dragDouble("Ring Width (ly)", st.params.ringWidthLy, 25.0, 1.0, 50'000.0, "%.0f");
      dirty |= dragDouble("Ring Power", st.params.ringPower, 0.05, 1.0, 8.0, "%.2f");

      ImGui::TextUnformatted("Warp / Flare");
      dirty |= dragDouble("Warp Amplitude (ly)", st.params.warpAmplitudeLy, 10.0, 0.0, 50'000.0, "%.0f");
      dirty |= dragDouble("Warp Start Radius (ly)", st.params.warpStartRadiusLy, 50.0, 0.0, 200'000.0, "%.0f");
      dirty |= dragDouble("Warp Power", st.params.warpPower, 0.05, 0.25, 8.0, "%.2f");
      dirty |= ImGui::DragInt("Warp Lobes", &st.params.warpLobes, 1.0f, 1, 8);
      dirty |= dragDouble("Warp Phase (deg)", st.params.warpPhaseDeg, 1.0, -180.0, 180.0, "%.1f");
      dirty |= dragDouble("Warp Noise Strength", st.params.warpNoiseStrength, 0.01, 0.0, 1.0, "%.2f");
      dirty |= dragDouble("Warp Noise Freq", st.params.warpNoiseFreq, 0.0001, 0.0, 0.05, "%.4f");

      dirty |= dragDouble("Flare Strength", st.params.flareStrength, 0.02, 0.0, 10.0, "%.2f");
      dirty |= dragDouble("Flare Power", st.params.flarePower, 0.05, 0.25, 8.0, "%.2f");

      ImGui::TextDisabled("Bar/Ring: density multipliers. Warp/Flare: midplane + thickness.");

      ImGui::Separator();
      ImGui::TextUnformatted("View");

      dirty |= dragDouble("Center X (ly)", st.centerLy.x, 5.0, -500'000.0, 500'000.0, "%.0f");
      dirty |= dragDouble("Center Y (ly)", st.centerLy.y, 5.0, -500'000.0, 500'000.0, "%.0f");
      dirty |= dragDouble("Center Z (ly)", st.centerLy.z, 5.0, -500'000.0, 500'000.0, "%.0f");
      dirty |= dragDouble("View Radius (ly)", st.viewRadiusLy, 5.0, 10.0, 50'000.0, "%.0f");
      dirty |= dragDouble("Z Half-Range (ly)", st.zHalfLy, 5.0, 0.0, 25'000.0, "%.0f");

      // Preview cap. (Keep this adjustable; the sampler tries to pick a good LOD, but the cap is the hard guardrail.)
      {
        int maxPreview = static_cast<int>(std::min<std::size_t>(st.maxStubs, 5000000ull));
        if (ImGui::DragInt("Max Preview Systems", &maxPreview, 500.0f, 500, 250000, "%d")) {
          st.maxStubs = static_cast<std::size_t>(std::max(1, maxPreview));
          dirty = true;
        }
      }

      ImGui::Separator();
      ImGui::TextUnformatted("Preview Performance");

      const bool cacheToggled = ImGui::Checkbox("Sector Cache (LRU)", &st.previewUseSectorCache);
      dirty |= cacheToggled;

      if (st.previewUseSectorCache) {
        ImGui::Indent();

        dirty |= ImGui::SliderInt("Cache Entries", &st.previewSectorCacheMaxEntries, 0, 200000);
        ImGui::SameLine();
        if (ImGui::Button("Clear Cache")) {
          st.previewSectorCache.clear();
          st.previewSectorCacheStubTotal = 0;
          if (toast) toast("Galaxy preview cache cleared", 2.0);
        }

        const std::size_t cap = st.previewSectorCache.capacity();
        if (cap == 0) {
          ImGui::TextDisabled("Cache: %zu entries (unbounded)", st.previewSectorCache.size());
        } else {
          ImGui::TextDisabled("Cache: %zu / %zu entries", st.previewSectorCache.size(), cap);
        }

        const double approxMb = (double)st.previewSectorCacheStubTotal * (double)sizeof(sim::SystemStub) / (1024.0 * 1024.0);
        ImGui::TextDisabled("Cached stubs: %zu (~%.2f MB)", st.previewSectorCacheStubTotal, approxMb);
        ImGui::TextDisabled("Hits: %d  Misses: %d  Evict: %d", st.previewSectorCacheHits, st.previewSectorCacheMisses, st.previewSectorCacheEvictions);
        ImGui::TextDisabled("Deterministic sectors are cached; ideal for pan/zoom exploration.");
        ImGui::Unindent();
      }


      // Parallel CPU preview generation. Useful when sampling many sectors (large view radius)
      // or when Max Preview Systems is high.
      const bool parallelToggled = ImGui::Checkbox("Parallel Preview", &st.previewParallel);
      dirty |= parallelToggled;
      if (st.previewParallel) {
        ImGui::Indent();

        int threads = st.previewParallelThreads;
        if (ImGui::DragInt("Threads (0=auto)", &threads, 0.25f, 0, 64)) {
          st.previewParallelThreads = threads;
          dirty = true;
        }

        ImGui::TextDisabled("Last rebuild used %d worker threads", st.previewParallelWorkersUsed);
        ImGui::Unindent();
      }
      ImGui::Separator();
      ImGui::TextUnformatted("Custom Rendering (Heatmap)");

      ImGui::Checkbox("GPU Heatmap", &st.heatmapEnabled);
      if (st.heatmapEnabled) {
        ImGui::Indent();

        const char* modeItems = "Density\0Spiral Arms\0Clusters\0Voids\0Morphology\0";
        ImGui::Combo("Mode", &st.heatmapMode, modeItems);

        int resIdx = 1; // half by default
        if (st.heatmapResolutionDiv == 1) resIdx = 0;
        else if (st.heatmapResolutionDiv == 2) resIdx = 1;
        else resIdx = 2;

        const char* resItems = "Full\0Half\0Quarter\0";
        if (ImGui::Combo("Resolution", &resIdx, resItems)) {
          st.heatmapResolutionDiv = (resIdx == 0) ? 1 : (resIdx == 1 ? 2 : 4);
        }

        ImGui::Checkbox("Animate", &st.heatmapAnimate);
        ImGui::SameLine();
        if (ImGui::Button("Force Re-render")) {
          st.heatmapLastHash = core::u64(-1);
          st.heatmapFrame = 0;
        }
        ImGui::TextDisabled("When Animate is off, the heatmap is cached and only redraws when inputs change.");

        ImGui::SliderFloat("Exposure", &st.heatmapExposure, 0.01f, 2.0f, "%.3f");
        ImGui::SliderFloat("Gamma", &st.heatmapGamma, 0.60f, 2.50f, "%.2f");

        ImGui::Checkbox("Contours", &st.heatmapContours);
        if (st.heatmapContours) {
          ImGui::SliderInt("Contour Count", &st.heatmapContourCount, 2, 64);
          ImGui::SliderFloat("Contour Width", &st.heatmapContourWidth, 0.005f, 0.08f, "%.3f");
        }

        ImGui::Separator();
        ImGui::TextUnformatted("Dynamic LOD (shader)");

        ImGui::Checkbox("Dynamic LOD", &st.heatmapDynamicLod);
        if (st.heatmapDynamicLod) {
          ImGui::Indent();
          ImGui::SliderFloat("LOD Px Lo", &st.heatmapLodPxLo, 0.10f, 4.0f, "%.2f");
          ImGui::SliderFloat("LOD Px Hi", &st.heatmapLodPxHi, 0.20f, 8.0f, "%.2f");
          ImGui::SliderFloat("Energy Preserve", &st.heatmapLodEnergy, 0.0f, 1.0f, "%.2f");
          ImGui::Checkbox("Skip tiny blob scans", &st.heatmapLodSkipTinyBlobs);
          ImGui::TextDisabled("Band-limits procedural FBM detail (reduces shimmer at distance).");
          ImGui::Unindent();
        }

        if (!st.heatmapError.empty()) {
          ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f), "Heatmap: %s", st.heatmapError.c_str());
        }

        ImGui::TextDisabled("Analytic GPU heatmap behind points (approx; may not match CPU exactly). ");
        ImGui::Unindent();
      }

      ImGui::Separator();
      ImGui::TextUnformatted("Hyperlanes (Navigation Graph)");

      // This overlay visualizes a deterministic, sparse navigation graph on top of the preview stubs.
      if (ImGui::Checkbox("Show Hyperlanes", &st.showHyperlanes)) {
        st.hyperlanesDirty = true;
      }

      if (st.showHyperlanes) {
        ImGui::Indent();
        ImGui::TextDisabled("Deterministic MST + kNN + extras. Great for corridors, trade, and route feel.");

        int maxNodes = (int)st.hyperlaneMaxNodes;
        if (maxNodes < 0) maxNodes = 0;
        if (ImGui::DragInt("Max Nodes (0=all)", &maxNodes, 250.0f, 0, 50000)) {
          st.hyperlaneMaxNodes = (std::size_t)std::max(0, maxNodes);
          st.hyperlanesDirty = true;
        }

        if (dragDouble("Neighbor Dist (ly)", st.hyperlaneParams.maxNeighborDistanceLy, 0.25, 1.0, 200.0, "%.2f")) st.hyperlanesDirty = true;

        int k = st.hyperlaneParams.neighborK;
        if (ImGui::DragInt("Neighbor K", &k, 1.0f, 1, 16)) {
          st.hyperlaneParams.neighborK = k;
          st.hyperlanesDirty = true;
        }

        int minDeg = st.hyperlaneParams.minDegree;
        if (ImGui::DragInt("Min Degree", &minDeg, 1.0f, 0, 8)) {
          st.hyperlaneParams.minDegree = minDeg;
          st.hyperlanesDirty = true;
        }

        if (ImGui::Checkbox("Force connected", &st.hyperlaneParams.forceConnected)) st.hyperlanesDirty = true;
        if (dragDouble("Extra Edge Chance", st.hyperlaneParams.extraEdgeChance, 0.01, 0.0, 1.0, "%.2f")) st.hyperlanesDirty = true;

        if (dragDouble("Region Cell Size (ly)", st.hyperlaneParams.regionCellSizeLy, 10.0, 0.0, 5000.0, "%.0f")) st.hyperlanesDirty = true;
        ImGui::TextDisabled("Set Region Cell Size <= 0 to disable region modulation.");

        if (ImGui::TreeNode("Hazard Modulation (nebula/storms)")) {
          if (ImGui::Checkbox("Enabled", &st.hyperlaneParams.hazards.enabled)) st.hyperlanesDirty = true;
          if (st.hyperlaneParams.hazards.enabled) {
            if (dragDouble("Time (days)", st.hyperlaneParams.hazards.timeDays, 0.25, 0.0, 500000.0, "%.2f")) st.hyperlanesDirty = true;
            if (dragDouble("Drift (ly/day)", st.hyperlaneParams.hazards.driftLyPerDay, 0.05, 0.0, 10.0, "%.2f")) st.hyperlanesDirty = true;
            if (dragDouble("Nebula Scale (ly)", st.hyperlaneParams.hazards.nebulaScaleLy, 10.0, 10.0, 5000.0, "%.0f")) st.hyperlanesDirty = true;
            if (dragDouble("Storm Scale (ly)", st.hyperlaneParams.hazards.stormScaleLy, 10.0, 10.0, 5000.0, "%.0f")) st.hyperlanesDirty = true;

            int es = st.hyperlaneParams.hazards.edgeSamples;
            if (ImGui::DragInt("Edge Samples", &es, 1.0f, 1, 11)) {
              st.hyperlaneParams.hazards.edgeSamples = es;
              st.hyperlanesDirty = true;
            }

            if (dragDouble("Hazard->Cost", st.hyperlaneParams.hazards.hazardToCost, 0.01, 0.0, 2.0, "%.2f")) st.hyperlanesDirty = true;
            if (dragDouble("Hazard->Risk", st.hyperlaneParams.hazards.hazardToRisk, 0.01, 0.0, 2.0, "%.2f")) st.hyperlanesDirty = true;
            if (dragDouble("Hazard->Bandwidth", st.hyperlaneParams.hazards.hazardToBandwidth, 0.01, 0.0, 0.95, "%.2f")) st.hyperlanesDirty = true;
          }
          ImGui::TreePop();
        }

        ImGui::Separator();
        ImGui::TextUnformatted("Draw LOD");

        int maxEdgesDraw = (int)st.hyperlaneMaxEdgesDraw;
        if (maxEdgesDraw < 0) maxEdgesDraw = 0;
        if (ImGui::DragInt("Max Edges (0=auto)", &maxEdgesDraw, 50.0f, 0, 200000)) {
          st.hyperlaneMaxEdgesDraw = (std::size_t)std::max(0, maxEdgesDraw);
          // If generation was previously truncated, we need a rebuild to expose more edges.
          st.hyperlanesDirty = true;
        }

        ImGui::SliderFloat("Opacity", &st.hyperlaneOpacity, 0.0f, 1.0f, "%.2f");
        ImGui::SliderFloat("Min Len (px)", &st.hyperlaneMinLenPx, 0.0f, 12.0f, "%.1f");
        ImGui::SliderFloat("Width (px)", &st.hyperlaneWidthPx, 0.5f, 4.0f, "%.1f");

        {
          const char* kColorItems = "Risk\0Bandwidth\0Traffic (betweenness)\0Criticality (bridges/articulation)\0";
          int mode = st.hyperlaneColorMode;
          if (ImGui::Combo("Color Mode", &mode, kColorItems)) {
            st.hyperlaneColorMode = mode;
            if (st.hyperlaneColorMode == 2) st.hyperlaneCentralityDirty = true;
            if (st.hyperlaneColorMode == 3) st.hyperlaneVulnerabilityDirty = true;
          }

          if (st.hyperlaneColorMode == 2) {
            ImGui::Indent();
            ImGui::TextDisabled("Traffic = betweenness centrality over shortest paths (approx).");

            int samples = st.hyperlaneTrafficSamples;
            if (ImGui::SliderInt("Traffic Samples (0=exact)", &samples, 0, 512)) {
              st.hyperlaneTrafficSamples = samples;
              st.hyperlaneCentralityDirty = true;
            }

            if (ImGui::Checkbox("Highlight chokepoints", &st.hyperlaneTrafficHighlightChokepoints)) {
              st.hyperlaneCentralityDirty = true;
            }
            if (st.hyperlaneTrafficHighlightChokepoints) {
              ImGui::Indent();
              if (ImGui::SliderInt("Top Chokepoints", &st.hyperlaneTrafficTopEdges, 1, 64)) {
                st.hyperlaneCentralityDirty = true;
              }
              ImGui::Unindent();
            }

            ImGui::SliderFloat("Traffic Width Boost", &st.hyperlaneTrafficWidthBoost, 0.25f, 3.0f, "%.2f");
            ImGui::SliderFloat("Traffic Alpha Boost", &st.hyperlaneTrafficAlphaBoost, 0.25f, 2.0f, "%.2f");

            if (ImGui::Button("Recompute Traffic")) {
              st.hyperlaneCentralityLastHash = core::u64(-1);
              st.hyperlaneCentralityDirty = true;
            }

            if (st.hyperlaneCentralityLastHash != core::u64(-1)) {
              if (st.hyperlaneCentralityMs > 0.0) {
                ImGui::TextDisabled("Traffic cached: %.2f ms (max edge %.3g, max node %.3g)",
                                  st.hyperlaneCentralityMs, st.hyperlaneCentralityMaxEdge, st.hyperlaneCentralityMaxNode);
              }
              if (!st.hyperlaneChokepointEdgeIdx.empty()) {
                ImGui::TextDisabled("Top chokepoint edges: %d", (int)st.hyperlaneChokepointEdgeIdx.size());
              }
            }

            ImGui::Unindent();
          }

          if (st.hyperlaneColorMode == 3) {
            ImGui::Indent();
            ImGui::TextDisabled("Criticality = structural fragility (bridges + articulation systems).");

            if (ImGui::Checkbox("Highlight bridges", &st.hyperlaneCriticalHighlightBridges)) {
              st.hyperlaneVulnerabilityDirty = true;
            }
            if (ImGui::Checkbox("Highlight articulation systems", &st.hyperlaneCriticalHighlightArticulation)) {
              st.hyperlaneVulnerabilityDirty = true;
            }

            ImGui::Indent();
            if (ImGui::SliderInt("Top Bridges", &st.hyperlaneCriticalTopBridges, 1, 64)) {
              st.hyperlaneVulnerabilityDirty = true;
            }
            if (ImGui::SliderInt("Top Articulation", &st.hyperlaneCriticalTopArticulation, 1, 64)) {
              st.hyperlaneVulnerabilityDirty = true;
            }
            ImGui::Unindent();

            ImGui::SliderFloat("Critical Width Boost", &st.hyperlaneCriticalWidthBoost, 0.25f, 3.0f, "%.2f");
            ImGui::SliderFloat("Critical Alpha Boost", &st.hyperlaneCriticalAlphaBoost, 0.25f, 2.0f, "%.2f");
            ImGui::SliderFloat("Articulation Ring Radius", &st.hyperlaneCriticalNodeRingRadiusPx, 2.0f, 18.0f, "%.1f");
            ImGui::SliderFloat("Articulation Ring Thickness", &st.hyperlaneCriticalNodeRingThicknessPx, 0.5f, 6.0f, "%.1f");

            if (ImGui::Button("Recompute Criticality")) {
              st.hyperlaneVulnerabilityLastHash = core::u64(-1);
              st.hyperlaneVulnerabilityDirty = true;
            }

            if (st.hyperlaneVulnerabilityLastHash != core::u64(-1)) {
              ImGui::TextDisabled("Criticality cached: %.2f ms (bridges %d, articulations %d, components %d)",
                                st.hyperlaneVulnerabilityMs,
                                st.hyperlaneVulnerability.bridgeCount,
                                st.hyperlaneVulnerability.articulationCount,
                                st.hyperlaneVulnerability.componentCount);
            }

            ImGui::Unindent();
          }
        }

        ImGui::Checkbox("Highlight hovered", &st.hyperlaneHighlightHovered);

        if (st.hyperlaneNodesUsed > 0 || st.hyperlaneEdgesUsed > 0) {
          ImGui::TextDisabled("Last build: %d nodes, %d edges (%.2f ms)", st.hyperlaneNodesUsed, st.hyperlaneEdgesUsed, st.hyperlaneGenMs);
        }
        if (ImGui::TreeNode("Routing / Pathfinding")) {
          ImGui::Checkbox("Enable Route Inspector", &st.hyperlaneRouteEnabled);
          if (st.hyperlaneRouteEnabled) {
            ImGui::Indent();
            ImGui::TextDisabled("Shift+Click = start, Ctrl+Click = end, Right-click = clear");

            auto showSel = [&](const char* label, sim::SystemId id) {
              if (id == 0) {
                ImGui::Text("%s: (none)", label);
                return;
              }
              const sim::SystemStub* ss = findStubById(st.stubs, id);
              if (ss) {
                ImGui::Text("%s: %s", label, ss->name.c_str());
              } else {
                ImGui::Text("%s: 0x%llx", label, (unsigned long long)id);
              }
            };

            showSel("Start", st.hyperlaneRouteFrom);
            showSel("End", st.hyperlaneRouteTo);

            if (ImGui::Button("Clear Route")) {
              st.hyperlaneRouteFrom = 0;
              st.hyperlaneRouteTo = 0;
              st.hyperlaneRouteDirty = true;
              if (toast) toast("Route cleared", 1.5);
            }
            ImGui::SameLine();
            if (ImGui::Button("Swap")) {
              std::swap(st.hyperlaneRouteFrom, st.hyperlaneRouteTo);
              st.hyperlaneRouteDirty = true;
            }

            ImGui::Separator();
            ImGui::TextUnformatted("Travel preference");

            if (dragDouble("Risk Weight", st.hyperlaneTravel.riskWeight, 0.05, 0.0, 3.0, "%.2f")) st.hyperlaneRouteDirty = true;
            if (dragDouble("Bandwidth Bias", st.hyperlaneTravel.bandwidthBias, 0.01, 0.0, 1.0, "%.2f")) st.hyperlaneRouteDirty = true;
            if (dragDouble("Min Bandwidth Factor", st.hyperlaneTravel.minBandwidthFactor, 0.01, 0.05, 1.0, "%.2f")) st.hyperlaneRouteDirty = true;
            ImGui::TextDisabled("Lower cost prefers high bandwidth / low risk edges.");

            ImGui::Separator();
            ImGui::TextUnformatted("Alternatives (K-shortest)");

            int kAlt = st.hyperlaneRouteK;
            if (ImGui::SliderInt("K Routes", &kAlt, 1, 8)) {
              st.hyperlaneRouteK = kAlt;
              st.hyperlaneRouteDirty = true;
            }

            ImGui::Checkbox("Draw alternatives", &st.hyperlaneRouteDrawAlternatives);
            if (st.hyperlaneRouteDrawAlternatives) {
              ImGui::Indent();
              ImGui::SliderFloat("Alt Opacity", &st.hyperlaneRouteAltOpacity, 0.0f, 1.0f, "%.2f");
              ImGui::SliderFloat("Alt Width Scale", &st.hyperlaneRouteAltWidthScale, 0.25f, 1.25f, "%.2f");
              ImGui::Unindent();
            }

            if (!st.hyperlaneKRoutes.empty()) {
              ImGui::TextDisabled("Found %d route(s). Click to preview:", (int)st.hyperlaneKRoutes.size());

              if (ImGui::BeginTable("##hyperlane_k_routes", 6,
                                    ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH |
                                      ImGuiTableFlags_SizingFixedFit)) {
                ImGui::TableSetupColumn("#");
                ImGui::TableSetupColumn("Hops");
                ImGui::TableSetupColumn("Cost");
                ImGui::TableSetupColumn("Dist");
                ImGui::TableSetupColumn("Risk");
                ImGui::TableSetupColumn("BW");
                ImGui::TableHeadersRow();

                for (int i = 0; i < (int)st.hyperlaneKRoutes.size(); ++i) {
                  const auto& r = st.hyperlaneKRoutes[(std::size_t)i];
                  ImGui::TableNextRow();

                  ImGui::TableNextColumn();
                  char lbl[16];
                  std::snprintf(lbl, sizeof(lbl), "R%d", i + 1);
                  const bool selected = (st.hyperlaneRouteSelected == i);
                  if (ImGui::Selectable(lbl, selected, ImGuiSelectableFlags_SpanAllColumns)) {
                    st.hyperlaneRouteSelected = i;
                  }

                  ImGui::TableNextColumn();
                  ImGui::Text("%d", r.metrics.hops);

                  ImGui::TableNextColumn();
                  ImGui::Text("%.2f", (float)r.metrics.costLy);

                  ImGui::TableNextColumn();
                  ImGui::Text("%.2f", (float)r.metrics.distanceLy);

                  ImGui::TableNextColumn();
                  ImGui::Text("%.2f", (float)r.metrics.risk01);

                  ImGui::TableNextColumn();
                  ImGui::Text("%.2f", (float)r.metrics.bottleneckBandwidth01);
                }

                ImGui::EndTable();
              }
            }


            ImGui::Separator();
            ImGui::TextUnformatted("Route rendering");
            ImGui::Checkbox("Draw Route Overlay", &st.hyperlaneRouteDraw);
            ImGui::SliderFloat("Route Opacity", &st.hyperlaneRouteOpacity, 0.0f, 1.0f, "%.2f");
            ImGui::SliderFloat("Route Width (px)", &st.hyperlaneRouteWidthPx, 1.0f, 10.0f, "%.1f");

            ImGui::Separator();
            if (!st.hyperlaneRouteError.empty()) {
              ImGui::TextColored(ImVec4(1.0f, 0.40f, 0.40f, 1.0f), "Route: %s", st.hyperlaneRouteError.c_str());
            }

            if (!st.hyperlaneRoutePath.empty()) {
              ImGui::Text("Path: %d hops", (int)std::max<int>(0, (int)st.hyperlaneRoutePath.size() - 1));
              ImGui::Text("Cost: %.2f   Distance: %.2f ly", (float)st.hyperlaneRouteMetrics.costLy, (float)st.hyperlaneRouteMetrics.distanceLy);
              ImGui::Text("Risk: %.2f   Bottleneck BW: %.2f", (float)st.hyperlaneRouteMetrics.risk01, (float)st.hyperlaneRouteMetrics.bottleneckBandwidth01);
              ImGui::TextDisabled("Last solve: %.2f ms", (float)st.hyperlaneRouteMs);

              if (ImGui::Button("Copy Route To Clipboard")) {
                std::string txt;
                txt.reserve(st.hyperlaneRoutePath.size() * 24);
                for (std::size_t i = 0; i < st.hyperlaneRoutePath.size(); ++i) {
                  const sim::SystemId sid = st.hyperlaneRoutePath[i];
                  const sim::SystemStub* ss = findStubById(st.stubs, sid);
                  if (ss) {
                    txt += ss->name;
                  } else {
                    char buf[32];
                    std::snprintf(buf, sizeof(buf), "0x%llx", (unsigned long long)sid);
                    txt += buf;
                  }
                  if (i + 1 < st.hyperlaneRoutePath.size()) txt += " -> ";
                }
                ImGui::SetClipboardText(txt.c_str());
                if (toast) toast("Copied route to clipboard", 1.5);
              }
            } else if (st.hyperlaneRouteFrom != 0 && st.hyperlaneRouteTo != 0 && st.hyperlaneRouteFrom != st.hyperlaneRouteTo) {
              ImGui::TextDisabled("(Route will update after next graph rebuild)");
            }

            ImGui::Unindent();
          }
          ImGui::TreePop();
        }

        ImGui::Unindent();
      }

      dirty |= ImGui::Checkbox("Auto regenerate", &st.autoRegenerate);

      const bool changedFaction = ImGui::Checkbox("Color: faction", &st.colorByFaction);
      ImGui::SameLine();
      const bool changedRegion = ImGui::Checkbox("region", &st.colorByRegion);
      ImGui::SameLine();
      const bool changedCluster = ImGui::Checkbox("cluster", &st.colorByCluster);
      ImGui::SameLine();
      const bool changedVoid = ImGui::Checkbox("void", &st.colorByVoid);
      ImGui::SameLine();
      const bool changedMorph = ImGui::Checkbox("morph", &st.colorByMorphology);
      ImGui::SameLine();
      dirty |= ImGui::Checkbox("Legend", &st.showLegend);

      const bool colorChanged = changedFaction || changedRegion || changedCluster || changedVoid || changedMorph;
      dirty |= colorChanged;

      // Make color modes mutually exclusive to avoid ambiguity.
      // IMPORTANT: pick the one that actually changed this frame.
      if (changedFaction && st.colorByFaction) {
        st.colorByRegion = false;
        st.colorByCluster = false;
        st.colorByVoid = false;
        st.colorByMorphology = false;
      }
      if (changedRegion && st.colorByRegion) {
        st.colorByFaction = false;
        st.colorByCluster = false;
        st.colorByVoid = false;
        st.colorByMorphology = false;
      }
      if (changedCluster && st.colorByCluster) {
        st.colorByFaction = false;
        st.colorByRegion = false;
        st.colorByVoid = false;
        st.colorByMorphology = false;
      }
      if (changedVoid && st.colorByVoid) {
        st.colorByFaction = false;
        st.colorByRegion = false;
        st.colorByCluster = false;
        st.colorByMorphology = false;
      }
      if (changedMorph && st.colorByMorphology) {
        st.colorByFaction = false;
        st.colorByRegion = false;
        st.colorByCluster = false;
        st.colorByVoid = false;
      }

      if (st.colorByRegion) {
        ImGui::Indent();
        dirty |= dragDouble("Region Cell Size (ly)", st.regionCellSizeLy, 10.0, 100.0, 5000.0, "%.0f");
        ImGui::Unindent();
      }

      ImGui::Checkbox("Arm guides", &st.showArmGuides);
      ImGui::SameLine();
      ImGui::Checkbox("Morph guides", &st.showMorphGuides);

      if (ImGui::Button("Regenerate")) {
        st.dirty = true;
      }

      ImGui::SameLine();
      if (ImGui::Button("Reset Params")) {
        st.params = proc::GalaxyParams{};
        st.dirty = true;
        if (toast) toast("Galaxy params reset", 2.0);
      }

      ImGui::Separator();
      ImGui::Text("Preview: %zu systems", st.stubs.size());
      ImGui::Text("Last gen: %.2f ms", st.lastGenMs);
      if (st.previewSectorsGenerated > 0) {
        ImGui::Text("Sectors sampled: %d  (stride XY: %d)", st.previewSectorsGenerated, st.previewSectorStrideXY);
      }
      if (st.previewCandidateStubs > 0) {
        ImGui::Text("Candidates: %d", st.previewCandidateStubs);
      }
      if (st.stubs.size() >= st.maxStubs) {
        ImGui::TextDisabled("(hit maxStubs cap: %zu)", st.maxStubs);
      }
    }

    ImGui::TableNextColumn();

    // ----- Map preview -----
    if (dirty) st.dirty = true;
    if (st.dirty && st.autoRegenerate) {
      rebuildPreview(st);
    }

    const ImVec2 canvasSize = ImGui::GetContentRegionAvail();
    ImGui::BeginChild("##galaxy_canvas", canvasSize, true, ImGuiWindowFlags_NoScrollWithMouse | ImGuiWindowFlags_NoScrollbar);

    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    const ImVec2 sz = ImGui::GetContentRegionAvail();
    const ImVec2 p1 = ImVec2(p0.x + sz.x, p0.y + sz.y);

    ImDrawList* dl = ImGui::GetWindowDrawList();

    // Optional GPU heatmap background (analytic approximation).
    if (st.heatmapEnabled) {
      renderGalaxyHeatmap(st, (int)sz.x, (int)sz.y, timeSec);
    }

    if (st.heatmapEnabled && st.heatmapInited && st.heatmapTarget.isInited() && st.heatmapTarget.color().handle() != 0) {
      const ImTextureID texId = (ImTextureID)(intptr_t)st.heatmapTarget.color().handle();
      dl->AddImage(texId, p0, p1, ImVec2(0, 1), ImVec2(1, 0));
      // Subtle dark overlay so point colors remain readable.
      dl->AddRectFilled(p0, p1, rgba(0.06f, 0.07f, 0.10f, 0.18f));
    } else {
      dl->AddRectFilled(p0, p1, rgba(0.06f, 0.07f, 0.10f));
    }

    dl->AddRect(p0, p1, rgba(0.30f, 0.30f, 0.35f));

    if (st.heatmapEnabled && !st.heatmapError.empty()) {
      dl->AddText(ImVec2(p0.x + 10.0f, p0.y + 8.0f), rgba(1.0f, 0.40f, 0.40f, 0.95f), "Heatmap error (see controls)");
    }

    const ImVec2 centerPx = ImVec2((p0.x + p1.x) * 0.5f, (p0.y + p1.y) * 0.5f);

    const double radiusLy = std::max(1.0, st.viewRadiusLy);
    const float scale = static_cast<float>(std::min(sz.x, sz.y) / (2.0 * radiusLy));

    const auto worldToScreen = [&](const math::Vec3d& w) -> ImVec2 {
      const double dx = w.x - st.centerLy.x;
      const double dy = w.y - st.centerLy.y;
      return ImVec2(centerPx.x + static_cast<float>(dx * scale),
                   centerPx.y - static_cast<float>(dy * scale));
    };

    // Hyperlane overlay caches (built lazily per-frame so tooltip/highlight can reuse them).
    std::unordered_map<sim::SystemId, ImVec2> hyperNodePos{};
    bool hyperNodePosReady = false;
    auto ensureHyperNodePos = [&]() {
      if (hyperNodePosReady) return;
      hyperNodePos.clear();
      hyperNodePos.reserve(st.hyperlaneNodes.size() * 2);
      for (const auto& n : st.hyperlaneNodes) {
        hyperNodePos.emplace(n.id, worldToScreen(n.posLy));
      }
      hyperNodePosReady = true;
    };

    // Axes + boundary.
    dl->AddLine(ImVec2(p0.x, centerPx.y), ImVec2(p1.x, centerPx.y), rgba(0.22f, 0.22f, 0.25f));
    dl->AddLine(ImVec2(centerPx.x, p0.y), ImVec2(centerPx.x, p1.y), rgba(0.22f, 0.22f, 0.25f));
    dl->AddCircle(centerPx, static_cast<float>(radiusLy * scale), rgba(0.25f, 0.25f, 0.28f));

    // Optional arm guides (only meaningful if you're looking near the galactic center).
    if (st.showArmGuides && st.params.spiralArmCount > 0 && st.params.spiralArmStrength > 0.0) {
      const int arms = st.params.spiralArmCount;
      const double pitchDeg = std::clamp(st.params.spiralPitchDeg, 1.0, 89.0);
      const double pitchRad = pitchDeg * (stellar::math::kPi / 180.0);
      const double k = 1.0 / std::tan(pitchRad);
      const double phaseRad = st.params.spiralArmPhaseDeg * (stellar::math::kPi / 180.0);
      const double twoPi = 2.0 * stellar::math::kPi;
      const double rRef = std::max(1.0, st.params.radiusLy * 0.02);

      // Draw each arm as a polyline. We sample r in log space.
      const int steps = 220;
      for (int a = 0; a < arms; ++a) {
        ImVec2 prev{};
        bool havePrev = false;
        for (int i = 0; i < steps; ++i) {
          const double t = static_cast<double>(i) / static_cast<double>(steps - 1);
          const double r = rRef * std::exp(std::log(std::max(1.0, radiusLy / rRef)) * t);
          const double lnTerm = std::log(std::max(1.0e-6, r / rRef));
          const double theta = k * lnTerm + phaseRad + twoPi * (static_cast<double>(a) / static_cast<double>(arms));

          const math::Vec3d w{r * std::cos(theta), r * std::sin(theta), st.centerLy.z};
          const ImVec2 p = worldToScreen(w);
          if (havePrev) {
            dl->AddLine(prev, p, rgba(0.16f, 0.19f, 0.28f, 0.65f));
          }
          prev = p;
          havePrev = true;
        }
      }
    }

    // Morphology guides overlay (bar + ring).
    if (st.showMorphGuides && (st.params.barStrength > 0.0 || st.params.ringStrength > 0.0)) {
      const double kPi = 3.14159265358979323846;

      // Bar: draw the ellipse in world space (galactic center at origin), rotated by Bar Angle.
      if (st.params.barStrength > 0.0 && st.params.barLengthLy > 0.0 && st.params.barWidthLy > 0.0) {
        const double a = st.params.barLengthLy;
        const double b = st.params.barWidthLy;
        const double ang = math::degToRad(st.params.barAngleDeg);
        const double ca = std::cos(ang);
        const double sa = std::sin(ang);

        ImVec2 prev{};
        bool havePrev = false;
        const int N = 128;
        for (int i = 0; i <= N; ++i) {
          const double t = (2.0 * kPi) * (static_cast<double>(i) / static_cast<double>(N));
          const double xLocal = a * std::cos(t);
          const double yLocal = b * std::sin(t);

          // Rotate by +ang into world coordinates.
          const double x = ca * xLocal - sa * yLocal;
          const double y = sa * xLocal + ca * yLocal;

          const ImVec2 p = worldToScreen(math::Vec3d{x, y, st.centerLy.z});
          if (havePrev) {
            dl->AddLine(prev, p, rgba(0.95f, 0.65f, 0.25f, 0.55f));
          }
          prev = p;
          havePrev = true;
        }
      }

      // Ring: draw a circle at Ring Radius.
      if (st.params.ringStrength > 0.0 && st.params.ringRadiusLy > 0.0) {
        const double r = st.params.ringRadiusLy;
        ImVec2 prev{};
        bool havePrev = false;
        const int N = 128;
        for (int i = 0; i <= N; ++i) {
          const double t = (2.0 * kPi) * (static_cast<double>(i) / static_cast<double>(N));
          const double x = r * std::cos(t);
          const double y = r * std::sin(t);

          const ImVec2 p = worldToScreen(math::Vec3d{x, y, st.centerLy.z});
          if (havePrev) {
            dl->AddLine(prev, p, rgba(0.35f, 0.85f, 0.75f, 0.55f));
          }
          prev = p;
          havePrev = true;
        }
      }
    }


    // Hyperlane overlay (procedural navigation graph) - drawn behind points.
    if (st.showHyperlanes) {
      rebuildHyperlanes(st);
      rebuildHyperlaneRoute(st);
      rebuildHyperlaneCentrality(st);
      rebuildHyperlaneVulnerability(st);

      if (!st.hyperlaneNet.edges.empty() && !st.hyperlaneNodes.empty()) {
        ensureHyperNodePos();

        std::size_t maxEdgesDraw = st.hyperlaneMaxEdgesDraw;
        if (maxEdgesDraw == 0) {
          const double area = (double)sz.x * (double)sz.y;
          maxEdgesDraw = (std::size_t)std::clamp(area / 45.0, 2500.0, 50000.0);
        }

        const float minLenPx = std::max(0.0f, st.hyperlaneMinLenPx);
        const float baseAlpha = std::clamp(st.hyperlaneOpacity, 0.0f, 1.0f);
        const float baseWidth = std::clamp(st.hyperlaneWidthPx, 0.5f, 8.0f);

        const auto smoothstep = [](float a, float b, float x) -> float {
          const float denom = std::max(1e-6f, b - a);
          const float t = std::clamp((x - a) / denom, 0.0f, 1.0f);
          return t * t * (3.0f - 2.0f * t);
        };

        const auto lerpF = [](float a, float b, float t) -> float { return a + (b - a) * t; };

        const auto laneColor = [&](float bw01, float risk01, float traffic01, float critical01, float alpha) -> ImU32 {
          float r = 0.20f, g = 0.85f, b = 1.00f;
          const int mode = st.hyperlaneColorMode;
          if (mode == 0) {
            const float t = std::clamp(risk01, 0.0f, 1.0f);
            // Low risk: cyan/blue; High risk: red/orange.
            r = lerpF(0.20f, 1.00f, t);
            g = lerpF(0.85f, 0.32f, t);
            b = lerpF(1.00f, 0.15f, t);
          } else if (mode == 1) {
            const float u = std::clamp(bw01, 0.0f, 1.0f);
            // Bandwidth gradient (low=dim, high=bright).
            r = lerpF(0.18f, 0.95f, u);
            g = lerpF(0.20f, 0.95f, u);
            b = lerpF(0.25f, 0.95f, u);
          } else if (mode == 2) {
            const float t = std::clamp(traffic01, 0.0f, 1.0f);
            // Traffic (betweenness): low=dim, high=hot (chokepoints).
            r = lerpF(0.15f, 1.00f, t);
            g = lerpF(0.20f, 0.92f, t);
            b = lerpF(0.35f, 0.20f, t);
          } else {
            const float t = std::clamp(critical01, 0.0f, 1.0f);
            // Structural criticality (bridges): low=stable, high=fragile.
            r = lerpF(0.18f, 1.00f, t);
            g = lerpF(0.55f, 0.35f, t);
            b = lerpF(0.95f, 0.12f, t);
          }
          return rgba(r, g, b, alpha);
        };

        // Draw lane network (LOD).
        const bool trafficMode = (st.hyperlaneColorMode == 2);
        const float trafficWidthBoost = std::clamp(st.hyperlaneTrafficWidthBoost, 0.05f, 10.0f);
        const float trafficAlphaBoost = std::clamp(st.hyperlaneTrafficAlphaBoost, 0.0f, 10.0f);

        const bool criticalMode = (st.hyperlaneColorMode == 3);
        const float criticalWidthBoost = std::clamp(st.hyperlaneCriticalWidthBoost, 0.05f, 10.0f);
        const float criticalAlphaBoost = std::clamp(st.hyperlaneCriticalAlphaBoost, 0.0f, 10.0f);

        std::size_t drawn = 0;
        for (std::size_t ei = 0; ei < st.hyperlaneNet.edges.size(); ++ei) {
          if (drawn >= maxEdgesDraw) break;
          const auto& e = st.hyperlaneNet.edges[ei];
          auto itA = hyperNodePos.find(e.a);
          auto itB = hyperNodePos.find(e.b);
          if (itA == hyperNodePos.end() || itB == hyperNodePos.end()) continue;

          const ImVec2 a = itA->second;
          const ImVec2 b2 = itB->second;
          const float dx = b2.x - a.x;
          const float dy = b2.y - a.y;
          const float len = std::sqrt(dx * dx + dy * dy);
          if (len < minLenPx) continue;

          const float bw = (float)std::clamp(e.bandwidth01, 0.0, 1.0);
          const float risk = (float)std::clamp(e.risk01, 0.0, 1.0);

          float traffic01 = 0.0f;
          if (trafficMode && ei < st.hyperlaneEdgeBetweenness.size() && st.hyperlaneCentralityMaxEdge > 0.0) {
            traffic01 = (float)(st.hyperlaneEdgeBetweenness[ei] / st.hyperlaneCentralityMaxEdge);
            traffic01 = std::clamp(traffic01, 0.0f, 1.0f);
          }

          bool isBridge = false;
          float critical01 = 0.0f;
          if (criticalMode && ei < st.hyperlaneVulnerability.edgeIsBridge.size() && st.hyperlaneVulnerability.edgeIsBridge[ei]) {
            isBridge = true;
            if (st.hyperlaneVulnerability.maxEdgeCut01 > 0.0f && ei < st.hyperlaneVulnerability.edgeCut01.size()) {
              critical01 = st.hyperlaneVulnerability.edgeCut01[ei] / st.hyperlaneVulnerability.maxEdgeCut01;
            } else if (st.hyperlaneVulnerability.maxEdgeCutSize > 0 && ei < st.hyperlaneVulnerability.edgeCutSize.size()) {
              critical01 = (float)st.hyperlaneVulnerability.edgeCutSize[ei] / (float)st.hyperlaneVulnerability.maxEdgeCutSize;
            }
            critical01 = std::clamp(critical01, 0.0f, 1.0f);
          }

          const float lenFade = smoothstep(minLenPx, minLenPx * 4.0f + 1.0f, len);

          float vis01 = bw;
          if (trafficMode) vis01 = traffic01;
          if (criticalMode) {
            if (isBridge) {
              vis01 = std::max(vis01 * 0.25f, 0.35f + 0.65f * critical01);
            } else {
              vis01 = vis01 * 0.25f;
            }
          }

          float alpha = baseAlpha * (0.25f + 0.75f * vis01) * lenFade;
          if (trafficMode) alpha *= trafficAlphaBoost;
          if (criticalMode) alpha *= criticalAlphaBoost;
          alpha = std::clamp(alpha, 0.0f, 1.0f);
          if (alpha <= 0.001f) continue;

          float widthMul = (0.75f + 1.75f * bw);
          if (trafficMode) widthMul = (0.65f + 2.35f * traffic01) * trafficWidthBoost;
          if (criticalMode && st.hyperlaneCriticalHighlightBridges && isBridge) {
            widthMul = (0.75f + 2.75f * critical01) * criticalWidthBoost;
          }
          const float width = baseWidth * widthMul;

          dl->AddLine(a, b2, laneColor(bw, risk, traffic01, critical01, alpha), width);
          drawn++;
        }

        // In traffic mode, optionally redraw the highest-betweenness edges on top so
        // chokepoints remain visible even under LOD edge caps.
        if (trafficMode && st.hyperlaneTrafficHighlightChokepoints && !st.hyperlaneChokepointEdgeIdx.empty()) {
          const float wBoost = trafficWidthBoost;
          const float aBoost = trafficAlphaBoost;
          for (int idx : st.hyperlaneChokepointEdgeIdx) {
            if (idx < 0 || idx >= (int)st.hyperlaneNet.edges.size()) continue;
            const auto& e = st.hyperlaneNet.edges[(std::size_t)idx];
            auto itA = hyperNodePos.find(e.a);
            auto itB = hyperNodePos.find(e.b);
            if (itA == hyperNodePos.end() || itB == hyperNodePos.end()) continue;

            const ImVec2 a = itA->second;
            const ImVec2 b2 = itB->second;
            const float dx = b2.x - a.x;
            const float dy = b2.y - a.y;
            const float len = std::sqrt(dx * dx + dy * dy);
            if (len < minLenPx) continue;
            const float lenFade = smoothstep(minLenPx, minLenPx * 4.0f + 1.0f, len);
            if (lenFade <= 0.0f) continue;

            const float bw = (float)std::clamp(e.bandwidth01, 0.0, 1.0);
            const float risk = (float)std::clamp(e.risk01, 0.0, 1.0);

            float traffic01 = 0.0f;
            if ((std::size_t)idx < st.hyperlaneEdgeBetweenness.size() && st.hyperlaneCentralityMaxEdge > 0.0) {
              traffic01 = (float)(st.hyperlaneEdgeBetweenness[(std::size_t)idx] / st.hyperlaneCentralityMaxEdge);
              traffic01 = std::clamp(traffic01, 0.0f, 1.0f);
            }

            const float alpha = std::clamp(baseAlpha * (0.45f + 0.55f * traffic01) * lenFade * 1.25f * aBoost, 0.0f, 1.0f);
            const float width = baseWidth * (1.25f + 3.50f * traffic01) * wBoost;
            dl->AddLine(a, b2, laneColor(bw, risk, traffic01, 0.0f, alpha), width);
          }
        }

        // In criticality mode, optionally redraw the most fragile bridges on top so
        // they remain visible even under LOD edge caps.
        if (criticalMode && st.hyperlaneCriticalHighlightBridges && !st.hyperlaneCriticalBridgeEdgeIdx.empty()) {
          const float wBoost = criticalWidthBoost;
          const float aBoost = criticalAlphaBoost;
          for (int idx : st.hyperlaneCriticalBridgeEdgeIdx) {
            if (idx < 0 || idx >= (int)st.hyperlaneNet.edges.size()) continue;
            const auto& e = st.hyperlaneNet.edges[(std::size_t)idx];
            auto itA = hyperNodePos.find(e.a);
            auto itB = hyperNodePos.find(e.b);
            if (itA == hyperNodePos.end() || itB == hyperNodePos.end()) continue;

            const ImVec2 a = itA->second;
            const ImVec2 b2 = itB->second;
            const float dx = b2.x - a.x;
            const float dy = b2.y - a.y;
            const float len = std::sqrt(dx * dx + dy * dy);
            if (len < minLenPx) continue;
            const float lenFade = smoothstep(minLenPx, minLenPx * 4.0f + 1.0f, len);
            if (lenFade <= 0.0f) continue;

            const float bw = (float)std::clamp(e.bandwidth01, 0.0, 1.0);
            const float risk = (float)std::clamp(e.risk01, 0.0, 1.0);

            float critical01 = 0.0f;
            if ((std::size_t)idx < st.hyperlaneVulnerability.edgeIsBridge.size() && st.hyperlaneVulnerability.edgeIsBridge[(std::size_t)idx]) {
              if (st.hyperlaneVulnerability.maxEdgeCut01 > 0.0f && (std::size_t)idx < st.hyperlaneVulnerability.edgeCut01.size()) {
                critical01 = st.hyperlaneVulnerability.edgeCut01[(std::size_t)idx] / st.hyperlaneVulnerability.maxEdgeCut01;
              } else if (st.hyperlaneVulnerability.maxEdgeCutSize > 0 && (std::size_t)idx < st.hyperlaneVulnerability.edgeCutSize.size()) {
                critical01 = (float)st.hyperlaneVulnerability.edgeCutSize[(std::size_t)idx] / (float)st.hyperlaneVulnerability.maxEdgeCutSize;
              }
              critical01 = std::clamp(critical01, 0.0f, 1.0f);
            }

            const float alpha = std::clamp(baseAlpha * (0.45f + 0.55f * critical01) * lenFade * 1.25f * aBoost, 0.0f, 1.0f);
            const float width = baseWidth * (1.25f + 3.50f * critical01) * wBoost;
            dl->AddLine(a, b2, laneColor(bw, risk, 0.0f, critical01, alpha), width);
          }
        }

        if (st.hyperlaneMaxEdgesDraw > 0 && st.hyperlaneNet.edges.size() > maxEdgesDraw) {
          dl->AddText(ImVec2(p0.x + 10.0f, p0.y + 26.0f), rgba(0.85f, 0.85f, 0.85f, 0.70f), "Hyperlanes LOD: drawing capped");
        }

        // Route overlay (when enabled).
        if (st.hyperlaneRouteEnabled && st.hyperlaneRouteDraw) {
          const float routeAlpha = std::clamp(st.hyperlaneRouteOpacity, 0.0f, 1.0f);
          const float routeWidthBase = std::clamp(st.hyperlaneRouteWidthPx, 1.0f, 16.0f);

          if (!st.hyperlaneRoutePath.empty() && routeAlpha > 0.001f) {
            // Build a quick lookup so per-segment coloring can use edge risk/bw.
            struct EdgeKey {
              sim::SystemId a{};
              sim::SystemId b{};
              bool operator==(const EdgeKey& o) const noexcept { return a == o.a && b == o.b; }
            };
            struct EdgeKeyHash {
              std::size_t operator()(const EdgeKey& k) const noexcept {
                core::u64 h = core::fnv1a64("hyperlaneEdgeKey");
                h = hashMix(h, (core::u64)k.a);
                h = hashMix(h, (core::u64)k.b);
                return (std::size_t)h;
              }
            };

            auto normKey = [](sim::SystemId x, sim::SystemId y) -> EdgeKey {
              if ((core::u64)x < (core::u64)y) return EdgeKey{x, y};
              return EdgeKey{y, x};
            };

            struct EdgeInfo {
              float bw{1.0f};
              float risk{0.0f};
              float traffic{0.0f};
              float critical{0.0f};
            };

            std::unordered_map<EdgeKey, EdgeInfo, EdgeKeyHash> lut;
            lut.reserve(st.hyperlaneNet.edges.size() * 2);
            for (std::size_t ei = 0; ei < st.hyperlaneNet.edges.size(); ++ei) {
              const auto& e = st.hyperlaneNet.edges[ei];
              const float bw = (float)std::clamp(e.bandwidth01, 0.0, 1.0);
              const float risk = (float)std::clamp(e.risk01, 0.0, 1.0);
              float traffic01 = 0.0f;
              if (st.hyperlaneColorMode == 2 && ei < st.hyperlaneEdgeBetweenness.size() && st.hyperlaneCentralityMaxEdge > 0.0) {
                traffic01 = (float)(st.hyperlaneEdgeBetweenness[ei] / st.hyperlaneCentralityMaxEdge);
                traffic01 = std::clamp(traffic01, 0.0f, 1.0f);
              }

              float critical01 = 0.0f;
              if (st.hyperlaneColorMode == 3 && ei < st.hyperlaneVulnerability.edgeIsBridge.size() && st.hyperlaneVulnerability.edgeIsBridge[ei]) {
                if (st.hyperlaneVulnerability.maxEdgeCut01 > 0.0f && ei < st.hyperlaneVulnerability.edgeCut01.size()) {
                  critical01 = st.hyperlaneVulnerability.edgeCut01[ei] / st.hyperlaneVulnerability.maxEdgeCut01;
                } else if (st.hyperlaneVulnerability.maxEdgeCutSize > 0 && ei < st.hyperlaneVulnerability.edgeCutSize.size()) {
                  critical01 = (float)st.hyperlaneVulnerability.edgeCutSize[ei] / (float)st.hyperlaneVulnerability.maxEdgeCutSize;
                }
                critical01 = std::clamp(critical01, 0.0f, 1.0f);
              }

              lut.emplace(normKey(e.a, e.b), EdgeInfo{bw, risk, traffic01, critical01});
            }

            const auto routeColor = [&](float bw01, float risk01, float traffic01, float critical01, float alpha) -> ImU32 {
              float r = 0.25f, g = 0.95f, b = 1.00f;
              const int mode = st.hyperlaneColorMode;
              if (mode == 0) {
                const float t = std::clamp(risk01, 0.0f, 1.0f);
                r = lerpF(0.12f, 1.00f, t);
                g = lerpF(0.95f, 0.25f, t);
                b = lerpF(1.00f, 0.10f, t);
              } else if (mode == 1) {
                const float u = std::clamp(bw01, 0.0f, 1.0f);
                r = lerpF(0.20f, 0.95f, u);
                g = lerpF(0.25f, 0.95f, u);
                b = lerpF(0.30f, 0.95f, u);
              } else if (mode == 2) {
                const float t = std::clamp(traffic01, 0.0f, 1.0f);
                r = lerpF(0.20f, 1.00f, t);
                g = lerpF(0.95f, 0.20f, t);
                b = lerpF(1.00f, 0.10f, t);
              } else {
                const float t = std::clamp(critical01, 0.0f, 1.0f);
                r = lerpF(0.18f, 1.00f, t);
                g = lerpF(0.55f, 0.35f, t);
                b = lerpF(0.95f, 0.12f, t);
              }

              // Pull toward white so it is clearly a highlight layer.
              const float w = 0.35f + 0.45f * std::clamp(bw01, 0.0f, 1.0f);
              r = lerpF(r, 1.0f, w);
              g = lerpF(g, 1.0f, w);
              b = lerpF(b, 1.0f, w);
              return rgba(r, g, b, alpha);
            };
            const std::size_t segCap = 12000;
            std::size_t segDrawn = 0;

            const auto drawPath = [&](const std::vector<sim::SystemId>& path, float alphaMul, float widthMul, bool markers) {
              if (path.empty()) return;

              for (std::size_t i = 1; i < path.size(); ++i) {
                if (segDrawn >= segCap) break;

                const sim::SystemId aId = path[i - 1];
                const sim::SystemId bId = path[i];

                auto itA = hyperNodePos.find(aId);
                auto itB = hyperNodePos.find(bId);
                if (itA == hyperNodePos.end() || itB == hyperNodePos.end()) continue;

                float bw = 1.0f;
                float risk = 0.0f;
                float traffic01 = 0.0f;
                float critical01 = 0.0f;
                if (auto itE = lut.find(normKey(aId, bId)); itE != lut.end()) {
                  bw = itE->second.bw;
                  risk = itE->second.risk;
                  traffic01 = itE->second.traffic;
                  critical01 = itE->second.critical;
                }

                float width = routeWidthBase * std::max(0.05f, widthMul) * (0.85f + 1.65f * bw);
                if (st.hyperlaneColorMode == 3) {
                  width *= (0.85f + 1.75f * std::clamp(critical01, 0.0f, 1.0f));
                }
                dl->AddLine(itA->second,
                            itB->second,
                            routeColor(bw, risk, traffic01, critical01, routeAlpha * std::clamp(alphaMul, 0.0f, 1.0f)),
                            width);
                segDrawn++;
              }

              if (markers) {
                // Start/end markers.
                auto itS = hyperNodePos.find(path.front());
                auto itT = hyperNodePos.find(path.back());
                if (itS != hyperNodePos.end()) {
                  dl->AddCircleFilled(itS->second, routeWidthBase * 1.25f + 3.0f, rgba(0.95f, 0.95f, 0.95f, routeAlpha));
                }
                if (itT != hyperNodePos.end()) {
                  dl->AddCircleFilled(itT->second, routeWidthBase * 1.25f + 3.0f, rgba(1.00f, 0.75f, 0.25f, routeAlpha));
                }
              }
            };

            // Draw alternative routes first so the selected route sits on top.
            if (st.hyperlaneRouteDrawAlternatives && st.hyperlaneKRoutes.size() > 1) {
              const float altA = std::clamp(st.hyperlaneRouteAltOpacity, 0.0f, 1.0f);
              const float altW = std::max(0.05f, st.hyperlaneRouteAltWidthScale);
              const int sel = st.hyperlaneRouteSelected;
              const int count = (int)st.hyperlaneKRoutes.size();

              for (int ri = 0; ri < count; ++ri) {
                if (ri == sel) continue;
                const auto& p = st.hyperlaneKRoutes[(std::size_t)ri].path;
                if (p.size() < 2) continue;

                // A gentle fade so later alternatives don't dominate the scene.
                const float fade = std::max(0.15f, 1.0f - 0.10f * (float)ri);
                drawPath(p, altA * fade, altW, false);
                if (segDrawn >= segCap) break;
              }
            }

            // Selected route (bold + markers).
            drawPath(st.hyperlaneRoutePath, 1.0f, 1.0f, true);

            if (segDrawn >= segCap) {
              dl->AddText(ImVec2(p0.x + 10.0f, p0.y + 44.0f), rgba(0.85f, 0.85f, 0.85f, 0.70f), "Route LOD: too many segments, drawing capped");
            }
          } else if (!st.hyperlaneRouteError.empty()) {
            dl->AddText(ImVec2(p0.x + 10.0f, p0.y + 44.0f), rgba(1.0f, 0.40f, 0.40f, 0.85f), st.hyperlaneRouteError.c_str());
          }
        }
      }
    }

    // Draw points.
    for (std::size_t i = 0; i < st.stubs.size(); ++i) {
      const auto& stub = st.stubs[i];

      const ImVec2 p = worldToScreen(stub.posLy);
      if (p.x < p0.x || p.x > p1.x || p.y < p0.y || p.y > p1.y) continue;

      ImU32 col = colorForStarClass(stub.primaryClass);
      if (st.colorByFaction) {
        col = colorForFaction(stub.factionId, st.factions);
      } else if (st.colorByRegion && i < st.stubRegionKind.size()) {
        col = colorForRegionKind(st.stubRegionKind[i]);
      } else if (st.colorByCluster && i < st.stubCluster01.size()) {
        col = colorForClusterInfluence(st.stubCluster01[i]);
      } else if (st.colorByVoid && i < st.stubVoid01.size()) {
        col = colorForVoidInfluence(st.stubVoid01[i]);
      } else if (st.colorByMorphology && i < st.stubMorph01.size()) {
        col = colorForMorphologyInfluence(st.stubMorph01[i]);
      }
      dl->AddCircleFilled(p, 2.0f, col);
    }

    // Criticality overlay: articulation rings (drawn on top of points).
    if (st.showHyperlanes && st.hyperlaneColorMode == 3 && st.hyperlaneCriticalHighlightArticulation &&
        !st.hyperlaneCriticalArticulationNodeIds.empty()) {
      ensureHyperNodePos();

      const float ringR = std::clamp(st.hyperlaneCriticalNodeRingRadiusPx, 2.0f, 64.0f);
      const float ringTh = std::clamp(st.hyperlaneCriticalNodeRingThicknessPx, 0.5f, 16.0f);
      const float max01 = std::max(0.0001f, st.hyperlaneVulnerability.maxNodeImpact01);

      for (std::size_t i = 0; i < st.hyperlaneCriticalArticulationNodeIds.size(); ++i) {
        const sim::SystemId id = st.hyperlaneCriticalArticulationNodeIds[i];
        auto it = hyperNodePos.find(id);
        if (it == hyperNodePos.end()) continue;

        float impact01 = 0.0f;
        if (i < st.hyperlaneCriticalArticulationNodeImpact01.size()) {
          impact01 = st.hyperlaneCriticalArticulationNodeImpact01[i];
        }
        const float t = std::clamp(impact01 / max01, 0.0f, 1.0f);

        const float r = lerpF(0.18f, 1.00f, t);
        const float g = lerpF(0.55f, 0.35f, t);
        const float b = lerpF(0.95f, 0.12f, t);
        const float a = std::clamp(0.35f + 0.55f * t, 0.0f, 1.0f);

        const float rr = ringR + 6.0f * t;
        dl->AddCircle(it->second, rr, rgba(r, g, b, a), 0, ringTh);
      }
    }

    // Hover tooltip.
    if (ImGui::IsMouseHoveringRect(p0, p1) && !st.stubs.empty()) {
      const ImVec2 m = ImGui::GetIO().MousePos;
      int bestIdx = -1;
      float bestD2 = 6.0f * 6.0f;

      for (int i = 0; i < static_cast<int>(st.stubs.size()); ++i) {
        const ImVec2 p = worldToScreen(st.stubs[static_cast<std::size_t>(i)].posLy);
        const float dx = m.x - p.x;
        const float dy = m.y - p.y;
        const float d2 = dx * dx + dy * dy;
        if (d2 < bestD2) {
          bestD2 = d2;
          bestIdx = i;
        }
      }

      if (bestIdx >= 0) {
        const auto& stub = st.stubs[static_cast<std::size_t>(bestIdx)];
        ImGui::BeginTooltip();
        ImGui::TextUnformatted(stub.name.c_str());
        ImGui::Separator();
        ImGui::Text("SystemId: 0x%llx", static_cast<unsigned long long>(stub.id));
        ImGui::Text("Pos (ly): [%.1f, %.1f, %.1f]", stub.posLy.x, stub.posLy.y, stub.posLy.z);
        ImGui::Text("Star: %c", static_cast<char>(stub.primaryClass));
        ImGui::Text("Planets: %u", stub.planetCount);
        ImGui::Text("Stations: %u", stub.stationCount);
        ImGui::Text("Faction: %u", stub.factionId);

        if (st.colorByRegion) {
          const auto reg = proc::sampleGalaxyRegion(st.seed, stub.posLy, st.regionCellSizeLy);
          ImGui::Separator();
          ImGui::Text("Region: %s", reg.name.c_str());
          ImGui::Text("Kind: %s", proc::galaxyRegionKindName(reg.kind));
          ImGui::Text("Edge: %.2f", reg.edge01);
        }
        if (st.colorByCluster || st.params.clusterStrength > 0.0) {
          proc::GalaxyClustersParams cp{};
          cp.cellSizeLy = st.params.clusterCellSizeLy;
          cp.chancePerCell = st.params.clusterChancePerCell;
          cp.radiusLy = st.params.clusterRadiusLy;
          cp.radiusJitter01 = st.params.clusterRadiusJitter;
          cp.strengthJitter01 = st.params.clusterStrengthJitter;
          cp.falloffPower = st.params.clusterFalloffPower;

          const auto cs = proc::sampleGalaxyClusters(st.seed, stub.posLy, cp);
          ImGui::Separator();
          ImGui::Text("Cluster influence: %.2f", cs.cluster01);
          if (cs.hasCluster) {
            ImGui::Text("ClusterId: 0x%llx", static_cast<unsigned long long>(cs.clusterId));
            ImGui::Text("Radius: %.0f ly", cs.radiusLy);
          }
        }

        if (st.colorByVoid || st.params.voidStrength > 0.0) {
          proc::GalaxyVoidsParams vp{};
          vp.cellSizeLy = st.params.voidCellSizeLy;
          vp.chancePerCell = st.params.voidChancePerCell;
          vp.radiusLy = st.params.voidRadiusLy;
          vp.radiusJitter01 = st.params.voidRadiusJitter;
          vp.strengthJitter01 = st.params.voidStrengthJitter;
          vp.falloffPower = st.params.voidFalloffPower;

          const auto vs = proc::sampleGalaxyVoids(st.seed, stub.posLy, vp);
          ImGui::Separator();
          ImGui::Text("Void influence: %.2f", vs.void01);
          if (vs.hasVoid) {
            ImGui::Text("VoidId: 0x%llx", static_cast<unsigned long long>(vs.voidId));
            ImGui::Text("Radius: %.0f ly", vs.radiusLy);
          }
        }

        if (st.colorByMorphology || st.params.barStrength > 0.0 || st.params.ringStrength > 0.0 ||
            std::abs(st.params.warpAmplitudeLy) > 0.0 || st.params.flareStrength > 0.0) {
          const auto ms = proc::sampleGalaxyMorphology(st.seed, st.params, stub.posLy);
          ImGui::Separator();
          ImGui::Text("Morph mul: %.2f", ms.densityMul);
          ImGui::Text("Bar01: %.2f  Ring01: %.2f", ms.bar01, ms.ring01);
          if (std::abs(st.params.warpAmplitudeLy) > 0.0 || st.params.flareStrength > 0.0) {
            ImGui::Text("Warp Z: %.0f ly", ms.warpZLy);
            ImGui::Text("Half-thickness: %.0f ly", ms.thicknessHalfLy);
            ImGui::Text("Zrel: %.0f ly", std::abs(stub.posLy.z - ms.warpZLy));
          }
        }


        // Hyperlane stats for this system (if the overlay is enabled).
        if (st.showHyperlanes && !st.hyperlaneNet.edges.empty()) {
          const sim::SystemId hid = stub.id;
          int deg = 0;
          double bwSum = 0.0;
          double riskSum = 0.0;
          for (const auto& e : st.hyperlaneNet.edges) {
            if (e.a == hid || e.b == hid) {
              deg++;
              bwSum += e.bandwidth01;
              riskSum += e.risk01;
            }
          }
          ImGui::Separator();
          if (deg > 0) {
            ImGui::Text("Hyperlanes: %d", deg);
            ImGui::Text("Avg BW: %.2f  Avg Risk: %.2f", (float)(bwSum / (double)deg), (float)(riskSum / (double)deg));
            if (st.hyperlaneColorMode == 2 && !st.hyperlaneCentralityNodeIds.empty() &&
                st.hyperlaneCentralityMaxNode > 0.0 && st.hyperlaneCentralityNodeIds.size() == st.hyperlaneNodeBetweenness.size()) {
              auto it = std::lower_bound(st.hyperlaneCentralityNodeIds.begin(), st.hyperlaneCentralityNodeIds.end(), hid);
              if (it != st.hyperlaneCentralityNodeIds.end() && *it == hid) {
                const std::size_t idx = (std::size_t)std::distance(st.hyperlaneCentralityNodeIds.begin(), it);
                const double raw = st.hyperlaneNodeBetweenness[idx];
                const float norm = (float)std::clamp(raw / st.hyperlaneCentralityMaxNode, 0.0, 1.0);
                ImGui::Text("Traffic: %.2f", norm);
              }
            }
          } else {
            ImGui::TextDisabled("Hyperlanes: none (not in lane subset)");
          }
        }

        ImGui::TextDisabled("Click to copy name");
        if (st.showHyperlanes && st.hyperlaneRouteEnabled) {
          ImGui::TextDisabled("Shift+Click: route start   Ctrl+Click: route end   Right-click: clear route");
        }
        ImGui::EndTooltip();

        // Optional highlight: redraw hyperlanes incident to the hovered system on top of points.
        if (st.showHyperlanes && st.hyperlaneHighlightHovered && !st.hyperlaneNet.edges.empty() && !st.hyperlaneNodes.empty()) {
          const sim::SystemId hid = stub.id;
          ensureHyperNodePos();
          if (hyperNodePos.find(hid) != hyperNodePos.end()) {
            for (const auto& e : st.hyperlaneNet.edges) {
              if (e.a != hid && e.b != hid) continue;
              auto itA = hyperNodePos.find(e.a);
              auto itB = hyperNodePos.find(e.b);
              if (itA == hyperNodePos.end() || itB == hyperNodePos.end()) continue;
              const float bw = (float)std::clamp(e.bandwidth01, 0.0, 1.0);
              const float alpha = std::clamp(st.hyperlaneOpacity * (0.75f + 0.75f * bw), 0.0f, 1.0f);
              const float width = std::clamp(st.hyperlaneWidthPx * (1.40f + 2.25f * bw), 0.5f, 10.0f);
              dl->AddLine(itA->second, itB->second, rgba(1.0f, 1.0f, 1.0f, alpha), width);
            }
          }
        }

        const bool routeGestures = (st.showHyperlanes && st.hyperlaneRouteEnabled);
        const ImGuiIO& io = ImGui::GetIO();

        if (ImGui::IsMouseClicked(1) && routeGestures) {
          st.hyperlaneRouteFrom = 0;
          st.hyperlaneRouteTo = 0;
          st.hyperlaneRouteDirty = true;
          if (toast) toast("Route cleared", 1.5);
        } else if (ImGui::IsMouseClicked(0)) {
          if (routeGestures && io.KeyShift) {
            st.hyperlaneRouteFrom = stub.id;
            st.hyperlaneRouteDirty = true;
            if (toast) toast("Route start set", 1.5);
          } else if (routeGestures && io.KeyCtrl) {
            st.hyperlaneRouteTo = stub.id;
            st.hyperlaneRouteDirty = true;
            if (toast) toast("Route end set", 1.5);
          } else {
            ImGui::SetClipboardText(stub.name.c_str());
            if (toast) toast("Copied system name to clipboard", 1.5);
          }
        }
      }
    }

    // Legend.
    if (st.showLegend) {
      ImGui::SetCursorScreenPos(ImVec2(p0.x + 12.0f, p0.y + 12.0f));
      ImGui::BeginChild("##legend", ImVec2(210, 160), true);

      const bool legendFaction = st.colorByFaction;
      const bool legendRegion = (!legendFaction && st.colorByRegion);
      const bool legendCluster = (!legendFaction && !legendRegion && st.colorByCluster);
      const bool legendVoid = (!legendFaction && !legendRegion && !legendCluster && st.colorByVoid);
      const bool legendMorph = (!legendFaction && !legendRegion && !legendCluster && !legendVoid && st.colorByMorphology);

      const char* legendTitle = legendFaction        ? "Legend (Faction)"
                                : legendRegion     ? "Legend (Region)"
                                : legendCluster    ? "Legend (Cluster Influence)"
                                : legendVoid       ? "Legend (Void Influence)"
                                : legendMorph      ? "Legend (Morph Influence)"
                                                   : "Legend (Star Class)";
      ImGui::TextUnformatted(legendTitle);
      ImGui::Separator();

      if (legendRegion) {
        const proc::GalaxyRegionKind kinds[] = {
            proc::GalaxyRegionKind::Core,
            proc::GalaxyRegionKind::InnerDisc,
            proc::GalaxyRegionKind::OuterRim,
            proc::GalaxyRegionKind::Nebula,
            proc::GalaxyRegionKind::Cluster,
            proc::GalaxyRegionKind::Rift,
        };
        for (auto k : kinds) {
          // Unique IDs: ImGui uses the label as the ID. These legend swatches
          // are purely decorative but must not collide.
          ImGui::PushID((int)k);
          ImGui::ColorButton("##r", ImGui::ColorConvertU32ToFloat4(colorForRegionKind(k)), ImGuiColorEditFlags_NoTooltip, ImVec2(16, 16));
          ImGui::PopID();
          ImGui::SameLine();
          ImGui::TextUnformatted(proc::galaxyRegionKindName(k));
        }
      } else if (legendCluster) {
        const float levels[] = {0.0f, 0.25f, 0.50f, 0.75f, 1.0f};
        for (int i = 0; i < 5; ++i) {
          const float t = levels[i];
          ImGui::PushID(i);
          ImGui::ColorButton("##cl", ImGui::ColorConvertU32ToFloat4(colorForClusterInfluence(t)), ImGuiColorEditFlags_NoTooltip, ImVec2(16, 16));
          ImGui::PopID();
          ImGui::SameLine();
          ImGui::Text("%.2f", t);
        }
      } else if (legendVoid) {
        const float levels[] = {0.0f, 0.25f, 0.50f, 0.75f, 1.0f};
        for (int i = 0; i < 5; ++i) {
          const float t = levels[i];
          ImGui::PushID(i);
          ImGui::ColorButton("##vd", ImGui::ColorConvertU32ToFloat4(colorForVoidInfluence(t)), ImGuiColorEditFlags_NoTooltip, ImVec2(16, 16));
          ImGui::PopID();
          ImGui::SameLine();
          ImGui::Text("%.2f", t);
        }
      } else if (legendMorph) {
        const float levels[] = {0.0f, 0.25f, 0.50f, 0.75f, 1.0f};
        for (int i = 0; i < 5; ++i) {
          const float t = levels[i];
          ImGui::PushID(i);
          ImGui::ColorButton("##m", ImGui::ColorConvertU32ToFloat4(colorForMorphologyInfluence(t)), ImGuiColorEditFlags_NoTooltip, ImVec2(16, 16));
          ImGui::PopID();
          ImGui::SameLine();
          ImGui::Text("%.2f", t);
        }
      } else if (!legendFaction) {
        const sim::StarClass classes[] = {sim::StarClass::O, sim::StarClass::B, sim::StarClass::A, sim::StarClass::F, sim::StarClass::G, sim::StarClass::K, sim::StarClass::M};
        for (auto c : classes) {
          ImGui::PushID((int)c);
          ImGui::ColorButton("##c", ImGui::ColorConvertU32ToFloat4(colorForStarClass(c)), ImGuiColorEditFlags_NoTooltip, ImVec2(16, 16));
          ImGui::PopID();
          ImGui::SameLine();
          ImGui::Text("%c", static_cast<char>(c));
        }
      } else {
        for (const auto& f : st.factions) {
          ImGui::PushID((int)f.id);
          ImGui::ColorButton("##f", ImVec4(static_cast<float>(f.colorRgb.x), static_cast<float>(f.colorRgb.y), static_cast<float>(f.colorRgb.z), 1.0f), ImGuiColorEditFlags_NoTooltip, ImVec2(16, 16));
          ImGui::PopID();
          ImGui::SameLine();
          ImGui::Text("%u %s", f.id, f.name.c_str());
        }
      }

      ImGui::EndChild();
    }

    ImGui::EndChild();

    ImGui::EndTable();
  }

  ImGui::End();
}

} // namespace stellar::game
