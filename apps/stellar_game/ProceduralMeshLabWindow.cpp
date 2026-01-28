#include "ProceduralMeshLabWindow.h"

#include "stellar/core/Log.h"
#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/math/Math.h"
#include "stellar/math/Mat4.h"
#include "stellar/math/Vec3.h"
#include "stellar/render/GltfExport.h"
#include "stellar/render/Gl.h"

#include "Screenshot.h"

#include <SDL_opengl.h>
#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_set>

namespace stellar::game {

namespace {

using namespace stellar;

// -----------------------------------------------------------------------------
// Undo / Redo history helpers
// -----------------------------------------------------------------------------
using HistorySnapshot = ProceduralMeshLabWindowState::HistorySnapshot;

static bool sdfNodeEquals(const render::SdfNode& a, const render::SdfNode& b) {
  if (a.op != b.op) return false;
  if (a.in0 != b.in0) return false;
  if (a.in1 != b.in1) return false;
  if (a.seed != b.seed) return false;
  if (a.p0 != b.p0) return false;
  if (a.p1 != b.p1) return false;
  if (a.p2 != b.p2) return false;
  if (a.p3 != b.p3) return false;
  if (a.p4 != b.p4) return false;
  if (a.p5 != b.p5) return false;
  if (a.p6 != b.p6) return false;
  if (a.p7 != b.p7) return false;
  return true;
}

static bool sdfGraphEquals(const render::SdfGraph& a, const render::SdfGraph& b) {
  if (a.seed != b.seed) return false;
  if (a.output != b.output) return false;
  if (a.nodes.size() != b.nodes.size()) return false;
  for (std::size_t i = 0; i < a.nodes.size(); ++i) {
    if (!sdfNodeEquals(a.nodes[i], b.nodes[i])) return false;
  }
  return true;
}

static bool procNodeEquals(const render::ProcNode& a, const render::ProcNode& b) {
  if (a.op != b.op) return false;
  if (a.in0 != b.in0) return false;
  if (a.in1 != b.in1) return false;
  if (a.in2 != b.in2) return false;
  if (a.seed != b.seed) return false;
  if (a.p0 != b.p0) return false;
  if (a.p1 != b.p1) return false;
  if (a.p2 != b.p2) return false;
  if (a.p3 != b.p3) return false;
  return true;
}

static bool procPaletteStopEquals(const render::ProcPaletteStop& a, const render::ProcPaletteStop& b) {
  if (a.pos != b.pos) return false;
  if (a.r != b.r) return false;
  if (a.g != b.g) return false;
  if (a.b != b.b) return false;
  return true;
}

static bool procGraphEquals(const render::ProcGraph& a, const render::ProcGraph& b) {
  if (a.seed != b.seed) return false;
  if (a.output != b.output) return false;
  if (a.usePalette != b.usePalette) return false;
  if (a.paletteCount != b.paletteCount) return false;
  for (int i = 0; i < a.paletteCount; ++i) {
    if (!procPaletteStopEquals(a.palette[i], b.palette[i])) return false;
  }
  if (a.nodes.size() != b.nodes.size()) return false;
  for (std::size_t i = 0; i < a.nodes.size(); ++i) {
    if (!procNodeEquals(a.nodes[i], b.nodes[i])) return false;
  }
  return true;
}

static bool sdfMesherParamsEquals(const render::SdfMesherParams& a, const render::SdfMesherParams& b) {
  return a.resolution == b.resolution &&
         a.bounds == b.bounds &&
         a.iso == b.iso &&
         a.computeNormalsFromField == b.computeNormalsFromField &&
         a.normalEps == b.normalEps &&
         a.fixWindingFromNormals == b.fixWindingFromNormals;
}

static bool raymarchSettingsEquals(const render::SdfRaymarchSettings& a, const render::SdfRaymarchSettings& b) {
  if (a.maxSteps != b.maxSteps) return false;
  if (a.epsilon != b.epsilon) return false;
  if (a.maxDistance != b.maxDistance) return false;
  if (a.useBoundsAabb != b.useBoundsAabb) return false;

  for (int i = 0; i < 3; ++i) {
    if (a.lightDir[i] != b.lightDir[i]) return false;
    if (a.baseColor[i] != b.baseColor[i]) return false;
    if (a.backgroundColor[i] != b.backgroundColor[i]) return false;
  }

  if (a.ambient != b.ambient) return false;
  if (a.diffuse != b.diffuse) return false;
  if (a.specular != b.specular) return false;
  if (a.shininess != b.shininess) return false;

  if (a.softShadows != b.softShadows) return false;
  if (a.shadowSteps != b.shadowSteps) return false;
  if (a.shadowMaxDistance != b.shadowMaxDistance) return false;
  if (a.shadowK != b.shadowK) return false;

  if (a.ambientOcclusion != b.ambientOcclusion) return false;
  if (a.aoSteps != b.aoSteps) return false;
  if (a.aoStepSize != b.aoStepSize) return false;
  if (a.aoStrength != b.aoStrength) return false;

  if (a.debug != b.debug) return false;
  return true;
}

static bool historySnapshotEquals(const HistorySnapshot& a, const HistorySnapshot& b) {
  if (a.preset != b.preset) return false;
  if (a.lockToPreset != b.lockToPreset) return false;
  if (a.seed != b.seed) return false;

  if (a.mesherType != b.mesherType) return false;
  if (!sdfMesherParamsEquals(a.mesher, b.mesher)) return false;
  if (a.dcQefRegularization != b.dcQefRegularization) return false;
  if (a.dcClampToCell != b.dcClampToCell) return false;
  if (a.dcProjectToIso != b.dcProjectToIso) return false;
  if (a.dcProjectIterations != b.dcProjectIterations) return false;

  if (a.buildLods != b.buildLods) return false;
  if (a.lodLevels != b.lodLevels) return false;
  if (a.lodRatioPerLevel != b.lodRatioPerLevel) return false;
  if (a.previewLod != b.previewLod) return false;
  if (a.exportLod != b.exportLod) return false;
  if (a.exportAllLods != b.exportAllLods) return false;

  if (a.previewRaymarch != b.previewRaymarch) return false;
  if (!raymarchSettingsEquals(a.raymarchSettings, b.raymarchSettings)) return false;

  if (a.useProceduralTexture != b.useProceduralTexture) return false;
  if (a.applyTextureToRaymarch != b.applyTextureToRaymarch) return false;
  if (a.texPreset != b.texPreset) return false;
  if (a.texLockToPreset != b.texLockToPreset) return false;
  if (a.texSeed != b.texSeed) return false;
  if (a.texResolution != b.texResolution) return false;
  if (a.texGenerateMips != b.texGenerateMips) return false;
  if (a.texDitherStrength != b.texDitherStrength) return false;
  if (a.texPackHeightInAlpha != b.texPackHeightInAlpha) return false;

  if (!sdfGraphEquals(a.graph, b.graph)) return false;
  if (!procGraphEquals(a.texGraph, b.texGraph)) return false;

  return true;
}

static HistorySnapshot makeHistorySnapshot(const ProceduralMeshLabWindowState& st) {
  HistorySnapshot s;
  s.graph = st.graph;
  s.preset = st.preset;
  s.lockToPreset = st.lockToPreset;
  s.seed = st.seed;

  s.mesherType = st.mesherType;
  s.mesher = st.mesher;
  s.dcQefRegularization = st.dcQefRegularization;
  s.dcClampToCell = st.dcClampToCell;
  s.dcProjectToIso = st.dcProjectToIso;
  s.dcProjectIterations = st.dcProjectIterations;

  s.buildLods = st.buildLods;
  s.lodLevels = st.lodLevels;
  s.lodRatioPerLevel = st.lodRatioPerLevel;
  s.previewLod = st.previewLod;
  s.exportLod = st.exportLod;
  s.exportAllLods = st.exportAllLods;

  s.previewRaymarch = st.previewRaymarch;
  s.raymarchSettings = st.raymarchSettings;

  s.useProceduralTexture = st.useProceduralTexture;
  s.applyTextureToRaymarch = st.applyTextureToRaymarch;
  s.texPreset = st.texPreset;
  s.texLockToPreset = st.texLockToPreset;
  s.texSeed = st.texSeed;
  s.texResolution = st.texResolution;
  s.texGenerateMips = st.texGenerateMips;
  s.texDitherStrength = st.texDitherStrength;
  s.texPackHeightInAlpha = st.texPackHeightInAlpha;
  s.texGraph = st.texGraph;

  return s;
}

static void applyHistorySnapshot(ProceduralMeshLabWindowState& st, const HistorySnapshot& s) {
  st.graph = s.graph;
  st.preset = s.preset;
  st.lockToPreset = s.lockToPreset;
  st.seed = s.seed;

  st.mesherType = s.mesherType;
  st.mesher = s.mesher;
  st.dcQefRegularization = s.dcQefRegularization;
  st.dcClampToCell = s.dcClampToCell;
  st.dcProjectToIso = s.dcProjectToIso;
  st.dcProjectIterations = s.dcProjectIterations;

  st.buildLods = s.buildLods;
  st.lodLevels = s.lodLevels;
  st.lodRatioPerLevel = s.lodRatioPerLevel;
  st.previewLod = s.previewLod;
  st.exportLod = s.exportLod;
  st.exportAllLods = s.exportAllLods;

  st.previewRaymarch = s.previewRaymarch;
  st.raymarchSettings = s.raymarchSettings;

  st.useProceduralTexture = s.useProceduralTexture;
  st.applyTextureToRaymarch = s.applyTextureToRaymarch;
  st.texPreset = s.texPreset;
  st.texLockToPreset = s.texLockToPreset;
  st.texSeed = s.texSeed;
  st.texResolution = s.texResolution;
  st.texGenerateMips = s.texGenerateMips;
  st.texDitherStrength = s.texDitherStrength;
  st.texPackHeightInAlpha = s.texPackHeightInAlpha;
  st.texGraph = s.texGraph;

  // Prevent "locked preset" auto-apply from clobbering a restored snapshot.
  st.texAppliedPreset = st.texPreset;
  st.texAppliedSeed = st.texSeed;

  // Force rebuilds after restoring.
  st.dirty = true;
  st.texDirty = true;
  st.hasMesh = false;
  st.hasLods = false;
  st.lastError.clear();
  st.lodError.clear();
}

static void historyTrim(std::vector<HistorySnapshot>& v, int maxEntries) {
  maxEntries = std::max(1, maxEntries);
  if ((int)v.size() > maxEntries) {
    const int drop = (int)v.size() - maxEntries;
    v.erase(v.begin(), v.begin() + drop);
  }
}

static void historyEnsureBaseline(ProceduralMeshLabWindowState& st) {
  if (!st.historyBaselineValid) {
    st.historyBaseline = makeHistorySnapshot(st);
    st.historyBaselineValid = true;
  }
}

static void historyClear(ProceduralMeshLabWindowState& st, float timeSec) {
  st.historyUndo.clear();
  st.historyRedo.clear();
  st.historyBaseline = makeHistorySnapshot(st);
  st.historyBaselineValid = true;
  st.historyPending = false;
  st.historyPendingSince = timeSec;
}

static void historyCommit(ProceduralMeshLabWindowState& st, const HistorySnapshot& cur, float timeSec) {
  historyEnsureBaseline(st);
  st.historyUndo.push_back(st.historyBaseline);
  historyTrim(st.historyUndo, st.historyMax);
  st.historyBaseline = cur;
  st.historyBaselineValid = true;
  st.historyPending = false;
  st.historyPendingSince = timeSec;
  st.historyRedo.clear();
}

static bool historyUndo(ProceduralMeshLabWindowState& st, float timeSec) {
  if (!st.historyEnabled) return false;
  historyEnsureBaseline(st);

  // If we have uncommitted edits, undo first cancels them back to baseline.
  if (st.historyPending) {
    const HistorySnapshot cur = makeHistorySnapshot(st);
    st.historyRedo.push_back(cur);
    historyTrim(st.historyRedo, st.historyMax);
    applyHistorySnapshot(st, st.historyBaseline);
    st.historyPending = false;
    st.historyPendingSince = timeSec;
    return true;
  }

  if (st.historyUndo.empty()) return false;
  const HistorySnapshot cur = st.historyBaseline;
  const HistorySnapshot prev = st.historyUndo.back();
  st.historyUndo.pop_back();

  st.historyRedo.push_back(cur);
  historyTrim(st.historyRedo, st.historyMax);

  applyHistorySnapshot(st, prev);
  st.historyBaseline = prev;
  st.historyBaselineValid = true;
  st.historyPending = false;
  st.historyPendingSince = timeSec;
  return true;
}

static bool historyRedo(ProceduralMeshLabWindowState& st, float timeSec) {
  if (!st.historyEnabled) return false;
  historyEnsureBaseline(st);
  if (st.historyRedo.empty()) return false;

  const HistorySnapshot cur = st.historyBaseline;
  const HistorySnapshot next = st.historyRedo.back();
  st.historyRedo.pop_back();

  st.historyUndo.push_back(cur);
  historyTrim(st.historyUndo, st.historyMax);

  applyHistorySnapshot(st, next);
  st.historyBaseline = next;
  st.historyBaselineValid = true;
  st.historyPending = false;
  st.historyPendingSince = timeSec;
  return true;
}

static void historyAutoCommit(ProceduralMeshLabWindowState& st, float timeSec) {
  historyEnsureBaseline(st);

  // If history is disabled, keep baseline synced and don't record.
  if (!st.historyEnabled) {
    st.historyBaseline = makeHistorySnapshot(st);
    st.historyBaselineValid = true;
    st.historyPending = false;
    st.historyPendingSince = timeSec;
    return;
  }

  const HistorySnapshot cur = makeHistorySnapshot(st);
  if (historySnapshotEquals(cur, st.historyBaseline)) {
    st.historyPending = false;
    return;
  }

  if (!st.historyPending) {
    st.historyPending = true;
    st.historyRedo.clear(); // new edits invalidate redo chain
  }

  // Debounce commits while dragging/typing.
  if (ImGui::IsAnyItemActive()) {
    st.historyPendingSince = timeSec;
  }

  const float delay = std::max(0.0f, st.historyCoalesceSec);
  if (!ImGui::IsAnyItemActive() && (timeSec - st.historyPendingSince) >= delay) {
    historyCommit(st, cur, timeSec);
  }
}

static const char* meshOpHelp(render::SdfNodeOp op) {
  switch (op) {
    case render::SdfNodeOp::Constant:
      return "Outputs a constant distance value (rarely useful alone; handy for debugging).";
    case render::SdfNodeOp::Sphere:
      return "Sphere SDF. Params: center (p0..p2), radius (p3).";
    case render::SdfNodeOp::Box:
      return "Axis-aligned box SDF. Params: center (p0..p2), halfSize (p3..p5).";
    case render::SdfNodeOp::Capsule:
      return "Capsule SDF. Params: a (p0..p2), b (p3..p5), radius (p6).";
    case render::SdfNodeOp::TorusY:
      return "Torus SDF around the Y axis. Params: center (p0..p2), majorR (p3), minorR (p4).";
    case render::SdfNodeOp::Plane:
      return "Plane SDF. Params: normal (p0..p2), offset (p3).";
    case render::SdfNodeOp::Union:
      return "CSG union. Output = min(A,B).";
    case render::SdfNodeOp::SmoothUnion:
      return "Smooth union (polynomial smin). Param: k (p0).";
    case render::SdfNodeOp::Intersect:
      return "CSG intersection. Output = max(A,B).";
    case render::SdfNodeOp::Subtract:
      return "CSG subtraction. Output = max(A,-B).";
    case render::SdfNodeOp::NoiseDisplace:
      return "Displaces SDF surface along normal by FBM noise. Params: freq (p0), amp (p1), octaves (p2), lacunarity (p3), gain (p4), offset (p5..p7).";
    case render::SdfNodeOp::NoiseDisplacePerlin:
      return "Displaces SDF surface along normal by FBM Perlin (gradient) noise. Same params as NoiseDisplace, but usually smoother / less grid-like.";
    case render::SdfNodeOp::Shell:
      return "Turns a solid into a thin shell. Output = abs(d) - thickness (p0).";

    case render::SdfNodeOp::Translate:
      return "Domain translate. Evaluates input at (p - t). Params: translate (p0..p2).";

    case render::SdfNodeOp::RotateX:
      return "Domain rotate about X. Evaluates input at rotX(p - pivot, -angle). Params: angleDeg (p0), pivot (p1..p3).";

    case render::SdfNodeOp::RotateY:
      return "Domain rotate about Y. Evaluates input at rotY(p - pivot, -angle). Params: angleDeg (p0), pivot (p1..p3).";

    case render::SdfNodeOp::RotateZ:
      return "Domain rotate about Z. Evaluates input at rotZ(p - pivot, -angle). Params: angleDeg (p0), pivot (p1..p3).";

    case render::SdfNodeOp::Scale:
      return "Uniform scale. Evaluates input at (p - pivot)/s + pivot and multiplies distance by s. Params: s (p0), pivot (p1..p3).";

    case render::SdfNodeOp::Repeat:
      return "Infinite grid repetition. Wraps p into [-period/2,period/2) per axis before evaluating input. Params: period (p0..p2), offset (p3..p5).";

    case render::SdfNodeOp::Mirror:
      return "Mirrors space by taking abs() across selected axes around a pivot. Params: mask (p0..p2 as 0/1), pivot (p3..p5).";

    case render::SdfNodeOp::TwistY:
      return "Twists space around Y: rotates XZ by angle=k*y. Params: k (deg/unit) (p0), pivot (p1..p3).";

    default:
      return "";
  }
}

static bool beginComboNodeIndex(const char* label,
                                int& idx,
                                const render::SdfGraph& g,
                                int maxIndex,
                                const char* noneLabel = "None") {
  // idx may be -1.
  std::string preview;
  if (idx < 0) {
    preview = noneLabel;
  } else if (idx >= 0 && idx <= maxIndex && idx < (int)g.nodes.size()) {
    preview = "#" + std::to_string(idx) + " " + render::sdfNodeOpName(g.nodes[(std::size_t)idx].op);
  } else {
    preview = "(invalid)";
  }

  bool changed = false;
  if (ImGui::BeginCombo(label, preview.c_str())) {
    // None
    if (ImGui::Selectable(noneLabel, idx < 0)) {
      idx = -1;
      changed = true;
    }

    // Nodes
    const int n = std::min(maxIndex + 1, (int)g.nodes.size());
    for (int i = 0; i < n; ++i) {
      std::string name = "#" + std::to_string(i) + " " + render::sdfNodeOpName(g.nodes[(std::size_t)i].op);
      const bool isSel = (idx == i);
      if (ImGui::Selectable(name.c_str(), isSel)) {
        idx = i;
        changed = true;
      }
      if (isSel) ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }
  return changed;
}

static bool beginComboOutputIndex(const char* label, int& idx, const render::SdfGraph& g) {
  const int maxIndex = (int)g.nodes.size() - 1;
  if (maxIndex < 0) {
    ImGui::TextUnformatted("(no nodes)");
    idx = -1;
    return false;
  }

  std::string preview;
  if (idx >= 0 && idx <= maxIndex) {
    preview = "#" + std::to_string(idx) + " " + render::sdfNodeOpName(g.nodes[(std::size_t)idx].op);
  } else {
    preview = "#" + std::to_string(maxIndex);
  }

  bool changed = false;
  if (ImGui::BeginCombo(label, preview.c_str())) {
    for (int i = 0; i <= maxIndex; ++i) {
      std::string name = "#" + std::to_string(i) + " " + render::sdfNodeOpName(g.nodes[(std::size_t)i].op);
      const bool isSel = (idx == i);
      if (ImGui::Selectable(name.c_str(), isSel)) {
        idx = i;
        changed = true;
      }
      if (isSel) ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }
  return changed;
}

static bool writeObj(const std::string& path, const render::SdfMeshData& mesh, std::string* outErr) {
  namespace fs = std::filesystem;

  if (path.empty()) {
    if (outErr) *outErr = "Empty output path.";
    return false;
  }

  if (mesh.vertices.empty() || mesh.indices.empty()) {
    if (outErr) *outErr = "No mesh data to export.";
    return false;
  }

  std::error_code ec;
  fs::path p(path);
  if (p.has_parent_path()) {
    fs::create_directories(p.parent_path(), ec);
    // If directory creation fails, we'll still try to open; the open will fail and report.
  }

  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f.is_open()) {
    if (outErr) *outErr = "Failed to open file for writing.";
    return false;
  }

  f << "# Stellar Procedural Mesh Lab export\n";
  f << "o procedural_mesh\n";

  // Positions
  for (const auto& v : mesh.vertices) {
    f << "v " << v.px << ' ' << v.py << ' ' << v.pz << "\n";
  }
  // UVs
  for (const auto& v : mesh.vertices) {
    f << "vt " << v.u << ' ' << v.v << "\n";
  }
  // Normals
  for (const auto& v : mesh.vertices) {
    f << "vn " << v.nx << ' ' << v.ny << ' ' << v.nz << "\n";
  }

  // Faces (triangles). We use 1-based indexing and keep v/vt/vn aligned.
  for (std::size_t i = 0; i + 2 < mesh.indices.size(); i += 3) {
    const std::uint32_t a = mesh.indices[i + 0] + 1u;
    const std::uint32_t b = mesh.indices[i + 1] + 1u;
    const std::uint32_t c = mesh.indices[i + 2] + 1u;
    f << "f "
      << a << '/' << a << '/' << a << ' '
      << b << '/' << b << '/' << b << ' '
      << c << '/' << c << '/' << c << "\n";
  }

  return true;
}

static bool writeObjWithMtl(const std::string& path,
                            const render::SdfMeshData& mesh,
                            const std::string& mtlFile,
                            const std::string& materialName,
                            std::string* outErr) {
  namespace fs = std::filesystem;

  if (path.empty()) {
    if (outErr) *outErr = "Empty output path.";
    return false;
  }

  if (mesh.vertices.empty() || mesh.indices.empty()) {
    if (outErr) *outErr = "No mesh data to export.";
    return false;
  }

  std::error_code ec;
  fs::path p(path);
  if (p.has_parent_path()) {
    fs::create_directories(p.parent_path(), ec);
  }

  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f.is_open()) {
    if (outErr) *outErr = "Failed to open file for writing.";
    return false;
  }

  f << "# Stellar Procedural Mesh Lab export\n";
  if (!mtlFile.empty()) f << "mtllib " << mtlFile << "\n";
  f << "o procedural_mesh\n";
  if (!materialName.empty()) f << "usemtl " << materialName << "\n";

  // Positions
  for (const auto& v : mesh.vertices) {
    f << "v " << v.px << ' ' << v.py << ' ' << v.pz << "\n";
  }
  // UVs
  for (const auto& v : mesh.vertices) {
    f << "vt " << v.u << ' ' << v.v << "\n";
  }
  // Normals
  for (const auto& v : mesh.vertices) {
    f << "vn " << v.nx << ' ' << v.ny << ' ' << v.nz << "\n";
  }

  // Faces (triangles). We use 1-based indexing and keep v/vt/vn aligned.
  for (std::size_t i = 0; i + 2 < mesh.indices.size(); i += 3) {
    const std::uint32_t a = mesh.indices[i + 0] + 1u;
    const std::uint32_t b = mesh.indices[i + 1] + 1u;
    const std::uint32_t c = mesh.indices[i + 2] + 1u;
    f << "f "
      << a << '/' << a << '/' << a << ' '
      << b << '/' << b << '/' << b << ' '
      << c << '/' << c << '/' << c << "\n";
  }

  return true;
}

static bool writeMtl(const std::string& path,
                     const std::string& materialName,
                     const std::string& albedoFile,
                     std::string* outErr) {
  namespace fs = std::filesystem;
  std::error_code ec;
  fs::path p(path);
  if (p.has_parent_path()) fs::create_directories(p.parent_path(), ec);

  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f.is_open()) {
    if (outErr) *outErr = "Failed to open MTL file for writing.";
    return false;
  }

  f << "# Stellar Procedural Mesh Lab material\n";
  f << "newmtl " << (materialName.empty() ? "mat0" : materialName) << "\n";
  f << "Ka 0.000 0.000 0.000\n";
  f << "Kd 1.000 1.000 1.000\n";
  f << "Ks 0.100 0.100 0.100\n";
  f << "Ns 32.000\n";
  f << "illum 2\n";
  if (!albedoFile.empty()) {
    f << "map_Kd " << albedoFile << "\n";
  }
  return true;
}

static core::u64 randomU64() {
  const auto t = (core::u64)std::chrono::high_resolution_clock::now().time_since_epoch().count();
  return core::hashCombine(0xC0FFEEULL, t);
}

// -----------------------------------------------------------------------------
// SDF Graph diagnostics + "Variation Studio" helpers
// -----------------------------------------------------------------------------

struct SdfSampleStats {
  float insideFraction{0.0f}; // [0,1] fraction of sampled points inside iso
  float boundingRadius{0.0f}; // max radius of inside points (approx)
  float minDistance{0.0f};
  float maxDistance{0.0f};
  int samples{0};
};

static core::u64 sdfGraphFingerprint(const render::SdfGraph& g) {
  // This is used for variant dedupe; it is not meant to be cryptographic.
  core::u64 h = 0x51F15EEDBADC0FFEuLL;
  h = core::hashCombine(h, g.seed);
  h = core::hashCombine(h, (core::u64)g.output);
  h = core::hashCombine(h, (core::u64)g.nodes.size());

  for (const auto& n : g.nodes) {
    h = core::hashCombine(h, (core::u64)n.op);
    h = core::hashCombine(h, (core::u64)(n.in0 + 2));
    h = core::hashCombine(h, (core::u64)(n.in1 + 2));
    h = core::hashCombine(h, n.seed);

    const float params[8] = {n.p0, n.p1, n.p2, n.p3, n.p4, n.p5, n.p6, n.p7};
    h = core::hashCombine(h, core::hashBytes(params, sizeof(params)));
  }

  return h;
}

static SdfSampleStats sampleSdfGraphStats(const render::SdfGraph& g, float bounds, float iso, int gridRes) {
  SdfSampleStats s{};

  if (bounds <= 0.0f) bounds = 1.0f;
  gridRes = std::clamp(gridRes, 4, 64);

  const int N = gridRes;
  const float step = (2.0f * bounds) / (float)(N - 1);

  int inside = 0;
  float minD = 1e9f;
  float maxD = -1e9f;
  float maxR = 0.0f;

  for (int z = 0; z < N; ++z) {
    const float pz = -bounds + step * (float)z;
    for (int y = 0; y < N; ++y) {
      const float py = -bounds + step * (float)y;
      for (int x = 0; x < N; ++x) {
        const float px = -bounds + step * (float)x;

        const float d = render::evalSdfGraph(g, px, py, pz);
        minD = std::min(minD, d);
        maxD = std::max(maxD, d);

        if (d <= iso) {
          inside++;
          const float r = std::sqrt(px * px + py * py + pz * pz);
          maxR = std::max(maxR, r);
        }
      }
    }
  }

  const int total = N * N * N;
  s.samples = total;
  s.insideFraction = (total > 0) ? ((float)inside / (float)total) : 0.0f;
  s.boundingRadius = maxR;
  s.minDistance = minD;
  s.maxDistance = maxD;
  return s;
}

static void markReachableRec(const render::SdfGraph& g, int idx, std::vector<unsigned char>& mark) {
  const int n = (int)g.nodes.size();
  if (idx < 0 || idx >= n) return;
  if (mark[(std::size_t)idx]) return;
  mark[(std::size_t)idx] = 1;

  const auto& node = g.nodes[(std::size_t)idx];
  markReachableRec(g, node.in0, mark);
  markReachableRec(g, node.in1, mark);
}

static bool sdfGraphHasCycle(const render::SdfGraph& g) {
  const int n = (int)g.nodes.size();
  if (n <= 0) return false;

  // 0=unvisited, 1=visiting, 2=done
  std::vector<unsigned char> state((std::size_t)n, 0);

  std::function<bool(int)> dfs = [&](int i) -> bool {
    if (i < 0 || i >= n) return false;
    const unsigned char st = state[(std::size_t)i];
    if (st == 1) return true;
    if (st == 2) return false;

    state[(std::size_t)i] = 1;
    const auto& node = g.nodes[(std::size_t)i];
    if (dfs(node.in0)) return true;
    if (dfs(node.in1)) return true;
    state[(std::size_t)i] = 2;
    return false;
  };

  for (int i = 0; i < n; ++i) {
    if (state[(std::size_t)i] == 0) {
      if (dfs(i)) return true;
    }
  }
  return false;
}

static void fixSdfGraphIndices(render::SdfGraph& g) {
  const int n = (int)g.nodes.size();
  if (n <= 0) {
    g.output = -1;
    return;
  }

  // Keep graphs DAG-ish (and matching editor behavior): inputs must refer to earlier nodes.
  for (int i = 0; i < n; ++i) {
    auto& node = g.nodes[(std::size_t)i];
    if (node.in0 < -1 || node.in0 >= i) node.in0 = -1;
    if (node.in1 < -1 || node.in1 >= i) node.in1 = -1;
  }

  if (g.output < 0 || g.output >= n) g.output = n - 1;
}

static render::SdfGraph pruneUnreachableNodes(const render::SdfGraph& g) {
  render::SdfGraph out = g;
  const int n = (int)g.nodes.size();
  if (n <= 0) return out;

  std::vector<unsigned char> reachable((std::size_t)n, 0);
  markReachableRec(g, g.output, reachable);

  // Build mapping old -> new
  std::vector<int> map((std::size_t)n, -1);
  std::vector<render::SdfNode> newNodes;
  newNodes.reserve((std::size_t)n);

  for (int i = 0; i < n; ++i) {
    if (!reachable[(std::size_t)i]) continue;
    map[(std::size_t)i] = (int)newNodes.size();
    newNodes.push_back(g.nodes[(std::size_t)i]);
  }

  out.nodes = std::move(newNodes);

  // Rewrite inputs
  const int m = (int)out.nodes.size();
  for (int i = 0; i < m; ++i) {
    auto& node = out.nodes[(std::size_t)i];
    const int in0 = node.in0;
    const int in1 = node.in1;
    node.in0 = (in0 >= 0 && in0 < n) ? map[(std::size_t)in0] : -1;
    node.in1 = (in1 >= 0 && in1 < n) ? map[(std::size_t)in1] : -1;

    // Keep DAG-ish order if possible; if pruning collapsed indices in a way that
    // creates forward refs, drop them to -1.
    if (node.in0 >= i) node.in0 = -1;
    if (node.in1 >= i) node.in1 = -1;
  }

  out.output = (g.output >= 0 && g.output < n) ? map[(std::size_t)g.output] : -1;
  if (out.output < 0 && !out.nodes.empty()) out.output = (int)out.nodes.size() - 1;

  fixSdfGraphIndices(out);
  return out;
}

enum class SdfOpGroup : core::u8 {
  Constant = 0,
  Primitive,
  Boolean,
  Modifier,
  Transform,
};

static SdfOpGroup sdfOpGroup(render::SdfNodeOp op) {
  using render::SdfNodeOp;
  switch (op) {
    case SdfNodeOp::Constant: return SdfOpGroup::Constant;

    case SdfNodeOp::Sphere:
    case SdfNodeOp::Box:
    case SdfNodeOp::Capsule:
    case SdfNodeOp::TorusY:
    case SdfNodeOp::Plane:
      return SdfOpGroup::Primitive;

    case SdfNodeOp::Union:
    case SdfNodeOp::SmoothUnion:
    case SdfNodeOp::Intersect:
    case SdfNodeOp::Subtract:
      return SdfOpGroup::Boolean;

    case SdfNodeOp::NoiseDisplace:
    case SdfNodeOp::NoiseDisplacePerlin:
    case SdfNodeOp::Shell:
      return SdfOpGroup::Modifier;

    default:
      return SdfOpGroup::Transform;
  }
}

static int sdfOpInputCount(render::SdfNodeOp op) {
  switch (sdfOpGroup(op)) {
    case SdfOpGroup::Boolean: return 2;
    case SdfOpGroup::Modifier: return 1;
    case SdfOpGroup::Transform: return 1;
    default: return 0;
  }
}

static render::SdfNodeOp randomOpInGroup(SdfOpGroup group, core::SplitMix64& rng, render::SdfNodeOp avoid) {
  using render::SdfNodeOp;
  const SdfNodeOp* ops = nullptr;
  int count = 0;

  static const SdfNodeOp kPrims[] = {SdfNodeOp::Sphere, SdfNodeOp::Box, SdfNodeOp::Capsule, SdfNodeOp::TorusY, SdfNodeOp::Plane};
  static const SdfNodeOp kBools[] = {SdfNodeOp::Union, SdfNodeOp::SmoothUnion, SdfNodeOp::Intersect, SdfNodeOp::Subtract};
  static const SdfNodeOp kMods[]  = {SdfNodeOp::NoiseDisplace, SdfNodeOp::NoiseDisplacePerlin, SdfNodeOp::Shell};
  static const SdfNodeOp kXforms[] = {SdfNodeOp::Translate, SdfNodeOp::RotateX, SdfNodeOp::RotateY, SdfNodeOp::RotateZ,
                                      SdfNodeOp::Scale, SdfNodeOp::Repeat, SdfNodeOp::Mirror, SdfNodeOp::TwistY};

  switch (group) {
    case SdfOpGroup::Primitive: ops = kPrims; count = (int)(sizeof(kPrims) / sizeof(kPrims[0])); break;
    case SdfOpGroup::Boolean: ops = kBools; count = (int)(sizeof(kBools) / sizeof(kBools[0])); break;
    case SdfOpGroup::Modifier: ops = kMods; count = (int)(sizeof(kMods) / sizeof(kMods[0])); break;
    case SdfOpGroup::Transform: ops = kXforms; count = (int)(sizeof(kXforms) / sizeof(kXforms[0])); break;
    default: ops = nullptr; count = 0; break;
  }

  if (!ops || count <= 0) return SdfNodeOp::Constant;
  if (count == 1) return ops[0];

  // Try to avoid choosing the same op to make mutation visible.
  for (int tries = 0; tries < 8; ++tries) {
    const SdfNodeOp pick = ops[rng.range<int>(0, count - 1)];
    if (pick != avoid) return pick;
  }
  return ops[rng.range<int>(0, count - 1)];
}

static void randomizeNodeParams(render::SdfNode& n, core::SplitMix64& rng, float bounds) {
  using render::SdfNodeOp;
  if (bounds <= 0.0f) bounds = 1.0f;

  // Helper: centered range.
  auto cr = [&](float s) -> float { return rng.range(-s, s); };

  // Start with zeros so ops that don't use some params are stable.
  n.p0 = n.p1 = n.p2 = n.p3 = n.p4 = n.p5 = n.p6 = n.p7 = 0.0f;

  switch (n.op) {
    case SdfNodeOp::Constant: {
      n.p0 = rng.range(-1.0f, 1.0f);
      break;
    }

    case SdfNodeOp::Sphere: {
      n.p0 = cr(0.35f * bounds);
      n.p1 = cr(0.35f * bounds);
      n.p2 = cr(0.35f * bounds);
      n.p3 = rng.range(0.08f * bounds, 0.85f * bounds); // radius
      break;
    }

    case SdfNodeOp::Box: {
      n.p0 = cr(0.35f * bounds);
      n.p1 = cr(0.35f * bounds);
      n.p2 = cr(0.35f * bounds);
      n.p3 = rng.range(0.05f * bounds, 0.75f * bounds);
      n.p4 = rng.range(0.05f * bounds, 0.75f * bounds);
      n.p5 = rng.range(0.05f * bounds, 0.75f * bounds);
      break;
    }

    case SdfNodeOp::Capsule: {
      n.p0 = cr(0.50f * bounds);
      n.p1 = cr(0.50f * bounds);
      n.p2 = cr(0.50f * bounds);
      n.p3 = cr(0.50f * bounds);
      n.p4 = cr(0.50f * bounds);
      n.p5 = cr(0.50f * bounds);
      n.p6 = rng.range(0.03f * bounds, 0.35f * bounds);
      break;
    }

    case SdfNodeOp::TorusY: {
      n.p0 = cr(0.25f * bounds);
      n.p1 = cr(0.25f * bounds);
      n.p2 = cr(0.25f * bounds);
      n.p3 = rng.range(0.10f * bounds, 0.85f * bounds); // major
      n.p4 = rng.range(0.03f * bounds, 0.35f * bounds); // minor
      break;
    }

    case SdfNodeOp::Plane: {
      // Random-ish normal (not uniformly distributed, but good enough for mutations).
      n.p0 = rng.range(-1.0f, 1.0f);
      n.p1 = rng.range(-1.0f, 1.0f);
      n.p2 = rng.range(-1.0f, 1.0f);
      n.p3 = cr(0.35f * bounds); // offset
      break;
    }

    case SdfNodeOp::SmoothUnion: {
      n.p0 = rng.range(0.01f, 0.35f * bounds); // k
      break;
    }

    case SdfNodeOp::NoiseDisplace:
    case SdfNodeOp::NoiseDisplacePerlin: {
      n.p0 = rng.range(0.5f, 10.0f); // freq
      n.p1 = rng.range(-0.35f, 0.35f); // amp
      n.p2 = (float)rng.range<int>(2, 7); // octaves (stored as float)
      n.p3 = rng.range(1.8f, 3.2f); // lacunarity
      n.p4 = rng.range(0.35f, 0.75f); // gain
      n.p5 = cr(50.0f);
      n.p6 = cr(50.0f);
      n.p7 = cr(50.0f);
      break;
    }

    case SdfNodeOp::Shell: {
      n.p0 = rng.range(0.01f, 0.25f * bounds); // thickness
      break;
    }

    case SdfNodeOp::Translate: {
      n.p0 = cr(0.6f * bounds);
      n.p1 = cr(0.6f * bounds);
      n.p2 = cr(0.6f * bounds);
      break;
    }

    case SdfNodeOp::RotateX:
    case SdfNodeOp::RotateY:
    case SdfNodeOp::RotateZ: {
      n.p0 = rng.range(-180.0f, 180.0f); // angle deg
      n.p1 = cr(0.35f * bounds);
      n.p2 = cr(0.35f * bounds);
      n.p3 = cr(0.35f * bounds);
      break;
    }

    case SdfNodeOp::Scale: {
      n.p0 = rng.range(0.35f, 2.75f);
      n.p1 = cr(0.25f * bounds);
      n.p2 = cr(0.25f * bounds);
      n.p3 = cr(0.25f * bounds);
      break;
    }

    case SdfNodeOp::Repeat: {
      n.p0 = rng.range(0.30f * bounds, 1.75f * bounds); // period x
      n.p1 = rng.range(0.30f * bounds, 1.75f * bounds); // period y
      n.p2 = rng.range(0.30f * bounds, 1.75f * bounds); // period z
      n.p3 = cr(0.30f * bounds);
      n.p4 = cr(0.30f * bounds);
      n.p5 = cr(0.30f * bounds);
      break;
    }

    case SdfNodeOp::Mirror: {
      n.p0 = (float)rng.range<int>(0, 1);
      n.p1 = (float)rng.range<int>(0, 1);
      n.p2 = (float)rng.range<int>(0, 1);
      n.p3 = cr(0.25f * bounds);
      n.p4 = cr(0.25f * bounds);
      n.p5 = cr(0.25f * bounds);
      break;
    }

    case SdfNodeOp::TwistY: {
      n.p0 = rng.range(-90.0f, 90.0f); // deg per unit
      n.p1 = cr(0.25f * bounds);
      n.p2 = cr(0.25f * bounds);
      n.p3 = cr(0.25f * bounds);
      break;
    }

    default:
      break;
  }
}

static void sanitizeNodeParams(render::SdfNode& n) {
  using render::SdfNodeOp;

  auto clampf = [](float v, float lo, float hi) -> float {
    return std::max(lo, std::min(v, hi));
  };

  switch (n.op) {
    case SdfNodeOp::Sphere:
      n.p3 = std::max(0.001f, n.p3);
      break;

    case SdfNodeOp::Box:
      n.p3 = std::max(0.001f, n.p3);
      n.p4 = std::max(0.001f, n.p4);
      n.p5 = std::max(0.001f, n.p5);
      break;

    case SdfNodeOp::Capsule:
      n.p6 = std::max(0.001f, n.p6);
      break;

    case SdfNodeOp::TorusY:
      n.p3 = std::max(0.001f, n.p3);
      n.p4 = std::max(0.001f, n.p4);
      break;

    case SdfNodeOp::Plane: {
      const float len = std::sqrt(n.p0 * n.p0 + n.p1 * n.p1 + n.p2 * n.p2);
      if (len > 1e-6f) {
        n.p0 /= len;
        n.p1 /= len;
        n.p2 /= len;
      } else {
        n.p0 = 0.0f;
        n.p1 = 1.0f;
        n.p2 = 0.0f;
      }
      break;
    }

    case SdfNodeOp::SmoothUnion:
      n.p0 = std::max(0.0001f, n.p0);
      break;

    case SdfNodeOp::NoiseDisplace:
    case SdfNodeOp::NoiseDisplacePerlin: {
      n.p0 = std::max(0.0f, n.p0);
      const int oct = std::max(1, std::min(12, (int)std::lround(n.p2)));
      n.p2 = (float)oct;
      if (n.p3 <= 0.0f) n.p3 = 2.0f;
      if (n.p4 <= 0.0f) n.p4 = 0.5f;
      break;
    }

    case SdfNodeOp::Shell:
      n.p0 = std::max(0.0f, n.p0);
      break;

    case SdfNodeOp::Scale:
      n.p0 = clampf(n.p0, 0.05f, 50.0f);
      break;

    case SdfNodeOp::Repeat:
      n.p0 = std::max(0.001f, n.p0);
      n.p1 = std::max(0.001f, n.p1);
      n.p2 = std::max(0.001f, n.p2);
      break;

    case SdfNodeOp::Mirror:
      n.p0 = (float)((int)std::lround(n.p0) ? 1 : 0);
      n.p1 = (float)((int)std::lround(n.p1) ? 1 : 0);
      n.p2 = (float)((int)std::lround(n.p2) ? 1 : 0);
      break;

    default:
      break;
  }
}

static render::SdfGraph makeSdfVariant(const render::SdfGraph& base,
                                      core::u64 variantSeed,
                                      float bounds,
                                      float paramJitter,
                                      float rewireChance,
                                      float opChance,
                                      float wrapChance,
                                      bool mutateNodeSeeds) {
  render::SdfGraph g = base;
  core::SplitMix64 rng(variantSeed);

  // Ensure the variant "feels" different even if only noise nodes exist.
  g.seed = core::hashCombine(base.seed, variantSeed);

  const int n = (int)g.nodes.size();
  if (n <= 0) {
    g.output = -1;
    return g;
  }

  for (int i = 0; i < n; ++i) {
    auto& node = g.nodes[(std::size_t)i];

    const SdfOpGroup grp = sdfOpGroup(node.op);

    bool opMutated = false;
    if (opChance > 0.0f && rng.chance(opChance)) {
      const auto newOp = randomOpInGroup(grp, rng, node.op);
      if (newOp != node.op) {
        node.op = newOp;
        opMutated = true;
        randomizeNodeParams(node, rng, bounds);
      }
    }

    // Mutate node-local noise seed tweaks (for NoiseDisplace*, etc).
    if (mutateNodeSeeds && rng.chance(0.25)) {
      node.seed = rng.nextU64();
    }

    // Rewire inputs (stay within earlier nodes to keep the graph DAG-ish).
    const int inputs = sdfOpInputCount(node.op);
    if (rewireChance > 0.0f && rng.chance(rewireChance) && i > 0) {
      if (inputs >= 1) node.in0 = rng.range(-1, i - 1);
      if (inputs >= 2) node.in1 = rng.range(-1, i - 1);
    }

    // If op changed, enforce the right number of inputs.
    if (opMutated) {
      if (inputs == 0) {
        node.in0 = -1;
        node.in1 = -1;
      } else if (inputs == 1) {
        if (node.in0 < -1 || node.in0 >= i) node.in0 = (i > 0) ? rng.range(-1, i - 1) : -1;
        node.in1 = -1;
      } else if (inputs == 2) {
        if (node.in0 < -1 || node.in0 >= i) node.in0 = (i > 0) ? rng.range(-1, i - 1) : -1;
        if (node.in1 < -1 || node.in1 >= i) node.in1 = (i > 0) ? rng.range(-1, i - 1) : -1;
      }
    }

    // Param jitter (relative-ish).
    if (paramJitter > 0.0f && rng.chance(0.92)) {
      auto jitter = [&](float& v) {
        const float scale = std::max(0.25f, std::fabs(v));
        v += rng.range(-1.0f, 1.0f) * paramJitter * scale;
      };

      jitter(node.p0);
      jitter(node.p1);
      jitter(node.p2);
      jitter(node.p3);
      jitter(node.p4);
      jitter(node.p5);
      jitter(node.p6);
      jitter(node.p7);
    }

    sanitizeNodeParams(node);
  }

  // Optionally wrap the output in a random domain op / modifier (adds complexity quickly).
  if (wrapChance > 0.0f && rng.chance(wrapChance) && (int)g.nodes.size() < render::kSdfGraphMaxNodes) {
    render::SdfNode wrap{};
    wrap.op = randomOpInGroup(SdfOpGroup::Transform, rng, render::SdfNodeOp::Translate);
    // Occasionally pick a modifier instead.
    if (rng.chance(0.35)) wrap.op = randomOpInGroup(SdfOpGroup::Modifier, rng, render::SdfNodeOp::Shell);

    wrap.in0 = g.output;
    wrap.in1 = -1;
    wrap.seed = rng.nextU64();
    randomizeNodeParams(wrap, rng, bounds);
    sanitizeNodeParams(wrap);

    g.nodes.push_back(wrap);
    g.output = (int)g.nodes.size() - 1;
  }

  fixSdfGraphIndices(g);
  return g;
}

static bool readTextureRgba8(const render::Texture2D& tex,
                             int width,
                             int height,
                             std::vector<unsigned char>& outPixels,
                             std::string* outError) {
  outPixels.assign((std::size_t)width * (std::size_t)height * 4u, 0u);
  glBindTexture(GL_TEXTURE_2D, (GLuint)tex.handle());
  glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, outPixels.data());
  const GLenum err = glGetError();
  if (err != GL_NO_ERROR) {
    if (outError) {
      std::ostringstream oss;
      oss << "glGetTexImage failed (GL error 0x" << std::hex << (unsigned int)err << ")";
      *outError = oss.str();
    }
    return false;
  }
  return true;
}

static bool ensurePreview(ProceduralMeshLabWindowState& st) {
  if (st.previewInited) return st.previewTarget.width() == st.previewResolution;

  st.previewInitError.clear();

  if (!st.previewTarget.init(st.previewResolution, st.previewResolution, &st.previewInitError)) {
    st.previewInited = false;
    return false;
  }

  // Basic checker texture for shading.
  st.checker.createChecker(128, 128, 16);

  std::string err;
  if (!st.previewRenderer.init(&err)) {
    st.previewInitError = err;
    st.previewInited = false;
    return false;
  }

  st.previewRenderer.setMesh(&st.previewMesh);
  st.previewRenderer.setTexture(&st.checker);
  st.previewRenderer.setNormalTexture(nullptr);
  st.previewRenderer.setSpecular(0.10f, 72.0f);
  st.previewRenderer.setUnlit(false);

  st.previewInited = true;
  return true;
}

static void ensurePreviewSize(ProceduralMeshLabWindowState& st) {
  if (!st.previewInited) return;
  st.previewTarget.ensureSize(st.previewResolution, st.previewResolution);
}

static void matToFloat(const math::Mat4d& m, float out16[16]) {
  for (int i = 0; i < 16; ++i) out16[i] = (float)m.m[i];
}

static void renderPreview(ProceduralMeshLabWindowState& st, float timeSec) {
  if (!st.previewInited || !st.hasMesh || st.previewMesh.indexCount() == 0) return;

  ensurePreviewSize(st);

  // Bind the current surface texture (procedural or checker) for this frame.
  if (st.useProceduralTexture && st.texBaker.isReady()) {
    st.previewRenderer.setTexture(&st.texBaker.texture());
  } else {
    st.previewRenderer.setTexture(&st.checker);
  }

  // Save GL state we touch.
  GLint prevFbo = 0;
  GLint prevViewport[4] = {0, 0, 0, 0};
  GLint prevProg = 0;
  GLint prevVao = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
  glGetIntegerv(GL_VIEWPORT, prevViewport);
  glGetIntegerv(GL_CURRENT_PROGRAM, &prevProg);
  glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &prevVao);

  const bool prevDepth = glIsEnabled(GL_DEPTH_TEST);
  const bool prevCull = glIsEnabled(GL_CULL_FACE);
  const bool prevBlend = glIsEnabled(GL_BLEND);

  GLint polyModeFront[2] = {GL_FILL, GL_FILL};
  glGetIntegerv(GL_POLYGON_MODE, polyModeFront);

  st.previewTarget.begin();

  glViewport(0, 0, st.previewTarget.width(), st.previewTarget.height());
  glDisable(GL_BLEND);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  glClearColor(0.02f, 0.02f, 0.03f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (st.wireframe) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  } else {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }

  // Orbit camera around origin.
  float yaw = math::degToRad((double)st.yawDeg);
  if (st.autoSpin) yaw += math::degToRad((double)(st.spinDegPerSec * timeSec));
  const float pitch = (float)math::degToRad((double)st.pitchDeg);

  const float cy = std::cos(yaw);
  const float sy = std::sin(yaw);
  const float cp = std::cos(pitch);
  const float sp = std::sin(pitch);

  const math::Vec3d eye{(double)(st.distance * cp * cy),
                       (double)(st.distance * sp),
                       (double)(st.distance * cp * sy)};

  const math::Mat4d view = math::Mat4d::lookAt(eye, math::Vec3d{0,0,0}, math::Vec3d{0,1,0});
  const double aspect = (double)st.previewTarget.width() / (double)st.previewTarget.height();
  const math::Mat4d proj = math::Mat4d::perspective(math::degToRad(55.0), aspect, 0.05, 200.0);

  float viewF[16];
  float projF[16];
  matToFloat(view, viewF);
  matToFloat(proj, projF);

  st.previewRenderer.setViewProj(viewF, projF);
  st.previewRenderer.setCameraPos((float)eye.x, (float)eye.y, (float)eye.z);
  // Light slightly above and to the side.
  st.previewRenderer.setLightPos(2.5f, 2.8f, 1.8f);

  render::InstanceData inst{};
  inst.px = 0.0f; inst.py = 0.0f; inst.pz = 0.0f;
  inst.sx = 1.0f; inst.sy = 1.0f; inst.sz = 1.0f;
  inst.qx = 0.0f; inst.qy = 0.0f; inst.qz = 0.0f; inst.qw = 1.0f;
  inst.cr = 0.85f; inst.cg = 0.90f; inst.cb = 1.0f;

  std::vector<render::InstanceData> insts;
  insts.push_back(inst);
  st.previewRenderer.drawInstances(insts);

  render::RenderTarget2D::end();

  // Restore state.
  if (!prevDepth) glDisable(GL_DEPTH_TEST);
  if (!prevCull) glDisable(GL_CULL_FACE);
  if (prevBlend) glEnable(GL_BLEND);

  render::gl::UseProgram((GLuint)prevProg);
  render::gl::BindVertexArray((GLuint)prevVao);
  render::gl::BindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
  glViewport(prevViewport[0], prevViewport[1], prevViewport[2], prevViewport[3]);

  // Restore polygon mode (OpenGL returns front/back packed on some drivers).
  glPolygonMode(GL_FRONT, polyModeFront[0]);
  glPolygonMode(GL_BACK, polyModeFront[1]);
}

static bool renderRaymarchPreview(ProceduralMeshLabWindowState& st, float timeSec) {
  st.raymarchError.clear();

  // Orbit camera around origin (same controls as raster preview).
  float yaw = (float)math::degToRad((double)st.yawDeg);
  if (st.autoSpin) yaw += (float)math::degToRad((double)(st.spinDegPerSec * timeSec));
  const float pitch = (float)math::degToRad((double)st.pitchDeg);

  const float cy = std::cos(yaw);
  const float sy = std::sin(yaw);
  const float cp = std::cos(pitch);
  const float sp = std::sin(pitch);

  const math::Vec3d eye{(double)(st.distance * cp * cy),
                       (double)(st.distance * sp),
                       (double)(st.distance * cp * sy)};

  const math::Vec3d fwd = (-eye).normalized();
  const math::Vec3d worldUp{0, 1, 0};
  math::Vec3d right = math::cross(fwd, worldUp).normalized();
  if (right.lengthSq() < 1e-10) right = math::Vec3d{1, 0, 0};
  const math::Vec3d up = math::cross(right, fwd).normalized();

  render::SdfRaymarchCamera cam{};
  cam.pos[0] = (float)eye.x; cam.pos[1] = (float)eye.y; cam.pos[2] = (float)eye.z;
  cam.forward[0] = (float)fwd.x; cam.forward[1] = (float)fwd.y; cam.forward[2] = (float)fwd.z;
  cam.right[0] = (float)right.x; cam.right[1] = (float)right.y; cam.right[2] = (float)right.z;
  cam.up[0] = (float)up.x; cam.up[1] = (float)up.y; cam.up[2] = (float)up.z;
  cam.fovYRadians = (float)math::degToRad((double)st.raymarchFovDeg);

  // Optional albedo texture from the 2D ProcGraph.
  const bool useTex = st.applyTextureToRaymarch && st.useProceduralTexture && st.texBaker.isReady();
  render::SdfRaymarchMaterial mat = st.raymarchMat;
  mat.albedoTex = useTex ? &st.texBaker.texture() : nullptr;
  const render::SdfRaymarchMaterial* matPtr = useTex ? &mat : nullptr;

  return st.raymarcher.render(st.graph,
                              st.previewResolution,
                              st.previewResolution,
                              timeSec,
                              st.mesher.bounds,
                              st.mesher.iso,
                              cam,
                              st.raymarchSettings,
                              &st.raymarchError,
                              matPtr);
}


static render::SdfMeshData buildMeshForSettings(const render::ScalarField3D& field,
                                                const render::SdfMesherParams& base,
                                                ProceduralMeshLabWindowState::MesherType type,
                                                float dcQefRegularization,
                                                bool dcClampToCell,
                                                bool dcProjectToIso,
                                                int dcProjectIterations) {
  if (type == ProceduralMeshLabWindowState::MesherType::DualContouring) {
    render::SdfDualContourParams p{};
    // Keep the common mesher settings unified between algorithms.
    p.resolution = base.resolution;
    p.bounds = base.bounds;
    p.iso = base.iso;
    p.computeNormalsFromField = base.computeNormalsFromField;
    p.normalEps = base.normalEps;
    p.fixWindingFromNormals = base.fixWindingFromNormals;

    // Dual Contouring specifics.
    p.qefRegularization = dcQefRegularization;
    p.clampToCell = dcClampToCell;
    p.projectToIso = dcProjectToIso;
    p.projectIterations = dcProjectIterations;

    return render::meshIsosurfaceDualContouring(field, p);
  }

  return render::meshIsosurfaceMarchingTetrahedra(field, base);
}

static void requestRemesh(ProceduralMeshLabWindowState& st) {
  st.lastError.clear();
  st.jobRunning = false;

  // Coalesce: bump job id and overwrite previous job.
  const int jobId = st.nextJobId++;
  st.latestJobId = jobId;

  const render::SdfGraph gCopy = st.graph;
  const render::SdfMesherParams params = st.mesher;

  // Mesher choice (captured for the worker thread).
  const ProceduralMeshLabWindowState::MesherType mesherType = st.mesherType;
  const float dcQef = st.dcQefRegularization;
  const bool dcClamp = st.dcClampToCell;
  const bool dcProject = st.dcProjectToIso;
  const int dcIters = st.dcProjectIterations;

  st.jobRunning = true;
  st.job = std::async(std::launch::async, [jobId, gCopy, params, mesherType, dcQef, dcClamp, dcProject, dcIters]() -> ProceduralMeshLabWindowState::MeshJobResult {
    ProceduralMeshLabWindowState::MeshJobResult out;
    out.id = jobId;

    const auto t0 = std::chrono::high_resolution_clock::now();

    try {
      const auto field = render::makeSdfField(gCopy);
      out.mesh = buildMeshForSettings(field, params, mesherType, dcQef, dcClamp, dcProject, dcIters);
      out.stats = render::measureSdfMesh(out.mesh);

      if (out.mesh.vertices.empty() || out.mesh.indices.empty()) {
        out.error = "Mesher produced an empty mesh (try increasing bounds or changing iso).";
      }
    } catch (const std::exception& e) {
      out.error = std::string("Meshing exception: ") + e.what();
    } catch (...) {
      out.error = "Meshing exception (unknown).";
    }

    const auto t1 = std::chrono::high_resolution_clock::now();
    out.ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return out;
  });
}

// Forward declarations (used before definition).
static void uploadPreviewMeshForCurrentLod(ProceduralMeshLabWindowState& st);
static void resetLodChainToBaseMesh(ProceduralMeshLabWindowState& st);
static void requestBuildLods(ProceduralMeshLabWindowState& st);

static void pollRemeshJob(ProceduralMeshLabWindowState& st, const ToastFn& toast) {
  if (!st.jobRunning) return;
  if (!st.job.valid()) {
    st.jobRunning = false;
    return;
  }

  using namespace std::chrono_literals;
  if (st.job.wait_for(0ms) != std::future_status::ready) return;

  ProceduralMeshLabWindowState::MeshJobResult res = st.job.get();
  st.jobRunning = false;

  // Stale result? discard.
  if (res.id != st.latestJobId) return;

  st.lastMeshMs = res.ms;
  st.lastError = res.error;

  if (res.error.empty()) {
    st.cpuMesh = std::move(res.mesh);
    st.cpuStats = res.stats;
    st.hasMesh = !st.cpuMesh.vertices.empty() && !st.cpuMesh.indices.empty();

    // Reset LOD chain to LOD0 = base mesh.
    // Also invalidate any in-flight LOD jobs from an older base mesh.
    st.latestLodJobId = st.nextLodJobId++;
    resetLodChainToBaseMesh(st);

    // Upload the selected preview LOD (defaults to LOD0).
    if (st.hasMesh) {
      uploadPreviewMeshForCurrentLod(st);
      toast("Procedural mesh rebuilt (" + std::to_string((int)st.cpuStats.vertexCount) + " verts)", 1.25);
    }

    // Auto-build LOD chain if enabled.
    if (st.buildLods && st.autoBuildLods && st.lodLevels > 1) {
      requestBuildLods(st);
    }
  } else {
    st.hasMesh = false;
    resetLodChainToBaseMesh(st);
  }
}

// ---- LOD chain / mesh simplification -------------------------------------

static void uploadPreviewMeshForCurrentLod(ProceduralMeshLabWindowState& st) {
  if (!st.previewInited) return;
  if (!st.hasMesh) return;

  const int maxLod = (int)st.lodMeshes.size(); // LOD0 + (lodMeshes)
  st.previewLod = std::clamp(st.previewLod, 0, std::max(0, maxLod));

  const render::SdfMeshData* m = &st.cpuMesh; // LOD0
  if (st.previewLod > 0) {
    const std::size_t idx = (std::size_t)(st.previewLod - 1);
    if (idx < st.lodMeshes.size()) m = &st.lodMeshes[idx];
  }

  if (m && !m->vertices.empty() && !m->indices.empty()) {
    st.previewMesh.upload(m->vertices, m->indices);
  }
}

static void resetLodChainToBaseMesh(ProceduralMeshLabWindowState& st) {
  st.lodMeshes.clear();
  st.lodStats.clear();
  st.lodSimplifyStats.clear();
  st.lodError.clear();
  st.lastLodMs = 0.0;
  st.hasLods = false;

  // LOD0 is always the base cpuMesh. Clearing the arrays means "LOD0 only".
  if (st.hasMesh) {
    st.previewLod = 0;
    st.exportLod = 0;
    st.exportAllLods = false;
  }
}

static void requestBuildLods(ProceduralMeshLabWindowState& st) {
  st.lodError.clear();

  if (!st.hasMesh || st.cpuMesh.vertices.empty() || st.cpuMesh.indices.empty()) {
    st.lodError = "No base mesh to simplify.";
    return;
  }

  const int jobId = st.nextLodJobId++;
  st.latestLodJobId = jobId;

  // Capture copies for the worker thread.
  render::SdfMeshData baseMesh = st.cpuMesh;
  const int levels = std::clamp(st.lodLevels, 1, 8);
  const float ratio = std::clamp(st.lodRatioPerLevel, 0.05f, 0.99f);

  st.lodJobRunning = true;
  st.lodJob = std::async(std::launch::async, [jobId, baseMesh = std::move(baseMesh), levels, ratio]() -> ProceduralMeshLabWindowState::LodJobResult {
    ProceduralMeshLabWindowState::LodJobResult out;
    out.id = jobId;

    const auto t0 = std::chrono::high_resolution_clock::now();

    try {
      out.meshes.clear();
      out.stats.clear();
      out.simplifyStats.clear();

      // LOD0 is the base cpuMesh (kept in state). This job only produces LOD1..LOD(N).
      const int extraLevels = std::max(0, levels - 1);
      out.meshes.reserve((std::size_t)extraLevels);
      out.stats.reserve((std::size_t)extraLevels);
      out.simplifyStats.reserve((std::size_t)extraLevels);

      render::SdfMeshData prev = baseMesh;
      for (int i = 1; i < levels; ++i) {
        render::MeshSimplifyParams p;
        p.targetTriangleRatio = ratio;
        p.maxIterations = 500000;
        p.recomputeNormals = true;
        p.recomputeSphericalUVs = true;

        render::MeshSimplifyStats simp;
        std::string err;
        render::SdfMeshData simplified = render::simplifyMeshQEM(prev, p, &simp, &err);
        if (!err.empty()) {
          out.error = err;
          break;
        }

        out.meshes.push_back(simplified);
        out.stats.push_back(render::measureSdfMesh(simplified));
        out.simplifyStats.push_back(simp);

        prev = std::move(simplified);
        if (prev.indices.size() / 3u <= 12u) break; // stop when it's already tiny
      }
    } catch (const std::exception& e) {
      out.error = std::string("LOD build exception: ") + e.what();
    } catch (...) {
      out.error = "LOD build exception (unknown).";
    }

    const auto t1 = std::chrono::high_resolution_clock::now();
    out.ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return out;
  });
}

static void pollLodJob(ProceduralMeshLabWindowState& st, const ToastFn& toast) {
  if (!st.lodJobRunning) return;
  if (!st.lodJob.valid()) {
    st.lodJobRunning = false;
    return;
  }

  using namespace std::chrono_literals;
  if (st.lodJob.wait_for(0ms) != std::future_status::ready) return;

  ProceduralMeshLabWindowState::LodJobResult res = st.lodJob.get();
  st.lodJobRunning = false;

  if (res.id != st.latestLodJobId) return; // stale

  st.lastLodMs = res.ms;
  st.lodError = res.error;

  // LOD0 lives in cpuMesh. The job returns LOD1.. (potentially empty).
  st.lodMeshes = std::move(res.meshes);
  st.lodStats = std::move(res.stats);
  st.lodSimplifyStats = std::move(res.simplifyStats);
  st.hasLods = !st.lodMeshes.empty();

  // Clamp selection and upload the chosen LOD for preview.
  const int maxLod = (int)st.lodMeshes.size();
  st.previewLod = std::clamp(st.previewLod, 0, std::max(0, maxLod));
  st.exportLod = std::clamp(st.exportLod, 0, std::max(0, maxLod));
  if (st.previewLod > maxLod) st.previewLod = maxLod;
  if (st.exportLod > maxLod) st.exportLod = maxLod;
  uploadPreviewMeshForCurrentLod(st);

  const int totalLods = 1 + (int)st.lodMeshes.size();
  if (st.lodError.empty()) {
    toast("Built " + std::to_string(totalLods) + " LODs", 1.25);
  } else {
    toast("LOD build finished with warnings", 2.0);
  }
}

static bool uiNodeOpCombo(const char* label, render::SdfNodeOp& op) {
  const render::SdfNodeOp ops[] = {
    // Primitives
    render::SdfNodeOp::Constant,
    render::SdfNodeOp::Sphere,
    render::SdfNodeOp::Box,
    render::SdfNodeOp::Capsule,
    render::SdfNodeOp::TorusY,
    render::SdfNodeOp::Plane,

    // CSG / blending
    render::SdfNodeOp::Union,
    render::SdfNodeOp::SmoothUnion,
    render::SdfNodeOp::Intersect,
    render::SdfNodeOp::Subtract,

    // Modifiers
    render::SdfNodeOp::NoiseDisplace,
    render::SdfNodeOp::NoiseDisplacePerlin,
    render::SdfNodeOp::Shell,

    // Domain transforms
    render::SdfNodeOp::Translate,
    render::SdfNodeOp::RotateX,
    render::SdfNodeOp::RotateY,
    render::SdfNodeOp::RotateZ,
    render::SdfNodeOp::Scale,
    render::SdfNodeOp::Repeat,
    render::SdfNodeOp::Mirror,
    render::SdfNodeOp::TwistY,
  };

  const char* preview = render::sdfNodeOpName(op);
  bool changed = false;
  if (ImGui::BeginCombo(label, preview)) {
    for (auto c : ops) {
      const bool sel = (op == c);
      if (ImGui::Selectable(render::sdfNodeOpName(c), sel)) {
        op = c;
        changed = true;
      }
      if (sel) ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }
  return changed;
}

static bool drawNodeEditor(render::SdfGraph& g, int i) {
  bool changed = false;
  if (i < 0 || i >= (int)g.nodes.size()) return false;

  auto& n = g.nodes[(std::size_t)i];

  ImGui::PushID(i);

  // Header
  {
    std::string title = "#" + std::to_string(i) + " " + render::sdfNodeOpName(n.op);
    const bool open = ImGui::CollapsingHeader(title.c_str(), ImGuiTreeNodeFlags_DefaultOpen);
    if (!open) {
      ImGui::PopID();
      return false;
    }
  }

  ImGui::Indent();

  changed |= uiNodeOpCombo("Op", n.op);

  const char* help = meshOpHelp(n.op);
  if (help && *help) {
    ImGui::TextDisabled("%s", help);
  }

  // Inputs
  switch (n.op) {
    case render::SdfNodeOp::Union:
    case render::SdfNodeOp::SmoothUnion:
    case render::SdfNodeOp::Intersect:
    case render::SdfNodeOp::Subtract:
      changed |= beginComboNodeIndex("A", n.in0, g, i - 1);
      changed |= beginComboNodeIndex("B", n.in1, g, i - 1);
      break;

    case render::SdfNodeOp::NoiseDisplace:
    case render::SdfNodeOp::NoiseDisplacePerlin:
    case render::SdfNodeOp::Shell:
    case render::SdfNodeOp::Translate:
    case render::SdfNodeOp::RotateX:
    case render::SdfNodeOp::RotateY:
    case render::SdfNodeOp::RotateZ:
    case render::SdfNodeOp::Scale:
    case render::SdfNodeOp::Repeat:
    case render::SdfNodeOp::Mirror:
    case render::SdfNodeOp::TwistY:
      changed |= beginComboNodeIndex("In", n.in0, g, i - 1);
      break;

    default:
      break;
  }

  // Params
  switch (n.op) {
    case render::SdfNodeOp::Constant:
      changed |= ImGui::DragFloat("Value", &n.p0, 0.01f);
      break;

    case render::SdfNodeOp::Sphere: {
      float c[3] = {n.p0, n.p1, n.p2};
      if (ImGui::DragFloat3("Center", c, 0.01f)) {
        n.p0 = c[0]; n.p1 = c[1]; n.p2 = c[2];
        changed = true;
      }
      changed |= ImGui::DragFloat("Radius", &n.p3, 0.01f, 0.001f, 10.0f);
    } break;

    case render::SdfNodeOp::Box: {
      float c[3] = {n.p0, n.p1, n.p2};
      float hs[3] = {n.p3, n.p4, n.p5};
      if (ImGui::DragFloat3("Center", c, 0.01f)) {
        n.p0 = c[0]; n.p1 = c[1]; n.p2 = c[2];
        changed = true;
      }
      if (ImGui::DragFloat3("HalfSize", hs, 0.01f, 0.001f, 10.0f)) {
        n.p3 = std::max(0.001f, hs[0]);
        n.p4 = std::max(0.001f, hs[1]);
        n.p5 = std::max(0.001f, hs[2]);
        changed = true;
      }
    } break;

    case render::SdfNodeOp::Capsule: {
      float a[3] = {n.p0, n.p1, n.p2};
      float b[3] = {n.p3, n.p4, n.p5};
      if (ImGui::DragFloat3("A", a, 0.01f)) {
        n.p0 = a[0]; n.p1 = a[1]; n.p2 = a[2];
        changed = true;
      }
      if (ImGui::DragFloat3("B", b, 0.01f)) {
        n.p3 = b[0]; n.p4 = b[1]; n.p5 = b[2];
        changed = true;
      }
      changed |= ImGui::DragFloat("Radius", &n.p6, 0.01f, 0.001f, 10.0f);
    } break;

    case render::SdfNodeOp::TorusY: {
      float c[3] = {n.p0, n.p1, n.p2};
      if (ImGui::DragFloat3("Center", c, 0.01f)) {
        n.p0 = c[0]; n.p1 = c[1]; n.p2 = c[2];
        changed = true;
      }
      changed |= ImGui::DragFloat("MajorR", &n.p3, 0.01f, 0.001f, 10.0f);
      changed |= ImGui::DragFloat("MinorR", &n.p4, 0.01f, 0.001f, 10.0f);
    } break;

    case render::SdfNodeOp::Plane: {
      float nn[3] = {n.p0, n.p1, n.p2};
      if (ImGui::DragFloat3("Normal", nn, 0.01f, -1.0f, 1.0f)) {
        n.p0 = nn[0]; n.p1 = nn[1]; n.p2 = nn[2];
        changed = true;
      }
      changed |= ImGui::DragFloat("Offset", &n.p3, 0.01f);
      if (ImGui::Button("Normalize normal")) {
        const float len = std::sqrt(n.p0*n.p0 + n.p1*n.p1 + n.p2*n.p2);
        if (len > 1e-6f) {
          n.p0 /= len; n.p1 /= len; n.p2 /= len;
          changed = true;
        }
      }
    } break;

    case render::SdfNodeOp::SmoothUnion:
      changed |= ImGui::DragFloat("k", &n.p0, 0.01f, 0.0001f, 2.0f);
      break;

    case render::SdfNodeOp::NoiseDisplace:
    case render::SdfNodeOp::NoiseDisplacePerlin: {
      changed |= ImGui::DragFloat("Freq", &n.p0, 0.01f, 0.0f, 64.0f);
      changed |= ImGui::DragFloat("Amp", &n.p1, 0.01f, 0.0f, 2.0f);
      int oct = (int)std::round(n.p2);
      if (ImGui::SliderInt("Octaves", &oct, 1, 10)) {
        n.p2 = (float)oct;
        changed = true;
      }
      changed |= ImGui::DragFloat("Lacunarity", &n.p3, 0.01f, 1.0f, 6.0f);
      changed |= ImGui::DragFloat("Gain", &n.p4, 0.01f, 0.0f, 1.0f);
      float off[3] = {n.p5, n.p6, n.p7};
      if (ImGui::DragFloat3("Offset", off, 0.01f)) {
        n.p5 = off[0]; n.p6 = off[1]; n.p7 = off[2];
        changed = true;
      }
      ImGui::TextDisabled("Node seed: %llu", (unsigned long long)n.seed);
      if (ImGui::Button("Randomize node seed")) {
        n.seed = core::hashCombine(n.seed ? n.seed : 1ull, core::fnv1a64("noise"));
        changed = true;
      }
    } break;

    case render::SdfNodeOp::Shell:
      changed |= ImGui::DragFloat("Thickness", &n.p0, 0.01f, 0.0f, 2.0f);
      break;

    case render::SdfNodeOp::Translate: {
      float t[3] = {n.p0, n.p1, n.p2};
      if (ImGui::DragFloat3("Translate", t, 0.01f)) {
        n.p0 = t[0]; n.p1 = t[1]; n.p2 = t[2];
        changed = true;
      }
    } break;

    case render::SdfNodeOp::RotateX:
    case render::SdfNodeOp::RotateY:
    case render::SdfNodeOp::RotateZ: {
      changed |= ImGui::DragFloat("Angle (deg)", &n.p0, 0.5f, -360.0f, 360.0f);
      float pivot[3] = {n.p1, n.p2, n.p3};
      if (ImGui::DragFloat3("Pivot", pivot, 0.01f)) {
        n.p1 = pivot[0]; n.p2 = pivot[1]; n.p3 = pivot[2];
        changed = true;
      }
      if (ImGui::Button("Pivot = origin")) {
        n.p1 = n.p2 = n.p3 = 0.0f;
        changed = true;
      }
    } break;

    case render::SdfNodeOp::Scale: {
      changed |= ImGui::DragFloat("Scale", &n.p0, 0.01f, 0.001f, 100.0f);
      n.p0 = std::max(0.001f, std::fabs(n.p0));

      float pivot[3] = {n.p1, n.p2, n.p3};
      if (ImGui::DragFloat3("Pivot", pivot, 0.01f)) {
        n.p1 = pivot[0]; n.p2 = pivot[1]; n.p3 = pivot[2];
        changed = true;
      }
      if (ImGui::Button("Pivot = origin")) {
        n.p1 = n.p2 = n.p3 = 0.0f;
        changed = true;
      }
    } break;

    case render::SdfNodeOp::Repeat: {
      float per[3] = {n.p0, n.p1, n.p2};
      if (ImGui::DragFloat3("Period", per, 0.01f, 0.001f, 100.0f)) {
        n.p0 = std::max(0.001f, per[0]);
        n.p1 = std::max(0.001f, per[1]);
        n.p2 = std::max(0.001f, per[2]);
        changed = true;
      }

      float off[3] = {n.p3, n.p4, n.p5};
      if (ImGui::DragFloat3("Offset", off, 0.01f)) {
        n.p3 = off[0]; n.p4 = off[1]; n.p5 = off[2];
        changed = true;
      }

      if (ImGui::Button("Offset = 0")) {
        n.p3 = n.p4 = n.p5 = 0.0f;
        changed = true;
      }
    } break;

    case render::SdfNodeOp::Mirror: {
      bool mx = n.p0 > 0.5f;
      bool my = n.p1 > 0.5f;
      bool mz = n.p2 > 0.5f;

      bool mChanged = false;
      mChanged |= ImGui::Checkbox("Mirror X", &mx);
      ImGui::SameLine();
      mChanged |= ImGui::Checkbox("Mirror Y", &my);
      ImGui::SameLine();
      mChanged |= ImGui::Checkbox("Mirror Z", &mz);
      if (mChanged) {
        n.p0 = mx ? 1.0f : 0.0f;
        n.p1 = my ? 1.0f : 0.0f;
        n.p2 = mz ? 1.0f : 0.0f;
        changed = true;
      }

      float pivot[3] = {n.p3, n.p4, n.p5};
      if (ImGui::DragFloat3("Pivot", pivot, 0.01f)) {
        n.p3 = pivot[0]; n.p4 = pivot[1]; n.p5 = pivot[2];
        changed = true;
      }
      if (ImGui::Button("Pivot = origin")) {
        n.p3 = n.p4 = n.p5 = 0.0f;
        changed = true;
      }
    } break;

    case render::SdfNodeOp::TwistY: {
      changed |= ImGui::DragFloat("Twist (deg/unit)", &n.p0, 0.25f, -720.0f, 720.0f);
      float pivot[3] = {n.p1, n.p2, n.p3};
      if (ImGui::DragFloat3("Pivot", pivot, 0.01f)) {
        n.p1 = pivot[0]; n.p2 = pivot[1]; n.p3 = pivot[2];
        changed = true;
      }
      if (ImGui::Button("Pivot = origin")) {
        n.p1 = n.p2 = n.p3 = 0.0f;
        changed = true;
      }
    } break;

    default:
      break;
  }

  ImGui::Unindent();
  ImGui::PopID();
  return changed;
}

static void addNode(render::SdfGraph& g, render::SdfNodeOp op) {
  if ((int)g.nodes.size() >= render::kSdfGraphMaxNodes) return;
  render::SdfNode n;
  n.op = op;

  // Give sane defaults.
  switch (op) {
    case render::SdfNodeOp::Sphere:
      n.p3 = 1.0f;
      break;
    case render::SdfNodeOp::Box:
      n.p3 = n.p4 = n.p5 = 0.5f;
      break;
    case render::SdfNodeOp::Capsule:
      n.p2 = -0.5f;
      n.p5 = 0.5f;
      n.p6 = 0.25f;
      break;
    case render::SdfNodeOp::TorusY:
      n.p3 = 0.75f;
      n.p4 = 0.25f;
      break;
    case render::SdfNodeOp::Plane:
      n.p1 = 1.0f;
      n.p3 = 0.0f;
      break;
    case render::SdfNodeOp::SmoothUnion:
      n.p0 = 0.25f;
      break;
    case render::SdfNodeOp::NoiseDisplace:
    case render::SdfNodeOp::NoiseDisplacePerlin:
      n.p0 = 3.0f;
      n.p1 = 0.15f;
      n.p2 = 5.0f;
      n.p3 = 2.0f;
      n.p4 = 0.5f;
      n.seed = 1337;
      break;
    case render::SdfNodeOp::Shell:
      n.p0 = 0.05f;
      break;

    case render::SdfNodeOp::Translate:
      // Translate vector defaults to 0.
      break;

    case render::SdfNodeOp::RotateX:
    case render::SdfNodeOp::RotateY:
    case render::SdfNodeOp::RotateZ:
      n.p0 = 45.0f; // deg
      break;

    case render::SdfNodeOp::Scale:
      n.p0 = 1.0f;
      break;

    case render::SdfNodeOp::Repeat:
      n.p0 = n.p1 = n.p2 = 2.0f; // period
      // offset p3..p5 = 0
      break;

    case render::SdfNodeOp::Mirror:
      n.p0 = 1.0f; // mirror X
      n.p1 = 0.0f;
      n.p2 = 0.0f;
      // pivot p3..p5 = 0
      break;

    case render::SdfNodeOp::TwistY:
      n.p0 = 45.0f; // deg/unit
      break;

    default:
      break;
  }

  g.nodes.push_back(n);
  g.output = (int)g.nodes.size() - 1;
}

static void removeNode(render::SdfGraph& g, int idx) {
  if (idx < 0 || idx >= (int)g.nodes.size()) return;
  g.nodes.erase(g.nodes.begin() + idx);

  // Fix up indices.
  for (auto& n : g.nodes) {
    auto fix = [&](int& in) {
      if (in == idx) in = -1;
      else if (in > idx) --in;
    };
    fix(n.in0);
    fix(n.in1);
  }
  if (g.output == idx) g.output = (int)g.nodes.size() - 1;
  else if (g.output > idx) --g.output;
}

} // namespace

void drawProceduralMeshLabWindow(ProceduralMeshLabWindowState& st, float timeSec, const ToastFn& toast) {
  if (!st.open) return;

  // One-time init defaults that are nicer than SdfMesherParams defaults.
  if (st.nextJobId == 1 && st.mesher.resolution == 32) {
    st.mesher.resolution = 64;
    st.mesher.bounds = 1.25f;
    st.mesher.iso = 0.0f;
    st.mesher.normalEps = 0.0025f;
    // New default: Dual Contouring tends to preserve sharper CSG edges than Marching Tets.
    st.mesherType = ProceduralMeshLabWindowState::MesherType::DualContouring;
  }

  if (!ensurePreview(st)) {
    // Keep the window usable even without preview.
  }

  pollRemeshJob(st, toast);
  pollLodJob(st, toast);

  ImGui::SetNextWindowSize(ImVec2(980, 720), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Procedural Mesh Lab", &st.open)) {
    ImGui::End();
    return;
  }

  ImGui::TextUnformatted("Signed-distance-field (SDF) node graph -> (Marching Tetrahedra | Dual Contouring) mesh -> OBJ export.");
  ImGui::Separator();

  // Tooling History (Undo / Redo)
  //
  // The mesh lab is an iteration-heavy tool; having a safe history makes it much
  // harder to accidentally lose a good shape while experimenting.
  {
    const bool focused = ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows);
    ImGuiIO& io = ImGui::GetIO();
    const bool wantText = io.WantTextInput;

    // Keyboard shortcuts (only when this window is focused and we're not typing into a field).
    if (focused && !wantText && st.historyEnabled) {
      const bool ctrl = io.KeyCtrl;
      const bool shift = io.KeyShift;

      if (ctrl && ImGui::IsKeyPressed(ImGuiKey_Z, false)) {
        if (shift) {
          if (historyRedo(st, timeSec)) toast("Redo", 1.0f);
        } else {
          if (historyUndo(st, timeSec)) toast("Undo", 1.0f);
        }
      } else if (ctrl && ImGui::IsKeyPressed(ImGuiKey_Y, false)) {
        if (historyRedo(st, timeSec)) toast("Redo", 1.0f);
      }
    }

    ImGui::SeparatorText("History");
    const bool canUndo = st.historyEnabled && (st.historyPending || !st.historyUndo.empty());
    const bool canRedo = st.historyEnabled && !st.historyRedo.empty();

    ImGui::BeginDisabled(!canUndo);
    if (ImGui::Button("Undo##mesh_lab")) {
      if (historyUndo(st, timeSec)) toast("Undo", 1.0f);
    }
    ImGui::EndDisabled();

    ImGui::SameLine();

    ImGui::BeginDisabled(!canRedo);
    if (ImGui::Button("Redo##mesh_lab")) {
      if (historyRedo(st, timeSec)) toast("Redo", 1.0f);
    }
    ImGui::EndDisabled();

    ImGui::SameLine();

    bool enabled = st.historyEnabled;
    if (ImGui::Checkbox("Enable##mesh_lab_history", &enabled)) {
      st.historyEnabled = enabled;
      // Reset stacks on toggle so behavior stays predictable.
      historyClear(st, timeSec);
      toast(st.historyEnabled ? "History enabled" : "History disabled", 1.0f);
    }

    ImGui::SameLine();
    ImGui::SetNextItemWidth(90.0f);
    if (ImGui::DragInt("Max##mesh_lab_history", &st.historyMax, 1.0f, 4, 256)) {
      st.historyMax = std::clamp(st.historyMax, 4, 256);
      historyTrim(st.historyUndo, st.historyMax);
      historyTrim(st.historyRedo, st.historyMax);
    }

    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    if (ImGui::DragFloat("Coalesce##mesh_lab_history", &st.historyCoalesceSec, 0.01f, 0.0f, 2.0f, "%.2f s")) {
      if (st.historyCoalesceSec < 0.0f) st.historyCoalesceSec = 0.0f;
    }

    ImGui::SameLine();
    if (ImGui::Button("Clear##mesh_lab_history")) {
      historyClear(st, timeSec);
      toast("History cleared", 1.0f);
    }

    ImGui::TextDisabled("Ctrl+Z undo, Ctrl+Y / Ctrl+Shift+Z redo. (Pending changes commit once you stop editing.)");
    ImGui::Separator();
  }

  bool needApplyPreset = false;

  // Preset + seed row
  {
    if (ImGui::BeginCombo("Preset", render::sdfGraphPresetName(st.preset))) {
      for (int i = 0; i <= (int)render::SdfGraphPreset::BooleanDemo; ++i) {
        auto p = (render::SdfGraphPreset)i;
        const bool sel = (p == st.preset);
        if (ImGui::Selectable(render::sdfGraphPresetName(p), sel)) {
          st.preset = p;
          needApplyPreset = true;
        }
        if (sel) ImGui::SetItemDefaultFocus();
      }
      ImGui::EndCombo();
    }

    ImGui::SameLine();
    ImGui::Checkbox("Lock", &st.lockToPreset);

    core::u64 seed = st.seed;
    if (ImGui::InputScalar("Seed", ImGuiDataType_U64, &seed)) {
      st.seed = seed;
      if (st.lockToPreset) needApplyPreset = true;
    }

    ImGui::SameLine();
    if (ImGui::Button("Random")) {
      st.seed = core::hashCombine(st.seed ? st.seed : 1ull, (core::u64)std::chrono::high_resolution_clock::now().time_since_epoch().count());
      if (st.lockToPreset) needApplyPreset = true;
      st.dirty = true;
    }

    if (!st.lockToPreset) {
      ImGui::SameLine();
      if (ImGui::Button("Reset to preset")) {
        needApplyPreset = true;
      }
    }
  }

  if (needApplyPreset) {
    st.graph = render::makeSdfGraphPreset(st.preset, st.seed);
    st.dirty = true;
  }

  // Graph save/load (asset pipeline)
  {
    ImGui::SeparatorText("Graph I/O");
    ImGui::InputText("SDF graph file", st.graphPath, sizeof(st.graphPath));

    if (ImGui::Button("Save SDF graph")) {
      std::string err;
      if (!render::saveSdfGraphToFile(st.graph, st.graphPath, &err)) {
        toast(std::string("Save failed: ") + err, 4.0);
      } else {
        toast(std::string("Saved ") + st.graphPath, 3.0);
      }
    }
    ImGui::SameLine();
    if (ImGui::Button("Load SDF graph")) {
      render::SdfGraph loaded;
      std::string err;
      if (!render::loadSdfGraphFromFile(st.graphPath, loaded, &err)) {
        toast(std::string("Load failed: ") + err, 4.0);
      } else {
        st.graph = std::move(loaded);
        st.seed = st.graph.seed;
        st.lockToPreset = false;
        st.dirty = true;
        toast(std::string("Loaded ") + st.graphPath, 3.0);
      }
    }
  }


  // -----------------------------------------------------------------------------
  // Graph diagnostics / maintenance
  // -----------------------------------------------------------------------------
  {
    if (ImGui::CollapsingHeader("Graph Diagnostics", ImGuiTreeNodeFlags_DefaultOpen)) {
      const int nodeCount = (int)st.graph.nodes.size();
      const bool outputValid = (st.graph.output >= 0 && st.graph.output < nodeCount);

      int invalidInputs = 0;
      int forwardRefs = 0;

      std::vector<unsigned char> reachable;
      reachable.assign((std::size_t)std::max(0, nodeCount), 0u);
      if (outputValid) {
        markReachableRec(st.graph, st.graph.output, reachable);
      }

      int reachableCount = 0;
      for (int i = 0; i < nodeCount; ++i) {
        if (reachable[(std::size_t)i]) reachableCount++;
        const auto& n = st.graph.nodes[(std::size_t)i];
        if (n.in0 < -1 || n.in0 >= nodeCount) invalidInputs++;
        if (n.in1 < -1 || n.in1 >= nodeCount) invalidInputs++;
        if (n.in0 >= i && n.in0 != -1) forwardRefs++;
        if (n.in1 >= i && n.in1 != -1) forwardRefs++;
      }

      const bool hasCycle = sdfGraphHasCycle(st.graph);

      ImGui::Text("Nodes: %d", nodeCount);
      ImGui::SameLine();
      if (outputValid) {
        ImGui::Text("Output: %d", st.graph.output);
      } else {
        ImGui::TextColored(ImVec4(1, 0.3f, 0.3f, 1), "Output: %d (INVALID)", st.graph.output);
      }

      if (outputValid) {
        ImGui::Text("Reachable: %d   Unreachable: %d", reachableCount, nodeCount - reachableCount);
      } else {
        ImGui::TextDisabled("(Reachability unavailable: invalid output index)");
      }

      if (invalidInputs > 0) {
        ImGui::TextColored(ImVec4(1, 0.7f, 0.2f, 1), "Invalid input indices: %d", invalidInputs);
      }
      if (forwardRefs > 0) {
        ImGui::TextColored(ImVec4(1, 0.7f, 0.2f, 1), "Forward references: %d (editor expects DAG)", forwardRefs);
      }
      if (hasCycle) {
        ImGui::TextColored(ImVec4(1, 0.5f, 0.2f, 1), "Cycle detected (eval has a safety guard, but results may be wrong)");
      }

      if (ImGui::Button("Fix indices")) {
        fixSdfGraphIndices(st.graph);
        st.seed = st.graph.seed;
        st.lockToPreset = false;
        st.dirty = true;
        toast("Graph indices fixed", 2.0);
      }

      ImGui::SameLine();
      if (ImGui::Button("Prune unreachable")) {
        render::SdfGraph pruned = pruneUnreachableNodes(st.graph);
        if (pruned.nodes.size() != st.graph.nodes.size()) {
          const int removed = (int)st.graph.nodes.size() - (int)pruned.nodes.size();
          st.graph = std::move(pruned);
          st.seed = st.graph.seed;
          st.lockToPreset = false;
          st.dirty = true;
          toast("Pruned " + std::to_string(removed) + " unreachable node(s)", 2.5);
        } else {
          toast("No unreachable nodes to prune", 1.2);
        }
      }

      ImGui::SameLine();
      if (ImGui::Button("Clear variants")) {
        st.variants.clear();
        st.selectedVariant = -1;
        st.variantsError.clear();
        toast("Variants cleared", 1.0);
      }
    }
  }

  // -----------------------------------------------------------------------------
  // Variation Studio (SDF ideation)
  // -----------------------------------------------------------------------------
  {
    // Use saved preference as the default open state (first time only).
    ImGui::SetNextItemOpen(st.showVariationStudio, ImGuiCond_FirstUseEver);
    const bool open = ImGui::CollapsingHeader("Variation Studio");
    st.showVariationStudio = open;

    if (open) {
      ImGui::TextDisabled("Generate mutated variants of the current SDF graph (fast ideation).");

      const float bounds = std::max(0.1f, st.mesher.bounds);
      const float iso = st.mesher.iso;
      const int statsRes = std::clamp(st.variationSampleRes, 4, 64);

      const SdfSampleStats baseStats = sampleSdfGraphStats(st.graph, bounds, iso, statsRes);
      ImGui::Text("Base stats: inside %.1f%%  R~%.2f  d[min %.2f, max %.2f]  (%d samples)",
                  baseStats.insideFraction * 100.0f,
                  baseStats.boundingRadius,
                  baseStats.minDistance,
                  baseStats.maxDistance,
                  baseStats.samples);

      bool regen = false;

      core::u64 vSeed = st.variationSeed;
      if (ImGui::InputScalar("Variation seed", ImGuiDataType_U64, &vSeed)) {
        st.variationSeed = vSeed;
        regen |= st.variationAutoRegenerate;
      }
      ImGui::SameLine();
      if (ImGui::Button("Random##var_seed")) {
        st.variationSeed = randomU64();
        regen |= true;
      }

      int count = st.variationCount;
      if (ImGui::SliderInt("Count", &count, 1, 24)) {
        st.variationCount = count;
        regen |= st.variationAutoRegenerate;
      }

      if (ImGui::DragFloat("Param jitter", &st.variationParamJitter, 0.01f, 0.0f, 1.5f, "%.2f")) {
        if (st.variationParamJitter < 0.0f) st.variationParamJitter = 0.0f;
        regen |= st.variationAutoRegenerate;
      }

      if (ImGui::DragFloat("Rewire chance", &st.variationRewireChance, 0.005f, 0.0f, 1.0f, "%.3f")) {
        st.variationRewireChance = std::clamp(st.variationRewireChance, 0.0f, 1.0f);
        regen |= st.variationAutoRegenerate;
      }

      if (ImGui::DragFloat("Op mutate chance", &st.variationOpChance, 0.005f, 0.0f, 1.0f, "%.3f")) {
        st.variationOpChance = std::clamp(st.variationOpChance, 0.0f, 1.0f);
        regen |= st.variationAutoRegenerate;
      }

      if (ImGui::DragFloat("Wrap output chance", &st.variationWrapChance, 0.005f, 0.0f, 1.0f, "%.3f")) {
        st.variationWrapChance = std::clamp(st.variationWrapChance, 0.0f, 1.0f);
        regen |= st.variationAutoRegenerate;
      }

      ImGui::Checkbox("Mutate node seeds", &st.variationMutateNodeSeeds);
      ImGui::SameLine();
      if (ImGui::Checkbox("Auto regenerate", &st.variationAutoRegenerate)) {
        // no-op (flag is used by the edits above)
      }

      int sr = st.variationSampleRes;
      if (ImGui::SliderInt("Stats grid", &sr, 4, 32)) {
        st.variationSampleRes = sr;
        regen |= st.variationAutoRegenerate;
      }

      const float maxTargetR = std::sqrt(3.0f) * bounds;
      if (ImGui::DragFloat("Target fill", &st.variationTargetFill, 0.01f, 0.0f, 1.0f, "%.2f")) {
        st.variationTargetFill = std::clamp(st.variationTargetFill, 0.0f, 1.0f);
        regen |= st.variationAutoRegenerate;
      }
      if (ImGui::DragFloat("Target radius", &st.variationTargetRadius, 0.01f, 0.0f, maxTargetR, "%.2f")) {
        st.variationTargetRadius = std::clamp(st.variationTargetRadius, 0.0f, maxTargetR);
        regen |= st.variationAutoRegenerate;
      }
      ImGui::Checkbox("Sort by score", &st.variationSortByScore);

      auto doGenerate = [&]() {
        st.variantsError.clear();
        st.selectedVariant = -1;
        st.variants.clear();

        if (st.graph.nodes.empty()) {
          st.variantsError = "Current graph is empty.";
          return;
        }

        std::unordered_set<core::u64> seen;
        seen.reserve((std::size_t)std::max(16, st.variationCount) * 2u);

        const core::u64 baseKey = sdfGraphFingerprint(st.graph);
        seen.insert(baseKey);

        const float targetFill = std::clamp(st.variationTargetFill, 0.0f, 1.0f);
        const float targetR = std::max(0.0f, st.variationTargetRadius);

        auto scoreVariant = [&](const SdfSampleStats& ss) -> float {
          // Heuristic score:
          //  - prefer fill ratio near targetFill (keeps "volume" in a useful range)
          //  - prefer a radius near targetR (keeps size consistent with your bounds)
          //  - require that the iso-surface exists in the sampled region (otherwise it's "empty" or "solid")
          const float fillDen = std::max(0.05f, targetFill);
          float fillScore = 1.0f - (std::fabs(ss.insideFraction - targetFill) / fillDen);
          fillScore = std::clamp(fillScore, 0.0f, 1.0f);

          float rScore = 0.0f;
          if (targetR > 0.0f) {
            const float rDen = std::max(0.05f, targetR);
            rScore = 1.0f - (std::fabs(ss.boundingRadius - targetR) / rDen);
            rScore = std::clamp(rScore, 0.0f, 1.0f);
          }

          const float surfaceScore = (ss.minDistance <= iso && ss.maxDistance >= iso) ? 1.0f : 0.0f;
          return 0.55f * fillScore + 0.35f * rScore + 0.10f * surfaceScore;
        };

        const int want = std::clamp(st.variationCount, 1, 64);
        const int maxTries = want * 6;

        int produced = 0;
        for (int i = 0; i < maxTries && produced < want; ++i) {
          const core::u64 seed = core::hashCombine(st.graph.seed, core::hashCombine(st.variationSeed, (core::u64)i));
          render::SdfGraph gVar = makeSdfVariant(st.graph,
                                                seed,
                                                bounds,
                                                st.variationParamJitter,
                                                st.variationRewireChance,
                                                st.variationOpChance,
                                                st.variationWrapChance,
                                                st.variationMutateNodeSeeds);

          const core::u64 key = sdfGraphFingerprint(gVar);
          if (seen.find(key) != seen.end()) continue;
          seen.insert(key);

          const SdfSampleStats ss = sampleSdfGraphStats(gVar, bounds, iso, statsRes);

          ProceduralMeshLabWindowState::SdfVariantCandidate cand;
          cand.graph = std::move(gVar);
          cand.key = key;
          cand.insideFraction = ss.insideFraction;
          cand.boundingRadius = ss.boundingRadius;
          cand.score = scoreVariant(ss);
          cand.minDistance = ss.minDistance;
          cand.maxDistance = ss.maxDistance;
          cand.samples = ss.samples;

          st.variants.push_back(std::move(cand));
          produced++;
        }

        if (st.variationSortByScore && !st.variants.empty()) {
          std::sort(st.variants.begin(),
                    st.variants.end(),
                    [](const ProceduralMeshLabWindowState::SdfVariantCandidate& a,
                       const ProceduralMeshLabWindowState::SdfVariantCandidate& b) { return a.score > b.score; });
        }

        if (st.variants.empty()) {
          st.variantsError = "No unique variants generated (try increasing jitter or seed).";
        }
      };

      if (ImGui::Button("Generate variants")) {
        doGenerate();
      }

      if (regen && st.variationAutoRegenerate) {
        doGenerate();
      }

      if (!st.variantsError.empty()) {
        ImGui::TextColored(ImVec4(1, 0.35f, 0.35f, 1), "Variant error: %s", st.variantsError.c_str());
      }

      if (!st.variants.empty()) {
        ImGui::SeparatorText("Candidates");

        const ImGuiTableFlags flags =
          ImGuiTableFlags_Borders |
          ImGuiTableFlags_RowBg |
          ImGuiTableFlags_Resizable |
          ImGuiTableFlags_ScrollY;

        const float tableH = std::min(300.0f, ImGui::GetTextLineHeightWithSpacing() * (float)(st.variants.size() + 2));
        if (ImGui::BeginTable("##sdf_variants", 8, flags, ImVec2(0, tableH))) {
          ImGui::TableSetupScrollFreeze(0, 1);
          ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 24.0f);
          ImGui::TableSetupColumn("Apply", ImGuiTableColumnFlags_WidthFixed, 58.0f);
          ImGui::TableSetupColumn("Nodes", ImGuiTableColumnFlags_WidthFixed, 52.0f);
          ImGui::TableSetupColumn("Inside", ImGuiTableColumnFlags_WidthFixed, 62.0f);
          ImGui::TableSetupColumn("R~", ImGuiTableColumnFlags_WidthFixed, 54.0f);
          ImGui::TableSetupColumn("Score", ImGuiTableColumnFlags_WidthFixed, 54.0f);
          ImGui::TableSetupColumn("d[min,max]", ImGuiTableColumnFlags_WidthStretch);
          ImGui::TableSetupColumn("Key", ImGuiTableColumnFlags_WidthFixed, 120.0f);
          ImGui::TableHeadersRow();

          for (int i = 0; i < (int)st.variants.size(); ++i) {
            auto& v = st.variants[(std::size_t)i];

            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%d", i);

            ImGui::TableSetColumnIndex(1);
            std::string btn = "Use##var_" + std::to_string(i);
            if (ImGui::SmallButton(btn.c_str())) {
              st.graph = v.graph;
              st.seed = st.graph.seed;
              st.lockToPreset = false;
              st.dirty = true;
              toast("Applied variant " + std::to_string(i), 1.5);
            }

            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%d", (int)v.graph.nodes.size());

            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%.1f%%", v.insideFraction * 100.0f);

            ImGui::TableSetColumnIndex(4);
            ImGui::Text("%.2f", v.boundingRadius);

            ImGui::TableSetColumnIndex(5);
            ImGui::Text("%.2f", v.score);

            ImGui::TableSetColumnIndex(6);
            ImGui::Text("[%.2f, %.2f]", v.minDistance, v.maxDistance);

            ImGui::TableSetColumnIndex(7);
            ImGui::Text("0x%llX", (unsigned long long)v.key);
          }

          ImGui::EndTable();
        }
      }
    }
  }

  // Mesher controls
  {
    ImGui::SeparatorText("Mesher");
    // Algorithm selection
    {
      const char* modes[] = {"Marching Tetrahedra", "Dual Contouring (QEF)"};
      int mode = (st.mesherType == ProceduralMeshLabWindowState::MesherType::DualContouring) ? 1 : 0;
      if (ImGui::Combo("Mesher algorithm", &mode, modes, IM_ARRAYSIZE(modes))) {
        st.mesherType = (mode == 1) ? ProceduralMeshLabWindowState::MesherType::DualContouring
                                 : ProceduralMeshLabWindowState::MesherType::MarchingTetrahedra;
        st.dirty = true;
      }
    }


    int res = st.mesher.resolution;
    if (ImGui::SliderInt("Resolution", &res, 16, 192)) {
      st.mesher.resolution = res;
      st.dirty = true;
    }

    if (ImGui::DragFloat("Bounds", &st.mesher.bounds, 0.01f, 0.1f, 10.0f)) {
      st.mesher.bounds = std::max(0.1f, st.mesher.bounds);
      st.dirty = true;
    }

    if (ImGui::DragFloat("Iso", &st.mesher.iso, 0.001f, -2.0f, 2.0f)) {
      st.dirty = true;
    }

    ImGui::Checkbox("Compute normals from field", &st.mesher.computeNormalsFromField);
    if (ImGui::IsItemDeactivatedAfterEdit()) st.dirty = true;

    if (st.mesher.computeNormalsFromField) {
      if (ImGui::DragFloat("Normal eps", &st.mesher.normalEps, 0.0001f, 1e-5f, 0.1f, "%.5f")) {
        st.mesher.normalEps = std::max(1e-6f, st.mesher.normalEps);
        st.dirty = true;
      }
    }

    ImGui::Checkbox("Fix winding from normals", &st.mesher.fixWindingFromNormals);
    if (ImGui::IsItemDeactivatedAfterEdit()) st.dirty = true;

    if (st.mesherType == ProceduralMeshLabWindowState::MesherType::DualContouring) {
      ImGui::SeparatorText("Dual Contouring");
      if (ImGui::DragFloat("QEF regularization", &st.dcQefRegularization, 1e-7f, 0.0f, 1e-2f, "%.7f")) {
        st.dcQefRegularization = std::max(0.0f, st.dcQefRegularization);
        st.dirty = true;
      }
      ImGui::Checkbox("Clamp vertex to cell", &st.dcClampToCell);
      if (ImGui::IsItemDeactivatedAfterEdit()) st.dirty = true;
      ImGui::Checkbox("Project vertex to iso", &st.dcProjectToIso);
      if (ImGui::IsItemDeactivatedAfterEdit()) st.dirty = true;
      if (st.dcProjectToIso) {
        int it = st.dcProjectIterations;
        if (ImGui::SliderInt("Projection iterations", &it, 0, 8)) {
          st.dcProjectIterations = std::max(0, it);
          st.dirty = true;
        }
      }
      ImGui::Spacing();
      ImGui::TextDisabled("Dual Contouring places 1 vertex per cell and solves it from Hermite samples (edge hits + normals) via a tiny QEF.");
    }

    ImGui::Checkbox("Auto-remesh", &st.autoRemesh);
    ImGui::SameLine();
    ImGui::Checkbox("Async", &st.asyncRemesh);

    if (ImGui::Button("Remesh now")) {
      st.dirty = true;
      st.autoRemesh = false;

      // Fire immediately.
      if (st.asyncRemesh) {
        requestRemesh(st);
      } else {
        // Sync path.
        const auto t0 = std::chrono::high_resolution_clock::now();
        st.cpuMesh = buildMeshForSettings(render::makeSdfField(st.graph),
                                            st.mesher,
                                            st.mesherType,
                                            st.dcQefRegularization,
                                            st.dcClampToCell,
                                            st.dcProjectToIso,
                                            st.dcProjectIterations);
        st.cpuStats = render::measureSdfMesh(st.cpuMesh);
        const auto t1 = std::chrono::high_resolution_clock::now();
        st.lastMeshMs = std::chrono::duration<double, std::milli>(t1 - t0).count();

        st.hasMesh = !st.cpuMesh.vertices.empty() && !st.cpuMesh.indices.empty();
        if (st.hasMesh) {
          // Reset LOD chain to LOD0 and (optionally) rebuild LODs.
          st.latestLodJobId = st.nextLodJobId++;
          resetLodChainToBaseMesh(st);
          uploadPreviewMeshForCurrentLod(st);
          toast("Procedural mesh rebuilt (sync)", 1.0);

          if (st.buildLods && st.autoBuildLods && st.lodLevels > 1) {
            requestBuildLods(st);
          }
        } else {
          st.lastError = "Mesher produced an empty mesh.";
        }
      }
    }

    if (st.jobRunning) {
      ImGui::SameLine();
      ImGui::TextUnformatted("Computing…");
    }

    if (!st.lastError.empty()) {
      ImGui::TextColored(ImVec4(1, 0.25f, 0.25f, 1), "Error: %s", st.lastError.c_str());
    }

    if (st.hasMesh) {
      ImGui::Text("Verts: %zu  Tris: %zu  Build: %.2f ms",
                  st.cpuStats.vertexCount,
                  st.cpuStats.indexCount / 3,
                  st.lastMeshMs);
    }
  }

  // Auto remesh trigger (coalesces while a job is running).
  if (st.autoRemesh && st.dirty && !st.jobRunning) {
    st.dirty = false;
    if (st.asyncRemesh) {
      requestRemesh(st);
    } else {
      const auto t0 = std::chrono::high_resolution_clock::now();
      st.cpuMesh = buildMeshForSettings(render::makeSdfField(st.graph),
                                          st.mesher,
                                          st.mesherType,
                                          st.dcQefRegularization,
                                          st.dcClampToCell,
                                          st.dcProjectToIso,
                                          st.dcProjectIterations);
      st.cpuStats = render::measureSdfMesh(st.cpuMesh);
      const auto t1 = std::chrono::high_resolution_clock::now();
      st.lastMeshMs = std::chrono::duration<double, std::milli>(t1 - t0).count();

      st.hasMesh = !st.cpuMesh.vertices.empty() && !st.cpuMesh.indices.empty();
      if (st.hasMesh) {
        st.latestLodJobId = st.nextLodJobId++;
        resetLodChainToBaseMesh(st);
        uploadPreviewMeshForCurrentLod(st);
        if (st.buildLods && st.autoBuildLods && st.lodLevels > 1) {
          requestBuildLods(st);
        }
      }
    }
  }

  // LOD / Simplify controls (QEM-based edge collapse)
  {
    ImGui::SeparatorText("LOD / Simplify");

    ImGui::Checkbox("Build LOD chain", &st.buildLods);
    ImGui::SameLine();
    ImGui::Checkbox("Auto-build", &st.autoBuildLods);

    if (!st.buildLods) {
      // If the user disables LODs, clear any existing LOD1+ meshes.
      if (!st.lodMeshes.empty()) {
        st.lodMeshes.clear();
        st.lodStats.clear();
        st.lodSimplifyStats.clear();
        st.hasLods = false;
        st.previewLod = 0;
        st.exportLod = 0;
        st.exportAllLods = false;
        uploadPreviewMeshForCurrentLod(st);
      }
      ImGui::TextDisabled("LOD chain disabled (LOD0 only).");
    } else {
      int levels = st.lodLevels;
      if (ImGui::SliderInt("LOD levels", &levels, 1, 6)) {
        st.lodLevels = levels;
      }
      if (ImGui::IsItemDeactivatedAfterEdit() && st.autoBuildLods && st.hasMesh && st.lodLevels > 1) {
        requestBuildLods(st);
      }

      float ratio = st.lodRatioPerLevel;
      if (ImGui::SliderFloat("Ratio per level", &ratio, 0.05f, 0.95f, "%.2f")) {
        st.lodRatioPerLevel = ratio;
      }
      if (ImGui::IsItemDeactivatedAfterEdit() && st.autoBuildLods && st.hasMesh && st.lodLevels > 1) {
        requestBuildLods(st);
      }

      if (ImGui::Button("Build LODs now")) {
        requestBuildLods(st);
      }
      if (st.lodJobRunning) {
        ImGui::SameLine();
        ImGui::TextUnformatted("Computing…");
      }
      if (st.lastLodMs > 0.0) {
        ImGui::SameLine();
        ImGui::TextDisabled("%.2f ms", st.lastLodMs);
      }

      if (!st.lodError.empty()) {
        ImGui::TextColored(ImVec4(1, 0.65f, 0.25f, 1), "LOD: %s", st.lodError.c_str());
      }

      if (st.hasMesh) {
        const int maxLod = (int)st.lodMeshes.size(); // 0 = LOD0 only
        if (maxLod > 0) {
          int pl = st.previewLod;
          if (ImGui::SliderInt("Preview LOD", &pl, 0, maxLod)) {
            st.previewLod = pl;
            uploadPreviewMeshForCurrentLod(st);
          }

          int el = st.exportLod;
          if (ImGui::SliderInt("Export LOD", &el, 0, maxLod)) {
            st.exportLod = el;
          }

          ImGui::Checkbox("Export all LODs", &st.exportAllLods);
        } else {
          st.previewLod = 0;
          st.exportLod = 0;
          st.exportAllLods = false;
          ImGui::TextDisabled("(Only LOD0 available — set LOD levels > 1 and rebuild.)");
        }

        if (ImGui::BeginTable("lod_table", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
          ImGui::TableSetupColumn("LOD");
          ImGui::TableSetupColumn("Verts");
          ImGui::TableSetupColumn("Tris");
          ImGui::TableSetupColumn("Simplify ms");
          ImGui::TableSetupColumn("Max err");
          ImGui::TableHeadersRow();

          const int total = 1 + (int)st.lodMeshes.size();
          for (int lod = 0; lod < total; ++lod) {
            const std::size_t idx = (lod == 0) ? 0u : (std::size_t)(lod - 1);
            render::SdfMeshStats s = (lod == 0)
              ? st.cpuStats
              : (idx < st.lodStats.size() ? st.lodStats[idx] : render::SdfMeshStats{});
            render::MeshSimplifyStats ss = (lod == 0)
              ? render::MeshSimplifyStats{}
              : (idx < st.lodSimplifyStats.size() ? st.lodSimplifyStats[idx] : render::MeshSimplifyStats{});

            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            if (lod == st.previewLod) {
              ImGui::TextColored(ImVec4(0.35f, 0.9f, 0.35f, 1), "%d", lod);
            } else {
              ImGui::Text("%d", lod);
            }
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%zu", s.vertexCount);
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%zu", s.indexCount / 3);
            ImGui::TableSetColumnIndex(3);
            if (lod == 0) {
              ImGui::TextUnformatted("-");
            } else {
              ImGui::Text("%.2f", ss.ms);
            }
            ImGui::TableSetColumnIndex(4);
            if (lod == 0) {
              ImGui::TextUnformatted("-");
            } else {
              ImGui::Text("%.3g", ss.maxAcceptedError);
            }
          }

          ImGui::EndTable();
        }
      } else {
        ImGui::TextDisabled("(Generate a mesh to build LODs.)");
      }
    }
  }

  // Surface material: bake a 2D ProcGraph to an albedo texture and apply it to both
  // the raster preview and the raymarch preview.
  {
    ImGui::SeparatorText("Surface Material");

    ImGui::Checkbox("Use procedural albedo", &st.useProceduralTexture);
    ImGui::SameLine();
    ImGui::Checkbox("Apply to raymarch", &st.applyTextureToRaymarch);

    if (st.useProceduralTexture) {
      // Preset/seed controls (optionally locked).
      if (ImGui::BeginCombo("Albedo preset", render::procGraphPresetName(st.texPreset))) {
        for (int i = 0; i <= (int)render::ProcGraphPreset::AlienCircuit; ++i) {
          auto p = (render::ProcGraphPreset)i;
          const bool sel = (p == st.texPreset);
          if (ImGui::Selectable(render::procGraphPresetName(p), sel)) {
            st.texPreset = p;
          }
          if (sel) ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
      }

      ImGui::SameLine();
      ImGui::Checkbox("Lock##texLock", &st.texLockToPreset);

      core::u64 texSeed = st.texSeed;
      if (ImGui::InputScalar("Tex seed", ImGuiDataType_U64, &texSeed)) {
        st.texSeed = texSeed;
        st.texDirty = true;
      }

      ImGui::SameLine();
      if (ImGui::Button("Random##tex")) {
        st.texSeed = randomU64();
        st.texDirty = true;
      }

      // Apply preset when locked (and on first use).
      if (st.texLockToPreset && (st.texPreset != st.texAppliedPreset || st.texSeed != st.texAppliedSeed)) {
        st.texGraph = render::makeProceduralGraphPreset(st.texPreset, st.texSeed);
        st.texAppliedPreset = st.texPreset;
        st.texAppliedSeed = st.texSeed;
        st.texDirty = true;
      }

      if (ImGui::Checkbox("Use palette##texPal", &st.texGraph.usePalette)) {
        st.texDirty = true;
      }

      if (ImGui::SliderInt("Tex resolution", &st.texResolution, 64, 2048)) {
        st.texResolution = std::max(64, st.texResolution);
        st.texDirty = true;
      }

      if (ImGui::Checkbox("Generate mips##texMips", &st.texGenerateMips)) {
        st.texDirty = true;
      }
      if (ImGui::SliderFloat("Dither##texDither", &st.texDitherStrength, 0.0f, 2.0f, "%.2f")) {
        st.texDirty = true;
      }

      if (ImGui::Checkbox("Pack height in alpha##texHeightA", &st.texPackHeightInAlpha)) {
        // Requires a re-bake to take effect.
        st.texDirty = true;
      }

      // Apply quality options to the texture baker (affects the next bake).
      st.texBaker.setGenerateMips(st.texGenerateMips);
      st.texBaker.setDitherStrength(st.texDitherStrength);
      st.texBaker.setPackHeightInAlpha(st.texPackHeightInAlpha);

      ImGui::Checkbox("Auto-bake##tex", &st.texAutoBake);
      ImGui::SameLine();
      const bool wantBake = ImGui::Button("Bake texture") || (st.texAutoBake && st.texDirty);
      if (wantBake) {
        std::string err;
        if (!st.texBaker.bake(st.texGraph, st.texResolution, st.texResolution, timeSec, &err)) {
          st.texError = err;
          st.texDirty = false;
        } else {
          st.texError.clear();
          st.texDirty = false;
        }
      }

      if (!st.texError.empty()) {
        ImGui::TextColored(ImVec4(1, 0.3f, 0.3f, 1), "Texture error: %s", st.texError.c_str());
      }

      if (st.texBaker.isReady()) {
        ImGui::SameLine();
        const auto& s = st.texBaker.stats();
        if (st.texGenerateMips) {
          ImGui::TextDisabled("(shader %.2f ms, draw %.2f ms, mips %.2f ms)", s.shaderBuildMs, s.drawMs, s.mipsGenerated ? s.mipGenMs : 0.0);
        } else {
          ImGui::TextDisabled("(shader %.2f ms, draw %.2f ms)", s.shaderBuildMs, s.drawMs);
        }
        ImGui::Image((ImTextureID)(intptr_t)st.texBaker.texture().handle(), ImVec2(128, 128), ImVec2(0, 1), ImVec2(1, 0));
      }

      // Optional UV transform for the raymarch preview.
      ImGui::DragFloat("RM UV scale", &st.raymarchMat.uvScale, 0.01f, 0.05f, 20.0f);
      ImGui::DragFloat("RM UV rotate", &st.raymarchMat.uvRotateDeg, 0.5f, -360.0f, 360.0f);
      ImGui::DragFloat2("RM UV offset", st.raymarchMat.uvOffset, 0.01f, -10.0f, 10.0f);

      ImGui::SeparatorText("Raymarch Mapping");
      ImGui::Checkbox("RM Triplanar blend##rmTri", &st.raymarchMat.triplanarBlend);
      ImGui::SliderFloat("RM Tri sharpness##rmTriSharp", &st.raymarchMat.triplanarSharpness, 1.0f, 16.0f, "%.1f");

      ImGui::SeparatorText("Raymarch Micro Normals");
      ImGui::Checkbox("RM Height from alpha##rmHeightA", &st.raymarchMat.heightFromAlpha);
      ImGui::SliderFloat("RM Micro normal strength##rmMicroN", &st.raymarchMat.microNormalStrength, 0.0f, 4.0f, "%.2f");
      ImGui::SliderFloat("RM Micro normal step (texels)##rmMicroStep", &st.raymarchMat.microNormalStepTexels, 0.25f, 4.0f, "%.2f");

      if (st.raymarchMat.microNormalStrength > 0.0f && (!st.texPackHeightInAlpha || !st.raymarchMat.heightFromAlpha)) {
        ImGui::TextDisabled("Tip: enable 'Pack height in alpha' and 'RM Height from alpha' for best results.");
      }

      // Graph save/load (texture graph)
      ImGui::InputText("Tex graph file", st.texGraphPath, sizeof(st.texGraphPath));
      if (ImGui::Button("Save tex graph")) {
        std::string err;
        if (!render::saveProcGraphToFile(st.texGraph, st.texGraphPath, &err)) {
          toast(std::string("Save failed: ") + err, 4.0);
        } else {
          toast(std::string("Saved ") + st.texGraphPath, 3.0);
        }
      }
      ImGui::SameLine();
      if (ImGui::Button("Load tex graph")) {
        render::ProcGraph loaded;
        std::string err;
        if (!render::loadProcGraphFromFile(st.texGraphPath, loaded, &err)) {
          toast(std::string("Load failed: ") + err, 4.0);
        } else {
          st.texGraph = std::move(loaded);
          st.texSeed = st.texGraph.seed;
          st.texLockToPreset = false;
          st.texDirty = true;
          toast(std::string("Loaded ") + st.texGraphPath, 3.0);
        }
      }

      ImGui::Checkbox("Show tex GLSL", &st.showTexShader);
      if (st.showTexShader && st.texBaker.isReady()) {
        ImGui::BeginChild("##tex_glsl", ImVec2(0, 160), true);
        ImGui::TextUnformatted(st.texBaker.lastFragmentSource().c_str());
        ImGui::EndChild();
      }
    } else {
      ImGui::TextDisabled("(Surface texture disabled: using checkerboard)");
    }
  }

  ImGui::SeparatorText("Preview");
  {
    // Mode selector.
    int pMode = st.previewRaymarch ? 1 : 0;
    ImGui::RadioButton("Mesh", &pMode, 0);
    ImGui::SameLine();
    ImGui::RadioButton("Raymarch (SDF)", &pMode, 1);
    st.previewRaymarch = (pMode == 1);

    if (ImGui::SliderInt("Preview res", &st.previewResolution, 256, 1024)) {
      if (st.previewInited) {
        st.previewTarget.ensureSize(st.previewResolution, st.previewResolution);
      }
    }

    // Shared orbit controls.
    ImGui::SameLine();
    ImGui::Checkbox("Auto spin", &st.autoSpin);
    if (st.autoSpin) {
      ImGui::SameLine();
      ImGui::DragFloat("Spin deg/s", &st.spinDegPerSec, 0.1f, -200.0f, 200.0f);
    }

    if (!st.previewRaymarch) {
      ImGui::Checkbox("Wireframe", &st.wireframe);
    }

    ImGui::DragFloat("Yaw", &st.yawDeg, 0.5f, -180.0f, 180.0f);
    ImGui::DragFloat("Pitch", &st.pitchDeg, 0.5f, -89.0f, 89.0f);
    ImGui::DragFloat("Distance", &st.distance, 0.01f, 0.2f, 50.0f);

    if (st.previewRaymarch) {
      ImGui::SeparatorText("Raymarch");

      ImGui::Checkbox("Auto render", &st.raymarchAutoRender);
      ImGui::SameLine();
      ImGui::DragFloat("FOV", &st.raymarchFovDeg, 0.25f, 15.0f, 120.0f, "%.2f deg");

      bool wantRender = false;
      if (st.raymarchAutoRender || !st.raymarcher.isReady()) {
        wantRender = true;
      } else {
        if (ImGui::Button("Render")) {
          wantRender = true;
        }
      }

      ImGui::SameLine();
      ImGui::Checkbox("Bounds AABB", &st.raymarchSettings.useBoundsAabb);

      ImGui::DragInt("Max steps", &st.raymarchSettings.maxSteps, 1.0f, 16, 512);
      ImGui::DragFloat("Epsilon", &st.raymarchSettings.epsilon, 0.0001f, 0.00005f, 0.02f, "%.5f");
      ImGui::DragFloat("Max dist", &st.raymarchSettings.maxDistance, 0.05f, 1.0f, 200.0f, "%.2f");

      // Lighting controls.
      ImGui::ColorEdit3("Base color", st.raymarchSettings.baseColor);
      ImGui::ColorEdit3("Background", st.raymarchSettings.backgroundColor);
      ImGui::DragFloat3("Light dir", st.raymarchSettings.lightDir, 0.01f, -1.0f, 1.0f);
      ImGui::DragFloat("Ambient", &st.raymarchSettings.ambient, 0.01f, 0.0f, 2.0f);
      ImGui::DragFloat("Diffuse", &st.raymarchSettings.diffuse, 0.01f, 0.0f, 2.0f);
      ImGui::DragFloat("Specular", &st.raymarchSettings.specular, 0.01f, 0.0f, 2.0f);
      ImGui::DragFloat("Shininess", &st.raymarchSettings.shininess, 1.0f, 1.0f, 512.0f);

      ImGui::Checkbox("Soft shadows", &st.raymarchSettings.softShadows);
      if (st.raymarchSettings.softShadows) {
        ImGui::SameLine();
        ImGui::DragInt("Shadow steps", &st.raymarchSettings.shadowSteps, 1.0f, 0, 64);
        ImGui::DragFloat("Shadow max", &st.raymarchSettings.shadowMaxDistance, 0.05f, 0.0f, 50.0f);
        ImGui::DragFloat("Shadow k", &st.raymarchSettings.shadowK, 0.25f, 0.0f, 64.0f);
      }

      ImGui::Checkbox("AO", &st.raymarchSettings.ambientOcclusion);
      if (st.raymarchSettings.ambientOcclusion) {
        ImGui::SameLine();
        ImGui::DragInt("AO steps", &st.raymarchSettings.aoSteps, 1.0f, 0, 16);
        ImGui::DragFloat("AO step", &st.raymarchSettings.aoStepSize, 0.01f, 0.0f, 1.0f);
        ImGui::DragFloat("AO strength", &st.raymarchSettings.aoStrength, 0.05f, 0.0f, 4.0f);
      }

      const char* debugNames[] = {"Shaded", "Normals", "Steps", "Distance"};
      int dbg = (int)st.raymarchSettings.debug;
      if (ImGui::Combo("Debug", &dbg, debugNames, IM_ARRAYSIZE(debugNames))) {
        st.raymarchSettings.debug = (render::SdfRaymarchDebug)dbg;
      }

      ImGui::Checkbox("Show shader", &st.showRaymarchShader);

      if (wantRender) {
        renderRaymarchPreview(st, timeSec);
      }

      if (!st.raymarchError.empty()) {
        ImGui::TextColored(ImVec4(1, 0.25f, 0.25f, 1), "Raymarch error: %s", st.raymarchError.c_str());
      }

      if (st.raymarcher.isReady()) {
        const ImVec2 size((float)st.raymarcher.texture().width(), (float)st.raymarcher.texture().height());
        ImGui::Image((ImTextureID)(intptr_t)st.raymarcher.texture().handle(), size, ImVec2(0, 1), ImVec2(1, 0));

        const auto& stats = st.raymarcher.stats();
        ImGui::Text("GPU: %.2f ms  Shader: %s (%.2f ms)",
                    stats.drawMs,
                    stats.shaderRebuilt ? "rebuilt" : "cached",
                    stats.shaderBuildMs);
      } else {
        ImGui::TextUnformatted("(raymarch preview not ready)");
      }

      if (st.showRaymarchShader) {
        ImGui::BeginChild("##raymarch_shader", ImVec2(0, 220), true);
        ImGui::TextUnformatted(st.raymarcher.lastFragmentSource().c_str());
        ImGui::EndChild();
      }
    } else {
      // Mesh preview.
      if (st.previewInited && st.hasMesh) {
        renderPreview(st, timeSec);

        const ImVec2 size((float)st.previewTarget.width(), (float)st.previewTarget.height());
        ImGui::Image((ImTextureID)(intptr_t)st.previewTarget.color().handle(), size, ImVec2(0, 1), ImVec2(1, 0));
      } else {
        if (!st.previewInitError.empty()) {
          ImGui::TextColored(ImVec4(1, 0.25f, 0.25f, 1), "Preview init failed: %s", st.previewInitError.c_str());
        } else {
          ImGui::TextUnformatted("(mesh preview not ready — generate a mesh first)");
        }
      }
    }
  }

  ImGui::SeparatorText("Graph");
  {
    // Output selection
    if (beginComboOutputIndex("Output", st.graph.output, st.graph)) {
      st.dirty = true;
    }

    // Quick add buttons
    if (ImGui::Button("+ Sphere")) { addNode(st.graph, render::SdfNodeOp::Sphere); st.dirty = true; }
    ImGui::SameLine();
    if (ImGui::Button("+ Box")) { addNode(st.graph, render::SdfNodeOp::Box); st.dirty = true; }
    ImGui::SameLine();
    if (ImGui::Button("+ Capsule")) { addNode(st.graph, render::SdfNodeOp::Capsule); st.dirty = true; }
    ImGui::SameLine();
    if (ImGui::Button("+ Torus")) { addNode(st.graph, render::SdfNodeOp::TorusY); st.dirty = true; }
    ImGui::SameLine();
    if (ImGui::Button("+ Noise")) { addNode(st.graph, render::SdfNodeOp::NoiseDisplace); st.dirty = true; }
    ImGui::SameLine();
    if (ImGui::Button("+ Perlin Noise")) { addNode(st.graph, render::SdfNodeOp::NoiseDisplacePerlin); st.dirty = true; }

    ImGui::SameLine();
    if (ImGui::Button("Clear")) {
      st.graph = render::SdfGraph::makeDefault();
      st.dirty = true;
    }

    ImGui::TextDisabled("Nodes: %zu / %d", st.graph.nodes.size(), render::kSdfGraphMaxNodes);

    ImGui::Separator();
    ImGui::TextDisabled("Domain transforms (wrap an input node):");
    if (ImGui::Button("+ Translate")) { addNode(st.graph, render::SdfNodeOp::Translate); st.dirty = true; }
    ImGui::SameLine();
    if (ImGui::Button("+ RotateY")) { addNode(st.graph, render::SdfNodeOp::RotateY); st.dirty = true; }
    ImGui::SameLine();
    if (ImGui::Button("+ Scale")) { addNode(st.graph, render::SdfNodeOp::Scale); st.dirty = true; }
    ImGui::SameLine();
    if (ImGui::Button("+ Repeat")) { addNode(st.graph, render::SdfNodeOp::Repeat); st.dirty = true; }
    ImGui::SameLine();
    if (ImGui::Button("+ Mirror")) { addNode(st.graph, render::SdfNodeOp::Mirror); st.dirty = true; }
    ImGui::SameLine();
    if (ImGui::Button("+ TwistY")) { addNode(st.graph, render::SdfNodeOp::TwistY); st.dirty = true; }

    // Node list
    ImGui::BeginChild("##sdf_nodes", ImVec2(0, 260), true);
    for (int i = 0; i < (int)st.graph.nodes.size(); ++i) {
      bool nodeChanged = drawNodeEditor(st.graph, i);
      if (nodeChanged) {
        st.dirty = true;
      }

      ImGui::PushID(i);
      if (ImGui::Button("Remove")) {
        removeNode(st.graph, i);
        st.dirty = true;
        ImGui::PopID();
        break; // list invalidated
      }
      ImGui::Separator();
      ImGui::PopID();
    }
    ImGui::EndChild();
  }

  ImGui::SeparatorText("Export");
  {
    ImGui::InputText("OBJ path", st.exportPath, sizeof(st.exportPath));
    ImGui::InputText("glTF path", st.exportGltfPath, sizeof(st.exportGltfPath));

    ImGui::Checkbox("Include material + albedo PNG", &st.exportWithMaterial);
    ImGui::SameLine();
    ImGui::Checkbox("Flip PNG Y", &st.exportTextureFlipY);

    ImGui::Checkbox("Include tangents (glTF TANGENT)", &st.exportGltfTangents);
    ImGui::SameLine();
    {
      ImGui::BeginDisabled(!st.exportAllLods);
      ImGui::Checkbox("Pack LODs into one glTF (MSFT_lod)", &st.exportGltfPackLodsMsft);
      ImGui::EndDisabled();
      if (!st.exportAllLods) {
        ImGui::TextDisabled("(Enable 'Export all LODs' to use MSFT_lod packaging)");
      }
    }

    ImGui::SliderFloat("PBR metallic", &st.exportPbrMetallic, 0.0f, 1.0f);
    ImGui::SliderFloat("PBR roughness", &st.exportPbrRoughness, 0.0f, 1.0f);

    ImGui::TextDisabled("glTF export writes .gltf + .bin (and optional PNG). glTF UV origin is upper-left; exporter flips V by default.\nOptionally exports TANGENT and can package LODs via MSFT_lod.");

    const bool clickedObj = ImGui::Button("Export OBJ");
    ImGui::SameLine();
    const bool clickedGltf = ImGui::Button("Export glTF");

    if (clickedObj) {
      namespace fs = std::filesystem;
      std::string err;

      if (!st.hasMesh) {
        st.lastError = "No mesh to export (generate a mesh first).";
        toast("Export failed", 2.0);
        ImGui::End();
        return;
      }

      // Choose export source.
      const int maxLod = (int)st.lodMeshes.size(); // 0 = LOD0 only
      st.exportLod = std::clamp(st.exportLod, 0, std::max(0, maxLod));
      const render::SdfMeshData* exportMesh = &st.cpuMesh; // LOD0
      if (st.exportLod > 0) {
        const std::size_t idx = (std::size_t)(st.exportLod - 1);
        if (idx < st.lodMeshes.size()) {
          exportMesh = &st.lodMeshes[idx];
        }
      }

      const bool canTex = st.useProceduralTexture && st.texBaker.isReady();
      const bool withMat = st.exportWithMaterial && canTex;

      if (withMat) {
        const fs::path objPath(st.exportPath);
        fs::path mtlPath = objPath;
        mtlPath.replace_extension(".mtl");

        const std::string baseName = objPath.stem().string();
        const fs::path texPath = objPath.parent_path() / (baseName + "_albedo.png");

        const std::string matName = "mat0";

        // 1) Export PNG
        std::vector<unsigned char> pixels;
        std::string texErr;
        const auto& tex = st.texBaker.texture();
        const int w = tex.width();
        const int h = tex.height();
        if (!readTextureRgba8(tex, w, h, pixels, &texErr)) {
          err = texErr;
        } else {
          const int stride = w * 4;
          if (!writePixelsToPng(texPath.string(), w, h, 4, pixels.data(), stride, st.exportTextureFlipY, &texErr)) {
            err = texErr;
          }
        }

        // 2) Write MTL
        if (err.empty()) {
          if (!writeMtl(mtlPath.string(), matName, texPath.filename().string(), &texErr)) {
            err = texErr;
          }
        }

        // 3) Write OBJ(s) that reference the MTL
        if (err.empty()) {
          if (st.exportAllLods && maxLod > 0) {
            // Batch export: foo.obj -> foo_lod0.obj, foo_lod1.obj, ...
            for (int li = 0; li <= maxLod && err.empty(); ++li) {
              const fs::path objPath(st.exportPath);
              const std::string ext = objPath.extension().string();
              const std::string baseName = objPath.stem().string();
              const fs::path outObj = objPath.parent_path() / (baseName + "_lod" + std::to_string(li) + ext);

              const render::SdfMeshData& m = (li == 0)
                ? st.cpuMesh
                : st.lodMeshes[(std::size_t)(li - 1)];
              if (!writeObjWithMtl(outObj.string(), m, mtlPath.filename().string(), matName, &texErr)) {
                err = texErr;
              }
            }
          } else {
            if (!writeObjWithMtl(st.exportPath, *exportMesh, mtlPath.filename().string(), matName, &texErr)) {
              err = texErr;
            }
          }
        }

        if (!err.empty()) {
          st.lastError = err;
          toast("Export failed", 3.0);
        } else {
          toast(std::string("Exported ") + st.exportPath + " (+ MTL/PNG)", 3.0);
        }
      } else {
        if (st.exportAllLods && maxLod > 0) {
          // Batch export without MTL.
          for (int li = 0; li <= maxLod && err.empty(); ++li) {
            const fs::path objPath(st.exportPath);
            const std::string ext = objPath.extension().string();
            const std::string baseName = objPath.stem().string();
            const fs::path outObj = objPath.parent_path() / (baseName + "_lod" + std::to_string(li) + ext);

            const render::SdfMeshData& m = (li == 0)
              ? st.cpuMesh
              : st.lodMeshes[(std::size_t)(li - 1)];
            if (!writeObj(outObj.string(), m, &err)) {
              // err already set
              break;
            }
          }
        } else {
          if (!writeObj(st.exportPath, *exportMesh, &err)) {
            // err already set
          }
        }

        if (!err.empty()) {
          st.lastError = err;
          toast("OBJ export failed", 2.5);
        } else {
          if (st.exportWithMaterial && !canTex) {
            toast(std::string("Exported ") + st.exportPath + " (no baked texture for MTL)", 3.0);
          } else if (st.exportAllLods && maxLod > 0) {
            toast("Exported LOD chain", 3.0);
          } else {
            toast(std::string("Exported ") + st.exportPath, 2.5);
          }
        }
      }
    }

    if (clickedGltf) {
      namespace fs = std::filesystem;
      std::string err;

      if (!st.hasMesh) {
        st.lastError = "No mesh to export (generate a mesh first).";
        toast("Export failed", 2.0);
        ImGui::End();
        return;
      }

      // Choose export source.
      const int maxLod = (int)st.lodMeshes.size(); // 0 = LOD0 only
      st.exportLod = std::clamp(st.exportLod, 0, std::max(0, maxLod));
      const render::SdfMeshData* exportMesh = &st.cpuMesh; // LOD0
      if (st.exportLod > 0) {
        const std::size_t idx = (std::size_t)(st.exportLod - 1);
        if (idx < st.lodMeshes.size()) {
          exportMesh = &st.lodMeshes[idx];
        }
      }

      const bool canTex = st.useProceduralTexture && st.texBaker.isReady();
      const bool withTex = st.exportWithMaterial && canTex;

      const fs::path rootPath(st.exportGltfPath);
      const fs::path rootDir = rootPath.parent_path();

      // Texture path derived from the *root* glTF path (so LOD batch exports can share a texture).
      const std::string rootBaseName = rootPath.stem().string();
      const fs::path sharedTexPath = (rootDir.empty() ? fs::path(".") : rootDir) / (rootBaseName + "_albedo.png");
      const std::string sharedTexUri = sharedTexPath.filename().string();

      // 1) Export texture (optional)
      if (withTex) {
        std::error_code ec;
        const fs::path texDir = sharedTexPath.parent_path();
        if (!texDir.empty() && !fs::exists(texDir, ec)) {
          fs::create_directories(texDir, ec);
        }

        std::vector<unsigned char> pixels;
        std::string texErr;
        const auto& tex = st.texBaker.texture();
        const int w = tex.width();
        const int h = tex.height();
        if (!readTextureRgba8(tex, w, h, pixels, &texErr)) {
          err = texErr;
        } else {
          const int stride = w * 4;
          if (!writePixelsToPng(sharedTexPath.string(), w, h, 4, pixels.data(), stride, st.exportTextureFlipY, &texErr)) {
            err = texErr;
          }
        }
      }

      auto buildOpt = [&]() {
        render::GltfExportOptions opt;
        opt.meshName = "procedural_mesh";
        opt.materialName = "mat0";
        opt.nodeName = "node0";
        opt.doubleSided = true;
        opt.baseColorTextureUri = withTex ? sharedTexUri : std::string{};
        opt.metallicFactor = std::clamp(st.exportPbrMetallic, 0.0f, 1.0f);
        opt.roughnessFactor = std::clamp(st.exportPbrRoughness, 0.0f, 1.0f);
        opt.exportTangents = st.exportGltfTangents;

        // Our engine's UVs are authored in OpenGL-style (origin bottom-left). glTF's UV origin is upper-left.
        // Flip V so the asset looks correct in standard glTF viewers.
        opt.flipTexcoordV = true;
        return opt;
      };

      auto exportOne = [&](const fs::path& outGltf, const render::SdfMeshData& m) {
        const render::GltfExportOptions opt = buildOpt();
        std::string gltfErr;
        if (!render::exportMeshToGltf(outGltf.string(), m, opt, &gltfErr)) {
          err = gltfErr;
        }
      };

      if (err.empty()) {
        if (st.exportAllLods && maxLod > 0) {
          if (st.exportGltfPackLodsMsft) {
            const render::GltfExportOptions opt = buildOpt();
            std::string gltfErr;
            const std::span<const render::SdfMeshData> extra(st.lodMeshes.data(), st.lodMeshes.size());
            if (!render::exportMeshLodsToGltf(rootPath.string(), st.cpuMesh, extra, opt, &gltfErr)) {
              err = gltfErr;
            }
          } else {
            // Batch export: foo.gltf -> foo_lod0.gltf, foo_lod1.gltf, ...
            for (int li = 0; li <= maxLod && err.empty(); ++li) {
              const render::SdfMeshData& m = (li == 0)
                ? st.cpuMesh
                : st.lodMeshes[(std::size_t)(li - 1)];

              fs::path out = rootPath;
              const std::string ext = out.extension().string();
              const std::string base = out.stem().string();
              out = out.parent_path() / (base + "_lod" + std::to_string(li) + ext);
              exportOne(out, m);
            }
          }
        } else {
          exportOne(rootPath, *exportMesh);
        }
      }

      if (!err.empty()) {
        st.lastError = err;
        toast("glTF export failed", 3.0);
      } else {
        if (withTex) {
          toast(std::string("Exported ") + st.exportGltfPath + " (+ BIN/PNG)", 3.0);
        } else {
          toast(std::string("Exported ") + st.exportGltfPath + " (+ BIN)", 3.0);
        }
      }
    }
  }


  // Auto-commit history changes (undo/redo snapshotting).
  historyAutoCommit(st, timeSec);

  ImGui::End();
}

} // namespace stellar::game
