#include "ProceduralMeshLabWindow.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Log.h"
#include "stellar/math/Math.h"
#include "stellar/math/Mat4.h"
#include "stellar/math/Vec3.h"
#include "stellar/render/Gl.h"
#include "stellar/render/GltfExport.h"

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

namespace stellar::game {

namespace {

using namespace stellar;

// Forward declarations for helpers used before their definitions.
static void uploadPreviewMeshForCurrentLod(ProceduralMeshLabWindowState& st);
static void resetLodChainToBaseMesh(ProceduralMeshLabWindowState& st);
static void requestBuildLods(ProceduralMeshLabWindowState& st);

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

  // Restore core GL object bindings.
  // NOTE: We intentionally use the engine's tiny GL loader wrappers here
  // instead of calling glUseProgram/glBindVertexArray/glBindFramebuffer
  // directly (those entry points are not exported from legacy Windows gl.h).
  if (render::gl::UseProgram) render::gl::UseProgram((GLuint)prevProg);
  if (render::gl::BindVertexArray) render::gl::BindVertexArray((GLuint)prevVao);
  if (render::gl::BindFramebuffer) render::gl::BindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
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
        ImGui::TextDisabled("(shader %.2f ms, draw %.2f ms)", st.texBaker.stats().shaderBuildMs, st.texBaker.stats().drawMs);
        ImGui::Image((ImTextureID)(intptr_t)st.texBaker.texture().handle(), ImVec2(128, 128), ImVec2(0, 1), ImVec2(1, 0));
      }

      // Optional UV transform for the raymarch preview.
      ImGui::DragFloat("RM UV scale", &st.raymarchMat.uvScale, 0.01f, 0.05f, 20.0f);
      ImGui::DragFloat("RM UV rotate", &st.raymarchMat.uvRotateDeg, 0.5f, -360.0f, 360.0f);
      ImGui::DragFloat2("RM UV offset", st.raymarchMat.uvOffset, 0.01f, -10.0f, 10.0f);

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

  ImGui::End();
}

} // namespace stellar::game
