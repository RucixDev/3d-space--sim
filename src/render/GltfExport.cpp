#include "stellar/render/GltfExport.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace stellar::render {

namespace fs = std::filesystem;

namespace {

static void setErr(std::string* outError, const std::string& msg) {
  if (outError) *outError = msg;
}

static void appendBytes(std::vector<std::uint8_t>& dst, const void* src, std::size_t n) {
  const auto* p = static_cast<const std::uint8_t*>(src);
  dst.insert(dst.end(), p, p + n);
}

template <class T>
static void appendPod(std::vector<std::uint8_t>& dst, const T& v) {
  static_assert(std::is_trivially_copyable_v<T>, "appendPod requires trivially copyable type");
  appendBytes(dst, &v, sizeof(T));
}

static void padTo(std::vector<std::uint8_t>& dst, std::size_t align) {
  if (align == 0) return;
  const std::size_t rem = dst.size() % align;
  if (rem == 0) return;
  const std::size_t pad = align - rem;
  dst.insert(dst.end(), pad, 0u);
}

static std::string jsonEscape(const std::string& s) {
  std::string out;
  out.reserve(s.size() + 8);
  for (unsigned char ch : s) {
    switch (ch) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\b': out += "\\b"; break;
      case '\f': out += "\\f"; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      default:
        if (ch < 0x20) {
          // Control characters must be escaped.
          char buf[8];
          std::snprintf(buf, sizeof(buf), "\\u%04x", (unsigned int)ch);
          out += buf;
        } else {
          out.push_back((char)ch);
        }
        break;
    }
  }
  return out;
}

struct ViewInfo {
  std::size_t offset{0};
  std::size_t length{0};
  int target{0};
};

// Simple internal vector math (double precision for stable tangent solves).
struct V2 {
  double x{0}, y{0};
};

struct V3 {
  double x{0}, y{0}, z{0};
};

static inline V3 operator+(const V3& a, const V3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
static inline V3 operator-(const V3& a, const V3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
static inline V3 operator*(const V3& a, double s) { return {a.x * s, a.y * s, a.z * s}; }

static inline double dot(const V3& a, const V3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

static inline V3 cross(const V3& a, const V3& b) {
  return {a.y * b.z - a.z * b.y,
          a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x};
}

static inline double lengthSq(const V3& v) { return dot(v, v); }

static inline V3 normalize(const V3& v, double eps = 1e-18) {
  const double l2 = lengthSq(v);
  if (l2 <= eps) return {0, 0, 0};
  const double invL = 1.0 / std::sqrt(l2);
  return v * invL;
}

static inline V3 vertexPos(const VertexPNUT& v) { return {v.px, v.py, v.pz}; }
static inline V3 vertexNrm(const VertexPNUT& v) { return {v.nx, v.ny, v.nz}; }

static inline V2 vertexUv(const VertexPNUT& v, bool flipV) {
  const double uu = (double)v.u;
  const double vv = flipV ? (1.0 - (double)v.v) : (double)v.v;
  return {uu, vv};
}

static void buildTangents4f(const SdfMeshData& mesh, bool flipV, std::vector<float>& outTangents) {
  const std::size_t n = mesh.vertices.size();
  outTangents.assign(n * 4, 0.0f);
  if (n == 0) return;

  std::vector<V3> tan1(n, {0, 0, 0});
  std::vector<V3> tan2(n, {0, 0, 0});

  const auto& idx = mesh.indices;
  for (std::size_t t = 0; t + 2 < idx.size(); t += 3) {
    const std::uint32_t i0 = idx[t + 0];
    const std::uint32_t i1 = idx[t + 1];
    const std::uint32_t i2 = idx[t + 2];
    if (i0 >= n || i1 >= n || i2 >= n) continue;

    const V3 p0 = vertexPos(mesh.vertices[i0]);
    const V3 p1 = vertexPos(mesh.vertices[i1]);
    const V3 p2 = vertexPos(mesh.vertices[i2]);

    const V2 w0 = vertexUv(mesh.vertices[i0], flipV);
    const V2 w1 = vertexUv(mesh.vertices[i1], flipV);
    const V2 w2 = vertexUv(mesh.vertices[i2], flipV);

    const V3 e1 = p1 - p0;
    const V3 e2 = p2 - p0;

    const double x1 = w1.x - w0.x;
    const double x2 = w2.x - w0.x;
    const double y1 = w1.y - w0.y;
    const double y2 = w2.y - w0.y;

    const double det = x1 * y2 - y1 * x2;
    // Skip degenerate UV triangles.
    if (std::fabs(det) < 1e-20) continue;

    const double r = 1.0 / det;
    const V3 sdir = (e1 * y2 - e2 * y1) * r;
    const V3 tdir = (e2 * x1 - e1 * x2) * r;

    tan1[i0] = tan1[i0] + sdir;
    tan1[i1] = tan1[i1] + sdir;
    tan1[i2] = tan1[i2] + sdir;

    tan2[i0] = tan2[i0] + tdir;
    tan2[i1] = tan2[i1] + tdir;
    tan2[i2] = tan2[i2] + tdir;
  }

  for (std::size_t i = 0; i < n; ++i) {
    V3 nrm = normalize(vertexNrm(mesh.vertices[i]), 1e-18);
    if (lengthSq(nrm) < 1e-18) nrm = {0, 0, 1};

    V3 t = tan1[i];
    // Gram-Schmidt orthogonalize.
    t = t - nrm * dot(nrm, t);

    double w = 1.0;
    if (lengthSq(t) < 1e-20) {
      // UVs are degenerate or missing; fall back to a stable basis from the normal.
      const V3 axis = (std::fabs(nrm.z) < 0.9) ? V3{0, 0, 1} : V3{0, 1, 0};
      t = normalize(cross(axis, nrm), 1e-18);
      w = 1.0;
    } else {
      t = normalize(t, 1e-18);
      const V3 b = tan2[i];
      // Handedness: if (n x t) points in opposite direction of b, flip.
      w = (dot(cross(nrm, t), b) < 0.0) ? -1.0 : 1.0;
    }

    outTangents[i * 4 + 0] = (float)t.x;
    outTangents[i * 4 + 1] = (float)t.y;
    outTangents[i * 4 + 2] = (float)t.z;
    outTangents[i * 4 + 3] = (float)w;
  }
}

struct AccessorInfo {
  int bufferView{0};
  std::size_t byteOffset{0};
  int componentType{5126};
  std::size_t count{0};
  const char* type{"SCALAR"};

  bool hasMinMax{false};
  int minMaxCount{0};
  float minv[4]{};
  float maxv[4]{};
};

struct PackedPrimitiveAccessors {
  int accPos{-1};
  int accNrm{-1};
  int accUv{-1};
  int accTan{-1};
  int accIdx{-1};

  bool hasTangents() const { return accTan >= 0; }
};

static bool packMeshToBuffer(const SdfMeshData& mesh,
                             const GltfExportOptions& opt,
                             std::vector<std::uint8_t>& bin,
                             std::vector<ViewInfo>& views,
                             std::vector<AccessorInfo>& accessors,
                             PackedPrimitiveAccessors& outPrim,
                             std::string* outError) {
  if (mesh.vertices.empty() || mesh.indices.empty()) {
    setErr(outError, "exportMeshToGltf: mesh is empty.");
    return false;
  }
  if (mesh.indices.size() % 3 != 0) {
    setErr(outError, "exportMeshToGltf: index count is not a multiple of 3.");
    return false;
  }

  const std::size_t vCount = mesh.vertices.size();
  const std::size_t iCount = mesh.indices.size();

  // Validate indices and decide u16/u32.
  std::uint32_t maxIndex = 0;
  for (std::uint32_t i : mesh.indices) {
    if (i >= (std::uint32_t)vCount) {
      setErr(outError, "exportMeshToGltf: index out of range for vertex buffer.");
      return false;
    }
    maxIndex = std::max(maxIndex, i);
  }

  const bool useU16 = (vCount <= 65535u) && (maxIndex <= 65535u);

  // Compute POSITION min/max for accessor.
  float minX = mesh.vertices[0].px;
  float minY = mesh.vertices[0].py;
  float minZ = mesh.vertices[0].pz;
  float maxX = mesh.vertices[0].px;
  float maxY = mesh.vertices[0].py;
  float maxZ = mesh.vertices[0].pz;
  for (const auto& v : mesh.vertices) {
    minX = std::min(minX, v.px);
    minY = std::min(minY, v.py);
    minZ = std::min(minZ, v.pz);
    maxX = std::max(maxX, v.px);
    maxY = std::max(maxY, v.py);
    maxZ = std::max(maxZ, v.pz);
  }

  std::vector<float> tangents;
  if (opt.exportTangents) {
    buildTangents4f(mesh, opt.flipTexcoordV, tangents);
  }

  // Positions (VEC3 float)
  {
    ViewInfo v{};
    padTo(bin, 4);
    v.offset = bin.size();
    for (const auto& vert : mesh.vertices) {
      appendPod(bin, vert.px);
      appendPod(bin, vert.py);
      appendPod(bin, vert.pz);
    }
    padTo(bin, 4);
    v.length = bin.size() - v.offset;
    v.target = 34962; // ARRAY_BUFFER

    const int viewIdx = (int)views.size();
    views.push_back(v);

    AccessorInfo a{};
    a.bufferView = viewIdx;
    a.byteOffset = 0;
    a.componentType = 5126; // FLOAT
    a.count = vCount;
    a.type = "VEC3";
    a.hasMinMax = true;
    a.minMaxCount = 3;
    a.minv[0] = minX;
    a.minv[1] = minY;
    a.minv[2] = minZ;
    a.maxv[0] = maxX;
    a.maxv[1] = maxY;
    a.maxv[2] = maxZ;

    outPrim.accPos = (int)accessors.size();
    accessors.push_back(a);
  }

  // Normals (VEC3 float)
  {
    ViewInfo v{};
    padTo(bin, 4);
    v.offset = bin.size();
    for (const auto& vert : mesh.vertices) {
      appendPod(bin, vert.nx);
      appendPod(bin, vert.ny);
      appendPod(bin, vert.nz);
    }
    padTo(bin, 4);
    v.length = bin.size() - v.offset;
    v.target = 34962;

    const int viewIdx = (int)views.size();
    views.push_back(v);

    AccessorInfo a{};
    a.bufferView = viewIdx;
    a.byteOffset = 0;
    a.componentType = 5126;
    a.count = vCount;
    a.type = "VEC3";

    outPrim.accNrm = (int)accessors.size();
    accessors.push_back(a);
  }

  // UVs (VEC2 float)
  {
    ViewInfo v{};
    padTo(bin, 4);
    v.offset = bin.size();
    for (const auto& vert : mesh.vertices) {
      const float uu = vert.u;
      const float vv = opt.flipTexcoordV ? (1.0f - vert.v) : vert.v;
      appendPod(bin, uu);
      appendPod(bin, vv);
    }
    padTo(bin, 4);
    v.length = bin.size() - v.offset;
    v.target = 34962;

    const int viewIdx = (int)views.size();
    views.push_back(v);

    AccessorInfo a{};
    a.bufferView = viewIdx;
    a.byteOffset = 0;
    a.componentType = 5126;
    a.count = vCount;
    a.type = "VEC2";

    outPrim.accUv = (int)accessors.size();
    accessors.push_back(a);
  }

  // Tangents (VEC4 float)
  if (opt.exportTangents) {
    ViewInfo v{};
    padTo(bin, 4);
    v.offset = bin.size();

    // tangents.size() == vCount*4
    for (float f : tangents) {
      appendPod(bin, f);
    }

    padTo(bin, 4);
    v.length = bin.size() - v.offset;
    v.target = 34962;

    const int viewIdx = (int)views.size();
    views.push_back(v);

    AccessorInfo a{};
    a.bufferView = viewIdx;
    a.byteOffset = 0;
    a.componentType = 5126;
    a.count = vCount;
    a.type = "VEC4";

    outPrim.accTan = (int)accessors.size();
    accessors.push_back(a);
  }

  // Indices (SCALAR u16/u32)
  {
    ViewInfo v{};
    padTo(bin, 4);
    v.offset = bin.size();

    if (useU16) {
      for (std::uint32_t i : mesh.indices) {
        const std::uint16_t ii = (std::uint16_t)i;
        appendPod(bin, ii);
      }
    } else {
      for (std::uint32_t i : mesh.indices) {
        appendPod(bin, i);
      }
    }

    padTo(bin, 4);
    v.length = bin.size() - v.offset;
    v.target = 34963; // ELEMENT_ARRAY_BUFFER

    const int viewIdx = (int)views.size();
    views.push_back(v);

    AccessorInfo a{};
    a.bufferView = viewIdx;
    a.byteOffset = 0;
    a.componentType = useU16 ? 5123 : 5125; // UNSIGNED_SHORT / UNSIGNED_INT
    a.count = iCount;
    a.type = "SCALAR";

    outPrim.accIdx = (int)accessors.size();
    accessors.push_back(a);
  }

  return true;
}

static bool ensureOutDir(const fs::path& outDir, std::string* outError) {
  std::error_code ec;
  if (!fs::exists(outDir, ec)) {
    fs::create_directories(outDir, ec);
  }
  if (ec) {
    setErr(outError, "exportMeshToGltf: failed to create directory: " + outDir.string());
    return false;
  }
  return true;
}

static bool writeBinaryFile(const fs::path& path, const std::vector<std::uint8_t>& bin, std::string* outError) {
  std::ofstream f(path, std::ios::binary | std::ios::out | std::ios::trunc);
  if (!f) {
    setErr(outError, "exportMeshToGltf: failed to open for writing: " + path.string());
    return false;
  }
  f.write(reinterpret_cast<const char*>(bin.data()), (std::streamsize)bin.size());
  if (!f) {
    setErr(outError, "exportMeshToGltf: failed to write: " + path.string());
    return false;
  }
  return true;
}

static bool writeTextFile(const fs::path& path, const std::string& text, std::string* outError) {
  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f) {
    setErr(outError, "exportMeshToGltf: failed to open for writing: " + path.string());
    return false;
  }
  f.write(text.c_str(), (std::streamsize)text.size());
  if (!f) {
    setErr(outError, "exportMeshToGltf: failed to write: " + path.string());
    return false;
  }
  return true;
}

static void emitBufferViews(std::ostringstream& js, const std::vector<ViewInfo>& views) {
  js << "  \"bufferViews\": [\n";
  for (std::size_t i = 0; i < views.size(); ++i) {
    const ViewInfo& v = views[i];
    js << "    { \"buffer\": 0, \"byteOffset\": " << v.offset
       << ", \"byteLength\": " << v.length;
    if (v.target != 0) {
      js << ", \"target\": " << v.target;
    }
    js << " }";
    js << ((i + 1 == views.size()) ? "\n" : ",\n");
  }
  js << "  ],\n";
}

static void emitAccessors(std::ostringstream& js, const std::vector<AccessorInfo>& accessors) {
  js << "  \"accessors\": [\n";
  for (std::size_t i = 0; i < accessors.size(); ++i) {
    const AccessorInfo& a = accessors[i];
    js << "    { \"bufferView\": " << a.bufferView
       << ", \"byteOffset\": " << a.byteOffset
       << ", \"componentType\": " << a.componentType
       << ", \"count\": " << a.count
       << ", \"type\": \"" << a.type << "\"";

    if (a.hasMinMax && a.minMaxCount > 0) {
      js << ", \"min\": [ ";
      for (int k = 0; k < a.minMaxCount; ++k) {
        js << a.minv[k];
        if (k + 1 < a.minMaxCount) js << ", ";
      }
      js << " ], \"max\": [ ";
      for (int k = 0; k < a.minMaxCount; ++k) {
        js << a.maxv[k];
        if (k + 1 < a.minMaxCount) js << ", ";
      }
      js << " ]";
    }

    js << " }";
    js << ((i + 1 == accessors.size()) ? "\n" : ",\n");
  }
  js << "  ]\n";
}

} // namespace

bool exportMeshToGltf(const std::string& gltfPath,
                      const SdfMeshData& mesh,
                      const GltfExportOptions& opt,
                      std::string* outError) {
  if (gltfPath.empty()) {
    setErr(outError, "exportMeshToGltf: empty glTF path.");
    return false;
  }

  fs::path outPath(gltfPath);
  fs::path outDir = outPath.parent_path();
  if (outDir.empty()) outDir = fs::path(".");

  if (!ensureOutDir(outDir, outError)) return false;

  fs::path binPath = outPath;
  binPath.replace_extension(".bin");
  const std::string binUri = binPath.filename().string();

  // ---- Pack binary buffer ----
  // We keep attributes in separate contiguous blocks for simplicity.
  // Ensure 4-byte alignment between blocks.
  std::vector<std::uint8_t> bin;
  bin.reserve(mesh.vertices.size() * (3 + 3 + 2 + (opt.exportTangents ? 4 : 0)) * sizeof(float)
              + mesh.indices.size() * sizeof(std::uint32_t) + 64);

  std::vector<ViewInfo> views;
  std::vector<AccessorInfo> accessors;
  PackedPrimitiveAccessors prim{};

  if (!packMeshToBuffer(mesh, opt, bin, views, accessors, prim, outError)) {
    return false;
  }

  if (!writeBinaryFile(binPath, bin, outError)) {
    return false;
  }

  // ---- Build JSON (.gltf) ----
  // We emit a small, valid subset of glTF 2.0.
  const bool hasTex = !opt.baseColorTextureUri.empty();

  std::ostringstream js;
  js.imbue(std::locale::classic());
  js << std::setprecision(9);

  js << "{\n";
  js << "  \"asset\": { \"version\": \"2.0\", \"generator\": \"stellar_forge\" },\n";
  js << "  \"scene\": 0,\n";
  js << "  \"scenes\": [ { \"nodes\": [ 0 ] } ],\n";
  js << "  \"nodes\": [ { \"name\": \"" << jsonEscape(opt.nodeName) << "\", \"mesh\": 0 } ],\n";

  // Mesh + primitive
  js << "  \"meshes\": [ { \"name\": \"" << jsonEscape(opt.meshName) << "\", \"primitives\": [ {\n";
  js << "    \"attributes\": { \"POSITION\": " << prim.accPos
     << ", \"NORMAL\": " << prim.accNrm
     << ", \"TEXCOORD_0\": " << prim.accUv;
  if (prim.hasTangents()) {
    js << ", \"TANGENT\": " << prim.accTan;
  }
  js << " },\n";
  js << "    \"indices\": " << prim.accIdx << ",\n";
  js << "    \"mode\": 4,\n"; // TRIANGLES
  js << "    \"material\": 0\n";
  js << "  } ] } ],\n";

  // Material
  js << "  \"materials\": [ {\n";
  js << "    \"name\": \"" << jsonEscape(opt.materialName) << "\",\n";
  js << "    \"doubleSided\": " << (opt.doubleSided ? "true" : "false") << ",\n";
  js << "    \"pbrMetallicRoughness\": {\n";
  js << "      \"baseColorFactor\": [ "
     << opt.baseColorFactorR << ", "
     << opt.baseColorFactorG << ", "
     << opt.baseColorFactorB << ", "
     << opt.baseColorFactorA << " ],\n";
  js << "      \"metallicFactor\": " << opt.metallicFactor << ",\n";
  js << "      \"roughnessFactor\": " << opt.roughnessFactor;
  if (hasTex) {
    js << ",\n      \"baseColorTexture\": { \"index\": 0 }\n";
  } else {
    js << "\n";
  }
  js << "    }\n";
  js << "  } ],\n";

  // Optional texture plumbing.
  if (hasTex) {
    js << "  \"samplers\": [ { \"magFilter\": 9729, \"minFilter\": 9987, \"wrapS\": 10497, \"wrapT\": 10497 } ],\n";
    js << "  \"images\": [ { \"uri\": \"" << jsonEscape(opt.baseColorTextureUri) << "\" } ],\n";
    js << "  \"textures\": [ { \"sampler\": 0, \"source\": 0 } ],\n";
  }

  // Buffers
  js << "  \"buffers\": [ { \"uri\": \"" << jsonEscape(binUri) << "\", \"byteLength\": " << bin.size() << " } ],\n";

  // BufferViews + Accessors
  emitBufferViews(js, views);
  emitAccessors(js, accessors);

  js << "}\n";

  if (!writeTextFile(outPath, js.str(), outError)) {
    return false;
  }

  return true;
}

bool exportMeshLodsToGltf(const std::string& gltfPath,
                          const SdfMeshData& lod0,
                          std::span<const SdfMeshData> lods,
                          const GltfExportOptions& opt,
                          std::string* outError) {
  // If no extra LODs were provided, fall back to the single-mesh exporter.
  if (lods.empty()) {
    return exportMeshToGltf(gltfPath, lod0, opt, outError);
  }

  if (gltfPath.empty()) {
    setErr(outError, "exportMeshLodsToGltf: empty glTF path.");
    return false;
  }

  fs::path outPath(gltfPath);
  fs::path outDir = outPath.parent_path();
  if (outDir.empty()) outDir = fs::path(".");

  if (!ensureOutDir(outDir, outError)) return false;

  fs::path binPath = outPath;
  binPath.replace_extension(".bin");
  const std::string binUri = binPath.filename().string();

  const std::size_t meshCount = 1 + lods.size();

  std::vector<std::uint8_t> bin;
  bin.reserve(1024);

  std::vector<ViewInfo> views;
  std::vector<AccessorInfo> accessors;
  std::vector<PackedPrimitiveAccessors> prims;
  prims.resize(meshCount);

  // Pack LOD0
  if (!packMeshToBuffer(lod0, opt, bin, views, accessors, prims[0], outError)) {
    return false;
  }

  // Pack LOD1..N
  for (std::size_t i = 0; i < lods.size(); ++i) {
    if (!packMeshToBuffer(lods[i], opt, bin, views, accessors, prims[i + 1], outError)) {
      return false;
    }
  }

  if (!writeBinaryFile(binPath, bin, outError)) {
    return false;
  }

  // ---- Build JSON (.gltf) ----
  const bool hasTex = !opt.baseColorTextureUri.empty();

  std::ostringstream js;
  js.imbue(std::locale::classic());
  js << std::setprecision(9);

  js << "{\n";
  js << "  \"asset\": { \"version\": \"2.0\", \"generator\": \"stellar_forge\" },\n";
  js << "  \"extensionsUsed\": [ \"MSFT_lod\" ],\n";
  js << "  \"scene\": 0,\n";
  js << "  \"scenes\": [ { \"nodes\": [ 0 ] } ],\n";

  // Nodes
  js << "  \"nodes\": [\n";
  for (std::size_t ni = 0; ni < meshCount; ++ni) {
    const bool isRoot = (ni == 0);
    std::string nodeName = opt.nodeName;
    if (!isRoot) {
      nodeName += "_lod" + std::to_string(ni);
    }

    js << "    { \"name\": \"" << jsonEscape(nodeName) << "\", \"mesh\": " << ni;

    if (isRoot) {
      // MSFT_lod: reference child node ids [1..meshCount-1].
      js << ", \"extensions\": { \"MSFT_lod\": { \"ids\": [ ";
      for (std::size_t li = 1; li < meshCount; ++li) {
        js << li;
        if (li + 1 < meshCount) js << ", ";
      }
      js << " ] } }";
    }

    js << " }";
    js << ((ni + 1 == meshCount) ? "\n" : ",\n");
  }
  js << "  ],\n";

  // Meshes
  js << "  \"meshes\": [\n";
  for (std::size_t mi = 0; mi < meshCount; ++mi) {
    std::string meshName = opt.meshName + "_lod" + std::to_string(mi);
    const auto& prim = prims[mi];

    js << "    { \"name\": \"" << jsonEscape(meshName) << "\", \"primitives\": [ {\n";
    js << "      \"attributes\": { \"POSITION\": " << prim.accPos
       << ", \"NORMAL\": " << prim.accNrm
       << ", \"TEXCOORD_0\": " << prim.accUv;
    if (prim.hasTangents()) {
      js << ", \"TANGENT\": " << prim.accTan;
    }
    js << " },\n";
    js << "      \"indices\": " << prim.accIdx << ",\n";
    js << "      \"mode\": 4,\n";
    js << "      \"material\": 0\n";
    js << "    } ] }";
    js << ((mi + 1 == meshCount) ? "\n" : ",\n");
  }
  js << "  ],\n";

  // Material (shared)
  js << "  \"materials\": [ {\n";
  js << "    \"name\": \"" << jsonEscape(opt.materialName) << "\",\n";
  js << "    \"doubleSided\": " << (opt.doubleSided ? "true" : "false") << ",\n";
  js << "    \"pbrMetallicRoughness\": {\n";
  js << "      \"baseColorFactor\": [ "
     << opt.baseColorFactorR << ", "
     << opt.baseColorFactorG << ", "
     << opt.baseColorFactorB << ", "
     << opt.baseColorFactorA << " ],\n";
  js << "      \"metallicFactor\": " << opt.metallicFactor << ",\n";
  js << "      \"roughnessFactor\": " << opt.roughnessFactor;
  if (hasTex) {
    js << ",\n      \"baseColorTexture\": { \"index\": 0 }\n";
  } else {
    js << "\n";
  }
  js << "    }\n";
  js << "  } ],\n";

  // Optional texture plumbing.
  if (hasTex) {
    js << "  \"samplers\": [ { \"magFilter\": 9729, \"minFilter\": 9987, \"wrapS\": 10497, \"wrapT\": 10497 } ],\n";
    js << "  \"images\": [ { \"uri\": \"" << jsonEscape(opt.baseColorTextureUri) << "\" } ],\n";
    js << "  \"textures\": [ { \"sampler\": 0, \"source\": 0 } ],\n";
  }

  // Buffers
  js << "  \"buffers\": [ { \"uri\": \"" << jsonEscape(binUri) << "\", \"byteLength\": " << bin.size() << " } ],\n";

  // BufferViews + Accessors
  emitBufferViews(js, views);
  emitAccessors(js, accessors);

  js << "}\n";

  if (!writeTextFile(outPath, js.str(), outError)) {
    return false;
  }

  return true;
}

} // namespace stellar::render
