#include "stellar/render/GltfExport.h"

#include <algorithm>
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

} // namespace

bool exportMeshToGltf(const std::string& gltfPath,
                      const SdfMeshData& mesh,
                      const GltfExportOptions& opt,
                      std::string* outError) {
  if (gltfPath.empty()) {
    setErr(outError, "exportMeshToGltf: empty glTF path.");
    return false;
  }
  if (mesh.vertices.empty() || mesh.indices.empty()) {
    setErr(outError, "exportMeshToGltf: mesh is empty.");
    return false;
  }
  if (mesh.indices.size() % 3 != 0) {
    setErr(outError, "exportMeshToGltf: index count is not a multiple of 3.");
    return false;
  }

  fs::path outPath(gltfPath);
  fs::path outDir = outPath.parent_path();
  if (outDir.empty()) outDir = fs::path(".");

  std::error_code ec;
  if (!fs::exists(outDir, ec)) {
    fs::create_directories(outDir, ec);
  }
  if (ec) {
    setErr(outError, "exportMeshToGltf: failed to create directory: " + outDir.string());
    return false;
  }

  fs::path binPath = outPath;
  binPath.replace_extension(".bin");
  const std::string binUri = binPath.filename().string();

  // ---- Pack binary buffer ----
  // We keep attributes in separate contiguous blocks for simplicity.
  // Ensure 4-byte alignment between blocks.
  std::vector<std::uint8_t> bin;
  bin.reserve(mesh.vertices.size() * (3 + 3 + 2) * sizeof(float) + mesh.indices.size() * sizeof(std::uint32_t) + 64);

  ViewInfo posView{}, nrmView{}, uvView{}, idxView{};

  // Positions (VEC3 float)
  padTo(bin, 4);
  posView.offset = bin.size();
  for (const auto& v : mesh.vertices) {
    appendPod(bin, v.px);
    appendPod(bin, v.py);
    appendPod(bin, v.pz);
  }
  padTo(bin, 4);
  posView.length = bin.size() - posView.offset;
  posView.target = 34962; // ARRAY_BUFFER

  // Normals (VEC3 float)
  padTo(bin, 4);
  nrmView.offset = bin.size();
  for (const auto& v : mesh.vertices) {
    appendPod(bin, v.nx);
    appendPod(bin, v.ny);
    appendPod(bin, v.nz);
  }
  padTo(bin, 4);
  nrmView.length = bin.size() - nrmView.offset;
  nrmView.target = 34962;

  // UVs (VEC2 float)
  padTo(bin, 4);
  uvView.offset = bin.size();
  for (const auto& v : mesh.vertices) {
    const float u = v.u;
    const float vv = opt.flipTexcoordV ? (1.0f - v.v) : v.v;
    appendPod(bin, u);
    appendPod(bin, vv);
  }
  padTo(bin, 4);
  uvView.length = bin.size() - uvView.offset;
  uvView.target = 34962;

  // Indices (SCALAR uint32)
  padTo(bin, 4);
  idxView.offset = bin.size();
  for (std::uint32_t i : mesh.indices) {
    appendPod(bin, i);
  }
  padTo(bin, 4);
  idxView.length = bin.size() - idxView.offset;
  idxView.target = 34963; // ELEMENT_ARRAY_BUFFER

  // Write .bin
  {
    std::ofstream f(binPath, std::ios::binary | std::ios::out | std::ios::trunc);
    if (!f) {
      setErr(outError, "exportMeshToGltf: failed to open for writing: " + binPath.string());
      return false;
    }
    f.write(reinterpret_cast<const char*>(bin.data()), (std::streamsize)bin.size());
    if (!f) {
      setErr(outError, "exportMeshToGltf: failed to write: " + binPath.string());
      return false;
    }
  }

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
  js << "    \"attributes\": { \"POSITION\": 0, \"NORMAL\": 1, \"TEXCOORD_0\": 2 },\n";
  js << "    \"indices\": 3,\n";
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

  // BufferViews
  js << "  \"bufferViews\": [\n";
  auto emitView = [&](const ViewInfo& v, bool last) {
    js << "    { \"buffer\": 0, \"byteOffset\": " << v.offset
       << ", \"byteLength\": " << v.length;
    if (v.target != 0) {
      js << ", \"target\": " << v.target;
    }
    js << " }";
    js << (last ? "\n" : ",\n");
  };
  emitView(posView, false);
  emitView(nrmView, false);
  emitView(uvView, false);
  emitView(idxView, true);
  js << "  ],\n";

  // Accessors
  const std::size_t vCount = mesh.vertices.size();
  const std::size_t iCount = mesh.indices.size();

  js << "  \"accessors\": [\n";

  // POSITION accessor (index 0)
  js << "    { \"bufferView\": 0, \"byteOffset\": 0, \"componentType\": 5126, \"count\": " << vCount
     << ", \"type\": \"VEC3\", \"min\": [ " << minX << ", " << minY << ", " << minZ
     << " ], \"max\": [ " << maxX << ", " << maxY << ", " << maxZ << " ] },\n";

  // NORMAL accessor (index 1)
  js << "    { \"bufferView\": 1, \"byteOffset\": 0, \"componentType\": 5126, \"count\": " << vCount
     << ", \"type\": \"VEC3\" },\n";

  // TEXCOORD_0 accessor (index 2)
  js << "    { \"bufferView\": 2, \"byteOffset\": 0, \"componentType\": 5126, \"count\": " << vCount
     << ", \"type\": \"VEC2\" },\n";

  // INDICES accessor (index 3)
  // NOTE: UNSIGNED_INT (5125) is allowed for indices. (Not valid for vertex attributes.)
  js << "    { \"bufferView\": 3, \"byteOffset\": 0, \"componentType\": 5125, \"count\": " << iCount
     << ", \"type\": \"SCALAR\" }\n";

  js << "  ]\n";
  js << "}\n";

  // Write .gltf
  {
    std::ofstream f(outPath, std::ios::out | std::ios::trunc);
    if (!f) {
      setErr(outError, "exportMeshToGltf: failed to open for writing: " + outPath.string());
      return false;
    }
    const std::string s = js.str();
    f.write(s.c_str(), (std::streamsize)s.size());
    if (!f) {
      setErr(outError, "exportMeshToGltf: failed to write: " + outPath.string());
      return false;
    }
  }

  return true;
}

} // namespace stellar::render
