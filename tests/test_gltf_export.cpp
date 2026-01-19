#include "test_harness.h"

#include "stellar/render/GltfExport.h"

#include <cmath>
#include <cstring>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace {

static std::string readTextFile(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::in);
  if (!f) return {};
  return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

static std::vector<std::uint8_t> readBinaryFile(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary | std::ios::in);
  if (!f) return {};
  f.seekg(0, std::ios::end);
  const std::streamsize sz = f.tellg();
  f.seekg(0, std::ios::beg);
  if (sz <= 0) return {};
  std::vector<std::uint8_t> buf((std::size_t)sz);
  f.read(reinterpret_cast<char*>(buf.data()), sz);
  if (!f) return {};
  return buf;
}

static float readF32(const std::vector<std::uint8_t>& buf, std::size_t byteOffset) {
  float v = 0.0f;
  if (byteOffset + sizeof(float) > buf.size()) return v;
  std::memcpy(&v, buf.data() + byteOffset, sizeof(float));
  return v;
}

static std::size_t align4(std::size_t x) {
  return (x + 3u) & ~3u;
}

static bool approx(float a, float b, float eps = 1e-3f) {
  return std::fabs(a - b) <= eps;
}

} // namespace

int test_gltf_export() {
  int failures = 0;

  namespace fs = std::filesystem;
  const fs::path outDir = fs::path("test_out");
  std::error_code ec;
  fs::create_directories(outDir, ec);

  // --- Build a simple quad (two triangles) ---
  stellar::render::SdfMeshData mesh;
  mesh.vertices.resize(4);
  mesh.vertices[0] = {0.0f, 0.0f, 0.0f,  0.0f, 0.0f, 1.0f,  0.0f, 0.0f};
  mesh.vertices[1] = {1.0f, 0.0f, 0.0f,  0.0f, 0.0f, 1.0f,  1.0f, 0.0f};
  mesh.vertices[2] = {1.0f, 1.0f, 0.0f,  0.0f, 0.0f, 1.0f,  1.0f, 1.0f};
  mesh.vertices[3] = {0.0f, 1.0f, 0.0f,  0.0f, 0.0f, 1.0f,  0.0f, 1.0f};
  mesh.indices = {0, 1, 2, 0, 2, 3};

  const fs::path gltfPath = outDir / "quad_tangent.gltf";
  const fs::path binPath = outDir / "quad_tangent.bin";

  stellar::render::GltfExportOptions opt;
  opt.meshName = "quad";
  opt.nodeName = "quad_node";
  opt.materialName = "mat";
  opt.doubleSided = true;
  opt.flipTexcoordV = false;
  opt.exportTangents = true;

  std::string err;
  CHECK(stellar::render::exportMeshToGltf(gltfPath.string(), mesh, opt, &err));
  CHECK(err.empty());

  const std::string gltfText = readTextFile(gltfPath);
  CHECK(!gltfText.empty());
  CHECK(gltfText.find("\"TANGENT\"") != std::string::npos);

  const auto bin = readBinaryFile(binPath);
  CHECK(!bin.empty());

  // Verify tangent data layout + values. The exporter packs:
  //   [pos][nrm][uv][tan?][indices]
  const std::size_t vCount = mesh.vertices.size();
  std::size_t off = 0;
  off = align4(off);
  off += vCount * 3 * sizeof(float); // pos
  off = align4(off);
  off += vCount * 3 * sizeof(float); // nrm
  off = align4(off);
  off += vCount * 2 * sizeof(float); // uv
  off = align4(off);
  const std::size_t tanOffset = off;

  // Each tangent is (1,0,0,1) for this quad.
  for (std::size_t i = 0; i < vCount; ++i) {
    const std::size_t base = tanOffset + i * 4 * sizeof(float);
    const float tx = readF32(bin, base + 0);
    const float ty = readF32(bin, base + 4);
    const float tz = readF32(bin, base + 8);
    const float tw = readF32(bin, base + 12);

    CHECK(approx(tx, 1.0f));
    CHECK(approx(ty, 0.0f));
    CHECK(approx(tz, 0.0f));
    CHECK(approx(tw, 1.0f));
  }

  // --- LOD export (MSFT_lod) ---
  stellar::render::SdfMeshData lod1;
  lod1.vertices.resize(3);
  lod1.vertices[0] = {0.0f, 0.0f, 0.0f,  0.0f, 0.0f, 1.0f,  0.0f, 0.0f};
  lod1.vertices[1] = {1.0f, 0.0f, 0.0f,  0.0f, 0.0f, 1.0f,  1.0f, 0.0f};
  lod1.vertices[2] = {0.0f, 1.0f, 0.0f,  0.0f, 0.0f, 1.0f,  0.0f, 1.0f};
  lod1.indices = {0, 1, 2};

  const fs::path lodGltfPath = outDir / "quad_lods.gltf";
  const fs::path lodBinPath = outDir / "quad_lods.bin";

  std::vector<stellar::render::SdfMeshData> lods;
  lods.push_back(lod1);

  std::string lodErr;
  CHECK(stellar::render::exportMeshLodsToGltf(lodGltfPath.string(), mesh,
                                             std::span<const stellar::render::SdfMeshData>(lods.data(), lods.size()),
                                             opt, &lodErr));
  CHECK(lodErr.empty());

  const std::string lodText = readTextFile(lodGltfPath);
  CHECK(!lodText.empty());
  CHECK(lodText.find("MSFT_lod") != std::string::npos);
  CHECK(lodText.find("\"extensionsUsed\"") != std::string::npos);
  CHECK(lodText.find("\"TANGENT\"") != std::string::npos);
  CHECK(fs::exists(lodBinPath));

  return failures;
}
