#pragma once

#include "stellar/render/SdfMesher.h" // SdfMeshData

#include <span>
#include <string>

namespace stellar::render {

// Minimal glTF 2.0 mesh exporter.
//
// This exporter is intentionally small and dependency-free:
//  - Writes .gltf (JSON) + external .bin buffer.
//  - Exports one mesh with POSITION/NORMAL/TEXCOORD_0 + triangle indices.
//  - Optionally references a baseColorTexture (PNG) for PBR metallic-roughness.
//
// NOTE:
//  - If you want normal mapping in glTF, you need tangents (TANGENT). Tangents
//    are optional and can be exported by enabling GltfExportOptions::exportTangents.

struct GltfExportOptions {
  // Optional. When non-empty, an image/texture is emitted and the material's
  // pbrMetallicRoughness.baseColorTexture refers to this URI.
  //
  // This should usually be a *file name*, not an absolute path, so the glTF is
  // portable (e.g. "my_mesh_albedo.png").
  std::string baseColorTextureUri{};

  // Material factors (multiplied against texture if present).
  float baseColorFactorR{1.0f};
  float baseColorFactorG{1.0f};
  float baseColorFactorB{1.0f};
  float baseColorFactorA{1.0f};
  float metallicFactor{0.0f};
  float roughnessFactor{1.0f};

  bool doubleSided{true};

  // If true, exports TEXCOORD_0 with v = (1 - v). This is useful when the mesh
  // uses OpenGL-style UVs but the texture is written with a top-left origin.
  bool flipTexcoordV{false};

  // If true, exports per-vertex tangents as TANGENT (VEC4 float) using a
  // standard UV-based tangent frame build.
  //
  // - This is recommended if you plan to use normal maps in downstream tools.
  // - When UVs are degenerate, the exporter falls back to an arbitrary stable
  //   tangent frame derived from the vertex normal.
  bool exportTangents{false};

  std::string meshName{"mesh0"};
  std::string materialName{"mat0"};
  std::string nodeName{"node0"};
};

// Export the given mesh as a glTF 2.0 asset.
//
// Output paths:
//  - Writes <gltfPath>
//  - Writes <gltfPath stem>.bin in the same directory
//
// The glTF will reference the .bin by file name (relative URI). If
// options.baseColorTextureUri is non-empty, it will be referenced as-is.
//
// Returns false on failure and optionally fills outError.
bool exportMeshToGltf(const std::string& gltfPath,
                      const SdfMeshData& mesh,
                      const GltfExportOptions& options,
                      std::string* outError = nullptr);

// Export a mesh LOD chain as a single glTF asset.
//
// This uses the `MSFT_lod` extension so glTF viewers that support it can pick
// an appropriate LOD automatically.
//
// - lod0 is required.
// - lods contains LOD1..LOD(N) (may be empty).
// - The exporter writes <gltfPath> and <gltfPath stem>.bin (same folder).
// - The same material is used for all LOD meshes.
bool exportMeshLodsToGltf(const std::string& gltfPath,
                          const SdfMeshData& lod0,
                          std::span<const SdfMeshData> lods,
                          const GltfExportOptions& options,
                          std::string* outError = nullptr);

} // namespace stellar::render
