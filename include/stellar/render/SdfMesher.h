#pragma once

#include "stellar/render/Mesh.h" // VertexPNUT (CPU-side)

#include <cstdint>
#include <functional>
#include <vector>

namespace stellar::render {

// Small, self-contained isosurface mesher for procedural geometry.
//
// - Works on any scalar field f(x,y,z).
// - Extracts an isosurface at a chosen iso value using Marching Tetrahedra.
// - Outputs CPU vertex/index arrays (upload via Mesh::upload()).
//
// Convention: by default, inside is f <= iso and outside is f > iso.
// For signed distance fields (SDFs), iso = 0 and f < 0 is "inside".

using ScalarField3D = std::function<float(float x, float y, float z)>;

struct SdfMeshData {
  std::vector<VertexPNUT> vertices;
  std::vector<std::uint32_t> indices;
};

struct SdfMesherParams {
  // Number of cells per axis in the cubic domain.
  // Total sampled grid points is (resolution+1)^3.
  int resolution{32};

  // Half-extent of the cubic domain: positions are sampled in [-bounds, +bounds].
  float bounds{1.25f};

  // Iso value threshold.
  float iso{0.0f};

  // Compute smooth per-vertex normals from the scalar field gradient (central differences).
  bool computeNormalsFromField{true};

  // Step size for gradient sampling (in world units).
  float normalEps{0.0025f};

  // Flip triangle winding so generated geometry is front-facing for
  // standard backface culling when the field is an SDF (positive outside).
  bool fixWindingFromNormals{true};
};

struct SdfMeshStats {
  std::size_t vertexCount{0};
  std::size_t indexCount{0};
  float minX{0}, minY{0}, minZ{0};
  float maxX{0}, maxY{0}, maxZ{0};
};

SdfMeshData meshIsosurfaceMarchingTetrahedra(const ScalarField3D& field,
                                             const SdfMesherParams& params = {});

SdfMeshStats measureSdfMesh(const SdfMeshData& mesh);

} // namespace stellar::render
