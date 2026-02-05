#pragma once

#include "stellar/math/Geometry.h"

#include <cstdint>
#include <vector>

namespace stellar::render {

struct VertexPNUT {
  float px, py, pz;
  float nx, ny, nz;
  float u, v;
};

class Mesh {
public:
  Mesh() = default;
  ~Mesh();

  Mesh(const Mesh&) = delete;
  Mesh& operator=(const Mesh&) = delete;

  Mesh(Mesh&&) noexcept;
  Mesh& operator=(Mesh&&) noexcept;

  void upload(const std::vector<VertexPNUT>& vertices, const std::vector<std::uint32_t>& indices);

  void bind() const;
  void draw() const;
  void drawInstanced(std::uint32_t instanceCount) const;

  std::uint32_t indexCount() const { return indexCount_; }

  // Local-space bounds for the uploaded mesh.
  //
  // These are computed at upload-time from the vertex positions.
  // For empty meshes, the bounds are non-finite.
  const stellar::math::Aabb3d& localBounds() const { return localBounds_; }
  bool hasLocalBounds() const { return localBounds_.isFinite(); }

  // Bounding sphere radius derived from the local AABB (half diagonal).
  // Intended for conservative culling/picking.
  double localBoundingRadius() const { return localBoundingRadius_; }

  static Mesh makeCube();
  static Mesh makeUvSphere(int slices = 32, int stacks = 16);

  // Create a flat annulus (ring) mesh in the XZ plane (Y=0).
  //
  // UV mapping:
  //   u = angle around the ring [0..1]
  //   v = radial coordinate [0..1] (inner->outer)
  //
  // When doubleSided is true, the mesh includes both top and bottom faces so
  // it is visible with backface culling enabled.
  static Mesh makeRing(int segments = 128,
                       float innerRadius = 0.65f,
                       float outerRadius = 1.25f,
                       bool doubleSided = true);

  // Lightweight procedural ship meshes (unit scale, forward +Z).
  //
  // These are intentionally simple, low-poly meshes used for ships in the
  // prototype so we don't rely on external model assets.
  static Mesh makeShipScout();
  static Mesh makeShipHauler();
  static Mesh makeShipFighter();

private:
  void destroy();

  unsigned int vao_{0};
  unsigned int vbo_{0};
  unsigned int ebo_{0};
  std::uint32_t indexCount_{0};

  // Cached local-space bounds.
  stellar::math::Aabb3d localBounds_{};
  double localBoundingRadius_{0.0};
};

} // namespace stellar::render
