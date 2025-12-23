#pragma once

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

  static Mesh makeCube();
  static Mesh makeUvSphere(int slices = 32, int stacks = 16);

private:
  void destroy();

  unsigned int vao_{0};
  unsigned int vbo_{0};
  unsigned int ebo_{0};
  std::uint32_t indexCount_{0};
};

} // namespace stellar::render
