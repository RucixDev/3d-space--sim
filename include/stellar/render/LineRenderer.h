#pragma once

#include "stellar/render/Shader.h"

#include <vector>
#include <cstddef>

namespace stellar::render {

struct LineVertex {
  float px, py, pz;
  float cr, cg, cb;
};

class LineRenderer {
public:
  LineRenderer() = default;
  ~LineRenderer();

  LineRenderer(const LineRenderer&) = delete;
  LineRenderer& operator=(const LineRenderer&) = delete;

  bool init(std::string* outError = nullptr);

  void setViewProj(const float* view, const float* proj);
  void drawLines(const std::vector<LineVertex>& vertices); // pairs

private:
  ShaderProgram shader_{};
  unsigned int vao_{0};
  unsigned int vbo_{0};

  // Capacity of the VBO allocation (bytes). Used to avoid reallocating the
  // data store every frame.
  std::size_t vboCapacityBytes_{0};

  float view_[16]{};
  float proj_[16]{};
};

} // namespace stellar::render
