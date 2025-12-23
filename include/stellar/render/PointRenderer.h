#pragma once

#include "stellar/render/Shader.h"

#include <vector>

namespace stellar::render {

struct PointVertex {
  float px, py, pz;
  float cr, cg, cb;
  float size;
};

class PointRenderer {
public:
  PointRenderer() = default;
  ~PointRenderer();

  PointRenderer(const PointRenderer&) = delete;
  PointRenderer& operator=(const PointRenderer&) = delete;

  bool init(std::string* outError = nullptr);

  void setViewProj(const float* view, const float* proj);
  void drawPoints(const std::vector<PointVertex>& points);

private:
  ShaderProgram shader_{};
  unsigned int vao_{0};
  unsigned int vbo_{0};

  float view_[16]{};
  float proj_[16]{};
};

} // namespace stellar::render
