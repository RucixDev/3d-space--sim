#pragma once

#include "stellar/render/Shader.h"
#include "stellar/render/Texture.h"

#include <vector>

namespace stellar::render {

enum class PointBlendMode {
  Alpha = 0,
  Additive,
};

struct PointVertex {
  float px, py, pz;
  float cr, cg, cb;
  float a;
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
  void drawPoints(const std::vector<PointVertex>& points, PointBlendMode blend = PointBlendMode::Alpha);

  // Draw points as textured point-sprites (samples the full sprite using gl_PointCoord).
  // Useful for stars/particles that want a smoother falloff than a hard-coded circle.
  void drawPointsSprite(const std::vector<PointVertex>& points,
                        const Texture2D& sprite,
                        PointBlendMode blend = PointBlendMode::Alpha);

private:
  ShaderProgram shader_{};
  unsigned int vao_{0};
  unsigned int vbo_{0};

  float view_[16]{};
  float proj_[16]{};
};

} // namespace stellar::render
