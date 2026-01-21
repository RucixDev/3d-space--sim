#pragma once

#include "stellar/render/PointRenderer.h"
#include "stellar/render/Shader.h"
#include "stellar/render/Starfield.h"
#include "stellar/render/Texture.h"

#include <cstddef>
#include <string>
#include <vector>

namespace stellar::render {

// GPU-driven starfield renderer.
//
// This eliminates the per-frame CPU expansion (cameraPos + dir*radius) and
// per-frame PointVertex streaming performed by Starfield::update().
// Instead, we upload the star directions + appearance parameters once after
// regenerate(), and compute position + twinkle in the vertex shader each frame.
//
// Intended as a CPU-side FPS win on drivers/CPUs that struggle with large
// per-frame dynamic VBO uploads (common on integrated GPUs).
class StarfieldGpuRenderer {
public:
  StarfieldGpuRenderer() = default;
  ~StarfieldGpuRenderer();

  StarfieldGpuRenderer(const StarfieldGpuRenderer&) = delete;
  StarfieldGpuRenderer& operator=(const StarfieldGpuRenderer&) = delete;

  bool init(std::string* outError = nullptr);

  void setViewProj(const float* view, const float* proj);
  void setCameraPos(float x, float y, float z);
  void setRadius(float r);
  void setTimeSeconds(float t);

  void upload(const std::vector<Starfield::GpuStar>& stars);
  void upload(const Starfield& starfield);

  // Draws the starfield as point sprites. If spriteTex is nullptr, we render a
  // soft circular sprite in the fragment shader.
  void draw(const Texture2D* spriteTex, PointBlendMode blend = PointBlendMode::Additive);

  std::size_t starCount() const { return starCount_; }

private:
  void ensureCapacityBytes(std::size_t bytes);

  ShaderProgram shader_;

  unsigned int vao_{0};
  unsigned int vbo_{0};
  std::size_t vboCapacityBytes_{0};
  std::size_t starCount_{0};

  float view_[16]{};
  float proj_[16]{};
  float camPos_[3]{0.0f, 0.0f, 0.0f};
  float radius_{16000.0f};
  float time_{0.0f};
};

} // namespace stellar::render
