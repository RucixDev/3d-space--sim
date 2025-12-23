#pragma once

#include "stellar/math/Mat4.h"
#include "stellar/render/Gl.h"
#include "stellar/render/Shader.h"

#include <span>
#include <string>

namespace stellar::render {

struct PointVertex {
  float px = 0.0f;
  float py = 0.0f;
  float pz = 0.0f;

  float r = 1.0f;
  float g = 1.0f;
  float b = 1.0f;

  float size = 4.0f;
};

class PointRenderer {
public:
  bool init(std::string* outError = nullptr);
  void shutdown();

  void resize(int width, int height);

  void beginFrame(float clearR, float clearG, float clearB, float clearA = 1.0f);
  void drawPoints(std::span<const PointVertex> pts, const stellar::math::Mat4f& viewProj);

private:
  ShaderProgram m_program;
  gl::GLuint m_vao = 0;
  gl::GLuint m_vbo = 0;
  gl::GLint m_uVP = -1;
  int m_width = 1;
  int m_height = 1;
};

} // namespace stellar::render
