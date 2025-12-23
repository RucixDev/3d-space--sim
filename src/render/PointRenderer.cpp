#include "stellar/render/PointRenderer.h"

#include "stellar/render/Gl.h"

#include <cstring>

namespace stellar::render {

PointRenderer::~PointRenderer() {
  if (vbo_) gl::DeleteBuffers(1, &vbo_);
  if (vao_) gl::DeleteVertexArrays(1, &vao_);
  vbo_ = vao_ = 0;
}

static const char* kVS = R"GLSL(
#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aColor;
layout(location=2) in float aSize;

uniform mat4 uView;
uniform mat4 uProj;

out vec3 vColor;

void main() {
  gl_Position = uProj * uView * vec4(aPos, 1.0);
  gl_PointSize = aSize;
  vColor = aColor;
}
)GLSL";

static const char* kFS = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 FragColor;

void main() {
  // soft circle
  vec2 p = gl_PointCoord * 2.0 - 1.0;
  float d = dot(p,p);
  float a = smoothstep(1.0, 0.7, d);
  FragColor = vec4(vColor, a);
}
)GLSL";

bool PointRenderer::init(std::string* outError) {
  if (!shader_.build(kVS, kFS, outError)) return false;

  gl::GenVertexArrays(1, &vao_);
  gl::GenBuffers(1, &vbo_);

  gl::BindVertexArray(vao_);
  gl::BindBuffer(GL_ARRAY_BUFFER, vbo_);
  gl::BufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);

  gl::EnableVertexAttribArray(0);
  gl::VertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(PointVertex), (void*)offsetof(PointVertex, px));
  gl::EnableVertexAttribArray(1);
  gl::VertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(PointVertex), (void*)offsetof(PointVertex, cr));
  gl::EnableVertexAttribArray(2);
  gl::VertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(PointVertex), (void*)offsetof(PointVertex, size));

  gl::BindVertexArray(0);

  for (int i = 0; i < 16; ++i) {
    view_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
    proj_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
  }

  return true;
}

void PointRenderer::setViewProj(const float* view, const float* proj) {
  std::memcpy(view_, view, sizeof(float) * 16);
  std::memcpy(proj_, proj, sizeof(float) * 16);
}

void PointRenderer::drawPoints(const std::vector<PointVertex>& points) {
  if (points.empty()) return;

  glEnable(GL_PROGRAM_POINT_SIZE);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  shader_.bind();
  shader_.setUniformMat4("uView", view_);
  shader_.setUniformMat4("uProj", proj_);

  gl::BindVertexArray(vao_);
  gl::BindBuffer(GL_ARRAY_BUFFER, vbo_);
  gl::BufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(points.size() * sizeof(PointVertex)),
                 points.data(),
                 GL_DYNAMIC_DRAW);

  glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(points.size()));

  gl::BindVertexArray(0);
}

} // namespace stellar::render
