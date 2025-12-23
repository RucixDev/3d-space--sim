#include "stellar/render/LineRenderer.h"

#include "stellar/render/Gl.h"

#include <cstring>

namespace stellar::render {

LineRenderer::~LineRenderer() {
  if (vbo_) gl::DeleteBuffers(1, &vbo_);
  if (vao_) gl::DeleteVertexArrays(1, &vao_);
  vbo_ = vao_ = 0;
}

static const char* kVS = R"GLSL(
#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aColor;

uniform mat4 uView;
uniform mat4 uProj;

out vec3 vColor;

void main() {
  gl_Position = uProj * uView * vec4(aPos, 1.0);
  vColor = aColor;
}
)GLSL";

static const char* kFS = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 FragColor;
void main() {
  FragColor = vec4(vColor, 1.0);
}
)GLSL";

bool LineRenderer::init(std::string* outError) {
  if (!shader_.build(kVS, kFS, outError)) return false;

  gl::GenVertexArrays(1, &vao_);
  gl::GenBuffers(1, &vbo_);

  gl::BindVertexArray(vao_);
  gl::BindBuffer(GL_ARRAY_BUFFER, vbo_);
  gl::BufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);

  gl::EnableVertexAttribArray(0);
  gl::VertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(LineVertex), (void*)offsetof(LineVertex, px));
  gl::EnableVertexAttribArray(1);
  gl::VertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(LineVertex), (void*)offsetof(LineVertex, cr));

  gl::BindVertexArray(0);

  for (int i = 0; i < 16; ++i) {
    view_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
    proj_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
  }

  return true;
}

void LineRenderer::setViewProj(const float* view, const float* proj) {
  std::memcpy(view_, view, sizeof(float) * 16);
  std::memcpy(proj_, proj, sizeof(float) * 16);
}

void LineRenderer::drawLines(const std::vector<LineVertex>& vertices) {
  if (vertices.empty()) return;

  shader_.bind();
  shader_.setUniformMat4("uView", view_);
  shader_.setUniformMat4("uProj", proj_);

  gl::BindVertexArray(vao_);
  gl::BindBuffer(GL_ARRAY_BUFFER, vbo_);
  gl::BufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(vertices.size() * sizeof(LineVertex)),
                 vertices.data(),
                 GL_DYNAMIC_DRAW);

  glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(vertices.size()));

  gl::BindVertexArray(0);
}

} // namespace stellar::render
