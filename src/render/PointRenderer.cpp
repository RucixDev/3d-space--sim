#include "stellar/render/PointRenderer.h"

#include <cstddef>

namespace stellar::render {
namespace {

constexpr const char* kVert = R"GLSL(
#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aColor;
layout(location = 2) in float aSize;

uniform mat4 uVP;

out vec3 vColor;

void main() {
  gl_Position = uVP * vec4(aPos, 1.0);
  vColor = aColor;
  gl_PointSize = aSize;
}
)GLSL";

constexpr const char* kFrag = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 FragColor;
void main() {
  FragColor = vec4(vColor, 1.0);
}
)GLSL";

}

bool PointRenderer::init(std::string* outError) {
  std::string err;
  if (!m_program.create(kVert, kFrag, &err)) {
    if (outError) *outError = err;
    return false;
  }

  m_uVP = m_program.uniformLocation("uVP");

  gl::GenVertexArrays(1, &m_vao);
  gl::BindVertexArray(m_vao);

  gl::GenBuffers(1, &m_vbo);
  gl::BindBuffer(gl::GL_ARRAY_BUFFER, m_vbo);
  gl::BufferData(gl::GL_ARRAY_BUFFER, 0, nullptr, gl::GL_DYNAMIC_DRAW);

  const gl::GLsizei stride = static_cast<gl::GLsizei>(sizeof(PointVertex));

  gl::EnableVertexAttribArray(0);
  gl::VertexAttribPointer(0, 3, gl::GL_FLOAT, gl::GL_FALSE, stride, reinterpret_cast<void*>(0));

  gl::EnableVertexAttribArray(1);
  gl::VertexAttribPointer(1, 3, gl::GL_FLOAT, gl::GL_FALSE, stride, reinterpret_cast<void*>(offsetof(PointVertex, r)));

  gl::EnableVertexAttribArray(2);
  gl::VertexAttribPointer(2, 1, gl::GL_FLOAT, gl::GL_FALSE, stride, reinterpret_cast<void*>(offsetof(PointVertex, size)));

  // Render state
  gl::Enable(gl::GL_PROGRAM_POINT_SIZE);
  gl::Enable(gl::GL_BLEND);
  gl::BlendFunc(gl::GL_SRC_ALPHA, gl::GL_ONE_MINUS_SRC_ALPHA);
  gl::Enable(gl::GL_DEPTH_TEST);

  gl::BindBuffer(gl::GL_ARRAY_BUFFER, 0);
  gl::BindVertexArray(0);
  return true;
}

void PointRenderer::shutdown() {
  m_program.destroy();

  if (m_vbo) {
    gl::DeleteBuffers(1, &m_vbo);
    m_vbo = 0;
  }
  if (m_vao) {
    gl::DeleteVertexArrays(1, &m_vao);
    m_vao = 0;
  }
}

void PointRenderer::resize(int width, int height) {
  m_width = (width <= 0) ? 1 : width;
  m_height = (height <= 0) ? 1 : height;
  gl::Viewport(0, 0, m_width, m_height);
}

void PointRenderer::beginFrame(float clearR, float clearG, float clearB, float clearA) {
  gl::ClearColor(clearR, clearG, clearB, clearA);
  gl::Clear(gl::GL_COLOR_BUFFER_BIT | gl::GL_DEPTH_BUFFER_BIT);
}

void PointRenderer::drawPoints(std::span<const PointVertex> pts, const stellar::math::Mat4f& viewProj) {
  m_program.bind();
  gl::UniformMatrix4fv(m_uVP, 1, gl::GL_FALSE, viewProj.m);

  gl::BindVertexArray(m_vao);
  gl::BindBuffer(gl::GL_ARRAY_BUFFER, m_vbo);

  const auto bytes = static_cast<gl::GLsizeiptr>(pts.size_bytes());
  gl::BufferData(gl::GL_ARRAY_BUFFER, bytes, pts.data(), gl::GL_DYNAMIC_DRAW);
  gl::DrawArrays(gl::GL_POINTS, 0, static_cast<gl::GLsizei>(pts.size()));

  gl::BindBuffer(gl::GL_ARRAY_BUFFER, 0);
  gl::BindVertexArray(0);
}

} // namespace stellar::render