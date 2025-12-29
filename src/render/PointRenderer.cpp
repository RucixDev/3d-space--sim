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
layout(location=2) in float aAlpha;
layout(location=3) in float aSize;

uniform mat4 uView;
uniform mat4 uProj;

out vec3 vColor;
out float vAlpha;

void main() {
  gl_Position = uProj * uView * vec4(aPos, 1.0);
  gl_PointSize = aSize;
  vColor = aColor;
  vAlpha = aAlpha;
}
)GLSL";

static const char* kFS = R"GLSL(
#version 330 core
in vec3 vColor;
in float vAlpha;
out vec4 FragColor;

uniform int uUseTex;
uniform sampler2D uTex;

void main() {
  // soft circle mask (keeps sprites round even if texture isn't)
  vec2 p = gl_PointCoord * 2.0 - 1.0;
  float d = dot(p,p);
  float circle = smoothstep(1.0, 0.7, d);

  if (uUseTex == 0) {
    float a = circle * vAlpha;
    FragColor = vec4(vColor, a);
  } else {
    vec4 t = texture(uTex, gl_PointCoord);
    float a = t.a * vAlpha * circle;
    vec3 rgb = t.rgb * vColor;
    FragColor = vec4(rgb, a);
  }
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
  gl::VertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(PointVertex), (void*)offsetof(PointVertex, a));
  gl::EnableVertexAttribArray(3);
  gl::VertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(PointVertex), (void*)offsetof(PointVertex, size));

  gl::BindVertexArray(0);

  for (int i = 0; i < 16; ++i) {
    view_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
    proj_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
  }

  // Default uniforms.
  shader_.bind();
  shader_.setUniform1i("uUseTex", 0);
  shader_.setUniform1i("uTex", 0);

  return true;
}

void PointRenderer::setViewProj(const float* view, const float* proj) {
  std::memcpy(view_, view, sizeof(float) * 16);
  std::memcpy(proj_, proj, sizeof(float) * 16);
}

void PointRenderer::drawPoints(const std::vector<PointVertex>& points, PointBlendMode blend) {
  if (points.empty()) return;

  glEnable(GL_PROGRAM_POINT_SIZE);
  glEnable(GL_BLEND);
  if (blend == PointBlendMode::Additive) {
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
  } else {
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }

  shader_.bind();
  shader_.setUniformMat4("uView", view_);
  shader_.setUniformMat4("uProj", proj_);
  shader_.setUniform1i("uUseTex", 0);

  gl::BindVertexArray(vao_);
  gl::BindBuffer(GL_ARRAY_BUFFER, vbo_);
  gl::BufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(points.size() * sizeof(PointVertex)),
                 points.data(),
                 GL_DYNAMIC_DRAW);

  glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(points.size()));

  gl::BindVertexArray(0);
}

void PointRenderer::drawPointsSprite(const std::vector<PointVertex>& points,
                                    const Texture2D& sprite,
                                    PointBlendMode blend) {
  if (points.empty()) return;

  glEnable(GL_PROGRAM_POINT_SIZE);
  glEnable(GL_BLEND);
  if (blend == PointBlendMode::Additive) {
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
  } else {
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }

  shader_.bind();
  shader_.setUniformMat4("uView", view_);
  shader_.setUniformMat4("uProj", proj_);
  shader_.setUniform1i("uUseTex", 1);
  shader_.setUniform1i("uTex", 0);

  sprite.bind(0);

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
