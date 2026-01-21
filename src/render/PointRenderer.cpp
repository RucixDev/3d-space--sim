#include "stellar/render/PointRenderer.h"

#include "stellar/render/Gl.h"

#include <cstring>

// Some platform/header combinations (notably Windows) ship legacy OpenGL
// headers that omit these tokens even though MapBufferRange is available.
#ifndef GL_MAP_WRITE_BIT
#define GL_MAP_WRITE_BIT 0x0002
#endif
#ifndef GL_MAP_INVALIDATE_RANGE_BIT
#define GL_MAP_INVALIDATE_RANGE_BIT 0x0004
#endif
#ifndef GL_MAP_INVALIDATE_BUFFER_BIT
#define GL_MAP_INVALIDATE_BUFFER_BIT 0x0008
#endif

namespace stellar::render {

PointRenderer::~PointRenderer() {
  for (int i = 0; i < kBufferRing; ++i) {
    if (vbo_[i]) gl::DeleteBuffers(1, &vbo_[i]);
    if (vao_[i]) gl::DeleteVertexArrays(1, &vao_[i]);
    vbo_[i] = vao_[i] = 0;
    vboCapacityBytes_[i] = 0;
  }
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

  // Create a small VAO/VBO ring so updates never stomp the buffer currently
  // in-flight on the GPU. This reduces implicit synchronization stalls on
  // some drivers when streaming many point sprites.
  for (int i = 0; i < kBufferRing; ++i) {
    gl::GenVertexArrays(1, &vao_[i]);
    gl::GenBuffers(1, &vbo_[i]);

    gl::BindVertexArray(vao_[i]);
    gl::BindBuffer(GL_ARRAY_BUFFER, vbo_[i]);
    gl::BufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
    vboCapacityBytes_[i] = 0;

    gl::EnableVertexAttribArray(0);
    gl::VertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(PointVertex), (void*)offsetof(PointVertex, px));
    gl::EnableVertexAttribArray(1);
    gl::VertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(PointVertex), (void*)offsetof(PointVertex, cr));
    gl::EnableVertexAttribArray(2);
    gl::VertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(PointVertex), (void*)offsetof(PointVertex, a));
    gl::EnableVertexAttribArray(3);
    gl::VertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(PointVertex), (void*)offsetof(PointVertex, size));
  }

  gl::BindVertexArray(0);
  ringIndex_ = 0;

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

  ringIndex_ = (ringIndex_ + 1) % kBufferRing;
  gl::BindVertexArray(vao_[ringIndex_]);
  gl::BindBuffer(GL_ARRAY_BUFFER, vbo_[ringIndex_]);

  const std::size_t bytes = points.size() * sizeof(PointVertex);
  if (bytes > vboCapacityBytes_[ringIndex_]) {
    gl::BufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(bytes), nullptr, GL_STREAM_DRAW);
    vboCapacityBytes_[ringIndex_] = bytes;
  }

  // Update buffer contents with a streaming-friendly path.
  // Prefer MapBufferRange (with invalidation) to avoid implicit reallocation.
  bool wrote = false;
  if (gl::MapBufferRange) {
    void* dst = gl::MapBufferRange(GL_ARRAY_BUFFER,
                                  0,
                                  static_cast<GLsizeiptr>(bytes),
                                  GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT);
    if (dst) {
      std::memcpy(dst, points.data(), bytes);
      if (gl::UnmapBuffer) {
        (void)gl::UnmapBuffer(GL_ARRAY_BUFFER);
      }
      wrote = true;
    }
  }

  if (!wrote) {
    gl::BufferSubData(GL_ARRAY_BUFFER,
                      0,
                      static_cast<GLsizeiptr>(bytes),
                      points.data());
  }

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

  ringIndex_ = (ringIndex_ + 1) % kBufferRing;
  gl::BindVertexArray(vao_[ringIndex_]);
  gl::BindBuffer(GL_ARRAY_BUFFER, vbo_[ringIndex_]);

  const std::size_t bytes = points.size() * sizeof(PointVertex);
  if (bytes > vboCapacityBytes_[ringIndex_]) {
    gl::BufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(bytes), nullptr, GL_STREAM_DRAW);
    vboCapacityBytes_[ringIndex_] = bytes;
  }

  bool wrote = false;
  if (gl::MapBufferRange) {
    void* dst = gl::MapBufferRange(GL_ARRAY_BUFFER,
                                  0,
                                  static_cast<GLsizeiptr>(bytes),
                                  GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT);
    if (dst) {
      std::memcpy(dst, points.data(), bytes);
      if (gl::UnmapBuffer) {
        (void)gl::UnmapBuffer(GL_ARRAY_BUFFER);
      }
      wrote = true;
    }
  }

  if (!wrote) {
    gl::BufferSubData(GL_ARRAY_BUFFER,
                      0,
                      static_cast<GLsizeiptr>(bytes),
                      points.data());
  }

  glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(points.size()));

  gl::BindVertexArray(0);
}

} // namespace stellar::render
