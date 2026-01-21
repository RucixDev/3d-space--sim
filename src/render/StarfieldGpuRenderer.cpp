#include "stellar/render/StarfieldGpuRenderer.h"

#include "stellar/render/Gl.h"

#include <algorithm>
#include <cstring>

namespace stellar::render {
namespace {

static constexpr const char* kStarfieldGpuVS = R"GLSL(
#version 330 core

layout(location=0) in vec3 aDir;
layout(location=1) in vec3 aColor;
layout(location=2) in float aBaseAlpha;
layout(location=3) in float aSize;
layout(location=4) in float aTwinkleSpeed;
layout(location=5) in float aPhase;

uniform mat4 uView;
uniform mat4 uProj;

uniform vec3 uCameraPos;
uniform float uRadius;
uniform float uTime;

out vec3 vColor;
out float vAlpha;

void main() {
  vec3 pos = uCameraPos + aDir * uRadius;

  gl_Position = uProj * uView * vec4(pos, 1.0);
  gl_PointSize = aSize;

  float tw = 0.85 + 0.15 * sin(uTime * aTwinkleSpeed + aPhase);
  vAlpha = clamp(aBaseAlpha * tw, 0.0, 1.0);
  vColor = aColor;
}
)GLSL";

static constexpr const char* kStarfieldGpuFS = R"GLSL(
#version 330 core

in vec3 vColor;
in float vAlpha;

uniform sampler2D uTex;
uniform int uUseTex;

out vec4 fragColor;

void main() {
  vec2 uv = gl_PointCoord;
  float a = vAlpha;

  if (uUseTex != 0) {
    vec4 t = texture(uTex, uv);
    a *= t.a;
  } else {
    // Soft circular mask.
    vec2 p = uv * 2.0 - 1.0;
    float r2 = dot(p, p);
    if (r2 > 1.0) discard;
    float edge = smoothstep(1.0, 0.0, r2);
    a *= edge;
  }

  fragColor = vec4(vColor, a);
}
)GLSL";

static void setIdentity(float* m) {
  std::fill(m, m + 16, 0.0f);
  m[0] = m[5] = m[10] = m[15] = 1.0f;
}

} // namespace

StarfieldGpuRenderer::~StarfieldGpuRenderer() {
  if (vbo_ != 0) {
    gl::DeleteBuffers(1, &vbo_);
    vbo_ = 0;
  }
  if (vao_ != 0) {
    gl::DeleteVertexArrays(1, &vao_);
    vao_ = 0;
  }
}

bool StarfieldGpuRenderer::init(std::string* outError) {
  if (!shader_.build(kStarfieldGpuVS, kStarfieldGpuFS, outError)) {
    return false;
  }

  setIdentity(view_);
  setIdentity(proj_);

  gl::GenVertexArrays(1, &vao_);
  gl::GenBuffers(1, &vbo_);

  gl::BindVertexArray(vao_);
  gl::BindBuffer(GL_ARRAY_BUFFER, vbo_);

  // Start with an empty buffer; we grow on first upload.
  gl::BufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);
  vboCapacityBytes_ = 0;

  const GLsizei stride = (GLsizei)sizeof(Starfield::GpuStar);

  // aDir
  gl::EnableVertexAttribArray(0);
  gl::VertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, (void*)offsetof(Starfield::GpuStar, dx));

  // aColor
  gl::EnableVertexAttribArray(1);
  gl::VertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, (void*)offsetof(Starfield::GpuStar, r));

  // aBaseAlpha
  gl::EnableVertexAttribArray(2);
  gl::VertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, stride, (void*)offsetof(Starfield::GpuStar, baseAlpha));

  // aSize
  gl::EnableVertexAttribArray(3);
  gl::VertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, stride, (void*)offsetof(Starfield::GpuStar, sizePx));

  // aTwinkleSpeed
  gl::EnableVertexAttribArray(4);
  gl::VertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, stride, (void*)offsetof(Starfield::GpuStar, twinkleSpeed));

  // aPhase
  gl::EnableVertexAttribArray(5);
  gl::VertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, stride, (void*)offsetof(Starfield::GpuStar, phase));

  gl::BindVertexArray(0);

  // One-time uniform setup.
  shader_.bind();
  shader_.setUniform1i("uUseTex", 0);
  shader_.setUniform1i("uTex", 0);

  return true;
}

void StarfieldGpuRenderer::setViewProj(const float* view, const float* proj) {
  if (view) std::memcpy(view_, view, sizeof(view_));
  if (proj) std::memcpy(proj_, proj, sizeof(proj_));
}

void StarfieldGpuRenderer::setCameraPos(float x, float y, float z) {
  camPos_[0] = x;
  camPos_[1] = y;
  camPos_[2] = z;
}

void StarfieldGpuRenderer::setRadius(float r) {
  radius_ = r;
}

void StarfieldGpuRenderer::setTimeSeconds(float t) {
  time_ = t;
}

void StarfieldGpuRenderer::ensureCapacityBytes(std::size_t bytes) {
  if (bytes <= vboCapacityBytes_) return;

  gl::BindBuffer(GL_ARRAY_BUFFER, vbo_);
  gl::BufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(bytes), nullptr, GL_STATIC_DRAW);
  vboCapacityBytes_ = bytes;
}

void StarfieldGpuRenderer::upload(const std::vector<Starfield::GpuStar>& stars) {
  starCount_ = stars.size();
  if (stars.empty()) return;

  const std::size_t bytes = stars.size() * sizeof(Starfield::GpuStar);
  ensureCapacityBytes(bytes);

  gl::BindBuffer(GL_ARRAY_BUFFER, vbo_);
  gl::BufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLsizeiptr>(bytes), stars.data());
}

void StarfieldGpuRenderer::upload(const Starfield& starfield) {
  std::vector<Starfield::GpuStar> tmp;
  starfield.exportGpuStars(tmp);
  upload(tmp);
  setRadius((float)starfield.radius());
}

void StarfieldGpuRenderer::draw(const Texture2D* spriteTex, PointBlendMode blend) {
  if (starCount_ == 0) return;

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
  shader_.setUniform3f("uCameraPos", camPos_[0], camPos_[1], camPos_[2]);
  shader_.setUniform1f("uRadius", radius_);
  shader_.setUniform1f("uTime", time_);

  if (spriteTex) {
    shader_.setUniform1i("uUseTex", 1);
    shader_.setUniform1i("uTex", 0);
    spriteTex->bind(0);
  } else {
    shader_.setUniform1i("uUseTex", 0);
  }

  gl::BindVertexArray(vao_);
  glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(starCount_));
  gl::BindVertexArray(0);
}

} // namespace stellar::render
