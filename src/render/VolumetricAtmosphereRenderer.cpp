#include "stellar/render/VolumetricAtmosphereRenderer.h"

#include "stellar/render/Gl.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace stellar::render {

VolumetricAtmosphereRenderer::~VolumetricAtmosphereRenderer() {
  if (instanceVbo_) {
    gl::DeleteBuffers(1, &instanceVbo_);
    instanceVbo_ = 0;
  }
}

static const char* kVS = R"GLSL(
#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aNrm;

layout(location=3) in vec3 iPos;
layout(location=4) in vec3 iScale;
layout(location=5) in vec4 iQuat;
layout(location=6) in vec3 iColor;
layout(location=7) in float iInnerRadius;
layout(location=8) in float iScaleHeight;
layout(location=9) in float iDensityMul;
layout(location=10) in float iMieMul;

uniform mat4 uView;
uniform mat4 uProj;

out vec3 vWorldPos;
out vec3 vNrm;
out vec3 vColor;
out vec3 vCenter;
out float vOuterRadius;
out float vInnerRadius;
out float vScaleHeight;
out float vDensityMul;
out float vMieMul;

vec3 quatRotate(vec4 q, vec3 v) {
  vec3 t = 2.0 * cross(q.xyz, v);
  return v + q.w * t + cross(q.xyz, t);
}

void main() {
  vec3 local = aPos * iScale;
  vec3 pos = quatRotate(iQuat, local) + iPos;
  gl_Position = uProj * uView * vec4(pos, 1.0);

  vWorldPos = pos;
  vNrm = quatRotate(iQuat, aNrm);
  vColor = iColor;
  vCenter = iPos;

  // We expect near-uniform scale for spheres; take X as radius.
  vOuterRadius = iScale.x;
  vInnerRadius = iInnerRadius;
  vScaleHeight = iScaleHeight;
  vDensityMul = iDensityMul;
  vMieMul = iMieMul;
}
)GLSL";

static const char* kFS = R"GLSL(
#version 330 core

in vec3 vWorldPos;
in vec3 vNrm;
in vec3 vColor;
in vec3 vCenter;
in float vOuterRadius;
in float vInnerRadius;
in float vScaleHeight;
in float vDensityMul;
in float vMieMul;

uniform vec3 uCamPos;
uniform vec3 uSunPos;

uniform float uIntensity;
uniform float uRayleighStrength;
uniform float uMieStrength;
uniform float uScaleHeightMul;
uniform float uDensityMul;

uniform int uViewSteps;
uniform int uSunSteps;
uniform float uDitherStrength;
uniform float uMieG;

uniform int uUsePhaseLut;
uniform sampler2D uPhaseLut;
uniform float uMiePhaseStrength;

// Analytic multiple scattering approximation.
uniform float uMsStrength;
uniform float uMsAlbedo;

uniform int uUseMsPhaseLut;
uniform sampler2D uMsPhaseLut;
uniform float uMsPhaseStrength;

out vec4 FragColor;

// --- Phase functions ---

// Rayleigh phase function (normalized over 4π): 3/(16π) * (1 + μ^2)
float phaseRayleigh(float mu) {
  const float k = 0.05968310365946075; // 3/(16*pi)
  return k * (1.0 + mu * mu);
}

// Henyey-Greenstein phase (normalized over 4π)
float phaseHG(float mu, float g) {
  g = clamp(g, -0.99, 0.99);
  float g2 = g * g;
  float denom = pow(max(1e-4, 1.0 + g2 - 2.0 * g * mu), 1.5);
  return (1.0 / (4.0 * 3.14159265)) * (1.0 - g2) / denom;
}

// Isotropic phase (normalized over 4π)
float phaseIsotropic() {
  return 1.0 / (4.0 * 3.14159265);
}

// --- Helpers ---

bool raySphere(vec3 ro, vec3 rd, vec3 c, float r, out float t0, out float t1) {
  vec3 oc = ro - c;
  float b = dot(oc, rd);
  float cc = dot(oc, oc) - r * r;
  float h = b * b - cc;
  if (h < 0.0) return false;
  float s = sqrt(h);
  t0 = -b - s;
  t1 = -b + s;
  return true;
}

float hash12(vec2 p) {
  // Cheap, deterministic hash for dithering.
  vec3 p3 = fract(vec3(p.xyx) * 0.1031);
  p3 += dot(p3, p3.yzx + 33.33);
  return fract((p3.x + p3.y) * p3.z);
}

float densityAt(vec3 p, vec3 center, float innerR, float scaleH) {
  float h = length(p - center) - innerR;
  h = max(0.0, h);
  scaleH = max(1e-5, scaleH);
  return exp(-h / scaleH);
}

float integrateOpticalDepth(vec3 startPos,
                            vec3 dir,
                            float tMax,
                            vec3 center,
                            float innerR,
                            float scaleH,
                            int steps) {
  // Uniform steps along [0, tMax].
  steps = clamp(steps, 2, 16);
  float dt = tMax / float(steps);
  float t = dt * 0.5;
  float od = 0.0;
  for (int i = 0; i < 16; ++i) {
    if (i >= steps) break;
    vec3 p = startPos + dir * t;
    float d = densityAt(p, center, innerR, scaleH);
    od += d * dt;
    t += dt;
  }
  return od;
}

void main() {
  // Compute the view ray for this fragment.
  vec3 ro = uCamPos;
  vec3 rd = normalize(vWorldPos - uCamPos);
  vec3 center = vCenter;

  float outerR = max(1e-4, vOuterRadius);
  float innerR = max(1e-4, vInnerRadius);

  // Find ray segment through the atmosphere shell.
  float t0, t1;
  if (!raySphere(ro, rd, center, outerR, t0, t1)) discard;
  if (t1 <= 0.0) discard;
  t0 = max(t0, 0.0);

  // Clamp to the planet surface if the ray hits the inner sphere.
  float ti0, ti1;
  if (innerR > 1e-4 && raySphere(ro, rd, center, innerR, ti0, ti1)) {
    if (ti0 > t0 && ti0 < t1) {
      t1 = ti0;
    }
  }
  if (t1 <= t0) discard;

  // Sampling controls.
  int viewSteps = clamp(uViewSteps, 2, 32);
  int sunSteps = clamp(uSunSteps, 2, 16);

  float scaleH = max(1e-5, vScaleHeight * max(0.01, uScaleHeightMul));
  float densMul = max(0.0, vDensityMul) * max(0.0, uDensityMul);

  // Wavelength dependence (relative): Rayleigh ~ λ^-4.
  // Use approximate visible wavelengths in micrometers to keep magnitudes stable.
  const vec3 lambda = vec3(0.680, 0.550, 0.440);
  vec3 invLambda4 = 1.0 / (lambda * lambda * lambda * lambda);

  // Scattering / extinction coefficients. These are in "game units"; tune via sliders.
  vec3 betaR = invLambda4 * (uRayleighStrength * densMul);
  vec3 betaM = vec3(uMieStrength * densMul * max(0.0, vMieMul));
  vec3 betaExt = betaR + betaM;

  // Dithered integration start (reduces banding).
  float jitter = (hash12(gl_FragCoord.xy) - 0.5) * clamp(uDitherStrength, 0.0, 1.0);
  float seg = (t1 - t0) / float(viewSteps);
  float t = t0 + seg * (0.5 + jitter);

  vec3 sum = vec3(0.0);
  float odView = 0.0;

  for (int i = 0; i < 32; ++i) {
    if (i >= viewSteps) break;

    vec3 p = ro + rd * t;
    // Local density sample (dimensionless). Overall thickness is controlled via
    // the scattering/extinction coefficients (densMul applied in betaR/betaM).
    float d = densityAt(p, center, innerR, scaleH);

    // Accumulate view optical depth up to this sample.
    odView += d * seg;

    // Sun ray from sample to the edge of atmosphere.
    vec3 l = normalize(uSunPos - p);
    float ts0, ts1;
    float odSun = 0.0;
    if (raySphere(p, l, center, outerR, ts0, ts1)) {
      float tExit = max(0.0, ts1);
      odSun = integrateOpticalDepth(p, l, tExit, center, innerR, scaleH, sunSteps);
    }

    // Beer-Lambert transmittance.
    vec3 transSun = exp(-betaExt * odSun);
    vec3 transView = exp(-betaExt * odView);
    vec3 trans = transSun * transView;

    // Scattering angle cosine μ between incoming light (-l) and view direction.
    vec3 vDir = normalize(uCamPos - p);
    float mu = clamp(dot(vDir, -l), -1.0, 1.0);

    float pr = phaseRayleigh(mu);

    vec3 pm = vec3(phaseHG(mu, uMieG));
    if (uUsePhaseLut != 0) {
      float tt = mu * 0.5 + 0.5;
      vec3 lut = texture(uPhaseLut, vec2(tt, 0.5)).rgb;
      pm = mix(pm, lut, clamp(uMiePhaseStrength, 0.0, 1.0));
    }

    // Single scattering contribution.
    vec3 scatter1 = (betaR * pr) + (betaM * pm);
    sum += d * seg * trans * scatter1;

    // --- Analytic multiple scattering (cheap) ---
    // We approximate higher orders with a geometric-series energy boost and a
    // broader phase function. This is not a full RT solution, but it adds the
    // missing "milkiness"/haze and a better horizon rolloff.
    float msStrength = max(0.0, uMsStrength);
    if (msStrength > 1e-5) {
      float betaExtAvg = (betaExt.r + betaExt.g + betaExt.b) * 0.3333333;
      float tauSun = betaExtAvg * odSun;
      float f = 1.0 - exp(-tauSun); // how much direct sun is removed on the way in

      float w = clamp(uMsAlbedo, 0.0, 0.99);
      float denom = max(0.02, 1.0 - w * f);
      float series = (w * f) / denom;
      series = min(series, 50.0);

      // Multiple scattering tends toward isotropic; Rayleigh is already broad.
      float prMs = phaseIsotropic();
      float gMs = clamp(uMieG * uMieG, 0.0, 0.99);
      vec3 pmMs = vec3(phaseHG(mu, gMs));
      if (uUseMsPhaseLut != 0) {
        float tt = mu * 0.5 + 0.5;
        vec3 lut = texture(uMsPhaseLut, vec2(tt, 0.5)).rgb;
        pmMs = mix(pmMs, lut, clamp(uMsPhaseStrength, 0.0, 1.0));
      }

      vec3 scatterMs = (betaR * prMs) + (betaM * pmMs);
      sum += d * seg * transView * scatterMs * (msStrength * series);
    }

    // Early out if we're basically fully attenuated.
    if (max(max(trans.r, trans.g), trans.b) < 1e-4) break;

    t += seg;
  }

  // Apply tint and intensity. vColor acts as an artistic bias.
  vec3 outRgb = sum * vColor * max(0.0, uIntensity);

  // Additive blend path expects alpha=1 (SRC_ALPHA, ONE).
  FragColor = vec4(outRgb, 1.0);
}
)GLSL";

bool VolumetricAtmosphereRenderer::init(std::string* outError) {
  if (!shader_.build(kVS, kFS, outError)) return false;

  gl::GenBuffers(1, &instanceVbo_);

  // Default identity matrices.
  for (int i = 0; i < 16; ++i) {
    view_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
    proj_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
  }

  return true;
}

void VolumetricAtmosphereRenderer::setViewProj(const float* view, const float* proj) {
  std::memcpy(view_, view, sizeof(float) * 16);
  std::memcpy(proj_, proj, sizeof(float) * 16);
}

void VolumetricAtmosphereRenderer::drawInstances(const std::vector<VolumetricAtmosphereInstance>& instances) {
  if (!mesh_ || instances.empty()) return;

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE);

  const bool hasLut = useMiePhaseLut_ && miePhaseLut_ && miePhaseLut_->handle() != 0;
  const bool hasMsLut = useMsPhaseLut_ && msPhaseLut_ && msPhaseLut_->handle() != 0;

  shader_.bind();
  shader_.setUniformMat4("uView", view_);
  shader_.setUniformMat4("uProj", proj_);
  shader_.setUniform3f("uCamPos", camPos_[0], camPos_[1], camPos_[2]);
  shader_.setUniform3f("uSunPos", sunPos_[0], sunPos_[1], sunPos_[2]);

  shader_.setUniform1f("uIntensity", intensity_);
  shader_.setUniform1f("uRayleighStrength", rayleighStrength_);
  shader_.setUniform1f("uMieStrength", mieStrength_);
  shader_.setUniform1f("uScaleHeightMul", scaleHeightMul_);
  shader_.setUniform1f("uDensityMul", densityMul_);
  shader_.setUniform1i("uViewSteps", viewSteps_);
  shader_.setUniform1i("uSunSteps", sunSteps_);
  shader_.setUniform1f("uDitherStrength", ditherStrength_);
  shader_.setUniform1f("uMieG", mieG_);

  // Multiple scattering approximation.
  shader_.setUniform1f("uMsStrength", msStrength_);
  shader_.setUniform1f("uMsAlbedo", msAlbedo_);
  shader_.setUniform1i("uUseMsPhaseLut", hasMsLut ? 1 : 0);
  shader_.setUniform1f("uMsPhaseStrength", msPhaseStrength_);

  shader_.setUniform1i("uUsePhaseLut", hasLut ? 1 : 0);
  shader_.setUniform1f("uMiePhaseStrength", miePhaseStrength_);
  if (hasLut) {
    constexpr int kPhaseUnit = 3;
    shader_.setUniform1i("uPhaseLut", kPhaseUnit);
    miePhaseLut_->bind(kPhaseUnit);
  }

  if (hasMsLut) {
    constexpr int kMsPhaseUnit = 4;
    shader_.setUniform1i("uMsPhaseLut", kMsPhaseUnit);
    msPhaseLut_->bind(kMsPhaseUnit);
  }

  mesh_->bind();

  gl::BindBuffer(GL_ARRAY_BUFFER, instanceVbo_);
  gl::BufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(instances.size() * sizeof(VolumetricAtmosphereInstance)),
                 instances.data(),
                 GL_DYNAMIC_DRAW);

  // Location 3..6 mirror MeshRenderer / AtmosphereRenderer.
  gl::EnableVertexAttribArray(3);
  gl::VertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(VolumetricAtmosphereInstance),
                          (void*)offsetof(VolumetricAtmosphereInstance, px));
  gl::VertexAttribDivisor(3, 1);

  gl::EnableVertexAttribArray(4);
  gl::VertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(VolumetricAtmosphereInstance),
                          (void*)offsetof(VolumetricAtmosphereInstance, sx));
  gl::VertexAttribDivisor(4, 1);

  gl::EnableVertexAttribArray(5);
  gl::VertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, sizeof(VolumetricAtmosphereInstance),
                          (void*)offsetof(VolumetricAtmosphereInstance, qx));
  gl::VertexAttribDivisor(5, 1);

  gl::EnableVertexAttribArray(6);
  gl::VertexAttribPointer(6, 3, GL_FLOAT, GL_FALSE, sizeof(VolumetricAtmosphereInstance),
                          (void*)offsetof(VolumetricAtmosphereInstance, cr));
  gl::VertexAttribDivisor(6, 1);

  // Extra per-instance parameters.
  gl::EnableVertexAttribArray(7);
  gl::VertexAttribPointer(7, 1, GL_FLOAT, GL_FALSE, sizeof(VolumetricAtmosphereInstance),
                          (void*)offsetof(VolumetricAtmosphereInstance, innerRadius));
  gl::VertexAttribDivisor(7, 1);

  gl::EnableVertexAttribArray(8);
  gl::VertexAttribPointer(8, 1, GL_FLOAT, GL_FALSE, sizeof(VolumetricAtmosphereInstance),
                          (void*)offsetof(VolumetricAtmosphereInstance, scaleHeight));
  gl::VertexAttribDivisor(8, 1);

  gl::EnableVertexAttribArray(9);
  gl::VertexAttribPointer(9, 1, GL_FLOAT, GL_FALSE, sizeof(VolumetricAtmosphereInstance),
                          (void*)offsetof(VolumetricAtmosphereInstance, densityMul));
  gl::VertexAttribDivisor(9, 1);

  gl::EnableVertexAttribArray(10);
  gl::VertexAttribPointer(10, 1, GL_FLOAT, GL_FALSE, sizeof(VolumetricAtmosphereInstance),
                          (void*)offsetof(VolumetricAtmosphereInstance, mieMul));
  gl::VertexAttribDivisor(10, 1);

  mesh_->drawInstanced(static_cast<std::uint32_t>(instances.size()));
}

} // namespace stellar::render
