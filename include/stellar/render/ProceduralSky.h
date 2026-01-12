#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/Shader.h"

#include <string>

namespace stellar::render {

// GPU procedural background sky renderer.
//
// Renders a full-screen triangle with a fragment shader that synthesizes:
//  - stars (grid-jittered, seeded)
//  - a cheap volumetric-ish nebula (few-step raymarch through seeded 3D value noise)
//
// Intended usage (caller controlled):
//  - call init() once after GL context creation
//  - each frame: with depth writes disabled (and preferably depth test disabled), call draw().
//
// The output is HDR-friendly and is meant to feed into PostFX tonemapping/bloom.
struct ProceduralSkySettings {
  // A 32-bit seed used inside the shader. (The game typically hashes a 64-bit world seed.)
  int seed{1};

  // Overall intensity multiplier (HDR).
  float intensity{1.0f};

  // ---- Stars ----
  float starIntensity{1.20f};
  // Approximate cell density in equirectangular UV space.
  float starDensity{900.0f};
  // Probability a given cell contains a star (0..1).
  float starProbability{0.055f};
  // Scales star radius.
  float starSize{1.0f};
  // 0 disables time-varying twinkle.
  float starTwinkle{0.20f};

  // ---- Nebula ----
  bool nebulaEnabled{true};
  float nebulaIntensity{1.15f};
  float nebulaFrequency{2.40f};
  float nebulaThreshold{0.48f};
  float nebulaSoftness{0.18f};
  float nebulaBandPower{1.8f};
  // 0 -> fixed in world (full parallax), 1 -> camera-locked (no parallax)
  float nebulaParallax{0.20f};
  // Number of raymarch steps (clamped in shader). Higher = smoother but slower.
  int nebulaSteps{10};
  // Scale applied to camera position before sampling noise (keeps numbers small).
  float nebulaWorldScale{1.0f / 12000.0f};
  // Time-based drift speed (0 disables).
  float nebulaDrift{0.03f};
};

class ProceduralSky {
public:
  ProceduralSky() = default;
  ~ProceduralSky();

  ProceduralSky(const ProceduralSky&) = delete;
  ProceduralSky& operator=(const ProceduralSky&) = delete;

  bool init(std::string* outError = nullptr);

  // Draw a procedural sky background.
  //
  // Camera basis vectors are in world/render space:
  //  - camRight:  (1,0,0) rotated by camera orientation
  //  - camUp:     (0,1,0) rotated by camera orientation
  //  - camForward:(0,0,1) rotated by camera orientation
  //
  // tanHalfFovY = tan(fovY/2)
  void draw(int w,
            int h,
            const ProceduralSkySettings& s,
            const float camRight[3],
            const float camUp[3],
            const float camForward[3],
            const float camPosU[3],
            float tanHalfFovY,
            float aspect,
            float timeSeconds) const;

private:
  ShaderProgram shader_{};
  unsigned int vao_{0};
};

} // namespace stellar::render
