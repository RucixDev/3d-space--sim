#pragma once

#include "stellar/render/Shader.h"

#include <cstdint>
#include <string>

namespace stellar::render {

// A GPU-resident particle field simulated via render-to-texture ping-pong.
//
// Notes:
//  - Requires OpenGL 3.3+.
//  - Simulation state lives in two RGBA32F textures for position and velocity
//    (each ping-ponged every step).
//  - Rendering uses GL_POINTS and gl_VertexID to fetch particle state with
//    texelFetch (no VBO uploads).
class GpuParticleField {
public:
  struct Settings {
    bool enabled{false};
    bool textured{true};

    // Simulation texture dimension (particle count = dim * dim).
    // Recommended: 64..512 (4096..262k particles).
    int dim{256};

    // Seed for deterministic init.
    int seed{1};

    // Simulation bounds in render units (a cube from -bounds..+bounds).
    float boundsU{45.0f};

    // Velocity damping (0..). Roughly: v *= exp(-drag * dt).
    float drag{0.18f};

    // Noise-driven acceleration field.
    float noiseFrequency{0.07f};
    float noiseStrength{1.10f};

    // Optional pull toward origin.
    float attractStrength{0.0f};

    // Cap on |v| (render units / second).
    float maxSpeed{18.0f};

    // Rendering
    float pointSizePx{2.25f};
    float intensity{1.0f};
    float alpha{0.65f};
  };

  GpuParticleField() = default;
  ~GpuParticleField();

  GpuParticleField(const GpuParticleField&) = delete;
  GpuParticleField& operator=(const GpuParticleField&) = delete;

  bool init(std::string* outError = nullptr);
  bool isInited() const { return inited_; }

  void setViewProj(const float* view, const float* proj);

  // Force a deterministic re-seed on the next update() based on Settings.
  void reset();

  // Step the simulation.
  // shiftDeltaU is a camera/floating-origin translation to apply to all
  // particle positions this frame (newPos = oldPos - shiftDeltaU).
  void update(double dtSec,
              const float shiftDeltaU[3],
              float timeSec,
              const Settings& s);

  // Draw particles. If spriteTexHandle==0 or Settings.textured==false,
  // a procedural soft circle is used instead of sampling a sprite.
  void draw(unsigned int spriteTexHandle, const Settings& s);

  int dim() const { return dim_; }
  std::int64_t particleCount() const { return (std::int64_t)dim_ * (std::int64_t)dim_; }

private:
  void destroyGL();
  bool ensureSize(int dim, std::string* outError);
  bool initState(const Settings& s, std::string* outError);

  // Helpers for ping-pong indices.
  int readIndex() const { return ping_ & 1; }
  int writeIndex() const { return (ping_ ^ 1) & 1; }

  bool inited_{false};
  bool needsInit_{true};

  int dim_{0};
  int lastSeed_{0};
  int ping_{0};

  unsigned int vao_{0};
  unsigned int fbo_{0};
  unsigned int posTex_[2]{0, 0};
  unsigned int velTex_[2]{0, 0};

  ShaderProgram initShader_{};
  ShaderProgram velShader_{};
  ShaderProgram posShader_{};
  ShaderProgram drawShader_{};

  float view_[16]{};
  float proj_[16]{};
};

} // namespace stellar::render
