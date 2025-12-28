#pragma once

#include "stellar/core/Random.h"
#include "stellar/math/Vec3.h"
#include "stellar/render/PointRenderer.h"

#include <cstddef>
#include <vector>

namespace stellar::render {

// A lightweight CPU particle system intended for small, high-impact VFX:
// - ship thruster plumes
// - impact sparks
// - small explosions
//
// It renders particles as point sprites (see PointRenderer). The goal is to
// provide a fun "feel" upgrade without pulling in a full GPU particle pipeline.
class ParticleSystem {
public:
  ParticleSystem() = default;

  void reseed(core::u64 seed) { rng_.reseed(seed); }

  void setMaxParticles(std::size_t max) { maxParticles_ = (max == 0) ? 1 : max; }
  std::size_t maxParticles() const { return maxParticles_; }
  std::size_t aliveCount() const { return particles_.size(); }

  void clear() { particles_.clear(); emitCarry_ = 0.0; }

  // Advance particle simulation by dt seconds.
  void update(double dtSeconds);

  // Append the current particles as PointVertex entries to be fed to PointRenderer.
  // The output vector is cleared and then filled.
  void buildPoints(std::vector<PointVertex>& out) const;

  // ---- Spawn helpers ----

  // Continuous-ish thruster plume emitter. Designed to be called once per frame.
  // - posU: emitter position in render units
  // - dirWorld: direction the exhaust travels (unit vec), in world space
  // - intensity: 0..1
  void emitThruster(const math::Vec3d& posU,
                    const math::Vec3d& dirWorld,
                    double intensity,
                    double dtSeconds,
                    bool boost = false);

  // A short, bright burst.
  void spawnExplosion(const math::Vec3d& posU,
                      const math::Vec3d& inheritVelU,
                      double energy,
                      int count = -1);

  // Small impact sparks that scatter along a hemisphere oriented around normal.
  void spawnSparks(const math::Vec3d& posU,
                   const math::Vec3d& normalWorld,
                   const math::Vec3d& inheritVelU,
                   double energy,
                   int count = -1);

private:
  struct Particle {
    math::Vec3d posU{0,0,0};
    math::Vec3d velU{0,0,0};
    float r{1}, g{1}, b{1};
    float baseAlpha{1.0f};
    float sizePx{4.0f};
    float life{1.0f};
    float maxLife{1.0f};
    float drag{0.0f};
  };

  void push(const Particle& p);

  std::vector<Particle> particles_;
  std::size_t maxParticles_{6000};
  core::SplitMix64 rng_{0};

  // Fractional emission accumulator used by emitThruster.
  double emitCarry_{0.0};
};

} // namespace stellar::render
