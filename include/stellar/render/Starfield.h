#pragma once

#include "stellar/core/Random.h"
#include "stellar/math/Vec3.h"
#include "stellar/render/PointRenderer.h"

#include <cstddef>
#include <vector>

namespace stellar::render {

// A simple deterministic procedural starfield meant as a cheap but effective
// background layer.
//
// The star positions are stored as direction vectors on the unit sphere and are
// expanded into world positions around the camera each frame.
class Starfield {
public:
  Starfield() = default;

  // Star distribution modes.
  enum class Distribution {
    // Legacy behavior: purely random uniform sampling on a sphere.
    UniformRandom = 0,

    // Deterministic quasi-uniform sampling using a Fibonacci lattice.
    Fibonacci,

    // Fibonacci lattice followed by a fast blue-noise style relaxation pass.
    // This reduces visible clumping/banding without an expensive Poisson solver.
    RelaxedFibonacci,
  };

  // Controls star placement + subtle appearance properties.
  //
  // NOTE: These settings are consumed by regenerate(). Starfield does not
  // automatically re-seed/regenerate when settings change.
  struct Settings {
    Distribution distribution{Distribution::RelaxedFibonacci};

    // 1 = uniform. Values >1 squeeze density toward the equator (a mild
    // "Milky Way" band). Values <1 bias toward the poles.
    double bandPower{1.0};

    // Apply a deterministic random rotation to hide lattice artifacts.
    bool randomRotation{true};

    // Small tangent-plane jitter (fraction of expected spacing). 0 disables.
    // Useful for hiding residual lattice patterns in Fibonacci mode.
    double jitter01{0.08};

    // Blue-noise relaxation parameters (only used in RelaxedFibonacci).
    int relaxIterations{4};

    // 0..1-ish step strength. Higher values converge faster but can overshoot.
    double relaxStrength{0.20};
  };

  // GPU-friendly packed star attributes.
  //
  // When rendering via StarfieldGpuRenderer, we upload these once after
  // regenerate() and then compute the final world position + twinkle on the GPU
  // (vertex shader) each frame. This avoids per-frame CPU expansion of stars
  // around the camera and dramatically reduces streaming bandwidth.
  struct GpuStar {
    // Unit direction vector on the sky sphere.
    float dx{0.0f}, dy{0.0f}, dz{1.0f};

    // Base color.
    float r{1.0f}, g{1.0f}, b{1.0f};

    float baseAlpha{1.0f};
    float sizePx{1.0f};
    float twinkleSpeed{1.0f};
    float phase{0.0f};
  };

  void setRadius(double radiusU) { radiusU_ = radiusU; }
  double radius() const { return radiusU_; }

  const Settings& settings() const { return settings_; }
  void setSettings(const Settings& s) { settings_ = s; }

  // Regenerate star distribution using current settings().
  void regenerate(core::u64 seed, int starCount);

  // Convenience: regenerate while temporarily overriding settings (also stores it).
  void regenerate(core::u64 seed, int starCount, const Settings& settings);

  // Rebuild the renderable point list for a given camera position.
  void update(const math::Vec3d& cameraPosU, double timeSeconds);

  const std::vector<PointVertex>& points() const { return points_; }
  std::size_t starCount() const { return stars_.size(); }

  // Export GPU-friendly packed star attributes.
  // Intended for use with StarfieldGpuRenderer.
  void exportGpuStars(std::vector<GpuStar>& out) const;

private:
  struct Star {
    math::Vec3d dir{0, 0, 1};
    float r{1}, g{1}, b{1};
    float baseAlpha{1.0f};
    float sizePx{1.0f};
    float twinkleSpeed{1.0f};
    float phase{0.0f};
  };

  double radiusU_{16000.0};
  core::u64 seed_{0};
  Settings settings_{};

  std::vector<Star> stars_;
  std::vector<PointVertex> points_;
};

} // namespace stellar::render
