#pragma once

#include "stellar/core/Random.h"
#include "stellar/math/Vec3.h"
#include "stellar/render/PointRenderer.h"

#include <vector>

namespace stellar::render {

// A cheap nebula/background gas layer built from large additive point sprites.
//
// This complements the Starfield (tiny points) by adding slow-moving, colored
// "cloud puffs" with optional banding (galactic plane) and mild parallax.
class NebulaField {
public:
  struct Settings {
    // Distribution of puffs in a spherical shell around an anchor point.
    double innerRadiusU{9000.0};
    double outerRadiusU{22000.0};

    // 0 -> fixed around origin (full parallax), 1 -> follows camera (no parallax).
    double parallax{0.25};

    // Overall visual tuning.
    float intensity{1.0f}; // multiplies RGB (HDR-friendly)
    float opacity{0.16f};  // base alpha
    float sizeMinPx{70.0f};
    float sizeMaxPx{240.0f};

    // Animation (0 disables)
    float turbulence{0.35f};       // alpha wobble amplitude
    float turbulenceSpeed{0.10f};  // Hz-ish
  };

  NebulaField() = default;

  // Regenerate the puff field.
  // bandPower > 1 squeezes puffs toward y=0 (galactic plane band).
  void regenerate(core::u64 seed, int puffCount, float bandPower = 1.8f);

  // Rebuild the renderable point list for a given camera position.
  void update(const math::Vec3d& cameraPosU, double timeSeconds, const Settings& s);

  const std::vector<PointVertex>& points() const { return points_; }
  std::size_t puffCount() const { return puffs_.size(); }

private:
  struct Puff {
    math::Vec3d dir{0,0,1};
    float r{1}, g{1}, b{1};
    float alpha{1.0f};
    float size01{0.5f};   // random in [0,1]
    float depth01{0.5f};  // random in [0,1] (maps to inner/outer radius)
    float twinkleSpeed{1.0f};
    float phase{0.0f};
  };

  core::u64 seed_{0};
  float bandPower_{1.8f};

  std::vector<Puff> puffs_;
  std::vector<PointVertex> points_;
};

} // namespace stellar::render
