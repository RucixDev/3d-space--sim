#pragma once

#include "stellar/core/Random.h"
#include "stellar/math/Vec3.h"
#include "stellar/render/PointRenderer.h"

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

  void setRadius(double radiusU) { radiusU_ = radiusU; }
  double radius() const { return radiusU_; }

  void regenerate(core::u64 seed, int starCount);

  // Rebuild the renderable point list for a given camera position.
  void update(const math::Vec3d& cameraPosU, double timeSeconds);

  const std::vector<PointVertex>& points() const { return points_; }
  std::size_t starCount() const { return stars_.size(); }

private:
  struct Star {
    math::Vec3d dir{0,0,1};
    float r{1}, g{1}, b{1};
    float baseAlpha{1.0f};
    float sizePx{1.0f};
    float twinkleSpeed{1.0f};
    float phase{0.0f};
  };

  double radiusU_{16000.0};
  core::u64 seed_{0};

  std::vector<Star> stars_;
  std::vector<PointVertex> points_;
};

} // namespace stellar::render
