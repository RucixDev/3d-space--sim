#include "stellar/sim/Occlusion.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

void OcclusionFieldKm::clear() {
  occluders_.clear();
  bvh_.clear();
}

bool OcclusionFieldKm::empty() const {
  return occluders_.empty() || bvh_.empty();
}

double OcclusionFieldKm::clamp01_(double v) {
  return std::clamp(v, 0.0, 1.0);
}

bool OcclusionFieldKm::segmentHitsSphereKm_(const math::Vec3d& aKm,
                                            const math::Vec3d& bKm,
                                            const math::Vec3d& centerKm,
                                            double radiusKm) {
  if (!math::isFinite(aKm) || !math::isFinite(bKm) || !math::isFinite(centerKm)) return false;
  if (!std::isfinite(radiusKm) || radiusKm <= 0.0) return false;

  const math::Vec3d ab = bKm - aKm;
  const double abLenSq = ab.lengthSq();
  if (!(abLenSq > 1e-24)) {
    // Degenerate segment.
    return (aKm - centerKm).lengthSq() <= radiusKm * radiusKm;
  }

  const double t = std::clamp(math::dot(centerKm - aKm, ab) / abLenSq, 0.0, 1.0);
  const math::Vec3d p = aKm + ab * t;
  return (p - centerKm).lengthSq() <= radiusKm * radiusKm;
}

void OcclusionFieldKm::build(std::vector<SphereOccluderKm> occluders, std::size_t leafSize) {
  clear();

  occluders_.reserve(occluders.size());
  for (const auto& o : occluders) {
    if (!math::isFinite(o.centerKm)) continue;
    if (!std::isfinite(o.radiusKm) || o.radiusKm <= 0.0) continue;
    const double s = clamp01_(o.strength01);
    if (!(s > 0.0)) continue;

    SphereOccluderKm out = o;
    out.strength01 = s;
    occluders_.push_back(out);
  }

  if (occluders_.empty()) return;

  std::vector<math::BvhItem3d> items;
  items.reserve(occluders_.size());

  for (std::size_t i = 0; i < occluders_.size(); ++i) {
    const auto& o = occluders_[i];
    math::BvhItem3d it{};
    it.bounds = math::Aabb3d::fromCenterExtents(o.centerKm, {o.radiusKm, o.radiusKm, o.radiusKm});
    it.id = static_cast<int>(i);
    items.push_back(it);
  }

  bvh_.build(std::move(items), leafSize);
}

double OcclusionFieldKm::occlusion01Segment(const math::Vec3d& aKm,
                                           const math::Vec3d& bKm,
                                           int maxSphereTests) const {
  if (empty()) return 0.0;
  if (!math::isFinite(aKm) || !math::isFinite(bKm)) return 0.0;

  maxSphereTests = std::clamp(maxSphereTests, 1, 256);

  int tests = 0;
  double occ = 0.0;

  bvh_.querySegment(aKm, bKm, [&](int id) -> bool {
    if (id < 0 || static_cast<std::size_t>(id) >= occluders_.size()) return true;

    // Cap refinement work in pathological dense scenes.
    if (++tests > maxSphereTests) return false;

    const auto& o = occluders_[static_cast<std::size_t>(id)];

    // If either endpoint is inside the occluder, ignore it to avoid
    // "always occluded" artifacts when the sensor sits near/within a field.
    const double rInside = o.radiusKm * 0.98;
    if ((aKm - o.centerKm).lengthSq() < rInside * rInside) return true;
    if ((bKm - o.centerKm).lengthSq() < rInside * rInside) return true;

    if (!segmentHitsSphereKm_(aKm, bKm, o.centerKm, o.radiusKm)) return true;

    const double s = clamp01_(o.strength01);
    occ = 1.0 - (1.0 - occ) * (1.0 - s);
    if (occ >= 0.985) return false;
    return true;
  });

  return clamp01_(occ);
}

} // namespace stellar::sim
