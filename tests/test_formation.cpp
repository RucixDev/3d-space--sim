#include "test_harness.h"

#include "stellar/sim/Formation.h"

#include <cmath>
#include <vector>

int test_formation() {
  int failures = 0;

  using namespace stellar;

  sim::FormationParams p{};
  p.radiusKm = 40000.0;
  p.lateralScale = 0.50;
  p.verticalScale = 0.25;
  p.behindScale = 1.00;
  p.rowSpacingScale = 0.60;
  p.jitterKm = 0.0;

  sim::FormationParams pj = p;
  pj.jitterKm = 1000.0;

  const core::u64 seed = 0x123456789abcdef0ull;

  auto finite = [](double v) { return std::isfinite(v); };
  auto vfinite = [&](const math::Vec3d& v) { return finite(v.x) && finite(v.y) && finite(v.z); };

  const sim::FormationPattern patterns[] = {
    sim::FormationPattern::Trail,
    sim::FormationPattern::Wedge,
    sim::FormationPattern::Diamond,
    sim::FormationPattern::Ring,
  };

  for (const auto pat : patterns) {
    for (int n = 1; n <= 10; ++n) {
      std::vector<math::Vec3d> offsets;
      offsets.reserve((std::size_t)n);

      for (int i = 0; i < n; ++i) {
        // Determinism should hold even with non-zero jitter.
        const math::Vec3d a = sim::formationOffsetLocalKm(pat, i, n, pj, seed);
        const math::Vec3d b = sim::formationOffsetLocalKm(pat, i, n, pj, seed);
        CHECK(vfinite(a));
        CHECK((a - b).length() < 1e-9);
        // Use jitter-free positions for spacing sanity (jitter is allowed to reduce separation).
        offsets.push_back(sim::formationOffsetLocalKm(pat, i, n, p, seed));
      }

      // Basic spacing sanity: we want escorts to be meaningfully separated.
      if (n > 1) {
        double minDist = 1e30;
        for (int i = 0; i < n; ++i) {
          for (int j = i + 1; j < n; ++j) {
            minDist = std::min(minDist, (offsets[(std::size_t)i] - offsets[(std::size_t)j]).length());
          }
        }
        CHECK(minDist > p.radiusKm * 0.15);
      }
    }
  }

  // World transform should behave like a basis transform.
  {
    const sim::Basis basis{math::Vec3d{1,0,0}, math::Vec3d{0,1,0}, math::Vec3d{0,0,1}};
    const math::Vec3d pos{10, 20, 30};
    const math::Vec3d off{1, 2, 3};
    const math::Vec3d tgt = sim::formationTargetWorldKm(pos, basis, off);
    CHECK(std::abs(tgt.x - 11.0) < 1e-12);
    CHECK(std::abs(tgt.y - 22.0) < 1e-12);
    CHECK(std::abs(tgt.z - 33.0) < 1e-12);
  }

  // Basis builder: should survive parallel forward/up inputs.
  {
    const sim::Basis b = sim::makeBasisFromForward(math::Vec3d{0, 1, 0}, math::Vec3d{0, 1, 0});
    CHECK(b.right.lengthSq() > 0.5);
    CHECK(b.up.lengthSq() > 0.5);
    CHECK(b.fwd.lengthSq() > 0.5);
    CHECK(std::abs(math::dot(b.right, b.up)) < 1e-9);
    CHECK(std::abs(math::dot(b.right, b.fwd)) < 1e-9);
    CHECK(std::abs(math::dot(b.up, b.fwd)) < 1e-9);
  }

  return failures;
}
