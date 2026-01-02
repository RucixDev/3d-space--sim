#include "stellar/sim/Universe.h"
#include "stellar/core/JobSystem.h"
#include "stellar/proc/GalaxyGenerator.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

struct Item {
  stellar::sim::SystemStub stub{};
  double d2{0.0};
};

static bool betterItem(const Item& a, const Item& b) {
  if (a.d2 != b.d2) return a.d2 < b.d2;
  return a.stub.id < b.stub.id;
}

std::vector<stellar::sim::SystemStub> bruteForceQueryNearby(const stellar::sim::Universe& u,
                                                            const stellar::math::Vec3d& posLy,
                                                            double radiusLy,
                                                            std::size_t maxResults) {
  std::vector<stellar::sim::SystemStub> out;
  if (radiusLy <= 0.0 || maxResults == 0) return out;

  const double r2 = radiusLy * radiusLy;
  const double s = u.galaxyParams().sectorSizeLy;

  const stellar::proc::SectorCoord minC{
    static_cast<stellar::core::i32>(std::floor((posLy.x - radiusLy) / s)),
    static_cast<stellar::core::i32>(std::floor((posLy.y - radiusLy) / s)),
    static_cast<stellar::core::i32>(std::floor((posLy.z - radiusLy) / s)),
  };
  const stellar::proc::SectorCoord maxC{
    static_cast<stellar::core::i32>(std::floor((posLy.x + radiusLy) / s)),
    static_cast<stellar::core::i32>(std::floor((posLy.y + radiusLy) / s)),
    static_cast<stellar::core::i32>(std::floor((posLy.z + radiusLy) / s)),
  };

  std::vector<Item> items;
  stellar::proc::GalaxyGenerator gen(u.seed(), u.galaxyParams());

  for (stellar::core::i32 x = minC.x; x <= maxC.x; ++x) {
    for (stellar::core::i32 y = minC.y; y <= maxC.y; ++y) {
      for (stellar::core::i32 z = minC.z; z <= maxC.z; ++z) {
        const stellar::proc::SectorCoord c{x, y, z};
        const auto sec = gen.generateSector(c, u.factions());

        for (const auto& stub : sec.systems) {
          const double dd = (stub.posLy - posLy).lengthSq();
          if (dd <= r2) items.push_back(Item{stub, dd});
        }
      }
    }
  }

  std::sort(items.begin(), items.end(), betterItem);
  if (items.size() > maxResults) items.resize(maxResults);

  out.reserve(items.size());
  for (auto& it : items) out.push_back(std::move(it.stub));
  return out;
}

} // namespace

int test_query_nearby() {
  int fails = 0;

  const stellar::core::u64 seed = 7777777;
  stellar::sim::Universe u(seed);

  // Use a fixed-size pool for deterministic, cross-platform tests.
  stellar::core::JobSystem jobs(4);

  struct Case {
    stellar::math::Vec3d pos{};
    double radiusLy{0.0};
    std::size_t maxResults{0};
  };

  const std::vector<Case> cases{
    {{0.0, 0.0, 0.0}, 60.0, 64},
    {{37.3, -12.7, 5.1}, 40.0, 32},
    {{9.9, 9.9, 9.9}, 25.0, 128},
    {{-15.0, 24.0, -3.3}, 55.0, 8},
    {{0.0, 0.0, 0.0}, 20.0, 512},
  };

  for (std::size_t ci = 0; ci < cases.size(); ++ci) {
    const auto& c = cases[ci];

    const auto got = u.queryNearby(c.pos, c.radiusLy, c.maxResults);
    const auto gotPar = u.queryNearbyParallel(jobs, c.pos, c.radiusLy, c.maxResults);
    const auto ref = bruteForceQueryNearby(u, c.pos, c.radiusLy, c.maxResults);

    const auto checkList = [&](const char* label,
                               const std::vector<stellar::sim::SystemStub>& list) {
      if (list.size() != ref.size()) {
        std::cerr << "[test_query_nearby] size mismatch " << label << " case=" << ci
                  << " got=" << list.size() << " ref=" << ref.size()
                  << " (radius=" << c.radiusLy << " maxResults=" << c.maxResults << ")\n";
        ++fails;
        return;
      }

      for (std::size_t i = 0; i < list.size(); ++i) {
        if (list[i].id != ref[i].id) {
          std::cerr << "[test_query_nearby] id mismatch " << label << " case=" << ci << " idx=" << i
                    << " got=" << list[i].id << " ref=" << ref[i].id << "\n";
          ++fails;
          break;
        }
      }
    };

    checkList("serial", got);
    checkList("parallel", gotPar);

  }

  if (fails == 0) std::cout << "[test_query_nearby] pass\n";
  return fails;
}
