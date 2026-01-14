#include "stellar/proc/GalaxyGenerator.h"
#include "stellar/sim/Faction.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace stellar;

static double dist2(const math::Vec3d& a, const math::Vec3d& b) {
  const double dx = a.x - b.x;
  const double dy = a.y - b.y;
  const double dz = a.z - b.z;
  return dx * dx + dy * dy + dz * dz;
}

int test_galaxy_minsep() {
  proc::GalaxyParams gp{};
  gp.sectorSizeLy = 10.0;
  gp.radiusLy = 500.0;
  gp.thicknessLy = 200.0;
  gp.radialScaleLengthLy = 150.0;
  gp.baseMeanSystemsPerSector = 80.0;

  // Enable the min-separation (blue-noise-ish) placement path.
  gp.minSystemSeparationLy = 1.0;

  proc::GalaxyGenerator gen(1337ull, gp);

  // Factions don't affect spacing, but this exercises the faction assignment path.
  std::vector<sim::Faction> factions = sim::generateFactions(1337ull, 4);

  const proc::SectorCoord c{0, 0, 0};
  const proc::Sector sec = gen.generateSector(c, factions);

  if (sec.systems.empty()) {
    std::cerr << "[test_galaxy_minsep] expected non-empty sector systems\n";
    return 1;
  }

  // The generator clamps separation to >= 0.25 for perf.
  const double r = std::max(0.25, gp.minSystemSeparationLy);
  const double r2 = r * r;

  for (std::size_t i = 0; i < sec.systems.size(); ++i) {
    for (std::size_t j = i + 1; j < sec.systems.size(); ++j) {
      const double d2 = dist2(sec.systems[i].posLy, sec.systems[j].posLy);
      if (d2 + 1e-12 < r2) {
        std::cerr << "[test_galaxy_minsep] systems too close: d=" << std::sqrt(d2)
                  << " < min=" << r << " (i=" << i << " j=" << j << ")\n";
        return 1;
      }
    }
  }

  // Determinism: regenerate and compare.
  const proc::Sector sec2 = gen.generateSector(c, factions);
  if (sec.systems.size() != sec2.systems.size()) {
    std::cerr << "[test_galaxy_minsep] non-deterministic system count\n";
    return 1;
  }

  for (std::size_t i = 0; i < sec.systems.size(); ++i) {
    const auto& a = sec.systems[i];
    const auto& b = sec2.systems[i];

    if (a.id != b.id) {
      std::cerr << "[test_galaxy_minsep] non-deterministic id at i=" << i << "\n";
      return 1;
    }
    if (dist2(a.posLy, b.posLy) > 1e-24) {
      std::cerr << "[test_galaxy_minsep] non-deterministic pos at i=" << i << "\n";
      return 1;
    }
  }

  std::cout << "[test_galaxy_minsep] PASS\n";
  return 0;
}
