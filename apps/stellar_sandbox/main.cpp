#include "stellar/core/Log.h"
#include "stellar/proc/GalaxyGenerator.h"
#include "stellar/sim/Universe.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

using namespace stellar;

static sim::StarClass parseStarClass(sim::StarClass c) { return c; }

static const char* starClassName(sim::StarClass c) {
  switch (c) {
    case sim::StarClass::O: return "O";
    case sim::StarClass::B: return "B";
    case sim::StarClass::A: return "A";
    case sim::StarClass::F: return "F";
    case sim::StarClass::G: return "G";
    case sim::StarClass::K: return "K";
    case sim::StarClass::M: return "M";
    default: return "?";
  }
}

static void printHelp() {
  std::cout << "stellar_sandbox\n"
            << "  --seed <u64>           Galaxy seed (default: 1337)\n"
            << "  --pos <x y z>          Query position in ly (default: 0 0 0)\n"
            << "  --radius <ly>          Query radius in ly (default: 50)\n"
            << "  --limit <n>            Max systems (default: 32)\n";
}

int main(int argc, char** argv) {
  core::setLogLevel(core::LogLevel::Info);

  core::u64 seed = 1337;
  math::Vec3d posLy{0,0,0};
  double radiusLy = 50.0;
  std::size_t limit = 32;

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--help" || a == "-h") {
      printHelp();
      return 0;
    } else if (a == "--seed" && i + 1 < argc) {
      seed = static_cast<core::u64>(std::strtoull(argv[++i], nullptr, 10));
    } else if (a == "--pos" && i + 3 < argc) {
      posLy.x = std::atof(argv[++i]);
      posLy.y = std::atof(argv[++i]);
      posLy.z = std::atof(argv[++i]);
    } else if (a == "--radius" && i + 1 < argc) {
      radiusLy = std::atof(argv[++i]);
    } else if (a == "--limit" && i + 1 < argc) {
      limit = static_cast<std::size_t>(std::strtoull(argv[++i], nullptr, 10));
    } else {
      std::cerr << "Unknown arg: " << a << "\n";
      printHelp();
      return 1;
    }
  }

  sim::Universe u(seed);

  const auto systems = u.queryNearby(posLy, radiusLy, limit);

  std::cout << "Seed: " << seed << "\n";
  std::cout << "Query @ (" << posLy.x << "," << posLy.y << "," << posLy.z << ") radius=" << radiusLy << " ly\n";
  std::cout << "Found " << systems.size() << " systems\n\n";

  for (const auto& s : systems) {
    const math::Vec3d d = s.posLy - posLy;
    const double dist = std::sqrt(d.lengthSq());

    std::cout << std::setw(14) << s.id
              << "  " << std::setw(14) << s.name
              << "  class=" << starClassName(s.primaryClass)
              << "  dist=" << std::fixed << std::setprecision(2) << dist << " ly"
              << "  planets=" << s.planetCount
              << "  stations=" << s.stationCount
              << "  faction=" << s.factionId
              << "\n";
  }

  if (!systems.empty()) {
    std::cout << "\n--- Example system detail: " << systems.front().name << " ---\n";
    const auto& sys = u.getSystem(systems.front().id, &systems.front());
    std::cout << "Star mass=" << sys.star.massSol << " Sol, lum=" << sys.star.luminositySol << " Sol\n";
    std::cout << "Planets:\n";
    for (const auto& p : sys.planets) {
      std::cout << "  - " << p.name
                << " a=" << p.orbit.semiMajorAxisAU << " AU"
                << " e=" << p.orbit.eccentricity
                << " period=" << p.orbit.periodDays << " days\n";
    }
    std::cout << "Stations:\n";
    for (const auto& st : sys.stations) {
      std::cout << "  - " << st.name << " (fee=" << st.feeRate << ")\n";
    }
  }

  return 0;
}
