#include "stellar/sim/Faction.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/proc/NameGenerator.h"

#include <array>
#include <cmath>

namespace stellar::sim {

std::vector<Faction> generateFactions(core::u64 seed, int count) {
  std::vector<Faction> out;
  out.reserve(static_cast<std::size_t>(count + 1));

  // id=0 reserved
  out.push_back(Faction{});

  static const std::array<const char*, 6> suffix = {"Union","League","Combine","Consortium","Republic","Directorate"};

  for (int i = 0; i < count; ++i) {
    const core::u32 id = static_cast<core::u32>(i + 1);
    core::SplitMix64 rng(core::hashCombine(seed, static_cast<core::u64>(id) * 0x9E3779B97F4A7C15ull));

    proc::NameGenerator ng(core::hashCombine(seed, id));
    std::string name = ng.systemName() + " " + suffix[rng.range<int>(0, static_cast<int>(suffix.size()-1))];

    // Place faction homes with a disc-ish distribution.
    const double r = 20000.0 * std::sqrt(rng.nextDouble()); // within 20k ly for now
    const double a = rng.nextDouble() * 6.283185307179586;
    const double z = (rng.nextDouble() - 0.5) * 600.0; // +-300 ly

    math::Vec3d home{ r * std::cos(a), r * std::sin(a), z };

    Faction f{};
    f.id = id;
    f.name = std::move(name);
    f.homePosLy = home;
    f.influenceRadiusLy = 800.0 + 1200.0 * rng.nextDouble();
    f.taxRate = 0.01 + 0.05 * rng.nextDouble();
    f.industryBias = rng.range(-1.0, 1.0);

    // simple bright-ish color
    f.color = {0.2 + 0.8*rng.nextDouble(), 0.2 + 0.8*rng.nextDouble(), 0.2 + 0.8*rng.nextDouble()};

    out.push_back(std::move(f));
  }

  return out;
}

} // namespace stellar::sim
