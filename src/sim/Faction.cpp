#include "stellar/sim/Faction.h"

#include "stellar/core/Random.h"

#include <array>
#include <span>
#include <string>

namespace stellar::sim {
namespace {

constexpr std::array<const char*, 8> kGovSuffix = {
  "Union",
  "Federation",
  "Republic",
  "Consortium",
  "Collective",
  "Directorate",
  "Kingdom",
  "Free State",
};

constexpr std::array<const char*, 6> kGovPrefix = {
  "United",
  "Outer",
  "Core",
  "New",
  "Greater",
  "Independent",
};

}

FactionGenerator::FactionGenerator(stellar::core::u64 universeSeed, std::size_t factionCount)
  : m_seed(universeSeed)
  , m_count(factionCount)
  , m_names(stellar::core::deriveSeed(universeSeed, "factions"))
{
  if (m_count == 0) m_count = 1;
}

Faction FactionGenerator::faction(std::size_t index) const {
  const auto base = stellar::core::deriveSeed(m_seed, "faction");
  const auto id = stellar::core::deriveSeed(base, static_cast<stellar::core::u64>(index));

  stellar::core::SplitMix64 rng(id);

  Faction f;
  f.id = id;

  const std::string coreName = m_names.makeName(rng, 2, 4);
  const std::string prefix = rng.chance(0.35) ? std::string(rng.pick(std::span(kGovPrefix))) + " " : std::string();
  const std::string suffix = std::string(rng.pick(std::span(kGovSuffix)));

  f.name = prefix + coreName + " " + suffix;

  // Lightly correlated attributes.
  f.techLevel = rng.uniform(0.15, 0.95);
  f.wealth = std::clamp(rng.uniform(0.10, 0.95) * (0.60 + 0.60 * f.techLevel), 0.0, 1.0);
  f.lawfulness = std::clamp(rng.uniform(0.05, 0.95) * (0.60 + 0.40 * f.wealth), 0.0, 1.0);
  return f;
}

std::size_t FactionGenerator::controllingFactionIndex(stellar::core::u64 systemId) const {
  const auto mapSeed = stellar::core::deriveSeed(m_seed, "faction_map");
  const auto h = stellar::core::hashCombine(systemId, mapSeed);
  return static_cast<std::size_t>(h % static_cast<stellar::core::u64>(m_count));
}

} // namespace stellar::sim
