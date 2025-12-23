#pragma once

#include "stellar/core/Types.h"
#include "stellar/proc/NameGenerator.h"

#include <cstddef>
#include <string>

namespace stellar::sim {

struct Faction {
  stellar::core::u64 id = 0;
  std::string name;

  // Simple 0..1 attributes used by procedural economy / game rules.
  double techLevel = 0.5;
  double lawfulness = 0.5;
  double wealth = 0.5;
};

// Deterministic faction generator / registry.
// Factions are generated purely from universeSeed + faction index, so they don't need persistence.
class FactionGenerator {
public:
  explicit FactionGenerator(stellar::core::u64 universeSeed, std::size_t factionCount = 32);

  std::size_t count() const { return m_count; }

  Faction faction(std::size_t index) const;

  // Maps a star system id to its controlling faction index.
  std::size_t controllingFactionIndex(stellar::core::u64 systemId) const;

  Faction controllingFaction(stellar::core::u64 systemId) const {
    return faction(controllingFactionIndex(systemId));
  }

private:
  stellar::core::u64 m_seed = 0;
  std::size_t m_count = 0;
  stellar::proc::NameGenerator m_names;
};

} // namespace stellar::sim
