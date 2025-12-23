#pragma once

#include "stellar/core/Random.h"

#include <string>

namespace stellar::proc {

// Deterministic, lightweight name generator (syllable-based).
class NameGenerator {
public:
  explicit NameGenerator(core::u64 seed = 0) : rng_(seed) {}

  void reseed(core::u64 seed) { rng_.reseed(seed); }

  std::string systemName();
  std::string planetName(const std::string& systemName, int index);
  std::string stationName(const std::string& systemName, int index);

private:
  core::SplitMix64 rng_;
};

} // namespace stellar::proc
