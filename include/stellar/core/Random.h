#pragma once

#include "stellar/core/Hash.h"
#include "stellar/core/Types.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string_view>
#include <type_traits>

namespace stellar::core {

// SplitMix64: fast, simple PRNG with 64-bit state.
// Great for procedural generation & seeding other RNGs.
// Not suitable for crypto.
class SplitMix64 {
public:
  explicit SplitMix64(u64 seed = 0) : state_(seed) {}

  void reseed(u64 seed) { state_ = seed; }

  u64 nextU64() {
    u64 z = (state_ += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }

  u32 nextU32() { return static_cast<u32>(nextU64() >> 32); }

  // [0,1)
  double nextDouble() {
    // 53 random bits -> double in [0,1).
    constexpr double inv = 1.0 / static_cast<double>(1ull << 53);
    return static_cast<double>(nextU64() >> 11) * inv;
  }

  // Alias kept for older gameplay/proc-gen code.
  // Returns a double in [0,1).
  double nextUnit() { return nextDouble(); }

  // Inclusive range for integers
  template <class Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0>
  Int range(Int minInclusive, Int maxInclusive) {
    if (maxInclusive < minInclusive) std::swap(minInclusive, maxInclusive);
    const u64 span = static_cast<u64>(maxInclusive) - static_cast<u64>(minInclusive) + 1ull;
    return static_cast<Int>(minInclusive + static_cast<Int>(nextU64() % span));
  }

  // Range for floating point
  template <class Float, std::enable_if_t<std::is_floating_point_v<Float>, int> = 0>
  Float range(Float minInclusive, Float maxInclusive) {
    if (maxInclusive < minInclusive) std::swap(minInclusive, maxInclusive);
    return static_cast<Float>(minInclusive + (maxInclusive - minInclusive) * nextDouble());
  }

  bool chance(double p) { return nextDouble() < p; }

private:
  u64 state_{0};
};

inline u64 seedFromText(std::string_view text) { return fnv1a64(text); }

} // namespace stellar::core
