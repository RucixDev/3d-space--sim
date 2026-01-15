#pragma once

#include <cmath>
#include <cstdint>

namespace stellar::core {

// Simple low-discrepancy sampling helpers for procedural generation.
//
// These are deterministic, lightweight, and header-only so they can be used
// across sim/proc/render without adding link dependencies.

static inline double frac(double x) {
  // Keep in [0,1) for any finite x.
  return x - std::floor(x);
}

// Halton sequence value for a given 1-based index (index >= 1) and base >= 2.
//
// Note: For procedural placement we only ever need small indices (tens/hundreds),
// so a straightforward implementation is fine.
static inline double halton(std::uint32_t index1Based, std::uint32_t base) {
  if (base < 2) base = 2;
  std::uint32_t i = index1Based;
  double f = 1.0;
  double r = 0.0;
  while (i > 0) {
    f /= static_cast<double>(base);
    r += f * static_cast<double>(i % base);
    i /= base;
  }
  return r;
}

struct Halton3 {
  double x{0.0};
  double y{0.0};
  double z{0.0};
};

static inline Halton3 halton3(std::uint32_t index1Based) {
  // Small coprime bases.
  return {halton(index1Based, 2), halton(index1Based, 3), halton(index1Based, 5)};
}

} // namespace stellar::core
