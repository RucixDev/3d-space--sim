#include "stellar/core/Hash.h"

#include <cstddef>

namespace stellar::core {

u64 fnv1a64(std::string_view text) {
  constexpr u64 offsetBasis = 14695981039346656037ull;
  constexpr u64 prime       = 1099511628211ull;

  u64 hash = offsetBasis;
  for (unsigned char c : text) {
    hash ^= static_cast<u64>(c);
    hash *= prime;
  }
  return hash;
}

u64 hashCombine(u64 a, u64 b) {
  // A common 64-bit combine (based on boost::hash_combine-ish mixing).
  // Not cryptographic, but good for procedural IDs/seeds.
  u64 x = a;
  x ^= b + 0x9E3779B97F4A7C15ull + (x << 6) + (x >> 2);
  // final avalanche
  x ^= x >> 33;
  x *= 0xff51afd7ed558ccdULL;
  x ^= x >> 33;
  x *= 0xc4ceb9fe1a85ec53ULL;
  x ^= x >> 33;
  return x;
}

u64 hashBytes(const void* data, std::size_t size) {
  const auto* bytes = static_cast<const unsigned char*>(data);
  constexpr u64 offsetBasis = 14695981039346656037ull;
  constexpr u64 prime       = 1099511628211ull;

  u64 hash = offsetBasis;
  for (std::size_t i = 0; i < size; ++i) {
    hash ^= static_cast<u64>(bytes[i]);
    hash *= prime;
  }
  return hash;
}

} // namespace stellar::core
