#pragma once

#include "stellar/core/Types.h"

#include <cmath>
#include <cstddef>
#include <cstring>
#include <string_view>

namespace stellar::core {

// A small, portable 64-bit *stable* hash builder.
//
// Designed for "regression signatures" in tests and headless tooling.
// The hash is stable across runs and platforms because:
//  - integers are serialized in explicit little-endian byte order
//  - strings are length-prefixed and hashed as raw bytes
//  - doubles can be hashed either by exact IEEE-754 bits or via quantization
//
// This is not cryptographic.
class StableHash64 {
public:
  static constexpr u64 kOffsetBasis = 14695981039346656037ull;
  static constexpr u64 kPrime       = 1099511628211ull;

  explicit StableHash64(u64 seed = kOffsetBasis) : h_(seed) {}

  void reset(u64 seed = kOffsetBasis) { h_ = seed; }
  u64 value() const { return h_; }

  void addByte(u8 b) {
    h_ ^= static_cast<u64>(b);
    h_ *= kPrime;
  }

  void addBytes(const void* data, std::size_t size) {
    if (!data || size == 0) return;
    const auto* p = static_cast<const u8*>(data);
    for (std::size_t i = 0; i < size; ++i) addByte(p[i]);
  }

  void addBool(bool v) { addByte(v ? 1u : 0u); }

  void addU16(u16 v) {
    const u8 b[2] = {
        static_cast<u8>((v >> 0) & 0xFFu),
        static_cast<u8>((v >> 8) & 0xFFu),
    };
    addBytes(b, 2);
  }

  void addU32(u32 v) {
    const u8 b[4] = {
        static_cast<u8>((v >> 0) & 0xFFu),
        static_cast<u8>((v >> 8) & 0xFFu),
        static_cast<u8>((v >> 16) & 0xFFu),
        static_cast<u8>((v >> 24) & 0xFFu),
    };
    addBytes(b, 4);
  }

  void addU64(u64 v) {
    const u8 b[8] = {
        static_cast<u8>((v >> 0) & 0xFFull),
        static_cast<u8>((v >> 8) & 0xFFull),
        static_cast<u8>((v >> 16) & 0xFFull),
        static_cast<u8>((v >> 24) & 0xFFull),
        static_cast<u8>((v >> 32) & 0xFFull),
        static_cast<u8>((v >> 40) & 0xFFull),
        static_cast<u8>((v >> 48) & 0xFFull),
        static_cast<u8>((v >> 56) & 0xFFull),
    };
    addBytes(b, 8);
  }

  void addI64(i64 v) { addU64(static_cast<u64>(v)); }
  void addInt(int v) { addI64(static_cast<i64>(v)); }

  // Hash a string with an explicit length prefix.
  // This prevents collisions like: ["ab","c"] vs ["a","bc"].
  void addString(std::string_view s) {
    addU64(static_cast<u64>(s.size()));
    addBytes(s.data(), s.size());
  }

  // Hash a double by its IEEE-754 bits.
  // Prefer this when you *want* exact bitwise stability.
  void addDoubleBits(double v) {
    u64 bits = 0;
    static_assert(sizeof(double) == sizeof(u64), "double must be 64-bit");
    std::memcpy(&bits, &v, sizeof(bits));
    addU64(bits);
  }

  // Hash a double by quantizing to an integer (default: 1e-9 units).
  // Prefer this when you want stability across minor floating differences.
  void addDoubleQ(double v, double scale = 1e9) {
    if (!std::isfinite(v) || !std::isfinite(scale) || scale <= 0.0) {
      // Encode "non-finite" distinctly.
      addU8(0xFFu);
      return;
    }
    const double scaled = v * scale;
    const i64 q = static_cast<i64>(std::llround(scaled));
    addI64(q);
  }

  void addU8(u8 v) { addByte(v); }

private:
  u64 h_{kOffsetBasis};
};

} // namespace stellar::core
