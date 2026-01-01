#pragma once

#include "stellar/core/Hash.h"
#include "stellar/core/Types.h"

namespace stellar::sim {

// Tag bit for deterministic/procedural world objects.
//
// This is used to distinguish procedurally generated, stable IDs from transient IDs
// allocated at runtime.
//
// IMPORTANT: This value must remain stable across versions, since it is used to
// persist state (e.g., asteroid depletion) keyed by world object id.
constexpr core::u64 kDeterministicWorldIdBit = 0x2000000000000000ull;

inline bool isDeterministicWorldId(core::u64 id) {
  return (id & kDeterministicWorldIdBit) != 0;
}

// Build a stable, deterministic world object id from two 64-bit inputs.
// The resulting id is tagged with kDeterministicWorldIdBit.
inline core::u64 makeDeterministicWorldId(core::u64 a, core::u64 b) {
  return (core::hashCombine(a, b) | kDeterministicWorldIdBit);
}

} // namespace stellar::sim
