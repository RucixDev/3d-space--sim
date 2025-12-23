#pragma once
#include "stellar/core/Types.h"
#include <string_view>

namespace stellar::core {

// 64-bit FNV-1a hash (stable, fast, good for IDs / seeds).
u64 fnv1a64(std::string_view text);

// Mix/combine two 64-bit hashes into one (order-sensitive).
u64 hashCombine(u64 a, u64 b);

// Hash an arbitrary POD-ish value by bytes (for seeds).
u64 hashBytes(const void* data, std::size_t size);

} // namespace stellar::core
