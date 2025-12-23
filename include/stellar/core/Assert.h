#pragma once

#include <cstdlib>
#include <string_view>

namespace stellar::core {

// Minimal assert helper. Kept small on purpose.
[[noreturn]] void panic(std::string_view message, const char* file, int line);

} // namespace stellar::core

#define STELLAR_ASSERT(expr) \
  do { \
    if (!(expr)) { \
      ::stellar::core::panic("Assertion failed: " #expr, __FILE__, __LINE__); \
    } \
  } while (0)

#define STELLAR_ASSERT_MSG(expr, msg) \
  do { \
    if (!(expr)) { \
      ::stellar::core::panic((msg), __FILE__, __LINE__); \
    } \
  } while (0)
