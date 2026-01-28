#pragma once

#include <iostream>
#include <utility>

// Minimal test helper to keep individual test files tidy.
//
// Usage:
//   int failures = 0;
//   CHECK(expr);
//   CHECK_EQ(a, b);
//
// The macro increments `failures` in the current scope and prints a message.

#ifndef CHECK
#define CHECK(expr)                                                                            \
  do {                                                                                         \
    if (!(expr)) {                                                                             \
      std::cerr << "[stellar_tests] CHECK failed: " #expr " (" << __FILE__ << ":" << __LINE__ \
                << ")\n";                                                                     \
      ++failures;                                                                              \
    }                                                                                          \
  } while (0)
#endif

#ifndef CHECK_EQ
#define CHECK_EQ(a, b)                                                                         \
  do {                                                                                         \
    auto _a = (a);                                                                             \
    auto _b = (b);                                                                             \
    if (!(_a == _b)) {                                                                         \
      std::cerr << "[stellar_tests] CHECK_EQ failed: " #a " == " #b " (" << __FILE__ << ":" \
                << __LINE__ << ")\n"                                                         \
                << "  left:  " << _a << "\n"                                                \
                << "  right: " << _b << "\n";                                               \
      ++failures;                                                                              \
    }                                                                                          \
  } while (0)
#endif
