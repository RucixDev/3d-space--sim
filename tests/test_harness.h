#pragma once

#include <iostream>

// Minimal test helper to keep individual test files tidy.
//
// Usage:
//   int failures = 0;
//   CHECK(expr);
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
