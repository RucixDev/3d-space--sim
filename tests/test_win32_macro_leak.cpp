#include "test_harness.h"

#include "stellar/core/AtomicWriteFile.h"

// Ensure Windows headers pulled in by AtomicWriteFile.h don't leak problematic macros.
//
// The Windows SDK historically defined `near`/`far` as macros (from 16-bit era).
// If these leak into user code they can break valid identifiers.
//
// This test is compile-time only on Windows. On non-Windows platforms it simply
// compiles and returns success.
int test_win32_macro_leak() {
  int failures = 0;

#if defined(_WIN32)
  #ifdef near
    #error "Windows macro 'near' leaked past AtomicWriteFile.h"
  #endif
  #ifdef far
    #error "Windows macro 'far' leaked past AtomicWriteFile.h"
  #endif
#endif

  CHECK(true);
  return failures;
}
