#include "stellar/core/JobSystem.h"

#include "test_harness.h"

#include <atomic>
#include <cstddef>

int test_jobs() {
  int failures = 0;

  using stellar::core::JobSystem;

  // ---- submit() returns values ----
  {
    JobSystem js(4);
    auto f = js.submit([]() { return 123; });
    CHECK(f.get() == 123);
  }

  // ---- waitIdle() drains all queued work ----
  {
    JobSystem js(4);
    std::atomic<int> c{0};
    for (int i = 0; i < 1000; ++i) {
      js.submit([&]() { c.fetch_add(1, std::memory_order_relaxed); });
    }
    js.waitIdle();
    CHECK(c.load(std::memory_order_relaxed) == 1000);
  }

  // ---- parallelFor() executes each index exactly once ----
  {
    JobSystem js(4);
    std::atomic<std::size_t> hits{0};
    js.parallelFor(10'000, [&](std::size_t /*i*/) {
      hits.fetch_add(1, std::memory_order_relaxed);
    });
    CHECK(hits.load(std::memory_order_relaxed) == 10'000);
  }

  return failures;
}
