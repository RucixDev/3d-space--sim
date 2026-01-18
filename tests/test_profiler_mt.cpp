#include "test_harness.h"

#include "stellar/core/Profiler.h"

#include <thread>
#include <unordered_set>
#include <vector>

int test_profiler_mt() {
  int failures = 0;

  stellar::core::Profiler profiler;
  profiler.setEnabled(true);
  stellar::core::setActiveProfiler(&profiler);

  const std::uint64_t mainTid = stellar::core::Profiler::threadIdHash();

  profiler.beginFrame();

  // Spawn a couple threads that emit nested scopes.
  std::vector<std::thread> threads;
  threads.reserve(2);

  for (int i = 0; i < 2; ++i) {
    threads.emplace_back([i]() {
      (void)i;
      STELLAR_PROFILE_SCOPE("worker_scope");
      {
        STELLAR_PROFILE_SCOPE("worker_nested");
      }
    });
  }

  {
    STELLAR_PROFILE_SCOPE("main_join");
    for (auto& t : threads) {
      t.join();
    }
  }

  profiler.endFrame();

  const stellar::core::ProfilerFrame* f = profiler.newest();
  CHECK(f != nullptr);
  if (f) {
    CHECK(f->endNs > f->startNs);
    CHECK(f->mainThreadId == mainTid);
    CHECK(!f->events.empty());
    CHECK(f->events.size() >= 5); // 2 threads * (scope+nested) + main_join

    std::unordered_set<std::uint64_t> tids;
    for (const auto& e : f->events) {
      CHECK(e.threadId != 0);
      tids.insert(e.threadId);
    }

    // We expect at least main + 1 worker.
    CHECK(tids.size() >= 2);
  }

  stellar::core::setActiveProfiler(nullptr);
  return failures;
}
