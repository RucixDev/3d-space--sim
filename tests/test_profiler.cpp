#include "stellar/core/Profiler.h"

#include "test_harness.h"

#include <chrono>
#include <cstring>
#include <thread>

int test_profiler() {
  int failures = 0;

  using stellar::core::Profiler;
  using stellar::core::ProfilerEvent;
  using stellar::core::setActiveProfiler;

  Profiler p;
  setActiveProfiler(&p);

  p.setEnabled(true);
  p.setMaxFrames(3);

  // ---- One frame with nesting
  p.beginFrame();
  {
    STELLAR_PROFILE_SCOPE("A");
    std::this_thread::sleep_for(std::chrono::microseconds(100));
    {
      STELLAR_PROFILE_SCOPE("B");
      std::this_thread::sleep_for(std::chrono::microseconds(50));
    }
  }
  p.endFrame();

  CHECK(p.frames().size() == 1);

  const auto* f0 = p.newest();
  CHECK(f0 != nullptr);
  CHECK(f0->endNs >= f0->startNs);
  CHECK(!f0->events.empty());

  bool foundA = false;
  bool foundB = false;
  for (const ProfilerEvent& e : f0->events) {
    CHECK(e.name != nullptr);
    CHECK(e.endNs >= e.startNs);

    if (std::strcmp(e.name, "A") == 0) {
      foundA = true;
      CHECK(e.depth == 0);
    }
    if (std::strcmp(e.name, "B") == 0) {
      foundB = true;
      CHECK(e.depth == 1);
    }
  }

  CHECK(foundA);
  CHECK(foundB);

  // ---- Ring buffer trimming
  for (int i = 0; i < 5; ++i) {
    p.beginFrame();
    { STELLAR_PROFILE_SCOPE("Frame"); }
    p.endFrame();
  }
  CHECK(p.frames().size() == 3);

  // ---- Disabled profiler should not record
  const std::size_t before = p.frames().size();
  p.setEnabled(false);
  p.beginFrame();
  { STELLAR_PROFILE_SCOPE("Nope"); }
  p.endFrame();
  CHECK(p.frames().size() == before);

  setActiveProfiler(nullptr);

  return failures;
}
