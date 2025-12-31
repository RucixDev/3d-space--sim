#include "stellar/core/Log.h"
#include "tests/test_harness.h"

#include <atomic>
#include <string>

using namespace stellar;

static void testSink(core::LogLevel /*level*/, std::string_view /*ts*/, std::string_view /*msg*/, void* user) {
  auto* count = reinterpret_cast<std::atomic<int>*>(user);
  if (count) count->fetch_add(1, std::memory_order_relaxed);
}

int test_log_sinks() {
  int failures = 0;

  const core::LogLevel prev = core::getLogLevel();
  core::setLogLevel(core::LogLevel::Trace);

  std::atomic<int> count{0};
  const core::LogSink sink{&testSink, &count};

  core::addLogSink(sink);
  core::log(core::LogLevel::Info, "hello");
  CHECK(count.load(std::memory_order_relaxed) == 1);

  // Removing should stop callbacks.
  core::removeLogSink(sink);
  core::log(core::LogLevel::Info, "world");
  CHECK(count.load(std::memory_order_relaxed) == 1);

  // Respect log-level filtering.
  core::addLogSink(sink);
  core::setLogLevel(core::LogLevel::Off);
  core::log(core::LogLevel::Error, "should_not_fire");
  CHECK(count.load(std::memory_order_relaxed) == 1);

  // Restore.
  core::removeLogSink(sink);
  core::setLogLevel(prev);
  return failures;
}
