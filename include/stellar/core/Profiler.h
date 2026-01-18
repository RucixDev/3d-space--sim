#pragma once

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <mutex>
#include <vector>

namespace stellar::core {

// -----------------------------------------------------------------------------
// Lightweight CPU profiler (headless)
// -----------------------------------------------------------------------------
//
// Records nested scope timings into a per-frame ring buffer.
//
// Notes:
//  - Event names are stored as `const char*`.
//    Prefer string literals (or other permanently alive strings).

struct ProfilerEvent {
  const char* name{nullptr};
  std::uint64_t startNs{0};
  std::uint64_t endNs{0};
  std::uint32_t depth{0};
  std::uint64_t threadId{0}; // hashed std::thread::id (stable within a run)

  std::uint64_t durationNs() const {
    return (endNs >= startNs) ? (endNs - startNs) : 0;
  }
};

struct ProfilerFrame {
  std::uint64_t startNs{0};
  std::uint64_t endNs{0};
  std::uint64_t mainThreadId{0}; // hashed std::thread::id that called beginFrame()
  std::vector<ProfilerEvent> events;

  std::uint64_t durationNs() const {
    return (endNs >= startNs) ? (endNs - startNs) : 0;
  }
};

class Profiler {
 public:
  Profiler() = default;

  void setEnabled(bool enabled) { enabled_ = enabled; }
  bool enabled() const { return enabled_.load(std::memory_order_relaxed); }

  // Returns a stable, hashed identifier for the calling thread.
  // This is primarily used for trace export (Chrome/Perfetto) and debug UIs.
  static std::uint64_t threadIdHash();

  // Called by the host once per frame.
  void beginFrame();
  void endFrame();

  void clear();

  void setMaxFrames(std::size_t maxFrames);
  std::size_t maxFrames() const { return maxFrames_; }

  const std::deque<ProfilerFrame>& frames() const { return frames_; }
  const ProfilerFrame* newest() const;

  // Internal: called by ProfilerScope.
  void record(const char* name,
              std::uint64_t startNs,
              std::uint64_t endNs,
              std::uint32_t depth);

  static std::uint64_t nowNs();

 private:
  std::atomic<bool> enabled_{false};
  std::size_t maxFrames_{240};
  std::size_t reserveEventsPerFrame_{256};
  std::atomic<bool> inFrame_{false};

  // Protects current_ + frames_ so record() can be called from multiple threads.
  // The profiler is designed for low overhead; this mutex is only contended when
  // profiling spans are emitted from multiple threads.
  mutable std::mutex mutex_{};

  ProfilerFrame current_{};
  std::deque<ProfilerFrame> frames_{};
};

// Active profiler pointer (optional). This keeps the scope macro callsite tiny.
Profiler* activeProfiler();
void setActiveProfiler(Profiler* p);

// RAII scope helper. Records an event into the active profiler (if any) when destroyed.
class ProfilerScope {
 public:
  explicit ProfilerScope(const char* name);
  ~ProfilerScope();

  ProfilerScope(const ProfilerScope&) = delete;
  ProfilerScope& operator=(const ProfilerScope&) = delete;

 private:
  const char* name_{nullptr};
  std::uint64_t startNs_{0};
  std::uint32_t depth_{0};
  bool active_{false};
};

} // namespace stellar::core

// Token-pasting helper for unique variable names.
#define STELLAR_PP_CONCAT_INNER(a, b) a##b
#define STELLAR_PP_CONCAT(a, b) STELLAR_PP_CONCAT_INNER(a, b)

// Profile the current scope (records to the active profiler if enabled).
#define STELLAR_PROFILE_SCOPE(name_literal) \
  ::stellar::core::ProfilerScope STELLAR_PP_CONCAT(_stellar_prof_scope_, __LINE__)(name_literal)
