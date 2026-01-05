#pragma once

#include <cstddef>
#include <cstdint>
#include <deque>
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

  std::uint64_t durationNs() const {
    return (endNs >= startNs) ? (endNs - startNs) : 0;
  }
};

struct ProfilerFrame {
  std::uint64_t startNs{0};
  std::uint64_t endNs{0};
  std::vector<ProfilerEvent> events;

  std::uint64_t durationNs() const {
    return (endNs >= startNs) ? (endNs - startNs) : 0;
  }
};

class Profiler {
 public:
  Profiler() = default;

  void setEnabled(bool enabled) { enabled_ = enabled; }
  bool enabled() const { return enabled_; }

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
  bool enabled_{false};
  std::size_t maxFrames_{240};
  std::size_t reserveEventsPerFrame_{256};
  bool inFrame_{false};

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
