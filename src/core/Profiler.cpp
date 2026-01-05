#include "stellar/core/Profiler.h"

#include <algorithm>
#include <atomic>
#include <chrono>

namespace stellar::core {

static std::atomic<Profiler*> gActiveProfiler{nullptr};

Profiler* activeProfiler() {
  return gActiveProfiler.load(std::memory_order_relaxed);
}

void setActiveProfiler(Profiler* p) {
  gActiveProfiler.store(p, std::memory_order_relaxed);
}

std::uint64_t Profiler::nowNs() {
  using clock = std::chrono::steady_clock;
  const auto now = clock::now().time_since_epoch();
  return (std::uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(now).count();
}

void Profiler::beginFrame() {
  if (!enabled_) return;

  inFrame_ = true;
  current_.startNs = nowNs();
  current_.endNs = 0;
  current_.events.clear();
  if (current_.events.capacity() < reserveEventsPerFrame_) {
    current_.events.reserve(reserveEventsPerFrame_);
  }
}

void Profiler::endFrame() {
  if (!enabled_) return;
  if (!inFrame_) return;

  // Mark as closed first so any scope destructors that happen after endFrame()
  // don't accidentally record into the next frame.
  inFrame_ = false;

  current_.endNs = nowNs();

  // Move the frame into history.
  frames_.push_back(std::move(current_));
  while (frames_.size() > maxFrames_) {
    frames_.pop_front();
  }

  // Prepare a fresh frame buffer.
  current_ = ProfilerFrame{};
  current_.events.reserve(reserveEventsPerFrame_);
}

void Profiler::clear() {
  frames_.clear();
  inFrame_ = false;

  current_ = ProfilerFrame{};
  current_.events.reserve(reserveEventsPerFrame_);
}

void Profiler::setMaxFrames(std::size_t maxFrames) {
  // Keep it sane.
  if (maxFrames < 1) maxFrames = 1;
  if (maxFrames > 4096) maxFrames = 4096;

  maxFrames_ = maxFrames;
  while (frames_.size() > maxFrames_) {
    frames_.pop_front();
  }
}

const ProfilerFrame* Profiler::newest() const {
  if (frames_.empty()) return nullptr;
  return &frames_.back();
}

void Profiler::record(const char* name,
                      std::uint64_t startNs,
                      std::uint64_t endNs,
                      std::uint32_t depth) {
  if (!enabled_) return;
  if (!inFrame_) return;

  ProfilerEvent e;
  e.name = name;
  e.startNs = startNs;
  e.endNs = endNs;
  e.depth = depth;
  current_.events.push_back(e);
}

static thread_local std::uint32_t gDepth = 0;

ProfilerScope::ProfilerScope(const char* name) : name_(name) {
  Profiler* p = activeProfiler();
  if (!p || !p->enabled()) {
    active_ = false;
    return;
  }

  active_ = true;
  depth_ = gDepth;
  ++gDepth;
  startNs_ = Profiler::nowNs();
}

ProfilerScope::~ProfilerScope() {
  if (!active_) return;

  // Always unwind depth even if profiling gets disabled mid-scope.
  if (gDepth > 0) {
    --gDepth;
  }

  Profiler* p = activeProfiler();
  if (!p || !p->enabled()) return;

  const std::uint64_t endNs = Profiler::nowNs();
  p->record(name_, startNs_, endNs, depth_);
}

} // namespace stellar::core
