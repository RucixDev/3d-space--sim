#pragma once

#include "stellar/core/Log.h"

#include <cstddef>
#include <cstdint>
#include <deque>
#include <mutex>
#include <string>
#include <string_view>
#include <vector>

namespace stellar::ui {

// Thread-safe ring buffer for log messages.
//
// Intended for tooling/UI integration (ImGui log viewer, on-screen diagnostics).
// Log callbacks can arrive from any thread (jobs, streaming, etc.), so this
// buffer uses a mutex to keep the data structure consistent.
//
// Notes:
//  - Entries are stored in insertion order.
//  - Each entry gets a monotonically increasing sequence id.
//  - When the capacity is exceeded, the oldest entries are dropped.
struct LogEntry {
  std::uint64_t seq{0};
  core::LogLevel level{core::LogLevel::Info};
  std::string timestamp;
  std::string message;
};

class LogBuffer {
public:
  explicit LogBuffer(std::size_t capacity = 4096);

  void setCapacity(std::size_t capacity);
  std::size_t capacity() const;

  void clear();

  std::size_t size() const;

  // Sequence values are stable across capacity changes.
  //
  // latestSeq() is the most recent seq currently stored, or 0 if empty.
  // earliestSeq() is the oldest seq currently stored, or (latestSeq()+1) if empty.
  std::uint64_t latestSeq() const;
  std::uint64_t earliestSeq() const;

  // Append a new entry.
  void push(core::LogLevel level, std::string_view timestamp, std::string_view message);

  // Copy all entries in insertion order.
  void copyAll(std::vector<LogEntry>& out) const;

  // Copy entries with seq > afterSeq in insertion order.
  // If maxCount is non-zero, at most that many entries are returned.
  void copySince(std::uint64_t afterSeq, std::vector<LogEntry>& out, std::size_t maxCount = 0) const;

  // Convenience: create a core::LogSink that writes into this buffer.
  core::LogSink makeSink();

  // Sink callback matching core::LogSink::Fn.
  static void sinkFn(core::LogLevel level,
                     std::string_view timestamp,
                     std::string_view message,
                     void* user);

private:
  void trimLocked();

  mutable std::mutex mutex_;
  std::deque<LogEntry> entries_;
  std::size_t capacity_{4096};
  std::uint64_t nextSeq_{1};
};

} // namespace stellar::ui
