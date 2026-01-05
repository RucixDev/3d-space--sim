#include "stellar/ui/LogBuffer.h"

#include <algorithm>

namespace stellar::ui {

LogBuffer::LogBuffer(std::size_t capacity) : capacity_(capacity) {
  if (capacity_ == 0) capacity_ = 1;
}

void LogBuffer::setCapacity(std::size_t capacity) {
  if (capacity == 0) capacity = 1;
  std::lock_guard<std::mutex> lock(mutex_);
  capacity_ = capacity;
  trimLocked();
}

std::size_t LogBuffer::capacity() const {
  std::lock_guard<std::mutex> lock(mutex_);
  return capacity_;
}

void LogBuffer::clear() {
  std::lock_guard<std::mutex> lock(mutex_);
  entries_.clear();
}

std::size_t LogBuffer::size() const {
  std::lock_guard<std::mutex> lock(mutex_);
  return entries_.size();
}

std::uint64_t LogBuffer::latestSeq() const {
  std::lock_guard<std::mutex> lock(mutex_);
  if (entries_.empty()) return 0;
  return entries_.back().seq;
}

std::uint64_t LogBuffer::earliestSeq() const {
  std::lock_guard<std::mutex> lock(mutex_);
  if (entries_.empty()) return nextSeq_; // convention: earliest = (latest+1)
  return entries_.front().seq;
}

void LogBuffer::push(core::LogLevel level, std::string_view timestamp, std::string_view message) {
  std::lock_guard<std::mutex> lock(mutex_);
  LogEntry e;
  e.seq = nextSeq_++;
  e.level = level;
  e.timestamp = std::string(timestamp);
  e.message = std::string(message);
  entries_.push_back(std::move(e));
  trimLocked();
}

void LogBuffer::copyAll(std::vector<LogEntry>& out) const {
  std::lock_guard<std::mutex> lock(mutex_);
  out.assign(entries_.begin(), entries_.end());
}

void LogBuffer::copySince(std::uint64_t afterSeq, std::vector<LogEntry>& out, std::size_t maxCount) const {
  std::lock_guard<std::mutex> lock(mutex_);

  out.clear();
  if (entries_.empty()) return;

  // Fast path: if cursor is ahead of our newest entry, nothing to return.
  if (afterSeq >= entries_.back().seq) return;

  // Clamp: if cursor is behind the oldest entry, return from the oldest.
  const std::uint64_t startSeq = std::max(afterSeq + 1, entries_.front().seq);

  // Linear scan (buffer sizes are small; this keeps it simple).
  for (const LogEntry& e : entries_) {
    if (e.seq < startSeq) continue;
    out.push_back(e);
    if (maxCount != 0 && out.size() >= maxCount) break;
  }
}

core::LogSink LogBuffer::makeSink() {
  return core::LogSink{&LogBuffer::sinkFn, this};
}

void LogBuffer::sinkFn(core::LogLevel level,
                       std::string_view timestamp,
                       std::string_view message,
                       void* user) {
  auto* buf = reinterpret_cast<LogBuffer*>(user);
  if (!buf) return;
  buf->push(level, timestamp, message);
}

void LogBuffer::trimLocked() {
  while (entries_.size() > capacity_) {
    entries_.pop_front();
  }
}

} // namespace stellar::ui
