#pragma once

#include "stellar/core/Profiler.h"

#include <cstdint>
#include <deque>
#include <string>
#include <vector>

namespace stellar::core {

// -----------------------------------------------------------------------------
// Chrome JSON Trace exporter (chrome://tracing / Perfetto legacy JSON)
// -----------------------------------------------------------------------------
//
// Writes a trace file compatible with the "Trace Event" / Chrome JSON trace format.
// This can be opened in chrome://tracing (legacy) and in the Perfetto UI.
//
// The exporter is intentionally lightweight:
//  - We emit Complete ("X") events for all recorded scopes.
//  - Optional frame-span events can be included.
//  - By default, events are split by the profiler's recorded thread ids and
//    exported on separate tids (so multi-thread work shows up on distinct tracks).
//    When splitThreads=false, all span events are collapsed onto opt.tid.

struct ChromeTraceWriteOptions {
  // If true, emit a "Frame" span event for each captured profiler frame.
  bool includeFrameEvents{true};

  // When true, the JSON is pretty printed (larger files, easier to read).
  bool pretty{false};

  // If true, emit spans on separate tids per ProfilerEvent::threadId.
  // This allows viewing worker-thread work in Perfetto / chrome://tracing.
  bool splitThreads{true};

  // Process/thread identifiers used by trace viewers.
  int pid{1};
  int tid{1};
};

// Write all frames in the provided deque.
// Returns false on failure and optionally sets err.
bool writeProfilerChromeTraceJson(const char* path,
                                 const std::deque<ProfilerFrame>& frames,
                                 const ChromeTraceWriteOptions& opt = {},
                                 std::string* err = nullptr);

// Convenience helper: write a single frame.
bool writeProfilerChromeTraceJson(const char* path,
                                 const ProfilerFrame& frame,
                                 const ChromeTraceWriteOptions& opt = {},
                                 std::string* err = nullptr);

// -----------------------------------------------------------------------------
// Counter events ("ph":"C") helper
// -----------------------------------------------------------------------------
//
// Chrome JSON counter events allow visualizing changing numeric values over time.
// We expose a simple table format: N timestamps * K keys, stored row-major.
// Each sample produces one counter event with args {key0: v0, key1: v1, ...}.

struct ChromeTraceCounterTable {
  std::string name{"counter"};
  std::string category{"counter"};

  // K keys (args object members).
  std::vector<std::string> keys;

  // N timestamps in microseconds relative to an arbitrary origin (usually 0).
  std::vector<std::uint64_t> tsUs;

  // Flattened row-major values (tsUs.size() * keys.size()).
  std::vector<double> values;
};

// Write a trace JSON containing only counter events.
// Returns false on failure and optionally sets err.
bool writeCounterChromeTraceJson(const char* path,
                                const ChromeTraceCounterTable& table,
                                const ChromeTraceWriteOptions& opt = {},
                                std::string* err = nullptr);

// -----------------------------------------------------------------------------
// Combined export helper
// -----------------------------------------------------------------------------
//
// Writes a single Chrome JSON trace that contains:
//  - CPU profiler span events ("X") from the provided frames
//  - one or more counter tables ("C")
//
// This is useful for overlaying CPU breakdowns with time-series metrics such as
// GPU pass timings.
bool writeProfilerChromeTraceJsonWithCounters(const char* path,
                                             const std::deque<ProfilerFrame>& frames,
                                             const std::vector<ChromeTraceCounterTable>& counters,
                                             const ChromeTraceWriteOptions& opt = {},
                                             std::string* err = nullptr);

} // namespace stellar::core
