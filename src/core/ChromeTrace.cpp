#include "stellar/core/ChromeTrace.h"

#include "stellar/core/JsonWriter.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

namespace stellar::core {

static std::uint64_t nsToUs(std::uint64_t ns) {
  return ns / 1000ull;
}

static void writeProcessThreadMetadata(JsonWriter& w, int pid, int tid) {
  // process_name metadata
  w.beginObject();
  w.key("name"); w.value("process_name");
  w.key("ph"); w.value("M");
  w.key("pid"); w.value(pid);
  w.key("tid"); w.value(tid);
  w.key("args");
  w.beginObject();
  w.key("name"); w.value("stellar");
  w.endObject();
  w.endObject();

  // thread_name metadata
  w.beginObject();
  w.key("name"); w.value("thread_name");
  w.key("ph"); w.value("M");
  w.key("pid"); w.value(pid);
  w.key("tid"); w.value(tid);
  w.key("args");
  w.beginObject();
  w.key("name"); w.value("Main");
  w.endObject();
  w.endObject();
}

static void writeSpanEvent(JsonWriter& w,
                           std::string_view name,
                           std::uint64_t tsUs,
                           std::uint64_t durUs,
                           int pid,
                           int tid,
                           std::uint32_t depth,
                           std::string_view cat) {
  w.beginObject();
  w.key("name"); w.value(name);
  w.key("cat"); w.value(cat);
  w.key("ph"); w.value("X");
  w.key("ts"); w.value((unsigned long long)tsUs);
  w.key("dur"); w.value((unsigned long long)durUs);
  w.key("pid"); w.value(pid);
  w.key("tid"); w.value(tid);
  w.key("args");
  w.beginObject();
  w.key("depth"); w.value((unsigned long long)depth);
  w.endObject();
  w.endObject();
}

static bool writeFramesImpl(std::ostream& out,
                            const std::deque<ProfilerFrame>& frames,
                            const ChromeTraceWriteOptions& opt,
                            std::string* err) {
  if (frames.empty()) {
    if (err) *err = "No frames to export.";
    return false;
  }

  const std::uint64_t baseNs = frames.front().startNs;

  JsonWriter w(out, opt.pretty);
  w.beginObject();
  w.key("displayTimeUnit");
  w.value("ms");

  w.key("traceEvents");
  w.beginArray();

  // Emit metadata first so viewers can label the track.
  writeProcessThreadMetadata(w, opt.pid, opt.tid);

  // A scratch vector used to sort events per-frame by start time.
  std::vector<const ProfilerEvent*> sorted;

  for (const auto& f : frames) {
    if (f.endNs <= f.startNs) continue;

    const std::uint64_t frameTsUs = nsToUs((f.startNs >= baseNs) ? (f.startNs - baseNs) : 0ull);
    const std::uint64_t frameDurUs = nsToUs(f.durationNs());

    if (opt.includeFrameEvents) {
      writeSpanEvent(w,
                     /*name=*/"Frame",
                     /*tsUs=*/frameTsUs,
                     /*durUs=*/frameDurUs,
                     /*pid=*/opt.pid,
                     /*tid=*/opt.tid,
                     /*depth=*/0,
                     /*cat=*/"frame");
    }

    if (f.events.empty()) continue;

    sorted.clear();
    sorted.reserve(f.events.size());
    for (const auto& e : f.events) {
      sorted.push_back(&e);
    }
    std::sort(sorted.begin(), sorted.end(), [](const ProfilerEvent* a, const ProfilerEvent* b) {
      if (a->startNs != b->startNs) return a->startNs < b->startNs;
      if (a->endNs != b->endNs) return a->endNs < b->endNs;
      return a->depth < b->depth;
    });

    for (const ProfilerEvent* e : sorted) {
      const char* n = e->name ? e->name : "<null>";
      const std::uint64_t tsUs = nsToUs((e->startNs >= baseNs) ? (e->startNs - baseNs) : 0ull);
      const std::uint64_t durUs = nsToUs(e->durationNs());
      writeSpanEvent(w,
                     /*name=*/n,
                     /*tsUs=*/tsUs,
                     /*durUs=*/durUs,
                     /*pid=*/opt.pid,
                     /*tid=*/opt.tid,
                     /*depth=*/e->depth,
                     /*cat=*/"cpu");
    }
  }

  w.endArray();
  w.endObject();
  return true;
}

bool writeProfilerChromeTraceJson(const char* path,
                                 const std::deque<ProfilerFrame>& frames,
                                 const ChromeTraceWriteOptions& opt,
                                 std::string* err) {
  if (!path || !*path) {
    if (err) *err = "Invalid path.";
    return false;
  }

  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f.is_open()) {
    if (err) *err = "Failed to open file.";
    return false;
  }

  return writeFramesImpl(f, frames, opt, err);
}

bool writeProfilerChromeTraceJson(const char* path,
                                 const ProfilerFrame& frame,
                                 const ChromeTraceWriteOptions& opt,
                                 std::string* err) {
  std::deque<ProfilerFrame> frames;
  frames.push_back(frame);
  return writeProfilerChromeTraceJson(path, frames, opt, err);
}



static void writeCounterEvent(JsonWriter& w,
                              std::string_view name,
                              std::string_view cat,
                              std::uint64_t tsUs,
                              int pid,
                              int tid,
                              const std::vector<std::string>& keys,
                              const double* row) {
  w.beginObject();
  w.key("name"); w.value(name);
  w.key("cat"); w.value(cat);
  w.key("ph"); w.value("C");
  w.key("ts"); w.value((unsigned long long)tsUs);
  w.key("pid"); w.value(pid);
  w.key("tid"); w.value(tid);
  w.key("args");
  w.beginObject();
  for (std::size_t k = 0; k < keys.size(); ++k) {
    w.key(keys[k]);
    w.value(row ? row[k] : 0.0);
  }
  w.endObject();
  w.endObject();
}

static bool writeCountersImpl(std::ostream& out,
                              const ChromeTraceCounterTable& table,
                              const ChromeTraceWriteOptions& opt,
                              std::string* err) {
  const std::size_t k = table.keys.size();
  const std::size_t n = table.tsUs.size();
  if (k == 0) {
    if (err) *err = "No counter keys provided.";
    return false;
  }
  if (n == 0) {
    if (err) *err = "No counter samples provided.";
    return false;
  }
  if (table.values.size() != n * k) {
    if (err) *err = "Counter table values size mismatch.";
    return false;
  }

  JsonWriter w(out, opt.pretty);
  w.beginObject();
  w.key("displayTimeUnit");
  w.value("ms");

  w.key("traceEvents");
  w.beginArray();

  writeProcessThreadMetadata(w, opt.pid, opt.tid);

  const std::string_view name = table.name.empty() ? std::string_view{"counter"} : std::string_view{table.name};
  const std::string_view cat = table.category.empty() ? std::string_view{"counter"} : std::string_view{table.category};

  for (std::size_t i = 0; i < n; ++i) {
    const double* row = &table.values[i * k];
    writeCounterEvent(w, name, cat, table.tsUs[i], opt.pid, opt.tid, table.keys, row);
  }

  w.endArray();
  w.endObject();
  return true;
}

bool writeCounterChromeTraceJson(const char* path,
                                const ChromeTraceCounterTable& table,
                                const ChromeTraceWriteOptions& opt,
                                std::string* err) {
  if (!path || !*path) {
    if (err) *err = "Invalid path.";
    return false;
  }

  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f.is_open()) {
    if (err) *err = "Failed to open file.";
    return false;
  }

  return writeCountersImpl(f, table, opt, err);
}

} // namespace stellar::core
