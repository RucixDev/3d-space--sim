#include "stellar/ui/LogBuffer.h"

#include "test_harness.h"

#include <vector>

using namespace stellar;

int test_log_buffer() {
  int failures = 0;

  ui::LogBuffer buf(/*capacity=*/2);

  buf.push(core::LogLevel::Info, "t1", "a");
  buf.push(core::LogLevel::Warn, "t2", "b");

  CHECK(buf.size() == 2);
  CHECK(buf.earliestSeq() == 1);
  CHECK(buf.latestSeq() == 2);

  // Capacity should drop the oldest.
  buf.push(core::LogLevel::Error, "t3", "c");
  CHECK(buf.size() == 2);
  CHECK(buf.earliestSeq() == 2);
  CHECK(buf.latestSeq() == 3);

  std::vector<ui::LogEntry> all;
  buf.copyAll(all);
  CHECK(all.size() == 2);
  CHECK(all[0].seq == 2);
  CHECK(all[0].message == "b");
  CHECK(all[1].seq == 3);
  CHECK(all[1].message == "c");

  // copySince should return seq > cursor.
  std::vector<ui::LogEntry> since;
  buf.copySince(2, since);
  CHECK(since.size() == 1);
  CHECK(since[0].seq == 3);

  // copySince with maxCount.
  buf.copySince(0, since, /*maxCount=*/1);
  CHECK(since.size() == 1);
  CHECK(since[0].seq == 2);

  // Sink callback should append like push().
  ui::LogBuffer::sinkFn(core::LogLevel::Info, "t4", "d", &buf);
  CHECK(buf.latestSeq() == 4);
  buf.copyAll(all);
  CHECK(all.size() == 2);
  CHECK(all.back().message == "d");

  return failures;
}
