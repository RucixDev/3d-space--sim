#include "test_harness.h"

#include "stellar/sim/Comms.h"

#include <string>

using namespace stellar;

int test_comms() {
  int failures = 0;

  sim::CommsLog log;
  sim::CommsLogParams params;
  params.maxMessages = 3;
  params.dedupe = true;

  sim::CommsMessage a;
  a.timeDays = 1.0;
  a.channel = sim::CommsChannel::System;
  a.from = "ATC";
  a.subject = "Docking clearance";
  a.body = "Proceed to pad 3.";

  const core::u64 idA1 = log.push(a, params);
  const core::u64 idA2 = log.push(a, params);

  CHECK(idA1 == idA2);
  CHECK(log.items().size() == 1);
  CHECK(log.unreadCount() == 1);

  sim::CommsMessage b;
  b.timeDays = 2.0;
  b.channel = sim::CommsChannel::Pirate;
  b.from = "Pirates";
  b.subject = "Pay up";
  b.body = "Drop cargo.";

  sim::CommsMessage c;
  c.timeDays = 3.0;
  c.channel = sim::CommsChannel::Mission;
  c.from = "Fixer";
  c.subject = "Contract";
  c.body = "Deliver package.";

  sim::CommsMessage d;
  d.timeDays = 4.0;
  d.channel = sim::CommsChannel::Security;
  d.from = "Security";
  d.subject = "Scan";
  d.body = "Submit.";

  log.push(b, params);
  log.push(c, params);
  log.push(d, params);
  CHECK(log.items().size() == 3);

  // Should have pruned the oldest duplicate (ATC) because max=3.
  bool foundATC = false;
  for (const auto& m : log.items()) {
    if (m.from == "ATC") foundATC = true;
  }
  CHECK(!foundATC);

  log.markAllRead();
  CHECK(log.unreadCount() == 0);

  // Preview should include markup tags and be truncated safely.
  const sim::CommsPreview prev = sim::makeCommsPreview(log.items().front(), 24);
  CHECK(!prev.titleMarkup.empty());
  CHECK(!prev.lineMarkup.empty());
  CHECK(prev.lineMarkup.find("[type") != std::string::npos);



  // MarkRead/MarkPinned APIs + pruning should preserve pinned messages.
  {
    sim::CommsLog log2;
    sim::CommsLogParams p;
    p.maxMessages = 2;
    p.dedupe = false;

    sim::CommsMessage m1;
    m1.timeDays = 1.0;
    m1.channel = sim::CommsChannel::System;
    m1.from = "A";
    m1.subject = "S1";
    m1.body = "B1";

    sim::CommsMessage m2;
    m2.timeDays = 2.0;
    m2.channel = sim::CommsChannel::System;
    m2.from = "B";
    m2.subject = "S2";
    m2.body = "B2";

    const core::u64 id1 = log2.push(m1, p);
    const core::u64 id2 = log2.push(m2, p);

    // markRead should flip the unread flag deterministically.
    CHECK(log2.markRead(id2, /*read=*/true));
    CHECK(log2.find(id2) && !log2.find(id2)->unread);
    CHECK(log2.markRead(id2, /*read=*/false));
    CHECK(log2.find(id2) && log2.find(id2)->unread);

    // markPinned should set pinned state.
    CHECK(log2.markPinned(id2, true));
    CHECK(log2.find(id2) && log2.find(id2)->pinned);

    // When capacity is exceeded, the oldest *non-pinned* message should be pruned.
    sim::CommsMessage m3;
    m3.timeDays = 3.0;
    m3.channel = sim::CommsChannel::System;
    m3.from = "C";
    m3.subject = "S3";
    m3.body = "B3";
    (void)log2.push(m3, p);

    CHECK(log2.size() == 2);
    CHECK(log2.find(id2) != nullptr);
    CHECK(log2.find(id1) == nullptr);

    // Unknown ids should return false.
    CHECK(!log2.markRead(999999999ull, true));
    CHECK(!log2.markPinned(999999999ull, true));
  }
  return failures;
}
