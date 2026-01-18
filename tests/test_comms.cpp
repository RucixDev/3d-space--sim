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

  return failures;
}
