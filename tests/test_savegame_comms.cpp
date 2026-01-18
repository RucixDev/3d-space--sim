#include "test_harness.h"

#include "stellar/sim/SaveGame.h"

#include <cstdio>
#include <string>

using namespace stellar;

static bool msgsEqual(const sim::CommsMessage& a, const sim::CommsMessage& b) {
  return a.id == b.id &&
         a.timeDays == b.timeDays &&
         a.channel == b.channel &&
         a.factionId == b.factionId &&
         a.systemId == b.systemId &&
         a.stationId == b.stationId &&
         a.relatedId == b.relatedId &&
         a.unread == b.unread &&
         a.pinned == b.pinned &&
         a.from == b.from &&
         a.subject == b.subject &&
         a.body == b.body;
}

int test_savegame_comms() {
  int failures = 0;

  sim::SaveGame s{};
  s.seed = 42ull;
  s.timeDays = 123.456;
  s.currentSystem = 77;
  s.dockedStation = 5;

  // Two messages: one with empty fields (exercises the "~" token), one with markup and newlines.
  {
    sim::CommsMessage m{};
    m.id = 1;
    m.timeDays = 100.0;
    m.channel = sim::CommsChannel::System;
    m.factionId = 0;
    m.systemId = 77;
    m.stationId = 0;
    m.relatedId = 0;
    m.unread = true;
    m.pinned = false;
    m.from = "";
    m.subject = "";
    m.body = "";
    s.comms.push_back(m);
  }

  {
    sim::CommsMessage m{};
    m.id = 2;
    m.timeDays = 123.0;
    m.channel = sim::CommsChannel::Pirate;
    m.factionId = 9;
    m.systemId = 77;
    m.stationId = 5;
    m.relatedId = 7777ull;
    m.unread = false;
    m.pinned = true;
    m.from = "[color #ff6666]Red Fang[/color]";
    m.subject = "[wave amp=3 freq=0.2 speed=1.8]PAY UP[/wave]";
    m.body = "Line1\nLine2\r\n[grad #ff00aa #00ccff]Nebula[/grad]";
    s.comms.push_back(m);
  }

  const std::string path = "__test_savegame_comms.sav";

  CHECK(sim::saveToFile(s, path));

  sim::SaveGame out{};
  CHECK(sim::loadFromFile(path, out));

  CHECK(out.comms.size() == s.comms.size());
  if (out.comms.size() == s.comms.size()) {
    for (std::size_t i = 0; i < s.comms.size(); ++i) {
      CHECK(msgsEqual(out.comms[i], s.comms[i]));
    }
  }

  // Cleanup (best-effort).
  std::remove(path.c_str());

  return failures;
}
