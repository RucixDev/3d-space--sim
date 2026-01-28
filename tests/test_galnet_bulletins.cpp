#include "test_harness.h"

#include "stellar/sim/GalNet.h"

#include <string>

int test_galnet_bulletins() {
  int failures = 0;

  // --- Bulletin composition -------------------------------------------------
  stellar::sim::StarSystem sys{};
  sys.stub.id = 42;
  sys.stub.name = "Arcadia";

  stellar::sim::SystemConditionsSnapshot snap{};
  snap.systemId = sys.stub.id;
  snap.base.controllingFactionId = 7;
  snap.base.security01 = 0.30;
  snap.base.piracy01 = 0.70;
  snap.base.traffic01 = 0.55;

  snap.effective = snap.base;
  snap.effective.security01 = 0.45;
  snap.effective.piracy01 = 0.60;
  snap.effective.traffic01 = 0.65;

  snap.event.active = true;
  snap.event.kind = stellar::sim::SystemEventKind::TradeBoom;
  snap.event.systemId = sys.stub.id;
  snap.event.startDay = 10.0;
  snap.event.endDay = 16.0;
  snap.event.severity01 = 0.80;

  {
    const double timeDays = 12.0;
    const auto r = stellar::sim::makeGalNetBulletin(sys,
                                                    snap,
                                                    timeDays,
                                                    "Zenith Combine",
                                                    "Manual bulletin",
                                                    /*allowWhenNoEvent=*/false);

    CHECK(r.ok);
    CHECK(r.bulletin.msg.channel == stellar::sim::CommsChannel::Trade);
    CHECK(r.bulletin.hasActiveEvent);
    CHECK(r.bulletin.importance01 > 0.79 && r.bulletin.importance01 < 0.81);

    CHECK(r.bulletin.msg.subject.find("GalNet: Trade Boom") != std::string::npos);
    CHECK(r.bulletin.msg.subject.find("Arcadia") != std::string::npos);

    CHECK(r.bulletin.msg.body.find("Manual bulletin") != std::string::npos);
    CHECK(r.bulletin.msg.body.find("System: Arcadia") != std::string::npos);
    CHECK(r.bulletin.msg.body.find("Controlling faction: Zenith Combine") != std::string::npos);
    CHECK(r.bulletin.msg.body.find("Severity: 80%") != std::string::npos);
    CHECK(r.bulletin.msg.body.find("Ends in: 4d") != std::string::npos);
    CHECK(r.bulletin.msg.body.find("Security: 45% (+15%)") != std::string::npos);
    CHECK(r.bulletin.msg.body.find("Piracy: 60% (-10%)") != std::string::npos);
    CHECK(r.bulletin.msg.body.find("Traffic: 65% (+10%)") != std::string::npos);
    CHECK(r.bulletin.msg.body.find("Market impact:") != std::string::npos);
  }

  // No-event bulletin (used for "event ended" / cycle status updates).
  {
    snap.event.active = false;
    snap.event.kind = stellar::sim::SystemEventKind::None;
    snap.event.startDay = 16.0;
    snap.event.endDay = 22.0;

    const auto r = stellar::sim::makeGalNetBulletin(sys,
                                                    snap,
                                                    /*timeDays=*/17.0,
                                                    "Zenith Combine",
                                                    /*contextTag=*/{},
                                                    /*allowWhenNoEvent=*/true);
    CHECK(r.ok);
    CHECK(!r.bulletin.hasActiveEvent);
    CHECK(r.bulletin.msg.subject.find("System status update") != std::string::npos);
    CHECK(r.bulletin.msg.body.find("Event: None") != std::string::npos);
  }

  // --- Auto-broadcast decision logic ---------------------------------------
  {
    stellar::sim::GalNetAnnounceState st{};

    stellar::sim::SystemEvent ev{};
    ev.startDay = 10.0;
    ev.endDay = 16.0;
    ev.active = true;
    ev.severity01 = 0.50;
    ev.kind = stellar::sim::SystemEventKind::PirateRaid;

    const auto d1 = stellar::sim::galNetMaybeAutoBroadcast(
        st, ev, /*minSeverity01=*/0.25, /*autoEnabled=*/true, /*broadcastEventEnds=*/true);
    CHECK(d1.shouldPublish);
    CHECK(!d1.allowWhenNoEvent);
    CHECK(st.lastActiveEventKind == stellar::sim::SystemEventKind::PirateRaid);
    CHECK(st.lastActiveEventSeverity01 > 0.49 && st.lastActiveEventSeverity01 < 0.51);

    // Same cycle: no publish.
    const auto d2 = stellar::sim::galNetMaybeAutoBroadcast(
        st, ev, /*minSeverity01=*/0.25, /*autoEnabled=*/true, /*broadcastEventEnds=*/true);
    CHECK(!d2.shouldPublish);

    // Event ends: emit a status update if broadcastEventEnds is enabled.
    stellar::sim::SystemEvent evEnd = ev;
    evEnd.active = false;
    evEnd.startDay = 16.0;
    evEnd.endDay = 22.0;

    const auto d3 = stellar::sim::galNetMaybeAutoBroadcast(
        st, evEnd, /*minSeverity01=*/0.25, /*autoEnabled=*/true, /*broadcastEventEnds=*/true);
    CHECK(d3.shouldPublish);
    CHECK(d3.allowWhenNoEvent);
    // Ended-cycle: should preserve the last active event kind for UI labeling.
    CHECK(st.lastActiveEventKind == stellar::sim::SystemEventKind::PirateRaid);

    // Auto disabled: should not publish, but should still advance state.
    stellar::sim::GalNetAnnounceState st2{};
    const auto d4 = stellar::sim::galNetMaybeAutoBroadcast(
        st2, ev, /*minSeverity01=*/0.25, /*autoEnabled=*/false, /*broadcastEventEnds=*/true);
    CHECK(!d4.shouldPublish);
    CHECK(st2.lastCycleStartDay == ev.startDay);
    CHECK(st2.lastActiveEventKind == stellar::sim::SystemEventKind::PirateRaid);

    // Enabling later should not retro-spam the same cycle.
    const auto d5 = stellar::sim::galNetMaybeAutoBroadcast(
        st2, ev, /*minSeverity01=*/0.25, /*autoEnabled=*/true, /*broadcastEventEnds=*/true);
    CHECK(!d5.shouldPublish);
  }

  return failures;
}
