#include "test_harness.h"

#include "stellar/sim/GalNet.h"
#include "stellar/sim/Universe.h"

#include <string>
#include <vector>

int test_galnet_digest() {
  int failures = 0;

  stellar::sim::Universe u(/*seed=*/1337);

  const auto stubs = u.queryNearby(stellar::math::Vec3d{0.0, 0.0, 0.0},
                                   /*radiusLy=*/200.0,
                                   /*maxResults=*/16);
  CHECK(stubs.size() >= 3);

  const auto id0 = stubs[0].id;
  const auto id1 = stubs[1].id;
  const auto id2 = stubs[2].id;

  // Intentionally include a duplicate and an invalid id to ensure the digest
  // sanitizes inputs.
  std::vector<stellar::sim::SystemId> ids = {id0, id1, id2, id0, stellar::sim::SystemId{0}};

  stellar::sim::GalNetDigestParams p{};
  p.maxItems = 2;
  p.minSeverity01 = 0.0;
  p.sortBySeverity = true;

  const double timeDays = 123.0;
  const auto r = stellar::sim::makeGalNetDigestBulletin(u, ids, timeDays, p, "Test digest");

  CHECK(r.ok);
  CHECK(r.bulletin.msg.from == "GalNet");
  CHECK(r.bulletin.msg.channel == stellar::sim::CommsChannel::System);
  CHECK(r.bulletin.msg.subject.find("GalNet Digest") != std::string::npos);
  CHECK(r.bulletin.msg.body.find("Test digest") != std::string::npos);
  CHECK(r.bulletin.msg.body.find("Watchlist digest") != std::string::npos);

  // Should mention at least one system name and should be truncated due to maxItems.
  CHECK(r.bulletin.msg.body.find(stubs[0].name) != std::string::npos);
  CHECK(r.bulletin.msg.body.find("not shown") != std::string::npos);

  // Importance should be within [0,1].
  CHECK(r.bulletin.importance01 >= 0.0 && r.bulletin.importance01 <= 1.0);

  return failures;
}
