#include "test_harness.h"

#include "stellar/ui/Bookmarks.h"

#include <fstream>

int test_bookmarks() {
  int failures = 0;

  const std::string path = "test_bookmarks_tmp.txt";

  // Roundtrip save/load.
  {
    stellar::ui::Bookmarks b = stellar::ui::makeDefaultBookmarks();
    stellar::ui::addOrUpdateSystemBookmark(b, (stellar::sim::SystemId)1234, "Sol");
    stellar::ui::addOrUpdateStationBookmark(b, (stellar::sim::SystemId)1234, (stellar::sim::StationId)77,
                                            "Sol / Galileo");

    CHECK(stellar::ui::saveBookmarksToFile(b, path));

    stellar::ui::Bookmarks out = stellar::ui::makeDefaultBookmarks();
    CHECK(stellar::ui::loadBookmarksFromFile(path, out));
    CHECK(out.items.size() == 2);
    CHECK(stellar::ui::hasAnyInSystem(out, (stellar::sim::SystemId)1234));
    CHECK(stellar::ui::findSystemBookmarkIndex(out, (stellar::sim::SystemId)1234).has_value());
    CHECK(stellar::ui::findStationBookmarkIndex(out, (stellar::sim::SystemId)1234, (stellar::sim::StationId)77).has_value());
  }

  // Parser resilience: unknown lines + empty label handling.
  {
    std::ofstream f(path, std::ios::out | std::ios::trunc);
    f << "StellarForgeBookmarks 1\n";
    f << "# comment\n";
    f << "unknown_token 1 2 3\n";
    f << "system 555\n"; // missing label -> should become (unnamed)
    f << "station 555 999   \n";
    f.close();

    stellar::ui::Bookmarks out = stellar::ui::makeDefaultBookmarks();
    CHECK(stellar::ui::loadBookmarksFromFile(path, out));
    CHECK(out.items.size() == 2);
    CHECK(out.items[0].label.size() > 0);
    CHECK(out.items[1].label.size() > 0);
  }

  return failures;
}
