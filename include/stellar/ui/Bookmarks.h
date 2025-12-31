#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/Celestial.h"

#include <optional>
#include <string>
#include <vector>

namespace stellar::ui {

// Simple persistent navigation shortcuts.
//
// Bookmarks are stored in the core library (no ImGui dependency) so they can be
// unit-tested and reused by different front-ends.

enum class BookmarkKind : core::u8 {
  System = 0,
  Station = 1,
};

struct Bookmark {
  BookmarkKind kind{BookmarkKind::System};
  sim::SystemId systemId{0};
  sim::StationId stationId{0}; // only meaningful when kind == Station
  std::string label;
};

struct Bookmarks {
  int version{1};
  std::vector<Bookmark> items;
};

// Default config file path (relative to the working directory).
std::string defaultBookmarksPath();

// Default (empty) data.
Bookmarks makeDefaultBookmarks();

// Persistence (simple text format).
bool saveBookmarksToFile(const Bookmarks& bookmarks, const std::string& path);
bool loadBookmarksFromFile(const std::string& path, Bookmarks& out);

// Helpers for managing bookmark sets (dedupe by ids).
std::optional<std::size_t> findSystemBookmarkIndex(const Bookmarks& b, sim::SystemId sysId);
std::optional<std::size_t> findStationBookmarkIndex(const Bookmarks& b, sim::SystemId sysId, sim::StationId stId);

// True if any bookmark targets the given system (either as a system bookmark or
// as a station bookmark within that system).
bool hasAnyInSystem(const Bookmarks& b, sim::SystemId sysId);

// Add/update.
void addOrUpdateSystemBookmark(Bookmarks& b, sim::SystemId sysId, std::string label);
void addOrUpdateStationBookmark(Bookmarks& b, sim::SystemId sysId, sim::StationId stId, std::string label);

// Remove. Returns true if something was removed.
bool removeSystemBookmark(Bookmarks& b, sim::SystemId sysId);
bool removeStationBookmark(Bookmarks& b, sim::SystemId sysId, sim::StationId stId);
bool removeAllInSystem(Bookmarks& b, sim::SystemId sysId);

} // namespace stellar::ui
