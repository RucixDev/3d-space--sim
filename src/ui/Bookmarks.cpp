#include "stellar/ui/Bookmarks.h"

#include "stellar/core/Log.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>

namespace stellar::ui {

static std::string lowerAscii(std::string s) {
  for (char& c : s) c = (char)std::tolower((unsigned char)c);
  return s;
}

static std::string trimAscii(std::string s) {
  std::size_t a = 0;
  while (a < s.size() && std::isspace((unsigned char)s[a])) ++a;
  std::size_t b = s.size();
  while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
  return s.substr(a, b - a);
}

static std::string sanitizeLabel(std::string s) {
  for (char& c : s) {
    if (c == '\n' || c == '\r' || c == '\t') c = ' ';
  }
  s = trimAscii(s);
  if (s.empty()) s = "(unnamed)";
  return s;
}

std::string defaultBookmarksPath() {
  return "bookmarks.txt";
}

Bookmarks makeDefaultBookmarks() {
  Bookmarks b{};
  b.version = 1;
  b.items.clear();
  return b;
}

bool saveBookmarksToFile(const Bookmarks& bookmarks, const std::string& path) {
  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Warn, "Bookmarks: failed to open for writing: " + path);
    return false;
  }

  f << "StellarForgeBookmarks " << bookmarks.version << "\n";
  for (const auto& bm : bookmarks.items) {
    if (bm.kind == BookmarkKind::System) {
      f << "system " << (unsigned long long)bm.systemId << " " << sanitizeLabel(bm.label) << "\n";
    } else {
      f << "station " << (unsigned long long)bm.systemId << " " << (unsigned long long)bm.stationId
        << " " << sanitizeLabel(bm.label) << "\n";
    }
  }
  return true;
}

static void dedupeKeepLast(Bookmarks& b) {
  // Keep the last entry for each (kind, ids) key.
  std::vector<Bookmark> out;
  out.reserve(b.items.size());

  auto keyEq = [](const Bookmark& a, const Bookmark& c) {
    return a.kind == c.kind && a.systemId == c.systemId && a.stationId == c.stationId;
  };

  for (const auto& bm : b.items) {
    bool replaced = false;
    for (auto& o : out) {
      if (keyEq(o, bm)) {
        o = bm;
        replaced = true;
        break;
      }
    }
    if (!replaced) out.push_back(bm);
  }
  b.items = std::move(out);
}

bool loadBookmarksFromFile(const std::string& path, Bookmarks& out) {
  std::ifstream f(path);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Debug, "Bookmarks: file not found: " + path);
    return false;
  }

  std::string header;
  int version = 0;
  if (!(f >> header >> version)) {
    stellar::core::log(stellar::core::LogLevel::Warn, "Bookmarks: failed to read header");
    return false;
  }
  if (header != "StellarForgeBookmarks") {
    stellar::core::log(stellar::core::LogLevel::Warn, "Bookmarks: bad header");
    return false;
  }

  Bookmarks loaded = makeDefaultBookmarks();
  loaded.version = version;

  std::string line;
  std::getline(f, line); // consume remainder
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream ss(line);
    std::string tok;
    if (!(ss >> tok)) continue;
    tok = lowerAscii(tok);

    if (tok == "system") {
      unsigned long long sysId = 0;
      if (!(ss >> sysId)) continue;
      std::string label;
      std::getline(ss, label);
      label = sanitizeLabel(label);

      if (sysId == 0ULL) continue;
      Bookmark bm;
      bm.kind = BookmarkKind::System;
      bm.systemId = (sim::SystemId)sysId;
      bm.stationId = 0;
      bm.label = std::move(label);
      loaded.items.push_back(std::move(bm));
      continue;
    }

    if (tok == "station") {
      unsigned long long sysId = 0;
      unsigned long long stId = 0;
      if (!(ss >> sysId >> stId)) continue;
      std::string label;
      std::getline(ss, label);
      label = sanitizeLabel(label);

      if (sysId == 0ULL || stId == 0ULL) continue;
      Bookmark bm;
      bm.kind = BookmarkKind::Station;
      bm.systemId = (sim::SystemId)sysId;
      bm.stationId = (sim::StationId)stId;
      bm.label = std::move(label);
      loaded.items.push_back(std::move(bm));
      continue;
    }

    // Unknown line type: ignore for forward-compat.
  }

  dedupeKeepLast(loaded);
  out = std::move(loaded);
  return true;
}

std::optional<std::size_t> findSystemBookmarkIndex(const Bookmarks& b, sim::SystemId sysId) {
  for (std::size_t i = 0; i < b.items.size(); ++i) {
    const auto& bm = b.items[i];
    if (bm.kind == BookmarkKind::System && bm.systemId == sysId) return i;
  }
  return std::nullopt;
}

std::optional<std::size_t> findStationBookmarkIndex(const Bookmarks& b, sim::SystemId sysId, sim::StationId stId) {
  for (std::size_t i = 0; i < b.items.size(); ++i) {
    const auto& bm = b.items[i];
    if (bm.kind == BookmarkKind::Station && bm.systemId == sysId && bm.stationId == stId) return i;
  }
  return std::nullopt;
}

bool hasAnyInSystem(const Bookmarks& b, sim::SystemId sysId) {
  for (const auto& bm : b.items) {
    if (bm.systemId == sysId) return true;
  }
  return false;
}

void addOrUpdateSystemBookmark(Bookmarks& b, sim::SystemId sysId, std::string label) {
  label = sanitizeLabel(std::move(label));

  for (auto& bm : b.items) {
    if (bm.kind == BookmarkKind::System && bm.systemId == sysId) {
      bm.label = std::move(label);
      return;
    }
  }

  Bookmark bm;
  bm.kind = BookmarkKind::System;
  bm.systemId = sysId;
  bm.stationId = 0;
  bm.label = std::move(label);
  b.items.push_back(std::move(bm));
}

void addOrUpdateStationBookmark(Bookmarks& b, sim::SystemId sysId, sim::StationId stId, std::string label) {
  label = sanitizeLabel(std::move(label));
  for (auto& bm : b.items) {
    if (bm.kind == BookmarkKind::Station && bm.systemId == sysId && bm.stationId == stId) {
      bm.label = std::move(label);
      return;
    }
  }

  Bookmark bm;
  bm.kind = BookmarkKind::Station;
  bm.systemId = sysId;
  bm.stationId = stId;
  bm.label = std::move(label);
  b.items.push_back(std::move(bm));
}

bool removeSystemBookmark(Bookmarks& b, sim::SystemId sysId) {
  const auto it = std::remove_if(b.items.begin(), b.items.end(), [&](const Bookmark& bm) {
    return bm.kind == BookmarkKind::System && bm.systemId == sysId;
  });
  const bool removed = (it != b.items.end());
  b.items.erase(it, b.items.end());
  return removed;
}

bool removeStationBookmark(Bookmarks& b, sim::SystemId sysId, sim::StationId stId) {
  const auto it = std::remove_if(b.items.begin(), b.items.end(), [&](const Bookmark& bm) {
    return bm.kind == BookmarkKind::Station && bm.systemId == sysId && bm.stationId == stId;
  });
  const bool removed = (it != b.items.end());
  b.items.erase(it, b.items.end());
  return removed;
}

bool removeAllInSystem(Bookmarks& b, sim::SystemId sysId) {
  const auto it = std::remove_if(b.items.begin(), b.items.end(), [&](const Bookmark& bm) {
    return bm.systemId == sysId;
  });
  const bool removed = (it != b.items.end());
  b.items.erase(it, b.items.end());
  return removed;
}

} // namespace stellar::ui
