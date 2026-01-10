#include "stellar/sim/Logbook.h"

#include <algorithm>

namespace stellar::sim {

const char* logbookEntryKindName(LogbookEntryKind k) {
  switch (k) {
    case LogbookEntryKind::StarScan: return "Star Scan";
    case LogbookEntryKind::PlanetScan: return "Planet Scan";
    case LogbookEntryKind::StationScan: return "Station Scan";
    case LogbookEntryKind::SignalScan: return "Signal Scan";
    case LogbookEntryKind::AsteroidProspect: return "Asteroid Prospect";
    case LogbookEntryKind::SystemSurveyBonus: return "System Survey Bonus";
    default: return "Unknown";
  }
}

double logbookUnsoldValueCr(const std::vector<LogbookEntry>& entries) {
  double sum = 0.0;
  for (const auto& e : entries) {
    if (e.sold) continue;
    if (e.valueCr <= 0.0) continue;
    sum += e.valueCr;
  }
  return sum;
}

double logbookMarkAllSold(std::vector<LogbookEntry>& entries) {
  double sold = 0.0;
  for (auto& e : entries) {
    if (e.sold) continue;
    if (e.valueCr <= 0.0) continue;
    e.sold = true;
    sold += e.valueCr;
  }
  return sold;
}

double logbookMarkSoldByKeys(std::vector<LogbookEntry>& entries,
                             const std::vector<core::u64>& keys) {
  if (keys.empty() || entries.empty()) return 0.0;

  double sold = 0.0;
  for (auto& e : entries) {
    if (e.sold) continue;
    if (e.valueCr <= 0.0) continue;
    if (std::find(keys.begin(), keys.end(), e.key) == keys.end()) continue;
    e.sold = true;
    sold += e.valueCr;
  }
  return sold;
}

void pruneLogbook(std::vector<LogbookEntry>& entries, std::size_t maxEntries) {
  if (maxEntries < 1) maxEntries = 1;
  if (entries.size() <= maxEntries) return;

  // Prefer pruning SOLD entries first (oldest first), then fall back to oldest overall.
  auto oldestFirst = [](const LogbookEntry& a, const LogbookEntry& b) {
    if (a.discoveredDay != b.discoveredDay) return a.discoveredDay < b.discoveredDay;
    return a.key < b.key;
  };

  // Partition: sold first.
  std::stable_sort(entries.begin(), entries.end(), [&](const LogbookEntry& a, const LogbookEntry& b) {
    if (a.sold != b.sold) return a.sold && !b.sold; // sold entries come first
    return oldestFirst(a, b);
  });

  // Remove excess from the front (oldest sold) while preserving relative order of the remainder.
  const std::size_t excess = entries.size() - maxEntries;
  entries.erase(entries.begin(), entries.begin() + (std::ptrdiff_t)excess);

  // Re-sort to a more natural append-order (oldest->newest) so future pruning remains stable.
  std::stable_sort(entries.begin(), entries.end(), oldestFirst);
}

} // namespace stellar::sim
