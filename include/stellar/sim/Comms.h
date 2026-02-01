#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/Celestial.h" // SystemId / StationId

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Comms / Transmissions
// -----------------------------------------------------------------------------
//
// This is a lightweight, UI-friendly message log designed to make the world feel
// more "alive" (pirate ultimatums, authority scans, mission contracts, etc.)
// without bloating the save format.
//
// Notes:
//  - Messages can be serialized as part of SaveGame (see SaveGame::comms).
//  - Strings may contain optional ui::textfx markup tags. The sim layer does not
//    interpret them; unknown tags remain literal.

enum class CommsChannel : core::u8 {
  System = 0,
  Mission,
  Security,
  Pirate,
  Trade,
  Distress,
  Custom,
};

const char* commsChannelName(CommsChannel c);

struct CommsMessage {
  core::u64 id{0};
  double timeDays{0.0};

  CommsChannel channel{CommsChannel::System};

  // Optional linkage back into the world.
  core::u32 factionId{0};
  SystemId systemId{0};
  StationId stationId{0};
  core::u64 relatedId{0}; // mission id / signal id / etc.

  bool unread{true};
  bool pinned{false};

  // UI-facing strings (may include textfx markup tags).
  std::string from;
  std::string subject;
  std::string body;
};

struct CommsLogParams {
  std::size_t maxMessages{256};
  bool dedupe{true};
};

class CommsLog {
 public:
  core::u64 push(CommsMessage msg, const CommsLogParams& params = {});
  void clear();

  // Replace the entire log contents (primarily for save/load).
  // Preserves existing message ids where possible and updates the internal id allocator.
  void replace(std::vector<CommsMessage> items);

  std::size_t size() const { return items_.size(); }
  std::size_t unreadCount() const;

  void markAllRead();
  bool markRead(core::u64 id, bool read = true);
  bool togglePinned(core::u64 id);

  // Explicitly set pinned state (UI convenience).
  bool markPinned(core::u64 id, bool pinned);

  CommsMessage* find(core::u64 id);
  const CommsMessage* find(core::u64 id) const;

  // Alias used by some UI code (matches naming used elsewhere in the project).
  CommsMessage* findMutable(core::u64 id) { return find(id); }

  std::vector<CommsMessage>& items() { return items_; }
  const std::vector<CommsMessage>& items() const { return items_; }

  // Alias used by some UI code (matches naming used elsewhere in the project).
  std::vector<CommsMessage>& itemsMutable() { return items_; }

 private:
  core::u64 nextId_{1};
  std::vector<CommsMessage> items_;
};

// A small, UI-oriented preview for transient overlays ("incoming transmission").
struct CommsPreview {
  std::string titleMarkup;
  std::string lineMarkup;
};

CommsPreview makeCommsPreview(const CommsMessage& msg, int maxPlainChars = 180);

// -----------------------------------------------------------------------------
// Message composers (pure helpers)
// -----------------------------------------------------------------------------

struct MissionBriefing;

CommsMessage makePirateDemandMessage(double timeDays,
                                    SystemId systemId,
                                    core::u64 pirateGroupId,
                                    std::string_view pirateName,
                                    double requiredValueCr,
                                    double untilDays);

CommsMessage makePoliceBountyDemandMessage(double timeDays,
                                          SystemId systemId,
                                          core::u32 factionId,
                                          std::string_view authorityName,
                                          double bountyCr,
                                          double untilDays);

CommsMessage makePoliceBribeOfferMessage(double timeDays,
                                        SystemId systemId,
                                        core::u32 factionId,
                                        std::string_view authorityName,
                                        double bribeCr,
                                        double fineCr,
                                        std::string_view detail,
                                        double untilDays);

CommsMessage makeMissionBriefingMessage(double timeDays,
                                       SystemId systemId,
                                       StationId stationId,
                                       core::u32 factionId,
                                       core::u64 missionId,
                                       const MissionBriefing& brief,
                                       bool accepted);

} // namespace stellar::sim
