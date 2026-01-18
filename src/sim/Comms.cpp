#include "stellar/sim/Comms.h"

#include "stellar/core/Hash.h"
#include "stellar/sim/MissionBriefing.h"
#include "stellar/ui/TextFx.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace stellar::sim {

const char* commsChannelName(CommsChannel c) {
  switch (c) {
    case CommsChannel::System: return "System";
    case CommsChannel::Mission: return "Mission";
    case CommsChannel::Security: return "Security";
    case CommsChannel::Pirate: return "Pirate";
    case CommsChannel::Trade: return "Trade";
    case CommsChannel::Distress: return "Distress";
    case CommsChannel::Custom: return "Custom";
  }
  return "Unknown";
}

static std::string channelColorTag(CommsChannel c) {
  // Hand-picked neon-ish palette, intentionally subtle.
  switch (c) {
    case CommsChannel::Pirate: return "#FF3B3B";
    case CommsChannel::Security: return "#FFB340";
    case CommsChannel::Mission: return "#67B7FF";
    case CommsChannel::Trade: return "#C77DFF";
    case CommsChannel::Distress: return "#55FF8A";
    case CommsChannel::System: return "#B8C0CC";
    case CommsChannel::Custom: return "#E0E0E0";
  }
  return "#E0E0E0";
}

std::size_t CommsLog::unreadCount() const {
  std::size_t n = 0;
  for (const auto& m : items_) {
    if (m.unread) ++n;
  }
  return n;
}

void CommsLog::markAllRead() {
  for (auto& m : items_) m.unread = false;
}

static bool isDuplicate(const CommsMessage& a, const CommsMessage& b) {
  return (a.channel == b.channel) &&
         (a.factionId == b.factionId) &&
         (a.systemId == b.systemId) &&
         (a.stationId == b.stationId) &&
         (a.relatedId == b.relatedId) &&
         (a.from == b.from) &&
         (a.subject == b.subject) &&
         (a.body == b.body);
}

core::u64 CommsLog::push(CommsMessage msg, const CommsLogParams& params) {
  if (params.dedupe && !items_.empty()) {
    const CommsMessage& last = items_.back();
    // If an event fires twice in the same frame we don't want spam.
    if (isDuplicate(last, msg)) {
      return last.id;
    }
  }

  msg.id = nextId_++;
  items_.push_back(std::move(msg));

  // Prune oldest, but try to preserve pinned messages.
  const std::size_t maxN = std::max<std::size_t>(1, params.maxMessages);
  if (items_.size() > maxN) {
    const std::size_t extra = items_.size() - maxN;

    // Remove from the front, skipping pinned messages when possible.
    std::size_t removed = 0;
    for (std::size_t i = 0; i < items_.size() && removed < extra; /*manual*/) {
      if (!items_[i].pinned) {
        items_.erase(items_.begin() + (std::ptrdiff_t)i);
        ++removed;
        continue;
      }
      ++i;
    }

    // If we couldn't remove enough because everything is pinned, fall back to
    // hard trimming from the front.
    if (items_.size() > maxN) {
      items_.erase(items_.begin(), items_.begin() + (std::ptrdiff_t)(items_.size() - maxN));
    }
  }

  return items_.back().id;
}

void CommsLog::clear() {
  items_.clear();
  nextId_ = 1;
}

void CommsLog::replace(std::vector<CommsMessage> items) {
  items_ = std::move(items);

  core::u64 maxId = 0;
  for (const auto& m : items_) {
    maxId = (m.id > maxId) ? m.id : maxId;
  }

  nextId_ = maxId + 1;
  if (nextId_ == 0) nextId_ = 1;
}


CommsMessage* CommsLog::find(core::u64 id) {
  for (auto& m : items_) {
    if (m.id == id) return &m;
  }
  return nullptr;
}

const CommsMessage* CommsLog::find(core::u64 id) const {
  for (const auto& m : items_) {
    if (m.id == id) return &m;
  }
  return nullptr;
}

bool CommsLog::markRead(core::u64 id, bool read) {
  if (CommsMessage* m = find(id)) {
    m->unread = !read;
    return true;
  }
  return false;
}

bool CommsLog::togglePinned(core::u64 id) {
  if (CommsMessage* m = find(id)) {
    m->pinned = !m->pinned;
    return true;
  }
  return false;
}

static std::string firstLine(const std::string& s) {
  const std::size_t n = s.find_first_of("\r\n");
  if (n == std::string::npos) return s;
  return s.substr(0, n);
}

static std::string truncateUtf8ByGlyphs(std::string_view plainUtf8, int maxGlyphs) {
  if (maxGlyphs <= 0) return {};
  // Conservative: iterate bytes and stop when we reach maxGlyphs glyphs.
  int glyphs = 0;
  std::size_t i = 0;
  while (i < plainUtf8.size() && glyphs < maxGlyphs) {
    unsigned char c = (unsigned char)plainUtf8[i];
    std::size_t adv = 1;
    if ((c & 0x80) == 0x00) adv = 1;
    else if ((c & 0xE0) == 0xC0) adv = 2;
    else if ((c & 0xF0) == 0xE0) adv = 3;
    else if ((c & 0xF8) == 0xF0) adv = 4;
    if (i + adv > plainUtf8.size()) adv = 1;
    i += adv;
    ++glyphs;
  }
  return std::string(plainUtf8.substr(0, i));
}

CommsPreview makeCommsPreview(const CommsMessage& msg, int maxPlainChars) {
  CommsPreview out;
  maxPlainChars = std::clamp(maxPlainChars, 32, 512);

  const std::string color = channelColorTag(msg.channel);
  const std::string subjPlain = ui::textfx::stripMarkup(msg.subject);
  const std::string fromPlain = ui::textfx::stripMarkup(msg.from);

  // Title: [CHANNEL] Subject
  {
    std::ostringstream ss;
    ss << "[color " << color << "]";
    ss << "[" << commsChannelName(msg.channel) << "]";
    ss << "[/color] ";
    ss << subjPlain;
    out.titleMarkup = ss.str();
  }

  // Line: From: (first line of message body), typewriter animated.
  const std::string bodyPlainFirst = firstLine(ui::textfx::stripMarkup(msg.body));
  std::string clipped = truncateUtf8ByGlyphs(bodyPlainFirst, maxPlainChars);
  if ((int)ui::textfx::glyphCountUtf8(bodyPlainFirst) > (int)ui::textfx::glyphCountUtf8(clipped)) {
    clipped += "…";
  }

  {
    std::ostringstream ss;
    ss << "[type cps=48 fade=0.06]";
    if (!fromPlain.empty()) {
      ss << "[color " << color << "]" << fromPlain << ": [/color]";
    }
    ss << clipped;
    ss << "[/type]";
    out.lineMarkup = ss.str();
  }

  return out;
}

CommsMessage makePirateDemandMessage(double timeDays,
                                    SystemId systemId,
                                    core::u64 pirateGroupId,
                                    std::string_view pirateName,
                                    double requiredValueCr,
                                    double untilDays) {
  (void)pirateGroupId;
  CommsMessage m;
  m.timeDays = timeDays;
  m.channel = CommsChannel::Pirate;
  m.systemId = systemId;
  m.relatedId = pirateGroupId;
  m.from = std::string(pirateName.empty() ? "Unknown Raiders" : pirateName);
  m.subject = "[shake amp=1.6 rate=18][color #FF3B3B]ULTIMATUM[/color][/shake]";

  const double secLeft = std::max(0.0, (untilDays - timeDays) * 24.0 * 3600.0);
  const int secInt = (int)std::llround(secLeft);

  std::ostringstream body;
  body << "[type cps=60 fade=0.05]";
  body << "We have you on scope. Jettison cargo worth ~" << (int)std::llround(std::max(0.0, requiredValueCr)) << " cr.\\n";
  body << "Time remaining: " << secInt << " s.\\n\\n";
  body << "Tip: Ship → Cargo → Jettison (or vent a few pods quickly).";
  body << "[/type]";
  m.body = body.str();
  return m;
}

CommsMessage makePoliceBountyDemandMessage(double timeDays,
                                          SystemId systemId,
                                          core::u32 factionId,
                                          std::string_view authorityName,
                                          double bountyCr,
                                          double untilDays) {
  CommsMessage m;
  m.timeDays = timeDays;
  m.channel = CommsChannel::Security;
  m.systemId = systemId;
  m.factionId = factionId;
  m.from = std::string(authorityName.empty() ? "System Authority" : authorityName);
  m.subject = "[color #FFB340]ENFORCEMENT ORDER[/color]";

  const double secLeft = std::max(0.0, (untilDays - timeDays) * 24.0 * 3600.0);
  const int secInt = (int)std::llround(secLeft);

  std::ostringstream body;
  body << "[type cps=56 fade=0.05]";
  body << "You are flagged with an outstanding bounty of " << (int)std::llround(std::max(0.0, bountyCr)) << " cr.\\n";
  body << "Submit to inspection immediately.\\n";
  body << "Deadline: " << secInt << " s.";
  body << "[/type]";
  m.body = body.str();
  return m;
}

CommsMessage makePoliceBribeOfferMessage(double timeDays,
                                        SystemId systemId,
                                        core::u32 factionId,
                                        std::string_view authorityName,
                                        double bribeCr,
                                        double fineCr,
                                        std::string_view detail,
                                        double untilDays) {
  CommsMessage m;
  m.timeDays = timeDays;
  m.channel = CommsChannel::Security;
  m.systemId = systemId;
  m.factionId = factionId;
  m.from = std::string(authorityName.empty() ? "System Authority" : authorityName);
  m.subject = "[color #FFB340]DISCRETIONARY RESOLUTION[/color]";

  const double secLeft = std::max(0.0, (untilDays - timeDays) * 24.0 * 3600.0);
  const int secInt = (int)std::llround(secLeft);

  std::ostringstream body;
  body << "[type cps=56 fade=0.05]";
  body << "Scan result: " << ui::textfx::stripMarkup(std::string(detail)).c_str() << "\\n\\n";
  body << "Option A: Pay bribe " << (int)std::llround(std::max(0.0, bribeCr)) << " cr.\\n";
  body << "Option B: Comply (fine ~" << (int)std::llround(std::max(0.0, fineCr)) << " cr, cargo may be seized).\\n\\n";
  body << "Decision window: " << secInt << " s.";
  body << "[/type]";
  m.body = body.str();
  return m;
}

CommsMessage makeMissionBriefingMessage(double timeDays,
                                       SystemId systemId,
                                       StationId stationId,
                                       core::u32 factionId,
                                       core::u64 missionId,
                                       const MissionBriefing& brief,
                                       bool accepted) {
  CommsMessage m;
  m.timeDays = timeDays;
  m.channel = CommsChannel::Mission;
  m.systemId = systemId;
  m.stationId = stationId;
  m.factionId = factionId;
  m.relatedId = missionId;

  const std::string contactPlain = ui::textfx::stripMarkup(brief.contactName);
  m.from = contactPlain.empty() ? "Contract Broker" : contactPlain;

  // Keep the subject short; titleMarkup can be long.
  if (accepted) {
    m.subject = "[color #67B7FF]CONTRACT ACCEPTED[/color]  " + ui::textfx::stripMarkup(brief.contractCode);
  } else {
    m.subject = "[color #67B7FF]CONTRACT BRIEFING[/color]  " + ui::textfx::stripMarkup(brief.contractCode);
  }

  std::ostringstream body;
  body << "[type cps=72 fade=0.04]";
  body << brief.titleMarkup << "\\n\\n";
  body << brief.synopsisMarkup << "\\n\\n";
  body << "[color #9BA6B4]Contact:[/color] " << brief.contactName << "  ";
  body << "[color #9BA6B4]Code:[/color] " << brief.contractCode << "\\n";
  body << "[color #9BA6B4]Risk:[/color] " << riskTierName(brief.risk.overall01) << " (";
  body << (int)std::llround(brief.risk.overall01 * 100.0) << "%)\\n\\n";
  body << "[color #9BA6B4]Details:[/color]\\n";
  for (const auto& line : brief.bulletsMarkup) {
    body << "• " << line << "\\n";
  }
  body << "[/type]";
  m.body = body.str();
  return m;
}

} // namespace stellar::sim
