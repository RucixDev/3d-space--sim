#include "CommsWindow.h"

#include "stellar/sim/Universe.h"
#include "stellar/ui/FuzzySearch.h"
#include "stellar/ui/TextFx.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>

#include <imgui.h>

namespace stellar::game {

namespace {
static ui::textfx::ProgramCache g_textFxCache(1024);
} // namespace

static std::string formatSimTime(double timeDays) {
  const double totalSec = std::max(0.0, timeDays) * 24.0 * 3600.0;
  const long long sec = (long long)std::llround(totalSec);
  const long long day = sec / 86400;
  const long long rem = sec % 86400;
  const int hh = (int)(rem / 3600);
  const int mm = (int)((rem % 3600) / 60);
  const int ss = (int)(rem % 60);
  char buf[64];
  std::snprintf(buf, sizeof(buf), "D%lld %02d:%02d:%02d", day, hh, mm, ss);
  return buf;
}

void enqueueCommsOverlay(CommsOverlayState& ov, core::u64 messageId) {
  if (messageId == 0) return;
  ov.queue.push_back(messageId);
}

void tickAndDrawCommsOverlay(CommsOverlayState& ov, const sim::CommsLog& log, double timeRealSec) {
  if (!ov.enabled) return;

  // Activate next queued message.
  if (ov.activeId == 0 && !ov.queue.empty()) {
    ov.activeId = ov.queue.front();
    ov.queue.erase(ov.queue.begin());
    ov.activeStartSec = timeRealSec;
    ov.activeUntilSec = timeRealSec + ov.baseHoldSec;
  }

  if (ov.activeId == 0) return;

  const sim::CommsMessage* msg = log.find(ov.activeId);
  if (!msg) {
    ov.activeId = 0;
    return;
  }

  // Recompute lifetime as a function of preview length (feels nicer).
  const sim::CommsPreview prev = sim::makeCommsPreview(*msg);
  const ui::textfx::Program& prevLineProg = g_textFxCache.get(prev.lineMarkup);
  const int glyphs = prevLineProg.glyphCount;
  const double dynamicHold = std::clamp(2.8 + 0.020 * (double)glyphs, 4.5, 9.0);
  ov.activeUntilSec = ov.activeStartSec + std::max(ov.baseHoldSec, dynamicHold);

  // Expire.
  if (timeRealSec > ov.activeUntilSec) {
    ov.activeId = 0;
    return;
  }

  // Fade-out.
  float alpha = 1.0f;
  if (ov.fadeOutSec > 0.01) {
    const double tLeft = ov.activeUntilSec - timeRealSec;
    if (tLeft < ov.fadeOutSec) {
      alpha = (float)std::clamp(tLeft / ov.fadeOutSec, 0.0, 1.0);
    }
  }

  ImGuiViewport* vp = ImGui::GetMainViewport();
  if (!vp) return;

  ui::textfx::DrawParams dp;
  dp.wrapWidthPx = std::max(10.0f, ov.widthPx - ov.padPx * 2.0f);

  const ui::textfx::Program& progTitle = g_textFxCache.get(prev.titleMarkup);
  const ui::textfx::Program& progLine = prevLineProg;

  const ImVec2 titleSz = ui::textfx::CalcSize(progTitle, dp);
  const ImVec2 lineSz  = ui::textfx::CalcSize(progLine, dp);

  const float w = ov.widthPx;
  const float h = ov.padPx * 2.0f + titleSz.y + 6.0f + lineSz.y;

  const ImVec2 p0(vp->WorkPos.x + vp->WorkSize.x - w - ov.marginPx,
                  vp->WorkPos.y + ov.marginPx);
  const ImVec2 p1(p0.x + w, p0.y + h);

  ImDrawList* draw = ImGui::GetForegroundDrawList(vp);
  if (!draw) return;

  const ImU32 bg = IM_COL32(0, 0, 0, (int)std::llround(165.0f * alpha));
  const ImU32 bd = IM_COL32(180, 200, 235, (int)std::llround(55.0f * alpha));

  draw->AddRectFilled(p0, p1, bg, 6.0f);
  draw->AddRect(p0, p1, bd, 6.0f, 0, 1.2f);

  const float tLocal = (float)(timeRealSec - ov.activeStartSec);

  dp.baseColor = IM_COL32(235, 240, 255, (int)std::llround(235.0f * alpha));
  ui::textfx::Draw(draw, ImVec2(p0.x + ov.padPx, p0.y + ov.padPx), progTitle, tLocal, dp);

  dp.baseColor = IM_COL32(235, 240, 255, (int)std::llround(230.0f * alpha));
  ui::textfx::Draw(draw, ImVec2(p0.x + ov.padPx, p0.y + ov.padPx + titleSz.y + 6.0f), progLine, tLocal, dp);
}

static std::string systemName(sim::Universe* u, sim::SystemId sysId) {
  if (!u || sysId == 0) return "(unknown)";
  return u->getSystem(sysId).stub.name;
}

static std::string stationName(sim::Universe* u, sim::SystemId sysId, sim::StationId stId) {
  if (!u || sysId == 0 || stId == 0) return "(unknown)";
  const auto& sys = u->getSystem(sysId);
  for (const auto& st : sys.stations) {
    if (st.id == stId) return st.name;
  }
  return "Station " + std::to_string(stId);
}

static void drawFxWrapped(std::string_view markup, float wrapWidthPx, float timeSec, ImU32 baseColor = 0, core::u64 seed = 0) {
  ui::textfx::DrawParams dp;
  dp.wrapWidthPx = wrapWidthPx;
  dp.baseColor = baseColor;
  dp.seed = seed;

  const ui::textfx::Program& prog = g_textFxCache.get(markup);
  ui::textfx::Draw(ImGui::GetWindowDrawList(), ImGui::GetCursorScreenPos(), prog, timeSec, dp);
  const ImVec2 sz = ui::textfx::CalcSize(prog, dp);
  ImGui::Dummy(sz);
}

void drawCommsWindow(CommsWindowState& st, CommsOverlayState& overlay, CommsWindowContext& ctx) {
  if (!st.open) return;
  if (!ctx.log) return;

  auto& log = *ctx.log;

  const std::size_t unread = log.unreadCount();
  std::string title = "Comms";
  if (unread > 0) title += "  (" + std::to_string(unread) + " unread)";
  title += "###Comms";

  ImGui::SetNextWindowSize(ImVec2(900, 560), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin(title.c_str(), &st.open)) {
    ImGui::End();
    return;
  }

  // Toolbar
  if (ImGui::Button("Mark all read")) {
    log.markAllRead();
  }
  ImGui::SameLine();
  ImGui::Checkbox("Unread only", &st.unreadOnly);
  ImGui::SameLine();
  ImGui::Checkbox("Pinned first", &st.pinnedFirst);
  ImGui::SameLine();
  ImGui::Checkbox("Newest first", &st.newestFirst);
  ImGui::SameLine();
  ImGui::Checkbox("Wrap", &st.wrapBody);
  ImGui::SameLine();
  ImGui::Checkbox("Overlay", &overlay.enabled);

  // Channel filter
  ImGui::SameLine();
  ImGui::SetNextItemWidth(140.0f);
  if (ImGui::BeginCombo("##channel", (st.channelFilter < 0) ? "All channels" : sim::commsChannelName((sim::CommsChannel)st.channelFilter))) {
    if (ImGui::Selectable("All channels", st.channelFilter < 0)) st.channelFilter = -1;
    for (int i = 0; i <= (int)sim::CommsChannel::Custom; ++i) {
      const auto c = (sim::CommsChannel)i;
      const bool sel = (st.channelFilter == i);
      if (ImGui::Selectable(sim::commsChannelName(c), sel)) st.channelFilter = i;
      if (sel) ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }

  // Text filter
  ImGui::SetNextItemWidth(-1);
  ImGui::InputTextWithHint("##filter", "filter (fuzzy)", st.filter, sizeof(st.filter));

  ImGui::Separator();

  // Two-pane layout
  const float availW = ImGui::GetContentRegionAvail().x;
  const float listW = std::max(300.0f, std::min(420.0f, availW * 0.42f));

  ImGui::BeginChild("##comms_list", ImVec2(listW, 0), true);

  const auto& items = log.items();
  std::vector<int> idx(items.size());
  std::iota(idx.begin(), idx.end(), 0);

  // Filter + score
  struct ScoredIdx { int i; int score; };
  std::vector<ScoredIdx> scored;
  scored.reserve(items.size());
  const std::string query(st.filter);

  for (int i : idx) {
    const auto& m = items[(std::size_t)i];
    if (st.unreadOnly && !m.unread) continue;
    if (st.channelFilter >= 0 && (int)m.channel != st.channelFilter) continue;

    int score = 0;
    if (!query.empty()) {
      const std::string hay = g_textFxCache.get(m.from).plain + " "
                            + g_textFxCache.get(m.subject).plain + " "
                            + g_textFxCache.get(m.body).plain;
      score = ui::fuzzyMatchScore(query, hay);
      if (score < 0) continue;
    }
    scored.push_back({i, score});
  }

  // Sort
  std::sort(scored.begin(), scored.end(), [&](const ScoredIdx& a, const ScoredIdx& b) {
    const auto& ma = items[(std::size_t)a.i];
    const auto& mb = items[(std::size_t)b.i];

    // Pinned first
    if (st.pinnedFirst && ma.pinned != mb.pinned) return ma.pinned > mb.pinned;

    // Query relevance
    if (!query.empty() && a.score != b.score) return a.score > b.score;

    // Time
    if (st.newestFirst) return ma.timeDays > mb.timeDays;
    return ma.timeDays < mb.timeDays;
  });

  // List header
  ImGui::TextDisabled("%zu messages", scored.size());
  ImGui::Separator();

  for (const auto& si : scored) {
    const auto& m = items[(std::size_t)si.i];

    const bool selected = (m.id == st.selectedId);
    const ui::textfx::Program& subjProg = g_textFxCache.get(m.subject);
    std::string label = subjProg.plain;
    if (label.empty()) label = "(no subject)";
    label += "##" + std::to_string(m.id);

    if (m.unread) {
      ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(255, 255, 255, 245));
    } else {
      ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(200, 210, 225, 210));
    }

    if (ImGui::Selectable(label.c_str(), selected)) {
      st.selectedId = m.id;
      st.selectedOpenedSec = ctx.timeRealSec;
      if (st.autoMarkReadOnSelect) {
        log.markRead(m.id, true);
      }
    }
    ImGui::PopStyleColor();

    // Mini meta line
    ImGui::SameLine();
    ImGui::TextDisabled("%s", m.pinned ? "[PIN]" : "");
  }

  ImGui::EndChild();
  ImGui::SameLine();

  ImGui::BeginChild("##comms_body", ImVec2(0, 0), true);

  sim::CommsMessage* sel = (st.selectedId != 0) ? log.find(st.selectedId) : nullptr;
  if (!sel) {
    ImGui::TextDisabled("Select a message");
    ImGui::EndChild();
    ImGui::End();
    return;
  }

  const auto& m = *sel;
  const ui::textfx::Program& subjProgSel = g_textFxCache.get(m.subject);
  ImGui::TextUnformatted(subjProgSel.plain.c_str());
  const ui::textfx::Program& fromProgSel = g_textFxCache.get(m.from);
  ImGui::TextDisabled("%s  ·  %s  ·  %s", fromProgSel.plain.c_str(),
                       sim::commsChannelName(m.channel),
                       formatSimTime(m.timeDays).c_str());
  if (m.systemId != 0) {
    ImGui::TextDisabled("Location: %s / %s",
                        systemName(ctx.universe, m.systemId).c_str(),
                        stationName(ctx.universe, m.systemId, m.stationId).c_str());
  }

  // Action bar
  if (ImGui::Button(m.unread ? "Mark read" : "Mark unread")) {
    log.markRead(m.id, m.unread);
  }
  ImGui::SameLine();
  if (ImGui::Button(m.pinned ? "Unpin" : "Pin")) {
    log.togglePinned(m.id);
  }
  ImGui::SameLine();
  if (ImGui::Button("Copy plain")) {
    const ui::textfx::Program& bodyPlain = g_textFxCache.get(m.body);
    ImGui::SetClipboardText(bodyPlain.plain.c_str());
    if (ctx.toast) ctx.toast("Copied message", 1.2);
  }
  ImGui::SameLine();
  if (ImGui::Button("Copy markup")) {
    ImGui::SetClipboardText(m.body.c_str());
    if (ctx.toast) ctx.toast("Copied markup", 1.2);
  }

  if (m.systemId != 0 && ctx.plotTo) {
    ImGui::SameLine();
    if (ImGui::Button("Plot route")) {
      ctx.plotTo(m.systemId, m.stationId);
    }
  }

  ImGui::Separator();

  // Body
  ImGui::BeginChild("##body_scroll", ImVec2(0, 0), false, ImGuiWindowFlags_HorizontalScrollbar);
  const float wrap = st.wrapBody ? ImGui::GetContentRegionAvail().x : 0.0f;
  const float tLocal = (float)std::max(0.0, ctx.timeRealSec - st.selectedOpenedSec);
  drawFxWrapped(m.body, wrap, tLocal);
  ImGui::EndChild();

  ImGui::EndChild();
  ImGui::End();
}

} // namespace stellar::game
