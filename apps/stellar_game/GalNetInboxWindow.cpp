#include "GalNetInboxWindow.h"

#include "stellar/sim/Universe.h"
#include "stellar/ui/TextFx.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace stellar::game {

namespace {

bool containsSystemId(const std::vector<stellar::sim::SystemId>& ids, stellar::sim::SystemId id) {
  for (const auto v : ids) {
    if (v == id) {
      return true;
    }
  }
  return false;
}

void toggleSystemId(std::vector<stellar::sim::SystemId>& ids, stellar::sim::SystemId id) {
  for (size_t i = 0; i < ids.size(); ++i) {
    if (ids[i] == id) {
      ids.erase(ids.begin() + static_cast<std::ptrdiff_t>(i));
      return;
    }
  }
  ids.push_back(id);
}

bool containsInsensitive(std::string_view haystack, std::string_view needle) {
  if (needle.empty()) {
    return true;
  }
  if (haystack.empty()) {
    return false;
  }

  auto lower = [](unsigned char c) { return static_cast<char>(std::tolower(c)); };

  for (size_t i = 0; i < haystack.size(); ++i) {
    size_t j = 0;
    while (i + j < haystack.size() && j < needle.size() && lower(haystack[i + j]) == lower(needle[j])) {
      ++j;
    }
    if (j == needle.size()) {
      return true;
    }
  }
  return false;
}

std::string ageShort(double nowDays, double msgTimeDays) {
  const double dtDays = std::max(0.0, nowDays - msgTimeDays);
  const double dtHours = dtDays * 24.0;
  if (dtHours < 1.0) {
    const int mins = static_cast<int>(std::round(dtHours * 60.0));
    return std::to_string(std::max(0, mins)) + "m";
  }
  if (dtHours < 48.0) {
    const int hrs = static_cast<int>(std::round(dtHours));
    return std::to_string(std::max(0, hrs)) + "h";
  }
  const int days = static_cast<int>(std::round(dtDays));
  return std::to_string(std::max(0, days)) + "d";
}

std::string systemNameForId(stellar::sim::Universe& universe, stellar::sim::SystemId id) {
  if (id == 0) {
    return {};
  }
  const auto& sys = universe.getSystem(id);
  return sys.stub.name;
}

bool messageHasActiveEventHint(const stellar::sim::CommsMessage& msg) {
  // Full GalNet bulletins include this sentence when posted without an active event.
  // Digest messages may contain both active and status lines; treat them as active
  // if at least one "[...]%" entry is present.
  if (msg.body.find("No active system event") != std::string::npos) {
    return false;
  }
  if (msg.subject.find("GalNet Digest") != std::string::npos) {
    return msg.body.find("%]") != std::string::npos;
  }
  return true;
}

struct Row {
  stellar::sim::CommsMessage* msg{};
  std::string systemName{};
};

} // namespace

void drawGalNetInboxWindow(GalNetInboxWindowState& st, const GalNetInboxWindowContext& ctx) {
  if (!st.open) {
    return;
  }

  if (!ImGui::Begin("GalNet Inbox", &st.open)) {
    ImGui::End();
    return;
  }

  ImGui::TextDisabled("Filtered view of your Comms inbox (GalNet only)");

  ImGui::Checkbox("Newest first", &st.newestFirst);
  ImGui::SameLine();
  ImGui::Checkbox("Watched only", &st.watchedOnly);
  ImGui::SameLine();
  ImGui::Checkbox("Active only", &st.activeOnly);

  ImGui::SetNextItemWidth(280.0f);
  ImGui::InputTextWithHint("##galnet_filter", "Filter: system / subject / body", st.filter, sizeof(st.filter));
  ImGui::SameLine();
  ImGui::Checkbox("Wrap body", &st.wrapBody);
  ImGui::SameLine();
  ImGui::Checkbox("Mark read on select", &st.autoMarkReadOnSelect);

  const std::string_view filterSv{st.filter};

  // Collect matching messages (pointers into CommsLog).
  std::vector<Row> rows;
  rows.reserve(ctx.comms.items().size());

  auto& items = ctx.comms.itemsMutable();
  for (auto& msg : items) {
    if (msg.from != "GalNet") {
      continue;
    }
    if (msg.systemId == 0) {
      continue;
    }

    if (st.watchedOnly && ctx.watchSystems) {
      if (!containsSystemId(*ctx.watchSystems, msg.systemId)) {
        continue;
      }
    }

    if (st.activeOnly) {
      if (!messageHasActiveEventHint(msg)) {
        continue;
      }
    }

    std::string systemName = systemNameForId(ctx.universe, msg.systemId);

    const bool matches = containsInsensitive(msg.subject, filterSv) || containsInsensitive(msg.body, filterSv) ||
                         containsInsensitive(systemName, filterSv);
    if (!matches) {
      continue;
    }

    rows.push_back(Row{&msg, std::move(systemName)});
  }

  std::stable_sort(rows.begin(), rows.end(), [&](const Row& a, const Row& b) {
    if (st.newestFirst) {
      return a.msg->timeDays > b.msg->timeDays;
    }
    return a.msg->timeDays < b.msg->timeDays;
  });

  // Two-pane layout.
  ImGui::BeginChild("##galnet_left", ImVec2(st.leftPaneWidth, 0), true);
  {
    const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_SizingFixedFit |
                                  ImGuiTableFlags_ScrollY;
    if (ImGui::BeginTable("##galnet_table", 4, flags)) {
      ImGui::TableSetupColumn("Age", ImGuiTableColumnFlags_WidthFixed, 46.0f);
      ImGui::TableSetupColumn("Sys", ImGuiTableColumnFlags_WidthFixed, 120.0f);
      ImGui::TableSetupColumn("Ch", ImGuiTableColumnFlags_WidthFixed, 72.0f);
      ImGui::TableSetupColumn("Subject", ImGuiTableColumnFlags_WidthStretch);
      ImGui::TableHeadersRow();

      for (auto& row : rows) {
        auto* msg = row.msg;
        ImGui::TableNextRow();

        ImGui::TableSetColumnIndex(0);
        const std::string age = ageShort(ctx.timeDays, msg->timeDays);
        ImGui::TextUnformatted(age.c_str());

        ImGui::TableSetColumnIndex(1);
        ImGui::TextUnformatted(row.systemName.c_str());

        ImGui::TableSetColumnIndex(2);
        ImGui::TextUnformatted(stellar::sim::commsChannelName(msg->channel));

        ImGui::TableSetColumnIndex(3);
        std::string label = stellar::ui::textfx::stripMarkup(msg->subject);
        if (msg->unread) {
          label = std::string("* ") + label;
        }
        if (msg->pinned) {
          label = std::string("[PIN] ") + label;
        }

        const bool selected = (st.selectedMsgId == msg->id);
        if (ImGui::Selectable(label.c_str(), selected, ImGuiSelectableFlags_SpanAllColumns)) {
          st.selectedMsgId = msg->id;
          if (st.autoMarkReadOnSelect) {
            ctx.comms.markRead(msg->id, /*read=*/true);
          }
        }
      }

      ImGui::EndTable();
    }
  }
  ImGui::EndChild();

  ImGui::SameLine();

  ImGui::BeginChild("##galnet_right", ImVec2(0, 0), true);
  {
    stellar::sim::CommsMessage* selectedMsg = ctx.comms.findMutable(st.selectedMsgId);
    if (!selectedMsg) {
      ImGui::TextDisabled("Select a bulletin from the list.");
      ImGui::EndChild();
      ImGui::End();
      return;
    }

    const std::string systemName = systemNameForId(ctx.universe, selectedMsg->systemId);

    ImGui::TextUnformatted(stellar::ui::textfx::stripMarkup(selectedMsg->subject).c_str());
    ImGui::Separator();

    // Quick actions.
    if (ImGui::Button(selectedMsg->unread ? "Mark read" : "Mark unread")) {
      const bool setRead = selectedMsg->unread;
      ctx.comms.markRead(selectedMsg->id, setRead);
    }
    ImGui::SameLine();
    if (ImGui::Button(selectedMsg->pinned ? "Unpin" : "Pin")) {
      ctx.comms.markPinned(selectedMsg->id, !selectedMsg->pinned);
    }

    if (ctx.watchSystems) {
      const bool watched = containsSystemId(*ctx.watchSystems, selectedMsg->systemId);
      ImGui::SameLine();
      if (ImGui::Button(watched ? "Unwatch" : "Watch")) {
        toggleSystemId(*ctx.watchSystems, selectedMsg->systemId);
        if (ctx.toast) {
          ctx.toast(watched ? "Removed from watchlist." : "Added to watchlist.", 1.6);
        }
      }
    }

    if (ctx.targetSystem) {
      ImGui::SameLine();
      if (ImGui::Button("Target system")) {
        ctx.targetSystem(selectedMsg->systemId);
      }
    }

    if (ctx.plotRouteToSystem) {
      ImGui::SameLine();
      if (ImGui::Button("Plot route")) {
        (void)ctx.plotRouteToSystem(selectedMsg->systemId);
      }
    }

    if (ctx.postGalNetBulletin) {
      ImGui::SameLine();
      if (ImGui::Button("Request update")) {
        ctx.postGalNetBulletin(selectedMsg->systemId, /*showOverlay=*/false, /*showToast=*/true);
      }
    }

    ImGui::Separator();

    ImGui::TextDisabled("From:");
    ImGui::SameLine();
    ImGui::TextUnformatted(selectedMsg->from.c_str());
    ImGui::SameLine();
    ImGui::TextDisabled("|");
    ImGui::SameLine();
    ImGui::TextDisabled("System:");
    ImGui::SameLine();
    ImGui::TextUnformatted(systemName.c_str());
    ImGui::SameLine();
    ImGui::TextDisabled("|");
    ImGui::SameLine();
    ImGui::TextDisabled("Age:");
    ImGui::SameLine();
    const std::string age = ageShort(ctx.timeDays, selectedMsg->timeDays);
    ImGui::TextUnformatted(age.c_str());

    ImGui::Separator();

    const std::string body = stellar::ui::textfx::stripMarkup(selectedMsg->body);
    if (st.wrapBody) {
      ImGui::TextWrapped("%s", body.c_str());
    } else {
      ImGui::TextUnformatted(body.c_str());
    }
  }
  ImGui::EndChild();

  ImGui::End();
}

} // namespace stellar::game
