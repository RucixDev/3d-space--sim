#include "CommandPalette.h"

#include "stellar/ui/FuzzySearch.h"

#include <imgui.h>

#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

namespace stellar::game {

static bool anyNonSpace(const char* s) {
  if (!s) return false;
  for (const unsigned char* p = (const unsigned char*)s; *p; ++p) {
    if (std::isspace(*p) == 0) return true;
  }
  return false;
}

static void drawHighlightedText(std::string_view text, const std::vector<int>& positions, ImU32 colText, ImU32 colHi) {
  if (positions.empty()) {
    ImGui::TextUnformatted(text.data(), text.data() + text.size());
    return;
  }

  // Build a mask for contiguous range coalescing.
  std::vector<unsigned char> mark(text.size(), 0);
  for (int p : positions) {
    if (p >= 0 && (std::size_t)p < text.size()) mark[(std::size_t)p] = 1;
  }

  ImDrawList* dl = ImGui::GetWindowDrawList();
  ImVec2 cur = ImGui::GetCursorScreenPos();
  // Colors are passed in so callers can render disabled text while still highlighting.
  const float lineH = ImGui::GetTextLineHeight();

  const char* base = text.data();

  std::size_t i = 0;
  while (i < text.size()) {
    const bool hi = mark[i] != 0;
    std::size_t j = i + 1;
    while (j < text.size() && (mark[j] != 0) == hi) ++j;

    const char* segB = base + i;
    const char* segE = base + j;
    dl->AddText(cur, hi ? colHi : colText, segB, segE);

    const ImVec2 sz = ImGui::CalcTextSize(segB, segE);
    cur.x += sz.x;
    i = j;
  }

  // Advance cursor.
  ImGui::Dummy(ImVec2(cur.x - ImGui::GetCursorScreenPos().x, lineH));
}

static void drawHighlightedText(std::string_view text, const std::vector<int>& positions) {
  drawHighlightedText(text, positions,
                      ImGui::GetColorU32(ImGuiCol_Text),
                      ImGui::GetColorU32(ImGuiCol_PlotHistogram));
}

void openCommandPalette(CommandPaletteState& st, bool resetQuery) {
  st.open = true;
  st.justOpened = true;
  st.focusInput = true;
  st.selected = 0;
  if (resetQuery) st.query[0] = '\0';
}

bool drawCommandPalette(CommandPaletteState& st, const std::vector<PaletteItem>& items) {
  bool executed = false;

  if (st.open && st.justOpened) {
    ImGui::OpenPopup("Command Palette");
    st.justOpened = false;
  }

  // Centered popup.
  const ImVec2 vp = ImGui::GetMainViewport()->GetWorkCenter();
  ImGui::SetNextWindowPos(vp, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
  ImGui::SetNextWindowSize(ImVec2(640.0f, 420.0f), ImGuiCond_Appearing);

  ImGuiWindowFlags flags = ImGuiWindowFlags_NoDocking;
  flags |= ImGuiWindowFlags_NoTitleBar;
  flags |= ImGuiWindowFlags_NoResize;
  flags |= ImGuiWindowFlags_NoMove;

  if (ImGui::BeginPopup("Command Palette", flags)) {
    // Close on Escape.
    if (ImGui::IsKeyPressed(ImGuiKey_Escape)) {
      ImGui::CloseCurrentPopup();
      st.open = false;
      st.focusInput = false;
      ImGui::EndPopup();
      return false;
    }

    ImGui::TextUnformatted("Command Palette");
    ImGui::SameLine();
    ImGui::TextDisabled("(Enter: run, Esc: close, ↑↓: navigate)");

    // Filter input.
    // (SetKeyboardFocusHere() must be called *before* the widget we want to focus.)
    if (st.focusInput) {
      ImGui::SetKeyboardFocusHere();
      st.focusInput = false;
    }

    ImGui::PushItemWidth(-1);
    if (ImGui::InputTextWithHint("##palette_query", "Type a command or search...", st.query, sizeof(st.query),
                                 ImGuiInputTextFlags_AutoSelectAll)) {
      st.selected = 0;
    }
    ImGui::PopItemWidth();

    // Build filtered index list.
    struct Match {
      int idx;
      int score;
      std::vector<int> labelPos;
      std::vector<int> detailPos;
    };
    std::vector<Match> matches;
    matches.reserve(items.size());

    const bool hasQuery = anyNonSpace(st.query);
    for (int i = 0; i < (int)items.size(); ++i) {
      const PaletteItem& it = items[(std::size_t)i];
      int score = 0;
      std::vector<int> labelPos;
      std::vector<int> detailPos;
      if (hasQuery) {
        // Match against label + detail to make discovery easier.
        std::string hay;
        hay.reserve(it.label.size() + it.detail.size() + 2);
        hay += it.label;
        if (!it.detail.empty()) {
          hay.push_back(' ');
          hay += it.detail;
        }
        const auto r = stellar::ui::fuzzyMatch(st.query, hay);
        score = r.score;
        if (score < 0) continue;

        const int split = (int)it.label.size();
        for (int p : r.positions) {
          if (p < 0) continue;
          if (p < split) {
            labelPos.push_back(p);
          } else if (!it.detail.empty() && p > split) {
            detailPos.push_back(p - split - 1);
          }
        }
      } else {
        score = it.priority;
      }
      matches.push_back({i, score, std::move(labelPos), std::move(detailPos)});
    }

    std::stable_sort(matches.begin(), matches.end(), [&](const Match& a, const Match& b) {
      if (a.score != b.score) return a.score > b.score;
      // Secondary: prefer higher priority even during matching.
      if (items[(std::size_t)a.idx].priority != items[(std::size_t)b.idx].priority)
        return items[(std::size_t)a.idx].priority > items[(std::size_t)b.idx].priority;
      return items[(std::size_t)a.idx].label < items[(std::size_t)b.idx].label;
    });

    const int maxShow = 80;
    if ((int)matches.size() > maxShow) matches.resize((std::size_t)maxShow);

    if (st.selected < 0) st.selected = 0;
    if (st.selected >= (int)matches.size()) st.selected = std::max(0, (int)matches.size() - 1);

    // Keyboard navigation.
    if (ImGui::IsKeyPressed(ImGuiKey_DownArrow) && !matches.empty()) {
      st.selected = std::min(st.selected + 1, (int)matches.size() - 1);
    }
    if (ImGui::IsKeyPressed(ImGuiKey_UpArrow) && !matches.empty()) {
      st.selected = std::max(st.selected - 1, 0);
    }

    // Run selected.
    const bool enter = ImGui::IsKeyPressed(ImGuiKey_Enter) || ImGui::IsKeyPressed(ImGuiKey_KeypadEnter);
    if (enter && !matches.empty()) {
      const PaletteItem& it = items[(std::size_t)matches[(std::size_t)st.selected].idx];
      if (it.action) {
        it.action();
        executed = true;
      }
      ImGui::CloseCurrentPopup();
      st.open = false;
      ImGui::EndPopup();
      return executed;
    }

    ImGui::Separator();

    // Result list.
    const float listHeight = 320.0f;
    if (ImGui::BeginChild("##palette_list", ImVec2(0, listHeight), true,
                          ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
      if (ImGui::BeginTable("##palette_table", 3,
                            ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_NoSavedSettings)) {
        ImGui::TableSetupColumn("Command", ImGuiTableColumnFlags_WidthStretch, 0.65f);
        ImGui::TableSetupColumn("Details", ImGuiTableColumnFlags_WidthStretch, 0.25f);
        ImGui::TableSetupColumn("Key", ImGuiTableColumnFlags_WidthStretch, 0.10f);

        for (int row = 0; row < (int)matches.size(); ++row) {
          const Match& m = matches[(std::size_t)row];
          const PaletteItem& it = items[(std::size_t)m.idx];
          ImGui::TableNextRow();
          ImGui::TableSetColumnIndex(0);

          const bool sel = (row == st.selected);
          const std::string selId = "##pal_sel_" + std::to_string(row);
          const ImVec2 startPos = ImGui::GetCursorScreenPos();
          const bool clicked = ImGui::Selectable(selId.c_str(), sel,
                                                ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowItemOverlap);
          const ImVec2 afterPos = ImGui::GetCursorScreenPos();

          // Draw highlighted label on top of the selectable.
          ImGui::SetCursorScreenPos(startPos);
          if (hasQuery) {
            drawHighlightedText(it.label, m.labelPos);
          } else {
            ImGui::TextUnformatted(it.label.c_str());
          }
          ImGui::SetCursorScreenPos(afterPos);

          if (clicked) {
            st.selected = row;
            if (it.action) {
              it.action();
              executed = true;
            }
            ImGui::CloseCurrentPopup();
            st.open = false;
            break;
          }

          ImGui::TableSetColumnIndex(1);
          if (!it.detail.empty()) {
            if (hasQuery) {
              drawHighlightedText(it.detail, m.detailPos,
                                  ImGui::GetColorU32(ImGuiCol_TextDisabled),
                                  ImGui::GetColorU32(ImGuiCol_PlotHistogram));
            } else {
              ImGui::TextDisabled("%s", it.detail.c_str());
            }
          }

          ImGui::TableSetColumnIndex(2);
          if (!it.shortcut.empty()) {
            ImGui::TextDisabled("%s", it.shortcut.c_str());
          }
        }
        ImGui::EndTable();
      }
    }
    ImGui::EndChild();

    if (matches.empty()) {
      ImGui::TextDisabled("No matches.");
    }

    ImGui::EndPopup();
  } else {
    // Popup not open anymore.
    st.open = false;
  }

  return executed;
}

} // namespace stellar::game
