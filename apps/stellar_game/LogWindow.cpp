#include "LogWindow.h"

#include "stellar/ui/FuzzySearch.h"

#include <imgui.h>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace stellar::game {

static const char* levelTag(core::LogLevel lvl) {
  switch (lvl) {
    case core::LogLevel::Trace: return "TRACE";
    case core::LogLevel::Debug: return "DEBUG";
    case core::LogLevel::Info:  return "INFO";
    case core::LogLevel::Warn:  return "WARN";
    case core::LogLevel::Error: return "ERROR";
    default: return "?";
  }
}

static bool levelEnabled(core::LogLevel lvl, const LogWindowState& st) {
  switch (lvl) {
    case core::LogLevel::Trace: return st.showTrace;
    case core::LogLevel::Debug: return st.showDebug;
    case core::LogLevel::Info:  return st.showInfo;
    case core::LogLevel::Warn:  return st.showWarn;
    case core::LogLevel::Error: return st.showError;
    default: return true;
  }
}

static bool anyNonSpace(const char* s) {
  if (!s) return false;
  for (const unsigned char* p = (const unsigned char*)s; *p; ++p) {
    if (std::isspace(*p) == 0) return true;
  }
  return false;
}

static void drawHighlighted(std::string_view text, const std::vector<int>& positions, ImU32 colText) {
  if (positions.empty()) {
    ImGui::TextUnformatted(text.data(), text.data() + text.size());
    return;
  }

  std::vector<unsigned char> mark(text.size(), 0);
  for (int p : positions) {
    if (p >= 0 && (std::size_t)p < text.size()) mark[(std::size_t)p] = 1;
  }

  ImDrawList* dl = ImGui::GetWindowDrawList();
  ImVec2 cur = ImGui::GetCursorScreenPos();
  const ImU32 colHi = ImGui::GetColorU32(ImGuiCol_PlotHistogram);
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

  // Advance cursor by the rendered width.
  ImGui::Dummy(ImVec2(cur.x - ImGui::GetCursorScreenPos().x, lineH));
}

static std::string formatLine(const ui::LogEntry& e) {
  std::ostringstream oss;
  oss << '[' << e.timestamp << "][" << levelTag(e.level) << "] " << e.message;
  return oss.str();
}

static bool saveLines(const char* path, const std::vector<std::string>& lines, std::string* outErr) {
  if (!path || !path[0]) {
    if (outErr) *outErr = "Empty path.";
    return false;
  }

  namespace fs = std::filesystem;

  std::error_code ec;
  fs::path p(path);
  if (p.has_parent_path()) {
    fs::create_directories(p.parent_path(), ec);
    // Ignore errors if the directory already exists; other failures will show when opening.
  }

  std::ofstream f(path, std::ios::binary);
  if (!f) {
    if (outErr) *outErr = "Failed to open file for writing.";
    return false;
  }

  for (const auto& l : lines) {
    f << l << "\n";
  }
  return true;
}

void drawLogWindow(LogWindowState& st,
                   ui::LogBuffer& buffer,
                   const ToastFn& toast) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(900.0f, 520.0f), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Log", &st.open)) {
    ImGui::End();
    return;
  }

  std::vector<ui::LogEntry> entries;
  buffer.copyAll(entries);

  ImGui::TextDisabled("%zu entries", entries.size());
  ImGui::SameLine();
  if (ImGui::Button("Clear")) {
    buffer.clear();
    st.selectedSeq = 0;
    if (toast) toast("Log cleared.", 1.2);
  }
  ImGui::SameLine();
  ImGui::Checkbox("Auto-scroll", &st.autoScroll);

  ImGui::Separator();

  // ---- Filters ----
  if (st.focusFilter) {
    ImGui::SetKeyboardFocusHere();
    st.focusFilter = false;
  }
  ImGui::InputTextWithHint("##log_filter", "Filter (fuzzy)...", st.filter, sizeof(st.filter));
  const bool hasFilter = anyNonSpace(st.filter);
  if (hasFilter) {
    ImGui::SameLine();
    ImGui::Checkbox("Sort by relevance", &st.sortByRelevance);
  }

  ImGui::TextDisabled("Levels:");
  ImGui::SameLine();
  ImGui::Checkbox("Trace", &st.showTrace);
  ImGui::SameLine();
  ImGui::Checkbox("Debug", &st.showDebug);
  ImGui::SameLine();
  ImGui::Checkbox("Info", &st.showInfo);
  ImGui::SameLine();
  ImGui::Checkbox("Warn", &st.showWarn);
  ImGui::SameLine();
  ImGui::Checkbox("Error", &st.showError);

  // ---- Selection actions ----
  const bool hasSelection = (st.selectedSeq != 0);
  if (!hasSelection) ImGui::BeginDisabled();
  if (ImGui::Button("Copy selected") && hasSelection) {
    for (const auto& e : entries) {
      if (e.seq == st.selectedSeq) {
        const std::string line = formatLine(e);
        ImGui::SetClipboardText(line.c_str());
        if (toast) toast("Copied selected log line.", 1.2);
        break;
      }
    }
  }
  if (!hasSelection) ImGui::EndDisabled();

  ImGui::SameLine();
  if (ImGui::Button("Copy filtered")) {
    std::ostringstream oss;
    int copied = 0;
    for (const auto& e : entries) {
      if (!levelEnabled(e.level, st)) continue;
      if (hasFilter) {
        if (ui::fuzzyMatchScore(st.filter, e.message) < 0) continue;
      }
      oss << formatLine(e) << "\n";
      ++copied;
    }
    const std::string txt = oss.str();
    ImGui::SetClipboardText(txt.c_str());
    if (toast) toast("Copied " + std::to_string(copied) + " log lines.", 1.4);
  }

  ImGui::Separator();

  // ---- Export ----
  ImGui::TextDisabled("Export");
  ImGui::SameLine();
  ImGui::SetNextItemWidth(-200.0f);
  ImGui::InputText("##log_export_path", st.exportPath, sizeof(st.exportPath));
  ImGui::SameLine();
  if (ImGui::Button("Save filtered")) {
    std::vector<std::string> lines;
    lines.reserve(entries.size());
    for (const auto& e : entries) {
      if (!levelEnabled(e.level, st)) continue;
      if (hasFilter) {
        if (ui::fuzzyMatchScore(st.filter, e.message) < 0) continue;
      }
      lines.push_back(formatLine(e));
    }
    std::string err;
    if (saveLines(st.exportPath, lines, &err)) {
      if (toast) toast(std::string("Saved ") + std::to_string(lines.size()) + " lines.", 1.8);
    } else {
      if (toast) toast(std::string("Failed to save log: ") + err, 2.4);
    }
  }

  ImGui::Separator();

  // ---- Table ----
  struct Row {
    const ui::LogEntry* e{nullptr};
    int score{0};
    std::vector<int> pos;
  };
  std::vector<Row> rows;
  rows.reserve(entries.size());

  for (const auto& e : entries) {
    if (!levelEnabled(e.level, st)) continue;
    if (!hasFilter) {
      rows.push_back(Row{&e, 0, {}});
      continue;
    }
    const auto r = ui::fuzzyMatch(st.filter, e.message);
    if (r.score < 0) continue;
    rows.push_back(Row{&e, r.score, r.positions});
  }

  if (hasFilter && st.sortByRelevance) {
    std::stable_sort(rows.begin(), rows.end(), [&](const Row& a, const Row& b) {
      if (a.score != b.score) return a.score > b.score;
      // Prefer newer entries when scores tie.
      return a.e->seq > b.e->seq;
    });
  }

  if (ImGui::BeginTable("##log_table", 3,
                        ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY
                        | ImGuiTableFlags_SizingStretchProp)) {
    ImGui::TableSetupColumn("Time", ImGuiTableColumnFlags_WidthFixed, 96.0f);
    ImGui::TableSetupColumn("Lvl", ImGuiTableColumnFlags_WidthFixed, 64.0f);
    ImGui::TableSetupColumn("Message", ImGuiTableColumnFlags_WidthStretch);
    ImGui::TableHeadersRow();

    for (int i = 0; i < (int)rows.size(); ++i) {
      const Row& row = rows[(std::size_t)i];
      const ui::LogEntry& e = *row.e;

      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      ImGui::TextDisabled("%s", e.timestamp.c_str());

      ImGui::TableSetColumnIndex(1);
      ImGui::TextUnformatted(levelTag(e.level));

      ImGui::TableSetColumnIndex(2);
      const bool sel = (st.selectedSeq == e.seq);
      const std::string id = "##log_sel_" + std::to_string((unsigned long long)e.seq);
      const ImVec2 startPos = ImGui::GetCursorScreenPos();
      const bool clicked = ImGui::Selectable(id.c_str(), sel,
                                            ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap);
      const ImVec2 afterPos = ImGui::GetCursorScreenPos();

      // Overlay text (so we can draw highlights while keeping full-row select).
      ImGui::SetCursorScreenPos(startPos);
      if (hasFilter) {
        drawHighlighted(e.message, row.pos, ImGui::GetColorU32(ImGuiCol_Text));
      } else {
        ImGui::TextUnformatted(e.message.c_str());
      }
      ImGui::SetCursorScreenPos(afterPos);

      if (clicked) st.selectedSeq = e.seq;
    }

    if (st.autoScroll && !rows.empty()) {
      ImGui::SetScrollHereY(1.0f);
    }

    ImGui::EndTable();
  }

  ImGui::End();
}

} // namespace stellar::game
