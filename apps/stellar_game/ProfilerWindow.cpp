#include "ProfilerWindow.h"

#include "stellar/core/ChromeTrace.h"
#include "stellar/core/Hash.h"

#include <imgui.h>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>

namespace stellar::game {

static double nsToMs(std::uint64_t ns) {
  return (double)ns / 1e6;
}

static bool icontains(std::string_view hay, std::string_view needle) {
  if (needle.empty()) return true;
  if (hay.empty()) return false;

  // naive ascii lowercase search; good enough for debug tooling.
  auto lower = [](char c) -> char {
    return (char)std::tolower((unsigned char)c);
  };

  for (std::size_t i = 0; i < hay.size(); ++i) {
    std::size_t j = 0;
    while (i + j < hay.size() && j < needle.size() &&
           lower(hay[i + j]) == lower(needle[j])) {
      ++j;
    }
    if (j == needle.size()) return true;
  }

  return false;
}

static ImU32 colorForName(std::string_view name) {
  const stellar::core::u64 h = stellar::core::fnv1a64(name);

  const float hue = (float)(h % 360) / 360.0f;
  const float sat = 0.55f + 0.25f * (float)((h >> 8) & 0xFF) / 255.0f;
  const float val = 0.55f + 0.25f * (float)((h >> 16) & 0xFF) / 255.0f;

  // Slight transparency so overlapping UI remains readable.
  return ImColor::HSV(hue, sat, val, 0.90f);
}

struct AggRow {
  std::string_view name;
  std::uint64_t totalNs{0};
  std::uint64_t maxNs{0};
  int count{0};
};

static void drawFramePlot(const std::deque<core::ProfilerFrame>& frames) {
  if (frames.empty()) {
    ImGui::TextUnformatted("No frames captured yet.");
    return;
  }

  static std::vector<float> values;
  values.clear();
  values.reserve(frames.size());

  float maxMs = 0.0f;
  float sumMs = 0.0f;
  for (const auto& f : frames) {
    const float ms = (float)nsToMs(f.durationNs());
    values.push_back(ms);
    sumMs += ms;
    maxMs = std::max(maxMs, ms);
  }

  const float avgMs = sumMs / (float)values.size();

  ImGui::Text("History: %zu frames | avg %.3f ms | max %.3f ms", values.size(), avgMs, maxMs);

  // PlotLines expects a pointer to a contiguous float buffer.
  ImGui::PlotLines("##frame_times",
                   values.data(),
                   (int)values.size(),
                   0,
                   nullptr,
                   0.0f,
                   maxMs * 1.10f,
                   ImVec2(0, 80));
}

static void drawFlameGraph(const core::ProfilerFrame& frame, const char* id) {
  const std::uint64_t frameDurNs = frame.durationNs();
  if (frameDurNs == 0 || frame.events.empty()) {
    ImGui::TextUnformatted("No events in selected frame.");
    return;
  }

  // Copy pointers & sort by start time to draw left-to-right.
  std::vector<const core::ProfilerEvent*> events;
  events.reserve(frame.events.size());
  std::uint32_t maxDepth = 0;
  for (const auto& e : frame.events) {
    events.push_back(&e);
    maxDepth = std::max(maxDepth, e.depth);
  }
  std::sort(events.begin(), events.end(), [](const auto* a, const auto* b) {
    if (a->startNs != b->startNs) return a->startNs < b->startNs;
    return a->depth < b->depth;
  });

  const float rowH = ImGui::GetTextLineHeight() + 4.0f;
  const float height = std::max(80.0f, rowH * (float)(maxDepth + 1) + 8.0f);

  ImGui::BeginChild(id, ImVec2(0, height), true, ImGuiWindowFlags_NoScrollWithMouse);

  const ImVec2 p0 = ImGui::GetCursorScreenPos();
  const ImVec2 avail = ImGui::GetContentRegionAvail();
  const ImVec2 size = ImVec2(avail.x, height - 8.0f);

  ImGui::InvisibleButton("##flame_canvas", size);
  const ImVec2 c0 = ImGui::GetItemRectMin();
  const ImVec2 c1 = ImGui::GetItemRectMax();

  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(c0, c1, ImGui::GetColorU32(ImGuiCol_FrameBg));
  dl->AddRect(c0, c1, ImGui::GetColorU32(ImGuiCol_Border));

  dl->PushClipRect(c0, c1, true);

  const float w = c1.x - c0.x;

  const core::ProfilerEvent* hovered = nullptr;

  for (const core::ProfilerEvent* e : events) {
    const std::uint64_t s = (e->startNs >= frame.startNs) ? (e->startNs - frame.startNs) : 0;
    const std::uint64_t t = (e->endNs >= frame.startNs) ? (e->endNs - frame.startNs) : 0;

    const float x0 = c0.x + (float)((double)s / (double)frameDurNs) * w;
    float x1 = c0.x + (float)((double)t / (double)frameDurNs) * w;
    if (x1 < x0 + 1.0f) x1 = x0 + 1.0f;

    const float y0 = c0.y + rowH * (float)e->depth;
    const float y1 = y0 + rowH - 1.0f;

    const ImVec2 r0(x0, y0);
    const ImVec2 r1(x1, y1);

    const char* name = e->name ? e->name : "<null>";
    const ImU32 col = colorForName(name);

    dl->AddRectFilled(r0, r1, col);
    dl->AddRect(r0, r1, ImGui::GetColorU32(ImGuiCol_Border));

    if ((x1 - x0) > 40.0f) {
      dl->AddText(ImVec2(x0 + 3.0f, y0 + 2.0f), ImGui::GetColorU32(ImGuiCol_Text), name);
    }

    if (ImGui::IsItemHovered() && ImGui::IsMouseHoveringRect(r0, r1, true)) {
      hovered = e;
    }
  }

  dl->PopClipRect();

  if (hovered) {
    const char* name = hovered->name ? hovered->name : "<null>";
    const double ms = nsToMs(hovered->durationNs());
    const double pct = 100.0 * (double)hovered->durationNs() / (double)frameDurNs;

    ImGui::BeginTooltip();
    ImGui::TextUnformatted(name);
    ImGui::Separator();
    ImGui::Text("%.3f ms (%.2f%%)", ms, pct);
    ImGui::Text("Depth: %u", hovered->depth);
    ImGui::TextUnformatted("Right-click to copy scope name");
    ImGui::EndTooltip();

    if (ImGui::IsMouseReleased(ImGuiMouseButton_Right)) {
      ImGui::SetClipboardText(name);
    }
  }

  ImGui::EndChild();
}

static void drawAggregates(const core::ProfilerFrame& frame,
                           std::string_view filter,
                           const ToastFn& toast) {
  if (frame.events.empty()) {
    ImGui::TextUnformatted("No events in selected frame.");
    return;
  }

  std::vector<AggRow> rows;
  rows.reserve(frame.events.size());

  auto findRow = [&](std::string_view name) -> AggRow* {
    for (auto& r : rows) {
      if (r.name == name) return &r;
    }
    return nullptr;
  };

  for (const auto& e : frame.events) {
    const char* n = e.name ? e.name : "<null>";
    const std::string_view name{n};
    const std::uint64_t d = e.durationNs();

    AggRow* r = findRow(name);
    if (!r) {
      rows.push_back(AggRow{name});
      r = &rows.back();
    }

    r->totalNs += d;
    r->maxNs = std::max(r->maxNs, d);
    r->count += 1;
  }

  // Filter & sort.
  std::vector<AggRow> filtered;
  filtered.reserve(rows.size());
  for (const auto& r : rows) {
    if (icontains(r.name, filter)) {
      filtered.push_back(r);
    }
  }

  std::sort(filtered.begin(), filtered.end(), [](const AggRow& a, const AggRow& b) {
    return a.totalNs > b.totalNs;
  });

  if (ImGui::BeginTable("##agg_table", 5, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV)) {
    ImGui::TableSetupColumn("Scope", ImGuiTableColumnFlags_WidthStretch);
    ImGui::TableSetupColumn("Total (ms)", ImGuiTableColumnFlags_WidthFixed);
    ImGui::TableSetupColumn("Avg (ms)", ImGuiTableColumnFlags_WidthFixed);
    ImGui::TableSetupColumn("Max (ms)", ImGuiTableColumnFlags_WidthFixed);
    ImGui::TableSetupColumn("Count", ImGuiTableColumnFlags_WidthFixed);
    ImGui::TableHeadersRow();

    int rowIdx = 0;
    for (const auto& r : filtered) {
      ImGui::TableNextRow();
      ImGui::TableNextColumn();

      // Color chip.
      const ImU32 col = colorForName(r.name);
      ImGui::ColorButton((std::string("##c") + std::to_string(rowIdx)).c_str(),
                         ImColor(col),
                         ImGuiColorEditFlags_NoTooltip,
                         ImVec2(10, 10));
      ImGui::SameLine();

      if (ImGui::Selectable(std::string(r.name).c_str(), false, ImGuiSelectableFlags_SpanAllColumns)) {
        ImGui::SetClipboardText(std::string(r.name).c_str());
        if (toast) {
          toast(std::string("Copied scope: ") + std::string(r.name), 1.2);
        }
      }

      ImGui::TableNextColumn();
      ImGui::Text("%.3f", nsToMs(r.totalNs));

      ImGui::TableNextColumn();
      ImGui::Text("%.3f", nsToMs(r.totalNs) / (double)std::max(1, r.count));

      ImGui::TableNextColumn();
      ImGui::Text("%.3f", nsToMs(r.maxNs));

      ImGui::TableNextColumn();
      ImGui::Text("%d", r.count);

      ++rowIdx;
    }

    ImGui::EndTable();
  }
}

void drawProfilerWindow(ProfilerWindowState& st,
                        core::Profiler& profiler,
                        const ToastFn& toast) {
  if (!st.open) return;

  if (!ImGui::Begin("Profiler", &st.open)) {
    ImGui::End();
    return;
  }

  // Controls
  {
    bool enabled = profiler.enabled();
    if (ImGui::Checkbox("Enabled", &enabled)) {
      profiler.setEnabled(enabled);
      if (toast) {
        toast(enabled ? "Profiler enabled" : "Profiler disabled", 1.2);
      }
    }

    ImGui::SameLine();
    if (ImGui::Button("Clear")) {
      profiler.clear();
      if (toast) toast("Profiler history cleared", 1.2);
    }

    ImGui::SameLine();
    ImGui::Checkbox("Pause", &st.pause);

    int maxFrames = (int)profiler.maxFrames();
    if (ImGui::SliderInt("History frames", &maxFrames, 30, 2000)) {
      profiler.setMaxFrames((std::size_t)maxFrames);
    }

    if (st.pause) {
      const int maxOffset = (int)std::max<std::size_t>(0, profiler.frames().size() ? profiler.frames().size() - 1 : 0);
      st.selectedFrameOffset = std::clamp(st.selectedFrameOffset, 0, maxOffset);
      ImGui::SliderInt("Frame offset", &st.selectedFrameOffset, 0, maxOffset);
    } else {
      st.selectedFrameOffset = 0;
    }
  }

  ImGui::Separator();

  const auto& frames = profiler.frames();
  if (st.showPlot) {
    drawFramePlot(frames);
    ImGui::Separator();
  }

  const core::ProfilerFrame* sel = nullptr;
  if (!frames.empty()) {
    const std::size_t offset = (std::size_t)std::clamp(st.selectedFrameOffset, 0, (int)frames.size() - 1);
    sel = &frames[frames.size() - 1 - offset];
  }

  if (!sel) {
    ImGui::TextUnformatted("No frames captured (enable the profiler and wait a frame).");
    ImGui::End();
    return;
  }

  ImGui::Text("Selected frame: %d from newest | %.3f ms | %zu events",
              st.selectedFrameOffset,
              nsToMs(sel->durationNs()),
              sel->events.size());

  if (ImGui::CollapsingHeader("Export Trace", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::TextWrapped("Writes a Chrome JSON trace (chrome://tracing / Perfetto legacy JSON). ");

    ImGui::InputTextWithHint("Path",
                            "profiler_trace.json",
                            st.exportPath,
                            sizeof(st.exportPath));

    ImGui::Checkbox("Export all history frames", &st.exportAllFrames);
    ImGui::SameLine();
    ImGui::Checkbox("Include frame spans", &st.exportIncludeFrameEvents);
    ImGui::SameLine();
    ImGui::Checkbox("Pretty JSON", &st.exportPretty);

    core::ChromeTraceWriteOptions opt;
    opt.includeFrameEvents = st.exportIncludeFrameEvents;
    opt.pretty = st.exportPretty;

    if (ImGui::Button("Export##trace")) {
      std::string err;
      bool ok = false;

      if (st.exportAllFrames) {
        ok = core::writeProfilerChromeTraceJson(st.exportPath, profiler.frames(), opt, &err);
      } else {
        ok = core::writeProfilerChromeTraceJson(st.exportPath, *sel, opt, &err);
      }

      if (toast) {
        if (ok) {
          toast(std::string("Wrote trace: ") + st.exportPath, 2.2);
        } else {
          toast(std::string("Trace export failed: ") + (err.empty() ? "unknown error" : err), 2.8);
        }
      }
    }

    ImGui::SameLine();
    if (ImGui::Button("Copy path##trace")) {
      ImGui::SetClipboardText(st.exportPath);
      if (toast) toast("Copied trace path to clipboard.", 1.2);
    }

    ImGui::TextDisabled("Open the JSON in chrome://tracing (legacy) or Perfetto UI.");
  }

  if (st.showFlameGraph) {
    if (ImGui::CollapsingHeader("Flame Graph", ImGuiTreeNodeFlags_DefaultOpen)) {
      drawFlameGraph(*sel, "##flame");
    }
  }

  if (st.showAggregates) {
    if (ImGui::CollapsingHeader("Aggregates", ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::InputTextWithHint("Filter", "substring...", st.filter, sizeof(st.filter));
      drawAggregates(*sel, std::string_view(st.filter), toast);
    }
  }

  ImGui::End();
}

} // namespace stellar::game
