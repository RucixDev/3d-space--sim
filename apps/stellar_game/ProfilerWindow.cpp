#include "ProfilerWindow.h"

#include "stellar/core/ChromeTrace.h"
#include "stellar/core/Hash.h"
#include "stellar/render/FrameGraph.h"

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

static void maybeCaptureGpuFrameGraphSample(ProfilerWindowState& st,
                                            const core::Profiler& profiler,
                                            const render::FrameGraph* fg) {
  if (!fg) return;
  const core::ProfilerFrame* newest = profiler.newest();
  if (!newest) return;

  // Capture exactly one GPU sample per completed CPU profiler frame.
  if (newest->startNs == st.gpuLastCapturedCpuFrameStartNs) return;
  st.gpuLastCapturedCpuFrameStartNs = newest->startNs;

  const auto& passes = fg->passes();
  const std::size_t passCount = passes.size();

  // Build stable counter keys: [0]=total GPU, then per-pass with an index prefix.
  std::vector<std::string> keys;
  keys.reserve(passCount + 1);
  keys.emplace_back("GPU Frame (ms)");
  for (std::size_t i = 0; i < passCount; ++i) {
    const std::string& name = passes[i].name;
    // Prefix with pass id to guarantee uniqueness.
    std::string k = "p";
    if (i < 10) k += "0";
    if (i < 100) {
      k += std::to_string(i);
    } else {
      k += std::to_string(i);
    }
    k += ": ";
    k += name.empty() ? std::string("<unnamed>") : name;
    keys.push_back(std::move(k));
  }

  if (st.gpuCounterKeys != keys) {
    st.gpuCounterKeys = std::move(keys);
    st.gpuHistoryRows.clear();
  }

  std::vector<double> row;
  row.assign(st.gpuCounterKeys.size(), 0.0);

  const bool canRead = fg->profilingSupported() && fg->profilingEnabled();
  if (canRead) {
    row[0] = fg->lastFrameTimeMs();
    const auto& times = fg->lastPassTimesMs();
    for (std::size_t i = 0; i < passCount && i < times.size(); ++i) {
      row[i + 1] = times[i];
    }
  }

  st.gpuHistoryRows.push_back(std::move(row));

  // Keep GPU history aligned to the CPU profiler ring buffer.
  const std::size_t target = profiler.frames().size();
  while (st.gpuHistoryRows.size() > target && !st.gpuHistoryRows.empty()) {
    st.gpuHistoryRows.pop_front();
  }
}

static void drawGpuFrameGraphPlot(const ProfilerWindowState& st) {
  if (st.gpuHistoryRows.empty()) {
    ImGui::TextUnformatted("No GPU timing samples captured yet.");
    return;
  }

  static std::vector<float> values;
  values.clear();
  values.reserve(st.gpuHistoryRows.size());

  float maxMs = 0.0f;
  float sumMs = 0.0f;
  for (const auto& r : st.gpuHistoryRows) {
    const float ms = (!r.empty()) ? (float)r[0] : 0.0f;
    values.push_back(ms);
    sumMs += ms;
    maxMs = std::max(maxMs, ms);
  }

  const float avgMs = sumMs / (float)std::max<std::size_t>(1, values.size());
  ImGui::Text("GPU History: %zu frames | avg %.3f ms | max %.3f ms",
              values.size(), avgMs, maxMs);

  ImGui::PlotLines("##gpu_frame_times",
                   values.data(),
                   (int)values.size(),
                   0,
                   nullptr,
                   0.0f,
                   maxMs * 1.10f,
                   ImVec2(0, 80));
}

struct GpuPassRow {
  std::string_view name;
  std::size_t passId{0};
  double lastMs{0.0};
  double avgMs{0.0};
  double maxMs{0.0};
};

static void drawGpuFrameGraphTable(const ProfilerWindowState& st,
                                   const render::FrameGraph& fg,
                                   std::string_view filter) {
  const auto& passes = fg.passes();
  const std::size_t passCount = passes.size();

  std::vector<GpuPassRow> rows;
  rows.reserve(passCount);

  const bool canRead = fg.profilingSupported() && fg.profilingEnabled();
  const auto& last = fg.lastPassTimesMs();

  const std::size_t histN = st.gpuHistoryRows.size();

  for (std::size_t i = 0; i < passCount; ++i) {
    const std::string_view name{passes[i].name};
    if (!icontains(name, filter)) continue;

    GpuPassRow r;
    r.name = name.empty() ? std::string_view{"<unnamed>"} : name;
    r.passId = i;
    r.lastMs = (canRead && i < last.size()) ? last[i] : 0.0;

    double sum = 0.0;
    double mx = 0.0;
    for (const auto& h : st.gpuHistoryRows) {
      const std::size_t idx = i + 1; // [0] is total
      const double v = (idx < h.size()) ? h[idx] : 0.0;
      sum += v;
      mx = std::max(mx, v);
    }
    r.avgMs = (histN > 0) ? (sum / (double)histN) : 0.0;
    r.maxMs = mx;

    rows.push_back(r);
  }

  if (st.gpuSortByLast) {
    std::sort(rows.begin(), rows.end(), [](const GpuPassRow& a, const GpuPassRow& b) {
      return a.lastMs > b.lastMs;
    });
  } else {
    std::sort(rows.begin(), rows.end(), [](const GpuPassRow& a, const GpuPassRow& b) {
      return a.avgMs > b.avgMs;
    });
  }

  const double frameMs = canRead ? fg.lastFrameTimeMs() : 0.0;

  if (ImGui::BeginTable("##gpu_fg_table", 6, ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV)) {
    ImGui::TableSetupColumn("Pass", ImGuiTableColumnFlags_WidthStretch);
    ImGui::TableSetupColumn("Last (ms)", ImGuiTableColumnFlags_WidthFixed);
    ImGui::TableSetupColumn("Avg (ms)", ImGuiTableColumnFlags_WidthFixed);
    ImGui::TableSetupColumn("Max (ms)", ImGuiTableColumnFlags_WidthFixed);
    ImGui::TableSetupColumn("% GPU", ImGuiTableColumnFlags_WidthFixed);
    ImGui::TableSetupColumn("ID", ImGuiTableColumnFlags_WidthFixed);
    ImGui::TableHeadersRow();

    for (const auto& r : rows) {
      ImGui::TableNextRow();

      ImGui::TableSetColumnIndex(0);
      const std::string label = std::string(r.name) + "##gpu_pass_" + std::to_string(r.passId);
      if (ImGui::Selectable(label.c_str(), false, ImGuiSelectableFlags_SpanAllColumns)) {
        ImGui::SetClipboardText(std::string(r.name).c_str());
      }

      ImGui::TableSetColumnIndex(1);
      ImGui::Text("%.3f", r.lastMs);

      ImGui::TableSetColumnIndex(2);
      ImGui::Text("%.3f", r.avgMs);

      ImGui::TableSetColumnIndex(3);
      ImGui::Text("%.3f", r.maxMs);

      ImGui::TableSetColumnIndex(4);
      const double pct = (frameMs > 1e-6) ? (100.0 * r.lastMs / frameMs) : 0.0;
      ImGui::Text("%.1f", pct);

      ImGui::TableSetColumnIndex(5);
      ImGui::Text("%d", (int)r.passId);

    }

    ImGui::EndTable();
  }
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
                        const ToastFn& toast,
                        const render::FrameGraph* frameGraph) {
  // Always keep GPU history in sync even if the window is closed.
  maybeCaptureGpuFrameGraphSample(st, profiler, frameGraph);

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

  if (st.showGpuFrameGraph) {
    if (ImGui::CollapsingHeader("GPU: FrameGraph", ImGuiTreeNodeFlags_DefaultOpen)) {
      if (!frameGraph) {
        ImGui::TextDisabled("(No FrameGraph provided — PostFX may be disabled)");
      } else {
        const bool supported = frameGraph->profilingSupported();
        const bool enabled = frameGraph->profilingEnabled();
        if (!supported) {
          ImGui::TextDisabled("Timer queries unavailable (GPU timings not supported). ");
        } else if (!enabled) {
          ImGui::TextDisabled("GPU timings are available but currently disabled (enable in FrameGraph debug). ");
        } else {
          ImGui::Text("GPU frame time: %.3f ms", frameGraph->lastFrameTimeMs());
        }

        ImGui::Checkbox("Plot", &st.gpuShowPlot);
        ImGui::SameLine();
        ImGui::Checkbox("Table", &st.gpuShowTable);
        ImGui::SameLine();
        ImGui::Checkbox("Sort by last", &st.gpuSortByLast);

        if (st.gpuShowPlot) {
          drawGpuFrameGraphPlot(st);
        }

        if (st.gpuShowTable) {
          ImGui::InputTextWithHint("GPU Filter", "substring...", st.gpuFilter, sizeof(st.gpuFilter));
          drawGpuFrameGraphTable(st, *frameGraph, std::string_view(st.gpuFilter));
          ImGui::TextDisabled("Tip: Click a pass to copy its name.");
        }
      }
    }
    ImGui::Separator();
  }

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

    if (frameGraph) {
      ImGui::Checkbox("Include GPU FrameGraph counters", &st.exportIncludeGpuCounters);
      ImGui::TextDisabled("GPU export uses timer query results (enable GPU pass timing in FrameGraph debug). ");
    } else {
      st.exportIncludeGpuCounters = false;
    }

    core::ChromeTraceWriteOptions opt;
    opt.includeFrameEvents = st.exportIncludeFrameEvents;
    opt.pretty = st.exportPretty;

    if (ImGui::Button("Export##trace")) {
      std::string err;
      bool ok = false;

      // Decide which CPU frames to export.
      std::deque<core::ProfilerFrame> framesOut;
      if (st.exportAllFrames) {
        framesOut = profiler.frames();
      } else {
        framesOut.push_back(*sel);
      }

      // Optional GPU counters.
      if (st.exportIncludeGpuCounters && frameGraph) {
        // Build a counter table aligned to the exported CPU frames.
        core::ChromeTraceCounterTable gpu;
        gpu.name = "FrameGraph GPU";
        gpu.category = "gpu";

        // Use the keys we captured from the live FrameGraph. If they are missing
        // (e.g., profiler enabled recently), fall back to a minimal table.
        if (!st.gpuCounterKeys.empty()) {
          gpu.keys = st.gpuCounterKeys;
        } else {
          gpu.keys = {"GPU Frame (ms)"};
        }

        const std::size_t k = gpu.keys.size();
        const std::size_t n = framesOut.size();
        gpu.tsUs.reserve(n);
        gpu.values.resize(n * k, 0.0);

        const std::uint64_t baseNs = framesOut.empty() ? 0ull : framesOut.front().startNs;

        // Align GPU history rows to the *current* profiler ring buffer.
        // If we have fewer GPU samples than frames, pad with zeros at the front.
        const std::size_t have = st.gpuHistoryRows.size();
        const std::size_t need = profiler.frames().size();
        const std::size_t pad = (need > have) ? (need - have) : 0;

        auto rowForProfilerIndex = [&](std::size_t idx) -> const std::vector<double>* {
          // idx is in [0, profiler.frames().size())
          if (idx < pad) return nullptr;
          const std::size_t hi = idx - pad;
          if (hi >= st.gpuHistoryRows.size()) return nullptr;
          return &st.gpuHistoryRows[hi];
        };

        // For each exported CPU frame, find its index in the current ring buffer.
        // We match by startNs (unique per frame).
        for (std::size_t i = 0; i < n; ++i) {
          const auto& f = framesOut[i];
          const std::uint64_t tsUs = (baseNs && f.startNs >= baseNs) ? ((f.startNs - baseNs) / 1000ull) : 0ull;
          gpu.tsUs.push_back(tsUs);

          // Locate frame in current profiler history.
          std::size_t profIdx = (std::size_t)-1;
          const auto& liveFrames = profiler.frames();
          for (std::size_t j = 0; j < liveFrames.size(); ++j) {
            if (liveFrames[j].startNs == f.startNs) {
              profIdx = j;
              break;
            }
          }

          const std::vector<double>* row = (profIdx != (std::size_t)-1) ? rowForProfilerIndex(profIdx) : nullptr;

          for (std::size_t col = 0; col < k; ++col) {
            const double v = (row && col < row->size()) ? (*row)[col] : 0.0;
            gpu.values[i * k + col] = v;
          }
        }

        ok = core::writeProfilerChromeTraceJsonWithCounters(st.exportPath, framesOut, {gpu}, opt, &err);
      } else {
        ok = core::writeProfilerChromeTraceJson(st.exportPath, framesOut, opt, &err);
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
