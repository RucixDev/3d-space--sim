#include "TextAnimationLabWindow.h"

#include "stellar/core/Hash.h"

#include <imgui.h>

#include <algorithm>
#include <cstdio>
#include <cctype>
#include <string>
#include <vector>

namespace stellar::game {

namespace {

static const char* kindName(ui::textfx::SpanKind k) {
  using ui::textfx::SpanKind;
  switch (k) {
    case SpanKind::Color: return "Color";
    case SpanKind::Gradient: return "Gradient";
    case SpanKind::Wave: return "Wave";
    case SpanKind::Shake: return "Shake";
    case SpanKind::Pulse: return "Pulse";
    case SpanKind::Rainbow: return "Rainbow";
    case SpanKind::Scramble: return "Scramble";
    case SpanKind::Typewriter: return "Typewriter";
    default: return "?";
  }
}

static void setPreset(TextAnimationLabWindowState& st, const char* text) {
  std::snprintf(st.markup, sizeof(st.markup), "%s", (text ? text : ""));
  st.lastHash = 0; // force recompile
}

} // namespace

void drawTextAnimationLabWindow(TextAnimationLabWindowState& st, float realTimeSec, const ToastFn& toast) {
  if (!st.open) return;

  // Update local animation clock.
  const float dt = (st.lastRealTimeSec > 0.0f) ? (realTimeSec - st.lastRealTimeSec) : 0.0f;
  st.lastRealTimeSec = realTimeSec;
  if (!st.paused) {
    st.timeSec += std::max(0.0f, dt) * std::clamp(st.timeScale, 0.0f, 8.0f);
  }

  ImGui::SetNextWindowSize(ImVec2(900, 680), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Text Animation Lab", &st.open)) {
    ImGui::End();
    return;
  }

  ImGui::TextUnformatted("Markup-driven text animation modifiers (wave / shake / pulse / gradient / rainbow / scramble / typewriter). ");
  ImGui::SameLine();
  ImGui::TextDisabled("Unknown [bracket] tokens are preserved literally.");

  // Presets
  if (ImGui::BeginTable("##presetRow", 6, ImGuiTableFlags_SizingFixedFit)) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    if (ImGui::Button("Wave + Gradient")) {
      setPreset(st,
        "[wave amp=7 freq=0.18 speed=1.9][grad #ff00aa #00ccff]STARFORGE SYSTEM ONLINE[/grad][/wave]\n"
        "[color #ffffff]Docking request received.[/color]");
    }
    ImGui::TableSetColumnIndex(1);
    if (ImGui::Button("Warning Pulse")) {
      setPreset(st,
        "[pulse min=0.15 max=1 speed=2.6][color #ff4444]WARNING:[/color] [color #ffffff]Heat critical[/color][/pulse]\n"
        "[color #999999]Route: [CONTRABAND] lanes detected[/color]");
    }
    ImGui::TableSetColumnIndex(2);
    if (ImGui::Button("Glitch Scramble")) {
      setPreset(st,
        "[scramble amount=0.9 rate=28 set=hex][color #44ff88]AUTH[/color] [color #ffaa44]CHALLENGE[/color] [color #ffffff]...[/color][/scramble]\n"
        "[color #aaaaaa]Handshake: [OK][/color]");
    }
    ImGui::TableSetColumnIndex(3);
    if (ImGui::Button("Typewriter")) {
      setPreset(st,
        "[type cps=28 fade=0.06]Incoming transmission...\n"
        "\"Keep your heat low. Scanner sweeps every 18 seconds.\"[/type]");
      st.timeSec = 0.0f;
    }
    ImGui::TableSetColumnIndex(4);
    if (ImGui::Button("Rainbow")) {
      setPreset(st,
        "[rainbow speed=0.45 freq=0.06]NEBULA TUNNEL CALIBRATION[/rainbow]\n"
        "[color #ffffff]Vector locked.\nProceed.[/color]");
    }
    ImGui::TableSetColumnIndex(5);
    if (ImGui::Button("Reset")) {
      setPreset(st,
        "[wave amp=6 freq=0.20 speed=1.8][grad #ff00aa #00ccff]STARFORGE[/grad][/wave]\n"
        "[pulse min=0.25 max=1 speed=2][color #ffffff]Incoming transmission...[/color][/pulse]\n"
        "[scramble amount=0.85 rate=24 set=hex][color #44ff88]AUTH[/color] [color #ffaa44]CHALLENGE[/color][/scramble]");
      st.timeSec = 0.0f;
    }
    ImGui::EndTable();
  }

  ImGui::Separator();

  // Editor
  bool edited = false;
  ImGui::TextUnformatted("Markup:");
  edited |= ImGui::InputTextMultiline("##markup", st.markup, sizeof(st.markup), ImVec2(-1, 140), ImGuiInputTextFlags_AllowTabInput);

  ImGui::Spacing();

  // Playback controls
  ImGui::Checkbox("Paused", &st.paused);
  ImGui::SameLine();
  ImGui::SetNextItemWidth(160);
  ImGui::SliderFloat("Time scale", &st.timeScale, 0.0f, 4.0f, "%.2fx");
  ImGui::SameLine();
  if (ImGui::Button("Reset time")) st.timeSec = 0.0f;
  ImGui::SameLine();
  if (ImGui::Button("Copy markup")) {
    ImGui::SetClipboardText(st.markup);
    if (toast) toast("Copied markup to clipboard.", 1.1);
  }

  ImGui::SetNextItemWidth(160);
  ImGui::DragFloat("Time (sec)", &st.timeSec, 0.01f, 0.0f, 1.0e6f, "%.2f");

  // Render settings
  ImGui::SeparatorText("Render");

  ImGui::ColorEdit4("Base color", st.baseColor, ImGuiColorEditFlags_NoInputs);

  ImGui::SetNextItemWidth(220);
  ImGui::InputScalar("Seed", ImGuiDataType_U64, &st.seed);

  ImGui::Checkbox("Wrap to preview", &st.wrapToPreview);
  ImGui::SameLine();
  ImGui::Checkbox("Show bounds", &st.showBounds);
  ImGui::SameLine();
  ImGui::Checkbox("Grid", &st.showGrid);

  if (!st.wrapToPreview) {
    ImGui::SetNextItemWidth(220);
    ImGui::DragFloat("Wrap width (px)", &st.wrapWidthPx, 1.0f, 0.0f, 2000.0f, "%.0f");
  }

  // Compile if needed.
  std::string_view sv(st.markup);
  const core::u64 h = core::fnv1a64(sv);
  if (h != st.lastHash || edited) {
    st.prog = ui::textfx::compile(sv);
    st.lastHash = h;
  }

  // Preview canvas
  ImGui::SeparatorText("Preview");
  const float previewH = 260.0f;
  if (ImGui::BeginChild("##textfx_preview", ImVec2(0, previewH), true, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 canvasPos = ImGui::GetCursorScreenPos();
    const ImVec2 canvasSize = ImGui::GetContentRegionAvail();

    const ImVec2 p0 = canvasPos;
    const ImVec2 p1 = ImVec2(canvasPos.x + std::max(1.0f, canvasSize.x), canvasPos.y + std::max(1.0f, canvasSize.y));

    // Background
    dl->AddRectFilled(p0, p1, IM_COL32(8, 9, 11, 255));
    dl->AddRect(p0, p1, IM_COL32(70, 80, 95, 70));

    if (st.showGrid) {
      const float step = 24.0f;
      const ImU32 col = IM_COL32(80, 90, 105, 30);
      for (float x = p0.x; x < p1.x; x += step) dl->AddLine(ImVec2(x, p0.y), ImVec2(x, p1.y), col);
      for (float y = p0.y; y < p1.y; y += step) dl->AddLine(ImVec2(p0.x, y), ImVec2(p1.x, y), col);
    }

    ui::textfx::DrawParams dp;
    dp.baseColor = ImGui::GetColorU32(ImVec4(st.baseColor[0], st.baseColor[1], st.baseColor[2], st.baseColor[3]));
    dp.seed = st.seed;
    dp.wrapWidthPx = 0.0f;

    const float pad = 10.0f;
    float wrapW = 0.0f;
    if (st.wrapToPreview) {
      wrapW = std::max(0.0f, canvasSize.x - pad * 2.0f);
      dp.wrapWidthPx = wrapW;
    } else {
      dp.wrapWidthPx = st.wrapWidthPx;
    }

    const ImVec2 textPos = ImVec2(p0.x + pad, p0.y + pad);

    ui::textfx::Draw(dl, textPos, st.prog, st.timeSec, dp);

    if (st.showBounds) {
      const ImVec2 sz = ui::textfx::CalcSize(st.prog, dp);
      dl->AddRect(textPos, ImVec2(textPos.x + sz.x, textPos.y + sz.y), IM_COL32(180, 210, 255, 120));
    }

    // Consume the remaining child area so the drawing stays behind an invisible item.
    ImGui::InvisibleButton("##canvas_btn", canvasSize);
  }
  ImGui::EndChild();

  // Warnings
  if (!st.prog.warn.empty()) {
    ImGui::SeparatorText("Parser warnings");
    for (const std::string& w : st.prog.warn) {
      ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.35f, 1.0f), "%s", w.c_str());
    }
  }

  // Span inspector
  ImGui::SeparatorText("Spans");
  ImGui::TextDisabled("Plain glyphs: %d    Spans: %zu", st.prog.glyphCount, st.prog.spans.size());

  if (ImGui::BeginTable("##spanTable", 5, ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders | ImGuiTableFlags_SizingFixedFit)) {
    ImGui::TableSetupColumn("#");
    ImGui::TableSetupColumn("Kind");
    ImGui::TableSetupColumn("Start");
    ImGui::TableSetupColumn("End");
    ImGui::TableSetupColumn("Params");
    ImGui::TableHeadersRow();

    for (std::size_t i = 0; i < st.prog.spans.size(); ++i) {
      const ui::textfx::Span& s = st.prog.spans[i];
      ImGui::TableNextRow();

      ImGui::TableSetColumnIndex(0);
      ImGui::Text("%zu", i);

      ImGui::TableSetColumnIndex(1);
      ImGui::TextUnformatted(kindName(s.kind));

      ImGui::TableSetColumnIndex(2);
      ImGui::Text("%d", s.startGlyph);

      ImGui::TableSetColumnIndex(3);
      ImGui::Text("%d", s.endGlyph);

      ImGui::TableSetColumnIndex(4);
      switch (s.kind) {
        case ui::textfx::SpanKind::Wave:
          ImGui::Text("amp=%.1f freq=%.2f speed=%.2f", s.wave.ampPx, s.wave.freq, s.wave.speed);
          break;
        case ui::textfx::SpanKind::Shake:
          ImGui::Text("amp=%.1f rate=%.1f", s.shake.ampPx, s.shake.rateHz);
          break;
        case ui::textfx::SpanKind::Pulse:
          ImGui::Text("min=%.2f max=%.2f speed=%.2f", s.pulse.minA, s.pulse.maxA, s.pulse.speed);
          break;
        case ui::textfx::SpanKind::Color:
          ImGui::Text("rgba=(%.2f,%.2f,%.2f,%.2f)", s.color.color.r, s.color.color.g, s.color.color.b, s.color.color.a);
          break;
        case ui::textfx::SpanKind::Gradient:
          ImGui::Text("from->to");
          break;
        case ui::textfx::SpanKind::Rainbow:
          ImGui::Text("freq=%.2f speed=%.2f", s.rainbow.freq, s.rainbow.speed);
          break;
        case ui::textfx::SpanKind::Scramble:
          ImGui::Text("p=%.2f rate=%.1f set=%s", s.scramble.amount, s.scramble.rateHz, s.scramble.charset.c_str());
          break;
        case ui::textfx::SpanKind::Typewriter:
          ImGui::Text("cps=%.1f start=%.2f fade=%.2f", s.typewriter.cps, s.typewriter.start, s.typewriter.fade);
          break;
        default:
          ImGui::TextUnformatted("-");
          break;
      }
    }

    ImGui::EndTable();
  }

  ImGui::End();
}

} // namespace stellar::game
