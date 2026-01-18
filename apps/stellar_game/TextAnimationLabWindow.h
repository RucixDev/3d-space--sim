#pragma once

#include "stellar/core/Types.h"
#include "stellar/ui/TextFx.h"

#include <functional>
#include <string>

namespace stellar::game {

struct TextAnimationLabWindowState {
  bool open{false};

  // Markup source (editable)
  char markup[4096]{
    "[wave amp=6 freq=0.20 speed=1.8][grad #ff00aa #00ccff]STARFORGE[/grad][/wave]\n"
    "[pulse min=0.25 max=1 speed=2][color #ffffff]Incoming transmission...[/color][/pulse]\n"
    "[scramble amount=0.85 rate=24 set=hex][color #44ff88]AUTH[/color] [color #ffaa44]CHALLENGE[/color][/scramble]"};

  // Timeline
  bool paused{false};
  float timeScale{1.0f};
  float timeSec{0.0f};
  float lastRealTimeSec{0.0f};

  // Presentation
  bool wrapToPreview{true};
  float wrapWidthPx{420.0f};
  bool showBounds{true};
  bool showGrid{true};

  float baseColor[4]{0.92f, 0.94f, 1.0f, 1.0f};
  core::u64 seed{0xC0FFEEu};

  // Internal cache
  core::u64 lastHash{0};
  ui::textfx::Program prog{};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

void drawTextAnimationLabWindow(TextAnimationLabWindowState& st, float realTimeSec, const ToastFn& toast);

} // namespace stellar::game
