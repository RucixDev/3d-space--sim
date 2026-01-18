#pragma once

#include "stellar/core/Types.h"
#include "stellar/ui/HudSettings.h" // Color4f

#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace stellar::ui::textfx {

// A tiny markup-driven text animation/modifier system.
//
// Design goals:
// - Lightweight (no external deps beyond the project)
// - Deterministic (seeded jitter/scramble)
// - Safe to use in gameplay UI strings that already contain bracketed text like "[CONTRABAND]"
//
// Markup syntax (subset):
//   [wave amp=6 freq=0.25 speed=2]Hello[/wave]
//   [shake amp=1.5 rate=18]ALERT[/shake]
//   [pulse min=0.2 max=1 speed=2]ping[/pulse]
//   [color #ff4040]DANGER[/color]
//   [grad #ff00aa #00ccff]NEBULA[/grad]
//   [rainbow speed=0.4 freq=0.08]STARFORGE[/rainbow]
//   [scramble amount=1 rate=25 set=hex]HACKING...[/scramble]
//   [type cps=28 start=0 fade=0.05]Incoming transmission...[/type]
//
// Unrecognized tags are treated as literal text (e.g. "[CONTRABAND]" remains visible).

enum class SpanKind : core::u8 {
  Color = 0,
  Gradient,
  Wave,
  Shake,
  Pulse,
  Rainbow,
  Scramble,
  Typewriter,
};

struct WaveParams {
  float ampPx{4.0f};
  float freq{0.25f};  // cycles per glyph
  float speed{1.5f};  // cycles per second
  float phase{0.0f};
};

struct ShakeParams {
  float ampPx{1.25f};
  float rateHz{18.0f};
  core::u64 seedOffset{0};
};

struct PulseParams {
  float minA{0.25f};
  float maxA{1.0f};
  float speed{1.5f};
  float phase{0.0f};
};

struct ColorParams {
  Color4f color{1, 1, 1, 1};
};

struct GradientParams {
  Color4f from{1, 1, 1, 1};
  Color4f to{1, 1, 1, 1};
};

struct RainbowParams {
  float freq{0.10f};   // hue cycles per glyph
  float speed{0.35f};  // hue cycles per second
  float sat{0.85f};
  float val{1.0f};
};

struct ScrambleParams {
  float amount{1.0f};  // probability [0..1] of scrambling a glyph
  float rateHz{22.0f}; // new glyph each 1/rateHz seconds
  std::string charset{"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"};
};

struct TypewriterParams {
  float cps{30.0f};    // characters per second
  float start{0.0f};   // seconds
  float fade{0.04f};   // seconds to fade-in per glyph (0 disables)
};

struct Span {
  SpanKind kind{SpanKind::Wave};
  int startGlyph{0}; // inclusive
  int endGlyph{0};   // exclusive

  // Param payload (only the matching member for 'kind' is used).
  WaveParams wave{};
  ShakeParams shake{};
  PulseParams pulse{};
  ColorParams color{};
  GradientParams gradient{};
  RainbowParams rainbow{};
  ScrambleParams scramble{};
  TypewriterParams typewriter{};
};

struct Program {
  std::string plain;              // UTF-8 plain text with markup removed.
  std::vector<Span> spans;        // modifier spans referencing glyph indices of 'plain'
  std::vector<std::string> warn;  // non-fatal parse warnings
  int glyphCount{0};              // UTF-8 codepoint count in `plain`
};

// A small LRU cache for compiled markup Programs.
//
// Motivation:
// - Many UI widgets redraw the same markup strings every frame.
// - Compiling markup allocates and parses tags/spans.
//
// This cache keeps a bounded set of compiled Programs so callers can request
// `get(markup)` each frame without repeatedly parsing on hits.
class ProgramCache {
 public:
  // `capacity` is the maximum number of distinct markup strings to keep.
  explicit ProgramCache(std::size_t capacity = 256);

  void setCapacity(std::size_t capacity);
  std::size_t capacity() const { return capacity_; }
  std::size_t size() const { return map_.size(); }

  void clear();

  // Get a compiled Program for the given markup string (compiles on miss).
  //
  // The returned reference remains valid until the cache evicts the entry or is cleared.
  const Program& get(std::string_view markup);

 private:
  struct Entry {
    Program program;
    core::u64 lastUse{0};
  };

  // Heterogeneous lookup so callers can query with std::string_view without allocating.
  struct Hash {
    using is_transparent = void;
    std::size_t operator()(std::string_view s) const noexcept;
    std::size_t operator()(const std::string& s) const noexcept;
  };

  struct Eq {
    using is_transparent = void;
    bool operator()(std::string_view a, std::string_view b) const noexcept;
    bool operator()(const std::string& a, std::string_view b) const noexcept;
    bool operator()(std::string_view a, const std::string& b) const noexcept;
    bool operator()(const std::string& a, const std::string& b) const noexcept;
  };

  void evictIfNeeded_();

  std::size_t capacity_{256};
  core::u64 tick_{0};
  std::unordered_map<std::string, Entry, Hash, Eq> map_;
};

// Compile markup into a plain string + spans.
Program compile(std::string_view markup);

// Convenience: remove markup (unknown tags are preserved literally).
std::string stripMarkup(std::string_view markup);

// Count UTF-8 codepoints.
int glyphCountUtf8(std::string_view plainUtf8);

struct EvalStyle {
  Color4f baseColor{1, 1, 1, 1};
  core::u64 seed{0};
};

struct GlyphVisual {
  float dx{0.0f};
  float dy{0.0f};
  Color4f color{1, 1, 1, 1};
  char32_t codepoint{0};
  bool visible{true};
};

// Evaluate all modifiers affecting a glyph.
//
// Notes:
// - 'x'/'y' are the base glyph pen position in pixels (can be used by effects)
// - 'glyphIndex' refers to the codepoint index in Program::plain
// - 'timeSec' is an arbitrary animation clock provided by the caller
GlyphVisual evalGlyph(const Program& prog,
                      int glyphIndex,
                      float x,
                      float y,
                      char32_t cp,
                      float timeSec,
                      const EvalStyle& style);

} // namespace stellar::ui::textfx

#if STELLAR_ENABLE_IMGUI

#include <imgui.h>

namespace stellar::ui::textfx {

struct DrawParams {
  ImFont* font{nullptr};
  float fontSizePx{0.0f};     // 0 => ImGui::GetFontSize()
  ImU32 baseColor{0};         // 0 => ImGuiCol_Text
  float wrapWidthPx{0.0f};    // 0 => no wrapping
  core::u64 seed{0};
};

// Draw compiled markup.
void Draw(ImDrawList* drawList, ImVec2 pos, const Program& prog, float timeSec, const DrawParams& params);

// Convenience: compile + draw.
void Draw(ImDrawList* drawList, ImVec2 pos, std::string_view markup, float timeSec, const DrawParams& params);

// Compute the pixel size of the rendered plain text (with wrapping disabled unless wrapWidthPx>0).
ImVec2 CalcSize(const Program& prog, const DrawParams& params);

// Convenience: compile + measure.
ImVec2 CalcSize(std::string_view markup, const DrawParams& params);

} // namespace stellar::ui::textfx

#endif // STELLAR_ENABLE_IMGUI
