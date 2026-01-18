#include "stellar/ui/TextFx.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <algorithm>
#include <functional>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <sstream>

namespace stellar::ui::textfx {

namespace {

static constexpr float kTwoPi = 6.28318530717958647692f;

// ---- UTF-8 helpers (minimal) ----

static int utf8DecodeOne(const char* s, const char* end, char32_t* outCp) {
  if (!s || s >= end) return 0;
  const unsigned char c0 = static_cast<unsigned char>(*s);
  if (c0 < 0x80) {
    if (outCp) *outCp = static_cast<char32_t>(c0);
    return 1;
  }

  // 2-byte
  if ((c0 & 0xE0u) == 0xC0u) {
    if (s + 1 >= end) return 0;
    const unsigned char c1 = static_cast<unsigned char>(s[1]);
    if ((c1 & 0xC0u) != 0x80u) return 0;
    const char32_t cp = ((char32_t)(c0 & 0x1Fu) << 6) | (char32_t)(c1 & 0x3Fu);
    if (outCp) *outCp = cp;
    return 2;
  }

  // 3-byte
  if ((c0 & 0xF0u) == 0xE0u) {
    if (s + 2 >= end) return 0;
    const unsigned char c1 = static_cast<unsigned char>(s[1]);
    const unsigned char c2 = static_cast<unsigned char>(s[2]);
    if ((c1 & 0xC0u) != 0x80u || (c2 & 0xC0u) != 0x80u) return 0;
    const char32_t cp = ((char32_t)(c0 & 0x0Fu) << 12) | ((char32_t)(c1 & 0x3Fu) << 6) | (char32_t)(c2 & 0x3Fu);
    if (outCp) *outCp = cp;
    return 3;
  }

  // 4-byte
  if ((c0 & 0xF8u) == 0xF0u) {
    if (s + 3 >= end) return 0;
    const unsigned char c1 = static_cast<unsigned char>(s[1]);
    const unsigned char c2 = static_cast<unsigned char>(s[2]);
    const unsigned char c3 = static_cast<unsigned char>(s[3]);
    if ((c1 & 0xC0u) != 0x80u || (c2 & 0xC0u) != 0x80u || (c3 & 0xC0u) != 0x80u) return 0;
    const char32_t cp = ((char32_t)(c0 & 0x07u) << 18) | ((char32_t)(c1 & 0x3Fu) << 12) | ((char32_t)(c2 & 0x3Fu) << 6) |
                        (char32_t)(c3 & 0x3Fu);
    if (outCp) *outCp = cp;
    return 4;
  }

  return 0;
}

static int utf8EncodeOne(char32_t cp, char out[8]) {
  if (!out) return 0;
  if (cp <= 0x7Fu) {
    out[0] = (char)cp;
    return 1;
  }
  if (cp <= 0x7FFu) {
    out[0] = (char)(0xC0u | ((cp >> 6) & 0x1Fu));
    out[1] = (char)(0x80u | (cp & 0x3Fu));
    return 2;
  }
  if (cp <= 0xFFFFu) {
    out[0] = (char)(0xE0u | ((cp >> 12) & 0x0Fu));
    out[1] = (char)(0x80u | ((cp >> 6) & 0x3Fu));
    out[2] = (char)(0x80u | (cp & 0x3Fu));
    return 3;
  }
  // 4-byte
  out[0] = (char)(0xF0u | ((cp >> 18) & 0x07u));
  out[1] = (char)(0x80u | ((cp >> 12) & 0x3Fu));
  out[2] = (char)(0x80u | ((cp >> 6) & 0x3Fu));
  out[3] = (char)(0x80u | (cp & 0x3Fu));
  return 4;
}

static int glyphCountForBytes(std::string_view s) {
  int n = 0;
  const char* p = s.data();
  const char* e = p + s.size();
  while (p < e) {
    char32_t cp = 0;
    int adv = utf8DecodeOne(p, e, &cp);
    if (adv <= 0) {
      // Invalid byte -> treat as a single glyph.
      ++p;
      ++n;
      continue;
    }
    p += adv;
    ++n;
  }
  return n;
}

static std::string toLowerAscii(std::string_view s) {
  std::string out;
  out.reserve(s.size());
  for (char c : s) {
    out.push_back((char)std::tolower((unsigned char)c));
  }
  return out;
}

static void trimAsciiInPlace(std::string& s) {
  std::size_t a = 0;
  while (a < s.size() && std::isspace((unsigned char)s[a])) ++a;
  std::size_t b = s.size();
  while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
  if (a == 0 && b == s.size()) return;
  s = s.substr(a, b - a);
}

struct ParsedTag {
  bool ok{false};
  bool closing{false};
  bool reset{false};
  std::string name;                 // lower-case
  std::vector<std::string> tokens;  // key=value or positional
};

static ParsedTag parseTagContent(std::string_view content) {
  ParsedTag t{};

  std::string s(content);
  trimAsciiInPlace(s);
  if (s.empty()) return t;

  if (s == "reset" || s == "/reset") {
    t.ok = true;
    t.reset = true;
    t.closing = false;
    t.name = "reset";
    return t;
  }

  if (!s.empty() && s[0] == '/') {
    t.closing = true;
    s.erase(s.begin());
    trimAsciiInPlace(s);
  }

  // Tokenize on ASCII whitespace.
  std::vector<std::string> toks;
  {
    std::string cur;
    for (std::size_t i = 0; i < s.size(); ++i) {
      const char c = s[i];
      if (std::isspace((unsigned char)c)) {
        if (!cur.empty()) {
          toks.push_back(cur);
          cur.clear();
        }
        continue;
      }
      cur.push_back(c);
    }
    if (!cur.empty()) toks.push_back(cur);
  }

  if (toks.empty()) return t;

  t.name = toLowerAscii(toks[0]);
  for (std::size_t i = 1; i < toks.size(); ++i) t.tokens.push_back(toks[i]);
  t.ok = true;
  return t;
}

static bool parseFloat(std::string_view s, float* out) {
  if (!out) return false;
  std::string tmp(s);
  trimAsciiInPlace(tmp);
  if (tmp.empty()) return false;
  char* endPtr = nullptr;
  const float v = std::strtof(tmp.c_str(), &endPtr);
  if (!endPtr || endPtr == tmp.c_str()) return false;
  *out = v;
  return true;
}

static bool parseU64(std::string_view s, core::u64* out) {
  if (!out) return false;
  std::string tmp(s);
  trimAsciiInPlace(tmp);
  if (tmp.empty()) return false;
  char* endPtr = nullptr;
  const unsigned long long v = std::strtoull(tmp.c_str(), &endPtr, 0);
  if (!endPtr || endPtr == tmp.c_str()) return false;
  *out = (core::u64)v;
  return true;
}

static bool parseHexByte(char c, int* outNibble) {
  if (!outNibble) return false;
  if (c >= '0' && c <= '9') { *outNibble = c - '0'; return true; }
  if (c >= 'a' && c <= 'f') { *outNibble = 10 + (c - 'a'); return true; }
  if (c >= 'A' && c <= 'F') { *outNibble = 10 + (c - 'A'); return true; }
  return false;
}

static bool parseHexColor(std::string_view s, Color4f* out) {
  if (!out) return false;
  std::string tmp(s);
  trimAsciiInPlace(tmp);
  if (tmp.empty()) return false;
  if (tmp[0] == '#') tmp.erase(tmp.begin());
  if (tmp.size() != 6 && tmp.size() != 8) return false;

  auto byteAt = [&](int idx, int* outByte) -> bool {
    int hi = 0, lo = 0;
    if (!parseHexByte(tmp[(std::size_t)idx * 2 + 0], &hi)) return false;
    if (!parseHexByte(tmp[(std::size_t)idx * 2 + 1], &lo)) return false;
    *outByte = (hi << 4) | lo;
    return true;
  };

  int r = 255, g = 255, b = 255, a = 255;
  if (!byteAt(0, &r) || !byteAt(1, &g) || !byteAt(2, &b)) return false;
  if (tmp.size() == 8) {
    if (!byteAt(3, &a)) return false;
  }

  out->r = (float)r / 255.0f;
  out->g = (float)g / 255.0f;
  out->b = (float)b / 255.0f;
  out->a = (float)a / 255.0f;
  return true;
}

static Color4f lerp(const Color4f& a, const Color4f& b, float t) {
  t = std::clamp(t, 0.0f, 1.0f);
  return Color4f{a.r + (b.r - a.r) * t,
                 a.g + (b.g - a.g) * t,
                 a.b + (b.b - a.b) * t,
                 a.a + (b.a - a.a) * t};
}

static Color4f hsvToRgb(float h01, float s, float v, float a) {
  // h01 in [0,1)
  h01 = std::fmod(h01, 1.0f);
  if (h01 < 0.0f) h01 += 1.0f;
  s = std::clamp(s, 0.0f, 1.0f);
  v = std::clamp(v, 0.0f, 1.0f);

  const float h = h01 * 6.0f;
  const int i = (int)std::floor(h);
  const float f = h - (float)i;
  const float p = v * (1.0f - s);
  const float q = v * (1.0f - s * f);
  const float t = v * (1.0f - s * (1.0f - f));

  float r = v, g = t, b = p;
  switch (i % 6) {
    case 0: r = v; g = t; b = p; break;
    case 1: r = q; g = v; b = p; break;
    case 2: r = p; g = v; b = t; break;
    case 3: r = p; g = q; b = v; break;
    case 4: r = t; g = p; b = v; break;
    case 5: r = v; g = p; b = q; break;
    default: break;
  }

  return Color4f{r, g, b, a};
}

static bool isRecognizedTag(std::string_view nameLower) {
  return nameLower == "color" || nameLower == "c" ||
         nameLower == "grad" || nameLower == "gradient" ||
         nameLower == "wave" || nameLower == "w" ||
         nameLower == "shake" || nameLower == "jitter" ||
         nameLower == "pulse" || nameLower == "p" ||
         nameLower == "rainbow" || nameLower == "rb" ||
         nameLower == "scramble" || nameLower == "glitch" ||
         nameLower == "type" || nameLower == "typewriter" || nameLower == "tw";
}

static SpanKind kindForTag(std::string_view nameLower) {
  if (nameLower == "color" || nameLower == "c") return SpanKind::Color;
  if (nameLower == "grad" || nameLower == "gradient") return SpanKind::Gradient;
  if (nameLower == "wave" || nameLower == "w") return SpanKind::Wave;
  if (nameLower == "shake" || nameLower == "jitter") return SpanKind::Shake;
  if (nameLower == "pulse" || nameLower == "p") return SpanKind::Pulse;
  if (nameLower == "rainbow" || nameLower == "rb") return SpanKind::Rainbow;
  if (nameLower == "scramble" || nameLower == "glitch") return SpanKind::Scramble;
  return SpanKind::Typewriter;
}

static void applyTokensToSpan(const ParsedTag& tag, Span& span, std::vector<std::string>& warn) {
  // Parse tokens with a minimal key=value approach.
  // Unknown tokens are ignored but recorded for debugging.

  auto getKv = [&](std::string_view tok, std::string& outK, std::string& outV) -> bool {
    const std::size_t eq = tok.find('=');
    if (eq == std::string_view::npos) return false;
    outK = toLowerAscii(tok.substr(0, eq));
    outV = std::string(tok.substr(eq + 1));
    return true;
  };

  std::vector<std::string> positional;
  for (const std::string& tok : tag.tokens) {
    std::string k, v;
    if (getKv(tok, k, v)) {
      // key=value
      if (span.kind == SpanKind::Wave) {
        if (k == "amp" || k == "a") parseFloat(v, &span.wave.ampPx);
        else if (k == "freq" || k == "f") parseFloat(v, &span.wave.freq);
        else if (k == "speed" || k == "s") parseFloat(v, &span.wave.speed);
        else if (k == "phase" || k == "ph") parseFloat(v, &span.wave.phase);
      } else if (span.kind == SpanKind::Shake) {
        if (k == "amp" || k == "a") parseFloat(v, &span.shake.ampPx);
        else if (k == "rate" || k == "hz") parseFloat(v, &span.shake.rateHz);
        else if (k == "seed") parseU64(v, &span.shake.seedOffset);
      } else if (span.kind == SpanKind::Pulse) {
        if (k == "min" || k == "mina") parseFloat(v, &span.pulse.minA);
        else if (k == "max" || k == "maxa") parseFloat(v, &span.pulse.maxA);
        else if (k == "speed" || k == "s") parseFloat(v, &span.pulse.speed);
        else if (k == "phase" || k == "ph") parseFloat(v, &span.pulse.phase);
      } else if (span.kind == SpanKind::Color) {
        if (k == "value" || k == "v") {
          Color4f c{};
          if (parseHexColor(v, &c)) span.color.color = c;
        } else if (k == "r") parseFloat(v, &span.color.color.r);
        else if (k == "g") parseFloat(v, &span.color.color.g);
        else if (k == "b") parseFloat(v, &span.color.color.b);
        else if (k == "a") parseFloat(v, &span.color.color.a);
      } else if (span.kind == SpanKind::Gradient) {
        if (k == "from" || k == "a" || k == "c0") {
          Color4f c{};
          if (parseHexColor(v, &c)) span.gradient.from = c;
        } else if (k == "to" || k == "b" || k == "c1") {
          Color4f c{};
          if (parseHexColor(v, &c)) span.gradient.to = c;
        }
      } else if (span.kind == SpanKind::Rainbow) {
        if (k == "freq" || k == "f") parseFloat(v, &span.rainbow.freq);
        else if (k == "speed" || k == "s") parseFloat(v, &span.rainbow.speed);
        else if (k == "sat") parseFloat(v, &span.rainbow.sat);
        else if (k == "val") parseFloat(v, &span.rainbow.val);
      } else if (span.kind == SpanKind::Scramble) {
        if (k == "amount" || k == "a") parseFloat(v, &span.scramble.amount);
        else if (k == "rate" || k == "hz") parseFloat(v, &span.scramble.rateHz);
        else if (k == "set" || k == "charset") {
          span.scramble.charset = v;
        } else if (k == "chars") {
          span.scramble.charset = v;
        }
      } else if (span.kind == SpanKind::Typewriter) {
        if (k == "cps" || k == "speed") parseFloat(v, &span.typewriter.cps);
        else if (k == "start") parseFloat(v, &span.typewriter.start);
        else if (k == "fade") parseFloat(v, &span.typewriter.fade);
      }
      continue;
    }

    // positional
    positional.push_back(tok);
  }

  // positional fallbacks for common tags
  if (span.kind == SpanKind::Color) {
    if (!positional.empty()) {
      Color4f c{};
      if (parseHexColor(positional[0], &c)) span.color.color = c;
    }
  } else if (span.kind == SpanKind::Gradient) {
    if (positional.size() >= 2) {
      Color4f a{}, b{};
      if (parseHexColor(positional[0], &a)) span.gradient.from = a;
      if (parseHexColor(positional[1], &b)) span.gradient.to = b;
    }
  } else if (span.kind == SpanKind::Scramble) {
    // 'set=hex' convenience
    std::string setLower = toLowerAscii(span.scramble.charset);
    if (setLower == "hex") {
      span.scramble.charset = "0123456789ABCDEF";
    } else if (setLower == "alnum") {
      span.scramble.charset = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
    } else if (setLower == "ascii") {
      span.scramble.charset = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
    }
  }

  // Clamp some safety ranges.
  span.wave.ampPx = std::clamp(span.wave.ampPx, 0.0f, 128.0f);
  span.wave.freq = std::clamp(span.wave.freq, 0.0f, 10.0f);
  span.wave.speed = std::clamp(span.wave.speed, 0.0f, 50.0f);

  span.shake.ampPx = std::clamp(span.shake.ampPx, 0.0f, 64.0f);
  span.shake.rateHz = std::clamp(span.shake.rateHz, 0.0f, 120.0f);

  span.pulse.minA = std::clamp(span.pulse.minA, 0.0f, 1.0f);
  span.pulse.maxA = std::clamp(span.pulse.maxA, 0.0f, 1.0f);
  span.pulse.speed = std::clamp(span.pulse.speed, 0.0f, 50.0f);

  span.scramble.amount = std::clamp(span.scramble.amount, 0.0f, 1.0f);
  span.scramble.rateHz = std::clamp(span.scramble.rateHz, 0.0f, 120.0f);
  if (span.scramble.charset.empty()) span.scramble.charset = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

  span.typewriter.cps = std::max(0.0f, span.typewriter.cps);
  span.typewriter.fade = std::max(0.0f, span.typewriter.fade);

  (void)warn;
}

struct OpenSpan {
  std::string nameLower;
  std::size_t spanIndex{0};
};

static void appendPlain(Program& out, std::string_view bytes, int* inOutGlyphIndex) {
  out.plain.append(bytes.data(), bytes.size());
  if (inOutGlyphIndex) *inOutGlyphIndex += glyphCountForBytes(bytes);
}

static bool isWhitespaceCp(char32_t cp) {
  return cp == ' ' || cp == '\t' || cp == '\n' || cp == '\r';
}

} // namespace

Program compile(std::string_view markup) {
  Program out;
  out.plain.reserve(markup.size());

  std::vector<OpenSpan> stack;
  int glyphIndex = 0;

  const char* s = markup.data();
  const char* e = s + markup.size();

  const char* p = s;
  while (p < e) {
    const char c = *p;
    if (c == '\\') {
      // Escape sequence: "\\[" => literal '['
      if (p + 1 < e && p[1] == '[') {
        appendPlain(out, std::string_view(p + 1, 1), &glyphIndex);
        p += 2;
        continue;
      }
    }

    if (c == '[') {
      // Try parse tag.
      const char* close = (const char*)std::memchr(p, ']', (std::size_t)(e - p));
      if (!close) {
        appendPlain(out, std::string_view(p, 1), &glyphIndex);
        ++p;
        continue;
      }

      const std::string_view inside(p + 1, (std::size_t)(close - (p + 1)));
      ParsedTag tag = parseTagContent(inside);
      if (!tag.ok) {
        // Not a tag, treat literally.
        appendPlain(out, std::string_view(p, (std::size_t)(close - p + 1)), &glyphIndex);
        p = close + 1;
        continue;
      }

      if (tag.reset) {
        // Close all open spans.
        for (auto& os : stack) {
          if (os.spanIndex < out.spans.size()) {
            out.spans[os.spanIndex].endGlyph = glyphIndex;
          }
        }
        stack.clear();
        p = close + 1;
        continue;
      }

      if (!isRecognizedTag(tag.name)) {
        // Preserve unknown tags literally so existing gameplay strings like "[CONTRABAND]" remain intact.
        appendPlain(out, std::string_view(p, (std::size_t)(close - p + 1)), &glyphIndex);
        p = close + 1;
        continue;
      }

      if (tag.closing) {
        // Close the most recent matching span.
        bool found = false;
        for (std::size_t si = stack.size(); si-- > 0;) {
          if (stack[si].nameLower == tag.name) {
            const std::size_t spanIdx = stack[si].spanIndex;
            if (spanIdx < out.spans.size()) {
              out.spans[spanIdx].endGlyph = glyphIndex;
            }
            stack.erase(stack.begin() + (std::ptrdiff_t)si);
            found = true;
            break;
          }
        }
        if (!found) {
          out.warn.push_back("TextFx: closing tag without opener: [/" + tag.name + "]");
          // Preserve literally so content isn't lost.
          appendPlain(out, std::string_view(p, (std::size_t)(close - p + 1)), &glyphIndex);
        }
        p = close + 1;
        continue;
      }

      // Opening tag -> create span.
      Span span;
      span.kind = kindForTag(tag.name);
      span.startGlyph = glyphIndex;
      span.endGlyph = glyphIndex; // temporary; fixed on close or end-of-string
      applyTokensToSpan(tag, span, out.warn);
      const std::size_t idx = out.spans.size();
      out.spans.push_back(std::move(span));
      stack.push_back(OpenSpan{tag.name, idx});

      p = close + 1;
      continue;
    }

    // Regular byte / UTF-8 sequence.
    char32_t cp = 0;
    int adv = utf8DecodeOne(p, e, &cp);
    if (adv <= 0) {
      appendPlain(out, std::string_view(p, 1), &glyphIndex);
      ++p;
    } else {
      appendPlain(out, std::string_view(p, (std::size_t)adv), &glyphIndex);
      p += adv;
    }
  }

  // Close any remaining open spans at end of string.
  for (const auto& os : stack) {
    if (os.spanIndex < out.spans.size()) {
      out.spans[os.spanIndex].endGlyph = glyphIndex;
      out.warn.push_back("TextFx: unclosed tag: [" + os.nameLower + "] (auto-closed at end)");
    }
  }

  out.glyphCount = glyphIndex;

  return out;
}

std::string stripMarkup(std::string_view markup) {
  return compile(markup).plain;
}



ProgramCache::ProgramCache(std::size_t capacity) {
  setCapacity(capacity);
}

void ProgramCache::setCapacity(std::size_t capacity) {
  capacity_ = (capacity < 1) ? 1 : capacity;
  evictIfNeeded_();
}

void ProgramCache::clear() {
  map_.clear();
  tick_ = 0;
}

std::size_t ProgramCache::Hash::operator()(std::string_view s) const noexcept {
  return std::hash<std::string_view>{}(s);
}

std::size_t ProgramCache::Hash::operator()(const std::string& s) const noexcept {
  return (*this)(std::string_view(s));
}

bool ProgramCache::Eq::operator()(std::string_view a, std::string_view b) const noexcept {
  return a == b;
}

bool ProgramCache::Eq::operator()(const std::string& a, std::string_view b) const noexcept {
  return std::string_view(a) == b;
}

bool ProgramCache::Eq::operator()(std::string_view a, const std::string& b) const noexcept {
  return a == std::string_view(b);
}

bool ProgramCache::Eq::operator()(const std::string& a, const std::string& b) const noexcept {
  return a == b;
}

const Program& ProgramCache::get(std::string_view markup) {
  // Heterogeneous lookup: avoid allocating on cache hits.
  auto it = map_.find(markup);
  if (it != map_.end()) {
    it->second.lastUse = ++tick_;
    return it->second.program;
  }

  Entry e;
  e.program = compile(markup);
  e.lastUse = ++tick_;

  auto [insIt, inserted] = map_.emplace(std::string(markup), std::move(e));
  (void)inserted;

  evictIfNeeded_();
  return insIt->second.program;
}

void ProgramCache::evictIfNeeded_() {
  // O(N) eviction is OK: this cache is intentionally small and UI-bound.
  while (map_.size() > capacity_) {
    auto victim = map_.begin();
    for (auto it = map_.begin(); it != map_.end(); ++it) {
      if (it->second.lastUse < victim->second.lastUse) victim = it;
    }
    map_.erase(victim);
  }
}
int glyphCountUtf8(std::string_view plainUtf8) {
  return glyphCountForBytes(plainUtf8);
}

static GlyphVisual evalGlyphIndexed_(const Program& prog,
                                 int glyphIndex,
                                 float x,
                                 float y,
                                 char32_t cp,
                                 float timeSec,
                                 const EvalStyle& style,
                                 const int* spanIdx,
                                 std::size_t spanCount) {
  (void)x;
  (void)y;

  GlyphVisual out;
  out.codepoint = cp;
  out.color = style.baseColor;

  // Always treat base alpha as a multiplier so external fade-outs (e.g. toast TTL)
  // apply even when a span overrides RGB.
  const float baseA = std::clamp(style.baseColor.a, 0.0f, 1.0f);

  float alphaMul = 1.0f;
  float dx = 0.0f;
  float dy = 0.0f;

  // Typewriter visibility may hide glyph entirely.
  bool visible = true;

  auto applySpan = [&](const Span& s) {
    if (glyphIndex < s.startGlyph || glyphIndex >= s.endGlyph) return;

    switch (s.kind) {
      case SpanKind::Color: {
        out.color.r = s.color.color.r;
        out.color.g = s.color.color.g;
        out.color.b = s.color.color.b;
        out.color.a = s.color.color.a * baseA;
      } break;
      case SpanKind::Gradient: {
        const int len = std::max(1, s.endGlyph - s.startGlyph);
        const float t = (len <= 1) ? 0.0f : (float)(glyphIndex - s.startGlyph) / (float)(len - 1);
        Color4f c = lerp(s.gradient.from, s.gradient.to, t);
        out.color.r = c.r;
        out.color.g = c.g;
        out.color.b = c.b;
        out.color.a = c.a * baseA;
      } break;
      case SpanKind::Rainbow: {
        const float hue = (float)(glyphIndex - s.startGlyph) * s.rainbow.freq + timeSec * s.rainbow.speed;
        Color4f c = hsvToRgb(hue, s.rainbow.sat, s.rainbow.val, 1.0f);
        out.color.r = c.r;
        out.color.g = c.g;
        out.color.b = c.b;
        out.color.a = baseA;
      } break;
      case SpanKind::Pulse: {
        const float ph = (timeSec * s.pulse.speed + s.pulse.phase) * kTwoPi;
        const float s01 = 0.5f + 0.5f * std::sin(ph);
        const float a = s.pulse.minA + (s.pulse.maxA - s.pulse.minA) * s01;
        alphaMul *= std::clamp(a, 0.0f, 1.0f);
      } break;
      case SpanKind::Wave: {
        const float ph = ((float)(glyphIndex - s.startGlyph) * s.wave.freq + timeSec * s.wave.speed + s.wave.phase) * kTwoPi;
        dy += std::sin(ph) * s.wave.ampPx;
      } break;
      case SpanKind::Shake: {
        const float rate = std::max(0.0f, s.shake.rateHz);
        const core::u64 tSlice = (rate > 0.0f) ? (core::u64)std::floor(std::max(0.0f, timeSec) * (double)rate) : 0ull;
        core::u64 seed = style.seed;
        seed = core::hashCombine(seed, (core::u64)glyphIndex);
        seed = core::hashCombine(seed, tSlice);
        seed = core::hashCombine(seed, s.shake.seedOffset);
        core::SplitMix64 rng(seed);
        const float a = std::max(0.0f, s.shake.ampPx);
        dx += (float)rng.range<double>(-a, a);
        dy += (float)rng.range<double>(-a, a);
      } break;
      case SpanKind::Scramble: {
        if (!isWhitespaceCp(cp) && !s.scramble.charset.empty()) {
          const float rate = std::max(0.0f, s.scramble.rateHz);
          const core::u64 tSlice = (rate > 0.0f) ? (core::u64)std::floor(std::max(0.0f, timeSec) * (double)rate) : 0ull;
          core::u64 seed = style.seed;
          seed = core::hashCombine(seed, (core::u64)glyphIndex);
          seed = core::hashCombine(seed, tSlice);
          core::SplitMix64 rng(seed);

          const float p = std::clamp(s.scramble.amount, 0.0f, 1.0f);
          if (rng.chance(p)) {
            const std::size_t idx = (std::size_t)rng.range<core::u32>(0u, (core::u32)(s.scramble.charset.size() - 1u));
            out.codepoint = (unsigned char)s.scramble.charset[idx];
          }
        }
      } break;
      case SpanKind::Typewriter: {
        const float cps = std::max(0.0f, s.typewriter.cps);
        if (cps <= 0.0f) break;

        const float tLocal = timeSec - s.typewriter.start;
        if (tLocal <= 0.0f) {
          visible = false;
          break;
        }

        const int reveal = (int)std::floor(tLocal * cps);
        const int localIdx = glyphIndex - s.startGlyph;
        if (localIdx > reveal) {
          visible = false;
          break;
        }

        // Fade-in for the last glyph entering.
        if (s.typewriter.fade > 0.0001f) {
          const float tGlyph = tLocal - (float)localIdx / cps;
          const float k = std::clamp(tGlyph / s.typewriter.fade, 0.0f, 1.0f);
          alphaMul *= k;
        }
      } break;
      default: break;
    }
  };

  if (spanIdx && spanCount > 0) {
    for (std::size_t i = 0; i < spanCount; ++i) {
      const int si = spanIdx[i];
      if (si < 0 || (std::size_t)si >= prog.spans.size()) continue;
      applySpan(prog.spans[(std::size_t)si]);
      if (!visible) break;
    }
  } else {
    for (const Span& s : prog.spans) {
      applySpan(s);
      if (!visible) break;
    }
  }

  out.dx = dx;
  out.dy = dy;
  out.visible = visible;

  out.color.a = std::clamp(out.color.a * alphaMul, 0.0f, 1.0f);

  return out;
}

GlyphVisual evalGlyph(const Program& prog,
                      int glyphIndex,
                      float x,
                      float y,
                      char32_t cp,
                      float timeSec,
                      const EvalStyle& style) {
  return evalGlyphIndexed_(prog, glyphIndex, x, y, cp, timeSec, style, nullptr, 0);
}

} // namespace stellar::ui::textfx

#if STELLAR_ENABLE_IMGUI

#include <imgui.h>

namespace stellar::ui::textfx {

namespace {

static ImU32 toU32(Color4f c) {
  auto toByte = [](float x) -> int {
    const int v = (int)std::llround(std::clamp(x, 0.0f, 1.0f) * 255.0f);
    return std::clamp(v, 0, 255);
  };
  return IM_COL32(toByte(c.r), toByte(c.g), toByte(c.b), toByte(c.a));
}

static ImVec2 calcPlainSize(ImFont* font, float fontSize, std::string_view plain, float wrapWidth) {
  if (!font) font = ImGui::GetFont();
  if (fontSize <= 0.0f) fontSize = ImGui::GetFontSize();

  if (wrapWidth <= 0.0f) {
    return ImGui::CalcTextSize(plain.data(), plain.data() + plain.size());
  }

  // Basic wrapping via ImGui helper.
  return ImGui::CalcTextSize(plain.data(), plain.data() + plain.size(), /*hide_text_after_double_hash=*/false, wrapWidth);
}

} // namespace

void Draw(ImDrawList* drawList, ImVec2 pos, const Program& prog, float timeSec, const DrawParams& params) {
  if (!drawList) return;

  ImFont* font = params.font ? params.font : ImGui::GetFont();
  float fontSize = (params.fontSizePx > 0.0f) ? params.fontSizePx : ImGui::GetFontSize();
  const float scale = (font && font->FontSize > 0.0f) ? (fontSize / font->FontSize) : 1.0f;

  // Base color: if not provided, use ImGui default.
  ImU32 baseColU32 = params.baseColor;
  if (baseColU32 == 0) baseColU32 = ImGui::GetColorU32(ImGuiCol_Text);
  ImVec4 baseColV4 = ImGui::ColorConvertU32ToFloat4(baseColU32);

  EvalStyle style;
  style.baseColor = Color4f{baseColV4.x, baseColV4.y, baseColV4.z, baseColV4.w};
  style.seed = params.seed;

  // Layout: simple left-to-right with explicit newlines.
  float penX = 0.0f;
  float penY = 0.0f;
  const float lineH = ImGui::GetTextLineHeightWithSpacing();

  int glyphIndex = 0;

  // Perf: keep a sweep-line active set of spans so we don't scan every span for every glyph.
  // (evalGlyph() is still available for standalone queries.)
  const auto& spans = prog.spans;
  std::vector<int> activeSpans;
  activeSpans.reserve(spans.size());
  std::size_t nextSpanStart = 0;

  const char* s = prog.plain.data();
  const char* e = s + prog.plain.size();

  const float wrapW = params.wrapWidthPx;

  // Cheap wrap: when enabled, we use CalcTextSize on words and wrap on spaces.
  // (Not perfect, but works well for UI microtext; full wrap would require ImGui internal helpers.)
  auto renderGlyph = [&](int idx, float x, float y, char32_t cp) {
    if (cp == '\n' || cp == '\r') return;

    const GlyphVisual gv = evalGlyphIndexed_(prog, idx, x, y, cp, timeSec, style,
                                            activeSpans.empty() ? nullptr : activeSpans.data(),
                                            activeSpans.size());
    if (!gv.visible) return;

    Color4f c = gv.color;
    const ImU32 col = toU32(c);

    // Render via AddText for correct baseline handling.
    char buf[8]{};
    const int n = utf8EncodeOne(gv.codepoint, buf);
    if (n <= 0) return;

    drawList->AddText(font, fontSize, ImVec2(pos.x + x + gv.dx, pos.y + y + gv.dy), col, buf, buf + n);
  };

  auto advanceFor = [&](char32_t cp) -> float {
    if (!font) return fontSize * 0.5f;
    // ImGui fonts are typically limited to 0..0xFFFF unless configured otherwise.
    const unsigned int c = (cp <= 0xFFFFu) ? (unsigned int)cp : (unsigned int)'?';
    const float adv = font->GetCharAdvance((ImWchar)c);
    return adv * scale;
  };

  auto flushNewline = [&]() {
    penX = 0.0f;
    penY += lineH;
  };

  // Wrapping helper: find next word boundary.
  auto nextWordLenBytes = [&](const char* start) -> int {
    const char* it = start;
    while (it < e) {
      char32_t cp = 0;
      int adv = utf8DecodeOne(it, e, &cp);
      if (adv <= 0) break;
      if (cp == '\n' || cp == '\r') break;
      if (cp == ' ' || cp == '\t') break;
      it += adv;
    }
    return (int)(it - start);
  };

  const float maxX = (wrapW > 0.0f) ? wrapW : std::numeric_limits<float>::infinity();

  auto syncActiveSpans = [&](int g) {
    // Add spans that start at or before this glyph index.
    while (nextSpanStart < spans.size() && spans[nextSpanStart].startGlyph <= g) {
      // Skip spans that are already over (can happen with zero-length spans).
      if (spans[nextSpanStart].endGlyph > g) {
        activeSpans.push_back((int)nextSpanStart);
      }
      ++nextSpanStart;
    }

    // Remove ended spans.
    activeSpans.erase(std::remove_if(activeSpans.begin(), activeSpans.end(), [&](int si) {
      if (si < 0 || (std::size_t)si >= spans.size()) return true;
      return spans[(std::size_t)si].endGlyph <= g;
    }), activeSpans.end());
  };

  const char* it = s;
  while (it < e) {
    char32_t cp = 0;
    int advB = utf8DecodeOne(it, e, &cp);
    if (advB <= 0) {
      cp = (unsigned char)*it;
      advB = 1;
    }

    // Update active spans for this glyph (even for newlines/spaces so typewriter timing stays correct).
    syncActiveSpans(glyphIndex);

    if (cp == '\n') {
      flushNewline();
      ++glyphIndex;
      it += advB;
      continue;
    }

    // Wrap check on spaces: look ahead to next word.
    if (wrapW > 0.0f && (cp == ' ' || cp == '\t')) {
      // Render the whitespace first.
      renderGlyph(glyphIndex, penX, penY, cp);
      penX += advanceFor(cp);
      ++glyphIndex;
      it += advB;

      // If the next word doesn't fit, newline.
      const int wordBytes = nextWordLenBytes(it);
      if (wordBytes > 0) {
        const ImVec2 wSz = ImGui::CalcTextSize(it, it + wordBytes);
        if (penX + wSz.x > maxX) flushNewline();
      }
      continue;
    }

    // If wrapping is enabled, do a crude check per-glyph as well.
    if (wrapW > 0.0f && penX > 0.0f && penX + advanceFor(cp) > maxX) {
      flushNewline();
    }

    renderGlyph(glyphIndex, penX, penY, cp);
    penX += advanceFor(cp);

    ++glyphIndex;
    it += advB;
  }
}

void Draw(ImDrawList* drawList, ImVec2 pos, std::string_view markup, float timeSec, const DrawParams& params) {
  Draw(drawList, pos, compile(markup), timeSec, params);
}

ImVec2 CalcSize(const Program& prog, const DrawParams& params) {
  ImFont* font = params.font ? params.font : ImGui::GetFont();
  float fontSize = (params.fontSizePx > 0.0f) ? params.fontSizePx : ImGui::GetFontSize();
  (void)fontSize;
  return calcPlainSize(font, fontSize, prog.plain, params.wrapWidthPx);
}

ImVec2 CalcSize(std::string_view markup, const DrawParams& params) {
  return CalcSize(compile(markup), params);
}

} // namespace stellar::ui::textfx

#endif // STELLAR_ENABLE_IMGUI
