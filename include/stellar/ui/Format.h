#pragma once

// UI-facing formatting helpers.
//
// Goals:
//  - keep dependencies light (no fmtlib)
//  - provide ergonomic conversions for common ID/seed types (u64) used in UI
//  - avoid accidental heap churn for per-frame UI labels via a small scratch formatter
//
// IMPORTANT: tmpFmt()/tmpToString() return pointers into a thread-local ring buffer.
// They are intended for immediate-mode UI calls (ImGui label strings, etc.) and are
// only valid until the next tmpFmt/tmpToString call on the same thread.

#include "stellar/core/Types.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <limits>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <vector>

namespace stellar::ui {

// -------------------------------
// toString (UI convenience)
// -------------------------------

inline std::string toString(bool v) {
  return v ? "true" : "false";
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool> && !std::is_enum_v<T>, std::string>
toString(T v) {
  // 64-bit signed range fits in 20 chars incl. sign.
  std::array<char, 32> buf{};
  auto res = std::to_chars(buf.data(), buf.data() + buf.size(), v);
  if (res.ec != std::errc{}) {
    return {};
  }
  return std::string(buf.data(), (std::size_t)(res.ptr - buf.data()));
}

template <typename T>
inline std::enable_if_t<std::is_floating_point_v<T> && !std::is_enum_v<T>, std::string>
toString(T v) {
  if (!std::isfinite((double)v)) {
    return "null";
  }

  // Use a compact, round-trippable-ish format.
  // - float: 9 significant digits
  // - double: 17 significant digits
  const char* fmt = (sizeof(T) <= sizeof(float)) ? "%.9g" : "%.17g";
  std::array<char, 128> buf{};
  const int n = std::snprintf(buf.data(), buf.size(), fmt, (double)v);
  if (n <= 0) {
    return {};
  }
  if ((std::size_t)n < buf.size()) {
    return std::string(buf.data(), (std::size_t)n);
  }

  // Rare: bigger than our stack buffer.
  std::string out;
  out.resize((std::size_t)n + 1);
  std::snprintf(out.data(), out.size(), fmt, (double)v);
  out.resize((std::size_t)n);
  return out;
}

// -------------------------------
// printf-style helpers
// -------------------------------

inline void appendf(std::string& out, const char* fmt, ...) {
  if (!fmt) {
    return;
  }

  va_list args;
  va_start(args, fmt);
  va_list args2;
  va_copy(args2, args);

  int n = 0;
#if defined(_MSC_VER)
  // MSVC's CRT historically prefers _vscprintf for sizing.
  n = _vscprintf(fmt, args2);
#else
  n = std::vsnprintf(nullptr, 0, fmt, args2);
#endif
  va_end(args2);

  if (n <= 0) {
    va_end(args);
    return;
  }

  const std::size_t start = out.size();
  // Reserve room for the formatted string + null terminator.
  out.resize(start + (std::size_t)n + 1);
  std::vsnprintf(out.data() + start, (std::size_t)n + 1, fmt, args);
  va_end(args);

  // Drop the null terminator we reserved.
  out.resize(start + (std::size_t)n);
}

inline std::string format(const char* fmt, ...) {
  std::string out;
  if (!fmt) {
    return out;
  }

  va_list args;
  va_start(args, fmt);
  va_list args2;
  va_copy(args2, args);

  int n = 0;
#if defined(_MSC_VER)
  n = _vscprintf(fmt, args2);
#else
  n = std::vsnprintf(nullptr, 0, fmt, args2);
#endif
  va_end(args2);

  if (n <= 0) {
    va_end(args);
    return out;
  }

  out.resize((std::size_t)n + 1);
  std::vsnprintf(out.data(), (std::size_t)n + 1, fmt, args);
  va_end(args);
  out.resize((std::size_t)n);
  return out;
}

// Small thread-local scratch formatter for ephemeral UI strings.
inline const char* tmpFmt(const char* fmt, ...) {
  thread_local std::array<std::array<char, 512>, 8> ring{};
  thread_local std::size_t slot = 0;

  auto& buf = ring[slot];
  slot = (slot + 1) % ring.size();
  buf[0] = '\0';

  if (!fmt) {
    return buf.data();
  }

  va_list args;
  va_start(args, fmt);
  std::vsnprintf(buf.data(), buf.size(), fmt, args);
  va_end(args);

  buf[buf.size() - 1] = '\0';
  return buf.data();
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool> && !std::is_enum_v<T>, const char*>
tmpToString(T v) {
  thread_local std::array<std::array<char, 48>, 8> ring{};
  thread_local std::size_t slot = 0;

  auto& buf = ring[slot];
  slot = (slot + 1) % ring.size();

  auto res = std::to_chars(buf.data(), buf.data() + buf.size(), v);
  if (res.ec != std::errc{}) {
    buf[0] = '\0';
    return buf.data();
  }
  *res.ptr = '\0';
  return buf.data();
}

// -------------------------------
// Minimal JSON writer
// -------------------------------

class JsonWriter {
public:
  explicit JsonWriter(std::string& out, int indentSpaces = 2)
      : out_(&out), indentSpaces_(indentSpaces) {
    out_->clear();
  }

  void beginObject() {
    beginValue_();
    out_->push_back('{');
    stack_.push_back(Frame{Kind::Object, true});
  }

  void endObject() {
    if (stack_.empty() || stack_.back().kind != Kind::Object) {
      return;
    }
    const bool hadMembers = !stack_.back().first;
    stack_.pop_back();
    afterKey_ = false;
    if (hadMembers) {
      newlineIndent_();
    }
    out_->push_back('}');
  }

  void beginArray() {
    beginValue_();
    out_->push_back('[');
    stack_.push_back(Frame{Kind::Array, true});
  }

  void endArray() {
    if (stack_.empty() || stack_.back().kind != Kind::Array) {
      return;
    }
    const bool hadElems = !stack_.back().first;
    stack_.pop_back();
    afterKey_ = false;
    if (hadElems) {
      newlineIndent_();
    }
    out_->push_back(']');
  }

  void key(std::string_view k) {
    if (stack_.empty() || stack_.back().kind != Kind::Object) {
      return;
    }

    Frame& f = stack_.back();
    if (!f.first) {
      out_->push_back(',');
    }
    out_->push_back('\n');
    indent_((int)stack_.size());
    appendJsonString_(*out_, k);
    out_->append(": ");
    f.first = false;
    afterKey_ = true;
  }

  void nullValue() {
    beginValue_();
    out_->append("null");
  }

  void value(std::string_view v) {
    beginValue_();
    appendJsonString_(*out_, v);
  }

  void value(const char* v) {
    value(std::string_view(v ? v : ""));
  }

  void value(bool v) {
    beginValue_();
    out_->append(v ? "true" : "false");
  }

  template <typename T>
  std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool> && !std::is_enum_v<T>, void>
  value(T v) {
    beginValue_();
    // JSON numbers must not be quoted.
    std::array<char, 32> buf{};
    auto res = std::to_chars(buf.data(), buf.data() + buf.size(), v);
    if (res.ec != std::errc{}) {
      out_->append("0");
      return;
    }
    out_->append(buf.data(), (std::size_t)(res.ptr - buf.data()));
  }

  template <typename T>
  std::enable_if_t<std::is_floating_point_v<T> && !std::is_enum_v<T>, void>
  value(T v) {
    if (!std::isfinite((double)v)) {
      nullValue();
      return;
    }
    beginValue_();
    // Use a compact representation.
    const char* fmt = (sizeof(T) <= sizeof(float)) ? "%.9g" : "%.17g";
    std::array<char, 128> buf{};
    const int n = std::snprintf(buf.data(), buf.size(), fmt, (double)v);
    if (n <= 0) {
      out_->append("0");
      return;
    }
    out_->append(buf.data(), (std::size_t)std::min<int>(n, (int)buf.size() - 1));
  }

  template <typename V>
  void member(std::string_view k, const V& v) {
    key(k);
    value(v);
  }

private:
  enum class Kind { Object, Array };

  struct Frame {
    Kind kind;
    bool first;
  };

  void beginValue_() {
    if (stack_.empty()) {
      afterKey_ = false;
      return;
    }

    Frame& f = stack_.back();
    if (f.kind == Kind::Array) {
      if (!f.first) {
        out_->push_back(',');
      }
      out_->push_back('\n');
      indent_((int)stack_.size());
      f.first = false;
    } else {
      // Objects are handled by key().
      // If a user forgets to call key(), we still emit a newline for readability.
      if (!afterKey_) {
        if (!f.first) {
          out_->push_back(',');
        }
        out_->push_back('\n');
        indent_((int)stack_.size());
        f.first = false;
      }
    }
    afterKey_ = false;
  }

  void newlineIndent_() {
    out_->push_back('\n');
    indent_((int)stack_.size());
  }

  void indent_(int depth) {
    if (indentSpaces_ <= 0) {
      return;
    }
    out_->append((std::size_t)depth * (std::size_t)indentSpaces_, ' ');
  }

  static void appendJsonString_(std::string& out, std::string_view sv) {
    out.push_back('"');
    for (unsigned char c : sv) {
      switch (c) {
        case '"': out.append("\\\""); break;
        case '\\': out.append("\\\\"); break;
        case '\b': out.append("\\b"); break;
        case '\f': out.append("\\f"); break;
        case '\n': out.append("\\n"); break;
        case '\r': out.append("\\r"); break;
        case '\t': out.append("\\t"); break;
        default:
          if (c < 0x20) {
            // Control character -> \u00XX
            static const char* kHex = "0123456789ABCDEF";
            out.append("\\u00");
            out.push_back(kHex[(c >> 4) & 0x0F]);
            out.push_back(kHex[c & 0x0F]);
          } else {
            out.push_back((char)c);
          }
          break;
      }
    }
    out.push_back('"');
  }

  std::string* out_{nullptr};
  int indentSpaces_{2};
  std::vector<Frame> stack_;
  bool afterKey_{false};
};

} // namespace stellar::ui
