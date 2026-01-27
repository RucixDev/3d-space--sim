#pragma once

#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <iosfwd>
#include <ostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace stellar::core {

// A tiny JSON writer designed for internal tools, save-files and traces.
//
// Goals:
//  - Zero third-party deps.
//  - Deterministic formatting.
//  - Can write either to a std::ostream or to an internal std::string buffer.
//
// NOTE: This is intentionally minimal and does not attempt to validate every
// misuse (like calling value() in an object without calling key() first).
class JsonWriter {
public:
  struct Options {
    bool pretty = false;
    int indentSpaces = 2;
  };

  JsonWriter() : JsonWriter(Options{}) {}

  // Writes into an internal string buffer.
  explicit JsonWriter(Options opt) : opt_(opt), outString_(&buffer_) {}

  // Writes into the given output stream.
  explicit JsonWriter(std::ostream& out, bool pretty = false) : outStream_(&out) {
    // NOTE: This project targets C++20, but MSVC support for some newer
    // initialization syntax can vary by version. Keep this ctor simple and
    // portable by avoiding designated initializers here.
    opt_.pretty = pretty;
    opt_.indentSpaces = 2;
  }

  JsonWriter(std::ostream& out, Options opt) : opt_(opt), outStream_(&out) {}

  JsonWriter(const JsonWriter&) = delete;
  JsonWriter& operator=(const JsonWriter&) = delete;

  std::string takeString() { return std::move(buffer_); }

  void beginObject() {
    beginValue();
    writeChar('{');
    stack_.push_back(Frame{Kind::Object, true});
    ++indentLevel_;
  }

  void endObject() {
    if (stack_.empty() || stack_.back().kind != Kind::Object) {
      // best effort
      writeChar('}');
      afterKey_ = false;
      return;
    }

    const bool hadMembers = !stack_.back().first;
    stack_.pop_back();
    --indentLevel_;

    if (opt_.pretty && hadMembers) {
      newline();
      indent();
    }
    writeChar('}');
    afterKey_ = false;
  }

  void beginArray() {
    beginValue();
    writeChar('[');
    stack_.push_back(Frame{Kind::Array, true});
    ++indentLevel_;
  }

  void endArray() {
    if (stack_.empty() || stack_.back().kind != Kind::Array) {
      // best effort
      writeChar(']');
      afterKey_ = false;
      return;
    }

    const bool hadElements = !stack_.back().first;
    stack_.pop_back();
    --indentLevel_;

    if (opt_.pretty && hadElements) {
      newline();
      indent();
    }
    writeChar(']');
    afterKey_ = false;
  }

  void key(std::string_view k) {
    // Only valid in objects.
    if (stack_.empty() || stack_.back().kind != Kind::Object) {
      return;
    }

    Frame& f = stack_.back();
    if (!f.first) {
      writeChar(',');
    }
    if (opt_.pretty) {
      newline();
      indent();
    }

    writeString(k);
    writeRaw(opt_.pretty ? ": " : ":");
    f.first = false;
    afterKey_ = true;
  }

  void nullValue() {
    beginValue();
    writeRaw("null");
  }

  void value(std::nullptr_t) { nullValue(); }

  void value(bool v) {
    beginValue();
    writeRaw(v ? "true" : "false");
  }

  void value(const char* s) { value(std::string_view{s ? s : ""}); }

  void value(std::string_view s) {
    beginValue();
    writeString(s);
  }

  void value(const std::string& s) { value(std::string_view{s}); }

  template <typename T, typename std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool>, int> = 0>
  void value(T v) {
    beginValue();
    // Avoid locale issues by using an ostringstream with classic formatting.
    std::ostringstream oss;
    oss.imbue(std::locale::classic());
    oss << v;
    writeRaw(oss.str());
  }

  template <typename T, typename std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
  void value(T v) {
    beginValue();
    std::ostringstream oss;
    oss.imbue(std::locale::classic());
    // JSON has no NaN/Inf. Clamp to null for non-finite.
    if (!std::isfinite(v)) {
      writeRaw("null");
      return;
    }
    oss << std::setprecision(std::numeric_limits<T>::max_digits10) << v;
    writeRaw(oss.str());
  }

  // Convenience helpers
  template <typename V, typename std::enable_if_t<!std::is_invocable_v<V>, int> = 0>
  void field(std::string_view k, const V& v) {
    key(k);
    value(v);
  }

  // For writing a field where the value is produced by custom calls (beginObject/array etc).
  template <typename Fn, typename std::enable_if_t<std::is_invocable_v<Fn>, int> = 0>
  void field(std::string_view k, Fn&& fn) {
    key(k);
    std::forward<Fn>(fn)();
  }

  template <typename Fn, typename std::enable_if_t<std::is_invocable_v<Fn>, int> = 0>
  void object(Fn&& fn) {
    beginObject();
    std::forward<Fn>(fn)();
    endObject();
  }

  template <typename Fn, typename std::enable_if_t<std::is_invocable_v<Fn>, int> = 0>
  void object(std::string_view k, Fn&& fn) {
    key(k);
    beginObject();
    std::forward<Fn>(fn)();
    endObject();
  }

  template <typename Fn, typename std::enable_if_t<std::is_invocable_v<Fn>, int> = 0>
  void array(Fn&& fn) {
    beginArray();
    std::forward<Fn>(fn)();
    endArray();
  }

  template <typename Fn, typename std::enable_if_t<std::is_invocable_v<Fn>, int> = 0>
  void array(std::string_view k, Fn&& fn) {
    key(k);
    beginArray();
    std::forward<Fn>(fn)();
    endArray();
  }

private:
  enum class Kind : uint8_t { Object, Array };

  struct Frame {
    Kind kind;
    bool first = true;
  };

  void beginValue() {
    // If we are writing a value for a key in an object, the key() call already
    // emitted separators/newlines.
    if (afterKey_) {
      afterKey_ = false;
      return;
    }

    if (stack_.empty()) {
      return;
    }

    Frame& f = stack_.back();
    if (f.kind == Kind::Array) {
      if (!f.first) {
        writeChar(',');
      }
      if (opt_.pretty) {
        newline();
        indent();
      }
      f.first = false;
    }
    // For objects, values should only happen after key().
  }

  void writeChar(char c) {
    if (outStream_) {
      outStream_->put(c);
      return;
    }
    if (outString_) {
      outString_->push_back(c);
    }
  }

  void writeRaw(std::string_view s) {
    if (outStream_) {
      outStream_->write(s.data(), static_cast<std::streamsize>(s.size()));
      return;
    }
    if (outString_) {
      outString_->append(s.data(), s.size());
    }
  }

  void newline() { writeChar('\n'); }

  void indent() {
    if (!opt_.pretty) {
      return;
    }
    const int spaces = std::max(0, opt_.indentSpaces);
    for (int i = 0; i < indentLevel_ * spaces; ++i) {
      writeChar(' ');
    }
  }

  void writeString(std::string_view s) {
    writeChar('"');
    for (unsigned char uc : s) {
      const char c = static_cast<char>(uc);
      switch (c) {
        case '"':
          writeRaw("\\\"");
          break;
        case '\\':
          writeRaw("\\\\");
          break;
        case '\b':
          writeRaw("\\b");
          break;
        case '\f':
          writeRaw("\\f");
          break;
        case '\n':
          writeRaw("\\n");
          break;
        case '\r':
          writeRaw("\\r");
          break;
        case '\t':
          writeRaw("\\t");
          break;
        default:
          if (uc < 0x20) {
            // Control chars -> \u00XX
            char buf[7]{};
            std::snprintf(buf, sizeof(buf), "\\u%04x", static_cast<unsigned int>(uc));
            writeRaw(buf);
          } else {
            writeChar(c);
          }
          break;
      }
    }
    writeChar('"');
  }

  Options opt_;
  std::ostream* outStream_ = nullptr;
  std::string* outString_ = nullptr;
  std::string buffer_;
  std::vector<Frame> stack_;
  int indentLevel_ = 0;
  bool afterKey_ = false;
};

}  // namespace stellar::core
