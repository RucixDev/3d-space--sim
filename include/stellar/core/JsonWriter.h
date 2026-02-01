#pragma once

#include "stellar/core/PathUtf8.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <system_error>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <iosfwd>
#include <ostream>
#include <limits>
#include <locale>
#include <map>
#include <tuple>
#include <set>
#include <unordered_set>
#include <span>
#include <optional>
#include <unordered_map>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
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

  void value(std::nullopt_t) { nullValue(); }

  template <typename T>
  void value(const std::optional<T>& v) {
    if (v) {
      value(*v);
    } else {
      nullValue();
    }
  }

  template <typename T>
  void value(const std::vector<T>& v) {
    beginArray();
    for (const auto& e : v) {
      value(e);
    }
    endArray();
  }

template <typename T, std::size_t N>
void value(const std::array<T, N>& v) {
  beginArray();
  for (const auto& e : v) {
    value(e);
  }
  endArray();
}

template <typename T, std::size_t Extent>
void value(std::span<T, Extent> v) {
  beginArray();
  for (const auto& e : v) {
    value(e);
  }
  endArray();
}

template <typename T, typename Compare, typename Alloc>
void value(const std::set<T, Compare, Alloc>& s) {
  beginArray();
  for (const auto& e : s) {
    value(e);
  }
  endArray();
}

template <typename T, typename Hash, typename KeyEq, typename Alloc>
void value(const std::unordered_set<T, Hash, KeyEq, Alloc>& s) {
  // Keep output deterministic when possible: sort by value if `T` is comparable.
  if constexpr (requires(const T& a, const T& b) { a < b; }) {
    std::vector<const T*> items;
    items.reserve(s.size());
    for (const auto& e : s) {
      items.push_back(&e);
    }
    std::sort(items.begin(), items.end(), [](const T* a, const T* b) { return *a < *b; });

    beginArray();
    for (const T* e : items) {
      value(*e);
    }
    endArray();
  } else {
    // Fallback: preserve the container's iteration order (may be non-deterministic).
    beginArray();
    for (const auto& e : s) {
      value(e);
    }
    endArray();
  }
}


  // Tuple-like containers as JSON arrays.
  template <typename A, typename B>
  void value(const std::pair<A, B>& p) {
    beginArray();
    value(p.first);
    value(p.second);
    endArray();
  }

  template <typename... Ts>
  void value(const std::tuple<Ts...>& t) {
    beginArray();
    std::apply([&](const auto&... elems) { (value(elems), ...); }, t);
    endArray();
  }

  // Variant-like containers: write the held alternative.
  void value(std::monostate) { nullValue(); }

  template <typename... Ts>
  void value(const std::variant<Ts...>& v) {
    if (v.valueless_by_exception()) {
      nullValue();
      return;
    }
    std::visit([&](const auto& alt) { value(alt); }, v);
  }


  // Map-like containers as JSON objects.
  //
  // NOTE: std::unordered_map has non-deterministic iteration order; we sort by
  // key to keep the output stable across runs (one of this writer's goals).
  template <typename V, typename Compare, typename Alloc>
  void value(const std::map<std::string, V, Compare, Alloc>& m) {
    beginObject();
    for (const auto& kv : m) {
      key(kv.first);
      value(kv.second);
    }
    endObject();
  }

  template <typename V, typename Hash, typename KeyEq, typename Alloc>
  void value(const std::unordered_map<std::string, V, Hash, KeyEq, Alloc>& m) {
    beginObject();
    using Map = std::unordered_map<std::string, V, Hash, KeyEq, Alloc>;
    using Pair = typename Map::value_type;
    std::vector<const Pair*> items;
    items.reserve(m.size());
    for (const auto& kv : m) {
      items.push_back(&kv);
    }
    std::sort(items.begin(), items.end(), [](const Pair* a, const Pair* b) {
      return a->first < b->first;
    });
    for (const Pair* kv : items) {
      key(kv->first);
      value(kv->second);
    }
    endObject();
  }


  void value(bool v) {
    beginValue();
    writeRaw(v ? "true" : "false");
  }

  void value(const char* s) {
    if (!s) {
      nullValue();
      return;
    }
    value(std::string_view{s});
  }


  void value(std::string_view s) {
    beginValue();
    writeString(s);
  }

  void value(const std::string& s) { value(std::string_view{s}); }

  void value(const std::filesystem::path& p) {
    value(pathToUtf8String(p));
  }

  template <typename T, typename std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool>, int> = 0>
  void value(T v) {
    beginValue();
    // Use to_chars for deterministic, locale-independent formatting.
    //
    // IMPORTANT: ostream insertion treats `char`, `signed char`, and `unsigned char`
    // as character types (not numbers). Using to_chars ensures 8-bit integer types
    // like uint8_t/int8_t serialize as valid JSON numbers.
    std::array<char, 32> buf{};
    auto res = std::to_chars(buf.data(), buf.data() + buf.size(), v);
    if (res.ec != std::errc{}) {
      writeRaw("0");
      return;
    }
    writeRaw(std::string_view(buf.data(), static_cast<std::size_t>(res.ptr - buf.data())));
  }

  template <typename T, typename std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
  void value(T v) {
    beginValue();

    // JSON has no NaN/Inf (RFC 8259). Clamp to null for non-finite.
    if (!std::isfinite(v)) {
      writeRaw("null");
      return;
    }

    // Prefer std::to_chars for locale-independent, deterministic formatting.
    // Some standard libraries still don't implement floating-point to_chars for
    // all T, so we use a requires-expression and fall back to iostreams.
    if constexpr (requires(char* first, char* last, T x) {
      std::to_chars(first, last, x, std::chars_format::general, std::numeric_limits<T>::max_digits10);
    }) {
      std::array<char, 128> buf{};
      const auto res = std::to_chars(buf.data(),
                                     buf.data() + buf.size(),
                                     v,
                                     std::chars_format::general,
                                     std::numeric_limits<T>::max_digits10);
      if (res.ec == std::errc{}) {
        writeRaw(std::string_view(buf.data(),
                                  static_cast<std::size_t>(res.ptr - buf.data())));
        return;
      }
    }

    // Fallback: iostreams with classic locale.
    std::ostringstream oss;
    oss.imbue(std::locale::classic());
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


  // Only emit a field when a value is present (omit instead of writing null).
  //
  // This is useful for "sparse" JSON objects where absence conveys meaning and
  // also keeps diffs smaller for traces/state dumps.
  template <typename V>
  void fieldIf(std::string_view k, const std::optional<V>& v) {
    if (v) {
      field(k, *v);
    }
  }

  // Pointer overload for non-char types (C-strings have their own overload).
  template <typename V,
            typename std::enable_if_t<!std::is_same_v<std::remove_cv_t<V>, char>, int> = 0>
  void fieldIf(std::string_view k, const V* v) {
    if (v) {
      field(k, *v);
    }
  }

  void fieldIf(std::string_view k, const char* v) {
    if (v) {
      field(k, v);
    }
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


  // -------------------------------
  // RAII scopes (exception-safe)
  // -------------------------------
  //
  // These helpers let callers write:
  //   {
  //     auto obj = w.scopedObject();
  //     w.field("a", 1);
  //   } // endObject() automatically
  class ObjectScope {
  public:
    explicit ObjectScope(JsonWriter& w) : w_(&w) { w_->beginObject(); }
    ObjectScope(JsonWriter& w, std::string_view key) : w_(&w) {
      w_->key(key);
      w_->beginObject();
    }

    ObjectScope(ObjectScope&& other) noexcept : w_(other.w_) { other.w_ = nullptr; }
    ObjectScope& operator=(ObjectScope&& other) noexcept {
      if (this == &other) return *this;
      if (w_) w_->endObject();
      w_ = other.w_;
      other.w_ = nullptr;
      return *this;
    }

    ~ObjectScope() {
      if (w_) w_->endObject();
    }

    ObjectScope(const ObjectScope&) = delete;
    ObjectScope& operator=(const ObjectScope&) = delete;

  private:
    JsonWriter* w_ = nullptr;
  };

  class ArrayScope {
  public:
    explicit ArrayScope(JsonWriter& w) : w_(&w) { w_->beginArray(); }
    ArrayScope(JsonWriter& w, std::string_view key) : w_(&w) {
      w_->key(key);
      w_->beginArray();
    }

    ArrayScope(ArrayScope&& other) noexcept : w_(other.w_) { other.w_ = nullptr; }
    ArrayScope& operator=(ArrayScope&& other) noexcept {
      if (this == &other) return *this;
      if (w_) w_->endArray();
      w_ = other.w_;
      other.w_ = nullptr;
      return *this;
    }

    ~ArrayScope() {
      if (w_) w_->endArray();
    }

    ArrayScope(const ArrayScope&) = delete;
    ArrayScope& operator=(const ArrayScope&) = delete;

  private:
    JsonWriter* w_ = nullptr;
  };

  [[nodiscard]] ObjectScope scopedObject() { return ObjectScope(*this); }
  [[nodiscard]] ObjectScope scopedObject(std::string_view k) { return ObjectScope(*this, k); }
  [[nodiscard]] ArrayScope scopedArray() { return ArrayScope(*this); }
  [[nodiscard]] ArrayScope scopedArray(std::string_view k) { return ArrayScope(*this, k); }


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
