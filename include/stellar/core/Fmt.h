#pragma once

// -----------------------------------------------------------------------------
// Minimal {fmt}-style formatting shim
// -----------------------------------------------------------------------------
//
// This project intentionally avoids depending on external formatting libraries.
// However, some game/UI code uses fmt::format("...{}...", args...) style calls.
//
// This header provides a tiny subset sufficient for in-repo usage:
//  - "{}" sequential placeholders
//  - "{:.Nf}" floating point precision (e.g. {:.2f}, {:.0f})
//  - escaped braces "{{" and "}}"
//
// It is NOT a full fmt replacement.

#include <iomanip>
#include <cmath>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>

namespace fmt {
namespace detail {

inline void appendString(std::string& out, std::string_view s) {
  out.append(s.data(), s.size());
}

template <typename T>
inline void appendValue(std::string& out, const T& v, int precision, bool fixed) {
  if constexpr (std::is_same_v<std::decay_t<T>, std::string>) {
    out += v;
  } else if constexpr (std::is_convertible_v<T, std::string_view>) {
    std::string_view sv = v;
    appendString(out, sv);
  } else if constexpr (std::is_same_v<std::decay_t<T>, const char*> ||
                       std::is_same_v<std::decay_t<T>, char*>) {
    out += (v ? v : "");
  } else {
    std::ostringstream ss;
    if (precision >= 0) {
      if (fixed) ss << std::fixed;
      ss << std::setprecision(precision);
    }
    ss << v;
    out += ss.str();
  }
}

template <std::size_t I = 0, typename Tuple, typename F>
inline bool visitTuple(std::size_t idx, Tuple&& t, F&& f) {
  constexpr std::size_t N = std::tuple_size_v<std::remove_reference_t<Tuple>>;
  if constexpr (I >= N) {
    return false;
  } else {
    if (idx == I) {
      f(std::get<I>(t));
      return true;
    }
    return visitTuple<I + 1>(idx, std::forward<Tuple>(t), std::forward<F>(f));
  }
}

inline void parseSpec(std::string_view spec, int& outPrecision, bool& outFixed) {
  // Supports: "" or ":.2f" or ":.0f".
  outPrecision = -1;
  outFixed = false;
  if (spec.empty()) return;

  if (spec.size() >= 2 && spec[0] == ':' && spec[1] == '.') {
    std::size_t i = 2;
    int p = 0;
    bool any = false;
    while (i < spec.size() && spec[i] >= '0' && spec[i] <= '9') {
      any = true;
      p = p * 10 + int(spec[i] - '0');
      ++i;
    }
    if (any) {
      outPrecision = p;
      outFixed = true;
    }
  }
}

} // namespace detail

template <typename... Args>
inline std::string format(std::string_view fmtStr, Args&&... args) {
  std::string out;
  out.reserve(fmtStr.size() + 32);

  auto tup = std::forward_as_tuple(std::forward<Args>(args)...);
  std::size_t argIndex = 0;

  for (std::size_t i = 0; i < fmtStr.size();) {
    const char c = fmtStr[i];
    if (c == '{') {
      // Escaped '{{'
      if (i + 1 < fmtStr.size() && fmtStr[i + 1] == '{') {
        out.push_back('{');
        i += 2;
        continue;
      }

      const std::size_t close = fmtStr.find('}', i + 1);
      if (close == std::string_view::npos) {
        // Unterminated; append the rest verbatim.
        detail::appendString(out, fmtStr.substr(i));
        break;
      }

      const std::string_view inside = fmtStr.substr(i + 1, close - (i + 1));
      int precision = -1;
      bool fixed = false;
      detail::parseSpec(inside, precision, fixed);

      bool ok = detail::visitTuple(argIndex, tup, [&](const auto& v) {
        detail::appendValue(out, v, precision, fixed);
      });
      if (!ok) {
        // Too few args; keep placeholder literal so it's obvious.
        out += "{";
        detail::appendString(out, inside);
        out += "}";
      }
      argIndex++;
      i = close + 1;
      continue;
    }

    if (c == '}') {
      // Escaped '}}'
      if (i + 1 < fmtStr.size() && fmtStr[i + 1] == '}') {
        out.push_back('}');
        i += 2;
        continue;
      }
      // Unmatched '}'
      out.push_back('}');
      i++;
      continue;
    }

    out.push_back(c);
    i++;
  }

  return out;
}

} // namespace fmt
