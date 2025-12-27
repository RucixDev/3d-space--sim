#pragma once

#include <cctype>
#include <cstdlib>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace stellar::core {

// Tiny, dependency-free argument parser for CLI tools.
//
// Supports:
//  - Flags:         --flag   -h
//  - KV args:       --key value   --key=value
//  - Positional:    everything else
//
// Notes:
//  - Unknown/duplicate keys are preserved (last() is convenient but values() exists).
//  - This is intentionally small; it's for tools like stellar_sandbox.
class Args {
public:
  Args() = default;
  Args(int argc, char** argv) { parse(argc, argv); }

  // Configure a long option to consume multiple subsequent values.
  // Example: setArity("pos", 3) allows: --pos 0 0 0
  void setArity(std::string_view key, int valueCount) {
    if (valueCount <= 0) return;
    arity_[std::string(key)] = valueCount;
  }

  void parse(int argc, char** argv) {
    program_.clear();
    kv_.clear();
    flags_.clear();
    positional_.clear();

    if (argc > 0 && argv && argv[0]) program_ = argv[0];

    for (int i = 1; i < argc; ++i) {
      const std::string a = argv[i] ? std::string(argv[i]) : std::string();
      if (a.empty()) continue;

      // End-of-options marker: everything after this is positional.
      // This is important for paths/strings that begin with '-' and for negative numbers.
      if (a == "--") {
        for (int j = i + 1; j < argc; ++j) {
          if (argv[j]) positional_.push_back(std::string(argv[j]));
        }
        break;
      }

      // Long options.
      if (startsWith(a, "--")) {
        const auto eq = a.find('=');
        if (eq != std::string::npos) {
          const std::string key = a.substr(2, eq - 2);
          const std::string val = a.substr(eq + 1);
          kv_[key].push_back(val);
          continue;
        }

        const std::string key = a.substr(2);

        const auto ar = arity_.find(key);
        const int need = (ar != arity_.end()) ? ar->second : 1;

        int took = 0;
        while (took < need && i + 1 < argc && argv[i + 1] && !isSwitch(argv[i + 1])) {
          kv_[key].push_back(std::string(argv[++i]));
          ++took;
        }

        // No value(s) -> treat as a flag.
        if (took == 0) flags_.push_back(key);
        continue;
      }

      // Short options (-h, -v, -abc, -o value, -o=value)
      if (startsWith(a, "-") && a.size() >= 2 && a[1] != '-') {
        // If this looks like a number (e.g. -1, -0.5, -.25), treat it as positional.
        // This matches common CLI behavior and avoids surprising "-1 => flags" bugs.
        if (looksLikeNumber(a.c_str())) {
          positional_.push_back(a);
          continue;
        }

        // Support a minimal "-k=value" form for single-letter options.
        // This is convenient for passing values that begin with '-' (e.g. "-o=-").
        const auto eq = a.find('=');
        if (eq == 2) {
          const std::string key = a.substr(1, 1);
          const std::string val = a.substr(eq + 1);
          kv_[key].push_back(val);
          continue;
        }

        // If a single-letter option is configured with an arity, allow "-k v1 v2 ...".
        // Otherwise fall back to grouped short flags (-abc).
        if (a.size() == 2) {
          const std::string key(1, a[1]);
          const auto ar = arity_.find(key);
          const int need = (ar != arity_.end()) ? ar->second : 0;

          if (need > 0) {
            int took = 0;
            while (took < need && i + 1 < argc && argv[i + 1] && !isSwitch(argv[i + 1])) {
              kv_[key].push_back(std::string(argv[++i]));
              ++took;
            }

            // No value(s) -> treat as a flag (consistent with long option behavior).
            if (took == 0) flags_.push_back(key);
            continue;
          }
        }

        for (std::size_t j = 1; j < a.size(); ++j) {
          const char c = a[j];
          if (std::isalnum((unsigned char)c) || c == '_') {
            flags_.push_back(std::string(1, c));
          }
        }
        continue;
      }

      // Positional.

      positional_.push_back(a);
    }
  }

  const std::string& program() const { return program_; }

  bool hasFlag(std::string_view key) const {
    for (const auto& f : flags_) {
      if (f == key) return true;
    }
    return false;
  }

  bool has(std::string_view key) const {
    if (hasFlag(key)) return true;
    return kv_.find(std::string(key)) != kv_.end();
  }

  std::optional<std::string> last(std::string_view key) const {
    const auto it = kv_.find(std::string(key));
    if (it == kv_.end() || it->second.empty()) return std::nullopt;
    return it->second.back();
  }

  std::vector<std::string> values(std::string_view key) const {
    const auto it = kv_.find(std::string(key));
    if (it == kv_.end()) return {};
    return it->second;
  }

  const std::vector<std::string>& flags() const { return flags_; }
  const std::vector<std::string>& positional() const { return positional_; }

  // Typed helpers (return true if provided & parsed).
  bool getU64(std::string_view key, unsigned long long& out) const {
    const auto v = last(key);
    if (!v) return false;
    char* end = nullptr;
    const auto val = std::strtoull(v->c_str(), &end, 10);
    if (end == v->c_str()) return false;
    out = val;
    return true;
  }

  bool getI64(std::string_view key, long long& out) const {
    const auto v = last(key);
    if (!v) return false;
    char* end = nullptr;
    const auto val = std::strtoll(v->c_str(), &end, 10);
    if (end == v->c_str()) return false;
    out = val;
    return true;
  }

  bool getInt(std::string_view key, int& out) const {
    long long v = 0;
    if (!getI64(key, v)) return false;
    out = static_cast<int>(v);
    return true;
  }

  bool getDouble(std::string_view key, double& out) const {
    const auto v = last(key);
    if (!v) return false;
    char* end = nullptr;
    const auto val = std::strtod(v->c_str(), &end);
    if (end == v->c_str()) return false;
    out = val;
    return true;
  }

  bool getString(std::string_view key, std::string& out) const {
    const auto v = last(key);
    if (!v) return false;
    out = *v;
    return true;
  }

private:
  static bool startsWith(const std::string& s, const char* prefix) {
    const std::size_t n = std::char_traits<char>::length(prefix);
    return s.size() >= n && s.compare(0, n, prefix) == 0;
  }

  // A small helper used to disambiguate negative numeric values from switches.
  // Accepts forms like: -1, -0.25, -.5, 1e-3, -2.0E+4
  static bool looksLikeNumber(const char* s) {
    if (!s || !*s) return false;

    int i = 0;
    if (s[i] == '+' || s[i] == '-') ++i;
    if (!s[i]) return false;

    bool anyDigit = false;
    bool anyDot = false;

    for (; s[i]; ++i) {
      const unsigned char c = (unsigned char)s[i];
      if (std::isdigit(c)) {
        anyDigit = true;
        continue;
      }
      if (c == '.' && !anyDot) {
        anyDot = true;
        continue;
      }

      // Scientific notation.
      if ((c == 'e' || c == 'E') && anyDigit) {
        ++i;
        if (s[i] == '+' || s[i] == '-') ++i;
        bool expDigit = false;
        for (; s[i]; ++i) {
          if (std::isdigit((unsigned char)s[i])) {
            expDigit = true;
            continue;
          }
          return false;
        }
        return expDigit;
      }

      return false;
    }

    return anyDigit;
  }

  static bool isSwitch(const char* s) {
    if (!s || !*s) return false;
    if (s[0] != '-') return false;

    // A single '-' is commonly used as a filename placeholder (stdin/stdout).
    // Treat it like a value/operand rather than an option introducer.
    if (s[1] == '\0') return false;

    // A leading '-' is ambiguous with negative numeric values.
    // Treat numeric-looking tokens as values.
    if (looksLikeNumber(s)) return false;
    return true;
  }

  std::string program_;
  std::unordered_map<std::string, int> arity_;
  std::unordered_map<std::string, std::vector<std::string>> kv_;
  std::vector<std::string> flags_;
  std::vector<std::string> positional_;
};

} // namespace stellar::core
