#include "stellar/core/CVar.h"

#include "stellar/core/Log.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace stellar::core {

static std::string_view trimView(std::string_view s) {
  std::size_t b = 0;
  while (b < s.size() && std::isspace((unsigned char)s[b])) ++b;
  std::size_t e = s.size();
  while (e > b && std::isspace((unsigned char)s[e - 1])) --e;
  return s.substr(b, e - b);
}

static std::string lowerAscii(std::string_view s) {
  std::string out;
  out.resize(s.size());
  for (std::size_t i = 0; i < s.size(); ++i) {
    out[i] = (char)std::tolower((unsigned char)s[i]);
  }
  return out;
}

static bool icontains(std::string_view hay, std::string_view needle) {
  if (needle.empty()) return true;
  const std::string h = lowerAscii(hay);
  const std::string n = lowerAscii(needle);
  return h.find(n) != std::string::npos;
}

static bool parseBool(std::string_view s, bool& out) {
  const std::string k = lowerAscii(trimView(s));
  if (k == "1" || k == "true" || k == "on" || k == "yes") { out = true; return true; }
  if (k == "0" || k == "false" || k == "off" || k == "no") { out = false; return true; }
  return false;
}

static bool parseInt(std::string_view s, std::int64_t& out) {
  s = trimView(s);
  if (s.empty()) return false;

  std::int64_t v = 0;
  const auto* begin = s.data();
  const auto* end = s.data() + s.size();
  const auto res = std::from_chars(begin, end, v, 10);
  if (res.ec != std::errc{} || res.ptr != end) return false;
  out = v;
  return true;
}

static bool parseFloat(std::string_view s, double& out) {
  s = trimView(s);
  if (s.empty()) return false;

  // strtod requires a null-terminated buffer.
  std::string tmp(s);
  char* end = nullptr;
  const double v = std::strtod(tmp.c_str(), &end);
  if (!end) return false;
  if ((std::size_t)(end - tmp.c_str()) != tmp.size()) return false;
  out = v;
  return true;
}

static std::string unquote(std::string_view s) {
  s = trimView(s);
  if (s.size() >= 2) {
    const char q0 = s.front();
    const char q1 = s.back();
    if ((q0 == '"' && q1 == '"') || (q0 == '\'' && q1 == '\'')) {
      s = s.substr(1, s.size() - 2);
    }
  }
  // Minimal escape handling: \" \\ \n \t
  std::string out;
  out.reserve(s.size());
  for (std::size_t i = 0; i < s.size(); ++i) {
    const char c = s[i];
    if (c == '\\' && i + 1 < s.size()) {
      const char n = s[i + 1];
      if (n == '\\' || n == '"' || n == '\'') { out.push_back(n); ++i; continue; }
      if (n == 'n') { out.push_back('\n'); ++i; continue; }
      if (n == 't') { out.push_back('\t'); ++i; continue; }
    }
    out.push_back(c);
  }
  return out;
}

static std::string quoteIfNeeded(std::string_view s) {
  bool needs = false;
  for (char c : s) {
    if (std::isspace((unsigned char)c) || c == '#' || c == '=' || c == '"' || c == '\\') {
      needs = true;
      break;
    }
  }
  if (!needs) return std::string(s);

  std::string out;
  out.reserve(s.size() + 8);
  out.push_back('"');
  for (char c : s) {
    switch (c) {
      case '\\': out += "\\\\"; break;
      case '"': out += "\\\""; break;
      case '\n': out += "\\n"; break;
      case '\t': out += "\\t"; break;
      default: out.push_back(c); break;
    }
  }
  out.push_back('"');
  return out;
}

const CVar* CVarRegistry::find(std::string_view name) const {
  std::lock_guard<std::mutex> lock(mutex_);
  const auto it = vars_.find(name);
  return (it != vars_.end()) ? &it->second : nullptr;
}

CVar* CVarRegistry::find(std::string_view name) {
  std::lock_guard<std::mutex> lock(mutex_);
  auto it = vars_.find(name);
  return (it != vars_.end()) ? &it->second : nullptr;
}

CVar* CVarRegistry::defineImpl(std::string_view name, CVarType type, CVarValue def,
                               std::uint32_t flags, std::string_view help) {
  std::lock_guard<std::mutex> lock(mutex_);

  auto it = vars_.find(name);
  if (it != vars_.end()) {
    if (it->second.type != type) {
      return nullptr;
    }
    // Keep existing current value, but allow help/flags/default to be updated.
    it->second.flags = flags;
    if (!help.empty()) it->second.help = std::string(help);
    it->second.defaultValue = def;

    applyPendingLocked(it->second);
    return &it->second;
  }

  CVar v;
  v.name = std::string(name);
  v.type = type;
  v.flags = flags;
  v.help = std::string(help);
  v.value = def;
  v.defaultValue = def;

  auto [insIt, ok] = vars_.emplace(v.name, std::move(v));
  (void)ok;
  applyPendingLocked(insIt->second);
  return &insIt->second;
}

CVar* CVarRegistry::defineBool(std::string_view name, bool defaultValue,
                               std::uint32_t flags, std::string_view help) {
  return defineImpl(name, CVarType::Bool, CVarValue{defaultValue}, flags, help);
}

CVar* CVarRegistry::defineInt(std::string_view name, std::int64_t defaultValue,
                              std::uint32_t flags, std::string_view help) {
  return defineImpl(name, CVarType::Int, CVarValue{defaultValue}, flags, help);
}

CVar* CVarRegistry::defineFloat(std::string_view name, double defaultValue,
                                std::uint32_t flags, std::string_view help) {
  return defineImpl(name, CVarType::Float, CVarValue{defaultValue}, flags, help);
}

CVar* CVarRegistry::defineString(std::string_view name, std::string defaultValue,
                                 std::uint32_t flags, std::string_view help) {
  return defineImpl(name, CVarType::String, CVarValue{std::move(defaultValue)}, flags, help);
}

bool CVarRegistry::setValueImpl(std::string_view name, const CVarValue& v, std::string* outError) {
  CVar* var = nullptr;
  std::vector<CVarListener> listeners;

  {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = vars_.find(name);
    if (it == vars_.end()) {
      if (outError) *outError = "Unknown cvar: " + std::string(name);
      return false;
    }
    if ((it->second.flags & CVar_ReadOnly) != 0u) {
      if (outError) *outError = "CVar is read-only: " + it->second.name;
      return false;
    }

    // Validate type.
    const bool ok =
      (it->second.type == CVarType::Bool   && std::holds_alternative<bool>(v)) ||
      (it->second.type == CVarType::Int    && std::holds_alternative<std::int64_t>(v)) ||
      (it->second.type == CVarType::Float  && std::holds_alternative<double>(v)) ||
      (it->second.type == CVarType::String && std::holds_alternative<std::string>(v));
    if (!ok) {
      if (outError) *outError = "Type mismatch for cvar: " + it->second.name;
      return false;
    }

    it->second.value = v;
    var = &it->second;
    listeners = it->second.listeners;
  }

  // Notify outside lock.
  if (var) {
    for (const auto& cb : listeners) {
      if (cb) cb(*var);
    }
  }
  return true;
}

bool CVarRegistry::setBool(std::string_view name, bool v, std::string* outError) {
  return setValueImpl(name, CVarValue{v}, outError);
}

bool CVarRegistry::setInt(std::string_view name, std::int64_t v, std::string* outError) {
  return setValueImpl(name, CVarValue{v}, outError);
}

bool CVarRegistry::setFloat(std::string_view name, double v, std::string* outError) {
  return setValueImpl(name, CVarValue{v}, outError);
}

bool CVarRegistry::setString(std::string_view name, std::string v, std::string* outError) {
  return setValueImpl(name, CVarValue{std::move(v)}, outError);
}

bool CVarRegistry::setFromString(std::string_view name, std::string_view value, std::string* outError) {
  // Type depends on existing variable.
  CVarType type;
  {
    std::lock_guard<std::mutex> lock(mutex_);
    const auto it = vars_.find(name);
    if (it == vars_.end()) {
      if (outError) *outError = "Unknown cvar: " + std::string(name);
      return false;
    }
    type = it->second.type;
  }

  switch (type) {
    case CVarType::Bool: {
      bool b = false;
      if (!parseBool(value, b)) {
        if (outError) *outError = "Invalid bool: " + std::string(value);
        return false;
      }
      return setBool(name, b, outError);
    }
    case CVarType::Int: {
      std::int64_t i = 0;
      if (!parseInt(value, i)) {
        if (outError) *outError = "Invalid int: " + std::string(value);
        return false;
      }
      return setInt(name, i, outError);
    }
    case CVarType::Float: {
      double f = 0.0;
      if (!parseFloat(value, f)) {
        if (outError) *outError = "Invalid float: " + std::string(value);
        return false;
      }
      return setFloat(name, f, outError);
    }
    case CVarType::String: {
      return setString(name, unquote(value), outError);
    }
  }
  if (outError) *outError = "Unknown cvar type.";
  return false;
}

bool CVarRegistry::reset(std::string_view name, std::string* outError) {
  CVarValue def;
  {
    std::lock_guard<std::mutex> lock(mutex_);
    const auto it = vars_.find(name);
    if (it == vars_.end()) {
      if (outError) *outError = "Unknown cvar: " + std::string(name);
      return false;
    }
    def = it->second.defaultValue;
  }
  return setValueImpl(name, def, outError);
}

bool CVarRegistry::addListener(std::string_view name, CVarListener cb, std::string* outError) {
  std::lock_guard<std::mutex> lock(mutex_);
  auto it = vars_.find(name);
  if (it == vars_.end()) {
    if (outError) *outError = "Unknown cvar: " + std::string(name);
    return false;
  }
  it->second.listeners.push_back(std::move(cb));
  return true;
}

bool CVarRegistry::getBool(std::string_view name, bool fallback) const {
  std::lock_guard<std::mutex> lock(mutex_);
  const auto it = vars_.find(name);
  if (it == vars_.end()) return fallback;
  if (!std::holds_alternative<bool>(it->second.value)) return fallback;
  return std::get<bool>(it->second.value);
}

std::int64_t CVarRegistry::getInt(std::string_view name, std::int64_t fallback) const {
  std::lock_guard<std::mutex> lock(mutex_);
  const auto it = vars_.find(name);
  if (it == vars_.end()) return fallback;
  if (!std::holds_alternative<std::int64_t>(it->second.value)) return fallback;
  return std::get<std::int64_t>(it->second.value);
}

double CVarRegistry::getFloat(std::string_view name, double fallback) const {
  std::lock_guard<std::mutex> lock(mutex_);
  const auto it = vars_.find(name);
  if (it == vars_.end()) return fallback;
  if (!std::holds_alternative<double>(it->second.value)) return fallback;
  return std::get<double>(it->second.value);
}

std::string CVarRegistry::getString(std::string_view name, std::string_view fallback) const {
  std::lock_guard<std::mutex> lock(mutex_);
  const auto it = vars_.find(name);
  if (it == vars_.end()) return std::string(fallback);
  if (!std::holds_alternative<std::string>(it->second.value)) return std::string(fallback);
  return std::get<std::string>(it->second.value);
}

std::vector<const CVar*> CVarRegistry::list(std::string_view filter) const {
  std::vector<const CVar*> out;
  std::lock_guard<std::mutex> lock(mutex_);
  out.reserve(vars_.size());
  for (const auto& kv : vars_) {
    if (!filter.empty() && !icontains(kv.first, filter)) continue;
    out.push_back(&kv.second);
  }
  return out;
}

const char* CVarRegistry::typeName(CVarType t) {
  switch (t) {
    case CVarType::Bool: return "bool";
    case CVarType::Int: return "int";
    case CVarType::Float: return "float";
    case CVarType::String: return "string";
  }
  return "?";
}

std::string CVarRegistry::valueToString(const CVar& v) {
  switch (v.type) {
    case CVarType::Bool:
      return std::get<bool>(v.value) ? "true" : "false";
    case CVarType::Int:
      return std::to_string(std::get<std::int64_t>(v.value));
    case CVarType::Float: {
      std::ostringstream oss;
      oss.setf(std::ios::fixed);
      oss.precision(6);
      oss << std::get<double>(v.value);
      return oss.str();
    }
    case CVarType::String:
      return std::get<std::string>(v.value);
  }
  return {};
}

void CVarRegistry::applyPendingLocked(CVar& var) {
  auto pit = pending_.find(var.name);
  if (pit == pending_.end()) return;

  const std::string pendingVal = pit->second;
  pending_.erase(pit);

  // Apply using the var's type.
  std::string err;
  // Avoid deadlock: setFromString takes the lock, so do parsing here while still locked.
  switch (var.type) {
    case CVarType::Bool: {
      bool b = false;
      if (parseBool(pendingVal, b)) var.value = b;
      break;
    }
    case CVarType::Int: {
      std::int64_t i = 0;
      if (parseInt(pendingVal, i)) var.value = i;
      break;
    }
    case CVarType::Float: {
      double f = 0.0;
      if (parseFloat(pendingVal, f)) var.value = f;
      break;
    }
    case CVarType::String:
      var.value = unquote(pendingVal);
      break;
  }
}

bool CVarRegistry::loadFile(const std::string& path, std::string* outError) {
  std::ifstream in(path);
  if (!in) {
    if (outError) *outError = "Failed to open cvar file: " + path;
    return false;
  }

  bool hadErrors = false;
  std::ostringstream errs;

  std::string line;
  int lineNo = 0;
  while (std::getline(in, line)) {
    ++lineNo;

    std::string_view sv(line);
    sv = trimView(sv);
    if (sv.empty()) continue;
    if (sv.rfind("#", 0) == 0) continue;
    if (sv.rfind("//", 0) == 0) continue;

    // Strip trailing comments: "x = y  # comment"
    {
      const std::size_t hash = sv.find('#');
      const std::size_t sl = sv.find("//");
      std::size_t cut = std::string_view::npos;
      if (hash != std::string_view::npos) cut = hash;
      if (sl != std::string_view::npos) cut = std::min(cut, sl);
      if (cut != std::string_view::npos) sv = trimView(sv.substr(0, cut));
    }

    if (sv.empty()) continue;

    std::string_view name;
    std::string_view val;

    const std::size_t eq = sv.find('=');
    if (eq != std::string_view::npos) {
      name = trimView(sv.substr(0, eq));
      val = trimView(sv.substr(eq + 1));
    } else {
      // Split at first whitespace.
      std::size_t sp = 0;
      while (sp < sv.size() && !std::isspace((unsigned char)sv[sp])) ++sp;
      name = trimView(sv.substr(0, sp));
      val = trimView(sp < sv.size() ? sv.substr(sp) : std::string_view{});
    }

    if (name.empty()) continue;

    // Apply or pend.
    bool known = false;
    {
      std::lock_guard<std::mutex> lock(mutex_);
      known = (vars_.find(name) != vars_.end());
      if (!known) {
        pending_[std::string(name)] = std::string(val);
      }
    }

    if (known) {
      std::string err;
      if (!setFromString(name, val, &err)) {
        hadErrors = true;
        errs << path << ":" << lineNo << ": " << err << "\n";
      }
    }
  }

  if (hadErrors && outError) {
    *outError = errs.str();
  }
  return !hadErrors;
}

bool CVarRegistry::saveFile(const std::string& path, std::string* outError) const {
  std::ofstream out(path);
  if (!out) {
    if (outError) *outError = "Failed to write cvar file: " + path;
    return false;
  }

  out << "# StellarForge CVars\n";
  out << "# Generated by cvar.save\n\n";

  std::lock_guard<std::mutex> lock(mutex_);

  for (const auto& kv : vars_) {
    const CVar& v = kv.second;
    if ((v.flags & CVar_Archive) == 0u) continue;

    out << v.name << " = ";
    if (v.type == CVarType::String) {
      out << quoteIfNeeded(std::get<std::string>(v.value));
    } else {
      out << valueToString(v);
    }
    out << "\n";
  }

  if (!pending_.empty()) {
    out << "\n# Pending (unknown at save time)\n";
    for (const auto& kv : pending_) {
      out << kv.first << " = " << quoteIfNeeded(kv.second) << "\n";
    }
  }

  return true;
}

bool CVarRegistry::hasPending(std::string_view name) const {
  std::lock_guard<std::mutex> lock(mutex_);
  return pending_.find(name) != pending_.end();
}

std::optional<std::string> CVarRegistry::pendingValue(std::string_view name) const {
  std::lock_guard<std::mutex> lock(mutex_);
  auto it = pending_.find(name);
  if (it == pending_.end()) return std::nullopt;
  return it->second;
}

CVarRegistry& cvars() {
  static CVarRegistry g;
  return g;
}

static bool parseLogLevelToken(std::string_view s, LogLevel& out) {
  const std::string k = lowerAscii(trimView(s));
  if (k == "trace") { out = LogLevel::Trace; return true; }
  if (k == "debug") { out = LogLevel::Debug; return true; }
  if (k == "info")  { out = LogLevel::Info;  return true; }
  if (k == "warn")  { out = LogLevel::Warn;  return true; }
  if (k == "error") { out = LogLevel::Error; return true; }
  if (k == "off")   { out = LogLevel::Off;   return true; }
  return false;
}

void installDefaultCVars() {
  static bool installed = false;
  if (installed) return;
  installed = true;

  // log.level: keep as a string for human-friendly console usage.
  {
    const auto cur = toString(getLogLevel());
    auto* v = cvars().defineString("log.level",
                                   std::string(cur),
                                   CVar_Archive,
                                   "Global log level: trace|debug|info|warn|error|off");
    if (v) {
      cvars().addListener("log.level", [](const CVar& cv) {
        if (!std::holds_alternative<std::string>(cv.value)) return;
        LogLevel lvl = LogLevel::Info;
        if (!parseLogLevelToken(std::get<std::string>(cv.value), lvl)) {
          // Invalid token: ignore.
          STELLAR_LOG_WARN("cvar log.level: invalid value (expected trace|debug|info|warn|error|off)");
          return;
        }
        setLogLevel(lvl);
      });
    }
  }

  // r.floating_origin.enabled: camera-relative rendering to reduce GPU float jitter
  // in huge coordinate spaces (space sims).
  cvars().defineBool("r.floating_origin.enabled",
                     true,
                     CVar_Archive,
                     "Enable floating origin / camera-relative rendering to reduce GPU float jitter at large distances.");

}

} // namespace stellar::core
