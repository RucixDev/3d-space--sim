#include "stellar/core/CVar.h"

#include "stellar/core/Log.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace stellar::core {


/* 
calling c lib functions like std::isspace in loops creates overhead
replaced std::isspace with direct comparisons, its much faster.
also added an early empty check to avoid unnecessary processing (;
and direct pointer access (cache) etc.

code may look less readable but its way faster
*/
static std::string_view trimView(std::string_view s) {
    if (s.empty()) return s;
    
    const char* data = s.data();
    std::size_t size = s.size();
    
    std::size_t b = 0;
    while (b < size && (data[b] == ' ' || data[b] == '\t' || data[b] == '\n' || 
                        data[b] == '\r' || data[b] == '\f' || data[b] == '\v')) {
        ++b;
    }
    
    if (b == size) return std::string_view();
    
    std::size_t e = size;
    while (e > b && (data[e - 1] == ' ' || data[e - 1] == '\t' || data[e - 1] == '\n' || 
                     data[e - 1] == '\r' || data[e - 1] == '\f' || data[e - 1] == '\v')) {
        --e;
    }
    
    return {data + b, e - b};
}

/*
cached size calc 
direct mem access
single alloc with init 
(im just getting lazy with the comment lol)
not a big change but its better (;
*/
static std::string lowerAscii(std::string_view s) {
    std::string out(s.size(), '\0');
    if (out.empty()) return out;
    
    char* dst = &out[0];
    const char* src = s.data();
    const char* end = src + s.size();
    
    while (src < end) {
        *dst++ = static_cast<char>(std::tolower(static_cast<unsigned char>(*src++)));
    }
    
    return out;
}

/*
early size check
removed expensive allocations i mean the original version allocated two full copies othe input strings
etc.
*/
static bool icontains(std::string_view hay, std::string_view needle) {
  if (needle.empty()) return true;
  if (needle.size() > hay.size()) return false;
  
  return std::search(
    hay.begin(), hay.end(),
    needle.begin(), needle.end(),
    [](unsigned char a, unsigned char b) {
      return std::tolower(a) == std::tolower(b);
    }
  ) != hay.end();
}

/*
hear me out, i know it looks BAD but its much faster bcs:
no memory allocations
no full string copies
early checks
only computes lowercase chars 

i mean thats the best i could do with how it work but tbh you never should even get to a stage where you need to use this function
*/

static bool parseBool(std::string_view s, bool& out) {
    auto t = trimView(s);
    size_t l = t.size();
    if (!l || l > 5) return false;
    
    auto c = [&t](size_t i) { return std::tolower(static_cast<unsigned char>(t[i])); };
    
    if (l == 1) {
        char x = c(0);
        if (x == '1' || x == 't') { out = true; return true; }
        if (x == '0' || x == 'f') { out = false; return true; }
        return false;
    }
    if (l == 2 && c(0)=='o' && c(1)=='n') { out = true; return true; }
    if (l == 2 && c(0)=='n' && c(1)=='o') { out = false; return true; }
    if (l == 3 && c(0)=='y' && c(1)=='e' && c(2)=='s') { out = true; return true; }
    if (l == 3 && c(0)=='o' && c(1)=='f' && c(2)=='f') { out = false; return true; }
    if (l == 4 && c(0)=='t' && c(1)=='r' && c(2)=='u' && c(3)=='e') { out = true; return true; }
    if (l == 5 && c(0)=='f' && c(1)=='a' && c(2)=='l' && c(3)=='s' && c(4)=='e') { out = false; return true; }
    return false;
}

/*
early check
stack buffer instead of from_chars
direct errno checking is faster than error_code comp
memcpy is very fast for small strings and null termination
*/

static bool parseInt(std::string_view s, std::int64_t& out) {
    s = trimView(s);
    if (s.empty() || s.size() > 20) return false; // max digits for int64_t is 19 + sign
    
    char buffer[21];
    std::memcpy(buffer, s.data(), s.size());
    buffer[s.size()] = '\0';
    
    char* end;
    errno = 0;
    const long long v = strtoll(buffer, &end, 10);
    
    if (errno || end != buffer + s.size()) return false;
    out = static_cast<std::int64_t>(v);
    return true;
}

/*
removed unnecessary std::string alloc, never copy the enire string to heap just to get null termination
*/

static bool parseFloat(std::string_view s, double& out) {
    s = trimView(s);
    if (s.empty() || s.size() >= 32) return false;
    
    char buf[32];
    std::memcpy(buf, s.data(), s.size());
    buf[s.size()] = '\0';
    
    char* end;
    errno = 0;
    out = std::strtod(buf, &end);
    
    return !errno && end == buf + s.size();
}

/*
fast path for strings without escapes, avoids full processing
removed redundant quote checks and simplified the escape handling wit switch etc.
*/

static std::string unquote(std::string_view s) {
    s = trimView(s);
    if (s.size() >= 2) {
        const char q = s.front();
        if ((q == '"' || q == '\'') && s.back() == q)
            s = s.substr(1, s.size() - 2);
    }
    
    if (s.find('\\') == std::string_view::npos)
        return std::string(s);
    
    std::string out;
    out.reserve(s.size());
    size_t i = 0;
    while (i < s.size()) {
        char c = s[i++];
        if (c == '\\' && i < s.size()) {
            switch (s[i++]) {
                case '\\': case '"': case '\'': out.push_back(s[i-1]); break;
                case 'n': out.push_back('\n'); break;
                case 't': out.push_back('\t'); break;
                default: out.push_back('\\'); out.push_back(s[i-1]); break;
            }
        } else {
            out.push_back(c);
        }
    }
    return out;
}

/*
again std::isspace

alr thats all im tired for today. 
*/

static std::string quoteIfNeeded(std::string_view s) {
  if (s.empty()) return {s};

  bool needs_quote = false;
  for (char c : s) {
    if (c <= ' ' || c == '#' || c == '=' || c == '"' || c == '\\') {
      needs_quote = true;
      break;
    }
  }
  if (!needs_quote) return {s};

  size_t escapes = 0;
  for (char c : s) {
    if (c == '\\' || c == '"' || c == '\n' || c == '\t') ++escapes;
  }

  std::string out;
  out.reserve(s.size() + 2 + escapes);
  out.push_back('"');
  for (char c : s) {
    switch (c) {
      case '\\': out.append("\\\\"); break;
      case '"':  out.append("\\\""); break;
      case '\n': out.append("\\n"); break;
      case '\t': out.append("\\t"); break;
      default:   out.push_back(c); break;
    }
  }
  out.push_back('"');
  return out;
}

static bool isTokenBoundary(char c) {
  return std::isspace((unsigned char)c) || c == '=';
}

static std::string_view stripTrailingComment(std::string_view sv) {
  // Remove trailing "#" or "//" comments, but only if they appear outside quotes.
  // This avoids breaking common string values such as URLs ("http://") or fragments
  // ("#") when they are quoted.
  bool inDouble = false;
  bool inSingle = false;
  bool escape = false;

  for (std::size_t i = 0; i < sv.size(); ++i) {
    const char c = sv[i];

    if (escape) {
      escape = false;
      continue;
    }

    if ((inDouble || inSingle) && c == '\\') {
      escape = true;
      continue;
    }

    // Only treat quotes as starting a quoted region if they appear at a token boundary.
    // This prevents apostrophes inside unquoted tokens (e.g. O'Reilly) from toggling
    // quote state and disabling comment stripping.
    if (!inSingle && c == '"') {
      if (inDouble) {
        inDouble = false;
        continue;
      }
      if (i == 0 || isTokenBoundary(sv[i - 1])) {
        inDouble = true;
        continue;
      }
    }
    if (!inDouble && c == '\'') {
      if (inSingle) {
        inSingle = false;
        continue;
      }
      if (i == 0 || isTokenBoundary(sv[i - 1])) {
        inSingle = true;
        continue;
      }
    }

    if (inDouble || inSingle) continue;

    if (c == '#') {
      return trimView(sv.substr(0, i));
    }

    if (c == '/' && i + 1 < sv.size() && sv[i + 1] == '/') {
      // Treat "//" as a comment start only if it looks like it begins a comment.
      // (i.e. at start of line or following whitespace). This prevents unquoted URLs
      // like http://example.com/a//b from being truncated.
      if (i == 0 || std::isspace((unsigned char)sv[i - 1])) {
        return trimView(sv.substr(0, i));
      }
    }
  }

  return sv;
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

    // Strip trailing comments: "x = y  # comment" (but ignore comment tokens inside quotes)
    sv = stripTrailingComment(sv);

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
      // Preserve pending values exactly as they appeared in the file.
      // (They may already be quoted/escaped, and double-quoting would corrupt them.)
      out << kv.first << " = " << kv.second << "\n";
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
