#include "stellar/render/ShaderToyParams.h"

#include "stellar/render/Shader.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <sstream>

namespace stellar::render {

namespace {

static std::string_view ltrim(std::string_view s) {
  size_t i = 0;
  while (i < s.size() && std::isspace((unsigned char)s[i])) ++i;
  return s.substr(i);
}

static std::string_view rtrim(std::string_view s) {
  size_t n = s.size();
  while (n > 0 && std::isspace((unsigned char)s[n - 1])) --n;
  return s.substr(0, n);
}

static std::string_view trim(std::string_view s) {
  return rtrim(ltrim(s));
}

static bool startsWith(std::string_view s, std::string_view pfx) {
  return s.size() >= pfx.size() && s.substr(0, pfx.size()) == pfx;
}

static std::string toUpper(std::string_view s) {
  std::string r;
  r.reserve(s.size());
  for (char c : s) r.push_back((char)std::toupper((unsigned char)c));
  return r;
}

static bool isIdentStart(char c) {
  return (c == '_') || std::isalpha((unsigned char)c);
}

static bool isIdentChar(char c) {
  return (c == '_') || std::isalnum((unsigned char)c);
}

static bool isValidIdentifier(std::string_view s) {
  if (s.empty()) return false;
  if (!isIdentStart(s[0])) return false;
  for (size_t i = 1; i < s.size(); ++i) {
    if (!isIdentChar(s[i])) return false;
  }
  return true;
}

// Tokenize a directive payload. Supports quoted strings.
static bool nextToken(std::string_view& s, std::string_view& outTok) {
  s = ltrim(s);
  if (s.empty()) return false;

  if (s.front() == '"') {
    // Quoted token.
    size_t end = 1;
    while (end < s.size() && s[end] != '"') ++end;
    if (end >= s.size()) {
      // Unterminated quote -> consume rest.
      outTok = s.substr(1);
      s = std::string_view{};
      return true;
    }
    outTok = s.substr(1, end - 1);
    s = s.substr(end + 1);
    return true;
  }

  // Unquoted.
  size_t end = 0;
  while (end < s.size() && !std::isspace((unsigned char)s[end])) ++end;
  outTok = s.substr(0, end);
  s = s.substr(end);
  return true;
}

static bool parseFloat(std::string_view s, float& outV) {
  s = trim(s);
  if (s.empty()) return false;

  // from_chars for float is C++17 but may not be implemented consistently on
  // older libstdc++; use strtof on a temp string.
  std::string tmp(s);
  char* end = nullptr;
  const float v = std::strtof(tmp.c_str(), &end);
  if (!end || end == tmp.c_str()) return false;
  outV = v;
  return true;
}

static bool parseInt(std::string_view s, int& outV) {
  s = trim(s);
  if (s.empty()) return false;
  int v = 0;
  auto res = std::from_chars(s.data(), s.data() + s.size(), v);
  if (res.ec != std::errc()) return false;
  outV = v;
  return true;
}

static void clampValueToDef(const ShaderToyParamDef& d, std::array<float, 4>& v) {
  const auto clamp1 = [&](int k) {
    const float lo = d.minValue[k];
    const float hi = d.maxValue[k];
    if (hi > lo) v[k] = std::clamp(v[k], lo, hi);
  };

  switch (d.type) {
    case ShaderToyParamType::Float:
    case ShaderToyParamType::Int:
    case ShaderToyParamType::Bool:
      clamp1(0);
      v[1] = 0.0f;
      v[2] = 0.0f;
      v[3] = 0.0f;
      break;
    case ShaderToyParamType::Vec2:
      clamp1(0);
      clamp1(1);
      v[2] = 0.0f;
      v[3] = 0.0f;
      break;
    case ShaderToyParamType::Vec3:
    case ShaderToyParamType::Color3:
      clamp1(0);
      clamp1(1);
      clamp1(2);
      v[3] = 0.0f;
      break;
  }
}

static std::unordered_map<std::string, ShaderToyParamType> schemaMap(const ShaderToyParamSet& s) {
  std::unordered_map<std::string, ShaderToyParamType> m;
  m.reserve(s.defs.size());
  for (const auto& d : s.defs) {
    m.emplace(d.name, d.type);
  }
  return m;
}

} // namespace

const char* shaderToyParamGlslType(ShaderToyParamType t) {
  switch (t) {
    case ShaderToyParamType::Float: return "float";
    case ShaderToyParamType::Int: return "int";
    case ShaderToyParamType::Bool: return "bool";
    case ShaderToyParamType::Vec2: return "vec2";
    case ShaderToyParamType::Vec3: return "vec3";
    case ShaderToyParamType::Color3: return "vec3";
  }
  return "float";
}

void ShaderToyParamSet::clear() {
  defs.clear();
  values.clear();
  indexByName.clear();
}

void ShaderToyParamSet::rebuildIndex() {
  indexByName.clear();
  indexByName.reserve(defs.size());
  for (int i = 0; i < (int)defs.size(); ++i) {
    indexByName[defs[i].name] = i;
  }
}

void ShaderToyParamSet::resetToDefaults() {
  values.clear();
  values.reserve(defs.size());
  for (const auto& d : defs) {
    std::array<float, 4> v = d.defaultValue;
    clampValueToDef(d, v);
    values.push_back(v);
  }
}

int ShaderToyParamSet::findIndex(std::string_view name) const {
  const auto it = indexByName.find(std::string(name));
  if (it == indexByName.end()) return -1;
  return it->second;
}

const ShaderToyParamDef* ShaderToyParamSet::findDef(std::string_view name) const {
  const int i = findIndex(name);
  if (i < 0 || i >= (int)defs.size()) return nullptr;
  return &defs[i];
}

bool ShaderToyParamSet::setValue(std::string_view name, const std::array<float, 4>& vIn) {
  const int i = findIndex(name);
  if (i < 0 || i >= (int)defs.size() || i >= (int)values.size()) return false;
  std::array<float, 4> v = vIn;
  clampValueToDef(defs[i], v);
  values[i] = v;
  return true;
}

std::string ShaderToyParamSet::buildUniformDecls() const {
  if (defs.empty()) return {};

  std::ostringstream ss;
  ss << "\n// ---- User parameters (Procedural Shader Lab) ----\n";
  for (const auto& d : defs) {
    ss << "uniform " << shaderToyParamGlslType(d.type) << " " << d.name << ";\n";
  }
  ss << "\n";
  return ss.str();
}

void ShaderToyParamSet::applyToShader(const ShaderProgram& shader) const {
  const int n = (int)defs.size();
  if (n == 0) return;
  if ((int)values.size() != n) return;

  for (int i = 0; i < n; ++i) {
    const ShaderToyParamDef& d = defs[i];
    const std::array<float, 4>& v = values[i];

    switch (d.type) {
      case ShaderToyParamType::Float:
        shader.setUniform1f(d.name.c_str(), v[0]);
        break;
      case ShaderToyParamType::Int: {
        const int iv = (int)std::lround(v[0]);
        shader.setUniform1i(d.name.c_str(), iv);
        break;
      }
      case ShaderToyParamType::Bool: {
        const int bv = (v[0] >= 0.5f) ? 1 : 0;
        shader.setUniform1i(d.name.c_str(), bv);
        break;
      }
      case ShaderToyParamType::Vec2:
        shader.setUniform2f(d.name.c_str(), v[0], v[1]);
        break;
      case ShaderToyParamType::Vec3:
      case ShaderToyParamType::Color3:
        shader.setUniform3f(d.name.c_str(), v[0], v[1], v[2]);
        break;
    }
  }
}

bool ShaderToyParamSet::schemaEquals(const ShaderToyParamSet& other) const {
  if (defs.size() != other.defs.size()) return false;
  const auto a = schemaMap(*this);
  const auto b = schemaMap(other);
  if (a.size() != b.size()) return false;
  for (const auto& [k, t] : a) {
    const auto it = b.find(k);
    if (it == b.end()) return false;
    if (it->second != t) return false;
  }
  return true;
}

ShaderToyParamSet parseShaderToyParamsFromSources(std::span<const std::string_view> sources,
                                                  const ShaderToyParamSet* preserveValuesFrom) {
  ShaderToyParamSet out;

  std::string currentGroup;
  out.defs.reserve(32);

  const auto addDefIfNew = [&](ShaderToyParamDef&& d) {
    if (!isValidIdentifier(d.name)) return;
    // First definition wins.
    if (out.indexByName.find(d.name) != out.indexByName.end()) return;

    const int idx = (int)out.defs.size();
    out.indexByName[d.name] = idx;
    out.defs.push_back(std::move(d));
  };

  for (std::string_view src : sources) {
    currentGroup.clear();

    while (!src.empty()) {
      // Consume a line.
      size_t end = src.find('\n');
      std::string_view line = (end == std::string_view::npos) ? src : src.substr(0, end);
      if (end == std::string_view::npos) src = std::string_view{};
      else src = src.substr(end + 1);

      line = rtrim(line);
      std::string_view t = ltrim(line);
      if (!startsWith(t, "//")) continue;

      t.remove_prefix(2);
      t = ltrim(t);
      if (t.empty()) continue;

      if (startsWith(t, "@group")) {
        t.remove_prefix(std::string_view("@group").size());
        t = trim(t);
        // Allow quoted group names.
        std::string_view tok;
        if (nextToken(t, tok)) currentGroup = std::string(tok);
        else currentGroup.clear();
        continue;
      }

      if (startsWith(t, "@endgroup")) {
        currentGroup.clear();
        continue;
      }

      if (!startsWith(t, "@param")) continue;
      t.remove_prefix(std::string_view("@param").size());
      t = trim(t);

      // Parse tokens.
      std::string_view tok;

      // type
      if (!nextToken(t, tok)) continue;
      const std::string typeU = toUpper(tok);
      ShaderToyParamType ty;
      if (typeU == "FLOAT") ty = ShaderToyParamType::Float;
      else if (typeU == "INT") ty = ShaderToyParamType::Int;
      else if (typeU == "BOOL" || typeU == "BOOLEAN") ty = ShaderToyParamType::Bool;
      else if (typeU == "VEC2") ty = ShaderToyParamType::Vec2;
      else if (typeU == "VEC3") ty = ShaderToyParamType::Vec3;
      else if (typeU == "COLOR" || typeU == "COLOR3" || typeU == "RGB") ty = ShaderToyParamType::Color3;
      else continue;

      // name
      if (!nextToken(t, tok)) continue;
      const std::string name = std::string(tok);

      // optional label (quoted)
      std::string label;
      {
        std::string_view probe = ltrim(t);
        if (!probe.empty() && probe.front() == '"') {
          std::string_view lt;
          if (nextToken(t, lt)) label = std::string(lt);
        }
      }
      if (label.empty()) label = name;

      ShaderToyParamDef d;
      d.type = ty;
      d.name = name;
      d.label = label;
      d.group = currentGroup;

      // Parse numeric payload based on type.
      auto needF = [&](float& outF) -> bool {
        std::string_view nt;
        if (!nextToken(t, nt)) return false;
        return parseFloat(nt, outF);
      };
      auto needI = [&](int& outI) -> bool {
        std::string_view nt;
        if (!nextToken(t, nt)) return false;
        return parseInt(nt, outI);
      };

      if (ty == ShaderToyParamType::Float) {
        float def = 0, mn = 0, mx = 1, step = 0.01f;
        if (!needF(def) || !needF(mn) || !needF(mx)) continue;
        std::string_view st;
        if (nextToken(t, st)) {
          float ss = 0;
          if (parseFloat(st, ss)) step = ss;
        }
        d.defaultValue = {def, 0, 0, 0};
        d.minValue = {mn, 0, 0, 0};
        d.maxValue = {mx, 0, 0, 0};
        d.step = step;
      } else if (ty == ShaderToyParamType::Int) {
        int def = 0, mn = 0, mx = 10;
        int stepI = 1;
        if (!needI(def) || !needI(mn) || !needI(mx)) continue;
        std::string_view st;
        if (nextToken(t, st)) {
          int si = 1;
          if (parseInt(st, si)) stepI = std::max(1, si);
        }
        d.defaultValue = {(float)def, 0, 0, 0};
        d.minValue = {(float)mn, 0, 0, 0};
        d.maxValue = {(float)mx, 0, 0, 0};
        d.step = (float)stepI;
      } else if (ty == ShaderToyParamType::Bool) {
        int def = 0;
        (void)needI(def); // default is optional
        def = (def != 0) ? 1 : 0;
        d.defaultValue = {(float)def, 0, 0, 0};
        d.minValue = {0, 0, 0, 0};
        d.maxValue = {1, 1, 1, 1};
        d.step = 1.0f;
      } else if (ty == ShaderToyParamType::Vec2) {
        float dx, dy, mnx, mny, mxx, mxy;
        float step = 0.01f;
        if (!needF(dx) || !needF(dy) || !needF(mnx) || !needF(mny) || !needF(mxx) || !needF(mxy)) continue;
        std::string_view st;
        if (nextToken(t, st)) {
          float ss = 0;
          if (parseFloat(st, ss)) step = ss;
        }
        d.defaultValue = {dx, dy, 0, 0};
        d.minValue = {mnx, mny, 0, 0};
        d.maxValue = {mxx, mxy, 0, 0};
        d.step = step;
      } else if (ty == ShaderToyParamType::Vec3) {
        float dx, dy, dz, mnx, mny, mnz, mxx, mxy, mxz;
        float step = 0.01f;
        if (!needF(dx) || !needF(dy) || !needF(dz)) continue;
        if (!needF(mnx) || !needF(mny) || !needF(mnz)) continue;
        if (!needF(mxx) || !needF(mxy) || !needF(mxz)) continue;
        std::string_view st;
        if (nextToken(t, st)) {
          float ss = 0;
          if (parseFloat(st, ss)) step = ss;
        }
        d.defaultValue = {dx, dy, dz, 0};
        d.minValue = {mnx, mny, mnz, 0};
        d.maxValue = {mxx, mxy, mxz, 0};
        d.step = step;
      } else if (ty == ShaderToyParamType::Color3) {
        float r, g, b;
        if (!needF(r) || !needF(g) || !needF(b)) continue;
        d.defaultValue = {r, g, b, 0};
        d.minValue = {0, 0, 0, 0};
        d.maxValue = {1, 1, 1, 1};
        d.step = 0.01f;
      }

      // Clamp default to range.
      clampValueToDef(d, d.defaultValue);

      addDefIfNew(std::move(d));
    }
  }

  out.resetToDefaults();
  out.rebuildIndex();

  // Preserve values from previous set when possible.
  if (preserveValuesFrom) {
    for (int i = 0; i < (int)out.defs.size(); ++i) {
      const auto& d = out.defs[i];
      // Look up matching param by name.
      const int prevIdx = preserveValuesFrom->findIndex(d.name);
      if (prevIdx < 0) continue;
      if (prevIdx >= (int)preserveValuesFrom->defs.size() || prevIdx >= (int)preserveValuesFrom->values.size()) continue;
      if (preserveValuesFrom->defs[prevIdx].type != d.type) continue;

      std::array<float, 4> v = preserveValuesFrom->values[prevIdx];
      clampValueToDef(d, v);
      out.values[i] = v;
    }
  }

  return out;
}

} // namespace stellar::render
