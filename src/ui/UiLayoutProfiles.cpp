#include "stellar/ui/UiLayoutProfiles.h"

#include "stellar/core/Log.h"

#include <algorithm>
#include <cctype>
#include <filesystem>

namespace stellar::ui {

std::string defaultUiLayoutsDir() {
  return "ui_layouts";
}

std::string sanitizeLayoutStem(std::string_view name) {
  std::string out;
  out.reserve(name.size());

  auto pushUnderscore = [&]() {
    if (out.empty()) return;
    if (out.back() == '_') return;
    out.push_back('_');
  };

  for (unsigned char uc : name) {
    const char c = (char)uc;
    if (std::isalnum(uc)) {
      out.push_back((char)std::tolower(uc));
      continue;
    }
    if (c == '_' || c == '-') {
      // Avoid leading underscores.
      if (!out.empty() || c == '-') out.push_back(c);
      continue;
    }
    if (std::isspace(uc)) {
      pushUnderscore();
      continue;
    }
    // Skip other characters.
  }

  // Trim leading/trailing underscores.
  while (!out.empty() && out.front() == '_') out.erase(out.begin());
  while (!out.empty() && out.back() == '_') out.pop_back();

  if (out.empty()) out = "profile";
  return out;
}

std::string layoutIniPath(std::string_view name, std::string_view dir) {
  namespace fs = std::filesystem;
  const std::string d = dir.empty() ? defaultUiLayoutsDir() : std::string(dir);
  const std::string stem = sanitizeLayoutStem(name);
  fs::path p(d);
  p /= stem + ".ini";
  return p.generic_string();
}

std::vector<UiLayoutProfile> listLayoutProfiles(const std::string& dir,
                                                const std::string& defaultIniPath) {
  namespace fs = std::filesystem;

  std::vector<UiLayoutProfile> out;
  out.push_back(UiLayoutProfile{"Default", defaultIniPath, true});

  std::error_code ec;
  if (!fs::exists(dir, ec) || !fs::is_directory(dir, ec)) {
    return out;
  }

  for (const auto& ent : fs::directory_iterator(fs::path(dir), ec)) {
    if (ec) break;
    if (!ent.is_regular_file(ec)) continue;
    const fs::path p = ent.path();
    if (p.extension() != ".ini") continue;
    UiLayoutProfile prof;
    prof.name = p.stem().string();
    prof.path = p.generic_string();
    prof.isDefault = false;
    out.push_back(std::move(prof));
  }

  if (out.size() > 1) {
    std::sort(out.begin() + 1, out.end(), [](const UiLayoutProfile& a, const UiLayoutProfile& b) {
      return a.name < b.name;
    });
  }

  return out;
}

bool ensureLayoutDir(const std::string& dir) {
  namespace fs = std::filesystem;
  std::error_code ec;
  if (dir.empty()) return false;
  if (fs::exists(dir, ec)) return true;
  fs::create_directories(dir, ec);
  if (ec) {
    stellar::core::log(stellar::core::LogLevel::Warn, "UiLayoutProfiles: failed to create dir: " + dir);
    return false;
  }
  return true;
}

bool isSafeLayoutIniPath(std::string_view dir, std::string_view path) {
  namespace fs = std::filesystem;
  if (dir.empty() || path.empty()) return false;
  fs::path p(path);
  fs::path d(dir);
  if (p.is_absolute()) return false;

  // Normalize lexical elements (.. etc) without touching filesystem.
  p = p.lexically_normal();
  d = d.lexically_normal();

  if (p.extension() != ".ini") return false;

  // Only allow direct children: <dir>/<file>.ini
  if (p.parent_path() != d) return false;
  const std::string filename = p.filename().string();
  if (filename.empty()) return false;
  if (filename == ".ini") return false;
  return true;
}

std::string uniqueLayoutIniPath(std::string_view desiredStem, const std::string& dir) {
  namespace fs = std::filesystem;
  const std::string baseStem = sanitizeLayoutStem(desiredStem);
  fs::path d(dir.empty() ? defaultUiLayoutsDir() : dir);

  // Try base name first.
  fs::path candidate = d / (baseStem + ".ini");
  std::error_code ec;
  if (!fs::exists(candidate, ec)) return candidate.generic_string();

  // Add suffixes.
  for (int k = 2; k < 1000; ++k) {
    fs::path c = d / (baseStem + "_" + std::to_string(k) + ".ini");
    ec.clear();
    if (!fs::exists(c, ec)) return c.generic_string();
  }

  // Fallback.
  fs::path c = d / (baseStem + "_copy.ini");
  return c.generic_string();
}

bool copyLayoutIniFile(const std::string& srcPath, const std::string& dstPath, bool overwrite) {
  namespace fs = std::filesystem;
  if (srcPath.empty() || dstPath.empty()) return false;
  std::error_code ec;
  fs::copy_options opt = overwrite ? fs::copy_options::overwrite_existing : fs::copy_options::none;
  fs::copy_file(fs::path(srcPath), fs::path(dstPath), opt, ec);
  return !ec;
}

bool renameLayoutIniFile(const std::string& srcPath, const std::string& dstPath, bool overwrite) {
  namespace fs = std::filesystem;
  if (srcPath.empty() || dstPath.empty()) return false;
  std::error_code ec;
  if (!overwrite && fs::exists(dstPath, ec)) return false;
  ec.clear();
  fs::rename(fs::path(srcPath), fs::path(dstPath), ec);
  return !ec;
}

bool deleteLayoutIniFile(const std::string& path) {
  namespace fs = std::filesystem;
  if (path.empty()) return false;
  std::error_code ec;
  fs::remove(fs::path(path), ec);
  return !ec;
}

} // namespace stellar::ui
