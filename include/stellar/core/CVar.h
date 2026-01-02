#pragma once

#include <cstdint>
#include <functional>
#include <map>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace stellar::core {

// Console variables ("CVars") are small runtime-tweakable settings.
//
// Goals:
//  - Simple, dependency-free (no JSON parser required).
//  - Deterministic iteration order (for stable console output / tests).
//  - Optional persistence via a line-based config file.
//
// This is intentionally lightweight; it's meant for a prototype / sandbox
// and can be expanded later (scopes, replication, change batching, etc.).

enum class CVarType : std::uint8_t {
  Bool   = 0,
  Int    = 1,
  Float  = 2,
  String = 3
};

enum CVarFlags : std::uint32_t {
  CVar_None     = 0u,
  CVar_Archive  = 1u << 0, // include in saveFile()
  CVar_ReadOnly = 1u << 1, // cannot be changed at runtime
  CVar_Cheat    = 1u << 2  // "cheat" var (UI can hide behind a toggle)
};

inline constexpr std::uint32_t operator|(CVarFlags a, CVarFlags b) {
  return (std::uint32_t)a | (std::uint32_t)b;
}

inline constexpr std::uint32_t operator&(std::uint32_t a, CVarFlags b) {
  return a & (std::uint32_t)b;
}

using CVarValue = std::variant<bool, std::int64_t, double, std::string>;
using CVarListener = std::function<void(const struct CVar&)>;

struct CVar {
  std::string name;
  std::string help;
  CVarType type{CVarType::String};
  std::uint32_t flags{CVar_None};

  CVarValue value{};
  CVarValue defaultValue{};

  // Optional change listeners. Called *after* a successful set/reset.
  std::vector<CVarListener> listeners;
};

class CVarRegistry {
public:
  CVarRegistry() = default;

  // ---- Lookup ----
  const CVar* find(std::string_view name) const;
  CVar*       find(std::string_view name);

  bool exists(std::string_view name) const { return find(name) != nullptr; }

  // ---- Definition (idempotent) ----
  // If the variable already exists with the same type, it is returned as-is.
  // If it exists with a different type, define* returns nullptr.
  CVar* defineBool(std::string_view name, bool defaultValue,
                   std::uint32_t flags = CVar_Archive, std::string_view help = {});
  CVar* defineInt(std::string_view name, std::int64_t defaultValue,
                  std::uint32_t flags = CVar_Archive, std::string_view help = {});
  CVar* defineFloat(std::string_view name, double defaultValue,
                    std::uint32_t flags = CVar_Archive, std::string_view help = {});
  CVar* defineString(std::string_view name, std::string defaultValue,
                     std::uint32_t flags = CVar_Archive, std::string_view help = {});

  // ---- Get (typed) ----
  bool        getBool(std::string_view name, bool fallback = false) const;
  std::int64_t getInt(std::string_view name, std::int64_t fallback = 0) const;
  double      getFloat(std::string_view name, double fallback = 0.0) const;
  std::string getString(std::string_view name, std::string_view fallback = {}) const;

  // ---- Set ----
  bool setBool(std::string_view name, bool v, std::string* outError = nullptr);
  bool setInt(std::string_view name, std::int64_t v, std::string* outError = nullptr);
  bool setFloat(std::string_view name, double v, std::string* outError = nullptr);
  bool setString(std::string_view name, std::string v, std::string* outError = nullptr);

  // Parse and set from a user-provided string.
  bool setFromString(std::string_view name, std::string_view value, std::string* outError = nullptr);

  // Reset to the default value.
  bool reset(std::string_view name, std::string* outError = nullptr);

  // ---- Observability ----
  bool addListener(std::string_view name, CVarListener cb, std::string* outError = nullptr);

  // ---- Introspection ----
  // List variables in deterministic (name-sorted) order. If `filter` is non-empty, only
  // vars whose name contains `filter` (case-insensitive) are returned.
  std::vector<const CVar*> list(std::string_view filter = {}) const;

  static const char* typeName(CVarType t);

  // Convert a variable's value to a user-displayable string.
  static std::string valueToString(const CVar& v);

  // ---- Persistence ----
  // Simple config file format:
  //   # comment
  //   some.bool = true
  //   some.int  = 42
  //   some.str  = "hello world"
  //
  // Unknown variables are stored as "pending" assignments, and applied automatically
  // when/if the variable is defined later.
  bool loadFile(const std::string& path, std::string* outError = nullptr);
  bool saveFile(const std::string& path, std::string* outError = nullptr) const;

  // Pending assignment introspection (mostly for tests/debug).
  bool hasPending(std::string_view name) const;
  std::optional<std::string> pendingValue(std::string_view name) const;

private:
  // Internal define helper that may apply pending assignments.
  CVar* defineImpl(std::string_view name, CVarType type, CVarValue def,
                   std::uint32_t flags, std::string_view help);

  bool setValueImpl(std::string_view name, const CVarValue& v, std::string* outError);

  // Called with mutex held.
  void applyPendingLocked(CVar& var);

  mutable std::mutex mutex_;

  // Keep deterministic ordering.
  std::map<std::string, CVar, std::less<>> vars_;

  // Unknown var assignments from loadFile().
  std::map<std::string, std::string, std::less<>> pending_;
};

// Global registry for convenience (especially for the in-game console).
CVarRegistry& cvars();

// Install a few core/default CVars (safe to call multiple times).
void installDefaultCVars();

} // namespace stellar::core
