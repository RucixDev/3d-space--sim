#pragma once

#include <filesystem>
#include <string>
#include <type_traits>

namespace stellar::core {

// Convert a filesystem path to a UTF-8 string suitable for logs / JSON.
//
// Implementation notes:
// - Uses `generic_u8string()` to get a stable, platform-independent representation
//   (forward-slash separators).
// - In C++20, `u8string()` and `generic_u8string()` return `std::u8string`
//   (char8_t). We explicitly byte-copy it into a std::string.
inline std::string pathToUtf8String(const std::filesystem::path& p) {
  const auto u8 = p.generic_u8string();
#if defined(__cpp_char8_t)
  using CharT = typename std::decay_t<decltype(u8)>::value_type;
  if constexpr (std::is_same_v<CharT, char8_t>) {
    return std::string(reinterpret_cast<const char*>(u8.data()), u8.size());
  } else {
    return std::string(u8.begin(), u8.end());
  }
#else
  return std::string(u8.begin(), u8.end());
#endif
}

} // namespace stellar::core
