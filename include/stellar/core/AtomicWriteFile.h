#pragma once

#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <system_error>
#include <type_traits>

namespace stellar::core {

namespace detail {

inline std::string pathToUtf8String(const std::filesystem::path& p) {
#if defined(_WIN32)
  // On Windows, prefer u8string to avoid mojibake when paths include non-ASCII.
  return p.generic_string();
#else
  // On POSIX, string() is typically UTF-8 already.
  return p.string();
#endif
}

inline std::filesystem::path makeTempSiblingPath(const std::filesystem::path& destPath) {
  static std::atomic<std::uint64_t> counter{0};

  const auto seq = counter.fetch_add(1, std::memory_order_relaxed);
  const auto now = static_cast<std::uint64_t>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());

  std::string tmpName = destPath.filename().generic_string();
  tmpName += ".tmp.";
  tmpName += std::to_string(now);
  tmpName += ".";
  tmpName += std::to_string(seq);

  return destPath.parent_path() / tmpName;
}

} // namespace detail

// WriteFn can be either:
//   - bool(std::ostream& out, std::string* outErr)
//   - void(std::ostream& out)
//
// Returns true on success. On failure, returns false and (optionally) fills outErr.
template <class WriteFn>
bool atomicWriteFile(const std::filesystem::path& destPath, WriteFn&& writeFn, std::string* outErr = nullptr) {
  if (destPath.empty()) {
    if (outErr) *outErr = "atomicWriteFile: destination path is empty";
    return false;
  }

  // Ensure parent directory exists.
  {
    const auto parent = destPath.parent_path();
    if (!parent.empty()) {
      std::error_code ec;
      std::filesystem::create_directories(parent, ec);
      if (ec) {
        if (outErr) {
          *outErr = "atomicWriteFile: failed to create parent directories for '" +
                    detail::pathToUtf8String(destPath) + "': " + ec.message();
        }
        return false;
      }
    }
  }

  const auto tmpPath = detail::makeTempSiblingPath(destPath);

  // Write to temp file.
  {
    std::ofstream out(tmpPath, std::ios::binary | std::ios::trunc);
    if (!out.is_open()) {
      if (outErr) {
        *outErr = "atomicWriteFile: failed to open temp file for write: '" +
                  detail::pathToUtf8String(tmpPath) + "'";
      }
      return false;
    }

    std::string writeErr;
    bool ok = true;

    if constexpr (std::is_invocable_r_v<bool, WriteFn, std::ostream&, std::string*>) {
      ok = writeFn(out, &writeErr);
    } else if constexpr (std::is_invocable_v<WriteFn, std::ostream&>) {
      writeFn(out);
      ok = true;
    } else {
      static_assert(std::is_invocable_v<WriteFn, std::ostream&> ||
                        std::is_invocable_r_v<bool, WriteFn, std::ostream&, std::string*>,
                    "atomicWriteFile: WriteFn must be invocable as void(std::ostream&) or "
                    "bool(std::ostream&, std::string*)");
    }

    out.flush();
    out.close();

    if (!ok) {
      std::error_code ec;
      std::filesystem::remove(tmpPath, ec);
      if (outErr) {
        if (!writeErr.empty()) {
          *outErr = writeErr;
        } else {
          *outErr = "atomicWriteFile: write callback returned failure";
        }
      }
      return false;
    }

    if (!writeErr.empty()) {
      std::error_code ec;
      std::filesystem::remove(tmpPath, ec);
      if (outErr) *outErr = writeErr;
      return false;
    }

    if (out.fail()) {
      std::error_code ec;
      std::filesystem::remove(tmpPath, ec);
      if (outErr) {
        *outErr = "atomicWriteFile: stream failure while writing temp file: '" +
                  detail::pathToUtf8String(tmpPath) + "'";
      }
      return false;
    }
  }

  // Rename temp file into place.
  {
    std::error_code ec;
    std::filesystem::rename(tmpPath, destPath, ec);

    // Windows: rename fails if destination exists; attempt best-effort replacement.
    if (ec) {
      std::error_code removeEc;
      std::filesystem::remove(destPath, removeEc);
      ec.clear();
      std::filesystem::rename(tmpPath, destPath, ec);
    }

    if (ec) {
      std::error_code cleanupEc;
      std::filesystem::remove(tmpPath, cleanupEc);
      if (outErr) {
        *outErr = "atomicWriteFile: failed to rename temp file into place ('" +
                  detail::pathToUtf8String(tmpPath) + "' -> '" +
                  detail::pathToUtf8String(destPath) + "'): " + ec.message();
      }
      return false;
    }
  }

  return true;
}

} // namespace stellar::core
