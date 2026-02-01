#pragma once

#include "stellar/core/PathUtf8.h"

#include <atomic>
#include <chrono>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>

#if defined(_WIN32)
  // Keep Windows headers as small as possible.
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  // Avoid <windows.h> defining min/max macros.
  #ifndef NOMINMAX
    #define NOMINMAX
  #endif
  #include <windows.h>

  // The Windows SDK can still define legacy `near`/`far` macros (from 16-bit era)
  // which can break otherwise valid code that uses these identifiers.
  // Undefine them here to keep this header from polluting downstream translation units.
  #ifdef near
    #undef near
  #endif
  #ifdef far
    #undef far
  #endif
#endif

namespace stellar::core {

namespace detail {

inline std::filesystem::path makeTempSiblingPath(const std::filesystem::path& destPath) {
  static std::atomic<std::uint64_t> counter{0};

  const auto seq = counter.fetch_add(1, std::memory_order_relaxed);
  const auto now = static_cast<std::uint64_t>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());

  // Build the temp filename using the path's native string type to avoid lossy
  // wide<->narrow conversions on Windows when destPath contains Unicode.
  using NativeStr = std::filesystem::path::string_type;
  NativeStr tmpName = destPath.filename().native();

  auto appendAscii = [](NativeStr& out, std::string_view ascii) {
    out.reserve(out.size() + ascii.size());
    for (unsigned char c : ascii) {
      out.push_back((NativeStr::value_type)c);
    }
  };

  appendAscii(tmpName, ".tmp.");
  appendAscii(tmpName, std::to_string(now));
  appendAscii(tmpName, ".");
  appendAscii(tmpName, std::to_string(seq));

  return destPath.parent_path() / std::filesystem::path(tmpName);
}

#if defined(_WIN32)
inline bool replaceFileWin32(const std::filesystem::path& srcTmp,
                             const std::filesystem::path& dest,
                             std::string* outErr) {
  // Use MoveFileExW with MOVEFILE_REPLACE_EXISTING to replace the destination
  // atomically. MOVEFILE_WRITE_THROUGH asks the OS to flush the operation
  // before returning (best-effort).
  //
  // Docs: MOVEFILE_REPLACE_EXISTING + MOVEFILE_WRITE_THROUGH.
  // https://learn.microsoft.com/windows/win32/api/winbase/nf-winbase-movefileexw
  const DWORD flags = MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH;

  if (MoveFileExW(srcTmp.c_str(), dest.c_str(), flags) != 0) {
    return true;
  }

  const DWORD e = GetLastError();
  std::error_code ec((int)e, std::system_category());
  if (outErr) {
    *outErr = "atomicWriteFile: MoveFileExW failed ('" + pathToUtf8String(srcTmp) +
              "' -> '" + pathToUtf8String(dest) + "'): " + ec.message();
  }
  return false;
}
#endif

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
                    pathToUtf8String(destPath) + "': " + ec.message();
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
                  pathToUtf8String(tmpPath) + "'";
      }
      return false;
    }

    std::string writeErr;
    bool ok = true;

    // Keep atomicWriteFile exception-safe: if the callback throws, we must
    // cleanup the temp file and report failure.
    try {
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
    } catch (const std::exception& e) {
      ok = false;
      writeErr = std::string("atomicWriteFile: exception in write callback: ") + e.what();
    } catch (...) {
      ok = false;
      writeErr = "atomicWriteFile: unknown exception in write callback";
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
                  pathToUtf8String(tmpPath) + "'";
      }
      return false;
    }
  }

  // Rename temp file into place.
  {
  #if defined(_WIN32)
    std::string renameErr;
    if (!detail::replaceFileWin32(tmpPath, destPath, &renameErr)) {
      std::error_code cleanupEc;
      std::filesystem::remove(tmpPath, cleanupEc);
      if (outErr) *outErr = renameErr;
      return false;
    }
  #else
    std::error_code ec;
    std::filesystem::rename(tmpPath, destPath, ec);
    if (ec) {
      std::error_code cleanupEc;
      std::filesystem::remove(tmpPath, cleanupEc);
      if (outErr) {
        *outErr = "atomicWriteFile: failed to rename temp file into place ('" +
                  pathToUtf8String(tmpPath) + "' -> '" +
                  pathToUtf8String(destPath) + "'): " + ec.message();
      }
      return false;
    }
  #endif
  }

  return true;
}

} // namespace stellar::core
