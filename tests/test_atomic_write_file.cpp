#include "stellar/core/AtomicWriteFile.h"

#include "test_harness.h"

#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>

static std::string readAllBytes(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

static bool hasTmpFiles(const std::filesystem::path& dir) {
  namespace fs = std::filesystem;

  // Avoid lossy/pathological wide<->narrow conversions on Windows when test files
  // include non-ASCII names (use native string type instead of .string()).
  const auto needle = fs::path(".tmp.").native();

  for (const auto& entry : fs::recursive_directory_iterator(dir)) {
    const auto name = entry.path().filename().native();
    if (name.find(needle) != decltype(name)::npos) return true;
  }
  return false;
}

int test_atomic_write_file() {
  int failures = 0;

  namespace fs = std::filesystem;

  std::error_code ec;
  const fs::path dir = fs::temp_directory_path() / "stellar_test_atomic_write";

  fs::remove_all(dir, ec);
  ec.clear();
  fs::create_directories(dir, ec);
  CHECK(!ec);

  const fs::path outPath = dir / "hello.txt";

  std::string err;
  const bool ok1 = stellar::core::atomicWriteFile(
      outPath,
      [&](std::ostream& out, std::string* outErr) -> bool {
        (void)outErr;
        out << "hello";
        return true;
      },
      &err);

  CHECK(ok1);
  CHECK(err.empty());
  CHECK(readAllBytes(outPath) == "hello");

  err.clear();
  const bool ok2 = stellar::core::atomicWriteFile(
      outPath,
      [&](std::ostream& out) {
        out << "world";
      },
      &err);

  CHECK(ok2);
  CHECK(err.empty());
  CHECK(readAllBytes(outPath) == "world");

  // ---------------------------------------------------------------------------
  // Failure path: callback returns false -> destination remains unchanged.
  // ---------------------------------------------------------------------------
  {
    err.clear();
    const bool okFail = stellar::core::atomicWriteFile(
        outPath,
        [&](std::ostream& out, std::string* outErr) -> bool {
          (void)outErr;
          out << "NEW_CONTENT_SHOULD_NOT_LAND";
          return false;
        },
        &err);

    CHECK(!okFail);
    CHECK(!err.empty());
    CHECK(readAllBytes(outPath) == "world");
    CHECK(!hasTmpFiles(dir));
  }

  // ---------------------------------------------------------------------------
  // Failure path: callback throws -> atomicWriteFile catches and cleans up.
  // ---------------------------------------------------------------------------
  {
    err.clear();
    const bool okThrow = stellar::core::atomicWriteFile(
        outPath,
        [&](std::ostream&) {
          throw std::runtime_error("boom");
        },
        &err);

    CHECK(!okThrow);
    CHECK(!err.empty());
    CHECK(readAllBytes(outPath) == "world");
    CHECK(!hasTmpFiles(dir));
  }

  // ---------------------------------------------------------------------------
  // Parent directory creation: atomicWriteFile should create missing folders.
  // ---------------------------------------------------------------------------
  {
    const fs::path deepPath = dir / "nested" / "deeper" / "ok.txt";
    err.clear();
    const bool okDeep = stellar::core::atomicWriteFile(
        deepPath,
        [&](std::ostream& out) {
          out << "deep";
        },
        &err);
    CHECK(okDeep);
    CHECK(err.empty());
    CHECK(readAllBytes(deepPath) == "deep");
  }


  // ---------------------------------------------------------------------------
  // Unicode filename: ensure atomicWriteFile works for non-ASCII paths.
  // ---------------------------------------------------------------------------
  {
#if defined(_WIN32)
    const fs::path uniPath = dir / fs::path(L"\u30e6\u30cb\u30b3\u30fc\u30c9.txt");
#else
    // "ユニコード.txt" as explicit UTF-8 bytes (avoids source-encoding surprises).
    const fs::path uniPath = dir / "\xE3\x83\xA6\xE3\x83\x8B\xE3\x82\xB3\xE3\x83\xBC\xE3\x83\x89.txt";
#endif

    err.clear();
    const bool okUni = stellar::core::atomicWriteFile(
        uniPath,
        [&](std::ostream& out) {
          out << "unicode";
        },
        &err);

    CHECK(okUni);
    CHECK(err.empty());
    CHECK(readAllBytes(uniPath) == "unicode");
    CHECK(!hasTmpFiles(dir));
  }

  // Ensure we didn't leave temp files behind.
  CHECK(!hasTmpFiles(dir));

  fs::remove_all(dir, ec);
  return failures;
}
