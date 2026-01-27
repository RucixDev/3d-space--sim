#include "stellar/core/AtomicWriteFile.h"

#include "test_harness.h"

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

static std::string readAllBytes(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
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

  // Ensure we didn't leave temp files behind.
  bool hasTmp = false;
  for (const auto& entry : fs::directory_iterator(dir)) {
    const std::string name = entry.path().filename().string();
    if (name.find(".tmp.") != std::string::npos) {
      hasTmp = true;
      break;
    }
  }
  CHECK(!hasTmp);

  fs::remove_all(dir, ec);
  return failures;
}
