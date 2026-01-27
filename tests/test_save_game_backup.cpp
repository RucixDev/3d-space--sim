#include "stellar/sim/SaveGame.h"

#include "test_harness.h"

#include <chrono>
#include <filesystem>

namespace {
std::filesystem::path uniqueTempPath(const char* stem) {
  const auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  auto dir = std::filesystem::temp_directory_path() / "stellar_save_game_tests";
  std::error_code ec;
  std::filesystem::create_directories(dir, ec);
  // Even if directory creation fails (permissions, sandbox), we'll still try to use it.
  return dir / (std::string(stem) + "_" + std::to_string(static_cast<long long>(now)) + ".sav");
}
}  // namespace

int test_save_game_backup() {
  int failures = 0;

  const auto savePath = uniqueTempPath("save_backup");
  auto bakPath = savePath;
  bakPath += ".bak";

  {
    std::error_code ec;
    std::filesystem::remove(savePath, ec);
    std::filesystem::remove(bakPath, ec);
  }

  stellar::sim::SaveGame s1;
  s1.seed = 111;
  CHECK(stellar::sim::saveToFile(s1, savePath.string()));

  stellar::sim::SaveGame s2;
  s2.seed = 222;
  CHECK(stellar::sim::saveToFile(s2, savePath.string()));

  CHECK(std::filesystem::exists(savePath));
  CHECK(std::filesystem::exists(bakPath));

  stellar::sim::SaveGame loadedMain;
  stellar::sim::SaveGame loadedBak;
  CHECK(stellar::sim::loadFromFile(savePath.string(), loadedMain));
  CHECK(stellar::sim::loadFromFile(bakPath.string(), loadedBak));

  CHECK(loadedMain.seed == 222);
  CHECK(loadedBak.seed == 111);

  return failures;
}
