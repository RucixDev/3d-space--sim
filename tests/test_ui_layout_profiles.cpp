#include "test_harness.h"

#include "stellar/ui/UiLayoutProfiles.h"

#include <filesystem>
#include <fstream>

using stellar::ui::UiLayoutProfile;

int test_ui_layout_profiles() {
  int failures = 0;

  namespace fs = std::filesystem;

  // sanitizeLayoutStem
  {
    CHECK(stellar::ui::sanitizeLayoutStem(" Trading 2 ") == "trading_2");
    CHECK(stellar::ui::sanitizeLayoutStem("Combat-Layout") == "combat-layout");
    CHECK(stellar::ui::sanitizeLayoutStem("!!!") == "profile");
  }

  const fs::path dir = fs::path("test_ui_layout_profiles_tmp");
  std::error_code ec;
  fs::remove_all(dir, ec);
  fs::create_directories(dir, ec);
  CHECK(!ec);

  // Create some dummy ini files.
  {
    std::ofstream(dir / "combat.ini").put('\n');
    std::ofstream(dir / "trading.ini").put('\n');
    std::ofstream(dir / "notes.txt").put('\n');
  }

  // listLayoutProfiles
  {
    const auto profs = stellar::ui::listLayoutProfiles(dir.string(), "imgui.ini");
    CHECK(!profs.empty());
    CHECK(profs.front().isDefault);
    CHECK(profs.front().name == "Default");
    CHECK(profs.front().path == "imgui.ini");
    bool sawCombat = false;
    bool sawTrading = false;
    for (const auto& p : profs) {
      if (p.path == (dir / "combat.ini").generic_string()) sawCombat = true;
      if (p.path == (dir / "trading.ini").generic_string()) sawTrading = true;
      // Should never include non-ini.
      CHECK(p.path.find("notes.txt") == std::string::npos);
    }
    CHECK(sawCombat);
    CHECK(sawTrading);
  }

  // isSafeLayoutIniPath
  {
    CHECK(stellar::ui::isSafeLayoutIniPath("ui_layouts", "ui_layouts/combat.ini"));
    CHECK(!stellar::ui::isSafeLayoutIniPath("ui_layouts", "../combat.ini"));
    CHECK(!stellar::ui::isSafeLayoutIniPath("ui_layouts", "ui_layouts/../combat.ini"));
    CHECK(!stellar::ui::isSafeLayoutIniPath("ui_layouts", "ui_layouts/sub/combat.ini"));
    CHECK(!stellar::ui::isSafeLayoutIniPath("ui_layouts", "ui_layouts/combat.txt"));
  }

  // uniqueLayoutIniPath
  {
    const std::string first = stellar::ui::uniqueLayoutIniPath("combat", dir.string());
    // combat.ini exists; should give a suffix.
    CHECK(first.find("combat_2.ini") != std::string::npos);
  }

  fs::remove_all(dir, ec);
  return failures;
}
