#include "test_harness.h"

#include "stellar/ui/HudSettings.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>

namespace fs = std::filesystem;

int test_hud_settings_threat_assist() {
  int failures = 0;

  const fs::path p = fs::temp_directory_path() / "stellar_hud_settings_threat_assist.txt";

  {
    std::ofstream f(p);
    CHECK(f.good());

    f << "Version 10\n";
    f << "ThreatAssistEnabled 0\n";
    f << "ThreatAssistHudIndicator 1\n";
    f << "ThreatAssistAutoApply 1\n";
    f << "ThreatAssistStrength 1.50\n";             // clamp -> 1.0
    f << "ThreatAssistMaxThrust01 0.05\n";          // clamp -> 0.1
    f << "ThreatAssistMissileEngageTtiSec 120.0\n"; // clamp -> 60.0
    f << "ThreatAssistPreferLateralJink 0\n";
    f << "ThreatAssistAllowBoost 0\n";
  }

  stellar::ui::HudSettings s = stellar::ui::makeDefaultHudSettings();
  CHECK(stellar::ui::loadFromFile(p.string(), s));
  CHECK_EQ(s.version, 10);

  CHECK(!s.threatAssistEnabled);
  CHECK(s.threatAssistHudIndicator);
  CHECK(s.threatAssistAutoApply);
  CHECK(std::abs(s.threatAssistStrength - 1.0) < 1e-9);
  CHECK(std::abs(s.threatAssistMaxThrust01 - 0.1) < 1e-9);
  CHECK(std::abs(s.threatAssistMissileEngageTtiSec - 60.0) < 1e-9);
  CHECK(!s.threatAssistPreferLateralJink);
  CHECK(!s.threatAssistAllowBoost);

  std::error_code ec;
  fs::remove(p, ec);
  CHECK(!ec);

  return failures;
}
