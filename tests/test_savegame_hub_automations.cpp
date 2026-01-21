#include "test_harness.h"

#include "stellar/sim/SaveGame.h"

#include <cmath>
#include <filesystem>
#include <string>

int test_savegame_hub_automations() {
  int failures = 0;

  sim::SaveGame s{};
  s.version = 32;
  s.seed = 123;
  s.hubAutomationsEnabled = false;

  sim::IntegrationHubAutomationRuleState r{};
  r.enabled = true;
  r.name = "Test rule";
  r.eventKind = 5;  // GameEventKind::Validation in apps/stellar_game
  r.tagMatch = 2;   // AutomationTagMatch::Contains
  r.tagText = "foo";
  r.cooldownSec = 1.25;

  sim::IntegrationHubAutomationActionState a{};
  a.kind = 10;      // GameActionKind::TransmitComms
  a.u64aSource = 1; // AutomationValueSource::EventU64a
  a.u64aConst = 42;
  a.u64bSource = 2; // AutomationValueSource::EventU64b
  a.u64bConst = 999;
  a.i32Const = -7;
  a.bConst = true;
  a.delaySec = 0.5;
  a.msgTemplate = "Hello|World";
  r.actions.push_back(a);

  s.hubAutomationRules.push_back(r);

  const std::string path = "__test_savegame_hub_automations.sav";
  CHECK(sim::saveToFile(s, path));

  sim::SaveGame loaded{};
  CHECK(sim::loadFromFile(path, loaded));

  CHECK(loaded.hubAutomationsEnabled == s.hubAutomationsEnabled);
  CHECK(loaded.hubAutomationRules.size() == 1);

  const auto& lr = loaded.hubAutomationRules[0];
  CHECK(lr.enabled == r.enabled);
  CHECK(lr.name == r.name);
  CHECK(lr.eventKind == r.eventKind);
  CHECK(lr.tagMatch == r.tagMatch);
  CHECK(lr.tagText == r.tagText);
  CHECK(std::fabs(lr.cooldownSec - r.cooldownSec) < 1e-9);

  CHECK(lr.actions.size() == 1);
  const auto& la = lr.actions[0];
  CHECK(la.kind == a.kind);
  CHECK(la.u64aSource == a.u64aSource);
  CHECK(la.u64aConst == a.u64aConst);
  CHECK(la.u64bSource == a.u64bSource);
  CHECK(la.u64bConst == a.u64bConst);
  CHECK(la.i32Const == a.i32Const);
  CHECK(la.bConst == a.bConst);
  CHECK(std::fabs(la.delaySec - a.delaySec) < 1e-9);
  CHECK(la.msgTemplate == a.msgTemplate);

  std::filesystem::remove(path);
  return failures;
}
