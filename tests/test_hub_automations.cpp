#include "test_harness.h"

#include "stellar/sim/HubAutomations.h"

#include <cmath>
#include <string>
#include <vector>

using namespace stellar;

int test_hub_automations() {
  int failures = 0;

  // Build a simple rule set:
  //   IF kind==5 AND tag contains "foo" THEN schedule one action after 0.5s.
  sim::IntegrationHubAutomationRuleState r;
  r.enabled = true;
  r.name = "RuleA";
  r.eventKind = 5;
  r.tagMatch = (core::u8)sim::HubAutomationTagMatch::Contains;
  r.tagText = "foo";
  r.cooldownSec = 1.0;

  sim::IntegrationHubAutomationActionState a;
  a.kind = 10;
  a.u64aSource = (core::u8)sim::HubAutomationValueSource::EventU64a;
  a.u64aConst = 111;
  a.u64bSource = (core::u8)sim::HubAutomationValueSource::EventU64b;
  a.u64bConst = 222;
  a.i32Const = -7;
  a.bConst = true;
  a.delaySec = 0.5;
  a.msgTemplate = "K={kind}|tag={tag}|a={u64a}|b={u64b}|rule={rule}|msg={msg}";
  r.actions.push_back(a);

  std::vector<sim::IntegrationHubAutomationRuleState> rules;
  rules.push_back(r);

  sim::HubAutomationRuntime rt;
  sim::hubAutomationInitRuntime(rt, rules);

  sim::HubAutomationEvent ev;
  ev.kind = 5;
  ev.kindName = "Validation";
  ev.tag = "zzz_foo_zzz";
  ev.msg = "hello";
  ev.u64a = 123;
  ev.u64b = 999;
  ev.tRealSec = 100.0;
  ev.tSimDays = 10.0;

  sim::hubAutomationApplyRules(rt, rules, ev);
  CHECK(rt.scheduled.size() == 1);

  // Not due yet.
  {
    std::vector<sim::HubAutomationAction> due;
    sim::hubAutomationTickScheduled(rt, 100.2, &due);
    CHECK(due.empty());
    CHECK(rt.scheduled.size() == 1);
  }

  // Due.
  {
    std::vector<sim::HubAutomationAction> due;
    sim::hubAutomationTickScheduled(rt, 100.6, &due);
    CHECK(due.size() == 1);
    CHECK(rt.scheduled.empty());
    const auto& out = due[0];
    CHECK(out.kind == 10);
    CHECK(out.u64a == 123);
    CHECK(out.u64b == 999);
    CHECK(out.i32a == -7);
    CHECK(out.b == true);
    CHECK(std::fabs(out.tRealSec - 100.5) < 1e-9);
    CHECK(std::fabs(out.tSimDays - (10.0 + 0.5 / 86400.0)) < 1e-12);

    CHECK(out.msg.find("K=Validation") != std::string::npos);
    CHECK(out.msg.find("tag=zzz_foo_zzz") != std::string::npos);
    CHECK(out.msg.find("a=123") != std::string::npos);
    CHECK(out.msg.find("b=999") != std::string::npos);
    CHECK(out.msg.find("rule=RuleA") != std::string::npos);
    CHECK(out.msg.find("msg=hello") != std::string::npos);
  }

  // Cooldown: firing again within 1.0s should produce no actions.
  {
    ev.tRealSec = 100.75;
    ev.tSimDays = 10.0 + (0.75 / 86400.0);
    sim::hubAutomationApplyRules(rt, rules, ev);
    CHECK(rt.scheduled.empty());
  }

  // After cooldown, should fire again.
  {
    ev.tRealSec = 101.25;
    ev.tSimDays = 10.0 + (1.25 / 86400.0);
    sim::hubAutomationApplyRules(rt, rules, ev);
    CHECK(rt.scheduled.size() == 1);
  }

  // Tag mismatch should not fire.
  {
    rt.scheduled.clear();
    ev.tRealSec = 200.0;
    ev.tag = "no_match";
    sim::hubAutomationApplyRules(rt, rules, ev);
    CHECK(rt.scheduled.empty());
  }

  return failures;
}
