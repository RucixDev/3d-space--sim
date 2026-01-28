#include "stellar/sim/HubAutomations.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <sstream>

namespace stellar::sim {

static std::string sanitizeForFilename(std::string s) {
  // Keep it conservative: avoid platform-specific illegal chars.
  for (char& c : s) {
    const unsigned char uc = (unsigned char)c;
    if (std::isalnum(uc) || c == '_' || c == '-' || c == '.') continue;
    c = '_';
  }
  return s;
}

bool hubAutomationTagMatches(core::u8 matchKind, std::string_view pattern, std::string_view tag) {
  const auto mk = static_cast<HubAutomationTagMatch>(matchKind);
  switch (mk) {
    case HubAutomationTagMatch::Any:
      return true;
    case HubAutomationTagMatch::Equals:
      return (!pattern.empty()) && (tag == pattern);
    case HubAutomationTagMatch::Contains:
      return (!pattern.empty()) && (tag.find(pattern) != std::string_view::npos);
    case HubAutomationTagMatch::Prefix:
      return (!pattern.empty()) && (tag.size() >= pattern.size()) &&
             (tag.substr(0, pattern.size()) == pattern);
    case HubAutomationTagMatch::Suffix:
      return (!pattern.empty()) && (tag.size() >= pattern.size()) &&
             (tag.substr(tag.size() - pattern.size()) == pattern);
    default:
      return false;
  }
}

core::u64 hubAutomationResolveU64(core::u8 sourceKind, core::u64 constant, const HubAutomationEvent& e) {
  const auto src = static_cast<HubAutomationValueSource>(sourceKind);
  switch (src) {
    case HubAutomationValueSource::Constant: return constant;
    case HubAutomationValueSource::EventU64a: return e.u64a;
    case HubAutomationValueSource::EventU64b: return e.u64b;
    default: return constant;
  }
}

std::string hubAutomationExpandTemplate(std::string_view tmpl, const HubAutomationEvent& e, std::string_view ruleName) {
  // Small, dependency-free template replacement:
  //   {tag}, {tag_safe}, {kind}, {msg}, {u64a}, {u64b}, {tRealMs}, {tRealSec}, {tSimDays}, {rule}
  //
  // If HubAutomationEvent::kindName is set, {kind} uses it. Otherwise it falls
  // back to the numeric value.
  std::string out;
  out.reserve(tmpl.size() + 32);

  auto appendU64 = [&](core::u64 v) {
    out += std::to_string((unsigned long long)v);
  };

  for (std::size_t i = 0; i < tmpl.size();) {
    const char c = tmpl[i];
    if (c != '{') {
      out.push_back(c);
      ++i;
      continue;
    }

    const std::size_t j = tmpl.find('}', i + 1);
    if (j == std::string_view::npos) {
      out.append(tmpl.substr(i));
      break;
    }

    const std::string_view key = tmpl.substr(i + 1, j - (i + 1));

    if (key == "tag") {
      out.append(e.tag);
    } else if (key == "tag_safe") {
      out.append(sanitizeForFilename(e.tag));
    } else if (key == "kind") {
      if (!e.kindName.empty()) {
        out.append(e.kindName);
      } else {
        appendU64((core::u64)e.kind);
      }
    } else if (key == "msg") {
      out.append(e.msg);
    } else if (key == "msg_safe") {
      out.append(sanitizeForFilename(e.msg));
    } else if (key == "u64a") {
      appendU64(e.u64a);
    } else if (key == "u64b") {
      appendU64(e.u64b);
    } else if (key == "tRealMs") {
      const long long ms = (long long)std::llround(e.tRealSec * 1000.0);
      out += std::to_string(ms);
    } else if (key == "tRealSec") {
      std::ostringstream ss;
      ss.setf(std::ios::fixed);
      ss.precision(3);
      ss << e.tRealSec;
      out += ss.str();
    } else if (key == "tSimDays") {
      std::ostringstream ss;
      ss.setf(std::ios::fixed);
      ss.precision(6);
      ss << e.tSimDays;
      out += ss.str();
    } else if (key == "rule") {
      out.append(ruleName);
    } else {
      // Unknown key: keep literal.
      out.push_back('{');
      out.append(key);
      out.push_back('}');
    }

    i = j + 1;
  }

  return out;
}

void hubAutomationInitRuntime(HubAutomationRuntime& rt, const std::vector<IntegrationHubAutomationRuleState>& rules) {
  // Preserve existing last-fired times when possible so callers can rebuild the
  // rules vector without resetting cooldown state every frame.
  if (rt.lastFiredRealSec.size() == rules.size()) return;

  std::vector<double> next;
  next.resize(rules.size(), -1e9);

  const std::size_t n = std::min(rt.lastFiredRealSec.size(), next.size());
  for (std::size_t i = 0; i < n; ++i) {
    const double v = rt.lastFiredRealSec[i];
    next[i] = std::isfinite(v) ? v : -1e9;
  }
  rt.lastFiredRealSec = std::move(next);
}

static void scheduleSorted(HubAutomationRuntime& rt, HubAutomationAction a) {
  // Insert so rt.scheduled remains sorted by tRealSec.
  const double t = a.tRealSec;
  auto it = std::upper_bound(rt.scheduled.begin(), rt.scheduled.end(), t,
                             [](double lhs, const HubAutomationAction& rhs) {
                               return lhs < rhs.tRealSec;
                             });
  rt.scheduled.insert(it, std::move(a));

  // Keep earliest actions; drop far-future overflow.
  while ((int)rt.scheduled.size() > rt.maxScheduledActions) {
    rt.scheduled.pop_back();
  }
}

void hubAutomationApplyRules(
    HubAutomationRuntime& rt,
    const std::vector<IntegrationHubAutomationRuleState>& rules,
    const HubAutomationEvent& ev) {
  if (!rt.enabled) return;
  if (rules.empty()) return;

  hubAutomationInitRuntime(rt, rules);

  // Hard recursion cap to prevent automation-driven event loops from causing
  // unbounded action scheduling.
  if (rt.recursionDepth > 3) return;
  struct DepthGuard {
    HubAutomationRuntime& rt;
    explicit DepthGuard(HubAutomationRuntime& r) : rt(r) { rt.recursionDepth++; }
    ~DepthGuard() { rt.recursionDepth--; }
  } guard(rt);

  int actionsFired = 0;

  for (std::size_t i = 0; i < rules.size(); ++i) {
    const auto& r = rules[i];
    if (!r.enabled) continue;
    if (r.eventKind != ev.kind) continue;
    if (!hubAutomationTagMatches(r.tagMatch, r.tagText, ev.tag)) continue;

    const double lastFired = (i < rt.lastFiredRealSec.size()) ? rt.lastFiredRealSec[i] : -1e9;
    const double dt = ev.tRealSec - lastFired;
    if (dt < r.cooldownSec) continue;

    if (i < rt.lastFiredRealSec.size()) {
      rt.lastFiredRealSec[i] = ev.tRealSec;
    }

    const std::string origin = std::string("Automation:") + r.name;

    for (const auto& t : r.actions) {
      if (actionsFired >= rt.maxActionsPerEvent) return;

      const double delaySec = std::max(0.0, t.delaySec);

      HubAutomationAction a;
      a.tRealSec = ev.tRealSec + delaySec;
      a.tSimDays = ev.tSimDays + (delaySec / 86400.0);
      a.origin = origin;

      a.kind = t.kind;
      a.u64a = hubAutomationResolveU64(t.u64aSource, t.u64aConst, ev);
      a.u64b = hubAutomationResolveU64(t.u64bSource, t.u64bConst, ev);
      a.i32a = t.i32Const;
      a.b = t.bConst;
      a.msg = hubAutomationExpandTemplate(t.msgTemplate, ev, r.name);

      scheduleSorted(rt, std::move(a));
      actionsFired++;
    }
  }
}

void hubAutomationTickScheduled(
    HubAutomationRuntime& rt,
    double nowRealSec,
    std::vector<HubAutomationAction>* outDue) {
  // Move due actions to outDue.
  while (!rt.scheduled.empty()) {
    const HubAutomationAction& a = rt.scheduled.front();
    if (a.tRealSec > nowRealSec + 1e-9) break;
    if (outDue) {
      outDue->push_back(a);
    }
    rt.scheduled.pop_front();
  }
}

} // namespace stellar::sim
