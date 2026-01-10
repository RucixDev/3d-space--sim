#include "LogbookWindow.h"

#include "stellar/econ/Commodity.h"
#include "stellar/sim/Signals.h"

#include <imgui.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>
#include <unordered_map>

namespace stellar::game {

namespace {

static bool containsCaseInsensitive(std::string_view haystack, std::string_view needle) {
  if (needle.empty()) return true;
  if (haystack.empty()) return false;

  // Naive search is fine at our expected sizes.
  const auto lower = [](unsigned char c) { return (unsigned char)std::tolower((int)c); };

  for (std::size_t i = 0; i + needle.size() <= haystack.size(); ++i) {
    bool ok = true;
    for (std::size_t j = 0; j < needle.size(); ++j) {
      if (lower((unsigned char)haystack[i + j]) != lower((unsigned char)needle[j])) {
        ok = false;
        break;
      }
    }
    if (ok) return true;
  }
  return false;
}

static std::string systemName(sim::Universe& u, sim::SystemId id) {
  if (id == 0) return std::string("(unknown)");
  const auto& sys = u.getSystem(id);
  return sys.stub.name;
}

static std::string stationNameInSystem(const sim::StarSystem& sys, sim::StationId stId) {
  if (stId == 0) return std::string();
  for (const auto& st : sys.stations) {
    if (st.id == stId) return st.name;
  }
  return std::string("Station#") + std::to_string((unsigned long long)stId);
}

static std::string planetNameInSystem(const sim::StarSystem& sys, std::size_t planetIndex) {
  if (planetIndex < sys.planets.size()) return sys.planets[planetIndex].name;
  return std::string("Planet#") + std::to_string((unsigned long long)planetIndex);
}

static std::string entryTargetLabel(sim::Universe& u, const sim::LogbookEntry& e) {
  switch (e.kind) {
    case sim::LogbookEntryKind::StarScan: {
      return "Primary Star";
    }
    case sim::LogbookEntryKind::PlanetScan: {
      const auto& sys = u.getSystem(e.systemId);
      return planetNameInSystem(sys, (std::size_t)e.objectId);
    }
    case sim::LogbookEntryKind::StationScan: {
      const auto& sys = u.getSystem(e.systemId);
      return stationNameInSystem(sys, e.stationId);
    }
    case sim::LogbookEntryKind::SignalScan: {
      const auto k = static_cast<sim::SignalKind>(e.subKind);
      std::string s = sim::signalKindName(k);
      s += " #";
      s += std::to_string((unsigned long long)e.objectId);
      return s;
    }
    case sim::LogbookEntryKind::AsteroidProspect: {
      const auto& cd = econ::commodityDef(e.commodity);
      std::string s = cd.name;
      if (e.units > 1e-6) {
        s += " (remaining ~";
        s += std::to_string((int)std::llround(e.units));
        s += "u)";
      }
      return s;
    }
    case sim::LogbookEntryKind::SystemSurveyBonus: {
      return "System survey complete";
    }
    default: return std::string();
  }
}

static double computeBrokerMultiplier(const ExplorationDataBrokerState& st,
                                      const sim::StarSystem& saleSystem,
                                      const sim::Station& saleStation,
                                      const sim::SystemStub& scanStub) {
  if (!st.enablePremium) return 1.0;

  double mult = 1.0;

  // Distance premium: linear ramp up to maxDistancePremium.
  const double distLy = (scanStub.posLy - saleSystem.stub.posLy).length();
  const double scale = std::max(1.0, (double)st.distanceScaleLy);
  const double t = std::clamp(distLy / scale, 0.0, 1.0);
  mult *= (1.0 + (double)st.maxDistancePremium * t);

  // Jurisdiction premium/penalty: faction alignment.
  const core::u32 saleF = saleStation.factionId;
  const core::u32 scanF = scanStub.factionId;
  if (saleF != 0 && scanF != 0) {
    if (saleF == scanF) mult *= (1.0 + (double)st.sameFactionBonus);
    else mult *= (1.0 - (double)st.otherFactionPenalty);
  }

  // Clamp to avoid extreme exploitation.
  mult = std::clamp(mult, 0.75, 1.50);
  return mult;
}

} // namespace

void drawExplorationDataBrokerPanel(ExplorationDataBrokerState& st,
                                    sim::Universe& universe,
                                    const sim::StarSystem& saleSystem,
                                    const sim::Station& saleStation,
                                    double /*timeDays*/,
                                    bool canSellHere,
                                    std::vector<sim::LogbookEntry>& logbook,
                                    double& explorationDataBankCr,
                                    double& creditsCr,
                                    std::function<void(std::string_view, double)> toast,
                                    std::function<void(core::u32, double)> addRep) {
  const double unsoldFromLog = sim::logbookUnsoldValueCr(logbook);

  ImGui::Text("Bank (base): %.0f cr", explorationDataBankCr);
  ImGui::TextDisabled("Unsold in log: %.0f cr", unsoldFromLog);

  const double diff = explorationDataBankCr - unsoldFromLog;
  if (std::fabs(diff) > 0.5) {
    ImGui::TextColored(ImVec4(1, 0.75f, 0.2f, 1),
                       "Note: bank/log mismatch (%.0f cr). Older saves won't have per-scan history.", diff);
  }

  if (!canSellHere) {
    ImGui::TextDisabled("Dock at a station to sell exploration data.");
    return;
  }

  struct Group {
    sim::SystemId systemId{0};
    core::u32 factionId{0};
    double distanceLy{0.0};
    double baseCr{0.0};
    double payoutCr{0.0};
    int entries{0};
    std::string name;
    double mult{1.0};
  };

  std::unordered_map<sim::SystemId, std::size_t> idx;
  std::vector<Group> groups;
  groups.reserve(64);

  for (const auto& e : logbook) {
    if (e.sold) continue;
    if (e.valueCr <= 0.0) continue;
    if (e.systemId == 0) continue;

    auto it = idx.find(e.systemId);
    if (it == idx.end()) {
      const std::size_t newIndex = groups.size();
      idx[e.systemId] = newIndex;
      groups.push_back(Group{});
      groups.back().systemId = e.systemId;
      it = idx.find(e.systemId);
    }
    Group& g = groups[it->second];
    g.baseCr += e.valueCr;
    g.entries += 1;
  }

  // Enrich groups with names/distances/multipliers.
  for (auto& g : groups) {
    const auto& sys = universe.getSystem(g.systemId);
    g.name = sys.stub.name;
    g.factionId = sys.stub.factionId;
    g.distanceLy = (sys.stub.posLy - saleSystem.stub.posLy).length();
    g.mult = computeBrokerMultiplier(st, saleSystem, saleStation, sys.stub);
    g.payoutCr = g.baseCr * g.mult;
  }

  // Sort by base value (descending).
  std::sort(groups.begin(), groups.end(), [&](const Group& a, const Group& b) {
    if (a.baseCr != b.baseCr) return a.baseCr > b.baseCr;
    return a.name < b.name;
  });

  // Premium settings.
  ImGui::Separator();
  if (ImGui::TreeNodeEx("Broker Premium Model", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Checkbox("Enable premium payout", &st.enablePremium);
    ImGui::SameLine();
    ImGui::Checkbox("Show details", &st.showDetails);

    ImGui::BeginDisabled(!st.enablePremium);
    ImGui::SliderFloat("Distance scale (ly)", &st.distanceScaleLy, 50.0f, 800.0f, "%.0f");
    ImGui::SliderFloat("Max distance premium", &st.maxDistancePremium, 0.0f, 0.75f, "+%.0f%%");
    ImGui::SliderFloat("Same-faction bonus", &st.sameFactionBonus, 0.0f, 0.50f, "+%.0f%%");
    ImGui::SliderFloat("Other-faction penalty", &st.otherFactionPenalty, 0.0f, 0.50f, "-%.0f%%");
    ImGui::EndDisabled();
    ImGui::TextDisabled("Payout = base * premium (clamped). Bank is tracked in base credits for compatibility.");
    ImGui::TreePop();
  }

  ImGui::Separator();

  if (groups.empty()) {
    ImGui::TextDisabled("No unsold scan data in the logbook.");
    return;
  }

  // Selection helpers.
  if (ImGui::SmallButton("Select all")) {
    st.selectedSystems.clear();
    for (const auto& g : groups) st.selectedSystems.insert(g.systemId);
  }
  ImGui::SameLine();
  if (ImGui::SmallButton("Clear selection")) {
    st.selectedSystems.clear();
  }

  double baseSel = 0.0;
  double payoutSel = 0.0;
  int entriesSel = 0;
  for (const auto& g : groups) {
    if (st.selectedSystems.find(g.systemId) == st.selectedSystems.end()) continue;
    baseSel += g.baseCr;
    payoutSel += g.payoutCr;
    entriesSel += g.entries;
  }

  const ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp;
  if (ImGui::BeginTable("##broker_table", st.showDetails ? 7 : 6, flags)) {
    ImGui::TableSetupColumn("Sell", ImGuiTableColumnFlags_WidthFixed, 38.0f);
    ImGui::TableSetupColumn("System");
    ImGui::TableSetupColumn("Dist (ly)", ImGuiTableColumnFlags_WidthFixed, 72.0f);
    ImGui::TableSetupColumn("Entries", ImGuiTableColumnFlags_WidthFixed, 58.0f);
    ImGui::TableSetupColumn("Base", ImGuiTableColumnFlags_WidthFixed, 90.0f);
    ImGui::TableSetupColumn("Payout", ImGuiTableColumnFlags_WidthFixed, 100.0f);
    if (st.showDetails) ImGui::TableSetupColumn("Mult", ImGuiTableColumnFlags_WidthFixed, 60.0f);
    ImGui::TableHeadersRow();

    for (const auto& g : groups) {
      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);

      bool sel = (st.selectedSystems.find(g.systemId) != st.selectedSystems.end());
      const std::string cbId = "##sell_sys_" + std::to_string((unsigned long long)g.systemId);
      if (ImGui::Checkbox(cbId.c_str(), &sel)) {
        if (sel) st.selectedSystems.insert(g.systemId);
        else st.selectedSystems.erase(g.systemId);
      }

      ImGui::TableSetColumnIndex(1);
      ImGui::TextUnformatted(g.name.c_str());

      ImGui::TableSetColumnIndex(2);
      ImGui::Text("%.0f", g.distanceLy);

      ImGui::TableSetColumnIndex(3);
      ImGui::Text("%d", g.entries);

      ImGui::TableSetColumnIndex(4);
      ImGui::Text("%.0f", g.baseCr);

      ImGui::TableSetColumnIndex(5);
      ImGui::Text("%.0f", g.payoutCr);

      if (st.showDetails) {
        ImGui::TableSetColumnIndex(6);
        ImGui::Text("x%.2f", g.mult);
      }
    }

    ImGui::EndTable();
  }

  ImGui::Separator();
  ImGui::Text("Selected: %d entries | base %.0f cr | payout %.0f cr", entriesSel, baseSel, payoutSel);

  auto doSell = [&](bool sellAll) {
    double baseSold = 0.0;
    double payoutSold = 0.0;
    int soldEntries = 0;

    // Map systemId -> multiplier to avoid repeated universe lookups.
    std::unordered_map<sim::SystemId, double> multCache;
    multCache.reserve(groups.size());
    for (const auto& g : groups) {
      multCache[g.systemId] = g.mult;
    }

    for (auto& e : logbook) {
      if (e.sold) continue;
      if (e.valueCr <= 0.0) continue;
      if (!sellAll) {
        if (st.selectedSystems.find(e.systemId) == st.selectedSystems.end()) continue;
      }
      const auto itM = multCache.find(e.systemId);
      const double mult = (itM != multCache.end()) ? itM->second : 1.0;
      baseSold += e.valueCr;
      payoutSold += e.valueCr * mult;
      e.sold = true;
      ++soldEntries;
    }

    if (baseSold <= 0.0 || soldEntries == 0) return;

    creditsCr += payoutSold;
    explorationDataBankCr = std::max(0.0, explorationDataBankCr - baseSold);
    sim::pruneLogbook(logbook);

    // Rep scales with the amount sold, capped.
    if (saleStation.factionId != 0) {
      const double repGain = std::clamp(baseSold / 5000.0, 0.05, 0.50);
      addRep(saleStation.factionId, repGain);
    }

    st.selectedSystems.clear();

    if (toast) {
      const int payoutI = (int)std::llround(payoutSold);
      const int baseI = (int)std::llround(baseSold);
      std::string msg = "Sold exploration data: ";
      msg += std::to_string(soldEntries);
      msg += " entries | base +";
      msg += std::to_string(baseI);
      msg += " cr | payout +";
      msg += std::to_string(payoutI);
      msg += " cr.";
      toast(msg, 3.0);
    }
  };

  ImGui::BeginDisabled(baseSel <= 0.0);
  if (ImGui::Button("Sell selected")) {
    doSell(/*sellAll=*/false);
  }
  ImGui::EndDisabled();
  ImGui::SameLine();

  const double baseAll = unsoldFromLog;
  ImGui::BeginDisabled(baseAll <= 0.0);
  if (ImGui::Button("Sell all")) {
    st.selectedSystems.clear();
    doSell(/*sellAll=*/true);
  }
  ImGui::EndDisabled();
}

void drawLogbookWindow(LogbookWindowState& st, const LogbookContext& ctx) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(780.0f, 520.0f), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Logbook", &st.open)) {
    ImGui::End();
    return;
  }

  const double unsold = sim::logbookUnsoldValueCr(ctx.logbook);
  ImGui::Text("Exploration data bank: %.0f cr", ctx.explorationDataBankCr);
  ImGui::SameLine();
  ImGui::TextDisabled("(unsold in log: %.0f cr)", unsold);

  const double diff = ctx.explorationDataBankCr - unsold;
  if (std::fabs(diff) > 0.5) {
    ImGui::TextColored(ImVec4(1, 0.75f, 0.2f, 1),
                       "Bank/log mismatch (%.0f cr). Older saves won't have per-scan history.", diff);
  }

  // Controls
  ImGui::Separator();
  ImGui::Checkbox("Show unsold", &st.showUnsold);
  ImGui::SameLine();
  ImGui::Checkbox("Show sold", &st.showSold);
  ImGui::SameLine();

  ImGui::SetNextItemWidth(220.0f);
  ImGui::InputTextWithHint("##log_filter", "Filter (system, station, commodity)...", st.filter, sizeof(st.filter));
  ImGui::SameLine();

  // Kind filter combo
  {
    const char* items[] = {"All", "Star Scan", "Planet Scan", "Station Scan", "Signal Scan", "Asteroid Prospect", "System Survey Bonus"};
    int idx = st.kindFilter;
    // Map: -1 => 0
    int comboIndex = (idx < 0) ? 0 : (idx + 1);
    ImGui::SetNextItemWidth(180.0f);
    if (ImGui::Combo("Kind", &comboIndex, items, (int)(sizeof(items) / sizeof(items[0])))) {
      st.kindFilter = (comboIndex == 0) ? -1 : (comboIndex - 1);
    }
  }

  ImGui::SameLine();
  {
    const char* sorts[] = {"Newest", "Highest value", "System name"};
    ImGui::SetNextItemWidth(150.0f);
    ImGui::Combo("Sort", &st.sortMode, sorts, (int)(sizeof(sorts) / sizeof(sorts[0])));
  }

  // Build filtered list.
  std::vector<const sim::LogbookEntry*> rows;
  rows.reserve(ctx.logbook.size());
  const std::string_view needle(st.filter);

  for (const auto& e : ctx.logbook) {
    if (!st.showSold && e.sold) continue;
    if (!st.showUnsold && !e.sold) continue;

    if (st.kindFilter >= 0) {
      if ((int)e.kind != st.kindFilter) continue;
    }

    if (!needle.empty()) {
      const auto& sys = ctx.universe.getSystem(e.systemId);
      const std::string sysName = sys.stub.name;
      const std::string target = entryTargetLabel(ctx.universe, e);
      if (!containsCaseInsensitive(sysName, needle) && !containsCaseInsensitive(target, needle)) {
        // For station scans, include station name too.
        if (e.kind == sim::LogbookEntryKind::StationScan && e.stationId != 0) {
          const std::string stName = stationNameInSystem(sys, e.stationId);
          if (!containsCaseInsensitive(stName, needle)) continue;
        } else {
          continue;
        }
      }
    }

    rows.push_back(&e);
  }

  // Sort.
  if (st.sortMode == 0) {
    std::sort(rows.begin(), rows.end(), [&](const sim::LogbookEntry* a, const sim::LogbookEntry* b) {
      if (a->discoveredDay != b->discoveredDay) return a->discoveredDay > b->discoveredDay;
      return a->key > b->key;
    });
  } else if (st.sortMode == 1) {
    std::sort(rows.begin(), rows.end(), [&](const sim::LogbookEntry* a, const sim::LogbookEntry* b) {
      if (a->valueCr != b->valueCr) return a->valueCr > b->valueCr;
      return a->discoveredDay > b->discoveredDay;
    });
  } else {
    std::sort(rows.begin(), rows.end(), [&](const sim::LogbookEntry* a, const sim::LogbookEntry* b) {
      const std::string an = systemName(ctx.universe, a->systemId);
      const std::string bn = systemName(ctx.universe, b->systemId);
      if (an != bn) return an < bn;
      return a->discoveredDay > b->discoveredDay;
    });
  }

  ImGui::Separator();

  const ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable |
                               ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_ScrollY;
  const float tableH = 280.0f;
  if (ImGui::BeginTable("##log_table", 6, flags, ImVec2(0, tableH))) {
    ImGui::TableSetupColumn("Day", ImGuiTableColumnFlags_WidthFixed, 60.0f);
    ImGui::TableSetupColumn("Kind", ImGuiTableColumnFlags_WidthFixed, 140.0f);
    ImGui::TableSetupColumn("System");
    ImGui::TableSetupColumn("Target");
    ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthFixed, 80.0f);
    ImGui::TableSetupColumn("Sold", ImGuiTableColumnFlags_WidthFixed, 45.0f);
    ImGui::TableHeadersRow();

    for (const auto* e : rows) {
      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      ImGui::Text("%.1f", e->discoveredDay);

      ImGui::TableSetColumnIndex(1);
      ImGui::TextUnformatted(sim::logbookEntryKindName(e->kind));

      ImGui::TableSetColumnIndex(2);
      const std::string sysName = systemName(ctx.universe, e->systemId);
      const bool selected = (st.selectedKey == e->key);
      if (ImGui::Selectable(sysName.c_str(), selected, ImGuiSelectableFlags_SpanAllColumns)) {
        st.selectedKey = e->key;
      }

      ImGui::TableSetColumnIndex(3);
      const std::string tgt = entryTargetLabel(ctx.universe, *e);
      ImGui::TextUnformatted(tgt.c_str());

      ImGui::TableSetColumnIndex(4);
      ImGui::Text("%.0f", e->valueCr);

      ImGui::TableSetColumnIndex(5);
      ImGui::TextUnformatted(e->sold ? "Y" : "");
    }

    ImGui::EndTable();
  }

  // Details panel
  ImGui::Separator();
  const sim::LogbookEntry* selected = nullptr;
  for (const auto& e : ctx.logbook) {
    if (e.key == st.selectedKey) {
      selected = &e;
      break;
    }
  }

  if (!selected) {
    ImGui::TextDisabled("Select an entry to see details.");
    ImGui::End();
    return;
  }

  const auto& sys = ctx.universe.getSystem(selected->systemId);
  ImGui::Text("%s", sys.stub.name.c_str());
  ImGui::TextDisabled("Kind: %s | Base value: %.0f cr | %s",
                      sim::logbookEntryKindName(selected->kind),
                      selected->valueCr,
                      selected->sold ? "SOLD" : "UNSOLD");

  if (selected->kind == sim::LogbookEntryKind::StationScan && selected->stationId != 0) {
    const std::string stName = stationNameInSystem(sys, selected->stationId);
    ImGui::Text("Station: %s", stName.c_str());
  }
  if (selected->kind == sim::LogbookEntryKind::AsteroidProspect) {
    const auto& cd = econ::commodityDef(selected->commodity);
    ImGui::Text("Yield: %s | Remaining ~%.0f u", cd.name, std::max(0.0, selected->units));
  }

  if (ctx.plotRouteToSystem) {
    if (ImGui::Button("Plot route to system")) {
      ctx.plotRouteToSystem(selected->systemId);
    }
  }

  if (ctx.targetStation && ctx.currentSystem && ctx.currentSystem->stub.id == selected->systemId
      && selected->stationId != 0) {
    ImGui::SameLine();
    if (ImGui::Button("Target station")) {
      ctx.targetStation(selected->stationId);
    }
  }

  ImGui::End();
}

} // namespace stellar::game
