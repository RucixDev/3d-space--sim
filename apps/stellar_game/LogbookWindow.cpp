#include "LogbookWindow.h"

#include "stellar/econ/Commodity.h"
#include "stellar/sim/ExplorationData.h"
#include "stellar/sim/Signals.h"

#include <imgui.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>
#include <unordered_map>
#include <utility>

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

static double brokerMultiplierForEntry(const ExplorationDataBrokerState& st,
                                      const sim::StarSystem& saleSystem,
                                      const sim::Station& saleStation,
                                      const sim::SystemStub& scanStub,
                                      const sim::LogbookEntry& e) {
  if (!st.enableMultipliers) return 1.0;
  return sim::explorationDataBrokerMultiplier(st.params,
                                              saleSystem.stub,
                                              saleStation,
                                              scanStub,
                                              e.kind,
                                              e.subKind);
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
  std::unordered_map<sim::SystemId, sim::SystemStub> stubCache;
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

    // Payout is entry-specific (station demand may depend on the entry kind).
    auto itStub = stubCache.find(e.systemId);
    if (itStub == stubCache.end()) {
      const auto& sys = universe.getSystem(e.systemId);
      itStub = stubCache.emplace(e.systemId, sys.stub).first;
    }
    const double mult = brokerMultiplierForEntry(st, saleSystem, saleStation, itStub->second, e);
    g.payoutCr += e.valueCr * mult;
  }

  // Enrich groups with names/distances/effective multiplier.
  for (auto& g : groups) {
    const auto itStub = stubCache.find(g.systemId);
    const sim::SystemStub& stub = (itStub != stubCache.end()) ? itStub->second : universe.getSystem(g.systemId).stub;
    g.name = stub.name;
    g.factionId = stub.factionId;
    g.distanceLy = (stub.posLy - saleSystem.stub.posLy).length();
    g.mult = (g.baseCr > 1e-6) ? (g.payoutCr / g.baseCr) : 1.0;
  }

  // Sort by base value (descending).
  std::sort(groups.begin(), groups.end(), [&](const Group& a, const Group& b) {
    if (a.baseCr != b.baseCr) return a.baseCr > b.baseCr;
    return a.name < b.name;
  });

  // Pricing model.
  ImGui::Separator();
  if (ImGui::TreeNodeEx("Broker Pricing Model", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Checkbox("Enable broker multipliers", &st.enableMultipliers);
    ImGui::SameLine();
    ImGui::Checkbox("Show details", &st.showDetails);

    ImGui::BeginDisabled(!st.enableMultipliers);

    ImGui::Checkbox("Distance premium", &st.params.enableDistancePremium);
    float distScaleLy = (float)st.params.distanceScaleLy;
    if (ImGui::SliderFloat("Distance scale (ly)", &distScaleLy, 50.0f, 800.0f, "%.0f")) {
      st.params.distanceScaleLy = (double)distScaleLy;
    }
    float maxDistPrem = (float)st.params.maxDistancePremium;
    if (ImGui::SliderFloat("Max distance premium", &maxDistPrem, 0.0f, 0.75f, "+%.0f%%")) {
      st.params.maxDistancePremium = (double)maxDistPrem;
    }

    float sameFaction = (float)st.params.sameFactionBonus;
    if (ImGui::SliderFloat("Same-faction bonus", &sameFaction, 0.0f, 0.50f, "+%.0f%%")) {
      st.params.sameFactionBonus = (double)sameFaction;
    }
    float otherFaction = (float)st.params.otherFactionPenalty;
    if (ImGui::SliderFloat("Other-faction penalty", &otherFaction, 0.0f, 0.50f, "-%.0f%%")) {
      st.params.otherFactionPenalty = (double)otherFaction;
    }

    ImGui::Separator();
    ImGui::Checkbox("Station demand shaping", &st.params.enableStationDemand);
    float demand = (float)st.params.demandStrength;
    if (ImGui::SliderFloat("Demand strength", &demand, 0.0f, 1.0f, "%.2f")) {
      st.params.demandStrength = (double)demand;
    }

    ImGui::EndDisabled();
    ImGui::TextDisabled("Payout = base × multiplier (distance, faction, demand; clamped). Bank stays in base credits.");
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
    for (auto& e : logbook) {
      if (e.sold) continue;
      if (e.valueCr <= 0.0) continue;
      if (!sellAll) {
        if (st.selectedSystems.find(e.systemId) == st.selectedSystems.end()) continue;
      }
      auto itStub = stubCache.find(e.systemId);
      if (itStub == stubCache.end()) {
        const auto& sys = universe.getSystem(e.systemId);
        itStub = stubCache.emplace(e.systemId, sys.stub).first;
      }
      const double mult = brokerMultiplierForEntry(st, saleSystem, saleStation, itStub->second, e);
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

  // ---------------------------------------------------------------------------
  // Sell Planner
  // ---------------------------------------------------------------------------
  // Ranks nearby stations by predicted payout using the same broker multiplier
  // as the station-side sale panel. This lets players plan a "data run".
  ImGui::Separator();
  if (ImGui::CollapsingHeader("Sell Planner", ImGuiTreeNodeFlags_DefaultOpen)) {
    auto& br = st.broker;

    if (!ctx.currentSystem) {
      ImGui::TextDisabled("No current system loaded yet.");
    } else {
      ImGui::Checkbox("Enable broker multipliers##planner", &br.enableMultipliers);
      ImGui::SameLine();
      ImGui::TextDisabled("(uses station demand + distance + faction)");

      ImGui::SetNextItemWidth(160.0f);
      ImGui::SliderFloat("Search radius (ly)", &br.plannerRadiusLy, 25.0f, 2000.0f, "%.0f");
      ImGui::SameLine();
      ImGui::SetNextItemWidth(120.0f);
      ImGui::SliderInt("Max systems", &br.plannerMaxSystems, 16, 256);
      ImGui::SameLine();
      ImGui::SetNextItemWidth(120.0f);
      ImGui::SliderInt("Max stations", &br.plannerMaxStations, 8, 256);

      const bool useSelection = !br.selectedSystems.empty();
      if (useSelection) {
        ImGui::TextDisabled("Using station-broker selection: %d systems", (int)br.selectedSystems.size());
      }

      struct BucketKey {
        sim::SystemId sys{0};
        sim::LogbookEntryKind kind{sim::LogbookEntryKind::StarScan};
        core::u8 sub{0};
        bool operator==(const BucketKey& o) const { return sys == o.sys && kind == o.kind && sub == o.sub; }
      };
      struct BucketKeyHash {
        std::size_t operator()(const BucketKey& k) const {
          std::size_t h = std::hash<sim::SystemId>{}(k.sys);
          h ^= (std::size_t)k.kind + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
          h ^= (std::size_t)k.sub + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
          return h;
        }
      };

      struct Bucket {
        BucketKey key{};
        double baseCr{0.0};
        int entries{0};
      };

      std::unordered_map<BucketKey, std::size_t, BucketKeyHash> bucketIdx;
      std::vector<Bucket> buckets;
      buckets.reserve(128);

      double baseTotal = 0.0;
      int entryCount = 0;
      for (const auto& e : ctx.logbook) {
        if (e.sold) continue;
        if (e.valueCr <= 0.0) continue;
        if (e.systemId == 0) continue;
        if (useSelection && br.selectedSystems.find(e.systemId) == br.selectedSystems.end()) continue;

        BucketKey k{e.systemId, e.kind, e.subKind};
        auto it = bucketIdx.find(k);
        if (it == bucketIdx.end()) {
          const std::size_t ni = buckets.size();
          bucketIdx.emplace(k, ni);
          buckets.push_back(Bucket{});
          buckets.back().key = k;
          it = bucketIdx.find(k);
        }
        Bucket& b = buckets[it->second];
        b.baseCr += e.valueCr;
        b.entries += 1;
        baseTotal += e.valueCr;
        entryCount += 1;
      }

      if (baseTotal <= 1e-6 || buckets.empty()) {
        ImGui::TextDisabled("No unsold exploration data to plan a sale.");
      } else {
        // Cache scan system stubs referenced by buckets.
        std::unordered_map<sim::SystemId, sim::SystemStub> scanStubCache;
        scanStubCache.reserve(buckets.size());
        for (const auto& b : buckets) {
          if (scanStubCache.find(b.key.sys) == scanStubCache.end()) {
            const auto& s = ctx.universe.getSystem(b.key.sys);
            scanStubCache.emplace(b.key.sys, s.stub);
          }
        }

        // Candidate station search.
        const auto nearby = ctx.universe.queryNearby(ctx.currentSystem->stub.posLy, (double)br.plannerRadiusLy,
                                                     (std::size_t)std::max(1, br.plannerMaxSystems));

        struct Candidate {
          sim::SystemId sysId{0};
          sim::StationId stationId{0};
          econ::StationType type{econ::StationType::Outpost};
          std::string systemName;
          std::string stationName;
          double travelLy{0.0};
          double payoutCr{0.0};
        };

        std::vector<Candidate> candidates;
        candidates.reserve((std::size_t)std::max(8, br.plannerMaxStations));

        auto stationTypeName = [](econ::StationType t) {
          switch (t) {
            case econ::StationType::Agricultural: return "Agricultural";
            case econ::StationType::Mining: return "Mining";
            case econ::StationType::Refinery: return "Refinery";
            case econ::StationType::Industrial: return "Industrial";
            case econ::StationType::Research: return "Research";
            case econ::StationType::TradeHub: return "Trade Hub";
            case econ::StationType::Shipyard: return "Shipyard";
            case econ::StationType::Outpost:
            default: return "Outpost";
          }
        };

        for (const auto& stub : nearby) {
          const auto& sys = ctx.universe.getSystem(stub.id, &stub);
          const double travelLy = (sys.stub.posLy - ctx.currentSystem->stub.posLy).length();
          for (const auto& stn : sys.stations) {
            double payout = 0.0;
            for (const auto& b : buckets) {
              const sim::SystemStub& scanStub = scanStubCache[b.key.sys];
              const double mult = br.enableMultipliers
                                    ? sim::explorationDataBrokerMultiplier(br.params, sys.stub, stn, scanStub, b.key.kind, b.key.sub)
                                    : 1.0;
              payout += b.baseCr * mult;
            }

            Candidate c;
            c.sysId = sys.stub.id;
            c.stationId = stn.id;
            c.type = stn.type;
            c.systemName = sys.stub.name;
            c.stationName = stn.name;
            c.travelLy = travelLy;
            c.payoutCr = payout;
            candidates.push_back(std::move(c));
          }
        }

        if (candidates.empty()) {
          ImGui::TextDisabled("No candidate stations found in radius.");
        } else {
          std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
            if (a.payoutCr != b.payoutCr) return a.payoutCr > b.payoutCr;
            if (a.travelLy != b.travelLy) return a.travelLy < b.travelLy;
            return a.stationName < b.stationName;
          });
          if ((int)candidates.size() > br.plannerMaxStations) {
            candidates.resize((std::size_t)std::max(1, br.plannerMaxStations));
          }

          const Candidate& best = candidates.front();
          ImGui::Text("Unsold: %d entries | base %.0f cr", entryCount, baseTotal);
          ImGui::Text("Best nearby buyer: %s (%s) | payout %.0f cr (x%.2f) | travel %.0f ly",
                      best.stationName.c_str(), best.systemName.c_str(), best.payoutCr,
                      (baseTotal > 1e-6) ? (best.payoutCr / baseTotal) : 1.0, best.travelLy);

          const ImGuiTableFlags tflags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp;
          if (ImGui::BeginTable("##sell_planner", 7, tflags)) {
            ImGui::TableSetupColumn("Rank", ImGuiTableColumnFlags_WidthFixed, 40.0f);
            ImGui::TableSetupColumn("Station");
            ImGui::TableSetupColumn("System");
            ImGui::TableSetupColumn("Type", ImGuiTableColumnFlags_WidthFixed, 86.0f);
            ImGui::TableSetupColumn("Payout", ImGuiTableColumnFlags_WidthFixed, 90.0f);
            ImGui::TableSetupColumn("Mult", ImGuiTableColumnFlags_WidthFixed, 60.0f);
            ImGui::TableSetupColumn("Go", ImGuiTableColumnFlags_WidthFixed, 92.0f);
            ImGui::TableHeadersRow();

            for (std::size_t i = 0; i < candidates.size(); ++i) {
              const auto& c = candidates[i];
              ImGui::TableNextRow();

              ImGui::TableSetColumnIndex(0);
              ImGui::Text("%d", (int)i + 1);

              ImGui::TableSetColumnIndex(1);
              ImGui::TextUnformatted(c.stationName.c_str());

              ImGui::TableSetColumnIndex(2);
              ImGui::TextUnformatted(c.systemName.c_str());

              ImGui::TableSetColumnIndex(3);
              ImGui::TextUnformatted(stationTypeName(c.type));

              ImGui::TableSetColumnIndex(4);
              ImGui::Text("%.0f", c.payoutCr);

              ImGui::TableSetColumnIndex(5);
              ImGui::Text("x%.2f", (baseTotal > 1e-6) ? (c.payoutCr / baseTotal) : 1.0);

              ImGui::TableSetColumnIndex(6);
              const std::string btn = "Plot##plan_" + std::to_string((unsigned long long)c.stationId);
              if (ImGui::SmallButton(btn.c_str())) {
                if (ctx.plotRouteToSystem) {
                  (void)ctx.plotRouteToSystem(c.sysId);
                }
                if (ctx.targetStation) {
                  ctx.targetStation(c.stationId);
                }
              }
              ImGui::SameLine();
              ImGui::TextDisabled("%.0f ly", c.travelLy);
            }

            ImGui::EndTable();
          }
        }
      }
    }
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
