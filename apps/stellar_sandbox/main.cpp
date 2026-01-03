#include "stellar/core/Args.h"
#include "stellar/core/JobSystem.h"
#include "stellar/core/JsonWriter.h"
#include "stellar/core/Log.h"
#include "stellar/econ/Commodity.h"
#include "stellar/econ/RoutePlanner.h"
#include "stellar/proc/GalaxyGenerator.h"
#include "stellar/sim/MissionLogic.h"
#include "stellar/sim/Contraband.h"
#include "stellar/sim/Law.h"
#include "stellar/sim/FactionProfile.h"
#include "stellar/sim/PoliceScan.h"
#include "stellar/sim/SecurityModel.h"
#include "stellar/sim/Industry.h"
#include "stellar/sim/Warehouse.h"
#include "stellar/sim/Reputation.h"
#include "stellar/sim/TradeScanner.h"
#include "stellar/sim/TradeLoopScanner.h"
#include "stellar/sim/TradeRunPlanner.h"
#include "stellar/sim/IndustryScanner.h"
#include "stellar/sim/NavRoute.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/Signature.h"
#include "stellar/sim/Signals.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <memory>
#include <string>
#include <vector>

using namespace stellar;

static const char* starClassName(sim::StarClass c) {
  switch (c) {
    case sim::StarClass::O: return "O";
    case sim::StarClass::B: return "B";
    case sim::StarClass::A: return "A";
    case sim::StarClass::F: return "F";
    case sim::StarClass::G: return "G";
    case sim::StarClass::K: return "K";
    case sim::StarClass::M: return "M";
    default: return "?";
  }
}

static const char* missionTypeName(sim::MissionType t) {
  using sim::MissionType;
  switch (t) {
    case MissionType::Courier: return "Courier";
    case MissionType::Delivery: return "Delivery";
    case MissionType::BountyScan: return "BountyScan";
    case MissionType::BountyKill: return "BountyKill";
    case MissionType::MultiDelivery: return "MultiDelivery";
    case MissionType::Passenger: return "Passenger";
    case MissionType::Smuggle: return "Smuggle";
    case MissionType::Salvage: return "Salvage";
    case MissionType::Escort: return "Escort";
    default: return "Unknown";
  }
}

static const char* stationTypeName(econ::StationType t) {
  using econ::StationType;
  switch (t) {
    case StationType::Outpost: return "Outpost";
    case StationType::Agricultural: return "Agricultural";
    case StationType::Mining: return "Mining";
    case StationType::Refinery: return "Refinery";
    case StationType::Industrial: return "Industrial";
    case StationType::Research: return "Research";
    case StationType::TradeHub: return "TradeHub";
    case StationType::Shipyard: return "Shipyard";
    default: return "?";
  }
}

static void printMission(const sim::Mission& m) {
  std::cout << "  [" << missionTypeName(m.type) << "] "
            << "toSystem=" << m.toSystem << " toStation=" << m.toStation;
  if (m.type == sim::MissionType::Delivery || m.type == sim::MissionType::MultiDelivery || m.type == sim::MissionType::Smuggle || m.type == sim::MissionType::Salvage || m.type == sim::MissionType::Escort) {
    std::cout << " cargo=" << econ::commodityCode(m.commodity) << " units=" << std::fixed << std::setprecision(0)
              << m.units;
  }
  if (m.type == sim::MissionType::Passenger) {
    std::cout << " party=" << std::fixed << std::setprecision(0) << m.units << " seats";
  }
  if (m.type == sim::MissionType::BountyScan || m.type == sim::MissionType::BountyKill || m.type == sim::MissionType::Salvage || m.type == sim::MissionType::Escort) {
    std::cout << " targetNpcId=" << m.targetNpcId;
  }
  if (m.type == sim::MissionType::Salvage) {
    std::cout << " siteVisited=" << (m.scanned ? "yes" : "no");
  }
  if (m.type == sim::MissionType::MultiDelivery && m.viaStation != 0) {
    std::cout << " via=" << m.viaSystem << "/" << m.viaStation << " leg=" << (int)m.leg;
  }
  std::cout << " reward=" << std::fixed << std::setprecision(0) << m.reward
            << " deadline=" << std::fixed << std::setprecision(1) << m.deadlineDay
            << (m.cargoProvided ? " (cargo-provided)" : " (source-cargo)")
            << "\n";
}

static void writeMissionJson(core::JsonWriter& j, const sim::Mission& m) {
  j.beginObject();
  j.key("id"); j.value((unsigned long long)m.id);
  j.key("type"); j.value(missionTypeName(m.type));
  j.key("factionId"); j.value((unsigned long long)m.factionId);
  j.key("fromSystem"); j.value((unsigned long long)m.fromSystem);
  j.key("fromStation"); j.value((unsigned long long)m.fromStation);
  j.key("toSystem"); j.value((unsigned long long)m.toSystem);
  j.key("toStation"); j.value((unsigned long long)m.toStation);
  j.key("reward"); j.value(m.reward);
  j.key("deadlineDay"); j.value(m.deadlineDay);
  j.key("cargoProvided"); j.value(m.cargoProvided);
  j.key("completed"); j.value(m.completed);
  j.key("failed"); j.value(m.failed);
  j.key("commodity"); j.value(econ::commodityCode(m.commodity));
  j.key("units"); j.value(m.units);
  j.key("viaSystem"); j.value((unsigned long long)m.viaSystem);
  j.key("viaStation"); j.value((unsigned long long)m.viaStation);
  j.key("leg"); j.value((int)m.leg);
  j.key("targetNpcId"); j.value((unsigned long long)m.targetNpcId);
  j.key("scanned"); j.value(m.scanned);
  j.endObject();
}

static void printHelp() {
  std::cout << "stellar_sandbox\n"
            << "  --seed <u64>           Galaxy seed (default: 1337)\n"
            << "  --pos <x y z>          Query position in ly (default: 0 0 0)\n"
            << "  --radius <ly>          Query radius in ly (default: 50)\n"
            << "  --limit <n>            Max systems (default: 32)\n"
            << "  --day <days>           Economy time (days) for trade quotes (default: 0)\n"
            << "  --threads <n>          Use a thread pool for bulk queries (0=auto; default: disabled)\n"
            << "  --sig                  Include a stable signature per system stub (determinism)\n"
            << "  --sysSig               Include a deep signature of full generated systems (slower)\n"
            << "  --json                 Emit machine-readable JSON to stdout (also works with --out)\n"
            << "  --out <path>           Write JSON output to a file instead of stdout ('-' means stdout)\n"
            << "\n"
            << "Trade scanning (headless route planner):\n"
            << "  --trade                Print best trade opportunities from a chosen station\n"
            << "  --tradeMix             Print best multi-commodity cargo manifests (trade mix)\n"
            << "  --tradeLoop            Print best closed trade loops (2- or 3-leg) starting at the origin station\n"
            << "  --tradeRun             Plan best multi-leg trade runs (A->B->C...) from the origin station\n"
            << "  --loopLegs <n>         Loop leg count (2 or 3) (default: 2)\n"
            << "  --loopLimit <n>        Max loops to print (default: 8)\n"
            << "  --loopLegCandidates <n> Max outgoing leg candidates expanded per stop (default: 16)\n"
            << "  --loopMinLegProfit <cr> Minimum net profit per leg (default: 0)\n"
            << "  --loopMinProfit <cr>   Minimum total net profit for the loop (default: 0)\n"
            << "  --runLegs <n>          Run leg count (1..4) (default: 2)\n"
            << "  --runLimit <n>         Max runs to print (default: 8)\n"
            << "  --runBeam <n>          Beam width (default: 32)\n"
            << "  --runLegCandidates <n> Max outgoing leg candidates expanded per stop (default: 16)\n"
            << "  --runMinLegProfit <cr> Minimum net profit per leg (default: 0)\n"
            << "  --runMinProfit <cr>    Minimum total net profit for the run (default: 0)\n"
            << "  --runScore <mode>      profit | profitLy | profitHop | profitCost (default: profit)\n"
            << "  --fromSys <idx>        Origin system index in the printed list (default: 0)\n"
            << "  --fromStation <idx>    Origin station index inside that system (default: 0)\n"
            << "  --cargoKg <kg>         Cargo capacity (kg) used for trip profit (default: 420)\n"
            << "  --tradeLimit <n>       Max trade opportunities to print (default: 12)\n"
            << "  --mixStepKg <kg>       Planner resolution in kg (default: 1.0)\n"
            << "  --mixNoImpact          Disable price-impact simulation (static quotes)\n"
            << "  --mixMaxSpend <cr>     Optional max buy budget (default: 0 = unlimited)\n"
            << "  --mixLines <n>         Max manifest lines to print per destination (default: 6)\n"
            << "  --commodity <c>        Optional filter (code or name), e.g. FOOD | H2O | water\n"
            << "  --minProfit <cr>       Min net trip profit filter (default: 0)\n"
            << "\n"
            << "Jump route planning (A*):\n"
            << "  --route                Plot a jump route between two systems in the printed list\n"
            << "  --toSys <idx>          Destination system index in the printed list (default: 0)\n"
            << "  --jr <ly>              Jump range in ly (default: 18)\n"
            << "  --kRoutes <n>         Print up to n alternative routes (default: 1)\n"
            << "  --routeCost <mode>     Route cost: hops | dist | fuel (default: hops)\n"
            << "  --fuelBase <u>         Base cost per jump for --routeCost fuel (default: 2.0)\n"
            << "  --fuelPerLy <u>        Cost per ly for --routeCost fuel (default: 0.5)\n"
            << "\n"
            << "Mission board (headless):\n"
            << "  --missions             Print mission offers for the chosen station\n"
            << "  --rep <r>              Player reputation used for offer count/rewards (default: 0)\n"
            << "  --credits <cr>         Starting credits when simulating acceptance (default: 10000)\n"
            << "  --seats <n>            Passenger seats when simulating acceptance (default: 2)\n"
            << "  --acceptOffer <idx>    Accept offer index and print resulting state\n"
            << "  --advanceDays <d>      Advance time by d days before auto-complete (default: 0)\n"
            << "  --autoComplete         After accepting, attempt to complete at destination\n"
            << "\n"
            << "Signals / space objects (headless):\n"
            << "  --signals             Print in-system signal sites (resource fields, derelicts, distress, mission salvage)\n"
            << "  --signalFields <n>    Number of persistent resource fields to generate (default: 1)\n"
            << "  --signalDistress <n>  Distress calls per day (default: 1)\n"
            << "  --signalDistressTtl <d> Distress time-to-live in days (default: 1.0)\n"
            << "  --signalNoDistress    Disable distress call generation\n"
            << "  --signalNoDerelict    Disable daily derelict generation\n"
            << "\n"
            << "Law / contraband (headless):\n"
            << "  --law                 Print law + contraband profile for the chosen station\n"
            << "  --faction <id>         Override faction id (default: station faction)\n"
            << "  --lawRep <r>           Player reputation used for bribe chance (default: 0)\n"
            << "  --heat <h>             Player heat used for bribe chance (default: 0)\n"
            << "  --illegalValue <cr>    Illegal cargo value used for fine/bribe examples (default: 5000)\n"
            << "  --smuggleHoldMk <mk>   Smuggle hold grade (0-3) for scan math (default: 0)\n"
            << "\n"
            << "Factions / diplomacy (tooling):\n"
            << "  --factions            Print all factions and their procedural profiles\n"
            << "  --diplomacy           With --factions, print top allies/hostiles per faction\n"
            << "\n"
            << "System security (tooling):\n"
            << "  --security            Include per-system security/piracy/traffic metrics\n"
            << "\n"
            << "Save/load (tooling):\n"
            << "  --load <path>          Load a save file before running missions\n"
            << "  --save <path>          Save the resulting state to a save file\n"
            << "\n"
            << "Industry (station processing):\n"
            << "  --industry             Print available industry recipes for the chosen station\n"
            << "  --industryTrade        Scan profitable process+haul routes (buy inputs @ origin, sell output elsewhere)\n"
            << "  --batches <n>           Batches to quote (default: 1)\n"
            << "\n"
            << "Warehouse (station storage):\n"
            << "  --warehouse            Quote storage fees for warehousing at the chosen station\n"
            << "  --storeCommodity <c>   Commodity to store for the quote (code or name; default: FOOD)\n"
            << "  --storeUnits <u>       Units to store for the quote (default: 10)\n";
}

int main(int argc, char** argv) {
  core::setLogLevel(core::LogLevel::Info);

  core::Args args;
  args.setArity("pos", 3);
  args.parse(argc, argv);

  if (args.hasFlag("help") || args.hasFlag("h")) {
    printHelp();
    return 0;
  }

  core::u64 seed = 1337;
  {
    unsigned long long s = (unsigned long long)seed;
    (void)args.getU64("seed", s);
    seed = (core::u64)s;
  }

  math::Vec3d posLy{0,0,0};
  {
    const auto p = args.values("pos");
    if (p.size() >= 3) {
      posLy.x = std::atof(p[p.size() - 3].c_str());
      posLy.y = std::atof(p[p.size() - 2].c_str());
      posLy.z = std::atof(p[p.size() - 1].c_str());
    }
  }

  double radiusLy = 50.0;
  (void)args.getDouble("radius", radiusLy);

  std::size_t limit = 32;
  {
    unsigned long long v = 0;
    if (args.getU64("limit", v)) limit = (std::size_t)v;
  }

  double timeDays = 0.0;
  (void)args.getDouble("day", timeDays);

  // Optional parallelism for bulk queries.
  bool wantThreads = false;
  std::size_t threads = 0;
  {
    unsigned long long v = 0;
    if (args.getU64("threads", v)) {
      wantThreads = true;
      threads = static_cast<std::size_t>(v);
    }
  }

  const bool json = args.hasFlag("json");
  std::string outPath;
  (void)args.getString("out", outPath);

  const bool emitStubSig = args.hasFlag("sig");
  const bool emitSysSig = args.hasFlag("sysSig");

  const bool doTrade = args.hasFlag("trade");
  const bool doSignals = args.hasFlag("signals");
  const bool doTradeMix = args.hasFlag("tradeMix");
  const bool doTradeLoop = args.hasFlag("tradeLoop");
  const bool doTradeRun = args.hasFlag("tradeRun");
  const bool doIndustry = args.hasFlag("industry");
  const bool doIndustryTrade = args.hasFlag("industryTrade");
  const bool doWarehouse = args.hasFlag("warehouse");
  const bool doFactions = args.hasFlag("factions");
  const bool doDiplomacy = args.hasFlag("diplomacy");
  const bool doSecurity = args.hasFlag("security");
  std::size_t fromSysIdx = 0;
  std::size_t fromStationIdx = 0;
  {
    unsigned long long v = 0;
    if (args.getU64("fromSys", v)) fromSysIdx = (std::size_t)v;
    if (args.getU64("fromStation", v)) fromStationIdx = (std::size_t)v;
  }

  // Signal generator settings.
  int signalFieldCount = 1;
  int signalDistressPerDay = 1;
  double signalDistressTtlDays = 1.0;
  {
    unsigned long long v = 0;
    if (args.getU64("signalFields", v)) signalFieldCount = (int)v;
    if (args.getU64("signalDistress", v)) signalDistressPerDay = (int)v;
  }
  (void)args.getDouble("signalDistressTtl", signalDistressTtlDays);
  if (!std::isfinite(signalDistressTtlDays)) signalDistressTtlDays = 1.0;
  signalDistressTtlDays = std::clamp(signalDistressTtlDays, 0.1, 10.0);
  const bool signalNoDistress = args.hasFlag("signalNoDistress");
  const bool signalNoDerelict = args.hasFlag("signalNoDerelict");
  double cargoCapacityKg = 420.0;
  (void)args.getDouble("cargoKg", cargoCapacityKg);
  std::size_t tradeLimit = 12;
  {
    unsigned long long v = 0;
    if (args.getU64("tradeLimit", v)) tradeLimit = (std::size_t)v;
  }

  // Trade-loop settings.
  std::size_t loopLimit = 8;
  std::size_t loopLegs = 2;
  std::size_t loopLegCandidates = 16;
  {
    unsigned long long v = 0;
    if (args.getU64("loopLimit", v)) loopLimit = (std::size_t)v;
    if (args.getU64("loopLegs", v)) loopLegs = (std::size_t)v;
    if (args.getU64("loopLegCandidates", v)) loopLegCandidates = (std::size_t)v;
  }
  double loopMinLegProfit = 0.0;
  double loopMinProfit = 0.0;
  (void)args.getDouble("loopMinLegProfit", loopMinLegProfit);
  (void)args.getDouble("loopMinProfit", loopMinProfit);
  if (!std::isfinite(loopMinLegProfit)) loopMinLegProfit = 0.0;
  if (!std::isfinite(loopMinProfit)) loopMinProfit = 0.0;
  loopMinLegProfit = std::max(0.0, loopMinLegProfit);
  loopMinProfit = std::max(0.0, loopMinProfit);


  // Trade-run settings.
  std::size_t runLimit = 8;
  std::size_t runLegs = 2;
  std::size_t runBeam = 32;
  std::size_t runLegCandidates = 16;
  {
    unsigned long long v = 0;
    if (args.getU64("runLimit", v)) runLimit = (std::size_t)v;
    if (args.getU64("runLegs", v)) runLegs = (std::size_t)v;
    if (args.getU64("runBeam", v)) runBeam = (std::size_t)v;
    if (args.getU64("runLegCandidates", v)) runLegCandidates = (std::size_t)v;
  }
  double runMinLegProfit = 0.0;
  double runMinProfit = 0.0;
  (void)args.getDouble("runMinLegProfit", runMinLegProfit);
  (void)args.getDouble("runMinProfit", runMinProfit);
  if (!std::isfinite(runMinLegProfit)) runMinLegProfit = 0.0;
  if (!std::isfinite(runMinProfit)) runMinProfit = 0.0;
  runMinLegProfit = std::max(0.0, runMinLegProfit);
  runMinProfit = std::max(0.0, runMinProfit);
  std::string runScore = "profit";
  (void)args.getString("runScore", runScore);


  // Trade-mix planner settings.
  double mixStepKg = 1.0;
  (void)args.getDouble("mixStepKg", mixStepKg);
  if (!std::isfinite(mixStepKg)) mixStepKg = 1.0;
  mixStepKg = std::max(0.05, mixStepKg);

  const bool mixNoImpact = args.hasFlag("mixNoImpact");

  double mixMaxSpend = 0.0;
  (void)args.getDouble("mixMaxSpend", mixMaxSpend);
  if (!std::isfinite(mixMaxSpend)) mixMaxSpend = 0.0;
  mixMaxSpend = std::max(0.0, mixMaxSpend);

  int mixLines = 6;
  (void)args.getInt("mixLines", mixLines);
  mixLines = std::max(0, mixLines);


  // Optional trade-scan filters.
  std::string commodityFilter;
  (void)args.getString("commodity", commodityFilter);

  double minProfit = 0.0;
  (void)args.getDouble("minProfit", minProfit);
  if (!std::isfinite(minProfit)) minProfit = 0.0;
  minProfit = std::max(0.0, minProfit);

  bool hasCommodityFilter = false;
  econ::CommodityId commodityFilterId = econ::CommodityId::Food;
  if (!commodityFilter.empty()) {
    hasCommodityFilter = econ::tryParseCommodity(commodityFilter, commodityFilterId);
    if (!hasCommodityFilter) {
      std::cerr << "Invalid --commodity: '" << commodityFilter << "'\n";
      std::cerr << "Valid codes:";
      for (const auto& def : econ::commodityTable()) {
        std::cerr << " " << def.code;
      }
      std::cerr << "\n";
      return 1;
    }
  }

  const bool doRoute = args.hasFlag("route");
  std::size_t toSysIdx = 0;
  {
    unsigned long long v = 0;
    if (args.getU64("toSys", v)) toSysIdx = (std::size_t)v;
  }
  double jumpRangeLy = 18.0;
  (void)args.getDouble("jr", jumpRangeLy);

  std::size_t kRoutes = 1;
  {
    unsigned long long v = 1;
    if (args.getU64("kRoutes", v)) kRoutes = (std::size_t)std::max<unsigned long long>(1ULL, v);
  }

  std::string routeCost = "hops";
  (void)args.getString("routeCost", routeCost);
  double fuelBase = 2.0;
  double fuelPerLy = 0.5;
  (void)args.getDouble("fuelBase", fuelBase);
  (void)args.getDouble("fuelPerLy", fuelPerLy);

  const bool doMissions = args.hasFlag("missions");
  double rep = 0.0;
  (void)args.getDouble("rep", rep);
  double credits = 10000.0;
  (void)args.getDouble("credits", credits);
  int seats = 2;
  (void)args.getInt("seats", seats);
  int acceptOffer = -1;
  (void)args.getInt("acceptOffer", acceptOffer);

  const bool doLaw = args.hasFlag("law");
  core::u32 lawFactionOverride = 0;
  {
    unsigned long long v = 0;
    if (args.getU64("faction", v)) lawFactionOverride = (core::u32)v;
  }
  double lawRep = 0.0;
  (void)args.getDouble("lawRep", lawRep);
  double lawHeat = 0.0;
  (void)args.getDouble("heat", lawHeat);
  double lawIllegalValueCr = 5000.0;
  (void)args.getDouble("illegalValue", lawIllegalValueCr);
  if (!std::isfinite(lawIllegalValueCr)) lawIllegalValueCr = 0.0;
  lawIllegalValueCr = std::max(0.0, lawIllegalValueCr);
  int smuggleHoldMk = 0;
  (void)args.getInt("smuggleHoldMk", smuggleHoldMk);
  smuggleHoldMk = std::clamp(smuggleHoldMk, 0, 3);

  double industryBatches = 1.0;
  (void)args.getDouble("batches", industryBatches);
  if (!std::isfinite(industryBatches)) industryBatches = 1.0;
  industryBatches = std::max(0.0, industryBatches);

  // Warehouse quote inputs. Uses --rep for reputation and --advanceDays for projection.
  std::string storeCommodityStr;
  (void)args.getString("storeCommodity", storeCommodityStr);
  double storeUnits = 10.0;
  (void)args.getDouble("storeUnits", storeUnits);
  if (!std::isfinite(storeUnits)) storeUnits = 0.0;
  storeUnits = std::max(0.0, storeUnits);

  econ::CommodityId storeCommodity = econ::CommodityId::Food;
  if (!storeCommodityStr.empty()) {
    if (!econ::tryParseCommodity(storeCommodityStr, storeCommodity)) {
      std::cerr << "Invalid --storeCommodity: '" << storeCommodityStr << "'\n";
      return 1;
    }
  }
  double advanceDays = 0.0;
  (void)args.getDouble("advanceDays", advanceDays);
  const bool autoComplete = args.hasFlag("autoComplete");
  std::string loadPath;
  std::string savePath;
  (void)args.getString("load", loadPath);
  (void)args.getString("save", savePath);

  std::unique_ptr<core::JobSystem> jobs;
  if (wantThreads && (threads == 0 || threads > 1)) {
    jobs = std::make_unique<core::JobSystem>(threads);
  }

  sim::Universe u(seed);

  const auto systems = jobs ? u.queryNearbyParallel(*jobs, posLy, radiusLy, limit)
                            : u.queryNearby(posLy, radiusLy, limit);

  // If JSON is requested, we build up a machine-readable object.
  // Otherwise, we keep the original human-readable console output.
  std::unique_ptr<std::ofstream> jsonFile;
  std::ostream* jsonStream = &std::cout;
  if (json && !outPath.empty() && outPath != "-") {
    jsonFile = std::make_unique<std::ofstream>(outPath, std::ios::out | std::ios::trunc);
    if (!*jsonFile) {
      std::cerr << "Failed to open --out file: " << outPath << "\n";
      return 1;
    }
    jsonStream = jsonFile.get();
  }

  core::JsonWriter j(*jsonStream, /*pretty=*/true);
  if (json) {
    j.beginObject();
    j.key("seed"); j.value((unsigned long long)seed);
    j.key("query");
    j.beginObject();
    j.key("posLy");
    j.beginArray();
    j.value(posLy.x); j.value(posLy.y); j.value(posLy.z);
    j.endArray();
    j.key("radiusLy"); j.value(radiusLy);
    j.key("limit"); j.value((unsigned long long)limit);
    j.key("day"); j.value(timeDays);
    if (jobs) {
      j.key("threads");
      j.value((unsigned long long)jobs->threadCount());
    }
    j.endObject();
    j.key("systems");
    j.beginArray();
    for (const auto& s : systems) {
      const double distLy = (s.posLy - posLy).length();
      j.beginObject();
      j.key("id"); j.value((unsigned long long)s.id);
      j.key("name"); j.value(s.name);
      j.key("class"); j.value(starClassName(s.primaryClass));
      j.key("distLy"); j.value(distLy);
      j.key("planetCount"); j.value((unsigned long long)s.planetCount);
      j.key("stationCount"); j.value((unsigned long long)s.stationCount);
      j.key("factionId"); j.value((unsigned long long)s.factionId);
      if (emitStubSig) {
        j.key("stubSig");
        j.value((unsigned long long)sim::signatureSystemStub(s));
      }
      const sim::StarSystem* sysPtr = nullptr;
      if (emitSysSig || doSecurity) {
        sysPtr = &u.getSystem(s.id, &s);
      }
      if (emitSysSig && sysPtr) {
        j.key("sysSig");
        j.value((unsigned long long)sim::signatureStarSystem(*sysPtr));
      }
      if (doSecurity && sysPtr) {
        const auto sp = sim::systemSecurityProfile(seed, *sysPtr);

        j.key("security");
        j.beginObject();
        j.key("controllingFactionId"); j.value((unsigned long long)sp.controllingFactionId);
        j.key("controlFrac"); j.value(sp.controlFrac);
        j.key("contest01"); j.value(sp.contest01);
        j.key("security01"); j.value(sp.security01);
        j.key("piracy01"); j.value(sp.piracy01);
        j.key("traffic01"); j.value(sp.traffic01);

        j.key("traits");
        j.beginObject();
        j.key("authority"); j.value(sp.traits.authority);
        j.key("corruption"); j.value(sp.traits.corruption);
        j.key("wealth"); j.value(sp.traits.wealth);
        j.key("stability"); j.value(sp.traits.stability);
        j.key("tech"); j.value(sp.traits.tech);
        j.key("militarism"); j.value(sp.traits.militarism);
        j.endObject();

        j.endObject();
      }
      j.endObject();
    }
    j.endArray();

    if (doFactions) {
      const auto& facs = u.factions();

      j.key("factions");
      j.beginArray();
      for (const auto& f : facs) {
        const auto p = sim::factionProfile(seed, f.id);

        j.beginObject();
        j.key("id"); j.value((unsigned long long)f.id);
        j.key("name"); j.value(f.name);
        j.key("taxRate"); j.value(f.taxRate);
        j.key("industryBias"); j.value(f.industryBias);
        j.key("influenceRadiusLy"); j.value(f.influenceRadiusLy);
        j.key("homePosLy");
        j.beginArray();
        j.value(f.homePosLy.x); j.value(f.homePosLy.y); j.value(f.homePosLy.z);
        j.endArray();

        j.key("profile");
        j.beginObject();
        j.key("authority"); j.value(p.authority);
        j.key("corruption"); j.value(p.corruption);
        j.key("wealth"); j.value(p.wealth);
        j.key("stability"); j.value(p.stability);
        j.key("tech"); j.value(p.tech);
        j.key("militarism"); j.value(p.militarism);
        j.endObject();

        if (doDiplomacy) {
          struct Rel { core::u32 id; double score; };
          std::vector<Rel> rels;
          rels.reserve(facs.size());
          for (const auto& o : facs) {
            if (o.id == f.id) continue;
            rels.push_back(Rel{o.id, sim::factionRelation(seed, f, o)});
          }
          std::sort(rels.begin(), rels.end(), [](const Rel& a, const Rel& b) { return a.score > b.score; });

          const auto writeTop = [&](const char* key, bool best) {
            j.key(key);
            j.beginArray();
            const std::size_t n = std::min<std::size_t>(3, rels.size());
            for (std::size_t i = 0; i < n; ++i) {
              const auto& r = best ? rels[i] : rels[rels.size() - 1 - i];
              j.beginObject();
              j.key("factionId"); j.value((unsigned long long)r.id);
              j.key("score"); j.value(r.score);
              const auto k = sim::classifyFactionRelation(r.score);
              j.key("kind"); j.value(sim::factionRelationKindName(k));
              j.endObject();
            }
            j.endArray();
          };

          writeTop("topAllies", true);
          writeTop("topHostiles", false);
        }

        j.endObject();
      }
      j.endArray();
    }
  } else {
    std::cout << "Seed: " << seed << "\n";
    std::cout << "Query @ (" << posLy.x << "," << posLy.y << "," << posLy.z << ") radius=" << radiusLy << " ly\n";
    std::cout << "Found " << systems.size() << " systems\n\n";

    if (doFactions) {
      const auto& facs = u.factions();
      std::cout << "Factions: " << facs.size() << "\n";
      for (const auto& f : facs) {
        const auto p = sim::factionProfile(seed, f.id);
        std::cout << "  [" << std::setw(2) << f.id << "] " << std::setw(18) << f.name
                  << "  tax=" << std::fixed << std::setprecision(3) << f.taxRate
                  << "  bias=" << std::fixed << std::setprecision(2) << f.industryBias
                  << "  infR=" << std::fixed << std::setprecision(0) << f.influenceRadiusLy
                  << "  home=(" << std::fixed << std::setprecision(0)
                  << f.homePosLy.x << "," << f.homePosLy.y << "," << f.homePosLy.z << ")"
                  << "\n";
        std::cout << "       profile: auth=" << std::fixed << std::setprecision(2) << p.authority
                  << "  corr=" << p.corruption
                  << "  wealth=" << p.wealth
                  << "  stab=" << p.stability
                  << "  tech=" << p.tech
                  << "  mil=" << p.militarism
                  << "\n";

        if (doDiplomacy) {
          struct Rel { core::u32 id; double score; };
          std::vector<Rel> rels;
          rels.reserve(facs.size());
          for (const auto& o : facs) {
            if (o.id == f.id) continue;
            rels.push_back(Rel{o.id, sim::factionRelation(seed, f, o)});
          }
          std::sort(rels.begin(), rels.end(), [](const Rel& a, const Rel& b) { return a.score > b.score; });

          const auto printTop = [&](const char* label, bool best) {
            const std::size_t n = std::min<std::size_t>(3, rels.size());
            std::cout << "       " << label << ":";
            for (std::size_t i = 0; i < n; ++i) {
              const auto& r = best ? rels[i] : rels[rels.size() - 1 - i];
              const auto k = sim::classifyFactionRelation(r.score);
              std::cout << "  " << r.id << "(" << std::fixed << std::setprecision(2) << r.score
                        << "," << sim::factionRelationKindName(k) << ")";
            }
            std::cout << "\n";
          };

          printTop("allies", true);
          printTop("hostiles", false);
        }
      }
      std::cout << "\n";
    }

    for (const auto& s : systems) {
      const math::Vec3d d = s.posLy - posLy;
      const double dist = std::sqrt(d.lengthSq());

      std::cout << std::setw(14) << s.id
                << "  " << std::setw(14) << s.name
                << "  class=" << starClassName(s.primaryClass)
                << "  dist=" << std::fixed << std::setprecision(2) << dist << " ly"
                << "  planets=" << s.planetCount
                << "  stations=" << s.stationCount
                << "  faction=" << s.factionId;

      if (emitStubSig) {
        const auto sig = sim::signatureSystemStub(s);
        std::cout << "  stubSig=" << (unsigned long long)sig;
      }
      const sim::StarSystem* sysPtr = nullptr;
      if (emitSysSig || doSecurity) {
        sysPtr = &u.getSystem(s.id, &s);
      }
      if (emitSysSig && sysPtr) {
        const auto sig = sim::signatureStarSystem(*sysPtr);
        std::cout << "  sysSig=" << (unsigned long long)sig;
      }
      if (doSecurity && sysPtr) {
        const auto sp = sim::systemSecurityProfile(seed, *sysPtr);
        std::cout << "  ctrl=" << sp.controllingFactionId
                  << "  sec=" << std::fixed << std::setprecision(2) << sp.security01
                  << "  piracy=" << sp.piracy01
                  << "  traffic=" << sp.traffic01
                  << "  contest=" << sp.contest01;
      }
      std::cout << "\n";
    }
  }

  if (!systems.empty() && !json) {
    std::cout << "\n--- Example system detail: " << systems.front().name << " ---\n";
    const auto& sys = u.getSystem(systems.front().id, &systems.front());
    std::cout << "Star mass=" << sys.star.massSol << " Sol, lum=" << sys.star.luminositySol << " Sol\n";
    std::cout << "Planets:\n";
    for (const auto& p : sys.planets) {
      std::cout << "  - " << p.name
                << " a=" << p.orbit.semiMajorAxisAU << " AU"
                << " e=" << p.orbit.eccentricity
                << " period=" << p.orbit.periodDays << " days\n";
    }
    std::cout << "Stations:\n";
    for (const auto& st : sys.stations) {
      std::cout << "  - " << st.name
                << " faction=" << st.factionId
                << " fee=" << st.feeRate
                << " illegal=["
                << sim::illegalCommodityListStringForStation(u.seed(), st.factionId, st.id, st.type)
                << "]\n";
    }
  }


  if (doSignals && !systems.empty()) {
    fromSysIdx = std::min(fromSysIdx, systems.size() - 1);
    const auto& stub = systems[fromSysIdx];
    const auto& sys = u.getSystem(stub.id, &stub);

    sim::SaveGame save{};
    if (!loadPath.empty()) {
      if (!sim::loadFromFile(loadPath, save)) {
        std::cerr << "Failed to load save from: " << loadPath << "\n";
        return 1;
      }
    }

    // Tool overrides keep output deterministic for the requested inputs.
    save.seed = seed;
    save.timeDays = timeDays;

    sim::SignalGenParams p{};
    p.resourceFieldCount = std::max(0, signalFieldCount);
    p.includeDailyDerelict = !signalNoDerelict;
    p.includeDistress = !signalNoDistress;
    p.distressPerDay = std::max(0, signalDistressPerDay);
    p.distressTtlDays = signalDistressTtlDays;

    const auto plan = sim::generateSystemSignals(seed, sys, timeDays, save.missions, save.resolvedSignalIds, p);

    if (json) {
      j.key("signals");
      j.beginObject();
      j.key("system");
      j.beginObject();
      j.key("systemId"); j.value((unsigned long long)stub.id);
      j.key("systemName"); j.value(stub.name);
      j.key("day"); j.value(timeDays);
      j.endObject();

      j.key("params");
      j.beginObject();
      j.key("resourceFieldCount"); j.value((unsigned long long)p.resourceFieldCount);
      j.key("includeDailyDerelict"); j.value(p.includeDailyDerelict);
      j.key("includeDistress"); j.value(p.includeDistress);
      j.key("distressPerDay"); j.value((unsigned long long)p.distressPerDay);
      j.key("distressTtlDays"); j.value(p.distressTtlDays);
      j.endObject();

      j.key("resourceFields");
      j.beginObject();
      j.key("fieldCount"); j.value((unsigned long long)plan.resourceFields.fields.size());
      j.key("asteroidCount"); j.value((unsigned long long)plan.resourceFields.asteroids.size());
      j.endObject();

      j.key("sites");
      j.beginArray();
      for (const auto& s : plan.sites) {
        j.beginObject();
        j.key("id"); j.value((unsigned long long)s.id);
        j.key("kind"); j.value(sim::signalKindName(s.kind));
        j.key("posKm");
        j.beginArray();
        j.value(s.posKm.x); j.value(s.posKm.y); j.value(s.posKm.z);
        j.endArray();
        j.key("expireDay"); j.value(s.expireDay);
        j.key("resolved"); j.value(s.resolved);
        if (s.kind == sim::SignalKind::ResourceField) {
          j.key("fieldKind"); j.value(sim::resourceFieldKindName(s.fieldKind));
        }
        if (s.kind == sim::SignalKind::Distress && s.hasDistressPlan) {
          j.key("distress");
          j.beginObject();
          j.key("scenario");
          switch (s.distress.scenario) {
            case sim::DistressScenario::Supplies: j.value("Supplies"); break;
            case sim::DistressScenario::Fuel: j.value("Fuel"); break;
            case sim::DistressScenario::Medical: j.value("Medical"); break;
            case sim::DistressScenario::Mechanical: j.value("Mechanical"); break;
            case sim::DistressScenario::Ambush: j.value("Ambush"); break;
            default: j.value("Unknown"); break;
          }
          j.key("needCommodity"); j.value(econ::commodityDef(s.distress.needCommodity).code);
          j.key("needUnits"); j.value(s.distress.needUnits);
          j.key("rewardCr"); j.value(s.distress.rewardCr);
          j.key("repReward"); j.value(s.distress.repReward);
          j.key("ambush"); j.value(s.distress.ambush);
          j.key("pirateCount"); j.value((unsigned long long)s.distress.pirateCount);
          j.key("risk"); j.value(s.distress.risk);
          j.endObject();
        }
        if (s.kind == sim::SignalKind::MissionSalvage) {
          j.key("missionId"); j.value((unsigned long long)s.missionId);
        }
        j.endObject();
      }
      j.endArray();

      j.endObject();
    } else {
      std::cout << "\n--- Signals @ " << stub.name << " (day=" << timeDays << ") ---\n";
      std::cout << "Resource fields: " << plan.resourceFields.fields.size()
                << " (asteroids=" << plan.resourceFields.asteroids.size() << ")\n";
      for (const auto& s : plan.sites) {
        std::cout << "  - [" << sim::signalKindName(s.kind) << "] id=" << (unsigned long long)s.id;
        std::cout << " posKm=(" << std::fixed << std::setprecision(1)
                  << s.posKm.x << "," << s.posKm.y << "," << s.posKm.z << ")";
        if (s.expireDay > 0.0) {
          std::cout << " expires@" << std::fixed << std::setprecision(2) << s.expireDay;
        }
        if (s.resolved) std::cout << " [resolved]";
        if (s.kind == sim::SignalKind::ResourceField) {
          std::cout << " kind=" << sim::resourceFieldKindName(s.fieldKind);
        }
        if (s.kind == sim::SignalKind::Distress && s.hasDistressPlan) {
          std::cout << " distress=";
          switch (s.distress.scenario) {
            case sim::DistressScenario::Supplies: std::cout << "Supplies"; break;
            case sim::DistressScenario::Fuel: std::cout << "Fuel"; break;
            case sim::DistressScenario::Medical: std::cout << "Medical"; break;
            case sim::DistressScenario::Mechanical: std::cout << "Mechanical"; break;
            case sim::DistressScenario::Ambush: std::cout << "Ambush"; break;
            default: std::cout << "Unknown"; break;
          }
          std::cout << " need=" << econ::commodityDef(s.distress.needCommodity).code << "x" << std::fixed << std::setprecision(0) << s.distress.needUnits;
          std::cout << " reward=" << std::fixed << std::setprecision(0) << s.distress.rewardCr;
          if (s.distress.ambush) std::cout << " ambush(pirates=" << s.distress.pirateCount << ")";
        }
        if (s.kind == sim::SignalKind::MissionSalvage) {
          std::cout << " missionId=" << (unsigned long long)s.missionId;
        }
        std::cout << "\n";
      }
    }
  }


  if (doLaw && !systems.empty()) {
    // Clamp indices so the tool remains easy to use from scripts.
    fromSysIdx = std::min(fromSysIdx, systems.size() - 1);

    const auto& fromStub = systems[fromSysIdx];
    const auto& fromSys = u.getSystem(fromStub.id, &fromStub);
    if (fromSys.stations.empty()) {
      std::cerr << "\n[law] system has no stations\n";
      return 1;
    }

    fromStationIdx = std::min(fromStationIdx, fromSys.stations.size() - 1);
    const auto& st = fromSys.stations[fromStationIdx];

    const core::u32 factionId = (lawFactionOverride != 0) ? lawFactionOverride : st.factionId;
    const sim::LawProfile law = sim::lawProfile(seed, factionId);
    const std::string illegalList = sim::illegalCommodityListString(seed, factionId);

    const double fineCr = law.contrabandFineCr(lawIllegalValueCr);
    const double bribeChance = sim::bribeOfferChance(law, lawRep, lawHeat, lawIllegalValueCr);
    const double bribeCr = sim::bribeAmountCrRounded(law, lawIllegalValueCr, 10.0);

    const double rateClean = sim::cargoScanStartRatePerSec(false, law, smuggleHoldMk);
    const double rateContraband = sim::cargoScanStartRatePerSec(true, law, smuggleHoldMk);

    const double durStationClean = sim::cargoScanDurationSecStation(false, smuggleHoldMk);
    const double durStationContraband = sim::cargoScanDurationSecStation(true, smuggleHoldMk);
    const double durPolice = sim::cargoScanDurationSecPolice(smuggleHoldMk);

    if (json) {
      j.key("law");
      j.beginObject();
      j.key("systemId"); j.value((unsigned long long)fromStub.id);
      j.key("stationId"); j.value((unsigned long long)st.id);
      j.key("systemName"); j.value(fromStub.name);
      j.key("stationName"); j.value(st.name);
      j.key("factionId"); j.value((unsigned long long)factionId);
      j.key("illegal"); j.value(illegalList);

      j.key("profile");
      j.beginObject();
      j.key("scanStrictness"); j.value(law.scanStrictness);
      j.key("corruption"); j.value(law.corruption);
      j.key("fineBaseCr"); j.value(law.fineBaseCr);
      j.key("fineRate"); j.value(law.fineRate);
      j.key("repBase"); j.value(law.repBase);
      j.key("repDiv"); j.value(law.repDiv);
      j.key("repMin"); j.value(law.repMin);
      j.key("repMax"); j.value(law.repMax);
      j.key("evadeRepMult"); j.value(law.evadeRepMult);
      j.endObject();

      j.key("smuggleHoldMk"); j.value((int)smuggleHoldMk);
      j.key("scan");
      j.beginObject();
      j.key("ratePerSecClean"); j.value(rateClean);
      j.key("ratePerSecContraband"); j.value(rateContraband);
      j.key("durationStationSecClean"); j.value(durStationClean);
      j.key("durationStationSecContraband"); j.value(durStationContraband);
      j.key("durationPoliceSec"); j.value(durPolice);
      j.endObject();

      j.key("example");
      j.beginObject();
      j.key("illegalValueCr"); j.value(lawIllegalValueCr);
      j.key("fineCr"); j.value(fineCr);
      j.key("bribeChance"); j.value(bribeChance);
      j.key("bribeCr"); j.value(bribeCr);
      j.key("rep"); j.value(lawRep);
      j.key("heat"); j.value(lawHeat);
      j.endObject();

      j.endObject();
    } else {
      std::cout << "\n--- Law / contraband profile @ " << fromStub.name << " / " << st.name << " ---\n";
      std::cout << "FactionId=" << factionId << "\n";
      std::cout << "Illegal: " << illegalList << "\n";
      std::cout << "Profile: scanStrictness=" << std::fixed << std::setprecision(2) << law.scanStrictness
                << " corruption=" << std::fixed << std::setprecision(2) << law.corruption
                << " fineBase=" << std::fixed << std::setprecision(0) << law.fineBaseCr
                << " fineRate=" << std::fixed << std::setprecision(2) << law.fineRate
                << " repBase=" << law.repBase
                << " repDiv=" << law.repDiv
                << " repMin=" << law.repMin
                << " repMax=" << law.repMax
                << " evadeMult=" << law.evadeRepMult
                << "\n";
      std::cout << "Scan: holdMk=" << smuggleHoldMk
                << " rateClean=" << std::fixed << std::setprecision(5) << rateClean << "/s"
                << " rateContraband=" << std::fixed << std::setprecision(5) << rateContraband << "/s"
                << " durStationClean=" << std::fixed << std::setprecision(2) << durStationClean << "s"
                << " durStationContraband=" << std::fixed << std::setprecision(2) << durStationContraband << "s"
                << " durPolice=" << std::fixed << std::setprecision(2) << durPolice << "s"
                << "\n";
      std::cout << "Example: illegalValue=" << std::fixed << std::setprecision(1) << lawIllegalValueCr
                << " fine=" << std::fixed << std::setprecision(1) << fineCr
                << " bribeChance=" << std::fixed << std::setprecision(3) << bribeChance
                << " bribe=" << std::fixed << std::setprecision(1) << bribeCr
                << " (rep=" << lawRep << " heat=" << lawHeat << ")\n";
    }
  }

  if (doWarehouse && !systems.empty()) {
    // Clamp indices so the tool remains easy to use from scripts.
    fromSysIdx = std::min(fromSysIdx, systems.size() - 1);

    const auto& fromStub = systems[fromSysIdx];
    const auto& fromSys = u.getSystem(fromStub.id, &fromStub);
    if (fromSys.stations.empty()) {
      std::cerr << "\n[warehouse] system has no stations\n";
      return 1;
    }

    fromStationIdx = std::min(fromStationIdx, fromSys.stations.size() - 1);
    const auto& st = fromSys.stations[fromStationIdx];
    const auto stType = st.economyModel.type;

    sim::StationStorage entry{};
    entry.stationId = st.id;
    entry.stationType = stType;
    entry.factionId = st.factionId;
    entry.feeRate = st.feeRate;
    entry.lastFeeDay = timeDays;
    entry.feesDueCr = 0.0;
    entry.cargo.fill(0.0);
    if (storeUnits > 0.0) {
      entry.cargo[(int)storeCommodity] = storeUnits;
    }

    const double massKg = sim::storageMassKg(entry);
    const double rate = sim::storageFeeRateCrPerKgPerDay(entry, rep);
    const double daily = sim::estimateStorageDailyFeeCr(entry, rep);

    sim::StationStorage future = entry;
    sim::accrueStorageFees(future, timeDays + std::max(0.0, advanceDays), rep);

    if (json) {
      j.key("warehouse");
      j.beginObject();
      j.key("systemId"); j.value((unsigned long long)fromStub.id);
      j.key("stationId"); j.value((unsigned long long)st.id);
      j.key("systemName"); j.value(fromStub.name);
      j.key("stationName"); j.value(st.name);
      j.key("stationType"); j.value(stationTypeName(stType));
      j.key("rep"); j.value(rep);
      j.key("commodity"); j.value(econ::commodityCode(storeCommodity));
      j.key("units"); j.value(storeUnits);
      j.key("massKg"); j.value(massKg);
      j.key("feeRateCrPerKgPerDay"); j.value(rate);
      j.key("dailyFeeCr"); j.value(daily);
      j.key("projectDays"); j.value(std::max(0.0, advanceDays));
      j.key("projectFeesDueCr"); j.value(future.feesDueCr);
      j.endObject();
    } else {
      std::cout << "\n--- Warehouse quote @ " << fromStub.name << " / " << st.name
                << " (" << stationTypeName(stType) << ") ---\n";
      std::cout << "Store: " << econ::commodityCode(storeCommodity)
                << " (" << econ::commodityName(storeCommodity) << ")"
                << " units=" << std::fixed << std::setprecision(2) << storeUnits
                << "  mass=" << std::fixed << std::setprecision(1) << massKg << " kg\n";
      std::cout << "Fee: rate=" << std::fixed << std::setprecision(5) << rate << " cr/kg/day"
                << "  estDaily=" << std::fixed << std::setprecision(1) << daily << " cr/day";
      if (advanceDays > 0.0) {
        std::cout << "  projected(" << std::fixed << std::setprecision(2) << advanceDays << "d)="
                  << std::fixed << std::setprecision(1) << future.feesDueCr << " cr";
      }
      std::cout << "\n";
    }
  }

  if (doIndustry && !systems.empty()) {
    // Clamp indices so the tool remains easy to use from scripts.
    fromSysIdx = std::min(fromSysIdx, systems.size() - 1);

    const auto& fromStub = systems[fromSysIdx];
    const auto& fromSys = u.getSystem(fromStub.id, &fromStub);
    if (fromSys.stations.empty()) {
      std::cerr << "\n[industry] system has no stations\n";
      return 1;
    }

    fromStationIdx = std::min(fromStationIdx, fromSys.stations.size() - 1);
    const auto& st = fromSys.stations[fromStationIdx];
    const auto stType = st.economyModel.type;
    const auto recipes = sim::availableIndustryRecipes(stType);

    if (json) {
      j.key("industry");
      j.beginObject();
      j.key("systemId"); j.value((unsigned long long)fromStub.id);
      j.key("stationId"); j.value((unsigned long long)st.id);
      j.key("systemName"); j.value(fromStub.name);
      j.key("stationName"); j.value(st.name);
      j.key("stationType"); j.value(stationTypeName(stType));
      j.key("batches"); j.value(industryBatches);
      j.key("rep"); j.value(rep);

      j.key("recipes");
      j.beginArray();
      for (const auto* r : recipes) {
        const auto q = sim::quoteIndustryOrder(*r, st.id, stType, industryBatches, st.feeRate, rep);
        j.beginObject();
        j.key("code"); j.value(r->code);
        j.key("name"); j.value(r->name);
        j.key("desc"); j.value(r->desc);

        j.key("inputA"); j.value(econ::commodityCode(q.inputA));
        j.key("inputAUnits"); j.value(q.inputAUnits);
        if (q.inputBUnits > 0.0) {
          j.key("inputB"); j.value(econ::commodityCode(q.inputB));
          j.key("inputBUnits"); j.value(q.inputBUnits);
        }
        j.key("output"); j.value(econ::commodityCode(q.output));
        j.key("outputUnits"); j.value(q.outputUnits);
        j.key("timeDays"); j.value(q.timeDays);
        j.key("feeCr"); j.value(q.serviceFeeCr);

        j.key("mods");
        j.beginObject();
        j.key("yieldMul"); j.value(q.mods.yieldMul);
        j.key("speedMul"); j.value(q.mods.speedMul);
        j.key("feeMul"); j.value(q.mods.feeMul);
        j.endObject();

        j.endObject();
      }
      j.endArray();
      j.endObject();
    } else {
      std::cout << "\n--- Industry @ " << fromStub.name << " / " << st.name
                << " (" << stationTypeName(stType) << ") ---\n";
      if (recipes.empty()) {
        std::cout << "No industry facilities here.\n";
      } else {
        std::cout << "Batches: " << std::fixed << std::setprecision(2) << industryBatches
                  << "  rep: " << std::fixed << std::setprecision(1) << rep
                  << "  feeRate: " << std::fixed << std::setprecision(3) << st.feeRate
                  << "\n\n";
        for (const auto* r : recipes) {
          const auto q = sim::quoteIndustryOrder(*r, st.id, stType, industryBatches, st.feeRate, rep);
          std::cout << "[" << r->code << "] " << r->name << "\n";
          std::cout << "  " << r->desc << "\n";
          std::cout << "  in:  " << econ::commodityCode(q.inputA) << " " << std::fixed << std::setprecision(2) << q.inputAUnits;
          if (q.inputBUnits > 0.0) {
            std::cout << " + " << econ::commodityCode(q.inputB) << " " << std::fixed << std::setprecision(2) << q.inputBUnits;
          }
          std::cout << "\n";
          std::cout << "  out: " << econ::commodityCode(q.output) << " " << std::fixed << std::setprecision(2) << q.outputUnits << "\n";
          std::cout << "  time: " << std::fixed << std::setprecision(2) << (q.timeDays * 24.0) << " h"
                    << "  fee: " << std::fixed << std::setprecision(0) << q.serviceFeeCr << " cr\n";
          std::cout << "  mods: yield=" << std::fixed << std::setprecision(3) << q.mods.yieldMul
                    << " speed=" << std::fixed << std::setprecision(3) << q.mods.speedMul
                    << " fee=" << std::fixed << std::setprecision(3) << q.mods.feeMul
                    << "\n\n";
        }
      }
    }
  }

  if (doIndustryTrade && !systems.empty()) {
    // Clamp indices so the tool remains easy to use from scripts.
    fromSysIdx = std::min(fromSysIdx, systems.size() - 1);

    const auto& fromStub = systems[fromSysIdx];
    const auto& fromSys = u.getSystem(fromStub.id, &fromStub);
    if (fromSys.stations.empty()) {
      std::cerr << "\n[industryTrade] system has no stations\n";
      return 1;
    }

    fromStationIdx = std::min(fromStationIdx, fromSys.stations.size() - 1);
    const auto& fromSt = fromSys.stations[fromStationIdx];

    const auto feeEff = [&](const sim::Station& st) {
      // Tooling mode: use a single rep value across all factions.
      return sim::applyReputationToFeeRate(st.feeRate, rep);
    };

    sim::IndustryTradeScanParams scan;
    scan.maxResults = tradeLimit;
    scan.perStationLimit = 1;
    scan.cargoCapacityKg = cargoCapacityKg;
    scan.cargoUsedKg = 0.0;
    scan.useFreeHold = true;
    scan.bidAskSpread = 0.10;
    scan.processingRep = rep;
    scan.minNetProfit = minProfit;
    scan.includeSameSystem = true;

    const auto ideas = sim::scanIndustryTradeOpportunities(u, fromStub, fromSt, timeDays, systems, scan, feeEff);

    if (json) {
      j.key("industryTrade");
      j.beginObject();
      j.key("fromSystemId"); j.value((unsigned long long)fromStub.id);
      j.key("fromStationId"); j.value((unsigned long long)fromSt.id);
      j.key("fromSystemName"); j.value(fromStub.name);
      j.key("fromStationName"); j.value(fromSt.name);
      j.key("fromStationType"); j.value(stationTypeName(fromSt.type));
      j.key("cargoKg"); j.value(cargoCapacityKg);
      j.key("rep"); j.value(rep);
      j.key("minProfit"); j.value(minProfit);

      j.key("opportunities");
      j.beginArray();
      for (const auto& t : ideas) {
        const auto* def = sim::findIndustryRecipe(t.recipe);
        j.beginObject();
        j.key("toSystemId"); j.value((unsigned long long)t.toSystem);
        j.key("toStationId"); j.value((unsigned long long)t.toStation);
        j.key("toSystemName"); j.value(t.toSystemName);
        j.key("toStationName"); j.value(t.toStationName);
        j.key("distanceLy"); j.value(t.distanceLy);
        j.key("recipe"); j.value(def ? def->code : "RECIPE");
        j.key("batches"); j.value(t.batches);
        j.key("inputA"); j.value(econ::commodityCode(t.inputA));
        j.key("inputAUnits"); j.value(t.inputAUnits);
        if (t.inputBUnits > 0.0) {
          j.key("inputB"); j.value(econ::commodityCode(t.inputB));
          j.key("inputBUnits"); j.value(t.inputBUnits);
        }
        j.key("output"); j.value(econ::commodityCode(t.output));
        j.key("outputUnits"); j.value(t.outputUnits);
        j.key("serviceFeeCr"); j.value(t.serviceFeeCr);
        j.key("timeDays"); j.value(t.timeDays);
        j.key("inputCostCr"); j.value(t.inputACostCr + t.inputBCostCr);
        j.key("outputRevenueCr"); j.value(t.outputRevenueCr);
        j.key("netProfitCr"); j.value(t.netProfitCr);
        j.key("netProfitPerKg"); j.value(t.netProfitPerKg);
        j.key("netProfitPerDay"); j.value(t.netProfitPerDay);
        j.endObject();
      }
      j.endArray();
      j.endObject();
    } else {
      std::cout << "\n--- Industry Trade @ " << fromStub.name << " / " << fromSt.name
                << " (" << stationTypeName(fromSt.type) << ") ---\n";
      std::cout << "Assumes buying recipe inputs at the origin station, then selling the output elsewhere.\n";
      if (ideas.empty()) {
        std::cout << "No profitable industrial routes found in the current query set.\n";
      } else {
        for (const auto& t : ideas) {
          const auto* def = sim::findIndustryRecipe(t.recipe);
          const char* code = def ? def->code : "RECIPE";

          std::cout << "[" << code << "] "
                    << "batches=" << std::fixed << std::setprecision(0) << t.batches
                    << "  out=" << econ::commodityCode(t.output) << " " << std::fixed << std::setprecision(1) << t.outputUnits
                    << "  profit=" << std::fixed << std::setprecision(0) << t.netProfitCr << " cr"
                    << " (" << std::fixed << std::setprecision(1) << t.netProfitPerKg << " cr/kg"
                    << ", " << std::fixed << std::setprecision(0) << t.netProfitPerDay << " cr/day)"
                    << "  -> " << t.toSystemName << " / " << t.toStationName
                    << "  dist=" << std::fixed << std::setprecision(1) << t.distanceLy << " ly\n";

          std::cout << "  in:  " << econ::commodityCode(t.inputA) << " " << std::fixed << std::setprecision(1) << t.inputAUnits
                    << " (ask=" << std::fixed << std::setprecision(1) << t.inputAAsk << ")";
          if (t.inputBUnits > 0.0) {
            std::cout << " + " << econ::commodityCode(t.inputB) << " " << std::fixed << std::setprecision(1) << t.inputBUnits
                      << " (ask=" << std::fixed << std::setprecision(1) << t.inputBAsk << ")";
          }
          std::cout << "\n";
          std::cout << "  fee: service=" << std::fixed << std::setprecision(0) << t.serviceFeeCr
                    << " cr  time=" << std::fixed << std::setprecision(1) << (t.timeDays * 24.0) << " h\n";
        }
      }
    }
  }

  if (doTrade && !systems.empty()) {
    // Clamp indices so the tool remains easy to use from scripts.
    fromSysIdx = std::min(fromSysIdx, systems.size() - 1);

    const auto& fromStub = systems[fromSysIdx];
    const auto& fromSys = u.getSystem(fromStub.id, &fromStub);
    if (fromSys.stations.empty()) {
      std::cerr << "\n[trade] origin system has no stations\n";
      return 1;
    }
    fromStationIdx = std::min(fromStationIdx, fromSys.stations.size() - 1);
    const auto& fromSt = fromSys.stations[fromStationIdx];

    const auto feeEff = [&](const sim::Station& st) {
      return sim::applyReputationToFeeRate(st.feeRate, rep);
    };

    sim::TradeScanParams scan;
    scan.maxResults = tradeLimit;
    scan.perStationLimit = 1;
    scan.cargoCapacityKg = cargoCapacityKg;
    scan.cargoUsedKg = 0.0;
    scan.useFreeHold = true;
    scan.bidAskSpread = 0.10;
    scan.minNetProfit = minProfit;
    scan.includeSameSystem = true;
    scan.commodityFilterEnabled = hasCommodityFilter;
    scan.commodityFilter = commodityFilterId;

    const auto ideas = sim::scanTradeOpportunities(u, fromStub, fromSt, timeDays, systems, scan, feeEff);

    if (json) {
      j.key("trade");
      j.beginObject();
      j.key("day"); j.value(timeDays);
      j.key("cargoKg"); j.value(cargoCapacityKg);
      j.key("commodityFilter");
      if (hasCommodityFilter) {
        j.value(econ::commodityCode(commodityFilterId));
      } else {
        j.nullValue();
      }
      j.key("minProfit"); j.value(minProfit);
      j.key("from");
      j.beginObject();
      j.key("systemId"); j.value((unsigned long long)fromStub.id);
      j.key("systemName"); j.value(fromStub.name);
      j.key("stationId"); j.value((unsigned long long)fromSt.id);
      j.key("stationName"); j.value(fromSt.name);
      j.key("feeRate"); j.value(fromSt.feeRate);
      j.key("feeRateEff"); j.value(feeEff(fromSt));
      j.endObject();
      j.key("ideas");
      j.beginArray();
      for (const auto& it : ideas) {
        j.beginObject();
        j.key("commodity"); j.value(econ::commodityCode(it.commodity));
        j.key("toSystemId"); j.value((unsigned long long)it.toSystem);
        j.key("toStationId"); j.value((unsigned long long)it.toStation);
        j.key("toSystemName"); j.value(it.toSystemName);
        j.key("toStationName"); j.value(it.toStationName);
        j.key("distLy"); j.value(it.distanceLy);
        j.key("unitsFrom"); j.value(it.unitsFrom);
        j.key("unitsToSpace"); j.value(it.unitsToSpace);
        j.key("unitsPossible"); j.value(it.unitsPossible);
        j.key("unitMassKg"); j.value(it.unitMassKg);
        j.key("feeFrom"); j.value(it.feeFrom);
        j.key("feeTo"); j.value(it.feeTo);
        j.key("netProfitPerUnit"); j.value(it.netProfitPerUnit);
        j.key("netTripProfit"); j.value(it.netProfitTotal);
        j.key("buyAsk"); j.value(it.buyPrice);
        j.key("sellBid"); j.value(it.sellPrice);
        j.endObject();
      }
      j.endArray();
      j.endObject();
    } else {
      std::cout << "\n--- Trade scan (day=" << timeDays << ", cargoKg=" << cargoCapacityKg << ") ---\n";
      std::cout << "From: " << fromStub.name << " / " << fromSt.name
                << "  (fee=" << fromSt.feeRate * 100.0
                << "%, eff=" << feeEff(fromSt) * 100.0 << "%)\n";
      if (hasCommodityFilter) {
        std::cout << "Filter: " << econ::commodityCode(commodityFilterId)
                  << " (" << econ::commodityName(commodityFilterId) << ")";
        if (minProfit > 0.0) std::cout << "  minProfit=" << std::fixed << std::setprecision(0) << minProfit;
        std::cout << "\n";
      } else if (minProfit > 0.0) {
        std::cout << "Min profit: " << std::fixed << std::setprecision(0) << minProfit << "\n";
      }
      if (ideas.empty()) {
        std::cout << "No profitable routes found inside the query set.\n";
      } else {
        for (const auto& it : ideas) {
          const auto code = econ::commodityCode(it.commodity);
          std::cout << std::setw(6) << code
                    << "  to=" << it.toSystemName << "/" << it.toStationName
                    << "  dist=" << std::fixed << std::setprecision(1) << it.distanceLy << " ly"
                    << "  units=" << std::fixed << std::setprecision(0) << it.unitsPossible
                    << "  net/unit=" << std::fixed << std::setprecision(2) << it.netProfitPerUnit
                    << "  net/trip=" << std::fixed << std::setprecision(0) << it.netProfitTotal
                    << "\n";
        }
        std::cout << "Tip: increase --radius (or move --pos) to scan a larger area.\n";
      }
    }
  }


  if (doTradeLoop && !systems.empty()) {
    // Clamp indices so the tool remains easy to use from scripts.
    fromSysIdx = std::min(fromSysIdx, systems.size() - 1);

    const auto& fromStub = systems[fromSysIdx];
    const auto& fromSys = u.getSystem(fromStub.id, &fromStub);
    if (fromSys.stations.empty()) {
      std::cerr << "\n[tradeLoop] origin system has no stations\n";
      return 1;
    }
    fromStationIdx = std::min(fromStationIdx, fromSys.stations.size() - 1);
    const auto& fromSt = fromSys.stations[fromStationIdx];

    const auto feeEff = [&](const sim::Station& st) {
      return sim::applyReputationToFeeRate(st.feeRate, rep);
    };

    sim::TradeLoopScanParams scan{};
    scan.manifest.cargoCapacityKg = cargoCapacityKg;
    scan.manifest.bidAskSpread = 0.10;
    scan.manifest.stepKg = mixStepKg;
    scan.manifest.maxBuyCreditsCr = mixMaxSpend;
    scan.manifest.simulatePriceImpact = !mixNoImpact;

    scan.minLegProfitCr = loopMinLegProfit;
    scan.minLoopProfitCr = loopMinProfit;
    scan.includeSameSystem = true;

    scan.legs = std::clamp<std::size_t>(loopLegs, 2, 3);
    scan.maxLegCandidates = std::max<std::size_t>(1, loopLegCandidates);
    scan.maxResults = std::max<std::size_t>(1, loopLimit);
    scan.maxStations = 256;

    const auto loops = jobs
      ? sim::scanTradeLoopsParallel(*jobs, u, fromStub, fromSt, timeDays, systems, scan, feeEff)
      : sim::scanTradeLoops(u, fromStub, fromSt, timeDays, systems, scan, feeEff);

    if (json) {
      j.key("tradeLoops");
      j.beginObject();
      j.key("day"); j.value(timeDays);
      j.key("cargoKg"); j.value(cargoCapacityKg);
      j.key("stepKg"); j.value(mixStepKg);
      j.key("simulatePriceImpact"); j.value(!mixNoImpact);
      j.key("maxBuyCreditsCr"); j.value(mixMaxSpend);

      j.key("legs"); j.value((unsigned long long)scan.legs);
      j.key("loopLimit"); j.value((unsigned long long)loopLimit);
      j.key("legCandidates"); j.value((unsigned long long)loopLegCandidates);
      j.key("minLegProfit"); j.value(loopMinLegProfit);
      j.key("minLoopProfit"); j.value(loopMinProfit);

      j.key("from");
      j.beginObject();
      j.key("systemId"); j.value((unsigned long long)fromStub.id);
      j.key("systemName"); j.value(fromStub.name);
      j.key("stationId"); j.value((unsigned long long)fromSt.id);
      j.key("stationName"); j.value(fromSt.name);
      j.key("feeRate"); j.value(fromSt.feeRate);
      j.key("feeRateEff"); j.value(feeEff(fromSt));
      j.endObject();

      j.key("loops");
      j.beginArray();
      for (const auto& lp : loops) {
        j.beginObject();
        j.key("totalProfitCr"); j.value(lp.totalProfitCr);
        j.key("totalDistanceLy"); j.value(lp.totalDistanceLy);
        j.key("profitPerLy"); j.value(lp.profitPerLy);

        j.key("legs");
        j.beginArray();
        for (const auto& lg : lp.legs) {
          j.beginObject();
          j.key("fromSystemId"); j.value((unsigned long long)lg.fromSystem);
          j.key("fromStationId"); j.value((unsigned long long)lg.fromStation);
          j.key("fromSystemName"); j.value(lg.fromSystemName);
          j.key("fromStationName"); j.value(lg.fromStationName);

          j.key("toSystemId"); j.value((unsigned long long)lg.toSystem);
          j.key("toStationId"); j.value((unsigned long long)lg.toStation);
          j.key("toSystemName"); j.value(lg.toSystemName);
          j.key("toStationName"); j.value(lg.toStationName);

          j.key("distanceLy"); j.value(lg.distanceLy);
          j.key("feeFrom"); j.value(lg.feeFrom);
          j.key("feeTo"); j.value(lg.feeTo);

          j.key("cargoFilledKg"); j.value(lg.manifest.cargoFilledKg);
          j.key("netBuyCr"); j.value(lg.manifest.netBuyCr);
          j.key("netSellCr"); j.value(lg.manifest.netSellCr);
          j.key("netProfitCr"); j.value(lg.manifest.netProfitCr);

          j.key("lines");
          j.beginArray();
          for (const auto& ln : lg.manifest.lines) {
            j.beginObject();
            j.key("commodity"); j.value(econ::commodityCode(ln.commodity));
            j.key("units"); j.value(ln.units);
            j.key("massKg"); j.value(ln.massKg);
            j.key("avgNetBuy"); j.value(ln.avgNetBuyPrice);
            j.key("avgNetSell"); j.value(ln.avgNetSellPrice);
            j.key("netProfit"); j.value(ln.netProfitCr);
            j.key("netProfitPerKg"); j.value(ln.netProfitPerKg);
            j.endObject();
          }
          j.endArray();

          j.endObject();
        }
        j.endArray();

        j.endObject();
      }
      j.endArray();
      j.endObject();
    } else {
      std::cout << "\n--- Trade loops (day=" << timeDays
                << ", cargoKg=" << cargoCapacityKg
                << ", legs=" << scan.legs
                << ", stepKg=" << mixStepKg
                << ", impact=" << (!mixNoImpact ? "on" : "off")
                << ") from " << fromStub.name << " / " << fromSt.name << " ---\n";

      if (loops.empty()) {
        std::cout << "No profitable loops found.\n";
        std::cout << "Tip: increase --radius (or move --pos) to scan a larger area.\n";
      } else {
        std::size_t rank = 0;
        for (const auto& lp : loops) {
          ++rank;
          std::cout << "\n[" << rank << "] totalProfit=" << std::fixed << std::setprecision(0) << lp.totalProfitCr
                    << " cr  dist=" << std::fixed << std::setprecision(2) << lp.totalDistanceLy
                    << " ly  profit/ly=" << std::fixed << std::setprecision(1) << lp.profitPerLy
                    << "\n";

          std::size_t legIdx = 0;
          for (const auto& lg : lp.legs) {
            ++legIdx;
            std::cout << "  " << legIdx << ") " << lg.fromSystemName << " / " << lg.fromStationName
                      << " -> " << lg.toSystemName << " / " << lg.toStationName
                      << "  dist=" << std::fixed << std::setprecision(2) << lg.distanceLy
                      << " ly  profit=" << std::fixed << std::setprecision(0) << lg.manifest.netProfitCr
                      << " cr\n";

            if (!lg.manifest.lines.empty()) {
              std::cout << "     cargo=" << std::fixed << std::setprecision(1) << lg.manifest.cargoFilledKg << " kg"
                        << "  buy=" << std::fixed << std::setprecision(0) << lg.manifest.netBuyCr << " cr"
                        << "  sell=" << std::fixed << std::setprecision(0) << lg.manifest.netSellCr << " cr\n";
              for (const auto& ln : lg.manifest.lines) {
                std::cout << "     - " << econ::commodityCode(ln.commodity)
                          << "  units=" << std::fixed << std::setprecision(0) << ln.units
                          << "  mass=" << std::fixed << std::setprecision(1) << ln.massKg << " kg"
                          << "  profit=" << std::fixed << std::setprecision(0) << ln.netProfitCr << " cr\n";
              }
            }
          }
        }
        std::cout << "\nTip: use --loopLegs 3 for triangular loops, or raise --loopLegCandidates for deeper search.\n";
      }
    }
  }


  if (doTradeRun && !systems.empty()) {
    // Clamp indices so the tool remains easy to use from scripts.
    fromSysIdx = std::min(fromSysIdx, systems.size() - 1);

    const auto& fromStub = systems[fromSysIdx];
    const auto& fromSys = u.getSystem(fromStub.id, &fromStub);
    if (fromSys.stations.empty()) {
      std::cerr << "\n[tradeRun] origin system has no stations\n";
      return 1;
    }
    fromStationIdx = std::min(fromStationIdx, fromSys.stations.size() - 1);
    const auto& fromSt = fromSys.stations[fromStationIdx];

    const auto feeEff = [&](const sim::Station& st) {
      return sim::applyReputationToFeeRate(st.feeRate, rep);
    };

    sim::TradeRunScanParams scan{};
    scan.manifest.cargoCapacityKg = cargoCapacityKg;
    scan.manifest.bidAskSpread = 0.10;
    scan.manifest.stepKg = mixStepKg;
    scan.manifest.maxBuyCreditsCr = mixMaxSpend;
    scan.manifest.simulatePriceImpact = !mixNoImpact;

    scan.minLegProfitCr = runMinLegProfit;
    scan.minRunProfitCr = runMinProfit;
    scan.includeSameSystem = true;
    scan.loopless = true;

    scan.legs = std::clamp<std::size_t>(runLegs, 1, 4);
    scan.beamWidth = std::max<std::size_t>(1, runBeam);
    scan.maxLegCandidates = std::max<std::size_t>(1, runLegCandidates);
    scan.maxResults = std::max<std::size_t>(1, runLimit);
    scan.maxStations = 256;

    // By default, do NOT constrain trade runs by jump-range unless --jr was explicitly provided.
    // (The sandbox's --jr has a default of 18 for the route tool, which would often yield no runs.)
    scan.jumpRangeLy = args.has("jr") ? jumpRangeLy : 0.0;

    // Reuse the route tool's cost params for ProfitPerCost ranking.
    scan.routeCostPerJump = fuelBase;
    scan.routeCostPerLy = fuelPerLy;

    {
      std::string rs = runScore;
      std::transform(rs.begin(), rs.end(), rs.begin(), [](unsigned char c) { return (char)std::tolower(c); });
      if (rs == "profitly" || rs == "profit/ly" || rs == "ly") {
        scan.scoreMode = sim::TradeRunScoreMode::ProfitPerLy;
      } else if (rs == "profithop" || rs == "profit/hop" || rs == "hop") {
        scan.scoreMode = sim::TradeRunScoreMode::ProfitPerHop;
      } else if (rs == "profitcost" || rs == "profit/cost" || rs == "cost") {
        scan.scoreMode = sim::TradeRunScoreMode::ProfitPerCost;
      } else {
        scan.scoreMode = sim::TradeRunScoreMode::TotalProfit;
      }
    }

    const auto runs = jobs ? sim::planTradeRunsParallel(*jobs, u, fromStub, fromSt, timeDays, systems, scan, feeEff)
                            : sim::planTradeRuns(u, fromStub, fromSt, timeDays, systems, scan, feeEff);

    if (json) {
      j.key("tradeRuns");
      j.beginObject();
      j.key("day"); j.value(timeDays);
      j.key("cargoKg"); j.value(cargoCapacityKg);
      j.key("stepKg"); j.value(mixStepKg);
      j.key("simulatePriceImpact"); j.value(!mixNoImpact);
      j.key("maxBuyCreditsCr"); j.value(mixMaxSpend);

      j.key("legs"); j.value((unsigned long long)scan.legs);
      j.key("runLimit"); j.value((unsigned long long)runLimit);
      j.key("beamWidth"); j.value((unsigned long long)runBeam);
      j.key("legCandidates"); j.value((unsigned long long)runLegCandidates);
      j.key("minLegProfit"); j.value(runMinLegProfit);
      j.key("minRunProfit"); j.value(runMinProfit);
      j.key("scoreMode"); j.value(runScore);
      j.key("jumpRangeLy"); j.value(scan.jumpRangeLy);
      j.key("routeCostPerJump"); j.value(scan.routeCostPerJump);
      j.key("routeCostPerLy"); j.value(scan.routeCostPerLy);

      j.key("from");
      j.beginObject();
      j.key("systemId"); j.value((unsigned long long)fromStub.id);
      j.key("systemName"); j.value(fromStub.name);
      j.key("stationId"); j.value((unsigned long long)fromSt.id);
      j.key("stationName"); j.value(fromSt.name);
      j.key("feeRate"); j.value(fromSt.feeRate);
      j.key("feeRateEff"); j.value(feeEff(fromSt));
      j.endObject();

      j.key("runs");
      j.beginArray();
      for (const auto& run : runs) {
        j.beginObject();
        j.key("totalProfitCr"); j.value(run.totalProfitCr);
        j.key("totalRouteDistanceLy"); j.value(run.totalRouteDistanceLy);
        j.key("totalRouteCost"); j.value(run.totalRouteCost);
        j.key("totalHops"); j.value((unsigned long long)run.totalHops);
        j.key("profitPerLy"); j.value(run.profitPerLy);
        j.key("profitPerHop"); j.value(run.profitPerHop);
        j.key("profitPerCost"); j.value(run.profitPerCost);

        j.key("legs");
        j.beginArray();
        for (const auto& lg : run.legs) {
          j.beginObject();
          j.key("fromSystemId"); j.value((unsigned long long)lg.fromSystem);
          j.key("fromStationId"); j.value((unsigned long long)lg.fromStation);
          j.key("fromSystemName"); j.value(lg.fromSystemName);
          j.key("fromStationName"); j.value(lg.fromStationName);
          j.key("toSystemId"); j.value((unsigned long long)lg.toSystem);
          j.key("toStationId"); j.value((unsigned long long)lg.toStation);
          j.key("toSystemName"); j.value(lg.toSystemName);
          j.key("toStationName"); j.value(lg.toStationName);

          j.key("routeHops"); j.value((unsigned long long)lg.routeHops);
          j.key("routeDistanceLy"); j.value(lg.routeDistanceLy);
          j.key("routeCost"); j.value(lg.routeCost);

          j.key("route");
          j.beginArray();
          for (const auto sysId : lg.route) {
            j.value((unsigned long long)sysId);
          }
          j.endArray();

          j.key("feeFrom"); j.value(lg.feeFrom);
          j.key("feeTo"); j.value(lg.feeTo);

          j.key("cargoFilledKg"); j.value(lg.manifest.cargoFilledKg);
          j.key("netBuyCr"); j.value(lg.manifest.netBuyCr);
          j.key("netSellCr"); j.value(lg.manifest.netSellCr);
          j.key("netProfitCr"); j.value(lg.manifest.netProfitCr);

          j.key("lines");
          j.beginArray();
          for (const auto& ln : lg.manifest.lines) {
            j.beginObject();
            j.key("commodity"); j.value(econ::commodityCode(ln.commodity));
            j.key("units"); j.value(ln.units);
            j.key("massKg"); j.value(ln.massKg);
            j.key("avgNetBuy"); j.value(ln.avgNetBuyPrice);
            j.key("avgNetSell"); j.value(ln.avgNetSellPrice);
            j.key("netProfit"); j.value(ln.netProfitCr);
            j.key("netProfitPerKg"); j.value(ln.netProfitPerKg);
            j.endObject();
          }
          j.endArray();

          j.endObject();
        }
        j.endArray();

        j.endObject();
      }
      j.endArray();
      j.endObject();
    } else {
      std::cout << "\n--- Trade runs (day=" << timeDays
                << ", cargoKg=" << cargoCapacityKg
                << ", legs=" << scan.legs
                << ", beam=" << scan.beamWidth
                << ", legCandidates=" << scan.maxLegCandidates
                << ", score=" << runScore
                << ", jr=" << (scan.jumpRangeLy > 0.0 ? scan.jumpRangeLy : 0.0)
                << ") from " << fromStub.name << " / " << fromSt.name << " ---\n";

      if (runs.empty()) {
        std::cout << "No profitable trade runs found.\n";
        std::cout << "Tip: increase --radius, raise --runBeam, or relax --runMinLegProfit/--runMinProfit.\n";
        std::cout << "Tip: add --jr <ly> to enforce jump-range reachability, or omit --jr to ignore it.\n";
      } else {
        std::size_t rank = 0;
        for (const auto& run : runs) {
          ++rank;
          std::cout << "\n[" << rank << "] totalProfit=" << std::fixed << std::setprecision(0) << run.totalProfitCr
                    << " cr  hops=" << run.totalHops
                    << "  dist=" << std::fixed << std::setprecision(2) << run.totalRouteDistanceLy
                    << " ly  cost=" << std::fixed << std::setprecision(2) << run.totalRouteCost
                    << "  profit/ly=" << std::fixed << std::setprecision(1) << run.profitPerLy
                    << "  profit/hop=" << std::fixed << std::setprecision(1) << run.profitPerHop
                    << "  profit/cost=" << std::fixed << std::setprecision(1) << run.profitPerCost
                    << "\n";

          std::size_t legIdx = 0;
          for (const auto& lg : run.legs) {
            ++legIdx;
            std::cout << "  " << legIdx << ") " << lg.fromSystemName << " / " << lg.fromStationName
                      << " -> " << lg.toSystemName << " / " << lg.toStationName
                      << "  hops=" << lg.routeHops
                      << "  dist=" << std::fixed << std::setprecision(2) << lg.routeDistanceLy << " ly"
                      << "  profit=" << std::fixed << std::setprecision(0) << lg.manifest.netProfitCr << " cr\n";

            if (!lg.manifest.lines.empty()) {
              std::cout << "     cargo=" << std::fixed << std::setprecision(1) << lg.manifest.cargoFilledKg << " kg"
                        << "  buy=" << std::fixed << std::setprecision(0) << lg.manifest.netBuyCr << " cr"
                        << "  sell=" << std::fixed << std::setprecision(0) << lg.manifest.netSellCr << " cr\n";

              const int take = std::min<int>(mixLines, (int)lg.manifest.lines.size());
              for (int i = 0; i < take; ++i) {
                const auto& ln = lg.manifest.lines[(std::size_t)i];
                std::cout << "     - " << econ::commodityCode(ln.commodity)
                          << "  units=" << std::fixed << std::setprecision(0) << ln.units
                          << "  mass=" << std::fixed << std::setprecision(1) << ln.massKg << " kg"
                          << "  profit=" << std::fixed << std::setprecision(0) << ln.netProfitCr << " cr\n";
              }
            }
          }
        }
        std::cout << "\nTip: use --runScore profitLy/profitHop/profitCost for more efficiency-oriented runs.\n";
        std::cout << "Tip: raise --runLegCandidates for deeper search (slower).\n";
      }
    }
  }


  if (doTradeMix && !systems.empty()) {
    // Clamp indices so the tool remains easy to use from scripts.
    fromSysIdx = std::min(fromSysIdx, systems.size() - 1);

    const auto& fromStub = systems[fromSysIdx];
    const auto& fromSys = u.getSystem(fromStub.id, &fromStub);
    if (fromSys.stations.empty()) {
      std::cerr << "\n[tradeMix] origin system has no stations\n";
      return 1;
    }
    fromStationIdx = std::min(fromStationIdx, fromSys.stations.size() - 1);
    const auto& fromSt = fromSys.stations[fromStationIdx];

    const auto feeEff = [&](const sim::Station& st) {
      return sim::applyReputationToFeeRate(st.feeRate, rep);
    };

    sim::TradeManifestScanParams scan;
    scan.maxResults = tradeLimit;
    scan.cargoCapacityKg = cargoCapacityKg;
    scan.cargoUsedKg = 0.0;
    scan.useFreeHold = true;
    scan.bidAskSpread = 0.10;
    scan.stepKg = mixStepKg;
    scan.maxBuyCreditsCr = mixMaxSpend;
    scan.simulatePriceImpact = !mixNoImpact;
    scan.minNetProfit = minProfit;
    scan.includeSameSystem = true;

    const auto ideas = sim::scanTradeManifests(u, fromStub, fromSt, timeDays, systems, scan, feeEff);

    if (json) {
      j.key("tradeMix");
      j.beginObject();
      j.key("day"); j.value(timeDays);
      j.key("cargoKg"); j.value(cargoCapacityKg);
      j.key("stepKg"); j.value(mixStepKg);
      j.key("simulatePriceImpact"); j.value(!mixNoImpact);
      j.key("maxBuyCreditsCr"); j.value(mixMaxSpend);
      j.key("minProfit"); j.value(minProfit);

      j.key("from");
      j.beginObject();
      j.key("systemId"); j.value((unsigned long long)fromStub.id);
      j.key("systemName"); j.value(fromStub.name);
      j.key("stationId"); j.value((unsigned long long)fromSt.id);
      j.key("stationName"); j.value(fromSt.name);
      j.key("feeRate"); j.value(fromSt.feeRate);
      j.key("feeRateEff"); j.value(feeEff(fromSt));
      j.endObject();

      j.key("ideas");
      j.beginArray();
      for (const auto& it : ideas) {
        j.beginObject();
        j.key("toSystemId"); j.value((unsigned long long)it.toSystem);
        j.key("toStationId"); j.value((unsigned long long)it.toStation);
        j.key("toSystemName"); j.value(it.toSystemName);
        j.key("toStationName"); j.value(it.toStationName);
        j.key("distLy"); j.value(it.distanceLy);

        j.key("cargoFilledKg"); j.value(it.cargoFilledKg);
        j.key("netBuyCr"); j.value(it.netBuyCr);
        j.key("netSellCr"); j.value(it.netSellCr);
        j.key("netTripProfit"); j.value(it.netProfitCr);

        j.key("lines");
        j.beginArray();
        for (const auto& ln : it.lines) {
          j.beginObject();
          j.key("commodity"); j.value(econ::commodityCode(ln.commodity));
          j.key("units"); j.value(ln.units);
          j.key("massKg"); j.value(ln.massKg);
          j.key("avgNetBuy"); j.value(ln.avgNetBuyPrice);
          j.key("avgNetSell"); j.value(ln.avgNetSellPrice);
          j.key("netProfit"); j.value(ln.netProfitCr);
          j.key("netProfitPerKg"); j.value(ln.netProfitPerKg);
          j.endObject();
        }
        j.endArray();

        j.endObject();
      }
      j.endArray();
      j.endObject();
    } else {
      std::cout << "\n--- Trade mix scan (day=" << timeDays
                << ", cargoKg=" << cargoCapacityKg
                << ", stepKg=" << mixStepKg
                << ", impact=" << (!mixNoImpact ? "yes" : "no")
                << ") ---\n";
      std::cout << "From: " << fromStub.name << " / " << fromSt.name
                << "  (fee=" << fromSt.feeRate * 100.0
                << "%, eff=" << feeEff(fromSt) * 100.0 << "%)\n";

      if (minProfit > 0.0) {
        std::cout << "Min profit: " << std::fixed << std::setprecision(0) << minProfit << "\n";
      }
      if (mixMaxSpend > 0.0) {
        std::cout << "Max spend: " << std::fixed << std::setprecision(0) << mixMaxSpend << " cr\n";
      }

      if (ideas.empty()) {
        std::cout << "No profitable manifests found inside the query set.\n";
      } else {
        for (const auto& it : ideas) {
          std::cout << "to=" << it.toSystemName << "/" << it.toStationName
                    << "  dist=" << std::fixed << std::setprecision(1) << it.distanceLy << " ly"
                    << "  fill=" << std::fixed << std::setprecision(1) << it.cargoFilledKg << " kg"
                    << "  net/trip=" << std::fixed << std::setprecision(0) << it.netProfitCr
                    << "\n";

          const int take = std::min<int>(mixLines, (int)it.lines.size());
          for (int i = 0; i < take; ++i) {
            const auto& ln = it.lines[(std::size_t)i];
            std::cout << "  - " << std::setw(6) << econ::commodityCode(ln.commodity)
                      << "  units=" << std::fixed << std::setprecision(1) << ln.units
                      << "  mass=" << std::fixed << std::setprecision(1) << ln.massKg << " kg"
                      << "  net/kg=" << std::fixed << std::setprecision(1) << ln.netProfitPerKg
                      << "  net=" << std::fixed << std::setprecision(0) << ln.netProfitCr
                      << "\n";
          }
        }
        std::cout << "Tip: decrease --mixStepKg for a more precise mix (slower), increase for speed.\n";
      }
    }
  }


  if (doRoute && !systems.empty()) {
    fromSysIdx = std::min(fromSysIdx, systems.size() - 1);
    toSysIdx = std::min(toSysIdx, systems.size() - 1);

    const auto& fromStub = systems[fromSysIdx];
    const auto& toStub = systems[toSysIdx];

    // Resolve cost model.
    std::string costModel = "hops";
    double costPerJump = 1.0;
    double costPerLy = 0.0;
    {
      std::string rc = routeCost;
      std::transform(rc.begin(), rc.end(), rc.begin(), [](unsigned char c) { return (char)std::tolower(c); });
      if (rc == "dist" || rc == "distance") {
        costModel = "dist";
        costPerJump = 0.0;
        costPerLy = 1.0;
      } else if (rc == "fuel") {
        costModel = "fuel";
        costPerJump = fuelBase;
        costPerLy = fuelPerLy;
      } else {
        costModel = "hops";
        costPerJump = 1.0;
        costPerLy = 0.0;
      }
    }

    // Always solve the best route once (gives expansion diagnostics).
    sim::RoutePlanStats bestStats{};
    std::vector<sim::SystemId> bestRoute;
    if (costModel == "dist") {
      bestRoute = sim::plotRouteAStarCost(systems, fromStub.id, toStub.id, jumpRangeLy, 0.0, 1.0, &bestStats);
    } else if (costModel == "fuel") {
      bestRoute = sim::plotRouteAStarCost(systems, fromStub.id, toStub.id, jumpRangeLy, fuelBase, fuelPerLy, &bestStats);
    } else {
      bestRoute = sim::plotRouteAStarHops(systems, fromStub.id, toStub.id, jumpRangeLy, &bestStats);
    }

    // If requested, also compute alternative loopless routes (Yen's algorithm).
    std::vector<sim::KRoute> routes;
    if (kRoutes <= 1) {
      if (!bestRoute.empty()) {
        routes.push_back(sim::KRoute{bestRoute, bestStats.hops, bestStats.distanceLy, bestStats.cost});
      }
    } else {
      if (costModel == "hops") {
        routes = sim::plotKRoutesAStarHops(systems, fromStub.id, toStub.id, jumpRangeLy, kRoutes);
      } else {
        routes = sim::plotKRoutesAStarCost(systems, fromStub.id, toStub.id, jumpRangeLy, costPerJump, costPerLy, kRoutes);
      }

      // In the unlikely event the k-shortest planner returns nothing but the best route exists,
      // keep the tool usable.
      if (routes.empty() && !bestRoute.empty()) {
        routes.push_back(sim::KRoute{bestRoute, bestStats.hops, bestStats.distanceLy, bestStats.cost});
      }
    }

    const double bestFuelEstimate = bestRoute.empty() ? 0.0 : (bestStats.hops * fuelBase + bestStats.distanceLy * fuelPerLy);

    if (json) {
      j.key("route");
      j.beginObject();
      j.key("fromSysIdx"); j.value((unsigned long long)fromSysIdx);
      j.key("toSysIdx"); j.value((unsigned long long)toSysIdx);
      j.key("fromSystemId"); j.value((unsigned long long)fromStub.id);
      j.key("toSystemId"); j.value((unsigned long long)toStub.id);
      j.key("jumpRangeLy"); j.value(jumpRangeLy);
      j.key("costModel"); j.value(costModel);
      j.key("cost"); j.value(bestStats.cost);
      j.key("kRoutes"); j.value((unsigned long long)kRoutes);
      j.key("fuelBase"); j.value(fuelBase);
      j.key("fuelPerLy"); j.value(fuelPerLy);
      j.key("fuelEstimate"); j.value(bestFuelEstimate);
      j.key("found"); j.value(!bestRoute.empty());
      j.key("hops"); j.value((unsigned long long)(bestRoute.size() > 1 ? bestRoute.size() - 1 : 0));
      j.key("distanceLy"); j.value(bestStats.distanceLy);
      j.key("visited"); j.value((unsigned long long)bestStats.visited);
      j.key("expansions"); j.value((unsigned long long)bestStats.expansions);
      j.key("path");
      j.beginArray();
      for (const auto id : bestRoute) j.value((unsigned long long)id);
      j.endArray();

      // Ranked route list (best first). Each entry mirrors the chosen cost model.
      j.key("routes");
      j.beginArray();
      for (std::size_t ri = 0; ri < routes.size(); ++ri) {
        const auto& r = routes[ri];
        const double fuelEst = r.path.empty() ? 0.0 : (r.hops * fuelBase + r.distanceLy * fuelPerLy);

        j.beginObject();
        j.key("rank"); j.value((unsigned long long)ri);
        j.key("hops"); j.value((unsigned long long)r.hops);
        j.key("distanceLy"); j.value(r.distanceLy);
        j.key("cost"); j.value(r.cost);
        j.key("fuelEstimate"); j.value(fuelEst);
        j.key("path");
        j.beginArray();
        for (const auto id : r.path) j.value((unsigned long long)id);
        j.endArray();
        j.endObject();
      }
      j.endArray();

      j.endObject();
    } else {
      std::cout << "\n--- Route plan (A*, cost=" << costModel
                << ", jr=" << std::fixed << std::setprecision(1) << jumpRangeLy
                << " ly, k=" << kRoutes << ") ---\n";
      std::cout << "From: [" << fromSysIdx << "] " << fromStub.name << " (" << fromStub.id << ")\n";
      std::cout << "To:   [" << toSysIdx << "] " << toStub.name << " (" << toStub.id << ")\n";

      if (bestRoute.empty()) {
        std::cout << "No route found inside the queried node set.\n";
        std::cout << "Tip: increase --radius (and/or --limit) so intermediate systems are available.\n";
      } else {
        std::unordered_map<sim::SystemId, const sim::SystemStub*> stubById;
        stubById.reserve(systems.size());
        for (const auto& s : systems) stubById[s.id] = &s;

        std::cout << "Best: hops=" << bestStats.hops
                  << "  dist=" << std::fixed << std::setprecision(2) << bestStats.distanceLy
                  << " ly"
                  << "  cost(" << costModel << ")=" << std::fixed << std::setprecision(2) << bestStats.cost
                  << "  fuelEst=" << std::fixed << std::setprecision(2) << bestFuelEstimate
                  << "  (visited=" << bestStats.visited << ")\n";

        auto printPath = [&](const std::vector<sim::SystemId>& path) {
          for (std::size_t i = 0; i < path.size(); ++i) {
            const auto id = path[i];
            const auto it = stubById.find(id);
            const auto* s = (it != stubById.end()) ? it->second : nullptr;

            std::cout << "  " << std::setw(2) << i << ": ";
            if (s) {
              std::cout << s->name << " (" << s->id << ")";
            } else {
              std::cout << "(unknown " << (unsigned long long)id << ")";
            }

            if (i > 0 && s) {
              const auto itPrev = stubById.find(path[i - 1]);
              const auto* p = (itPrev != stubById.end()) ? itPrev->second : nullptr;
              if (p) {
                const double d = (s->posLy - p->posLy).length();
                std::cout << "  +" << std::fixed << std::setprecision(2) << d << " ly";
              }
            }
            std::cout << "\n";
          }
        };

        // Print best path in full.
        printPath(bestRoute);

        if (kRoutes > 1) {
          if (routes.size() <= 1) {
            std::cout << "(No alternative loopless routes found.)\n";
          } else {
            std::cout << "\nAlternatives (ranked):\n";
            for (std::size_t ri = 0; ri < routes.size(); ++ri) {
              const auto& r = routes[ri];
              const double fuelEst = r.path.empty() ? 0.0 : (r.hops * fuelBase + r.distanceLy * fuelPerLy);
              std::cout << "#" << ri
                        << " hops=" << r.hops
                        << " dist=" << std::fixed << std::setprecision(2) << r.distanceLy << " ly"
                        << " cost=" << std::fixed << std::setprecision(2) << r.cost
                        << " fuelEst=" << std::fixed << std::setprecision(2) << fuelEst
                        << "\n";
            }
            std::cout << "Tip: set --routeCost dist for shortest-distance alternatives, or fuel for fuel-like trade routes.\n";
          }
        }
      }
    }
  }

  if (doMissions && !systems.empty()) {
    fromSysIdx = std::min(fromSysIdx, systems.size() - 1);
    const auto& fromStub = systems[fromSysIdx];
    const auto& fromSys = u.getSystem(fromStub.id, &fromStub);
    if (fromSys.stations.empty()) {
      std::cerr << "\n[missions] chosen system has no stations\n";
      if (json) {
        j.key("missions");
        j.beginObject();
        j.key("error"); j.value("chosen system has no stations");
        j.endObject();
      }
    } else {
      fromStationIdx = std::min(fromStationIdx, fromSys.stations.size() - 1);
      const auto& dockedStation = fromSys.stations[fromStationIdx];

      sim::SaveGame save{};
      if (!loadPath.empty()) {
        if (!sim::loadFromFile(loadPath, save)) {
          std::cerr << "Failed to load save from: " << loadPath << "\n";
          return 1;
        }
      }

      // Tool defaults / overrides.
      save.seed = seed;
      save.timeDays = timeDays;
      save.currentSystem = fromStub.id;
      save.dockedStation = dockedStation.id;
      save.credits = credits;
      save.cargoCapacityKg = cargoCapacityKg;
      save.passengerSeats = std::max(0, seats);

      // Ensure offers exist + are deterministic.
      sim::refreshMissionOffers(u, fromSys, dockedStation, timeDays, rep, save);

      if (json) {
        j.key("missions");
        j.beginObject();
        j.key("station");
        j.beginObject();
        j.key("systemId"); j.value((unsigned long long)fromStub.id);
        j.key("systemName"); j.value(fromStub.name);
        j.key("stationId"); j.value((unsigned long long)dockedStation.id);
        j.key("stationName"); j.value(dockedStation.name);
        j.key("rep"); j.value(rep);
        j.endObject();

        // Expose the contextual tuning inputs that shape the mission board.
        {
          const auto sec = sim::systemSecurityProfile(u.seed(), fromSys);
          const auto issuer = sim::factionProfile(u.seed(), dockedStation.factionId);
          const auto w = sim::computeMissionTypeWeights(sim::MissionBoardParams{}, sec, issuer, rep);

          j.key("context");
          j.beginObject();
          j.key("controllingFactionId"); j.value((unsigned long long)sec.controllingFactionId);
          j.key("security01"); j.value(sec.security01);
          j.key("piracy01"); j.value(sec.piracy01);
          j.key("traffic01"); j.value(sec.traffic01);
          j.key("contest01"); j.value(sec.contest01);
          j.key("issuerTraits");
          j.beginObject();
          j.key("authority"); j.value(issuer.authority);
          j.key("corruption"); j.value(issuer.corruption);
          j.key("wealth"); j.value(issuer.wealth);
          j.key("stability"); j.value(issuer.stability);
          j.key("tech"); j.value(issuer.tech);
          j.key("militarism"); j.value(issuer.militarism);
          j.endObject();
          j.endObject();

          j.key("weights");
          j.beginObject();
          j.key("courier"); j.value(w.wCourier);
          j.key("delivery"); j.value(w.wDelivery);
          j.key("multiDelivery"); j.value(w.wMultiDelivery);
          j.key("escort"); j.value(w.wEscort);
          j.key("salvage"); j.value(w.wSalvage);
          j.key("passenger"); j.value(w.wPassenger);
          j.key("smuggle"); j.value(w.wSmuggle);
          j.key("bountyScan"); j.value(w.wBountyScan);
          j.endObject();
        }
        j.key("offers");
        j.beginArray();
        for (const auto& m : save.missionOffers) writeMissionJson(j, m);
        j.endArray();
      } else {
        std::cout << "\n--- Mission board (day=" << timeDays << ", rep=" << rep << ") ---\n";
        std::cout << "At: " << fromStub.name << " / " << dockedStation.name << "\n";

        {
          const auto sec = sim::systemSecurityProfile(u.seed(), fromSys);
          const auto issuer = sim::factionProfile(u.seed(), dockedStation.factionId);
          const auto w = sim::computeMissionTypeWeights(sim::MissionBoardParams{}, sec, issuer, rep);

          std::cout << "Local context: security=" << sec.security01
                    << " piracy=" << sec.piracy01
                    << " traffic=" << sec.traffic01
                    << " contest=" << sec.contest01
                    << " ctrlFaction=" << sec.controllingFactionId
                    << "\n";
          std::cout << "Weights: courier=" << w.wCourier
                    << " delivery=" << w.wDelivery
                    << " multi=" << w.wMultiDelivery
                    << " escort=" << w.wEscort
                    << " salvage=" << w.wSalvage
                    << " passenger=" << w.wPassenger
                    << " smuggle=" << w.wSmuggle
                    << " bountyScan=" << w.wBountyScan
                    << "\n";
        }
        if (save.missionOffers.empty()) {
          std::cout << "No mission offers.\n";
        } else {
          for (std::size_t i = 0; i < save.missionOffers.size(); ++i) {
            std::cout << "Offer #" << i << ":\n";
            printMission(save.missionOffers[i]);
          }
        }
      }

      std::string acceptError;
      bool accepted = false;
      sim::Mission acceptedMission{};

      if (acceptOffer >= 0) {
        const std::size_t idx = (std::size_t)acceptOffer;
        if (idx >= save.missionOffers.size()) {
          acceptError = "acceptOffer index out of range";
        } else {
          acceptedMission = save.missionOffers[idx];
          accepted = sim::acceptMissionOffer(u, dockedStation, timeDays, rep, save, idx, &acceptError);
        }
      }

      if (json) {
        if (acceptOffer >= 0) {
          j.key("accept");
          j.beginObject();
          j.key("ok"); j.value(accepted);
          if (!acceptError.empty()) {
            j.key("error"); j.value(acceptError);
          }
          j.key("credits"); j.value(save.credits);
          j.key("cargoCapacityKg"); j.value(save.cargoCapacityKg);
          j.key("passengerSeats"); j.value(save.passengerSeats);
          j.key("activeMissions");
          j.beginArray();
          for (const auto& m : save.missions) writeMissionJson(j, m);
          j.endArray();
          j.endObject();
        }

        if (accepted && autoComplete) {
          // Best-effort completion simulation.
          sim::MissionDockResult r{};
          save.timeDays += advanceDays;

          const auto completeAt = [&](sim::SystemId sysId, sim::StationId stId) {
            const auto& sys = u.getSystem(sysId, nullptr);
            const sim::Station* st = nullptr;
            for (const auto& s : sys.stations) {
              if (s.id == stId) { st = &s; break; }
            }
            if (!st) return;
            save.currentSystem = sysId;
            save.dockedStation = stId;
            r = sim::tryCompleteMissionsAtDock(u, sys, *st, save.timeDays, save);
          };

          if (acceptedMission.type == sim::MissionType::MultiDelivery && acceptedMission.viaStation != 0) {
            completeAt(acceptedMission.viaSystem, acceptedMission.viaStation);
            completeAt(acceptedMission.toSystem, acceptedMission.toStation);
          } else {
            completeAt(acceptedMission.toSystem, acceptedMission.toStation);
          }

          j.key("autoComplete");
          j.beginObject();
          j.key("advancedDays"); j.value(advanceDays);
          j.key("completed"); j.value(r.completed);
          j.key("progressedMultiLeg"); j.value(r.progressedMultiLeg);
          j.key("credits"); j.value(save.credits);
          j.endObject();
        }

        j.endObject(); // missions
      } else {
        if (acceptOffer >= 0) {
          if (!accepted) {
            std::cout << "\n[accept] failed: " << acceptError << "\n";
          } else {
            std::cout << "\n[accept] ok. credits=" << std::fixed << std::setprecision(0) << save.credits
                      << " activeMissions=" << save.missions.size() << "\n";
          }
        }

        if (accepted && autoComplete) {
          save.timeDays += advanceDays;
          sim::MissionDockResult r{};
          const auto completeAt = [&](sim::SystemId sysId, sim::StationId stId) {
            const auto& sys = u.getSystem(sysId, nullptr);
            const sim::Station* st = nullptr;
            for (const auto& s : sys.stations) {
              if (s.id == stId) { st = &s; break; }
            }
            if (!st) {
              std::cout << "[autoComplete] station not found: " << stId << "\n";
              return;
            }
            save.currentSystem = sysId;
            save.dockedStation = stId;
            r = sim::tryCompleteMissionsAtDock(u, sys, *st, save.timeDays, save);
          };

          if (acceptedMission.type == sim::MissionType::MultiDelivery && acceptedMission.viaStation != 0) {
            completeAt(acceptedMission.viaSystem, acceptedMission.viaStation);
            completeAt(acceptedMission.toSystem, acceptedMission.toStation);
          } else {
            completeAt(acceptedMission.toSystem, acceptedMission.toStation);
          }

          std::cout << "[autoComplete] advancedDays=" << advanceDays
                    << " completed=" << r.completed
                    << " progressedMultiLeg=" << r.progressedMultiLeg
                    << " credits=" << std::fixed << std::setprecision(0) << save.credits
                    << "\n";
        }
      }

      if (!savePath.empty()) {
        if (!sim::saveToFile(save, savePath)) {
          std::cerr << "Failed to write save to: " << savePath << "\n";
          return 1;
        }
        if (!json) {
          std::cout << "[save] wrote " << savePath << "\n";
        }
      }
    }
  }

  if (json) {
    j.endObject();
    *jsonStream << "\n";
  }

  return 0;
}
