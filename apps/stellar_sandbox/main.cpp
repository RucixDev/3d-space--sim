#include "stellar/core/Args.h"
#include "stellar/core/JsonWriter.h"
#include "stellar/core/Log.h"
#include "stellar/econ/Commodity.h"
#include "stellar/econ/RoutePlanner.h"
#include "stellar/proc/GalaxyGenerator.h"
#include "stellar/sim/MissionLogic.h"
#include "stellar/sim/Contraband.h"
#include "stellar/sim/Law.h"
#include "stellar/sim/PoliceScan.h"
#include "stellar/sim/Industry.h"
#include "stellar/sim/Warehouse.h"
#include "stellar/sim/Reputation.h"
#include "stellar/sim/TradeScanner.h"
#include "stellar/sim/IndustryScanner.h"
#include "stellar/sim/NavRoute.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/Signature.h"
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
            << "  --sig                  Include a stable signature per system stub (determinism)\n"
            << "  --sysSig               Include a deep signature of full generated systems (slower)\n"
            << "  --json                 Emit machine-readable JSON to stdout (also works with --out)\n"
            << "  --out <path>           Write JSON output to a file instead of stdout ('-' means stdout)\n"
            << "\n"
            << "Trade scanning (headless route planner):\n"
            << "  --trade                Print best trade opportunities from a chosen station\n"
            << "  --tradeMix             Print best multi-commodity cargo manifests (trade mix)\n"
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
            << "Law / contraband (headless):\n"
            << "  --law                 Print law + contraband profile for the chosen station\n"
            << "  --faction <id>         Override faction id (default: station faction)\n"
            << "  --lawRep <r>           Player reputation used for bribe chance (default: 0)\n"
            << "  --heat <h>             Player heat used for bribe chance (default: 0)\n"
            << "  --illegalValue <cr>    Illegal cargo value used for fine/bribe examples (default: 5000)\n"
            << "  --smuggleHoldMk <mk>   Smuggle hold grade (0-3) for scan math (default: 0)\n"
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

  const bool json = args.hasFlag("json");
  std::string outPath;
  (void)args.getString("out", outPath);

  const bool emitStubSig = args.hasFlag("sig");
  const bool emitSysSig = args.hasFlag("sysSig");

  const bool doTrade = args.hasFlag("trade");
  const bool doTradeMix = args.hasFlag("tradeMix");
  const bool doIndustry = args.hasFlag("industry");
  const bool doIndustryTrade = args.hasFlag("industryTrade");
  const bool doWarehouse = args.hasFlag("warehouse");
  std::size_t fromSysIdx = 0;
  std::size_t fromStationIdx = 0;
  {
    unsigned long long v = 0;
    if (args.getU64("fromSys", v)) fromSysIdx = (std::size_t)v;
    if (args.getU64("fromStation", v)) fromStationIdx = (std::size_t)v;
  }
  double cargoCapacityKg = 420.0;
  (void)args.getDouble("cargoKg", cargoCapacityKg);
  std::size_t tradeLimit = 12;
  {
    unsigned long long v = 0;
    if (args.getU64("tradeLimit", v)) tradeLimit = (std::size_t)v;
  }


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

  sim::Universe u(seed);

  const auto systems = u.queryNearby(posLy, radiusLy, limit);

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
      if (emitSysSig) {
        const auto& sys = u.getSystem(s.id, &s);
        j.key("sysSig");
        j.value((unsigned long long)sim::signatureStarSystem(sys));
      }
      j.endObject();
    }
    j.endArray();
  } else {
    std::cout << "Seed: " << seed << "\n";
    std::cout << "Query @ (" << posLy.x << "," << posLy.y << "," << posLy.z << ") radius=" << radiusLy << " ly\n";
    std::cout << "Found " << systems.size() << " systems\n\n";

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
      if (emitSysSig) {
        const auto& sys = u.getSystem(s.id, &s);
        const auto sig = sim::signatureStarSystem(sys);
        std::cout << "  sysSig=" << (unsigned long long)sig;
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
      std::cout << "  - " << st.name << " (fee=" << st.feeRate << ")\n";
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

    sim::RoutePlanStats stats{};
    std::vector<sim::SystemId> route;
    std::string costModel = "hops";
    {
      std::string rc = routeCost;
      std::transform(rc.begin(), rc.end(), rc.begin(), [](unsigned char c) { return (char)std::tolower(c); });
      if (rc == "dist" || rc == "distance") {
        costModel = "dist";
        route = sim::plotRouteAStarCost(systems, fromStub.id, toStub.id, jumpRangeLy, 0.0, 1.0, &stats);
      } else if (rc == "fuel") {
        costModel = "fuel";
        route = sim::plotRouteAStarCost(systems, fromStub.id, toStub.id, jumpRangeLy, fuelBase, fuelPerLy, &stats);
      } else {
        costModel = "hops";
        route = sim::plotRouteAStarHops(systems, fromStub.id, toStub.id, jumpRangeLy, &stats);
      }
    }
    const double fuelEstimate = route.empty() ? 0.0 : (stats.hops * fuelBase + stats.distanceLy * fuelPerLy);

    if (json) {
      j.key("route");
      j.beginObject();
      j.key("fromSysIdx"); j.value((unsigned long long)fromSysIdx);
      j.key("toSysIdx"); j.value((unsigned long long)toSysIdx);
      j.key("fromSystemId"); j.value((unsigned long long)fromStub.id);
      j.key("toSystemId"); j.value((unsigned long long)toStub.id);
      j.key("jumpRangeLy"); j.value(jumpRangeLy);
      j.key("costModel"); j.value(costModel);
      j.key("cost"); j.value(stats.cost);
      j.key("fuelBase"); j.value(fuelBase);
      j.key("fuelPerLy"); j.value(fuelPerLy);
      j.key("fuelEstimate"); j.value(fuelEstimate);
      j.key("found"); j.value(!route.empty());
      j.key("hops"); j.value((unsigned long long)(route.size() > 1 ? route.size() - 1 : 0));
      j.key("distanceLy"); j.value(stats.distanceLy);
      j.key("visited"); j.value((unsigned long long)stats.visited);
      j.key("expansions"); j.value((unsigned long long)stats.expansions);
      j.key("path");
      j.beginArray();
      for (const auto id : route) j.value((unsigned long long)id);
      j.endArray();
      j.endObject();
    } else {
      std::cout << "\n--- Route plan (A*, cost=" << costModel << ", jr=" << std::fixed << std::setprecision(1) << jumpRangeLy << " ly) ---\n";
      std::cout << "From: [" << fromSysIdx << "] " << fromStub.name << " (" << fromStub.id << ")\n";
      std::cout << "To:   [" << toSysIdx << "] " << toStub.name << " (" << toStub.id << ")\n";

      if (route.empty()) {
        std::cout << "No route found inside the queried node set.\n";
        std::cout << "Tip: increase --radius (and/or --limit) so intermediate systems are available.\n";
      } else {
        std::unordered_map<sim::SystemId, const sim::SystemStub*> stubById;
        stubById.reserve(systems.size());
        for (const auto& s : systems) stubById[s.id] = &s;

        std::cout << "Hops: " << stats.hops
                  << "  Dist: " << std::fixed << std::setprecision(2) << stats.distanceLy
                  << " ly"
                  << "  Cost(" << costModel << "): " << std::fixed << std::setprecision(2) << stats.cost
                  << "  FuelEst: " << std::fixed << std::setprecision(2) << fuelEstimate
                  << "  (visited=" << stats.visited << ")\n";

        for (std::size_t i = 0; i < route.size(); ++i) {
          const auto id = route[i];
          const auto it = stubById.find(id);
          const auto* s = (it != stubById.end()) ? it->second : nullptr;

          std::cout << "  " << std::setw(2) << i << ": ";
          if (s) {
            std::cout << s->name << " (" << s->id << ")";
          } else {
            std::cout << "(unknown " << (unsigned long long)id << ")";
          }

          if (i > 0 && s) {
            const auto itPrev = stubById.find(route[i - 1]);
            const auto* p = (itPrev != stubById.end()) ? itPrev->second : nullptr;
            if (p) {
              const double d = (s->posLy - p->posLy).length();
              std::cout << "  +" << std::fixed << std::setprecision(2) << d << " ly";
            }
          }
          std::cout << "\n";
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
        j.key("offers");
        j.beginArray();
        for (const auto& m : save.missionOffers) writeMissionJson(j, m);
        j.endArray();
      } else {
        std::cout << "\n--- Mission board (day=" << timeDays << ", rep=" << rep << ") ---\n";
        std::cout << "At: " << fromStub.name << " / " << dockedStation.name << "\n";
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
