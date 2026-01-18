#include "stellar/sim/MissionBriefing.h"

#include "stellar/core/Clamp.h"
#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/econ/Commodity.h"
#include "stellar/sim/BlackMarket.h"
#include "stellar/sim/Contraband.h"
#include "stellar/sim/Law.h"
#include "stellar/sim/SystemConditions.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace stellar::sim {

static double clamp01(double x) {
  if (!std::isfinite(x)) return 0.0;
  return std::clamp(x, 0.0, 1.0);
}

static const SystemSecurityDeltaState* findDelta(std::span<const SystemSecurityDeltaState> deltas, SystemId sysId) {
  if (sysId == 0) return nullptr;
  for (const auto& d : deltas) {
    if (d.systemId == sysId) return &d;
  }
  return nullptr;
}

static const Station* findStationById(const StarSystem& sys, StationId stId) {
  if (stId == 0) return nullptr;
  for (const auto& st : sys.stations) {
    if (st.id == stId) return &st;
  }
  return nullptr;
}

static std::string factionNameById(const Universe& u, core::u32 factionId) {
  if (factionId == 0) return "Independent";
  for (const auto& f : u.factions()) {
    if (f.id == factionId) return f.name;
  }
  return std::string("Faction #") + std::to_string((core::u64)factionId);
}

static std::pair<SystemId, StationId> missionNextStop(const Mission& m) {
  const bool hasViaStop = (m.viaSystem != 0 && m.viaStation != 0 && m.leg == 0
                           && (m.type == MissionType::MultiDelivery || m.type == MissionType::Passenger));
  if (hasViaStop) {
    return {m.viaSystem, m.viaStation};
  }

  // Salvage jobs: before the site is visited, the objective is the mission signal
  // (not a station). We still treat the destination system as m.toSystem.
  if (m.type == MissionType::Salvage && !m.scanned) {
    return {m.toSystem, 0};
  }

  // Bounty missions can provide a station hint where the target tends to lurk.
  if (m.type == MissionType::BountyScan || m.type == MissionType::BountyKill) {
    return {m.toSystem, m.toStation};
  }

  return {m.toSystem, m.toStation};
}

static double missionTypeDangerBias(const Mission& m) {
  switch (m.type) {
    case MissionType::Courier: return 0.02;
    case MissionType::Delivery: return 0.05;
    case MissionType::MultiDelivery: return 0.08;
    case MissionType::Passenger: return 0.06;
    case MissionType::Smuggle: return 0.08;
    case MissionType::Salvage: return 0.14;
    case MissionType::Escort: return 0.18;
    case MissionType::BountyScan: return 0.20;
    case MissionType::BountyKill: return 0.28;
    default: return 0.05;
  }
}

static std::string contractCodeFromSeed(core::u64 seed) {
  core::SplitMix64 r(seed);

  auto letter = [&](int i) -> char {
    return (char)('A' + (r.nextU32() + (core::u32)i) % 26u);
  };

  const core::u32 n = r.nextU32();
  const core::u32 m = r.nextU32();

  std::ostringstream oss;
  oss << letter(0) << letter(1) << '-';

  // Keep it readable: 4 hex digits + 2 base-36 digits.
  const core::u32 h = (n ^ (m * 0x9E3779B9u));
  const core::u32 h16 = h & 0xFFFFu;
  oss.setf(std::ios::uppercase);
  oss << std::hex;
  oss.width(4);
  oss.fill('0');
  oss << h16;

  auto base36 = [](core::u32 v) -> char {
    v %= 36u;
    if (v < 10u) return (char)('0' + v);
    return (char)('A' + (v - 10u));
  };

  oss << base36((h >> 16) & 0xFFu);
  oss << base36((h >> 24) & 0xFFu);

  return oss.str();
}

static std::string contactNameFromSeed(core::u64 seed, MissionType type) {
  core::SplitMix64 r(seed);

  static const std::array<const char*, 12> kNeutralCalls{
      "Kestrel", "Orchid", "Vega", "Atlas", "Meridian", "Helios",
      "Juniper", "Sable", "Ardent", "Pioneer", "Cobalt", "Lumen",
  };
  static const std::array<const char*, 10> kNeutralOrgs{
      "Logistics", "Courier", "Freight", "Charters", "Consortium",
      "Exchange", "Transit", "Holdings", "Brokers", "Dispatch",
  };

  static const std::array<const char*, 12> kCovertCalls{
      "Ghost", "Noir", "Cipher", "Vesper", "Wraith", "Umbra",
      "Nyx", "Raven", "Specter", "Hush", "Violet", "Ash",
  };
  static const std::array<const char*, 10> kCovertOrgs{
      "Syndicate", "Market", "Ledger", "Circuit", "Backroom",
      "Fence", "Brokerage", "Cartel", "Network", "Drop",
  };

  const bool covert = (type == MissionType::Smuggle);

  const auto& A = covert ? kCovertCalls : kNeutralCalls;
  const auto& B = covert ? kCovertOrgs : kNeutralOrgs;

  const char* a = A[(std::size_t)(r.nextU32() % (core::u32)A.size())];
  const char* b = B[(std::size_t)(r.nextU32() % (core::u32)B.size())];

  std::string out = a;
  out += ' ';
  out += b;
  return out;
}

static std::string riskTierMarkup(double risk01, bool useMarkup) {
  const char* tier = riskTierName(risk01);
  if (!useMarkup) return tier;

  // Subtle UI flair: keep tags simple and readable.
  if (risk01 >= 0.75) return std::string("[shake][color #ff4a4a]") + tier + "[/color][/shake]";
  if (risk01 >= 0.50) return std::string("[color #ff9a4a]") + tier + "[/color]";
  if (risk01 >= 0.25) return std::string("[color #ffd36a]") + tier + "[/color]";
  return std::string("[color #7cffb2]") + tier + "[/color]";
}

static std::string titleForType(MissionType t, bool useMarkup) {
  const char* plain = "Contract";
  const char* gradA = "#ffffff";
  const char* gradB = "#ffffff";

  switch (t) {
    case MissionType::Courier:      plain = "Courier Contract";       gradA = "#7cffb2"; gradB = "#5ad1ff"; break;
    case MissionType::Delivery:     plain = "Freight Delivery";       gradA = "#ffd36a"; gradB = "#ff9a4a"; break;
    case MissionType::MultiDelivery:plain = "Multi-hop Logistics";    gradA = "#ffd36a"; gradB = "#ff7c32"; break;
    case MissionType::Passenger:    plain = "Passenger Charter";      gradA = "#5ad1ff"; gradB = "#b86aff"; break;
    case MissionType::Smuggle:      plain = "Covert Delivery";        gradA = "#ff4a4a"; gradB = "#b86aff"; break;
    case MissionType::Salvage:      plain = "Salvage Retrieval";      gradA = "#ff9a4a"; gradB = "#ffd36a"; break;
    case MissionType::Escort:       plain = "Convoy Escort";          gradA = "#ffd36a"; gradB = "#ffffff"; break;
    case MissionType::BountyScan:   plain = "Bounty Scan";            gradA = "#ff4a4a"; gradB = "#ff9a4a"; break;
    case MissionType::BountyKill:   plain = "Bounty Hunt";            gradA = "#ff4a4a"; gradB = "#ffd36a"; break;
    default:                        plain = "Contract";               gradA = "#ffffff"; gradB = "#ffffff"; break;
  }

  if (!useMarkup) return plain;
  return std::string("[grad ") + gradA + " " + gradB + "]" + plain + "[/grad]";
}

const char* riskTierName(double risk01) {
  const double r = clamp01(risk01);
  if (r < 0.25) return "LOW";
  if (r < 0.50) return "MED";
  if (r < 0.75) return "HIGH";
  return "EXTREME";
}

MissionRisk computeMissionRisk(Universe& universe,
                               const StarSystem& originSystem,
                               const Station& originStation,
                               double timeDays,
                               double playerRepWithIssuerFaction,
                               const Mission& mission,
                               std::span<const SystemSecurityDeltaState> securityDeltas,
                               const MissionBriefingParams& params) {
  MissionRisk r{};

  if (!std::isfinite(timeDays)) return r;

  const auto [dstSysId, dstStId] = missionNextStop(mission);
  if (dstSysId == 0) return r;

  const StarSystem* dstSysPtr = &originSystem;
  if (dstSysId != originSystem.stub.id) {
    dstSysPtr = &universe.getSystem(dstSysId);
  }
  const StarSystem& dstSys = *dstSysPtr;

  const math::Vec3d dpos = dstSys.stub.posLy - originSystem.stub.posLy;
  r.distanceLy = std::isfinite(dpos.length()) ? dpos.length() : 0.0;

  const SystemSecurityDeltaState* delta = nullptr;
  if (params.applySecurityDeltas) {
    delta = findDelta(securityDeltas, dstSysId);
  }

  const SystemSecurityProfile sec = effectiveSystemSecurityProfile(universe.seed(),
                                                                   dstSys,
                                                                   timeDays,
                                                                   delta,
                                                                   params.dynamicsParams,
                                                                   params.eventParams,
                                                                   nullptr);

  r.security01 = clamp01(sec.security01);
  r.piracy01 = clamp01(sec.piracy01);
  r.traffic01 = clamp01(sec.traffic01);
  r.contest01 = clamp01(sec.contest01);

  const double dist01 = clamp01(r.distanceLy / 420.0);

  // Baseline danger: piracy pressure dominates, contestedness adds volatility,
  // and low security increases expected time-to-response.
  double danger = 0.0;
  danger += 0.56 * r.piracy01;
  danger += 0.20 * r.contest01;
  danger += 0.16 * (1.0 - r.security01);
  danger += 0.08 * dist01;

  danger += missionTypeDangerBias(mission);
  danger = clamp01(danger);
  r.danger01 = danger;

  // Legality/customs risk (primarily for smuggling).
  double lawRisk = 0.0;
  if (mission.type == MissionType::Smuggle) {
    const Station* dstSt = findStationById(dstSys, (dstStId != 0 ? dstStId : mission.toStation));

    if (dstSt) {
      const LawProfile law = lawProfile(universe.seed(), dstSt->factionId);

      const double strict01 = clamp01((law.scanStrictness - 0.70) / 0.95);
      const double corruption01 = clamp01(law.corruption);

      const bool illegal = isIllegalCommodityAtStation(universe.seed(),
                                                       dstSt->factionId,
                                                       dstSt->id,
                                                       dstSt->type,
                                                       mission.commodity);
      r.contrabandIllegalAtDestination = illegal;

      const BlackMarketProfile bm = blackMarketProfile(universe.seed(),
                                                       dstSys.stub.id,
                                                       *dstSt,
                                                       sec,
                                                       law,
                                                       timeDays,
                                                       playerRepWithIssuerFaction);
      r.blackMarketAccess01 = clamp01(bm.access01);
      r.stingChance01 = clamp01(bm.stingChance);

      // Expected fine magnitude for the mission payload.
      const double units = std::max(0.0, mission.units);
      const double unitMid = std::max(0.0, econ::commodityDef(mission.commodity).basePrice);
      const double illegalValueCr = units * unitMid;
      r.expectedFineCr = std::max(0.0, law.contrabandFineCr(illegalValueCr));

      // Convert to a 0..1 risk proxy.
      // - strict factions push up scans
      // - low corruption reduces bribe odds
      // - higher system security increases enforcement effectiveness
      // - higher stingChance captures "this feels like a trap" black market profile
      lawRisk = 0.60 * (0.65 * strict01 + 0.35 * (1.0 - corruption01));
      lawRisk += 0.40 * clamp01(bm.risk01);
      lawRisk *= (0.70 + 0.70 * r.security01);

      // Higher volume adds a little extra attention.
      const double vol01 = clamp01(units / 60.0);
      lawRisk = clamp01(lawRisk + 0.10 * vol01);
    }
  }

  r.lawRisk01 = clamp01(lawRisk);

  // Overall: for most missions danger dominates; for smuggling the law component
  // matters much more.
  const double wLaw = (mission.type == MissionType::Smuggle) ? 0.42 : 0.18;
  r.overall01 = clamp01((1.0 - wLaw) * r.danger01 + wLaw * r.lawRisk01);

  return r;
}

MissionBriefing generateMissionBriefing(Universe& universe,
                                        const StarSystem& originSystem,
                                        const Station& originStation,
                                        double timeDays,
                                        double playerRepWithIssuerFaction,
                                        const Mission& mission,
                                        std::span<const SystemSecurityDeltaState> securityDeltas,
                                        const MissionBriefingParams& params) {
  MissionBriefing b{};

  const core::u64 seedBase = core::hashCombine(core::seedFromText("mission_briefing_v1"), (core::u64)mission.id);
  const core::u64 seed = core::hashCombine(seedBase,
                                          core::hashCombine((core::u64)(core::i64)mission.type,
                                                            core::hashCombine((core::u64)mission.factionId,
                                                                              (core::u64)mission.targetNpcId)));

  b.contractCode = contractCodeFromSeed(seed);
  b.contactName = contactNameFromSeed(seed, mission.type);

  b.risk = computeMissionRisk(universe,
                              originSystem,
                              originStation,
                              timeDays,
                              playerRepWithIssuerFaction,
                              mission,
                              securityDeltas,
                              params);

  const auto [nextSysId, nextStId] = missionNextStop(mission);

  const StarSystem* nextSysPtr = &originSystem;
  if (nextSysId != originSystem.stub.id) {
    nextSysPtr = &universe.getSystem(nextSysId);
  }
  const StarSystem& nextSys = *nextSysPtr;

  const Station* nextSt = findStationById(nextSys, nextStId);

  const std::string sysName = nextSys.stub.name;
  const std::string stName = nextSt ? nextSt->name : std::string("(site)");

  b.titleMarkup = titleForType(mission.type, params.useMarkup);

  // Short synopsis
  {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss.precision(0);

    switch (mission.type) {
      case MissionType::Courier:
        oss << "Deliver encrypted data to " << stName << " (" << sysName << ").";
        break;
      case MissionType::Delivery:
        oss << "Deliver " << (int)std::llround(mission.units) << " units of "
            << econ::commodityDef(mission.commodity).name << " to "
            << stName << " (" << sysName << ").";
        break;
      case MissionType::MultiDelivery:
        oss << "Multi-hop delivery: " << econ::commodityDef(mission.commodity).name
            << " x" << (int)std::llround(mission.units) << ".";
        break;
      case MissionType::Passenger:
        oss << "Charter: transport party of " << (int)std::llround(mission.units) << " passengers.";
        break;
      case MissionType::Smuggle:
        oss << "Move " << (int)std::llround(mission.units) << " units of "
            << econ::commodityDef(mission.commodity).name << " to "
            << stName << " without attracting attention.";
        break;
      case MissionType::Salvage:
        oss << "Recover loose cargo from a derelict signal in " << sysName << " and return.";
        break;
      case MissionType::Escort:
        oss << "Escort a convoy to " << stName << " (" << sysName << ") and keep it alive.";
        break;
      case MissionType::BountyScan:
        oss << "Locate and scan a wanted pilot in " << sysName << ".";
        break;
      case MissionType::BountyKill:
        oss << "Eliminate a wanted pilot in " << sysName << ".";
        break;
      default:
        oss << "Complete the contract objectives.";
        break;
    }

    b.synopsisMarkup = oss.str();
  }

  // Bullets (objective + risk + logistics)
  {
    std::vector<std::string> lines;
    lines.reserve(12);

    const std::string issuerFaction = factionNameById(universe, mission.factionId);

    lines.push_back(std::string("Contract ID: ") + b.contractCode);
    lines.push_back(std::string("Contact: ") + b.contactName + "  (" + issuerFaction + ")");

    {
      std::ostringstream oss;
      oss.setf(std::ios::fixed);
      oss.precision(0);
      oss << "Reward: " << mission.reward << " cr";
      lines.push_back(oss.str());
    }

    if (mission.deadlineDay > 0.0 && std::isfinite(mission.deadlineDay)) {
      const double hrsLeft = std::max(0.0, (mission.deadlineDay - timeDays) * 24.0);
      std::ostringstream oss;
      oss.setf(std::ios::fixed);
      oss.precision(1);
      oss << "Deadline: " << hrsLeft << " h";
      lines.push_back(oss.str());
    }

    {
      std::ostringstream oss;
      oss.setf(std::ios::fixed);
      oss.precision(1);
      oss << "Route: ~" << b.risk.distanceLy << " ly straight-line";
      lines.push_back(oss.str());
    }

    // Objective-specific details.
    switch (mission.type) {
      case MissionType::Courier: {
        lines.push_back(std::string("Delivery point: ") + stName + " (" + sysName + ")");
        lines.push_back("Payload: encrypted data-core (non-cargo)." );
      } break;
      case MissionType::Delivery:
      case MissionType::MultiDelivery:
      case MissionType::Smuggle:
      case MissionType::Salvage:
      case MissionType::Escort: {
        const double units = std::max(0.0, mission.units);
        const auto& def = econ::commodityDef(mission.commodity);
        const double massKg = std::max(0.0, units * std::max(0.0, def.massKg));

        std::ostringstream oss;
        oss.setf(std::ios::fixed);
        oss.precision(0);
        oss << "Cargo: " << def.name << " x" << (int)std::llround(units) << "  (" << massKg << " kg)";
        lines.push_back(oss.str());

        if (mission.type == MissionType::Smuggle) {
          lines.push_back(std::string("Status: ") + (params.useMarkup ? "[shake][color #ff4a4a]CONTRABAND[/color][/shake]" : "CONTRABAND"));
        } else if (mission.type == MissionType::Salvage) {
          // Salvage uses targetNpcId as a stable signal id (high bit set).
          std::ostringstream site;
          site << "Mission site: derelict signal 0x";
          site.setf(std::ios::uppercase);
          site << std::hex << (mission.targetNpcId & 0xFFFFFFFFull);
          lines.push_back(site.str());
        }

        if (mission.type == MissionType::Delivery || mission.type == MissionType::MultiDelivery) {
          lines.push_back(std::string("Acquisition: ") + (mission.cargoProvided ? "cargo provided" : "procure on market"));
        }
        if (mission.type == MissionType::Smuggle) {
          lines.push_back("Acquisition: cargo provided on acceptance." );
        }
      } break;
      case MissionType::Passenger: {
        std::ostringstream oss;
        oss << "Passengers: " << (int)std::llround(std::max(0.0, mission.units)) << " (requires seats)";
        lines.push_back(oss.str());
        if (mission.viaSystem != 0 && mission.viaStation != 0 && mission.leg == 0) {
          const auto& viaSys = universe.getSystem(mission.viaSystem);
          const Station* viaSt = findStationById(viaSys, mission.viaStation);
          lines.push_back(std::string("Stopover: ") + (viaSt ? viaSt->name : std::string("(station)")) + " (" + viaSys.stub.name + ")");
        }
        lines.push_back("Guidance: avoid combat and aggressive thermal loads." );
      } break;
      case MissionType::BountyScan:
      case MissionType::BountyKill: {
        std::ostringstream oss;
        oss << "Target tag: " << (core::u64)mission.targetNpcId;
        lines.push_back(oss.str());
        if (mission.type == MissionType::BountyScan) {
          lines.push_back("Objective: perform a close-range scan." );
        } else {
          lines.push_back("Objective: confirm elimination." );
        }
      } break;
      default:
        break;
    }

    // Risk / conditions.
    {
      std::ostringstream oss;
      oss.setf(std::ios::fixed);
      oss.precision(2);
      oss << "Threat: " << riskTierMarkup(b.risk.overall01, params.useMarkup)
          << "  (danger " << b.risk.danger01
          << ", law " << b.risk.lawRisk01
          << ", piracy " << b.risk.piracy01
          << ", security " << b.risk.security01 << ")";
      lines.push_back(oss.str());
    }

    if (mission.type == MissionType::Smuggle) {
      std::ostringstream oss;
      oss.setf(std::ios::fixed);
      oss.precision(0);

      const bool illegal = b.risk.contrabandIllegalAtDestination;
      oss << "Customs: " << (illegal ? "illegal" : "legal")
          << " | black market access " << (int)std::llround(b.risk.blackMarketAccess01 * 100.0) << "%"
          << " | sting chance " << (int)std::llround(b.risk.stingChance01 * 100.0) << "%";
      lines.push_back(oss.str());

      if (b.risk.expectedFineCr > 1.0) {
        std::ostringstream fine;
        fine.setf(std::ios::fixed);
        fine.precision(0);
        fine << "If scanned: expected fine ~" << b.risk.expectedFineCr << " cr (before rep penalties).";
        lines.push_back(fine.str());
      }

      lines.push_back("Recommendation: approach cold, avoid loitering in comms range." );
    } else {
      // Generic recommendations (based on danger).
      if (b.risk.overall01 >= 0.60) {
        lines.push_back("Recommendation: bring shields, countermeasures, and enough fuel for reroutes." );
      } else if (b.risk.overall01 >= 0.35) {
        lines.push_back("Recommendation: stay alert; plot a safe refuel stop if needed." );
      }
    }

    b.bulletsMarkup = std::move(lines);
  }

  return b;
}

} // namespace stellar::sim
