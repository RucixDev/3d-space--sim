#include "stellar/sim/MissionLogic.h"

#include "stellar/core/Random.h"
#include "stellar/econ/Cargo.h"
#include "stellar/econ/Market.h"
#include "stellar/sim/Reputation.h"
#include "stellar/sim/Contraband.h"
#include "stellar/sim/FactionProfile.h"
#include "stellar/sim/Orbit.h"
#include "stellar/sim/SecurityModel.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace stellar::sim {

struct CargoMissionSpec {
  econ::CommodityId commodity{econ::CommodityId::Food};
  int units{0};
  econ::MarketQuote originQ{};
  econ::MarketQuote destQ{};
};

// Pick a commodity + unit count for a cargo-style mission.
//
// Goals:
//  - Prefer goods that are in higher demand at the destination (dest mid > origin mid).
//  - Ensure the origin station has enough inventory for the chosen unit count.
//  - Avoid generating missions that would be contradictory to contraband rules.
static bool pickCargoMission(core::SplitMix64& rng,
                             const econ::StationEconomyState& originState,
                             const econ::StationEconomyModel& originModel,
                             const econ::StationEconomyState& destState,
                             const econ::StationEconomyModel& destModel,
                             core::u32 destIllegalMask,
                             bool requireLegalAtDest,
                             bool cargoProvided,
                             double maxCargoKg,
                             int minUnits,
                             int maxUnitsHard,
                             CargoMissionSpec& out) {
  struct Cand {
    econ::CommodityId id{econ::CommodityId::Food};
    double score{0.0};
    econ::MarketQuote o{};
    econ::MarketQuote d{};
    int maxUnits{0};
  };

  minUnits = std::max(1, minUnits);
  maxUnitsHard = std::max(minUnits, maxUnitsHard);
  maxCargoKg = std::max(0.2, maxCargoKg);

  // Provided cargo should be limited to avoid draining the station economy too aggressively.
  const double invFrac = cargoProvided ? 0.28 : 0.80;

  std::vector<Cand> cands;
  cands.reserve(econ::kCommodityCount);

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const auto cid = static_cast<econ::CommodityId>(i);

    if (requireLegalAtDest && ((destIllegalMask & commodityBit(cid)) != 0u)) continue;

    const auto oq = econ::quote(originState, originModel, cid, 0.10);
    const auto dq = econ::quote(destState, destModel, cid, 0.10);

    if (oq.inventory < 1.0) continue;

    const double massKg = std::max(1e-6, econ::commodityDef(cid).massKg);
    const int maxByKg = (int)std::floor(maxCargoKg / massKg);
    const int maxByInv = (int)std::floor(std::max(0.0, oq.inventory) * invFrac);
    const int maxUnits = std::min(maxUnitsHard, std::min(maxByKg, maxByInv));

    if (maxUnits < minUnits) continue;

    // Score cargo by demand: destination mid-price relative to origin mid-price.
    // Also bias slightly towards higher-value goods so boards aren't always dominated by
    // ultra-bulky low-value commodities.
    const double ratio = dq.mid / std::max(1e-6, oq.mid);
    const double valueBias = 0.65 + 0.35 * std::min(1.0, dq.mid / 120.0);
    const double score = ratio * valueBias;

    cands.push_back(Cand{cid, score, oq, dq, maxUnits});
  }

  if (cands.empty()) return false;

  std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) {
    if (a.score != b.score) return a.score > b.score;
    return (int)a.id < (int)b.id;
  });

  const int top = std::min<int>(3, (int)cands.size());
  const int pick = (top <= 1) ? 0 : rng.range(0, top - 1);
  const Cand& c = cands[(std::size_t)pick];

  const int units = rng.range(minUnits, c.maxUnits);
  out.commodity = c.id;
  out.units = units;
  out.originQ = c.o;
  out.destQ = c.d;
  return true;
}

static int activePassengerCount(const SaveGame& s) {
  int n = 0;
  for (const auto& m : s.missions) {
    if (m.completed || m.failed) continue;
    if (m.type != MissionType::Passenger) continue;
    n += std::max(0, (int)std::llround(m.units));
  }
  return n;
}

static core::u64 missionBoardSeed(StationId stationId, int dayStamp, core::u32 factionId) {
  // Similar to the in-game seeding approach: mix stable identifiers.
  // Uses large odd constants to improve bit diffusion.
  return (core::u64)stationId * 1469598103934665603ull ^ (core::u64)dayStamp * 1099511628211ull ^ (core::u64)factionId;
}

MissionTypeWeights computeMissionTypeWeights(const MissionBoardParams& params,
                                             const SystemSecurityProfile& local,
                                             const FactionProfile& issuer,
                                             double playerRep) {
  (void)playerRep; // reserved for future tuning (e.g., "trust" gating)

  MissionTypeWeights w{};
  w.wCourier = params.wCourier;
  w.wDelivery = params.wDelivery;
  w.wMultiDelivery = params.wMultiDelivery;
  w.wEscort = params.wEscort;
  w.wSalvage = params.wSalvage;
  w.wPassenger = params.wPassenger;
  w.wSmuggle = params.wSmuggle;
  w.wBountyScan = params.wBountyScan;

  auto bias01 = [](double x) {
    // Map [0,1] -> [-1,1] around baseline 0.5.
    return std::clamp((x - 0.5) * 2.0, -1.0, 1.0);
  };

  const double piracyB = bias01(local.piracy01);
  const double trafficB = bias01(local.traffic01);
  const double securityB = bias01(local.security01);
  const double contest = std::clamp(local.contest01, 0.0, 1.0);

  const double authorityB = bias01(issuer.authority);
  const double corruptionB = bias01(issuer.corruption);
  const double wealthB = bias01(issuer.wealth);
  const double stabilityB = bias01(issuer.stability);
  const double techB = bias01(issuer.tech);
  const double militarismB = bias01(issuer.militarism);

  // Apply small, bounded multipliers. These are tuned so that defaults remain
  // close to MissionBoardParams while still producing noticeable shifts for
  // extreme systems/factions.
  w.wCourier = std::max(0.0, w.wCourier * (1.0 + 0.22 * trafficB - 0.10 * piracyB));
  w.wDelivery = std::max(0.0, w.wDelivery * (1.0 + 0.18 * trafficB + 0.10 * wealthB));
  w.wMultiDelivery = std::max(0.0, w.wMultiDelivery * (1.0 + 0.25 * trafficB + 0.10 * techB));

  // Escort/salvage/bounty scan become more common in contested, pirate-heavy space.
  w.wEscort = std::max(0.0, w.wEscort * (1.0 + 0.55 * piracyB + 0.30 * contest + 0.18 * militarismB));
  w.wSalvage = std::max(0.0, w.wSalvage * (1.0 + 0.25 * piracyB + 0.25 * contest - 0.10 * stabilityB));
  w.wBountyScan = std::max(0.0, w.wBountyScan * (1.0 + 0.55 * piracyB + 0.35 * contest + 0.15 * militarismB));

  // Passenger traffic is suppressed by piracy pressure.
  w.wPassenger = std::max(0.0, w.wPassenger * (1.0 + 0.45 * trafficB - 0.25 * piracyB + 0.18 * wealthB));

  // Smuggling boards skew towards authoritarian/corrupt contexts (more contraband,
  // more opportunity, higher risk premiums).
  w.wSmuggle = std::max(0.0, w.wSmuggle * (1.0 + 0.30 * securityB + 0.25 * authorityB + 0.35 * corruptionB));

  // Keep the "remainder is BountyKill" semantic stable by capping the non-kill
  // weights. This prevents extreme systems from starving out BountyKill offers.
  const double cap = 0.92;
  const double sum = w.sum();
  if (sum > cap && sum > 1e-9) {
    const double s = cap / sum;
    w.wCourier *= s;
    w.wDelivery *= s;
    w.wMultiDelivery *= s;
    w.wEscort *= s;
    w.wSalvage *= s;
    w.wPassenger *= s;
    w.wSmuggle *= s;
    w.wBountyScan *= s;
  }

  return w;
}

void refreshMissionOffers(Universe& universe,
                          const StarSystem& currentSystem,
                          const Station& dockedStation,
                          double timeDays,
                          double playerRep,
                          SaveGame& ioSave,
                          const MissionBoardParams& params) {
  const int dayStamp = (int)std::floor(timeDays);
  if (ioSave.missionOffersStationId == dockedStation.id && ioSave.missionOffersDayStamp == dayStamp) return;

  ioSave.missionOffersStationId = dockedStation.id;
  ioSave.missionOffersDayStamp = dayStamp;
  ioSave.missionOffers.clear();

  core::SplitMix64 mrng(missionBoardSeed(dockedStation.id, dayStamp, dockedStation.factionId));

  auto candidates = universe.queryNearby(currentSystem.stub.posLy, params.searchRadiusLy, params.maxCandidates);
  std::vector<SystemStub> dests;
  dests.reserve(candidates.size());
  for (const auto& s : candidates) {
    if (s.id == currentSystem.stub.id) continue;
    if (s.stationCount <= 0) continue;
    dests.push_back(s);
  }

  const int baseCount = params.baseOfferCount
                      + (playerRep >= params.repThreshold1 ? 1 : 0)
                      + (playerRep >= params.repThreshold2 ? 1 : 0);
  const int offerCount = std::clamp(baseCount, params.minOfferCount, params.maxOfferCount);

  const double repScale = 1.0 + std::clamp(playerRep, 0.0, 100.0) / 100.0 * params.repRewardBonus;

  // Used for pricing delivery missions that require purchasing cargo.
  const double feeEff = applyReputationToFeeRate(dockedStation.feeRate, playerRep);

  // Cache origin economy snapshot for this refresh.
  auto& originEcon = universe.stationEconomy(dockedStation, timeDays);

  // Ship-fit: tailor offer payload sizes to the player's current ship capacity.
  // For cargo missions that grant cargo immediately on acceptance, we use the *free*
  // cargo space so the board avoids generating offers the player can't accept right now.
  const double cargoFreeKg = std::max(0.0, ioSave.cargoCapacityKg - econ::cargoMassKg(ioSave.cargo));
  const double cargoCapKg = std::max(0.0, ioSave.cargoCapacityKg);
  const int freeSeats = std::max(0, ioSave.passengerSeats - activePassengerCount(ioSave));

  // Context: local system security + issuing faction traits.
  //
  // This lets the board "feel" different in pirate-heavy frontier space vs.
  // stable core systems, without requiring any non-deterministic simulation.
  const SystemSecurityProfile localSec = systemSecurityProfile(universe.seed(), currentSystem);
  const FactionProfile issuerTraits = factionProfile(universe.seed(), dockedStation.factionId);

  // Precompute cumulative weight thresholds once (the RNG is deterministic per board seed).
  const MissionTypeWeights w = computeMissionTypeWeights(params, localSec, issuerTraits, playerRep);
  const double tCourier = w.wCourier;
  const double tDelivery = tCourier + w.wDelivery;
  const double tMulti = tDelivery + w.wMultiDelivery;
  const double tEscort = tMulti + w.wEscort;
  const double tSalvage = tEscort + w.wSalvage;
  const double tPassenger = tSalvage + w.wPassenger;
  const double tSmuggle = tPassenger + w.wSmuggle;
  const double tBountyScan = tSmuggle + w.wBountyScan;

  // Helpers for destination scoring/selection.
  const auto& factions = universe.factions();

  auto findFactionById = [&](core::u32 id) -> const Faction* {
    if (id == 0) return nullptr;
    for (const auto& f : factions) {
      if (f.id == id) return &f;
    }
    return nullptr;
  };

  auto relationScore = [&](core::u32 aId, core::u32 bId) -> double {
    if (aId == 0 || bId == 0) return 0.0;
    if (aId == bId) return 1.0;
    const Faction* a = findFactionById(aId);
    const Faction* b = findFactionById(bId);
    if (!a || !b) return 0.0;
    return std::clamp(factionRelation(universe.seed(), *a, *b), -1.0, 1.0);
  };

  auto scoreDestination = [&](MissionType type,
                              const SystemSecurityProfile& sp,
                              const FactionProfile& destTraits,
                              double rel,
                              bool hasContraband,
                              double distLy) -> double {
    const double rel01 = (rel + 1.0) * 0.5;          // 0..1
    const double hostility01 = std::max(0.0, -rel);  // 0..1
    const double dist01 = params.searchRadiusLy > 1e-6 ? std::clamp(distLy / params.searchRadiusLy, 0.0, 1.0) : 0.0;

    // "Legit" travel jobs prefer friendly/secure trade corridors.
    if (type == MissionType::Courier || type == MissionType::Delivery || type == MissionType::MultiDelivery || type == MissionType::Passenger) {
      double s = 0.58 * rel01
               + 0.18 * sp.traffic01
               + 0.10 * sp.security01
               + 0.10 * destTraits.wealth
               - 0.24 * sp.piracy01;

      // Avoid deep-hostile destinations for standard jobs.
      if (rel < -0.25 && sp.controllingFactionId != 0) s -= 0.45;

      // Slight preference for leaving the local cluster so the board doesn't
      // collapse into ultra-short hops.
      s += 0.04 * dist01;
      return s;
    }

    // Smuggling prefers strict, wealthy targets with actual contraband laws.
    if (type == MissionType::Smuggle) {
      if (!hasContraband) return -1e9;
      const double strict01 = std::max(0.0, sp.security01 - 0.5) * 2.0;
      double s = 0.28 * destTraits.wealth
               + 0.26 * destTraits.authority
               + 0.22 * strict01
               + 0.14 * (1.0 - destTraits.corruption)
               + 0.10 * hostility01
               - 0.08 * sp.piracy01;
      s += 0.03 * dist01;
      return s;
    }

    // Bounties prefer pirate-heavy / contested destinations.
    if (type == MissionType::BountyKill || type == MissionType::BountyScan) {
      double s = 0.56 * sp.piracy01
               + 0.25 * sp.contest01
               + 0.12 * (1.0 - sp.security01)
               + 0.07 * hostility01
               + 0.05 * destTraits.militarism;
      s += 0.03 * dist01;
      return s;
    }

    return 0.0;
  };

  auto pickDestination = [&](MissionType type,
                             SystemStub& outStub,
                             const Station*& outSt,
                             double& outDistLy,
                             SystemSecurityProfile& outSysSec,
                             double& outRel) -> bool {
    if (dests.empty()) return false;

    int tries = 14;
    if (type == MissionType::Smuggle) tries = 22;
    if (type == MissionType::BountyKill || type == MissionType::BountyScan) tries = 22;

    bool found = false;
    double bestScore = -1e9;
    SystemStub bestStub{};
    const Station* bestSt = nullptr;
    double bestDist = 0.0;
    SystemSecurityProfile bestSec{};
    double bestRel = 0.0;

    for (int t = 0; t < tries; ++t) {
      const auto& stub = dests[(std::size_t)(mrng.nextU32() % (core::u32)dests.size())];
      const auto& sys = universe.getSystem(stub.id, &stub);
      if (sys.stations.empty()) continue;

      const auto* st = &sys.stations[(std::size_t)(mrng.nextU32() % (core::u32)sys.stations.size())];
      const double distLy = (stub.posLy - currentSystem.stub.posLy).length();

      const SystemSecurityProfile sp = systemSecurityProfile(universe.seed(), sys);
      const FactionProfile stTraits = factionProfile(universe.seed(), st->factionId);

      const double rel = relationScore(dockedStation.factionId, st->factionId);

      bool hasContraband = true;
      if (type == MissionType::Smuggle) {
        const core::u32 mask = illegalCommodityMaskForStation(universe.seed(), st->factionId, st->id, st->type);
        hasContraband = (mask != 0u);
      }

      double score = scoreDestination(type, sp, stTraits, rel, hasContraband, distLy);
      // Add a tiny bit of noise so repeated offers don't always pick the same "best".
      score += mrng.range(-0.03, 0.03);

      if (!found || score > bestScore) {
        found = true;
        bestScore = score;
        bestStub = stub;
        bestSt = st;
        bestDist = distLy;
        bestSec = sp;
        bestRel = rel;
      }
    }

    if (!found || !bestSt) return false;
    outStub = bestStub;
    outSt = bestSt;
    outDistLy = bestDist;
    outSysSec = bestSec;
    outRel = bestRel;
    return true;
  };

  for (int i = 0; i < offerCount; ++i) {
    const double r = mrng.nextUnit();

    MissionType type = MissionType::BountyKill;
    if (r < tCourier) type = MissionType::Courier;
    else if (r < tDelivery) type = MissionType::Delivery;
    else if (r < tMulti) type = MissionType::MultiDelivery;
    else if (r < tEscort) type = MissionType::Escort;
    else if (r < tSalvage) type = MissionType::Salvage;
    else if (r < tPassenger) type = MissionType::Passenger;
    else if (r < tSmuggle) type = MissionType::Smuggle;
    else if (r < tBountyScan) type = MissionType::BountyScan;
    else type = MissionType::BountyKill;

    // Multi-hop deliveries need at least two distinct candidate systems.
    if (type == MissionType::MultiDelivery && dests.size() < 2) {
      type = MissionType::Delivery;
    }

    // If we don't have any destination candidates, we can still offer local Salvage/Escort jobs.
    if (type != MissionType::Salvage && type != MissionType::Escort && dests.empty()) {
      type = MissionType::Salvage;
    }

    Mission m{};
    m.id = 0; // offers are not persisted as "accepted" missions
    m.factionId = dockedStation.factionId;
    m.fromSystem = currentSystem.stub.id;
    m.fromStation = dockedStation.id;

    // Default to a local job that returns to the issuing station.
    m.toSystem = currentSystem.stub.id;
    m.toStation = dockedStation.id;

    // Default deadline (overridden per type)
    m.deadlineDay = timeDays + 1.0;

    double distLy = 0.0;
    const Station* destSt = nullptr;
    sim::SystemStub destStub{};

    // When the mission targets another system, we also capture the destination system's
    // security profile so rewards can incorporate risk/strictness.
    SystemSecurityProfile destSysSec{};
    double destRel = 0.0;
    bool haveDestSysSec = false;

    if (type != MissionType::Salvage && type != MissionType::Escort) {
      // Pick a destination that fits the mission "flavor" (trade corridors, pirate space, strict targets, ...).
      if (!pickDestination(type, destStub, destSt, distLy, destSysSec, destRel)) {
        // If we can't find any suitable destination, keep the board usable.
        type = MissionType::Salvage;
      } else {
        m.toSystem = destStub.id;
        m.toStation = destSt->id;
        m.deadlineDay = timeDays + 1.0 + distLy / 20.0;
        haveDestSysSec = true;
      }
    }

    if (type == MissionType::Salvage) {
      // Salvage is a local in-system job: make it snappy.
      m.deadlineDay = timeDays + 0.35 + mrng.range(0.05, 0.20);
    } else if (type == MissionType::Escort) {
      // Escort is a local in-system job that moves between *stations* in the current system.
      if (currentSystem.stations.size() < 2) {
        type = MissionType::Salvage;
        m.deadlineDay = timeDays + 0.35 + mrng.range(0.05, 0.20);
      } else {
        // Pick a destination station different from the issuing station.
        const Station* picked = nullptr;
        for (int tries = 0; tries < 8; ++tries) {
          const auto& cand = currentSystem.stations[(std::size_t)(mrng.nextU32() % (core::u32)currentSystem.stations.size())];
          if (cand.id != dockedStation.id) { picked = &cand; break; }
        }
        if (!picked) {
          // Deterministic fallback: use the first non-issuing station.
          for (const auto& cand : currentSystem.stations) {
            if (cand.id != dockedStation.id) { picked = &cand; break; }
          }
        }

        if (!picked) {
          type = MissionType::Salvage;
          m.deadlineDay = timeDays + 0.35 + mrng.range(0.05, 0.20);
        } else {
          destStub = currentSystem.stub;
          destSt = picked;
          m.toSystem = currentSystem.stub.id;
          m.toStation = destSt->id;

          // In-system escorts should be short compared to interstellar courier jobs.
          // Scale the deadline loosely with orbital separation in AU.
          const double distAU = (orbitPosition3DAU(dockedStation.orbit, timeDays) - orbitPosition3DAU(destSt->orbit, timeDays)).length();
          m.deadlineDay = timeDays + 0.35 + std::clamp(distAU * 0.10, 0.0, 0.65) + mrng.range(0.02, 0.12);
        }
      }
    } else {
      // Cross-system mission should have a resolved destination station.
      if (!destSt) continue;
    }

    if (type == MissionType::Courier) {
      m.type = MissionType::Courier;
      m.reward = 350.0 + distLy * 110.0;
    } else if (type == MissionType::Delivery) {
      m.type = MissionType::Delivery;
      m.cargoProvided = (mrng.nextUnit() < 0.25);

      // Ship-fit: cargo missions grant cargo immediately on acceptance, so only generate
      // jobs that fit into the player's current free cargo space.
      const double freeKg = cargoFreeKg;
      if (freeKg + 1e-9 < 0.2) {
        // Hold is effectively full (or too small to carry any commodity); keep the board usable.
        m.type = MissionType::Courier;
        m.reward = 330.0 + distLy * 115.0;
      } else {
        const double maxCargoKg = std::min(260.0, freeKg * 0.95);
        const int minUnits = std::clamp((int)std::floor(freeKg / 30.0), 1, 12);
        const int maxUnitsHard = std::clamp((int)std::floor(freeKg / 1.2), 6, 180);

        // Economy-aware cargo selection. Ensure we don't ask the player to deliver contraband
        // as a "legal" delivery job.
        CargoMissionSpec spec{};
        bool ok = false;
        if (destSt) {
          auto& destEcon = universe.stationEconomy(*destSt, timeDays);
          const core::u32 destIllegalMask =
              illegalCommodityMaskForStation(universe.seed(), destSt->factionId, destSt->id, destSt->type);
          ok = pickCargoMission(mrng,
                                originEcon,
                                dockedStation.economyModel,
                                destEcon,
                                destSt->economyModel,
                                destIllegalMask,
                                /*requireLegalAtDest=*/true,
                                /*cargoProvided=*/m.cargoProvided,
                                /*maxCargoKg=*/maxCargoKg,
                                /*minUnits=*/minUnits,
                                /*maxUnitsHard=*/maxUnitsHard,
                                spec);
        }

        if (!ok) {
          // Fallback: if the economy can't support a cargo job right now, degrade to courier.
          m.type = MissionType::Courier;
          m.reward = 330.0 + distLy * 115.0;
        } else {
          m.commodity = spec.commodity;
          m.units = (double)spec.units;

          const double base = 250.0 + distLy * 58.0 + mrng.range(0.0, 140.0);

          if (m.cargoProvided) {
            // Provided cargo: payout is more "service" oriented; scale with destination price.
            m.reward = base + spec.destQ.bid * (double)spec.units * 0.55;
          } else {
            // Purchased cargo: ensure missions are meaningfully profitable even for high-value goods.
            const double cost = spec.originQ.ask * (double)spec.units * (1.0 + feeEff);
            const double margin = 0.18 + mrng.range(0.0, 0.08);
            m.reward = base + cost * (1.0 + margin);
          }
        }
      }
    } else if (type == MissionType::MultiDelivery) {
      m.type = MissionType::MultiDelivery;
      m.cargoProvided = (mrng.nextUnit() < 0.35);

      // Multi-hop deliveries grant cargo immediately on acceptance; size them to the player's
      // current free cargo capacity to avoid impossible offers.
      const double freeKg = cargoFreeKg;
      if (freeKg + 1e-9 < 0.2) {
        m.type = MissionType::Courier;
        m.reward = 360.0 + distLy * 120.0;
      } else {
        const double maxCargoKg = std::min(240.0, freeKg * 0.90);
        const int minUnits = std::clamp((int)std::floor(freeKg / 35.0), 1, 10);
        const int maxUnitsHard = std::clamp((int)std::floor(freeKg / 1.4), 6, 160);

        // Pick a cargo load that makes sense for the final destination's market.
        CargoMissionSpec spec{};
        bool ok = false;
        if (destSt) {
          auto& destEcon = universe.stationEconomy(*destSt, timeDays);
          const core::u32 destIllegalMask =
              illegalCommodityMaskForStation(universe.seed(), destSt->factionId, destSt->id, destSt->type);
          ok = pickCargoMission(mrng,
                                originEcon,
                                dockedStation.economyModel,
                                destEcon,
                                destSt->economyModel,
                                destIllegalMask,
                                /*requireLegalAtDest=*/true,
                                /*cargoProvided=*/m.cargoProvided,
                                /*maxCargoKg=*/maxCargoKg,
                                /*minUnits=*/minUnits,
                                /*maxUnitsHard=*/maxUnitsHard,
                                spec);
        }

        if (!ok) {
          // Fallback: keep the board populated.
          m.type = MissionType::Courier;
          m.reward = 360.0 + distLy * 120.0;
        } else {
          m.commodity = spec.commodity;
          m.units = (double)spec.units;

          // Pick a via system/station. Prefer routes that aren't extreme detours.
          SystemStub viaStub{};
          const Station* viaSt = nullptr;
          double routeLy = 0.0;

          const double directLy = std::max(1e-6, distLy);
          for (int tries = 0; tries < 18; ++tries) {
            const auto candStub = dests[(std::size_t)(mrng.nextU32() % (core::u32)dests.size())];
            if (candStub.id == destStub.id) continue;

            const double d1 = (candStub.posLy - currentSystem.stub.posLy).length();
            const double d2 = (destStub.posLy - candStub.posLy).length();
            const double route = d1 + d2;
            const double detour = route / directLy;
            if (detour > 2.25) continue;

            const auto& candSys = universe.getSystem(candStub.id, &candStub);
            if (candSys.stations.empty()) continue;
            const auto& candSt = candSys.stations[(std::size_t)(mrng.nextU32() % (core::u32)candSys.stations.size())];

            if (!viaSt || route < routeLy) {
              viaStub = candStub;
              viaSt = &candSt;
              routeLy = route;
            }
          }

          if (!viaSt) {
            // Couldn't find a reasonable intermediate stop; downgrade to a normal delivery.
            m.type = MissionType::Delivery;

            const double base = 270.0 + distLy * 60.0 + mrng.range(0.0, 160.0);
            if (m.cargoProvided) {
              m.reward = base + spec.destQ.bid * (double)spec.units * 0.55;
            } else {
              const double cost = spec.originQ.ask * (double)spec.units * (1.0 + feeEff);
              const double margin = 0.18 + mrng.range(0.0, 0.08);
              m.reward = base + cost * (1.0 + margin);
            }

            // Keep the direct-leg deadline picked above.
          } else {
            m.viaSystem = viaStub.id;
            m.viaStation = viaSt->id;

            // Reward/deadline scale with the actual route distance (origin->via + via->dest),
            // plus a small "extra stop" premium.
            const double base = 380.0 + routeLy * 70.0 + mrng.range(0.0, 220.0);

            if (m.cargoProvided) {
              m.reward = base + spec.destQ.bid * (double)spec.units * 0.65;
            } else {
              const double cost = spec.originQ.ask * (double)spec.units * (1.0 + feeEff);
              const double margin = 0.22 + mrng.range(0.0, 0.10);
              m.reward = base + cost * (1.0 + margin);
            }

            m.reward += 140.0 + mrng.range(0.0, 70.0);
            m.deadlineDay = timeDays + 2.05 + routeLy / 18.0 + mrng.range(0.0, 0.25);
          }
        }
      }
    } else if (type == MissionType::Escort) {
      m.type = MissionType::Escort;

      if (!destSt || destSt->id == dockedStation.id) {
        // If we can't resolve a destination station, keep the board usable.
        m.type = MissionType::Courier;
        m.reward = 360.0;
      } else {
        // Convoy payload: reserve real station inventory at accept time.
        m.cargoProvided = true;

        // Stable convoy id for the game layer.
        m.targetNpcId = (mrng.nextU64() | 0x4000000000000000ull);

        // Pick a legal cargo load that is in demand at the destination station.
        CargoMissionSpec spec{};
        bool ok = false;
        {
          auto& destEcon = universe.stationEconomy(*destSt, timeDays);

          // Use a fixed NPC "freighter" capacity. This is separate from the player's ship hold.
          const double maxCargoKg = 320.0;
          const int minUnits = 6;
          const int maxUnitsHard = 220;

          const core::u32 destIllegalMask =
              illegalCommodityMaskForStation(universe.seed(), destSt->factionId, destSt->id, destSt->type);
          ok = pickCargoMission(mrng,
                                originEcon,
                                dockedStation.economyModel,
                                destEcon,
                                destSt->economyModel,
                                destIllegalMask,
                                /*requireLegalAtDest=*/true,
                                /*cargoProvided=*/true,
                                /*maxCargoKg=*/maxCargoKg,
                                /*minUnits=*/minUnits,
                                /*maxUnitsHard=*/maxUnitsHard,
                                spec);
        }

        if (!ok) {
          // If the economy can't support a convoy job right now, downgrade to courier.
          m.type = MissionType::Courier;
          m.reward = 360.0;
          m.cargoProvided = false;
          m.targetNpcId = 0;
        } else {
          m.commodity = spec.commodity;
          m.units = (double)spec.units;

          const double distAU = (orbitPosition3DAU(dockedStation.orbit, timeDays) - orbitPosition3DAU(destSt->orbit, timeDays)).length();

          // Reward: distance service component + small share of cargo value.
          const double cargoValue = spec.originQ.ask * (double)spec.units;
          const double base = 780.0 + distAU * 210.0 + mrng.range(0.0, 160.0);
          m.reward = base + cargoValue * (0.10 + mrng.range(0.0, 0.04));

          // Deadline already set in the local-destination selection above.
        }
      }
    } else if (type == MissionType::Salvage) {
      m.type = MissionType::Salvage;

      // Salvage jobs are local: recover loose goods from a mission derelict signal and return here.
      static const std::array<econ::CommodityId, 4> kSalvage = {
        econ::CommodityId::Machinery,
        econ::CommodityId::Electronics,
        econ::CommodityId::Metals,
        econ::CommodityId::Luxury,
      };

      // Ensure the requested salvage can fit in the player's cargo hold (at completion time).
      const double capKg = cargoCapKg;
      bool ok = false;

      for (int tries = 0; tries < 8; ++tries) {
        m.commodity = kSalvage[(std::size_t)(mrng.nextU32() % (core::u32)kSalvage.size())];
        const double massKg = std::max(1e-6, econ::commodityDef(m.commodity).massKg);
        const int maxByKg = (int)std::floor(capKg / massKg + 1e-9);
        if (maxByKg >= 1) {
          const int maxUnits = std::clamp(maxByKg, 1, 18);
          const int minUnits = std::min(3, maxUnits);
          m.units = (double)mrng.range(minUnits, maxUnits);
          ok = true;
          break;
        }
      }

      if (!ok) {
        // If the ship can't carry any salvage at all, keep the board usable.
        m.type = MissionType::Courier;
        m.reward = 360.0;
      } else {
        const double base = econ::commodityDef(m.commodity).basePrice;
        m.reward = 520.0 + (double)m.units * base * 1.05 + mrng.range(0.0, 180.0);

        // Reuse targetNpcId as a stable "signal id" so the game can spawn/track the mission site.
        m.targetNpcId = (mrng.nextU64() | 0x8000000000000000ull);
        m.cargoProvided = false;
      }
    } else if (type == MissionType::Passenger) {
      // Ship-fit: only offer passenger party sizes that the player can accept right now.
      const int maxParty = std::min(6, freeSeats);
      if (maxParty <= 0) {
        m.type = MissionType::Courier;
        m.reward = 340.0 + distLy * 118.0;
      } else {
        m.type = MissionType::Passenger;
        // units is interpreted as "passenger count" for Passenger missions.
        m.units = (double)mrng.range(1, maxParty);
        // Reward scales with distance and party size.
        m.reward = 450.0 + distLy * 125.0 + (double)m.units * 85.0;
        // Slightly tighter deadline than courier to encourage routing decisions.
        m.deadlineDay = timeDays + 0.9 + distLy / 24.0;
        m.cargoProvided = false;
      }
    } else if (type == MissionType::Smuggle) {
      m.type = MissionType::Smuggle;

      // Smuggling jobs should only exist where there is actual contraband.
      const core::u32 destFaction = destSt ? destSt->factionId : 0;
      const core::u32 mask =
          destSt ? illegalCommodityMaskForStation(universe.seed(), destSt->factionId, destSt->id, destSt->type) : 0u;

      std::vector<econ::CommodityId> illegal;
      illegal.reserve(econ::kCommodityCount);
      for (std::size_t c = 0; c < econ::kCommodityCount; ++c) {
        if ((mask & ((core::u32)1u << (core::u32)c)) != 0u) illegal.push_back((econ::CommodityId)c);
      }

      if (destFaction == 0 || illegal.empty()) {
        // No contraband rules here; downgrade to a normal courier job.
        m.type = MissionType::Courier;
        m.reward = 340.0 + distLy * 118.0;
      } else {
        // Ship-fit: smuggling grants cargo immediately on acceptance, so filter to commodities that
        // can actually fit in the player's *current* free hold.
        const double freeKg = cargoFreeKg;

        std::vector<econ::CommodityId> fitIllegal;
        fitIllegal.reserve(illegal.size());
        for (const auto cid : illegal) {
          const double massKg = std::max(1e-6, econ::commodityDef(cid).massKg);
          if (massKg <= freeKg * 0.95 + 1e-9) fitIllegal.push_back(cid);
        }

        if (fitIllegal.empty()) {
          // Not enough free capacity to accept even 1 unit of any contraband.
          m.type = MissionType::Courier;
          m.reward = 340.0 + distLy * 118.0;
        } else {
          m.commodity = fitIllegal[(std::size_t)(mrng.nextU32() % (core::u32)fitIllegal.size())];

          const double massKg = std::max(1e-6, econ::commodityDef(m.commodity).massKg);
          const int maxByKg = (int)std::floor((freeKg * 0.95) / massKg + 1e-9);
          const int maxUnits = std::clamp(maxByKg, 1, 55);

          const int minUnits = std::min(8, maxUnits);
          const int units = mrng.range(minUnits, maxUnits);
          m.units = (double)units;

          // Smuggling cargo is provided by the contact (not drawn from open-market inventory).
          m.cargoProvided = true;

          // Reward scales with destination price (risk/interest) + distance.
          double destMid = econ::commodityDef(m.commodity).basePrice;
          if (destSt) {
            auto& destEcon = universe.stationEconomy(*destSt, timeDays);
            destMid = econ::quote(destEcon, destSt->economyModel, m.commodity, 0.10).mid;
          }

          m.reward = 720.0
                   + distLy * 170.0
                   + destMid * (double)units * 0.28
                   + mrng.range(0.0, 220.0);
          m.deadlineDay = timeDays + 1.15 + distLy / 20.0;
        }
      }
    } else if (type == MissionType::BountyScan) {
      m.type = MissionType::BountyScan;
      m.targetNpcId = std::max<core::u64>(1, mrng.nextU64());
      m.reward = 900.0 + distLy * 80.0;
      m.deadlineDay = timeDays + 1.5 + distLy / 22.0;
    } else {
      m.type = MissionType::BountyKill;
      m.targetNpcId = std::max<core::u64>(1, mrng.nextU64());
      m.reward = 1400.0 + distLy * 90.0;
      m.deadlineDay = timeDays + 1.8 + distLy / 20.0;
    }

    // Risk/strictness bonus: tie payout to the same security model used by
    // destination selection, so "dangerous" or tightly policed space pays more.
    {
      const SystemSecurityProfile& sec = haveDestSysSec ? destSysSec : localSec;

      const double piracyRisk01 = std::max(0.0, sec.piracy01 - 0.5) * 2.0; // 0..1
      const double contest01 = std::clamp(sec.contest01, 0.0, 1.0);
      const double strict01 = std::max(0.0, sec.security01 - 0.5) * 2.0; // 0..1

      // Prefer destination faction traits when available.
      const FactionProfile destTraits = destSt ? factionProfile(universe.seed(), destSt->factionId) : issuerTraits;

      double mult = 1.0;
      switch (m.type) {
        case MissionType::Courier:
          mult = 1.0 + 0.10 * piracyRisk01 + 0.08 * contest01;
          break;
        case MissionType::Passenger:
          mult = 1.0 + 0.08 * piracyRisk01 + 0.06 * contest01;
          break;
        case MissionType::Delivery:
        case MissionType::MultiDelivery:
          mult = 1.0 + 0.12 * piracyRisk01 + 0.10 * contest01;
          break;
        case MissionType::Escort:
          mult = 1.0 + 0.18 * piracyRisk01 + 0.14 * contest01;
          break;
        case MissionType::Salvage:
          mult = 1.0 + 0.16 * piracyRisk01 + 0.18 * contest01;
          break;
        case MissionType::Smuggle:
          // Smuggling pays mostly for strictness and authority, reduced by corruption.
          mult = 1.0
               + 0.22 * strict01
               + 0.10 * std::clamp(destTraits.authority, 0.0, 1.0)
               + 0.08 * (1.0 - std::clamp(destTraits.corruption, 0.0, 1.0));
          break;
        case MissionType::BountyScan:
        case MissionType::BountyKill:
          mult = 1.0
               + 0.28 * piracyRisk01
               + 0.20 * contest01
               + 0.10 * (1.0 - std::clamp(sec.security01, 0.0, 1.0));
          break;
      }

      m.reward *= std::clamp(mult, 1.0, 1.65);
    }

    m.reward *= repScale;
    ioSave.missionOffers.push_back(std::move(m));
  }
}

bool acceptMissionOffer(Universe& universe,
                        const Station& dockedStation,
                        double timeDays,
                        double playerRep,
                        SaveGame& ioSave,
                        std::size_t offerIndex,
                        std::string* outError) {
  if (offerIndex >= ioSave.missionOffers.size()) {
    if (outError) *outError = "Offer index out of range.";
    return false;
  }

  Mission m = ioSave.missionOffers[offerIndex];

  const double feeEff = applyReputationToFeeRate(dockedStation.feeRate, playerRep);
  auto& stEcon = universe.stationEconomy(dockedStation, timeDays);

  // Delivery-like missions need cargo acquisition.
  const bool isDelivery = (m.type == MissionType::Delivery || m.type == MissionType::MultiDelivery || m.type == MissionType::Smuggle);
  if (isDelivery && m.units > 0.0) {
    const econ::CommodityId cid = m.commodity;
    const double massKg = std::max(1e-6, econ::commodityDef(cid).massKg);
    const double needKg = massKg * (double)m.units;
    if (econ::cargoMassKg(ioSave.cargo) + needKg > ioSave.cargoCapacityKg + 1e-6) {
      if (outError) *outError = "Not enough cargo capacity for this delivery.";
      return false;
    }

    if (m.type == MissionType::Smuggle) {
      // Smuggling cargo is provided by the job contact (not pulled from the legal station market).
      ioSave.cargo[(std::size_t)cid] += (double)m.units;
    } else if (m.cargoProvided) {
      // Provided cargo should come from real station inventory.
      const double taken = econ::takeInventory(stEcon, dockedStation.economyModel, cid, (double)m.units);
      if (taken + 1e-6 < (double)m.units) {
        // Restore any partial transfer.
        econ::addInventory(stEcon, dockedStation.economyModel, cid, taken);
        if (outError) *outError = "Station can't provide that much cargo right now.";
        return false;
      }
      ioSave.cargo[(std::size_t)cid] += taken;
    } else {
      const auto q = econ::quote(stEcon, dockedStation.economyModel, cid, 0.10);
      const double totalCost = q.ask * (double)m.units * (1.0 + feeEff);
      if (q.inventory + 1e-6 < (double)m.units) {
        if (outError) *outError = "Station doesn't have enough inventory to stock this delivery.";
        return false;
      }
      if (ioSave.credits + 1e-6 < totalCost) {
        if (outError) *outError = "Not enough credits to buy the delivery cargo.";
        return false;
      }
      const auto tr = econ::buy(stEcon, dockedStation.economyModel, cid, (double)m.units, ioSave.credits, 0.10, feeEff);
      if (!tr.ok) {
        if (outError) *outError = tr.reason ? tr.reason : "Trade failed.";
        return false;
      }
      ioSave.cargo[(std::size_t)cid] += tr.unitsDelta;
    }
  }

  // Escort missions reserve station inventory for an NPC convoy. The player doesn't receive cargo,
  // but the origin market should still be debited so the economy reflects that the convoy was
  // loaded.
  if (m.type == MissionType::Escort && m.units > 0.0) {
    const econ::CommodityId cid = m.commodity;
    const double taken = econ::takeInventory(stEcon, dockedStation.economyModel, cid, (double)m.units);
    if (taken + 1e-6 < (double)m.units) {
      // Restore any partial transfer.
      econ::addInventory(stEcon, dockedStation.economyModel, cid, taken);
      if (outError) *outError = "Station can't load that convoy right now.";
      return false;
    }
    m.cargoProvided = true;
  }

  // Salvage missions don't grant cargo up-front, but we should still ensure the requested salvage will
  // fit in the player's hold (otherwise the job is impossible to complete).
  if (m.type == MissionType::Salvage && m.units > 0.0) {
    const econ::CommodityId cid = m.commodity;
    const double massKg = std::max(1e-6, econ::commodityDef(cid).massKg);
    const double needKg = massKg * (double)m.units;
    if (needKg > ioSave.cargoCapacityKg + 1e-6) {
      if (outError) *outError = "This salvage won't fit in your cargo hold.";
      return false;
    }
  }

  // Passenger missions require cabin capacity.
  if (m.type == MissionType::Passenger) {
    const int party = std::max(0, (int)std::llround(m.units));
    const int used = activePassengerCount(ioSave);
    const int cap = std::max(0, ioSave.passengerSeats);
    if (party <= 0) {
      if (outError) *outError = "Passenger party size is invalid.";
      return false;
    }
    if (used + party > cap) {
      if (outError) *outError = "Not enough passenger seats for this job.";
      return false;
    }
  }

  m.id = ioSave.nextMissionId++;
  ioSave.missions.push_back(m);
  ioSave.missionOffers.erase(ioSave.missionOffers.begin() + (std::ptrdiff_t)offerIndex);
  return true;
}

MissionTickResult tickMissionDeadlines(SaveGame& ioSave, double timeDays, double repPenaltyOnFail) {
  MissionTickResult r{};
  for (auto& m : ioSave.missions) {
    if (m.completed || m.failed) continue;
    if (m.deadlineDay > 0.0 && timeDays > m.deadlineDay) {
      m.failed = true;
      addReputation(ioSave, m.factionId, repPenaltyOnFail);
      ++r.failed;
    }
  }
  return r;
}

MissionDockResult tryCompleteMissionsAtDock(Universe& universe,
                                            const StarSystem& currentSystem,
                                            const Station& dockedStation,
                                            double timeDays,
                                            SaveGame& ioSave,
                                            double repRewardOnComplete) {
  MissionDockResult r{};
  const SystemId sysId = currentSystem.stub.id;
  const StationId here = dockedStation.id;

  for (auto& m : ioSave.missions) {
    if (m.completed || m.failed) continue;

    const bool atFinal = (sysId == m.toSystem && here == m.toStation);
    const bool atVia = (m.viaSystem != 0 && sysId == m.viaSystem && here == m.viaStation);

    if (m.type == MissionType::Courier) {
      if (atFinal) {
        m.completed = true;
        ioSave.credits += m.reward;
        addReputation(ioSave, m.factionId, repRewardOnComplete);
        ++r.completed;
      }
    } else if (m.type == MissionType::Passenger) {
      if (atFinal) {
        m.completed = true;
        ioSave.credits += m.reward;
        addReputation(ioSave, m.factionId, repRewardOnComplete);
        ++r.completed;
      }
    } else if (m.type == MissionType::Escort) {
      // Escort jobs are completed by docking at the destination after the convoy arrives.
      // The game layer sets m.scanned=true when the NPC convoy reaches the target station.
      if (atFinal && m.scanned) {
        if (m.units > 0.0) {
          auto& econHere = universe.stationEconomy(dockedStation, timeDays);
          econ::addInventory(econHere, dockedStation.economyModel, m.commodity, m.units);
        }

        m.completed = true;
        ioSave.credits += m.reward;
        addReputation(ioSave, m.factionId, repRewardOnComplete);
        ++r.completed;
      }
    } else if (m.type == MissionType::Salvage) {
      // Salvage jobs require that the player has actually visited the mission site (m.scanned),
      // then returned with the requested goods.
      if (atFinal && m.scanned) {
        const econ::CommodityId cid = m.commodity;
        const double have = ioSave.cargo[(std::size_t)cid];
        if (have + 1e-6 >= m.units) {
          ioSave.cargo[(std::size_t)cid] -= m.units;
          if (ioSave.cargo[(std::size_t)cid] < 1e-6) ioSave.cargo[(std::size_t)cid] = 0.0;

          // Feed recovered salvage into the local market (best-effort).
          auto& econHere = universe.stationEconomy(dockedStation, timeDays);
          econ::addInventory(econHere, dockedStation.economyModel, cid, m.units);

          m.completed = true;
          ioSave.credits += m.reward;
          addReputation(ioSave, m.factionId, repRewardOnComplete);
          ++r.completed;
        }
      }
    } else if (m.type == MissionType::Delivery || m.type == MissionType::MultiDelivery || m.type == MissionType::Smuggle) {
      if (m.type == MissionType::MultiDelivery && m.viaSystem != 0 && m.leg == 0 && atVia) {
        m.leg = 1;
        ++r.progressedMultiLeg;
      } else if (atFinal && (m.viaSystem == 0 || m.leg >= 1)) {
        const econ::CommodityId cid = m.commodity;
        const double have = ioSave.cargo[(std::size_t)cid];
        if (have + 1e-6 >= m.units) {
          ioSave.cargo[(std::size_t)cid] -= m.units;
          if (ioSave.cargo[(std::size_t)cid] < 1e-6) ioSave.cargo[(std::size_t)cid] = 0.0;

          if (m.type != MissionType::Smuggle) {
            // Feed delivered cargo back into the destination market (best-effort).
            auto& econHere = universe.stationEconomy(dockedStation, timeDays);
            econ::addInventory(econHere, dockedStation.economyModel, cid, m.units);
          }

          m.completed = true;
          ioSave.credits += m.reward;
          addReputation(ioSave, m.factionId, repRewardOnComplete);
          ++r.completed;
        }
      }
    }
  }

  return r;
}


MissionEventResult tryCompleteBountyScan(SaveGame& ioSave,
                                         SystemId currentSystemId,
                                         core::u64 targetNpcId,
                                         double repRewardOnComplete) {
  MissionEventResult r{};
  if (currentSystemId == 0 || targetNpcId == 0) return r;

  for (auto& m : ioSave.missions) {
    if (m.completed || m.failed) continue;
    if (m.type != MissionType::BountyScan) continue;
    if (m.toSystem != currentSystemId) continue;
    if (m.targetNpcId != targetNpcId) continue;

    m.scanned = true;
    m.completed = true;
    ioSave.credits += m.reward;
    addReputation(ioSave, m.factionId, repRewardOnComplete);
    ++r.completed;
    r.rewardCr += m.reward;
  }

  return r;
}

MissionEventResult tryCompleteBountyKill(SaveGame& ioSave,
                                         SystemId currentSystemId,
                                         core::u64 targetNpcId,
                                         double repRewardOnComplete) {
  MissionEventResult r{};
  if (currentSystemId == 0 || targetNpcId == 0) return r;

  for (auto& m : ioSave.missions) {
    if (m.completed || m.failed) continue;
    if (m.type != MissionType::BountyKill) continue;
    if (m.toSystem != currentSystemId) continue;
    if (m.targetNpcId != targetNpcId) continue;

    m.completed = true;
    ioSave.credits += m.reward;
    addReputation(ioSave, m.factionId, repRewardOnComplete);
    ++r.completed;
    r.rewardCr += m.reward;
  }

  return r;
}

} // namespace stellar::sim
