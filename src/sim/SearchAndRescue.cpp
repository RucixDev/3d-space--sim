#include "stellar/sim/SearchAndRescue.h"

#include "stellar/sim/Mission.h"
#include "stellar/sim/Reputation.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace stellar::sim {

static bool podExpired(const RescuedPod& p, double timeDays) {
  if (!std::isfinite(timeDays)) return false;
  if (!std::isfinite(p.lifeSupportEndDay) || p.lifeSupportEndDay <= 0.0) return false;
  return timeDays >= p.lifeSupportEndDay;
}

static double urgencyFactor01(const RescuedPod& p, double timeDays) {
  // Returns 0..1-ish where 1 means freshly recovered and 0 means at/after expiry.
  if (!std::isfinite(timeDays)) return 0.5;
  if (!std::isfinite(p.lifeSupportEndDay) || p.lifeSupportEndDay <= 0.0) return 0.5;
  if (!std::isfinite(p.recoveredDay)) return 0.5;

  const double duration = std::max(1e-6, p.lifeSupportEndDay - p.recoveredDay);
  const double remain = p.lifeSupportEndDay - timeDays;
  const double frac = remain / duration;
  return std::clamp(frac, 0.0, 1.0);
}

RescuePodSummary summarizeRescuedPods(const SaveGame& s, double timeDays) {
  RescuePodSummary sum{};
  sum.total = (int)s.rescuedPods.size();
  double soonest = 0.0;

  for (const auto& p : s.rescuedPods) {
    const bool expired = podExpired(p, timeDays);
    const bool rewardable = !p.fromPlayerKill;

    if (!rewardable) ++sum.noReward;

    if (expired) {
      ++sum.expired;
      continue;
    }

    ++sum.alive;

    if (std::isfinite(p.lifeSupportEndDay) && p.lifeSupportEndDay > 0.0) {
      if (soonest == 0.0 || p.lifeSupportEndDay < soonest) soonest = p.lifeSupportEndDay;
    }
  }

  sum.soonestLifeSupportEndDay = soonest;
  return sum;
}

int passengerMissionSeatUsage(const SaveGame& s) {
  int used = 0;
  for (const auto& m : s.missions) {
    if (m.completed || m.failed) continue;
    if (m.type != MissionType::Passenger) continue;
    used += std::max(0, (int)std::llround(m.units));
  }
  return used;
}

int rescuedPodSeatUsage(const SaveGame& s) {
  return (int)s.rescuedPods.size();
}

int occupiedPassengerSeatsTotal(const SaveGame& s) {
  return passengerMissionSeatUsage(s) + rescuedPodSeatUsage(s);
}

int freePassengerSeats(const SaveGame& s) {
  const int cap = std::max(0, s.passengerSeats);
  const int used = std::max(0, occupiedPassengerSeatsTotal(s));
  return std::max(0, cap - used);
}

bool addRescuedPod(SaveGame& ioSave, const RescuedPod& pod, double timeDays, const char** outReason) {
  if (outReason) *outReason = nullptr;

  if (pod.id == 0) {
    if (outReason) *outReason = "Invalid pod id.";
    return false;
  }

  for (const auto& existing : ioSave.rescuedPods) {
    if (existing.id == pod.id) {
      if (outReason) *outReason = "Pod already onboard.";
      return false;
    }
  }

  if (freePassengerSeats(ioSave) <= 0) {
    if (outReason) *outReason = "No free passenger seats.";
    return false;
  }

  RescuedPod p = pod;

  if (!std::isfinite(p.recoveredDay) || p.recoveredDay <= 0.0) {
    p.recoveredDay = std::isfinite(timeDays) ? timeDays : 0.0;
  }

  if (!std::isfinite(p.lifeSupportEndDay) || p.lifeSupportEndDay <= 0.0) {
    p.lifeSupportEndDay = p.recoveredDay + kRescuePodDefaultLifeSupportDays;
  }

  if (p.lifeSupportEndDay < p.recoveredDay) {
    p.lifeSupportEndDay = p.recoveredDay;
  }

  ioSave.rescuedPods.push_back(p);
  return true;
}

static std::size_t getOrCreateBreakdown(std::vector<RescuePodTurnInBreakdown>& v,
                                       std::unordered_map<core::u32, std::size_t>& idx,
                                       core::u32 factionId) {
  auto it = idx.find(factionId);
  if (it != idx.end()) return it->second;
  const std::size_t i = v.size();
  v.push_back(RescuePodTurnInBreakdown{});
  v[i].factionId = factionId;
  idx[factionId] = i;
  return i;
}

RescuePodTurnInQuote quoteRescuePodTurnIn(const SaveGame& s,
                                         double timeDays,
                                         core::u32 stationFactionId,
                                         bool stationHasAuthority) {
  RescuePodTurnInQuote q{};
  q.totalPods = (int)s.rescuedPods.size();

  if (q.totalPods <= 0) {
    q.ok = false;
    q.reason = "No rescued pods onboard.";
    return q;
  }

  q.ok = true;

  // If there's no authority here, we still allow disembarking the pods but pay nothing.
  if (!stationHasAuthority) {
    q.reason = "No local authority: pods will be disembarked without reward.";
    return q;
  }

  std::unordered_map<core::u32, std::size_t> idx;
  idx.reserve(std::min<std::size_t>(s.rescuedPods.size(), 32));

  for (const auto& p : s.rescuedPods) {
    const bool expired = podExpired(p, timeDays);
    const bool rewardable = !p.fromPlayerKill;

    const core::u32 effFaction = (p.registryFactionId != 0) ? p.registryFactionId : stationFactionId;
    const bool factionMatch = (effFaction == 0) || (stationFactionId != 0 && stationFactionId == effFaction);

    auto& b = q.byFaction[getOrCreateBreakdown(q.byFaction, idx, effFaction)];
    ++b.pods;

    if (expired) {
      ++b.expired;
      ++q.expiredPods;
    } else {
      ++b.alive;
    }

    if (!rewardable) {
      ++b.noReward;
      ++q.noRewardPods;
      continue;
    }

    ++q.rewardablePods;

    const double frac = urgencyFactor01(p, timeDays);

    // Small urgency bonus for quick recoveries (bounded).
    const double urgencyCredits = 0.85 + 0.30 * frac; // [0.85, 1.15]
    const double urgencyRep = 0.80 + 0.40 * frac;     // [0.80, 1.20]

    double credits = kRescuePodBasePayoutCr * urgencyCredits;
    double rep = kRescuePodBaseRep * urgencyRep;

    if (expired) {
      credits *= kRescuePodExpiredPayoutMul;
      rep *= kRescuePodExpiredRepMul;
    }

    if (!factionMatch) {
      credits *= kRescuePodOffFactionPayoutMul;
      rep *= kRescuePodOffFactionRepMul;
    }

    if (effFaction == 0) {
      rep = 0.0;
    }

    if (!std::isfinite(credits) || credits < 0.0) credits = 0.0;
    if (!std::isfinite(rep) || rep < 0.0) rep = 0.0;

    b.credits += credits;
    b.rep += rep;

    q.creditsTotal += credits;
    q.repTotal += rep;
  }

  return q;
}

RescuePodTurnInResult applyRescuePodTurnIn(SaveGame& ioSave,
                                          double timeDays,
                                          core::u32 stationFactionId,
                                          bool stationHasAuthority) {
  RescuePodTurnInResult r{};
  RescuePodTurnInQuote q = quoteRescuePodTurnIn(ioSave, timeDays, stationFactionId, stationHasAuthority);
  if (!q.ok) {
    r.ok = false;
    r.reason = q.reason;
    return r;
  }

  // Payouts.
  if (stationHasAuthority) {
    ioSave.credits += q.creditsTotal;
    for (const auto& b : q.byFaction) {
      if (b.factionId == 0) continue;
      if (b.rep == 0.0) continue;
      addReputation(ioSave, b.factionId, b.rep);
    }
  }

  r.ok = true;
  r.reason = q.reason;
  r.turnedIn = q.totalPods;
  r.creditsPaid = stationHasAuthority ? q.creditsTotal : 0.0;
  r.repAwarded = stationHasAuthority ? q.repTotal : 0.0;

  ioSave.rescuedPods.clear();
  return r;
}

} // namespace stellar::sim
