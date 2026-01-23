#include "stellar/sim/StationSecurity.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

const char* stationOffenseKindName(StationOffenseKind k) {
  switch (k) {
    case StationOffenseKind::None:            return "None";
    case StationOffenseKind::WeaponDischarge: return "WeaponDischarge";
    case StationOffenseKind::Speeding:        return "Speeding";
    case StationOffenseKind::Trespass:        return "Trespass";
    default:                                 return "?";
  }
}

const char* stationOffenseActionName(StationOffenseAction a) {
  switch (a) {
    case StationOffenseAction::None:    return "None";
    case StationOffenseAction::Warning: return "Warning";
    case StationOffenseAction::Fine:    return "Fine";
    case StationOffenseAction::Bounty:  return "Bounty";
    default:                            return "?";
  }
}

double stationNoFireZoneRadiusKm(const Station& st, const StationSecurityParams& params) {
  return std::max(0.0, st.radiusKm) * std::max(0.0, params.noFireZoneRadiusFactor);
}

double stationSpeedZoneRadiusKm(const Station& st, const StationSecurityParams& params) {
  return std::max(0.0, st.radiusKm) * std::max(0.0, params.speedZoneRadiusFactor);
}

static double smoothstep01(double t) {
  t = std::clamp(t, 0.0, 1.0);
  return t * t * (3.0 - 2.0 * t);
}

double stationSpeedLimitKmS(const Station& st, const StationSecurityParams& params, double distanceKm) {
  const double base = std::max(1e-6, st.maxApproachSpeedKmS);
  const double rInner = std::max(0.0, st.radiusKm) * std::max(0.0, params.speedZoneInnerFactor);
  const double rOuter = stationSpeedZoneRadiusKm(st, params);

  // If the zone is degenerate, fall back to a simple multiplier.
  if (rOuter <= rInner + 1e-9) {
    const double mult = std::max(params.speedLimitNearMult, params.speedLimitFarMult);
    return base * std::max(1.0, mult);
  }

  const double t = (distanceKm - rInner) / (rOuter - rInner);
  const double s = smoothstep01(t);
  const double mult = params.speedLimitNearMult + (params.speedLimitFarMult - params.speedLimitNearMult) * s;
  return base * std::max(1.0, mult);
}

bool insideStationSlotTunnel(const Station& st, const math::Vec3d& relLocalKm) {
  // Mirror the geometry used in Docking.cpp for the hull cutout.
  const double wz = st.radiusKm * 1.10;

  const double slotHalfW = st.slotWidthKm * 0.5;
  const double slotHalfH = st.slotHeightKm * 0.5;

  const double zEntrance = wz;
  const double zMin = zEntrance - st.slotDepthKm;

  const bool insideTunnel =
    (std::abs(relLocalKm.x) < slotHalfW) &&
    (std::abs(relLocalKm.y) < slotHalfH) &&
    (relLocalKm.z <= zEntrance) &&
    (relLocalKm.z >= zMin);

  return insideTunnel;
}

static void resetStrikesIfStale(StationOffenseTracker& t, double timeDays, double resetDays) {
  if (t.strikes <= 0) return;
  if ((timeDays - t.lastStrikeDays) > resetDays) {
    t.strikes = 0;
    t.cooldownUntilDays = std::min(t.cooldownUntilDays, timeDays);
  }
}

static void bumpStrike(StationOffenseTracker& t, double timeDays) {
  t.strikes = std::min(9, t.strikes + 1);
  t.lastStrikeDays = timeDays;
}

static double fineMultiplier(const LawProfile& law) {
  const double strict = std::clamp(law.scanStrictness, 0.50, 2.00);
  const double antiCorrupt = std::clamp(1.08 - 0.38 * std::clamp(law.corruption, 0.0, 1.0), 0.65, 1.15);
  return std::clamp(0.65 + 0.55 * strict, 0.60, 1.75) * antiCorrupt;
}

static double repStrictness(const LawProfile& law) {
  return std::clamp(law.scanStrictness, 0.50, 2.00);
}

static double computeFineCr(StationOffenseKind kind,
                            const StationSecurityParams& params,
                            const LawProfile& law,
                            double severity01) {
  const double sev = std::clamp(severity01, 0.0, 1.0);
  const double m = fineMultiplier(law);
  const double fb = std::max(0.0, law.fineBaseCr);

  if (kind == StationOffenseKind::WeaponDischarge) {
    const double base = params.baseFineWeaponCr + 1.05 * fb;
    return base * m * (1.0 + 0.65 * sev);
  }
  if (kind == StationOffenseKind::Trespass) {
    const double base = params.baseFineTrespassCr + 0.95 * fb;
    return base * m * (1.0 + 0.55 * sev);
  }
  // Speeding
  const double base = params.baseFineSpeedingCr + 0.70 * fb;
  return base * m * (1.0 + 1.75 * sev);
}

static double computeBountyCr(StationOffenseKind kind,
                              const StationSecurityParams& params,
                              const LawProfile& law,
                              double severity01) {
  const double sev = std::clamp(severity01, 0.0, 1.0);
  const double strict = repStrictness(law);
  const double antiCorrupt = std::clamp(1.05 - 0.30 * std::clamp(law.corruption, 0.0, 1.0), 0.70, 1.10);
  const double fb = std::max(0.0, law.fineBaseCr);

  if (kind == StationOffenseKind::WeaponDischarge) {
    const double base = params.baseBountyWeaponCr + 1.25 * fb;
    return base * (0.90 + 0.35 * strict) * (1.0 + 0.35 * sev) * antiCorrupt;
  }

  // Trespass
  const double base = params.baseBountyTrespassCr + 1.05 * fb;
  return base * (0.88 + 0.30 * strict) * (1.0 + 0.25 * sev) * antiCorrupt;
}

static double computeRepPenalty(StationOffenseKind kind,
                                const LawProfile& law,
                                double severity01,
                                StationOffenseAction action) {
  const double sev = std::clamp(severity01, 0.0, 1.0);
  const double strict = repStrictness(law);

  // Warnings are informational by default.
  if (action == StationOffenseAction::Warning) {
    return 0.0;
  }

  double base = 0.3;
  double mult = 1.0;

  switch (kind) {
    case StationOffenseKind::WeaponDischarge:
      base = 0.85;
      mult = 1.0 + 0.55 * sev;
      break;
    case StationOffenseKind::Trespass:
      base = 0.55;
      mult = 1.0 + 0.45 * sev;
      break;
    case StationOffenseKind::Speeding:
      base = 0.25;
      mult = 1.0 + 0.75 * sev;
      break;
    default:
      break;
  }

  const double strictAdj = 1.0 + 0.35 * std::max(0.0, strict - 1.0);
  double rep = -std::clamp(base * mult * strictAdj, 0.10, 4.0);

  // Bounty escalation is a more serious mark.
  if (action == StationOffenseAction::Bounty) {
    rep *= 1.45;
  }
  return rep;
}

static StationOffenseAction actionFromStrikes(StationOffenseKind kind,
                                              int strikes,
                                              const StationSecurityParams& params,
                                              double severity01) {
  const bool canBounty = (kind == StationOffenseKind::WeaponDischarge) || (kind == StationOffenseKind::Trespass);

  if (strikes < params.strikesToFine) return StationOffenseAction::Warning;

  if (canBounty && strikes >= params.strikesToBounty) {
    return StationOffenseAction::Bounty;
  }

  // Speeding never becomes a bounty by default, but very reckless speeding can.
  if (!canBounty && strikes >= params.strikesToBounty && severity01 > 0.92) {
    return StationOffenseAction::Bounty;
  }

  return StationOffenseAction::Fine;
}

std::optional<StationOffenseEvent> updateStationSecurity(StationSecurityState& ioState,
                                                        const StationSecurityParams& params,
                                                        const Station& st,
                                                        const LawProfile& law,
                                                        const math::Vec3d& stationPosKm,
                                                        const math::Vec3d& stationVelKmS,
                                                        const math::Quatd& stationOrient,
                                                        const math::Vec3d& shipPosKm,
                                                        const math::Vec3d& shipVelKmS,
                                                        bool hasClearance,
                                                        bool shipDocked,
                                                        bool shipFiredWeaponThisTick,
                                                        double timeDays,
                                                        double dtSimSec) {
  // Housekeeping: reset strikes after a sufficiently long "cool off" period.
  const double resetDays = std::max(0.0, params.strikeResetSec) / 86400.0;
  resetStrikesIfStale(ioState.weapons, timeDays, resetDays);
  resetStrikesIfStale(ioState.speeding, timeDays, resetDays);
  resetStrikesIfStale(ioState.trespass, timeDays, resetDays);

  if (shipDocked) {
    ioState.speedingAccumSec = 0.0;
    ioState.trespassAccumSec = 0.0;
    return std::nullopt;
  }

  const math::Vec3d rel = shipPosKm - stationPosKm;
  const double distKm = rel.length();

  const double noFireKm = stationNoFireZoneRadiusKm(st, params);
  const double speedKm = stationSpeedZoneRadiusKm(st, params);

  const bool inNoFireZone = distKm <= noFireKm;
  const bool inSpeedZone = distKm <= speedKm;

  // Decay accumulators when out of their relevant zones.
  const double decayFast = std::max(0.0, dtSimSec) * 2.0;
  const double decaySlow = std::max(0.0, dtSimSec) * 1.25;

  if (!inSpeedZone) {
    ioState.speedingAccumSec = std::max(0.0, ioState.speedingAccumSec - decayFast);
  }
  if (hasClearance || !inNoFireZone) {
    ioState.trespassAccumSec = std::max(0.0, ioState.trespassAccumSec - decayFast);
  }

  // Priority #1: weapons discharge in no-fire zone.
  if (shipFiredWeaponThisTick && inNoFireZone) {
    auto& tr = ioState.weapons;
    if (timeDays >= tr.cooldownUntilDays) {
      bumpStrike(tr, timeDays);
      tr.cooldownUntilDays = timeDays + std::max(0.0, params.weaponEventCooldownSec) / 86400.0;

      StationOffenseEvent ev{};
      ev.kind = StationOffenseKind::WeaponDischarge;
      ev.severity01 = 1.0;
      ev.strikes = tr.strikes;
      ev.action = actionFromStrikes(ev.kind, ev.strikes, params, ev.severity01);
      if (ev.action == StationOffenseAction::Fine) {
        ev.fineCr = computeFineCr(ev.kind, params, law, ev.severity01);
        ev.repPenalty = computeRepPenalty(ev.kind, law, ev.severity01, ev.action);
      } else if (ev.action == StationOffenseAction::Bounty) {
        ev.bountyCr = computeBountyCr(ev.kind, params, law, ev.severity01);
        ev.repPenalty = computeRepPenalty(ev.kind, law, ev.severity01, ev.action);
      }
      return ev;
    }
  }

  // Priority #2: trespass (slot tunnel entry without clearance).
  bool trespassNow = false;
  if (!hasClearance && inNoFireZone) {
    const math::Vec3d relLocal = stationOrient.conjugate().rotate(shipPosKm - stationPosKm);
    trespassNow = insideStationSlotTunnel(st, relLocal);
  }

  if (trespassNow) {
    ioState.trespassAccumSec += std::max(0.0, dtSimSec);
  } else {
    ioState.trespassAccumSec = std::max(0.0, ioState.trespassAccumSec - decaySlow);
  }

  if (trespassNow && ioState.trespassAccumSec >= std::max(0.0, params.trespassTriggerSec)) {
    auto& tr = ioState.trespass;
    if (timeDays >= tr.cooldownUntilDays) {
      bumpStrike(tr, timeDays);
      tr.cooldownUntilDays = timeDays + std::max(0.0, params.genericCooldownSec) / 86400.0;
      ioState.trespassAccumSec = 0.0;

      StationOffenseEvent ev{};
      ev.kind = StationOffenseKind::Trespass;
      ev.severity01 = 0.90;
      ev.strikes = tr.strikes;
      ev.action = actionFromStrikes(ev.kind, ev.strikes, params, ev.severity01);
      if (ev.action == StationOffenseAction::Fine) {
        ev.fineCr = computeFineCr(ev.kind, params, law, ev.severity01);
        ev.repPenalty = computeRepPenalty(ev.kind, law, ev.severity01, ev.action);
      } else if (ev.action == StationOffenseAction::Bounty) {
        ev.bountyCr = computeBountyCr(ev.kind, params, law, ev.severity01);
        ev.repPenalty = computeRepPenalty(ev.kind, law, ev.severity01, ev.action);
      }
      return ev;
    }
  }

  // Priority #3: speeding (soft envelope).
  if (inSpeedZone) {
    const double speedLimit = stationSpeedLimitKmS(st, params, distKm);
    const double relSpeed = (shipVelKmS - stationVelKmS).length();

    const double tol = std::max(0.0, params.speedToleranceFrac);
    const bool overspeed = relSpeed > speedLimit * (1.0 + tol);

    double severity = 0.0;
    if (speedLimit > 1e-9) {
      const double excess = (relSpeed / speedLimit) - 1.0;
      severity = std::clamp((excess - tol) / std::max(1e-9, (1.0 - tol)), 0.0, 1.0);
    }

    if (overspeed) {
      ioState.speedingAccumSec += std::max(0.0, dtSimSec);
    } else {
      ioState.speedingAccumSec = std::max(0.0, ioState.speedingAccumSec - decaySlow);
    }

    if (overspeed && ioState.speedingAccumSec >= std::max(0.0, params.speedingTriggerSec)) {
      auto& tr = ioState.speeding;
      if (timeDays >= tr.cooldownUntilDays) {
        bumpStrike(tr, timeDays);
        tr.cooldownUntilDays = timeDays + std::max(0.0, params.genericCooldownSec) / 86400.0;
        ioState.speedingAccumSec = 0.0;

        StationOffenseEvent ev{};
        ev.kind = StationOffenseKind::Speeding;
        ev.severity01 = severity;
        ev.strikes = tr.strikes;
        ev.speedLimitKmS = speedLimit;
        ev.measuredSpeedKmS = relSpeed;
        ev.action = actionFromStrikes(ev.kind, ev.strikes, params, ev.severity01);

        if (ev.action == StationOffenseAction::Fine) {
          ev.fineCr = computeFineCr(ev.kind, params, law, ev.severity01);
          ev.repPenalty = computeRepPenalty(ev.kind, law, ev.severity01, ev.action);
        } else if (ev.action == StationOffenseAction::Bounty) {
          // Extremely reckless speeding can become a bounty.
          ev.bountyCr = computeBountyCr(StationOffenseKind::Trespass, params, law, ev.severity01);
          ev.repPenalty = computeRepPenalty(ev.kind, law, ev.severity01, ev.action);
        }

        return ev;
      }
    }
  }

  return std::nullopt;
}

} // namespace stellar::sim
