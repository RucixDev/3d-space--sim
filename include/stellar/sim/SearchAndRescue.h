#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/SaveGame.h"

#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Search & Rescue (SAR)
// -----------------------------------------------------------------------------
//
// Escape pods recovered in space are stored in SaveGame::rescuedPods.
// Each pod consumes one passenger seat while onboard (similar to passenger
// missions). This module provides headless (renderer-independent) business
// logic for:
//   - seat accounting
//   - adding a rescued pod record when scooped
//   - turning in rescued pods at a station authority for credits + reputation
//
// The prototype intentionally treats pods as abstract records (no NPC entities).
// Gameplay/UI layers can spawn visuals and call into this module.

// Default life support window used by addRescuedPod() when the caller does not
// provide a lifeSupportEndDay. Expressed in SaveGame::timeDays units.
inline constexpr double kRescuePodDefaultLifeSupportDays = 1.25; // ~30 hours

// Baseline rewards per pod (before modifiers).
inline constexpr double kRescuePodBasePayoutCr = 1850.0;
inline constexpr double kRescuePodBaseRep = 0.12;

// Multipliers applied when the pod has expired.
inline constexpr double kRescuePodExpiredPayoutMul = 0.25;
inline constexpr double kRescuePodExpiredRepMul = 0.35;

// Multipliers applied when turning in at a station that does not match the pod's
// registry faction.
inline constexpr double kRescuePodOffFactionPayoutMul = 0.70;
inline constexpr double kRescuePodOffFactionRepMul = 0.50;

struct RescuePodSummary {
  int total{0};
  int alive{0};
  int expired{0};
  int noReward{0};

  // Soonest life support deadline among non-expired pods.
  // 0 means none.
  double soonestLifeSupportEndDay{0.0};
};

RescuePodSummary summarizeRescuedPods(const SaveGame& s, double timeDays);

// -----------------------------------------------------------------------------
// Passenger seat accounting
// -----------------------------------------------------------------------------

// Seats currently occupied by active passenger missions.
int passengerMissionSeatUsage(const SaveGame& s);

// Seats occupied by rescued pods (== rescuedPods.size()).
int rescuedPodSeatUsage(const SaveGame& s);

// Total occupied seats (passenger missions + rescued pods).
int occupiedPassengerSeatsTotal(const SaveGame& s);

// Free seats available for new passenger cargo (>= 0).
int freePassengerSeats(const SaveGame& s);

// Add a new rescued pod record to SaveGame.
//
// Enforces that there is at least one free passenger seat.
// Returns true on success.
// On failure returns false and optionally writes a reason string.
bool addRescuedPod(SaveGame& ioSave, const RescuedPod& pod, double timeDays, const char** outReason = nullptr);

// -----------------------------------------------------------------------------
// Station turn-in
// -----------------------------------------------------------------------------

struct RescuePodTurnInBreakdown {
  // Registry faction for the pods in this bucket.
  // If 0, the stationFactionId is used as the effective faction.
  core::u32 factionId{0};

  int pods{0};
  int alive{0};
  int expired{0};
  int noReward{0};

  double credits{0.0};
  double rep{0.0};
};

struct RescuePodTurnInQuote {
  bool ok{false};
  const char* reason{nullptr};

  int totalPods{0};
  int rewardablePods{0};
  int expiredPods{0};
  int noRewardPods{0};

  double creditsTotal{0.0};
  double repTotal{0.0};

  std::vector<RescuePodTurnInBreakdown> byFaction{};
};

// Quote turning in all rescued pods currently onboard.
//
// stationFactionId: faction id of the station authority.
// stationHasAuthority: if false, the quote will succeed but pay 0 (pods are simply disembarked).
RescuePodTurnInQuote quoteRescuePodTurnIn(const SaveGame& s,
                                         double timeDays,
                                         core::u32 stationFactionId,
                                         bool stationHasAuthority);

struct RescuePodTurnInResult {
  bool ok{false};
  const char* reason{nullptr};

  int turnedIn{0};
  double creditsPaid{0.0};
  double repAwarded{0.0};
};

// Apply turn-in of all rescued pods currently onboard.
// Mutates credits and reputation, then clears SaveGame::rescuedPods.
// If stationHasAuthority is false, pods are still removed but no rewards are granted.
RescuePodTurnInResult applyRescuePodTurnIn(SaveGame& ioSave,
                                          double timeDays,
                                          core::u32 stationFactionId,
                                          bool stationHasAuthority);

} // namespace stellar::sim
