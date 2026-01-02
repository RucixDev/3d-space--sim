#pragma once

#include "stellar/sim/SaveGame.h"

#include <cstddef>
#include <string>

namespace stellar::sim {

class Universe;
struct StarSystem;
struct Station;

struct MissionBoardParams {
  double searchRadiusLy{160.0};
  std::size_t maxCandidates{128};

  int baseOfferCount{6};
  int minOfferCount{4};
  int maxOfferCount{9};

  // Rep thresholds that add one extra offer each.
  double repThreshold1{25.0};
  double repThreshold2{50.0};

  // Mission type weights (must sum <= 1.0; remainder goes to BountyKill).
  // Keep the default sum <= 0.92 so the "remainder goes to BountyKill" rule
  // yields a meaningful baseline presence without requiring contextual capping.
  double wCourier{0.24};
  double wDelivery{0.24};
  double wMultiDelivery{0.12};
  double wEscort{0.05};
  double wSalvage{0.08};
  double wPassenger{0.10};
  double wSmuggle{0.04};
  double wBountyScan{0.05};

  // Reward scaling (positive rep only): reward *= 1 + (rep/100)*repRewardBonus
  double repRewardBonus{0.10};
};

// Forward decls (keep MissionLogic.h light).
struct SystemSecurityProfile;
struct FactionProfile;

// Contextual mission distribution weights.
//
// refreshMissionOffers() adapts mission frequencies based on the local system's
// security profile and the issuing faction's traits. This struct exposes the
// resulting weights for tools/tests/UI.
//
// NOTE: These weights follow the same semantics as MissionBoardParams:
//  - They must sum to <= 1.0
//  - The remainder is assigned to BountyKill
struct MissionTypeWeights {
  double wCourier{0.0};
  double wDelivery{0.0};
  double wMultiDelivery{0.0};
  double wEscort{0.0};
  double wSalvage{0.0};
  double wPassenger{0.0};
  double wSmuggle{0.0};
  double wBountyScan{0.0};

  double sum() const {
    return wCourier + wDelivery + wMultiDelivery + wEscort + wSalvage + wPassenger + wSmuggle + wBountyScan;
  }
};

// Compute contextual mission weights from base params.
//
// Deterministic: callers should feed deterministic inputs (SystemSecurityProfile,
// FactionProfile). refreshMissionOffers() uses this internally.
MissionTypeWeights computeMissionTypeWeights(const MissionBoardParams& params,
                                             const SystemSecurityProfile& local,
                                             const FactionProfile& issuer,
                                             double playerRep = 0.0);

struct MissionTickResult {
  int failed{0};
};

struct MissionDockResult {
  int completed{0};
  int progressedMultiLeg{0};
};

struct MissionEventResult {
  int completed{0};
  double rewardCr{0.0};
};

// Ensure SaveGame::missionOffers is populated for the given docked station and day.
// Offers are deterministic per (station id, dayStamp, faction id) and persisted in the save.
void refreshMissionOffers(Universe& universe,
                          const StarSystem& currentSystem,
                          const Station& dockedStation,
                          double timeDays,
                          double playerRep,
                          SaveGame& ioSave,
                          const MissionBoardParams& params = {});

// Accept an offer from SaveGame::missionOffers.
// On success:
//  - consumes station inventory / player credits as required
//  - adds cargo to player if needed
//  - moves offer into SaveGame::missions and removes it from missionOffers
bool acceptMissionOffer(Universe& universe,
                        const Station& dockedStation,
                        double timeDays,
                        double playerRep,
                        SaveGame& ioSave,
                        std::size_t offerIndex,
                        std::string* outError = nullptr);

// Mark missions as failed when deadlines pass (applies a small rep hit).
MissionTickResult tickMissionDeadlines(SaveGame& ioSave,
                                       double timeDays,
                                       double repPenaltyOnFail = -4.0);

// When docked, attempt to complete Courier/Delivery/MultiDelivery missions.
// Applies rewards, rep gains, and returns delivered cargo to the destination market.
MissionDockResult tryCompleteMissionsAtDock(Universe& universe,
                                            const StarSystem& currentSystem,
                                            const Station& dockedStation,
                                            double timeDays,
                                            SaveGame& ioSave,
                                            double repRewardOnComplete = +2.0);

// When you scan a bounty target in-system (prototype scanner), complete matching BountyScan missions.
MissionEventResult tryCompleteBountyScan(SaveGame& ioSave,
                                         SystemId currentSystemId,
                                         core::u64 targetNpcId,
                                         double repRewardOnComplete = +2.0);

// When you destroy a bounty target in-system, complete matching BountyKill missions.
MissionEventResult tryCompleteBountyKill(SaveGame& ioSave,
                                         SystemId currentSystemId,
                                         core::u64 targetNpcId,
                                         double repRewardOnComplete = +2.0);

} // namespace stellar::sim
