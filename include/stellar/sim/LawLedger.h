#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/SaveGame.h"

#include <unordered_map>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Law / Crime Ledger (core)
// -----------------------------------------------------------------------------
//
// The SDL prototype historically tracked per-faction bounties, fines, and bounty
// vouchers directly inside apps/stellar_game/main.cpp as ad-hoc unordered_maps.
//
// This module lifts that bookkeeping into the core library so:
//  - the logic is shared across apps (stellar_game, sandbox, servers, tools)
//  - save/load round-tripping is centralized and less error-prone
//  - tests can lock down edge cases deterministically
//
// The ledger intentionally stores only the *accounting* state. Game-specific
// "response" state (police alert timers, local heat, etc.) remains in the app.

struct LawFineEntry {
  double amountCr{0.0};
  double dueDay{0.0};
};

class LawLedger {
public:
  // Clear all tracked state.
  void clear();

  // ---- Bounties -------------------------------------------------------------
  double bounty(core::u32 factionId) const;
  void addBounty(core::u32 factionId, double deltaCr);
  void clearBounty(core::u32 factionId);

  // ---- Bounty vouchers ------------------------------------------------------
  double voucher(core::u32 factionId) const;
  void addVoucher(core::u32 factionId, double deltaCr);
  void clearVoucher(core::u32 factionId);

  // Redeem vouchers for the current faction into credits.
  // Returns the amount redeemed (0 if none).
  double redeemVoucher(core::u32 factionId, double& credits);

  // ---- Fines ----------------------------------------------------------------
  double fine(core::u32 factionId) const;
  double fineDueDay(core::u32 factionId) const;

  // Add to an outstanding fine ledger for this faction (unpaid minor offenses).
  // dueDay is interpreted as the "earliest due" date across outstanding fines.
  void addFine(core::u32 factionId, double deltaCr, double dueDay);

  void clearFine(core::u32 factionId);

  // Pay down a fine (returns amount actually paid).
  double payFine(core::u32 factionId, double& credits, double desiredCr);

  struct OverdueConversion {
    core::u32 factionId{0};
    double fineCr{0.0};
    double lateFeeCr{0.0};
    double bountyAddedCr{0.0};
  };

  // Convert overdue fines into bounties (warrants).
  // Returns a list of conversions applied (for UI/telemetry).
  std::vector<OverdueConversion> processOverdueFines(double timeDays, double lateFeeRate = 0.25);

  struct InterstellarFactorsResult {
    double baseFinesCr{0.0};
    double serviceFeeCr{0.0};
    double totalCostCr{0.0};
    double paidCr{0.0};
    int clearedFactions{0};
  };

  // Interstellar Factors: pay fines from other jurisdictions at a surcharge.
  //
  // On success:
  //  - credits are decreased by (baseFines * (1 + feeRate))
  //  - all fines where factionId != homeFactionId are cleared
  InterstellarFactorsResult payAllOtherFines(core::u32 homeFactionId,
                                            double& credits,
                                            double feeRate = 0.20);

  // ---- Persistence ----------------------------------------------------------
  void importFromSave(const SaveGame& s);
  void exportToSave(SaveGame& s) const;

  // ---- Accessors ------------------------------------------------------------
  const std::unordered_map<core::u32, double>& bounties() const { return bounties_; }
  const std::unordered_map<core::u32, double>& vouchers() const { return vouchers_; }
  const std::unordered_map<core::u32, LawFineEntry>& fines() const { return fines_; }

private:
  std::unordered_map<core::u32, double> bounties_;
  std::unordered_map<core::u32, double> vouchers_;
  std::unordered_map<core::u32, LawFineEntry> fines_;
};

} // namespace stellar::sim
