#include "stellar/sim/LawLedger.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

namespace {

static bool isNearlyZero(double v, double eps = 1e-6) {
  return std::abs(v) <= eps;
}

static double nonNeg(double v) {
  return (v > 0.0) ? v : 0.0;
}

} // namespace

void LawLedger::clear() {
  bounties_.clear();
  vouchers_.clear();
  fines_.clear();
}

// ---- Bounties --------------------------------------------------------------

double LawLedger::bounty(core::u32 factionId) const {
  if (factionId == 0) return 0.0;
  const auto it = bounties_.find(factionId);
  if (it == bounties_.end()) return 0.0;
  return nonNeg(it->second);
}

void LawLedger::addBounty(core::u32 factionId, double deltaCr) {
  if (factionId == 0) return;
  if (std::abs(deltaCr) <= 1e-12) return;

  const double next = std::max(0.0, bounty(factionId) + deltaCr);
  if (isNearlyZero(next)) {
    bounties_.erase(factionId);
  } else {
    bounties_[factionId] = next;
  }
}

void LawLedger::clearBounty(core::u32 factionId) {
  if (factionId == 0) return;
  bounties_.erase(factionId);
}

// ---- Vouchers --------------------------------------------------------------

double LawLedger::voucher(core::u32 factionId) const {
  if (factionId == 0) return 0.0;
  const auto it = vouchers_.find(factionId);
  if (it == vouchers_.end()) return 0.0;
  return nonNeg(it->second);
}

void LawLedger::addVoucher(core::u32 factionId, double deltaCr) {
  if (factionId == 0) return;
  if (std::abs(deltaCr) <= 1e-12) return;

  const double next = std::max(0.0, voucher(factionId) + deltaCr);
  if (isNearlyZero(next)) {
    vouchers_.erase(factionId);
  } else {
    vouchers_[factionId] = next;
  }
}

void LawLedger::clearVoucher(core::u32 factionId) {
  if (factionId == 0) return;
  vouchers_.erase(factionId);
}

double LawLedger::redeemVoucher(core::u32 factionId, double& credits) {
  if (factionId == 0) return 0.0;
  const auto it = vouchers_.find(factionId);
  if (it == vouchers_.end()) return 0.0;

  const double amount = std::max(0.0, it->second);
  vouchers_.erase(it);

  if (amount > 0.0) credits += amount;
  return amount;
}

// ---- Fines -----------------------------------------------------------------

double LawLedger::fine(core::u32 factionId) const {
  if (factionId == 0) return 0.0;
  const auto it = fines_.find(factionId);
  if (it == fines_.end()) return 0.0;
  return nonNeg(it->second.amountCr);
}

double LawLedger::fineDueDay(core::u32 factionId) const {
  if (factionId == 0) return 0.0;
  const auto it = fines_.find(factionId);
  if (it == fines_.end()) return 0.0;
  return it->second.dueDay;
}

void LawLedger::addFine(core::u32 factionId, double deltaCr, double dueDay) {
  if (factionId == 0) return;
  if (!(deltaCr > 0.0)) return;

  LawFineEntry& e = fines_[factionId];
  e.amountCr = std::max(0.0, e.amountCr + deltaCr);

  // Earliest due date wins.
  if (dueDay > 0.0) {
    if (!(e.dueDay > 0.0)) {
      e.dueDay = dueDay;
    } else {
      e.dueDay = std::min(e.dueDay, dueDay);
    }
  }

  if (isNearlyZero(e.amountCr)) {
    fines_.erase(factionId);
  }
}

void LawLedger::clearFine(core::u32 factionId) {
  if (factionId == 0) return;
  fines_.erase(factionId);
}

double LawLedger::payFine(core::u32 factionId, double& credits, double desiredCr) {
  if (factionId == 0) return 0.0;
  if (!(desiredCr > 0.0)) return 0.0;
  if (!(credits > 0.0)) return 0.0;

  auto it = fines_.find(factionId);
  if (it == fines_.end()) return 0.0;

  const double owed = std::max(0.0, it->second.amountCr);
  if (!(owed > 0.0)) {
    fines_.erase(it);
    return 0.0;
  }

  const double pay = std::min(desiredCr, std::min(credits, owed));
  if (pay <= 0.0) return 0.0;

  credits = std::max(0.0, credits - pay);
  it->second.amountCr = std::max(0.0, it->second.amountCr - pay);

  if (isNearlyZero(it->second.amountCr)) {
    fines_.erase(it);
  }

  return pay;
}

std::vector<LawLedger::OverdueConversion> LawLedger::processOverdueFines(double timeDays, double lateFeeRate) {
  std::vector<OverdueConversion> out;

  const double rate = std::max(0.0, lateFeeRate);
  std::vector<core::u32> overdueIds;
  overdueIds.reserve(fines_.size());

  for (const auto& kv : fines_) {
    const core::u32 fid = kv.first;
    const LawFineEntry& e = kv.second;
    if (fid == 0) continue;
    if (!(e.amountCr > 1e-6)) continue;
    if (!(e.dueDay > 0.0)) continue;
    if (timeDays > e.dueDay) overdueIds.push_back(fid);
  }

  for (core::u32 fid : overdueIds) {
    auto it = fines_.find(fid);
    if (it == fines_.end()) continue;

    const double fineCr = std::max(0.0, it->second.amountCr);
    const double lateFee = fineCr * rate;
    const double bountyAdd = fineCr + lateFee;

    fines_.erase(it);
    addBounty(fid, bountyAdd);

    OverdueConversion c{};
    c.factionId = fid;
    c.fineCr = fineCr;
    c.lateFeeCr = lateFee;
    c.bountyAddedCr = bountyAdd;
    out.push_back(c);
  }

  return out;
}

LawLedger::InterstellarFactorsResult LawLedger::payAllOtherFines(core::u32 homeFactionId,
                                                                double& credits,
                                                                double feeRate) {
  InterstellarFactorsResult r{};

  const double fee = std::max(0.0, feeRate);
  std::vector<core::u32> clearIds;
  clearIds.reserve(fines_.size());

  for (const auto& kv : fines_) {
    const core::u32 fid = kv.first;
    const LawFineEntry& e = kv.second;
    if (fid == 0) continue;
    if (fid == homeFactionId) continue;
    if (!(e.amountCr > 1e-6)) continue;

    r.baseFinesCr += e.amountCr;
    clearIds.push_back(fid);
  }

  if (!(r.baseFinesCr > 0.0)) return r;

  r.serviceFeeCr = r.baseFinesCr * fee;
  r.totalCostCr = r.baseFinesCr + r.serviceFeeCr;

  if (credits + 1e-6 < r.totalCostCr) {
    return r;
  }

  credits -= r.totalCostCr;
  r.paidCr = r.totalCostCr;

  for (core::u32 fid : clearIds) {
    fines_.erase(fid);
    r.clearedFactions++;
  }

  return r;
}

// ---- Persistence -----------------------------------------------------------

void LawLedger::importFromSave(const SaveGame& s) {
  clear();

  for (const auto& b : s.bounties) {
    if (b.factionId == 0) continue;
    if (b.bountyCr <= 0.0) continue;
    addBounty(b.factionId, b.bountyCr);
  }

  for (const auto& v : s.bountyVouchers) {
    if (v.factionId == 0) continue;
    if (v.bountyCr <= 0.0) continue;
    addVoucher(v.factionId, v.bountyCr);
  }

  for (const auto& f : s.fines) {
    if (f.factionId == 0) continue;
    if (f.fineCr <= 0.0) continue;
    addFine(f.factionId, f.fineCr, f.dueDay);
  }
}

void LawLedger::exportToSave(SaveGame& s) const {
  s.bounties.clear();
  s.fines.clear();
  s.bountyVouchers.clear();

  s.bounties.reserve(bounties_.size());
  for (const auto& kv : bounties_) {
    const core::u32 fid = kv.first;
    const double amount = std::max(0.0, kv.second);
    if (fid == 0 || amount <= 1e-6) continue;
    s.bounties.push_back(FactionBounty{fid, amount});
  }
  std::sort(s.bounties.begin(), s.bounties.end(), [](const FactionBounty& a, const FactionBounty& b) {
    return a.factionId < b.factionId;
  });

  s.bountyVouchers.reserve(vouchers_.size());
  for (const auto& kv : vouchers_) {
    const core::u32 fid = kv.first;
    const double amount = std::max(0.0, kv.second);
    if (fid == 0 || amount <= 1e-6) continue;
    s.bountyVouchers.push_back(FactionBounty{fid, amount});
  }
  std::sort(s.bountyVouchers.begin(), s.bountyVouchers.end(), [](const FactionBounty& a, const FactionBounty& b) {
    return a.factionId < b.factionId;
  });

  s.fines.reserve(fines_.size());
  for (const auto& kv : fines_) {
    const core::u32 fid = kv.first;
    const LawFineEntry& e = kv.second;
    const double amount = std::max(0.0, e.amountCr);
    if (fid == 0 || amount <= 1e-6) continue;
    s.fines.push_back(FactionFine{fid, amount, e.dueDay});
  }
  std::sort(s.fines.begin(), s.fines.end(), [](const FactionFine& a, const FactionFine& b) {
    return a.factionId < b.factionId;
  });
}

} // namespace stellar::sim
