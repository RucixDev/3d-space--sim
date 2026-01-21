#include "stellar/sim/FormationPlanner.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace stellar::sim {
namespace {

using CostT = long long;

static CostT safeQuantizeMeters(double meters, double qMeters) {
  if (!(meters >= 0.0)) meters = 0.0;
  if (!(qMeters > 0.0)) qMeters = 1.0;
  const double v = meters / qMeters;
  // llround is deterministic for finite inputs.
  if (!std::isfinite(v)) return std::numeric_limits<CostT>::max() / 4;
  return (CostT)std::llround(v);
}

static CostT quantizedDistanceCost(const math::Vec3d& aKm,
                                  const math::Vec3d& bKm,
                                  double qMeters) {
  const double distKm = (bKm - aKm).length();
  const double meters = std::max(0.0, distKm) * 1000.0;
  return safeQuantizeMeters(meters, qMeters);
}

// Hungarian / Kuhn–Munkres algorithm (min-cost assignment).
//
// Returns a row->col assignment for a square cost matrix of size n.
// Input matrix is 0-indexed (n x n). Output values are 0-indexed.
static std::vector<int> hungarianMinCost(const std::vector<std::vector<CostT>>& cost) {
  const int n = (int)cost.size();
  std::vector<int> rowToCol;
  if (n <= 0) return rowToCol;

  // We use the classic O(n^3) potentials implementation.
  // 1-index everything internally for clarity with the standard reference.
  std::vector<CostT> u(n + 1, 0), v(n + 1, 0);
  std::vector<int> p(n + 1, 0), way(n + 1, 0);

  for (int i = 1; i <= n; ++i) {
    p[0] = i;
    int j0 = 0;
    std::vector<CostT> minv(n + 1, std::numeric_limits<CostT>::max() / 4);
    std::vector<unsigned char> used(n + 1, 0);

    do {
      used[j0] = 1;
      const int i0 = p[j0];
      int j1 = 0;
      CostT delta = std::numeric_limits<CostT>::max() / 4;

      for (int j = 1; j <= n; ++j) {
        if (used[j]) continue;
        // Reduced cost.
        const CostT cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
        if (cur < minv[j]) {
          minv[j] = cur;
          way[j] = j0;
        }
        // Deterministic tie-break: prefer lower slot index.
        if (minv[j] < delta || (minv[j] == delta && j < j1)) {
          delta = minv[j];
          j1 = j;
        }
      }

      for (int j = 0; j <= n; ++j) {
        if (used[j]) {
          u[p[j]] += delta;
          v[j] -= delta;
        } else {
          minv[j] -= delta;
        }
      }
      j0 = j1;
    } while (p[j0] != 0);

    // Augmenting.
    do {
      const int j1 = way[j0];
      p[j0] = p[j1];
      j0 = j1;
    } while (j0 != 0);
  }

  // Extract row->col.
  rowToCol.assign(n, -1);
  for (int j = 1; j <= n; ++j) {
    if (p[j] > 0) {
      const int i = p[j] - 1;
      const int col = j - 1;
      rowToCol[i] = col;
    }
  }

  // Defensive: fill any gaps (shouldn't happen for a valid square matrix).
  for (int i = 0; i < n; ++i) {
    if (rowToCol[i] >= 0) continue;
    rowToCol[i] = 0;
  }

  return rowToCol;
}

static int lookupPriorSlot(core::u64 id, std::span<const FormationSlotHint> hints) {
  for (const auto& h : hints) {
    if (h.memberId == id) return h.priorSlot;
  }
  return -1;
}

} // namespace

FormationPlan planFormation(const math::Vec3d& leaderPosKm,
                            const math::Vec3d& leaderVelKmS,
                            const math::Vec3d& leaderForwardWorld,
                            const math::Vec3d& leaderUpWorld,
                            std::span<const FormationMember> members,
                            core::u64 seed,
                            const FormationPlanParams& params,
                            std::span<const FormationSlotHint> priorSlots) {
  FormationPlan out{};

  const int n = (int)members.size();
  if (n <= 0) {
    out.pattern = FormationPattern::Trail;
    out.basis = makeBasisFromForward(leaderForwardWorld, leaderUpWorld);
    return out;
  }

  // Build a stable basis.
  math::Vec3d fwd = leaderVelKmS;
  if (fwd.lengthSq() <= 1e-9) fwd = leaderForwardWorld;
  if (fwd.lengthSq() <= 1e-12) fwd = {0, 0, 1};
  out.basis = makeBasisFromForward(fwd, leaderUpWorld);

  out.pattern = chooseEscortFormationPattern(n);

  // Generate slot positions.
  std::vector<math::Vec3d> slotPosKm;
  slotPosKm.reserve((std::size_t)n);
  for (int slot = 0; slot < n; ++slot) {
    const math::Vec3d offL = formationOffsetLocalKm(out.pattern, slot, n, params.formation, seed);
    slotPosKm.push_back(formationTargetWorldKm(leaderPosKm, out.basis, offL));
  }

  // Sort members by id for stable behaviour.
  std::vector<FormationMember> sorted;
  sorted.reserve((std::size_t)n);
  for (const auto& m : members) sorted.push_back(m);
  std::sort(sorted.begin(), sorted.end(), [](const FormationMember& a, const FormationMember& b) {
    return a.id < b.id;
  });

  const double qMeters = (params.costQuantizeMeters > 0.0) ? params.costQuantizeMeters : 1.0;
  const CostT stickyCost = (params.stickySlotPenaltyMeters > 0.0)
    ? safeQuantizeMeters(params.stickySlotPenaltyMeters, qMeters)
    : 0;

  // Build cost matrix.
  std::vector<std::vector<CostT>> cost;
  cost.resize((std::size_t)n);
  for (int i = 0; i < n; ++i) {
    cost[i].resize((std::size_t)n);
    const int prior = lookupPriorSlot(sorted[i].id, priorSlots);

    for (int j = 0; j < n; ++j) {
      CostT c = quantizedDistanceCost(sorted[i].posKm, slotPosKm[(std::size_t)j], qMeters);
      if (stickyCost > 0 && prior >= 0 && prior != j) c += stickyCost;
      cost[i][j] = c;
    }
  }

  const std::vector<int> rowToCol = hungarianMinCost(cost);

  out.assignments.clear();
  out.assignments.reserve((std::size_t)n);
  for (int i = 0; i < n; ++i) {
    const int slot = (i >= 0 && i < (int)rowToCol.size()) ? rowToCol[(std::size_t)i] : 0;

    FormationAssignment a{};
    a.memberId = sorted[(std::size_t)i].id;
    a.slotIndex = std::clamp(slot, 0, n - 1);
    a.slotPosKm = slotPosKm[(std::size_t)a.slotIndex];
    out.assignments.push_back(a);
  }

  return out;
}

} // namespace stellar::sim
