#include "stellar/proc/HyperlaneRouter.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>

namespace stellar::proc {

static double clamp01(double x) {
  return std::clamp(x, 0.0, 1.0);
}

static double edgeCostLy(double distanceLy, double bandwidth01, double risk01, const HyperlaneTravelParams& p) {
  const double dist = std::max(0.0, distanceLy);

  const double risk = clamp01(risk01);
  const double riskMul = 1.0 + std::max(0.0, p.riskWeight) * risk;

  const double bw = clamp01(bandwidth01);
  const double minF = std::clamp(p.minBandwidthFactor, 0.05, 1.0);
  const double bwF = minF + (1.0 - minF) * bw;
  const double bias = clamp01(p.bandwidthBias);
  const double denom = (1.0 - bias) + bias * bwF;

  const double c = (dist * riskMul) / std::max(0.05, denom);
  if (!std::isfinite(c)) return dist;
  return std::max(0.0, c);
}

HyperlaneRouter::HyperlaneRouter(const HyperlaneNetwork& net, const std::vector<sim::SystemStub>& stubs) {
  reset(net, stubs);
}

void HyperlaneRouter::reset(const HyperlaneNetwork& net, const std::vector<sim::SystemStub>& stubs) {
  origin_ = 0;
  computed_ = false;

  ids_.clear();
  idToIndex_.clear();
  adj_.clear();

  // Collect nodes from stubs.
  ids_.reserve(stubs.size());
  for (const auto& s : stubs) {
    if (s.id != 0) ids_.push_back(s.id);
  }

  // Also include any edge endpoints not present in the stub list.
  for (const auto& e : net.edges) {
    if (e.a != 0) ids_.push_back(e.a);
    if (e.b != 0) ids_.push_back(e.b);
  }

  // Make node ordering stable regardless of input ordering.
  std::sort(ids_.begin(), ids_.end());
  ids_.erase(std::unique(ids_.begin(), ids_.end()), ids_.end());

  idToIndex_.reserve(ids_.size() * 2);
  for (std::size_t i = 0; i < ids_.size(); ++i) {
    idToIndex_[ids_[i]] = (int)i;
  }

  adj_.assign(ids_.size(), {});

  // Build adjacency.
  for (const auto& he : net.edges) {
    auto ita = idToIndex_.find(he.a);
    auto itb = idToIndex_.find(he.b);
    if (ita == idToIndex_.end() || itb == idToIndex_.end()) continue;
    const int ia = ita->second;
    const int ib = itb->second;
    if (ia == ib) continue;

    Adj a2b;
    a2b.to = ib;
    a2b.distanceLy = he.distanceLy;
    a2b.bandwidth01 = he.bandwidth01;
    a2b.risk01 = he.risk01;
    adj_[ia].push_back(a2b);

    Adj b2a = a2b;
    b2a.to = ia;
    adj_[ib].push_back(b2a);
  }

  // Deterministic neighbor iteration.
  for (auto& v : adj_) {
    std::sort(v.begin(), v.end(), [&](const Adj& a, const Adj& b) {
      if (a.to != b.to) return ids_[a.to] < ids_[b.to];
      if (a.distanceLy != b.distanceLy) return a.distanceLy < b.distanceLy;
      if (a.bandwidth01 != b.bandwidth01) return a.bandwidth01 > b.bandwidth01;
      return a.risk01 < b.risk01;
    });
  }

  bestCost_.assign(ids_.size(), std::numeric_limits<double>::infinity());
  bestDistance_.assign(ids_.size(), std::numeric_limits<double>::infinity());
  bestRiskNot_.assign(ids_.size(), 0.0);
  bestBottleneckBw_.assign(ids_.size(), 0.0);
  bestHops_.assign(ids_.size(), 0);
  bestPrev_.assign(ids_.size(), -1);
}

bool HyperlaneRouter::compute(sim::SystemId originId, const HyperlaneTravelParams& params) {
  computed_ = false;
  origin_ = originId;

  auto itO = idToIndex_.find(originId);
  if (itO == idToIndex_.end()) return false;
  const int o = itO->second;

  const std::size_t n = ids_.size();
  if (n == 0) return false;

  // Reset result arrays.
  std::fill(bestCost_.begin(), bestCost_.end(), std::numeric_limits<double>::infinity());
  std::fill(bestDistance_.begin(), bestDistance_.end(), std::numeric_limits<double>::infinity());
  std::fill(bestRiskNot_.begin(), bestRiskNot_.end(), 0.0);
  std::fill(bestBottleneckBw_.begin(), bestBottleneckBw_.end(), 0.0);
  std::fill(bestHops_.begin(), bestHops_.end(), 0);
  std::fill(bestPrev_.begin(), bestPrev_.end(), -1);

  bestCost_[o] = 0.0;
  bestDistance_[o] = 0.0;
  bestRiskNot_[o] = 1.0;
  bestBottleneckBw_[o] = 1.0;
  bestHops_[o] = 0;
  bestPrev_[o] = -1;

  struct Q {
    double cost{0.0};
    double dist{0.0};
    int hops{0};
    int idx{0};
  };

  struct Cmp {
    bool operator()(const Q& a, const Q& b) const {
      if (a.cost != b.cost) return a.cost > b.cost;
      if (a.dist != b.dist) return a.dist > b.dist;
      if (a.hops != b.hops) return a.hops > b.hops;
      return a.idx > b.idx;
    }
  };

  std::priority_queue<Q, std::vector<Q>, Cmp> pq;
  pq.push(Q{0.0, 0.0, 0, o});

  auto better = [&](double nc, double nd, int nh, double nbw, double nrNot, int nprev,
                    double bc, double bd, int bh, double bbw, double brNot, int bprev) -> bool {
    const double eps = 1e-12;
    if (nc < bc - eps) return true;
    if (nc > bc + eps) return false;
    if (nd < bd - eps) return true;
    if (nd > bd + eps) return false;
    if (nh < bh) return true;
    if (nh > bh) return false;
    // Prefer higher bottleneck bandwidth.
    if (nbw > bbw + 1e-12) return true;
    if (nbw < bbw - 1e-12) return false;
    // Prefer lower risk (higher riskNot).
    if (nrNot > brNot + 1e-12) return true;
    if (nrNot < brNot - 1e-12) return false;
    // Final tie-breaker: predecessor index.
    return nprev < bprev;
  };

  while (!pq.empty()) {
    const Q q = pq.top();
    pq.pop();

    const int u = q.idx;
    if ((std::size_t)u >= n) continue;

    // Skip stale queue entries.
    if (q.cost > bestCost_[u] + 1e-12) continue;

    for (const Adj& e : adj_[u]) {
      const int v = e.to;
      if ((std::size_t)v >= n) continue;

      const double ec = edgeCostLy(e.distanceLy, e.bandwidth01, e.risk01, params);
      const double nc = bestCost_[u] + ec;
      const double nd = bestDistance_[u] + std::max(0.0, e.distanceLy);
      const int nh = bestHops_[u] + 1;
      const double nbw = std::min(bestBottleneckBw_[u], clamp01(e.bandwidth01));
      const double nrNot = bestRiskNot_[u] * (1.0 - clamp01(e.risk01));

      if (!std::isfinite(nc) || !std::isfinite(nd)) continue;

      if (better(nc, nd, nh, nbw, nrNot, u,
                 bestCost_[v], bestDistance_[v], bestHops_[v], bestBottleneckBw_[v], bestRiskNot_[v], bestPrev_[v])) {
        bestCost_[v] = nc;
        bestDistance_[v] = nd;
        bestHops_[v] = nh;
        bestBottleneckBw_[v] = nbw;
        bestRiskNot_[v] = nrNot;
        bestPrev_[v] = u;
        pq.push(Q{nc, nd, nh, v});
      }
    }
  }

  computed_ = true;
  return true;
}

HyperlanePathMetrics HyperlaneRouter::metricsTo(sim::SystemId targetId) const {
  HyperlanePathMetrics m{};
  if (!computed_) return m;
  if (targetId == 0) return m;

  auto it = idToIndex_.find(targetId);
  if (it == idToIndex_.end()) return m;
  const int i = it->second;
  if (i < 0 || (std::size_t)i >= bestCost_.size()) return m;

  const double c = bestCost_[i];
  if (!std::isfinite(c) || c == std::numeric_limits<double>::infinity()) return m;

  m.reachable = true;
  m.costLy = std::max(0.0, c);
  m.distanceLy = std::max(0.0, bestDistance_[i]);
  m.hops = std::max(0, bestHops_[i]);
  m.bottleneckBandwidth01 = clamp01(bestBottleneckBw_[i]);
  const double riskNot = std::clamp(bestRiskNot_[i], 0.0, 1.0);
  m.risk01 = clamp01(1.0 - riskNot);
  return m;
}

bool HyperlaneRouter::buildPathTo(sim::SystemId targetId, std::vector<sim::SystemId>& outPath) const {
  outPath.clear();
  if (!computed_) return false;
  if (origin_ == 0 || targetId == 0) return false;

  auto it = idToIndex_.find(targetId);
  if (it == idToIndex_.end()) return false;
  const int t = it->second;
  if (t < 0 || (std::size_t)t >= ids_.size()) return false;

  // Must be reachable.
  const double c = bestCost_[t];
  if (!std::isfinite(c) || c == std::numeric_limits<double>::infinity()) return false;

  // Backtrack via predecessor chain.
  const std::size_t guardMax = ids_.size() + 4;
  int cur = t;
  std::size_t guard = 0;
  while (cur >= 0 && (std::size_t)cur < ids_.size() && guard++ < guardMax) {
    outPath.push_back(ids_[(std::size_t)cur]);
    const int prev = bestPrev_[(std::size_t)cur];
    if (prev < 0) break;
    cur = prev;
  }

  if (outPath.empty()) return false;
  std::reverse(outPath.begin(), outPath.end());

  // Sanity: path should start at origin.
  if (outPath.front() != origin_) {
    outPath.clear();
    return false;
  }
  return true;
}

} // namespace stellar::proc
